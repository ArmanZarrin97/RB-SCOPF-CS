# Load in the same module as Esev_cal_exact and Rcas_cal_RS, or pass callbacks.
# No FigConobj definition is needed: the result is a NamedTuple.
# All matrices use [time, contingency]. Stages are zero-based; report rows
# and optstage[2] are one-based. Report rows must represent consecutive stages.

_cr_value(x::Number) = Float64(x)
_cr_value(x) = Float64(JuMP.value(x))

function _cr_key(key, n)
    k = key isa CartesianIndex ? Tuple(key) : key
    k = k isa Tuple ? k : (k,)
    if length(k) == 1 && (only(k) isa Tuple || only(k) isa CartesianIndex)
        return _cr_key(only(k), n)
    end
    length(k) == n && all(x -> x isa Integer, k) ||
        throw(ArgumentError("Expected $n integer indices; received $key"))
    return ntuple(i -> Int(k[i]), n)
end

# Public JuMP rowtable API preserves actual axis labels, including a single
# axis of tuple keys. Never align different containers by storage position.
function _cr_entries(container, n)
    entries = if container isa AbstractDict
        [(_cr_key(k, n), _cr_value(v)) for (k, v) in pairs(container)]
    else
        [(_cr_key(Tuple(r)[1:end-1], n), Float64(last(Tuple(r))))
         for r in JuMP.Containers.rowtable(_cr_value, container)]
    end
    length(unique(first.(entries))) == length(entries) ||
        throw(ArgumentError("Duplicate normalized container indices"))
    return entries
end

function _cr_binary(v, tol)
    isfinite(v) || throw(ArgumentError("Nonfinite binary solution value"))
    abs(v) <= tol && return false
    abs(v - 1) <= tol && return true
    throw(ArgumentError("Expected a binary solution value, received $v"))
end

function _cr_nonnegative(v, label; tol=1e-8)
    x = Float64(v)
    isfinite(x) && x >= -tol || throw(ArgumentError("Invalid $label: $v"))
    return max(x, 0.0)
end

# Missing/NaN report statistics are displayed as unavailable, never overwritten.
function _cr_stat(row, name)
    hasproperty(row, Symbol(name)) || return NaN
    x = getproperty(row, Symbol(name))
    return x isa Real && isfinite(x) ? Float64(x) : NaN
end

function _cr_vstats(row, suffix)
    lo, hi = "Vmin_" * suffix, "Vmax_" * suffix
    # Original pre-contingency maxima use Vmin_Abs and Vmax_Abs.
    abslo = suffix == "N" ? "Vmin_Abs" : lo * "_Abs"
    abshi = suffix == "N" ? "Vmax_Abs" : hi * "_Abs"
    nlo, nhi = _cr_stat(row, lo * "_N"), _cr_stat(row, hi * "_N")
    n = nlo + nhi
    weighted(m, count) = count == 0 ? 0.0 : m * count
    avg = n == 0 ? 0.0 :
          (weighted(_cr_stat(row, lo * "_M"), nlo) +
           weighted(_cr_stat(row, hi * "_M"), nhi)) / n
    return (max(_cr_stat(row, abslo), _cr_stat(row, abshi)), avg, n)
end

function _cr_summary(row)
    row === nothing && return nothing
    flow(s) = (_cr_stat(row, "Smax_" * s),
               _cr_stat(row, "Smax_" * s * "_M"),
               _cr_stat(row, "Smax_" * s * "_N"))
    return [("Voltage (p.u.)", (_cr_vstats(row, "N")..., _cr_vstats(row, "C")...)),
            ("Line flow (p.u.)", (flow("N")..., flow("C")...))]
end

"""
    optstage_out, data = Cas_risk_process(Glob, casreps, iter, model=nothing; ...)

Required report columns: Sev, LL. VM is needed only to calculate comparison
inaccuracy. Each report must have at least one row (stage zero).

Preserved stopping rules:
* RSS: pbd <= rss_tol selects stage zero; otherwise the final stage.
* SS: bin_Cas == 1 selects stage zero; zero selects the final stage.
* MS: the unique bin_Cas[t,c,s] == 1 selects report row s.
* Slack fallback: slk_cas >= slack_tol selects the final stage.

Severity always comes from Sev at the selected row. Expected severity comes
from exact_risk(report) in simulation or rcas_OP in RSS optimization.
Other optimization modes have no expected-severity matrix.

Pairs absent from a stopping container retain their final simulated stage.
Absent RSS risk entries default to zero, matching the original aggregate;
set missing_risk=:error to reject incomplete risk coverage instead. The
modeled_risk mask identifies which RSS pairs actually have model values.

`aggregate_risk` receives a [time, contingency] matrix, without additional
weighting here. Your Rcas_cal_RS must perform the intended aggregation.
"""
function Cas_risk_process_new(Glob, casreps, iter, model=nothing;
        reprow=nothing, optstage=nothing, h_cas=nothing,
        exact_risk=Esev_cal_exact, aggregate_risk=Rcas_cal_RS,
        binary_tol=1e-5, rss_tol=1e-6, slack_tol=1e-4,
        missing_risk=:zero)
    nTP, nCon = Int(Glob["nTL"]), size(Glob["SysData"]["ConL"], 1)
    nTP > 0 && nCon > 0 || throw(ArgumentError("Empty time/contingency domain"))
    size(casreps) == (nTP, nCon) || throw(DimensionMismatch("casreps must be ($nTP, $nCon)"))
    0 <= binary_tol < 0.5 || throw(ArgumentError("binary_tol must be in [0, 0.5)"))
    rss_tol >= 0 && slack_tol >= 0 || throw(ArgumentError("Thresholds must be nonnegative"))
    missing_risk in (:zero, :error) || throw(ArgumentError("missing_risk must be :zero or :error"))
    checkpair(tc) = checkbounds(Bool, casreps, tc...) ? nothing :
        throw(ArgumentError("Pair $tc is outside the report matrix"))
    for report in casreps
        size(report, 1) > 0 || throw(ArgumentError("Empty cascade report"))
        all(s -> s in propertynames(report), (:Sev, :LL)) ||
            throw(ArgumentError("Cascade report requires Sev and LL columns"))
    end

    ism = model === nothing
    mode = get(Glob, "Risk_cascade", "")
    stop_stage = [size(casreps[t,c], 1) - 1 for t in 1:nTP, c in 1:nCon]
    selected = Dict{Tuple{Int,Int},Int}() # pair => one-based report row
    risk = nothing
    modeled_risk = falses(nTP, nCon)
    tr = nothing

    if !ism
        JuMP.has_values(model) || throw(ArgumentError("Model has no solution values"))
        objects = JuMP.object_dictionary(model)
        if mode == "RSS"
            for (tc, v) in _cr_entries(model[:pbd], 2)
                checkpair(tc)
                isfinite(v) || throw(ArgumentError("Nonfinite pbd at $tc"))
                selected[tc] = v <= rss_tol ? 1 : size(casreps[tc...], 1)
            end
            # risk = zeros(nTP, nCon)
            risk = [Float64(casreps[t, c][end, "Sev"]) for t in 1:nTP, c in 1:nCon]
            for (tc, v) in _cr_entries(model[:rcas_OP], 2)
                checkpair(tc)
                risk[tc...] = _cr_nonnegative(v, "rcas_OP at $tc")
                modeled_risk[tc...] = true
            end
            missing_risk == :error && !all(modeled_risk) &&
                throw(ArgumentError("rcas_OP does not cover every (time, contingency) pair"))
        elseif haskey(objects, :bin_Cas)
            if mode == "MS"
                groups = Dict{Tuple{Int,Int},Vector{Int}}()
                for ((t,c,s), v) in _cr_entries(model[:bin_Cas], 3)
                    tc = (t,c)
                    checkpair(tc)
                    1 <= s <= size(casreps[t,c], 1) ||
                        throw(ArgumentError("Invalid MS report row $s at $tc"))
                    rows = get!(groups, tc, Int[])
                    _cr_binary(v, binary_tol) && push!(rows, s)
                end
                for (tc, rows) in groups
                    length(rows) == 1 || throw(ArgumentError("MS requires exactly one selected row at $tc; got $rows"))
                    selected[tc] = only(rows)
                end
            elseif mode == "SS"
                for (tc, v) in _cr_entries(model[:bin_Cas], 2)
                    checkpair(tc)
                    selected[tc] = _cr_binary(v, binary_tol) ? 1 : size(casreps[tc...], 1)
                end
            else
                throw(ArgumentError("bin_Cas supplied with unsupported mode $mode"))
            end
        elseif haskey(objects, :slk_cas)
            for (tc, v) in _cr_entries(model[:slk_cas], 2)
                checkpair(tc)
                isfinite(v) || throw(ArgumentError("Nonfinite slk_cas at $tc"))
                selected[tc] = v >= slack_tol ? size(casreps[tc...], 1) : 1
            end
        else
            throw(ArgumentError("Expected pbd, bin_Cas, or slk_cas"))
        end
        if haskey(objects, :TR)
            # Preserve original association: TR[2] voltage; TR[1] line flow.
            tr = (voltage=_cr_value(model[:TR][2]), lineflow=_cr_value(model[:TR][1]))
        end
    else
        risk = [_cr_nonnegative(exact_risk(casreps[t,c]), "expected severity at ($t,$c)")
                for t in 1:nTP, c in 1:nCon]
        modeled_risk .= true
    end

    for (tc, row) in selected
        stop_stage[tc...] = row - 1
    end
    severity = [_cr_nonnegative(casreps[t,c][stop_stage[t,c]+1, "Sev"], "severity at ($t,$c)")
                for t in 1:nTP, c in 1:nCon]
    optstage_out = if ism
        nothing
    else
        pairs = sort!(collect(keys(selected)))
        (pairs, [selected[tc] for tc in pairs],
         Any[deepcopy(casreps[tc...][selected[tc], "LL"]) for tc in pairs])
    end

    greenbar = Dict{Tuple{Int,Int},Int}()
    redbar = Dict{Tuple{Int,Int},Int}()
    tc_new = Set{Tuple{Int,Int}}()
    inaccuracy = nothing
    if optstage !== nothing
        pairs, rows, elements = optstage
        length(pairs) == length(rows) == length(elements) ||
            throw(DimensionMismatch("optstage vectors must have equal lengths"))
        seen = Set{Tuple{Int,Int}}()
        for i in eachindex(pairs)
            tc = _cr_key(pairs[i], 2)
            checkpair(tc)
            tc in seen && throw(ArgumentError("Duplicate optstage pair $tc"))
            push!(seen, tc)
            row = Int(rows[i])
            row >= 1 || throw(ArgumentError("optstage rows are one-based"))
            report = casreps[tc...]
            if row >= size(report, 1)
                greenbar[tc] = row
            elseif issubset(report[row, "LL"], elements[i])
                redbar[tc] = row
            end
        end
        history = h_cas !== nothing && hasproperty(h_cas, :tup) ? Set(h_cas.tup) : nothing
        for t in 1:nTP, c in 1:nCon
            tc, report = (t,c), casreps[t,c]
            active = size(report, 1) > 1 || report[1, "Sev"] > 0
            (!active || haskey(greenbar, tc) || haskey(redbar, tc)) && continue
            initial = report[1, "LL"]
            # Empty initial LL is not evidence of a previously handled overload.
            if history !== nothing && !isempty(initial) && all(e -> (t,e,c) in history, initial)
                redbar[tc] = 1
            else
                push!(tc_new, tc)
            end
        end
        if isempty(redbar)
            inaccuracy = 0.0
        elseif haskey(Glob, "I_relax") && all(tc -> :VM in propertynames(casreps[tc...]), keys(redbar))
            inaccuracy = 0.0
            for (tc, row) in redbar
                vm = casreps[tc...][row, "VM"]
                peak = isempty(vm) ? 0.0 : maximum(vm)
                isfinite(peak) || throw(ArgumentError("Nonfinite VM at $tc"))
                inaccuracy = max(inaccuracy, peak - Glob["I_relax"])
            end
        end
    end
    #
    # sysrisk = risk === nothing ? nothing : Float64(aggregate_risk(copy(risk)))
    # maxrisk = risk === nothing ? nothing : maximum(risk)
    aggregate_values = if risk === nothing nothing
    elseif !ism && mode == "RSS" ifelse.(modeled_risk, risk, 0.0) else copy(risk) end
    sysrisk = aggregate_values === nothing ? nothing : Float64(aggregate_risk(aggregate_values))
    maxrisk = aggregate_values === nothing ? nothing : maximum(aggregate_values)
    #
    title = (ism ? "Simulation" : "Optimization") * " — iteration $iter"
    data = (; title, severity, risk, stop_stage, modeled_risk,
             Sys_Risk=sysrisk, Max_OP_Risk=maxrisk, ism, optstage=deepcopy(optstage),
             greenbar, redbar, tc_new, inaccuracy, tr, summary=_cr_summary(reprow))
    return optstage_out, data
end

_cr_fmt(x; count=false) = x === nothing || !isfinite(x) ? "—" :
    count ? string(round(Int, x)) : @sprintf("%.3f", x)

"""
    fig = plot_risk(Glob, data; metric=:auto, annotations=:auto, ...)

metric=:auto plots expected severity if available, otherwise severity.
Use metric=:severity or :risk explicitly for consistent figure comparisons.
annotations=:auto restores the original system-dependent annotations:
pass test_sys="39bus/" for detailed cells, otherwise values only.
Use :compact, :detailed, or :none to override this choice.
Default figure sizing, typography, footer geometry, and color scale follow
 the original plot. units labels the existing numbers; it does not convert them.
Units are labels only: no implicit MW conversion or base-power scaling occurs.
Set colorrange=(0, upper) to share an identical scale across figures.
"""
function plot_risk_new(Glob, data; metric=:auto, annotations=:auto,
        test_sys=get(Glob, "Test_Sys", ""), units="p.u.", colorrange=nothing, figure_size=nothing,
        fontsize=nothing, cell_fontsize=nothing, zero_tol=1e-10)
    metric in (:auto, :severity, :risk) || throw(ArgumentError("Invalid metric"))
    annotations in (:auto, :compact, :detailed, :none) || throw(ArgumentError("Invalid annotation mode"))
    zero_tol >= 0 || throw(ArgumentError("zero_tol must be nonnegative"))
    chosen = metric == :auto ? (data.risk === nothing ? :severity : :risk) : metric
    values = chosen == :risk ? data.risk : data.severity
    values === nothing && throw(ArgumentError("Expected severity unavailable for this optimization mode"))
    nTP, nCon = size(values)
    size(data.stop_stage) == size(values) || throw(DimensionMismatch("Stage/value matrices differ"))
    (nTP, nCon) == (Int(Glob["nTL"]), size(Glob["SysData"]["ConL"], 1)) ||
        throw(DimensionMismatch("Glob and processed data differ"))
    # A zero expected value alone does not establish that no event occurred.
    nc = (abs.(data.severity) .<= zero_tol) .& (data.stop_stage .== 0) .& (abs.(values) .<= zero_tol)
    z = copy(values)
    z[nc] .= NaN
    # Original system-capacity color scale (explicit colorrange still overrides).
    sys = Glob["SysData"]
    ls_max = sum(sys["Load"][!, "Pd (pu)"]; init=0.0) +
             sum(sys["Gen"][!, "P_max (pu)"]; init=0.0) +
             sum(sys["RES"][!, "Cap"]; init=0.0)
    limits = colorrange === nothing ? (0.0, max(Float64(ls_max), eps(Float64))) : colorrange
    length(limits) == 2 && all(isfinite, limits) && limits[1] < limits[2] ||
        throw(ArgumentError("colorrange must contain two finite increasing bounds"))

    # Original adaptive sizing, including its compact bottom row.
    cell_w = min(2000, 2000/nTP)
    cell_h = min(cell_w, 600/nCon)
    fig_w = min(2200, nTP*cell_w + 200)
    fig_h = min(750, nCon*cell_h + 150)
    footer_h = data.summary === nothing ? 60 : 80
    dims = figure_size === nothing ? (round(Int, fig_w), round(Int, fig_h + footer_h)) : figure_size
    fig_w = dims[1]
    cell_size = min(cell_w, cell_h)
    fsize_cell = cell_fontsize === nothing ? clamp(round(Int, cell_size/9), 18, 28) : cell_fontsize
    fsize_ax = fontsize === nothing ? clamp(round(Int, cell_size/8), 18, 28) : fontsize
    font_table = fsize_ax
    off2 = clamp(cell_h/300, 0.10, 0.20)
    linespace = 0.22
    detailed = annotations == :detailed || (annotations == :auto && test_sys == "39bus/")
    iso = data.optstage === nothing
    label = chosen == :risk ? "Expected severity" : "Sev"
    fig = Figure(size=dims, fontsize=fsize_ax)
    metric_title = chosen == :risk ?
    "Expected severity (ESEV)" :
    "Worst severity (SEV_Max)"

    ax = Axis(fig[1,1]; title="$(data.title) — $metric_title",
              xlabel="Time Period", ylabel="Contingency",
              xticks=(1:nTP, string.(1:nTP)), yticks=(1:nCon, string.(1:nCon)),
              xticklabelrotation=pi/4)
    hm = CairoMakie.heatmap!(ax, 1:nTP, 1:nCon, z;
        colormap=[:mistyrose, :red, :darkred], colorrange=limits, nan_color=:lightgrey)
    for t in 1:nTP, c in 1:nCon
        tc, val, stage = (t,c), values[t,c], data.stop_stage[t,c]
        color = val >= (limits[1]+limits[2])/2 ? :white : :black
        annotations == :none && continue
        if nc[t,c]
            text!(ax, t, c; text="NC", align=(:center,:center), fontsize=fsize_cell, color=:grey)
        elseif !detailed
            text!(ax, t, c; text=@sprintf("%.2f", val), align=(:center,:center),
                fontsize=fsize_cell-2, color=color)
        elseif stage == 0
            text!(ax, t, c+off2; text=@sprintf("%.2f", val), align=(:center,:center), fontsize=fsize_cell, color=color)
            text!(ax, t, c-off2; text=units, align=(:center,:center), fontsize=fsize_cell, color=color)
        else
            text!(ax, t, c+0.4; text=@sprintf("%.2f", val), align=(:center,:center), fontsize=fsize_cell, color=color)
            text!(ax, t, c+0.4-linespace; text=units, align=(:center,:center), fontsize=fsize_cell, color=color)
            sublabel = get(Glob, "Risk_cascade", "") == "RSS" && data.risk !== nothing ?
                @sprintf("r: %.1e", data.risk[t,c]) : "s: $(stage),"
            text!(ax, t, c-0.35+linespace; text=sublabel,
                align=(:center,:center), fontsize=fsize_cell, color=color)
            if tc in data.tc_new
                text!(ax, t, c-0.4; text="*NEW*", align=(:center,:center), fontsize=fsize_cell-3, color=color)
            elseif haskey(data.redbar, tc) || haskey(data.greenbar, tc)
                barcolor = haskey(data.redbar, tc) ? :black : :green
                y = c-0.375+linespace/2
                lines!(ax, [t-0.3,t+0.3], [y,y]; color=barcolor, linewidth=1.5)
            end
        end
    end
    Colorbar(fig[1,2], hm; label="$label ($units)", ticks=LinearTicks(5), flipaxis=true)
    bottom = GridLayout(fig[2,1:2])
    if data.summary !== nothing
        ax_t = Axis(bottom[1,1])
        hidedecorations!(ax_t)
        hidespines!(ax_t)
            # Grid definition
            if data.ism && iso
                cw = [0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
            else
                cw = [0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
            end
            # cw = [0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
            rh = [0.10, 0.10, 0.10, 0.10]
            # CairoMakie.xlims!(ax_t, 0, sum(cw))
            # CairoMakie.ylims!(ax_t, 0, sum(rh))

            x0_2=sum(cw)+0.02; x1_2=sum(cw)+0.06; y0_2=0; y1_2=sum(rh)
            CairoMakie.xlims!(ax_t, 0, x1_2)
            CairoMakie.ylims!(ax_t, y0_2, y1_2)

            cx = cumsum([0; cw])
            ry = cumsum([0; rh])
            flip(y) = sum(rh) - y

            cell_center_x(col) = (cx[col+1] + cx[col+2]) / 2
            cell_center_y(row) = flip((ry[row+1] + ry[row+2]) / 2)

            # Header backgrounds
            for (span, clr) in [((0, 1:3), (:green, 0.15)), ((0, 4:6), (:darkorange, 0.15))]
                row, cols = span
                x0, x1 = cx[first(cols)+1], cx[last(cols)+2]
                y0, y1 = flip(ry[row+2]), flip(ry[row+1])
                poly!(ax_t, Point2f[(x0,y0),(x1,y0),(x1,y1),(x0,y1)];
                    color = clr, strokecolor = :lightgrey, strokewidth = 0.5)
            end

            poly!(ax_t, Point2f[(x0_2,y0_2),(x1_2,y0_2),(x1_2,y1_2),(x0_2,y1_2)];
                    color = (:grey, 0.15), strokecolor = :lightgrey, strokewidth = 0.5)

            # poly!(ax_t, Point2f[(x0,y0),(x1,y0),(x1,y1),(x0,y1)];
            #         color = clr, strokecolor = :lightgrey, strokewidth = 0.5)

            if !data.ism || !iso
                row, cols = 0 , 7:7
                clr = (:steelblue, 0.15)
                x0, x1 = cx[first(cols)+1], cx[last(cols)+2]
                y0, y1 = flip(ry[row+3]), flip(ry[row+1])
                poly!(ax_t, Point2f[(x0,y0),(x1,y0),(x1,y1),(x0,y1)];
                    color = clr, strokecolor = :lightgrey, strokewidth = 0.5)
            end


            # Grid lines
            for (i, x) in enumerate(cx)
                i > 1 && lines!(ax_t, [x, x], [0, sum(rh)-rh[end]]; color = :lightgrey, linewidth = 0.6)
            end
            lines!(ax_t, [cx[5], cx[5]], [0, sum(rh)]; color = :grey, linewidth = 1.2)
            lines!(ax_t, [cx[8], cx[8]], [0, sum(rh)]; color = :grey, linewidth = 1.2)

            # Row 0: Pre / Post headers
            text!(ax_t, (cx[2]+cx[5])/2, cell_center_y(0);
                text="Pre-Contingency",  align=(:center,:center), fontsize = font_table, color=:green,  font=:bold)
            text!(ax_t, (cx[5]+cx[8])/2, cell_center_y(0);
                text="Post-Contingency", align=(:center,:center), fontsize = font_table, color=:darkorange, font=:bold)
            if !data.ism
                text!(ax_t, (cx[8]+cx[9])/2, cell_center_y(0);
                    text="TR", align=(:center,:center), fontsize = font_table, color=:steelblue, font=:bold)
            elseif (!iso)
                text!(ax_t, (cx[8]+cx[9])/2, (cell_center_y(0)+cell_center_y(1))/2;
                    text="inacc", align=(:center,:center), fontsize = font_table, color=:steelblue, font=:bold)
            end

            text!(ax_t, (x0_2+x1_2)/2, (cell_center_y(0));
                    text="Max OP Risk", align=(:center,:center), fontsize = font_table, color=:steelblue, font=:bold)
            text!(ax_t, (x0_2+x1_2)/2, (cell_center_y(2));
                    text="System Risk", align=(:center,:center), fontsize = font_table, color=:steelblue, font=:bold)
            # Row 1: column headers
            for (col, lbl) in enumerate(["", "Max", "Mean", "N", "Max", "Mean", "N"])
                text!(ax_t, cell_center_x(col-1), cell_center_y(1);
                    text=lbl, align=(:center,:center), fontsize = font_table, color=:black, font=:bold)
            end

            table_rows = data.summary
            for (r, (label, vals)) in enumerate(table_rows)
                row = r + 1
                text!(ax_t, cx[2] - 0.01, cell_center_y(row);
                    text=label, align=(:right,:center), fontsize=font_table, color=:black, font=:bold)
                for (col, val) in enumerate(vals)
                    text!(ax_t, cell_center_x(col), cell_center_y(row);
                        text=_cr_fmt(val; count=col in (3,6)),
                        align=(:center,:center), fontsize=font_table)
                end
            end
            if !data.ism
                text!(ax_t, cell_center_x(7), cell_center_y(1);
                    text="size", align=(:center,:center), fontsize=font_table, color=:steelblue)
                for (row, val) in ((2, data.tr === nothing ? nothing : data.tr.voltage),
                                   (3, data.tr === nothing ? nothing : data.tr.lineflow))
                    text!(ax_t, cell_center_x(7), cell_center_y(row);
                        text=_cr_fmt(val), align=(:center,:center), fontsize=font_table)
                end
            elseif !iso
                text!(ax_t, cell_center_x(7), cell_center_y(2);
                    text="-", align=(:center,:center), fontsize=font_table)
                text!(ax_t, cell_center_x(7), cell_center_y(3);
                    text=_cr_fmt(data.inaccuracy), align=(:center,:center), fontsize=font_table)
            end
            riskfmt(v) = v === nothing || !isfinite(v) ? "—" : @sprintf("%.3e", v)
            text!(ax_t, (x0_2+x1_2)/2, cell_center_y(1);
                text=riskfmt(data.Max_OP_Risk), align=(:center,:center), fontsize=font_table)
            text!(ax_t, (x0_2+x1_2)/2, cell_center_y(3);
                text=riskfmt(data.Sys_Risk), align=(:center,:center), fontsize=font_table)
            colsize!(bottom, 1, fig_w/2 * ((data.ism && iso) ? 1 : 35/31))
    end
    rowsize!(fig.layout, 2, footer_h)
    Legend(bottom[1, data.summary === nothing ? (1:2) : 2],
        [PolyElement(color=:lightgrey)], ["NC — no cascade and zero severity"],
        orientation=:vertical, tellwidth=false, tellheight=false,
        valign=:center, halign=:right)
    rowgap!(fig.layout, 1, 6)
    return fig
end

# EXAMPLES (run after loading your application and this file):
# optstage, optdata = Cas_risk_process(Glob, casreps, iter, model; reprow=reprow)
# fig = plot_risk(Glob, optdata; annotations=:detailed)
# display(fig)
# save("optimization_risk.png", fig)
#
# _, simdata = Cas_risk_process(Glob, new_casreps, iter;
#     optstage=optstage, h_cas=h_cas, reprow=new_reprow)
# display(plot_risk(Glob, simdata; metric=:risk))
#
# NOTE: this replaces the two original methods. It does not use FigConobj,
# LS/LB/Risk_map, or hidden Test_Sys/model/RS/casreps globals in plotting.
# `data.severity`, `data.risk`, and `data.stop_stage` all use [t,c].