function make_cascade_stopping_stage_plot4(Glob, casreps, iter, model=nothing; reprow=nothing, stoppoint=nothing)

    nCon = size(Glob["SysData"]["ConL"], 1)
    nTP  = Glob["nTL"]

    # --- Adaptive sizing ---
    cell_w    = min(50,500/nTP)                              # pixels per time period
    cell_h    = min(20,150/nCon)                             # pixels per contingency
    fig_w     = min(700, nTP  * cell_w + 200)  # +200 for colorbar + margins
    fig_h     = min(300, nCon * cell_h + 150)  # +150 for title + bottom row
    fig_h_tot = fig_h + (reprow === nothing ? 60 : 80)

    # Font and offset scale with cell size
    cell_size   = min(cell_w, cell_h)
    fsize_cell  = clamp(round(Int, cell_size / 9),  7, 12)   # annotation font
    fsize_ax    = clamp(round(Int, cell_size / 8),  9, 14)   # axis font
    off3        = clamp(cell_h / 220, 0.15, 0.32)            # 3-line vertical offset
    off2        = clamp(cell_h / 300, 0.10, 0.20)            # 2-line vertical offset
    # line_h = fsize_cell * (nCon / (fig_h - 150))   # axis units per line of text
    linespace = fsize_cell / cell_h *0.5
    # padding = cell_h / fsize_cell

    if model === nothing
        # --- Simulation ---
        tc_stop   = findall(x -> size(x, 1) > 1, casreps)
        tc_stopat = getfield.(tc_stop, :I)
        stopat    = map(x -> x[end, 1], casreps[tc_stop])
        # casv=0
        # casvtc=0
        # for d in tc_stop
        #     stagevio=maximum.(casreps[d]."VM",init=-10)
        #     casvio=minimum(filter(x-> x>=0,stagevio))
        #     if casvio>casv
        #         casv=casvio
        #         casvtc=d
        #     end
        # end
        LNS_stop  = [casreps[CartesianIndex(tc_stopat[x])][stopat[x]+1, "LNS_T"]
                     for x in 1:length(stopat)]
        labels_map = Dict(tc_stopat[x] => stopat[x]  for x in 1:length(stopat))
        title      = "Simulation — Total Load Shed per (Time, Contingency),  iter=$(iter)"

    else
        # --- Optimization ---
        bin_cas  = model[:bin_Cas]
        indices  = getfield.(getfield.(keys(bin_cas), :I), 1)
        tc_all   = Set(map(((t, c, _),) -> (t, c), indices))
        vals     = value.(bin_cas).data
        stages   = findall(x -> x == 1, vals)
        @assert length(tc_all) == length(stages)
        stopat    = getfield.(indices[stages], 3) .- 1
        tc_stopat = map(((t, c, _),) -> (t, c), indices[stages])
        LNS_stop  = [casreps[CartesianIndex(tc_stopat[x])][stopat[x]+1, "LNS_T"]
                     for x in 1:length(stages)]
        labels_map = Dict(tc_stopat[x] => stopat[x]  for x in 1:length(stages))
        title      = "Optimization — Load Shed at Stopping Stage,  iter=$(iter)"
    end

    # --- Build matrices ---
    loadshed_map = Dict(tc_stopat[x] => LNS_stop[x] for x in 1:length(stopat))

    LS = [get(loadshed_map, (t, c), casreps[t, c][end, "LNS_T"]) for c in 1:nCon, t in 1:nTP]
    LB = [get(labels_map,   (t, c), 0)                           for c in 1:nCon, t in 1:nTP]

    LS_plot = Float64.(LS)
    LS_plot[LS .== 0.0 .&& LB .== 0] .= NaN

    ls_max = sum(Glob["SysData"]["Load"][!, "Pd (pu)"]; init=0)

    # --- Figure ---
    fig = Figure(size = (fig_w, fig_h_tot), fontsize = fsize_ax)

    ax = Axis(fig[1, 1];
        title              = title,
        xlabel             = "Time Period",
        ylabel             = "Contingency",
        xticks             = (1:nTP,  string.(1:nTP)),
        yticks             = (1:nCon, string.(1:nCon)),
        xticklabelrotation = π/4
    )

    hm = CairoMakie.heatmap!(ax, 1:nTP, 1:nCon, LS_plot';
        colormap   = cgrad([:mistyrose, :red, :darkred]),
        colorrange = (0, ls_max),
        nan_color  = :lightgrey
    )

    for t in 1:nTP, c in 1:nCon
        ls        = LS[c, t]
        lbl       = LB[c, t]
        no_event  = ls == 0.0 && lbl == 0
        fontcolor = ls >= 0.6 * ls_max ? :white : :black

        if no_event
            text!(ax, t, c;
                text = "NC", align = (:center, :center), fontsize = fsize_cell, color = :grey)
        elseif lbl == 0
            # Load shed but no cascade — two lines
            text!(ax, t, c + off2;
                text = @sprintf("%.2f", ls), align = (:center, :center), fontsize = fsize_cell, color = fontcolor)
            text!(ax, t, c - off2;
                text = "MW",                 align = (:center, :center), fontsize = fsize_cell, color = fontcolor)
        else
            # Load shed + cascade — three lines
            text!(ax, t, c + 0.4;
                text = @sprintf("%.2f", ls), align = (:center, :center), fontsize = fsize_cell, color = fontcolor)
            text!(ax, t, c+ 0.4 - linespace;
                text = "MW",                 align = (:center, :center), fontsize = fsize_cell, color = fontcolor)
            text!(ax, t, c - 0.4+linespace;
                text = "s: $(lbl)",          align = (:center, :center), fontsize = fsize_cell, color = fontcolor)
        end
    end

    #   text!(ax, 8, 5;
    #             text = @sprintf("%.2f", ls), align = (:center, :center), fontsize = 9, color = fontcolor)


    Colorbar(fig[1, 2], hm;
        label    = "Load Shed (MW)",
        ticks    = LinearTicks(5),
        flipaxis = true
    )

    # --- Bottom row: table (left) + legend (right) ---
    bottom = GridLayout(fig[2, 1:2])

    # ── Table (left) ──────────────────────────────────────────
    if reprow !== nothing
        ax_t = Axis(bottom[1, 1])
        hidedecorations!(ax_t)
        hidespines!(ax_t)

        # Grid definition
        if isnothing(model)
            cw = [0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
        else
            cw = [0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
        end
        # cw = [0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
        rh = [0.10, 0.10, 0.10, 0.10]
        CairoMakie.xlims!(ax_t, 0, sum(cw))
        CairoMakie.ylims!(ax_t, 0, sum(rh))

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

        if !isnothing(model)
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
            text="Pre-Contingency",  align=(:center,:center), fontsize=8, color=:green,  font=:bold)
        text!(ax_t, (cx[5]+cx[8])/2, cell_center_y(0);
            text="Post-Contingency", align=(:center,:center), fontsize=8, color=:darkorange, font=:bold)
        if !isnothing(model)
            text!(ax_t, (cx[8]+cx[9])/2, cell_center_y(0);
                text="TR", align=(:center,:center), fontsize=8, color=:steelblue, font=:bold)
        end
        # Row 1: column headers
        for (col, lbl) in enumerate(["", "Max", "Mean", "N", "Max", "Mean", "N"])
            text!(ax_t, cell_center_x(col-1), cell_center_y(1);
                text=lbl, align=(:center,:center), fontsize=8, color=:black, font=:bold)
        end

        # Rows 2-3: data
        foreach(col -> reprow[col] isa AbstractFloat && isnan(reprow[col]) && (reprow[col] = 0.0), names(reprow))
        XX = (
            voltage  = (pre  = (max  = max(reprow."Vmin_Abs", reprow."Vmax_Abs"),
                                mean = (reprow."Vmin_N_M" * reprow."Vmin_N_N" + reprow."Vmax_N_M" * reprow."Vmax_N_N") / (reprow."Vmin_N_N" + reprow."Vmax_N_N"),
                                n    = sum(reprow."Vmin_N_N" + reprow."Vmax_N_N")),
                        post = (max  = max(reprow."Vmin_C_Abs", reprow."Vmax_C_Abs"),
                                mean = (reprow."Vmin_C_M" * reprow."Vmin_C_N" + reprow."Vmax_C_M" * reprow."Vmax_C_N") / (reprow."Vmin_C_N" + reprow."Vmax_C_N"),
                                n    = sum(reprow."Vmin_C_N" + reprow."Vmax_C_N"))),
            lineflow = (pre  = (max  = reprow."Smax_N",   mean = reprow."Smax_N_M", n = reprow."Smax_N_N"),
                        post = (max  = reprow."Smax_C",   mean = reprow."Smax_C_M", n = reprow."Smax_C_N")),
        )
        table_rows = [
            ("Voltage (p.u)",   [XX.voltage.pre.max,   XX.voltage.pre.mean,   XX.voltage.pre.n,
                                  XX.voltage.post.max,  XX.voltage.post.mean,  XX.voltage.post.n]),
            ("Line Flow (p.u)", [XX.lineflow.pre.max,  XX.lineflow.pre.mean,  XX.lineflow.pre.n,
                                  XX.lineflow.post.max, XX.lineflow.post.mean, XX.lineflow.post.n]),
        ]

        for (r, (label, vals)) in enumerate(table_rows)
            row = r + 1
            text!(ax_t, cx[2] - 0.01, cell_center_y(row);
                text=label, align=(:right,:center), fontsize=9, color=:black, font=:bold)
            for (col, val) in enumerate(vals)
                text!(ax_t, cell_center_x(col), cell_center_y(row);
                    text  = col in [3,6] ? string(round(Int, val)) : @sprintf("%.3f", val),
                    align = (:center, :center),
                    fontsize = 9)
            end
        end
        if !isnothing(model)
            text!(ax_t, cell_center_x(7), cell_center_y(1);text  = "size",align = (:center, :center),fontsize = 9,color=:steelblue)
            text!(ax_t, cell_center_x(7), cell_center_y(2);text  = @sprintf("%.3f", model[:TR][2]),align = (:center, :center),fontsize = 9)
            text!(ax_t, cell_center_x(7), cell_center_y(3);text  = @sprintf("%.3f", model[:TR][1]),align = (:center, :center),fontsize = 9)
        end

        # colsize!(bottom, 1, Relative(0.5))
        colsize!(bottom, 1, 350)

    end

    rowsize!(fig.layout, 2, 50)

    # ── Legend (right) ────────────────────────────────────────
    Legend(bottom[1, reprow === nothing ? (1:2) : 2],
        [PolyElement(color = :lightgrey)],
        ["NC — no cascade and no load shed"],
        orientation = :vertical,
        tellwidth   = false,
        tellheight  = false,
        valign      = :center,
        halign      = :right
    )

    rowgap!(fig.layout, 1, 6)
    return fig
end
