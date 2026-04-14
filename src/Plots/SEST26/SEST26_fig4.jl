using CairoMakie

case_labels = ["Mild", "Moderate", "Severe"]
Rcs = [(3,3) (4,3) (5,3)]
n = 49
Plots.plot(MC_Results[(4,3)]."Time")
MC_Results[(4,3)]."Time"[3,:]

MC_Results[(4,3)]=MC_Results[(4,3)][Not(3),:]
# ── Data ──────────────────────────────────────────────────────────────────────
time_data = vec([Float64.(MC_Results[Rc]."Time") for Rc in Rcs])
iter_data = vec([Float64.(MC_Results[Rc]."iter") for Rc in Rcs])
cost_data = vec([Float64.(MC_Results[Rc]."Cost") for Rc in Rcs])

# ── IEEE two-column dimensions ────────────────────────────────────────────────
pt_per_mm  = 2.83465
fig_width  = 190 * pt_per_mm
fig_height = 120 * pt_per_mm

# ── Palette: muted, print-safe ────────────────────────────────────────────────
palette = [
    RGBf(0.267, 0.537, 0.729),   # muted blue    – Mild
    RGBf(0.820, 0.400, 0.180),   # muted orange  – Moderate
    RGBf(0.290, 0.627, 0.310),   # muted green   – Severe
]

# ── Helper ────────────────────────────────────────────────────────────────────
function add_boxplot_panel!(ax, datasets, ylabel; integer_ticks=false)
    for (i, data) in enumerate(datasets)
        CairoMakie.boxplot!(ax, fill(i, length(data)), data;
            color            = (palette[i], 0.55),
            strokecolor      = palette[i],
            strokewidth      = 1.2,
            mediancolor      = :white,
            medianlinewidth  = 2.0,
            whiskerwidth     = 0.4,
            outliercolor     = (palette[i], 0.7),
            outlierstrokecolor = palette[i],
            outlierstrokewidth = 0.8,
            markersize       = 5,
        )
    end

    ax.xticks            = (1:3, case_labels)
    ax.ylabel            = ylabel
    ax.ylabelsize        = 10
    ax.xlabelsize        = 10
    ax.xticklabelsize    = 10
    ax.yticklabelsize    = 9
    ax.xgridvisible      = false
    ax.ygridvisible      = true
    ax.ygridcolor        = (:black, 0.10)
    ax.ygridwidth        = 0.5
    ax.yminorgridvisible = true
    ax.yminorgridcolor   = (:black, 0.04)
    ax.yminorgridwidth   = 0.4
    ax.spinewidth        = 0.8
    ax.xtickwidth        = 0.8
    ax.ytickwidth        = 0.8
    ax.xticksize         = 4
    ax.yticksize         = 4
    ax.topspinevisible   = false
    ax.rightspinevisible = false
    ax.bottomspinevisible = true
    ax.leftspinevisible  = true
    ax.xticksmirrored    = false
    ax.yticksmirrored    = false

    # x-axis limits with small padding so boxes don't touch the spines
    CairoMakie.xlims!(ax, 0.4, 3.6)

    if integer_ticks
        lo, hi = extrema(vcat(datasets...))
        step = max(1, round(Int, (hi - lo) / 5 / n) * n)
        nice_ticks = range(floor(Int, lo / step) * step,
                           ceil(Int,  hi / step) * step; step=step)
        ax.yticks = nice_ticks
    end
end

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(
    size           = (fig_width, fig_height),
    fontsize       = 10,
    figure_padding = (6, 14, 6, 6),
    fonts          = (; regular = "Times New Roman",
                        bold    = "Times New Roman Bold"),
)

# Thin horizontal line separating content from any future caption space
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[1, 3])

add_boxplot_panel!(ax1, time_data, "Solution time (s)")
add_boxplot_panel!(ax2, iter_data, "No. of iterations"; integer_ticks=true)
add_boxplot_panel!(ax3, cost_data, "Cost (€)")


# Sub-figure labels (a) (b) (c) — placed inside top-left of each panel
for (ax, lbl) in zip([ax1, ax2, ax3], ["(a)", "(b)", "(c)"])
    text!(ax, 0.03, 0.97;
        text     = lbl,
        space    = :relative,
        align    = (:left, :top),
        fontsize = 10,
        font     = "Times New Roman Bold",
        color    = :black,
    )
end

colgap!(fig.layout, 6)
rowgap!(fig.layout, 0)

# ── Save ──────────────────────────────────────────────────────────────────────
Pathplot = Glob["Path"] * "../Result/Iter_Report/toplot/"
save(Pathplot * "Boxplots.pdf", fig)
fig