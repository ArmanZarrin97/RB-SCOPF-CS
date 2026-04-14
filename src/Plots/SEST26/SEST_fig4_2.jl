using CairoMakie

case_labels = ["Mild", "Moderate", "Severe"]
Rcs = [(3, 3), (4, 3), (5, 3)]

# ── Data ──────────────────────────────────────────────────────────────────────
time_data = [Float64.(MC_Results[Rc]."Time") for Rc in Rcs]
iter_data = [Float64.(MC_Results[Rc]."iter") for Rc in Rcs]
cost_data = [Float64.(MC_Results[Rc]."Cost") for Rc in Rcs]

# ── IEEE two-column dimensions ────────────────────────────────────────────────
pt_per_mm = 2.83465
fig_width  = 190 * pt_per_mm
fig_height = 100 * pt_per_mm  # Slightly taller for legend space

# ── Consistent palette ────────────────────────────────────────────────────────
palette = [
    RGBf(0.267, 0.537, 0.729),  # Mild
    RGBf(0.820, 0.400, 0.180),  # Moderate
    RGBf(0.290, 0.627, 0.310),  # Severe
]

function style_axis!(ax, ylabel, datasets; integer_ticks=false)
    ax.xlabel = ""
    ax.xticklabelsvisible = false
    ax.xticksvisible = false
    ax.ylabel = ylabel
    ax.ylabelsize = 14
    ax.yticklabelsize = 11
    ax.ygridvisible = true
    ax.ygridcolor = (:black, 0.10)
    ax.ygridwidth = 0.5
    ax.yminorgridvisible = true
    ax.yminorgridcolor = (:black, 0.04)
    ax.yminorgridwidth = 0.4
    ax.spinewidth = 0.8
    ax.ytickwidth = 0.8
    ax.yticksize = 4
    ax.topspinevisible = false
    ax.rightspinevisible = false
    ax.bottomspinevisible = true
    ax.leftspinevisible = true
    CairoMakie.xlims!(ax, 0.4, 3.6)

    if integer_ticks
        all_data = vcat(map(d -> round.(Int, d), datasets)...)
        lo, hi = extrema(all_data)
        step = max(1, round(Int, (hi - lo) / 5))
        ax.yticks = lo:step:hi
    end
end

function add_boxplot_panel!(ax, datasets, ylabel; integer_ticks=false)
    for i in eachindex(datasets)
        data = datasets[i]
        CairoMakie.boxplot!(ax, fill(i, length(data)), data;
            color = (palette[i], 0.55),
            strokecolor = palette[i],
            strokewidth = 1.2,
            mediancolor = :white,
            medianlinewidth = 2.0,
            whiskerwidth = 0.4,
            outliercolor = (palette[i], 0.7),
            outlierstrokecolor = palette[i],
            outlierstrokewidth = 0.8,
            markersize = 5,
        )
    end
    style_axis!(ax, ylabel, datasets; integer_ticks=integer_ticks)
end

fig = Figure(
    size = (fig_width, fig_height),
    fontsize = 12,
    figure_padding = (6, 10, 2, 6),  # Less bottom padding
    fonts = (; regular = "Times New Roman Bold", bold = "Times New Roman Bold"),
)

# Subplots in row 1
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[1, 3])

add_boxplot_panel!(ax1, time_data, "Solution time (s)")
add_boxplot_panel!(ax2, iter_data, "No. of iterations"; integer_ticks=true)
add_boxplot_panel!(ax3, cost_data, "Cost (€)")

# Subfigure labels
for (ax, lbl) in zip([ax1, ax2, ax3], ["(a)", "(b)", "(c)"])
    text!(ax, 0.03, 0.97;
        text = lbl,
        space = :relative,
        align = (:left, :top),
        fontsize = 10,
        font = "Times New Roman Bold",
        color = :black,
    )
end

# Legend BELOW - spans all 3 columns in row 2
legend_elements = [
    PolyElement(polycolor = (palette[i], 0.55), strokecolor = palette[i])
    for i in eachindex(case_labels)
]
Legend(fig[2, 1:3], legend_elements, case_labels;
    orientation = :horizontal,
    framevisible = false,
    # nbanks = 3,  # One per column
    patchsize = (12, 10),
    labelsize = 14,
    titlesize = 10,
    rowgap = 2,
    halign = :center,
)

colgap!(fig.layout, 6)
rowgap!(fig.layout, 8)  # Space between plots and legend

# ── Save ──────────────────────────────────────────────────────────────────────
Pathplot = Glob["Path"] * "../Result/Iter_Report/toplot/"
save(Pathplot * "Boxplots.pdf", fig)
fig