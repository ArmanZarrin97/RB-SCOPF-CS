using CairoMakie

case_labels = ["Mild", "Moderate", "Severe"]
Rcs = [(3,3) (4,3) (5,3)]
n = 50

# ── Data loading (vec() handles the 1×3 matrix from comprehension over Rcs) ──
time_data = vec([Float64.(MC_Results[Rc]."Time") for Rc in Rcs])
iter_data = vec([Float64.(MC_Results[Rc]."iter") for Rc in Rcs])
cost_data = vec([Float64.(MC_Results[Rc]."Cost") for Rc in Rcs])

# ── IEEE two-column figure dimensions ────────────────────────────────────────
width_mm  = 190
height_mm = 130   # extra height for legend row

pt_per_mm  = 2.83465
fig_width  = width_mm  * pt_per_mm
fig_height = height_mm * pt_per_mm

# ── Colour palette (colour-blind friendly, works in B&W) ─────────────────────
palette = [
    RGBf(0.337, 0.706, 0.914),   # steel blue   – Mild
    RGBf(0.933, 0.522, 0.290),   # burnt orange – Moderate
    RGBf(0.416, 0.800, 0.396),   # sage green   – Severe
]

# ── Helper: draw one density sub-panel ───────────────────────────────────────
function add_density_panel!(ax, datasets, xlabel)
    for (i, data) in enumerate(datasets)
        CairoMakie.density!(ax, data;
            color       = (palette[i], 0.25),
            strokecolor = palette[i],
            strokewidth = 1.2,
            label       = case_labels[i],
        )
    end

    ax.xlabel          = xlabel
    ax.ylabel          = "Density"
    ax.xlabelsize      = 11
    ax.ylabelsize      = 11
    ax.xticklabelsize  = 10
    ax.yticklabelsize  = 10
    ax.xgridvisible    = false
    ax.ygridvisible    = true
    ax.ygridcolor      = (:black, 0.08)
    ax.ygridwidth      = 0.5
    ax.spinewidth      = 0.6
    ax.xtickwidth      = 0.6
    ax.ytickwidth      = 0.6
    ax.topspinevisible   = false
    ax.rightspinevisible = false
end

# ── Build figure ─────────────────────────────────────────────────────────────
fig = Figure(
    size           = (fig_width, fig_height),
    fontsize       = 7,
    figure_padding = (4, 4, 6, 10),
    fonts          = (; regular = "Times New Roman",
                        bold    = "Times New Roman Bold"),
)

ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
ax3 = Axis(fig[3, 1])

add_density_panel!(ax1, time_data, "Solution time (s)")
add_density_panel!(ax2, iter_data, "Iterations")
add_density_panel!(ax3, cost_data, "Cost")

# Sub-figure labels (a), (b), (c) – IEEE style
for (ax, lbl) in zip([ax1, ax2, ax3], ["(a)", "(b)", "(c)"])
    text!(ax, 0.5, -0.18;
        text     = lbl,
        space    = :relative,
        align    = (:center, :top),
        fontsize = 7,
        font     = "Times New Roman Bold",
        color    = :black,
    )
end

# Shared legend below all panels
Legend(fig[1:3, 2], ax1;
    orientation  = :vertical,
    framevisible = false,
    labelsize    = 11,
    patchsize    = (12, 6),
)

colgap!(fig.layout, 8)
rowgap!(fig.layout, 4)

# ── Save ─────────────────────────────────────────────────────────────────────
Pathplot = Glob["Path"] * "../Result/Iter_Report/toplot/"
save(Pathplot * "PDF.pdf", fig)
fig