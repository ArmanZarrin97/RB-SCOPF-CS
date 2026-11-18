using DelimitedFiles

# Example: If your text file is a delimited table or CSV
# Alternatively, you can construct the matrices directly from your parsed variables:
nCon = 33  # Number of contingencies (rows)
nTP = 24   # Number of time periods (columns)

# Initialize matrices for severity values and text annotations
LS_matrix = zeros(Float64, nCon, nTP)
Text_matrix = fill("NC", nCon, nTP)

# Populate your matrices here based on your loaded data file.
# For example, setting specific cells:
# LS_matrix[18, 1] = 212.52
# Text_matrix[18, 1] = "212.52"

using CSV
using DataFrames

# Load CSV file
df = CSV.read("contingency_grid.csv", DataFrame)

# Extract data columns (excluding the first 'Contingency' column)
data_matrix_raw = Matrix(df[:, 2:end])

# Create Text_matrix (strings) and LS_matrix (Float64 with 0.0 for "NC")
nCon, nTP = size(data_matrix_raw)
Text_matrix = string.(data_matrix_raw)
LS_matrix = zeros(Float64, nCon, nTP)

for c in 1:nCon, t in 1:nTP
    val_str = string(data_matrix_raw[c, t])
    if val_str == "NC" || val_str == "missing"
        LS_matrix[c, t] = 0.0
        Text_matrix[c, t] = "NC"
    else
        LS_matrix[c, t] = parse(Float64, val_str)
        Text_matrix[c, t] = val_str
    end
end

using CairoMakie
using Printf

function plot_cascade_manual(LS_matrix, Text_matrix; iter=5, ls_max=nothing)
    nCon, nTP = size(LS_matrix)

    # --- Adaptive sizing (from your original function) ---
    cell_w    = min(2000, 2000 / nTP)
    cell_h    = min(cell_w, 600 / nCon)
    fig_w     = min(2200, nTP * cell_w + 200)
    fig_h     = min(750, nCon * cell_h + 150)

    cell_size   = min(cell_w, cell_h)
    fsize_cell  = clamp(round(Int, cell_size / 9), 14, 26)
    fsize_ax    = clamp(round(Int, cell_size / 8), 16, 26)

    # Determine max scale if not provided
    if isnothing(ls_max)
        valid_vals = filter(x -> !isnan(x) && x > 0, LS_matrix)
        ls_max = isempty(valid_vals) ? 1.0 : maximum(valid_vals)
    end

    # Prepare plot matrix (treat 0.0 or explicit NC markers as NaN for gray background)
    LS_plot = copy(LS_matrix)
    LS_plot[LS_matrix .== 0.0] .= NaN

    # --- Figure Setup ---
    fig = Figure(size = (fig_w, fig_h), fontsize = fsize_ax)
    ax = Axis(fig[1, 1];
        title              = "Simulation — Severity per (Time, Contingency), iter=$(iter)",
        xlabel             = "Time Period",
        ylabel             = "Contingency",
        xticks             = (1:nTP, string.(1:nTP)),
        yticks             = (1:nCon, string.(1:nCon)),
        xticklabelrotation = π/4
    )

    hm = CairoMakie.heatmap!(ax, 1:nTP, 1:nCon, LS_plot';
        colormap   = cgrad([:mistyrose, :red, :darkred]),
        colorrange = (0, ls_max),
        nan_color  = :lightgrey
    )

    # --- Text Overlays ---
    for t in 1:nTP, c in 1:nCon
        val  = LS_matrix[c, t]
        txt  = Text_matrix[c, t]
        is_nc = (txt == "NC" || val == 0.0 || isnan(val))
        fontcolor = val >= 0.5 * ls_max ? :white : :black

        if is_nc
            text!(ax, t, c;
                text = "NC", align = (:center, :center), fontsize = fsize_cell, color = :grey)
        else
            text!(ax, t, c;
                text = isa(txt, Number) ? @sprintf("%.2f", txt) : string(txt),
                align = (:center, :center), fontsize = fsize_cell, color = fontcolor)
        end
    end

    Colorbar(fig[1, 2], hm;
        label    = "Sev (p.u.)",
        ticks    = LinearTicks(5),
        flipaxis = true
    )

    return fig
end

fig = plot_cascade_manual(LS_matrix, Text_matrix; iter=5, ls_max=215.0)
display(fig)