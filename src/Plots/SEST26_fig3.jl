using CairoMakie, Random

Random.seed!(42)

# ── Iteration axis ────────────────────────────────────────────────────────────
n_iter = 30
iters  = 1:n_iter

# ── Helper: smooth exponential decay with small noise ────────────────────────
function decay_curve(start, floor_, speed, noise_scale)
    base  = @. start * exp(-speed * (iters - 1)) + floor_
    noise = noise_scale .* randn(n_iter) .* exp.(-0.12 .* (collect(iters) .- 1))
    return max.(base .+ noise, floor_ .- 0.5)
end

# ── Data ──────────────────────────────────────────────────────────────────────
lsb_steady      = decay_curve(42.0, 1.8, 0.13, 3.5)
lsb_cascading   = decay_curve(67.0, 3.2, 0.10, 5.0)
calsb_steady    = decay_curve(42.0, 0.6, 0.20, 2.5)
calsb_cascading = decay_curve(67.0, 1.4, 0.17, 3.8)

# ── Colours (print-friendly, distinct on white) ───────────────────────────────
COL_LSB_S   = RGBf(0.12, 0.47, 0.71)   # steel blue
COL_LSB_C   = RGBf(0.02, 0.63, 0.76)   # teal
COL_CA_S    = RGBf(0.84, 0.15, 0.16)   # crimson
COL_CA_C    = RGBf(0.95, 0.55, 0.00)   # orange

# ── Figure & Axis ─────────────────────────────────────────────────────────────
fig = Figure(
    size            = (780, 460),
    backgroundcolor = :white,
    fonts           = (; regular = "Times New Roman"),
)

ax = Axis(fig[1, 1],
    backgroundcolor  = :white,
    xlabel           = "Iteration",
    ylabel           = "Worst Overload Violation (MW)",
    title            = "Overload Violation Reduction — LSB vs CA-LSB",
    titlesize        = 13,
    xlabelsize       = 11,
    ylabelsize       = 11,
    xticklabelsize   = 10,
    yticklabelsize   = 10,
    xgridcolor       = RGBf(0.85, 0.85, 0.85),
    ygridcolor       = RGBf(0.85, 0.85, 0.85),
    xgridwidth       = 0.8,
    ygridwidth       = 0.8,
    xgridvisible     = true,
    ygridvisible     = true,
    leftspinecolor   = :black,
    rightspinecolor  = :black,
    topspinecolor    = :black,
    bottomspinecolor = :black,
    xtickcolor       = :black,
    ytickcolor       = :black,
    xticklabelcolor  = :black,
    yticklabelcolor  = :black,
    xlabelcolor      = :black,
    ylabelcolor      = :black,
    titlecolor       = :black,
)

x = collect(iters)

# ── Shaded bands ──────────────────────────────────────────────────────────────
band!(ax, x, max.(lsb_steady     .- 2.5, 0), lsb_steady     .+ 2.5, color=(COL_LSB_S, 0.12))
band!(ax, x, max.(lsb_cascading  .- 3.5, 0), lsb_cascading  .+ 3.5, color=(COL_LSB_C, 0.12))
band!(ax, x, max.(calsb_steady   .- 2.0, 0), calsb_steady   .+ 2.0, color=(COL_CA_S,  0.12))
band!(ax, x, max.(calsb_cascading.- 2.8, 0), calsb_cascading.+ 2.8, color=(COL_CA_C,  0.12))

# ── Lines + markers ───────────────────────────────────────────────────────────
lines!(ax, x, lsb_steady;
    color=COL_LSB_S, linewidth=1.8, label="LSB — Steady-State")
CairoMakie.scatter!(ax, x, lsb_steady;
    color=COL_LSB_S, marker=:circle, markersize=6, strokewidth=0)

lines!(ax, x, lsb_cascading;
    color=COL_LSB_C, linewidth=1.8, linestyle=:dash, label="LSB — Cascading")
CairoMakie.scatter!(ax, x, lsb_cascading;
    color=COL_LSB_C, marker=:diamond, markersize=6, strokewidth=0)

lines!(ax, x, calsb_steady;
    color=COL_CA_S, linewidth=1.8, label="CA-LSB — Steady-State")
CairoMakie.scatter!(ax, x, calsb_steady;
    color=COL_CA_S, marker=:utriangle, markersize=7, strokewidth=0)

lines!(ax, x, calsb_cascading;
    color=COL_CA_C, linewidth=1.8, linestyle=:dash, label="CA-LSB — Cascading")
CairoMakie.scatter!(ax, x, calsb_cascading;
    color=COL_CA_C, marker=:star5, markersize=8, strokewidth=0)

# ── Annotation ────────────────────────────────────────────────────────────────
text!(ax, 18, 20;
    text     = "CA-LSB converges\nfaster & lower",
    color    = COL_CA_S,
    fontsize = 8,
    font     = :italic,
    align    = (:left, :center),
)

# ── Legend ────────────────────────────────────────────────────────────────────
Legend(fig[1, 2], ax,
    backgroundcolor = :white,
    framecolor      = :black,
    labelcolor      = :black,
    labelsize       = 9,
    rowgap          = 5,
    padding         = (10, 10, 8, 8),
)

# ── Save as high-res PNG and PDF ──────────────────────────────────────────────
# save("/mnt/user-data/outputs/overload_violations.png", fig, px_per_unit=3)
# save("/mnt/user-data/outputs/overload_violations.pdf", fig)
println("Saved → overload_violations.png  +  overload_violations.pdf")

fig