using CairoMakie, Random

#Figure 3: Overload Violation Reduction — LSB vs CA-LSB




# ── Data ──────────────────────────────────────────────────────────────────────
Pathplot=Glob["Path"]*"../Result/Iter_Report/toplot/" 
initp=[3,4,5,6]
method_code=[3]
if Test_Sys=="Polish/"
    obj_digits=1e6
    unit="M"
else
    obj_digits=1e3
    unit="K"
end
LSB_R=read_result(Pathplot,"SCOPF")
CALSB_R=read_result(Pathplot,"CA_SCOPF")
NLP_R=read_NLresult(Pathplot,"SCOPF")
NLP_R=NLP_R[(3,3)]

LSB_CR=read_result(Pathplot*"Cascade_report/","SCOPF","CR")
CALSB_CR=read_result(Pathplot*"Cascade_report/","CA_SCOPF","CR")
NLP_CR=read_NLresult(Pathplot*"Cascade_report/","SCOPF","CR")
NLP_CR=NLP_CR[(3,3)]

i=5



lsb_steady      = LSB_R[(i,3)]."Smax_C" .*100
lsb_cascading   = LSB_CR[(i,3)]."Worst_Cas_MW" .*100
calsb_steady    = CALSB_R[(i,3)]."Smax_C" .*100
calsb_cascading = CALSB_CR[(i,3)]."Worst_Cas_MW" .*100

lsb_time      = LSB_R[(i,3)]."Res_Time"
calsb_time    = CALSB_R[(i,3)]."Res_Time"

# ── Colours: colorblind-safe, print-friendly ─────────────────────────────────
COL_LSB_S  = RGBf(0.00, 0.45, 0.70)   # blue
COL_LSB_C  = RGBf(0.34, 0.71, 0.91)   # sky blue
COL_CA_S   = RGBf(0.84, 0.13, 0.13)   # red
COL_CA_C   = RGBf(0.95, 0.58, 0.10)   # amber

BAND_ALPHA  = 0.10
LINE_W      = 2.2
MARK_S      = 12

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(
    size            = (820, 480),
    backgroundcolor = :white,
    fonts           = (; regular = "Times New Roman", bold = "Times New Roman Bold"),
    figure_padding  = (16, 20, 12, 12),
)

ax = Axis(fig[1, 1],
    backgroundcolor  = :white,
    xlabel           = "Computation Time (s)",
    ylabel           = "Worst Overload Violation (MW)",
    title            = "Overload Violation Reduction — LSB vs CA-LSB",
    titlesize        = 14,
    titlefont        = :bold,
    titlegap         = 10,
    xlabelsize       = 12,
    ylabelsize       = 12,
    xticklabelsize   = 10,
    yticklabelsize   = 10,
    xgridcolor       = RGBf(0.88, 0.88, 0.88),
    ygridcolor       = RGBf(0.88, 0.88, 0.88),
    xgridwidth       = 0.7,
    ygridwidth       = 0.7,
    xgridvisible     = true,
    ygridvisible     = true,
    xminorgridvisible = false,
    yminorgridvisible = false,
    # Clean spine: only bottom + left
    rightspinevisible  = false,
    topspinevisible    = false,
    leftspinecolor     = RGBf(0.3, 0.3, 0.3),
    bottomspinecolor   = RGBf(0.3, 0.3, 0.3),
    xtickcolor         = RGBf(0.3, 0.3, 0.3),
    ytickcolor         = RGBf(0.3, 0.3, 0.3),
    xticklabelcolor    = :black,
    yticklabelcolor    = :black,
    xlabelcolor        = :black,
    ylabelcolor        = :black,
    titlecolor         = :black,
    xtickalign         = 1,
    ytickalign         = 1,
    ylabelpadding      = 8,
    xlabelpadding      = 8,
)

# ── Shaded uncertainty bands ──────────────────────────────────────────────────
band!(ax, lsb_time,   max.(lsb_steady     .- 2.5, 0), lsb_steady     .+ 2.5, color=(COL_LSB_S, BAND_ALPHA))
band!(ax, lsb_time,   max.(lsb_cascading  .- 3.5, 0), lsb_cascading  .+ 3.5, color=(COL_LSB_C, BAND_ALPHA))
band!(ax, calsb_time, max.(calsb_steady   .- 2.0, 0), calsb_steady   .+ 2.0, color=(COL_CA_S,  BAND_ALPHA))
band!(ax, calsb_time, max.(calsb_cascading.- 2.8, 0), calsb_cascading.+ 2.8, color=(COL_CA_C,  BAND_ALPHA))

# ── Lines + markers ───────────────────────────────────────────────────────────
lines!(ax, lsb_time, lsb_steady;
    color=COL_LSB_S, linewidth=LINE_W, label="LSB — Steady-State")
CairoMakie.scatter!(ax, lsb_time, lsb_steady;
    color=COL_LSB_S, marker=:diamond, markersize=MARK_S, strokewidth=0.8, strokecolor=:white)

lines!(ax, lsb_time, lsb_cascading;
    color=COL_LSB_C, linewidth=LINE_W, linestyle=:dash, label="LSB — Cascading")
CairoMakie.scatter!(ax, lsb_time, lsb_cascading;
    color=COL_LSB_C, marker=:diamond, markersize=MARK_S, strokewidth=0.8, strokecolor=:white)

lines!(ax, calsb_time, calsb_steady;
    color=COL_CA_S, linewidth=LINE_W, label="CA-LSB — Steady-State")
CairoMakie.scatter!(ax, calsb_time, calsb_steady;
    color=COL_CA_S, marker=:star5, markersize=MARK_S+1, strokewidth=0.8, strokecolor=:white)

lines!(ax, calsb_time, calsb_cascading;
    color=COL_CA_C, linewidth=LINE_W, linestyle=:dash, label="CA-LSB — Cascading")
CairoMakie.scatter!(ax, calsb_time, calsb_cascading;
    color=COL_CA_C, marker=:star5, markersize=MARK_S+2, strokewidth=0.8, strokecolor=:white)

# ── Legend (inside plot, top-right) ──────────────────────────────────────────
axislegend(ax,
    position        = :rt,          # inside: right-top
    backgroundcolor = (:white, 0.88),
    framecolor      = RGBf(0.7, 0.7, 0.7),
    framewidth      = 0.8,
    labelcolor      = :black,
    labelsize       = 9.5,
    rowgap          = 4,
    padding         = (10, 10, 8, 8),
    patchsize       = (22, 12),
)

# ── Save ──────────────────────────────────────────────────────────────────────
# save(Pathplot * "overload_violations.png", fig, px_per_unit=3)
# save(Pathplot * "overload_violations.pdf", fig)
# println("Saved → overload_violations.png  +  overload_violations.pdf")

fig