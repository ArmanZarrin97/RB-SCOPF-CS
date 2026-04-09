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

scale_factor = 0.9  # adjust as needed
BAND_ALPHA  = 0.10* scale_factor
LINE_W      = 2.2* scale_factor
MARK_S      = 12* scale_factor
font_scale=1.2;
# ── Figure: 2 rows, shared x-axis ────────────────────────────────────────────
fig = Figure(
    size            = (500, 600) .* scale_factor,   # taller to fit two panels
    backgroundcolor = :white,
    fontsize = 12 * scale_factor*font_scale,
    # fonts           = (; regular = "Times New Roman", bold = "Times New Roman Bold"),
    figure_padding  = (16, 20, 12, 12) .* scale_factor,
)

ax1 = Axis(fig[1, 1],
    backgroundcolor  = :white,
    ylabel           = "Worst Overload Violation (MW)",
    # title            = "Overload Violation Reduction — LSB vs CA-LSB",
    # titlesize        = 14,
    # titlefont        = :bold,
    # titlegap         = 10,
    xlabelsize       = 12* scale_factor*font_scale,
    ylabelsize       = 12* scale_factor*font_scale,
    xticklabelsize   = 10* scale_factor*font_scale,
    yticklabelsize   = 10* scale_factor*font_scale,
    xgridcolor       = RGBf(0.88, 0.88, 0.88),
    ygridcolor       = RGBf(0.88, 0.88, 0.88),
    xgridwidth       = 0.7* scale_factor,
    ygridwidth       = 0.7* scale_factor,
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
    ylabelpadding      = 8* scale_factor,
    xticklabelsvisible = false,   # hide x labels on top panel
)

ax2 = Axis(fig[2, 1],
    backgroundcolor  = :white,
    xlabel           = "Computation Time (s)",
    ylabel           = "Max Voltage Violation (p.u.)",
    titlesize        = 14* scale_factor,
    titlefont        = :bold,
    titlegap         = 10* scale_factor*font_scale,
    xlabelsize       = 12* scale_factor*font_scale,
    ylabelsize       = 12* scale_factor*font_scale,
    xticklabelsize   = 10* scale_factor*font_scale,
    yticklabelsize   = 10* scale_factor*font_scale,
    xgridcolor       = RGBf(0.88, 0.88, 0.88),
    ygridcolor       = RGBf(0.88, 0.88, 0.88),
    xgridwidth       = 0.7* scale_factor,
    ygridwidth       = 0.7* scale_factor,
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
    ylabelpadding      = 8* scale_factor,
    xlabelpadding      = 8* scale_factor,
)

ax1.xticks = range(0, 10; step=1)
ax2.xticks = range(0, 10; step=1)
# ax2.xticks = range(0, 10; length=10)
# ax2.xticks = (vmin, vmax) -> range(vmin, vmax; length=9)

# ── Shaded uncertainty bands ──────────────────────────────────────────────────
# band!(ax1, lsb_time,   max.(lsb_steady     .- 2.5, 0), lsb_steady     .+ 2.5, color=(COL_LSB_S, BAND_ALPHA))
# band!(ax1, lsb_time,   max.(lsb_cascading  .- 3.5, 0), lsb_cascading  .+ 3.5, color=(COL_LSB_C, BAND_ALPHA))
# band!(ax1, calsb_time, max.(calsb_steady   .- 2.0, 0), calsb_steady   .+ 2.0, color=(COL_CA_S,  BAND_ALPHA))
# band!(ax1, calsb_time, max.(calsb_cascading.- 2.8, 0), calsb_cascading.+ 2.8, color=(COL_CA_C,  BAND_ALPHA))

# ── Lines + markers ───────────────────────────────────────────────────────────
lines!(ax1, lsb_time, lsb_steady;
    color=COL_LSB_S, linewidth=LINE_W, label="SCOPF — Steady-State")
CairoMakie.scatter!(ax1, lsb_time, lsb_steady;
    color=COL_LSB_S, marker=:diamond, markersize=MARK_S, strokewidth=0.8* scale_factor, strokecolor=:white)

lines!(ax1, lsb_time, lsb_cascading;
    color=COL_LSB_C, linewidth=LINE_W, linestyle=:dash, label="SCOPF — Cascading")
CairoMakie.scatter!(ax1, lsb_time, lsb_cascading;
    color=COL_LSB_C, marker=:diamond, markersize=MARK_S, strokewidth=0.8* scale_factor, strokecolor=:white)

lines!(ax1, calsb_time, calsb_steady;
    color=COL_CA_S, linewidth=LINE_W, label="CA-SCOPF — Steady-State")
CairoMakie.scatter!(ax1, calsb_time, calsb_steady;
    color=COL_CA_S, marker=:star5, markersize=MARK_S+1* scale_factor, strokewidth=0.8* scale_factor, strokecolor=:white)

lines!(ax1, calsb_time, calsb_cascading;
    color=COL_CA_C, linewidth=LINE_W, linestyle=:dash, label="CA-SCOPF — Cascading")
CairoMakie.scatter!(ax1, calsb_time, calsb_cascading;
    color=COL_CA_C, marker=:star5, markersize=MARK_S+2* scale_factor, strokewidth=0.8* scale_factor, strokecolor=:white)

# ── Legend (inside plot, top-right) ──────────────────────────────────────────
axislegend(ax1,
    position        = :rt,          # inside: right-top
    backgroundcolor = (:white, 0.88),
    framecolor      = RGBf(0.7, 0.7, 0.7),
    framewidth      = 0.8* scale_factor,
    labelcolor      = :black,
    labelsize       = 9.5* scale_factor*font_scale,
    rowgap          = 4* scale_factor,
    padding         = (10, 10, 8, 8) .* scale_factor,
    patchsize       = (22, 12) .* scale_factor,
)

# ── Voltage data (adjust field names to match your result struct) ─────────────
lsb_voltage   = max.(LSB_R[(i,3)]."Vmax_C",LSB_R[(i,3)]."Vmin_C")    # or whatever your field is called
calsb_voltage = max.(CALSB_R[(i,3)]."Vmax_C",CALSB_R[(i,3)]."Vmin_C")

# ── Voltage panel (ax2) ───────────────────────────────────────────────────────
# band!(ax2, lsb_time,   max.(lsb_voltage   .- 0.005, 0), lsb_voltage   .+ 0.005, color=(COL_LSB_S, BAND_ALPHA))
# band!(ax2, calsb_time, max.(calsb_voltage  .- 0.004, 0), calsb_voltage .+ 0.004, color=(COL_CA_S,  BAND_ALPHA))

lines!(ax2, lsb_time, lsb_voltage;
    color=COL_LSB_S, linewidth=LINE_W, label="SCOPF")
CairoMakie.scatter!(ax2, lsb_time, lsb_voltage;
    color=COL_LSB_S, marker=:diamond, markersize=MARK_S, strokewidth=0.8* scale_factor, strokecolor=:white)

lines!(ax2, calsb_time, calsb_voltage;
    color=COL_CA_S, linewidth=LINE_W, label="CA-SCOPF")
CairoMakie.scatter!(ax2, calsb_time, calsb_voltage;
    color=COL_CA_S, marker=:star5, markersize=MARK_S+1* scale_factor, strokewidth=0.8* scale_factor, strokecolor=:white)

axislegend(ax2, position=:rt,
    backgroundcolor = (:white, 0.88),
    framecolor      = RGBf(0.7, 0.7, 0.7),
    framewidth      = 0.8* scale_factor,
    labelcolor      = :black,
    labelsize       = 9.5* scale_factor*font_scale,
    rowgap          = 4* scale_factor,
    padding         = (10, 10, 8, 8) .* scale_factor,
    patchsize       = (22, 12) .* scale_factor,
)
# Link x-axes so zoom/pan stays in sync
linkxaxes!(ax1, ax2)

# Reduce gap between the two panels
rowgap!(fig.layout, 1, 6* scale_factor)


# ── Save ──────────────────────────────────────────────────────────────────────
save(Pathplot * "violations.png", fig)
save(Pathplot * "violations.pdf", fig)
# println("Saved → overload_violations.png  +  overload_violations.pdf")

fig
