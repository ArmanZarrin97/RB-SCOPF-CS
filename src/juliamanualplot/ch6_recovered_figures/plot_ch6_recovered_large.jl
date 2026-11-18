using CairoMakie
using Printf


"""
    plot_ch6_figure(
        data;
        outpath=data.suggested_filename,
        size=(2600, 1700),
        ...
    )

Recreates one of the recovered Chapter 6 heatmap figures.

The recovered data use:
    severity[c, t]

with:
    c = contingency 1:33
    t = time period 1:24
    NaN = "NC"

Compared with the earlier version:
  * the detailed bottom security table is removed;
  * only the two risk/severity indicators are retained;
  * they are displayed horizontally;
  * their displayed names are:
        "Maximum simplified Sev"
        "Overall system risk"
  * axis labels, tick labels, and cell values are larger.
"""
function plot_ch6_figure(
    data;
    outpath = data.suggested_filename,

    # Larger and less shallow than the original thesis figure.
    size = (2600, 1700),

    # Main heatmap typography.
    cell_font = 27,
    tick_font = 28,
    axis_font = 36,
    title_font = 36,

    # Colorbar typography.
    colorbar_font = 36,

    # Bottom risk indicators.
    risk_label_font = 36,
    risk_value_font = 36,
)

    # ================================================================
    # Data
    # ================================================================

    sev = Matrix{Float64}(data.severity)

    # Base.size is intentional because "size" is also a keyword above.
    nCon, nTP = Base.size(sev)

    nCon == 33 || @warn "Expected 33 contingencies; found $nCon."
    nTP  == 24 || @warn "Expected 24 time periods; found $nTP."


    # ================================================================
    # Figure
    # ================================================================

    fig = Figure(
        size = size,
        fontsize = axis_font
    )


    # ================================================================
    # Main heatmap
    # ================================================================

    ax = Axis(
        fig[1, 1];

        title = data.title,
        xlabel = "Time Period",
        ylabel = "Contingency",

        titlesize = title_font,
        xlabelsize = axis_font,
        ylabelsize = axis_font,

        xticks = (
            1:nTP,
            string.(1:nTP)
        ),

        yticks = (
            1:nCon,
            string.(1:nCon)
        ),

        xticklabelsize = tick_font,
        yticklabelsize = tick_font,

        xticklabelrotation = π / 4,

        # Set limits directly in Axis to avoid xlims!/ylims! conflicts.
        limits = (
            0.5,
            nTP + 0.5,
            0.5,
            nCon + 0.5
        ),
    )


    hm = CairoMakie.heatmap!(
        ax,
        1:nTP,
        1:nCon,
        sev';

        colormap = cgrad([
            :mistyrose,
            :red,
            :darkred
        ]),

        colorrange = (
            0.0,
            data.color_max
        ),

        nan_color = :lightgrey,
    )


    # ================================================================
    # Cell annotations
    # ================================================================

    for t in 1:nTP, c in 1:nCon

        v = sev[c, t]

        if isnan(v)

            text!(
                ax,
                t,
                c;

                text = "NC",
                align = (:center, :center),
                fontsize = cell_font,
                color = :grey,
            )

        else

            fontcolor =
                v >= 0.5 * data.color_max ?
                :white :
                :black

            text!(
                ax,
                t,
                c;

                text = @sprintf("%.2f", v),
                align = (:center, :center),
                fontsize = cell_font,
                color = fontcolor,
            )

        end
    end


    # ================================================================
    # Colorbar
    # ================================================================

    Colorbar(
        fig[1, 2],
        hm;

        label = "Sev (p.u.)",
        ticks = 0:50:200,
        flipaxis = true,

        labelsize = colorbar_font,
        ticklabelsize = colorbar_font,
    )


    # ================================================================
    # Bottom row: only the two requested indicators
    # ================================================================

    s = data.summary

    riskrow = GridLayout(
        fig[2, 1:2]
    )


    # Maximum simplified Sev
    Label(
        riskrow[1, 1],
        "Maximum ESev:";
        fontsize = risk_label_font,
        font = :bold,
        color = :steelblue,
        halign = :right,
        valign = :center,
    )

    Label(
        riskrow[1, 2],
        @sprintf("%.3f", s.max_op_risk);
        fontsize = risk_value_font,
        halign = :left,
        valign = :center,
    )


    # Overall system risk
    Label(
        riskrow[1, 3],
        "Overall Risk:";
        fontsize = risk_label_font,
        font = :bold,
        color = :steelblue,
        halign = :right,
        valign = :center,
    )

    Label(
        riskrow[1, 4],
        @sprintf("%.3f", s.system_risk);
        fontsize = risk_value_font,
        halign = :left,
        valign = :center,
    )


    # Keep the summary compact so most of the page is used by the heatmap.
    rowsize!(
        fig.layout,
        2,
        105
    )

    rowgap!(
        fig.layout,
        1,
        16
    )


    # ================================================================
    # Save
    # ================================================================

    save(
        outpath,
        fig
    )

    return fig
end



# ====================================================================
# Plot all four recovered Chapter 6 figures
# ====================================================================

function plot_all_ch6_figures(
    ;
    outdir = "Chapter6_Recovered_Plots",
    kwargs...
)

    mkpath(outdir)

    for data in (
        FIG6_3,
        FIG6_4,
        FIG6_5,
        FIG6_6
    )

        outpath = joinpath(
            outdir,
            data.suggested_filename
        )

        plot_ch6_figure(
            data;
            outpath = outpath,
            kwargs...
        )

        println(
            "Saved: ",
            outpath
        )

    end

    return nothing
end
