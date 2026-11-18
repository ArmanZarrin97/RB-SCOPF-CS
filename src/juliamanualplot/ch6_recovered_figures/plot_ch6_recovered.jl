using CairoMakie
using Printf


"""
    plot_ch6_figure(
        data;
        outpath=data.suggested_filename,
        size=(2400, 1500),
        ...
    )

Recreates one of the recovered Chapter 6 heatmap figures.

The data files store:
  * `severity[c,t]` for contingency c=1:33 and time period t=1:24
  * `NaN` for cells displayed as "NC"
  * the bottom conventional-security summary
  * the optional `inacc` or `TR size` column
  * Max OP Risk and System Risk

The default aspect ratio is intentionally less shallow than the original
thesis figure so that the 33 contingency rows and their annotations are
much easier to read after the PDF is inserted as a sideways figure in LaTeX.
"""
function plot_ch6_figure(
    data;
    outpath = data.suggested_filename,
    size = (2400, 1500),
    cell_font = 22,
    tick_font = 22,
    axis_font = 25,
    title_font = 28,
    table_font = 20,
    legend_font = 20,
    colorbar_font = 22,
)

    # ================================================================
    # Data
    # ================================================================

    sev = Matrix{Float64}(data.severity)

    # Base.size is intentional because the keyword argument above
    # is also called "size".
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

        # No xlims! / ylims! required.
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

        colormap = cgrad(
            [
                :mistyrose,
                :red,
                :darkred
            ]
        ),

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
    # Bottom summary
    # ================================================================

    bottom = GridLayout(
        fig[2, 1:2]
    )


    s = data.summary

    has_extra =
        s.extra_label !== nothing


    # ================================================================
    # Table geometry
    # ================================================================

    cw = if has_extra

        [
            0.12,
            0.075,
            0.075,
            0.075,
            0.075,
            0.075,
            0.075,
            0.075
        ]

    else

        [
            0.12,
            0.075,
            0.075,
            0.075,
            0.075,
            0.075,
            0.075
        ]

    end


    rh = [
        0.10,
        0.10,
        0.10,
        0.10
    ]


    cx = cumsum(
        [0.0; cw]
    )

    ry = cumsum(
        [0.0; rh]
    )


    H = sum(rh)


    flipy(y) =
        H - y


    cell_center_x(col) =
        (
            cx[col + 1] +
            cx[col + 2]
        ) / 2


    cell_center_y(row) =
        flipy(
            (
                ry[row + 1] +
                ry[row + 2]
            ) / 2
        )


    # ================================================================
    # Space for risk box
    # ================================================================

    risk_x0 =
        sum(cw) + 0.025

    risk_x1 =
        risk_x0 + 0.105


    # ================================================================
    # Bottom table Axis
    #
    # IMPORTANT:
    # We create it here because risk_x1 and H are now known.
    # This avoids xlims! and ylims! completely.
    # ================================================================

    ax_t = Axis(
        bottom[1, 1];

        limits = (
            0.0,
            risk_x1,
            0.0,
            H
        )
    )


    hidedecorations!(
        ax_t
    )

    hidespines!(
        ax_t
    )


    # ================================================================
    # Header backgrounds
    # ================================================================

    for (cols, clr) in [

        (
            1:3,
            (:green, 0.15)
        ),

        (
            4:6,
            (:darkorange, 0.15)
        )

    ]

        x0 =
            cx[first(cols) + 1]

        x1 =
            cx[last(cols) + 2]

        y0 =
            flipy(
                ry[2]
            )

        y1 =
            flipy(
                ry[1]
            )


        poly!(
            ax_t,

            Point2f[
                (x0, y0),
                (x1, y0),
                (x1, y1),
                (x0, y1)
            ];

            color = clr,
            strokecolor = :lightgrey,
            strokewidth = 0.7,
        )

    end


    # ================================================================
    # Extra column background: inacc / TR
    # ================================================================

    if has_extra

        x0, x1 =
            cx[8],
            cx[9]

        y0, y1 =
            flipy(ry[3]),
            flipy(ry[1])


        poly!(
            ax_t,

            Point2f[
                (x0, y0),
                (x1, y0),
                (x1, y1),
                (x0, y1)
            ];

            color = (
                :steelblue,
                0.15
            ),

            strokecolor = :lightgrey,
            strokewidth = 0.7,
        )

    end


    # ================================================================
    # Risk box background
    # ================================================================

    poly!(
        ax_t,

        Point2f[
            (risk_x0, 0),
            (risk_x1, 0),
            (risk_x1, H),
            (risk_x0, H)
        ];

        color = (
            :grey,
            0.15
        ),

        strokecolor = :lightgrey,
        strokewidth = 0.7,
    )


    # ================================================================
    # Vertical grid lines
    # ================================================================

    for x in cx[2:end]

        lines!(
            ax_t,
            [x, x],
            [0, H - rh[end]];

            color = :lightgrey,
            linewidth = 0.7
        )

    end


    lines!(
        ax_t,

        [
            cx[5],
            cx[5]
        ],

        [
            0,
            H
        ];

        color = :grey,
        linewidth = 1.2
    )


    lines!(
        ax_t,

        [
            cx[8],
            cx[8]
        ],

        [
            0,
            H
        ];

        color = :grey,
        linewidth = 1.2
    )


    # ================================================================
    # Pre-contingency header
    # ================================================================

    text!(
        ax_t,

        (
            cx[2] +
            cx[5]
        ) / 2,

        cell_center_y(0);

        text = "Pre-Contingency",
        align = (:center, :center),
        fontsize = table_font,
        color = :green,
        font = :bold,
    )


    # ================================================================
    # Post-contingency header
    # ================================================================

    text!(
        ax_t,

        (
            cx[5] +
            cx[8]
        ) / 2,

        cell_center_y(0);

        text = "Post-Contingency",
        align = (:center, :center),
        fontsize = table_font,
        color = :darkorange,
        font = :bold,
    )


    # ================================================================
    # Extra column header
    # ================================================================

    if has_extra

        text!(
            ax_t,

            (
                cx[8] +
                cx[9]
            ) / 2,

            (
                cell_center_y(0) +
                cell_center_y(1)
            ) / 2;

            text = s.extra_label,
            align = (:center, :center),
            fontsize = table_font,
            color = :steelblue,
            font = :bold,
        )

    end


    # ================================================================
    # Max / Mean / N column headers
    # ================================================================

    for (col, lbl) in enumerate(
        [
            "",
            "Max",
            "Mean",
            "N",
            "Max",
            "Mean",
            "N"
        ]
    )

        text!(
            ax_t,

            cell_center_x(
                col - 1
            ),

            cell_center_y(1);

            text = lbl,
            align = (:center, :center),
            fontsize = table_font,
            color = :black,
            font = :bold,
        )

    end


    # ================================================================
    # Formatting helper
    # ================================================================

    function fmt_mean(v)

        isnan(v) ?
        "NaN" :
        @sprintf("%.3f", v)

    end


    # ================================================================
    # Voltage table values
    # ================================================================

    voltage_vals = [

        @sprintf(
            "%.3f",
            s.voltage_pre.max
        ),

        fmt_mean(
            s.voltage_pre.mean
        ),

        string(
            s.voltage_pre.n
        ),

        @sprintf(
            "%.3f",
            s.voltage_post.max
        ),

        fmt_mean(
            s.voltage_post.mean
        ),

        string(
            s.voltage_post.n
        ),

    ]


    # ================================================================
    # Line-flow table values
    # ================================================================

    line_vals = [

        @sprintf(
            "%.3f",
            s.line_pre.max
        ),

        fmt_mean(
            s.line_pre.mean
        ),

        string(
            s.line_pre.n
        ),

        @sprintf(
            "%.3f",
            s.line_post.max
        ),

        fmt_mean(
            s.line_post.mean
        ),

        string(
            s.line_post.n
        ),

    ]


    # ================================================================
    # Draw table values
    # ================================================================

    for (r, (label, vals)) in enumerate(
        [

            (
                "Voltage (p.u.)",
                voltage_vals
            ),

            (
                "Line Flow (p.u.)",
                line_vals
            )

        ]
    )

        row =
            r + 1


        text!(
            ax_t,

            cx[2] - 0.012,

            cell_center_y(row);

            text = label,
            align = (:right, :center),
            fontsize = table_font,
            color = :black,
            font = :bold,
        )


        for (col, val) in enumerate(vals)

            text!(
                ax_t,

                cell_center_x(col),

                cell_center_y(row);

                text = val,
                align = (:center, :center),
                fontsize = table_font,
            )

        end

    end


    # ================================================================
    # Extra column values
    # ================================================================

    if has_extra

        ev =
            s.extra_voltage === nothing ?
            "-" :
            @sprintf(
                "%.3f",
                s.extra_voltage
            )


        el =
            s.extra_line === nothing ?
            "-" :
            @sprintf(
                "%.3f",
                s.extra_line
            )


        text!(
            ax_t,

            cell_center_x(7),

            cell_center_y(2);

            text = ev,
            align = (:center, :center),
            fontsize = table_font,
        )


        text!(
            ax_t,

            cell_center_x(7),

            cell_center_y(3);

            text = el,
            align = (:center, :center),
            fontsize = table_font,
        )

    end


    # ================================================================
    # Risk block
    # ================================================================

    risk_x =
        (
            risk_x0 +
            risk_x1
        ) / 2


    # Max OP Risk label

    text!(
        ax_t,

        risk_x,

        cell_center_y(0);

        text = "Max OP Risk",
        align = (:center, :center),
        fontsize = table_font,
        color = :steelblue,
        font = :bold,
    )


    # Max OP Risk value

    text!(
        ax_t,

        risk_x,

        cell_center_y(1);

        text = @sprintf(
            "%.3e",
            s.max_op_risk
        ),

        align = (:center, :center),
        fontsize = table_font,
    )


    # System Risk label

    text!(
        ax_t,

        risk_x,

        cell_center_y(2);

        text = "System Risk",
        align = (:center, :center),
        fontsize = table_font,
        color = :steelblue,
        font = :bold,
    )


    # System Risk value

    text!(
        ax_t,

        risk_x,

        cell_center_y(3);

        text = @sprintf(
            "%.3e",
            s.system_risk
        ),

        align = (:center, :center),
        fontsize = table_font,
    )


    # ================================================================
    # NC legend
    # ================================================================

    Legend(
        bottom[1, 2],

        [
            PolyElement(
                color = :lightgrey
            )
        ],

        [
            "NC — no cascade and no load shed"
        ];

        orientation = :vertical,

        tellwidth = false,
        tellheight = false,

        valign = :center,
        halign = :right,

        labelsize = legend_font,

        framevisible = false,
    )


    # ================================================================
    # Layout
    # ================================================================

    rowsize!(
        fig.layout,
        2,
        250
    )


    rowgap!(
        fig.layout,
        1,
        12
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

    mkpath(
        outdir
    )


    for data in (
        FIG6_3,
        FIG6_4,
        FIG6_5,
        FIG6_6
    )

        outpath =
            joinpath(
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