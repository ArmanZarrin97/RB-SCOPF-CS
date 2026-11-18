# Recovered directly from thesis Figure 6.3.
# Matrix orientation:
#   rows    = contingencies 1:33
#   columns = time periods 1:24
# NaN denotes an "NC" cell (no cascading propagation and no load shedding).

const FIG6_3 = (
    figure_id = "6.3",
    title = "QSS Simulation — Simplified Severity by Time Period and Contingency, iter=5",
    suggested_filename = "SCOPF.pdf",

    # Common color scale used by the four Chapter 6 figures.
    color_max = 214.0,

    severity = [
        7.95  8.14  7.84  8.04  9.41  8.97  9.93  9.93  9.95  8.20  9.08  10.18  9.06  9.06  8.55  8.06  1.21  1.21  7.10  7.01  NaN  NaN  6.90  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  6.67  7.29  7.02  7.02  NaN  1.21  7.35  6.16  6.15  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  1.21  1.21  1.21  1.21  1.21  1.61  1.61  1.21  1.21  1.21  1.61  1.21  1.21  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  1.21  1.21  1.21  1.21  1.21  1.61  1.61  1.61  1.21  1.21  1.61  1.21  1.21  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  NaN  1.21  1.21  214.19  213.12  213.12  213.12  1.21  1.21  213.12  214.33  214.33  214.08  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  214.19  213.12  213.12  213.12  1.21  1.21  213.12  214.33  214.33  214.08  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  213.12  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        212.52  NaN  214.03  212.82  1.21  214.19  213.12  213.12  213.12  212.82  1.21  213.12  214.33  214.33  214.08  214.03  1.21  1.21  213.73  212.52  NaN  NaN  213.73  214.03;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  1.21  1.21  1.21  NaN  1.21  2.71  1.21  1.21  1.21  2.71  NaN  NaN  1.21  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  NaN  1.21  1.21  212.98  210.83  204.20  209.22  1.21  1.21  201.56  206.04  209.56  212.87  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  210.83  204.20  209.22  1.21  1.21  201.56  206.04  209.56  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  NaN  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  NaN  1.21  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        2.71  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN
    ],

    summary = (
        voltage_pre  = (max = 0.000, mean = NaN, n = 0),
        voltage_post = (max = 0.000, mean = NaN, n = 0),
        line_pre     = (max = 0.000, mean = 0.000, n = 0),
        line_post    = (max = 0.000, mean = 0.000, n = 0),

        extra_label   = nothing,
        extra_voltage = nothing,
        extra_line    = nothing,

        max_op_risk = 214.33,
        system_risk = 5.619
    )
)
