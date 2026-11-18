# Recovered directly from thesis Figure 6.4.
# Matrix orientation:
#   rows    = contingencies 1:33
#   columns = time periods 1:24
# NaN denotes an "NC" cell (no cascading propagation and no load shedding).

const FIG6_4 = (
    figure_id = "6.4",
    title = "QSS Simulation — Simplified Severity by Time Period and Contingency, iter=7",
    suggested_filename = "CASCOPF.pdf",

    # Common color scale used by the four Chapter 6 figures.
    color_max = 214.0,

    severity = [
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  1.21  1.21  1.21  1.61  1.21  2.01  2.47  2.01  1.61  1.61  2.47  1.21  1.21  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  1.21  1.21  1.21  1.61  1.21  2.47  4.64  2.47  1.61  1.61  4.64  1.21  1.21  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  2.07  2.47  2.07  1.21  1.21  2.47  1.21  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.07  1.21  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.07  1.21  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.07  1.21  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        2.71  NaN  1.21  1.21  1.21  NaN  3.57  2.07  2.07  1.21  1.21  3.57  2.71  1.50  NaN  NaN  NaN  1.50  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.07  1.21  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.07  1.21  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  2.61  NaN  3.02  3.48  3.02  2.61  2.61  3.48  1.41  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.07  1.21  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  1.21  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.07  1.21  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  2.61  NaN  3.48  3.48  3.48  2.61  2.61  3.48  1.21  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  3.57  3.11  2.71  2.71  3.57  NaN  NaN  NaN  NaN  NaN  1.50  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
    ],

    summary = (
        voltage_pre  = (max = 0.000, mean = NaN, n = 0),
        voltage_post = (max = 0.000, mean = NaN, n = 0),
        line_pre     = (max = 0.000, mean = 0.000, n = 0),
        line_post    = (max = 0.000, mean = 0.000, n = 0),

        extra_label   = "inacc",
        extra_voltage = nothing,
        extra_line    = 0.0,

        max_op_risk = 4.639,
        system_risk = 0.7174
    )
)
