# Recovered directly from thesis Figure 6.5.
# Matrix orientation:
#   rows    = contingencies 1:33
#   columns = time periods 1:24
# NaN denotes an "NC" cell (no cascading propagation and no load shedding).

const FIG6_5 = (
    figure_id = "6.5",
    title = "QSS Simulation — Simplified Severity by Time Period and Contingency, iter=5",
    suggested_filename = "RBSCOPF.pdf",

    # Common color scale used by the four Chapter 6 figures.
    color_max = 214.0,

    severity = [
        1.21  0.20  0.80  1.76  3.52  NaN  3.51  5.75  4.82  3.43  3.44  6.80  1.03  2.41  0.09  0.36  NaN  NaN  NaN  0.01  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  0.49  1.93  2.35  1.93  1.21  1.21  4.53  0.00  0.00  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  1.21  1.21  1.21  1.21  1.21  1.61  2.01  1.61  1.21  1.21  2.47  1.21  1.21  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  1.21  1.21  1.21  1.21  1.21  2.47  2.47  2.47  1.21  1.61  4.64  1.21  1.21  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.47  NaN  NaN  NaN  NaN  1.21  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  10.38  20.48  19.91  16.60  1.21  1.21  20.10  7.33  NaN  7.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  5.98  1.61  1.21  1.21  6.05  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  9.84  2.07  1.61  1.21  1.21  21.23  NaN  NaN  NaN  NaN  NaN  NaN  4.10  1.21  NaN  NaN  1.01  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        2.71  NaN  1.21  1.21  1.21  NaN  3.11  3.57  3.11  1.21  1.21  3.57  1.50  NaN  NaN  NaN  1.21  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  21.52  3.64  19.00  1.21  1.21  3.15  20.02  20.45  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  7.91  3.10  6.07  1.21  1.21  2.73  3.62  3.76  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  3.02  2.61  1.21  1.21  3.48  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.61  3.48  3.02  1.21  1.21  3.48  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  3.57  NaN  NaN  NaN  NaN  NaN  1.50  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
    ],

    summary = (
        voltage_pre  = (max = 0.000, mean = NaN, n = 0),
        voltage_post = (max = 0.000, mean = NaN, n = 0),
        line_pre     = (max = 0.000, mean = 0.000, n = 0),
        line_post    = (max = 0.000, mean = 0.000, n = 0),

        extra_label   = "inacc",
        extra_voltage = nothing,
        extra_line    = 0.0,

        max_op_risk = 21.517,
        system_risk = 1.017
    )
)
