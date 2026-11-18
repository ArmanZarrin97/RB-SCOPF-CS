# Recovered directly from thesis Figure 6.6.
# Matrix orientation:
#   rows    = contingencies 1:33
#   columns = time periods 1:24
# NaN denotes an "NC" cell (no cascading propagation and no load shedding).

const FIG6_6 = (
    figure_id = "6.6",
    title = "RBCA-SCOPF — Estimated Simplified Severity by Time Period and Contingency, iter=5",
    suggested_filename = "RBSCOPF_O.pdf",

    # Common color scale used by the four Chapter 6 figures.
    color_max = 214.0,

    severity = [
        1.21  0.65  1.96  3.33  5.11  NaN  5.17  7.02  6.24  4.95  4.96  7.68  2.41  4.07  0.34  1.04  NaN  NaN  NaN  0.08  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  1.25  2.77  2.78  2.99  1.21  1.21  5.48  0.04  0.04  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  1.21  1.21  1.21  1.21  1.21  2.01  2.47  1.61  1.21  1.21  2.47  1.21  1.21  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  1.21  1.21  1.21  1.21  1.21  2.47  2.47  2.47  1.21  1.61  4.64  1.21  1.21  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  2.07  2.07  2.07  1.21  1.21  2.47  NaN  NaN  NaN  NaN  1.21  NaN  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  29.84  42.63  48.85  44.72  1.21  1.21  48.12  22.74  NaN  21.68  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  16.61  1.61  1.21  1.21  15.59  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  29.09  2.07  1.61  1.21  1.21  50.00  NaN  NaN  0.00  NaN  NaN  NaN  13.67  1.21  NaN  NaN  4.51  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.61  2.07  1.61  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        2.71  NaN  1.21  1.21  1.21  NaN  3.11  3.57  3.11  1.21  1.21  3.57  1.50  NaN  NaN  NaN  1.21  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  50.00  5.03  46.54  1.21  1.21  4.28  50.00  50.00  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  20.91  4.34  16.25  1.21  1.21  3.65  12.44  12.90  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  3.02  2.61  1.21  1.21  3.48  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.61  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
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

        extra_label   = "TR size",
        extra_voltage = 0.016,
        extra_line    = 0.025,

        max_op_risk = 49.998,
        system_risk = 0.9992
    )
)
