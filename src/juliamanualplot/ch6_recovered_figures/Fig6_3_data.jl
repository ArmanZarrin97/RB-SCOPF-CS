# Recovered from _6, 3_.pdf, page 6.
# Selected figure: SCOPF
#
# Matrix convention:
#   severity[c, t]
#   c = contingency 1:33
#   t = time period 1:24
#   NaN = "NC" (no cascade and no load shed)
#
# The common color scale is retained at 214.0 p.u., matching the source figures.

const FIG6_3 = (
    figure_id = "6.3",
    title = "SCOPF — Simulation-Based Cascading Expected Severity",
    suggested_filename = "SCOPF.pdf",
    color_max = 214.0,

    severity = [
        1.28  1.79  0.82  1.81  3.32  3.36  4.15  4.39  4.31  2.40  2.65  7.14  3.26  3.27  3.07  1.62  1.21  1.21  0.31  0.21  NaN  NaN  0.04  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  4.37  3.19  1.61  1.57  NaN  1.21  3.19  1.15  1.11  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  1.21  1.21  1.21  1.21  1.21  1.61  1.61  1.21  1.21  1.21  1.61  1.21  1.21  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  1.21  1.21  1.21  1.21  1.21  1.61  1.61  1.61  1.21  1.21  1.61  1.21  1.21  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  1.21  1.21;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  NaN  1.21  1.21  76.98  83.29  88.89  88.47  1.21  1.21  116.84  101.09  100.06  40.78  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  57.38  63.52  70.57  70.46  1.21  1.21  100.69  82.38  81.07  22.72  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  3.51  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        13.77  NaN  14.37  10.21  1.21  35.36  139.54  92.82  50.37  48.78  1.21  209.09  60.07  57.58  8.88  34.37  1.21  1.21  30.91  6.62  NaN  NaN  22.39  1.98;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  NaN  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  1.21  1.21  1.21  1.21  1.21  1.21  NaN  NaN  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  1.21  1.21  1.21  NaN  1.21  2.71  1.21  1.21  1.21  2.71  NaN  NaN  1.21  NaN  1.21  1.21  NaN  1.21  NaN  NaN  1.21  NaN;
        1.21  NaN  NaN  1.21  1.21  1.16  93.50  45.68  116.54  1.21  1.21  21.01  33.62  86.68  1.86  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  72.98  39.32  98.67  1.21  1.21  18.45  28.91  73.54  NaN  NaN  1.21  1.21  NaN  1.21  NaN  NaN  NaN  NaN;
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
        voltage_pre  = (max = 0, mean = NaN, n = 0),
        voltage_post = (max = 0, mean = NaN, n = 0),
        line_pre     = (max = 0, mean = 0, n = 0),
        line_post    = (max = 0, mean = 0, n = 0),

        extra_label   = nothing,
        extra_voltage = nothing,
        extra_line    = nothing,

        # In the simplified plotting function these are displayed as:
        #   Maximum simplified Sev
        #   Overall system risk
        max_op_risk = 209.1,
        system_risk = 4.477
    )
)
