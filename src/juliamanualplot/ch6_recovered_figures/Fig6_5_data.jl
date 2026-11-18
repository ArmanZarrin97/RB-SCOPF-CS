# Recovered from RB_6, 3_.pdf, page 11.
# Selected figure: RBSCOPF
#
# Matrix convention:
#   severity[c, t]
#   c = contingency 1:33
#   t = time period 1:24
#   NaN = "NC" (no cascade and no load shed)
#
# The common color scale is retained at 214.0 p.u., matching the source figures.

const FIG6_5 = (
    figure_id = "6.5",
    title = "RBCA-SCOPF — Simulation-Based Cascading Expected Severity",
    suggested_filename = "RBSCOPF.pdf",
    color_max = 214.0,

    severity = [
        1.21  0.08  0.37  1.09  3.30  NaN  3.01  5.75  4.82  2.72  2.67  6.80  0.66  1.99  0.04  0.15  NaN  NaN  NaN  0.00  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  NaN  1.21  2.07  1.21  1.21  1.21  2.07  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  NaN  1.21  0.13  1.58  2.24  1.69  1.21  1.21  4.52  0.00  0.00  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN;
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
        1.21  NaN  NaN  1.21  1.21  NaN  7.83  3.64  6.64  1.21  1.21  3.15  4.61  4.79  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
        1.21  NaN  NaN  1.21  1.21  NaN  3.58  3.10  2.90  1.21  1.21  2.73  0.83  0.88  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
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
        voltage_pre  = (max = 0, mean = NaN, n = 0),
        voltage_post = (max = 0, mean = NaN, n = 0),
        line_pre     = (max = 0, mean = 0, n = 0),
        line_post    = (max = 0, mean = 0, n = 0),

        extra_label   = "inacc",
        extra_voltage = nothing,
        extra_line    = 0,

        # In the simplified plotting function these are displayed as:
        #   Maximum simplified Sev
        #   Overall system risk
        max_op_risk = 21.227,
        system_risk = 0.9145
    )
)
