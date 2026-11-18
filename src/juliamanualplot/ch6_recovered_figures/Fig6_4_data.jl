# Recovered from CA_6, 3_.pdf, page 8.
# Selected figure: CASCOPF
#
# Matrix convention:
#   severity[c, t]
#   c = contingency 1:33
#   t = time period 1:24
#   NaN = "NC" (no cascade and no load shed)
#
# The common color scale is retained at 214.0 p.u., matching the source figures.

const FIG6_4 = (
    figure_id = "6.4",
    title = "CA-SCOPF — Simulation-Based Cascading Expected Severity",
    suggested_filename = "CASCOPF.pdf",
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
        2.71  NaN  1.21  1.21  1.21  NaN  3.57  2.07  2.07  1.21  1.21  3.57  2.71  1.50  NaN  NaN  NaN  NaN  NaN  1.21  NaN  NaN  NaN  NaN;
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
        max_op_risk = 4.639,
        system_risk = 0.717
    )
)
