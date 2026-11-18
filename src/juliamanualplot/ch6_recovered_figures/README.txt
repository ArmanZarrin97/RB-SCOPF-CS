Recovered Chapter 6 figure data
================================

Files:
  Fig6_3_data.jl  - thesis Figure 6.3 (SCOPF)
  Fig6_4_data.jl  - thesis Figure 6.4 (CA-SCOPF)
  Fig6_5_data.jl  - thesis Figure 6.5 (RBCA-SCOPF)
  Fig6_6_data.jl  - thesis Figure 6.6 (RBCA-SCOPF optimization estimate)
  plot_ch6_recovered.jl - CairoMakie plotting function

Matrix convention:
  severity[c,t]
  c = contingency 1:33
  t = time period 1:24
  NaN = "NC"

Quick use in Julia:

  include("Fig6_3_data.jl")
  include("Fig6_4_data.jl")
  include("Fig6_5_data.jl")
  include("Fig6_6_data.jl")
  include("plot_ch6_recovered.jl")

  plot_all_ch6_figures()

Or one figure:

  plot_ch6_figure(
      FIG6_5;
      outpath="RBSCOPF_large.pdf",
      size=(2400, 1500),
      cell_font=22
  )

Sanity checks against the thesis:
  Figure 6.3 max = 214.33 at contingencies 13, 15, 18 and periods 13, 14.
  Figure 6.4 max = 4.64 at contingency 6 and periods 8, 12.
  Figure 6.5 max = 21.52 at contingency 22 and period 7.
  Figure 6.6 max = 50.00.
