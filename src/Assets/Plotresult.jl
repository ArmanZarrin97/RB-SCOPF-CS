# # For plotting Sensitivities
#     fig = Figure(resolution=(900, 600))

#     Busnumber=800
#     Linenumber=2227
#     # Primary axis (for Max Violation)
#     ax1 = Axis(fig[1, 1],
#             xlabel="index of Bus", ylabel="Value of Sensitivity",
#             title="Sensitivity of Voltage to 1p.u. active power injections",
#             backgroundcolor=:white,
#             xticklabelsize=12, yticklabelsize=12,
#             xgridstyle=:dash, ygridstyle=:dash,
#             xgridcolor=:gray, ygridcolor=:gray)

#     ax2 = Axis(fig[1, 2],
#             xlabel="index of Bus", ylabel="Value of Sensitivity",
#             title="Sensitivity of Voltage to 1p.u. Reactive power injections",
#             backgroundcolor=:white,
#             xticklabelsize=12, yticklabelsize=12,
#             xgridstyle=:dash, ygridstyle=:dash,
#             xgridcolor=:gray, ygridcolor=:gray)

#     ax3 = Axis(fig[2, 1],
#             xlabel="index of Bus", ylabel="Value of Sensitivity",
#             title="Sensitivity of Line Flow to 1p.u. Active power injections",
#             backgroundcolor=:white,
#             xticklabelsize=12, yticklabelsize=12,
#             xgridstyle=:dash, ygridstyle=:dash,
#             xgridcolor=:gray, ygridcolor=:gray)

#     ax4 = Axis(fig[2, 2],
#             xlabel="index of Bus", ylabel="Value of Sensitivity",
#             title="Sensitivity of Line Flow to 1p.u. Reactive power injections",
#             backgroundcolor=:white,
#             xticklabelsize=12, yticklabelsize=12,
#             xgridstyle=:dash, ygridstyle=:dash,
#             xgridcolor=:gray, ygridcolor=:gray)


#     line1 = scatterlines!(ax1, S_VP[1,1][Busnumber,:], label="Max Violation (Method 1)", color=:blue, alpha=0.3, markersize=0.5)
#     line2 = scatterlines!(ax2, S_VQ[1,1][Busnumber,:], label="Max Violation (Method 1)", color=:orange, alpha=0.3,markersize=0.5)

#     line1 = scatterlines!(ax3, S_IP[1,1][Linenumber,:], label="Max Violation (Method 1)", color=:blue, alpha=0.3, markersize=0.5)
#     line2 = scatterlines!(ax4, S_IQ[1,1][Linenumber,:], label="Max Violation (Method 1)", color=:orange, alpha=0.3,markersize=0.5)


#     S_VP[1,1][600,:]/100
#     vlines!(ax, t, color=:gray, linestyle=:dash, linewidth=1, alpha=0.3)
#     vlines!(ax, t, color=:gray, linestyle=:dash, linewidth=1, alpha=0.3,markersize=0.5)


#     save("Sensistivities.png", fig)


# # Plot at-risk constraints

# fig = Figure(resolution=(900, 600))
# # Primary axis (for Max Violation)
# ax1 = Axis(fig[1:7, 1:3],
#         xlabel="percent of constraints", ylabel="potential violation",
#         title="Constraint Filtering trade-off",
#         backgroundcolor=:white,
#         xticklabelsize=12, yticklabelsize=12,
#         xgridstyle=:dash, ygridstyle=:dash,
#         xgridcolor=:gray, ygridcolor=:gray)


# line1 = scatterlines!(ax1, report_TR[3,"plotdata"][1], label=string(TR), alpha=0.3, markersize=0.5)
# remaining_vio1=round(report_TR[3,"y_acc_Q"],digits=5)
# TR1= round(report_TR[3,"TR"],digits=4)

# line2 = scatterlines!(ax1, report_TR[11,"plotdata"], label=string(TR),color=:green, alpha=0.3, markersize=0.5)
# remaining_vio2=round(report_TR[11,"y_acc"],digits=5)
# TR2= round(report_TR[11,"TR"],digits=4)

# line3 = scatterlines!(ax1, report_TR[2,"plotdata"], label=string(TR), alpha=0.3, markersize=0.5)
# remaining_vio3=round(report_TR[2,"y_acc"],digits=5)
# TR3= round(report_TR[2,"TR"],digits=4)

# vlines!(ax1, zmax, color=:red, linestyle=:dot,alpha=0.9)


Legend(fig[8, 2],
[line1,line2,line3],
["TR= $TR1 and remaining violation is $remaining_vio1",
"TR= $TR2 and remaining violation is $remaining_vio2",
"TR= $TR3 and remaining violation is $remaining_vio3"]) #patchsize = (10, 10), rowgap = 3


save("Potential_violation.png", fig)


# plot results
methods_=[1,3,5]
StartP=6
global Pathplot=Path*"../Result/Iter_Report/toplot/" 
for StartP in [ 3,4,5,6]
        code_=[StartP,methods_]
        fig=plot_Iter_rep2(Problem_Type,code_;scale_="on")
        save(Pathplot*Problem_Type*string(code_)*"scaled.png", fig)
end

# for Runcode in [(3,3) (4,3) (5,3) (6,3)]
#     fig=plot_Iter_rep(read_iter_rep(Problem_Type,Runcode))
#     save(Path*"../Result/Iter_Report/"*Problem_Type*string(Runcode)*".png", fig)
# end
