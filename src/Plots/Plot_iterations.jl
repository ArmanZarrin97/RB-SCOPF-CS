include("../Funcs/IO/Packages.jl") # load necessary packages
using CairoMakie


global Pathplot=Path*"../Result/Iter_Report/toplot/" 
datadict=Dict()
rowsum=Dict()
# initp=[6],methods=[3]

runcodes=["(6, 3)","(6, 3)_A","(6, 3)_B","(6, 3)_C"]



for rr in runcodes
    try
        raw=read_iter_rep(Problem_Type,rr)
        fi=findfirst(x->x=="end",raw[:,1])
        datadict[rr]=coalesce.(raw[2:fi-1,:],0)
        # datadict[i,m]=read_iter_rep(Problem_Type,(i,m))[2:end-1,:]
    catch
        continue
    end
end

# Example time vector (10 points, e.g., hours)
# time = 1:10  # replace with actual time values, e.g., 0:1:9 or DateTime array

methods = ["LSB", "V_I", "V_II","V_III"]
colors = [:green, :red, :blue, :orange]

# # Random example values
# voltage_normal = [rand(10) for _ in 1:3]
line=Vector{Any}(undef, 4);voltage=Vector{Any}(undef, 4);voltage_normal=Vector{Any}(undef, 4);voltage_contingency=Vector{Any}(undef, 4);

for i in 1:4
    raw=datadict[runcodes[i]]
    voltage_normal[i]=mapslices(maximum, Matrix(raw[:,["Vmin_Abs","Vmax_Abs"]]), dims=2)[:,1]
    voltage_contingency[i]=mapslices(maximum, Matrix(raw[:,["Vmin_C_Abs","Vmax_C_Abs"]]), dims=2)[:,1]
    voltage[i]=max.(voltage_normal[i],voltage_contingency[i])
    line[i]=max.(raw[:,"Smax_Abs"],raw[:,"Smax_C_Abs"])
end
line*=100
# voltage=max.(voltage_normal,voltage_contingency)
# voltage_contingency = [rand(10) for _ in 1:3]
# line_normal = [rand(10) for _ in 1:3]
# line_contingency = [rand(10) for _ in 1:3]

# Create figure with 2 subplots (vertical layout)
fig = Figure(resolution = (600, 600),fontsize=16)

ax1 = Axis(fig[1, 1], title = "(a)", xlabel = "Time (seconds)", ylabel = "Voltage violation (p.u.)")
ax2 = Axis(fig[2, 1], title = "(b)", xlabel = "Time (seconds)", ylabel = "Line flow violation (MW)")

legend_lines = []
legend_labels = []
linestylevec = [:dash, :dash, :dash, :dash]
markersvec = [:circle, :rect, :cross, :diamond]
# Plot voltage violations
for i in 1:4
    time=datadict[runcodes[i]][:,"Res_Time"]
    l1=scatterlines!(ax1, time, voltage[i], color = colors[i], label = "$(methods[i])", linestyle = linestylevec[i], marker= markersvec[i])

    # l1=lines!(ax1, time, voltage_normal[i], color = colors[i], label = "$(methods[i]) Normal")
    # l2=lines!(ax1, time, voltage_contingency[i], color = colors[i], linestyle = :dash, label = "$(methods[i]) Contingency")
    # CairoMakie.scatter!(ax1, time, voltage_normal[i], color = colors[i])
    # CairoMakie.scatter!(ax1, time, voltage_contingency[i], color = colors[i])

    # push!(legend_lines, l1)
    # push!(legend_labels, "$(methods[i]) Normal")
    # push!(legend_lines, l2)
    # push!(legend_labels, "$(methods[i]) Contingency")
end

# Plot line flow violations
for i in 1:4
    scatterlines!(ax2, datadict[runcodes[i]][:,"Res_Time"], line[i], color = colors[i], label = "$(methods[i])", linestyle = linestylevec[i], marker= markersvec[i])

    # scatterlines!(ax2, time, line_normal[i], color = colors[i], label = "$(methods[i]) Normal")
    # scatterlines!(ax2, time, line_contingency[i], color = colors[i], linestyle = :dash, label = "$(methods[i]) Contingency")
    # CairoMakie.scatter!(ax2, time, line_normal[i], color = colors[i])
    # CairoMakie.scatter!(ax2, time, line_contingency[i], color = colors[i])
end

# # Add legends
# axislegend(ax1, position = :rt,nbanks=2)
# axislegend(ax2, position = :rt,fontsize=12,nbanks=2)

axislegend(ax1, position = :rt,orientation = :horizontal)
axislegend(ax2, position = :rt,fontsize=12,orientation= :horizontal)
# legend = Legend(fig, legend_lines, legend_labels, orientation = :horizontal)
# fig[3, 1] = legend

fig


using Printf

using Printf  # needed for @sprintf

# === Example data (replace with yours) ===
iterations = 1:10
methods = ["LSB", "V_I", "V_II","V_III"]

# TR  = [1 1 1; 1 2 1; 2 2 2; 2 2 2; 3 3 2; 3 3 3; 3 3 3; 2 2 3; 2 2 3; 1 1 2]
# Time = [5 4 6; 10 8 11; 18 15 20; 26 22 29; 35 31 38; 46 40 50; 58 51 63; 70 62 75; 85 75 90; 100 88 108]
# FI  = rand(10,3) .* 1e-2  # Feasibility index


TR=Vector{Any}(undef, 4);obj=Vector{Any}(undef, 4);N_I=Vector{Any}(undef, 4);N_V=Vector{Any}(undef, 4);

# Vector{Vector{Any}{zero,10}}(undef, 4)
# Vector{Any}(zeros, 4)
for i in 1:4
    raw=datadict[runcodes[i]]
    TR[i]=eval.(Meta.parse.(raw[:,"TR"]))
    obj[i]=raw[:,"real_Obj"]
    N_I[i]=eval.(Meta.parse.(raw[:,"Num_C_I"]))
    N_V[i]=eval.(Meta.parse.(raw[:,"Num_C_V"]))
end

# for i in 1:length(iterations)
#     tr_P=zeros(4,length(iterations));tr_Q=zeros(4,length(iterations))
#     for m in 1:4
#         tr_P[m,i]=TR[m][i][1]
#         tr_Q[m,i]=TR[m][i][2]

#         latex_table *= @sprintf(
#         "%d & %d & %.0f & %.2e & %d & %.0f & %.2e & %d & %.0f & %.2e \\\\\n",
#         iterations[i],
#         tr_P[m,i], obj[1][i], N_I[1][i][1],
#         TR[2][i][1], obj[2][i], N_I[2][i][1],
#         TR[3][i][1], obj[3][i], N_I[3][i][1],
#         TR[4][i][1], obj[4][i], N_I[4][i][1]
#     )
# end


if Test_Sys=="Polish/"
    obj_digits=1e6
    unit="M"
else
    obj_digits=1e3
    unit="K"
end

maxiter=7
dataframe2 = DataFrame(
    Iter=1:maxiter,
    trp=[vcat([ length(TR[m])>=i+1 ? round(TR[m][i+1][1]*100,digits=1) : "--" for m in 1:length(methods)]) for i=1:maxiter],
    trq=[vcat([ length(TR[m])>=i+1 ? round(TR[m][i+1][2]*100,digits=1) : "--" for m in 1:length(methods)]) for i=1:maxiter],
    obj=[vcat([ length(obj[m])>=i+1 ? round(obj[m][i+1]/obj_digits,digits=2) : "--" for m in 1:length(methods)]) for i=1:maxiter],
    volt=[vcat([ length(N_V[m])>=i+1 ? round((N_V[m][i+1][1]/N_V[m][i+1][2])*100,digits=2) : "--" for m in 1:length(methods)]) for i=1:maxiter],
    line=[vcat([ length(N_I[m])>=i+1 ? round((N_I[m][i+1][1]/N_I[m][i+1][2])*100,digits=2) : "--" for m in 1:length(methods)]) for i=1:maxiter]
)



latex_table = raw"""
\begin{table*}[ht]
\centering
\caption{Comparison of TRs, objective value, and retained constraints for four methods}
\resizebox{\textwidth}{!}{%
\begin{tabular}{c|cccc|cccc|cccc|cccc|cccc}
\toprule
\multirow{4}{*}{\textbf{Iter}} &
\multicolumn{8}{c|}{\textbf{TRs}} &
\multicolumn{4}{c|}{\multirow{2}{*}{\textbf{dispatch cost (M€)}}} &
\multicolumn{8}{c}{\textbf{Retained Constraints (\%)}} \\
\cmidrule(lr){2-9} \cmidrule(lr){14-21}
 & \multicolumn{4}{c|}{$\mathbf{tr_P}$ (MW)} & \multicolumn{4}{c|}{$\mathbf{tr_Q}$ (MVar)} & \multicolumn{4}{c|}{} & 
   \multicolumn{4}{c|}{Voltage} & \multicolumn{4}{c}{Line Flow} \\
\cmidrule(lr){2-5} \cmidrule(lr){6-9} \cmidrule(lr){10-13} \cmidrule(lr){14-17} \cmidrule(lr){18-21}"""
latex_table *="
 & $(methods[1]) & $(methods[2]) & $(methods[3]) & $(methods[4]) & $(methods[1]) & $(methods[2]) & $(methods[3]) & $(methods[4]) & $(methods[1]) & $(methods[2]) & $(methods[3]) & $(methods[4]) & $(methods[1]) & $(methods[2]) & $(methods[3]) & $(methods[4]) & $(methods[1]) & $(methods[2]) & $(methods[3]) & $(methods[4]) \\\\\n
\\midrule
"

for i in 1:nrow(dataframe2)
    row = dataframe2[i, :]
    latex_table *= "$(row.Iter) & $(row.trp[1]) & $(row.trp[2]) & $(row.trp[3]) & $(row.trp[4]) & $(row.trq[1]) & $(row.trq[2]) & $(row.trq[3]) & $(row.trq[4]) & $(row.obj[1]) & $(row.obj[2]) & $(row.obj[3]) & $(row.obj[4]) & $(row.volt[1]) & $(row.volt[2]) & $(row.volt[3]) & $(row.volt[4]) & $(row.line[1]) & $(row.line[2]) & $(row.line[3]) & $(row.line[4]) \\\\\n"
end



latex_table *= raw"""
\bottomrule
\end{tabular}}
\label{tab:CF-iterative}
\end{table*}
"""

println(latex_table)


