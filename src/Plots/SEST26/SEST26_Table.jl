using DataFrames

Pathplot=Glob["Path"]*"../Result/Iter_Report/toplot/" 
initp=[3,4,5,6]
method_code=[3]
if Test_Sys=="Polish/"
    obj_digits=1e6
    unit="M"
else
    obj_digits=1e3
    unit="K"
end
LSB_R=read_result(Pathplot,"SCOPF")
CALSB_R=read_result(Pathplot,"CA_SCOPF")
NLP_R=read_NLresult(Pathplot,"SCOPF")
NLP_R=NLP_R[(3,3)]

LSB_CR=read_result(Pathplot*"Cascade_report/","SCOPF","CR")
CALSB_CR=read_result(Pathplot*"Cascade_report/","CA_SCOPF","CR")
NLP_CR=read_NLresult(Pathplot*"Cascade_report/","SCOPF","CR")
NLP_CR=NLP_CR[(3,3)]


result0=[LSB_R, CALSB_R]
result0_C=[LSB_CR, CALSB_CR]
method=["SCOPF (LSB)","CA-SCOPF (LSB)"]
data = DataFrame(
    Initial = ["Mild",["" for _ in 1:length(method)]..., "Moderate",["" for _ in 1:length(method)]..., "Severe",["" for _ in 1:length(method)]...,"Extreme",["" for _ in 1:length(method)]...],
    Method=repeat(["Initial point","SCOPF (LSB)","CA-SCOPF (LSB)"], outer=length(initp)),
    # Time = ["--", summary[rowsum[3,3],"Res_Time"], summary[rowsum[3,1],"Res_Time"], summary[rowsum[3,5],"Res_Time"],
    #         "--", summary[rowsum[4,3],"Res_Time"], summary[rowsum[4,1],"Res_Time"], summary[rowsum[4,5],"Res_Time"],
    #         "--", summary[rowsum[5,3],"Res_Time"], summary[rowsum[5,1],"Res_Time"], summary[rowsum[5,5],"Res_Time"],
    #         "--", summary[rowsum[6,3],"Res_Time"], summary[rowsum[6,1],"Res_Time"], summary[rowsum[6,5],"Res_Time"]],

    # Time=vcat([["--",[round(summary[rowsum[r,m],"Res_Time"], digits=2) for m in methods]...] for r in initp]...),
    Time=vcat([["--",[round(result0[m][(r,3)][end,"Res_Time"], digits=2) for m in 1:length(method)]...] for r in initp]...),
  
    # Iterations = [0,summary[rowsum[3,3],"Iter"],summary[rowsum[3,1],"Iter"],summary[rowsum[3,5],"Iter"], 0,summary[rowsum[4,3],"Iter"],summary[rowsum[4,1],"Iter"],summary[rowsum[4,5],"Iter"], 0,summary[rowsum[5,3],"Iter"],summary[rowsum[5,1],"Iter"],summary[rowsum[5,5],"Iter"],0,summary[rowsum[6,3],"Iter"],summary[rowsum[6,1],"Iter"],summary[rowsum[6,5],"Iter"]]
    # 
    # Iterations = vcat([[0,[summary[rowsum[r,m],"Iter"] for m in methods]...] for r in initp]...),
    Iterations = vcat([[0,[result0[m][(r,3)][end,"Iter"] for m in 1:length(method)]...] for r in initp]...),

    Objective = round.((vcat([[result0[1][(r,3)][1,"real_Obj"],[result0[m][(r,3)][end,"real_Obj"] for m in 1:length(method)]...] for r in initp]...))./obj_digits,digits=3),
    LNS_cas= round.((vcat([[result0[1][(r,3)][1,"LNS_cas"],[result0[m][(r,3)][end,"LNS_cas"] for m in 1:length(method)]...] for r in initp]...)),digits=2),
    OptGap= round.((vcat([[result0[1][(r,3)][1,"OptGap"],[result0[m][(r,3)][end,"OptGap"] for m in 1:length(method)]...] for r in initp]...)),digits=4),
    Max_Normal_V= round.((vcat([[max(result0[1][(r,3)][1,"Vmin_N"],result0[1][(r,3)][1,"Vmax_N"]),[ max(result0[m][(r,3)][end,"Vmin_N"],result0[m][(r,3)][end,"Vmax_N"]) for m in 1:length(method)]...] for r in initp]...)),digits=4),
    Max_Cont_V= round.((vcat([[max(result0[1][(r,3)][1,"Vmin_C"],result0[1][(r,3)][1,"Vmax_C"]),[ max(result0[m][(r,3)][end,"Vmin_C"],result0[m][(r,3)][end,"Vmax_C"]) for m in 1:length(method)]...] for r in initp]...)),digits=4),
    Nb_Normal_V= vcat([[(result0[1][(r,3)][1,"Vmin_N_N"]+result0[1][(r,3)][1,"Vmax_N_N"]),[(result0[m][(r,3)][end,"Vmin_N_N"]+result0[m][(r,3)][end,"Vmax_N_N"]) for m in 1:length(method)]...] for r in initp]...),
    Nb_Cont_V= vcat([[(result0[1][(r,3)][1,"Vmin_C_N"]+result0[1][(r,3)][1,"Vmax_C_N"]),[(result0[m][(r,3)][end,"Vmin_C_N"]+result0[m][(r,3)][end,"Vmax_C_N"]) for m in 1:length(method)]...] for r in initp]...),
    Mean_Normal_V= round.((vcat([[ (result0[1][(r,3)][1,"Vmin_N_M"]*result0[1][(r,3)][1,"Vmin_N_N"]+result0[1][(r,3)][1,"Vmax_N_M"]*result0[1][(r,3)][1,"Vmax_N_N"])/(result0[1][(r,3)][1,"Vmin_N_N"]+result0[1][(r,3)][1,"Vmax_N_N"]),[ (result0[m][(r,3)][end,"Vmin_N_M"]*result0[m][(r,3)][end,"Vmin_N_N"]+result0[m][(r,3)][end,"Vmax_N_M"]*result0[m][(r,3)][end,"Vmax_N_N"])/(result0[m][(r,3)][end,"Vmin_N_N"]+result0[m][(r,3)][end,"Vmax_N_N"]) for m in 1:length(method)]...] for r in initp]...)),digits=4),
    Mean_Cont_V= round.((vcat([[ (result0[1][(r,3)][1,"Vmin_C_M"]*result0[1][(r,3)][1,"Vmin_C_N"]+result0[1][(r,3)][1,"Vmax_C_M"]*result0[1][(r,3)][1,"Vmax_C_N"])/(result0[1][(r,3)][1,"Vmin_C_N"]+result0[1][(r,3)][1,"Vmax_C_N"]),[ (result0[m][(r,3)][end,"Vmin_C_M"]*result0[m][(r,3)][end,"Vmin_C_N"]+result0[m][(r,3)][end,"Vmax_C_M"]*result0[m][(r,3)][end,"Vmax_C_N"])/(result0[m][(r,3)][end,"Vmin_C_N"]+result0[m][(r,3)][end,"Vmax_C_N"]) for m in 1:length(method)]...] for r in initp]...)),digits=4),
    
    Max_Normal_I= round.((vcat([[result0[1][(r,3)][1,"Smax_N"],[result0[m][(r,3)][end,"Smax_N"] for m in 1:length(method)]...] for r in initp]...)).*100,digits=2),
    Max_Cont_I= round.((vcat([[result0[1][(r,3)][1,"Smax_C"],[result0[m][(r,3)][end,"Smax_C"] for m in 1:length(method)]...] for r in initp]...)).*100,digits=2),
    Nb_Normal_I= vcat([[(result0[1][(r,3)][1,"Smax_N_N"]),[(result0[m][(r,3)][end,"Smax_N_N"]) for m in 1:length(method)]...] for r in initp]...),
    Nb_Cont_I= vcat([[(result0[1][(r,3)][1,"Smax_C_N"]),[(result0[m][(r,3)][end,"Smax_C_N"]) for m in 1:length(method)]...] for r in initp]...),
    Mean_Normal_I= round.((vcat([[result0[1][(r,3)][1,"Smax_N_M"],[result0[m][(r,3)][end,"Smax_N_M"] for m in 1:length(method)]...] for r in initp]...)).*100,digits=2),
    Mean_Cont_I= round.((vcat([[result0[1][(r,3)][1,"Smax_C_M"],[result0[m][(r,3)][end,"Smax_C_M"] for m in 1:length(method)]...] for r in initp]...)).*100,digits=2),

    Max_Cas_I=round.((vcat([[result0_C[1][(r,3)][1,"Worst_Cas_MW"],[result0_C[m][(r,3)][end,"Worst_Cas_MW"] for m in 1:length(method)]...] for r in initp]...)).*100,digits=2),
    Nb_Cas_I=vcat([[(result0_C[1][(r,3)][1,"CasNum"]),[(result0_C[m][(r,3)][end,"CasNum"]) for m in 1:length(method)]...] for r in initp]...),
    Mean_Cas_I=round.((vcat([[result0_C[1][(r,3)][1,"Mean_Cas_MW"],[result0_C[m][(r,3)][end,"Mean_Cas_MW"] for m in 1:length(method)]...] for r in initp]...)).*100,digits=2),
    Nbcont_Cas_I=vcat([[(result0_C[1][(r,3)][1,"CasNum"]),[(result0_C[m][(r,3)][end,"CasNum"]) for m in 1:length(method)]...] for r in initp]...),
)

for (i,m) in enumerate(initp)
    # insert!(data,2+(i-1)*(length(methods)+2),["" "exact NLP"  round(summary[1,"Time"],digits=1) "--" round(summary[1,"Obj"]/obj_digits,digits=2) Int(round.(100 .*nonlinear_redispatch[initp[i]][1], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][2], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][3], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][4], digits=0)) "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--"],promote=true)
    insert!(data,2+(i-1)*(length(method)+2),["" "SCOPF (NLP)"  round(NLP_R[1,"Time"],digits=1) "--" round(NLP_R[1,"Optimizer_Obj"]/obj_digits,digits=3) round(NLP_R[1,"LNS_cas"],digits=2) round(NLP_R[1,"OptGap"],digits=4) 0 0 0 0 "--" "--" 0 0 0 0 "--" "--" round(100 .*NLP_CR[1,"Worst_Cas_MW"],digits=2) NLP_CR[1,"CasNum"] round(100 .*NLP_CR[1,"Mean_Cas_MW"],digits=2) NLP_CR[1,"CasNum"]],promote=true)
end


for c in names(data)
    data[!, c] = [ismissing(x) ? "0.0" : x for x in data[!, c]]
    data[!, c] = [(x isa AbstractFloat && isnan(x)) ? "-" : x for x in data[!, c]]
end

latex_table = raw"""
\begin{table*}[th!]
\centering
\caption{Performance of Methods to solve SCOPF problem on the modified Nordic 60-bus System}
\label{tab:nordic_results}
\resizebox{2\columnwidth}{!}{%
\begin{tabular}{cccccccccccccccccccc}
\toprule
\multirow{2}{*}{Stress} & \multirow{3}{*}{Method} & \multirow{3}{*}{Time} & \multirow{2}{*}{Number of} & \multirow{2}{*}{Dispatch} 
& \multicolumn{6}{c}{Voltage Violations ($p.u.$)} 
& \multicolumn{9}{c}{Line Overload (MVA)} \\ 
\cmidrule(lr){6-11} \cmidrule(lr){12-20}
level & & & Iterations & cost 
& \multicolumn{3}{c}{Normal state} & \multicolumn{3}{c}{Contingency states}
& \multicolumn{3}{c}{Normal state} & \multicolumn{3}{c}{Contingency states} & \multicolumn{3}{c}{Cascading} \\
\cmidrule(lr){6-8} \cmidrule(lr){9-11} \cmidrule(lr){12-14} \cmidrule(lr){15-17} \cmidrule(lr){18-20}
& & (s) & & (K\euro{}) & Max & Mean & Nb & Max & Mean & Nb & Max & Mean & Nb & Max & Mean & Nb & Max & Mean & Nb \\
\midrule
"""
# - & exact non-linear & $(round(summary[1,"Time"],digits=1)) & - & $(round(summary[1,"Obj"]/obj_digits,digits=2)) & -- & -- & - & - & - & - & - & - & - & - & - & - & - & - \\\\

for i in 1:nrow(data)
    row = data[i, :]
    # Light gray for initial points
    color = (row.Method == "Initial point") ? "\\rowcolor{lightgray} " : ""
    # star=(row.Method == "Initial point") ? "^{*}" : ""
    # phantom=(row.Method == "Initial point") ? "\\phantom{*}" : ""
    star=""
    phantom=""
    latex_table *= "$color $(row.Initial) & $(row.Method) & $(row.Time) & $(row.Iterations) & $(row.Objective) & $(row.Max_Normal_V) & $(row.Mean_Normal_V) & $(row.Nb_Normal_V) & $(row.Max_Cont_V) & $(row.Mean_Cont_V) & $(row.Nb_Cont_V) & $(row.Max_Normal_I) & $(row.Mean_Normal_I) & $(row.Nb_Normal_I) & $(row.Max_Cont_I) & $(row.Mean_Cont_I) & $(row.Nb_Cont_I) & $(row.Max_Cas_I) & $(row.Mean_Cas_I)  & $(row.Nb_Cas_I) \\\\\n" #& $(row.Nbcont_Cas_I)

    # print("$(row.Normal)$star")
end

# End of table
latex_table *= """
\\bottomrule
\\end{tabular}
}
\\end{table*}
"""

println(latex_table)