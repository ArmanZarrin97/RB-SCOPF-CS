using DataFrames

global Pathplot=Glob["Path"]*"../Result/Iter_Report/toplot/"
initp=[3,4,5,6]
method_code=[3]
if Test_Sys=="Polish/"
    obj_digits=1e6
    unit="M"
else
    obj_digits=1e3
    unit="K"
end
LSB_R=read_result(Pathplot,"_SCOPF")
CALSB_R=read_result(Pathplot,"CA_SCOPF")
RBLSB_R=read_result(Pathplot,"RB_SCOPF")
NLP_R=read_NLresult(Pathplot,"_SCOPF")
NLP_R=NLP_R[(3,3)]

LSB_CR=read_result(Pathplot*"Cascade_report/","_SCOPF","CR")
CALSB_CR=read_result(Pathplot*"Cascade_report/","CA_SCOPF","CR")
RBLSB_CR=read_result(Pathplot*"Cascade_report/","RB_SCOPF","CR")
NLP_CR=read_NLresult(Pathplot*"Cascade_report/","_SCOPF","CR")
NLP_CR=NLP_CR[(3,3)]


result0=[LSB_R, CALSB_R, RBLSB_R]
result0_C=[LSB_CR, CALSB_CR, RBLSB_CR]
method=["SCOPF (LSB)","CA-SCOPF (LSB)","RB-SCOPF (LSB)"]

data_2 = DataFrame(
    Initial = ["Mild",["" for _ in 1:length(method)]..., "Moderate",["" for _ in 1:length(method)]..., "Severe",["" for _ in 1:length(method)]...,"Extreme",["" for _ in 1:length(method)]...],
    Method=repeat(["Initial point","SCOPF","CA-SCOPF","RBCA-SCOPF"], outer=length(initp)),
    Risk= round.((vcat([[result0_C[1][(r,3)][1,"R"],[result0_C[m][(r,3)][end,"R"] for m in 1:length(method)]...] for r in initp]...)),digits=3),
    Risk_S= round.((vcat([[result0_C[1][(r,3)][1,"R_S"],[result0_C[m][(r,3)][end,"R_S"] for m in 1:length(method)]...] for r in initp]...)),digits=3),
    Risk_OV= round.((vcat([[result0_C[1][(r,3)][1,"R_0"],[result0_C[m][(r,3)][end,"R_0"] for m in 1:length(method)]...] for r in initp]...)),digits=3),

    Sevmax= round.((vcat([[result0_C[1][(r,3)][1,"ROP"],[result0_C[m][(r,3)][end,"ROP"] for m in 1:length(method)]...] for r in initp]...)),digits=3),
    Sevmax_S= round.((vcat([[result0_C[1][(r,3)][1,"ROP_S"],[result0_C[m][(r,3)][end,"ROP_S"] for m in 1:length(method)]...] for r in initp]...)),digits=3),
    Sevmax_OV= round.((vcat([[result0_C[1][(r,3)][1,"ROP_0"],[result0_C[m][(r,3)][end,"ROP_0"] for m in 1:length(method)]...] for r in initp]...)),digits=3),

    Maxrating= round.((vcat([[result0_C[1][(r,3)][1,"Worst_Cas_per"],[result0_C[m][(r,3)][end,"Worst_Cas_per"] for m in 1:length(method)]...] for r in initp]...))*100,digits=1),
    MaxMW= round.((vcat([[result0_C[1][(r,3)][1,"Worst_Cas_MW"],[result0_C[m][(r,3)][end,"Worst_Cas_MW"] for m in 1:length(method)]...] for r in initp]...))*100,digits=1),
)

for (i,m) in enumerate(initp)
    # insert!(data,2+(i-1)*(length(methods)+2),["" "exact NLP"  round(summary[1,"Time"],digits=1) "--" round(summary[1,"Obj"]/obj_digits,digits=2) Int(round.(100 .*nonlinear_redispatch[initp[i]][1], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][2], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][3], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][4], digits=0)) "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--"],promote=true)
    insert!(data_2,2+(i-1)*(length(method)+2),["" "SCOPF (NLP)"  round(NLP_CR[1,"R"],digits=3) round(NLP_CR[1,"R_S"],digits=3) round(NLP_CR[1,"R_0"],digits=3) round(NLP_CR[1,"ROP"],digits=3) round(NLP_CR[1,"ROP_S"],digits=3) round(NLP_CR[1,"ROP_0"],digits=3) round(NLP_CR[1,"Worst_Cas_per"]*100,digits=1) round(NLP_CR[1,"Worst_Cas_MW"]*100,digits=1)],promote=true)
end

for c in names(data_2)
    data_2[!, c] = [ismissing(x) ? "0.0" : x for x in data_2[!, c]]
    data_2[!, c] = [(x isa AbstractFloat && isnan(x)) ? "-" : x for x in data_2[!, c]]
end

latex_table = raw"""
\begin{table*}[th!]
\centering
\caption{Comparison of Cascaing Risk and maximum Severity for different methods.}
\label{tab:risk_results}
\resizebox{\columnwidth}{!}{%
\begin{tabular}{ccccccccccccccccccccc}
\toprule
\multirow{2}{*}{Stress} & \multirow{3}{*}{Method} & \multicolumn{3}{c}{Risk} & \multicolumn{3}{c}{Severity Metrics} & \multicolumn{2}{c}{Maximum Cascading Overload} \\
\cmidrule(lr){3-5} \cmidrule(lr){6-8} \cmidrule(lr){9-10}
level & & Risk & Risk^{0} & Risk^{OV}  & {Sev} & {Sev}^{0} & {Sev}^{OV} & {rating} & {amount} \\
& & ($p.u.$) & ($p.u.$) & - & ($p.u.$) & ($p.u.$) & - & (\%) & (MVA)\\
\midrule
"""
# - & exact non-linear & $(round(summary[1,"Time"],digits=1)) & - & $(round(summary[1,"Obj"]/obj_digits,digits=2)) & -- & -- & - & - & - & - & - & - & - & - & - & - & - & - \\\\

for i in 1:nrow(data_2)
    row = data_2[i, :]
    # Light gray for initial points
    color = (row.Method == "Initial point") ? "\\rowcolor{lightgray} " : ""
    # star=(row.Method == "Initial point") ? "^{*}" : ""
    # phantom=(row.Method == "Initial point") ? "\\phantom{*}" : ""
    star=""
    phantom=""
    latex_table *= "$color $(row.Initial) & $(row.Method) & $(row.Risk) & $(row.Risk_S) & $(row.Risk_OV) & $(row.Sevmax) & $(row.Sevmax_S) & $(row.Sevmax_OV) & $(row.Maxrating) & $(row.MaxMW) \\\\\n" #& $(row.Nbcont_Cas_I)

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