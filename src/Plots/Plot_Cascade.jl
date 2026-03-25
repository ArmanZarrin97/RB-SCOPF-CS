using LaTeXStrings
using DataFrames

# Define your data as a DataFrame
Pathplot=Glob["Path"]*"../Result/Iter_Report/toplot/Cascade/" 
if Test_Sys=="Polish/"
    obj_digits=1e6
    unit="M"
else
    obj_digits=1e3
    unit="K"
end
# datadict=Dict()
# rowsum=Dict()
method_code=[3]
initp=[3,4,5,6]
# summary=ods_read(Pathplot*"LIN_Sum_SCOPF",sheetName="General", retType= "DataFrame")[1:end-1,:]
function read_result(path)
    # path=Pathplot*"CA_LSB/"
    datadict=Dict()
    for i in initp, m in method_code
        try
            raw=read_iter_rep(path,Problem_Type,(i,m))
            fi=findfirst(x->x=="end",raw[:,1])
            datadict[i,m]=coalesce.(raw[2:fi-1,:],0)
            # datadict[i,m]=read_iter_rep(Problem_Type,(i,m))[2:end-1,:]
        catch
            continue
        end
    end
    return datadict
end
function read_NLresult(path)
    # path=Pathplot*"CA_LSB/"
    datadict=Dict()
    for i in initp, m in method_code
        try
            raw=read_iter_rep(path,Problem_Type,(i,m))
            # fi=findfirst(x->x=="end",raw[:,1])
            datadict[i,m]=coalesce.(raw[1:1,:],0)
            # datadict[i,m]=read_iter_rep(Problem_Type,(i,m))[2:end-1,:]
        catch
            continue
        end
    end
    return datadict
end

LSB_result=read_result(Pathplot*"LSB/")
CA_LSB_result=read_result(Pathplot*"CA_LSB/")
RB_LSB_result=read_result(Pathplot*"RB_LSB/")
result0=[LSB_result,CA_LSB_result,RB_LSB_result]
NL_result=read_NLresult(Pathplot*"LSB/")[(3,3)]
methods=["LSB","CA_LSB","RB_LSB"]
# \$({CA}_{LSB})\$
nTP=nTL

data = DataFrame(
    Initial = ["Mild",["" for _ in 1:length(methods)]..., "Moderate",["" for _ in 1:length(methods)]..., "Severe",["" for _ in 1:length(methods)]...,"Extreme",["" for _ in 1:length(methods)]...],
    Method=repeat(["Initial point", "LSB","CA-LSB","RB-LSB"], outer=length(initp)),
    # Time = ["--", summary[rowsum[3,3],"Res_Time"], summary[rowsum[3,1],"Res_Time"], summary[rowsum[3,5],"Res_Time"],
    #         "--", summary[rowsum[4,3],"Res_Time"], summary[rowsum[4,1],"Res_Time"], summary[rowsum[4,5],"Res_Time"],
    #         "--", summary[rowsum[5,3],"Res_Time"], summary[rowsum[5,1],"Res_Time"], summary[rowsum[5,5],"Res_Time"],
    #         "--", summary[rowsum[6,3],"Res_Time"], summary[rowsum[6,1],"Res_Time"], summary[rowsum[6,5],"Res_Time"]],

    # Time=vcat([["--",[round(summary[rowsum[r,m],"Res_Time"], digits=2) for m in methods]...] for r in initp]...),
    Time=vcat([["--",[round(result0[m][(r,3)][end,"Res_Time"], digits=2) for m in 1:length(methods)]...] for r in initp]...),
  
    # Iterations = [0,summary[rowsum[3,3],"Iter"],summary[rowsum[3,1],"Iter"],summary[rowsum[3,5],"Iter"], 0,summary[rowsum[4,3],"Iter"],summary[rowsum[4,1],"Iter"],summary[rowsum[4,5],"Iter"], 0,summary[rowsum[5,3],"Iter"],summary[rowsum[5,1],"Iter"],summary[rowsum[5,5],"Iter"],0,summary[rowsum[6,3],"Iter"],summary[rowsum[6,1],"Iter"],summary[rowsum[6,5],"Iter"]]
    # 
    # Iterations = vcat([[0,[summary[rowsum[r,m],"Iter"] for m in methods]...] for r in initp]...),
    Iterations = vcat([[0,[result0[m][(r,3)][end,"Iter"] for m in 1:length(methods)]...] for r in initp]...),

    Objective = round.((vcat([[result0[m][(r,3)][1,"real_Obj"],[result0[m][(r,3)][end,"real_Obj"] for m in 1:length(methods)]...] for r in initp]...))./obj_digits,digits=3),
    LNS_cas= round.((vcat([[result0[m][(r,3)][1,"LNS_cas"],[result0[m][(r,3)][end,"LNS_cas"] for m in 1:length(methods)]...] for r in initp]...)),digits=2),
    OptGap= round.((vcat([[result0[m][(r,3)][1,"OptGap"],[result0[m][(r,3)][end,"OptGap"] for m in 1:length(methods)]...] for r in initp]...)),digits=4),
    # Objective = [datadict[3,3][1,"real_Obj"],summary[rowsum[3,3],"Obj"],summary[rowsum[3,1],"Obj"],summary[rowsum[3,5],"Obj"],
    #              datadict[4,3][1,"real_Obj"],summary[rowsum[4,3],"Obj"],summary[rowsum[4,1],"Obj"],summary[rowsum[4,5],"Obj"],
    #              datadict[5,3][1,"real_Obj"],summary[rowsum[5,3],"Obj"],summary[rowsum[5,1],"Obj"],summary[rowsum[5,5],"Obj"]],
    # Normal=vcat([["--",[Int(round.(100 .*summary[rowsum[r,m],"Max_N_redis"], digits=0)) for m in methods]...] for r in initp]...),
    # Contingency=vcat([["--",[Int(round.(100 .*summary[rowsum[r,m],"Max_C_redis"], digits=0)) for m in methods]...] for r in initp]...),
    # NormalQ=vcat([["--",[Int(round.(100 .*summary[rowsum[r,m],"Max_N_redisQ"], digits=0)) for m in methods]...] for r in initp]...),
    # ContingencyQ=vcat([["--",[Int(round.(100 .*summary[rowsum[r,m],"Max_C_redisQ"], digits=0)) for m in methods]...] for r in initp]...),
    # # Active = ["--",0.0,1.2,0.8,"--",0.0,1.1,0.9,"--",0.0,1.0,0.7],
    # # Reactive = ["--",0.0,"Feas.","Feas.","Infeas.",0.0,"Feas.","Feas.","Infeas.",0.0,"Feas.","Feas."],
    # Max_Normal_Line=round.(100 .*vcat([[datadict[r,methods[1]][1,"Smax_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_N"] for m in methods]...] for r in initp]...),digits=1),
    # Mean_Normal_Line=round.(100 .*vcat([[datadict[r,methods[1]][1,"Smax_N_M"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_N_M"] for m in methods]...] for r in initp]...),digits=1),
    # Nb_Normal_Line=vcat([[datadict[r,methods[1]][1,"Smax_N_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_N_N"] for m in methods]...] for r in initp]...),
    # # Max_Normal_Line = ["9.6","-","8.1","9.4","18.9","-","8.1","9.4","16.1","-","8.2","9.4"],
    # # Mean_Normal_Line = ["3.8","-","4.0","4.3","5.7","-","3.9","4.3","5.8","-","4.1","4.3"],
    # # Nb_Normal_Line = ["19","0","29","29","26","0","29","29","18","0","28","29"],
    # Max_Cont_Line=round.(100 .*vcat([[datadict[r,methods[1]][1,"Smax_C"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_C"] for m in methods]...] for r in initp]...),digits=1),
    # Mean_Cont_Line=round.(100 .*vcat([[datadict[r,methods[1]][1,"Smax_C_M"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_C_M"] for m in methods]...] for r in initp]...),digits=1),
    # Nb_Cont_Line=vcat([[datadict[r,methods[1]][1,"Smax_C_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_C_N"] for m in methods]...] for r in initp]...),
    # # Max_Cont_Line = ["32.6","0.2","50.4","45.5","60.5","0.2","55.0","45.5","140.2","1.3","47.5","45.5"],
    # # Mean_Cont_Line = ["4.6","0.2","4.6","5.0","7.5","0.2","4.8","5.0","8.3","0.6","4.7","5.0"],
    # # Nb_Cont_Line = ["609","2","869","753","670","2","782","753","560","10","809","753"],
    # Max_Normal_V=round.(vcat([[max(datadict[r,methods[1]][1,"Vmax_N"],datadict[r,methods[1]][1,"Vmin_N"]),[max(datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N"],datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N"]) for m in methods]...] for r in initp]...),digits=3),
    # Mean_Normal_V=round.(vcat([[ (datadict[r,methods[1]][1,"Vmax_N_M"]*datadict[r,methods[1]][1,"Vmax_N_N"]+datadict[r,methods[1]][1,"Vmin_N_M"]*datadict[r,methods[1]][1,"Vmin_N_N"])/(datadict[r,methods[1]][1,"Vmax_N_N"]+datadict[r,methods[1]][1,"Vmin_N_N"]),[ (datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N_M"]*datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N_M"]*datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N_N"])/(datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N_N"]) for m in methods]...] for r in initp]...),digits=3),
    # Nb_Normal_V=vcat([[datadict[r,methods[1]][1,"Vmax_N_N"]+datadict[r,methods[1]][1,"Vmin_N_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N_N"] for m in methods]...] for r in initp]...),
    # Max_Cont_V=round.(vcat([[max(datadict[r,methods[1]][1,"Vmax_C"],datadict[r,methods[1]][1,"Vmin_C"]),[max(datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C"],datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C"]) for m in methods]...] for r in initp]...),digits=3),
    # Mean_Cont_V=round.(vcat([[ (datadict[r,methods[1]][1,"Vmax_C_M"]*datadict[r,methods[1]][1,"Vmax_C_N"]+datadict[r,methods[1]][1,"Vmin_C_M"]*datadict[r,methods[1]][1,"Vmin_C_N"])/(datadict[r,methods[1]][1,"Vmax_C_N"]+datadict[r,methods[1]][1,"Vmin_C_N"]),[ (datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C_M"]*datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C_M"]*datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C_N"])/(datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C_N"]) for m in methods]...] for r in initp]...),digits=3),
    # Nb_Cont_V=vcat([[datadict[r,methods[1]][1,"Vmax_C_N"]+datadict[r,methods[1]][1,"Vmin_C_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C_N"] for m in methods]...] for r in initp]...)
   
    # # Max_Normal_V = ["0.003","-","0.005","0.005","0.004","-","0.005","0.005","0.007","-","0.005","0.005"],
    # # Mean_Normal_V = ["0.002","-","0.002","0.003","0.002","-","0.002","0.003","0.003","-","0.003","0.003"],
    # # Nb_Normal_V = ["15","0","55","65","25","0","52","65","60","0","65","65"],
    # Max_Cont_V = ["0.006","0.002","0.006","0.006","0.013","0.002","0.006","0.006","0.017","0.015","0.005","0.006"],
    # Mean_Cont_V = ["0.003","0.002","0.002","0.002","0.004","0.002","0.002","0.002","0.005","0.009","0.003","0.002"],
    # Nb_Cont_V = ["318","7","300","169","599","9","235","169","1546","18","134","169"]
)
for (i,m) in enumerate(initp)
    # insert!(data,2+(i-1)*(length(methods)+2),["" "exact NLP"  round(summary[1,"Time"],digits=1) "--" round(summary[1,"Obj"]/obj_digits,digits=2) Int(round.(100 .*nonlinear_redispatch[initp[i]][1], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][2], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][3], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][4], digits=0)) "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--"],promote=true)
    insert!(data,2+(i-1)*(length(methods)+2),["" "exact NLP"  round(NL_result[1,"Time"],digits=1) "--" round(NL_result[1,"Optimizer_Obj"]/obj_digits,digits=3) round(NL_result[1,"LNS_cas"],digits=2) round(NL_result[1,"OptGap"],digits=4) ],promote=true)

end
    # insert!(data,2+,["" "exact NLP"  round(summary[1,"Time"],digits=1) "--" round(summary[1,"Obj"]/obj_digits,digits=2) Int(round.(100 .*nonlinear_redispatch[initp[2]][1], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[2]][2], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[2]][3], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[2]][4], digits=0)) "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--"],promote=true)
    # insert!(data,2,["" "exact NLP"  round(summary[1,"Time"],digits=1) "--" round(summary[1,"Obj"]/obj_digits,digits=2) Int(round.(100 .*nonlinear_redispatch[initp[3]][1], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[3]][2], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[3]][3], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[3]][4], digits=0)) "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--"],promote=true)
    # insert!(data,2,["" "exact NLP"  round(summary[1,"Time"],digits=1) "--" round(summary[1,"Obj"]/obj_digits,digits=2) Int(round.(100 .*nonlinear_redispatch[initp[4]][1], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[4]][2], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[4]][3], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[4]][4], digits=0)) "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--"],promote=true)

for c in names(data)
    data[!, c] = [ismissing(x) ? "0.0" : x for x in data[!, c]]
    data[!, c] = [(x isa AbstractFloat && isnan(x)) ? "-" : x for x in data[!, c]]
end

# Begin LaTeX table
latex_table = """
\\begin{table}[th!]
\\centering
\\caption{Performance of Methods to solve the redispatch problem with 24 look-ahead time periods on the New England 39-bus System}
\\label{tab:nordic_results_short}
\\resizebox{1\\columnwidth}{!}{%
\\begin{tabular}{ccccccc}
\\toprule
Stress level & Method & Time (s) & Iterations & Dispatch cost (K€) & Deviation (\\%) & \${LNS}_{Cas}\$ \\\\
\\midrule
"""

# - & exact non-linear & $(round(summary[1,"Time"],digits=1)) & - & $(round(summary[1,"Obj"]/obj_digits,digits=2)) & -- & -- & - & - & - & - & - & - & - & - & - & - & - & - \\\\
# \\midrule

# \\multicolumn{2}{c}{Active Power(MW)} & \\multicolumn{2}{c}{Reactive Power(MVA)}

# Function to add rows
for i in 1:nrow(data)
    row = data[i, :]
    # Light gray for initial points
    color = (row.Method == "Initial point") ? "\\rowcolor{lightgray} " : ""
    # star=(row.Method == "Initial point") ? "^{*}" : ""
    # phantom=(row.Method == "Initial point") ? "\\phantom{*}" : ""
    star=""
    phantom=""
    latex_table *= "$color $(row.Initial) & $(row.Method) & $(row.Time) & $(row.Iterations) & $(row.Objective) & + $(row.OptGap) & $(row.LNS_cas) \\\\\n"

    # print("$(row.Normal)$star")
end

# End of table
latex_table *= """
\\bottomrule
\\end{tabular}
}
\\end{table}
"""

println(latex_table)