using LaTeXStrings
using DataFrames

# Define your data as a DataFrame
global Pathplot=Path*"../Result/Iter_Report/toplot/" 
if Test_Sys=="Polish/"
    obj_digits=1e6
    unit="M"
else
    obj_digits=1e3
    unit="K"
end
datadict=Dict()
rowsum=Dict()
methods=[1,3,5]
initp=[3,4,5,6]
summary=ods_read(Pathplot*"LIN_Sum_SCOPF",sheetName="General", retType= "DataFrame")[1:end-1,:]
for i in initp, m in methods
    try
        raw=read_iter_rep(Problem_Type,(i,m))
        fi=findfirst(x->x=="end",raw[:,1])
        datadict[i,m]=coalesce.(raw[2:fi-1,:],0)
        # datadict[i,m]=read_iter_rep(Problem_Type,(i,m))[2:end-1,:]
    catch
        continue
    end
    try
        rowsum[i,m]=findall(x->x=="($(i), $(m))",summary[!,"Runcode"])[1]
    catch
        continue
    end
end

datadict[3,3][1,"real_Obj"]
# summary[rowsum[],"Res_Time"]
nTP=nTL
nonlinear_redispatch = Dict(
    r => find_total_redispatch!(OP_NLP, Init_OP[Linearization_INIT[r]]) 
    for r in initp
)
data = DataFrame(
    Initial = ["Mild",["" for _ in 1:length(methods)]..., "Moderate",["" for _ in 1:length(methods)]..., "Severe",["" for _ in 1:length(methods)]...,"Extreme",["" for _ in 1:length(methods)]...],
    Method=repeat(["Initial point", "\\cite{yangMWonly}", "LSB", "Hybrid LSB"], outer=length(initp)),
    # Time = ["--", summary[rowsum[3,3],"Res_Time"], summary[rowsum[3,1],"Res_Time"], summary[rowsum[3,5],"Res_Time"],
    #         "--", summary[rowsum[4,3],"Res_Time"], summary[rowsum[4,1],"Res_Time"], summary[rowsum[4,5],"Res_Time"],
    #         "--", summary[rowsum[5,3],"Res_Time"], summary[rowsum[5,1],"Res_Time"], summary[rowsum[5,5],"Res_Time"],
    #         "--", summary[rowsum[6,3],"Res_Time"], summary[rowsum[6,1],"Res_Time"], summary[rowsum[6,5],"Res_Time"]],

    Time=vcat([["--",[round(summary[rowsum[r,m],"Res_Time"], digits=2) for m in methods]...] for r in initp]...),

    # Iterations = [0,summary[rowsum[3,3],"Iter"],summary[rowsum[3,1],"Iter"],summary[rowsum[3,5],"Iter"], 0,summary[rowsum[4,3],"Iter"],summary[rowsum[4,1],"Iter"],summary[rowsum[4,5],"Iter"], 0,summary[rowsum[5,3],"Iter"],summary[rowsum[5,1],"Iter"],summary[rowsum[5,5],"Iter"],0,summary[rowsum[6,3],"Iter"],summary[rowsum[6,1],"Iter"],summary[rowsum[6,5],"Iter"]]
    # 
    Iterations = vcat([[0,[summary[rowsum[r,m],"Iter"] for m in methods]...] for r in initp]...),

    Objective = round.((vcat([[datadict[r,methods[1]][1,"real_Obj"],[summary[rowsum[r,m],"Obj"] for m in methods]...] for r in initp]...))./obj_digits,digits=2),
    # Objective = [datadict[3,3][1,"real_Obj"],summary[rowsum[3,3],"Obj"],summary[rowsum[3,1],"Obj"],summary[rowsum[3,5],"Obj"],
    #              datadict[4,3][1,"real_Obj"],summary[rowsum[4,3],"Obj"],summary[rowsum[4,1],"Obj"],summary[rowsum[4,5],"Obj"],
    #              datadict[5,3][1,"real_Obj"],summary[rowsum[5,3],"Obj"],summary[rowsum[5,1],"Obj"],summary[rowsum[5,5],"Obj"]],
    Normal=vcat([["--",[Int(round.(100 .*summary[rowsum[r,m],"Max_N_redis"], digits=0)) for m in methods]...] for r in initp]...),
    Contingency=vcat([["--",[Int(round.(100 .*summary[rowsum[r,m],"Max_C_redis"], digits=0)) for m in methods]...] for r in initp]...),
    NormalQ=vcat([["--",[Int(round.(100 .*summary[rowsum[r,m],"Max_N_redisQ"], digits=0)) for m in methods]...] for r in initp]...),
    ContingencyQ=vcat([["--",[Int(round.(100 .*summary[rowsum[r,m],"Max_C_redisQ"], digits=0)) for m in methods]...] for r in initp]...),
    # Active = ["--",0.0,1.2,0.8,"--",0.0,1.1,0.9,"--",0.0,1.0,0.7],
    # Reactive = ["--",0.0,"Feas.","Feas.","Infeas.",0.0,"Feas.","Feas.","Infeas.",0.0,"Feas.","Feas."],
    Max_Normal_Line=round.(100 .*vcat([[datadict[r,methods[1]][1,"Smax_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_N"] for m in methods]...] for r in initp]...),digits=1),
    Mean_Normal_Line=round.(100 .*vcat([[datadict[r,methods[1]][1,"Smax_N_M"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_N_M"] for m in methods]...] for r in initp]...),digits=1),
    Nb_Normal_Line=vcat([[datadict[r,methods[1]][1,"Smax_N_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_N_N"] for m in methods]...] for r in initp]...),
    # Max_Normal_Line = ["9.6","-","8.1","9.4","18.9","-","8.1","9.4","16.1","-","8.2","9.4"],
    # Mean_Normal_Line = ["3.8","-","4.0","4.3","5.7","-","3.9","4.3","5.8","-","4.1","4.3"],
    # Nb_Normal_Line = ["19","0","29","29","26","0","29","29","18","0","28","29"],
    Max_Cont_Line=round.(100 .*vcat([[datadict[r,methods[1]][1,"Smax_C"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_C"] for m in methods]...] for r in initp]...),digits=1),
    Mean_Cont_Line=round.(100 .*vcat([[datadict[r,methods[1]][1,"Smax_C_M"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_C_M"] for m in methods]...] for r in initp]...),digits=1),
    Nb_Cont_Line=vcat([[datadict[r,methods[1]][1,"Smax_C_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Smax_C_N"] for m in methods]...] for r in initp]...),
    # Max_Cont_Line = ["32.6","0.2","50.4","45.5","60.5","0.2","55.0","45.5","140.2","1.3","47.5","45.5"],
    # Mean_Cont_Line = ["4.6","0.2","4.6","5.0","7.5","0.2","4.8","5.0","8.3","0.6","4.7","5.0"],
    # Nb_Cont_Line = ["609","2","869","753","670","2","782","753","560","10","809","753"],
    Max_Normal_V=round.(vcat([[max(datadict[r,methods[1]][1,"Vmax_N"],datadict[r,methods[1]][1,"Vmin_N"]),[max(datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N"],datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N"]) for m in methods]...] for r in initp]...),digits=3),
    Mean_Normal_V=round.(vcat([[ (datadict[r,methods[1]][1,"Vmax_N_M"]*datadict[r,methods[1]][1,"Vmax_N_N"]+datadict[r,methods[1]][1,"Vmin_N_M"]*datadict[r,methods[1]][1,"Vmin_N_N"])/(datadict[r,methods[1]][1,"Vmax_N_N"]+datadict[r,methods[1]][1,"Vmin_N_N"]),[ (datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N_M"]*datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N_M"]*datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N_N"])/(datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N_N"]) for m in methods]...] for r in initp]...),digits=3),
    Nb_Normal_V=vcat([[datadict[r,methods[1]][1,"Vmax_N_N"]+datadict[r,methods[1]][1,"Vmin_N_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_N_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_N_N"] for m in methods]...] for r in initp]...),
    Max_Cont_V=round.(vcat([[max(datadict[r,methods[1]][1,"Vmax_C"],datadict[r,methods[1]][1,"Vmin_C"]),[max(datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C"],datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C"]) for m in methods]...] for r in initp]...),digits=3),
    Mean_Cont_V=round.(vcat([[ (datadict[r,methods[1]][1,"Vmax_C_M"]*datadict[r,methods[1]][1,"Vmax_C_N"]+datadict[r,methods[1]][1,"Vmin_C_M"]*datadict[r,methods[1]][1,"Vmin_C_N"])/(datadict[r,methods[1]][1,"Vmax_C_N"]+datadict[r,methods[1]][1,"Vmin_C_N"]),[ (datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C_M"]*datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C_M"]*datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C_N"])/(datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C_N"]) for m in methods]...] for r in initp]...),digits=3),
    Nb_Cont_V=vcat([[datadict[r,methods[1]][1,"Vmax_C_N"]+datadict[r,methods[1]][1,"Vmin_C_N"],[datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmax_C_N"]+datadict[r,m][summary[rowsum[r,m],"Iter"]+1,"Vmin_C_N"] for m in methods]...] for r in initp]...)
   
    # # Max_Normal_V = ["0.003","-","0.005","0.005","0.004","-","0.005","0.005","0.007","-","0.005","0.005"],
    # # Mean_Normal_V = ["0.002","-","0.002","0.003","0.002","-","0.002","0.003","0.003","-","0.003","0.003"],
    # # Nb_Normal_V = ["15","0","55","65","25","0","52","65","60","0","65","65"],
    # Max_Cont_V = ["0.006","0.002","0.006","0.006","0.013","0.002","0.006","0.006","0.017","0.015","0.005","0.006"],
    # Mean_Cont_V = ["0.003","0.002","0.002","0.002","0.004","0.002","0.002","0.002","0.005","0.009","0.003","0.002"],
    # Nb_Cont_V = ["318","7","300","169","599","9","235","169","1546","18","134","169"]
)
for (i,m) in enumerate(initp)
    insert!(data,2+(i-1)*(length(methods)+2),["" "exact NLP"  round(summary[1,"Time"],digits=1) "--" round(summary[1,"Obj"]/obj_digits,digits=2) Int(round.(100 .*nonlinear_redispatch[initp[i]][1], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][2], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][3], digits=0)) Int(round.(100 .*nonlinear_redispatch[initp[i]][4], digits=0)) "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--"],promote=true)
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
\\begin{table*}[th!]
\\centering
\\caption{Performance of Methods on the Polish System using Different Initial set-points}
\\label{tab:nordic_results}
\\resizebox{2\\columnwidth}{!}{%
\\begin{tabular}{ccccccccccccccccccccc}
\\toprule
\\multirow{2}{*}{Required} & \\multirow{3}{*}{Method} & \\multirow{3}{*}{Time} & \\multirow{2}{*}{Number of} & \\multirow{2}{*}{Objective} 
& \\multicolumn{4}{c}{Overall Re-dispatch} 
& \\multicolumn{6}{c}{Line Overload (MVA)} & \\multicolumn{6}{c}{Voltage Violations (\$p.u.\$)} \\\\ 
\\cmidrule(lr){6-9} \\cmidrule(lr){10-15} \\cmidrule(lr){16-21}
redispatch & & & Iterations & value &  \\multicolumn{2}{c}{Normal state} & \\multicolumn{2}{c}{Contingency states} & \\multicolumn{3}{c}{Normal state} & \\multicolumn{3}{c}{Contingency states} 
& \\multicolumn{3}{c}{Normal state} & \\multicolumn{3}{c}{Contingency states} \\\\
\\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-12} \\cmidrule(lr){13-15} \\cmidrule(lr){16-18} \\cmidrule(lr){19-21}
 & & (s) & & ($(unit)€) & (MW) & (MVar) & (MW) & (MVar) & Max & Mean & Nb & Max & Mean & Nb & Max & Mean & Nb & Max & Mean & Nb \\\\
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
    latex_table *= "$color $(row.Initial) & $(row.Method) & $(row.Time) & $(row.Iterations) & $(row.Objective) & $phantom$(row.Normal)$star & $phantom$(row.NormalQ)$star & $phantom$(row.Contingency)$star & $phantom$(row.ContingencyQ)$star & $(row.Max_Normal_Line) & $(row.Mean_Normal_Line) & $(row.Nb_Normal_Line) & $(row.Max_Cont_Line) & $(row.Mean_Cont_Line) & $(row.Nb_Cont_Line) & $(row.Max_Normal_V) & $(row.Mean_Normal_V) & $(row.Nb_Normal_V) & $(row.Max_Cont_V) & $(row.Mean_Cont_V) & $(row.Nb_Cont_V) \\\\\n"

    # print("$(row.Normal)$star")
end

# End of table
latex_table *= """
\\bottomrule
\\end{tabular}
}
\\end{table*}
"""

# \\\\[0.5mm]
# \\multicolumn{19}{l}{\\footnotesize{\\textbf{Note:} "Max" and "Mean" denote the worst and average violation magnitudes, respectively, and "Nb" indicates the number of violations.}}\\\\
# \\multicolumn{19}{l}{\\footnotesize{\\phantom{Note: }* In the Initial point rows, redispatch values correspond to the non-linear solution.}}

# Print LaTeX code
println(latex_table)




# ddd_extreme=[(1716, 198848)	(17575, 351488)	(122304, 1684032)	(13183304, 21524160)	;
# (1502, 198848)	(17575, 351488)	(114102, 1298304)	(9742721, 16254882)	;
# (785, 198848)	(17575, 351488)	(214110, 623070)	(9982689, 16550478)	;
# (329, 198848)	(17575, 351488)	(134055, 259854)	(10639438, 16714698)	;
# (343, 198848)	(17575, 351488)	(126546, 260820)	(10821914, 15848196)	;
# (283, 198848)	(17575, 351488)	(100490, 234030)	(10135374, 14157510)	;
# (396, 198848)	(1660, 351488)	(157818, 270758)	(951392, 1187736)	;
# (303, 198848)	(929, 351488)	(127198, 208980)	(505515, 681210)	;
# (420, 198848)	(961, 351488)	(174578, 280376)	(517887, 702152)	;
# (299, 198848)	(17575, 351488)	(102643, 202808)	(9359879, 12780944)	]

# ddd_severe= [(1773, 198848)	(17575, 351488)	(148532, 1344200)	(8691023, 16961736)	;
# (432, 198848)	(17575, 351488)	(160361, 351094)	(10003783, 16521538)	;
# (362, 198848)	(17575, 351488)	(116163, 299222)	(9363670, 15566110)	;
# (362, 198848)	(15050, 351488)	(111379, 270756)	(7776764, 11521620)	;
# (392, 198848)	(1602, 351488)	(123448, 253216)	(776227, 1119400)	;
# (379, 198848)	(1025, 351488)	(132455, 234850)	(489823, 733040)	
# ]

# ddd_moderate=[(1819, 198848)	(17575, 351488)	(142711, 921328)	(5738387, 11609160)	;
# (834, 198848)	(10841, 351488)	(125500, 388732)	(3723404, 5742432)	;
# (403, 198848)	(1949, 351488)	(103803, 204462)	(698695, 1013100)	;
# (390, 198848)	(4124, 351488)	(84282, 203184)	(1503449, 2101608)	;
# (407, 198848)	(1034, 351488)	(91506, 210368)	(346505, 584896)	
# ]

# ddd_mild=[(352, 198848)	(17575, 351488)	(71443, 153170)	(5745689, 7727930)	;
# (505, 198848)	(3527, 351488)	(98886, 204050)	(1043295, 1670030)	;
# (439, 198848)	(2961, 351488)	(79471, 188150)	(830102, 1331890)	;
# (485, 198848)	(4271, 351488)	(109468, 209880)	(1303517, 2112050)	
# ]


ddd_extreme=[(246, 4896)	(563, 11232)	(3430, 12366)	(19758, 29160);
(190, 4896)	(563, 11232)	(5010, 9396)	(21594, 28134);
(190, 4896)	(563, 11232)	(4919, 9018)	(19883, 29646);
(201, 4896)	(563, 11232)	(5197, 9936)	(21105, 28566);
(195, 4896)	(563, 11232)	(5579, 9450)	(19996, 29160);
(192, 4896)	(563, 11232)	(4913, 9072)	(21187, 28944);
(193, 4896)	(563, 11232)	(5990, 9450)	(19667, 29106);
(195, 4896)	(563, 11232)	(4328, 9234)	(26047, 28296);
(185, 4896)	(563, 11232)	(5372, 8910)	(22202, 28674);
(202, 4896)	(563, 11232)	(6339, 8856)	(20319, 29376);
(193, 4896)	(563, 11232)	(5079, 9288)	(22331, 30402);
(188, 4896)	(563, 11232)	(5451, 9126)	(20357, 30132);
(195, 4896)	(563, 11232)	(5339, 9396)	(19213, 29646);
(204, 4896)	(563, 11232)	(5598, 10368)	(19306, 29268);
(198, 4896)	(563, 11232)	(5176, 9396)	(19006, 29538);
]
ddd_severe= [(246, 4896)	(563, 11232)	(2654, 7412)	(14936, 18292);
(189, 4896)	(563, 11232)	(3497, 5344)	(13899, 16064);
(193, 4896)	(563, 11232)	(3596, 5344)	(14318, 17408);
(193, 4896)	(563, 11232)	(3070, 3912)	(10648, 13200);
(193, 4896)	(563, 11232)	(2913, 3648)	(11345, 12720);
(192, 4896)	(563, 11232)	(3067, 4032)	(11347, 12600);
(193, 4896)	(563, 11232)	(3122, 3718)	(9788, 12210);
(193, 4896)	(563, 11232)	(3036, 3608)	(9473, 11550);
(192, 4896)	(563, 11232)	(2754, 3234)	(9098, 10758);
(192, 4896)	(563, 11232)	(2833, 3454)	(9960, 11946);
(192, 4896)	(563, 11232)	(2689, 3278)	(9463, 11418);
(192, 4896)	(563, 11232)	(3090, 3500)	(9291, 11000);
(192, 4896)	(563, 11232)	(3146, 3540)	(9663, 10720);
(192, 4896)	(563, 11232)	(2928, 3200)	(9647, 11260);
(192, 4896)	(563, 11232)	(2918, 3380)	(9457, 10980)
]
ddd_moderate=[(246, 4896)	(563, 11232)	(2507, 5112)	(11986, 12960);
(192, 4896)	(563, 11232)	(2825, 3652)	(10635, 11704)
]
ddd_mild=[(207, 4896)	(563, 11232)	(2042, 3600)	(10220, 10480)	;
(192, 4896)	(7, 11232)	(1790, 3380)	(68, 80)	;
(192, 4896)	(563, 11232)	(2022, 3300)	(10185, 10380)	
]

using Statistics

datadict[(5,3)][2:end,end-3:end]
# Example data for four initial points

data_dict = Dict(
    "Mild" => datadict[(3,3)][2:end,end-3:end],
    "Moderate" => datadict[(4,3)][2:end,end-3:end],
    "Severe" => datadict[(5,3)][2:end,end-3:end],
    "Extreme" => datadict[(6,3)][2:end,end-3:end]
)

# Compute averages and problem size reduction
function compute_avg(ddd)
    # Filtered constraints
    filtered_line = mean([(1-ddd[i,1][1]/ddd[i,1][2])*100 for i in 1:size(ddd,1)])
    filtered_voltage = mean([(1-ddd[i,2][1]/ddd[i,2][2])*100 for i in 1:size(ddd,1)])
    filtered_overall = mean([(1-(ddd[i,1][1]+ddd[i,2][1]) / (ddd[i,1][2]+ddd[i,2][2])) * 100 for i in 1:size(ddd,1)])

    # Suppressed control actions
    suppressed_line = mean([(1-(ddd[i,3][1]/ddd[i,3][2]))*100 for i in 1:size(ddd,1)])
    suppressed_voltage = mean([(1-(ddd[i,4][1]/ddd[i,4][2]))*100 for i in 1:size(ddd,1)])
    suppressed_overall = mean([(1-(ddd[i,3][1]+ddd[i,4][1]) / (ddd[i,3][2]+ddd[i,4][2])) * 100 for i in 1:size(ddd,1)])

    # Problem size reduction as product of overall percentages
    ps_reduction = (1-((1-filtered_overall/100) * (1-suppressed_overall/100)))*100 # convert back to %

    return [filtered_line, filtered_voltage, filtered_overall,
            suppressed_line, suppressed_voltage, suppressed_overall,
            ps_reduction]
end

function compute_avg2(ddd)
    # ddd=data_dict["Mild"]
    # Filtered constraints
    # parse.(Int,ddd[1,1])
    # tuple(ddd[1,1])
    filtered_line_min = minimum([(ddd[i,1][1]/ddd[i,1][2])*100 for i in 1:size(ddd,1)])
    filtered_voltage_min = minimum([(ddd[i,2][1]/ddd[i,2][2])*100 for i in 1:size(ddd,1)])
    filtered_overall_min = minimum([((ddd[i,1][1]+ddd[i,2][1]) / (ddd[i,1][2]+ddd[i,2][2])) * 100 for i in 1:size(ddd,1)])

    filtered_line_max = maximum([(ddd[i,1][1]/ddd[i,1][2])*100 for i in 1:size(ddd,1)])
    filtered_voltage_max = maximum([(ddd[i,2][1]/ddd[i,2][2])*100 for i in 1:size(ddd,1)])
    filtered_overall_max = maximum([((ddd[i,1][1]+ddd[i,2][1]) / (ddd[i,1][2]+ddd[i,2][2])) * 100 for i in 1:size(ddd,1)])

    # Suppressed control actions
    # suppressed_line = mean([((ddd[i,3][1]/ddd[i,3][2]))*100 for i in 1:size(ddd,1)])
    # suppressed_voltage = mean([((ddd[i,4][1]/ddd[i,4][2]))*100 for i in 1:size(ddd,1)])
    # suppressed_overall = mean([((ddd[i,3][1]+ddd[i,4][1]) / (ddd[i,3][2]+ddd[i,4][2])) * 100 for i in 1:size(ddd,1)])

    # Problem size reduction as product of overall percentages
    ps_reduction_max = (100-filtered_overall_max) # convert back to %
    ps_reduction_min = (100-filtered_overall_min) # convert back to %

    return [(filtered_line_min,filtered_line_max), (filtered_voltage_min,filtered_voltage_max), (filtered_overall_min,filtered_overall_max),
            (ps_reduction_min,ps_reduction_max)]
end
# compute_avg2(data_dict["Mild"])

# eval.(Meta.parse.(ddd))

# parse.(Int, split.(strip.(ddd, ['(', ')']), ","))
# scenario_avgs = Dict(k => compute_avg(v) for (k,v) in data_dict)
scenario_avgs = Dict(k => compute_avg2(eval.(Meta.parse.(v))) for (k,v) in data_dict)



# Build LaTeX table string
latex_table = """
\\begin{table}[ht]
\\centering
\\caption{Average Problem Size Reduction over Iterations on the 2746-bus Test System (Proposed Method)}
\\label{tab:constraint_summary}
\\resizebox{\\columnwidth}{!}{%
\\begin{tabular}{|c|cc|c|}
\\hline
\\multirow{2}{*}{\\textbf{Initial Point}} & 
\\multicolumn{2}{c|}{\\begin{tabular}{c}\\textbf{Filtered}\\\\\\textbf{Constraints (\\%)}\\end{tabular}} & 
\\multirow{2}{*}{\\begin{tabular}{c}\\textbf{Problem Size}\\\\\\textbf{Reduction (\\%)}\\end{tabular}} \\\\
\\cline{2-3}
 & \\textbf{Line Flow} & \\textbf{Voltage} & \\\\
\\hline
"""

for label in ["Mild", "Moderate", "Severe", "Extreme"]
    avg = scenario_avgs[label]
    latex_table *= "$(label) & $(round(avg[1], digits=1)) & $(round(avg[2], digits=1)) & $(round(avg[7], digits=1)) \\\\\n"
end

latex_table *= """
\\hline
\\end{tabular}}
\\end{table}
"""

println(latex_table)



methods=[1,3,5]
initp=[3,4,5,6]
scenario_avgs=Dict()
Names=Dict("Polish/"=>"Polish 2746-bus","Nordic60/"=>"Nordic 60-bus","39bus/"=>"IEEE 39-bus")
for TS in ["Polish/","Nordic60/","39bus/"]
    # Path= "Test Systems/"*Test_Sys*"Data/"
    # TS="Nordic60/"
    global Pathplot="Test Systems/"*TS*"Data/"*"../Result/Iter_Report/toplot/" 
    datadict=Dict()
    summary=ods_read(Pathplot*"LIN_Sum_SCOPF",sheetName="General", retType= "DataFrame")[1:end-1,:]
    for i in initp, m in methods
        try
            raw=read_iter_rep(Problem_Type,(i,m))
            fi=findfirst(x->x=="end",raw[:,1])
            datadict[i,m]=coalesce.(raw[2:fi-1,:],0)
            # datadict[i,m]=read_iter_rep(Problem_Type,(i,m))[2:end-1,:]
        catch
            continue
        end
        # try
        #     rowsum[i,m]=findall(x->x=="($(i), $(m))",summary[!,"Runcode"])[1]
        # catch
        #     continue
        # end
    end
    data_dict2 = Dict(
    "Mild" => datadict[(3,3)][2:end,end-3:end],
    "Moderate" => datadict[(4,3)][2:end,end-3:end],
    "Severe" => datadict[(5,3)][2:end,end-3:end],
    "Extreme" => datadict[(6,3)][2:end,end-3:end]
    )
    v=data_dict2["Moderate"]
    scenario_avgs[TS] = Dict(k => compute_avg2(eval.(Meta.parse.(v))) for (k,v) in data_dict2)
end


latex_table = """
\\begin{table}[ht]
\\centering
\\caption{Average Problem Size Reduction over Iterations for Three Test Systems (Proposed Method)}
\\label{tab:constraint_summary_all}
\\resizebox{\\columnwidth}{!}{%
\\begin{tabular}{|c|cc|c|}
\\hline
\\multirow{2}{*}{\\textbf{Initial Point}} & 
\\multicolumn{2}{c|}{\\begin{tabular}{c}\\textbf{Retained}\\\\\\textbf{Constraints (\\%)}\\end{tabular}} & 
\\multirow{2}{*}{\\begin{tabular}{c}\\textbf{Problem Size}\\\\\\textbf{Reduction (\\%)}\\end{tabular}} \\\\
\\cline{2-3}
& \\textbf{Line Flow} & \\textbf{Voltage} & \\\\
\\hline
"""

for TS in ["39bus/","Nordic60/","Polish/"]
    latex_table *="\\multicolumn{4}{|c|}{\\textbf{$(Names[TS])}} \\\\
    \\hline\n"
    for label in ["Mild", "Moderate", "Severe", "Extreme"]
        avg = scenario_avgs[TS][label]
        latex_table *= "$(label) & $(round(avg[1][1], digits=1))-$(round(avg[1][2], digits=1)) & $(round(avg[2][1], digits=1))-$(round(avg[2][2], digits=1)) & $(round(avg[4][2], digits=1))-$(round(avg[4][1], digits=1)) \\\\\n"
    end
    latex_table *="\\hline"
end
latex_table *= """
\\end{tabular}}
\\end{table}
"""

print(latex_table)