include("Funcs/IO/Packages.jl") # load necessary packages
include("Funcs/IO/utils.jl") # load general functions
include("Assets/L-MPSCOPF.jl") # load the multi-period security-constrained optimal power flow model
# include("Assets/MPSCOPF.jl") # load the multi-period security-constrained optimal power flow model

# export JULIA_NUM_THREADS=8    #run this code in terminal for multi-threads

#JuMP: Modeling Optimization Problem #Ipopt: Optimization solver package # OdsIO: Input & output for .ods files # DataFrames: !!
include("Assets/ACPF-V1_4.jl") # load power flow model
include("Assets/linearization.jl") # load power flow model
include("Funcs/IO/Result_generation.jl") # save MPSCOPF results in Output file
include("Assets/L-Problems.jl") # load the different models of problems
include("Assets/ACPF_Check.jl") # load the power flow check functions
include("Funcs/Cascading/Cascading_tool.jl")


# global Problem_Type="SCOPF";global Ptype_="DA";
# #---------------------Reading Input files
# global Test_Sys="39bus/";global nTL=24; global Report="feasiblity"# temporary values
# global ERORR=Dict(); ERORR["DA"]=[0.00 0.00];ERORR["RT"]=[0.00 0.00];new_profile=[true true];global printlevel=[1 0]; #[terminal print, ods print];
# global loadbusinjection="false";global errortype="Normal"

# Set_global_parameters(Test_Sys)

Problem_Type="SCOPF";Ptype_="DA";
#---------------------Reading Input files
Test_Sys="Nordic60/";nTL=24;TS=1; Report="feasiblity"# temporary values
ERORR=Dict(); ERORR["DA"]=[0.00 0.00];ERORR["RT"]=[0.00 0.00];new_profile=[true true];printlevel=[1 0]; #[terminal print, ods print];
Glob=Set_global_parameters1(Test_Sys)

#----------------------Set Initial values
Linearization_INIT=Dict(1=>"Flat",2=> "Exact", 3=> "DA5",4=> "DA10",5=> "DA20", 6=>"DA30" ) #
Linearization_TYPE=Dict(1=>"Yang",2=> "SOCP", 3=> "SENS", 4=> "None", 5=> "Hybrid")

if Test_Sys=="Polish/"
    Q_eff=2; P_eff=4;
    Q_con_percent=0.05; P_con_percent=0.01;
    Q_con_percent_coupled=0.01; P_con_percent_coupled=0.01;
    Pen_S=0.0001 *24/nTL; Pen_W=0.1 *24/nTL # Penalty factors in the objective function
    genprox=50
elseif Test_Sys=="39bus/"
    Q_eff=2; P_eff=4;
    Q_con_percent=0.08; P_con_percent=0.05;
    Q_con_percent_coupled=0.02; P_con_percent_coupled=0.01;
    Pen_S=0.0001 *24/nTL; Pen_W=0.1 *24/nTL # Penalty factors in the objective function
    genprox=10
    Glob["I_relax"]=1.03
    Glob["I_relaxmax"]=1.1
elseif Test_Sys=="Nordic60/"
    Q_eff=2; P_eff=4;
    Q_con_percent=0.7; P_con_percent=0.2;
    Q_con_percent_coupled=0.3; P_con_percent_coupled=0.2;
    Pen_S=0.01 *24/nTL; Pen_W=0.1 *24/nTL
    # OP_init0=Init_OP[Linearization_INIT[Glob["Runcode"][1]]]# Penalty factors in the objective function
    # SysData["Bus"][14,5:6].=[1.1,0.9]
    genprox=10

elseif Test_Sys=="Florin5/"
    Q_eff=2; P_eff=4;
    Q_con_percent=1; P_con_percent=0.5;
    Q_con_percent_coupled=1; P_con_percent_coupled=0.5;
    Pen_S=0.0001 *24/nTL; Pen_W=0.1 *24/nTL # Penalty factors in the objective function
end
relaxpermit=true
Glob["relaxpermit"]=true
Glob["monte_carlo"]=false
include("Assets/Risk.jl")

# DA_out,DA_model,init_out,init_model=Initial_Set_points()
# SLP_Report_Prepare(init_model,init_out)
# Init_OP=Dict()
# for Init in collect(values(Linearization_INIT))
#     Init_OP[Init]=Initialize(Init,init_out,init_model,DA_out,DA_model)
# end
# BSON.@save Path*"../Result/NL_Solution.bson" Init_OP SLP_Rep
# #case2: ramp_rates are relaxed, reactive power range is a 50% relaxed

BSON.@load Glob["Path"]*"../Result/NL_Solution.bson" Init_OP SLP_Rep
insertcols!(SLP_Rep, :Max_N_redis, :Max_N_redisQ => ["-"], after=true)
insertcols!(SLP_Rep, :Max_C_redis, :Max_C_redisQ => ["-"], after=true)
# BSON.@load Path*"../Result/NL_Solution_Ramp3.bson" Init_OP SLP_Rep

# V_1,θ_1,PG_1,QG_1,cost1=Matrix(init_out."mag_v"[1]),Matrix(init_out."ang_v"[1]),Matrix(init_out."PG"[1]),Matrix(init_out."QG"[1]),init_out."OFValue_WP"[1]
# if Problem_Type=="SCOPF"
#     V_c_Line_1,θ_c_Line_1,PG_c_Line_1,QG_c_Line_1=sqrt.(value.(init_model[:e_c_Line].^2 +init_model[:f_c_Line].^2)),atan.(value.(<init_model[:f_c_Line]),value.(init_model[:e_c_Line])),value.(init_model[:PG_c_Line]),value.(init_model[:QG_c_Line])
# end
Snapshots=Dict()
Runcode=(3,3)
# @load "NL_Solution.jld2"

OP_NLP=Init_OP["Exact"]
# ``for Runcode in [ (6,3) (3,3) (4,3) (5,3) ] # (4,3) (5,3) (2,1) (2,3) (1,1) (2,1) (3,1) (4,1) (5,1) (6,1) (7,1) $#(2,3) (3,3) (4,3) (1,1) (1,3) (5,1) (5,3) (2,1) (3,1) (4,1) (5,1)
# for Runcode in [ (3,1) (4,1) (5,1) (6,1) (3,3) (6,3) (4,3) (5,3)  (6,5) (5,5) (4,5) (3,5) ] #(3,1) (4,1) (5,1) (6,1) (4,3) (5,3) (2,1) (2,3) (1,1) (2,1) (3,1) (4,1) (5,1) (6,1) (7,1) $#(2,3) (3,3) (4,3) (1,1) (1,3) (5,1) (5,3) (2,1) (3,1) (4,1) (5,1) 
for Runcode in [(3,3) (6,3) (4,3) (5,3)] #(3,1) (4,1) (5,1) (6,1) (4,3) (5,3) (2,1) (2,3) (1,1) (2,1) (3,1) (4,1) (5,1) (6,1) (7,1) $#(2,3) (3,3) (4,3) (1,1) (1,3) (5,1) (5,3) (2,1) (3,1) (4,1) (5,1) 
    #(6,5) (5,5) (4,5) (3,5)

    Glob["Runcode"]=Runcode
    OP_init0=Init_OP[Linearization_INIT[Glob["Runcode"][1]]]
    OP_init0["acc"]=[0.0 0.0]

    if nTL==1
        for key in keys(OP_init0)
            # key= first(keys(OP_init0))
            if size(OP_init0[key], 1) == 24
                OP_init0[key]=OP_init0[key][1:nTL,ntuple(d -> :, ndims(OP_init0[key])-1)...]
            end
        end

        for key in keys(OP_NLP)
            # key= first(keys(OP_init0))
            if size(OP_NLP[key], 1) == 24
                OP_NLP[key]=OP_NLP[key][1:nTL,ntuple(d -> :, ndims(OP_NLP[key])-1)...]
            end
        end
    end

    # mincost,Best_solution,iter,min_infmax_N,min_infmax_C,minOptGap,Best_vio_general_N,Best_vio_sum_N,Best_vio_general_C,Best_vio_sum_C,sol_time,restoration_time,Best_Cas_MW, Best_Cas_MW_mean,Cascade_snapshots= SLA_clean(Glob,OP_init0,OP_NLP;BItermax=7)
    mincost,Best_solution,iter,min_infmax_N,min_infmax_C,minOptGap,Best_vio_general_N,Best_vio_sum_N,Best_vio_general_C,Best_vio_sum_C,sol_time,restoration_time,Best_Cas_MW, Best_Cas_MW_mean,Cascade_snapshots= SLA_RCAS(Glob,OP_init0,OP_NLP;BItermax=7)

    Snapshots[Runcode]=Cascade_snapshots
    if Test_Sys in ["39bus/", "Nordic60/"]
        makecasfig(Runcode,Cascade_snapshots,Glob["Path"])
    end
    Max_N_redis,Max_C_redis,Max_N_redisQ,Max_C_redisQ=find_total_redispatch!(Glob, Best_solution,OP_init0)

    SLP_Rep=Add_SLP_report!(SLP_Rep,Runcode,Best_solution,mincost,minOptGap,min_infmax_N,min_infmax_C,Best_vio_general_N,Best_vio_general_C,iter,sol_time,restoration_time,Max_N_redis,Max_C_redis,Max_N_redisQ,Max_C_redisQ)
end

# BSON.@save Glob["Path"]*"../Result/Cascade_Snapshots.bson" Snapshots
push!(SLP_Rep,repeat(["end"],size(SLP_Rep,2)))
ods_write(Glob["Path"]*"../Result/LIN_Sum_"*Problem_Type,Dict(("General",1,1)=>SLP_Rep))


# BSON.@load Glob["Path"]*"../Result/Cascade_Snapshots.bson" Snapshots

# for g in 1:length(Cascade_snapshots)
#     display(Cascade_snapshots[g])
# end
# snaps_old=Snapshots[Runcode]
# # @save "SLP_Rep_Yang.jld2" SLP_Rep
# # @load "SLP_Rep_Yang.jld2" SLP_Rep

# Snapshots[(4, 3)][1]
