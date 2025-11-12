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


global Problem_Type="SCOPF";global Ptype_="DA";
#---------------------Reading Input files
global Test_Sys="Polish/";global nTL=36; global Report="feasiblity"# temporary values
global ERORR=Dict(); ERORR["DA"]=[0.00 0.00];ERORR["RT"]=[0.00 0.00];new_profile=[true true];global printlevel=[1 0]; #[terminal print, ods print];
global loadbusinjection="false";global errortype="Normal"

Set_global_parameters(Test_Sys)

#----------------------Set Initial values
Linearization_INIT=Dict(0 => "DC", 1=>"Flat",2=> "Exact" ) #
Linearization_TYPE=Dict(1=>"Yang",2=> "SOCP", 3=> "SENS", 4=> "None", 5=> "Hybrid")

if Test_Sys=="Polish/"
    global Q_eff=2; global P_eff=4;
    Q_con_percent=0.05; P_con_percent=0.01;
    Q_con_percent_coupled=0.01; P_con_percent_coupled=0.01;
    global Pen_S=0.0001; Pen_W=0.1 # Penalty factors in the objective function
    global genprox=50
elseif Test_Sys=="39bus/"
    global Q_eff=2; global P_eff=4; 
    Q_con_percent=0.08; P_con_percent=0.05;
    Q_con_percent_coupled=0.02; P_con_percent_coupled=0.01;
    global Pen_S=0.0001; Pen_W=0.1 # Penalty factors in the objective function
    global genprox=10
elseif Test_Sys=="Nordic60/"
    global Q_eff=2; global P_eff=4; 
    Q_con_percent=0.7; P_con_percent=0.2;
    Q_con_percent_coupled=0.3; P_con_percent_coupled=0.2;
    global Pen_S=0.01; Pen_W=0.1 # Penalty factors in the objective function
    # SysData["Bus"][14,5:6].=[1.1,0.9]
    global genprox=10
elseif Test_Sys=="Florin5/"
    global Q_eff=2; global P_eff=4; 
    Q_con_percent=1; P_con_percent=0.5;
    Q_con_percent_coupled=1; P_con_percent_coupled=0.5;
    global Pen_S=0.0001; Pen_W=0.1 # Penalty factors in the objective function
end

# DA_out,DA_model,init_out,init_model=Initial_Set_points()
# SLP_Report_Prepare(init_model,init_out)
# Init_OP=Dict()
# for Init in collect(values(Linearization_INIT))
#     Init_OP[Init]=Initialize(Init,init_out,init_model,DA_out,DA_model)
# end
# BSON.@save Path*"../Result/NL_Solution.bson" Init_OP SLP_Rep
# #case2: ramp_rates are relaxed, reactive power range is a 50% relaxed
Init00=Dict{String,Any}()
Init00["V"]=ones(nTL,size(SysData["Bus"],1))
Init00["V_c_Line"]=ones(nTL,size(SysData["Line"],1),size(SysData["ConL"],1))
Init00["θ"]=zeros(nTL,size(SysData["Bus"],1))
Init00["θ_c_Line"]=zeros(nTL,size(SysData["Line"],1),size(SysData["ConL"],1))
Init00["Pshed"]=zeros(nTL,size(SysData["Load"],1))
Init00["PG"]=zeros(nTL,size(SysData["Gen"],1))
Init00["QG"]=zeros(nTL,size(SysData["Gen"],1))
Init00["warmstart_yang"]=0;Init00["TR"]=[1.0 0.001];Init00["TRB"]=[0.0001];Init00["FilterIndex"]=[0.7]

# OP_NLP["warmstart_yang"]
# OP_NLP["Pshed"]
# OP_NLP["θ_c_Line"]
OP,out,model=Optimizer(Linearization_TYPE[1],deepcopy(Init00))
                # model[:e_c_Line]= sqrt.(model[:V2_c_Line]).*cos.(model[:θ_c_Line])
                # model[:f_c_Line]= sqrt.(model[:V2_c_Line]).*sin.(model[:θ_c_Line])
SLP_Report_Prepare(model,out)


Runcode=(0,3)
# @load "NL_Solution.jld2"

OP_NLP=OP

# for Runcode in [ (6,3) (3,3) (4,3) (5,3) ] # (4,3) (5,3) (2,1) (2,3) (1,1) (2,1) (3,1) (4,1) (5,1) (6,1) (7,1) $#(2,3) (3,3) (4,3) (1,1) (1,3) (5,1) (5,3) (2,1) (3,1) (4,1) (5,1)
# for Runcode in [  (3,3) (6,3) (4,3) (5,3) (6,5) (5,5) (4,5) (3,5) (3,1) (4,1) (5,1) (6,1) ] # (4,3) (5,3) (2,1) (2,3) (1,1) (2,1) (3,1) (4,1) (5,1) (6,1) (7,1) $#(2,3) (3,3) (4,3) (1,1) (1,3) (5,1) (5,3) (2,1) (3,1) (4,1) (5,1)
#(6,5) (5,5) (4,5) (3,5)
    global relaxpermit=true
    OP_init0=deepcopy(OP_NLP)
    OP_init0["acc"]=[0.0 0.0]
    mincost,Best_solution,iter,min_infmax_N,min_infmax_C,minOptGap,Best_vio_general_N,Best_vio_sum_N,Best_vio_general_C,Best_vio_sum_C,sol_time,restoration_time= SLA_clean(Runcode,OP_init0,OP_NLP;BItermax=7)
    Max_N_redis,Max_C_redis,Max_N_redisQ,Max_C_redisQ=find_total_redispatch!(Best_solution,OP_init0)

    SLP_Rep=Add_SLP_report!(SLP_Rep,Runcode,Best_solution,mincost,minOptGap,min_infmax_N,min_infmax_C,Best_vio_general_N,Best_vio_general_C,iter,sol_time,restoration_time,Max_N_redis,Max_C_redis,Max_N_redisQ,Max_C_redisQ)
# end

push!(SLP_Rep,repeat(["end"],size(SLP_Rep,2)))
ods_write(Path*"../Result/LIN_Sum_"*Problem_Type,Dict(("General",1,1)=>SLP_Rep))


# @save "SLP_Rep_Yang.jld2" SLP_Rep
# @load "SLP_Rep_Yang.jld2" SLP_Rep