include("Funcs/IO/Packages.jl") # load necessary packages
include("Funcs/IO/utils.jl") # load general functions
include("Assets/L-MPSCOPF.jl") # load the multi-period security-constrained optimal power flow model
include("Funcs/IO/Result_generation.jl") # save MPSCOPF results in Output file
include("Assets/ACPF.jl") # load power flow model
include("Assets/L-Problems.jl") # load the different models of problems
include("Assets/linearization.jl") # load linearization functions
include("Assets/L-Simulations.jl")

##########################################################################################################################################################
#***#***#***##***#***#***##***#***#***##***#***#***##***#***#**# Part I: Problem Setting #***#***#***##***#***#***##***#***#***##***#***#***##***#***#***#
##########################################################################################################################################################

################## Select the Test System and the type of  Problem #################
global Test_Sys="Florin5/";global nTL=24; global FullReport=0; LinMethod="SENS"# temporary values
Initialization=true # Initialization: true=>Run DA_SCOPF before the simulation, false=>No need to run DA_SCOPF
Simulation_Type=["MPC-OPF"] # Simulation Type: ["MPC-OPF","SP-SCOPF","MPC-SCOPF"]
global ERORR=Dict(); ERORR["DA"]=[0.1 0.1];ERORR["RT"]=[0.05 0.05];##default values
#-->||| Day-ahead OPF=>DA_OPF||| Real-time OPF=>RT_OPF||| Day-ahead SCOPF=>DA_SCOPF||| Real-time SCOPF=>RT_SCOPF|||
global printlevel=[1 1]; ## To control terminal prints
new_profile=[true true]; ## [true true]=>Generate new profiles, [false false]=>Read profiles from the file, [true same]=>Generate new profiles for DA and use same for RT
######################----------------------------------------######################
ERORR["DA"]=[0.00 0.00];ERORR["RT"]=[0.00 0.00] #temporary values
##############################------------------------##############################

###################### Set global Vars and Read Initial Points #####################
Set_global_parameters(Test_Sys)


##########################################################################################################################################################
#***#***#***##***#***#***##***#***#***##***#***#***##***#***#**# Part II: Initialization #***#***#***##***#***#***##***#***#***##***#***#***##***#***#***#
##########################################################################################################################################################

# Run Day-Ahead multi-period SCOPF and the result will be used as Market cleared set-points
if Initialization==true
    init_out,init_model=DA_SCOPF()
    global Init_Data=InitData_Reader(Path,["Zero","DA_SCOPF"])
end

##########################################################################################################################################################
#***#***#***##***#***#***##***#***#***##***#***#***##***#***#****# Part III: Simulations #***#***#***##***#***#***##***#***#***##***#***#***##***#***#***#
##########################################################################################################################################################
global TD=CCT; global RT_result=DataFrame(case1=0.00,case2=0.00,case3=0.00);

################################### 1.Real-time OPF #################################
global OPFres=DataFrame(t=1:nTL*nTS,cost=zeros(),EConcost=zeros(),State=DataFrame(),Before=Vector(undef,nTL*nTS),After=Vector(undef,nTL*nTS),time=zeros(),InfRep=Vector(undef,nTL*nTS))

if "MPC-OPF" in Simulation_Type
    Problem_Type="OPF"
    Lininit_0=Initialize("Exact",init_out,init_model,0,0)
    MPC_OPF_Simulation(LinMethod,Lininit_0)
    @save "OPFres.jld2" OPFres RT_result
end

################################ 2.Real-time SP-SCOPF ##############################
global SP_SCOPFres=DataFrame(t=1:nTL*nTS,cost=zeros(),cost_pre=zeros(),EConcost=zeros(),State=DataFrame(),Before=Vector(undef,nTL*nTS),After=Vector(undef,nTL*nTS),time=zeros(),InfRep=Vector(undef,nTL*nTS))

if "SP-SCOPF" in Simulation_Type
    SP_SCOPF_Simulation()
    @save "SP-SCOPFres.jld2" SP_SCOPFres RT_result
end

################################ 3.Real-time MP-SCOPF ##############################
global MP_SCOPFres=DataFrame(t=1:nTL*nTS,cost=zeros(),cost_pre=zeros(),EConcost=zeros(),State=DataFrame(),Before=Vector(undef,nTL*nTS),After=Vector(undef,nTL*nTS),iteration=zeros(),time=zeros(),InfRep=Vector(undef,nTL*nTS))

if "MPC-SCOPF" in Simulation_Type
    MPC_SCOPF_Simulation()
    @save "MPC_SCOPFres.jld2" MP_SCOPFres RT_result
end