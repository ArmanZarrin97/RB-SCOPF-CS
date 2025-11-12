using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn))
DIR="../../../src/"
include(DIR*"Funcs/IO/Packages.jl") # load necessary packages
include(DIR*"Funcs/IO/utils.jl") # load general functions
include(DIR*"Assets/MPSCOPF.jl") # load the multi-period security-constrained optimal power flow model
include(DIR*"Funcs/IO/Result_generation.jl") # save MPSCOPF results in Output file
include(DIR*"Assets/ACPF.jl") # load power flow model
include(DIR*"Assets/Problems.jl") # load the different models of problems

##########################################################################################################################################################
#***#***#***##***#***#***##***#***#***##***#***#***##***#***#**# Part I: Problem Setting #***#***#***##***#***#***##***#***#***##***#***#***##***#***#***#
##########################################################################################################################################################

################## Select the Test System and the type of  Problem #################
global Test_Sys="Florin5/"
Initialization=true # Initialization: true=>Run DA_SCOPF before the simulation, false=>No need to run DA_SCOPF
Simulation_Type="all" # Simulation_Type: MP-OPF=>Multi-period OPF, MP-SCOPF=>Multi-period SCOPF, SP-SCOPF=>Single-period SCOPF, all=>All
ERORR=Dict(); ERORR["DA"]=[0.1 0.1];ERORR["RT"]=[0.05 0.05];##default values
#-->||| Day-ahead OPF=>DA_OPF||| Real-time OPF=>RT_OPF||| Day-ahead SCOPF=>DA_SCOPF||| Real-time SCOPF=>RT_SCOPF|||
global printlevel=0; ## To control terminal prints
new_profile=[true true]; ## [true true]=>Generate new profiles, [false false]=>Read profiles from the file, [true same]=>Generate new profiles for both DA and use same for RT
##############################------------------------##############################

###################### Set global Vars and Read Initial Points #####################
Set_global_parameters(Test_Sys)
ERORR["DA"]=[0 0]  #temporary
global Init_Data=InitData_Reader(Path,["Zero"])
######################----------------------------------------######################

################################### Profile Reader #################################

True_Profile = True_Profile_Reader() # Read Real profiles
## new Profile generation ##
if new_profile[1]
    DA_Profiles=DA_Profile_Generator(True_Profile;Data_=[("Load",ERORR["DA"][1]),("Wind",ERORR["DA"][2])])
end
if new_profile[2]
    RT_Profiles=RT_Profile_Generator(True_Profile;Data=[("Load",ERORR["RT"][1]),("Wind",ERORR["RT"][2])])
elseif new_profile[2]=="same"
    RT_Profiles=RT_Profile_Generator(DA_Profiles;Data=[("Load",0),("Wind",0)])
end

global Profile_Data= Profile_Reader(Path)
######################----------------------------------------######################

Ptype_="DA";init_="Zero";nTP=nTL;TD=Td

raw,initraw,profraw= try Data_Reader(;Ptype=Ptype_ , init=init_)
catch e
    @warn "Problem with running Data_Reader function"
    rethrow(e)
end
glob_input=(nTP,TD,raw,initraw,initraw,Mode,Inc_coef,Dec_coef,CCT,LSpermit)

Prof=Load_Data(nTP,raw,profraw) # Prepare Load and Shunt data and generation setpoints

SCOPF,N_Var,C_Var=try
    MPSCOPF(glob_input,Prof;fixed_SOC_=1,conprob=conprob_,OPFon=1)
catch e
    @warn "Error in MPSCOPF function: $e"
    global debug_hint=glob_input,Prof
    rethrow(e)
end

OPFout,Initials=try  Get_Result(SCOPF,raw,initraw;OPFon=1)
catch e
    @warn "Error in Get_Result function: $e"
    global debug_hint= SCOPF,N_Var,C_Var
    rethrow(e)
end

@assert round(OPFout."OFValue"[1])==1.465013e6
@assert round.(Matrix(OPFout."PG"[1]),digits=3)==repeat([8.470 6.370 1.500],nTL)
global_logger(ConsoleLogger(stderr, Logging.Info))
@info "MP-DA-OPF Test passed!"