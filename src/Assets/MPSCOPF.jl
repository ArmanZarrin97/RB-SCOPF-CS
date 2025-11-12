
    include("../Funcs/SCOPF/Model_Build.jl")
    include("../Funcs/SCOPF/PowerFlow.jl")

    export MPSCOPF
    function MPSCOPF(glob_input,Proffin; fixed_SOC_=1,conprob=0.05,OPFon=0,Winit=0,t0=1,WP=1) # WP: without Penalty #OPFon: 0=>SCOPF, 1=>OPF #Winit: 0=>No initial value, 1=>Initial value #t0: first time period # fixed_SOC_: 0=>No fixed SOC, 1=>Fixed SOC # conprob: contingency probability summation # Proffin: Profile information
        ##-----------------SETTING INITIAL PARAMETERS---------------------------
        nTP,TD,raw,initraw,initial,Mode,Inc_coef,Dec_coef,CCT,LSpermit=glob_input
        fixed_SOC=fixed_SOC_;
        sl_F,sl_V=0,0
        ##------------------ Model Intialization --------------------------------
        SCOPF =Initialize_Model()
        ##------------------ Normal Operation Model -----------------------------
        SCOPF,N_Var=Normal_PB(SCOPF,nTP,raw,initraw,fixed_SOC,Proffin,Winit,initial)
        C_Var= 0,0,0,0,0,0,0,0,0
        ##------------------ Contingency Operation Model ------------------------
        if (OPFon==0 && raw["nConL"]!=0) SCOPF,C_Var=Contingency_L_PB(SCOPF,nTP,raw,initraw,fixed_SOC,N_Var,Proffin,CCT,initial,LSpermit) end
        ##------------------ Setting Objective Function ---------------------
        SCOPF=Objective_function(SCOPF,nTP,raw,initraw,Proffin,N_Var,C_Var,sl_F,sl_V,t0;lambda=conprob,OPFon=OPFon,WP=WP)
        #----------------Solve and Results --------------------------------------
        @info "Model is ready"
        SCOPF=Solve_P(SCOPF)
        @info "Solve compeleted!"
        return SCOPF,N_Var,C_Var
    end