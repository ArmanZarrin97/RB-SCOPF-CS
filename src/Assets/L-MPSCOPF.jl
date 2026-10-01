
include("../Funcs/SCOPF/L_Model_Build.jl")
include("../Funcs/SCOPF/L_PowerFlow.jl")

# export L_MPSCOPF
function L_MPSCOPF(Glob,glob_input,Proffin,casreps,h_cas=DataFrame(),comcas=DataFrame(); fixed_SOC_=1,conprob=0.05,OPFon=0,Winit=0,t0=1,WP=1,LinMethod="None",Lininit=Dict{String,Any}(),relaxconfac=0.0,optimal_tcs=nothing) #Array{Union{Float64,Array{Float64}}} # WP: without Penalty #OPFon: 0=>SCOPF, 1=>OPF #Winit: 0=>No initial value, 1=>Initial value #t0: first time period # fixed_SOC_: 0=>No fixed SOC, 1=>Fixed SOC # conprob: contingency probability summation # Proffin: Profile information
    # casreps=nothing
    #conprob=conprob_;LinMethod=LinMethod;Lininit=OP_init
    #Winit=0;t0=1;WP=1;
    #Winit=0;t0=1;WP=1;fixed_SOC_=1;OPFon=0;
    if LinMethod!="None" @assert !isempty(Lininit) end #"Linearization method is selected but no initial value is provided!"
    ##-----------------SETTING INITIAL PARAMETERS---------------------------
    nTP,TD,raw,initraw,initial,Mode,Inc_coef,Dec_coef,CCT,LSpermit,ctdall=glob_input
    fixed_SOC=fixed_SOC_;
    sl_F,sl_V=0,0
    ##------------------ Model Intialization --------------------------------
    LSCOPF =Initialize_Model(LinMethod)

    if LinMethod == "SENS"
        Rampvio_N,Rampvio_C=find_ramp_rate_violations(Glob,Lininit,raw,OPFon)
        LSCOPF[:Rampvio]=[Rampvio_N,Rampvio_C]
    else
        LSCOPF[:Rampvio]=[0,0]
    end
    # Rampvio_C[1,15,3]
    # findmax(Rampvio_C)
    ##------------------ Normal Operation Model -----------------------------
    LSCOPF,N_Var=Normal_PB(Glob,TD, LSCOPF,nTP,raw,initraw,ctdall,fixed_SOC,Proffin,Winit,initial,LinMethod,Lininit);

    # ------------------- Cascade Model -----------------------------
    if Glob["Cascade_Opt"] && (OPFon==0 && raw["nConL"]!=0 && (!isnothing(casreps) || !isempty(h_cas) || !isnothing(comcas)) && LinMethod == "SENS")
        LSCOPF, CascadeInfo=Cascade_L_PB(LSCOPF,nTP,raw,casreps,h_cas,comcas,Lininit,Glob,optimal_tcs);
        # CascadeInfo["optimal_tcs"]=optimal_tcs
        # CascadeInfo=Dict()
        # LSCOPF[:CascadeInfo]=CascadeInfo

        # LSCOPF[:CNcascade_LP][1]
        # LSCOPF[:slack_Cascade_LP]
    else
        CascadeInfo=Dict()
    end
    # CascadeInfo=Dict()

    C_Var= 0,0,0,0,0,0,0,0,0,[0,0],[0],Dict()
    ##------------------ Contingency Operation Model ------------------------
    if (OPFon==0 && raw["nConL"]!=0)
        LSCOPF,C_Var=Contingency_L_PB(Glob, LSCOPF,nTP,raw,initraw,ctdall,fixed_SOC,N_Var,Proffin,CCT,initial,LSpermit,LinMethod,Lininit,relaxconfac);
    end
    ##------------------ Setting SENS Constraints -------------------------------
    if LinMethod =="SENS"
        LSCOPF= Sensitivity_constraints(Glob,LSCOPF,CascadeInfo,N_Var,C_Var,OPFon,raw,Lininit,ctdall)
        # LSCOPF= Sensitivity_constraints(Glob,LSCOPF,Dict(),N_Var,C_Var,OPFon,raw,Lininit,ctdall)
    else
        LSCOPF[:modTR]=[0.01 0.01]
        LSCOPF[:acc]=[0.0 0.0]
    end
    ##------------------ Setting Objective Function ---------------------
    LSCOPF=Objective_function(Glob, TD, LSCOPF,CascadeInfo,nTP,raw,initraw,Proffin,N_Var,C_Var,sl_F,sl_V,t0,Lininit,LinMethod;lambda=conprob,OPFon=OPFon,WP=WP)
    # #----------------Solve and Results --------------------------------------
    # # @constraint(LSCOPF,C_acc, LSCOPF[:slk_Cont]<=0.9)
    # @info "Model is ready"
    # optimize!(LSCOPF)
    # SCOPF=LSCOPF
    # println("Status:  ",termination_status(SCOPF))
    # println("Objective value:  ", value(SCOPF[:real_objective]))
    # println("Objective value+Cas:  ", value(SCOPF[:real_objective_cas]))
    # println("Solver Time:   ", JuMP.solve_time(SCOPF))
    # println("slk_normal_P:",value.(LSCOPF[:slk_Norm]))
    # println("slk_cont_P:",value.(LSCOPF[:slk_Cont]))
    # println("slk_normal_Q:",value.(LSCOPF[:slk_Norm1]))
    # println("slk_cont_Q:",value.(LSCOPF[:slk_Cont1]))
    # if Cascade_Opt
    #     println("slk_cascade_LP:",value.(LSCOPF[:cascade_penalty]))
    # end
    
    # println("cost_pen:",value.(LSCOPF[:cost_TR_relaxC]))
    # objective_value(SCOPF)
    
    # resu=DataFrame(obj=objective_value(SCOPF),real_objective_cas=value(SCOPF[:real_objective_cas]),real_objective=value(SCOPF[:real_objective]),cost_TR_relaxN=value(SCOPF[:cost_TR_relaxN]),cost_TR_PenN=value(SCOPF[:cost_TR_PenN]),cost_TR_relaxC=value(SCOPF[:cost_TR_relaxC]),cascade_penalty=value(SCOPF[:cascade_penalty]))
    # push!(resu, (objective_value(SCOPF),value(SCOPF[:real_objective_cas]),value(SCOPF[:real_objective]),value(SCOPF[:cost_TR_relaxN]),value(SCOPF[:cost_TR_PenN]),value(SCOPF[:cost_TR_relaxC]),value(SCOPF[:cascade_penalty])))

    # fix(LSCOPF[:bin_Cas][(18,4,1)],1;force=true)
    LSCOPF=Solve_P(LSCOPF)
    # value.(LSCOPF[:CV2magexp])[1,38,1]
    # LSCOPF[:CV2magexp][1,38,1]
    # LSCOPF[:slkDataFrame]["slkVmax"]

    # LSCOPF[:TR]
    # # Plots.plot((value.(LSCOPF[:slk_cas])).data)
    # sum(SCOPF[:slk_cas].*CascadeInfo["LNS_tc"]).*coef6

    # value.(LSCOPF[:cascade_penalty])
    # value.(LSCOPF[:slk_cas])



    @info "Solve compeleted!"
    if LinMethod != "None"
        if LinMethod=="Yang"
            LSCOPF[:e]=sqrt.(LSCOPF[:V2]).*cos.(LSCOPF[:θ])
            LSCOPF[:f]=sqrt.(LSCOPF[:V2]).*sin.(LSCOPF[:θ])
            if (OPFon==0 && raw["nConL"]!=0)
                LSCOPF[:e_c_Line]= sqrt.(LSCOPF[:V2_c_Line]).*cos.(LSCOPF[:θ_c_Line])
                LSCOPF[:f_c_Line]= sqrt.(LSCOPF[:V2_c_Line]).*sin.(LSCOPF[:θ_c_Line])
            end
        elseif LinMethod=="SOCP"
            @error "SOCP is not implemented yet!"
            LSCOPF=Find_θ_SOCP(LSCOPF)
            LSCOPF[:e]=value.(LSCOPF[:V2]).*cos.(value.(LSCOPF[:θ]))
            LSCOPF[:f]=value.(LSCOPF[:V2]).*sin.(value.(LSCOPF[:θ]))
        elseif LinMethod == "SENS"
            θ,V=find_state_variables(nTP, LSCOPF,Lininit,raw)
            LSCOPF[:Nθexp]=θ
            LSCOPF[:NV2magexporg]=LSCOPF[:NV2magexp]
            LSCOPF[:NV2magexp]=V
            LSCOPF[:e]= V[:,:,1].*cos.(θ)
            LSCOPF[:f]= V[:,:,1].*sin.(θ)
            LSCOPF[:I]= value.(LSCOPF[:Line_flow])
            LSCOPF[:I_cal]= calculate_I_sol( LSCOPF[:e],LSCOPF[:f],raw)
            if (OPFon==0 && raw["nConL"]!=0)
                θ_c_Line,V_c_Line=find_state_variables(nTP, LSCOPF,Lininit,raw,"Cont")
                LSCOPF[:Cθexp]=θ_c_Line
                LSCOPF[:CV2magexporg]=LSCOPF[:CV2magexp]
                LSCOPF[:CV2magexp]=V_c_Line
                LSCOPF[:e_c_Line]= V_c_Line.*cos.(θ_c_Line)
                LSCOPF[:f_c_Line]= V_c_Line.*sin.(θ_c_Line)
                LSCOPF[:I_c_Line]= value.(LSCOPF[:Line_flow_c_Line])
                LSCOPF[:I_cal_c_Line]= calculate_I_sol( LSCOPF[:e_c_Line],LSCOPF[:f_c_Line],raw)
            end
        end
    end
    return LSCOPF,N_Var,C_Var
end