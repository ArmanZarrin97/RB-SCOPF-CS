
include("../Funcs/SCOPF/L_Model_Build.jl")
include("../Funcs/SCOPF/L_PowerFlow.jl")

export L_MPSCOPF
function L_MPSCOPF(glob_input,Proffin; fixed_SOC_=1,conprob=0.05,OPFon=0,Winit=0,t0=1,WP=1,LinMethod="None",Lininit=Dict{String,Any}(),relaxconfac=0.0) #Array{Union{Float64,Array{Float64}}} # WP: without Penalty #OPFon: 0=>SCOPF, 1=>OPF #Winit: 0=>No initial value, 1=>Initial value #t0: first time period # fixed_SOC_: 0=>No fixed SOC, 1=>Fixed SOC # conprob: contingency probability summation # Proffin: Profile information
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
        Rampvio_N,Rampvio_C=find_ramp_rate_violations(Lininit,raw,OPFon)
        LSCOPF[:Rampvio]=[Rampvio_N,Rampvio_C]
    else
        LSCOPF[:Rampvio]=[0,0]
    end
    ##------------------ Normal Operation Model -----------------------------
    LSCOPF,N_Var=Normal_PB(LSCOPF,nTP,raw,initraw,ctdall,fixed_SOC,Proffin,Winit,initial,LinMethod,Lininit);

    # findall(N_Var[15][1].<N_Var[15][2])
    # findall(N_Var[15][3].<N_Var[15][4])

    # @assert begin
    #     sum(N_Var[15][1].< -1e-5)==0
    #     sum(N_Var[15][2].> 1e-5)==0
    #     sum(N_Var[15][3].< -1e-5)==0
    #     sum(N_Var[15][4].> 1e-5)==0
    #         end

    C_Var= 0,0,0,0,0,0,0,0,0,[0,0],[0],Dict()
    ##------------------ Contingency Operation Model ------------------------
    if (OPFon==0 && raw["nConL"]!=0)
        LSCOPF,C_Var=Contingency_L_PB(LSCOPF,nTP,raw,initraw,ctdall,fixed_SOC,N_Var,Proffin,CCT,initial,LSpermit,LinMethod,Lininit,relaxconfac);
    end
    # C_Var[12]
    # findall(C_Var[12][1].<C_Var[12][2].-1e-5)

    # findall(C_Var[12][3].<C_Var[12][4])
    ##------------------ Setting SENS Constraints -------------------------------
    if LinMethod =="SENS"
        LSCOPF= Sensitivity_constraints(LSCOPF,N_Var,C_Var,OPFon,raw,Lininit,ctdall)
    else
        LSCOPF[:modTR]=[0.01 0.01]
        LSCOPF[:acc]=[0.0 0.0]
    end

    ##------------------ Setting Objective Function ---------------------
    LSCOPF=Objective_function(LSCOPF,nTP,raw,initraw,Proffin,N_Var,C_Var,sl_F,sl_V,t0,Lininit,LinMethod;lambda=conprob,OPFon=OPFon,WP=WP)
    #----------------Solve and Results --------------------------------------
    # fix.(LSCOPF[:shed], 0.0; force=true)
    # fix.(LSCOPF[:ΔPinjP], 0.0; force=true)
    # fix.(LSCOPF[:ΔPinjN], 0.0; force=true)
    # fix.(LSCOPF[:ΔQinjP], 0.0; force=true)
    # fix.(LSCOPF[:ΔQinjN], 0.0; force=true)
    # fix.(LSCOPF[:ΔPinjP_c_Line], 0.0; force=true)
    # fix.(LSCOPF[:ΔPinjN_c_Line], 0.0; force=true)
    # fix.(LSCOPF[:ΔQinjP_c_Line], 0.0; force=true)
    # fix.(LSCOPF[:ΔQinjN_c_Line], 0.0; force=true)
    # @objective(LSCOPF, Min,  0)
    @info "Model is ready"
    LSCOPF=Solve_P(LSCOPF)

    # vars = all_variables(LSCOPF)
    # v = vars[4710]<
    # ctdall[2572,:]
    # # println("Lower bound: ", lower_bound(v))
    # # println("Upper bound: ", upper_bound(v))
    # Lininit["Pshed"][:,1872]
    # raw["PActBus"][589]
    # raw["Load"][1872,:]
    # con_ref=LSCOPF[:NTrust_RegionP2][6,589]
    # delete(LSCOPF, con_ref)
    # con_ref=LSCOPF[:NTrust_RegionP1][6,589]

    # v=LSCOPF[:shed][6,1872]

    # DataFrame(real_obj=value.(LSCOPF[:real_objective]),
    #          TR_relax=value.(LSCOPF[:cost_TR_relaxN]),
    #          TR_relaxC=value.(LSCOPF[:cost_TR_relaxC]),
    #          TR_PenN=value.(LSCOPF[:cost_TR_PenN]),
    #          TR_PenC=value.(LSCOPF[:cost_TR_PenC])
    #          )
    # termination_status(LSCOPF)
    # map=relax_with_penalty!(LSCOPF)
    # LSCOPF[:CEActive_power_balance][8,3,2]
    # OP_init["Pshed_c_Line"][8,3,2]
    # OP_init["Pinc_c_Line"][8,1,2]
    # for (con, penalty) in map
    #     violation = value(penalty)
    #     if violation > 1e-1
    #         println("Constraint `$(name(con))` is violated by $violation")
    #     end
    # end
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
            θ,V=find_state_variables(LSCOPF,Lininit,raw)
            LSCOPF[:Nθexp]=θ
            LSCOPF[:NV2magexporg]=LSCOPF[:NV2magexp]
            LSCOPF[:NV2magexp]=V
            LSCOPF[:e]= V[:,:,1].*cos.(θ)
            LSCOPF[:f]= V[:,:,1].*sin.(θ)
            LSCOPF[:I]= value.(LSCOPF[:Line_flow])
            LSCOPF[:I_cal]= calculate_I_sol( LSCOPF[:e],LSCOPF[:f],raw)
            if (OPFon==0 && raw["nConL"]!=0)
                θ_c_Line,V_c_Line=find_state_variables(LSCOPF,Lininit,raw,"Cont")
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