function Initialize_Model(LinMethod)
    if LinMethod=="None"
        model = Model(Ipopt.Optimizer)
        set_optimizer_attributes(model,"tol"=>1e-04)
        set_attribute(model, "hsllib", HSL_jll.libhsl_path)
        set_attribute(model, "linear_solver", "ma97")
        if(printlevel==0) set_optimizer_attribute(model, "print_level",0) end
    else
        model = Model(Ipopt.Optimizer)
        # set_optimizer_attributes(model,"tol"=>1e-04)
        # set_attribute(model, "hsllib", HSL_jll.libhsl_path)
        # set_attribute(model, "linear_solver", "ma97")

        model = Model(CPLEX.Optimizer)
        model =direct_model(CPLEX.Optimizer())
        set_attribute(model, "CPXPARAM_Simplex_Tolerances_Feasibility", 1.5e-5)
        set_attribute(model, "CPX_PARAM_LPMETHOD",4)
        set_optimizer_attribute(model, "CPX_PARAM_SOLUTIONTYPE", 2)
        set_optimizer_attribute(model, "CPXPARAM_Barrier_StartAlg", 4)
        set_optimizer_attribute(model, "CPXPARAM_Barrier_Display", 2)
        # set_optimizer_attribute(model, "CPXPARAM_Threads", 1)
        # set_optimizer_attribute(model, "CPX_PARAM_BARALG", 2)

        set_optimizer_attribute(model, CPLEX.PassNames(), true)
        set_attribute(model, "CPXPARAM_Network_Tolerances_Feasibility", 1.5e-5)

        set_attribute(model,"CPXPARAM_MIP_Tolerances_Integrality",1.5e-5)
        set_attribute(model,"CPXPARAM_Barrier_QCPConvergeTol",1.5e-5)
        set_attribute(model,"CPXPARAM_Barrier_ConvergeTol",1.5e-5)
        set_attribute(model,"CPXPARAM_Barrier_QCPConvergeTol",1.5e-5)
        # set_attribute(model,"CPXPARAM_Barrier_QCPConvergeTol",1e-5)
        set_attribute(model, "CPXPARAM_Simplex_Tolerances_Feasibility", 1.5e-5)
        set_attribute(model, "CPXPARAM_MIP_Tolerances_MIPGap", 1.5e-5)
        set_attribute(model, "CPXPARAM_Simplex_Tolerances_Optimality", 1.5e-5)
        set_attribute(model, "CPXPARAM_MIP_Tolerances_Linearization", 1.5e-5)
        set_attribute(model, "CPXPARAM_Network_Tolerances_Optimality", 1.5e-5)

        


        # set_attribute(model,"CPXPARAM_Barrier_Limits_ObjRange",1e30)
        # set_optimizer_attributes(model,"acceptable_constr_viol_tol"=>1e-04)
        # set_attribute(model, "CPX_PARAM_EPINT", 1e-1)
    end
    # set_attribute(model, "linear_solver", "ma97")
    # set_optimizer_attributes(model,"tol"=>1e-05)
    # # set_optimizer_attributes(model,"constr_viol_tol"=>1e-05)
    # # # set_optimizer_attributes(model,"acceptable_constr_viol_tol"=>1e-05)
    # # # set_optimizer_attributes(model,"acceptable_dual_inf_tol"=>1e-05)
    # set_optimizer_attributes(model,"constr_viol_tol"=>1e-05)
    # set_optimizer_attributes(model,"dual_inf_tol"=>1e-05)
    # set_optimizer_attributes(model,"compl_inf_tol"=>1e-05)
    # set_optimizer_attributes(model,"bound_relax_factor"=>0.0)
    # set_optimizer_attributes(model,"mu_target"=>1e-02)
    # set_optimizer_attributes(model,"mu_strategy"=>"adaptive")
    # set_optimizer_attributes(model,"mumps_mem_percent"=>100)
    # set_optimizer_attributes(model,"max_iter"=>100)
    return model
end
function Solve_P(SCOPF)
    optimize!(SCOPF)
    # println("Objective penalized: ", value(SCOPF[:real_objective]))
    println("Status:  ",termination_status(SCOPF))
    println("Objective value:  ", value(SCOPF[:real_objective]))
    println("Solver Time:   ", JuMP.solve_time(SCOPF))
    return SCOPF
end

function Objective_function(SCOPF,nTP,raw,initraw,Proffin,N_Var,C_Var,sl_F,sl_V,t0,Lininit,LinMethod="None",Lccost=1e4,Rccost=0;lambda=0.05,OPFon=0,WP=1)
    # SCOPF,nTP,raw,initraw,Proffin,N_Var,C_Var,sl_F,sl_V,t0=LSCOPF,nTP,raw,initraw,Proffin,N_Var,C_Var,sl_F,sl_V,t0
singleconprob= (raw["nConL"]==0 ? 0 : 1/raw["nConL"])
#####################################
#####################################
Lccost=raw["LScost"]*raw["sbase"]
Rccost=raw["RCcost"]*raw["sbase"]

FactorT=max(1,500*(t0-132))
λG=0.005*FactorT*raw["Gen"]."cost_b"./(2)
λS=0.005*FactorT*((144)^2)*raw["ESS"]."cost_b"./(2*TD)
λF=0.005*FactorT*((144)^2)*raw["Load"]."ShCost"[raw["FL"]."iFL"]./(2*TD)
    #---------------- convert pu to MW for cost function
#real cost of Lc and Rc = value/sbase

    #---------------- convert pu to MW for cost function
PG,gar,P_ch,P_dis,P_inc,P_dec,gar,shed,curt,curtc,soc,UE,Vars =N_Var
PG_c_Line,gar,P_ch_c_Line,P_dis_c_Line,P_inc_c_Line,P_dec_c_Line,shed_c_Line,curt_c_Line,curtc_c_Line =C_Var
sl_F_cost,sl_V_cost=1e5,1e5
Prob_c= zeros(raw["nConL"]).+singleconprob
#Preventive
@expression(SCOPF,cost_ESS[t=1:nTP,i=1:raw["nESS"]], raw["ESS"]."cost_c"[i]+ raw["ESS"]."cost_b"[i]*raw["sbase"]*(P_dis[t,i]+P_ch[t,i])) #+raw["ESS"]."cost_a"[i]*raw["sbase"]^2*(P_dis[t,i]+P_ch[t,i])*(P_dis[t,i]+P_ch[t,i])
@expression(SCOPF,cost_Gen[t=1:nTP,i=1:raw["nGen"]], raw["Gen"]."cost_c"[i]+ raw["Gen"]."cost_b"[i]*raw["sbase"]*(PG[t,i])) #+raw["Gen"]."cost_a"[i]*(raw["sbase"]*(PG[t,i]))^2
@expression(SCOPF,cost_FL[t=1:nTP,i=1:raw["nFL"]], raw["Load"][raw["FL"]."iFL"[i],"ShCost"]*(P_inc[t,i]+P_dec[t,i])*raw["sbase"])
            #Sumup

@expression(SCOPF,cost_Flex[t=1:nTP], (sum(cost_ESS[t,i] for i=1:raw["nESS"]) +sum(cost_Gen[t,i] for i=1:raw["nGen"]) +sum(cost_FL[t,i] for i=1:raw["nFL"])))
@expression(SCOPF,cost_Shed[t=1:nTP], Lccost*sum(shed[t,n] for n=1:raw["nLoad"] ))
@expression(SCOPF,cost_Curt[t=1:nTP], Rccost*sum(curt[t,n] for n=1:raw["nRES"] ))
@expression(SCOPF,cost_Curtc[t=1:nTP], Lccost*sum(curtc[t,n] for n=1:raw["nGen"] ))
if WP==1
    @expression(SCOPF,cost_Pen, 0)
elseif WP=="SP"
    @expression(SCOPF,FPG[i=1:raw["nGen"]], λG[i]*((PG[nTP,i])^2))
    # @expression(SCOPF,FPS[i=1:raw["nESS"]], λS[i]*((Init_Data["DA_SCOPF"]["soc"][div(t0-1,6)+1,i+1]-soc[nTP,i])/(1/TD -mod(t0-1,6)))^2)
    # @expression(SCOPF,FPF[i=1:raw["nFL"]], λF[i]*((Init_Data["DA_SCOPF"]["UE"][div(t0-1,6)+1,i+1]-UE[nTP,i])/(1/TD -mod(t0-1,6)))^2)
    @expression(SCOPF,FPS[i=1:raw["nESS"]], λS[i]*(soc[nTP,i]-initraw["soc"][nTP,i])^2)
    @expression(SCOPF,FPF[i=1:raw["nFL"]], λF[i]*(UE[nTP,i]-initraw["UE"][nTP,i])^2)
    @expression(SCOPF,cost_Pen, sum(FPG)+sum(FPS)+sum(FPF))
    @constraint(SCOPF,socmaxd[i=1:raw["nESS"],t=1:nTP], Init_Data["DA_SCOPF"]["soc"][nTL,i+1]-soc[nTP,i]<= (nTL*nTS-t0)* raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"]*TD)
    @constraint(SCOPF,UEmaxd[i=1:raw["nFL"],t=1:nTP],  Init_Data["DA_SCOPF"]["UE"][nTL,i+1]-UE[nTP,i]<= (nTL*nTS-t0)*Inc_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]]*TD)
else
    @expression(SCOPF,FPG[i=1:raw["nGen"]], λG[i]*((PG[nTP,i])^2))
    @expression(SCOPF,FPS[i=1:raw["nESS"]], λS[i]*(soc[nTP,i]-initraw["soc"][nTP,i])^2)
    @expression(SCOPF,FPF[i=1:raw["nFL"]], λF[i]*(UE[nTP,i]-initraw["UE"][nTP,i])^2)
    tempvar=sum(FPG;init=0)+sum(FPS;init=0)+sum(FPF;init=0)
    @expression(SCOPF,cost_Pen, tempvar)
end
# @expression(SCOPF,cost_Pre, TD* sum(cost_Flex[t]+cost_Curt[t]+cost_Shed[t]  for t=1:nTP))
@expression(SCOPF,cost_Pre_T[t=1:nTP], TD*(cost_Flex[t]+cost_Curt[t]+cost_Curtc[t]+cost_Shed[t]))
# @expression(SCOPF,cost_Pre_T[t=1:nTP], TD*(cost_Flex[t]))
# @expression(SCOPF,cost_Pre_T[t=1:nTP], 0)
#Corrective
if OPFon==0
    @expression(SCOPF,cost_Gen_C[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]], raw["Gen"]."cost_c"[i]+ (raw["Gen"]."cost_b"[i]*raw["sbase"]*(PG_c_Line[t,i,cL]))-cost_Gen[t,i]) #+(raw["Gen"]."cost_a"[i]*(raw["sbase"]*(PG_c_Line[t,i,cL]))^2)
    @expression(SCOPF,cost_ESS_C[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]], raw["ESS"]."cost_c"[i]+ raw["ESS"]."cost_b"[i]*raw["sbase"]*(P_dis_c_Line[t,i,cL]+P_ch_c_Line[t,i,cL])-cost_ESS[t,i]) #+raw["ESS"]."cost_a"[i]*raw["sbase"]^2*(P_dis_c_Line[t,i,cL]+P_ch_c_Line[t,i,cL])^2
    @expression(SCOPF,cost_FL_C[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]], raw["Load"][raw["FL"]."iFL"[i],"ShCost"]*(P_inc_c_Line[t,i,cL]+P_dec_c_Line[t,i,cL])*raw["sbase"]-cost_FL[t,i])
    @expression(SCOPF,cost_Flex_C[t=1:nTP,cL=1:raw["nConL"]], (sum(cost_ESS_C[t,i,cL] for i=1:raw["nESS"]) +sum(cost_Gen_C[t,i,cL] for i=1:raw["nGen"]) +sum(cost_FL_C[t,i,cL] for i=1:raw["nFL"])))

    @expression(SCOPF,cost_Shed_C[t=1:nTP,cL=1:raw["nConL"]], Lccost*sum((shed_c_Line[t,n,cL]) for n=1:raw["nLoad"])-cost_Shed[t])
    @expression(SCOPF,cost_Curt_C[t=1:nTP,cL=1:raw["nConL"]], Rccost*sum((curt_c_Line[t,n,cL]) for n=1:raw["nRES"] )-cost_Curt[t])
    @expression(SCOPF,cost_Curtc_C[t=1:nTP,cL=1:raw["nConL"]], Lccost*sum((curtc_c_Line[t,n,cL]) for n=1:raw["nGen"] )-cost_Curtc[t])
    # @constraint(SCOPF,positive_Post_cost[t=1:nTP,cL=1:raw["nConL"]], cost_Shed_C[t,cL]+cost_Curt_C[t,cL]+cost_Curtc_C[t,cL]+cost_Flex_C[t,cL]>=0 )


    @expression(SCOPF,cost_CS_C[cL=1:raw["nConL"]], sum(cost_Shed_C[t,cL]+cost_Curt_C[t,cL]+cost_Curtc_C[t,cL] for t=1:nTP))

    @expression(SCOPF,cost_Cor_C[cL=1:raw["nConL"]], sum(cost_Flex_C[t,cL]+cost_Shed_C[t,cL]+cost_Curt_C[t,cL]+cost_Curtc_C[t,cL] for t=1:nTP))
    # @constraint(SCOPF,positive_Post_cost[cL=1:raw["nConL"]], cost_Cor_C[cL]>=0 )


    # @expression(SCOPF,cost_Cor_C[cL=1:raw["nConL"]], cost_CS_C[cL])
    @expression(SCOPF,cost_Cor_T, TD*sum(Prob_c*(cost_Cor_C[cL]) for cL=1:raw["nConL"]))
    # @expression(SCOPF,cost_Cor_T_OF, TD*sum(Prob_c*(cost_Cor_C[cL]) for cL=1:raw["nConL"]))
    @expression(SCOPF,cost_Post_T[t=1:nTP], TD*(sum(Prob_c[cL]*(cost_Shed_C[t,cL]+cost_Curt_C[t,cL]+cost_Curtc_C[t,cL]+cost_Flex_C[t,cL]) for cL=1:raw["nConL"])))
    # @expression(SCOPF,cost_SL_F[t=1:nTP,i=1:raw["nLine"],cL=1:raw["nConL"]], sl_F[t,i,cL]*sl_F_cost)
    # @expression(SCOPF,cost_SL_V[t=1:nTP,i=1:raw["nBus"],cL=1:raw["nConL"]], sl_V[t,i,cL]*sl_V_cost)
    # #         #Sumup
    # @expression(SCOPF,cost_SL[cL=1:raw["nConL"]], sum(cost_SL_V[t,i,cL] for t=1:nTP,i=1:raw["nBus"])+sum(cost_SL_F[t,i,cL] for t=1:nTP,i=1:raw["nLine"]))
else
    # @expression(SCOPF,cost_Cor_T_OF, 0)
    @expression(SCOPF,cost_Cor_T, 0)
    @expression(SCOPF,cost_Post_T[t=1:nTP], 0)
end

@expression(SCOPF,cost_T[t=1:nTP], cost_Pre_T[t] +(lambda*cost_Post_T[t]))
@expression(SCOPF,real_objective, sum(cost_T[t] for t=1:nTP)+cost_Pen)
if LinMethod !="SENS"
    @objective(SCOPF, Min,  real_objective)
elseif LinMethod=="SENS"
    # coef1=2e5;coef2=4e5;coef3=1e5;coef4=1000; coef5=mean(SysData["Gen"]."cost_b")*0.5*raw["sbase"]
    # coef4=Lininit["obj"]*0.05;coef1=coef4/1000;coef2=coef4/1000;coef3=1.1*Lccost; coef5=mean(SysData["Gen"]."cost_b")*0.1*raw["sbase"]/100
    # coef4=Lininit["obj"]*0.2;coef1=0.3*Lccost;coef2= 0.3*Lccost;coef3=1.1*Lccost;
    # cof=0.01*Lininit["obj"]/sqrt(length(SCOPF[:slk_Imax])+length(SCOPF[:slk_Vmax])+length(SCOPF[:slk_Vmin])+length(SCOPF[:slk_Imax_c])+length(SCOPF[:slk_Vmax_c])+length(SCOPF[:slk_Vmin_c]))
    coef4=Lininit["obj"]*Pen_W ;coef1=Lininit["obj"]*Pen_S*200;coef2= Lininit["obj"]*Pen_S;coef3=Lininit["obj"]*Pen_S
    # coef5=mean(SysData["Gen"]."cost_b")*0.1*raw["sbase"]
    # coef4=Lininit["obj"]*0.1 ;coef1=cof;coef2= cof;coef3=cof

    # coef1=Lininit["obj"]*0.05;coef2=Lininit["obj"]*0.05;#coef5=Lininit["obj"]*0.05
    # coef4=0
    # coef1=0;coef2=0;coef5=0
    # @expression(SCOPF,cost_TR_PenN, coef5* (sumJuMPi(SCOPF[:ΔQinjP].+SCOPF[:ΔQinjN])+ sumJuMPi(SCOPF[:ΔPinjP].+SCOPF[:ΔPinjN])))
    @expression(SCOPF,cost_TR_PenN,0)
    if relaxpermit==true
        # GBusindex=LinearIndices(zeros(nTP,raw["nBus"]))[:,raw["Gen"].Bus]
        # GBusindex_Vmax=filter(x-> x in GBusindex, axes(SCOPF[:slk_Vmax],1))
        # GBusindex_Vmin=filter(x-> x in GBusindex, axes(SCOPF[:slk_Vmin],1))
        # @constraint(SCOPF,worst_rampN1,SCOPF[:W_ramp].>=SCOPF[:slk_NNramp])
        # @constraint(SCOPF,worst_rampN2,SCOPF[:slk_Norm].>=SCOPF[:slk_NNramp])
        @variable(SCOPF,slk_Norm>=0)
        @variable(SCOPF,slk_Norm1>=0)
        @constraint(SCOPF,Nworst_slk1,slk_Norm.>=SCOPF[:slk_Imax])
        @constraint(SCOPF,Nworst_slk2,slk_Norm1.>=SCOPF[:slk_Vmax])
        @constraint(SCOPF,Nworst_slk3,slk_Norm1.>=SCOPF[:slk_Vmin])
        @constraint(SCOPF,Nworst_slk4,slk_Norm.>=SCOPF[:slk_NNramp])
        @expression(SCOPF,I_relaxN, (coef2) * (sumJuMPi(SCOPF[:slk_Imax])))#/max(1e-6,SCOPF[:W_VIO][1])
        @expression(SCOPF,V_relaxN, (coef1) * ( sumJuMPi(SCOPF[:slk_Vmax]) + sumJuMPi(SCOPF[:slk_Vmin])))#+100*sumJuMPi(SCOPF[:slk_Vmax][GBusindex_Vmax])+100*sumJuMPi(SCOPF[:slk_Vmin][GBusindex_Vmin])))/max(1e-2,SCOPF[:W_VIO][2])
        @expression(SCOPF,ramp_relaxN,(coef3)* sumJuMPi(SCOPF[:slk_NNramp]))
        # @expression(SCOPF,Worst_relaxN, (coef4/max(1e-6,maximum(SCOPF[:W_VIO][1:2])))*SCOPF[:slk_Norm])#
        @expression(SCOPF,Worst_relaxN, coef4 * (slk_Norm+0.5*slk_Norm1))
        @expression(SCOPF,cost_TR_relaxN, I_relaxN+V_relaxN+ramp_relaxN+Worst_relaxN)
        # @expression(SCOPF,cost_TR_relaxN, coef1 * ( sumJuMPi(SCOPF[:slk_Vmax]) + sumJuMPi(SCOPF[:slk_Vmin])) + coef2 * (sumJuMPi(SCOPF[:slk_Imax]))+ coef3* sum(SCOPF[:slk_NNramp1].+SCOPF[:slk_NNramp2]))#+ coef3* (sum(SCOPF[:N_TRfreeP])+ sum(SCOPF[:N_TRfreeQ])+sum(SCOPF[:N_SBfreeP])+sum(SCOPF[:N_SBfreeQ]))
    else
        @expression(SCOPF,cost_TR_relaxN, 0)
    end
    if OPFon==0
        # @expression(SCOPF,cost_TR_PenC, coef5* (sumJuMPi(SCOPF[:ΔQinjP_c_Line].+SCOPF[:ΔQinjN_c_Line])+sumJuMPi(SCOPF[:ΔPinjP_c_Line].+SCOPF[:ΔPinjN_c_Line])))
        @expression(SCOPF,cost_TR_PenC, 0)
        if relaxpermit==true
            # GBusindex_c=LinearIndices(zeros(nTP,raw["nBus"],raw["nConL"]))[:,raw["Gen"].Bus,:]
            # GBusindex_Vmax_c=filter(x-> x in GBusindex_c, axes(SCOPF[:slk_Vmax_c],1))
            # GBusindex_Vmin_c=filter(x-> x in GBusindex_c, axes(SCOPF[:slk_Vmin_c],1))
            # @constraint(SCOPF,worst_rampC1,SCOPF[:slk_Cont].>=SCOPF[:slk_CNramp])
            # @constraint(SCOPF,worst_rampC1,SCOPF[:W_ramp].>=SCOPF[:slk_CNramp])
            @variable(SCOPF,slk_Cont>=0)
            @variable(SCOPF,slk_Cont1>=0)

            @constraint(SCOPF,Cworst_slk1,slk_Cont.>=SCOPF[:slk_Imax_c])
            @constraint(SCOPF,Cworst_slk2,slk_Cont1.>=SCOPF[:slk_Vmax_c])
            @constraint(SCOPF,Cworst_slk3,slk_Cont1.>=SCOPF[:slk_Vmin_c])
            @constraint(SCOPF,Cworst_slk4,slk_Cont.>=SCOPF[:slk_CNramp])
            # @constraint(SCOPF,worst_rampC2,SCOPF[:slk_Cont].>=SCOPF[:slk_CNramp])
            @expression(SCOPF,I_relaxC,(coef2) * (sumJuMPi(SCOPF[:slk_Imax_c])))#/min(1e-2,SCOPF[:W_VIO][1])
            @expression(SCOPF,V_relaxC,(coef1) * (sumJuMPi(SCOPF[:slk_Vmax_c])+ sumJuMPi(SCOPF[:slk_Vmin_c])))#+100* sumJuMPi(SCOPF[:slk_Vmax_c][GBusindex_Vmax_c])+100*sumJuMPi(SCOPF[:slk_Vmin_c][GBusindex_Vmin_c])))/min(1e-2,SCOPF[:W_VIO][2])
            @expression(SCOPF,ramp_relaxC,(coef3)* sumJuMPi(SCOPF[:slk_CNramp]))
            # @expression(SCOPF,ramp_relaxC,0)
            # @expression(SCOPF,Worst_relaxC,(coef4/max(1e-6,maximum(SCOPF[:W_VIO][1:2])))*SCOPF[:slk_Cont])#
            @expression(SCOPF,Worst_relaxC,coef4*(slk_Cont+0.5*slk_Cont1))#
            @expression(SCOPF,cost_TR_relaxC, lambda*(I_relaxC + V_relaxC+ramp_relaxC)+0.9*Worst_relaxC)
            # @expression(SCOPF,cost_TR_relaxC, coef1 * (sumJuMPi(SCOPF[:slk_Vmax_c])+ sumJuMPi(SCOPF[:slk_Vmin_c])) + coef2 * sumJuMPi(SCOPF[:slk_Imax_c]) + coef3 * sumJuMPi(SCOPF[:slk_CNramp1].+SCOPF[:slk_CNramp2]))#+ coef3 *(sumJuMPi(SCOPF[:C_TRfreeP])+ sumJuMPi(SCOPF[:C_TRfreeQ])+sumJuMPi(SCOPF[:C_SBfreeP])+sumJuMPi(SCOPF[:C_SBfreeQ]))+ coef4 * sum(SCOPF[:slkfree_Vmax_c].+SCOPF[:slkfree_Vmin_c]))# 
        else
            @expression(SCOPF,cost_TR_relaxC, 0)
        end
    else
        @expression(SCOPF,cost_TR_PenC, 0)
        @expression(SCOPF,cost_TR_relaxC, 0)
    end
    @objective(SCOPF, Min,  real_objective +cost_TR_relaxN+cost_TR_PenN+cost_TR_relaxC+ lambda*singleconprob*cost_TR_PenC)
    # @objective(SCOPF, Min,  0)
end
# @objective(SCOPF, Min,  0)
return SCOPF
end


function Normal_PB(SCOPF::Model,nTP::Int64,raw::Dict{String, Any},initraw::Dict{String, DataFrame},ctdall,fixed_SOC::Int64,Proffin::Dict{String, Array{Float64}},Winit::Int64,initial::Dict{String, DataFrame},LinMethod::String,Lininit::Dict{String,Any})
    # SCOPF,nTP,raw,initraw,fixed_SOC,Proffin,Winit,initial,LinMethod,Lininit=LSCOPF,nTP,raw,initraw,fixed_SOC,Proffin,Winit,initial,LinMethod,Lininit
    #----------------------OPF State Variables
    SCOPF,Vars=Define_Variable_Normal(SCOPF,nTP,raw,LinMethod,Lininit)
    #----------------------Other Variables
    @variables(
        SCOPF,
        begin
            PG[t = 1:nTP, i = 1:raw["nGen"]] 
            QG[t = 1:nTP, i = 1:raw["nGen"]] 
            P_ch[t = 1:nTP,i = 1:raw["nESS"]] 
            P_dis[t = 1:nTP,i = 1:raw["nESS"]] 
            soc[t = 1:nTP,i = 1:raw["nESS"]]
            UE[t = 1:nTP, i = 1:raw["nFL"]]
            P_inc[t = 1:nTP, i = 1:raw["nFL"]]
            P_dec[t = 1:nTP, i = 1:raw["nFL"]]
            shed[t = 1:nTP, i = 1:raw["nLoad"]]
            curt[t = 1:nTP, i = 1:raw["nRES"]]
            curtc[t = 1:nTP, i = 1:raw["nGen"]]
            # shunt[t = 1:nTP, i = 1:raw["nShunt"]]
        end
    )
    # for  t=1:nTP, i=1:raw["nShunt"]
    #     set_start_value(shunt[t,i], raw["Shunt"]."bsh0"[i])                         # Setting start values
    #     set_lower_bound(shunt[t,i], raw["Shunt"]."bshmin"[i])
    #     set_upper_bound(shunt[t,i], raw["Shunt"]."bshmax"[i])
    # end
    shunt=raw["Shunt"]."bsh0"
    if raw["nRES"]>0
        if isempty(Lininit)
            set_start_value.(curt, 0) else set_start_value.(curt, Lininit["Pcurt"])  end                     # Setting start values
        set_lower_bound.(curt, -1 .* initraw["Pcurt"])
        set_upper_bound.(curt, Proffin["Wind"] .-initraw["Pcurt"])
    end

    if isempty(Lininit)
        set_start_value.(shed, 0) else set_start_value.(shed, Lininit["Pshed"])  end                     # Setting start values
    set_lower_bound.(shed, -1 .* initraw["Pshed"])
    set_upper_bound.(shed, Proffin["PL"] .- initraw["Pshed"])
    # fix.(shed,0; force=true)

    set_lower_bound.(PG, (raw["Gen"]."P_min (pu)")' .- initraw["PG"])
    set_lower_bound.(PG, (raw["Gen"]."P_min (pu)")' .- initraw["PG"])
    set_upper_bound.(PG, (raw["Gen"]."P_max (pu)")' .- initraw["PG"])
    if isempty(Lininit)
        set_start_value.(PG, 0.5*(raw["Gen"]."P_max (pu)")') else set_start_value.(PG, Lininit["PG"])  end
    set_lower_bound.(QG, (raw["Gen"]."Q_min (pu)")' .- initraw["QG"])
    set_upper_bound.(QG, (raw["Gen"]."Q_max (pu)")' .- initraw["QG"])
    if isempty(Lininit)
        set_start_value.(QG, 0.5*(raw["Gen"]."Q_max (pu)")') else set_start_value.(QG, Lininit["QG"])  end
    # set_lower_bound.(curtc, -1* initraw["Pcurtc"])
    fix.(curtc,0)

    # @constraint(SCOPF,NEcurtclim[t=1:nTP,i=1:raw["nGen"]], curtc[t,i]<=PG[t,i]+initraw["PG"][t,i]-initraw["Pcurtc"][t,i])
    @constraint(SCOPF,NEcurtclim, curtc .<=PG .+ initraw["PG"] .- initraw["Pcurtc"])
    if raw["nESS"]>0
        if isempty(Lininit)
            set_start_value.(P_ch, 0.5*(raw["ESS"]."Charge_rating (MW)")'/raw["sbase"]) else set_start_value.(P_ch, Lininit["Pch"])  end
        set_lower_bound.(P_ch, -1 .* initraw["Pch"]) ;set_upper_bound.(P_ch, (raw["ESS"]."Charge_rating (MW)")' ./ raw["sbase"] .- initraw["Pch"]);
        if isempty(Lininit)
            set_start_value.(P_dis, 0.5*(raw["ESS"]."Discharge_rating (MW)")' /raw["sbase"]) else set_start_value.(P_dis, Lininit["Pdis"])  end;
        set_lower_bound.(P_dis, -1 .* initraw["Pdis"]) ;set_upper_bound.(P_dis, (raw["ESS"]."Discharge_rating (MW)")' ./raw["sbase"] .- initraw["Pdis"]);
        if isempty(Lininit)
            set_start_value.(soc, 0.5*(raw["ESS"]."E_rating (MWh)")' ./raw["sbase"]) else set_start_value.(soc, Lininit["soc"]) end;
        set_lower_bound.(soc, (raw["ESS"]."E_rating_min (MWh)")' ./raw["sbase"]);set_upper_bound.(soc, (raw["ESS"]."E_rating (MWh)")'/raw["sbase"]);
    end

    # set_start_value(P_inc[t,i], 0.5*load_inc_prct*prof_ploads[i,t])                         # Setting start values
    # set_start_value(P_dec[t,i], 0.5*load_dec_prct*prof_ploads[i,t])

    if raw["nFL"]==0 || isempty(Lininit)
            set_start_value.(UE, 0);
            set_start_value.(P_inc, 0);
            set_start_value.(P_dec, 0)
    else
        set_start_value.(UE, Lininit["UE"])
        set_start_value.(P_inc, Lininit["Pinc"])
        set_start_value.(P_dec, Lininit["Pdec"])
    end

    if raw["nFL"]>0
        set_lower_bound.(P_inc, -1 .* initraw["Pinc"])
        set_lower_bound.(P_dec, -1 .* initraw["Pdec"])

        # for t=1:nTP, i= 1:raw["nFL"]
        #     set_upper_bound.(P_inc[t,i], Inc_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]]-initraw["Pinc"][t,i])
        #     set_upper_bound.(P_dec[t,i], Dec_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]]-initraw["Pdec"][t,i])
        # end

        set_upper_bound.(P_inc, Inc_coef*Proffin["PL"][:,raw["FL"]."iFL"] .- initraw["Pinc"])
        set_upper_bound.(P_dec, Dec_coef*Proffin["PL"][:,raw["FL"]."iFL"] .- initraw["Pdec"])
    end

    #-------------------------------Power balance constraint #Oc
    #calculating injected power from the bus to grid based on input and output power considering line outage and networc topology
    # if LinMethod !="SENS"
    #     SENSDict=Dict{String, Matrix{Matrix{Float64}}}()
    # elseif LinMethod=="SENS"
    #     @time SENSDict =Cal_SENS_new(Lininit,ctdall,nTP,raw;C="N")
    # end
    #-----------------Line flow Modeling  #OK
    # SENSDict["S_VQ"][1,1]
    # SENSDict["S_VP"][1,1]

    Npij_L, Nqij_L, Npji_L, Nqji_L, Parlist=Flow_Expr_Simple(LinMethod,Lininit,nTP,raw,Vars,'L')
    #
    # @expression(SCOPF, Npij_Lout[t=1:nTP,i=1:raw["nLine"]], Npij_L[t,i])
    # @expression(SCOPF, Nqij_Lout[t=1:nTP,i=1:raw["nLine"]], Nqij_L[t,i])
    # @expression(SCOPF, Npji_Lout[t=1:nTP,i=1:raw["nLine"]], Npji_L[t,i])
    # @expression(SCOPF, Nqji_Lout[t=1:nTP,i=1:raw["nLine"]], Nqji_L[t,i])
    #
    # Npij_L1, Nqij_L1, Npji_L1, Nqji_L1, Parlist1=Flow_Expr(LinMethod,Lininit,nTP,raw,Vars,'L';SENSDict=SensModel)
    #---------------------- Transformer Modelling ---------------------
    Npij_T, Nqij_T, Npji_T, Nqji_T=Flow_Expr_Simple(LinMethod,Lininit,nTP,raw,Vars,'T')
    # Npij_T1, Nqij_T1, Npji_T1, Nqji_T1=Flow_Expr(LinMethod,Lininit,nTP,raw,Vars,'T';SENSDict=SensModel)
    Np_inj=Power_Inj(nTP,raw,ctdall,Npij_L,Npji_L,Npij_T,Npji_T,LinMethod);
    Nq_inj=Power_Inj(nTP,raw,ctdall,Nqij_L,Nqji_L,Nqij_T,Nqji_T,LinMethod);
    # Np_inj1=Power_Inj(nTP,raw,Npij_L1,Npji_L1,Npij_T1,Npji_T1);
    # Nq_inj1=Power_Inj(nTP,raw,Nqij_L1,Nqji_L1,Nqij_T1,Nqji_T1);
    ####################################################
    PBusset=1:raw["nBus"];QBusset=1:raw["nBus"]
    ####################################################
    maxActive_g = zeros(nTP, length(PBusset))
    minActive_g = zeros(nTP, length(PBusset))
    maxReactive_g = zeros(nTP, length(QBusset))
    minReactive_g = zeros(nTP, length(QBusset))
    if LinMethod=="SENS"
        maxG= ((raw["Gen"]."P_max (pu)")')[ones(Int8,nTP),:]+Matrix(initraw["PG"])
        minG= ((raw["Gen"]."P_min (pu)")')[ones(Int8,nTP),:]+Matrix(initraw["PG"])
        maxQ= ((raw["Gen"]."Q_max (pu)")')[ones(Int8,nTP),:]+Matrix(initraw["QG"])
        minQ= ((raw["Gen"]."Q_min (pu)")')[ones(Int8,nTP),:]+Matrix(initraw["QG"])
        maxfl= isempty(P_dec) ? zeros(nTP,raw["nFL"]) : upper_bound.(P_dec)-lower_bound.(P_inc)+Matrix(initraw["Pdec"])-Matrix(initraw["Pinc"])
        minfl= isempty(P_dec) ? zeros(nTP,raw["nFL"]) : lower_bound.(P_dec)-upper_bound.(P_inc)+Matrix(initraw["Pdec"])-Matrix(initraw["Pinc"])
        maxPess= isempty(P_dis) ? zeros(nTP,raw["nESS"]) : upper_bound.(P_dis)-lower_bound.(P_ch)+Matrix(initraw["Pdis"])-Matrix(initraw["Pch"])
        minPess= isempty(P_dis) ? zeros(nTP,raw["nESS"]) : lower_bound.(P_dis)-upper_bound.(P_ch)+Matrix(initraw["Pdis"])-Matrix(initraw["Pch"])
        maxshed_P= isempty(shed) ? zeros(nTP,raw["nLoad"]) : 0.01*(Proffin["PL"]-Matrix(initraw["Pshed"]).-Lininit["Pshed"])
        maxcurt_P= isempty(curt) ? zeros(nTP,raw["nRES"]) : 0.01*(Proffin["Wind"].-Matrix(initraw["Pcurt"]).-Lininit["Pcurt"])
        maxshed_restore_P= isempty(shed) ? zeros(nTP,raw["nLoad"]) : Matrix(initraw["Pshed"]).+Lininit["Pshed"]
        maxcurt_restore_P= isempty(curt) ? zeros(nTP,raw["nRES"]) : Matrix(initraw["Pcurt"]).+Lininit["Pcurt"]
        maxshed_Q=maxshed_P.*reshape(raw["Load"]."PF" , 1, :)
        maxshed_restore_Q=maxshed_restore_P.*reshape(raw["Load"]."PF" , 1, :)

        # reshape(raw["Load"]."PF" , 1, :)[findall(reshape(raw["Load"]."PF" , 1, :).<0)]
        # minEess=
        # maxEess=
        R_0=  isempty(curt) ? zeros(nTP,raw["nRES"]) : (Proffin["Wind"]-Matrix(initraw["Pcurt"])).-Lininit["Pcurt"]
        ld_0= isempty(shed) ? zeros(nTP,raw["nLoad"]) : Proffin["PL"]-Matrix(initraw["Pshed"]).-Lininit["Pshed"]
        Shunt_0= raw["Shunt"]."bsh0"
        Q_ld_0=isempty(shed) ? -shed : Proffin["QL"]-(Matrix(initraw["Pshed"]).+Lininit["Pshed"]).* reshape(raw["Load"]."PF" , 1, :)
        Np_inj0=value.(Np_inj)
        Nq_inj0=value.(Nq_inj)
        for t in 1:nTP
            for i in PBusset
                maxActive_g[t,i]= sum(maxcurt_restore_P[t,ctdall[i,"res"]];init=0)+sum(maxshed_P[t,ctdall[i,"load"]];init=0)+sum(maxG[t,ctdall[i,"gen"]];init=0)+sum(R_0[t,ctdall[i,"res"]];init=0)+sum(maxfl[t,ctdall[i,"fl"]];init=0)+sum(maxPess[t,ctdall[i,"ess"]];init=0)-sum(ld_0[t,ctdall[i,"load"]];init=0)-Np_inj0[t,i]
                minActive_g[t,i]= -sum(maxcurt_P[t,ctdall[i,"res"]];init=0)-sum(maxshed_restore_P[t,ctdall[i,"load"]];init=0)+sum(minG[t,ctdall[i,"gen"]];init=0)+sum(R_0[t,ctdall[i,"res"]];init=0)+sum(minfl[t,ctdall[i,"fl"]];init=0)+sum(minPess[t,ctdall[i,"ess"]];init=0)-sum(ld_0[t,ctdall[i,"load"]];init=0)-Np_inj0[t,i]
            end
            for i in QBusset
                maxReactive_g[t,i]= sum(maxshed_Q[t,ctdall[i,"load"]];init=0)+sum(maxQ[t,ctdall[i,"gen"]];init=0)-sum(Q_ld_0[t,ctdall[i,"load"]];init=0)+sum(Shunt_0[ctdall[i,"shunt"]];init=0)-Nq_inj0[t,i]
                minReactive_g[t,i]= -sum(maxshed_restore_Q[t,ctdall[i,"load"]];init=0)+sum(minQ[t,ctdall[i,"gen"]];init=0)-sum(Q_ld_0[t,ctdall[i,"load"]];init=0)+sum(Shunt_0[ctdall[i,"shunt"]];init=0)-Nq_inj0[t,i]
            end
        end
    end

    if LinMethod=="SENS"
        Np_inj[:,raw["PActBus"],:]=Np_inj[:,raw["PActBus"],:]+Vars[1]-Vars[2]
        Nq_inj[:,raw["QActBus"],:]=Nq_inj[:,raw["QActBus"],:]+Vars[3]-Vars[4]
    end
    # map_coefficients_inplace!.(x->round(x,digits=7),Np_inj1)
    # map_coefficients_inplace!.(x->round(x,digits=7),Nq_inj1)
    # drop_zeros!.(Np_inj1)
    # drop_zeros!.(Nq_inj1)


    #---------------------- Define Network Constraints-----------------
    if LinMethod!="SENS"
        SCOPF=ConstraintPars(LinMethod,SCOPF,Vars,Parlist,nTP,raw)
    end
    # if LinMethod=="SENS"
    #     PBusset=raw["PActBus"];QBusset=raw["QActBus"]
    # end
    G= isempty(PG) ? PG : (PG+Matrix(initraw["PG"])-curtc)
    R=  isempty(curt) ? -curt : (Proffin["Wind"]-curt-Matrix(initraw["Pcurt"]))
    ess= isempty(P_dis) ? P_dis-P_ch : (P_dis+Matrix(initraw["Pdis"])-P_ch-Matrix(initraw["Pch"]))
    fl= isempty(P_dec) ? P_dec-P_inc : (P_dec+Matrix(initraw["Pdec"])-P_inc-Matrix(initraw["Pinc"]))
    ld= isempty(shed) ? -shed : Proffin["PL"]-shed-Matrix(initraw["Pshed"])
    Shunt= raw["Shunt"]."bsh0"
    Q_ld=isempty(shed) ? -shed : Proffin["QL"]-(shed+Matrix(initraw["Pshed"])).* reshape(raw["Load"]."PF" , 1, :)
    # reshape(raw["Load"]."PF" , 1, :) .* ld
    Q_G=isempty(QG) ? QG : (QG+Matrix(initraw["QG"]))
    if LinMethod!="None"
        Active_g = zeros(AffExpr,nTP, length(PBusset))
        Reactive_g = zeros(AffExpr,nTP, length(QBusset))
    else
        Active_g = zeros(QuadExpr,nTP, length(PBusset))
        Reactive_g = zeros(QuadExpr,nTP, length(QBusset))
    end
    for t in 1:nTP
        for i in PBusset
            Active_g[t,i]= sum(G[t,ctdall[i,"gen"]];init=0)+sum(R[t,ctdall[i,"res"]];init=0)+sum(ess[t,ctdall[i,"ess"]];init=0)+sum(fl[t,ctdall[i,"fl"]];init=0)-sum(ld[t,ctdall[i,"load"]];init=0)-Np_inj[t,i]
        end
        for i in QBusset
            Reactive_g[t,i]= sum(Q_G[t,ctdall[i,"gen"]];init=0)-sum(Q_ld[t,ctdall[i,"load"]];init=0)+sum(Shunt[ctdall[i,"shunt"]];init=0)-Nq_inj[t,i]
        end
    end
    if LinMethod=="SENS"
        Active_g=.!isempty.(getfield.(Active_g, :terms)) .*Active_g
        Reactive_g=.!isempty.(getfield.(Reactive_g, :terms)) .*Reactive_g
    end
    # !isempty.(getfield.(Active_g, :terms))
    # !isempty.(getfield.(Active_g, :constant))
    # index=findall(x->abs(x)>1e-4,getfield.(Active_g, :constant))
    # getfield.(Active_g, :constant)[index]
    # Active_g[18]
    # OP_init0["PG"]
    # Nq_inj
    # Active_g=filter(x->x!=0,Active_g)
    # Reactive_g=filter(x->x!=0,Reactive_g)
    @constraint(SCOPF, NEactive_power_balance, Active_g.==0)
    @constraint(SCOPF, NEreactive_power_balance, Reactive_g.==0)
    # set_normalized_rhs.(NEreactive_power_balance,round.(normalized_rhs.(NEreactive_power_balance),digits=7))
    # set_normalized_rhs.(NEactive_power_balance,round.(normalized_rhs.(NEactive_power_balance),digits=7))
    #-----------Battery constraints
    @constraint(SCOPF,NNbiess[t=1:nTP, i=1:raw["nESS"]], P_ch[t,i]/(raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"])
    + P_dis[t,i]/(raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"])<=1)

    #-----------Coupling constraints
    if Winit==1
        # @constraint(SCOPF,NNramp1[t=1:1,i=1:raw["nGen"]],(PG[t,i])^2<=(raw["Gen"]."Ramp"[i]*TD)^2)
        @constraint(SCOPF,NNrampi1[t=1:1,i=1:raw["nGen"]],PG[t,i]+(initraw["PG"][t,i]-initial["PG"][1,i])<=raw["Gen"]."Ramp"[i]*TD)
        @constraint(SCOPF,NNrampi2[t=1:1,i=1:raw["nGen"]],-raw["Gen"]."Ramp"[i]*TD<=PG[t,i]+(initraw["PG"][t,i]-initial["PG"][1,i]))

    end
        #ESS
    @constraint(SCOPF, NEsoc0[t=1:1, i=1:raw["nESS"]], soc[t,i]==initial["soc"][1,i]
                                            + TD* ((P_ch[t,i]+initraw["Pch"][t,i])*raw["ESS"]."Charge_efficiency"[i] - (P_dis[t,i]+initraw["Pdis"][t,i])/(raw["ESS"]."Discharge_efficiency"[i]) ))
    @constraint(SCOPF, NEsocn[t=2:nTP, i=1:raw["nESS"]], soc[t,i]==soc[t-1,i]
                                            + TD* ((P_ch[t,i]+initraw["Pch"][t,i])*raw["ESS"]."Charge_efficiency"[i] - (P_dis[t,i]+initraw["Pdis"][t,i])/(raw["ESS"]."Discharge_efficiency"[i]) ))
        #Gens
    if relaxpermit==true && LinMethod=="SENS"
        # @variable(SCOPF, slk_NNramp1[t=1:nTP,i=1:raw["nGen"]] >= 0)
        # @variable(SCOPF, slk_NNramp2[t=1:nTP,i=1:raw["nGen"]] >= 0)
        @variable(SCOPF, 1>=slk_NNramp[t=1:nTP,i=1:raw["nGen"]] >= 0)
        # @constraint(SCOPF, slk_NNramp .>= (slk_NNramp1))
        # @constraint(SCOPF, slk_NNramp .>= (slk_NNramp2))
    end
    # slk_NNramp.-(slk_NNramp1)
       # -----------Coupling constraints
    @constraint(SCOPF,NNramp1[t=2:nTP,i=1:raw["nGen"]],PG[t,i]-PG[t-1,i]+(initraw["PG"][t,i]-initraw["PG"][t-1,i])<=raw["Gen"]."Ramp"[i]*TD+((relaxpermit==true && LinMethod=="SENS") ? (slk_NNramp[t,i]* max.(0,SCOPF[:Rampvio][1][t,i])) : (0)))
    @constraint(SCOPF,NNramp2[t=2:nTP,i=1:raw["nGen"]],-raw["Gen"]."Ramp"[i]*TD-((relaxpermit==true && LinMethod=="SENS") ? (slk_NNramp[t,i]* max.(0,-SCOPF[:Rampvio][1][t,i])) : (0))<=PG[t,i]-PG[t-1,i]+(initraw["PG"][t,i]-initraw["PG"][t-1,i]))

    # @constraint(SCOPF,NNramp1[t=2:nTP,i=1:raw["nGen"]],PG[t,i]-PG[t-1,i]+(initraw["PG"][t,i]-initraw["PG"][t-1,i])<=raw["Gen"]."Ramp"[i]*TD)
    # @constraint(SCOPF,NNramp2[t=2:nTP,i=1:raw["nGen"]],-raw["Gen"]."Ramp"[i]*TD<=PG[t,i]-PG[t-1,i]+(initraw["PG"][t,i]-initraw["PG"][t-1,i]))
        ##fixed Energy all_constraints
    if fixed_SOC==1
        @constraint(SCOPF, NEfixSoc1[i=1:raw["nESS"]], soc[nTP,i]==initial["soc"][1,i] )
        @constraint(SCOPF,NEfl1[i=1:raw["nFL"]],sum(P_inc[t,i] for t in 1:nTP)==sum(P_dec[t,i] for t in 1:nTP))
    end
#     elseif fixed_SOC==2
#         @constraint(SCOPF, NEfixSoc2[i=1:raw["nESS"]], sum(((P_ch[t,i])*raw["ESS"]."Charge_efficiency"[i] - (P_dis[t,i])/(raw["ESS"]."Discharge_efficiency"[i])) for t=1:nTP)==0)
#         @constraint(SCOPF,NEfl2[i=1:raw["nFL"]],sum((P_inc[t,i]-P_dec[t,i]) for t in 1:nTP)==0)
#     elseif fixed_SOC==3
#         # @constraint(SCOPF, NEfixSoc2[i=1:raw["nESS"]], sum(((P_ch[t,i])*raw["ESS"]."Charge_efficiency"[i] - (P_dis[t,i])/(raw["ESS"]."Discharge_efficiency"[i])) for t=1:nTP)==0)
#         @constraint(SCOPF, NEfixSoc3[i=1:raw["nESS"]], soc[nTP,i]<=initraw["soc"][nTP,i]+ 0.1(raw["ESS"]."E_rating (MWh)"[i]/raw["sbase"]-initraw["soc"][nTP,i]))
#         @constraint(SCOPF, NEfixSoc4[i=1:raw["nESS"]], soc[nTP,i]>=initraw["soc"][nTP,i]- 0.1(initraw["soc"][nTP,i]-raw["ESS"]."E_rating_min (MWh)"[i]/raw["sbase"]))
# ### I don't know what I am doing here
#         @constraint(SCOPF,NEfl2[i=1:raw["nFL"]],sum((P_inc[t,i]-P_dec[t,i]) for t in 1:nTP)<= (sum(initraw["Pinc"][t,i]-initraw["Pdec"][t,i] for t=1:nTP)))
#         @constraint(SCOPF,NEfl3[i=1:raw["nFL"]],sum((P_inc[t,i]-P_dec[t,i]) for t in 1:nTP)>= -1*(sum(initraw["Pinc"][t,i]-initraw["Pdec"][t,i] for t=1:nTP)))
#     end


    #FL
    @constraint(SCOPF, NEUE0[t=1:1, i=1:raw["nFL"]], UE[t,i]== initial["UE"][1,i]
                                            + TD* ( (P_inc[t,i]+initraw["Pinc"][t,i]) - (P_dec[t,i]+initraw["Pdec"][t,i]) ))
    @constraint(SCOPF, NEUEn[t=2:nTP, i=1:raw["nFL"]], UE[t,i]==UE[t-1,i]
                                            + TD* ((P_inc[t,i]+initraw["Pinc"][t,i]) - (P_dec[t,i]+initraw["Pdec"][t,i])))

    if Inc_coef>=0.01 && Dec_coef>=0.01
        @constraint(SCOPF,NNfl[t=1:nTP,i=1:raw["nFL"]],P_inc[t,i]/(Inc_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]])+P_dec[t,i]/(Dec_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]]) <=1)
    end
     #-----------Coupling constraints
    N_Var= PG,QG,P_ch,P_dis,P_inc,P_dec,shunt,shed,curt,curtc,soc,UE,Vars,Parlist,(maxActive_g,minActive_g,maxReactive_g,minReactive_g)
    return SCOPF,N_Var
end


function Contingency_L_PB( SCOPF::Model,nTP::Int64,raw::Dict{String, Any},initraw::Dict{String, DataFrame},ctdall,fixed_SOC::Int64,N_Var::Tuple,Proffin::Dict{String, Array{Float64}},CCT::Float64,initial::Dict{String, DataFrame},LSpermit::Int64,LinMethod::String,Lininit::Dict{String,Any},relaxconfac=0.0) #Array{Union{Float64,Array{Float64}}}

    # SCOPF,nTP,raw,initraw,fixed_SOC,N_Var,Proffin,CCT,initial,LSpermit,LinMethod,Lininit=LSCOPF,nTP,raw,initraw,fixed_SOC,N_Var,Proffin,CCT,initial,LSpermit,LinMethod,Lininit
    #----------------------Initializing Model----------------------
    PG,QG,P_ch,P_dis,P_inc,P_dec,NNNN,shed,NNNN,curtc,soc,UE,N_IVars=N_Var
    SCOPF,Vars=Define_Variable_Contingency_L(SCOPF,nTP,raw,LinMethod,Lininit)
    @variables(
        SCOPF,
        begin
            PG_c_Line[t=1:nTP, i = 1:raw["nGen"],cL = 1:raw["nConL"]]
            QG_c_Line[t=1:nTP, i = 1:raw["nGen"],cL = 1:raw["nConL"]]
            P_ch_c_Line[t=1:nTP,i= 1:raw["nESS"],cL = 1:raw["nConL"]]
            P_dis_c_Line[t=1:nTP,i= 1:raw["nESS"],cL = 1:raw["nConL"]]
            soc_c_Line[t=1:nTP,i= 1:raw["nESS"],cL = 1:raw["nConL"]]
            P_inc_c_Line[t = 1:nTP, i = 1:raw["nFL"],cL = 1:raw["nConL"]]
            P_dec_c_Line[t = 1:nTP, i = 1:raw["nFL"],cL = 1:raw["nConL"]]
            shed_c_Line[t = 1:nTP, i = 1:raw["nLoad"],cL = 1:raw["nConL"]]
            curt_c_Line[t = 1:nTP, i = 1:raw["nRES"],cL = 1:raw["nConL"]]
            curtc_c_Line[t = 1:nTP, i = 1:raw["nGen"],cL = 1:raw["nConL"]]
        end
    )
    # for  t=1:nTP, i=1:raw["nRES"],cL=1:raw["nConL"]     fix(curt_c_Line[t,i,cL], 0)  end
    # for  t=1:nTP, i=1:raw["nLoad"],cL=1:raw["nConL"]    fix(shed_c_Line[t,i,cL], 0) end
    if LSpermit==1
        @constraint(SCOPF,CNcurtlim1[t=1:nTP,i=1:raw["nRES"],cL=1:raw["nConL"]],curt_c_Line[t,i,cL]<=Proffin["Wind"][t,i]-initraw["Pcurt"][t,i])
        @constraint(SCOPF,CNcurtlim2[t=1:nTP,i=1:raw["nRES"],cL=1:raw["nConL"]],-initraw["Pcurt"][t,i] <=curt_c_Line[t,i,cL])
        @constraint(SCOPF,CNshedlim1[t=1:nTP,i=1:raw["nLoad"],cL=1:raw["nConL"]],shed_c_Line[t,i,cL]<=Proffin["PL"][t,i]-initraw["Pshed"][t,i])
        @constraint(SCOPF,CNshedlim2[t=1:nTP,i=1:raw["nLoad"],cL=1:raw["nConL"]],-initraw["Pshed"][t,i] <=shed_c_Line[t,i,cL])
        @constraint(SCOPF,CNcurtclim1[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],-initraw["Pcurtc"][t,i]<=curtc_c_Line[t,i,cL])
        @constraint(SCOPF,CNcurtclim2[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],curtc_c_Line[t,i,cL]<=PG_c_Line[t,i,cL]+initraw["PG"][t,i]-initraw["Pcurtc"][t,i])
    else
        @constraint(SCOPF,CNcurtlim1[t=1:nTP,i=1:raw["nRES"],cL=1:raw["nConL"]],curt_c_Line[t,i,cL]<=Proffin["Wind"][t,i]-initraw["Pcurt"][t,i])
        @constraint(SCOPF,CNcurtlim2[t=1:nTP,i=1:raw["nRES"],cL=1:raw["nConL"]],-initraw["Pcurt"][t,i] <=curt_c_Line[t,i,cL])
        # @constraint(SCOPF,CNcurtlim[t=1:nTP,i=1:raw["nRES"],cL=1:raw["nConL"]],curt_c_Line[t,i,cL]<=0)
        @constraint(SCOPF,CNshedlim[t=1:nTP,i=1:raw["nLoad"],cL=1:raw["nConL"]],shed_c_Line[t,i,cL]==shed[t,i])
        @constraint(SCOPF,CNcurtclim[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],curtc_c_Line[t,i,cL]==curtc[t,i])

    end

    for  t=1:nTP, i=1:raw["nRES"],cL=1:raw["nConL"]     set_start_value(curt_c_Line[t,i,cL], 0)  end                       # Setting start values
    for  t=1:nTP, i=1:raw["nLoad"],cL=1:raw["nConL"]    set_start_value(shed_c_Line[t,i,cL], 0) end

    # for t=1:nTP, i=1:raw["nBus"],cL=1:raw["nConL"]
    #     # fix(sl_V[t,i,cL], 0.0)
    #     # set_lower_bound(sl_V[t,i,cL], 0.0)
    #     # set_start_value(sl_V[t,i,cL], 0.0)
    # end
    # for t=1:nTP, i=1:raw["nLine"],cL=1:raw["nConL"]
    #     # fix(sl_F[t,i,cL], 0.0)
    #     # set_lower_bound(sl_F[t,i,cL], 0.0)
    #     # set_start_value(sl_F[t,i,cL], 0.0)
    # end

    @constraint(SCOPF,CNplim1[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],PG_c_Line[t,i,cL]+initraw["PG"][t,i]<=raw["Gen"]."P_max (pu)"[i])
    @constraint(SCOPF,CNplim2[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],raw["Gen"]."P_min (pu)"[i]<=PG_c_Line[t,i,cL]+initraw["PG"][t,i])
    @constraint(SCOPF,CNqlim1[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],QG_c_Line[t,i,cL]+initraw["QG"][t,i]<=raw["Gen"]."Q_max (pu)"[i])
    @constraint(SCOPF,CNqlim2[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],raw["Gen"]."Q_min (pu)"[i]<=QG_c_Line[t,i,cL]+initraw["QG"][t,i])

    for t=1:nTP, i=1:raw["nGen"],cL=1:raw["nConL"]
        set_start_value(PG_c_Line[t,i,cL], 0)
        set_start_value(QG_c_Line[t,i,cL], 0)
        fix(curtc_c_Line[t,i,cL], 0)
    end

    @constraint(SCOPF,CNbchlim1[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],P_ch_c_Line[t,i,cL]<=raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"])
    @constraint(SCOPF,CNbchlim2[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],0<=P_ch_c_Line[t,i,cL])
    @constraint(SCOPF,CNbdislim1[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],P_dis_c_Line[t,i,cL]<=raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"])
    @constraint(SCOPF,CNbdislim2[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],0<=P_dis_c_Line[t,i,cL])
    for  t=1:nTP, i=1:raw["nESS"], cL=1: raw["nConL"]
        set_start_value(P_ch_c_Line[t,i,cL], 0.5*raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"])
        set_start_value(P_dis_c_Line[t,i,cL], 0.5*raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"])
        set_start_value(soc_c_Line[t,i,cL], 0.5*raw["ESS"]."E_rating (MWh)"[i]/raw["sbase"]);set_lower_bound(soc_c_Line[t,i,cL], raw["ESS"]."E_rating_min (MWh)"[i]/raw["sbase"]);set_upper_bound(soc_c_Line[t,i,cL], raw["ESS"]."E_rating (MWh)"[i]/raw["sbase"]);
    end

    @constraint(SCOPF,CNfilim1[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],P_inc_c_Line[t,i,cL]+initraw["Pinc"][t,i] <=Inc_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]])
    @constraint(SCOPF,CNfilim2[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],0<=P_inc_c_Line[t,i,cL]+initraw["Pinc"][t,i])
    @constraint(SCOPF,CNfdlim1[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],P_dec_c_Line[t,i,cL]+initraw["Pdec"][t,i] <=Dec_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]])
    @constraint(SCOPF,CNfdlim2[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],0<=P_dec_c_Line[t,i,cL]+initraw["Pdec"][t,i])
    for t=1:nTP, i= 1:raw["nFL"],cL=1:raw["nConL"]
        set_start_value(P_inc_c_Line[t,i,cL], 0)                         # Setting start values
        set_start_value(P_dec_c_Line[t,i,cL], 0)
    end

    # if LinMethod !="SENS"
    #     SENSDict=Dict{String, Matrix{Matrix{Float64}}}()
    # elseif LinMethod =="SENS"
    #    SENSDict=Cal_SENS_new(Lininit,ctdall,nTP,raw;C="CL");
    # end
    #-----------------Line flow Modeling  #OK
     CLpij_L, CLqij_L, CLpji_L, CLqji_L, Parlist=Flow_Expr_Simple(LinMethod,Lininit,nTP,raw,Vars,'L';Cont=raw["ConL"]);
    #---------------------- Transformer Modelling ---------------------
     CLpij_T, CLqij_T, CLpji_T, CLqji_T=Flow_Expr_Simple(LinMethod,Lininit,nTP,raw,Vars,'T';Cont=raw["ConL"]);
     CLp_inj=Power_Inj(nTP,raw,ctdall,CLpij_L,CLpji_L,CLpij_T,CLpji_T,LinMethod;C="CL");
     CLq_inj=Power_Inj(nTP,raw,ctdall,CLqij_L,CLqji_L,CLqij_T,CLqji_T,LinMethod;C="CL");

     PBusset= 1:raw["nBus"];QBusset= 1:raw["nBus"]

     maxActive_g = zeros(nTP, length(PBusset),raw["nConL"])
     minActive_g = zeros(nTP, length(PBusset),raw["nConL"])
     maxReactive_g = zeros(nTP, length(QBusset),raw["nConL"])
     minReactive_g = zeros(nTP, length(QBusset),raw["nConL"])
    if LinMethod =="SENS"
        maxG= ((raw["Gen"]."P_max (pu)")')[ones(Int8,nTP),:]+Matrix(initraw["PG"])
        minG= ((raw["Gen"]."P_min (pu)")')[ones(Int8,nTP),:]+Matrix(initraw["PG"])
        maxQ= ((raw["Gen"]."Q_max (pu)")')[ones(Int8,nTP),:]+Matrix(initraw["QG"])
        minQ= ((raw["Gen"]."Q_min (pu)")')[ones(Int8,nTP),:]+Matrix(initraw["QG"])
        # maxfl= isempty(P_dec_c_Line) ? zeros(nTP,raw["nFL"],raw["nConL"]) : upper_bound.(P_dec_c_Line)-lower_bound.(P_inc_c_Line)+Matrix(initraw["Pdec"])-Matrix(initraw["Pinc"])
        # minfl= isempty(P_dec_c_Line) ? zeros(nTP,raw["nFL"],raw["nConL"]) : lower_bound.(P_dec_c_Line)-upper_bound.(P_inc_c_Line)+Matrix(initraw["Pdec"])-Matrix(initraw["Pinc"])
        maxfl= isempty(P_dec_c_Line) ? zeros(nTP,raw["nFL"]) : Dec_coef*Proffin["PL"][:,raw["FL"]."iFL"] .- Matrix(initraw["Pdec"]) - (-1 .* Matrix(initraw["Pinc"])) +Matrix(initraw["Pdec"])-Matrix(initraw["Pinc"])
        minfl= isempty(P_dec_c_Line) ? zeros(nTP,raw["nFL"]) : -1 .* Matrix(initraw["Pdec"]) - (Inc_coef*Proffin["PL"][:,raw["FL"]."iFL"] .- Matrix(initraw["Pinc"])) +Matrix(initraw["Pdec"])-Matrix(initraw["Pinc"])
         
        # minPess= isempty(P_dis_c_Line) ? zeros(nTP,raw["nESS"],raw["nConL"]) : upper_bound.(P_dis_c_Line)-lower_bound.(P_ch_c_Line)+Matrix(initraw["Pdis"])-Matrix(initraw["Pch"])
        # maxPess= isempty(P_dis_c_Line) ? zeros(nTP,raw["nESS"],raw["nConL"]) : lower_bound.(P_dis_c_Line)-upper_bound.(P_ch_c_Line)+Matrix(initraw["Pdis"])-Matrix(initraw["Pch"])
        maxPess= isempty(P_dis_c_Line) ? zeros(nTP,raw["nESS"]) : (raw["ESS"]."Discharge_rating (MW)")' ./raw["sbase"] .- Matrix(initraw["Pdis"])-(-1 .* Matrix(initraw["Pch"]))+Matrix(initraw["Pdis"])-Matrix(initraw["Pch"])
        minPess= isempty(P_dis_c_Line) ? zeros(nTP,raw["nESS"]) : -1 .* Matrix(initraw["Pdis"]).-(raw["ESS"]."Charge_rating (MW)")' ./ raw["sbase"] .- Matrix(initraw["Pch"]) +Matrix(initraw["Pdis"])-Matrix(initraw["Pch"])

        maxshed_P= isempty(shed_c_Line) ? zeros(nTP,raw["nLoad"],raw["nConL"]) : 0.01*(Proffin["PL"].-Matrix(initraw["Pshed"]).-Lininit["Pshed_c_Line"])
        maxcurt_P= isempty(curt_c_Line) ? zeros(nTP,raw["nRES"],raw["nConL"]) :  0.01*(Proffin["Wind"].-Matrix(initraw["Pcurt"]).-Lininit["Pcurt_c_Line"])
        maxshed_restore_P= isempty(shed_c_Line) ? zeros(nTP,raw["nLoad"],raw["nConL"]) : Matrix(initraw["Pshed"]).+Lininit["Pshed_c_Line"]
        maxcurt_restore_P= isempty(curt_c_Line) ? zeros(nTP,raw["nRES"],raw["nConL"]) : Matrix(initraw["Pcurt"]).+Lininit["Pcurt_c_Line"]
        maxshed_Q=maxshed_P.*reshape(raw["Load"]."PF" , 1, :)
        maxshed_restore_Q= maxshed_restore_P.*reshape(raw["Load"]."PF" , 1, :)
        # minEess=
        # maxEess=
        R_0=  isempty(curt_c_Line) ? zeros(nTP,raw["nRES"],raw["nConL"]) : (Proffin["Wind"]-Matrix(initraw["Pcurt"])).-Lininit["Pcurt_c_Line"]
        ld_0= isempty(shed_c_Line) ? zeros(nTP,raw["nLoad"],raw["nConL"]) : Proffin["PL"]-Matrix(initraw["Pshed"]).-Lininit["Pshed_c_Line"]
        Shunt_0= raw["Shunt"]."bsh0"
        # Q_ld_0=reshape(raw["Load"]."PF" , 1, :) .* ld_0
        Q_ld_0=isempty(shed_c_Line) ? -shed : Proffin["QL"].-(Matrix(initraw["Pshed"]).+Lininit["Pshed_c_Line"]).* reshape(raw["Load"]."PF" , 1, :)


        CLp_inj0=value.(CLp_inj)
        CLq_inj0=value.(CLq_inj)

        for t in 1:nTP,cL in 1:raw["nConL"]
            for i in PBusset
                maxActive_g[t,i,cL]= sum(maxcurt_restore_P[t,ctdall[i,"res"],cL];init=0)+sum(maxshed_P[t,ctdall[i,"load"],cL];init=0)+sum(maxG[t,ctdall[i,"gen"]];init=0)+sum(R_0[t,ctdall[i,"res"],cL];init=0)+sum(maxPess[t,ctdall[i,"ess"]];init=0)+sum(maxfl[t,ctdall[i,"fl"]];init=0)-sum(ld_0[t,ctdall[i,"load"],cL];init=0)-CLp_inj0[t,i,cL] 
                minActive_g[t,i,cL]= -sum(maxcurt_P[t,ctdall[i,"res"],cL];init=0)-sum(maxshed_restore_P[t,ctdall[i,"load"],cL];init=0)+sum(minG[t,ctdall[i,"gen"]];init=0)+sum(R_0[t,ctdall[i,"res"],cL];init=0)+sum(minPess[t,ctdall[i,"ess"]];init=0)+sum(minfl[t,ctdall[i,"fl"]];init=0)-sum(ld_0[t,ctdall[i,"load"],cL];init=0)-CLp_inj0[t,i,cL]
            end
            for i in QBusset
                maxReactive_g[t,i,cL]= sum(maxshed_Q[t,ctdall[i,"load"],cL];init=0)+sum(maxQ[t,ctdall[i,"gen"]];init=0)-sum(Q_ld_0[t,ctdall[i,"load"],cL];init=0)+sum(Shunt_0[ctdall[i,"shunt"]];init=0)-CLq_inj0[t,i,cL]
                minReactive_g[t,i,cL]= -sum(maxshed_restore_Q[t,ctdall[i,"load"],cL];init=0)+sum(minQ[t,ctdall[i,"gen"]];init=0)-sum(Q_ld_0[t,ctdall[i,"load"],cL];init=0)+sum(Shunt_0[ctdall[i,"shunt"]];init=0)-CLq_inj0[t,i,cL]
            end
        end
    end


     if LinMethod=="SENS"
        CLp_inj[:,raw["PActBus"],:]=CLp_inj[:,raw["PActBus"],:]+Vars[1]-Vars[2]
        CLq_inj[:,raw["QActBus"],:]=CLq_inj[:,raw["QActBus"],:]+Vars[3]-Vars[4]
     end

    #  map_coefficients_inplace!.(x->round(x,digits=7),CLp_inj)
    #  map_coefficients_inplace!.(x->round(x,digits=7),CLq_inj)
    #  drop_zeros!.(CLp_inj)
    #  drop_zeros!.(CLq_inj)
    #-------------------------- Network constraints--------------------
    if LinMethod!="SENS"
        SCOPF = ConstraintPars(LinMethod,SCOPF,Vars,Parlist,nTP,raw,relaxconfac;States="Cont");
    end
    # if LinMethod=="SENS"
    #     PBusset=raw["PActBus"];QBusset=raw["QActBus"]
    # end
    G= isempty(PG_c_Line) ? PG_c_Line : (PG_c_Line.+Matrix(initraw["PG"])-curtc_c_Line)
    R=  isempty(curt_c_Line) ? -curt_c_Line : (Proffin["Wind"].-curt_c_Line.-Matrix(initraw["Pcurt"]))
    ess= isempty(P_dis_c_Line) ? P_dis_c_Line-P_ch_c_Line : (P_dis_c_Line.+Matrix(initraw["Pdis"])-P_ch_c_Line.-Matrix(initraw["Pch"]))
    fl= isempty(P_dec_c_Line) ? P_dec_c_Line-P_inc_c_Line : (P_dec_c_Line.+Matrix(initraw["Pdec"])-P_inc_c_Line.-Matrix(initraw["Pinc"]))
    ld= isempty(shed_c_Line) ? -shed_c_Line : Proffin["PL"].-shed_c_Line.-Matrix(initraw["Pshed"])
    Shunt= raw["Shunt"]."bsh0"
    Q_ld=isempty(shed_c_Line) ? -shed_c_Line : Proffin["QL"].-(shed_c_Line.+Matrix(initraw["Pshed"])).* reshape(raw["Load"]."PF" , 1, :)
    # Q_ld=reshape(raw["Load"]."PF" , 1, :) .* ld
    Q_G=isempty(QG_c_Line) ? QG_c_Line : (QG_c_Line.+Matrix(initraw["QG"]))

    if LinMethod!="None"
        Active_g = zeros(AffExpr,nTP, length(PBusset),raw["nConL"])
        Reactive_g = zeros(AffExpr,nTP, length(QBusset),raw["nConL"])
    else
        Active_g = zeros(QuadExpr,nTP, length(PBusset),raw["nConL"])
        Reactive_g = zeros(QuadExpr,nTP, length(QBusset),raw["nConL"])
    end
    for t in 1:nTP,cL in 1:raw["nConL"]
        for i in PBusset
            Active_g[t,i,cL]= sum(G[t,ctdall[i,"gen"],cL];init=0)+sum(R[t,ctdall[i,"res"],cL];init=0)+sum(ess[t,ctdall[i,"ess"],cL];init=0)+sum(fl[t,ctdall[i,"fl"],cL];init=0)-sum(ld[t,ctdall[i,"load"],cL];init=0)-CLp_inj[t,i,cL]
        end
        for i in QBusset
            Reactive_g[t,i,cL]= sum(Q_G[t,ctdall[i,"gen"],cL];init=0)-sum(Q_ld[t,ctdall[i,"load"],cL];init=0)+sum(Shunt[ctdall[i,"shunt"]];init=0)-CLq_inj[t,i,cL]
        end
    end
    if LinMethod=="SENS"
        Active_g=.!isempty.(getfield.(Active_g, :terms)) .*Active_g
        Reactive_g=.!isempty.(getfield.(Reactive_g, :terms)) .*Reactive_g
    end
    @constraint(SCOPF, CEActive_power_balance, Active_g.==0)
    @constraint(SCOPF, CEReactive_power_balance, Reactive_g.==0)
    # fix.(SCOPF[:ΔQinj_c_Line][10,13,14],0)

    @assert(Mode in ["Preventive","Corrective"])
    if Mode=="Preventive"
        @constraint(SCOPF,CEfixp[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[raw["Gen"]."Bus"[i]] in [1,2]],PG_c_Line[t,i,cL]==PG[t,i])
        # @constraint(SCOPF,CEfixQ[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[raw["Gen"]."Bus"[i]]!=3],QG_c_Line[t,i,cL]==QG[t,i])
        @constraint(SCOPF,CEfixESS1[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],P_dis_c_Line[t,i,cL]==P_dis[t,i])
        @constraint(SCOPF,CEfixESS2[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],P_ch_c_Line[t,i,cL]==P_ch[t,i])
        @constraint(SCOPF,CEfixFL1[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],P_inc_c_Line[t,i,cL]==P_inc[t,i])
        @constraint(SCOPF,CEfixFL2[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],P_dec_c_Line[t,i,cL]==P_dec[t,i])
        @constraint(SCOPF,CEfixedvoltage[t=1:nTP,i=1:raw["nBus"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[i] in [2,3]],SCOPF[:CV2magexp][t,i,cL]==SCOPF[:NV2magexp][t,i])
    elseif Mode =="Corrective"
        if relaxpermit==true && LinMethod=="SENS"
            # @variable(SCOPF, slk_CNramp1[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]] >= 0)
            # @variable(SCOPF, slk_CNramp2[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]] >= 0)
            @variable(SCOPF, slk_CNramp[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]] >= 0)
            # reshapedslk=reshape(slk_CNramp, size(slk_CNramp,1), 1, size(slk_CNramp,2)) 
            # @constraint(SCOPF, reshapedslk .>= (slk_CNramp1))
            # @constraint(SCOPF, reshapedslk .>= (slk_CNramp2))
        end
        #-----------Coupling constraints
        @constraint(SCOPF,CNramp1[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],(PG_c_Line[t,i,cL]-PG[t,i])<=(raw["Gen"]."Ramp"[i]*CCT) +((relaxpermit==true && LinMethod=="SENS") ? (slk_CNramp[t,i,cL]* max.(0,SCOPF[:Rampvio][2][t,i,cL])) : (0)) )
        @constraint(SCOPF,CNramp2[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],(PG_c_Line[t,i,cL]-PG[t,i])>=(-raw["Gen"]."Ramp"[i]*CCT) - ((relaxpermit==true && LinMethod=="SENS") ? (slk_CNramp[t,i,cL]* max.(0,-SCOPF[:Rampvio][2][t,i,cL])) : (0)) )
        # @constraint(SCOPF,CNramp1[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],(PG_c_Line[t,i,cL]-PG[t,i])<=(raw["Gen"]."Ramp"[i]*CCT))
        # @constraint(SCOPF,CNramp2[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],(PG_c_Line[t,i,cL]-PG[t,i])>=(-raw["Gen"]."Ramp"[i]*CCT))

        # # @constraint(SCOPF,CNramp[t=2:nTP,i=1:raw["nGen"],cL=1: raw["nConL"]],-raw["Gen"]."Ramp"[i]*TD<=PG_c_Line[t,i,cL]-PG_c_Line[t-1,i,cL]+initraw["PG"][t,i]-initraw["PG"][t-1,i]<=raw["Gen"]."Ramp"[i]*TD)
        #-----------Battery constraints
        @constraint(SCOPF,CNbiess[t=1:nTP, i=1:raw["nESS"],cL=1:raw["nConL"]], (P_ch_c_Line[t,i,cL])/(raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"])
        + (P_dis_c_Line[t,i,cL])/(raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"])<=1)
        @constraint(SCOPF, CEsoc0[t=1:1, i=1:raw["nESS"],cL=1:raw["nConL"]], soc_c_Line[t,i,cL]==initial["soc"][1,i]
                + CCT* ((P_ch_c_Line[t,i,cL])*raw["ESS"]."Charge_efficiency"[i] - (P_dis_c_Line[t,i,cL])/(raw["ESS"]."Discharge_efficiency"[i])) )
        @constraint(SCOPF, CEsocn[t=2:nTP, i=1:raw["nESS"],cL=1:raw["nConL"]], soc_c_Line[t,i,cL]==soc[t-1,i]
                + CCT* ((P_ch_c_Line[t,i,cL])*raw["ESS"]."Charge_efficiency"[i] - (P_dis_c_Line[t,i,cL])/(raw["ESS"]."Discharge_efficiency"[i])) )
        # if fixed_SOC==1 @constraint(SCOPF, CEfixsoc[i=1:raw["nESS"],cL=1:raw["nConL"]], soc_c_Line[nTP,i,cL]==initial["soc"][1,i] ) end
        #FL
        # @constraint(SCOPF,CEfl[i=1:raw["nFL"],cL=1:raw["nConL"]],sum(P_inc_c_Line[t,i,cL] for t in 1:nTP)==sum(P_dec_c_Line[t,i,cL] for t in 1:nTP))
        if Inc_coef>=0.01 && Dec_coef>=0.01
            @constraint(SCOPF,CNfl[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],(P_inc_c_Line[t,i,cL])/(Inc_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]])+(P_dec_c_Line[t,i,cL])/(Dec_coef*Proffin["PL"][t,raw["FL"]."iFL"[i]]) <=1)
        end
    end

    C_Var= PG_c_Line,QG_c_Line,P_ch_c_Line,P_dis_c_Line,P_inc_c_Line,P_dec_c_Line,shed_c_Line,curt_c_Line,curtc_c_Line,Vars,Parlist,(maxActive_g,minActive_g,maxReactive_g,minReactive_g)
    return SCOPF,C_Var
end

function find_ramp_rate_violations(Lininit,raw,OPFon)
    if size(Lininit["PG"],1)>1
        PreT = vcat(Lininit["PG"][2:end,:], Lininit["PG"][end:end,:])  # Append row of zeros
        nextT=vcat( Lininit["PG"][1:1,:],Lininit["PG"][1:end-1,:])

        adiff=Lininit["PG"]-nextT
        # bdiff=Lininit["PG"]-PreT
        Base_TR_N=round.( sign.(adiff).*max.(0, abs.(adiff).- Td*(raw["Gen"].Ramp)'),digits=7)
        # Base_TR_N=max.(abs.(adiff).- Td*(raw["Gen"].Ramp)',abs.(bdiff).- Td*(raw["Gen"].Ramp)',abs.(adiff+bdiff).-(2*Td*raw["Gen"].Ramp)',0)
        if OPFon==0
            cdiff=Lininit["PG_c_Line"].-Lininit["PG"]

            Base_TR_C=round.(sign.(cdiff).*max.(0, abs.(cdiff).-CCT*(raw["Gen"].Ramp)'),digits=7)
        else
            Base_TR_C=zeros(size(Lininit["PG"]))
        end
    else
        Base_TR_N=zeros(size(Lininit["PG"]))
        Base_TR_C=zeros(size(Lininit["PG"]))
    end
    return Base_TR_N,Base_TR_C
end
