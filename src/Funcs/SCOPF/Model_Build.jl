function Initialize_Model()

    model = Model(Ipopt.Optimizer)
    if(printlevel==0) set_optimizer_attribute(model, "print_level",0) end
    # set_attribute(model, "linear_solver", "ma97")
    # set_optimizer_attributes(model,"tol"=>1e-04)
    # set_optimizer_attributes(model,"mu_strategy"=>"adaptive")
    # set_optimizer_attributes(model,"mumps_mem_percent"=>100)
    return model
end
function Solve_P(SCOPF)
    optimize!(SCOPF)
    println("Objective value:  ", JuMP.objective_value(SCOPF))
    println("Solver Time:   ", JuMP.solve_time(SCOPF))
    return SCOPF
end

function Objective_function(SCOPF,nTP,raw,initraw,Proffin,N_Var,C_Var,sl_F,sl_V,t0,Lccost=1e4,Rccost=0;lambda=0.05,OPFon=0,singleconprob=1/33,WP=1)
#####################################
#####################################
Lccost=raw["LScost"]*raw["sbase"]
Rccost=raw["RCcost"]*raw["sbase"]

FactorT=max(1,500*(t0-132))
λG=0.005*FactorT*raw["Gen"]."cost_b"./(2)
λS=0.005*FactorT*((144)^2)*raw["ESS"]."cost_b"./(2*TD)
λF=0.005*FactorT*((144)^2)*raw["Load"]."ShCost"[raw["FL"]."Bus"]./(2*TD)
    #---------------- convert pu to MW for cost function
#real cost of Lc and Rc = value/sbase

    #---------------- convert pu to MW for cost function
PG,gar,P_ch,P_dis,P_inc,P_dec,gar,shed,curt,curtc,soc,UE =N_Var
PG_c_Line,gar,P_ch_c_Line,P_dis_c_Line,P_inc_c_Line,P_dec_c_Line,shed_c_Line,curt_c_Line,curtc_c_Line =C_Var
sl_F_cost,sl_V_cost=1e5,1e5
Prob_c= zeros(raw["nConL"]).+singleconprob
#Preventive
@expression(SCOPF,cost_Gen[t=1:nTP,i=1:raw["nGen"]], raw["Gen"]."cost_c"[i]+ raw["Gen"]."cost_b"[i]*raw["sbase"]*(PG[t,i])+raw["Gen"]."cost_a"[i]*(raw["sbase"]*(PG[t,i]))^2)
@expression(SCOPF,cost_ESS[t=1:nTP,i=1:raw["nESS"]], raw["ESS"]."cost_c"[i]+ raw["ESS"]."cost_b"[i]*raw["sbase"]*(P_dis[t,i]+P_ch[t,i])+raw["ESS"]."cost_a"[i]*raw["sbase"]^2*(P_dis[t,i]+P_ch[t,i])*(P_dis[t,i]+P_ch[t,i]))
@expression(SCOPF,cost_FL[t=1:nTP,i=1:raw["nFL"]], raw["Load"][raw["FL"]."Bus"[i],"ShCost"]*(P_inc[t,i]+P_dec[t,i])*raw["sbase"])
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
    @constraint(SCOPF,UEmaxd[i=1:raw["nFL"],t=1:nTP],  Init_Data["DA_SCOPF"]["UE"][nTL,i+1]-UE[nTP,i]<= (nTL*nTS-t0)*Inc_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]]*TD)
else
    @expression(SCOPF,FPG[i=1:raw["nGen"]], λG[i]*((PG[nTP,i])^2))
    @expression(SCOPF,FPS[i=1:raw["nESS"]], λS[i]*(soc[nTP,i]-initraw["soc"][nTP,i])^2)
    @expression(SCOPF,FPF[i=1:raw["nFL"]], λF[i]*(UE[nTP,i]-initraw["UE"][nTP,i])^2)
    tempvar=sum(FPG;init=0)+sum(FPS;init=0)+sum(FPF;init=0)
    @expression(SCOPF,cost_Pen, tempvar)
end
# @expression(SCOPF,cost_Pre, TD* sum(cost_Flex[t]+cost_Curt[t]+cost_Shed[t]  for t=1:nTP))
@expression(SCOPF,cost_Pre_T[t=1:nTP], TD*(cost_Flex[t]+cost_Curt[t]+cost_Curtc[t]+cost_Shed[t]))
#Corrective
if OPFon==0
    @expression(SCOPF,cost_Gen_C[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]], raw["Gen"]."cost_c"[i]+ (raw["Gen"]."cost_b"[i]*raw["sbase"]*(PG_c_Line[t,i,cL]))+(raw["Gen"]."cost_a"[i]*(raw["sbase"]*(PG_c_Line[t,i,cL]))^2)-cost_Gen[t,i])
    @expression(SCOPF,cost_ESS_C[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]], raw["ESS"]."cost_c"[i]+ raw["ESS"]."cost_b"[i]*raw["sbase"]*(P_dis_c_Line[t,i,cL]+P_ch_c_Line[t,i,cL])+raw["ESS"]."cost_a"[i]*raw["sbase"]^2*(P_dis_c_Line[t,i,cL]+P_ch_c_Line[t,i,cL])^2-cost_ESS[t,i])
    @expression(SCOPF,cost_FL_C[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]], raw["Load"][raw["FL"]."Bus"[i],"ShCost"]*(P_inc_c_Line[t,i,cL]+P_dec_c_Line[t,i,cL])*raw["sbase"]-cost_FL[t,i])
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
@objective(SCOPF, Min,  sum(cost_T[t] for t=1:nTP)+cost_Pen)
return SCOPF
end


function Normal_PB(SCOPF,nTP,raw,initraw,fixed_SOC,Proffin,Winit,initial)
    #----------------------Initializing Model
    @variables(
        SCOPF,
        begin
            e[t = 1:nTP, i = 1:raw["nBus"]]
            f[t = 1:nTP, i = 1:raw["nBus"]]
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
    for  t=1:nTP, i=1:raw["nRES"]
        set_start_value(curt[t,i], 0)                         # Setting start values
        set_lower_bound(curt[t,i], -1* initraw["Pcurt"][t,i])
        set_upper_bound(curt[t,i], Proffin["Wind"][t,i]-initraw["Pcurt"][t,i])
    end
    for  t=1:nTP, i=1:raw["nLoad"]
        set_start_value(shed[t,i], 0)
        set_lower_bound(shed[t,i], -initraw["Pshed"][t,i])
        set_upper_bound(shed[t,i], Proffin["PL"][t,i]-initraw["Pshed"][t,i])
    end
    for  t=1:nTP, i=1:raw["nBus"]
        set_start_value(e[t,i], 1.00)
        set_start_value(f[t,i], 0.00)
    end
    for t=1:nTP, i=1:raw["nGen"]
        set_lower_bound(PG[t,i], raw["Gen"]."P_min (pu)"[i]-initraw["PG"][t,i])
        set_upper_bound(PG[t,i], raw["Gen"]."P_max (pu)"[i]-initraw["PG"][t,i])
        set_start_value(PG[t,i], 0)
        set_lower_bound(QG[t,i], raw["Gen"]."Q_min (pu)"[i]-initraw["QG"][t,i])
        set_upper_bound(QG[t,i], raw["Gen"]."Q_max (pu)"[i]-initraw["QG"][t,i])
        set_start_value(QG[t,i], 0)
        # set_lower_bound(curtc[t,i], -1* initraw["Pcurtc"][t,i])
        fix(curtc[t,i],0)
    end
    @constraint(SCOPF,NEcurtclim[t=1:nTP,i=1:raw["nGen"]], curtc[t,i]<=PG[t,i]+initraw["PG"][t,i]-initraw["Pcurtc"][t,i])
     for  t=1:nTP, i=1:raw["nESS"]
        set_start_value(P_ch[t,i], 0.5*raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"]);set_lower_bound(P_ch[t,i], -initraw["Pch"][t,i]) ;set_upper_bound(P_ch[t,i], raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"]-initraw["Pch"][t,i])
        set_start_value(P_dis[t,i], 0.5*raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"]);set_lower_bound(P_dis[t,i], -initraw["Pdis"][t,i]) ;set_upper_bound(P_dis[t,i], raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"]-initraw["Pdis"][t,i])
        set_start_value(soc[t,i], 0.5*raw["ESS"]."E_rating (MWh)"[i]/raw["sbase"]);set_lower_bound(soc[t,i], raw["ESS"]."E_rating_min (MWh)"[i]/raw["sbase"]);set_upper_bound(soc[t,i], raw["ESS"]."E_rating (MWh)"[i]/raw["sbase"]);
    end
    for t=1:nTP, i= 1:raw["nFL"]
        # set_start_value(P_inc[t,i], 0.5*load_inc_prct*prof_ploads[i,t])                         # Setting start values
        # set_start_value(P_dec[t,i], 0.5*load_dec_prct*prof_ploads[i,t])
        set_start_value(P_inc[t,i], 0)                         # Setting start values
        set_start_value(P_dec[t,i], 0)
        set_lower_bound(P_inc[t,i], -initraw["Pinc"][t,i])
        set_lower_bound(P_dec[t,i], -initraw["Pdec"][t,i])
        set_upper_bound(P_inc[t,i], Inc_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]]-initraw["Pinc"][t,i])
        set_upper_bound(P_dec[t,i], Dec_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]]-initraw["Pdec"][t,i])
    end

    #-----------------Line flow Modeling  #OK
    Npij_L, Nqij_L, Npji_L, Nqji_L, NIflow2=Flow_Expr(nTP,raw["Line"],e,f,'L')

    @expression(SCOPF,Line_flow[t=1:nTP, n=1:raw["nLine"]],NIflow2[t,n])

    #---------------------- Transformer Modelling ---------------------
    Npij_T, Nqij_T, Npji_T, Nqji_T=Flow_Expr(nTP,raw["Trans"],e,f,'T')

    #-------------------------------Power balance constraint #Oc
    #calculating injected power from the bus to grid based on input and output power considering line outage and networc topology

    Np_inj=Power_Inj(nTP,raw,initraw,P_ch,P_dis, P_inc, P_dec, shed,curt,Npij_L,Npji_L,Npij_T,Npji_T;xPQ='P');
    Nq_inj=Power_Inj(nTP,raw,initraw,P_ch,P_dis, P_inc, P_dec, shed,curt,Nqij_L,Nqji_L,Nqij_T,Nqji_T;xPQ='Q');

    @constraint(SCOPF, NEactive_power_balance[t=1:nTP,i=1:raw["nBus"]], sum(PG[t,m]+initraw["PG"][t,m]-curtc[t,m] for m in 1:raw["nGen"] if raw["Gen"]."Bus"[m]==i)
                                                                      +sum(Proffin["Wind"][t,m] for m in 1:raw["nRES"] if (raw["RES"]."Bus"[m]==i);init=0)
                                                                      ==Proffin["PD"][t,i]+Np_inj[t,i])
                                                                      ## For code clear purpose, we have included flexibility assets in p_inj
    @constraint(SCOPF, NEreactive_power_balance[t=1:nTP,i=1:raw["nBus"]], sum(QG[t,m]+initraw["QG"][t,m] for m in 1:raw["nGen"] if raw["Gen"]."Bus"[m]==i)==Proffin["QD"][t,i]+Nq_inj[t,i])
    #-------------------------- Network constraints
    #@constraint(SCOPF,NPCC[i=1:raw["nBus"];raw["Bus"]."Bus_Type"[i]==3],f[i]==0)
    for t in 1:nTP, i in 1:raw["nBus"] if raw["Bus"]."Bus_Type"[i]==3 fix(f[t,i],0) end end#fixing the slack bus voltage angle to zero

    @constraint(SCOPF,NNvconstraint[t=1:nTP, i=1:raw["nBus"]], raw["Bus"]."V_min (pu)"[i]^2<=e[t,i]^2+f[t,i]^2<=raw["Bus"]."V_max (pu)"[i]^2)

    #---!! we should ignore this constraint for the lines which are already out of service
    @constraint(SCOPF,NNlflow[t=1:nTP, n=1:raw["nLine"]; raw["Line"]."RATE_A"[n]!=0],NIflow2[t,n] <= raw["Line"]."RATE_A"[n]^2)

    #@NLconstraint(SCOPF,NTflow[n=1:raw["nTrans"]],Npji_T[n]^2 + Nqji_T[n]^2 <= raw["Trans"]."Snom"[n]^2)

    #-----------Battery constraints
    @constraint(SCOPF,NNbiess[t=1:nTP, i=1:raw["nESS"]], P_ch[t,i]/(raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"])
    + P_dis[t,i]/(raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"])<=1)

    #-----------Coupling constraints
    if Winit==1
        # @constraint(SCOPF,NNramp1[t=1:1,i=1:raw["nGen"]],(PG[t,i])^2<=(raw["Gen"]."Ramp"[i]*TD)^2)
        @constraint(SCOPF,NNramp1[t=1:1,i=1:raw["nGen"]],-raw["Gen"]."Ramp"[i]*TD<=PG[t,i]+(initraw["PG"][t,i]-initial["PG"][1,i])<=raw["Gen"]."Ramp"[i]*TD)

    end
        #ESS
    @constraint(SCOPF, NEsoc0[t=1:1, i=1:raw["nESS"]], soc[t,i]==initial["soc"][1,i]
                                            + TD* ((P_ch[t,i]+initraw["Pch"][t,i])*raw["ESS"]."Charge_efficiency"[i] - (P_dis[t,i]+initraw["Pdis"][t,i])/(raw["ESS"]."Discharge_efficiency"[i]) ))
    @constraint(SCOPF, NEsocn[t=2:nTP, i=1:raw["nESS"]], soc[t,i]==soc[t-1,i]
                                            + TD* ((P_ch[t,i]+initraw["Pch"][t,i])*raw["ESS"]."Charge_efficiency"[i] - (P_dis[t,i]+initraw["Pdis"][t,i])/(raw["ESS"]."Discharge_efficiency"[i]) ))
        #Gens
    @constraint(SCOPF,NNramp[t=2:nTP,i=1:raw["nGen"]],-raw["Gen"]."Ramp"[i]*TD<=PG[t,i]-PG[t-1,i]+(initraw["PG"][t,i]-initraw["PG"][t-1,i])<=raw["Gen"]."Ramp"[i]*TD)
#     ##fixed Energy all_constraints
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
        @constraint(SCOPF,NNfl[t=1:nTP,i=1:raw["nFL"]],P_inc[t,i]/(Inc_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]])+P_dec[t,i]/(Dec_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]]) <=1)
    end
     #-----------Coupling constraints
    N_Var= PG,QG,P_ch,P_dis,P_inc,P_dec,shunt,shed,curt,curtc,soc,UE,e,f
    return SCOPF,N_Var
end


function Contingency_L_PB( SCOPF,nTP,raw,initraw,fixed_SOC,N_Var,Proffin,CCT,initial,LSpermit)
    #----------------------PARAMETERS----------------------
    v_relax_factor_Line_min=0;
    v_relax_factor_Line_max = 0;
    #----------------------Initializing Model----------------------
    PG,QG,P_ch,P_dis,P_inc,P_dec,NNNN,shed,NNNN,curtc,soc,UE,e,f=N_Var
    @variables(
        SCOPF,
        begin
            e_c_Line[t=1:nTP, i = 1:raw["nBus"],cL = 1:raw["nConL"]]
            f_c_Line[t=1:nTP, i = 1:raw["nBus"],cL = 1:raw["nConL"]]
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
            # sl_F[t = 1:nTP, i = 1:raw["nLine"],cL = 1:raw["nConL"]]
            # sl_V[t = 1:nTP, i = 1:raw["nBus"],cL = 1:raw["nConL"]]
        end
    )
    # for  t=1:nTP, i=1:raw["nRES"],cL=1:raw["nConL"]     fix(curt_c_Line[t,i,cL], 0)  end
    # for  t=1:nTP, i=1:raw["nLoad"],cL=1:raw["nConL"]    fix(shed_c_Line[t,i,cL], 0) end
    if LSpermit==1
        @constraint(SCOPF,CNcurtlim[t=1:nTP,i=1:raw["nRES"],cL=1:raw["nConL"]],-initraw["Pcurt"][t,i] <=curt_c_Line[t,i,cL]<=Proffin["Wind"][t,i]-initraw["Pcurt"][t,i])
        @constraint(SCOPF,CNshedlim[t=1:nTP,i=1:raw["nLoad"],cL=1:raw["nConL"]],-initraw["Pshed"][t,i] <=shed_c_Line[t,i,cL]<=Proffin["PL"][t,i]-initraw["Pshed"][t,i])
        @constraint(SCOPF,CNcurtclim1[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],-initraw["Pcurtc"][t,i]<=curtc_c_Line[t,i,cL])
        @constraint(SCOPF,CNcurtclim2[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],curtc_c_Line[t,i,cL]<=PG_c_Line[t,i,cL]+initraw["PG"][t,i]-initraw["Pcurtc"][t,i])

    else
        @constraint(SCOPF,CNcurtlim[t=1:nTP,i=1:raw["nRES"],cL=1:raw["nConL"]],-initraw["Pcurt"][t,i] <=curt_c_Line[t,i,cL]<=Proffin["Wind"][t,i]-initraw["Pcurt"][t,i])
        # @constraint(SCOPF,CNcurtlim[t=1:nTP,i=1:raw["nRES"],cL=1:raw["nConL"]],curt_c_Line[t,i,cL]<=0)
        @constraint(SCOPF,CNshedlim[t=1:nTP,i=1:raw["nLoad"],cL=1:raw["nConL"]],shed_c_Line[t,i,cL]==shed[t,i])
        @constraint(SCOPF,CNcurtclim[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],curtc_c_Line[t,i,cL]==curtc[t,i])

    end

    for  t=1:nTP, i=1:raw["nRES"],cL=1:raw["nConL"]     set_start_value(curt_c_Line[t,i,cL], 0)  end                       # Setting start values
    for  t=1:nTP, i=1:raw["nLoad"],cL=1:raw["nConL"]    set_start_value(shed_c_Line[t,i,cL], 0) end

    for t=1:nTP, i=1:raw["nBus"],cL=1:raw["nConL"]
        set_start_value(e_c_Line[t,i,cL], 1.00)                         # Setting start values
        set_start_value(f_c_Line[t,i,cL], 0.00)
        # fix(sl_V[t,i,cL], 0.0)
        # set_lower_bound(sl_V[t,i,cL], 0.0)
        # set_start_value(sl_V[t,i,cL], 0.0)
    end
    # for t=1:nTP, i=1:raw["nLine"],cL=1:raw["nConL"]
    #     # fix(sl_F[t,i,cL], 0.0)
    #     # set_lower_bound(sl_F[t,i,cL], 0.0)
    #     # set_start_value(sl_F[t,i,cL], 0.0)
    # end

    @constraint(SCOPF,CNplim[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],raw["Gen"]."P_min (pu)"[i]<=PG_c_Line[t,i,cL]+initraw["PG"][t,i]<=raw["Gen"]."P_max (pu)"[i])
    @constraint(SCOPF,CNqlim[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],raw["Gen"]."Q_min (pu)"[i]<=QG_c_Line[t,i,cL]+initraw["QG"][t,i]<=raw["Gen"]."Q_max (pu)"[i])

    for t=1:nTP, i=1:raw["nGen"],cL=1:raw["nConL"]
        set_start_value(PG_c_Line[t,i,cL], 0)
        set_start_value(QG_c_Line[t,i,cL], 0)
        fix(curtc_c_Line[t,i,cL], 0)
    end

    @constraint(SCOPF,CNbchlim[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],0<=P_ch_c_Line[t,i,cL]<=raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"])
    @constraint(SCOPF,CNbdislim[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],0<=P_dis_c_Line[t,i,cL]<=raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"])
    for  t=1:nTP, i=1:raw["nESS"], cL=1: raw["nConL"]
        set_start_value(P_ch_c_Line[t,i,cL], 0.5*raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"])
        set_start_value(P_dis_c_Line[t,i,cL], 0.5*raw["ESS"]."Discharge_rating (MW)"[i]/raw["sbase"])
        set_start_value(soc_c_Line[t,i,cL], 0.5*raw["ESS"]."E_rating (MWh)"[i]/raw["sbase"]);set_lower_bound(soc_c_Line[t,i,cL], raw["ESS"]."E_rating_min (MWh)"[i]/raw["sbase"]);set_upper_bound(soc_c_Line[t,i,cL], raw["ESS"]."E_rating (MWh)"[i]/raw["sbase"]);
    end

    @constraint(SCOPF,CNfilim[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],0<=P_inc_c_Line[t,i,cL]+initraw["Pinc"][t,i] <=Inc_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]])
    @constraint(SCOPF,CNfdlim[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],0<=P_dec_c_Line[t,i,cL]+initraw["Pdec"][t,i] <=Dec_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]])
    for t=1:nTP, i= 1:raw["nFL"],cL=1:raw["nConL"]
        set_start_value(P_inc_c_Line[t,i,cL], 0)                         # Setting start values
        set_start_value(P_dec_c_Line[t,i,cL], 0)
    end


    #-----------------Line flow Modeling  #OK
    CLpij_L, CLqij_L, CLpji_L, CLqji_L, CLIflow2=Flow_Expr(nTP,raw["Line"],e_c_Line,f_c_Line,'L',raw["ConL"])
    @expression(SCOPF,Line_flow_c_Line[t=1:nTP, n=1:raw["nLine"],cL=1:raw["nConL"]],CLIflow2[t,n,cL])

    #---------------------- Transformer Modelling ---------------------
    CLpij_T, CLqij_T, CLpji_T, CLqji_T=Flow_Expr(nTP,raw["Trans"],e_c_Line,f_c_Line,'T',raw["ConL"])


    #-------------------------------Power balance constraint #Oc
    # calculating injected power from the bus to grid based on iCLNput and output power considering line outage and networc topology

    CLp_inj=Power_Inj(nTP,raw,initraw,P_ch_c_Line,P_dis_c_Line,P_inc_c_Line,P_dec_c_Line,shed_c_Line,curt_c_Line,CLpij_L,CLpji_L,CLpij_T,CLpji_T;xPQ='P',C="CL")
    CLq_inj=Power_Inj(nTP,raw,initraw,P_ch_c_Line,P_dis_c_Line,P_inc_c_Line,P_dec_c_Line,shed_c_Line,curt_c_Line,CLqij_L,CLqji_L,CLqij_T,CLqji_T;xPQ='Q',C="CL")



    @constraint(SCOPF, CEActive_power_balance[t=1:nTP, i=1:raw["nBus"],cL=1:raw["nConL"]], sum(PG_c_Line[t,m,cL]+initraw["PG"][t,m]-curtc_c_Line[t,m,cL] for m in 1:raw["nGen"] if raw["Gen"]."Bus"[m]==i)
                                                                                          +sum(Proffin["Wind"][t,m] for m in 1:raw["nRES"] if (raw["RES"]."Bus"[m]==i);init=0)
                                                                                          -Proffin["PD"][t,i]==CLp_inj[t,i,cL])
    @constraint(SCOPF, CEReactive_power_balance[t=1:nTP, i=1:raw["nBus"],cL=1:raw["nConL"]], sum(QG_c_Line[t,m,cL]+initraw["QG"][t,m] for m in 1:raw["nGen"] if raw["Gen"]."Bus"[m]==i)-Proffin["QD"][t,i]==CLq_inj[t,i,cL])
    #-------------------------- Network constraints
    # @constraint(SCOPF,CLpCC[i=1:raw["nBus"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[i]==3],f_c_Line[i,cL]==0)
    for t in 1:nTP, i in 1:raw["nBus"],cL=1:raw["nConL"]
        if raw["Bus"]."Bus_Type"[i]==3
            fix(f_c_Line[t,i,cL],0)
            # set_lower_bound(e_c_Line[t,i,cL],raw["Bus"]."V_min (pu)"[i])
            # set_upper_bound(e_c_Line[t,i,cL],raw["Bus"]."V_max (pu)"[i])
        end
    end#fixing the slack bus voltage angle to zero

    # @constraint(SCOPF,CLVConstraintMIN[t=1:nTP, i=1:raw["nBus"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[i]!=3 ], e_c_Line[t,i,cL]^2+f_c_Line[t,i,cL]^2>=((1-sl_V[t,i,cL])*(raw["Bus"]."V_min (pu)"[i]-v_relax_factor_Line_min))^2)
    # @constraint(SCOPF,CLVConstraintMAX[t=1:nTP, i=1:raw["nBus"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[i]!=3 ], e_c_Line[t,i,cL]^2+f_c_Line[t,i,cL]^2<=((1+sl_V[t,i,cL])*(raw["Bus"]."V_max (pu)"[i]+v_relax_factor_Line_max))^2)
    @constraint(SCOPF,CNvconstraint[t=1:nTP, i=1:raw["nBus"],cL=1:raw["nConL"] ], (raw["Bus"]."V_min (pu)"[i]-v_relax_factor_Line_min)^2<=e_c_Line[t,i,cL]^2+f_c_Line[t,i,cL]^2<=(raw["Bus"]."V_max (pu)"[i]+v_relax_factor_Line_max)^2)
    # @constraint(SCOPF,CLflow[t=1:nTP,n=1:raw["nLine"],cL=1:raw["nConL"]],CLIflow2[t,n,cL] <= ((1+sl_F[t,n,cL])*raw["Line"]."RATE_A"[n])^2)
    @constraint(SCOPF,CNlflow[t=1:nTP,n=1:raw["nLine"],cL=1:raw["nConL"]],CLIflow2[t,n,cL] <= raw["Line"]."RATE_A"[n]^2)

    @assert(Mode in ["Preventive","Corrective"])
    if Mode=="Preventive"
        @constraint(SCOPF,CEfixp[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[raw["Gen"]."Bus"[i]] in [1,2]],PG_c_Line[t,i,cL]==PG[t,i])
        # @constraint(SCOPF,CEfixQ[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[raw["Gen"]."Bus"[i]]!=3],QG_c_Line[t,i,cL]==QG[t,i])
        @constraint(SCOPF,CEfixESS1[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],P_dis_c_Line[t,i,cL]==P_dis[t,i])
        @constraint(SCOPF,CEfixESS2[t=1:nTP,i=1:raw["nESS"],cL=1:raw["nConL"]],P_ch_c_Line[t,i,cL]==P_ch[t,i])
        @constraint(SCOPF,CEfixFL1[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],P_inc_c_Line[t,i,cL]==P_inc[t,i])
        @constraint(SCOPF,CEfixFL2[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],P_dec_c_Line[t,i,cL]==P_dec[t,i])
        @constraint(SCOPF,CEfixedvoltage[t=1:nTP,i=1:raw["nBus"],cL=1:raw["nConL"];raw["Bus"]."Bus_Type"[i] in [2,3]],e_c_Line[t,i,cL]^2+f_c_Line[t,i,cL]^2==e[t,i]^2+f[t,i]^2)
    elseif Mode =="Corrective"
        #-----------Coupling constraints
        @constraint(SCOPF,CNramp1[t=1:nTP,i=1:raw["nGen"],cL=1:raw["nConL"]],(PG_c_Line[t,i,cL]-PG[t,i])^2<=(raw["Gen"]."Ramp"[i]*CCT)^2 )
        # @constraint(SCOPF,CNramp[t=2:nTP,i=1:raw["nGen"],cL=1: raw["nConL"]],-raw["Gen"]."Ramp"[i]*TD<=PG_c_Line[t,i,cL]-PG_c_Line[t-1,i,cL]+initraw["PG"][t,i]-initraw["PG"][t-1,i]<=raw["Gen"]."Ramp"[i]*TD)
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
            @constraint(SCOPF,CNfl[t=1:nTP,i=1:raw["nFL"],cL=1:raw["nConL"]],(P_inc_c_Line[t,i,cL])/(Inc_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]])+(P_dec_c_Line[t,i,cL])/(Dec_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]]) <=1)
        end
    end

    C_Var= PG_c_Line,QG_c_Line,P_ch_c_Line,P_dis_c_Line,P_inc_c_Line,P_dec_c_Line,shed_c_Line,curt_c_Line,curtc_c_Line
    return SCOPF,C_Var #,sl_F,sl_V
end