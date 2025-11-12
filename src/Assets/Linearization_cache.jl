function Lin_Yang(init)
    V_0,θ_0,PG_0=init
    #-------------------------Optimizer
    LACOPF = Model(Ipopt.Optimizer)
    @variables(
        LACOPF,
        begin
            V2[k = 1:size(SysData["Bus"],1)]
            θ[k = 1:size(SysData["Bus"],1)]
            PG[k = 1:size(SysData["Gen"],1)]
            QG[k = 1:size(SysData["Gen"],1)]
        end
    )
    @expression(LACOPF, pflow_ij_sh[k=1:size(SysData["Line"],1)],(SysData["Line"]."g_sh (pu)"[k]/2)*(V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, pflow_ij_sr1[k=1:size(SysData["Line"],1)],(SysData["Line"]."g (pu)"[k])*(V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, pflow_ij_sr2[k=1:size(SysData["Line"],1)],-(SysData["Line"]."g (pu)"[k]/2)*(V2[SysData["Line"]."From"[k]]+V2[SysData["Line"]."To"[k]])-(SysData["Line"]."b (pu)"[k])*(θ[SysData["Line"]."From"[k]]-θ[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, ploss_ij_V[k=1:size(SysData["Line"],1)],(SysData["Line"]."g (pu)"[k]/2)*(2*((V_0[SysData["Line"]."From"[k]]-V_0[SysData["Line"]."To"[k]])/(V_0[SysData["Line"]."From"[k]]+V_0[SysData["Line"]."To"[k]]))*(V2[SysData["Line"]."From"[k]]-V2[SysData["Line"]."To"[k]])-(V_0[SysData["Line"]."From"[k]]-V_0[SysData["Line"]."To"[k]])^2))
    @expression(LACOPF, ploss_ij_θ[k=1:size(SysData["Line"],1)], (SysData["Line"]."g (pu)"[k])*((θ_0[SysData["Line"]."From"[k]]-θ_0[SysData["Line"]."To"[k]])*(θ[SysData["Line"]."From"[k]]-θ[SysData["Line"]."To"[k]])-((θ_0[SysData["Line"]."From"[k]]-θ_0[SysData["Line"]."To"[k]])^2)/2))
    @expression(LACOPF, ploss_ij, ploss_ij_V + ploss_ij_θ)
    @expression(LACOPF, pflow_ij_sr3,pflow_ij_sr2+ ploss_ij) #Line flow expression
    @expression(LACOPF, pij, pflow_ij_sh + pflow_ij_sr1 + pflow_ij_sr2+ ploss_ij)
    @constraint(LACOPF, Plosscon[k=1:size(SysData["Line"],1)], ploss_ij_V[k]>=0)
    #---------------------------------------------------
    @expression(LACOPF, pflow_ji_sh[k=1:size(SysData["Line"],1)],(SysData["Line"]."g_sh (pu)"[k]/2)*(V2[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, pflow_ji_sr1[k=1:size(SysData["Line"],1)],(SysData["Line"]."g (pu)"[k]/2)*(V2[SysData["Line"]."To"[k]]-V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, pflow_ji_sr2[k=1:size(SysData["Line"],1)],-(SysData["Line"]."b (pu)"[k])*(θ[SysData["Line"]."To"[k]]-θ[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, ploss_ji_V[k=1:size(SysData["Line"],1)],(SysData["Line"]."g (pu)"[k]/2)*(2*((V_0[SysData["Line"]."To"[k]]-V_0[SysData["Line"]."From"[k]])/(V_0[SysData["Line"]."To"[k]]+V_0[SysData["Line"]."From"[k]]))*(V2[SysData["Line"]."To"[k]]-V2[SysData["Line"]."From"[k]])-(V_0[SysData["Line"]."To"[k]]-V_0[SysData["Line"]."From"[k]])^2))
    @expression(LACOPF, ploss_ji_θ[k=1:size(SysData["Line"],1)], (SysData["Line"]."g (pu)"[k])*((θ_0[SysData["Line"]."To"[k]]-θ_0[SysData["Line"]."From"[k]])*(θ[SysData["Line"]."To"[k]]-θ[SysData["Line"]."From"[k]])-((θ_0[SysData["Line"]."To"[k]]-θ_0[SysData["Line"]."From"[k]])^2)/2))
    @expression(LACOPF, ploss_ji, ploss_ji_V + ploss_ji_θ)
    @expression(LACOPF, pji, pflow_ji_sh + pflow_ji_sr1 + pflow_ji_sr2+ ploss_ji)
    @constraint(LACOPF, Plosscon1[k=1:size(SysData["Line"],1)], ploss_ji_V[k]>=0)
    #---------------------------------------------------
    @expression(LACOPF, qflow_ij_sh[k=1:size(SysData["Line"],1)],(-SysData["Line"]."b_sh (pu)"[k]/2)*(V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, qflow_ij_sr1[k=1:size(SysData["Line"],1)],(-SysData["Line"]."b (pu)"[k]/2)*(V2[SysData["Line"]."From"[k]]-V2[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, qflow_ij_sr2[k=1:size(SysData["Line"],1)],(-SysData["Line"]."g (pu)"[k])*(θ[SysData["Line"]."From"[k]]-θ[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, qloss_ij_V[k=1:size(SysData["Line"],1)], (-SysData["Line"]."b (pu)"[k]/2)*(2*((V_0[SysData["Line"]."From"[k]]-V_0[SysData["Line"]."To"[k]])/(V_0[SysData["Line"]."From"[k]]+V_0[SysData["Line"]."To"[k]]))*(V2[SysData["Line"]."From"[k]]-V2[SysData["Line"]."To"[k]])-(V_0[SysData["Line"]."From"[k]]-V_0[SysData["Line"]."To"[k]])^2))
    @expression(LACOPF, qloss_ij_θ[k=1:size(SysData["Line"],1)], (-SysData["Line"]."b (pu)"[k])*((θ_0[SysData["Line"]."From"[k]]-θ_0[SysData["Line"]."To"[k]])*(θ[SysData["Line"]."From"[k]]-θ[SysData["Line"]."To"[k]])-((θ_0[SysData["Line"]."From"[k]]-θ_0[SysData["Line"]."To"[k]])^2)/2))
    @expression(LACOPF, qloss_ij, qloss_ij_V + qloss_ij_θ)
    @expression(LACOPF, qij, qflow_ij_sh + qflow_ij_sr1 + qflow_ij_sr2+ qloss_ij)
    #---------------------------------------------------
    @expression(LACOPF, qflow_ji_sh[k=1:size(SysData["Line"],1)],(-SysData["Line"]."b_sh (pu)"[k]/2)*(V2[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, qflow_ji_sr1[k=1:size(SysData["Line"],1)],(-SysData["Line"]."b (pu)"[k]/2)*(V2[SysData["Line"]."To"[k]]-V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, qflow_ji_sr2[k=1:size(SysData["Line"],1)],(-SysData["Line"]."g (pu)"[k])*(θ[SysData["Line"]."To"[k]]-θ[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, qloss_ji_V[k=1:size(SysData["Line"],1)], (-SysData["Line"]."b (pu)"[k]/2)*(2*((V_0[SysData["Line"]."To"[k]]-V_0[SysData["Line"]."From"[k]])/(V_0[SysData["Line"]."To"[k]]+V_0[SysData["Line"]."From"[k]]))*(V2[SysData["Line"]."To"[k]]-V2[SysData["Line"]."From"[k]])-(V_0[SysData["Line"]."To"[k]]-V_0[SysData["Line"]."From"[k]])^2))
    @expression(LACOPF, qloss_ji_θ[k=1:size(SysData["Line"],1)], (-SysData["Line"]."b (pu)"[k])*((θ_0[SysData["Line"]."To"[k]]-θ_0[SysData["Line"]."From"[k]])*(θ[SysData["Line"]."To"[k]]-θ[SysData["Line"]."From"[k]])-((θ_0[SysData["Line"]."To"[k]]-θ_0[SysData["Line"]."From"[k]])^2)/2))
    @expression(LACOPF, qloss_ji, qloss_ji_V + qloss_ji_θ)
    @expression(LACOPF, qji, qflow_ji_sh + qflow_ji_sr1 + qflow_ji_sr2+ qloss_ji)
    #-------------------------------Power balance constraint #OK
    @expression(LACOPF, p_inj[i=1:size(SysData["Bus"],1)], sum(pij[k] for k in 1:size(SysData["Line"],1) if SysData["Line"]."From"[k]==i)+sum(pji[k] for k in 1:size(SysData["Line"],1) if SysData["Line"]."To"[k]==i))
    @expression(LACOPF, q_inj[i=1:size(SysData["Bus"],1)], sum(qij[k] for k in 1:size(SysData["Line"],1) if SysData["Line"]."From"[k]==i)+sum(qji[k] for k in 1:size(SysData["Line"],1) if SysData["Line"]."To"[k]==i))
    @expression(LACOPF, PD[i=1:size(SysData["Bus"],1)], sum(SysData["Load"]."Pd (pu)"[j] for j in 1:size(SysData["Load"],1) if SysData["Load"]."Bus_N"[j]==i))
    @expression(LACOPF, QD[i=1:size(SysData["Bus"],1)], sum(SysData["Load"]."Qd (pu)"[j] for j in 1:size(SysData["Load"],1) if SysData["Load"]."Bus_N"[j]==i))
    @constraint(LACOPF, Active_power_balance[i=1:size(SysData["Bus"],1)], sum(PG[j] for j in 1:size(SysData["Gen"],1) if SysData["Gen"]."Bus"[j]==i)-PD[i]==p_inj[i])
    @constraint(LACOPF, Reactive_power_balance[i=1:size(SysData["Bus"],1)], sum(QG[j] for j in 1:size(SysData["Gen"],1) if SysData["Gen"]."Bus"[j]==i)-QD[i]==q_inj[i])
    #-------------------------- Network constraints
    @constraint(LACOPF,PCC[i=1:size(SysData["Bus"],1);SysData["Bus"]."Bus_Type"[i]==3],θ[i]==0)
    # @constraint(LACOPF,PCC1[i=1:size(SysData["Bus"],1);SysData["Bus"]."Bus_Type"[i]==3],V2[i]<=1.1)
    # @constraint(LACOPF,PCC_V[i=1:size(SysData["Bus"],1);SysData["Bus"]."Bus_Type"[i]==3], SysData["Bus"]."V_min (pu)"[i]<=e[i]<=SysData["Bus"]."V_max (pu)"[i])
    @constraint(LACOPF,VConstraint[i=1:size(SysData["Bus"],1)], SysData["Bus"]."V_min (pu)"[i]^2<=V2[i]<=SysData["Bus"]."V_max (pu)"[i] ^2)
    @constraint(LACOPF,Gen_Limit_P[i=1:size(SysData["Gen"],1)],SysData["Gen"]."P_min (pu)"[i] <= PG[i] <= SysData["Gen"]."P_max (pu)"[i])
    @constraint(LACOPF,Gen_Limit_Q[i=1:size(SysData["Gen"],1)],SysData["Gen"]."Q_min (pu)"[i] <= QG[i] <= SysData["Gen"]."Q_max (pu)"[i])
    @constraint(LACOPF, flow1[k=1:size(SysData["Line"],1)], pij[k]+tan(π/6)*qij[k]<= SysData["Line"]."RATE_A"[k])
    @constraint(LACOPF, flow2[k=1:size(SysData["Line"],1)], pij[k]-tan(π/6)*qij[k]<= SysData["Line"]."RATE_A"[k])
    @constraint(LACOPF, flow3[k=1:size(SysData["Line"],1)], pij[k]+tan(π/6)*qij[k]>= -SysData["Line"]."RATE_A"[k])
    @constraint(LACOPF, flow4[k=1:size(SysData["Line"],1)], pij[k]-tan(π/6)*qij[k]>= -SysData["Line"]."RATE_A"[k])
    #----------------Objective function #OK

    @expression(LACOPF,cost_gen[i=1:size(SysData["Gen"],1)], SysData["Gen"]."cost_c"[i]+ SysData["Gen"]."cost_b"[i]*SysData["Base_MVA"]."base_MVA"[1]*PG[i]+SysData["Gen"]."cost_a"[i]*SysData["Base_MVA"]."base_MVA"[1]^2*PG[i]*PG[i])
    @objective(LACOPF, Min, sum(cost_gen[i] for i=1:size(SysData["Gen"],1))) # Setting objective function

    #----------------Solve and Results
    optimize!(LACOPF)
    println("Objective value:  ", JuMP.objective_value(LACOPF))
    println("Solver Time:   ", JuMP.solve_time(LACOPF))

    θ1=value.(LACOPF[:θ])
    θ1[5]=0
    V1=sqrt.(value.(LACOPF[:V2]))
    PG1=value.(LACOPF[:PG])

    return V1,θ1,PG1,LACOPF
end

function Lin_SOCP(init)
    V_0,θ_0,PG_0=init
    V2_0=V_0.^2
    A_0=zeros(size(SysData["Line"],1))
    B_0=zeros(size(SysData["Line"],1))
    for k in 1:size(SysData["Line"],1)
        A_0[k]=V_0[SysData["Line"]."From"[k]]*V_0[SysData["Line"]."To"[k]]*cos(θ_0[SysData["Line"]."From"[k]]-θ_0[SysData["Line"]."To"[k]])
        B_0[k]=V_0[SysData["Line"]."From"[k]]*V_0[SysData["Line"]."To"[k]]*sin(θ_0[SysData["Line"]."From"[k]]-θ_0[SysData["Line"]."To"[k]])
    end
    #-------------------------Optimizer
    LACOPF = Model(Ipopt.Optimizer)
    @variables(
        LACOPF,
        begin
            V2[k = 1:size(SysData["Bus"],1)]
            A[k = 1:size(SysData["Line"],1)]
            B[k = 1:size(SysData["Line"],1)]
            PG[k = 1:size(SysData["Gen"],1)]
            QG[k = 1:size(SysData["Gen"],1)]
        end
    )
    #-----------------Line flow Modeling #OK
    @expression(LACOPF, pflow_ij_sh[k=1:size(SysData["Line"],1)],(SysData["Line"]."g_sh (pu)"[k]/2)*(V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, pflow_ij_sr1[k=1:size(SysData["Line"],1)],(SysData["Line"]."g (pu)"[k])*(V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, pflow_ij_sr2[k=1:size(SysData["Line"],1)],-(SysData["Line"]."g (pu)"[k])*(A[k])-(SysData["Line"]."b (pu)"[k])*(B[k])) #Line flow expression 
    @expression(LACOPF, pij, pflow_ij_sh + pflow_ij_sr1 + pflow_ij_sr2)
    #-------------------
    @expression(LACOPF, pflow_ji_sh[k=1:size(SysData["Line"],1)],(SysData["Line"]."g_sh (pu)"[k]/2)*(V2[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, pflow_ji_sr1[k=1:size(SysData["Line"],1)],(SysData["Line"]."g (pu)"[k])*(V2[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, pflow_ji_sr2[k=1:size(SysData["Line"],1)],-(SysData["Line"]."g (pu)"[k])*(A[k])+(SysData["Line"]."b (pu)"[k])*(B[k])) #Line flow expression 
    @expression(LACOPF, pji, pflow_ji_sh + pflow_ji_sr1 + pflow_ji_sr2)
    #-------------------
    @expression(LACOPF, qflow_ij_sh[k=1:size(SysData["Line"],1)],(-SysData["Line"]."b_sh (pu)"[k]/2)*(V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, qflow_ij_sr1[k=1:size(SysData["Line"],1)],(-SysData["Line"]."b (pu)"[k])*(V2[SysData["Line"]."From"[k]])) #Line flow expression
    @expression(LACOPF, qflow_ij_sr2[k=1:size(SysData["Line"],1)],(SysData["Line"]."b (pu)"[k])*(A[k])-(SysData["Line"]."g (pu)"[k])*(B[k])) #Line flow expression 
    @expression(LACOPF, qij, qflow_ij_sh + qflow_ij_sr1 + qflow_ij_sr2)
    #-------------------
    @expression(LACOPF, qflow_ji_sh[k=1:size(SysData["Line"],1)],(-SysData["Line"]."b_sh (pu)"[k]/2)*(V2[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, qflow_ji_sr1[k=1:size(SysData["Line"],1)],(-SysData["Line"]."b (pu)"[k])*(V2[SysData["Line"]."To"[k]])) #Line flow expression
    @expression(LACOPF, qflow_ji_sr2[k=1:size(SysData["Line"],1)],(SysData["Line"]."b (pu)"[k])*(A[k])+(SysData["Line"]."g (pu)"[k])*(B[k])) #Line flow expression 
    @expression(LACOPF, qji, qflow_ji_sh + qflow_ji_sr1 + qflow_ji_sr2)
    #-------------------------------Power balance constraint #OK
    @expression(LACOPF, p_inj[i=1:size(SysData["Bus"],1)], sum(pij[k] for k in 1:size(SysData["Line"],1) if SysData["Line"]."From"[k]==i)+sum(pji[k] for k in 1:size(SysData["Line"],1) if SysData["Line"]."To"[k]==i))
    @expression(LACOPF, q_inj[i=1:size(SysData["Bus"],1)], sum(qij[k] for k in 1:size(SysData["Line"],1) if SysData["Line"]."From"[k]==i)+sum(qji[k] for k in 1:size(SysData["Line"],1) if SysData["Line"]."To"[k]==i))
    @expression(LACOPF, PD[i=1:size(SysData["Bus"],1)], sum(SysData["Load"]."Pd (pu)"[j] for j in 1:size(SysData["Load"],1) if SysData["Load"]."Bus_N"[j]==i))
    @expression(LACOPF, QD[i=1:size(SysData["Bus"],1)], sum(SysData["Load"]."Qd (pu)"[j] for j in 1:size(SysData["Load"],1) if SysData["Load"]."Bus_N"[j]==i))
    #-------------------
    @constraint(LACOPF, Active_power_balance[i=1:size(SysData["Bus"],1)], sum(PG[j] for j in 1:size(SysData["Gen"],1) if SysData["Gen"]."Bus"[j]==i)-PD[i]==p_inj[i])
    @constraint(LACOPF, Reactive_power_balance[i=1:size(SysData["Bus"],1)], sum(QG[j] for j in 1:size(SysData["Gen"],1) if SysData["Gen"]."Bus"[j]==i)-QD[i]==q_inj[i])
    #-------------------------- Network constraints
    # @constraint(LACOPF, LPCC[k=1:size(SysData["Line"],1);SysData["Bus"]."Bus_Type"[SysData["Line"]."To"[k]]==3], V2[SysData["Line"]."From"[k]]*V2[SysData["Line"]."To"[k]]==((A[k]-V2[SysData["Line"]."From"[k]]-V2[SysData["Line"]."To"[k]])/2)^2) # Voltage magnitude constraint
    @constraint(LACOPF, LPCC[k=1:size(SysData["Line"],1)], V2[SysData["Line"]."From"[k]]*V2_0[SysData["Line"]."To"[k]]+V2[SysData["Line"]."To"[k]]*V2_0[SysData["Line"]."From"[k]]==2*A[k]*A_0[k]+2*B[k]*B_0[k]) # Voltage magnitude constraint
    # @constraint(LACOPF,PCC[i=1:size(SysData["Bus"],1);SysData["Bus"]."Bus_Type"[i]==3],V2[i]>=1.0968)
    # @constraint(LACOPF,PCC1[i=1:size(SysData["Bus"],1);SysData["Bus"]."Bus_Type"[i]==3],V2[i]<=1.1)
    # @constraint(LACOPF,PCC_V[i=1:size(SysData["Bus"],1);SysData["Bus"]."Bus_Type"[i]==3], SysData["Bus"]."V_min (pu)"[i]<=es[i]<=SysData["Bus"]."V_max (pu)"[i])
    @constraint(LACOPF,VConstraint[i=1:size(SysData["Bus"],1)], SysData["Bus"]."V_min (pu)"[i]^2<=V2[i]<=SysData["Bus"]."V_max (pu)"[i] ^2)
    @constraint(LACOPF,Gen_Limit_P[i=1:size(SysData["Gen"],1)],SysData["Gen"]."P_min (pu)"[i] <= PG[i] <= SysData["Gen"]."P_max (pu)"[i])
    @constraint(LACOPF,Gen_Limit_Q[i=1:size(SysData["Gen"],1)],SysData["Gen"]."Q_min (pu)"[i] <= QG[i] <= SysData["Gen"]."Q_max (pu)"[i])
    @expression(LACOPF,Iflow2_ij[k=1:size(SysData["Line"],1)], (SysData["Line"]."g (pu)"[k]^2 + SysData["Line"]."b (pu)"[k]^2)*(V2[SysData["Line"]."From"[k]]+V2[SysData["Line"]."To"[k]]-2*A[k]))# lineflow current limitation
    @constraint(LACOPF,flow[k=1:size(SysData["Line"],1)],Iflow2_ij[k] <= SysData["Line"]."RATE_A"[k]^2)
    #----------------Objective function #OK
    @expression(LACOPF,cost_gen[i=1:size(SysData["Gen"],1)], SysData["Gen"]."cost_c"[i]+ SysData["Gen"]."cost_b"[i]*SysData["Base_MVA"]."base_MVA"[1]*PG[i]+SysData["Gen"]."cost_a"[i]*SysData["Base_MVA"]."base_MVA"[1]^2*PG[i]*PG[i])
    @objective(LACOPF, Min, sum(cost_gen[i] for i=1:size(SysData["Gen"],1))) # Setting objective function
    #----------------Solve and Results
    optimize!(LACOPF)
    println("Objective value:  ", JuMP.objective_value(LACOPF))
    println("Solver Time:   ", JuMP.solve_time(LACOPF))
    # approx_error=zeros(size(SysData["Line"],1))
    # for k=1:size(SysData["Line"],1)
    #     approx_error[k]=value.(V2[SysData["Line"]."From"[k]]*V2[SysData["Line"]."To"[k]]-A[k]^2-B[k]^2)
    # end

    ############## Model Recovery ##############################
    A_res=value.(A)
    B_res=value.(B)
    V2_res=value.(V2)

    Recover = Model(Ipopt.Optimizer)
    @variables(
        Recover,
        begin
            es[k = 1:size(SysData["Bus"],1)]
            fs[k = 1:size(SysData["Bus"],1)]
            ϵ1[k=1:size(SysData["Line"],1)]
            ϵ2[k=1:size(SysData["Line"],1)]
            ϵ3[k=1:size(SysData["Bus"],1)]
        end
    )

    @constraint(Recover, first[k=1:size(SysData["Line"],1)], es[SysData["Line"]."From"[k]]*es[SysData["Line"]."To"[k]]+fs[SysData["Line"]."From"[k]]*fs[SysData["Line"]."To"[k]] == A_res[k]+ϵ1[k])
    @constraint(Recover, second[k=1:size(SysData["Line"],1)], fs[SysData["Line"]."From"[k]]*es[SysData["Line"]."To"[k]]-es[SysData["Line"]."From"[k]]*fs[SysData["Line"]."To"[k]] == B_res[k]+ϵ2[k])
    @constraint(Recover, Third[k=1:size(SysData["Bus"],1)], es[k]^2+fs[k]^2 == V2_res[k]+ϵ3[k])
    @constraint(Recover, PCC[k=1:size(SysData["Bus"],1) ; SysData["Bus"]."Bus_Type"[k]==3],fs[k]==0 )
    @constraint(Recover,limit[k=1:size(SysData["Bus"],1)], es[k]>=0)

    @objective(Recover,Min, sum(ϵ1.^2)+sum(ϵ2.^2)+sum(ϵ3.^2))
    optimize!(Recover)
    println("Objective value:  ", JuMP.objective_value(Recover))
    if JuMP.objective_value(Recover)>=1e-4
        @error("Recovery error is higher than 1e-4")
        @info("Recovery error is $(JuMP.objective_value(Recover))")
    end


    θ1=atan.(value.(Recover[:fs])./value.(Recover[:es])) 
    θ1[5]=0
    V1=sqrt.(value.(LACOPF[:V2])) 
    PG1=value.(LACOPF[:PG])

    return V1,θ1,PG1
end
function ACOPF_temp()
#---------------- Set Parameter
nGens=size(SysData["Gen"],1); nBus=size(SysData["Bus"],1);nload=size(SysData["Load"],1);nline=size(SysData["Line"],1);sbase=SysData["Base_MVA"][1,"base_MVA"]

#-------------------Start Modeling #OK
ACOPF = Model(Ipopt.Optimizer)

@variables(
    ACOPF,
    begin
        e[k = 1:nBus]
        f[k = 1:nBus]
        PG[k = 1:nGens]
        QG[k = 1:nGens]
    end
)

for i = 1:nBus
    set_start_value(e[i], 1.00)                         # Setting start values
    set_start_value(f[i], 0.00)
end

#-----------------Line flow Modeling #OK

@expression(ACOPF, pflow_ij_sh[k=1:nline],(SysData["Line"]."g_sh (pu)"[k]/2)*(e[SysData["Line"]."From"[k]]^2+f[SysData["Line"]."From"[k]]^2)) #Line flow expression 
@expression(ACOPF, pflow_ij_sr1[k=1:nline],(SysData["Line"]."g (pu)"[k])*(e[SysData["Line"]."From"[k]]^2+f[SysData["Line"]."From"[k]]^2)) #Line flow expression 
@expression(ACOPF, pflow_ij_sr2[k=1:nline],-(SysData["Line"]."g (pu)"[k])*(e[SysData["Line"]."From"[k]]*e[SysData["Line"]."To"[k]]+f[SysData["Line"]."From"[k]]*f[SysData["Line"]."To"[k]])-(SysData["Line"]."b (pu)"[k])*(f[SysData["Line"]."From"[k]]*e[SysData["Line"]."To"[k]]-e[SysData["Line"]."From"[k]]*f[SysData["Line"]."To"[k]])) #Line flow expression 
@expression(ACOPF, pij, pflow_ij_sh + pflow_ij_sr1 + pflow_ij_sr2)

@expression(ACOPF, pflow_ji_sh[k=1:nline],(SysData["Line"]."g_sh (pu)"[k]/2)*(e[SysData["Line"]."To"[k]]^2+f[SysData["Line"]."To"[k]]^2)) #Line flow expression 
@expression(ACOPF, pflow_ji_sr1[k=1:nline],(SysData["Line"]."g (pu)"[k])*(e[SysData["Line"]."To"[k]]^2+f[SysData["Line"]."To"[k]]^2)) #Line flow expression 
@expression(ACOPF, pflow_ji_sr2[k=1:nline],-(SysData["Line"]."g (pu)"[k])*(e[SysData["Line"]."To"[k]]*e[SysData["Line"]."From"[k]]+f[SysData["Line"]."To"[k]]*f[SysData["Line"]."From"[k]])-(SysData["Line"]."b (pu)"[k])*(f[SysData["Line"]."To"[k]]*e[SysData["Line"]."From"[k]]-e[SysData["Line"]."To"[k]]*f[SysData["Line"]."From"[k]])) #Line flow expression 
@expression(ACOPF, pji, pflow_ji_sh + pflow_ji_sr1 + pflow_ji_sr2)

@expression(ACOPF, qflow_ij_sh[k=1:nline],(-SysData["Line"]."b_sh (pu)"[k]/2)*(e[SysData["Line"]."From"[k]]^2+f[SysData["Line"]."From"[k]]^2)) #Line flow expression 
@expression(ACOPF, qflow_ij_sr1[k=1:nline],(-SysData["Line"]."b (pu)"[k])*(e[SysData["Line"]."From"[k]]^2+f[SysData["Line"]."From"[k]]^2)) #Line flow expression 
@expression(ACOPF, qflow_ij_sr2[k=1:nline],(SysData["Line"]."b (pu)"[k])*(e[SysData["Line"]."From"[k]]*e[SysData["Line"]."To"[k]]+f[SysData["Line"]."From"[k]]*f[SysData["Line"]."To"[k]])-(SysData["Line"]."g (pu)"[k])*(f[SysData["Line"]."From"[k]]*e[SysData["Line"]."To"[k]]-e[SysData["Line"]."From"[k]]*f[SysData["Line"]."To"[k]])) #Line flow expression 
@expression(ACOPF, qij, qflow_ij_sh + qflow_ij_sr1 + qflow_ij_sr2)

@expression(ACOPF, qflow_ji_sh[k=1:nline],(-SysData["Line"]."b_sh (pu)"[k]/2)*(e[SysData["Line"]."To"[k]]^2+f[SysData["Line"]."To"[k]]^2)) #Line flow expression 
@expression(ACOPF, qflow_ji_sr1[k=1:nline],(-SysData["Line"]."b (pu)"[k])*(e[SysData["Line"]."To"[k]]^2+f[SysData["Line"]."To"[k]]^2)) #Line flow expression 
@expression(ACOPF, qflow_ji_sr2[k=1:nline],(SysData["Line"]."b (pu)"[k])*(e[SysData["Line"]."To"[k]]*e[SysData["Line"]."From"[k]]+f[SysData["Line"]."To"[k]]*f[SysData["Line"]."From"[k]])-(SysData["Line"]."g (pu)"[k])*(f[SysData["Line"]."To"[k]]*e[SysData["Line"]."From"[k]]-e[SysData["Line"]."To"[k]]*f[SysData["Line"]."From"[k]])) #Line flow expression 
@expression(ACOPF, qji, qflow_ji_sh + qflow_ji_sr1 + qflow_ji_sr2)


#-------------------------------Power balance constraint #OK
@expression(ACOPF, p_inj[i=1:nBus], sum(pij[k] for k in 1:nline if SysData["Line"]."From"[k]==i)+sum(pji[k] for k in 1:nline if SysData["Line"]."To"[k]==i))
@expression(ACOPF, q_inj[i=1:nBus], sum(qij[k] for k in 1:nline if SysData["Line"]."From"[k]==i)+sum(qji[k] for k in 1:nline if SysData["Line"]."To"[k]==i))
@expression(ACOPF, PD[i=1:nBus], sum(SysData["Load"]."Pd (pu)"[j] for j in 1:nload if SysData["Load"]."Bus_N"[j]==i))
@expression(ACOPF, QD[i=1:nBus], sum(SysData["Load"]."Qd (pu)"[j] for j in 1:nload if SysData["Load"]."Bus_N"[j]==i))

@constraint(ACOPF, Active_power_balance[i=1:nBus], sum(PG[j] for j in 1:nGens if SysData["Gen"]."Bus"[j]==i)-PD[i]==p_inj[i])
@constraint(ACOPF, Reactive_power_balance[i=1:nBus], sum(QG[j] for j in 1:nGens if SysData["Gen"]."Bus"[j]==i)-QD[i]==q_inj[i])

#-------------------------- Network constraints
@constraint(ACOPF,PCC[i=1:nBus;SysData["Bus"]."Bus_Type"[i]==3],f[i]==0)
@constraint(ACOPF,PCC_V[i=1:nBus;SysData["Bus"]."Bus_Type"[i]==3], SysData["Bus"]."V_min (pu)"[i]<=e[i]<=SysData["Bus"]."V_max (pu)"[i])
@constraint(ACOPF,VConstraint[i=1:nBus;SysData["Bus"]."Bus_Type"[i]!=3], SysData["Bus"]."V_min (pu)"[i]^2<=e[i]^2+f[i]^2<=SysData["Bus"]."V_max (pu)"[i] ^2)
@constraint(ACOPF,Gen_Limit_P[i=1:nGens],SysData["Gen"]."P_min (pu)"[i] <= PG[i] <= SysData["Gen"]."P_max (pu)"[i])
@constraint(ACOPF,Gen_Limit_Q[i=1:nGens],SysData["Gen"]."Q_min (pu)"[i] <= QG[i] <= SysData["Gen"]."Q_max (pu)"[i])
@expression(ACOPF,Iflow2_ij[k=1:nline], (SysData["Line"]."g (pu)"[k]^2 + SysData["Line"]."b (pu)"[k]^2)*((e[SysData["Line"]."From"[k]]-e[SysData["Line"]."To"[k]])^2+(f[SysData["Line"]."From"[k]]-f[SysData["Line"]."To"[k]])^2))# lineflow current limitation
@constraint(ACOPF,flow[k=1:nline],Iflow2_ij[k] <= SysData["Line"]."RATE_A"[k]^2)



#----------------Objective function #OK

@expression(ACOPF,cost_gen[i=1:nGens], SysData["Gen"]."cost_c"[i]+ SysData["Gen"]."cost_b"[i]*sbase*PG[i] + SysData["Gen"]."cost_a"[i]*sbase^2*PG[i]*PG[i])
@objective(ACOPF, Min, sum(cost_gen[i] for i=1:nGens)) # Setting objective function


#----------------Solve and Results

optimize!(ACOPF)
println("Objective value:  ", JuMP.objective_value(ACOPF))
println("Solver Time:   ", JuMP.solve_time(ACOPF))
value.(ACOPF[:PG])

e1=value.(ACOPF[:e])
f1=value.(ACOPF[:f])
V1=zeros(nBus)
V1=sqrt.(e1.^2+f1.^2)

θ1=zeros(nBus)
θ1=atan.(f1./e1)

return V1,θ1,value.(ACOPF[:PG]),JuMP.objective_value(ACOPF)
end


function Find_θ_SOCP(SCOPF,maxerror=1e-4)
    SCOPF[:A]
    cL=1;t=1
    for t=1:nTP, c=cL
        A_res=value.(SCOPF[:A][t,:,c])
        B_res=value.(SCOPF[:B][t,:,c])
        V2_res=value.(SCOPF[:V2][t,:,c])

        Recover = Model(Ipopt.Optimizer)
        @variables(
            Recover,
            begin
                es[k = 1:size(SysData["Bus"],1)]
                fs[k = 1:size(SysData["Bus"],1)]
                ϵ1[k=1:size(SysData["Line"],1)]
                ϵ2[k=1:size(SysData["Line"],1)]
                ϵ3[k=1:size(SysData["Bus"],1)]
            end
        )

        @constraint(Recover, first[k=1:size(SysData["Line"],1)], es[SysData["Line"]."From"[k]]*es[SysData["Line"]."To"[k]]+fs[SysData["Line"]."From"[k]]*fs[SysData["Line"]."To"[k]] == A_res[k]+ϵ1[k])
        @constraint(Recover, second[k=1:size(SysData["Line"],1)], fs[SysData["Line"]."From"[k]]*es[SysData["Line"]."To"[k]]-es[SysData["Line"]."From"[k]]*fs[SysData["Line"]."To"[k]] == B_res[k]+ϵ2[k])
        @constraint(Recover, Third[k=1:size(SysData["Bus"],1)], es[k]^2+fs[k]^2 == V2_res[k]+ϵ3[k])
        @constraint(Recover, PCC[k=1:size(SysData["Bus"],1) ; SysData["Bus"]."Bus_Type"[k]==3],fs[k]==0 )
        @constraint(Recover,limit[k=1:size(SysData["Bus"],1)], es[k]>=0)

        @objective(Recover,Min, sum(ϵ1.^2)+sum(ϵ2.^2)+sum(ϵ3.^2))
        optimize!(Recover)
        println("Objective value:  ", JuMP.objective_value(Recover))
        if JuMP.objective_value(Recover)>=maxerror
            converge=[false,JuMP.objective_value(Recover)]
        else
            converge=[true,JuMP.objective_value(Recover)]
        end

        SCOPF[:θ]=atan.(value.(Recover[:fs])./value.(Recover[:es]))
        SCOPF[:θ][5]=0
        V1=sqrt.(value.(LACOPF[:V2]))
        PG1=value.(LACOPF[:PG])
    end

    return V1,θ1,PG1
end

