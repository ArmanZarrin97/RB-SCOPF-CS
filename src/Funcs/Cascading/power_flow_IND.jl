
function multislack_acpf_order_pmaxqmax3(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order, ctd)
    
    eps_g=1e-4
    nBus = size(B, 1)
    ign = G[:, 1]
    Gtemp = deepcopy(G)
    exit = 0
    Gvio = 1
    ish = Shu[:, 1]
    QshuntBUS = zeros(Float64, nBus)
    QshuntBUS[ish]  = Shu[:, 2]

    (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i, vior_i, Iij, percentIij, is_slack, Jcb) = multislack_acpf_order_V2(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order, ctd)
    while exit == 0 && cg_t==1
        # G[:, 2] = p_i[ign] + p_V[ign]
        # G[:, 3] = q_i[ign] + q_V[ign] - QshuntBUS[ign]
        vio_g1 = find_beyond(G[:, 2], G[:, 7], eps_g)
        vio_g2 = find_below(G[:, 2], G[:, 8], eps_g)
        vio_g = union(vio_g1, vio_g2)
        # G[vio_g, :]
        # sum(G[:, 2])
        # sum(G[:, 7])
        # OP["QG_c_Line"][13]
        # G[13, 2]
        # non_vio_g = setdiff(collect(1:size(G, 1)), vio_g)
        non_vio_g = setdiff(findall(==(1),G.slack), vio_g)
        Gvio = length(vio_g)

        if cg_t == 1 && Gvio == 0
            exit = 1
            # println("generataor limits are satisfied")
            break
        else
            println("P are not in range")
            G[vio_g1, 2] = G[vio_g1, 7]
            G[vio_g2, 2] = G[vio_g2, 8]
            G[non_vio_g, 2] = Gtemp[non_vio_g, 2]
            G[vio_g, 12] = zeros(length(vio_g))
        end
        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i, vior_i, Iij, percentIij, is_slack, Jcb) = multislack_acpf_order_V2(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order, ctd)
    end

    return vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i, vior_i, Iij, percentIij, is_slack, Jcb
end

function multislack_acpf_order_V2(y_bus,
    B, L, G, Ld, Shu, REN, BAT, order, ctd;print=false)

    # y_bus=y_bus_i
    

    #----------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 1: Settings--[-][-][-][-][-][-][-]---------
    #----------------------------------------------------------------------------------------------
    nBus = size(B, 1)
    PgenBUS = zeros(Float64, nBus)
    QgenBUS = zeros(Float64, nBus)
    PloadBUS = zeros(Float64, nBus)
    QloadBUS = zeros(Float64, nBus)
    QshuntBUS = zeros(Float64, nBus)
    QgenBUSmin = zeros(Float64, nBus)
    QgenBUSmax = zeros(Float64, nBus)
    p_inj = zeros(Float64, nBus)
    q_inj = zeros(Float64, nBus)
    p_V =  zeros(Float64, nBus)
    q_V =  zeros(Float64, nBus)


    dp  = zeros(Float64, nBus, 3)
    dq  = zeros(Float64, nBus, 3)
    kip = zeros(Float64, nBus)
    kiq = zeros(Float64, nBus)
    Gs = ones(Float64, nBus)
    Bs = ones(Float64, nBus)
    Gs0 = ones(Float64, nBus)
    Bs0 = ones(Float64, nBus)
    Vi0 = ones(Float64, nBus)
    ind = zeros(Float64, nBus)
    diff_Tf  = zeros(Float64, nBus)
    diff_Tfm = zeros(Float64, nBus)
    diff_TV  = zeros(Float64, nBus)
    diff_Gf  = zeros(Float64, nBus)
    diff_Gfm = zeros(Float64, nBus)
    diff_Bf  = zeros(Float64, nBus)
    diff_Bfm = zeros(Float64, nBus)

    ign = G[:, 1]
    isG=0
    ipq = findall(x -> x == 1, B[:, 2])
    ipv = findall(x -> x != 1, B[:, 2])
    isl = findall(x -> x == 3, B[:, 2])

     
    if isl ==[]
        (xxx,isG)=findmax(G[:, 7])
        isl=G[isG,1]
    else
        isl=isl[1]
    end
    
    iload = Ld[:, 1]
    idu = findall(x -> x != 0, Ld[:, 23]) 
    Ld_ind=Ld[idu, 15:20]
    idu = iload[idu] # index of buses that induction machine are connected
    ish = Shu[:, 1]
    slk = deepcopy(G[:, 12])
   
    global slk1=deepcopy(slk)

    PgenBUS[ign] = G[:, 2]
   



    QgenBUSmax[ign] = G[:, 4]
    QgenBUSmin[ign] = G[:, 5]
    PloadBUS[iload] = Ld[:, 2]
    QloadBUS[iload] = Ld[:, 3]
    QshuntBUS[ish]  = Shu[:, 2]
    Vi0[iload]      = Ld[:, 13]
    ind[iload]      = Ld[:, 23]
    Gs0[iload]      = Ld[:, 21]
    Bs0[iload]      = Ld[:, 22]
    kip[iload]      = Ld[:, 24]
    kiq[iload]      = Ld[:, 25]

    dp[iload, 1:3].=Ld[:, 7:9]
    dq[iload, 1:3].=Ld[:, 10:12]

    PloadBUS[REN[:,2]] -= REN[:, 5]
    PloadBUS[BAT[:,1]] -= BAT[:, 26]
    QloadBUS[REN[:,2]] -= REN[:, 6]
    QloadBUS[BAT[:,1]] -= BAT[:, 27]
    dp[REN[:,2], 3].=1
    dq[REN[:,2], 3].=1
    dp[BAT[:,1], 3].=1
    dq[BAT[:,1], 3].=1


    npv = length(ipv)
    npq = length(ipq)

    max_mismatch_pf = 1     #
    epsilon = 1e-5        #convergence criteria
    epsilon_suf = 1e-5
    epsilon_suf_p = 1e-5
    itr_max = 100       # max number of iterations
    maxiter = []


    v_mag_pf = 1 * ones(Float64, nBus)
    v_angle_pf = zeros(Float64, nBus)
    v_mag_pf[ign] = G[:, 6]

   
    #.............................................................................
    # Note: In general, all (ipv=ipv+isl) buses can be considered as the slack bus
    #.............................................................................
    df = 0.00
    ki = zeros(Float64, nBus + npq + length(idu), 1)
    ki = spc_slack2(ki, G)
   
    #.............................................................................
    # Intitalize motro slips, calculate motor admitance
    #.............................................................................
    #.............................................................................
   
   
   sm = 0.001 * ones(length(idu))
   fm = (1 .-sm) * (1 + df)
    #.............................................................................
    max_mismatch_pf = 1
    iteration = 0
    mdf_idx = zeros(Int64, nBus)
    mdf=1
    convergence_indicator = 0
    delta_mismatch_inj = zeros(Float64, nBus + npq + length(idu))
    delta_mismatch=[]
    # Jacob1 = zeros(Float64, nBus + npq + length(idu), nBus + npq + length(idu))
    #---------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 2: Calcultae Scheduled Powers   [-][-][-][-][-][-][-]
    #---------------------------------------------------------------------------------------------------
    (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)

    #---------------------------------------------------------------------------------------------------


    #---------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 3: Main Loop of the Powr Flow--[-][-][-][-][-][-][-][-]
    #---------------------------------------------------------------------------------------------------
    # while (max_mismatch_pf > epsilon && iteration <= itr_max) || mdf!=0
    while max_mismatch_pf > epsilon && iteration <= itr_max


        # .......Step 3-1: Mismatch between  Scheduled and Calculated Powers............
        #...............................................................................
        (delta_mismatch_inj, p_inj, q_inj, diff_Tf1, diff_Tfm1, diff_TV1, diff_Gf1, diff_Gfm1, diff_Bf1, diff_Bfm1, p_V, q_V) = mismatch_vectorV_orderINDA2(nBus, isl, ipv, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd)
        diff_Tf[idu]  = diff_Tf1
        diff_Tfm[idu] = diff_Tfm1
        diff_TV[idu]  = diff_TV1
        diff_Gf[idu]  = diff_Gf1
        diff_Gfm[idu] = diff_Gfm1
        diff_Bf[idu]  = diff_Bf1
        diff_Bfm[idu] = diff_Bfm1
       
        delta_mismatch_inj =  convert(Vector{Float64}, delta_mismatch_inj)
        # .......Step 3-2: Jacobian Matrix.............................................
        #..............................................................................
        Jacob1 = jacobianV_orderINDA(nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, df, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd)
        
        # .......Step 3-3: Inverse of Jacobian Matrix...................................
        #...............................................................................
        Jcb = lu(Jacob1, check=false)
        # dense_Jacob1 = Array(Jacob1)
        # conditioning=cond(dense_Jacob1)
        # printstyled("conditioning factor is: $conditioning:"; color=:red, blink=true)
        if !issuccess(Jcb)
            if print
                println("")
                show("Sigularity observed for iteration= $iteration")
                println("")
            end
            break
        end

        # .......Step 3-4: Update State Variables (V, delta,Pim)..........................
        #.................................................................................
        delta_x_pf = Jcb \ delta_mismatch_inj
        (v_angle_sv, v_mag_sv, df1, fm1) = solution_voltageIND(nBus, isl, ipq, ipv, v_angle_pf, v_mag_pf, df, fm, delta_x_pf)
        v_mag_pf = v_mag_sv                                                  # sv = solution voltage (Here, voltage is updated)
        v_angle_pf = v_angle_sv
        df = df1
        fm = fm1

        (Gm, Bm) = y_motor(Ld_ind, df, fm)
        Gs[idu] = Gm
        Bs[idu] = Bm

        # .......Step 3-5: Check if PV-PQ Switching required ............................
        #.................................................................................
        max_mismatch_pf = maximum(abs.(delta_mismatch_inj))
        max_mismatch_pf_P = maximum(abs.(delta_mismatch_inj[1:nBus]))
        max_mismatch_pf_Q = maximum(abs.(delta_mismatch_inj[nBus+1:nBus+npq]))
        # max_mismatch_pf_T = maximum(abs.(delta_mismatch_inj[nBus+length(ipq)+1: end]))

        if max_mismatch_pf < epsilon
            (mdf_idx, QgenBUS, ipv, ipq, G) = pv_pq_switching1_order(epsilon_suf, y_bus, q_inj, v_mag_pf, v_angle_pf, df, ipv, ipq, mdf_idx, QgenBUS, kiq, QloadBUS, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G, order, Bs, Bs0, Vi0, ind)
            mdf = sum(mdf_idx)
            if mdf != 0
                # println("pq switching occured in iteration $iteration")
                npv = length(ipv)
                npq = length(ipq)
                ki = zeros(Float64, nBus + npq + length(idu), 1)
                ki = spc_slack2(ki, G)
                (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)
            end
        end

        
        iteration = iteration + 1
    #    push!(delta_mismatch,maximum(delta_mismatch_inj'))
    #    delta_mismatch_inj_dis = round.(delta_mismatch_inj', digits=4)
    #    println(delta_mismatch_inj_dis)
    # println(maximum(abs.(delta_mismatch_inj')))

    end # End while loop

    JacobOut=[]
    if !(max_mismatch_pf > epsilon && iteration <= itr_max)
        JacobOut =  jacobianV_orderINDA(nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, df, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd; flag=1)
        # Jacob=jacobianV_orderIND(nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd)
    end
    #--------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 4: Power Flow Results     [-][-][-][-][-][-][-][-][-]
    #--------------------------------------------------------------------------------------------------

    # .......Step 4-1: Report the Convergence Status..................................
    #.................................................................................
    push!(maxiter, iteration)
    if max_mismatch_pf < epsilon && iteration < itr_max
        convergence_indicator = 1
    else
        if print
            println("")
            show("The power flow has NOT converged in the required number of iteration")
            println("")
        end
    end

                                                  
    (delta_mismatch_inj, p_inj, q_inj, diff_Tf1, diff_Tfm1, diff_TV1, diff_Gf1, diff_Gfm1, diff_Bf1, diff_Bfm1, p_V, q_V) = mismatch_vectorV_orderINDA2(nBus, isl, ipv, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd)
   
    (sij, violated_lines, atrisk_lines, Iij, percentIij) = line_flow2(L, y_bus, v_mag_pf, v_angle_pf)
   
   
    pim_pf = - sum(ki) * df

    sm=(1 + df .-fm) / (1 + df)

    is_slack  = findall(x -> x == 1, G[:, 12])

    G[:, 12]=deepcopy(slk)

    return v_mag_pf, v_angle_pf, delta_mismatch_inj, delta_mismatch, p_inj, q_inj, p_V, q_V, convergence_indicator, pim_pf, df, fm, sm, Gs, Bs, iteration,sij, violated_lines, atrisk_lines, Iij, percentIij, is_slack, JacobOut

end


function multislack_acpf_order2(y_bus,
    B, L, G, Ld, Shu, order, ctd)
    

    #----------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 1: Settings--[-][-][-][-][-][-][-]---------
    #----------------------------------------------------------------------------------------------
    nBus = size(B, 1)
    PgenBUS = zeros(Float64, nBus)
    QgenBUS = zeros(Float64, nBus)
    PloadBUS = zeros(Float64, nBus)
    QloadBUS = zeros(Float64, nBus)
    QshuntBUS = zeros(Float64, nBus)
    QgenBUSmin = zeros(Float64, nBus)
    QgenBUSmax = zeros(Float64, nBus)
    p_inj = zeros(Float64, nBus)
    q_inj = zeros(Float64, nBus)
    p_V =  zeros(Float64, nBus)
    q_V =  zeros(Float64, nBus)


    dp  = zeros(Float64, nBus, 3)
    dq  = zeros(Float64, nBus, 3)
    kip = zeros(Float64, nBus)
    kiq = zeros(Float64, nBus)
    Gs = ones(Float64, nBus)
    Bs = ones(Float64, nBus)
    Gs0 = ones(Float64, nBus)
    Bs0 = ones(Float64, nBus)
    Vi0 = ones(Float64, nBus)
    ind = zeros(Float64, nBus)
    diff_Tf  = zeros(Float64, nBus)
    diff_Tfm = zeros(Float64, nBus)
    diff_TV  = zeros(Float64, nBus)
    diff_Gf  = zeros(Float64, nBus)
    diff_Gfm = zeros(Float64, nBus)
    diff_Bf  = zeros(Float64, nBus)
    diff_Bfm = zeros(Float64, nBus)

    ign = G[:, 1]
    isG=0
    ipq = findall(x -> x == 1, B[:, 2])
    ipv = findall(x -> x != 1, B[:, 2])
    isl = findall(x -> x == 3, B[:, 2])

     
    if isl ==[]
        (xxx,isG)=findmax(G[:, 7])
        isl=G[isG,1]
    else
        isl=isl[1]
    end
    
    iload = Ld[:, 1]
    idu = findall(x -> x != 0, Ld[:, 23]) 
    Ld_ind=Ld[idu, 15:20]
    idu = iload[idu] # index of buses that induction machine are connected
    ish = Shu[:, 1]
    slk = deepcopy(G[:, 12])
   
    global slk1=deepcopy(slk)

    PgenBUS[ign] = G[:, 2]
   



    QgenBUSmax[ign] = G[:, 4]
    QgenBUSmin[ign] = G[:, 5]
    PloadBUS[iload] = Ld[:, 2]
    QloadBUS[iload] = Ld[:, 3]
    QshuntBUS[ish]  = Shu[:, 2]
    Vi0[iload]      = Ld[:, 13]
    ind[iload]      = Ld[:, 23]
    Gs0[iload]      = Ld[:, 21]
    Bs0[iload]      = Ld[:, 22]
    kip[iload]      = Ld[:, 24]
    kiq[iload]      = Ld[:, 25]

    dp[iload, 1:3].=Ld[:, 7:9]
    dq[iload, 1:3].=Ld[:, 10:12]


    npv = length(ipv)
    npq = length(ipq)

    max_mismatch_pf = 1     #
    epsilon = 1e-4        #convergence criteria
    epsilon_suf = 1e-3
    epsilon_suf_p = 1e-3
    itr_max = 100       # max number of iterations
    maxiter = []


    v_mag_pf = 1 * ones(Float64, nBus)
    v_angle_pf = zeros(Float64, nBus)
    v_mag_pf[ign] = G[:, 6]

   
    #.............................................................................
    # Note: In general, all (ipv=ipv+isl) buses can be considered as the slack bus
    #.............................................................................
    df = 0.00
    ki = zeros(Float64, nBus + npq + length(idu), 1)
    ki = spc_slack2(ki, G)
   
    #.............................................................................
    # Intitalize motro slips, calculate motor admitance
    #.............................................................................
    #.............................................................................
   
   
   sm = 0.001 * ones(length(idu))
   fm = (1 .-sm) * (1 + df)
    #.............................................................................
    max_mismatch_pf = 1
    iteration = 0
    mdf_idx = zeros(Int64, nBus)
    mdf=1
    convergence_indicator = 0
    delta_mismatch_inj = zeros(Float64, nBus + npq + length(idu))
    delta_mismatch=[]
    # Jacob1 = zeros(Float64, nBus + npq + length(idu), nBus + npq + length(idu))
    #---------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 2: Calcultae Scheduled Powers   [-][-][-][-][-][-][-]
    #---------------------------------------------------------------------------------------------------
    (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)

    #---------------------------------------------------------------------------------------------------


    #---------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 3: Main Loop of the Powr Flow--[-][-][-][-][-][-][-][-]
    #---------------------------------------------------------------------------------------------------
    # while (max_mismatch_pf > epsilon && iteration <= itr_max) || mdf!=0
    while max_mismatch_pf > epsilon && iteration <= itr_max


        # .......Step 3-1: Mismatch between  Scheduled and Calculated Powers............
        #...............................................................................
        (delta_mismatch_inj, p_inj, q_inj, diff_Tf1, diff_Tfm1, diff_TV1, diff_Gf1, diff_Gfm1, diff_Bf1, diff_Bfm1, p_V, q_V) = mismatch_vectorV_orderINDA2(nBus, isl, ipv, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd)
        diff_Tf[idu]  = diff_Tf1
        diff_Tfm[idu] = diff_Tfm1
        diff_TV[idu]  = diff_TV1
        diff_Gf[idu]  = diff_Gf1
        diff_Gfm[idu] = diff_Gfm1
        diff_Bf[idu]  = diff_Bf1
        diff_Bfm[idu] = diff_Bfm1
       
        delta_mismatch_inj =  convert(Vector{Float64}, delta_mismatch_inj)
        # .......Step 3-2: Jacobian Matrix.............................................
        #..............................................................................
        Jacob1 = jacobianV_orderINDA(nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, df, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd)

        # .......Step 3-3: Inverse of Jacobian Matrix...................................
        #...............................................................................
        Jcb = lu(Jacob1, check=false)
        # dense_Jacob1 = Array(Jacob1)
        # conditioning=cond(dense_Jacob1)
        # printstyled("conditioning factor is: $conditioning:"; color=:red, blink=true)
        if !issuccess(Jcb)
            println("")
            show("Sigularity observed for iteration= $iteration")
            println("")
            # break
        end

        # .......Step 3-4: Update State Variables (V, delta,Pim)..........................
        #.................................................................................
        delta_x_pf = Jcb \ delta_mismatch_inj
        (v_angle_sv, v_mag_sv, df1, fm1) = solution_voltageIND(nBus, isl, ipq, ipv, v_angle_pf, v_mag_pf, df, fm, delta_x_pf)
        v_mag_pf = v_mag_sv                                                  # sv = solution voltage (Here, voltage is updated)
        v_angle_pf = v_angle_sv
        df = df1
        fm = fm1

        (Gm, Bm) = y_motor(Ld_ind, df, fm)
        Gs[idu] = Gm
        Bs[idu] = Bm

        # .......Step 3-5: Check if PV-PQ Switching required ............................
        #.................................................................................
        max_mismatch_pf = maximum(abs.(delta_mismatch_inj))
        max_mismatch_pf_P = maximum(abs.(delta_mismatch_inj[1:nBus]))
        max_mismatch_pf_Q = maximum(abs.(delta_mismatch_inj[nBus+1:nBus+npq]))
        # max_mismatch_pf_T = maximum(abs.(delta_mismatch_inj[nBus+length(ipq)+1: end]))

        if max_mismatch_pf_P < epsilon_pq
            (mdf_idx, QgenBUS, ipv, ipq, G) = pv_pq_switching1_order(epsilon_suf, y_bus, q_inj, v_mag_pf, v_angle_pf, df, ipv, ipq, mdf_idx, QgenBUS, kiq, QloadBUS, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G, order, Bs, Bs0, Vi0, ind)
            mdf = sum(mdf_idx)
            if mdf != 0
                # println("pq switching occured in iteration $iteration")
                npv = length(ipv)
                npq = length(ipq)
                ki = zeros(Float64, nBus + npq + length(idu), 1)
                ki = spc_slack2(ki, G)
                (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)
            end
        end

        
        iteration = iteration + 1
    #    push!(delta_mismatch,maximum(delta_mismatch_inj'))
    #    delta_mismatch_inj_dis = round.(delta_mismatch_inj', digits=4)
    #    println(delta_mismatch_inj_dis)
    # println(maximum(abs.(delta_mismatch_inj')))

    end # End while loop
    
    if !(max_mismatch_pf > epsilon && iteration <= itr_max)
        JacobOut =  jacobianV_orderINDA(nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, df, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd; flag=1)
        # Jacob=jacobianV_orderIND(nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd)
    end

    #--------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 4: Power Flow Results     [-][-][-][-][-][-][-][-][-]
    #--------------------------------------------------------------------------------------------------

    # .......Step 4-1: Report the Convergence Status..................................
    #.................................................................................
    push!(maxiter, iteration)
    if max_mismatch_pf < epsilon && iteration < itr_max
        convergence_indicator = 1
    else
        println("")
        show("The power flow has NOT converged in the required number of iteration")
        println("")
    end

                                                  
    (delta_mismatch_inj, p_inj, q_inj, diff_Tf1, diff_Tfm1, diff_TV1, diff_Gf1, diff_Gfm1, diff_Bf1, diff_Bfm1, p_V, q_V) = mismatch_vectorV_orderINDA2(nBus, isl, ipv, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd)
   
    (sij, violated_lines, Iij, percentIij) = line_flow2(L, y_bus, v_mag_pf, v_angle_pf)
   
   
    pim_pf = - sum(ki) * df

    sm=(1 + df .-fm) / (1 + df)

    is_slack  = findall(x -> x == 1, G[:, 12])

    G[:, 12]=deepcopy(slk)

    return v_mag_pf, v_angle_pf, delta_mismatch_inj, delta_mismatch, p_inj, q_inj, p_V, q_V, convergence_indicator, pim_pf, df, fm, sm, Gs, Bs, iteration,sij, violated_lines, Iij, percentIij, is_slack, JacobOut

end












