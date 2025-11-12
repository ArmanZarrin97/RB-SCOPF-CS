function multislack_acpf_extended(y_bus,
    B, L, G, Ld, Shu, order, ctd)

    #----------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 1: Settings--[-][-][-][-][-][-][-]---------
    #----------------------------------------------------------------------------------------------
    nBus = size(B, 1)
    PgenBUS = zeros(Float64, nBus)
    QgenBUS = zeros(Float64, nBus)
    kv = zeros(Float64, nBus)
    v_set = ones(Float64, nBus)
    X1d = ones(Float64, nBus)
    q_set = zeros(Float64, nBus)
    PloadBUS = zeros(Float64, nBus)
    QloadBUS = zeros(Float64, nBus)
    QshuntBUS = zeros(Float64, nBus)
    QgenBUSmin = zeros(Float64, nBus)
    QgenBUSmax = zeros(Float64, nBus)
    IgenBUSmax = zeros(Float64, nBus)
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
    ign_t = G[:, 27]
    isG=0
    ipq = findall(x -> x == 1, B[:, 2])
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
   



    QgenBUSmax[ign_t] = G[:, 4]
    QgenBUSmin[ign_t] = G[:, 5]
    IgenBUSmax[ign_t] = G[:, 31]
    PloadBUS[iload] = Ld[:, 2]
    QloadBUS[iload] = Ld[:, 3]
    QshuntBUS[ish]  = Shu[:, 2]
    Vi0[iload]      = Ld[:, 13]
    ind[iload]      = Ld[:, 23]
    Gs0[iload]      = Ld[:, 21]
    Bs0[iload]      = Ld[:, 22]
    kip[iload]      = Ld[:, 24]
    kiq[iload]      = Ld[:, 25]

    dp[iload, 1:3]=Ld[:, 7:9]
    dq[iload, 1:3]=Ld[:, 10:12]
    
    Sbm = G[:, 24]
    Sb = 100
    X1d[ign_t] = G[:, 18] .* (Sb ./ Sbm)

    ish_fx_idx = findall(x -> x == 1, Shu[:, 5])
    ish_fx = Shu[ish_fx_idx, 1]

    # igm_idx = findall(x -> x == 2, G[:, 29])
    # igl_idx = findall(x -> x == 3, G[:, 29])
    # igm_t = G[igm_idx, 27]
    # igl_t = G[igl_idx, 27]
    # igm   =  igm_t .+ length(ign)
    # igl   =  igl_t .+ length(ign)

   
    #.............................................................................
    # Note: In general, all generation buses can be considered as frequency and/or voltage control buses
    #.............................................................................
    kv[ign_t]    = G[:, 28]
    q_set[ign_t] = G[:, 3]
    v_set[ign_t] = G[:, 30]

    fx_idx = findall(x -> x == 1, G[:, 29])
    vr_idx = findall(x -> x != 1, G[:, 29])
    ign_fx = G[fx_idx, 1]
    ign_vr = G[vr_idx, 1]
    IPQ= sort(union(ign_vr,ipq))
    IPQ = convert(Vector{Int64}, IPQ)
    NPQ = length(IPQ)

    ki = zeros(Float64, nBus + NPQ + length(idu), 1)
    ki = spc_slack(ki, G)
   
     #.............................................................................
    # Intitalize voltages, angles, and frequency
    #.............................................................................
    #.............................................................................
    v_mag_pf = 1 * ones(Float64, nBus)
    v_angle_pf = 0.0001 * rand(nBus)

    v_mag_pf[ign_fx] = G[fx_idx, 6]
    v_mag_pf[ish_fx] = Shu[ish_fx_idx, 6]

    df = 0.0001

    #.............................................................................
    # Intitalize motor slips, calculate motor admitance
    #.............................................................................
    #.............................................................................
    sm = 0.001 * ones(length(idu))
    fm = (1 .-sm) * (1 + df)
    
    #.............................................................................
    # Power flow settings
    #.............................................................................
    #.............................................................................
    max_mismatch_pf = 1     #
    epsilon = 1e-4        #convergence criteria
    epsilon_suf = 1e-3
    epsilon_suf_p = 1e-3
    itr_max = 100       # max number of iterations
    maxiter = []

    max_mismatch_pf = 1
    iteration = 0
    mdf_idx = zeros(Int64, nBus)
    mdf=0
    convergence_indicator = 0
    delta_mismatch_inj = zeros(Float64, nBus + NPQ + length(idu))
    delta_mismatch=[]

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
        (delta_mismatch_inj, p_inj, q_inj, Pt, Qt, diff_Tf1, diff_Tfm1, diff_TV1, diff_Gf1, diff_Gfm1, diff_Bf1, diff_Bfm1, p_V, q_V) = mismatch_extended(nBus, ign, ign_vr, ign_t, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kv, v_set, q_set, X1d, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd)
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
        Jacob1 = jacobianV_extended(nBus, ign, ign_vr, isl, IPQ, ipq, idu, v_angle_pf, v_mag_pf, df, y_bus, p_inj, q_inj, Pt, Qt, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kv, X1d, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd)

        # .......Step 3-3: Inverse of Jacobian Matrix...................................
        #...............................................................................
        Jcb = lu(Jacob1, check=false)
        if !issuccess(Jcb)
            println("")
            show("Sigularity observed for iteration= $iteration")
            println("")
            break
        end

        # .......Step 3-4: Update State Variables (V, delta,df, fm)..........................
        #.................................................................................
        delta_x_pf = Jcb \ delta_mismatch_inj
        (v_angle_sv, v_mag_sv, df1, fm1) = solution_voltage_extended(nBus, isl, IPQ, ipq, v_angle_pf, v_mag_pf, df, fm, delta_x_pf)
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
        max_mismatch_pf_Q = maximum(abs.(delta_mismatch_inj[nBus+1:nBus+NPQ]))
        # max_mismatch_pf_T = maximum(abs.(delta_mismatch_inj[nBus+length(ipq)+1: end]))

        if max_mismatch_pf_P < epsilon_pq
        # if max_mismatch_pf < epsilon
            (mdf_idx, G, v_mag_pf) = pv_pq_switching_extended(epsilon_suf, mdf_idx, ign, ign_t, v_mag_pf, v_angle_pf, QgenBUSmax, QgenBUSmin, IgenBUSmax, G, X1d)
            mdf = sum(mdf_idx)
            if mdf != 0
                println("pq switching occured in iteration $iteration")
                kv[ign_t]    = G[:, 28]
                q_set[ign_t] = G[:, 3]
                v_set[ign_t] = G[:, 30]
                fx_idx = findall(x -> x == 1, G[:, 29])
                vr_idx = findall(x -> x != 1, G[:, 29])
                ign_fx = G[fx_idx, 1]
                ign_vr = G[vr_idx, 1]
                IPQ= sort(union(ign_vr,ipq))
                IPQ = convert(Vector{Int64}, IPQ)
                NPQ = length(IPQ)
                ki = zeros(Float64, nBus + NPQ + length(idu), 1)
                ki = spc_slack(ki, G)
                # (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)
            end
            mdf_idx = zeros(Int64, nBus)
        end

        
        iteration = iteration + 1
    #    push!(delta_mismatch,maximum(delta_mismatch_inj'))
    #    delta_mismatch_inj_dis = round.(delta_mismatch_inj', digits=4)
    #    println(delta_mismatch_inj_dis)
    println(maximum(abs.(delta_mismatch_inj')))

    end # End while loop

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

                                                  
    (delta_mismatch_inj, p_inj, q_inj, Pt, Qt, diff_Tf1, diff_Tfm1, diff_TV1, diff_Gf1, diff_Gfm1, diff_Bf1, diff_Bfm1, p_V, q_V) = mismatch_extended(nBus, ign, ign_vr, ign_t, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kv, v_set, q_set, X1d, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd)
   
    (sij, violated_lines) = line_flow1(L, y_bus, v_mag_pf, v_angle_pf)
   
   
    pim_pf = - sum(ki) * df

    sm=(1 + df .-fm) / (1 + df)

    G[:, 12]=deepcopy(slk)

    return v_mag_pf, v_angle_pf, delta_mismatch_inj, delta_mismatch, p_inj, q_inj, Qt, p_V, q_V, convergence_indicator, pim_pf, df, fm, sm, Gs, Bs, iteration,sij, violated_lines

end













