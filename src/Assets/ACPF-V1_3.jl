function multislack_acpf_order(y_bus,
    B, L, G, Ld, Shu, REN, BAT, order, ctd)


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

    PgenBUS[ipv] = G[:, 2]




    QgenBUSmax[ipv] = G[:, 4]
    QgenBUSmin[ipv] = G[:, 5]
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
    epsilon = 1e-4        #convergence criteria
    epsilon_suf = 1e-3
    epsilon_suf_p = 1e-3
    itr_max = 100       # max number of iterations
    maxiter = []


    v_mag_pf = 1 * ones(Float64, nBus)
    v_angle_pf = zeros(Float64, nBus)
    v_mag_pf[ipv] = G[:, 6]


    #.............................................................................
    # Note: In general, all (ipv=ipv+isl) buses can be considered as the slack bus
    #.............................................................................
    df = 0.00
    ki = zeros(Float64, npv + npq + npq + length(idu), 1)
    (ki, beti) = spc_slack(ki, G)

    #.............................................................................
    # Intitalize motro slips, calculate motor admitance
    #.............................................................................
    #.............................................................................


   sm = 0.001 * ones(length(idu))
   fm = (1 .-sm) * (1 + df)
    #.............................................................................
    max_mismatch_pf = 1
    iteration = 0
    mdf=0
    mdf_before=0
    mdf_idx = zeros(Int64, nBus)
    convergence_indicator = 0
    delta_mismatch_inj = zeros(Float64, nBus + npq + length(idu))
    delta_mismatch=[]

    #---------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 2: Calcultae Scheduled Powers   [-][-][-][-][-][-][-]
    #---------------------------------------------------------------------------------------------------
    (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)

    #---------------------------------------------------------------------------------------------------


    #---------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 3: Main Loop of the Powr Flow--[-][-][-][-][-][-][-][-]
    #---------------------------------------------------------------------------------------------------
    JacobOut=0;
    while (max_mismatch_pf > epsilon && iteration <= itr_max) || mdf>mdf_before

        # .......Step 3-1: Mismatch between  Scheduled and Calculated Powers............
        #...............................................................................
        (delta_mismatch_inj, p_inj, q_inj, diff_Tf1, diff_Tfm1, diff_TV1, diff_Gf1, diff_Gfm1, diff_Bf1, diff_Bfm1, p_V, q_V) = mismatch_vectorV_orderINDA(nBus, isl, ipv, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd)
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
        Jacob1 = jacobianV_orderIND(nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd)

        # .......Step 3-3: Inverse of Jacobian Matrix...................................
        #...............................................................................
        Jcb = lu(Jacob1, check=false)
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
        mdf_before=mdf
        max_mismatch_pf = maximum(abs.(delta_mismatch_inj))
        max_mismatch_pf_Q = maximum(abs.(delta_mismatch_inj[nBus+1:nBus+npq]))
        max_mismatch_pf_P = maximum(abs.(delta_mismatch_inj[1:nBus]))

        # max_mismatch_pf_T = maximum(abs.(delta_mismatch_inj[nBus+length(ipq)+1: end]))

        if max_mismatch_pf_P < epsilon
            (mdf_idx, QgenBUS, ipv, ipq, G) = pv_pq_switching1_order(epsilon_suf, y_bus, q_inj, v_mag_pf, v_angle_pf, ipv, ipq, mdf_idx, QgenBUS, QloadBUS, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G, order, Bs, Bs0, Vi0, ind)
            mdf = sum(mdf_idx)
            if mdf != 0
                println("pq switching occured in iteration $iteration")
                npv = length(ipv)
                npq = length(ipq)
                ki = zeros(npv + npq + npq + length(idu), 1)
                (ki, beti) = spc_slack(ki, G)
                (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)
            end
        end


        iteration = iteration + 1
        #    push!(delta_mismatch,maximum(delta_mismatch_inj'))
        #    delta_mismatch_inj_dis = round.(delta_mismatch_inj', digits=4)
        #    println(delta_mismatch_inj_dis)
        println(maximum(abs.(delta_mismatch_inj')))
     
    end # End while loop

    if !(max_mismatch_pf > epsilon && iteration <= itr_max)
        JacobOut = jacobianV_orderIND(nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd;flag=1)
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


    (delta_mismatch_inj, p_inj, q_inj, diff_Tf1, diff_Tfm1, diff_TV1, diff_Gf1, diff_Gfm1, diff_Bf1, diff_Bfm1, p_V, q_V) = mismatch_vectorV_orderINDA(nBus, isl, ipv, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd)

    (sij, violated_lines, Iij, percentIij) = line_flow1(L, y_bus, v_mag_pf, v_angle_pf)


    pim_pf = - sum(beti) * df

    sm=(1 + df .-fm) / (1 + df)

    G[:, 12]=deepcopy(slk)
    return v_mag_pf, v_angle_pf, delta_mismatch_inj, delta_mismatch, p_inj, q_inj, p_V, q_V, convergence_indicator, pim_pf, df, fm, sm, Gs, Bs, iteration,sij, violated_lines, Iij, percentIij, JacobOut

end
function Y_bus1(nBus, L, T, ctd)
    yij = zeros(ComplexF64, nBus, nBus)
    ybus = zeros(ComplexF64, nBus, nBus)
    fb = L[:, 1]
    tb = L[:, 2]
    g = L[:, 3]
    b = L[:, 4]
    bsh = L[:, 6]
    y = g .+ im * b

    # yij matrix
    for l in 1:length(fb)
        yij[fb[l], tb[l]] = yij[fb[l], tb[l]] + y[l]
        yij[tb[l], fb[l]] = yij[tb[l], fb[l]] + y[l]

        yij[fb[l], fb[l]] = yij[fb[l], fb[l]] + im * bsh[l] / 2
        yij[tb[l], tb[l]] = yij[tb[l], tb[l]] + im * bsh[l] / 2
    end


    for i in 1:nBus
        ctd_i = ctd[i]
        for j in ctd_i
            ybus[i, j] = -yij[i, j]
        end
        ybus[i, i] = sum(yij[i, ctd_i])
    end

    # ybus matrix with transformers
    fb_T = T[:, 1]
    tb_T = T[:, 2]
    g_T = T[:, 3]
    b_T = T[:, 4]
    y_T = g_T .+ im * b_T
    Tratio = 1.00 * T[:, 6]
    Tratio = 1 ./ Tratio # becasue tap is on the primary side

    for tr in 1:length(fb_T)
        ybus[fb_T[tr], fb_T[tr]] = ybus[fb_T[tr], fb_T[tr]] + abs(Tratio[tr])^2 * y_T[tr]
        ybus[fb_T[tr], tb_T[tr]] = ybus[fb_T[tr], tb_T[tr]] - conj(Tratio[tr]) * y_T[tr]
        ybus[tb_T[tr], fb_T[tr]] = ybus[tb_T[tr], fb_T[tr]] - Tratio[tr] * y_T[tr]
        ybus[tb_T[tr], tb_T[tr]] = ybus[tb_T[tr], tb_T[tr]] + 1 * y_T[tr]
    end

    return ybus
end




function A_net(nBus, L, T)
    A = spzeros(Int64, nBus, nBus)

    fb = L[:, 1]
    tb = L[:, 2]
    fb_T = T[:, 1]
    tb_T = T[:, 2]

    # aij matrix
    for l in 1:length(fb)
        A[fb[l], tb[l]] = A[fb[l], tb[l]] + 1
        A[tb[l], fb[l]] = A[tb[l], fb[l]] + 1
    end

    for i in 1:nBus
        A[i, i] = 1
    end

    # other connections due to transformers
    for tr in 1:length(fb_T)
        A[fb_T[tr], tb_T[tr]] = A[fb_T[tr], tb_T[tr]] + 1
        A[tb_T[tr], fb_T[tr]] = A[tb_T[tr], fb_T[tr]] + 1
    end

    return A
end




function connections(u)
    k = [1]
    uc = u[1]
    for i in 2:length(u)

        if length(intersect(u[1], u[i])) != 0
            uc = union(uc, u[i])
            k = union(k, i)
        end
    end

    return uc, k
end


function island6(Atemp)
    n = size(Atemp, 1)
    sum_n = 0.5 * n * (n + 1)
    sum_nu = 0.0
    Atemp0 = deepcopy(Atemp)
    flag_c = 0
    u = []
    unew = []
    h = 1

    # 1. Creat one set (initial clustering) corresponding to each row
    for i in 1:size(Atemp0, 1)
        ui = findall3(x -> x == 1, Atemp[i, :])
        push!(u, ui)
    end

    # 2. While there is no clustering (Sets have intersections)
    while flag_c != 1

        if length(u) == 1
            uc = u[1]
            k = 1
        else
            (uc, k) = connections(u)
        end

        # in any way
        push!(unew, uc)
        kres = setdiff(collect(1:size(u, 1)), k)
        u = u[kres]

        # if u is epmty then start with unew
        if length(u) == 0
            u = unew
            unew = []
            sum_nu = 0.0
        end

        # check the stopping criteria
        for i in 1:length(u)
            sum_nu += sum(u[i])
        end
        # println(sum_nu)

        if sum_nu == sum_n
            flag_c = 1
        end
        h += 1
    end

    return u, unew, h
end



function build_i(bi, B, L, T, G, Ld, Shu, REN, BAT, y_bus)
    y_bus_i = y_bus[bi, bi]

    # Cut required parts
    B = B[bi, :]
    L = L[find_avb(bi, L[:, 1]), :]
    T = T[find_avb(bi, T[:, 1]), :]
    G = G[find_avb(bi, G[:, 1]), :]
    Ld = Ld[find_avb(bi, Ld[:, 1]), :]
    Shu = Shu[find_avb(bi, Shu[:, 1]), :]

    # Rename bus numbers
    nBus = length(bi)
    B[:, 1] = re_number(bi, B[:, 1])

    L[:, 1] = re_number(bi, L[:, 1])
    L[:, 2] = re_number(bi, L[:, 2])

    T[:, 1] = re_number(bi, T[:, 1])
    T[:, 2] = re_number(bi, T[:, 2])

    G[:, 1] = re_number(bi, G[:, 1])
    Ld[:, 1] = re_number(bi, Ld[:, 1])
    Shu[:, 1] = re_number(bi, Shu[:, 1])


    return B, L, T, G, Ld, Shu, REN, BAT, y_bus_i

end


function re_number(r1, r2)
    # r1 is indeed bi
    # r2 is a new vector to be re-named
    r2_new = zeros(Int64, length(r2))
    for i in 1:length(r2)
        rr = findall3(x -> x == r2[i], r1)
        r2_new[i] = rr[1]
    end

    return r2_new

end

function assessV3(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
    # B = B_i
    # L = L_i
    # T = T_i
    # G = G_i
    # Ld = Ld_i
    # Shu = Shu_i

    Vsh = Vth * ones(size(Ld, 1))
    iload = Ld[:, 1]
    vmag_t = zeros(size(Ld, 1))
    Vvio_t = 0
    uflag = 0
    cg_t = 0
    cg_s = 0
    order_max = 2
    order_step = 1
    order_min = 0
    order = deepcopy(order_max)
    p_V = zeros(size(Ld, 1))
    q_V = zeros(size(Ld, 1))

    while true

        if cg_s == 0
            Ld = load_order(Ld, order)
            (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order)

            if cg_t == 0
                if order == order_max
                    println(" No short-term equilibrium")
                    kmain = Ld[:, 14]
                    LdI[kmain] = 2 * ones(length(kmain))
                    Ldsh_v[kmain] = Ldsh_v[kmain] + Ld[:, 2]
                    break
                else #~~~
                    priority = vmag_t
                    priority = sort(priority)
                    for i in 1:1
                        k = findall3(x -> x == priority[i], vmag_t)
                        k = k[1]
                        printstyled("There is a voltage vioaltion load $(Ld[k, 14]), {$(vmag_t[k])}:"; color=:blue, blink=true)

                        (Ld, uflag) = SV_alarm(Ld, vmag_t, k)
                        if uflag == 1
                            println(" No available load for shedding")
                            break
                        end
                        df_i = fr1(G, Ld)
                        (G, Ld) = updateG(G, Ld, df_i)
                        Ld = load_order(Ld, order)
                        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order)
                        if cg_t == 1
                            println(" GPFz has been converged, after $i single-load-shedding in short-term ")
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~
                            vmag_t = deepcopy(vmag_i[iload])
                            #~~~~~~~~
                            break
                        end
                    end
                end #~~~
            else
                println(" GPFz has been converged")
                vmag_dis = round.(vmag_i, digits=4)
                println(vmag_dis)

            end
            vmag_t = deepcopy(vmag_i[iload])
            Vvio_t = length(find_below(vmag_t, Vsh, eps_v))
            # bvio_t = find_below(vmag_t, Vsh, eps_v)
        end
        (p_V, q_V) = consumed_power(Ld, vmag_t, order)
        p_V_dis = round.(p_V, digits=4)
        q_V_dis = round.(q_V, digits=4)

        println(Ld[:, 14])
        println(p_V_dis)
        println(q_V_dis)

        if Vvio_t != 0
            println(" There is short-term voltage violation")
            order = 2
            alarm = 0
            (Ld, uflag) = SV(Ld, vmag_t, alarm)
            # (G, Ld) = updateGV(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
            df_i = fr1(G, Ld)
            (G, Ld) = updateG(G, Ld, df_i)
            if uflag == 1
                println(" No available load for shedding")
                break
            end

        else
            println(" No short-term voltage violation in order = $order")

            order_s = 0
            Ld = load_order(Ld, order_s)
            (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_s, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order_s)
            if cg_s == 0
                if order - order_step > order_min
                    println(" GPFp has been diverged, reduce order, and go back to short-term")
                    order = order - order_step
                    println(order)
                else
                    println(" GPFp has been diverged after convergence in last step")
                    alarm = 1
                    # order = order + order_step
                    priority = vmag_t
                    priority = sort(priority)
                    for i in 1:length(vmag_t)
                        k = findall3(x -> x == priority[i], vmag_t)
                        k = k[1]
                        printstyled("There is a voltage vioaltion load $(Ld[k, 14]), {$(vmag_t[k])}:"; color=:blue, blink=true)

                        (Ld, uflag) = SV_alarm(Ld, vmag_t, k)
                        if uflag == 1
                            println(" No available load for shedding")
                            break
                        end
                        df_i = fr1(G, Ld)
                        (G, Ld) = updateG(G, Ld, df_i)
                        order_s = 0
                        Ld = load_order(Ld, order_s)
                        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_s, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order_s)
                        if cg_s == 0
                            println(" GPFp has been converged, after $i single-load-shedding in the last short-term ")
                            break
                        end
                    end
                end
            else
                println(" GPFp has been converged,  check voltage violation")
                vmag_dis = round.(vmag_i, digits=4)
                println(vmag_dis)

                vmag_s = deepcopy(vmag_i[iload])
                Vvio_s = length(find_below(vmag_s, Vsh, eps_v))
                if Vvio_s != 0
                    println(" There is long-term voltage violation")
                    alarm = 0
                    (Ld, uflag) = SV(Ld, vmag_s, alarm)
                    Ld = load_s(Ld)
                    # (G, Ld) = updateGV(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
                    df_i = fr1(G, Ld)
                    (G, Ld) = updateG(G, Ld, df_i)
                    if uflag == 1
                        println(" No available load for shedding")
                        break
                    end
                else
                    println(" No long-term voltage violation")
                    break
                end
            end
        end
    end

    df_i = fr1(G, Ld)
    println(" after uvls,   the frequency is :       $df_i pu")

    return B, L, T, G, Ld, Shu
end







function assessV2(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
    # B = B_i
    # L = L_i
    # T = T_i
    # G = G_i
    # Ld = Ld_i
    # Shu = Shu_i

    Vsh = Vth * ones(size(Ld, 1))
    iload = Ld[:, 1]
    vmag_t = zeros(size(Ld, 1))
    Vvio_t = 0
    uflag = 0
    cg_t = 0
    cg_s = 0
    order_step = 1
    order = 2
    p_V = zeros(size(Ld, 1))
    q_V = zeros(size(Ld, 1))

    while true

        if cg_s == 0
            Ld = load_order(Ld, order)
            (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order)


            if cg_t == 0
                println(" No short-term equilibrium")
                kmain = Ld[:, 14]
                LdI[kmain] = 2 * ones(length(kmain))
                Ldsh_v[kmain] = Ldsh_v[kmain] + Ld[:, 2]
                break
            else
                println(" GPFz has been converged")
                vmag_dis = round.(vmag_i, digits=4)
                println(vmag_dis)

            end
            vmag_t = deepcopy(vmag_i[iload])
            Vvio_t = length(find_below(vmag_t, Vsh, eps_v))
            # bvio_t = find_below(vmag_t, Vsh, eps_v)
        end
        (p_V, q_V) = consumed_power(Ld, vmag_t, order)
        p_V_dis = round.(p_V, digits=4)
        q_V_dis = round.(q_V, digits=4)

        println(Ld[:, 14])
        println(p_V_dis)
        println(q_V_dis)

        if Vvio_t != 0
            println(" There is short-term voltage violation")
            order = 2
            alarm = 0
            (Ld, uflag) = SV(Ld, vmag_t, alarm)
            # (G, Ld) = updateGV(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
            df_i = fr1(G, Ld)
            (G, Ld) = updateG(G, Ld, df_i)
            if uflag == 1
                println(" No available load for shedding")
                break
            end

        else
            println(" No short-term voltage violation in order = $order")

            order_s = 0
            Ld = load_order(Ld, order_s)
            (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_s, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order_s)
            if cg_s == 0
                println(" GPFp has been diverged, reduce order, and go back to short-term")
                order = order - order_step
                println(order)
            else
                println(" GPFp has been converged,  check voltage violation")
                vmag_dis = round.(vmag_i, digits=4)
                println(vmag_dis)

                vmag_s = deepcopy(vmag_i[iload])
                Vvio_s = length(find_below(vmag_s, Vsh, eps_v))
                if Vvio_s != 0
                    println(" There is long-term voltage violation")
                    alarm = 0
                    (Ld, uflag) = SV(Ld, vmag_s, alarm)
                    Ld = load_s(Ld)
                    # (G, Ld) = updateGV(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
                    df_i = fr1(G, Ld)
                    (G, Ld) = updateG(G, Ld, df_i)
                    if uflag == 1
                        println(" No available load for shedding")
                        break
                    end
                else
                    println(" No long-term voltage violation")
                    break
                end
            end
        end
    end

    df_i = fr1(G, Ld)
    println(" after uvls,   the frequency is :       $df_i pu")

    return B, L, T, G, Ld, Shu
end



function assessV(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
    # B = B_i
    # L = L_i
    # T = T_i
    # G = G_i
    # Ld = Ld_i
    # Shu = Shu_i

    Vsh = Vth * ones(size(Ld, 1))
    iload = Ld[:, 1]
    vmag_t = zeros(size(Ld, 1))
    vmag_s = zeros(size(Ld, 1))
    vmag_i_t = zeros(size(B, 1))
    vmag_i_s = zeros(size(B, 1))
    vmag_i_temp = zeros(20, size(B, 1))
    temp = 1

    Vvio_t = 0
    Vvio_s = 0
    uflag = 0
    cg_t = 0
    cg_s = 0
    bvio_t = []
    while true

        if cg_s == 0
            order = 2
            Ld = load_order(Ld, order)
            (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order)
            println(sm_i')

            vmag_i_t = vmag_i
            if cg_t == 0
                println(" No short-term equilibrium")
                kmain = Ld[:, 14]
                LdI[kmain] = 2 * ones(length(kmain))
                Ldsh_v[kmain] = Ldsh_v[kmain] + Ld[:, 2]
                break
            else
                println(" GPFz has been converged")
                vmag_dis = round.(vmag_i, digits=4)
                println(vmag_dis)
            end
            vmag_t = deepcopy(vmag_i[iload])
            Vvio_t = length(find_below(vmag_t, Vsh, eps_v))
        end

        if Vvio_t != 0
            println(" There is short-term voltage violation")
            alarm = 0
            (Ld, uflag) = SV(Ld, vmag_t, alarm)
            # (G, Ld) = updateGV(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
            df_i = fr1(G, Ld)
            (G, Ld) = updateG(G, Ld, df_i)
            if uflag == 1
                println(" No available load for shedding")
                break
            end

        else
            println(" No short-term voltage violation")
            order = 0
            Ld = load_order(Ld, order)
            (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_s, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order)
            vmag_i_s = vmag_i
            if cg_s == 0
                println(" GPFp has been diverged,  go back to short-term")
                alarm = 1
                priority = vmag_t
                priority = sort(priority)
                # for i in 1:Vvio_t
                for i in 1:6

                    k = findall3(x -> x == priority[i], vmag_t)
                    k = k[1]
                    printstyled("There is a voltage vioaltion load $(Ld[k, 14]), {$(vmag_t[k])}:"; color=:blue, blink=true)

                    (Ld, uflag) = SV_alarm(Ld, vmag_t, k)
                    if uflag == 1
                        println(" No available load for shedding")
                        break
                    end
                    df_i = fr1(G, Ld)
                    (G, Ld) = updateG(G, Ld, df_i)
                    order = 0
                    Ld = load_order(Ld, order)
                    (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_s, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order)
                    #    #~~~~~~~~~~~~~~~~~~~~~~~~~~
                    #     vmag_t = deepcopy(vmag_i[iload])
                    #    #~~~~~~~~
                    if cg_s == 1
                        println(" GPFp has been converged, after $i more short-term shedding ")
                        break
                    end
                end
            else
                println(" GPFp has been converged,  check voltage violation")
                vmag_dis = round.(vmag_i, digits=4)
                println(vmag_dis)

                Vth = 0.9
                Vsh = Vth * ones(size(Ld, 1))
                vmag_s = deepcopy(vmag_i[iload])
                Vvio_s = length(find_below(vmag_s, Vsh, eps_v))
                if Vvio_s != 0
                    println(" There is long-term voltage violation")
                    alarm = 0
                    (Ld, uflag) = SV(Ld, vmag_s, alarm)
                    Ld = load_s(Ld)
                    # (G, Ld) = updateGV(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
                    df_i = fr1(G, Ld)
                    (G, Ld) = updateG(G, Ld, df_i)
                    if uflag == 1
                        println(" No available load for shedding")
                        break
                    end
                else
                    println(" No long-term voltage violation")
                    break
                end
            end
        end
    end

    df_i = fr1(G, Ld)
    println(" after uvls,   the frequency is :       $df_i pu")

    return B, L, T, G, Ld, Shu
end


function SV(Ld, vmag, alarm)
    Ld_ini = sum(Ld[:, 2])
    uflag = 0
    if alarm == 1
        k = findmin(vmag)
        k = k[2]
        Ld0v = Ld[k, 2]
        load_block = Ld[k, 5]
        min_load = Ld[k, 6]
        q_p = Ld[k, 3] / Ld[k, 2]
        if (Ld[k, 2] - 1 * load_block) >= min_load
            Ld[k, 2] = Ld[k, 2] - 1 * load_block
        else
            Ld[k, 2] = 0
        end
        Ld[k, 3] = q_p * Ld[k, 2]
        kmain = Ld[k, 14]
        LdI[kmain] = 2
        Ldsh_v[kmain] = Ldsh_v[kmain] + Ld0v - Ld[k, 2]
    else
        for k in 1:length(vmag)
            V = vmag[k]
            dVsh = Vth - V
            if dVsh > eps_v
                printstyled("There is a voltage vioaltion load $(Ld[k, 14]), {$V}:"; color=:blue, blink=true)
                println()
                Ld = shedv(Ld, dVsh, k)
            end
        end
    end
    if sum(Ld[:, 2]) == Ld_ini
        uflag = 1
    end

    return Ld, uflag
end


function SV_alarm(Ld, vmag, k)
    Ld_ini = sum(Ld[:, 2])
    uflag = 0


    Ld0v = Ld[k, 2]
    load_block = 2 * Ld[k, 5]
    min_load = Ld[k, 6]
    # if Ld[k, 2] == 0
    #     q_p = 0
    # else
    q_p = Ld[k, 3] / Ld[k, 2]
    # end

    if (Ld[k, 2] - 1 * load_block) >= min_load
        Ld[k, 2] = Ld[k, 2] - 1 * load_block
    else
        Ld[k, 2] = 0
    end
    Ld[k, 3] = q_p * Ld[k, 2]
    kmain = Ld[k, 14]
    LdI[kmain] = 2
    Ldsh_v[kmain] = Ldsh_v[kmain] + Ld0v - Ld[k, 2]

    if sum(Ld[:, 2]) == Ld_ini
        uflag = 1
    end

    return Ld, uflag
end




function shedv(Ld, dVsh, k)

    Ld0v = Ld[k, 2]
    load_block = 1 * Ld[k, 5]
    min_load = Ld[k, 6]
    q_p = Ld[k, 3] / Ld[k, 2]


    if (0.00 + eps_v < dVsh) && (dVsh <= 0.05 + eps_v)
        if (Ld[k, 2] - 1 * load_block) >= min_load
            Ld[k, 2] = (Ld[k, 2] - 3 * load_block)
        else
            Ld[k, 2] = 0
        end
    elseif (0.05 + eps_v < dVsh) && (dVsh <= 0.1 + eps_v)
        if (Ld[k, 2] - 2 * load_block) >= min_load
            Ld[k, 2] = Ld[k, 2] - 4 * load_block
        else
            Ld[k, 2] = 0
        end
    elseif (0.1 + eps_v < dVsh) && (dVsh <= 0.15 + eps_v)
        if (Ld[k, 2] - 3 * load_block) >= min_load
            Ld[k, 2] = Ld[k, 2] - 5 * load_block
        else
            Ld[k, 2] = 0
        end
    else
        Ld[k, 2] = 0
    end
    Ld[k, 3] = q_p * Ld[k, 2] # q

    kmain = Ld[k, 14]
    LdI[kmain] = 2
    Ldsh_v[kmain] = Ldsh_v[kmain] + Ld0v - Ld[k, 2]

    return Ld
end



function uls_simple(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
    eps_d = 0.001
    cg_i = 0
    # load_block = 0.01 * Ld[:, 2]
    load_block = 1 * Ld[:, 5]
    qp = Ld[:, 3] ./ Ld[:, 2]
    h = 0
    vmag_i = zeros(size(B, 1))
    iload = Ld[:, 1]
    Ldt0 = sum(Ld[:, 2])
    Ldini = deepcopy(Ld[:, 2])


    while true
        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_i, pim_i, df_i, itr_i, sij_i, vio_i) = multislack_acpf(y_bus_i, B, L, G, Ld, Shu)
        if cg_i == 1
            println(" power flow converged after $h run of ULS")
            break
        end
        # if sum(Ld[:, 2]) == 0
        if sum(Ld[:, 2]) < Ldt0 / 5
            kmain = Ld[:, 14]
            LdI[kmain] = 2 * ones(length(kmain))
            Ldsh_v[kmain] = Ldini[:, 1]
            Ld[:, 2] .= 0
            println(" blackout of the island")
            break
        end
        # ld1 = sum(Ld[:, 2])

        h += 1
        Ld[:, 2] = Ld[:, 2] - 1 * load_block
        vio_d = find_below(Ld[:, 2], Ld[:, 6], eps_d)
        Ld[vio_d, 2] = zeros(length(vio_d))
        Ld[:, 3] = qp .* Ld[:, 2]

        # ld2 = sum(Ld[:, 2])
        # LdI[1] = 2
        # Ldsh_v[1] = Ldsh_v[1] + ld1 - ld2

        kmain = Ld[:, 14]
        LdI[kmain] = 2 * ones(length(kmain))
        Ldsh_v[kmain] = Ldsh_v[kmain] + load_block


        df_i = fr1(G, Ld)
        (G, Ld) = updateG(G, Ld, df_i)
    end
    (B, L, T, G, Ld, Shu) = L1_r18(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)


    # vmag_t = deepcopy(vmag_i[iload])
    # (Ld, uflag) = SV(Ld, vmag_t, alarm)



    return B, L, T, G, Ld, Shu, REN, BAT, vmag_i
end




function uvls_r18(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)

    # B = B_i
    # L = L_i
    # T = T_i
    # G = G_i
    # Ld = Ld_i
    # Shu = Shu_i

    nBus = size(B, 1)
    vmag_i = zeros(size(B, 1))
    iload = Ld[:, 1]
    CBId = zeros(length(Ld[:, 1]))
    uflag = 0
    blk = 0
    (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_i, pim_i, df_i, itr_i, sij_i, vio_i, J, ipv, ipq, isl) = multislack_acpf_J(y_bus_i, B, L, G, Ld, Shu)
    # vmag_i_c=vmag_i
    if cg_i == 1
        println(" PF has been converged")
        (B, L, T, G, Ld, Shu) = L1_r18(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
    else
        println(" PF has been diverged")
        Ld_temp = deepcopy(Ld)
        G_temp = deepcopy(G)
        (B, L, T, G_temp, Ld_temp, Shu, vmag_i, blk) = uls_s(B, L, T, G_temp, Ld_temp, Shu, y_bus_i)
        if blk == 0
            (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_i, pim_i, df_i, itr_i, sij_i, vio_i, J, ipv, ipq, isl) = multislack_acpf_J(y_bus_i, B, L, G_temp, Ld_temp, Shu)
            CBI = zeros(size(B, 1))
            CBI1 = zeros(size(B, 1))
            CBI2 = zeros(size(B, 1))
            Pbus = zeros(size(B, 1))
            Jt = J'
            ires = setdiff(collect(1:size(Jt, 1)), isl)
            Jt = Jt[ires, :]
            J = Jt'
            ires = setdiff(collect(1:size(J, 1)), isl)
            J = J[ires, :]
            J11 = J[1:nBus-1, 1:nBus-1]
            J12 = J[1:nBus-1, nBus:end]
            J21 = J[nBus:end, 1:nBus-1]
            J22 = J[nBus:end, nBus:end]
            # cancel abs(V) which is already multiplied with elements of J12 and J22
            vmag = vmag_i[ipq]
            for i in 1:size(J12, 1)
                J12[i, :] = J12[i, :] ./ vmag
            end
            for i in 1:size(J22, 1)
                J22[i, :] = J22[i, :] ./ vmag
            end
            JR = J22 - J21 * inv(J11) * J12
            # JR=[0 -0.1092 -0.1236 0 0 0;376.99 0 0 0 0 0;0 -0.1938 -0.4229 -27.3172 0 27.3172;0 -7.3125 20.8391 -50 0 0 ; 0 -1.0372 -1.1738 0 -0.7143 0;0 -4.8404 -5.4777 0 26.9697 -30.3030]
            eig_JR = eigvals(JR)
            Landa = diagm(eig_JR)
            U = eigvecs(JR)
            W = inv(U')
            P = abs.(U .* W)
            imode = findmin(abs.(real.(eig_JR)))
            imode = imode[2]
            # iBUS = findmax(P[:, imode])
            # iBUS = iBUS[2]
            Pmain = P[:, imode]
            Pbus[ipq] = Pmain
            CBI1 = (Vth * ones(length(vmag_i)) - vmag_i) / maximum(abs.(Vth * ones(length(vmag_i)) - vmag_i))
            CBI2 = Pbus / maximum(Pbus)
            CBI = CBI1 + CBI2
            CBId = CBI[iload]

            vmag = vmag_i[iload]
            Vsh = Vth * ones(size(Ld, 1))
            Vvio_s = length(find_below(vmag, Vsh, eps_v))

            if Vvio_s != 0
                while true
                    (Ld, uflag) = SV2_r18(Ld, vmag, CBId)
                    if uflag == 1
                        println(" No available load for shedding, black out")
                        break
                    else
                        df_i = fr1(G, Ld)
                        (G, Ld) = updateG(G, Ld, df_i)
                        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_i, pim_i, df_i, itr_i, sij_i, vio_i, J, ipv, ipq, isl) = multislack_acpf_J(y_bus_i, B, L, G, Ld, Shu)
                        if cg_i == 1
                            println(" PF has been converged after uvls18")
                            vmag_dis = round.(vmag_i, digits=4)
                            println(vmag_dis)
                            (B, L, T, G, Ld, Shu) = L1_r18(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
                            break
                        end
                    end
                end
            else
                println(" No long-term voltage violation, initially")
                h = 1
                while true
                    (Ld, uflag) = SV4_r18(Ld, vmag, CBId)
                    if uflag == 1
                        println(" No available load for shedding, black out")
                        break
                    else
                        df_i = fr1(G, Ld)
                        (G, Ld) = updateG(G, Ld, df_i)
                        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_i, pim_i, df_i, itr_i, sij_i, vio_i, J, ipv, ipq, isl) = multislack_acpf_J(y_bus_i, B, L, G, Ld, Shu)
                        if cg_i == 1
                            println(" PF has been converged after $h uvls18")
                            (B, L, T, G, Ld, Shu) = L1_r18(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
                            break
                        end
                    end
                    h += 1
                end
            end
        end
    end
    return B, L, T, G, Ld, Shu
end



function uls_s(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
    eps_d = 0.001
    cg_i = 0
    load_block = 1 * Ld[:, 5]
    qp = Ld[:, 3] ./ Ld[:, 2]
    h = 0
    vmag_i = zeros(size(B, 1))
    blk = 0
    Ldt0 = sum(Ld[:, 2])
    Ldini = deepcopy(Ld[:, 2])

    while true
        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_i, pim_i, df_i, itr_i, sij_i, vio_i) = multislack_acpf(y_bus_i, B, L, G, Ld, Shu)
        if cg_i == 1
            println(" power flow converged after $h run of ULS")
            break
        end

        if sum(Ld[:, 2]) < Ldt0 / 5
            kmain = Ld[:, 14]
            LdI[kmain] = 2 * ones(length(kmain))
            Ldsh_v[kmain] = Ldsh_v[kmain] + Ldini[:, 1]
            Ld[:, 2] .= 0
            blk = 1
            println(" blackout of the island")
            break
        end

        h += 1
        Ld[:, 2] = Ld[:, 2] - 1 * load_block
        vio_d = find_below(Ld[:, 2], Ld[:, 6], eps_d)
        # Ld[vio_d, 2] = Ld[vio_d, 6]
        Ld[vio_d, 2] = zeros(length(vio_d))
        Ld[:, 3] = qp .* Ld[:, 2]

        df_i = fr1(G, Ld)
        (G, Ld) = updateG(G, Ld, df_i)

    end

    return B, L, T, G, Ld, Shu, REN, BAT, vmag_i, blk
end


function L1_r18(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i)
    iload = Ld[:, 1]
    while true
        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_i, pim_i, df_i, itr_i, sij_i, vio_i, J, ipv, ipq, isl) = multislack_acpf_J(y_bus_i, B, L, G, Ld, Shu)
        vmag = deepcopy(vmag_i[iload])
        Vsh = Vth * ones(size(Ld, 1))
        Vvio = length(find_below(vmag, Vsh, eps_v))
        if Vvio != 0
            println(" There is voltage violation")

            alarm = 0
            # (Ld, uflag) = SV(Ld, vmag, alarm)


            (Ld, uflag) = SV1_r18(Ld, vmag)
            if uflag == 1
                println(" No available load for shedding")
                break
            else
                df_i = fr1(G, Ld)
                (G, Ld) = updateG(G, Ld, df_i)
            end
        else
            println(" No voltage violation")
            break
        end
    end

    return B, L, T, G, Ld, Shu
end



function SV1_r18(Ld, vmag)
    Ld_ini = sum(Ld[:, 2])
    uflag = 0

    for k in 1:length(vmag)
        V = vmag[k]
        dVsh = Vth - V
        if dVsh > eps_v
            printstyled("There is a voltage vioaltion load $(Ld[k, 14]), {$V}:"; color=:blue, blink=true)
            println()
            Ld = shedv_r18(Ld, dVsh, k)
            # Ld = shedv(Ld, dVsh, k)
            # Ld = shedv(Ld, dVsh, k)
        end
    end

    if sum(Ld[:, 2]) == Ld_ini
        uflag = 1
    end

    return Ld, uflag
end


function SV2_r18(Ld, vmag, CBId)
    Ld_ini = sum(Ld[:, 2])
    uflag = 0
    priority = CBId
    priority = sort(priority, rev=true)
    i = 1
    dVsh = []

    for i in 1:length(CBId)

        k = findall3(x -> x == priority[i], CBId)
        V = vmag[k]
        dVsh = Vth - V[1]
        if dVsh > eps_v
            printstyled("There is a voltage vioaltion load $(Ld[k, 14]), {$V}:"; color=:blue, blink=true)
            println()
            Ld = shedv_r18(Ld, dVsh, k)
        else
            i += 1
        end
    end

    if sum(Ld[:, 2]) == Ld_ini
        uflag = 1
    end

    return Ld, uflag
end



function SV3_r18(Ld, vmag, CBId)
    Ld_ini = sum(Ld[:, 2])
    uflag = 0
    dVsh = []
    for k in 1:length(vmag)
        Ld = shedv_r18(Ld, dVsh, k)
    end

    if sum(Ld[:, 2]) == Ld_ini
        uflag = 1
    end

    return Ld, uflag
end



function SV4_r18(Ld, vmag, CBId)
    Ld_ini = sum(Ld[:, 2])
    uflag = 0
    # priority = vmag
    # priority = sort(priority)

    priority = CBId
    priority = sort(priority, rev=true)

    i = 1
    dVsh = []
    println(length(CBId))
    println(priority')

    for i in 1:length(CBId)
        # for i in 1:11
        k = findall3(x -> x == priority[i], CBId)
        # k = findall3(x -> x == priority[i], vmag)

        printstyled("Voltage ofload $(Ld[k, 14]) is $(vmag[k])}:"; color=:blue, blink=true)
        println()
        Ld = shedv_r18(Ld, dVsh, k)
    end

    if sum(Ld[:, 2]) == Ld_ini
        uflag = 1
    end

    return Ld, uflag
end



function shedv_r18(Ld, dVsh, k)
    Ld_ini = sum(Ld[:, 2])

    Ld0v = Ld[k, 2]
    load_block = Ld[k, 5]
    min_load = Ld[k, 6]
    q_p = Ld[k, 3] ./ Ld[k, 2]

    #**********
    if (Ld[k, 2] - 1 * load_block) >= min_load
        Ld[k, 2] = Ld[k, 2] - 1 * load_block
    end
    #**********

    Ld[k, 3] = q_p .* Ld[k, 2] # q


    kmain = Ld[k, 14]
    kmain = kmain[1]
    LdI[kmain] = 2
    Ldsh_v[kmain] = Ldsh_v[kmain] + sum(Ld0v - Ld[k, 2])


    return Ld
end


function updateGV(B, L, T, G, Ld, Shu, REN, BAT, y_bus_i, order)
    Gtemp = deepcopy(G)
    slk0 = deepcopy(G[:, 12])
    iload = Ld[:, 1]
    igen = G[:, 1]
    nBus = size(B, 1)
    p0 = zeros(nBus)
    dp = zeros(nBus, 3)

    p0[iload] = Ld[:, 2]
    dp[iload, :] = Ld[:, 7:9]

    while true
        # Step 1: update G for slack generators

        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_t, pim_i, df_i, itr_i, sij_i, vio_i) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order)

        if order != 0
            p_V = p0 .* (dp[:, 1] .* vmag_i .^ order .+ dp[:, 2] .* vmag_i .+ dp[:, 3])
        else
            p_V = p0
        end

        G[:, 2] = p_i[igen] + p_V[igen]

        # Step 2: check and fix violations
        vio_g1 = find_beyond(G[:, 2], G[:, 7], eps_g)
        vio_g2 = find_below(G[:, 2], G[:, 8], eps_g)
        vio_g = union(vio_g1, vio_g2)
        non_vio_g = setdiff(collect(1:size(G, 1)), vio_g)
        Gvio = length(vio_g)

        if Gvio == 0
            break
        else
            G[vio_g1, 2] = G[vio_g1, 7]
            G[vio_g2, 2] = G[vio_g2, 8]
            G[non_vio_g, 2] = Gtemp[non_vio_g, 2]
            G[vio_g, 12] = zeros(length(vio_g))
        end
    end

    G[:, 12] = deepcopy(slk0)

    return G, Ld
end



function pre_uvls(y_bus_i, B, L, G, Ld, Shu)

    # B = B_i
    # L = L_i
    # T = T_i
    # G = G_i
    # Ld = Ld_i
    # Shu = Shu_i


    slk0 = deepcopy(G[:, 12])
    iload = Ld[:, 1]
    igen = G[:, 1]
    nBus = size(B, 1)
    p0 = zeros(nBus)
    dp = zeros(nBus, 3)

    p0[iload] = Ld[:, 2]
    dp[iload, :] = Ld[:, 7:9]
    blkout = 0
    cg_t = 0
    Gvio = 0

    while true
        (G, Ld, df_i, blkout) = floop(G, Ld)

        # ONLY IF power flow converges, check G and f. based more accurately based on power flow
        Gtemp = deepcopy(G)
        Ld = load_t(Ld)
        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, cg_t, pim_i, df_i, itr_i, sij_i, vio_i) = multislack_acpf(y_bus_i, B, L, G, Ld, Shu) # GPF_p0
        if cg_t == 1
            p_V = p0 .* (dp[:, 1] .* vmag_i .^ 2 .+ dp[:, 2] .* vmag_i .+ dp[:, 3])
            G[:, 2] = p_i[igen] + p_V[igen]

            vio_g1 = find_beyond(G[:, 2], G[:, 7], eps_g)
            vio_g2 = find_below(G[:, 2], G[:, 8], eps_g)
            vio_g = union(vio_g1, vio_g2)
            non_vio_g = setdiff(collect(1:size(G, 1)), vio_g)
            Gvio = length(vio_g)

            df_i = fr1(G, Ld)
        end




        if (df_i <= df_max) && (df_i >= -df_min) && Gvio == 0 && cg_t == 1
            println(" considering g_limit, ploss, and v, final frequency change is :       $df_i pu")
            break
        elseif (df_i <= df_max) && (df_i >= -df_min) && Gvio == 0 && cg_t == 0
            println(" Without considering g_limit, ploss, and v, final frequency change is :       $df_i pu")
            break
        elseif (df_i <= df_max) && (df_i >= -df_min) && Gvio != 0
            G[vio_g1, 2] = G[vio_g1, 7]
            G[vio_g2, 2] = G[vio_g2, 8]
            G[non_vio_g, 2] = Gtemp[non_vio_g, 2]
            G[vio_g, 12] = zeros(length(vio_g))
        end
    end

    G[:, 12] = deepcopy(slk0)
    return G, Ld, blkout
end









function floop(G, Ld)
    #  B = B_i
    #     L = L_i
    #     T = T_i
    #     G = G_i
    #     Ld = Ld_i
    #     Shu = Shu_i


    blkout = 0
    df_i = fr1(G, Ld)
    println(" by no action.................  the frequency is :       $df_i pu")

    if df_i > df_max + df_eps
        # (G, Ld, df_i) = ofgt(G, Ld, df_i)
        (G, Ld, df_i) = ofgt(G, Ld, df_i)
        println(" by ofgt (no network effect),   the frequency is :       $df_i pu")
    elseif df_i < -(df_min + df_eps)
        (G, Ld, df_i, blkout) = ufls(G, Ld, df_i)
        if blkout == 1
            println("blackout of island, due to frequency issue")
        else
            println(" by ufls (no network effect),   the frequency is :       $df_i pu")
        end
    else
        println(" no need to ofgt or ufls  ")
    end

    return G, Ld, df_i, blkout
end



function ofgt(G, Ld, df_i)

    priority = G[:, 14]
    priority = sort(priority)
    i = 1
    bg = G[:, 15]

    while df_i > df_max + df_eps

        if df_i < -(df_min + df_eps)
            println("ufls is also done in ofgt")
            (G, Ld, df_i, blkout) = ufls(G, Ld, df_i)
            break
        end


        k = findall3(x -> x == priority[i], G[:, 14])
        if bg[k] == 0
            i += 1
            G[k, 12] = 0
        else
            G[k, 2] = G[k, 2] - G[k, 2] / sum(bg[k])
            G[k, 4] = G[k, 4] - G[k, 4] / sum(bg[k])
            G[k, 5] = G[k, 5] - G[k, 5] / sum(bg[k])
            G[k, 7] = G[k, 7] - G[k, 7] / sum(bg[k])
            G[k, 8] = G[k, 8] - G[k, 8] / sum(bg[k])

            bg[k] = bg[k] - [1]

            kmain = G[k, 16]
            Gt[kmain] = Gt[kmain] + G[k, 2] / sum(bg[k])
            df_i = fr1(G, Ld)
            println(df_i)
            (G, Ld) = updateG(G, Ld, df_i)
            df_i = fr1(G, Ld)
            println(df_i)
        end
    end


    return G, Ld, df_i
end



function ufls(G, Ld, df_i)

    priority = Ld[:, 4]
    priority = sort(priority)

    i = 1
    blkout = 0

    while df_i < -(df_min + df_eps)

        Load = sum(Ld[:, 2])
        minLd = sum(Ld[:, 6])

        if (Load <= minLd)
            println("some  loads need to be totally shed")
            (G, Ld, df_i, blkout) = ufls_special(G, Ld, df_i)
            break
        else
            k = findall3(x -> x == priority[i], Ld[:, 4])
            load_block = Ld[k, 5]
            min_load = Ld[k, 6]
            if (Ld[k, 2] - 1 * load_block) < min_load
                i += 1
            else
                load_block = Ld[k, 5]
                q_p = Ld[k, 3] ./ Ld[k, 2]
                q_p = q_p[1]
                Ld[k, 2] = Ld[k, 2] - load_block
                Ld[k, 3] = q_p * Ld[k, 2] # q

                kmain = Ld[k, 14]
                LdI[kmain] = ones(length(kmain))
                Ldsh_f[kmain] = Ldsh_f[kmain] + load_block
                printstyled("Load $(Ld[k, 14]) is shed, {$load_block}:"; color=:blue, blink=true)
                println()
                df_i = fr1(G, Ld)
                println(df_i)
                (G, Ld) = updateG(G, Ld, df_i)
                df_i = fr1(G, Ld)
                println(df_i)
            end

        end
    end

    return G, Ld, df_i, blkout
end



function ufls_special(G, Ld, df_i)

    priority = Ld[:, 4]
    priority = sort(priority)
    Ld0 = Ld[:, 2]

    i = 1
    blkout = 0

    while df_i < -(df_min + df_eps)
        k = findall3(x -> x == priority[i], Ld[:, 4])
        if Ld[k, 2] == 0
            i += 1
        else
            Ld[k, 2] .= 0
            Ld[k, 3] .= 0
            kmain = Ld[k, 14]
            LdI[kmain] = ones(length(kmain))
            Ldsh_f[kmain] = Ldsh_f[kmain] + Ld0[kmain]

            df_i = fr1(G, Ld)

            (G, Ld) = updateG(G, Ld, df_i)

        end
    end


    return G, Ld, df_i, blkout
end



function updateG(G, Ld, df_i)
    Gtemp = deepcopy(G)
    slk0 = deepcopy(G[:, 12])

    while true
        # Step 1: update G for slack generators
        slk = findall3(x -> x == 1, G[:, 12])
        Ri = G[slk, 13]
        beti = ones(length(Ri)) ./ Ri
        G[slk, 2] = G[slk, 2] - df_i * beti

        # Step 2: check and fix violations
        vio_g1 = find_beyond(G[:, 2], G[:, 7], eps_g)
        vio_g2 = find_below(G[:, 2], G[:, 8], eps_g)
        vio_g = union(vio_g1, vio_g2)
        non_vio_g = setdiff(collect(1:size(G, 1)), vio_g)
        Gvio = length(vio_g)

        if Gvio == 0
            break
        else
            G[vio_g1, 2] = G[vio_g1, 7]
            G[vio_g2, 2] = G[vio_g2, 8]
            G[non_vio_g, 2] = Gtemp[non_vio_g, 2]
            G[vio_g, 12] = zeros(length(vio_g))
        end
    end

    G[:, 12] = deepcopy(slk0)

    return G, Ld
end


function fr1(G, Ld)

    Pd = sum(Ld[:, 2])
    Pg = sum(G[:, 2])
    ploss = 0.0

    slk = findall3(x -> x == 1, G[:, 12])
    Ri = G[slk, 13]
    beti = ones(length(Ri)) ./ Ri

    Pmm = Pd + ploss - Pg

    # this imbalance is compensated by the 'slack' units not all
    df_i = -Pmm / sum(beti)

    return df_i
end




function load_order(Ld, order)

    V0 = Ld[:, 13]
    if order != 0
        Ld[:, 7] = ones(size(Ld, 1)) ./ (V0 .^ order)    # d2p=1/(V0^2)
        Ld[:, 10] = ones(size(Ld, 1)) ./ (V0 .^ order)   # d2q=1/(V0^2)
        Ld[:, 8:9] = zeros(size(Ld, 1), 2)           # d1p=d0p=0
        Ld[:, 11:12] = zeros(size(Ld, 1), 2)         # d1q=d0q=0
    else
        Ld[:, 7:8] = zeros(size(Ld, 1), 2)      # d2p=d1p=0
        Ld[:, 9] = ones(size(Ld, 1))    # d0p=1
        Ld[:, 10:11] = zeros(size(Ld, 1), 2)      # d2q=d1q=0
        Ld[:, 12] = ones(size(Ld, 1))     # d0q=1
    end


    return Ld
end

function load_t(Ld)

    V0 = Ld[:, 13]

    Ld[:, 7] = ones(size(Ld, 1)) ./ (V0 .^ 2)    # d2p=1/(V0^2)
    Ld[:, 10] = ones(size(Ld, 1)) ./ (V0 .^ 2)   # d2q=1/(V0^2)
    Ld[:, 8:9] = zeros(size(Ld, 1), 2)           # d1p=d0p=0
    Ld[:, 11:12] = zeros(size(Ld, 1), 2)         # d1q=d0q=0

    return Ld
end




function load_s(Ld)

    Ld[:, 7:8] = zeros(size(Ld, 1), 2)      # d2p=d1p=0
    Ld[:, 9] = ones(size(Ld, 1), 1)      # d0p=1
    Ld[:, 10:11] = zeros(size(Ld, 1), 2)      # d2q=d1q=0
    Ld[:, 12] = ones(size(Ld, 1), 1)      # d0q=1

    return Ld
end


function multislack_acpf_order_pmaxqmax(B, L, T, G, Ld, Shu, REN, BAT, Lout)
    # B, L, T, G, Ld, Shu,REN,BAT, Lout=NSysData["Bus"], NSysData["Line"], NSysData["Trans"], NSysData["Gen"], NSysData["Load"], NSysData["Shunt"],NSysData["RES"],NSysData["ESS"],[]

    eps_g=0.01
    ires = setdiff(collect(1:size(L, 1)), Lout)
    L = L[ires, :]

    nBus = size(B, 1)
    #.............................................................................
    # A matrix
    #.............................................................................
    A = spzeros(Int64, nBus, nBus)
    A = A_net(nBus, L, T)
    for i in 1:size(A, 1)
        for j in 1:size(A, 1)
            if A[i, j] != 0
                A[i, j] = 1
            end
        end
    end
    # A = convert(Matrix{Int64}, A)
    ctd = find_ones(A)

    iload = Ld[:, 1]
    ign = G[:, 1]
    order = 0
    Ld = load_order(Ld, order)
    y_bus_i = Y_bus1(nBus, L, T, ctd)
    Gtemp = G
    exit = 0
    Gvio = 1

    (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order, ctd)

    while exit == 0
        G[:, 2] = p_i[ign] + p_V[ign]
        vio_g1 = find_beyond(G[:, 2], G[:, 7], eps_g)
        vio_g2 = find_below(G[:, 2], G[:, 8], eps_g)
        vio_g = union(vio_g1, vio_g2)
        non_vio_g = setdiff(collect(1:size(G, 1)), vio_g)
        Gvio = length(vio_g)

        if cg_t == 1 && Gvio == 0
            exit = 1
            println("generataor limits are satisfied")
            break
        else
            G[vio_g1, 2] = G[vio_g1, 7]
            G[vio_g2, 2] = G[vio_g2, 8]
            G[non_vio_g, 2] = Gtemp[non_vio_g, 2]
            G[vio_g, 12] = zeros(length(vio_g))
        end
        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order, ctd)
    end

    PGGG=deepcopy(G[:,2])
    QGGG=deepcopy(G[:,3])
    # return vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i
    return vmag_i, vang_i, p_i, q_i, cg_t,Iij,percentIij, PGGG, QGGG,JacobOut
end




function initiate_IND_gfm(B, L, G, Ld, Shu)
    nBus = size(B, 1)
    iload = Ld[:, 1]
    ign = G[:, 1]

    order = 0

    Ld = load_order(Ld, order)
    y_bus_i = Y_bus1(nBus, L, T)

    for i in 1:1
        (vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_t, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i,Iij,percentIij,JacobOut) = multislack_acpf_order(y_bus_i, B, L, G, Ld, Shu, REN, BAT, order)

        Ld[:, 2] = p_V[iload]
        Ld[:, 3] = q_V[iload]
        G[:, 2] = p_i[ign] + p_V[ign]
        G[:, 3] = q_i[ign] + q_V[ign]

        Ld[:, 13] = vmag_i[iload]
        Ld[:, 21] = Gs[iload]
        Ld[:, 22] = Bs[iload]

        print_multislack_acpf(vmag_i, vang_i, p_i, q_i, p_V, q_V, cg_t, df_i, sm_i, Gs, Bs, itr_i)
    end

    return Ld, G, y_bus_i
end



function modif_Y(B, L, T, Ld, G, y_bus)
    nBus = size(B, 1)
    y_ld = zeros(ComplexF64, nBus)
    y_sg = zeros(ComplexF64, nBus)
    iload = Ld[:, 1]
    ign = G[:, 1]
    ign = convert(Vector{Int64}, ign)


    p0 = Ld[:, 2]
    q0 = Ld[:, 3]
    V0 = Ld[:, 13]
    Sb = 100
    Sbm = G[:, 24]
    X1d = G[:, 18] .* (Sb ./ Sbm)
    Vb = 345
    Zb = (Vb^2) / Sb
    Zf = 0.1 + im * 0# Ohm
    Zf /= Zb

    y_ld[iload] = (p0 .- (im) * q0) ./ (V0 .^ 2)
    y_sg[ign] = 1 ./ (im * X1d)
    y_f = 1 ./ Zf
    idx_f = 4


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    isl = findall(x -> x == 3, B[:, 2])
    if isl == []
        (xxx, isG) = findmax(G[:, 7])
        isl = G[isG, 1]
    else
        isl = isl[1]
    end

    idx = findall3(x -> x == isl, ign)
    isl = nBus + idx[1]

    B = ones(Int64, nBus + length(ign), 2)
    B[:, 1] = collect(1:nBus+length(ign))
    B[nBus+1:end, 2] .= 2
    B[isl, 2] = 3

    nL = size(L, 1)
    L = vcat(L, zeros(length(ign), size(L, 2)))
    L[nL+1:end, 1] = ign
    L[nL+1:end, 2] = collect(nBus+1:nBus+length(ign))
    L[nL+1:end, 4] = -1 ./ abs.(X1d)

    G[:, 1] = collect(nBus+1:nBus+length(ign))

    nBus_new = size(B, 1)
    y_bus = Y_bus1(nBus_new, L, T)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Add (load, synchronous machine, fault) impedances to ybus matrix
    for i in 1:nBus
        #~~~ for isgg
        # y_bus[i, i] += y_ld[i] + y_sg[i]
        # y_bus[i, i] +=  y_sg[i]
        #~~~ for ind
        # y_bus[i, i] += y_ld[i]
    end
    y_bus[idx_f, idx_f] += y_f

    # ires = []
    # Ld = Ld[ires, :]

    return B, L, Ld, G, y_bus
end



















function find_below(x, x_min, eps)
    i_vio = Int64[]
    for i in 1:length(x)
        if x[i] <= x_min[i] - eps
            push!(i_vio, i)
        end
    end
    return i_vio

end


function find_beyond(x, x_max, eps)
    i_vio = Int64[]
    for i in 1:length(x)
        if x[i] >= x_max[i] + eps
            push!(i_vio, i)
        end
    end
    return i_vio

end


#----------------------------------------------------------------------------------------------
# [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
# [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
# [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
# [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
# [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
# [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
# [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
#----------------------------------------------------------------------------------------------


function schedule_powerV(nBus, PgenBUS, QgenBUS, QshuntBUS)
    p_sch_pf = zeros(nBus)
    q_sch_pf = zeros(nBus)

    p_sch_pf = PgenBUS
    q_sch_pf = QgenBUS + QshuntBUS

    return p_sch_pf, q_sch_pf
end



function spc_slack(ki, G)

    is_slack = findall(x -> x == 1, G[:, 12])

    if length(is_slack) != 0
        beti = zeros(length(is_slack))
        Ri = G[is_slack, 13]
        beti = ones(length(Ri)) ./ Ri
        indices_slacks = G[is_slack, 1]
        ki[indices_slacks] = beti
    else
        beti = 100000
    end


    ki = convert(Array{Float64}, ki)
    return ki, beti
end

function y_motor(Ld_ind, df, fm)

    R1 = Ld_ind[:, 1]
    X1 = Ld_ind[:, 2]
    Xm = Ld_ind[:, 3]
    R2 = Ld_ind[:, 4]
    X2 = Ld_ind[:, 5]

    A = R2 .* Xm .^ 2
    B = (Xm + X2) .^ 2
    C = (R2 .^ 2) .* Xm
    D = X2 .* Xm .* (Xm + X2)

    Rs = R1 + (A * (1 + df) .* (1 + df .- fm)) ./ (R2 .^ 2 + B .* ((1 + df .- fm) .^ 2))
    Xs = (1 + df) * (X1 + (C + D .* ((1 + df .- fm) .^ 2)) ./ (R2 .^ 2 + B .* (1 + df .- fm) .^ 2))

    Gm = (Rs ./ (Rs .^ 2 + Xs .^ 2))
    Bm = (Xs ./ (Rs .^ 2 + Xs .^ 2))

    return Gm, Bm
end

function torque_motor(v_mag, df, fm, Ld_ind)

    R1 = Ld_ind[:, 1]
    X1 = Ld_ind[:, 2]
    Xm = Ld_ind[:, 3]
    R2 = Ld_ind[:, 4]
    X2 = Ld_ind[:, 5]
    am = Ld_ind[:, 6]

    A = R2 .* Xm .^ 2
    B = (Xm + X2) .^ 2
    C = (R2 .^ 2) .* Xm
    D = X2 .* Xm .* (Xm + X2)

    #    A1 = A * (1 + df) ^ 2
    #    B1 = B * (1 + df) ^ 2
    #    C1 = C * (1 + df) ^ 1
    #    D1 = D * (1 + df) ^ 3

    # sm=(1+ df .-fm) / (1 + df)
    # Rs = R1 + (A1 .* sm) ./ (R2 .^ 2 + B1 .* (sm .^ 2))
    # Xs = X1 + (C1 + D1 .* (sm .^ 2)) ./ (R2 .^ 2 + B1 .* sm .^ 2)


    Rs = R1 + (A * (1 + df) .* (1 + df .- fm)) ./ (R2 .^ 2 + B .* ((1 + df .- fm) .^ 2))
    Xs = (1 + df) * (X1 + (C + D .* ((1 + df .- fm) .^ 2)) ./ (R2 .^ 2 + B .* (1 + df .- fm) .^ 2))


    Gs = (Rs ./ (Rs .^ 2 + Xs .^ 2))
    Bs = (Xs ./ (Rs .^ 2 + Xs .^ 2))

    # g1 = (sm .* R2 .* Xm .^ 2) ./ (R2 .^ 2 + B .* sm .^ 2)
    g1 = (1 + df) * (A .* (1 + df .- fm)) ./ (R2 .^ 2 + B .* (1 + df .- fm) .^ 2)
    g2 = (Gs .^ 2 + Bs .^ 2) .* (v_mag .^ 2)

    Te = g1 .* g2
    km = 1.00
    Tm0 = 1
    Tm1 = 0
    # Tm2=0.7888+1.23
    Tm2 = 0
    Tm = km .* (Tm0 .+ Tm1 .* fm .+ Tm2 .* fm .^ 2)
    delta_tq = Te - Tm

    # ~ diff___fm
    diff_Rfm = -(1 + df) * A .* (R2 .^ 2 - B .* (1 + df .- fm) .^ 2) ./ ((R2 .^ 2 + B .* (1 + df .- fm) .^ 2) .^ 2)
    diff_Xfm = -(1 + df) * (2 * D .* (1 + df .- fm) .* R2 .^ 2 - 2 * B .* C .* (1 + df .- fm)) ./ ((R2 .^ 2 + B .* (1 + df .- fm) .^ 2) .^ 2)
    diff_Gfm = (diff_Rfm .* (Rs .^ 2 + Xs .^ 2) - Rs .* (2 * Rs .* diff_Rfm + 2 * Xs .* diff_Xfm)) ./ ((Rs .^ 2 + Xs .^ 2) .^ 2)
    diff_Bfm = (diff_Xfm .* (Rs .^ 2 + Xs .^ 2) - Xs .* (2 * Rs .* diff_Rfm + 2 * Xs .* diff_Xfm)) ./ ((Rs .^ 2 + Xs .^ 2) .^ 2)
    diff_g1fm = -(1 + df) * A .* (R2 .^ 2 - B .* (1 + df .- fm) .^ 2) ./ ((R2 .^ 2 + B .* (1 + df .- fm) .^ 2) .^ 2)
    diff_g2fm = (v_mag .^ 2) .* (2 .* Gs .* diff_Gfm + 2 .* Bs .* diff_Bfm)
    diff_Tefm = g1 .* diff_g2fm + g2 .* diff_g1fm
    diff_Tmfm = km .* (Tm1 .+ 2 .* Tm2 .* fm)
    diff_Tfm = diff_Tefm - diff_Tmfm

    # ~ diff___df
    diff_Rf = (A .* R2 .^ 2 .* (2 + 2 * df .- fm) - A .* B .* fm .* (1 + df .- fm) .^ 2) ./ ((R2 .^ 2 + B .* (1 + df .- fm) .^ 2) .^ 2)
    diff_Xf = X1 + ((2 * D .* R2 .^ 2 - 2 * B .* C) .* (1 + df .- fm)) ./ ((R2 .^ 2 + B .* (1 + df .- fm) .^ 2) .^ 2)
    diff_Gf = (diff_Rf .* (Rs .^ 2 + Xs .^ 2) - Rs .* (2 * Rs .* diff_Rf + 2 * Xs .* diff_Xf)) ./ ((Rs .^ 2 + Xs .^ 2) .^ 2)
    diff_Bf = (diff_Xf .* (Rs .^ 2 + Xs .^ 2) - Xs .* (2 * Rs .* diff_Rf + 2 * Xs .* diff_Xf)) ./ ((Rs .^ 2 + Xs .^ 2) .^ 2)
    diff_g1f = ((A .* R2 .^ 2) .* (2 + 2 * df .- fm) - (A .* B .* fm) .* (1 + df .- fm) .^ 2) ./ ((R2 .^ 2 + B .* (1 + df .- fm) .^ 2) .^ 2)
    diff_g2f = (v_mag .^ 2) .* (2 .* Gs .* diff_Gf + 2 .* Bs .* diff_Bf)
    diff_Tef = g1 .* diff_g2f + g2 .* diff_g1f
    diff_Tmf = 0.00
    diff_Tf = diff_Tef .- diff_Tmf

    # ~ diff___V
    diff_TV = 2 .* g1 .* (Gs .^ 2 + Bs .^ 2) .* (v_mag)

    return delta_tq, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm
end


function consumed_power(Ld, v_mag, order)

    p0 = Ld[:, 2]
    q0 = Ld[:, 3]
    dp = Ld[:, 7:9]
    dq = Ld[:, 10:12]
    p_V = zeros(size(Ld, 1))
    q_V = zeros(size(Ld, 1))

    if order != 0
        p_V = p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3])
        q_V = q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3])
    else
        p_V = p0
        q_V = q0
    end

    return p_V, q_V
end




function jacobian1(nBus, isl, ipq, v_angle, v_mag, y_bus, p_inj, q_inj, ki)
    J11 = zeros(nBus, nBus)
    J12 = zeros(nBus, nBus)
    J21 = zeros(nBus, nBus)
    J22 = zeros(nBus, nBus)

    for i in 1:nBus
        for j in 1:nBus
            if i == j                                             # Diagonal elements of Jacobian sub-matrices

                J11[j, j] = -q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j])
                J21[j, j] = p_inj[j, 1] - v_mag[j, 1]^2 * real(y_bus[j, j])
                J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j])
                J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j])

            else                                                  # Off-diagonal elements of Jacobian sub-matrices

                theta_diff_ij = v_angle[i, 1] - v_angle[j, 1]
                real_val = real(y_bus[i, j])
                imag_val = imag(y_bus[i, j])
                sin_val = sin(theta_diff_ij)
                cos_val = cos(theta_diff_ij)
                partial_p_inj = +1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * sin_val - imag_val * cos_val)
                partial_q_inj = -1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * cos_val + imag_val * sin_val)
                J11[i, j] = partial_p_inj
                J21[i, j] = partial_q_inj
                J12[i, j] = -partial_q_inj
                J22[i, j] = partial_p_inj
            end
        end
    end


    J11_tmp = J11[1:end, 1:end]
    J12_tmp = J12[1:end, ipq]
    J21_tmp = J21[ipq, 1:end]
    J22_tmp = J22[ipq, ipq]
    J = [J11_tmp J12_tmp; J21_tmp J22_tmp]
    J[:, isl] = -ki

    return J
end




function mismatch_vector1(nBus, isl, ipv, ipq, v_angle, v_mag, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki)
    p_inj = zeros(nBus)
    q_inj = zeros(nBus)
    p_calc = zeros(nBus)
    q_calc = zeros(nBus)

    p_plus = -ki * pim_pf

    for i in 1:nBus  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses

        theta_diff = v_angle[i, 1] .- v_angle
        real_val = real(y_bus[i, :])
        imag_val = imag(y_bus[i, :])
        sin_val = sin.(theta_diff)
        cos_val = cos.(theta_diff)

        p_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* cos_val + (imag_val) .* sin_val), dims=1)
        q_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* sin_val - (imag_val) .* cos_val), dims=1)


        p_inj[i] = p_inj_1[1]
        q_inj[i] = q_inj_1[1]


        p_calc[i] = p_inj[i] + p_plus[i]
        q_calc[i] = q_inj[i]
    end


    delta_p = p_sch_pf .- p_calc    # Mismatch active power
    delta_q = q_sch_pf .- q_calc    # Mismatch reactive power

    if isempty(ipv)                                                    # No PV bus in the network
        delta_mismatch_inj = vcat(delta_p, delta_q)
    else                                                                         # Existence of PV buses in a system
        # Remove the corresponding reactive power mismatch entry from the delta_Y vector

        delta_p_tmp = delta_p[1:end]
        delta_q_tmp = delta_q[ipq]
        delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp)
    end
    return delta_mismatch_inj, p_inj, q_inj
end




function jacobianV(nBus, isl, ipq, v_angle, v_mag, y_bus, p_inj, q_inj, ki, p0, q0, dp, dq;flag=0)
    J11 = zeros(nBus, nBus)
    J12 = zeros(nBus, nBus)
    J21 = zeros(nBus, nBus)
    J22 = zeros(nBus, nBus)

    for i in 1:nBus
        for j in 1:nBus
            if i == j                                             # Diagonal elements of Jacobian sub-matrices

                J11[j, j] = -q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j])
                J21[j, j] = p_inj[j, 1] - v_mag[j, 1]^2 * real(y_bus[j, j])
                J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j]) + p0[j] * (2 * dp[j, 1] * v_mag[j, 1] + dp[j, 2])
                J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) + q0[j] * (2 * dq[j, 1] * v_mag[j, 1] + dq[j, 2])

            else                                                  # Off-diagonal elements of Jacobian sub-matrices

                theta_diff_ij = v_angle[i, 1] - v_angle[j, 1]
                real_val = real(y_bus[i, j])
                imag_val = imag(y_bus[i, j])
                sin_val = sin(theta_diff_ij)
                cos_val = cos(theta_diff_ij)
                partial_p_inj = +1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * sin_val - imag_val * cos_val)
                partial_q_inj = -1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * cos_val + imag_val * sin_val)
                J11[i, j] = partial_p_inj
                J21[i, j] = partial_q_inj
                J12[i, j] = -partial_q_inj
                J22[i, j] = partial_p_inj
            end
        end
    end


    J11_tmp = J11[1:end, 1:end]
    J12_tmp = J12[1:end, ipq]
    J21_tmp = J21[ipq, 1:end]
    J22_tmp = J22[ipq, ipq]
    J = [J11_tmp J12_tmp; J21_tmp J22_tmp]
    J[:, isl] = -ki

    return J
end


function jacobianV_order(nBus, isl, ipq, v_angle, v_mag, y_bus, p_inj, q_inj, ki, p0, q0, dp, dq, order)
    J11 = zeros(nBus, nBus)
    J12 = zeros(nBus, nBus)
    J21 = zeros(nBus, nBus)
    J22 = zeros(nBus, nBus)

    for i in 1:nBus
        for j in 1:nBus
            if i == j                                             # Diagonal elements of Jacobian sub-matrices

                J11[j, j] = -q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j])
                J21[j, j] = p_inj[j, 1] - v_mag[j, 1]^2 * real(y_bus[j, j])
                if order != 0
                    J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j]) + p0[j] * (order * dp[j, 1] * v_mag[j, 1] .^ (order - 1))
                    J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) + q0[j] * (order * dq[j, 1] * v_mag[j, 1] .^ (order - 1))
                else
                    J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j])
                    J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j])
                end

            else                                                  # Off-diagonal elements of Jacobian sub-matrices

                theta_diff_ij = v_angle[i, 1] - v_angle[j, 1]
                real_val = real(y_bus[i, j])
                imag_val = imag(y_bus[i, j])
                sin_val = sin(theta_diff_ij)
                cos_val = cos(theta_diff_ij)
                partial_p_inj = +1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * sin_val - imag_val * cos_val)
                partial_q_inj = -1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * cos_val + imag_val * sin_val)
                J11[i, j] = partial_p_inj
                J21[i, j] = partial_q_inj
                J12[i, j] = -partial_q_inj
                J22[i, j] = partial_p_inj
            end
        end
    end


    J11_tmp = J11[1:end, 1:end]
    J12_tmp = J12[1:end, ipq]
    J21_tmp = J21[ipq, 1:end]
    J22_tmp = J22[ipq, ipq]
    J = [J11_tmp J12_tmp; J21_tmp J22_tmp]
    J[:, isl] = -ki

    return J
end


function jacobianV_orderIND(nBus, isl, ipq, idu, v_angle, v_mag, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, p0, q0, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind, ctd;flag=0)
    # nBus, isl, ipq, idu, v_angle, v_mag, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, p0, q0, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind, ctd=nBus, isl, ipq, idu, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, ctd

    J11 = spzeros(nBus, nBus)
    J12 = spzeros(nBus, nBus)
    J21 = spzeros(nBus, nBus)
    J22 = spzeros(nBus, nBus)

    J13 = spzeros(nBus, nBus)
    J23 = spzeros(nBus, nBus)
    J32 = spzeros(nBus, nBus)
    J33 = spzeros(nBus, nBus)

    for i in 1:nBus
        ctd_i = ctd[i]
        for j in ctd_i                                      # Off-diagonal elements of Jacobian sub-matrices
            theta_diff_ij = v_angle[i, 1] - v_angle[j, 1]
            real_val = real(y_bus[i, j])
            imag_val = imag(y_bus[i, j])
            sin_val = sin(theta_diff_ij)
            cos_val = cos(theta_diff_ij)
            partial_p_inj = +1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * sin_val - imag_val * cos_val)
            partial_q_inj = -1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * cos_val + imag_val * sin_val)
            J11[i, j] = partial_p_inj
            J21[i, j] = partial_q_inj
        end
    end
    J12 = -1 * J21
    J22 = +1 * J11

    for j in 1:nBus                                          # Diagonal elements of Jacobian sub-matrices
        if order != 0
            J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j]) + (1 - ind[j]) * p0[j] * (order * dp[j, 1] * v_mag[j, 1] .^ (order - 1)) + ind[j] * p0[j] * (Gs[j] ./ Gs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
            J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) + (1 - ind[j]) * q0[j] * (order * dq[j, 1] * v_mag[j, 1] .^ (order - 1)) + ind[j] * q0[j] * (Bs[j] ./ Bs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
        else
            J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j]) + ind[j] * p0[j] * (Gs[j] ./ Gs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
            J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) + ind[j] * q0[j] * (Bs[j] ./ Bs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
        end
        J11[j, j] = -q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j])
        J21[j, j] = p_inj[j, 1] - v_mag[j, 1]^2 * real(y_bus[j, j])

        J13[j, j] = ind[j] * p0[j] * (diff_Gfm[j] ./ Gs0[j]) .* (v_mag[j, 1] ./ v_mag0[j]) .^ 2
        J23[j, j] = ind[j] * q0[j] * (diff_Bfm[j] ./ Bs0[j]) .* (v_mag[j, 1] ./ v_mag0[j]) .^ 2
        J32[j, j] = -diff_TV[j]
        J33[j, j] = -diff_Tfm[j]
    end

    if flag==1
        J = [J11 J12; J21 J22]
        # J[:,isl].=1
        ind_AR=findall(x->x!=0,ki[:,1])
        val_AR=filter(x->x!=0,ki[:,1])
        J[ind_AR, isl] = val_AR
    else
        J11_tmp = J11[1:end, 1:end]
        J12_tmp = J12[1:end, ipq]
        J21_tmp = J21[ipq, 1:end]
        J22_tmp = J22[ipq, ipq]

        J13_tmp = J13[1:end, idu]
        J23_tmp = J23[ipq, idu]
        J31_tmp = zeros(length(idu), nBus)
        J32_tmp = J32[idu, ipq]
        J33_tmp = J33[idu, idu]

        J = [J11_tmp J12_tmp J13_tmp; J21_tmp J22_tmp J23_tmp; J31_tmp J32_tmp J33_tmp]
        J[:, isl] = ki

        diff_pf = (1 .- ind) .* p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3]) .* kip + ind .* p0 .* (diff_Gf ./ Gs0) .* (v_mag ./ v_mag0) .^ 2
        diff_qf = (1 .- ind) .* q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3]) .* kiq + ind .* q0 .* (diff_Bf ./ Bs0) .* (v_mag ./ v_mag0) .^ 2

        J[1:nBus, isl] = J[1:nBus, isl] + diff_pf
        J[nBus+1:nBus+length(ipq), isl] = J[nBus+1:nBus+length(ipq), isl] + diff_qf[ipq]
        J[nBus+length(ipq)+1:end, isl] = J[nBus+length(ipq)+1:end, isl] - diff_Tf[idu]
    end

    return J
end



function jacobianV_isGG(nBus, isl, ign, ipq, idu, v_angle, v_mag, y_bus, p_inj, q_inj, Pnor, Qnor, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, p0, q0, e1m, e1a, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind, sg, X1d)
    # v_angle=v_angle_pf
    # v_mag= v_mag_pf
    # p0=PloadBUS
    # q0=QloadBUS
    # v_mag0=Vi0

    J11 = zeros(nBus, nBus)
    J12 = zeros(nBus, nBus)
    J21 = zeros(nBus, nBus)
    J22 = zeros(nBus, nBus)

    J13 = zeros(nBus, nBus)
    J23 = zeros(nBus, nBus)
    J31 = zeros(nBus, nBus)
    J32 = zeros(nBus, nBus)
    J33 = zeros(nBus, nBus)

    J41 = zeros(nBus, nBus)
    J42 = zeros(nBus, nBus)
    J43 = zeros(nBus, nBus)
    J14 = zeros(nBus, nBus)
    J24 = zeros(nBus, nBus)
    J34 = zeros(nBus, nBus)
    J44 = zeros(nBus, nBus)

    for i in 1:nBus
        for j in 1:nBus
            if i == j                                             # Diagonal elements of Jacobian sub-matrices
                # J11[j, j] = -q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j])
                J11[j, j] = -q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) + Qnor[j]
                J21[j, j] = p_inj[j, 1] - v_mag[j, 1]^2 * real(y_bus[j, j]) - Pnor[j]
                if order != 0
                    # J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j])  + (1 - ind[j]) * p0[j] * (order * dp[j, 1] * v_mag[j, 1] .^ (order - 1)) + ind[j] * p0[j] * (Gs[j] ./ Gs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                    J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j]) - Pnor[j] + (1 - ind[j]) * p0[j] * (order * dp[j, 1] * v_mag[j, 1] .^ (order - 1)) + ind[j] * p0[j] * (Gs[j] ./ Gs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                    J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) - Qnor[j] + (1 - ind[j]) * q0[j] * (order * dq[j, 1] * v_mag[j, 1] .^ (order - 1)) + ind[j] * q0[j] * (Bs[j] ./ Bs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                else
                    # J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j])  + ind[j] * p0[j] * (Gs[j] ./ Gs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                    J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j]) - Pnor[j] + ind[j] * p0[j] * (Gs[j] ./ Gs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                    J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) - Qnor[j] + ind[j] * q0[j] * (Bs[j] ./ Bs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                end
                J13[j, j] = ind[j] * p0[j] * (diff_Gfm[j] ./ Gs0[j]) .* (v_mag[j, 1] ./ v_mag0[j]) .^ 2
                J23[j, j] = ind[j] * q0[j] * (diff_Bfm[j] ./ Bs0[j]) .* (v_mag[j, 1] ./ v_mag0[j]) .^ 2

                J32[j, j] = -diff_TV[j]
                J33[j, j] = -diff_Tfm[j]

                J41[j, j] = Qnor[j]
                J42[j, j] = -Pnor[j]
                J24[j, j] = Pnor[j]
                J44[j, j] = -Qnor[j]



                # J14[j, j] = + sg[j] * ((v_mag[j, 1]  * e1m[j]) / X1d[j]) * cos(e1a[j] - v_angle[j])

                # J41[j, j] = -((v_mag[j, 1] * e1m[j]) / X1d[j]) * cos(e1a[j] - v_angle[j])
                # J42[j, j] =             (1 * e1m[j]  / X1d[j]) * sin(e1a[j] - v_angle[j]) * (v_mag[j, 1])
                # J44[j, j] = ((v_mag[j, 1]  * e1m[j]) / X1d[j]) * cos(e1a[j] - v_angle[j])

                J41[j, j] = -J41[j, j]
                J42[j, j] = -J42[j, j]
                J44[j, j] = -J44[j, j]
            else                                                  # Off-diagonal elements of Jacobian sub-matrices

                theta_diff_ij = v_angle[i, 1] - v_angle[j, 1]
                real_val = real(y_bus[i, j])
                imag_val = imag(y_bus[i, j])
                sin_val = sin(theta_diff_ij)
                cos_val = cos(theta_diff_ij)
                partial_p_inj = +1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * sin_val - imag_val * cos_val)
                partial_q_inj = -1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * cos_val + imag_val * sin_val)
                J11[i, j] = partial_p_inj
                J21[i, j] = partial_q_inj
                J12[i, j] = -partial_q_inj
                J22[i, j] = partial_p_inj
            end
        end
    end



    J11_tmp = J11[1:end, 1:end]
    J12_tmp = J12[1:end, 1:end]
    J21_tmp = J21[1:end, 1:end]
    J22_tmp = J22[1:end, 1:end]

    J13_tmp = J13[1:end, idu]
    J23_tmp = J23[1:end, idu]
    J31_tmp = J31[idu, 1:end]
    J32_tmp = J32[idu, 1:end]
    J33_tmp = J33[idu, idu]

    J41_tmp = J41[ign, 1:end]
    J42_tmp = J42[ign, 1:end]
    J43_tmp = J43[ign, idu]
    J44_tmp = J44[ign, ign]
    J14_tmp = J14[1:end, ign]
    J24_tmp = J24[1:end, ign]
    J34_tmp = J34[idu, ign]

    J = [J11_tmp J12_tmp J13_tmp J14_tmp; J21_tmp J22_tmp J23_tmp J24_tmp; J31_tmp J32_tmp J33_tmp J34_tmp; J41_tmp J42_tmp J43_tmp J44_tmp]
    # J = [J11_tmp J12_tmp J13_tmp; J21_tmp J22_tmp J23_tmp ; J31_tmp J32_tmp J33_tmp ]

    # ires = setdiff(collect(1:size(J, 2)), isl)
    # J = J[:, ires]
    # ires = setdiff(collect(1:size(J, 1)), isl)
    # J = J[ires, :]

    J[:, isl] = zeros(2 * nBus + length(idu) + length(ign), 1)

    diff_pf = (1 .- ind) .* p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3]) .* kip + ind .* p0 .* (diff_Gf ./ Gs0) .* (v_mag ./ v_mag0) .^ 2
    diff_qf = (1 .- ind) .* q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3]) .* kiq + ind .* q0 .* (diff_Bf ./ Bs0) .* (v_mag ./ v_mag0) .^ 2

    J[1:nBus, isl] = J[1:nBus, isl] + diff_pf
    J[nBus+1:2*nBus, isl] = J[nBus+1:2*nBus, isl] + diff_qf
    J[2*nBus+1:2*nBus+length(idu), isl] = J[2*nBus+1:2*nBus+length(idu), isl] - diff_Tf[idu]
    J[2*nBus+length(idu)+1:end, isl] = J[2*nBus+length(idu)+1:end, isl] - ki[ign]

    return J
end











function jacobianV_IND_GFM(nBus, isl, ign, ipq, idu, v_angle, v_mag, y_bus, p_inj, q_inj, diff_Tf, diff_Tfm, diff_TV, diff_Gf, diff_Gfm, diff_Bf, diff_Bfm, ki, kip, kiq, p0, q0, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind)
    # v_angle=v_angle_pf
    # v_mag= v_mag_pf
    # p0=PloadBUS
    # q0=QloadBUS
    # v_mag0=Vi0


    J11 = zeros(nBus, nBus)
    J12 = zeros(nBus, nBus)
    J21 = zeros(nBus, nBus)
    J22 = zeros(nBus, nBus)

    J13 = zeros(nBus, nBus)
    J23 = zeros(nBus, nBus)
    J32 = zeros(nBus, nBus)
    J33 = zeros(nBus, nBus)

    J41 = zeros(nBus, nBus)
    J42 = zeros(nBus, nBus)
    J43 = zeros(nBus, nBus)
    J51 = zeros(nBus, nBus)
    J52 = zeros(nBus, nBus)
    J53 = zeros(nBus, nBus)

    for i in 1:nBus
        for j in 1:nBus
            if i == j                                             # Diagonal elements of Jacobian sub-matrices
                J11[j, j] = -q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j])
                J21[j, j] = p_inj[j, 1] - v_mag[j, 1]^2 * real(y_bus[j, j])
                if order != 0
                    J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j]) + (1 - ind[j]) * p0[j] * (order * dp[j, 1] * v_mag[j, 1] .^ (order - 1)) + ind[j] * p0[j] * (Gs[j] ./ Gs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                    J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) + (1 - ind[j]) * q0[j] * (order * dq[j, 1] * v_mag[j, 1] .^ (order - 1)) + ind[j] * q0[j] * (Bs[j] ./ Bs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                else
                    J12[j, j] = p_inj[j, 1] + v_mag[j, 1]^2 * real(y_bus[j, j]) + ind[j] * p0[j] * (Gs[j] ./ Gs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                    J22[j, j] = q_inj[j, 1] - v_mag[j, 1]^2 * imag(y_bus[j, j]) + ind[j] * q0[j] * (Bs[j] ./ Bs0[j]) * (2 * (1 / (v_mag0[j]^2)) * v_mag[j, 1])
                end
                J13[j, j] = ind[j] * p0[j] * (diff_Gfm[j] ./ Gs0[j]) .* (v_mag[j, 1] ./ v_mag0[j]) .^ 2
                J23[j, j] = ind[j] * q0[j] * (diff_Bfm[j] ./ Bs0[j]) .* (v_mag[j, 1] ./ v_mag0[j]) .^ 2

                J32[j, j] = -diff_TV[j]
                J33[j, j] = -diff_Tfm[j]

                # J41[j, j] =  (2 * J11[j, j] * p_inj[j, 1] + 2 * J21[j, j] * q_inj[j, 1]) / (v_mag[j, 1] ^ 2)
                # J42[j, j] = ((2 * J12[j, j] * p_inj[j, 1] + 2 * J22[j, j] * q_inj[j, 1]) * (v_mag[j, 1] ^ 2) - 2 * v_mag[j, 1] * (p_inj[j, 1] ^ 2 + q_inj[j, 1] ^ 2))/(v_mag[j, 1] ^ 4)
                # J51[j, j] =  -   (J21[j, j] * p_inj[j, 1]     - J11[j, j] * q_inj[j, 1]) / (p_inj[j, 1] ^ 2 + q_inj[j, 1] ^ 2) + 1
                # J52[j, j] =  -   (J22[j, j] * p_inj[j, 1]     - J12[j, j] * q_inj[j, 1]) / (p_inj[j, 1] ^ 2 + q_inj[j, 1] ^ 2)

                J41[j, j] = J11[j, j] * (cos(v_angle[j, 1]) / v_mag[j, 1]) + p_inj[j, 1] * (-sin(v_angle[i, 1]) / v_mag[i, 1]) + J21[j, j] * (sin(v_angle[j, 1]) / v_mag[j, 1]) + q_inj[j, 1] * (cos(v_angle[i, 1]) / v_mag[i, 1])
                J42[j, j] = J12[j, j] * (cos(v_angle[j, 1]) / v_mag[j, 1]) + p_inj[j, 1] * (-cos(v_angle[i, 1]) / v_mag[i, 1]) + J22[j, j] * (sin(v_angle[j, 1]) / v_mag[j, 1]) + q_inj[j, 1] * (-sin(v_angle[i, 1]) / v_mag[i, 1])
                J51[j, j] = J11[j, j] * (-sin(v_angle[j, 1]) / v_mag[j, 1]) + p_inj[j, 1] * (-cos(v_angle[i, 1]) / v_mag[i, 1]) + J21[j, j] * (-cos(v_angle[j, 1]) / v_mag[j, 1]) + q_inj[j, 1] * (sin(v_angle[i, 1]) / v_mag[i, 1])
                J52[j, j] = J12[j, j] * (-sin(v_angle[j, 1]) / v_mag[j, 1]) + p_inj[j, 1] * (sin(v_angle[i, 1]) / v_mag[i, 1]) + J22[j, j] * (-cos(v_angle[j, 1]) / v_mag[j, 1]) + q_inj[j, 1] * (cos(v_angle[i, 1]) / v_mag[i, 1])




                # J41[j, j] =  - J41[j, j]
                # J42[j, j] =  - J42[j, j]
                # J51[j, j] =  - J51[j, j]
                # J52[j, j] =  - J52[j, j]
            else                                                  # Off-diagonal elements of Jacobian sub-matrices

                theta_diff_ij = v_angle[i, 1] - v_angle[j, 1]
                real_val = real(y_bus[i, j])
                imag_val = imag(y_bus[i, j])
                sin_val = sin(theta_diff_ij)
                cos_val = cos(theta_diff_ij)
                partial_p_inj = +1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * sin_val - imag_val * cos_val)
                partial_q_inj = -1 * (v_mag[i, 1] * v_mag[j, 1]) * (real_val * cos_val + imag_val * sin_val)
                J11[i, j] = partial_p_inj
                J21[i, j] = partial_q_inj
                J12[i, j] = -partial_q_inj
                J22[i, j] = partial_p_inj

                # J41[i, j] =  (2 * J11[i, j] * p_inj[i, 1] + 2 * J21[i, j] * q_inj[i, 1]) / (v_mag[i, 1] ^ 2)
                # J42[i, j] =  (2 * J12[i, j] * p_inj[i, 1] + 2 * J22[i, j] * q_inj[i, 1]) / (v_mag[i, 1] ^ 2)
                # J51[i, j] =  -   (J21[i, j] * p_inj[i, 1]     - J11[i, j] * q_inj[i, 1]) / (p_inj[i, 1] ^ 2 + q_inj[i, 1] ^ 2)
                # J52[i, j] =  -   (J22[i, j] * p_inj[i, 1]     - J12[i, j] * q_inj[i, 1]) / (p_inj[i, 1] ^ 2 + q_inj[i, 1] ^ 2)

                J41[i, j] = J11[i, j] * (cos(v_angle[i, 1]) / v_mag[i, 1]) + J21[i, j] * (sin(v_angle[i, 1]) / v_mag[i, 1])
                J42[i, j] = J12[i, j] * (cos(v_angle[i, 1]) / v_mag[i, 1]) + J22[i, j] * (sin(v_angle[i, 1]) / v_mag[i, 1])
                J51[i, j] = J11[i, j] * (-sin(v_angle[i, 1]) / v_mag[i, 1]) + J21[i, j] * (-cos(v_angle[i, 1]) / v_mag[i, 1])
                J52[i, j] = J12[i, j] * (-sin(v_angle[i, 1]) / v_mag[i, 1]) + J22[i, j] * (-cos(v_angle[i, 1]) / v_mag[i, 1])

                # J41[i, j] =  - J41[i, j]
                # J42[i, j] =  - J42[i, j]
                # J51[i, j] =  - J51[i, j]
                # J52[i, j] =  - J52[i, j]

            end
        end
    end


    J11_tmp = J11[ipq, 1:end]
    J12_tmp = J12[ipq, 1:end]
    J21_tmp = J21[ipq, 1:end]
    J22_tmp = J22[ipq, 1:end]

    J13_tmp = J13[ipq, idu]
    J23_tmp = J23[ipq, idu]
    J31_tmp = zeros(length(idu), nBus)
    J32_tmp = J32[idu, 1:end]
    J33_tmp = J33[idu, idu]






    ires = findall3(x -> x != isl, ign)
    ign1 = ign[ires, :]
    ign1 = vec(ign1)

    J41_tmp = J41[ign, 1:end]
    J42_tmp = J42[ign, 1:end]
    J43_tmp = J43[ign, idu]

    J51_tmp = J51[ign1, 1:end]
    J52_tmp = J52[ign1, 1:end]
    J53_tmp = J53[ign1, idu]




    J = [J11_tmp J12_tmp J13_tmp; J21_tmp J22_tmp J23_tmp; J31_tmp J32_tmp J33_tmp; J41_tmp J42_tmp J43_tmp; J51_tmp J52_tmp J53_tmp]
    ires = setdiff(collect(1:size(J, 2)), isl)
    J = J[:, ires]


    # J[:, isl] = ki

    # diff_pf = (1 .- ind) .* p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3]) .* kip + ind .* p0 .* (diff_Gf ./ Gs0) .* (v_mag ./ v_mag0) .^ 2
    # diff_qf = (1 .- ind) .* q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3]) .* kiq + ind .* q0 .* (diff_Bf ./ Bs0) .* (v_mag ./ v_mag0) .^ 2

    # J[1:nBus, isl] = J[1:nBus, isl] + diff_pf
    # J[nBus+1:nBus+length(ipq), isl] = J[nBus+1:nBus+length(ipq), isl] + diff_qf[ipq]
    # J[nBus+length(ipq)+1:end, isl] = J[nBus+length(ipq)+1:end, isl] - diff_Tf[idu]


    return J
end



function mismatch_vectorV(nBus, isl, ipv, ipq, v_angle, v_mag, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, p0, q0, dp, dq)
    p_inj = zeros(nBus)
    q_inj = zeros(nBus)
    p_calc = zeros(nBus)
    q_calc = zeros(nBus)

    p_plus = -ki * pim_pf
    p_V = p0 .* (dp[:, 1] .* v_mag .^ 2 .+ dp[:, 2] .* v_mag .+ dp[:, 3])
    q_V = q0 .* (dq[:, 1] .* v_mag .^ 2 .+ dq[:, 2] .* v_mag .+ dq[:, 3])

    for i in 1:nBus  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses

        theta_diff = v_angle[i, 1] .- v_angle
        real_val = real(y_bus[i, :])
        imag_val = imag(y_bus[i, :])
        sin_val = sin.(theta_diff)
        cos_val = cos.(theta_diff)

        p_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* cos_val + (imag_val) .* sin_val), dims=1)
        q_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* sin_val - (imag_val) .* cos_val), dims=1)


        p_inj[i] = p_inj_1[1]
        q_inj[i] = q_inj_1[1]


        p_calc[i] = p_inj[i] + p_plus[i] + p_V[i]
        q_calc[i] = q_inj[i] + q_V[i]
    end


    delta_p = p_sch_pf .- p_calc    # Mismatch active power
    delta_q = q_sch_pf .- q_calc    # Mismatch reactive power

    if isempty(ipv)                                                    # No PV bus in the network
        delta_mismatch_inj = vcat(delta_p, delta_q)
    else                                                                         # Existence of PV buses in a system
        # Remove the corresponding reactive power mismatch entry from the delta_Y vector

        delta_p_tmp = delta_p[1:end]
        delta_q_tmp = delta_q[ipq]
        delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp)
    end
    return delta_mismatch_inj, p_inj, q_inj
end




function mismatch_vectorV_order(nBus, isl, ipv, ipq, v_angle, v_mag, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, p0, q0, dp, dq, order)
    p_inj = zeros(nBus)
    q_inj = zeros(nBus)
    p_calc = zeros(nBus)
    q_calc = zeros(nBus)

    p_plus = -ki * pim_pf
    if order != 0
        p_V = p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3])
        q_V = q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3])
    else
        p_V = p0
        q_V = q0
    end

    for i in 1:nBus  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses

        theta_diff = v_angle[i, 1] .- v_angle
        real_val = real(y_bus[i, :])
        imag_val = imag(y_bus[i, :])
        sin_val = sin.(theta_diff)
        cos_val = cos.(theta_diff)

        p_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* cos_val + (imag_val) .* sin_val), dims=1)
        q_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* sin_val - (imag_val) .* cos_val), dims=1)


        p_inj[i] = p_inj_1[1]
        q_inj[i] = q_inj_1[1]


        p_calc[i] = p_inj[i] + p_plus[i] + p_V[i]
        q_calc[i] = q_inj[i] + q_V[i]
    end


    delta_p = p_sch_pf .- p_calc    # Mismatch active power
    delta_q = q_sch_pf .- q_calc    # Mismatch reactive power

    if isempty(ipv)                                                    # No PV bus in the network
        delta_mismatch_inj = vcat(delta_p, delta_q)
    else                                                                         # Existence of PV buses in a system
        # Remove the corresponding reactive power mismatch entry from the delta_Y vector

        delta_p_tmp = delta_p[1:end]
        delta_q_tmp = delta_q[ipq]
        delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp)
    end
    return delta_mismatch_inj, p_inj, q_inj
end

function mismatch_vector(rdata_buses, v_angle, v_mag, y_bus, p_node, q_node, idx_dPV_buses, p_sch_pf, q_sch_pf)
    idx = 0
    p_injection = Dict()
    q_injection = Dict()
    @inbounds for i in 1:size(rdata_buses, 1)                                                  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses
        bus_type = rdata_buses[i, 2]
        nd_num = rdata_buses[i, 1]

        if bus_type == 3
            idx = idx + 1
        elseif bus_type != 3                                                       # Injection for PQ and PV buses
            # idx_nd_bsheet = findall3(x->x==nd_num,rdata_buses[:,1])
            idx_nd_bsheet = convert(Int64, nd_num)
            idx_nd_bsheet = idx_nd_bsheet[1, 1] - idx

            y_bus_row = @view y_bus[idx_nd_bsheet+idx, :]

            connected_indices = Vector{Int}(undef, length(node_data_new[idx_nd_bsheet+idx].node_cnode) + 1)
            real_val = Vector{Float64}(undef, length(connected_indices))
            imag_val = Vector{Float64}(undef, length(connected_indices))
            sin_val = Vector{Float64}(undef, length(connected_indices))
            cos_val = Vector{Float64}(undef, length(connected_indices))

            # View of y_bus row to avoid copying
            y_bus_row = @view y_bus[idx_nd_bsheet+idx, :]

            # Get the connected indices and append the current index
            connected_indices[1:end-1] .= node_data_new[idx_nd_bsheet+idx].node_cnode
            connected_indices[end] = idx_nd_bsheet + idx

            # Extract real and imaginary parts of y_bus_row at connected indices
            real_val .= real(y_bus_row[connected_indices])
            imag_val .= imag(y_bus_row[connected_indices])

            # Calculate sin and cos values
            theta_diff = v_angle[idx_nd_bsheet+idx, 1] .- v_angle[connected_indices]
            sin_val .= sin.(theta_diff)
            cos_val .= cos.(theta_diff)

            # Magnitude of the current voltage
            v_mag_current = v_mag[idx_nd_bsheet+idx, 1]

            # Compute p_inj and q_inj
            v_mag_connected = v_mag[connected_indices]
            p_inj = v_mag_current * sum(v_mag_connected .* (real_val .* cos_val .+ imag_val .* sin_val))
            q_inj = v_mag_current * sum(v_mag_connected .* (real_val .* sin_val .- imag_val .* cos_val))


            push!(p_injection, [idx_nd_bsheet + idx] => p_inj)
            push!(q_injection, [idx_nd_bsheet + idx] => q_inj)
            p_node[idx_nd_bsheet, 1] = p_inj[1, 1]
            q_node[idx_nd_bsheet, 1] = q_inj[1, 1]
        end
    end
    delta_p = p_sch_pf .- p_node    # Mismatch active power
    delta_q = q_sch_pf .- q_node    # Mismatch reactive power

    if isempty(idx_dPV_buses)                                                    # No PV bus in the network
        delta_mismatch_inj = vcat(delta_p, delta_q)
    else                                                                         # Existence of PV buses in a system
        delta_q_tmp = deleteat!(delta_q, idx_dPV_buses)                           # Remove the corresponding reactive power mismatch entry from the delta_Y vector
        delta_mismatch_inj = vcat(delta_p, delta_q_tmp)
    end
    return delta_mismatch_inj, p_injection, q_injection
end






function mismatch_vectorV_orderIND(nBus, isl, ipv, ipq, idu, v_angle, v_mag, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, p0, q0, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind, Ld_ind)
    p_inj = zeros(nBus)
    q_inj = zeros(nBus)
    p_calc = zeros(nBus)
    q_calc = zeros(nBus)
    delta_tq_tmp = zeros(length(idu))
    diff_Tf_tmp = zeros(length(idu))
    diff_Tfm_tmp = zeros(length(idu))
    diff_TV_tmp = zeros(length(idu))
    diff_Gf_tmp = zeros(length(idu))
    diff_Gfm_tmp = zeros(length(idu))
    diff_Bf_tmp = zeros(length(idu))
    diff_Bfm_tmp = zeros(length(idu))

    p_plus = ki * df

    p_V = (1 .- ind) .* p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3]) .* (1 .+ kip * df) + ind .* p0 .* (Gs ./ Gs0) .* (v_mag ./ v_mag0) .^ 2
    q_V = (1 .- ind) .* q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3]) .* (1 .+ kiq * df) + ind .* q0 .* (Bs ./ Bs0) .* (v_mag ./ v_mag0) .^ 2



    for i in 1:nBus  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses

        theta_diff = v_angle[i, 1] .- v_angle
        real_val = real(y_bus[i, :])
        imag_val = imag(y_bus[i, :])
        sin_val = sin.(theta_diff)
        cos_val = cos.(theta_diff)

        p_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* cos_val + (imag_val) .* sin_val), dims=1)
        q_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* sin_val - (imag_val) .* cos_val), dims=1)


        p_inj[i] = p_inj_1[1]
        q_inj[i] = q_inj_1[1]


        p_calc[i] = p_inj[i] + p_plus[i] + p_V[i]
        q_calc[i] = q_inj[i] + q_V[i]
    end


    delta_p = p_sch_pf .- p_calc    # Mismatch active power
    delta_q = q_sch_pf .- q_calc    # Mismatch reactive power

    if isempty(ipv)                                                    # No PV bus in the network
        delta_mismatch_inj = vcat(delta_p, delta_q)
    else                                                                         # Existence of PV buses in a system
        # Remove the corresponding reactive power mismatch entry from the delta_Y vector

        delta_p_tmp = delta_p[1:end]
        delta_q_tmp = delta_q[ipq]

        (delta_tq_tmp, diff_Tf_tmp, diff_Tfm_tmp, diff_TV_tmp, diff_Gf_tmp, diff_Gfm_tmp, diff_Bf_tmp, diff_Bfm_tmp) = torque_motor(v_mag[idu], df, fm, Ld_ind)

        delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp, delta_tq_tmp)
        # delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp)

    end
    return delta_mismatch_inj, p_inj, q_inj, diff_Tf_tmp, diff_Tfm_tmp, diff_TV_tmp, diff_Gf_tmp, diff_Gfm_tmp, diff_Bf_tmp, diff_Bfm_tmp, p_V, q_V
end


function mismatch_vectorV_orderINDA(nBus, isl, ipv, ipq, idu, v_angle, v_mag, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, p0, q0, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind, Ld_ind, ctd)
    # nBus, isl, ipv, ipq, idu, v_angle, v_mag, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, p0, q0, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind, Ld_ind, ctd=nBus, isl, ipv, ipq, idu, v_angle_pf, v_mag_pf, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, PloadBUS, QloadBUS, dp, dq, order, Gs, Bs, Gs0, Bs0, Vi0, ind, Ld_ind, ctd
    p_inj = zeros(nBus)
    q_inj = zeros(nBus)
    p_calc = zeros(nBus)
    q_calc = zeros(nBus)
    delta_tq_tmp = zeros(length(idu))
    diff_Tf_tmp = zeros(length(idu))
    diff_Tfm_tmp = zeros(length(idu))
    diff_TV_tmp = zeros(length(idu))
    diff_Gf_tmp = zeros(length(idu))
    diff_Gfm_tmp = zeros(length(idu))
    diff_Bf_tmp = zeros(length(idu))
    diff_Bfm_tmp = zeros(length(idu))

    p_plus = ki * df

    p_V = (1 .- ind) .* p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3]) .* (1 .+ kip * df) + ind .* p0 .* (Gs ./ Gs0) .* (v_mag ./ v_mag0) .^ 2
    q_V = (1 .- ind) .* q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3]) .* (1 .+ kiq * df) + ind .* q0 .* (Bs ./ Bs0) .* (v_mag ./ v_mag0) .^ 2



    for i in 1:nBus  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses
        ctd_i = ctd[i]

        theta_diff = v_angle[i, 1] .- v_angle[ctd_i]
        real_val = real(y_bus[i, ctd_i])
        imag_val = imag(y_bus[i, ctd_i])
        sin_val = sin.(theta_diff)
        cos_val = cos.(theta_diff)

        p_inj_1 = v_mag[i, 1] * sum(v_mag[ctd_i] .* ((real_val) .* cos_val + (imag_val) .* sin_val), dims=1)
        q_inj_1 = v_mag[i, 1] * sum(v_mag[ctd_i] .* ((real_val) .* sin_val - (imag_val) .* cos_val), dims=1)


        p_inj[i] = p_inj_1[1]
        q_inj[i] = q_inj_1[1]


        p_calc[i] = p_inj[i] + p_plus[i] + p_V[i]
        q_calc[i] = q_inj[i] + q_V[i]
    end


    delta_p = p_sch_pf .- p_calc    # Mismatch active power
    delta_q = q_sch_pf .- q_calc    # Mismatch reactive power

    if isempty(ipv)                                                    # No PV bus in the network
        delta_mismatch_inj = vcat(delta_p, delta_q)
    else                                                                         # Existence of PV buses in a system
        # Remove the corresponding reactive power mismatch entry from the delta_Y vector

        delta_p_tmp = delta_p[1:end]
        delta_q_tmp = delta_q[ipq]

        (delta_tq_tmp, diff_Tf_tmp, diff_Tfm_tmp, diff_TV_tmp, diff_Gf_tmp, diff_Gfm_tmp, diff_Bf_tmp, diff_Bfm_tmp) = torque_motor(v_mag[idu], df, fm, Ld_ind)

        delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp, delta_tq_tmp)
        # delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp)

    end
    return delta_mismatch_inj, p_inj, q_inj, diff_Tf_tmp, diff_Tfm_tmp, diff_TV_tmp, diff_Gf_tmp, diff_Gfm_tmp, diff_Bf_tmp, diff_Bfm_tmp, p_V, q_V
end




function mismatch_vectorV_isGG(nBus, isl, ign, ipv, ipq, idu, v_angle, v_mag, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, p0, q0, e1m, e1a, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind, Ld_ind, sg, X1d)
    p_inj = zeros(nBus)
    q_inj = zeros(nBus)
    p_calc = zeros(nBus)
    q_calc = zeros(nBus)
    q_plus = zeros(nBus)
    delta_tq_tmp = zeros(length(idu))
    diff_Tf_tmp = zeros(length(idu))
    diff_Tfm_tmp = zeros(length(idu))
    diff_TV_tmp = zeros(length(idu))
    diff_Gf_tmp = zeros(length(idu))
    diff_Gfm_tmp = zeros(length(idu))
    diff_Bf_tmp = zeros(length(idu))
    diff_Bfm_tmp = zeros(length(idu))

    p_plus = ki * df
    q_plus = sg .* (v_mag ./ X1d) .* (e1m .* cos.(e1a .- v_angle) .- v_mag)

    p_V = (1 .- ind) .* p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3]) .* (1 .+ kip * df) + ind .* p0 .* (Gs ./ Gs0) .* (v_mag ./ v_mag0) .^ 2
    q_V = (1 .- ind) .* q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3]) .* (1 .+ kiq * df) + ind .* q0 .* (Bs ./ Bs0) .* (v_mag ./ v_mag0) .^ 2

    Pe = ((e1m .* v_mag) ./ X1d) .* sin.(e1a .- v_angle)

    Pnor = v_mag .* e1m .* cos.(e1a .- v_angle)
    Qnor = -v_mag .* e1m .* sin.(e1a .- v_angle)




    for i in 1:nBus  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses

        theta_diff = v_angle[i, 1] .- v_angle
        real_val = real(y_bus[i, :])
        imag_val = imag(y_bus[i, :])
        sin_val = sin.(theta_diff)
        cos_val = cos.(theta_diff)

        p_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* cos_val + (imag_val) .* sin_val), dims=1)
        q_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* sin_val - (imag_val) .* cos_val), dims=1)


        p_inj[i] = p_inj_1[1]
        q_inj[i] = q_inj_1[1]


        # p_calc[i] = p_inj[i] + p_plus[i] + p_V[i]
        # q_calc[i] = q_inj[i] - q_plus[i] + q_V[i]

        p_calc[i] = p_inj[i] + p_V[i]
        q_calc[i] = q_inj[i] + q_V[i]
    end


    # delta_p = p_sch_pf .- p_calc    # Mismatch active power
    # delta_q = q_sch_pf .- q_calc    # Mismatch reactive power

    # delta_p = Pe .- p_calc    # Mismatch active power
    # delta_q = q_sch_pf .- q_calc    # Mismatch reactive power


    delta_p = Pnor .- p_calc    # Mismatch active power
    delta_q = Qnor .- q_calc    # Mismatch reactive power
    delta_ps = Pnor .- p_sch_pf .+ p_plus[1:nBus]

    if isempty(ipv)                                                    # No PV bus in the network
        delta_mismatch_inj = vcat(delta_p, delta_q)
    else                                                                         # Existence of PV buses in a system
        # Remove the corresponding reactive power mismatch entry from the delta_Y vector
        ires = findall3(x -> x != isl, collect(1:nBus))
        delta_p_tmp = delta_p[1:end]
        delta_q_tmp = delta_q[1:end]
        delta_ps_tmp = delta_ps[ign]

        # Pm = p_sch_pf .- p_plus[1:nBus]
        # Pe = ((e1m .* v_mag) ./ X1d) .* sin.(e1a .- v_angle)
        # dP = Pm - Pe
        # delta_P_tmp =  dP[ign]

        (delta_tq_tmp, diff_Tf_tmp, diff_Tfm_tmp, diff_TV_tmp, diff_Gf_tmp, diff_Gfm_tmp, diff_Bf_tmp, diff_Bfm_tmp) = torque_motor(v_mag[idu], df, fm, Ld_ind)

        delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp, delta_tq_tmp, delta_ps_tmp)
        # delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp, delta_tq_tmp)
    end
    return delta_mismatch_inj, p_inj, q_inj, Pnor, Qnor, diff_Tf_tmp, diff_Tfm_tmp, diff_TV_tmp, diff_Gf_tmp, diff_Gfm_tmp, diff_Bf_tmp, diff_Bfm_tmp, p_V, q_V
end



















function mismatch_vectorV_IND_GFM(nBus, isl, ign, ipq, idu, v_angle, v_mag, df, fm, y_bus, p_sch_pf, q_sch_pf, ki, kip, kiq, p0, q0, i20, ia0, dp, dq, order, Gs, Bs, Gs0, Bs0, v_mag0, ind, Ld_ind)
    p_inj = zeros(nBus)
    q_inj = zeros(nBus)
    p_calc = zeros(nBus)
    q_calc = zeros(nBus)
    i2_calc = zeros(nBus)
    ia_calc = zeros(nBus)

    delta_tq_tmp = zeros(length(idu))
    diff_Tf_tmp = zeros(length(idu))
    diff_Tfm_tmp = zeros(length(idu))
    diff_TV_tmp = zeros(length(idu))
    diff_Gf_tmp = zeros(length(idu))
    diff_Gfm_tmp = zeros(length(idu))
    diff_Bf_tmp = zeros(length(idu))
    diff_Bfm_tmp = zeros(length(idu))

    p_plus = ki * df

    p_V = (1 .- ind) .* p0 .* (dp[:, 1] .* v_mag .^ order .+ dp[:, 2] .* v_mag .+ dp[:, 3]) .* (1 .+ kip * df) + ind .* p0 .* (Gs ./ Gs0) .* (v_mag ./ v_mag0) .^ 2
    q_V = (1 .- ind) .* q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3]) .* (1 .+ kiq * df) + ind .* q0 .* (Bs ./ Bs0) .* (v_mag ./ v_mag0) .^ 2



    for i in 1:nBus  # Power Injection is calculated for both PV and PQ buses although there is no need for the calculation of Q injection for PV buses

        theta_diff = v_angle[i, 1] .- v_angle
        real_val = real(y_bus[i, :])
        imag_val = imag(y_bus[i, :])
        sin_val = sin.(theta_diff)
        cos_val = cos.(theta_diff)

        p_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* cos_val + (imag_val) .* sin_val), dims=1)
        q_inj_1 = v_mag[i, 1] * sum(v_mag .* ((real_val) .* sin_val - (imag_val) .* cos_val), dims=1)


        p_inj[i] = p_inj_1[1]
        q_inj[i] = q_inj_1[1]


        # p_calc[i] = p_inj[i] + p_plus[i] + p_V[i]
        p_calc[i] = p_inj[i] + p_V[i]
        q_calc[i] = q_inj[i] + q_V[i]

        # i2_calc[i] = (p_inj[i]^2 + q_inj[i]^2) / (v_mag[i, 1]^2)
        # ia_calc[i] = atan(q_inj[i]^2 / (-p_inj[i]^2)) + v_angle[i, 1]

        i2_calc[i] = p_inj[i] * (cos(v_angle[i, 1]) / v_mag[i, 1]) + q_inj[i] * (sin(v_angle[i, 1]) / v_mag[i, 1])
        ia_calc[i] = p_inj[i] * (-sin(v_angle[i, 1]) / v_mag[i, 1]) + q_inj[i] * (-cos(v_angle[i, 1]) / v_mag[i, 1])
    end


    delta_p = p_sch_pf .- p_calc    # Mismatch active power
    delta_q = q_sch_pf .- q_calc    # Mismatch reactive power
    delta_i2 = i20 .- i2_calc        # Mismatch current squred
    delta_ia = ia0 .- ia_calc        # Mismatch current angle


    # Remove the corresponding reactive power mismatch entry from the delta_Y vector

    # delta_p_tmp = delta_p[1:end]
    delta_p_tmp = delta_p[ipq]
    delta_q_tmp = delta_q[ipq]

    (delta_tq_tmp, diff_Tf_tmp, diff_Tfm_tmp, diff_TV_tmp, diff_Gf_tmp, diff_Gfm_tmp, diff_Bf_tmp, diff_Bfm_tmp) = torque_motor(v_mag[idu], df, fm, Ld_ind)

    ires = findall3(x -> x != isl, ign)
    ign1 = ign[ires, :]

    delta_i2_tmp = delta_i2[ign]
    delta_ia_tmp = delta_ia[ign1]


    delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp, delta_tq_tmp, delta_i2_tmp, delta_ia_tmp)
    # delta_mismatch_inj = vcat(delta_p_tmp, delta_q_tmp)

    return delta_mismatch_inj, p_inj, q_inj, diff_Tf_tmp, diff_Tfm_tmp, diff_TV_tmp, diff_Gf_tmp, diff_Gfm_tmp, diff_Bf_tmp, diff_Bfm_tmp, p_V, q_V
end

function solution_voltageIND(nBus, isl, ipq, ipv, v_angle, v_mag, df, fm, delta_x)

    if isempty(ipv)  # All buses are PQ buses and there is no PV bus

    else

        v_angle[ipv] = v_angle[ipv] + delta_x[ipv]
        v_angle[ipq] = v_angle[ipq] + delta_x[ipq]
        v_angle[isl] = 0.0

        v_mag[ipq] = v_mag[ipq] .* (1.0 .+ delta_x[nBus+1:nBus+length(ipq)])

        df = df .+ delta_x[isl]
        fm = fm .+ delta_x[nBus+length(ipq)+1:end]
    end
    return v_angle, v_mag, df, fm
end


function solution_voltage_isGG(nBus, isl, ipq, ipv, idu, ign, v_angle, v_mag, df, fm, e1_angle, delta_x)

    if isempty(ipv)  # All buses are PQ buses and there is no PV bus

    else

        # v_angle[ipv] = v_angle[ipv] + delta_x[ipv]
        # v_angle[ipq] = v_angle[ipq] + delta_x[ipq]

        # delta_x = insert!(delta_x, isl, 0.0)

        v_angle = delta_x[1:nBus]
        # v_angle[isl] = 0.0

        # v_mag[ipq] = v_mag[ipq] .* (1.0 .+ delta_x[nBus+1:nBus+length(ipq)])
        v_mag = v_mag .* (1.0 .+ delta_x[nBus+1:2*nBus])

        df = df .+ delta_x[isl]
        fm = fm .+ delta_x[2*nBus+1:2*nBus+length(idu)]
        e1_angle[ign] = e1_angle[ign] .+ delta_x[2*nBus+length(idu)+1:end]
        e1_angle[isl] = 0
    end
    return v_angle, v_mag, df, fm, del
end



function solution_voltageIND_GFM(nBus, isl, ipq, ipv, v_angle, v_mag, df, fm, delta_x)

    if isempty(ipv)  # All buses are PQ buses and there is no PV bus

    else
        delta_x = insert!(delta_x, isl, 0.0)

        v_angle[ipv] = v_angle[ipv] + delta_x[ipv]
        v_angle[ipq] = v_angle[ipq] + delta_x[ipq]
        # v_angle = delta_x[1:nBus, 1]
        v_angle[isl] = 0.0

        # v_mag[ipq] = v_mag[ipq] .* (1.0 .+ delta_x[nBus+1:nBus+length(ipq)])
        v_mag = v_mag .* (1.0 .+ delta_x[nBus+1:2*nBus])

        df = df .+ delta_x[isl]
        # fm = fm .+ delta_x[nBus+length(ipq)+1:end]
        fm = fm .+ delta_x[2*nBus+1:end]

    end
    return v_angle, v_mag, df, fm
end



function solution_voltage1(nBus, isl, ipq, ipv, v_angle, v_mag, pim_pf, delta_x)

    if isempty(ipv)  # All buses are PQ buses and there is no PV bus

    else

        v_angle[ipv] = v_angle[ipv] + delta_x[ipv]
        v_angle[ipq] = v_angle[ipq] + delta_x[ipq]
        v_angle[isl] = 0.0

        v_mag[ipq] = v_mag[ipq] .* (1.0 .+ delta_x[nBus+1:end])
        pim_pf = pim_pf .+ delta_x[isl]
    end
    return v_angle, v_mag, pim_pf
end





## -------------- Determining reactive generation at PV buses -------------####
# and switching PV bus to PQ bus in case Q-Gen hits the minimum or maximum limit#

function pv_pq_switching1(epsilon_suf, y_bus, v_mag, v_angle, ipv, ipq, mdf_idx, QgenBUS, q0, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G)
    q_gen = zeros(1)
    q_V = q0 .* (dq[:, 1] .* v_mag .^ 2 .+ dq[:, 2] .* v_mag .+ dq[:, 3])

    for i in ipv
        q_gen_max = QgenBUSmax[i, 1]
        q_gen_min = QgenBUSmin[i, 1]

        theta_diff = v_angle[i, 1] .- v_angle
        q_inj = v_mag[i, 1] * sum(v_mag .* ((real(y_bus[i, :])) .* sin.(theta_diff) - (imag(y_bus[i, :])) .* cos.(theta_diff)), dims=1)
        q_gen = q_inj + [q_V[i] - QshuntBUS[i]]

        if q_gen[1, 1] > q_gen_max + epsilon_suf                                     # Generated power becomes greator than the upper reactive power generation limit
            ipq = [ipq; [i]]  # Switch to PQ bus
            ipq = sort(ipq)
            idx_ipv = findall3(x -> x == i, ipv)
            ipv = ipv[1:end.!=idx_ipv]
            QgenBUS[i] = q_gen_max
            mdf_idx[i] = 1

            idx_G = findall3(x -> x == i, G[:, 1])
            idx_G = idx_G[1]
            G[idx_G, 12] = 0

        elseif q_gen[1, 1] < q_gen_min - epsilon_suf                                     # Generated power becomes greator than the upper reactive power generation limit
            ipq = [ipq; [i]] # Switch to PQ bus
            ipq = sort(ipq)
            idx_ipv = findall3(x -> x == i, ipv)
            ipv = ipv[1:end.!=idx_ipv]
            QgenBUS[i] = q_gen_min
            mdf_idx[i] = 1

            idx_G = findall3(x -> x == i, G[:, 1])
            idx_G = idx_G[1]
            G[idx_G, 12] = 0

        elseif q_gen_min - epsilon_suf <= q_gen[1, 1] <= q_gen_max + epsilon_suf && mdf_idx[i] == 1   # Qgen is within limits but the bus has switched from PV to PQ in the previous iteration
            mdf_idx[i] = 0
            ipv = [ipv; [i]] # Again Switch to PV bus
            ipv = sort(ipv)
            idx_ipq = findall3(x -> x == i, ipq)
            ipq = ipq[1:end.!=idx_ipq]

            idx_G = findall3(x -> x == i, G[:, 1])
            idx_G = idx_G[1]
            G[idx_G, 12] = 1
        end
    end



    return mdf_idx, QgenBUS, ipv, ipq, G
end


function pv_pq_switching1_order(epsilon_suf, y_bus, q_inj, v_mag, v_angle, ipv, ipq, mdf_idx, QgenBUS, q0, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G, order, Bs, Bs0, v_mag0, ind)
    # epsilon_suf, y_bus, q_inj, v_mag, v_angle, ipv, ipq, mdf_idx, QgenBUS, q0, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G, order, Bs, Bs0, v_mag0, ind=epsilon_suf, y_bus, q_inj, v_mag_pf, v_angle_pf, ipv, ipq, mdf_idx, QgenBUS, QloadBUS, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G, order, Bs, Bs0, Vi0, ind
    q_gen = zeros(1)

    if order != 0
        q_V = (1 .- ind) .* q0 .* (dq[:, 1] .* v_mag .^ order .+ dq[:, 2] .* v_mag .+ dq[:, 3]) + ind .* q0 .* (Bs ./ Bs0) .* (v_mag ./ v_mag0) .^ 2
    else
        q_V = (1 .- ind) .* q0 + ind .* q0 .* (Bs ./ Bs0) .* (v_mag ./ v_mag0) .^ 2
    end


    for i in ipv
        q_gen_max = QgenBUSmax[i, 1]
        q_gen_min = QgenBUSmin[i, 1]

        # theta_diff = v_angle[i, 1] .- v_angle
        # q_inj = v_mag[i, 1] * sum(v_mag .* ((real(y_bus[i, :])) .* sin.(theta_diff) - (imag(y_bus[i, :])) .* cos.(theta_diff)), dims=1)
        q_gen = q_inj[i, 1] + q_V[i] - QshuntBUS[i]

        if q_gen[1, 1] > q_gen_max + epsilon_suf                                     # Generated power becomes greator than the upper reactive power generation limit
            ipq = [ipq; [i]]  # Switch to PQ bus
            ipq = sort(ipq)
            idx_ipv = findall3(x -> x == i, ipv)
            ipv = ipv[1:end.!=idx_ipv]
            QgenBUS[i] = q_gen_max
            mdf_idx[i] = 1

            idx_G = findall3(x -> x == i, G[:, 1])
            idx_G = idx_G[1]
            G[idx_G, 12] = 0

        elseif q_gen[1, 1] < q_gen_min - epsilon_suf                                     # Generated power becomes greator than the upper reactive power generation limit
            ipq = [ipq; [i]] # Switch to PQ bus
            ipq = sort(ipq)
            idx_ipv = findall3(x -> x == i, ipv)
            ipv = ipv[1:end.!=idx_ipv]
            QgenBUS[i] = q_gen_min
            mdf_idx[i] = 1

            idx_G = findall3(x -> x == i, G[:, 1])
            idx_G = idx_G[1]
            G[idx_G, 12] = 0

        elseif q_gen_min - epsilon_suf <= q_gen[1, 1] <= q_gen_max + epsilon_suf && mdf_idx[i] == 1   # Qgen is within limits but the bus has switched from PV to PQ in the previous iteration
            mdf_idx[i] = 0
            ipv = [ipv; [i]] # Again Switch to PV bus
            ipv = sort(ipv)
            idx_ipq = findall3(x -> x == i, ipq)
            ipq = ipq[1:end.!=idx_ipq]

            idx_G = findall3(x -> x == i, G[:, 1])
            idx_G = idx_G[1]
            G[idx_G, 12] = 1
        end
    end



    return mdf_idx, QgenBUS, ipv, ipq, G
end

function findall3(f, a::Array{T,N}) where {T,N}
    j = 1
    b = Vector{Int}(undef, length(a))
    @inbounds for i in eachindex(a)
        @inbounds if f(a[i])
            b[j] = i
            j += 1
        end
    end
    resize!(b, j - 1)
    sizehint!(b, length(b))
    return b
end


function find_avb(x1, x2)
    # x1=bi
    # x2=L[:,1]

    x2_av = []

    for i in 1:length(x2)
        if length(intersect(x2[i], x1)) != 0
            push!(x2_av, i)
        end
    end

    return x2_av
end


function find_ones(M::SparseMatrixCSC{Int64,Int64})
    return [findall(x -> x == 1, row) for row in eachrow(M)]
end




function line_flow1(L, y_bus, V, del)
    fb = L[:, 1]
    tb = L[:, 2]
    g = L[:, 3]
    b = L[:, 4]
    bsh = L[:, 6]
    smax = L[:, 7]
    Imax=L[:, 7]

    y = g .+ im * b
    # V = vol_nodes_mag_pf
    # del = vol_nodes_theta_pf



    sij_s = zeros(length(fb))
    sij_r = zeros(length(fb))
    sij = zeros(length(fb))
    Iij = zeros(length(fb))
    percentIij = zeros(length(fb))



    for l in 1:length(fb)
        Vs = V[fb[l]] .* exp(im * del[fb[l]])
        Vr = V[tb[l]] .* exp(im * del[tb[l]])

        sij_s[l] = abs(Vs * (Vs - Vr) * y[l])
        sij_r[l] = abs(Vr * (Vs - Vr) * y[l])
        sij[l] = maximum([sij_s[l], sij_r[l]]) / smax[l]
        Iij[l] = abs((Vs - Vr) * y[l])
    end
    percentIij = Iij./Imax
    violated_lines = findall3(x -> x >= 1, sij)
    return sij, violated_lines, Iij, percentIij
end



function line_fault(L, G, Ld, y_bus, V, del)
    y_bus = y_bus_i
    V = vmag_i
    del = vang_i



    nBus = length(V)
    fb = L[:, 1]
    tb = L[:, 2]
    g = L[:, 3]
    b = L[:, 4]

    l = 8
    d = 0
    zb = (345^2) / 100
    ib = 100 / (sqrt(3) * 345)
    Rf = 0 / zb




    Inet = zeros(ComplexF64, nBus)
    V0 = zeros(ComplexF64, nBus)
    V_fault = zeros(ComplexF64, nBus)

    z = zeros(ComplexF64, length(fb))
    Iij_s = zeros(ComplexF64, length(fb))
    I_fault = zeros(ComplexF64, length(fb))

    y = g .+ im * b
    z = ones(length(fb)) ./ y

    V0 = V .* exp.(im * del)




    for i in 1:size(G, 1)
        ig = G[i, 1]
        zg = im * G[i, 20]
        # zg=im*0.2

        y_bus[ig, ig] = y_bus[ig, ig] + 1 / (zg)
    end

    for i in 1:size(Ld, 1)
        id = Ld[i, 1]
        pd = Ld[i, 2]
        qd = Ld[i, 3]
        vd = V0[id]

        yd = -(pd - im * qd) / (abs(vd)^2)
        y_bus[id, id] = y_bus[id, id] + yd
    end


    z_bus = inv(y_bus)


    # pre_fault
    Vs = V0[fb[l]]
    Vr = V0[tb[l]]
    Iij_s[l] = (Vs - Vr) / z[l]

    Vf_pre = Vs - d * z[l] * Iij_s[l]


    IF = Vf_pre / (Rf + z_bus[fb[l], fb[l]])

    Inet[fb[l]] = -IF

    V_fault = V0 + z_bus * Inet

    V_fault[14] = V0[14] + z_bus[14, 4] * IF

    for l in 1:length(fb)

        # pre_fault
        Vs = V_fault[fb[l]]
        Vr = V_fault[tb[l]]
        I_fault[l] = (Vs - Vr) / z[l]
    end



    # y1 = lu(y_bus, check=false)
    # eye_n=Matrix{Float64}(I, nBus, nBus)
    # z_bus1 = y1 \eye_n






    for l in 1:length(fb)

        # pre_fault
        Vs = V[fb[l]] .* exp(im * del[fb[l]])
        Vr = V[tb[l]] .* exp(im * del[tb[l]])
        Iij_s[l] = (Vs - Vr) / z[l]


        Vf_pre = Vs - d * z[l] * Iij_s[l]


        IF = Vf_pre / (Rf + z_bus[fb[l], fb[l]])

        Inet[fb[l]] = -IF

        V_fault = V0 + z_bus * Inet









        # pure_fault

        Zii = z_bus[fb[l], fb[l]]
        Zjj = z_bus[tb[l], tb[l]]
        k1 = (z[l]^2) * Zjj
        k2 = -z[l]^2 + 2 * z[l] * Zjj
        k3 = Zjj - Zii - z[l]

        rts = roots(Polynomial([k1, k2, k3]))
        ZH = rts[1]
        ZG = (Zii * z[l] + Zii * ZH) / (z[l] + ZH - Zii)


        # ZG1=(Zjj*z[l]+Zjj*ZH-ZH*z[l])/(ZH-Zjj)
        # println(ZG-ZG1)
        # Zii_vrf=(ZG*(z[l]+ZH))/(ZG+z[l]+ZH)
        # Zjj_vrf=(ZH*(z[l]+ZG))/(ZG+z[l]+ZH)
        # println(Zii_vrf-Zii)
        # println(abs(Zjj_vrf-Zjj))


        Zth = (ZG + d * z[l]) * (ZH + (1 - d) * z[l]) / (ZG + ZH + z[l])
        IF = Vf_pre / (Rf + Zth)
        d_Vf = -Vf_pre + Rf * IF
        d_IG = -d_Vf / (ZG + d * z[l])

        # fault
        Vf = Vf_pre + d_Vf
        IG = Iij_s[l] + d_IG

        VG = Vf + d * z[l] * IG



    end

    return Iij_s
end



function print_multislack_acpf(vmag_i, vang_i, p_i, q_i, p_V, q_V, cg_i, df_i, sm_i, Gs, Bs, itr_i)

    if cg_i == 1
        println(" Success!")

        XLSX.openxlsx("output_data\\R_ResultsE1.xlsx", mode="w") do xf
            Sheet = xf[1]
            XLSX.rename!(Sheet, "Voltage_ph")
            Sheet["A1"] = ["V"]
            Sheet["A2", dim=1] = vmag_i

            Sheet["B1"] = ["delt"]
            Sheet["B2", dim=1] = vang_i * 180 / pi

            Sheet["C1"] = ["Pnet"]
            Sheet["C2", dim=1] = p_i

            Sheet["D1"] = ["Qnet"]
            Sheet["D2", dim=1] = q_i

            Sheet["E1"] = ["sm"]
            Sheet["E2", dim=1] = sm_i

            Sheet["F1"] = ["p_V"]
            Sheet["F2", dim=1] = p_V

            Sheet["G1"] = ["q_V"]
            Sheet["G2", dim=1] = q_V

            Sheet["H1"] = ["Gs"]
            Sheet["H2", dim=1] = Gs

            Sheet["I1"] = ["Bs"]
            Sheet["I2", dim=1] = Bs

        end

        println(" power flow has been convereged in $itr_i iterations")
        println(" the frequency will be changed as :        $df_i pu")

    else
        println(" power flow has been diverged!")
    end
    println("[-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]")

end


function print_v(vmag_i_d)
    XLSX.openxlsx("output_data\\R_ResultsV.xlsx", mode="w") do xf
        Sheet = xf[1]
        XLSX.rename!(Sheet, "Voltage_d")
        # Sheet["A1", dim=1] = vmag_i_d
        Sheet["A1"] = vmag_i_d

    end
end




function statsL(Btemp, Ltemp, Ttemp, Gtemp, Ldtemp, Shutemp, Lout_stat)

    nscn = size(Lout_stat, 1)
    nout = size(Lout_stat, 2)
    LNS_stat = zeros(Float64, nscn)
    # Lout_stat = zeros(Int, nscn, nout)

    i = 1
    while i <= nscn
        Btemp = deepcopy(Bs)
        Ltemp = deepcopy(Ls)
        Ttemp = deepcopy(Ts)
        Gtemp = deepcopy(Gs)
        Ldtemp = deepcopy(Lds)
        Shutemp = deepcopy(Shus)
        println(sum(Ldtemp[:, 2]))
        println(i)


        Lout = Lout_stat[i, :]
        # Lout = rand(1:size(Ltemp, 1), nout)
        # Lout_stat[i, :] = Lout

        (Rn0, H0, u0, B0, L0, T0, G0, Ld0, Shu0) = CSC(Btemp, Ltemp, Ttemp, Gtemp, Ldtemp, Shutemp, Lout)
        LNS_v = sum(Ldsh_v) * sb
        LNS_stat[i, :] .= LNS_v

        Ldsh_v[1:size(Lds, 1)] = zeros(size(Lds, 1))
        i += 1
    end

    return LNS_stat
end



function statsB(Btemp, Ltemp, Ttemp, Gtemp, Ldtemp, Shutemp, Bout_stat)

    nscn = size(Bout_stat, 1)
    nout = size(Bout_stat, 2)
    LNS_stat = zeros(Float64, nscn)
    Lnew_stat = zeros(Int64, nscn)
    Ltot_stat = zeros(Int64, nscn)
    Nlout_stat = zeros(Int, nscn)

    i = 1
    while i <= nscn
        Btemp = deepcopy(Bs)
        Ltemp = deepcopy(Ls)
        Ttemp = deepcopy(Ts)
        Gtemp = deepcopy(Gs)
        Ldtemp = deepcopy(Lds)
        Shutemp = deepcopy(Shus)
        println(sum(Ldtemp[:, 2]))
        println(i)


        bout = Bout_stat[i, :]
        Lout = B2Lout(bout, Ltemp)
        Nlout_stat[i] = length(Lout)

        (Rn0, H0, u0, B0, L0, T0, G0, Ld0, Shu0, new_L, tot_L) = CSC(Btemp, Ltemp, Ttemp, Gtemp, Ldtemp, Shutemp, Lout)
        LNS_v = sum(Ldsh_v) * sb
        LNS_stat[i, :] .= LNS_v
        Lnew_stat[i, :] .= new_L
        Ltot_stat[i, :] .= tot_L

        Ldsh_v[1:size(Lds, 1)] = zeros(size(Lds, 1))
        i += 1
    end

    return LNS_stat, Nlout_stat, Lnew_stat, Ltot_stat
end




function B2Lout(bout, L)
    fb = L[:, 1]
    tb = L[:, 2]
    Lout = []
    for i in 1:length(bout)
        for j in 1:length(fb)
            if fb[j] == bout[i] || tb[j] == bout[i]
                Lout = union(Lout, j)
            end
        end
    end

    return Lout
end