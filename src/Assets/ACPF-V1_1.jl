
function PF(raw_temp,CSP,Prof,t)
    raw_temp["Gen"]."Ramp"=convert.(Float64,raw_temp["Gen"]."Ramp")
    raw_temp["Gen"]."Ramp".=0.02; #Mousa selections: ask him for details
    raw_temp["Trans"]=select(SysData["Trans"],Not("g_sh (pu)"));
    raw_temp["Trans"]."g_sh (pu)"=SysData["Trans"]."g_sh (pu)" #bringing shunt conductance to last column

    PG=Matrix(CSP["PG"])[t,:];QG=Matrix(CSP["QG"])[t,:];
    if isempty(CSP["Pshed"]) Shed= CSP["Pshed"] else Shed=Matrix(CSP["Pshed"])[t,:] end;
    if isempty(CSP["Pcurt"]) Curt= CSP["Pcurt"] else Curt=Matrix(CSP["Pcurt"])[t,:] end;
    if isempty(CSP["Pinc"]) Pinc= CSP["Pinc"] else Pinc=Matrix(CSP["Pinc"])[t,:] end;
    if isempty(CSP["Pdec"]) Pdec= CSP["Pdec"] else Pdec=Matrix(CSP["Pdec"])[t,:] end;
    if isempty(CSP["Pdis"]) Pdis= CSP["Pdis"] else Pdis=Matrix(CSP["Pdis"])[t,:] end;
    if isempty(CSP["Pch"]) Pch= CSP["Pch"] else Pch=Matrix(CSP["Pch"])[t,:] end;

    raw_temp["RES"]."PG"=convert.(Float64,raw_temp["RES"]."PG")
    raw_temp["ESS"]."PG"=convert.(Float64,raw_temp["ESS"]."PG")
    if isempty(Prof["Wind"]) raw_temp["RES"]."PG".=0.0 else raw_temp["RES"]."PG".= Matrix(Prof["Wind"])[t,:] end;
    if isempty(Prof["PL"]) raw_temp["Load"]."Pd (pu)".=0.0 else raw_temp["Load"]."Pd (pu)".= Matrix(Prof["PL"])[t,:] end;
    if isempty(Prof["QL"]) raw_temp["Load"]."Qd (pu)".=0.0 else raw_temp["Load"]."Qd (pu)".= Matrix(Prof["QL"])[t,:] end;
    raw_temp["Gen"]."V_set (pu)".=Matrix(CSP["Vmag"])[t,:][SysData["Gen"]."Bus"];


    raw_temp["Gen"]."Pg (pu)".=PG;raw_temp["Gen"]."Qg (pu)".=QG;
    if !isempty(Pinc) raw_temp["Load"]."Pd (pu)"[SysData["FL"]."Bus"].+= ((Pinc)+ (Pdec*-1)) end;
    if !isempty(Pdis) raw_temp["ESS"]."PG".= Pdis+ (-1* Pch) end;
    if !isempty(Curt) raw_temp["RES"]."PG".+=(-1*Curt); end
    if !isempty(Shed) raw_temp["Load"]."Pd (pu)".+= -1*Shed; end


    y_bus = Y_bus1(raw_temp["nBus"],raw_temp["Line"],raw_temp["Trans"])    #(nBus, L, T)
    vmag_i, vang_i, ind_converge,Iij,percent,PGGG = multislack_acpf(y_bus, raw_temp["Bus"], raw_temp["Line"], raw_temp["Gen"], raw_temp["Load"], raw_temp["Shunt"],raw_temp["RES"],raw_temp["ESS"]) #(y_bus_i, B , L, G, Ld, Shu)

    # voltage mag, voltage angle, mismatch injection,
    ov=findall(x-> x>1.05+ϵ,vmag_i)
    uv=findall(x-> x<0.95-ϵ,vmag_i)
    lo=raw_temp["Line"]."Line_N"[findall(x-> x>1+ϵ   ,percent)]
    lo_index=findall(x-> x>1+ϵ,percent)

    return ov,uv,lo,lo_index,vmag_i,vang_i,ind_converge,Iij,percent,PGGG
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

function multislack_acpf(y_bus,
    B, L, G, Ld, Shu,REN,BAT)

    # B=SysData["Bus"]; L=SysData["Line"]; G=SysData["Gen"]; Ld=SysData["Load"]; Shu=SysData["Shunt"];REN=SysData["RES"];BAT=SysData["ESS"];
    # B=deepcopy(SysData["Bus"]); L=deepcopy(SysData["Line"]); G=deepcopy(SysData["Gen"]); Ld=deepcopy(SysData["Load"]); Shu=deepcopy(SysData["Shunt"]);REN=deepcopy(SysData["RES"]);BAT=deepcopy(SysData["ESS"]);
    G."Ramp"=convert.(Float64,G."Ramp")
    G."Ramp".=0.02; #Mousa selections: ask him for details


    #----------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 1: Settings--[-][-][-][-][-][-][-]---------
    #----------------------------------------------------------------------------------------------
    nBus = size(B, 1)
    PgenBUS = zeros(nBus)
    QgenBUS = zeros(nBus)
    PloadBUS = zeros(nBus)
    QloadBUS = zeros(nBus)
    QshuntBUS = zeros(nBus)
    QgenBUSmin = zeros(nBus)
    QgenBUSmax = zeros(nBus)
    PgenBUSmin = zeros(nBus)
    PgenBUSmax = zeros(nBus)
    dp= zeros(nBus, 3)
    dq= zeros(nBus, 3)


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
    ish = Shu[:, 1]
    slk = deepcopy(G[:, 12])

    global slk1=deepcopy(slk)

    PgenBUS[ipv] = G[:, 2]


    QgenBUSmax[ipv] = G[:, 4]
    QgenBUSmin[ipv] = G[:, 5]
    PgenBUSmax[ipv] = G[:, 7]
    PgenBUSmin[ipv] = G[:, 8]

    PloadBUS[iload] = Ld[:, 2]
    PloadBUS[REN[:,2]] -= REN[:, 5]
    PloadBUS[BAT[:,1]] -= BAT[:, 26]
    QloadBUS[iload] = Ld[:, 3]
    QloadBUS[REN[:,2]] -= REN[:, 6]
    QloadBUS[BAT[:,1]] -= BAT[:, 27]
    QshuntBUS[ish] = Shu[:, 2]

    dp[iload, 3].=1
    dq[iload, 3].=1
    dp[REN[:,2], 3].=1
    dq[REN[:,2], 3].=1
    dp[BAT[:,1], 3].=1
    dq[BAT[:,1], 3].=1

    npv = length(ipv)
    npq = length(ipq)





    max_mismatch_pf = 1     #
    epsilon = 1e-9       #convergence criteria
    epsilon_suf = 1e-6
    epsilon_suf_p = 1e-6
    itr_max = 100       # max number of iterations
    maxiter = []


    v_mag_pf = ones(nBus)
    v_angle_pf = zeros(nBus)
    v_mag_pf[ipv] = G[:, 6]





    #.............................................................................
    # Note: In general, all (ipv=ipv+isl) buses can be considered as the slack bus
    #.............................................................................
    pim_pf = zeros(1)

    ki = zeros(npv + npq + npq, 1)
    (ki, beti) = spc_slack(ki, G)


    #.............................................................................


    vol_nodes_mag_pf = zeros(nBus)
    vol_nodes_theta_pf = zeros(nBus)
    max_mismatch_pf = 1
    iteration = 0
    mdf_idx = zeros(nBus)
    mdf_idxG = zeros(nBus)
    convergence_indicator = 0
    delta_mismatch_inj = zeros(nBus + npq)
    delta_mismatch=[]

    #---------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 2: Calcultae Scheduled Powers   [-][-][-][-][-][-][-]
    #---------------------------------------------------------------------------------------------------
    (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)

    #---------------------------------------------------------------------------------------------------


    #---------------------------------------------------------------------------------------------------
    # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 3: Main Loop of the Powr Flow--[-][-][-][-][-][-][-][-]
    #---------------------------------------------------------------------------------------------------
    Jacob1=0;p_inj=0;q_inj=0;
    while max_mismatch_pf > epsilon && iteration <= itr_max
        # .......Step 3-1: Mismatch between  Scheduled and Calculated Powers............
        #...............................................................................
        (delta_mismatch_inj, p_inj, q_inj) = mismatch_vectorV(nBus, isl, ipv, ipq, v_angle_pf, v_mag_pf, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, PloadBUS, QloadBUS, dp, dq)



        # .......Step 3-2: Jacobian Matrix.............................................
        #..............................................................................
         Jacob1 = jacobianV(nBus, isl, ipq, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, ki, PloadBUS, QloadBUS, dp, dq)




        # .......Step 3-3: Inverse of Jacobian Matrix...................................
        #...............................................................................
        Jcb = lu(Jacob1, check=false)
        if !issuccess(Jcb)
            println("")
            show("Singularity observed for iteration= $iteration")
            println("")
            break
        end



        # .......Step 3-4: Update State Variables (V, delta,Pim)..........................
        #.................................................................................
        delta_x_pf = Jcb \ delta_mismatch_inj
        (v_angle_sv, v_mag_sv, pim_pf1) = solution_voltage1(nBus, isl, ipq, ipv, v_angle_pf, v_mag_pf, pim_pf, delta_x_pf)
        v_mag_pf = v_mag_sv                                                  # sv = solution voltage (Here, voltage is updated)
        v_angle_pf = v_angle_sv
        pim_pf = pim_pf1


        # .......Step 3-5: Check if PV-PQ Switching required ............................
        #.................................................................................
        max_mismatch_pf = maximum(abs.(delta_mismatch_inj))
        max_mismatch_pf_P = maximum(abs.(delta_mismatch_inj[1:nBus]))
        max_mismatch_pf_Q = maximum(abs.(delta_mismatch_inj[nBus+1:end]))
        if max_mismatch_pf_P < 0.05
            (mdf_idx, QgenBUS, ipv, ipq, G) = pv_pq_switching1(epsilon_suf, y_bus, v_mag_pf, v_angle_pf, ipv, ipq, mdf_idx, QgenBUS, QloadBUS, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G)

            if sum(mdf_idx) != 0
                println("pq switching occured in iteration $iteration")
                npv = length(ipv)
                npq = length(ipq)
                ki = zeros(npv + npq + npq, 1)
                (ki, beti) = spc_slack(ki, G)
                (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)
            end
        end


        iteration = iteration + 1
       push!(delta_mismatch,maximum(delta_mismatch_inj'))
    end # End while loop

    #Arman
    fullJacob = jacobianV(nBus, isl, ipq, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, ki, PloadBUS, QloadBUS, dp, dq;flag=1)
    JacobArman=deepcopy(fullJacob)
    JacobArman[:,nBus+1:end]=fullJacob[:,nBus+1:end]./repeat(v_mag_pf',2*nBus)
    JcbArman = lu(JacobArman, check=false)
    Sensivities=JcbArman\Matrix(1.0I,2*nBus,2*nBus)
    Sensivities[isl,:].=0
    #Arman

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


    for i in 1:nBus
        vol_nodes_mag_pf[i] = v_mag_pf[i]
        vol_nodes_theta_pf[i] = (v_angle_pf[i])                         # Angles are in Radian Because input of Cosine and Sine functions are in radian (i.e., Julia consideres them as radian)
    end

    (delta_mismatch_inj, p_inj, q_inj) = mismatch_vectorV(nBus, isl, ipv, ipq, v_angle_pf, v_mag_pf, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, PloadBUS, QloadBUS, dp, dq)


    (Iij, violated_lines, percent) = line_flow1(L, vol_nodes_mag_pf, vol_nodes_theta_pf)



    df = -pim_pf/sum(beti)
    df=df[1]



    G[:, 12]=deepcopy(slk)


    PGGG=(p_inj+PloadBUS)[findall(x -> x != 1, B[:, 2])]
    QGGG=(q_inj+QloadBUS)[findall(x -> x != 1, B[:, 2])]
    # PGGG=0
    # return vol_nodes_mag_pf, vol_nodes_theta_pf, delta_mismatch_inj, delta_mismatch, p_inj, q_inj, convergence_indicator, pim_pf, df, iteration,sij, violated_lines
    return p_inj,vol_nodes_mag_pf, vol_nodes_theta_pf, convergence_indicator,Iij,percent,PGGG,QGGG,Sensivities

end


function Y_bus1(nBus::Int64, L::DataFrame, T::DataFrame)
    T=select(T,Not("g_sh (pu)"));
    T."g_sh (pu)"=SysData["Trans"]."g_sh (pu)" #bringing shunt conductance to last column
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

    # ybus matrix without transformers
    for i in 1:nBus
        for j in 1:nBus

            if i == j
                ybus[i, i] = sum(yij[i, :])
            elseif i != j
                ybus[i, j] = -yij[i, j]

            end

        end
    end

    # ybus matrix with transformers
    fb_T = T[:, 1]
    tb_T = T[:, 2]
    g_T = T[:, 3]
    b_T = T[:, 4]
    y_T = g_T .+ im * b_T
    Tratio = 1.00 * T[:, 6]
    # Tratio = 1 ./ Tratio # becasue tap is on the primary side

    for tr in 1:length(fb_T)
        ybus[fb_T[tr], fb_T[tr]] = ybus[fb_T[tr], fb_T[tr]] + abs(Tratio[tr])^2 * y_T[tr]
        ybus[fb_T[tr], tb_T[tr]] = ybus[fb_T[tr], tb_T[tr]] - conj(Tratio[tr]) * y_T[tr]
        ybus[tb_T[tr], fb_T[tr]] = ybus[tb_T[tr], fb_T[tr]] - Tratio[tr] * y_T[tr]
        ybus[tb_T[tr], tb_T[tr]] = ybus[tb_T[tr], tb_T[tr]] + 1 * y_T[tr]
    end

    return ybus
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
        ki[indices_slacks] = beti / sum(beti)
    else
        beti = 100000
    end


    ki = convert(Array{Float64}, ki)
    return ki, beti
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

    if flag==1
        J = [J11 J12; J21 J22]
        J[:, isl] = -[ki... zeros(length(setdiff(1:nBus,ipq)))...]
    else
        J = [J11_tmp J12_tmp; J21_tmp J22_tmp]
        J[:, isl] = -ki
    end

    return J
end




function mismatch_vectorV(nBus, isl, ipv, ipq, v_angle, v_mag, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, p0, q0, dp, dq)
    # nBus, isl, ipv, ipq, v_angle, v_mag, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, p0, q0, dp, dq=nBus, isl, ipv, ipq, v_angle_pf, v_mag_pf, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, PloadBUS, QloadBUS, dp, dq
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
    return delta_mismatch_inj, p_inj, q_inj,p_calc,q_calc
end





function solution_voltage1(nBus, isl, ipq, ipv, v_angle, v_mag, pim_pf, delta_x)
    #nBus, isl, ipq, ipv, v_angle, v_mag, pim_pf, delta_x=nBus, isl, ipq, ipv, v_angle_pf, v_mag_pf, pim_pf, delta_x_pf
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


function line_flow1(L, V, del)
    fb = L[:, 1]
    tb = L[:, 2]
    g = L[:, 3]
    b = L[:, 4]
    bsh = L[:, 6]
    bsh.=0
    Imax = L[:, 7]

    y = g .+ im * b
    # V = vol_nodes_mag_pf
    # del = vol_nodes_theta_pf



    # sij_s = zeros(length(fb))
    # sij_r = zeros(length(fb))
    Iij = zeros(length(fb))
    percent = zeros(length(fb))


    for l in 1:length(fb)
        Vs = V[fb[l]] .* exp(im * del[fb[l]])
        Vr = V[tb[l]] .* exp(im * del[tb[l]])

        # sij_s[l] = abs(Vs * (Vs - Vr) * y[l])
        # sij_r[l] = abs(Vr * (Vs - Vr) * y[l])
        Iij[l] = abs((Vs - Vr) * y[l])
        # sij[l] = sqrt(abs(y[l])^2 * (V[fb[l]]^2+V[tb[l]]^2-(2*V[fb[l]]*V[tb[l]]*cos(del[fb[l]]-del[tb[l]]))))
        # sij[l] = maximum([sij_s[l], sij_r[l]])
        percent[l] = Iij[l] / Imax[l]
    end
    violated_lines = findall3(x -> x >= 1, percent)
    return Iij, violated_lines, percent
end


function print_multislack_acpf(vmag_i, vang_i, p_i, q_i, cg_i, df_i, itr_i)

    if cg_i == 1
        println(" Success!")

        XLSX.openxlsx("output_data\\R_ResultsE.xlsx", mode="w") do xf
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

        end

        println(" power flow has been convereged in $itr_i iterations")
        println(" Exactly, ( + network effect), the frequency will be changed :        $df_i pu")

    else
        println(" power flow has been diverged!")
    end
    println("[-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]")

end


function multislack_acpf_pglimits(y_bus,
    B, L, G, Ld, Shu,REN,BAT)

    # B, L, G, Ld, Shu,REN,BAT=deepcopy(NSysData["Bus"]), deepcopy(NSysData["Line"]), deepcopy(NSysData["Gen"]), deepcopy(NSysData["Load"]), deepcopy(NSysData["Shunt"]),deepcopy(NSysData["RES"]),deepcopy(NSysData["ESS"])

    G."Ramp"=convert.(Float64,G."Ramp")
    G."Ramp".=0.02; #Mousa selections: ask him for details
    Gvio=0;Gtemp = copy(G);eps_g=10e-8 # Min-Max Gen Limit
    PloadBUS,QloadBUS,p_calc,q_calc,q_inj,p_inj,vol_nodes_mag_pf, vol_nodes_theta_pf, convergence_indicator,Iij,percent,PGGG,QGGG,Jacob1=zeros(14)
    igen = G[:, 1]
    Sensivities=0
    Singularity=0
    while true
        # B=SysData["Bus"]; L=SysData["Line"]; G=SysData["Gen"]; Ld=SysData["Load"]; Shu=SysData["Shunt"];REN=SysData["RES"];BAT=SysData["ESS"];
        #----------------------------------------------------------------------------------------------
        # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 1: Settings--[-][-][-][-][-][-][-]---------
        #----------------------------------------------------------------------------------------------
        nBus = size(B, 1)
        PgenBUS = zeros(nBus)
        QgenBUS = zeros(nBus)
        PloadBUS = zeros(nBus)
        QloadBUS = zeros(nBus)
        QshuntBUS = zeros(nBus)
        QgenBUSmin = zeros(nBus)
        QgenBUSmax = zeros(nBus)
        PgenBUSmin = zeros(nBus)
        PgenBUSmax = zeros(nBus)
        dp= zeros(nBus, 3)
        dq= zeros(nBus, 3)


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
        ish = Shu[:, 1]
        slk = deepcopy(G[:, 12])

        global slk1=deepcopy(slk)

        PgenBUS[ipv] = G[:, 2]


        QgenBUSmax[ipv] = G[:, 4]
        QgenBUSmin[ipv] = G[:, 5]
        PgenBUSmax[ipv] = G[:, 7]
        PgenBUSmin[ipv] = G[:, 8]

        PloadBUS[iload] = Ld[:, 2]
        PloadBUS[REN[:,2]] -= REN[:, 5]
        PloadBUS[BAT[:,1]] -= BAT[:, 26]
        QloadBUS[iload] = Ld[:, 3]
        QloadBUS[REN[:,2]] -= REN[:, 6]
        QloadBUS[BAT[:,1]] -= BAT[:, 27]
        QshuntBUS[ish] = Shu[:, 2]

        dp[iload, 3].=1
        dq[iload, 3].=1
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


        v_mag_pf = ones(nBus)
        v_angle_pf = zeros(nBus)
        v_mag_pf[ipv] = G[:, 6]





        #.............................................................................
        # Note: In general, all (ipv=ipv+isl) buses can be considered as the slack bus
        #.............................................................................
        pim_pf = zeros(1)

        ki = zeros(npv + npq + npq, 1)
        (ki, beti) = spc_slack(ki, G)


        #.............................................................................


        vol_nodes_mag_pf = zeros(nBus)
        vol_nodes_theta_pf = zeros(nBus)
        max_mismatch_pf = 1
        iteration = 0
        mdf_idx = zeros(nBus)
        mdf_idxG = zeros(nBus)
        convergence_indicator = 0
        delta_mismatch_inj = zeros(nBus + npq)
        delta_mismatch=[]

        #---------------------------------------------------------------------------------------------------
        # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 2: Calcultae Scheduled Powers   [-][-][-][-][-][-][-]
        #---------------------------------------------------------------------------------------------------
        (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)

        #---------------------------------------------------------------------------------------------------


        #---------------------------------------------------------------------------------------------------
        # [-][-][-][-][-][-][-][-][-][-][-][-]  STEP 3: Main Loop of the Powr Flow--[-][-][-][-][-][-][-][-]
        #---------------------------------------------------------------------------------------------------
        Jacob1=0;p_inj;q_inj;p_calc;q_calc
        while max_mismatch_pf > epsilon && iteration <= itr_max
            # .......Step 3-1: Mismatch between  Scheduled and Calculated Powers............
            #...............................................................................
            (delta_mismatch_inj, p_inj, q_inj,p_calc,q_calc) = mismatch_vectorV(nBus, isl, ipv, ipq, v_angle_pf, v_mag_pf, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, PloadBUS, QloadBUS, dp, dq)



            # .......Step 3-2: Jacobian Matrix.............................................
            #..............................................................................
            Jacob1 = jacobianV(nBus, isl, ipq, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, ki, PloadBUS, QloadBUS, dp, dq)




            # .......Step 3-3: Inverse of Jacobian Matrix...................................
            #...............................................................................
            Jcb = lu(Jacob1, check=false)
            if !issuccess(Jcb)
                println("")
                show("Singularity observed for iteration= $iteration")
                println("")
                Singularity=1
                # @error("Stop")
                break
            end



            # .......Step 3-4: Update State Variables (V, delta,Pim)..........................
            #.................................................................................
            delta_x_pf = Jcb \ delta_mismatch_inj
            (v_angle_sv, v_mag_sv, pim_pf1) = solution_voltage1(nBus, isl, ipq, ipv, v_angle_pf, v_mag_pf, pim_pf, delta_x_pf)
            v_mag_pf = v_mag_sv                                                  # sv = solution voltage (Here, voltage is updated)
            v_angle_pf = v_angle_sv
            pim_pf = pim_pf1


            # .......Step 3-5: Check if PV-PQ Switching required ............................
            #.................................................................................
            max_mismatch_pf = maximum(abs.(delta_mismatch_inj))
            max_mismatch_pf_P = maximum(abs.(delta_mismatch_inj[1:nBus]))
            max_mismatch_pf_Q = maximum(abs.(delta_mismatch_inj[nBus+1:end]))
            if max_mismatch_pf_P < epsilon
                (mdf_idx, QgenBUS, ipv, ipq, G) = pv_pq_switching1(epsilon_suf, y_bus, v_mag_pf, v_angle_pf, ipv, ipq, mdf_idx, QgenBUS, QloadBUS, dq, QshuntBUS, QgenBUSmax, QgenBUSmin, G)

                if sum(mdf_idx) != 0
                    println("pq switching occured in iteration $iteration")
                    npv = length(ipv)
                    npq = length(ipq)
                    ki = zeros(npv + npq + npq, 1)
                    (ki, beti) = spc_slack(ki, G)
                    (p_sch_pf, q_sch_pf) = schedule_powerV(nBus,PgenBUS, QgenBUS, QshuntBUS)
                end
            end


            iteration = iteration + 1
        push!(delta_mismatch,maximum(delta_mismatch_inj'))
        end # End while loop

        #Arman
        # JacobArman[:,nBus+1:end]./repeat(v_mag_pf',2*nBus)
        # Jcb=lu(Jacob1, check=false)
        # Sens=Jcb\Matrix(1.0I,7,7)

        fullJacob = jacobianV(nBus, isl, ipq, v_angle_pf, v_mag_pf, y_bus, p_inj, q_inj, ki, PloadBUS, QloadBUS, dp, dq;flag=1)
        JacobArman=deepcopy(fullJacob)
        JacobArman[:,nBus+1:end]=fullJacob[:,nBus+1:end]./repeat(v_mag_pf',2*nBus)
        JcbArman = lu(JacobArman, check=false)
        Sensivities=JcbArman\Matrix(1.0I,2*nBus,2*nBus)
        Sensivities[isl,:].=0
        #Arman

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


        for i in 1:nBus
            vol_nodes_mag_pf[i] = v_mag_pf[i]
            vol_nodes_theta_pf[i] = (v_angle_pf[i])                         # Angles are in Radian Because input of Cosine and Sine functions are in radian (i.e., Julia consideres them as radian)
        end

        (delta_mismatch_inj, p_inj, q_inj,p_calc,q_calc) = mismatch_vectorV(nBus, isl, ipv, ipq, v_angle_pf, v_mag_pf, pim_pf, y_bus, p_sch_pf, q_sch_pf, ki, PloadBUS, QloadBUS, dp, dq)


        (Iij, violated_lines, percent) = line_flow1(L, vol_nodes_mag_pf, vol_nodes_theta_pf)



        df = -pim_pf/sum(beti)
        df=df[1]

        G[:, 2] = p_inj[igen] + PloadBUS[igen]

        vio_g1 = find_beyond(G[:, 2], G[:, 7], eps_g)
        vio_g2 = find_below(G[:, 2], G[:, 8], eps_g)
        vio_g = union(vio_g1, vio_g2)
        non_vio_g = setdiff(collect(1:size(G, 1)), vio_g)
        Gvio = length(vio_g)

        if Gvio == 0 || Singularity==1
            break
        else
            if length(G[:,12])==length(vio_g)
                println("All generators are violated")
                @warn("All generators are violated")
                break
            else
                G[vio_g1, 2] = G[vio_g1, 7]
                G[vio_g2, 2] = G[vio_g2, 8]
                G[non_vio_g, 2] = Gtemp[non_vio_g, 2]
                G[vio_g, 12] = zeros(length(vio_g))
                # println(" g_limit are considered")
            end
        end
    end
    G=Gtemp
    PGGG=(p_inj+PloadBUS)[findall(x -> x != 1, B[:, 2])]
    QGGG=(q_inj+QloadBUS)[findall(x -> x != 1, B[:, 2])]
    # PGGG=0
    # return vol_nodes_mag_pf, vol_nodes_theta_pf, delta_mismatch_inj, delta_mismatch, p_inj, q_inj, convergence_indicator, pim_pf, df, iteration,sij, violated_lines
    return p_inj,q_inj,vol_nodes_mag_pf, vol_nodes_theta_pf, convergence_indicator,Iij,percent,PGGG,QGGG,Sensivities

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