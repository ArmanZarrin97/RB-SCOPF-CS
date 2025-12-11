function PF_Valid_N(OP::Dict{String,Any}, Runcode) #Array{Union{Float64,Array{Float64}}}
    #OP=OP_init0
    #OP=Lininit_0_BV
    OP_Valid= deepcopy(OP) #Array{Union{Float64,Array{Float64}}}::Dict{String,Any} 
    nT = size(OP["PG"], 1)
    nBus = size(OP["V"], 2)
    TSysData = deepcopy(SysData::Dict{String,DataFrame})
    N_cost = zeros(nT)
    # OP_Valid["NSens"] = Array{Matrix{Float64},1}(undef, nT)
    OP_Valid["NJacob"]=[spzeros(Float64,2*nBus, 2*nBus) for _ in 1:nT]
    OP_Valid["Pinj"] = Array{Float64,2}(undef, nT, nBus)
    OP_Valid["Qinj"] = Array{Float64,2}(undef, nT, nBus)
    Vio_general = [DataFrame() for _ in 1:nT]
    convergence_indicator = 1
    OP_Valid=pre_prepration!(OP_Valid,TSysData;nT=nT)
    # OP_Valid_fix=deepcopy(OP_Valid)
    Threads.@threads for t = 1:nT #
        print(t)
        ASysData= Prepare_Dataset(OP_Valid, TSysData, t)
        #power injections before PF:
        # p_inj00 = zeros(size(ASysData["Bus"], 1))
        # q_inj00 = zeros(size(ASysData["Bus"], 1))
        # p_inj00[ASysData["Gen"]."Bus"] .= OP["PG"][t, :]
        # q_inj00[ASysData["Gen"]."Bus"] .= OP["QG"][t, :]
        # p_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."PG"
        # q_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."QG"
        # p_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Pd (pu)"
        # q_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Qd (pu)"
        # q_inj00[ASysData["Shunt"]."Bus"] .+= ASysData["Shunt"]."bsh0"
        NSysData = deepcopy(ASysData)
        V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob = multislack_acpf_order_pmaxqmax(NSysData["Bus"], NSysData["Line"], NSysData["Trans"], NSysData["Gen"], NSysData["Load"], NSysData["Shunt"], NSysData["RES"], NSysData["ESS"], [])
        # for i in [maximum(abs.(PGGG-OP_Valid["PG"][t,:])) maximum(abs.(QGGG-OP_Valid["QG"][t,:])) maximum(abs.(V-OP_Valid["V"][t,:])) maximum(abs.(θ-OP_Valid["θ"][t,:])) maximum(abs.(p_inj-p_inj00)) maximum(abs.(q_inj-q_inj00))]
        #     println(i)
        # end
        # Jcb=lu(Jacob)
        # Se1=Jcb \ I
        # round.(Se1, digits=4)

        # Jacob
        # isl
        # newJacob=Jacob[setdiff(1:end, isl),setdiff(1:end, isl)]
        # newJcb=lu(newJacob)
        # Se=newJcb \ I

        # sparse(round.(Se, digits=4))
        # maximum(abs.(PGGG-OP["PG"][t,:]))
        # maximum(abs.(QGGG-OP["QG"][t,:]))
        # maximum(abs.(V-OP["V"][t,:]))
        # maximum(abs.(θ-OP["θ"][t,:]))
        # maximum(abs.(p_inj-p_inj00))
        # maximum(abs.(q_inj-q_inj00))
        # maximum(abs.(v_mag-OP["V"][t,:]))
        # maximum(abs.(v_angle-OP["θ"][t,:]))
        #q_inj and QGGG are wrong
        if Linearization_TYPE[Runcode[2]] in ["SENS", "Hybrid"]
            isl=findall(x -> x == 3, NSysData["Bus"]."Bus_Type")
            nBus=size(NSysData["Bus"], 1)
            # OP_Valid["NSens"][t].=Partials(Jacob, V,isl,nBus)
            OP_Valid["NJacob"][t].=Jacob
        end

        if convergence_indicator != 1
            @error("PF not converged in normal operation and time step $t")
            break
        end
        Vio_general[t] = make_vio_general(V, percent, Iij)
        OP_Valid["V"][t, :] = V
        OP_Valid["θ"][t, :] = θ
        OP_Valid["PG"][t, :] = PGGG
        OP_Valid["QG"][t, :] = QGGG
        N_cost[t] = obj_cal(PGGG, OP_Valid, t)
        OP_Valid["Pinj"][t, :] = p_inj
        OP_Valid["Qinj"][t, :] = q_inj
        GC.gc()
    end
    Vio_Sum = make_vio_sum(Vio_general)

    return N_cost, Vio_Sum, OP_Valid, Vio_general, convergence_indicator
end
function PF_Valid_C(OP_NValid::Dict{String,Any}, Runcode) #Array{Union{Float64,Array{Float64}}}
    #OP_NValid=OP_init_NV
    OP_Valid = deepcopy(OP_NValid)
    nBus = size(OP_NValid["V_c_Line"], 2)
    nCon = size(OP_NValid["PG_c_Line"], 3)
    nT = size(OP_NValid["PG_c_Line"], 1)
    TSysData = deepcopy(SysData::Dict{String,DataFrame})
    C_cost = zeros(nT, nCon)
    # OP_Valid["CSens"] = Array{Matrix{Float64},2}(undef, nT, nCon)
    OP_Valid["CJacob"]=[spzeros(Float64,2*nBus, 2*nBus) for _ in 1:nT, _ in 1:nCon]
    OP_Valid["Pinj_c_Line"] = Array{Float64,3}(undef, nT, nBus, nCon)
    OP_Valid["Qinj_c_Line"] = Array{Float64,3}(undef, nT, nBus, nCon)
    Vio_general = [DataFrame() for _ in 1:nT, _ in 1:nCon]
    OP_Valid=pre_prepration!(OP_Valid,TSysData;nCon=nCon,nT=nT,C="ConL")
    Threads.@threads for t = 1:nT
        for i = 1:nCon
            print(t)
            ASysData= Prepare_Dataset(OP_Valid, TSysData, t; i=i, C="ConL")
            #power injections before PF:
            # p_inj00 = zeros(size(ASysData["Bus"], 1))
            # q_inj00 = zeros(size(ASysData["Bus"], 1))
            # p_inj00[ASysData["Gen"]."Bus"] .= OP_Valid["PG_c_Line"][t, :, i]
            # q_inj00[ASysData["Gen"]."Bus"] .= OP_Valid["QG_c_Line"][t, :, i]
            # p_inj00[ASysData["ESS"]."Bus"] .+= ASysData["ESS"]."PG" 
            # p_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."PG"
            # q_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."QG"
            # p_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Pd (pu)"
            # q_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Qd (pu)"
            # q_inj00[ASysData["Shunt"]."Bus"] .+= ASysData["Shunt"]."bsh0"
            CSysData = deepcopy(ASysData)
             V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob = multislack_acpf_order_pmaxqmax(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"], CSysData["RES"], CSysData["ESS"], [])
            # V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob, LNS_stage = cascading_sim(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[];printcascade=true);
            # V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob, LNS_stage, casrep = cascading_sim(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[]);
            # casrep=DataFrame()
            # casrep.Stage[1]=1
            # casrep.Jcbs[3]
            # casreps[t,i]=casrep
            # if sum(LNS_stage)!=0
            #     println("time '$(t)' and con '$i'")
            #     println("$(length(LNS_stage)) stage of Cascade cause $(sum(LNS_stage)) p.u. LNS ")
            # end
            if Linearization_TYPE[Runcode[2]] in ["SENS", "Hybrid"]
                isl=findall(x -> x == 3, CSysData["Bus"]."Bus_Type")
                nBus=size(CSysData["Bus"], 1)
                # OP_Valid["CSens"][t, i] = Partials(Jacob, V,isl,nBus)
                OP_Valid["CJacob"][t, i] .= Jacob
            end
            if convergence_indicator != 1
                @error("PF not converged for contingency number $i and time step $t")
            end
            percentall = zeros(size(TSysData["Line"], 1))
            Iijall = zeros(size(TSysData["Line"], 1))
            Iijall[ASysData["Line"]."Line_N"] = Iij
            percentall[ASysData["Line"]."Line_N"] = percent
            # for j=1:size(TSysData["Line"],1)
            #     if j<ASysData["ConL"]."#"[i] percentall[j]=percent[j] elseif j==TSysData["ConL"]."#"[i] percentall[j]=0 elseif j>TSysData["ConL"]."#"[i] percentall[j]=percent[j-1] end
            # end
            Vio_general[t, i] = make_vio_general(V, percentall, Iijall)
            OP_Valid["V_c_Line"][t, :, i] = V
            OP_Valid["θ_c_Line"][t, :, i] = θ
            OP_Valid["PG_c_Line"][t, :, i] = PGGG
            OP_Valid["QG_c_Line"][t, :, i] = QGGG
            C_cost[t, i] = obj_cal(PGGG, OP_Valid, t; i=i, C="ConL")
            # OP_Valid["CSens"][t, i] = Sens
            OP_Valid["Pinj_c_Line"][t, :, i] = p_inj
            OP_Valid["Qinj_c_Line"][t, :, i] = q_inj
        end
    end
    Vio_Sum = make_vio_sum(Vio_general)
    Avgcost = sum(C_cost, dims=2) ./ nCon
    return Avgcost, Vio_Sum, OP_Valid, Vio_general
end
function ST_PF_Valid_C(OP_NValid::Dict{String,Any}, Runcode) #Array{Union{Float64,Array{Float64}}}
    #OP_NValid=OP_init_NV
    OP_Valid = deepcopy(OP_NValid)
    nBus = size(OP_NValid["V_c_Line"], 2)
    nCon = size(OP_NValid["PG_c_Line"], 3)
    nT = size(OP_NValid["PG_c_Line"], 1)
    TSysData = deepcopy(SysData::Dict{String,DataFrame})
    C_cost = zeros(nT, nCon)
    # OP_Valid["CSens"] = Array{Matrix{Float64},2}(undef, nT, nCon)
    OP_Valid["CJacob"]=[spzeros(Float64,2*nBus, 2*nBus) for _ in 1:nT, _ in 1:nCon]
    OP_Valid["Pinj_c_Line"] = Array{Float64,3}(undef, nT, nBus, nCon)
    OP_Valid["Qinj_c_Line"] = Array{Float64,3}(undef, nT, nBus, nCon)
    Vio_general = [DataFrame() for _ in 1:nT, _ in 1:nCon]
    OP_Valid=pre_prepration!(OP_Valid,TSysData;nCon=nCon,nT=nT,C="ConL")
    Threads.@threads for t = 1:nT
        for i = 1:nCon
            print(t)
            ASysData= Prepare_Dataset(OP_Valid, TSysData, t; i=i, C="CASCADE")
            #power injections before PF:
            # p_inj00 = zeros(size(ASysData["Bus"], 1))
            # q_inj00 = zeros(size(ASysData["Bus"], 1))
            # p_inj00[ASysData["Gen"]."Bus"] .= OP_Valid["PG_c_Line"][t, :, i]
            # q_inj00[ASysData["Gen"]."Bus"] .= OP_Valid["QG_c_Line"][t, :, i]
            # p_inj00[ASysData["ESS"]."Bus"] .+= ASysData["ESS"]."PG" 
            # p_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."PG"
            # q_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."QG"
            # p_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Pd (pu)"
            # q_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Qd (pu)"
            # q_inj00[ASysData["Shunt"]."Bus"] .+= ASysData["Shunt"]."bsh0"
            CSysData = deepcopy(ASysData)
            V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob = multislack_acpf_order_pmaxqmax(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"], CSysData["RES"], CSysData["ESS"], [])
            # V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob, LNS_stage = cascading_sim(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[];printcascade=true);
            # V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob, LNS_stage, casrep = cascading_sim(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[]);
            # casrep=DataFrame()
            # casrep.Stage[1]=1
            # casrep.Jcbs[3]
            # casreps[t,i]=casrep
            # if sum(LNS_stage)!=0
            #     println("time '$(t)' and con '$i'")
            #     println("$(length(LNS_stage)) stage of Cascade cause $(sum(LNS_stage)) p.u. LNS ")
            # end
            if Linearization_TYPE[Runcode[2]] in ["SENS", "Hybrid"]
                isl=findall(x -> x == 3, CSysData["Bus"]."Bus_Type")
                nBus=size(CSysData["Bus"], 1)
                # OP_Valid["CSens"][t, i] = Partials(Jacob, V,isl,nBus)
                OP_Valid["CJacob"][t, i] .= Jacob
            end
            if convergence_indicator != 1
                @error("PF not converged for contingency number $i and time step $t")
            end
            percentall = zeros(size(TSysData["Line"], 1))
            Iijall = zeros(size(TSysData["Line"], 1))
            Iijall[ASysData["Line"]."Line_N"] = Iij
            percentall[ASysData["Line"]."Line_N"] = percent
            # for j=1:size(TSysData["Line"],1)
            #     if j<ASysData["ConL"]."#"[i] percentall[j]=percent[j] elseif j==TSysData["ConL"]."#"[i] percentall[j]=0 elseif j>TSysData["ConL"]."#"[i] percentall[j]=percent[j-1] end
            # end
            Vio_general[t, i] = make_vio_general(V, percentall, Iijall)
            OP_Valid["V_c_Line"][t, :, i] = V
            OP_Valid["θ_c_Line"][t, :, i] = θ
            OP_Valid["PG_c_Line"][t, :, i] = PGGG
            OP_Valid["QG_c_Line"][t, :, i] = QGGG
            C_cost[t, i] = obj_cal(PGGG, OP_Valid, t; i=i, C="ConL")
            # OP_Valid["CSens"][t, i] = Sens
            OP_Valid["Pinj_c_Line"][t, :, i] = p_inj
            OP_Valid["Qinj_c_Line"][t, :, i] = q_inj
        end
    end
    Vio_Sum = make_vio_sum(Vio_general)
    Avgcost = sum(C_cost, dims=2) ./ nCon
    return Avgcost, Vio_Sum, OP_Valid, Vio_general
end

function Cascade_Sim(OP_NValid::Dict{String,Any}, Runcode) #Array{Union{Float64,Array{Float64}}}
    #OP_NValid=OP_init_NV
    OP_Valid = deepcopy(OP_NValid)
    I_relax=1.037
    nBus = size(OP_NValid["V_c_Line"], 2)
    nCon = size(OP_NValid["PG_c_Line"], 3)
    nT = size(OP_NValid["PG_c_Line"], 1)
    TSysData = deepcopy(SysData::Dict{String,DataFrame})
    C_cost = zeros(nT, nCon)
    # OP_Valid["CSens"] = Array{Matrix{Float64},2}(undef, nT, nCon)
    OP_Valid["CJacob"]=[spzeros(Float64,2*nBus, 2*nBus) for _ in 1:nT, _ in 1:nCon]
    OP_Valid["Pinj_c_Line"] = Array{Float64,3}(undef, nT, nBus, nCon)
    OP_Valid["Qinj_c_Line"] = Array{Float64,3}(undef, nT, nBus, nCon)
    Vio_general = [DataFrame() for _ in 1:nT, _ in 1:nCon]
    casreps=[DataFrame() for _ in 1:nT, _ in 1:nCon]
    OP_Valid=pre_prepration!(OP_Valid,TSysData;nCon=nCon,nT=nT,C="ConL")
    Threads.@threads for t = 1:nT
        for i = 1:nCon
            print(t)
            ASysData= Prepare_Dataset(OP_Valid, TSysData, t; i=i, C="CASCADE", I_relax=I_relax)
            #power injections before PF:
            # p_inj00 = zeros(size(ASysData["Bus"], 1))
            # q_inj00 = zeros(size(ASysData["Bus"], 1))
            # p_inj00[ASysData["Gen"]."Bus"] .= OP_Valid["PG_c_Line"][t, :, i]
            # q_inj00[ASysData["Gen"]."Bus"] .= OP_Valid["QG_c_Line"][t, :, i]
            # p_inj00[ASysData["ESS"]."Bus"] .+= ASysData["ESS"]."PG" 
            # p_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."PG"
            # q_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."QG"
            # p_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Pd (pu)"
            # q_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Qd (pu)"
            # q_inj00[ASysData["Shunt"]."Bus"] .+= ASysData["Shunt"]."bsh0"
            CSysData = deepcopy(ASysData)
            # V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob = multislack_acpf_order_pmaxqmax(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"], CSysData["RES"], CSysData["ESS"], [])
            # V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob, LNS_stage = cascading_sim(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[];printcascade=true);
            V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob, LNS_stage, casrep = cascading_sim(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[]);
            
            @. casrep."VM" = casrep."VM".*I_relax
            # casrep=DataFrame()
            # casrep.Stage[1]=1
            # casrep.Jcbs[3]
            casreps[t,i]=casrep
            # if sum(LNS_stage)!=0
            #     println("time '$(t)' and con '$i'")
            #     println("$(length(LNS_stage)) stage of Cascade cause $(sum(LNS_stage)) p.u. LNS ")
            # end
            if Linearization_TYPE[Runcode[2]] in ["SENS", "Hybrid"]
                isl=findall(x -> x == 3, CSysData["Bus"]."Bus_Type")
                nBus=size(CSysData["Bus"], 1)
                # OP_Valid["CSens"][t, i] = Partials(Jacob, V,isl,nBus)
                OP_Valid["CJacob"][t, i] .= Jacob
            end
            if convergence_indicator != 1
                @error("PF not converged for contingency number $i and time step $t")
            end
            percentall = zeros(size(TSysData["Line"], 1))
            Iijall = zeros(size(TSysData["Line"], 1))
            Iijall[ASysData["Line"]."Line_N"] = Iij
            percentall[ASysData["Line"]."Line_N"] = percent
            # for j=1:size(TSysData["Line"],1)
            #     if j<ASysData["ConL"]."#"[i] percentall[j]=percent[j] elseif j==TSysData["ConL"]."#"[i] percentall[j]=0 elseif j>TSysData["ConL"]."#"[i] percentall[j]=percent[j-1] end
            # end
            Vio_general[t, i] = make_vio_general(V, percentall, Iijall)
            OP_Valid["V_c_Line"][t, :, i] = V
            OP_Valid["θ_c_Line"][t, :, i] = θ
            OP_Valid["PG_c_Line"][t, :, i] = PGGG
            OP_Valid["QG_c_Line"][t, :, i] = QGGG
            C_cost[t, i] = obj_cal(PGGG, OP_Valid, t; i=i, C="ConL")
            # OP_Valid["CSens"][t, i] = Sens
            OP_Valid["Pinj_c_Line"][t, :, i] = p_inj
            OP_Valid["Qinj_c_Line"][t, :, i] = q_inj
        end
    end
    Vio_Sum = make_vio_sum(Vio_general)
    Avgcost = sum(C_cost, dims=2) ./ nCon
    return Avgcost, Vio_Sum, OP_Valid, Vio_general,casreps
end

function PF_Valid_C2(OP_NValid::Dict{String,Any}, Runcode) #Array{Union{Float64,Array{Float64}}}
    #OP_NValid=OP_init_NV
    OP_Valid = deepcopy(OP_NValid)
    nBus = size(OP_NValid["V_c_Line"], 2)
    nCon = size(OP_NValid["PG_c_Line"], 3)
    nT = size(OP_NValid["PG_c_Line"], 1)
    TSysData = deepcopy(SysData::Dict{String,DataFrame})
    C_cost = zeros(nT, nCon)
    # OP_Valid["CSens"] = Array{Matrix{Float64},2}(undef, nT, nCon)
    OP_Valid["CJacob"]=[spzeros(Float64,2*nBus, 2*nBus) for _ in 1:nT, _ in 1:nCon]
    OP_Valid["Pinj_c_Line"] = Array{Float64,3}(undef, nT, nBus, nCon)
    OP_Valid["Qinj_c_Line"] = Array{Float64,3}(undef, nT, nBus, nCon)
    Vio_general = [DataFrame() for _ in 1:nT, _ in 1:nCon]
    casreps=[DataFrame() for _ in 1:nT, _ in 1:nCon]

    OP_Valid=pre_prepration!(OP_Valid,TSysData;nCon=nCon,nT=nT,C="ConL")
    Threads.@threads for t = 1:nT
        for i = 1:nCon
            print(t)
            ASysData= Prepare_Dataset(OP_Valid, TSysData, t; i=i, C="ConL")
            #power injections before PF:
            # p_inj00 = zeros(size(ASysData["Bus"], 1))
            # q_inj00 = zeros(size(ASysData["Bus"], 1))
            # p_inj00[ASysData["Gen"]."Bus"] .= OP_Valid["PG_c_Line"][t, :, i]
            # q_inj00[ASysData["Gen"]."Bus"] .= OP_Valid["QG_c_Line"][t, :, i]
            # p_inj00[ASysData["ESS"]."Bus"] .+= ASysData["ESS"]."PG" 
            # p_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."PG"
            # q_inj00[ASysData["RES"]."Bus"] += ASysData["RES"]."QG"
            # p_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Pd (pu)"
            # q_inj00[ASysData["Load"]."Bus_N"] .-= ASysData["Load"]."Qd (pu)"
            # q_inj00[ASysData["Shunt"]."Bus"] .+= ASysData["Shunt"]."bsh0"
            CSysData = deepcopy(ASysData)
            # V1, θ1, p_inj1, q_inj1, convergence_indicator1, Iij1, percent1, PGGG1, QGGG1, Jacob1 = multislack_acpf_order_pmaxqmax(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"], CSysData["RES"], CSysData["ESS"], [])
            # V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob, LNS_stage = cascading_sim(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[];printcascade=true);
            V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob, LNS_stage, casrep = cascading_sim(CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[]);
            # casrep=DataFrame()
            # casrep.Stage[1]=1
            # casrep.Jcbs[3]
            casreps[t,i]=casrep
            if sum(LNS_stage)!=0
                println("time '$(t)' and con '$i'")
                println("$(length(LNS_stage)) stage of Cascade cause $(sum(LNS_stage)) p.u. LNS ")
            end
            if Linearization_TYPE[Runcode[2]] in ["SENS", "Hybrid"]
                isl=findall(x -> x == 3, CSysData["Bus"]."Bus_Type")
                nBus=size(CSysData["Bus"], 1)
                # OP_Valid["CSens"][t, i] = Partials(Jacob, V,isl,nBus)
                OP_Valid["CJacob"][t, i] .= Jacob
            end
            if convergence_indicator != 1
                @error("PF not converged for contingency number $i and time step $t")
            end
            percentall = zeros(size(TSysData["Line"], 1))
            Iijall = zeros(size(TSysData["Line"], 1))
            Iijall[ASysData["Line"]."Line_N"] = Iij
            percentall[ASysData["Line"]."Line_N"] = percent
            # for j=1:size(TSysData["Line"],1)
            #     if j<ASysData["ConL"]."#"[i] percentall[j]=percent[j] elseif j==TSysData["ConL"]."#"[i] percentall[j]=0 elseif j>TSysData["ConL"]."#"[i] percentall[j]=percent[j-1] end
            # end
            Vio_general[t, i] = make_vio_general(V, percentall, Iijall)
            OP_Valid["V_c_Line"][t, :, i] = V
            OP_Valid["θ_c_Line"][t, :, i] = θ
            OP_Valid["PG_c_Line"][t, :, i] = PGGG
            OP_Valid["QG_c_Line"][t, :, i] = QGGG
            C_cost[t, i] = obj_cal(PGGG, OP_Valid, t; i=i, C="ConL")
            # OP_Valid["CSens"][t, i] = Sens
            OP_Valid["Pinj_c_Line"][t, :, i] = p_inj
            OP_Valid["Qinj_c_Line"][t, :, i] = q_inj
        end
    end
    Vio_Sum = make_vio_sum(Vio_general)
    Avgcost = sum(C_cost, dims=2) ./ nCon
    return Avgcost, Vio_Sum, OP_Valid, Vio_general,casreps
end

function make_vio_general(V::Vector{Float64}, percentall::Vector{Float64}, Iij::Vector{Float64})
    # # SysData["Line"][221,:]
    # Iij[221]
    VS = DataFrame(Type=["Vmax"; "Vmin"; "Smax"], N=[0; 0; 0], max=[0.0; 0.0; 0.0], mean=[0.0; 0.0; 0.0], Element=[Int64[], Int64[], Int64[]], percent=[Float64[], Float64[], Float64[]], abs=[Float64[], Float64[], Float64[]])
    Vio_type = [any(V .> SysData["Bus"]."V_max (pu)" .+ ϵ) any(V .< SysData["Bus"]."V_min (pu)" .- ϵ) any(percentall .> 1 + ϵ)]
    if any(Vio_type)
        if Vio_type[1]
            VS[1, "Element"] = findall(V .> SysData["Bus"]."V_max (pu)" .+ ϵ)
            VS[1, "N"] = size(VS[1, "Element"], 1)
            VS[1, "abs"] = V[VS[1, "Element"]] - SysData["Bus"]."V_max (pu)"[VS[1, "Element"]]
            VS[1, "percent"] = 100 * ((V[VS[1, "Element"]] - SysData["Bus"]."V_max (pu)"[VS[1, "Element"]]) ./ SysData["Bus"]."V_max (pu)"[VS[1, "Element"]])
            VS[1, "max"] = maximum(VS[1, "abs"])
            VS[1, "mean"] = sum(VS[1, "abs"]; init=0) / VS[1, "N"]
        end
        if Vio_type[2]
            VS[2, "Element"] = findall(V .< SysData["Bus"]."V_min (pu)" .- ϵ)
            VS[2, "N"] = size(VS[2, "Element"], 1)
            VS[2, "abs"] = SysData["Bus"]."V_min (pu)"[VS[2, "Element"]] - V[VS[2, "Element"]]
            VS[2, "percent"] = 100 * ((SysData["Bus"]."V_min (pu)"[VS[2, "Element"]] - V[VS[2, "Element"]]) ./ SysData["Bus"]."V_min (pu)"[VS[2, "Element"]])
            VS[2, "max"] = maximum(VS[2, "abs"])
            VS[2, "mean"] = sum(VS[2, "abs"]; init=0) / VS[2, "N"]
        end
        if Vio_type[3]
            VS[3, "Element"] = findall(percentall .> 1 + ϵ)
            VS[3, "N"] = size(VS[3, "Element"], 1)
            VS[3, "abs"] = Iij[VS[3, "Element"]] - SysData["Line"]."RATE_A"[VS[3, "Element"]]
            VS[3, "percent"] = 100 * (percentall[VS[3, "Element"]] .- 1)
            VS[3, "max"] = maximum(VS[3, "abs"])
            VS[3, "mean"] = sum(VS[3, "abs"]; init=0) / VS[3, "N"]
        end
    end
    return VS
end
function make_vio_sum(Vio_general::Array{DataFrame})
    nCon = size(Vio_general, 2)
    nT = size(Vio_general, 1)
    Vio_name = ["Vmax"; "Vmin"; "Smax"]
    Nins = 0
    NVio = 0
    Reg = 0
    max_ = 0
    type = 0
    for t = 1:nT, i = 1:nCon
        if sum(Vio_general[t, i][:, "N"]) > 0
            Nins += 1
            for j = 1:3
                Reg += Vio_general[t, i][j, "N"] * Vio_general[t, i][j, "mean"]
                NVio += Vio_general[t, i][j, "N"]
                if Vio_general[t, i][j, "max"] > max_
                    max_ = Vio_general[t, i][j, "max"]
                    type = Vio_name[j] * "-t" * string(t) * "-c" * string(i)
                end
            end
        end
    end
    Vio_Sum = DataFrame(Nins=Nins, NVio=NVio, Mean=Reg / NVio, Max=max_, Type=type)
    return Vio_Sum
end
function obj_cal(PG, OP, t; i=0, C="N")
    #PG=PGGG;OP=OP_Valid;

    costG = SysData["Gen"]."cost_c" + SysData["Gen"]."cost_b" * SysData["Base_MVA"]."base_MVA"[1] .* PG #+ SysData["Gen"]."cost_a"*SysData["Base_MVA"]."base_MVA"[1]^2 .*(PG.^2)

    if C == "N"
        if !isempty(OP["Pdis"])
            PBat = (OP["Pdis"]+OP["Pch"])[t, :]
            costESS = (SysData["ESS"]."cost_a" .* PBat .^ 2 + SysData["ESS"]."cost_b" .* PBat + SysData["ESS"]."cost_c") * SysData["Base_MVA"]."base_MVA"[1]
        else
            costESS = 0
        end

        if !isempty(OP["Pinc"])
            PFL = (OP["Pinc"]+OP["Pdec"])[t, :]
            costFL = SysData["FL"]."ShCost" * SysData["Base_MVA"]."base_MVA"[1] .* PFL
        else
            costFL = 0
        end

        if !isempty(OP["Pcurt"])
            costcurt = sum(OP["Pcurt"][t, :]) * SysData["Base_MVA"]."Rccost"[1] * SysData["Base_MVA"]."base_MVA"[1]
        else
            costcurt = 0
        end

        costshed = sum(OP["Pshed"][t, :]) * SysData["Base_MVA"]."Lscost"[1] * SysData["Base_MVA"]."base_MVA"[1]
        costcurtc = sum(OP["Pcurtc"][t, :]) * SysData["Base_MVA"]."Lscost"[1] * SysData["Base_MVA"]."base_MVA"[1]

        # Lccost=raw["LScost"]*raw["sbase"]
        # Rccost=raw["RCcost"]*raw["sbase"]

        # FactorT=max(1,500*(t0-132))
        # λG=0.005*FactorT*raw["Gen"]."cost_b"./(2)
        # λS=0.005*FactorT*((144)^2)*raw["ESS"]."cost_b"./(2*TD)
        # λF=0.005*FactorT*((144)^2)*raw["Load"]."ShCost"[raw["FL"]."Bus"]./(2*TD)

        # if WP==1
        #     @expression(SCOPF,cost_Pen, 0)
        # elseif WP=="SP"
        #     @expression(SCOPF,FPG[i=1:raw["nGen"]], λG[i]*((PG[nTP,i])^2))
        #     # @expression(SCOPF,FPS[i=1:raw["nESS"]], λS[i]*((Init_Data["DA_SCOPF"]["soc"][div(t0-1,6)+1,i+1]-soc[nTP,i])/(1/TD -mod(t0-1,6)))^2)
        #     # @expression(SCOPF,FPF[i=1:raw["nFL"]], λF[i]*((Init_Data["DA_SCOPF"]["UE"][div(t0-1,6)+1,i+1]-UE[nTP,i])/(1/TD -mod(t0-1,6)))^2)
        #     @expression(SCOPF,FPS[i=1:raw["nESS"]], λS[i]*(soc[nTP,i]-initraw["soc"][nTP,i])^2)
        #     @expression(SCOPF,FPF[i=1:raw["nFL"]], λF[i]*(UE[nTP,i]-initraw["UE"][nTP,i])^2)
        #     @expression(SCOPF,cost_Pen, sum(FPG)+sum(FPS)+sum(FPF))
        #     @constraint(SCOPF,socmaxd[i=1:raw["nESS"],t=1:nTP], Init_Data["DA_SCOPF"]["soc"][nTL,i+1]-soc[nTP,i]<= (nTL*nTS-t0)* raw["ESS"]."Charge_rating (MW)"[i]/raw["sbase"]*TD)
        #     @constraint(SCOPF,UEmaxd[i=1:raw["nFL"],t=1:nTP],  Init_Data["DA_SCOPF"]["UE"][nTL,i+1]-UE[nTP,i]<= (nTL*nTS-t0)*Inc_coef*Proffin["PL"][t,raw["FL"]."Bus"[i]]*TD)
        # else
        #     @expression(SCOPF,FPG[i=1:raw["nGen"]], λG[i]*((PG[nTP,i])^2))
        #     @expression(SCOPF,FPS[i=1:raw["nESS"]], λS[i]*(soc[nTP,i]-initraw["soc"][nTP,i])^2)
        #     @expression(SCOPF,FPF[i=1:raw["nFL"]], λF[i]*(UE[nTP,i]-initraw["UE"][nTP,i])^2)
        #     tempvar=sum(FPG;init=0)+sum(FPS;init=0)+sum(FPF;init=0)
        #     @expression(SCOPF,cost_Pen, tempvar)
        # end
    elseif C == "ConL"
        if !isempty(OP["Pdis_c_Line"])
            PBat = (OP["Pdis_c_Line"]+OP["Pch_c_Line"])[t, :, i]
            costESS = (SysData["ESS"]."cost_a" .* PBat .^ 2 + SysData["ESS"]."cost_b" .* PBat + SysData["ESS"]."cost_c") * SysData["Base_MVA"]."base_MVA"[1]
        else
            costESS = 0
        end

        if !isempty(OP["Pinc_c_Line"])
            PFL = (OP["Pinc_c_Line"]+OP["Pdec_c_Line"])[t, :, i]
            costFL = SysData["FL"]."ShCost" * SysData["Base_MVA"]."base_MVA"[1] .* PFL
        else
            costFL = 0
        end

        if !isempty(OP["Pcurt_c_Line"])
            costcurt = sum(OP["Pcurt_c_Line"][t, :, i]) * SysData["Base_MVA"]."Rccost"[1] * SysData["Base_MVA"]."base_MVA"[1]
        else
            costcurt = 0
        end

        costshed = sum(OP["Pshed_c_Line"][t, :, i]; init=0) * SysData["Base_MVA"]."Lscost"[1] * SysData["Base_MVA"]."base_MVA"[1]
        costcurtc = sum(OP["Pcurtc_c_Line"][t, :, i]; init=0) * SysData["Base_MVA"]."Lscost"[1] * SysData["Base_MVA"]."base_MVA"[1]
    end

    return sum(costG; init=0) + sum(costESS; init=0) + sum(costFL; init=0) + costshed + costcurt + costcurtc
end
function AC_PFCheck(OP, Runcode)
    # OP=OP_init0
    casreps=[]
    if !isempty(OP["Pinc"])
        OP = fix_FL!(OP) 
    end
    Ncost, Vio_Sum_N, OP_init_NV, Vio_general_N =PF_Valid_N(OP, Runcode)
    GC.gc()
    if Problem_Type in ["SCOPF","RBOPF"]
        # AVG_C_cost, Vio_Sum_C, OP_initV, Vio_general_C = ST_PF_Valid_C(OP_init_NV, Runcode)
        AVG_C_cost, Vio_Sum_C, OP_initV, Vio_general_C = PF_Valid_C(OP_init_NV, Runcode)
        # AVG_C_cost, Vio_Sum_C, OP_initV, Vio_general_C, casreps = PF_Valid_C2(OP_init_NV, Runcode)
        _, _, _, _, casreps = Cascade_Sim(OP_initV, Runcode)
        # Cascade_Sim(OP_init_NV, Runcode);
        if !isempty(findall(x-> sum(x."LNS_T")>0.0, casreps))
            @warn "There is new Short-term security problem"
        end
        GC.gc()
    elseif Problem_Type == "OPF"
        OP_initV = OP_init_NV
        Vio_Sum_C = DataFrame(Nins=0, NVio=0, Mean=0, Max=0, Type="")
        AVG_C_cost = Ncost
        Vio_general_C = 0
        casreps=[]
    end
    Vio_general_N, Vio_general_C = coupling_violation_check!(OP_initV, Vio_general_N, Vio_general_C)
    return Ncost, Vio_Sum_N, Vio_general_N, AVG_C_cost, Vio_Sum_C, OP_initV, Vio_general_C, casreps
end


function coupling_violation_check!(OP_initV, Vio_general_N, Vio_general_C)

    ramp = abs.(OP_initV["PG"][2:end, :] - OP_initV["PG"][1:end-1, :])
    Act_RampN = [zeros(size(ramp, 2))'; ramp]
    if Problem_Type == "SCOPF"
        Act_RampC = abs.(OP_initV["PG_c_Line"] .- OP_initV["PG"])
    end

    for t = 1:nTL
        if maximum(Act_RampN[t, :] .- (SysData["Gen"]."Ramp")') < ϵ
            new_row = ["RampN", 0, 0.0, 0.0, Int64[], Float64[], Float64[]]
            push!(Vio_general_N[t], new_row)
        else
            maximum(Act_RampN[t, :] .- SysData["Gen"]."Ramp") > ϵ
            mag = (Act_RampN[t, :] .- SysData["Gen"]."Ramp")
            percent = 100 * (mag ./ (SysData["Gen"]."Ramp"))
            elements = findall(x -> x > ϵ, mag)
            average = size(elements, 1) == 0 ? 0.0 : sum(mag[elements]) / size(elements, 1)
            new_row = ["RampN", size(elements, 1), maximum(mag[elements], init=0), average, elements, percent[elements], mag[elements]]
            push!(Vio_general_N[t], new_row)
        end
        if Problem_Type == "SCOPF"
            for c = 1:size(SysData["ConL"], 1)
                if maximum(Act_RampC[t, :, c] .- SysData["Gen"]."Ramp" * CCT) < ϵ
                    new_row = ["RampC", 0, 0.0, 0.0, Int64[], Float64[], Float64[]]
                    push!(Vio_general_C[t, c], new_row)
                else
                    maximum(Act_RampC[t, :, c] .- SysData["Gen"]."Ramp" * CCT) > ϵ
                    mag = (Act_RampC[t, :, c] .- SysData["Gen"]."Ramp" * CCT)
                    percent = 100 * (mag ./ (SysData["Gen"]."Ramp" * CCT))
                    elements = findall(x -> x > ϵ, mag)
                    average = size(elements, 1) == 0 ? 0.0 : sum(mag[elements]) / size(elements, 1)
                    new_row = ["RampC", size(elements, 1), maximum(mag[elements], init=0), average, elements, percent[elements], mag[elements]]
                    push!(Vio_general_C[t, c], new_row)
                end
            end
        end
    end
    return Vio_general_N, Vio_general_C
end
# function find_min_Trust_Region(Vio,raw,SEN,Lininit,t;c=1)
#     model=Model(CPLEX.Optimizer)
#     MOI.set(model, MOI.Silent(), true)
#     set_optimizer_attribute(model, CPLEX.PassNames(), true)
#     @variable(model,ΔPinjinc[i = 1:raw["nBus"]]>=0)
#     @variable(model,ΔQinjinc[i = 1:raw["nBus"]]>=0)
#     @variable(model,ΔPinjdec[i = 1:raw["nBus"]]>=0)
#     @variable(model,ΔQinjdec[i = 1:raw["nBus"]]>=0)

#     @expression(model,ΔPinj[i = 1:raw["nBus"]],ΔPinjinc[i]-ΔPinjdec[i])
#     @expression(model,ΔQinj[i = 1:raw["nBus"]],ΔQinjinc[i]-ΔQinjdec[i])

#     @expression(model, ΔSmagexp[n=1:raw["nLine"]], sum(SEN["S_IP"][n,i]*(ΔPinj[i])+ SEN["S_IQ"][n,i]*(ΔQinj[i]) for i=1:raw["nBus"]))
#     for n=1:raw["nLine"]
#         minthreshold=0.01*(raw["Line"]."RATE_A"[n]^2)/(TRfree+1)
#         map_coefficients_inplace!.(x->ifelse.(abs(x) .> minthreshold, x, 0.0),ΔSmagexp[n])
#     end
#     drop_zeros!.(ΔSmagexp)

#     @expression(model, ΔVmagexp[n=1:raw["nBus"]], sum(SEN["S_VP"][n,i]*(ΔPinj[i])+ SEN["S_VQ"][n,i]*(ΔQinj[i]) for i=1:raw["nBus"]))
#     minthreshold=0.001*raw["Bus"]."V_min (pu)"[1]/(0.1*(TRfree+1))
#     map_coefficients_inplace!.(x->ifelse.(abs(x) .> minthreshold, x, 0.0),ΔVmagexp)
#     drop_zeros!.(ΔVmagexp)


#     Vmaxviodict=Dict(Vio[1,"Element"] .=> Vio[1,"abs"])
#     @constraint(model,Vmaxvio[i in keys(Vmaxviodict)], ΔVmagexp[i]<=-Vmaxviodict[i])

#     Vminviodict=Dict(Vio[2,"Element"] .=> Vio[2,"abs"])
#     @constraint(model,Vminvio[i in keys(Vminviodict)], ΔVmagexp[i]>=Vminviodict[i])

#     Smaxviodict=Dict(Vio[3,"Element"] .=> Vio[3,"abs"])
#     @constraint(model,Smaxvio[i in keys(Smaxviodict)], ΔSmagexp[i]<=-Smaxviodict[i])

#     maxPupadj=zeros(raw["nBus"]);maxPdwadj=zeros(raw["nBus"]);maxQupadj=zeros(raw["nBus"]);maxQdwadj=zeros(raw["nBus"])

#     for i in 1:raw["nBus"]
#         maxPupadj[i]=sum(raw["Gen"]."P_max (pu)"[j]-Lininit["PG"][t,j,c] for j in 1:raw["nGen"] if raw["Gen"]."Bus"[j]==i;init=0)
#         maxPdwadj[i]=sum(Lininit["PG"][t,j,c]-raw["Gen"]."P_min (pu)"[j] for j in 1:raw["nGen"] if raw["Gen"]."Bus"[j]==i;init=0)
#         maxQupadj[i]=sum(raw["Gen"]."Q_max (pu)"[j]-Lininit["QG"][t,j,c] for j in 1:raw["nGen"] if raw["Gen"]."Bus"[j]==i;init=0)
#         maxQdwadj[i]=sum(Lininit["QG"][t,j,c]-raw["Gen"]."Q_min (pu)"[j] for j in 1:raw["nGen"] if raw["Gen"]."Bus"[j]==i;init=0)
#     end

#     @constraint(model,Genlimit1[i = 1:raw["nBus"]],ΔPinjinc[i]<=maxPupadj[i])
#     @constraint(model,Genlimit2[i = 1:raw["nBus"]],ΔPinjdec[i]<=maxPdwadj[i])
#     @constraint(model,Genlimit3[i = 1:raw["nBus"]],ΔQinjinc[i]<=maxQupadj[i])
#     @constraint(model,Genlimit4[i = 1:raw["nBus"]],ΔQinjdec[i]<=maxQdwadj[i])

#     @constraint(model,sum(ΔPinj)<=min(0.01*(TRfree),1))
#     @constraint(model,-min(0.01*(TRfree),1)<=sum(ΔPinj))
#     @constraint(model,sum(ΔQinj)<=min(0.01*(TRfree),1))
#     @constraint(model,-min(0.01*(TRfree),1)<=sum(ΔQinj))

#     @objective(model,Min, sum(ΔPinjdec.+ΔPinjinc.+ΔQinjdec.+ΔQinjinc))
#     optimize!(model)
#     return DataFrame(ΔPinjinc=value.(ΔPinjinc),ΔPinjdec=value.(ΔPinjdec),ΔQinjinc=value.(ΔQinjinc),ΔQinjdec=value.(ΔQinjdec))
# end