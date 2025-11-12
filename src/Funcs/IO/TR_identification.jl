function findTR_Q_coupled5(V_0,Vmax_VP,Vmin_VP,VmaxIndex_P,VminIndex_P,Vmaxlimit0,Vminlimit0,S_VQ,QBound,raw,cl;Improve_Rate=5)
    
    z0=sum((V_0-Vmaxlimit0).>=0)+sum(Vminlimit0-V_0.>=0)
    zmax=round(size(Vmaxlimit0)[1]*size(Vmaxlimit0)[2]*size(Vmaxlimit0)[3]*Q_con_percent)+1
    TR=0;TR_max=0.5;TR_min=0.0
    iter=0
    Itermax=10
    Vmaxinrisk=falses(size(V_0));Vmininrisk=falses(size(V_0));
    Vmaxinrisk_P=falses(size(V_0));Vmininrisk_P=falses(size(V_0));
    # Vmaxinrisk_Q=falses(size(V_0));Vmaxinrisk_Q=falses(size(V_0));
    Vmaxviolated=0
    Vminviolated=0
    # worst_vio=max(maximum(Vmax_VP-Vmaxlimit0),maximum(-(Vmin_VP-Vminlimit0)))
    worst_vio=max(maximum(V_0-Vmaxlimit0),maximum(-(V_0-Vminlimit0)))

    # target_acc=(worst_vio/Improve_Rate)
    adj_acc=0
    acc=0
    best_acc=deepcopy(worst_vio)
    target_acc=0
    Vmaxlimit_new=deepcopy(Vmaxlimit0)
    genbus=raw["Gen"].Bus
    genbus_s = Set(genbus)  # O(1) lookup
    allbus=1:raw["nBus"]
    allbus_s = Set(allbus)  # O(1) lookup

    V0_length= length(V_0)
    nBus      = raw["nBus"]
    nQActBus  = raw["nQActBus"]
    ovr_Sens_VQ_up = zeros(nTP, raw["nBus"], cl);ovr_Sens_VQ_dn = zeros(nTP, raw["nBus"], cl)
    ovr_Sens_VQ_up_org = zeros(nTP, raw["nBus"], cl);ovr_Sens_VQ_dn_org = zeros(nTP, raw["nBus"], cl)

    ovr_Sens_VQ_upmat = zeros(Float64, nBus, nQActBus, nTP, cl)
    ovr_Sens_VQ_dnmat = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_up             = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_dn             = zeros(Float64, nBus, nQActBus, nTP, cl)
    # upelem            = zeros(Float64, nBus, nQActBus, nTP, cl)
    # dnelem            = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_up1      =falses(nBus, nQActBus, nTP, cl)
    VQ_dn1      =falses(nBus, nQActBus, nTP, cl)
    V_max_VPQ   = similar(ovr_Sens_VQ_up);V_min_VPQ   = similar(ovr_Sens_VQ_dn)
    V_max_VQ   = similar(ovr_Sens_VQ_up);V_min_VQ   = similar(ovr_Sens_VQ_dn)
    vals_VP=zeros(Float64,2* V0_length);vals_VQ=zeros(Float64,2* V0_length)

    DV_up=Vmax_VP.- Vmaxlimit0
    DV_dn=Vmin_VP.- Vminlimit0

    DV_up_4D = zeros(Float64, nBus, 1, nTP, cl)
    DV_dn_4D = zeros(Float64, nBus, 1, nTP, cl)
    S_VQ_4D = zeros(Float64, nBus, nQActBus, nTP, cl)
    Threads.@threads for t in 1:nTP
        for c in 1:cl
            @views S_VQ_4D[:, :, t, c] .= S_VQ[t, c]
            @views DV_up_4D[:, 1, t, c] .= DV_up[t, :, c]  # direct assignment
            @views DV_dn_4D[:, 1, t, c] .= DV_dn[t, :, c]  # direct assignment
        end
    end
    S_VQ1=keep_K_percent_twoside_4D(S_VQ_4D,0.9)
    Qinj_up,Qinj_dn=0,0
    y=0
    z=0
    TR=0.01
    Vmaxlimit_new=Vmaxlimit0
    Vminlimit_new=Vminlimit0
    report_TR= DataFrame(worst_vio=0.0,acc=0.0,TR=0.0,z=z0,y=z0,z_acc_Q=0.0,y_acc_Q=0.0,z_acc_P=0.0,y_acc_P=0.0,plotdata=[0.0])
    deleteat!(report_TR,1)
    z_acc_P=0;z_acc_Q=0
    y_acc_P=0;y_acc_Q=0
    S_VQ_pos = max.(S_VQ_4D, 0.0)
    S_VQ_neg = min.(S_VQ_4D, 0.0)

    TR_best=0.0
    while TR_max - TR_min > 1e-7  # Stop when the range is small enough
        VQ_up1.=0
        VQ_dn1.=0
        # Update TR using bisection
        if iter==0
            nothing
        elseif iter>=Itermax
            TR=TR_best
        else
            TR = (TR_min + TR_max) / 2
        end
        Qinj_up=min.(TR,QBound[1])
        Qinj_dn=max.(-TR,QBound[2])
        Threads.@threads for t in 1:nTP
            for c in 1:cl
            # @views S_VQ_tc = S_VQ[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
            @views up = Qinj_up[t, :, c]  # 1 × nPActBus
            @views dn = Qinj_dn[t, :, c]

            @views @. ovr_Sens_VQ_upmat[:, :, t, c] = S_VQ_pos[:,:,t, c]* up' + S_VQ_neg[:,:,t, c]* dn'
            @views @. ovr_Sens_VQ_dnmat[:, :, t, c] = S_VQ_pos[:,:,t, c]* dn' + S_VQ_neg[:,:,t, c]* up'


            @views ovr_Sens_VQ_up[t, :, c] .= sum(ovr_Sens_VQ_upmat[:, :, t, c], dims=2)
            @views ovr_Sens_VQ_dn[t, :, c] .= sum(ovr_Sens_VQ_dnmat[:, :, t, c], dims=2)
            # @views ovr_Sens_VQ_upmat[:, :, t, c] .= up_max_V
            # @views ovr_Sens_VQ_dnmat[:, :, t, c] .= dn_max_V
            end
        end
        @views @. ovr_Sens_VQ_up_org = ovr_Sens_VQ_up + 1e-6
        @views @. ovr_Sens_VQ_dn_org = ovr_Sens_VQ_dn + 1e-6


        @views @. VQ_up = ovr_Sens_VQ_upmat + DV_up_4D
        @views @. VQ_dn = -ovr_Sens_VQ_dnmat - DV_dn_4D
        keep_K_percent_new_4D_V2!(VQ_up1,VQ_up,0.0)
        keep_K_percent_new_4D_V2!(VQ_dn1,VQ_dn,0.0)

        keep_top_4D_V2!(VQ_up1,VQ_up,10)
        keep_top_4D_V2!(VQ_dn1,VQ_dn,10)
        @. VQ_up1 |=S_VQ1
        @. VQ_dn1 |=S_VQ1
        @. ovr_Sens_VQ_upmat *= VQ_up1 
        @. ovr_Sens_VQ_dnmat *= VQ_dn1

        @views ovr_Sens_VQ_up.= permutedims(sum(ovr_Sens_VQ_upmat,dims=2)[:,1,:,:], (2, 1, 3)) 
        @views ovr_Sens_VQ_dn.= permutedims(sum(ovr_Sens_VQ_dnmat,dims=2)[:,1,:,:], (2, 1, 3)) 
        @. V_max_VPQ=(Vmax_VP+ovr_Sens_VQ_up/1 .-Vmaxlimit_new)
        @. V_max_VQ=(V_0+ovr_Sens_VQ_up/1 .-Vmaxlimit_new)
        @. V_max_VQ[VmaxIndex_P].=-10
        @. V_max_VPQ[VmaxIndex_P].=-10

        @. V_min_VPQ=-(Vmin_VP+ovr_Sens_VQ_dn/1 .-Vminlimit_new)
        @. V_min_VQ=-(V_0+ovr_Sens_VQ_dn/1 .-Vminlimit_new)
        @. V_min_VQ[VminIndex_P].=-10
        @. V_min_VPQ[VminIndex_P].=-10
        @views @. vals_VP[1:V0_length] = V_max_VPQ[:]
        @views @. vals_VP[1+V0_length:end] = V_min_VPQ[:]
        k_VP = min(length(vals_VP), Int(zmax))
        z_acc_P = partialsort(vals_VP, 1:k_VP; rev=true)[end]
        vals_VQ = vcat( vec(V_max_VQ),vec((V_min_VQ)))

        @views @. vals_VQ[1:V0_length] = V_max_VQ[:]
        @views @. vals_VQ[1+V0_length:end] = V_min_VQ[:]

        k_VQ = min(length(vals_VQ), Int(zmax))
        z_acc_Q = partialsort(vals_VQ, 1:k_VQ; rev=true)[end]
        y_acc_Q=max(maximum(V_0+(ovr_Sens_VQ_dn_org/Q_eff).-Vmaxlimit_new),maximum(-(V_0+(ovr_Sens_VQ_up_org/Q_eff).-Vminlimit_new)),0)

        y_acc_P=max(maximum(Vmax_VP+(ovr_Sens_VQ_dn_org/Q_eff).-Vmaxlimit_new),maximum(-(Vmin_VP+(ovr_Sens_VQ_up_org/Q_eff).-Vminlimit_new)),0)

        acc=max(-0.0,y_acc_Q)
        @. Vmaxinrisk_P = Vmax_VP+ovr_Sens_VQ_up/1 >= Vmaxlimit_new+max(z_acc_P,1e-4)
        @. Vmininrisk_P = Vmin_VP+ovr_Sens_VQ_dn/1 <= Vminlimit_new-max(z_acc_P,1e-4)

        # @. Vmaxinrisk_Q = V_0+ovr_Sens_VQ_up/1 .>= Vmaxlimit_new.+max(z_acc_Q,1e-4)
        # @. Vmininrisk_Q = V_0+ovr_Sens_VQ_dn/1 .<= Vminlimit_new.-max(z_acc_Q,1e-4)

        @. Vmaxinrisk= Vmaxinrisk_P #.|| Vmaxinrisk_Q
        @. Vmininrisk= Vmininrisk_P #.|| Vmininrisk_Q
        z= sum(Vmaxinrisk) + sum(Vmininrisk)
        y = sum(Vmax_VP+(ovr_Sens_VQ_dn/Q_eff) .> Vmaxlimit_new.+y_acc_P) + sum(Vmin_VP+(ovr_Sens_VQ_up/Q_eff) .<Vminlimit_new.-y_acc_P)
        push!(report_TR,[worst_vio,acc,TR,z,y,z_acc_Q,y_acc_Q,z_acc_P,y_acc_P,0.0])
        if iter>=Itermax
            break
        end

        if acc<=best_acc && z_acc_Q <= y_acc_Q && z_acc_P<=worst_vio  #&& max(z_acc_Q,z_acc_P) <= max(0.05,y_acc_Q)
            best_acc=acc
            TR_best=TR
        end
        if z_acc_P>=worst_vio  || ((y_acc_P<=max(0.001,maximum(Vmax_VP),maximum(-Vmin_VP))-0.001) && max(0,z_acc_Q)>=y_acc_Q) #|| z_acc_Q >= worst_vio #|| max(z_acc_P,z_acc_Q)>= max(0.05,y_acc_Q)
            TR_max = TR  # Reduce TR
        else
            TR_min = TR  # Increase TR
        end
        iter+=1
    end
    VmaxIndex=findall(Vmaxinrisk)
    VminIndex=findall(Vmininrisk)
    display(report_TR)
    return TR,VmaxIndex,VminIndex,max(0.0,z_acc_Q),[Qinj_up,Qinj_dn],worst_vio,report_TR
end

function findTR_P_coupled5(Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,V_0,Vmaxlimit0,Vminlimit0,S_VP,PBound,raw,cl,acc_V0)

    # acc_V0=0.05
    Imaxlimit0_sqrt=sqrt.(Imaxlimit0)
    z0=0
    zmax=round(size(Imaxlimit0)[1]*size(Imaxlimit0)[2]*size(Imaxlimit0)[3]*P_con_percent)+1
    zmax_V=round(size(Vmaxlimit0)[1]*size(Vmaxlimit0)[2]*size(Vmaxlimit0)[3]*Q_con_percent_coupled)+1
    y=0
    z=0
    TR=0;TR_min=0.0
    TR_max=5;
    V0_length= length(V_0)
    acc_V=min(acc_V0,max(maximum(V_0-Vmaxlimit0),maximum(-(V_0-Vminlimit0))))
    iter=0
    Iter_max=10
    # atrisk=0
    atrisk= falses(size(I2_0));
    # Vmaxinrisk=0;Vmininrisk=0;
    Vmaxinrisk=falses(size(V_0));Vmininrisk=falses(size(V_0));
    # Vmax_VP=0;Vmin_VP=0;Imax_IP=0;
    Vmax_VP=zeros(size(V_0));Vmin_VP=zeros(size(V_0)); Imax_IP=zeros(size(I2_0));
    violated=0
    worst_vio=maximum(I2_0-Imaxlimit0)
    worst_vio_sqrt=max.(maximum(sqrt.(I2_0)-Imaxlimit0_sqrt),maximum(Base_TR_N),maximum(Base_TR_C))
    best_acc=deepcopy(worst_vio_sqrt)
    acc=0
    z_acc_VP=0
    adj_acc=0
    target_acc_sqrt=0
    Imaxlimit_new=deepcopy(Imaxlimit0)
    genbus=raw["Gen"].Bus
    genbus_s = Set(genbus)  # O(1) lookup
    allbus_s= Set(1:raw["nBus"])
    ovr_Sens_IP_up = zeros(nTP, raw["nLine"], cl);ovr_Sens_IP_dn = zeros(nTP, raw["nLine"], cl)
    ovr_Sens_VP_up = zeros(nTP, raw["nBus"], cl);ovr_Sens_VP_dn = zeros(nTP, raw["nBus"], cl)
    nBus      = raw["nBus"]
    nPActBus  = raw["nPActBus"]

    # Allocate once, contiguous in memory
    ovr_Sens_VP_upmat = zeros(Float64, nBus, nPActBus, nTP, cl)
    ovr_Sens_VP_dnmat = zeros(Float64, nBus, nPActBus, nTP, cl)
    ovr_Sens_IP_upmat = zeros(Float64, raw["nLine"]+raw["nTrans"], nPActBus, nTP, cl)
    ovr_Sens_IP_dnmat = zeros(Float64, raw["nLine"]+raw["nTrans"], nPActBus, nTP, cl)
    VP_up             = similar(ovr_Sens_VP_upmat)
    VP_dn             = similar(ovr_Sens_VP_upmat)
    # upelem            = similar(ovr_Sens_VP_upmat)
    # dnelem            = similar(ovr_Sens_VP_upmat)
    VP_up1      =falses(nBus, nPActBus, nTP, cl)
    VP_dn1      =similar(VP_up1)
    vals_IP= zeros(Float64, length(I2_0))
    vals_VP= zeros(Float64, 2* length(V_0))

    DV_up = V_0.- (acc_V .+ Vmaxlimit0)
    DV_dn = V_0.- (-acc_V .+ Vminlimit0)
    DV_up_4D = zeros(Float64, nBus, 1, nTP, cl)
    DV_dn_4D = zeros(Float64, nBus, 1, nTP, cl)
    S_VP_4D = zeros(Float64, nBus, nPActBus, nTP, cl)
    S_IP_4D = zeros(Float64, raw["nLine"]+raw["nTrans"], nPActBus, nTP, cl)
    Threads.@threads for t in 1:nTP
        for c in 1:cl
            @views S_VP_4D[:, :, t, c] .= S_VP[t, c]
            @views S_IP_4D[:, :, t, c] .= S_IP[t, c]
            @views DV_up_4D[:, 1, t, c] .= DV_up[t, :, c]  # direct assignment
            @views DV_dn_4D[:, 1, t, c] .= DV_dn[t, :, c]  # direct assignment
        end
    end
    S_VP1=keep_K_percent_twoside_4D(S_VP_4D,0.99)
    Pinj_up,Pinj_dn=0,0
    Gen_ind=indexin(raw["Gen"]."Bus", raw["PActBus"])
    report_TR=DataFrame(worst_vio_sqrt=0.0,acc=0.0,TR=0.0,z=z0,y=z0,z_acc=0.0,y_acc=0.0,ramp_acc=0.0,z_V=0,acc_VP=0.0)
    deleteat!(report_TR,1)
    TR=1
    z_acc=0
    y_acc=0
    ramp_acc=0
    TR_best=0.0
    S_VP_pos = max.(S_VP_4D, 0.0)
    S_VP_neg = min.(S_VP_4D, 0.0)
    S_IP_pos = max.(S_IP_4D, 0.0)
    S_IP_neg = min.(S_IP_4D, 0.0)
        while TR_max - TR_min > 1e-7# Stop when the range is small enough
            if iter==0
                nothing
            elseif iter>=Iter_max
                TR=TR_best
            elseif (z_acc)/max(y_acc,ramp_acc)>5 && y==0
                TR=max(TR/((z_acc)/max(y_acc)),1.2 *TR_min)
            else
                TR = (TR_min + TR_max) / 2
            end

            VP_up1.=0
            VP_dn1.=0
            Pinj_up=min.(TR,PBound[1])
            Pinj_dn=max.(-TR,PBound[2])
            Threads.@threads for t in 1:nTP
                for c in 1:cl
                    # @views S_IP_tc = S_IP[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
                    # @views S_VP_tc = S_VP[t,c]
                    @views up = Pinj_up[t, :, c]  # 1 × nPActBus
                    @views dn = Pinj_dn[t, :, c]
                    
                    @views @. ovr_Sens_IP_upmat[:, :, t, c] = S_IP_pos[:,:,t, c]* up' + S_IP_neg[:,:,t, c]* dn'
                    @views @. ovr_Sens_IP_dnmat[:, :, t, c] = S_IP_pos[:,:,t, c]* dn' + S_IP_neg[:,:,t, c]* up'
                    @views ovr_Sens_IP_up[t, :, c] .= sum(ovr_Sens_IP_upmat[:, :, t, c], dims=2)[1:raw["nLine"]] .+ 1e-6
                    @views ovr_Sens_IP_dn[t, :, c] .= sum(ovr_Sens_IP_dnmat[:, :, t, c], dims=2)[1:raw["nLine"]] .+ 1e-6

                    @views @. ovr_Sens_VP_upmat[:, :, t, c] .= S_VP_pos[:,:,t, c].* up' .+ S_VP_neg[:,:,t, c].* dn'
                    @views @. ovr_Sens_VP_dnmat[:, :, t, c] .= S_VP_pos[:,:,t, c].* dn' .+ S_VP_neg[:,:,t, c].* up'
                    @views ovr_Sens_VP_up[t, :, c] .= sum(ovr_Sens_VP_upmat[:, :, t, c], dims=2) .+ 1e-6
                    @views ovr_Sens_VP_dn[t, :, c] .= sum(ovr_Sens_VP_dnmat[:, :, t, c], dims=2) .+ 1e-6
                end
            end

            @views @. VP_up = ovr_Sens_VP_upmat + DV_up_4D
            @views @. VP_dn = -ovr_Sens_VP_dnmat - DV_dn_4D

            keep_K_percent_new_4D_V2!(VP_up1,VP_up,0.99)
            keep_K_percent_new_4D_V2!(VP_dn1,VP_dn,0.99)

            VP_up1 .|= S_VP1
            VP_dn1 .|= S_VP1
            @. ovr_Sens_VP_upmat *=VP_up1 #||S_VP1 
            @. ovr_Sens_VP_dnmat *=VP_dn1 #||S_VP1

            @views ovr_Sens_VP_up .= permutedims(sum(ovr_Sens_VP_upmat,dims=2)[:,1,:,:], (2, 1, 3)) 
            @views ovr_Sens_VP_dn .= permutedims(sum(ovr_Sens_VP_dnmat,dims=2)[:,1,:,:], (2, 1, 3)) 

            @views @. vals_IP = sqrt(max(0.0, I2_0[:] + ovr_Sens_IP_up[:])) - Imaxlimit0_sqrt[:]
            k_IP = min(length(vals_IP), Int(zmax))
            z_acc = partialsort(vals_IP, 1:k_IP; rev=true)[end]

            # @views vals_VP .= vcat(vec(V_0 .+ ovr_Sens_VP_up .- Vmaxlimit0),
            #    vec(-(V_0 .+ ovr_Sens_VP_dn .- Vminlimit0)))
            
             @views @. vals_VP[1:V0_length] = (V_0 + ovr_Sens_VP_up - Vmaxlimit0)[:]
             @views @. vals_VP[V0_length+1:end] = (-(V_0 .+ ovr_Sens_VP_dn .- Vminlimit0))[:]
        
            k_VP = min(length(vals_VP), Int(zmax_V))
            z_acc_VP = partialsort(vals_VP, 1:k_VP; rev=true)[end]

            @views Vmax_VP .= V_0 .+ ovr_Sens_VP_up
            @views Vmin_VP .= V_0 .+ ovr_Sens_VP_dn
            Vmaxinrisk .= Vmax_VP .>= Vmaxlimit0.+max(z_acc_VP,1e-4)
            Vmininrisk .= Vmin_VP .<= Vminlimit0.-max(z_acc_VP,1e-4)

        

            ramp_acc=0
            y_acc=max(maximum(sqrt.(max.(0,I2_0+ovr_Sens_IP_dn./P_eff)).-Imaxlimit0_sqrt),-0.1)
            acc=max(z_acc,y_acc)
            @. Imax_IP = I2_0 + ovr_Sens_IP_up/1
            @. atrisk = Imax_IP >= (Imaxlimit0_sqrt + max(z_acc,1e-4))^2

            z = sum(atrisk)
            z_V = sum(Vmaxinrisk)+sum(Vmininrisk)
            y = sum(@. (I2_0.+ovr_Sens_IP_dn./P_eff .>= (Imaxlimit0_sqrt .+ y_acc).^2))

            push!(report_TR,[worst_vio_sqrt,acc,TR,z,y,z_acc,y_acc,ramp_acc,z_V,z_acc_VP])
            if iter>=Iter_max
                break
            end

            if acc<=best_acc && z_acc_VP <= acc_V
                best_acc=max(0,acc)
                TR_best=TR
            end
            if z_acc>=max(y_acc) || y_acc<=0 || z_acc_VP >= acc_V
                TR_max = TR  # Reduce TR
            else
                TR_min = TR  # Increase TR
            end
            iter+=1
        end
    display(report_TR)
    Vmaxindex_P = findall(Vmaxinrisk)
    Vminindex_P = findall(Vmininrisk)
    println("Number of voltage in risk due to P adjustment: ", length(Vmaxindex_P)+length(Vminindex_P))
    println("voltage acc left out due to P adjustment: ", z_acc_VP)
    ImaxIndex = findall(atrisk)
    # ImaxIndex = findall(.~isnan.(I2_0));    z_acc=0;
    return TR,ImaxIndex,Vmaxindex_P,Vminindex_P,max(0.0,z_acc),[Pinj_up,Pinj_dn],worst_vio,worst_vio_sqrt,report_TR,Vmax_VP,Vmin_VP,Imax_IP
end


