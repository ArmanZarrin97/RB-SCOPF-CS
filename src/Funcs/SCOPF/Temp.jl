
        block=zeros(AffExpr,raw["nBus"])
        ΔV=fill(block,nTP,cl);Δθ=fill(block,nTP,cl);
        V=zeros(AffExpr,nTP,raw["nBus"],cl);θ=zeros(AffExpr,nTP,raw["nBus"],cl)
        V.=V_0;θ.=θ_0
        # abs.(SENSDict["S_VQ"][t,c]).>(0.0001/TR2)
        # abs.(SENSDict["S_VP"][t,c]).>(0.0001/TR1)
        # raw["PActBus"]
        # raw["QActBus"]
        # findall(x-> x!=1,raw["Bus"]."Bus_Type")
        for t=1:nTP,c=1:cl
            Qcoef_V=ifelse.(abs.(SENSDict["S_VQ"][t,c]).>(0.0001/TR2),SENSDict["S_VQ"][t,c],0.0)
            Pcoef_V=ifelse.(abs.(SENSDict["S_VP"][t,c]).>(0.0001/TR1),SENSDict["S_VP"][t,c],0.0)
            Pcoef_V[GBus,:]=SENSDict["S_VP"][t,c][GBus,:]
            Qcoef_V[GBus,:]=SENSDict["S_VQ"][t,c][GBus,:]
            # Pcoef_V=SENSDict["S_VP"][t,c]
            # Qcoef_V=SENSDict["S_VQ"][t,c]
            ΔV[t,c]=Pcoef_V*xΔPinj[t,:,c]+ Qcoef_V*xΔQinj[t,:,c]
            V[t,:,c]=V_0[t,:,c]+ΔV[t,c]
            Δθ[t,c]=SENSDict["S_θP"][t,c]*xΔPinj[t,:,c]+ SENSDict["S_θQ"][t,c]*xΔQinj[t,:,c]
            θ[t,:,c]=θ_0[t,:,c]+Δθ[t,c]
        end
        @expression(SCOPF,CV2magexp, V)
        @expression(SCOPF,Cθexp, θ)
        # minthreshold=0.001*raw["Bus"]."V_min (pu)"[1]/(0.1*(TRfree+1))
        # map_coefficients_inplace!.(x->ifelse.(abs(x) .> minthreshold, x, 0.0),NV2magexp)
        # drop_zeros!.(NV2magexp)
        # map_coefficients_inplace!.(x->round(x,digits=4),NV2magexp)

        # Vminindex=findall(x->x==1,V_0.<=(raw["Bus"]."V_min (pu)".+0.04)'[ones(Int8,nTP),:,ones(Int8,cl)])
        # setdiff!(Vminindex,findall(x->x==0,V-V_0))
        # Vmaxindex=findall(x->x==1,V_0.>=(raw["Bus"]."V_max (pu)".-0.04)'[ones(Int8,nTP),:,ones(Int8,cl)])
        # setdiff!(Vmaxindex,findall(x->x==0,V-V_0))
        # Vminlimit0=(raw["Bus"]."V_min (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
        # Vmaxlimit0=(raw["Bus"]."V_max (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
        Vminlimit0=(raw["Bus"]."V_min (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
        Vmaxlimit0=(raw["Bus"]."V_max (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
        Vmaxindex1= findall(x->x==1,Vmaxlimit0.-V_0 .<=ovr_Sens_V_Q*TR2)
        Vmaxindex2= findall(x->x==1,V_0.>=(raw["Bus"]."V_max (pu)".-0.00)'[ones(Int8,nTP),:,ones(Int8,cl)])
        Vmaxindex= union([Vmaxindex1;Vmaxindex2])
        Vminindex1= findall(x->x==1,V_0.- Vminlimit0.<=1*ovr_Sens_V_Q*TR2)
        Vminindex2= findall(x->x==1,V_0.<=(raw["Bus"]."V_min (pu)".+0.00)'[ones(Int8,nTP),:,ones(Int8,cl)])
        Vminindex= union([Vminindex1;Vminindex2])
        linear_Vminindex = [LinearIndices(CV2magexp)[ind] for ind in Vminindex]
        linear_Vmaxindex = [LinearIndices(CV2magexp)[ind] for ind in Vmaxindex]

        if relaxpermit==true
            @variable(SCOPF, slk_Vmax_c[ind in linear_Vmaxindex]>=0)
            @variable(SCOPF, slk_Vmin_c[ind in linear_Vminindex]>=0)
            # @expression(SCOPF,CVminlimit, Vminlimit0.*(1 .- slk_Vmin_c))
            # @expression(SCOPF,CVmaxlimit, Vmaxlimit0.*(1 .+ slk_Vmax_c))
            CVminlimit=Vminlimit0[linear_Vminindex].*(1 .- slk_Vmin_c)
            CVmaxlimit=Vmaxlimit0[linear_Vmaxindex].*(1 .+ slk_Vmax_c)
        else
            # @expression(SCOPF,CVminlimit, Vminlimit0)
            # @expression(SCOPF,CVmaxlimit, Vmaxlimit0)
            CVminlimit=Vminlimit0
            CVmaxlimit=Vmaxlimit0
        end
        @constraint(SCOPF,CV2mag1[ind in linear_Vminindex], CVminlimit[ind]<= CV2magexp[ind])
        @constraint(SCOPF,CV2mag2[ind in linear_Vmaxindex], CV2magexp[ind]<=CVmaxlimit[ind])
        # @constraint(SCOPF,CslackVangle, Cθexp[:,slackbus,:].==0)
    # Trust region and modelling constraints
    #---------------------------------------
        # Dont forget to uncomment the following lines
        if false #relaxpermit==true
            @variable(SCOPF, C_TRfreeP[t=1:nTP,i=1:raw["nPActBus"],c=1:cl]>=0,start=0,upper_bound=10)
            @variable(SCOPF, C_TRfreeQ[t=1:nTP,n=1:raw["nQActBus"],c=1:cl]>=0,start=0,upper_bound=10)
            @variable(SCOPF, C_SBfreeP[t=1:nTP,c=1:cl]>=0,start=0,upper_bound=100)
            @variable(SCOPF, C_SBfreeQ[t=1:nTP,c=1:cl]>=0,start=0,upper_bound=100)
        else
            @expression(SCOPF, C_TRfreeP, zeros(nTP,raw["nPActBus"],cl))
            @expression(SCOPF, C_TRfreeQ, zeros(nTP,raw["nQActBus"],cl))
            @expression(SCOPF, C_SBfreeP, zeros(nTP,cl))
            @expression(SCOPF, C_SBfreeQ, zeros(nTP,cl))
        end

        @constraint(SCOPF,CTrust_RegionP1,xΔPinj.<=TR1.*(C_TRfreeP.+1));
        @constraint(SCOPF,CTrust_RegionP2,-TR1.*(C_TRfreeP.+1).<=xΔPinj);
        @constraint(SCOPF,CTrust_RegionQ1,xΔQinj.<=TR2.*(C_TRfreeQ.+1));
        @constraint(SCOPF,CTrust_RegionQ2,-TR2.*(C_TRfreeQ.+1).<=xΔQinj);


        @constraint(SCOPF,CSensitivity_balance_P1[t=1:nTP,c=1:cl],sum(xΔPinj[t,:,c])<=Lininit["TRB"][1]*(C_SBfreeP[t,c]+1));
        @constraint(SCOPF,CSensitivity_balance_Q1[t=1:nTP,c=1:cl],sum(xΔQinj[t,:,c])<=300*Lininit["TRB"][1]*(C_SBfreeQ[t,c]+1));
        @constraint(SCOPF,CSensitivity_balance_P2[t=1:nTP,c=1:cl],-Lininit["TRB"][1]*(C_SBfreeP[t,c]+1)<=sum(xΔPinj[t,:,c]));
        @constraint(SCOPF,CSensitivity_balance_Q2[t=1:nTP,c=1:cl],-300*Lininit["TRB"][1]*(C_SBfreeQ[t,c]+1)<=sum(xΔQinj[t,:,c]));

        # fix.(xΔPinj,0;force=true)
        # fix.(xΔQinj,0;force=true)

        # @constraint(SCOPF, CTRLimit1[t=1:nTP,c=1:cl], sum(xΔPinj[t,:,c].^2)<=36)
        # @constraint(SCOPF, CTRLimit2[t=1:nTP,c=1:cl], sum(xΔQinj[t,:,c].^2)<=36)
    end