function Non_Linear_PB(nTP,nl,cl,exc,TR,df,Vars,TorL,inv='t')

    gsh=df."g_sh (pu)";bsh=df."b_sh (pu)";g=df."g (pu)";b=df."b (pu)";
    xe,xf=Vars
    temp=Model(); # Create a temporary model to remove the needless expression
        @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (gsh[n]/2)*(xe[t,df."From"[n],c]^2+xf[t,df."From"[n],c]^2))
        @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (g[n])*(xe[t,df."From"[n],c]^2+xf[t,df."From"[n],c]^2))
        @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : (TR[n])*(
                                            -(g[n])*(xe[t,df."From"[n],c]*xe[t,df."To"[n],c]+xf[t,df."From"[n],c]*xf[t,df."To"[n],c])
                                            -(b[n])*(xf[t,df."From"[n],c]*xe[t,df."To"[n],c]-xe[t,df."From"[n],c]*xf[t,df."To"[n],c])))
        @expression(temp, Npij, sh.+sr1.+sr2)
    temp=Model()
        @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-bsh[n]/2)*(xe[t,df."From"[n],c]^2+xf[t,df."From"[n],c]^2))
        @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-b[n])*(xe[t,df."From"[n],c]^2+xf[t,df."From"[n],c]^2))
        @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : TR[n]*(
                                            (b[n])*(xe[t,df."From"[n],c]*xe[t,df."To"[n],c]+xf[t,df."From"[n],c]*xf[t,df."To"[n],c])
                                            -(g[n])*(xf[t,df."From"[n],c]*xe[t,df."To"[n],c]-xe[t,df."From"[n],c]*xf[t,df."To"[n],c])))
        @expression(temp, Nqij, sh.+sr1.+sr2)

    if(inv=='t')
        temp=Model()
            @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(gsh[n]/2)*(xe[t,df."To"[n],c]^2+xf[t,df."To"[n],c]^2))
            @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(g[n])*(xe[t,df."To"[n],c]^2+xf[t,df."To"[n],c]^2))
            @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                                    -(g[n])*(xe[t,df."From"[n],c]*xe[t,df."To"[n],c]+xf[t,df."From"[n],c]*xf[t,df."To"[n],c])
                                                    +(b[n])*(xf[t,df."From"[n],c]*xe[t,df."To"[n],c]-xe[t,df."From"[n],c]*xf[t,df."To"[n],c])))
            @expression(temp, Npji, sh.+sr1.+sr2)



        temp=Model()
            @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-bsh[n]/2)*(xe[t,df."To"[n],c]^2+xf[t,df."To"[n],c]^2))
            @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-b[n])*(xe[t,df."To"[n],c]^2+xf[t,df."To"[n],c]^2))
            @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : TR[n]*(
                                                    (b[n])*(xe[t,df."From"[n],c]*xe[t,df."To"[n],c]+xf[t,df."From"[n],c]*xf[t,df."To"[n],c])
                                                +(g[n])*(xf[t,df."From"[n],c]*xe[t,df."To"[n],c]-xe[t,df."From"[n],c]*xf[t,df."To"[n],c])))
            @expression(temp, Nqji, sh.+sr1.+sr2)

    end
    if TorL == 'L'
        @expression(temp,NIflow2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : (g[n]^2 + b[n]^2)
                                            *((xe[t,df."From"[n],c]-xe[t,df."To"[n],c])^2+(xf[t,df."From"[n],c]-xf[t,df."To"[n],c])^2))# lineflow current limitation
    end
    if TorL == 'L'
        return (Npij, Nqij, Npji, Nqji, NIflow2)
    elseif TorL == 'T'
        return (Npij, Nqij, Npji, Nqji,0)
    end
end
function Linear_Yang_PB(Lininit::Dict{String,Array{Union{Float64,Array{Float64}}}},nTP,nl,cl,exc,TR,df,Vars,TorL,State,inv='t')
    if State=="Norm"
        V_0=Lininit["V"];θ_0=Lininit["θ"];
    elseif State=="Cont"
        V_0=Lininit["V_c_Line"];θ_0=Lininit["θ_c_Line"];
    end
    gsh=df."g_sh (pu)";bsh=df."b_sh (pu)";g=df."g (pu)";b=df."b (pu)";
    xV2,xθ=Vars
    temp=Model(); # Create a temporary model to remove the needless expression
        @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (gsh[n]/2)*(xV2[t,df."From"[n],c])) #Line flow expression
        @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (g[n])*(xV2[t,df."From"[n],c])) #Line flow expression
        @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
                                    -(g[n]/2)*(xV2[t,df."From"[n],c]+xV2[t,df."To"[n],c])
                                    -(b[n])*(xθ[t,df."From"[n],c]-xθ[t,df."To"[n],c]))) #Line flow expression

        @expression(temp, lossNpij_V[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                (g[n]/2)*(2*((V_0[t,df."From"[n],c]-V_0[t,df."To"[n],c])/(V_0[t,df."From"[n],c]+V_0[t,df."To"[n],c]))*(xV2[t,df."From"[n],c]-xV2[t,df."To"[n],c])-(V_0[t,df."From"[n],c]-V_0[t,df."To"[n],c])^2)))
        @expression(temp, loss_θ[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
                                (g[n])*((θ_0[t,df."From"[n],c]-θ_0[t,df."To"[n],c])*(xθ[t,df."From"[n],c]-xθ[t,df."To"[n],c])-((θ_0[t,df."From"[n],c]-θ_0[t,df."To"[n],c])^2)/2)))
        @expression(temp, lossij, lossNpij_V .+ loss_θ)
        @expression(temp, sr3, sr2.+ lossij) #Line flow expression
        @expression(temp, Npij, sh.+sr1.+sr2.+lossij)
    temp=Model()
        @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-bsh[n]/2)*(xV2[t,df."From"[n],c])) #Line flow expression
        @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-b[n])*(xV2[t,df."From"[n],c])) #Line flow expression

        @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                    (b[n]/2)*(xV2[t,df."From"[n],c]+xV2[t,df."To"[n],c])
                                    -(g[n])*(xθ[t,df."From"[n],c]-xθ[t,df."To"[n],c]))) #Line flow expression
        @expression(temp, loss_V[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
                                (-b[n]/2)*(2*((V_0[t,df."From"[n],c]-V_0[t,df."To"[n],c])/(V_0[t,df."From"[n],c]+V_0[t,df."To"[n],c]))*(xV2[t,df."From"[n],c]-xV2[t,df."To"[n],c])-(V_0[t,df."From"[n],c]-V_0[t,df."To"[n],c])^2)))
        @expression(temp, loss_θ[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
                                (-b[n])*((θ_0[t,df."From"[n],c]-θ_0[t,df."To"[n],c])*(xθ[t,df."From"[n],c]-xθ[t,df."To"[n],c])-((θ_0[t,df."From"[n],c]-θ_0[t,df."To"[n],c])^2)/2)))
        @expression(temp, lossij, loss_V .+ loss_θ)
        @expression(temp, Nqij, sh .+ sr1 .+ sr2 .+ lossij)

    if(inv=='t')
        temp=Model()
            @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(gsh[n]/2)*(xV2[t,df."To"[n],c])) #Line flow expression
            @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(g[n])*(xV2[t,df."To"[n],c])) #Line flow expression

            @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                        -(g[n]/2)*(xV2[t,df."To"[n],c]+xV2[t,df."From"[n],c])
                                        -(b[n])*(xθ[t,df."To"[n],c]-xθ[t,df."From"[n],c])))#Line flow expression

            @expression(temp, lossNpji_V[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                    (g[n]/2)*(2*((V_0[t,df."To"[n],c]-V_0[t,df."From"[n],c])/(V_0[t,df."To"[n],c]+V_0[t,df."From"[n],c]))*(xV2[t,df."To"[n],c]-xV2[t,df."From"[n],c])-(V_0[t,df."To"[n],c]-V_0[t,df."From"[n],c])^2)))
            @expression(temp, loss_θ[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
                                    (g[n])*((θ_0[t,df."To"[n],c]-θ_0[t,df."From"[n],c])*(xθ[t,df."To"[n],c]-xθ[t,df."From"[n],c])-((θ_0[t,df."To"[n],c]-θ_0[t,df."From"[n],c])^2)/2)))
            @expression(temp, lossji, lossNpji_V .+ loss_θ)
            @expression(temp, Npji, sh .+ sr1 .+ sr2 .+ lossji)


        temp=Model()
            @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(-bsh[n]/2)*(xV2[t,df."To"[n],c])) #Line flow expression
            @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(-b[n])*(xV2[t,df."To"[n],c])) #Line flow expression

            @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                        (b[n]/2)*(xV2[t,df."To"[n],c]+xV2[t,df."From"[n],c])
                                        -(g[n])*(xθ[t,df."To"[n],c]-xθ[t,df."From"[n],c]))) #Line flow expression

            @expression(temp, loss_V[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
            (-b[n]/2)*(2*((V_0[t,df."To"[n],c]-V_0[t,df."From"[n],c])/(V_0[t,df."To"[n],c]+V_0[t,df."From"[n],c]))*(xV2[t,df."To"[n],c]-xV2[t,df."From"[n],c])-(V_0[t,df."To"[n],c]-V_0[t,df."From"[n],c])^2)))
            @expression(temp, loss_θ[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
            (-b[n])*((θ_0[t,df."To"[n],c]-θ_0[t,df."From"[n],c])*(xθ[t,df."To"[n],c]-xθ[t,df."From"[n],c])-((θ_0[t,df."To"[n],c]-θ_0[t,df."From"[n],c])^2)/2)))
            @expression(temp, lossji, loss_V .+ loss_θ)
            @expression(temp, Nqji, sh .+ sr1 .+ sr2 .+ lossji)

    end

    Parlist=[Npij, Nqij,lossNpij_V,lossNpji_V]
    if TorL == 'L'
        return (Npij, Nqij, Npji, Nqji, Parlist)
    elseif TorL == 'T'
        return (Npij, Nqij, Npji, Nqji)
    end

end
function Linear_SOCP_PB(Lininit,nTP,nl,cl,exc,TR,df,Vars,TorL,State,inv='t')
    if State=="Norm"
        V_0=Lininit["V"];θ_0=Lininit["θ"];
    elseif State=="Cont"
        V_0=Lininit["V_c_Line"];θ_0=Lininit["θ_c_Line"];
    end
    gsh=df."g_sh (pu)";bsh=df."b_sh (pu)";g=df."g (pu)";b=df."b (pu)";
    xV2,xA,xB=Vars
    temp=Model(); # Create a temporary model to remove the needless expression
        @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (gsh[n]/2)*(xV2[t,df."From"[n],c])) #Line flow expression
        @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (g[n])*(xV2[t,df."From"[n],c])) #Line flow expression
        @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
                                    -(g[n])*(xA[t,n,c])
                                    -(b[n])*(xB[t,n,c]))) #Line flow expression

        @expression(temp, Npij, sh.+sr1.+sr2)
    temp=Model()
        @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-bsh[n]/2)*(xV2[t,df."From"[n],c])) #Line flow expression
        @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-b[n])*(xV2[t,df."From"[n],c])) #Line flow expression

        @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                                                            +(b[n])*(xA[t,n,c])
                                                                            -(g[n])*(xB[t,n,c]))) #Line flow expression
        @expression(temp, Nqij, sh .+ sr1 .+ sr2)

    if(inv=='t')
        temp=Model()
            @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(gsh[n]/2)*(xV2[t,df."To"[n],c])) #Line flow expression
            @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(g[n])*(xV2[t,df."To"[n],c])) #Line flow expression

            @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                                                            -(g[n])*(xA[t,n,c])
                                                                            +(b[n])*(xB[t,n,c])))#Line flow expression
            @expression(temp, Npji, sh .+ sr1 .+ sr2)


        temp=Model()
            @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(-bsh[n]/2)*(xV2[t,df."To"[n],c])) #Line flow expression
            @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(-b[n])*(xV2[t,df."To"[n],c])) #Line flow expression

            @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                                                            +(g[n])*(xA[t,n,c])
                                                                            +(b[n])*(xB[t,n,c])))#Line flow expression

            @expression(temp, Nqji, sh .+ sr1 .+ sr2)

    end

    if TorL == 'L'
        @expression(temp,NIflow2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : (g[n]^2 + b[n]^2)
                                            *(xV2[t,df."From"[n],c]+xV2[t,df."To"[n],c]-2*xA[t,n,c]))# lineflow current limitation
    end

    if TorL == 'L'
        return (Npij, Nqij, Npji, Nqji, NIflow2)
    elseif TorL == 'T'
        return (Npij, Nqij, Npji, Nqji)
    end

end
function Linear_SENS_PB(Lininit::Dict{String,Array{Union{Float64,Array{Float64}}}},nTP::Int64,nl,cl,exc,Vars,TorL,State,SENSDict::Dict{String, Array{Float64, 4}},raw::Dict{String, Any},Cont,inv='t')
    xΔPinj,xΔQinj=Vars
    if State=="Norm"
        V_0,θ_0=Lininit["V"],Lininit["θ"];
    elseif State=="Cont"
        V_0,θ_0=Lininit["V_c_Line"],Lininit["θ_c_Line"];
    end
    initialVar=(V_0.*cos.(θ_0), V_0.*sin.(θ_0));
    pij, qij, pji, qji,I2 =Flow_Expr("None",Dict{String,Array{Union{Float64,Array{Float64}}}}(),nTP,raw,initialVar,TorL;Cont=Cont);
    if TorL=='L' n_1=0 elseif TorL=='T' n_1=raw["nLine"]  end

    Npij=[AffExpr(0.0) for t=1:nTP,n=1:nl,c=1:cl]
    Nqij=[AffExpr(0.0) for t=1:nTP,n=1:nl,c=1:cl]
    if TorL == 'L'
        NIflow2=[AffExpr(0.0) for t=1:nTP,n=1:nl,c=1:cl]
    end
    if(inv=='t')
        Npji=[AffExpr(0.0) for t=1:nTP,n=1:nl,c=1:cl]
        Nqji=[AffExpr(0.0) for t=1:nTP,n=1:nl,c=1:cl]
    end
    for t=1:nTP,n=1:nl,c=1:cl
        if exc[c] != n+n_1
            add_to_expression!(Npij[t,n,c],pij[t,n,c])
            add_to_expression!(Nqij[t,n,c],qij[t,n,c])
            if inv=='t'
                add_to_expression!(Npji[t,n,c],pji[t,n,c])
                add_to_expression!(Nqji[t,n,c],qji[t,n,c])
            end
            if TorL == 'L'
                add_to_expression!(NIflow2[t,n,c],I2[t,n,c])
            end
            for i=1:raw["nBus"]
                add_to_expression!(Npij[t,n,c],SENSDict["S_PsP"][n+n_1,i,t,c]*xΔPinj[t,i,c])
                add_to_expression!(Npij[t,n,c],SENSDict["S_PsQ"][n+n_1,i,t,c]*xΔQinj[t,i,c])
                add_to_expression!(Nqij[t,n,c],SENSDict["S_QsP"][n+n_1,i,t,c]*xΔPinj[t,i,c])
                add_to_expression!(Nqij[t,n,c],SENSDict["S_QsQ"][n+n_1,i,t,c]*xΔQinj[t,i,c])
                if inv=='t'
                    add_to_expression!(Npji[t,n,c],SENSDict["S_PrP"][n+n_1,i,t,c]*xΔPinj[t,i,c])
                    add_to_expression!(Npji[t,n,c],SENSDict["S_PrQ"][n+n_1,i,t,c]*xΔQinj[t,i,c])
                    add_to_expression!(Nqji[t,n,c],SENSDict["S_QrP"][n+n_1,i,t,c]*xΔPinj[t,i,c])
                    add_to_expression!(Nqji[t,n,c],SENSDict["S_QrQ"][n+n_1,i,t,c]*xΔQinj[t,i,c])
                end
                if TorL == 'L'
                    add_to_expression!(NIflow2[t,n,c],SENSDict["S_IP"][n+n_1,i,t,c]*xΔPinj[t,i,c])
                    add_to_expression!(NIflow2[t,n,c],SENSDict["S_IQ"][n+n_1,i,t,c]*xΔQinj[t,i,c])
                end
            end
        end
    end

    if TorL == 'L'
        Parlist= NIflow2,V_0,θ_0
        return (Npij, Nqij, Npji, Nqji, Parlist)
    elseif TorL == 'T'
        return (Npij, Nqij, Npji, Nqji)
    end
end