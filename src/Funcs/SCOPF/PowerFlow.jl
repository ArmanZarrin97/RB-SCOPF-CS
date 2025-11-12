
function Flow_Expr(nTP,raw::DataFrame,xe,xf, TorL,Cont=[]; inv='t')
    df=raw
    nl=size(df,1)
    if nl==0 return (0,0,0,0); end
    cls=size(Cont,1)
    if cls==0
        cl=1;exc=0;
    else
        cl=cls;exc=Cont."#";
    end
    if TorL=='L'
        tr=ones(nl)
    elseif TorL=='T'
        exc=zeros(cl);tr=raw."ratio"
    end
    gsh=df."g_sh (pu)";bsh=df."b_sh (pu)";g=df."g (pu)";b=df."b (pu)";


    temp=Model(); # Create a temporary model to remove the needless expression
    @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (gsh[n]/2)*(xe[t,df."From"[n],c]^2+xf[t,df."From"[n],c]^2))
    @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (g[n])*(xe[t,df."From"[n],c]^2+xf[t,df."From"[n],c]^2))
    @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : (tr[n])*(
                                        -(g[n])*(xe[t,df."From"[n],c]*xe[t,df."To"[n],c]+xf[t,df."From"[n],c]*xf[t,df."To"[n],c])
                                        -(b[n])*(xf[t,df."From"[n],c]*xe[t,df."To"[n],c]-xe[t,df."From"[n],c]*xf[t,df."To"[n],c])))

    @expression(temp, Npij, sh.+sr1.+sr2)


    temp=Model()
    @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-bsh[n]/2)*(xe[t,df."From"[n],c]^2+xf[t,df."From"[n],c]^2))
    @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-b[n])*(xe[t,df."From"[n],c]^2+xf[t,df."From"[n],c]^2))
    @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : tr[n]*(
                                        (b[n])*(xe[t,df."From"[n],c]*xe[t,df."To"[n],c]+xf[t,df."From"[n],c]*xf[t,df."To"[n],c])
                                        -(g[n])*(xf[t,df."From"[n],c]*xe[t,df."To"[n],c]-xe[t,df."From"[n],c]*xf[t,df."To"[n],c])))
    @expression(temp, Nqij, sh.+sr1.+sr2)



    if(inv=='t')
        temp=Model()
        @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (tr[n]^2)*(gsh[n]/2)*(xe[t,df."To"[n],c]^2+xf[t,df."To"[n],c]^2))
        @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (tr[n]^2)*(g[n])*(xe[t,df."To"[n],c]^2+xf[t,df."To"[n],c]^2))
        @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (tr[n])*(
                                                -(g[n])*(xe[t,df."From"[n],c]*xe[t,df."To"[n],c]+xf[t,df."From"[n],c]*xf[t,df."To"[n],c])
                                                +(b[n])*(xf[t,df."From"[n],c]*xe[t,df."To"[n],c]-xe[t,df."From"[n],c]*xf[t,df."To"[n],c])))
        @expression(temp, Npji, sh.+sr1.+sr2)



        temp=Model()
        @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-bsh[n]/2)*(xe[t,df."To"[n],c]^2+xf[t,df."To"[n],c]^2))
        @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-b[n])*(xe[t,df."To"[n],c]^2+xf[t,df."To"[n],c]^2))
        @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : tr[n]*(
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
        return (Npij, Nqij, Npji, Nqji)
    end
end

function Power_Inj(nTP,raw,initraw,x_ch,x_dis, x_inc, x_dec, xshed,xcurt,xij,xji,xijT,xjiT;xPQ='P',C="N")
    temp=Model()
    if C=="N"
        CC=1;
    elseif C=="CL"
        CC=raw["nConL"];
    end
    @expression(temp, xinj[t=1:nTP, i=1:raw["nBus"],c=1:CC], sum(xij[t,n,c] for n in 1:raw["nLine"] if (raw["Line"]."From"[n]==i))
                                                + sum(xji[t,n,c] for n in 1:raw["nLine"] if (raw["Line"]."To"[n]==i))
                                                + sum(xijT[t,n,c] for n in 1:raw["nTrans"] if (raw["Trans"]."From"[n]==i))
                                                + sum(xjiT[t,n,c] for n in 1:raw["nTrans"] if (raw["Trans"]."To"[n]==i))
                                                - (xPQ=='P' ? (sum(x_dis[t,n,c]+initraw["Pdis"][t,n]- x_ch[t,n,c]-initraw["Pch"][t,n] for n in 1:raw["nESS"] if (raw["ESS"]."Bus"[n]==i);init=0)) : 0)
# for coordinating with mousa voltage is 0      - (xPQ=='Q' ? (xe[t,i,c]^2+xf[t,i,c]^2)*(sum(raw["Shunt"]."bsh0"[n] for n in 1:raw["nShunt"] if raw["Shunt"]."Bus"[n]==i;init=0)) : 0)
                                                - (xPQ=='Q' ? (sum(raw["Shunt"]."bsh0"[n] for n in 1:raw["nShunt"] if raw["Shunt"]."Bus"[n]==i;init=0)) : 0)
                                                - (xPQ=='Q' ? (sum((xshed[t,n,c]+initraw["Pshed"][t,n])*raw["Load"]."PF"[n] for n in 1:raw["nLoad"] if (raw["Load"]."Bus_N"[n]==i);init=0)) : 0)
                                                - (xPQ=='P' ? (sum(x_dec[t,n,c]+initraw["Pdec"][t,n]- x_inc[t,n,c]-initraw["Pinc"][t,n] for n in 1:raw["nFL"] if (raw["FL"]."Bus"[n]==i);init=0)) : 0)
                                                - (xPQ=='P' ? (sum(xshed[t,n,c]+initraw["Pshed"][t,n] for n in 1:raw["nLoad"] if (raw["Load"]."Bus_N"[n]==i);init=0)) : 0)
                                                + (xPQ=='P' ? (sum(xcurt[t,n,c]+initraw["Pcurt"][t,n] for n in 1:raw["nRES"] if (raw["RES"]."Bus"[n]==i);init=0)) : 0)
                                                )
return xinj
end
