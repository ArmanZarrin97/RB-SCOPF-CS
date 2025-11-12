# filtered_matrix = map(x -> x >= 1 ? x : 0, matrix)
function filter_matrix!(Mat::Array{Float64},threshold;type="Upper")
    if type=="Upper"
        Mat = map(x -> x >= threshold ? x : 0, Mat)
    elseif type=="Lower"
        Mat = map(x -> x <= threshold ? x : 0, Mat)
    end
    return Mat
end

SensModel["S_IQ"][:,:,1,1]=filter_matrix!(SensModel["S_IQ"][:,:,1,1],1.5)

filter(x-> x<1.5,abs.(SensModel["S_IQ"][:,:,1,1]))

median(abs.(SensModel["S_IP"][:,:,1,1]))
SensModel

function Calculate_I(OP,raw)
    V=OP["V"];θ=OP["θ"];df=raw["Line"]
    I2=zeros(nTP,raw["nLine"])
    for t=1:nTP,n=1:raw["nLine"]
        I2[t,n]=(df."g (pu)"[n]^2 + df."b (pu)"[n]^2)*(V[t,df."From"[n]]^2+V[t,df."To"[n]]^2-2*V[t,df."From"[n]]*V[t,df."To"[n]]*cos(θ[t,df."From"[n]]-θ[t,df."To"[n]]))
    end
    return I2
end
function Calculate_PQflow(OP,raw)
    V=OP["V"];θ=OP["θ"];df=raw["Line"]
    Pnm,Qnm,Pmn,Qmn=zeros(nTP,raw["nLine"]),zeros(nTP,raw["nLine"]),zeros(nTP,raw["nLine"]),zeros(nTP,raw["nLine"])
    for t=1:nTP,n=1:raw["nLine"]
        Pnm[t,n]=   (   df."g_sh (pu)"[n]/2+df."g (pu)"[n])*(V[t,df."From"[n]]^2)   -df."g (pu)"[n]*(V[t,df."From"[n]]*V[t,df."To"[n]]*cos(θ[t,df."From"[n]]-θ[t,df."To"[n]]))      -df."b (pu)"[n]*(V[t,df."From"[n]]*V[t,df."To"[n]]*sin(θ[t,df."From"[n]]-θ[t,df."To"[n]]))
        Qnm[t,n]=   (  -df."b_sh (pu)"[n]/2-df."b (pu)"[n])*(V[t,df."From"[n]]^2)   +df."b (pu)"[n]*(V[t,df."From"[n]]*V[t,df."To"[n]]*cos(θ[t,df."From"[n]]-θ[t,df."To"[n]]))      -df."g (pu)"[n]*(V[t,df."From"[n]]*V[t,df."To"[n]]*sin(θ[t,df."From"[n]]-θ[t,df."To"[n]]))
        Pmn[t,n]=   (   df."g_sh (pu)"[n]/2+df."g (pu)"[n])*(V[t,df."To"[n]]^2)     -df."g (pu)"[n]*(V[t,df."To"[n]]*V[t,df."From"[n]]*cos(θ[t,df."To"[n]]-θ[t,df."From"[n]]))      +df."b (pu)"[n]*(V[t,df."From"[n]]*V[t,df."To"[n]]*sin(θ[t,df."From"[n]]-θ[t,df."To"[n]]))
        Qmn[t,n]=   (  -df."b_sh (pu)"[n]/2-df."b (pu)"[n])*(V[t,df."To"[n]]^2)     +df."b (pu)"[n]*(V[t,df."To"[n]]*V[t,df."From"[n]]*cos(θ[t,df."To"[n]]-θ[t,df."From"[n]]))      +df."g (pu)"[n]*(V[t,df."From"[n]]*V[t,df."To"[n]]*sin(θ[t,df."From"[n]]-θ[t,df."To"[n]]))
    end
    return Pnm,Qnm,Pmn,Qmn
end

# function SENSfilter(SensModel,OP_init,raw)
#     for i in ["S_IQ","S_IP"]
#         if sum(abs.(SensModel["S_IQ"][:,:,1,1][1,:]))
#         SensModel[i]=filter_matrix!(SENSmodel[i],1.5)
#     end

#     sum(abs.(SensModel["S_IQ"][:,:,1,1][1,:]))
#     raw["Line"]."RATE_A"
#     (raw["Line"]."RATE_A").*raw["sbase"]-Calculate_I(OP_init,raw)[1,:]
#     Calculate_I(OP_init,raw)
#     OP_init
# end
# OP_init=OP_init0

function SENSfilter(SensModel,OP_init,raw)
    i="S_PsQ"
    for i in ["S_IQ","S_IP"]
        # SensModel[i]=filter_matrix!(SENSmodel[i],)
        for l in 1:size(SensModel[i],1)
            SensModel[i][l,:,:,:] = map(x -> x >= 1.5 ? x : 0, SensModel[i][l,:,:])
        end
        SensModel[i] = map(x -> x >= Mat[3] ? x : 0, SensModel[i])
    end

    sum(abs.(SensModel["S_IQ"][:,:,1,1][1,:]))
    raw["Line"]."RATE_A"
    (raw["Line"]."RATE_A").*raw["sbase"]-Calculate_I(OP_init,raw)[1,:]
    Calculate_I(OP_init,raw)
    OP_init
end
OP_init=OP_init0

Mat = map(x -> x >= Mat[3] ? x : 0, Mat)
function filter_matrix!(Mat::Array{Float64},threshold;type="Upper")
    if type=="Upper"
        Mat = map(x -> x >= threshold ? x : 0, Mat)
    elseif type=="Lower"
        Mat = map(x -> x <= threshold ? x : 0, Mat)
    end
    return Mat
end