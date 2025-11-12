Runcode=(2,3)

# OP_init0=Initialize(Linearization_INIT[Runcode[1]],init_out,init_model,DA_out,DA_model)
# Ncost0,Vio_Sum_N0,OP_NValid0,Vio_general_N0,convergence_indicator0=PF_Valid_N(OP_init0)
OP_init0=Initialize(Linearization_INIT[Runcode[1]],init_out,init_model,DA_out,DA_model)
Ncost0,Vio_Sum_N0,Vio_general_N0,AVG_C_cost0,Vio_Sum_C0,OP_init0V,Vio_general_C0=AC_PFCheck(OP_init0,Runcode)



Head=["Change";"Change p.u.";"-";["ΔQinj$i" for i in 1:size(SysData["Gen"],1)];["Δθ$i" for i in 1:size(SysData["Bus"],1)];["ΔPinj$i" for i in 1:size(SysData["Gen"],1)];["ΔV$i" for i in 1:size(SysData["Bus"],1)];["ΔI2$i" for i in 1:size(SysData["Line"],1)];["ΔPnm$i" for i in 1:size(SysData["Line"],1)];["ΔQnm$i" for i in 1:size(SysData["Line"],1)];["ΔPmn$i" for i in 1:size(SysData["Line"],1)];["ΔQmn$i" for i in 1:size(SysData["Line"],1)]]
Head=["Change";"Change p.u.";"-";"ΔQinj" ;"Δθ";"ΔPinj";"ΔV";"ΔI2";"ΔPnm";"ΔQnm";"ΔPmn";"ΔQmn"]

global Sensivity_errror=DataFrame(repeat([[]],size(Head,1)), Head)

OP_init=OP_init0V
Ptype_="DA";init_="Zero";
global nTP=nTL;
global TD=Td;
raw,initraw,profraw= try Data_Reader(;Ptype=Ptype_ , init=init_)
catch e
    @warn "Problem with running Data_Reader function"
    rethrow(e)
end
ctdall=ctdMaker(raw)
raw["Bus"][!,:PActive].=0;raw["Bus"][!,:QActive].=0
# raw["PActBus"]= filter(x-> !isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
# raw["QActBus"]= filter(x-> !isempty(x."shunt") || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
# raw["PActBus"]= 1:raw["nBus"]
# raw["QActBus"]= 1:raw["nBus"]
if LinMethod=="SENS"
    loadsheddingactivity=maximum(OP_init["Pshed"],dims=1)[1,:]+ ((Problem_Type=="OPF") ? (zeros(raw["nLoad"])) : maximum(OP_init["Pshed_c_Line"],dims=(1,3))[1,:,1])
    Activeloadbus=findall(x->x>1e-4,loadsheddingactivity)
else
    Activeloadbus=1:raw["nLoad"]
end
if loadbusinjection=="true"
    raw["PActBus"]= filter(x-> !isempty(x."load")||!isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
    raw["QActBus"]= filter(x-> !isempty(x."load")||!isempty(x."shunt") || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
else
    raw["PActBus"]= filter(x-> !isempty(intersect(x."load",Activeloadbus))||!isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
    raw["QActBus"]= filter(x-> !isempty(intersect(x."load",Activeloadbus)) || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus] #||!isempty(x."shunt")
    # raw["PActBus"]= filter(x-> !isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
    # raw["QActBus"]= filter(x-> !isempty(x."shunt") || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
end
raw["nPActBus"]=size(raw["PActBus"])[1];raw["nQActBus"]=size(raw["QActBus"])[1]
raw["Bus"][raw["PActBus"],"PActive"].=1;raw["Bus"][raw["QActBus"],"QActive"].=1
Lininit=OP_init


SENSDict = Cal_SENS_Simple(Lininit,ctdall,nTP,raw;C="N")

minimum(raw["Gen"]."cost_b")
incb=filter(x->x."cost_b"==minimum(raw["Gen"]."cost_b"),raw["Gen"])[1,"Bus"]
decb=filter(x->x."cost_b"==maximum(raw["Gen"]."cost_b"),raw["Gen"])[1,"Bus"]
namerow=Dict("P"=>"P$incb+,P$decb-", "Q"=>"Q$incb+,Q$decb-")
for change in [1 0.1 0.01], Var in ["P","Q"]
    DeltaQ=zeros(raw["nBus"]);DeltaP=zeros(raw["nBus"]);
    if Var=="P" DeltaP[incb]=change;DeltaP[decb]=-change
    elseif Var=="Q" DeltaQ[incb]=change;DeltaQ[decb]=-change end
    expectedΔQG=DeltaQ[raw["Gen"]."Bus"]
    expectedΔPG=DeltaP[raw["Gen"]."Bus"]
    if Var=="P"
        expectedΔV=(SensModel[:S_VP][:,:,1,1]*DeltaP)'
        expectedΔθ=(SensModel[:S_θP][:,:,1,1]*DeltaP)'
        expectedΔI2=(SensModel[:S_IP][:,:,1,1]*DeltaP)[1:raw["nLine"]]'
        expectedΔPnm=(SensModel[:S_PsP][:,:,1,1]*DeltaP)[1:raw["nLine"]]'
        expectedΔQnm=(SensModel[:S_QsP][:,:,1,1]*DeltaP)[1:raw["nLine"]]'
        expectedΔPmn=(SensModel[:S_PrP][:,:,1,1]*DeltaP)[1:raw["nLine"]]'
        expectedΔQmn=(SensModel[:S_QrP][:,:,1,1]*DeltaP)[1:raw["nLine"]]'
    elseif Var=="Q"
        expectedΔV=(SensModel[:S_VQ][:,:,1,1]*DeltaQ)'
        expectedΔθ=(SensModel[:S_θQ][:,:,1,1]*DeltaQ)'
        expectedΔI2=(SensModel[:S_IQ][:,:,1,1]*DeltaQ)[1:raw["nLine"]]'
        expectedΔPnm=(SensModel[:S_PsQ][:,:,1,1]*DeltaQ)[1:raw["nLine"]]'
        expectedΔQnm=(SensModel[:S_QsQ][:,:,1,1]*DeltaQ)[1:raw["nLine"]]'
        expectedΔPmn=(SensModel[:S_PrQ][:,:,1,1]*DeltaQ)[1:raw["nLine"]]'
        expectedΔQmn=(SensModel[:S_QrQ][:,:,1,1]*DeltaQ)[1:raw["nLine"]]'
    end
    # Vestimated=OP_init["V"].+expectedΔV# θestimated=OP_init["θ"].+expectedΔθ# PGestimated=OP_init["PG"].+Delta[3:end]# QGestimated=OP_init["QG"]

    NewOP=deepcopy(OP_init)
    if Var=="P" NewOP["PG"][filter(x-> x."Bus"==incb, raw["Gen"])[1,"G_N"]]+=change;NewOP["PG"][filter(x-> x."Bus"==decb, raw["Gen"])[1,"G_N"]]-=change end
    NewOP["V"]+=expectedΔV
    Ncost0,Vio_Sum_N,NewOP_V,Vio_general_N,convergence_indicator=PF_Valid_N(NewOP)
    realΔV=NewOP_V["V"]-OP_init["V"]
    realΔθ=NewOP_V["θ"]-OP_init["θ"]
    realΔQG=NewOP_V["QG"]-OP_init["QG"]
    realΔPG=NewOP_V["PG"]-OP_init["PG"]
    realΔI2=Calculate_I(NewOP_V,raw)-Calculate_I(OP_init,raw)
    realΔPnm,realΔQnm,realΔPmn,realΔQmn=Calculate_PQflow(NewOP_V,raw).-Calculate_PQflow(OP_init,raw)

    # errorθ=(θnew-θestimated).^2# meanerrorθ=sqrt.(sum(errorθ))/raw["nBus"]# maxerrorθ=sqrt(maximum(errorθ))
    # errorQG=(NewOP_V["QG"].-OP_init["QG"]).^2# meanerrorQG=sqrt.(sum(errorQG))/raw["nGen"]# maxerrorQG=sqrt(maximum(errorQG))
    # newrow=[namerow[Var] change "real" realΔQG realΔθ realΔPG realΔV realΔI2 realΔPnm realΔQnm realΔPmn realΔQmn]
    newrow=[namerow[Var] change "real" 0 0 0 0 0 0 0 0 0]
    push!(Sensivity_errror,newrow[1,:])
    # newrow=[namerow[Var] change "expected" expectedΔQG' expectedΔθ expectedΔPG' expectedΔV expectedΔI2 expectedΔPnm expectedΔQnm expectedΔPmn expectedΔQmn]
    # newrow=[namerow[Var] change "real" mean(sqrt.((realΔQG-expectedΔQG').^2)) mean(sqrt.((realΔθ-expectedΔθ).^2)) mean(sqrt.((realΔPG-expectedΔPG').^2)) mean(sqrt.((realΔV-expectedΔV).^2)) mean(sqrt.((realΔI2-expectedΔI2).^2)) mean(sqrt.((realΔPnm-expectedΔPnm).^2)) mean(sqrt.((realΔQnm-expectedΔQnm).^2)) mean(sqrt.((realΔPmn-expectedΔPmn).^2)) mean(sqrt.((realΔQmn-expectedΔQmn).^2))]
    newrow=[namerow[Var] change "real" mean(abs.(realΔQG-expectedΔQG')) mean(abs.(realΔθ-expectedΔθ)) mean(abs.(realΔPG-expectedΔPG')) mean(abs.(realΔV-expectedΔV)) mean(abs.(realΔI2-expectedΔI2)) mean(abs.(realΔPnm-expectedΔPnm)) mean(abs.(realΔQnm-expectedΔQnm)) mean(abs.(realΔPmn-expectedΔPmn)) mean(abs.(realΔQmn-expectedΔQmn))]
    push!(Sensivity_errror,newrow[1,:])
end
ods_write(Path*"../Result/Sensivities Error",Dict(("General",2,1)=>Sensivity_errror))

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


