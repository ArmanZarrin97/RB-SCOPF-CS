# include("Funcs/IO/Packages.jl") # load necessary packages
# include("Funcs/IO/utils.jl") # load general functions
# include("Assets/L-MPSCOPF.jl") # load the multi-period security-constrained optimal power flow model
# #JuMP: Modeling Optimization Problem #Ipopt: Optimization solver package # OdsIO: Input & output for .ods files # DataFrames: !!
# include("Assets/ACPF.jl") # load power flow model
# include("Assets/linearization.jl") # load power flow model
# include("Funcs/IO/Result_generation.jl") # save MPSCOPF results in Output file
# include("Assets/L-Problems.jl") # load the different models of problems

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


Runcode=(2,3)

OP_init0=Initialize(Linearization_INIT[Runcode[1]],init_out,init_model,DA_out,DA_model)
Ncost0,Vio_Sum_N,OP_NValid,Vio_general_N=PF_Valid_N(OP_init0,Runcode)
if Problem_Type=="SCOPF"
    AVG_C_cost0,Vio_Sum_C0,OP_init,Vio_general_C0=PF_Valid_C(OP_NValid,Runcode)
else
    OP_init=OP_NValid
end

# Head=["Change";"Change p.u.";"-";["ΔQinj$i" for i in 1:size(SysData["Gen"],1)];["Δθ$i" for i in 1:size(SysData["Bus"],1)];["ΔPinj$i" for i in 1:size(SysData["Gen"],1)];["ΔV$i" for i in 1:size(SysData["Bus"],1)];["ΔI2$i" for i in 1:size(SysData["Line"],1)];["ΔPnm$i" for i in 1:size(SysData["Line"],1)];["ΔQnm$i" for i in 1:size(SysData["Line"],1)];["ΔPmn$i" for i in 1:size(SysData["Line"],1)];["ΔQmn$i" for i in 1:size(SysData["Line"],1)]]
Head=["Change";"Change p.u.";"-";["ΔQinj$i" for i in 1:size(SysData["Gen"],1)];["Δθ$i" for i in 1:size(SysData["Bus"],1)];["ΔPinj$i" for i in 1:size(SysData["Gen"],1)];["ΔV$i" for i in 1:size(SysData["Bus"],1)];["ΔI2$i" for i in 1:size(SysData["Line"],1)]]

global Sensivity_errror=DataFrame(repeat([[]],size(Head,1)), Head)

# OP_init=OP_NValid
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
raw["PActBus"]= filter(x-> !isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
raw["QActBus"]= filter(x-> !isempty(x."shunt") || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
raw["nPActBus"]=size(raw["PActBus"])[1];raw["nQActBus"]=size(raw["QActBus"])[1]
raw["Bus"][raw["PActBus"],"PActive"].=1;raw["Bus"][raw["QActBus"],"QActive"].=1
Lininit=OP_init



# Lininit["NSens"][1][5,:].=0
# Lininit["CSens"][1][5,:].=0
SensModel=Cal_SENS_Simple(Lininit,nTP,raw;C="N")
SensModel_old,Pa_IV_old,Pa_Iθ_old=Cal_SENS_Simple_old(Lininit,nTP,raw;C="N")
# SensModel_old=Cal_SENS(Lininit,nTP,raw;C="N")

# SensModel["S_θQ"][1,1]

# SensModel_old["S_θQ"][:,3:5,:,:]

# Pa_Iθ[1]
# Pa_Iθ_old[:,:,1,1]

# key="S_PsP"
# sizeof(SensModel[key])
# SensModel[key]
# SENSDict=deepcopy(SensModel)
# threshold=0.1
# SENSDict[key] .= ifelse.(abs.(SENSDict[key]) .> threshold, SENSDict[key], 0.0)
# sizeof(SENSDict[key])
# sparse(SENSDict[key][:,:,1,1])
# sparse(SensModel[key][:,:,1,1])

namerow=Dict("P"=>"P5+,P3-", "Q"=>"Q5+,Q3-")
for change in [1 0.1 0.01], Var in ["P","Q"]
    DeltaQ=zeros(5);DeltaP=zeros(5);
    if Var=="P" DeltaP[5]=change;DeltaP[3]=-change
    elseif Var=="Q" DeltaQ[5]=change;DeltaQ[3]=-change end
    DeltaQA=DeltaQ[raw["QActBus"]]
    DeltaPA=DeltaP[raw["PActBus"]]
    # DeltaQA=DeltaQ
    # DeltaPA=DeltaP
    expectedΔQG=DeltaQ[3:5]
    expectedΔPG=DeltaP[3:5]
    if Var=="P"
        expectedΔV=(SensModel["S_VP"][1,1]*DeltaPA)'
        expectedΔθ=(SensModel["S_θP"][1,1]*DeltaPA)'
        expectedΔI2=(SensModel["S_IP"][1,1]*DeltaPA)'
        # expectedΔV=(SensModel_old["S_VP"][:,:,1,1]*DeltaPA)'
        # expectedΔθ=(SensModel_old["S_θP"][:,:,1,1]*DeltaPA)'
        # expectedΔI2=(SensModel_old["S_IP"][:,:,1,1]*DeltaPA)'
        # expectedΔPnm=(SensModel["S_PsP"][1,1]*DeltaP)'
        # expectedΔQnm=(SensModel["S_QsP"][1,1]*DeltaP)'
        # expectedΔPmn=(SensModel["S_PrP"][1,1]*DeltaP)'
        # expectedΔQmn=(SensModel["S_QrP"][1,1]*DeltaP)'
    elseif Var=="Q"
        expectedΔV=(SensModel["S_VQ"][1,1]*DeltaQA)'
        expectedΔθ=(SensModel["S_θQ"][1,1]*DeltaQA)'
        expectedΔI2=(SensModel["S_IQ"][1,1]*DeltaQA)'
        # expectedΔV=(SensModel_old["S_VQ"][:,:,1,1]*DeltaQA)'
        # expectedΔθ=(SensModel_old["S_θQ"][:,:,1,1]*DeltaQA)'
        # expectedΔI2=(SensModel_old["S_IQ"][:,:,1,1]*DeltaQA)'
        # expectedΔPnm=(SensModel["S_PsQ"][1,1]*DeltaQ)'
        # expectedΔQnm=(SensModel["S_QsQ"][1,1]*DeltaQ)'
        # expectedΔPmn=(SensModel["S_PrQ"][1,1]*DeltaQ)'
        # expectedΔQmn=(SensModel["S_QrQ"][1,1]*DeltaQ)'
    end
    # Vestimated=OP_init["V"].+expectedΔV# θestimated=OP_init["θ"].+expectedΔθ# PGestimated=OP_init["PG"].+Delta[3:end]# QGestimated=OP_init["QG"]

    NewOP=deepcopy(OP_init)
    if Var=="P" NewOP["PG"][3]+=change;NewOP["PG"][1]-=change end
    NewOP["V"]+=expectedΔV
    Ncost,Vio_Sum,NewOP_V,Vio_general,convergence_indicator=PF_Valid_N(NewOP,Runcode)
    realΔV=NewOP_V["V"]-OP_init["V"]
    realΔθ=NewOP_V["θ"]-OP_init["θ"]
    realΔQG=NewOP_V["QG"]-OP_init["QG"]
    realΔPG=NewOP_V["PG"]-OP_init["PG"]
    realΔI2=Calculate_I(NewOP_V,raw)-Calculate_I(OP_init,raw)
    # realΔPnm,realΔQnm,realΔPmn,realΔQmn=Calculate_PQflow(NewOP_V,raw).-Calculate_PQflow(OP_init,raw)

    # errorθ=(θnew-θestimated).^2# meanerrorθ=sqrt.(sum(errorθ))/raw["nBus"]# maxerrorθ=sqrt(maximum(errorθ))
    # errorQG=(NewOP_V["QG"].-OP_init["QG"]).^2# meanerrorQG=sqrt.(sum(errorQG))/raw["nGen"]# maxerrorQG=sqrt(maximum(errorQG))
    # newrow=[namerow[Var] change "real" realΔQG realΔθ realΔPG realΔV realΔI2 realΔPnm realΔQnm realΔPmn realΔQmn]
    newrow=[namerow[Var] change "real" realΔQG realΔθ realΔPG realΔV realΔI2]
    push!(Sensivity_errror,newrow[1,:])
    # newrow=[namerow[Var] change "expected" expectedΔQG' expectedΔθ expectedΔPG' expectedΔV expectedΔI2 expectedΔPnm expectedΔQnm expectedΔPmn expectedΔQmn]
    newrow=[namerow[Var] change "expected" expectedΔQG' expectedΔθ expectedΔPG' expectedΔV expectedΔI2]
    push!(Sensivity_errror,newrow[1,:])
end
ods_write(Path*"../Result/Sensivities Error",Dict(("General",2,1)=>Sensivity_errror))


# ### First check active power sensitivity
# change=1
# Delta=zeros(5);Delta[5]=change;Delta[4]=-change

# expectedV=(SensModel["S_VQ"][:,:,1,1]*Delta)'
# expectedθ=(SensModel["S_θQ"][:,:,1,1]*Delta)'
# NewOP=deepcopy(OP_init)
# NewOP["PG"][3]+=change
# NewOP["PG"][2]-=change
# NewOP["V"]+=expectedV
# Ncost,Vio_Sum,NewOP_V,Vio_general,convergence_indicator=PF_Valid_N(NewOP)
# NewOP_V["V"].-OP_init["V"]
# NewOP_V["θ"].-OP_init["θ"]
# NewOP_V["QG"].-OP_init["QG"]

# ### Second check reactive power sensitivity
# change=1
# Delta=zeros(5);Delta[3]=change;Delta[4]=-Delta[3]

# expectedV=(SensModel["S_VQ"][:,:,1,1]*Delta)'
# expectedθ=(SensModel["S_θQ"][:,:,1,1]*Delta)'
# NewOP=deepcopy(OP_init)
# NewOP["V"]+=expectedV
# Ncost,Vio_Sum,NewOP_V,Vio_general,convergence_indicator=PF_Valid_N(NewOP)
# NewOP_V["V"].-OP_init["V"]
# NewOP_V["θ"].-OP_init["θ"]
# NewOP_V["QG"].-OP_init["QG"]


function Cal_SENS_Simple_old(Lininit::Dict{String,Array{Union{Float64,Array{Float64}}}},nTP::Int64,raw::Dict{String, Any};C::String="N")
    nbus::Int64=raw["nBus"];nline::Int64=raw["nLine"];ntrans::Int64=raw["nTrans"];
    if C=="N"
        CC=1;
        V_0::Array{Float64}=Lininit["V"];θ_0::Array{Float64}=Lininit["θ"];
        exc=[0]
        SENS::Array{Array{Float64}}=Lininit["NSens"];
    elseif C=="CL"
        CC::Int64=raw["nConL"];
        V_0=Lininit["V_c_Line"];θ_0=Lininit["θ_c_Line"];#Union{Array{Float64, 3}, Matrix{Float64},Float64}
        exc::Vector{Int64}=raw["ConL"]."#";
        SENS=Lininit["CSens"];#::Union{Matrix{Matrix{Float64}},Vector{Matrix{Float64}}}
    end
    TR=ones(Float64,nline+ntrans)
    if ntrans!=0 TR[nline+1:end]= raw["Trans"]."ratio" end

    Pa_IV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Iθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);
    # Pa_PsV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Psθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);
    # Pa_PrV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Prθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);
    # Pa_QsV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Qsθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);
    # Pa_QrV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Qrθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);

    Names=["From","To","g (pu)","b (pu)","g_sh (pu)","b_sh (pu)"]
    df::DataFrame=vcat(raw["Line"][:,Names],raw["Trans"][:,Names])
    # df=raw["Line"];
    g_sh::Vector{Int64}=df."g_sh (pu)";b_sh::Vector{Float64}=df."b_sh (pu)";g::Vector{Float64}=df."g (pu)";b::Vector{Float64}=df."b (pu)";
    for l::Int64=1:nline+ntrans,i::Int64=1:nbus,t::Int64=1:nTP,c::Int64=1:CC
        Δθ::Float64=θ_0[t,df."From"[l],c]-θ_0[t,df."To"[l],c]
        Vn::Float64=V_0[t,df."From"[l],c];Vm::Float64=V_0[t,df."To"[l],c]
        if df[l,"From"]==i
            Pa_IV[l,i,t,c]= (g[l]^2+b[l]^2)*(2*Vn-2*Vm*cos(Δθ))
            Pa_Iθ[l,i,t,c]= (g[l]^2+b[l]^2)*(2*Vn*Vm*sin(Δθ))
            # Pa_PsV[l,i,t,c]=    (g_sh[l]/2+g[l])*Vn*2+  TR[l]*( -g[l]*Vm*cos(Δθ)    -b[l]*Vm*sin(Δθ))
            # Pa_Psθ[l,i,t,c]=                            TR[l]*( +g[l]*Vn*Vm*sin(Δθ) -b[l]*Vn*Vm*cos(Δθ))
            # Pa_QsV[l,i,t,c]=    (-b_sh[l]/2-b[l])*Vn*2+ TR[l]*( +b[l]*Vm*cos(Δθ)    -g[l]*Vm*sin(Δθ))
            # Pa_Qsθ[l,i,t,c]=                            TR[l]*( -b[l]*Vn*Vm*sin(Δθ) -g[l]*Vn*Vm*cos(Δθ))
            # Pa_PrV[l,i,t,c]=                            TR[l]*( -g[l]*Vm*cos(Δθ)    +b[l]*Vm*sin(Δθ))
            # Pa_Prθ[l,i,t,c]=                            TR[l]*( +g[l]*Vn*Vm*sin(Δθ) +b[l]*Vn*Vm*cos(Δθ))
            # Pa_QrV[l,i,t,c]=                            TR[l]*( +b[l]*Vm*cos(Δθ)    +g[l]*Vm*sin(Δθ))
            # Pa_Qrθ[l,i,t,c]=                            TR[l]*( -b[l]*Vn*Vm*sin(Δθ) +g[l]*Vn*Vm*cos(Δθ))
        elseif df[l,"To"]==i
            Pa_IV[l,i,t,c]= (g[l]^2+b[l]^2)*(2*Vm-2*Vn*cos(Δθ))
            Pa_Iθ[l,i,t,c]= (g[l]^2+b[l]^2)*(-2*Vn*Vm*sin(Δθ))
            # Pa_PsV[l,i,t,c]=                                    TR[l]*( -g[l]*Vn*cos(Δθ)    -b[l]*Vn*sin(Δθ))
            # Pa_Psθ[l,i,t,c]=                                    TR[l]*( -g[l]*Vn*Vm*sin(Δθ) +b[l]*Vn*Vm*cos(Δθ))
            # Pa_QsV[l,i,t,c]=                                    TR[l]*( +b[l]*Vn*cos(Δθ)    -g[l]*Vn*sin(Δθ))
            # Pa_Qsθ[l,i,t,c]=                                    TR[l]*( +b[l]*Vn*Vm*sin(Δθ) +g[l]*Vn*Vm*cos(Δθ))
            # Pa_PrV[l,i,t,c]=(TR[l]^2)*(g_sh[l]/2+g[l])*Vm*2+    TR[l]*( -g[l]*Vn*cos(Δθ)    +b[l]*Vn*sin(Δθ))
            # Pa_Prθ[l,i,t,c]=                                    TR[l]*( -g[l]*Vn*Vm*sin(Δθ) -b[l]*Vn*Vm*cos(Δθ))
            # Pa_QrV[l,i,t,c]= (TR[l]^2)*(-b_sh[l]/2-b[l])*Vm*2+  TR[l]*( +b[l]*Vm*cos(Δθ)    +g[l]*Vm*sin(Δθ))
            # Pa_Qrθ[l,i,t,c]=                                    TR[l]*( +b[l]*Vn*Vm*sin(Δθ) -g[l]*Vn*Vm*cos(Δθ))
        end
    end
    SENSDict = Dict{String, Array{Float64, 4}}(
        "S_IP" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_IQ" => zeros(nline + ntrans, nbus, nTP, CC),
        # "S_PrP" => zeros(nline + ntrans, nbus, nTP, CC),
        # "S_PsP" => zeros(nline + ntrans, nbus, nTP, CC),
        # "S_QrP" => zeros(nline + ntrans, nbus, nTP, CC),
        # "S_QsP" => zeros(nline + ntrans, nbus, nTP, CC),
        # "S_PrQ" => zeros(nline + ntrans, nbus, nTP, CC),
        # "S_PsQ" => zeros(nline + ntrans, nbus, nTP, CC),
        # "S_QrQ" => zeros(nline + ntrans, nbus, nTP, CC),
        # "S_QsQ" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_VP" => zeros(nbus, nbus, nTP, CC),
        "S_VQ" => zeros(nbus, nbus, nTP, CC),
        "S_θP" => zeros(nbus, nbus, nTP, CC),
        "S_θQ" => zeros(nbus, nbus, nTP, CC)
    )
    for i::Int64=1:nbus,t::Int64=1:nTP,c::Int64=1:CC
        SENStemp=SENS[t,c]
        for l::Int64=1:nline+ntrans
            if exc[c]==l
                for text in ["S_IP","S_IQ"] #,"S_PrP","S_PsP","S_QrP","S_QsP","S_PrQ","S_PsQ","S_QrQ","S_QsQ"
                    SENSDict[text][l,i,t,c]=0
                end
            else
                sb::Int64=df."From"[l];rb::Int64=df."To"[l];
                SENSDict["S_IP"][l,i,t,c]= (Pa_IV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_IV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Iθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Iθ[l,rb,t,c]*SENStemp[rb,i]);
                SENSDict["S_IQ"][l,i,t,c]= (Pa_IV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_IV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Iθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Iθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
                # SENSDict["S_PrP"][l,i,t,c]= (Pa_PrV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_PrV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Prθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Prθ[l,rb,t,c]*SENStemp[rb,i]);
                # SENSDict["S_PsP"][l,i,t,c]= (Pa_PsV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_PsV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Psθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Psθ[l,rb,t,c]*SENStemp[rb,i]);
                # SENSDict["S_QrP"][l,i,t,c]= (Pa_QrV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_QrV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Qrθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Qrθ[l,rb,t,c]*SENStemp[rb,i]);
                # SENSDict["S_QsP"][l,i,t,c]= (Pa_QsV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_QsV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Qsθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Qsθ[l,rb,t,c]*SENStemp[rb,i]);
                # SENSDict["S_PrQ"][l,i,t,c]= (Pa_PrV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_PrV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Prθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Prθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
                # SENSDict["S_PsQ"][l,i,t,c]= (Pa_PsV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_PsV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Psθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Psθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
                # SENSDict["S_QrQ"][l,i,t,c]= (Pa_QrV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_QrV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Qrθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Qrθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
                # SENSDict["S_QsQ"][l,i,t,c]= (Pa_QsV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_QsV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Qsθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Qsθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
            end
        end
        for j::Int64=1:nbus
            SENSDict["S_VP"][j,i,t,c]= SENStemp[nbus+j,i]
            SENSDict["S_VQ"][j,i,t,c]= SENStemp[nbus+j,nbus+i]
            SENSDict["S_θP"][j,i,t,c]= SENStemp[j,i]
            SENSDict["S_θQ"][j,i,t,c]= SENStemp[j,nbus+i]
        end
    end

    threshold = 0.0001
    for key in ["S_VP" "S_VQ" "S_IP" "S_IQ"]#collect(keys(SENSDict))
        # Broadcast comparison and set elements less than threshold to 0
        SENSDict[key] .= ifelse.(abs.(SENSDict[key]) .> threshold, SENSDict[key], 0.0)
    end
    # threshold = 0.001
    # for key in [ "S_PrP" "S_PsP" "S_PrQ" "S_PsQ"]#collect(keys(SENSDict)) "S_VP" "S_VQ"
    #     # Broadcast comparison and set elements less than threshold to 0
    #     SENSDict[key] .= ifelse.(abs.(SENSDict[key]) .> threshold, SENSDict[key], 0.0)
    # end
    # threshold = 0.001
    # for key in [ "S_QrP" "S_QsP" "S_QrQ" "S_QsQ"]#collect(keys(SENSDict)) "S_VP" "S_VQ"
    #     # Broadcast comparison and set elements less than threshold to 0
    #     SENSDict[key] .= ifelse.(abs.(SENSDict[key]) .> threshold, SENSDict[key], 0.0)
    # end
    return SENSDict,Pa_IV,Pa_Iθ
end