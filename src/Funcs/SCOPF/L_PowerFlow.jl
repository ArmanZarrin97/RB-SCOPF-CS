function Cal_SENS(Lininit::Dict{String,Any},nTP::Int64,raw::Dict{String, Any};C::String="N")
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
    Pa_PsV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Psθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);
    Pa_PrV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Prθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);
    Pa_QsV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Qsθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);
    Pa_QrV=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);Pa_Qrθ=Array{Float64,4}(undef, nline+ntrans,nbus,nTP,CC);

    Names=["From","To","g (pu)","b (pu)","g_sh (pu)","b_sh (pu)"]
    df::DataFrame=vcat(raw["Line"][:,Names],raw["Trans"][:,Names])
    # df=raw["Line"];
    g_sh::Vector{Int64}=df."g_sh (pu)";b_sh::Vector{Float64}=df."b_sh (pu)";g::Vector{Float64}=df."g (pu)";b::Vector{Float64}=df."b (pu)";
    for l::Int64=1:nline+ntrans,i::Int64=1:nbus,t::Int64=1:nTP,c::Int64=1:CC
        Δθ::Float64=θ_0[t,fb[l],c]-θ_0[t,tb[l],c]
        Vn::Float64=V_0[t,fb[l],c];Vm::Float64=V_0[t,tb[l],c]
        if df[l,"From"]==i
            Pa_IV[l,i,t,c]= (g[l]^2+b[l]^2)*(2*Vn-2*Vm*cos(Δθ))
            Pa_Iθ[l,i,t,c]= (g[l]^2+b[l]^2)*(2*Vn*Vm*sin(Δθ))
            Pa_PsV[l,i,t,c]=    (g_sh[l]/2+g[l])*Vn*2+  TR[l]*( -g[l]*Vm*cos(Δθ)    -b[l]*Vm*sin(Δθ))
            Pa_Psθ[l,i,t,c]=                            TR[l]*( +g[l]*Vn*Vm*sin(Δθ) -b[l]*Vn*Vm*cos(Δθ))
            Pa_QsV[l,i,t,c]=    (-b_sh[l]/2-b[l])*Vn*2+ TR[l]*( +b[l]*Vm*cos(Δθ)    -g[l]*Vm*sin(Δθ))
            Pa_Qsθ[l,i,t,c]=                            TR[l]*( -b[l]*Vn*Vm*sin(Δθ) -g[l]*Vn*Vm*cos(Δθ))
            Pa_PrV[l,i,t,c]=                            TR[l]*( -g[l]*Vm*cos(Δθ)    +b[l]*Vm*sin(Δθ))
            Pa_Prθ[l,i,t,c]=                            TR[l]*( +g[l]*Vn*Vm*sin(Δθ) +b[l]*Vn*Vm*cos(Δθ))
            Pa_QrV[l,i,t,c]=                            TR[l]*( +b[l]*Vm*cos(Δθ)    +g[l]*Vm*sin(Δθ))
            Pa_Qrθ[l,i,t,c]=                            TR[l]*( -b[l]*Vn*Vm*sin(Δθ) +g[l]*Vn*Vm*cos(Δθ))


        elseif df[l,"To"]==i

            Pa_IV[l,i,t,c]= (g[l]^2+b[l]^2)*(2*Vm-2*Vn*cos(Δθ))
            Pa_Iθ[l,i,t,c]= (g[l]^2+b[l]^2)*(-2*Vn*Vm*sin(Δθ))
            Pa_PsV[l,i,t,c]=                                    TR[l]*( -g[l]*Vn*cos(Δθ)    -b[l]*Vn*sin(Δθ))
            Pa_Psθ[l,i,t,c]=                                    TR[l]*( -g[l]*Vn*Vm*sin(Δθ) +b[l]*Vn*Vm*cos(Δθ))
            Pa_QsV[l,i,t,c]=                                    TR[l]*( +b[l]*Vn*cos(Δθ)    -g[l]*Vn*sin(Δθ))
            Pa_Qsθ[l,i,t,c]=                                    TR[l]*( +b[l]*Vn*Vm*sin(Δθ) +g[l]*Vn*Vm*cos(Δθ))
            Pa_PrV[l,i,t,c]=(TR[l]^2)*(g_sh[l]/2+g[l])*Vm*2+    TR[l]*( -g[l]*Vn*cos(Δθ)    +b[l]*Vn*sin(Δθ))
            Pa_Prθ[l,i,t,c]=                                    TR[l]*( -g[l]*Vn*Vm*sin(Δθ) -b[l]*Vn*Vm*cos(Δθ))
            Pa_QrV[l,i,t,c]= (TR[l]^2)*(-b_sh[l]/2-b[l])*Vm*2+  TR[l]*( +b[l]*Vm*cos(Δθ)    +g[l]*Vm*sin(Δθ))
            Pa_Qrθ[l,i,t,c]=                                    TR[l]*( +b[l]*Vn*Vm*sin(Δθ) -g[l]*Vn*Vm*cos(Δθ))
        end
    end

    SENSDict = Dict{String, Array{Float64, 4}}(
        "S_IP" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_IQ" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_PrP" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_PsP" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_QrP" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_QsP" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_PrQ" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_PsQ" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_QrQ" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_QsQ" => zeros(nline + ntrans, nbus, nTP, CC),
        "S_VP" => zeros(nbus, nbus, nTP, CC),
        "S_VQ" => zeros(nbus, nbus, nTP, CC),
        "S_θP" => zeros(nbus, nbus, nTP, CC),
        "S_θQ" => zeros(nbus, nbus, nTP, CC)
    )

    for i::Int64=1:nbus,t::Int64=1:nTP,c::Int64=1:CC
        SENStemp=SENS[t,c]
        for l::Int64=1:nline+ntrans
            if exc[c]==l
                for text in ["S_IP","S_IQ","S_PrP","S_PsP","S_QrP","S_QsP","S_PrQ","S_PsQ","S_QrQ","S_QsQ"]
                    SENSDict[text][l,i,t,c]=0
                end
            else
                sb::Int64=fb[l];rb::Int64=tb[l];
                SENSDict["S_IP"][l,i,t,c]= (Pa_IV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_IV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Iθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Iθ[l,rb,t,c]*SENStemp[rb,i]);
                SENSDict["S_IQ"][l,i,t,c]= (Pa_IV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_IV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Iθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Iθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
                SENSDict["S_PrP"][l,i,t,c]= (Pa_PrV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_PrV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Prθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Prθ[l,rb,t,c]*SENStemp[rb,i]);
                SENSDict["S_PsP"][l,i,t,c]= (Pa_PsV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_PsV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Psθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Psθ[l,rb,t,c]*SENStemp[rb,i]);
                SENSDict["S_QrP"][l,i,t,c]= (Pa_QrV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_QrV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Qrθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Qrθ[l,rb,t,c]*SENStemp[rb,i]);
                SENSDict["S_QsP"][l,i,t,c]= (Pa_QsV[l,sb,t,c]*SENStemp[nbus+sb,i]+ Pa_QsV[l,rb,t,c]*SENStemp[nbus+rb,i]+ Pa_Qsθ[l,sb,t,c]*SENStemp[sb,i]+ Pa_Qsθ[l,rb,t,c]*SENStemp[rb,i]);
                SENSDict["S_PrQ"][l,i,t,c]= (Pa_PrV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_PrV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Prθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Prθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
                SENSDict["S_PsQ"][l,i,t,c]= (Pa_PsV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_PsV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Psθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Psθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
                SENSDict["S_QrQ"][l,i,t,c]= (Pa_QrV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_QrV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Qrθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Qrθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
                SENSDict["S_QsQ"][l,i,t,c]= (Pa_QsV[l,sb,t,c]*SENStemp[nbus+sb,i+nbus]+ Pa_QsV[l,rb,t,c]*SENStemp[nbus+rb,i+nbus]+ Pa_Qsθ[l,sb,t,c]*SENStemp[sb,nbus+i]+ Pa_Qsθ[l,rb,t,c]*SENStemp[rb,nbus+i]);
            end
        end
        for j::Int64=1:nbus
            SENSDict["S_VP"][j,i,t,c]= SENStemp[nbus+j,i]
            SENSDict["S_VQ"][j,i,t,c]= SENStemp[nbus+j,nbus+i]
            SENSDict["S_θP"][j,i,t,c]= SENStemp[j,i]
            SENSDict["S_θQ"][j,i,t,c]= SENStemp[j,nbus+i]
        end
    end

    # threshold = 0.01
    # for key in ["S_VP" "S_VQ"]#collect(keys(SENSDict))
    #     # Broadcast comparison and set elements less than threshold to 0
    #     SENSDict[key] .= ifelse.(abs.(SENSDict[key]) .> threshold, SENSDict[key], 0.0)
    # end
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
    return SENSDict
end
function Cal_SENS_Simple(Lininit::Dict{String,Any},ctdall,nTP::Int64,raw::Dict{String, Any};C::String="N")
    nbus::Int64=raw["nBus"];nline::Int64=raw["nLine"];ntrans::Int64=raw["nTrans"];nPbus::Int64=raw["nPActBus"];nQbus::Int64=raw["nQActBus"];
    icondens=50;
    Pbus=raw["PActBus"];Qbus=raw["QActBus"];
    if C=="N"
        CC=1;
        V_0::Array{Float64}=Lininit["V"];θ_0::Array{Float64}=Lininit["θ"];
        exc=[[]]
        SENS::Array{Array{Float64}}=Lininit["NSens"];
    elseif C=="CL"
        CC::Int64=raw["nConL"];
        V_0=Lininit["V_c_Line"];θ_0=Lininit["θ_c_Line"];#Union{Array{Float64, 3}, Matrix{Float64},Float64}
        exc=raw["ConL"]."#";
        SENS=Lininit["CSens"];#::Union{Matrix{Matrix{Float64}},Vector{Matrix{Float64}}}
    end
    TR=ones(Float64,nline+ntrans)
    if ntrans!=0 TR[nline+1:end]= raw["Trans"]."ratio" end
    block=spzeros(Float64,nline+ntrans,nbus)
    pa_block=fill(block,nTP,CC)
    block2=spzeros(Float64,nbus,nbus)
    pa_block2=fill(block2,nTP,CC)

    Pa_IV=deepcopy(pa_block);Pa_Iθ=deepcopy(pa_block);

    Names=["From","To","g (pu)","b (pu)","g_sh (pu)","b_sh (pu)"]
    df::DataFrame=vcat(raw["Line"][:,Names],raw["Trans"][:,Names])
    # df=raw["Line"];
    g::Vector{Float64}=df."g (pu)";b::Vector{Float64}=df."b (pu)";
    fb::Vector{Int64}=df."From";tb::Vector{Int64}=df."To";
    Δθ=@. θ_0[:,fb,:]-θ_0[:,tb,:];
    Vn=V_0[:,fb,:];Vm=V_0[:,tb,:];
    term=@. 2*Vn-2*Vm*cos(Δθ)
    Line_admittance= @. g^2+b^2
    Line_admittancemat=(Line_admittance')[ones(Int8,nTP),:,ones(Int8,CC)]
    indexfb=CartesianIndex.(1:nline+ntrans, fb)
    indextb=CartesianIndex.(1:nline+ntrans, tb)

    blockP=spzeros(Float64,nline+ntrans,nPbus)
    pa_blockP=fill(copy(blockP),nTP,CC)
    block2P=spzeros(Float64,nPbus,nPbus)
    pa_block2P=fill(copy(block2P),nTP,CC)

    blockQ=spzeros(Float64,nline+ntrans,nQbus)
    pa_blockQ=fill(copy(blockQ),nTP,CC)
    block2Q=spzeros(Float64,nQbus,nQbus)
    pa_block2Q=fill(copy(block2Q),nTP,CC)

    SENSDict = Dict(
        "S_IP" => copy(pa_blockP),
        "S_IQ" => copy(pa_blockQ),
        "S_VP" => copy(pa_block2P),
        "S_VQ" => copy(pa_block2Q),
        "S_θP" => copy(pa_block2P),
        "S_θQ" => copy(pa_block2Q)
    )
    thresholdpercent = Lininit["FilterIndex"][1]
    for t=1:nTP,c=1:CC
        Pa_IV[t,c][indexfb]=(@. Line_admittancemat*((2*Vn)-2*Vm.*cos.(Δθ)))[t,:,c]
        Pa_IV[t,c][indextb]=(@. Line_admittancemat*((2*Vm)-2*Vn*cos(Δθ)))[t,:,c]
        term= (@. Line_admittancemat*(2*Vn*Vm*sin(Δθ)))[t,:,c]
        Pa_Iθ[t,c][indexfb]=term
        Pa_Iθ[t,c][indextb]=-term
        SENStemp=SENS[t,c]
        # SENStempV=copy(SENStemp)
        # @. SENSDict["S_IP"][t,c]= (Pa_IV[t,c].*SENStemp[nbus.+ (1:nbus),Pbus]+ Pa_IV[t,c][:,Pbus].*SENStemp[nbus.+tb,Pbus]+ Pa_Iθ[t,c][:,Pbus]*SENStemp[fb,Pbus]+ Pa_Iθ[t,c][:,Pbus]*SENStemp[tb,Pbus]);
        # @. SENSDict["S_IQ"][t,c]= (Pa_IV[t,c][:,Qbus]*SENStemp[nbus.+fb,nbus.+Qbus]+ Pa_IV[t,c][:,Qbus]*SENStemp[nbus.+tb,nbus.+Qbus]+ Pa_Iθ[t,c][:,Qbus]*SENStemp[fb,nbus.+Qbus]+ Pa_Iθ[t,c][:,Qbus]*SENStemp[tb,nbus.+Qbus]);
        # @. SENSDict["S_IP"][t,c][exc[c],:]=0
        # @. SENSDict["S_IQ"][t,c][exc[c],:]=0

        SENSDict["S_IP"][t,c]= (Pa_IV[t,c]*SENStemp[nbus+1:end,Pbus]+ Pa_Iθ[t,c]*SENStemp[1:nbus,Pbus]);
        SENSDict["S_IQ"][t,c]= (Pa_IV[t,c]*SENStemp[nbus+1:end,nbus.+Qbus]+ Pa_Iθ[t,c]*SENStemp[1:nbus,nbus.+Qbus]);
        @. SENSDict["S_IP"][t,c][exc[c],:]=0
        @. SENSDict["S_IQ"][t,c][exc[c],:]=0

#         # SENSDict["S_IP"][t,c][l,i]= (Pa_IV[t,c][l,sb]*SENStemp[nbus+sb,i]+ Pa_IV[t,c][l,rb]*SENStemp[nbus+rb,i]+ Pa_Iθ[t,c][l,sb]*SENStemp[sb,i]+ Pa_Iθ[t,c][l,rb]*SENStemp[rb,i]);
#         # SENSDict["S_IQ"][t,c][l,i]= (Pa_IV[t,c][l,sb]*SENStemp[nbus+sb,i+nbus]+ Pa_IV[t,c][l,rb]*SENStemp[nbus+rb,i+nbus]+ Pa_Iθ[t,c][l,sb]*SENStemp[sb,nbus+i]+ Pa_Iθ[t,c][l,rb]*SENStemp[rb,nbus+i]);
#         # sb=fb[5]
#         # rb=tb[5]
#         # l=5
#         # i=3

#         # lo_thre=[ quantile(SENStempV[i,:],thresholdpercent) for i in 1:size(SENStempV,1)]
#         # up_thre=[ quantile(SENStempV[i,:],1-thresholdpercent) for i in 1:size(SENStempV,1)]
#         # mean_thre=[ mean(abs.(SENStempV[i,:])) for i in 1:size(SENStempV,1)]
#         # histogram(SENStempV[60,:])
#         # lo_thre=ifelse.((lo_thre.<0), lo_thre, 0.0)
#         # up_thre=ifelse.((up_thre.>0), up_thre, 0.0)
#         lo_thre=zeros(size(SENStempV,1))
#         up_thre=zeros(size(SENStempV,1))

#         for i in 1:size(SENStempV,1)
#             if isempty([ctdall;ctdall][i,"gen"])
#                 SENStempVpos=filter(x-> x>0,SENStempV[i,:])
#                 SENStempVneg=filter(x-> x<0,SENStempV[i,:])
#                 if length(SENStempVpos)>7
#                     # up_thre[i]=mean(SENStempVpos)
#                     up_thre[i]=quantile(SENStempVpos,1-thresholdpercent)
#                 end
#                 if length(SENStempVneg)>7
#                     # lo_thre[i]=mean(SENStempVneg)
#                     lo_thre[i]=quantile(SENStempVneg,thresholdpercent)
#                 end
#             else
#                 lo_thre[i]=0
#                 up_thre[i]=0
#             end
#         end
#          SENStempV .= ifelse.((SENStempV .<= lo_thre .|| SENStempV .>= up_thre), SENStempV, 0.0)

#         # SENStempV .= ifelse.(abs.(SENStempV) .> threshold, SENStempV, 0.0)

#         # lo_thre[icondens]=0;up_thre[icondens]=0;
#         # lo_thre[icondens+nbus]=0;up_thre[icondens+nbus]=0;
# #
# #
#         # SENStempV[50,:]
#         # SENStempV[49,:]
#         # SENStempV[109,:]
#         # SENStempV[110,:]
#         # SENStempV2=copy(SENStempV)
#         # SENStempV .= ifelse.((SENStempV .<= lo_thre .|| SENStempV .>= up_thre), SENStempV, 0.0)

#         # SENStempV .= ifelse.((SENStempV .<= quantile(lo_thre,0.9) .|| SENStempV .>= quantile(up_thre,0.1) .|| SENStempV .<= lo_thre .||SENStempV .>= up_thre), SENStempV, 0.0)

#         # findall(x-> x==1,SENStempV .<= lo_thre .|| SENStempV .>= up_thre)
#         # findall(x-> x==0,SENStempV2)
#         # findall(x-> x==0,SENStempV)
#         # findall(x-> x==0,lo_thre)

        SENSDict["S_VP"][t,c]= SENStemp[nbus+1:end,Pbus]
        SENSDict["S_VQ"][t,c]= SENStemp[nbus+1:end,nbus.+Qbus]
        SENSDict["S_θP"][t,c]= SENStemp[1:nbus,Pbus]
        SENSDict["S_θQ"][t,c]= SENStemp[1:nbus,nbus.+Qbus]
        # for key in collect(keys(SENSDict)) #["S_VP" "S_VQ" "S_θP" "S_θQ"]
        #     # Broadcast comparison and set elements less than threshold to 0
        #     SENSDict[key][t,c] .= ifelse.(abs.(SENSDict[key][t,c]) .> threshold, SENSDict[key][t,c], 0.0)
        # end
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
    return SENSDict

end
function Cal_SENS_new(Lininit::Dict{String,Any},ctdall,nTP::Int64,raw::Dict{String, Any};C::String="N")
    nbus::Int64=raw["nBus"];nline::Int64=raw["nLine"];ntrans::Int64=raw["nTrans"];nPbus::Int64=raw["nPActBus"];nQbus::Int64=raw["nQActBus"];
    icondens=50;
    Pbus=raw["PActBus"];Qbus=raw["QActBus"];
    if C=="N"
        CC=1;
        V_0::Array{Float64}=Lininit["V"];θ_0::Array{Float64}=Lininit["θ"];
        # V_0=Lininit["V"];θ_0=Lininit["θ"];
        exc=[[]]
        Jacob::Array{SparseMatrixCSC{Float64, Int64}}=Lininit["NJacob"];
    elseif C=="CL"
        CC::Int64=raw["nConL"];
        V_0=Lininit["V_c_Line"];θ_0=Lininit["θ_c_Line"];#Union{Array{Float64, 3}, Matrix{Float64},Float64}
        exc=raw["ConL"]."#";
        Jacob=Lininit["CJacob"];#::Union{Matrix{Matrix{Float64}},Vector{Matrix{Float64}}}
    end
    TR=ones(Float64,nline+ntrans)
    if ntrans!=0 TR[nline+1:end]= raw["Trans"]."ratio" end
    pa_block=[spzeros(Float64, nline + ntrans, nbus) for _ in 1:nTP, _ in 1:CC]
    Pa_IV=deepcopy(pa_block);Pa_Iθ=deepcopy(pa_block);

    Names=["From","To","g (pu)","b (pu)","g_sh (pu)","b_sh (pu)"]
    df::DataFrame=vcat(raw["Line"][:,Names],raw["Trans"][:,Names])
    # df=raw["Line"];
    g::Vector{Float64}=df."g (pu)";b::Vector{Float64}=df."b (pu)";
    fb::Vector{Int64}=df."From";tb::Vector{Int64}=df."To";
    Δθ=@. θ_0[:,fb,:]-θ_0[:,tb,:];
    Vn=V_0[:,fb,:];Vm=V_0[:,tb,:];
    term=@. 2*Vn-2*Vm*cos(Δθ)
    Line_admittance= @. g^2+b^2
    Line_admittancemat=(Line_admittance')[ones(Int8,nTP),:,ones(Int8,CC)]
    indexfb=CartesianIndex.(1:nline+ntrans, fb)
    indextb=CartesianIndex.(1:nline+ntrans, tb)

    pa_block=[spzeros(Float64, nline + ntrans, nbus) for _ in 1:nTP, _ in 1:CC]


    pa_blockP=[zeros(Float64, nline + ntrans, nPbus) for _ in 1:nTP, _ in 1:CC]
    pa_block2P=[zeros(Float64, nPbus, nPbus) for _ in 1:nTP, _ in 1:CC]
    pa_blockQ=[zeros(Float64, nline + ntrans, nQbus) for _ in 1:nTP, _ in 1:CC]
    pa_block2Q=[zeros(Float64, nQbus, nQbus) for _ in 1:nTP, _ in 1:CC]

    SENSDict = Dict(
        "S_IP" => [zeros(Float64, nline + ntrans, nPbus) for _ in 1:nTP, _ in 1:CC],
        "S_IQ" => [zeros(Float64, nline + ntrans, nQbus) for _ in 1:nTP, _ in 1:CC],
        "S_VP" => [zeros(Float64, nPbus, nPbus) for _ in 1:nTP, _ in 1:CC],
        "S_VQ" => [zeros(Float64, nQbus, nQbus) for _ in 1:nTP, _ in 1:CC],
        "S_θP" => [zeros(Float64, nPbus, nPbus) for _ in 1:nTP, _ in 1:CC],
        "S_θQ" => [zeros(Float64, nQbus, nQbus) for _ in 1:nTP, _ in 1:CC]
    )
    # thresholdpercent = Lininit["FilterIndex"][1]

    Threads.@threads for t=1:nTP
        for c=1:CC

            Pa_IV[t,c][indexfb]=(@. Line_admittancemat*((2*Vn)-2*Vm.*cos.(Δθ)))[t,:,c]
            Pa_IV[t,c][indextb]=(@. Line_admittancemat*((2*Vm)-2*Vn*cos(Δθ)))[t,:,c]
            term= (@. Line_admittancemat*(2*Vn*Vm*sin(Δθ)))[t,:,c]
            Pa_Iθ[t,c][indexfb]=term
            Pa_Iθ[t,c][indextb]=-term

            Jcb=Jacob[t,c]
            fJcb=lu(Jcb)
            n=size(Jcb, 1)
            I_sub_Pbus = I(n)[:, Pbus]
            I_sub_Qbus = I(n)[:, nbus.+Qbus]
            SENStemp_Pbus = fJcb \ I_sub_Pbus
            SENStemp_Qbus = fJcb \ I_sub_Qbus

            # SENSDict["S_IP"][t,c]= (Pa_IV[t,c]*SENStemp[nbus+1:end,Pbus]+ Pa_Iθ[t,c]*SENStemp[1:nbus,Pbus]);
            mul!(SENSDict["S_IP"][t,c], Pa_IV[t,c], SENStemp_Pbus[nbus+1:end,:],1,0)
            mul!(SENSDict["S_IP"][t,c], Pa_Iθ[t,c], SENStemp_Pbus[1:nbus,:], 1,1)
            mul!(SENSDict["S_IQ"][t,c], Pa_IV[t,c], SENStemp_Qbus[nbus+1:end,:],1,0)
            mul!(SENSDict["S_IQ"][t,c], Pa_Iθ[t,c], SENStemp_Qbus[1:nbus,:], 1,1)

            # SENSDict["S_IQ"][t,c]= (Pa_IV[t,c]*SENStemp[nbus+1:end,nbus.+Qbus]+ Pa_Iθ[t,c]*SENStemp[1:nbus,nbus.+Qbus]);
            @. SENSDict["S_IP"][t,c][exc[c],:]=0
            @. SENSDict["S_IQ"][t,c][exc[c],:]=0

            SENSDict["S_VP"][t,c]= SENStemp_Pbus[nbus+1:end,:]
            SENSDict["S_VQ"][t,c]= SENStemp_Qbus[nbus+1:end,:]
            SENSDict["S_θP"][t,c]= SENStemp_Pbus[1:nbus,:]
            SENSDict["S_θQ"][t,c]= SENStemp_Qbus[1:nbus,:]
        end
    end
    return SENSDict

end

function Cal_SENS_new2(Lininit::Dict{String,Any}, ctdall, nTP::Int64, raw::Dict{String, Any}; C::String="N")
    # --- unpack / init ---
    nbus::Int64 = raw["nBus"]
    nline::Int64 = raw["nLine"]
    ntrans::Int64 = raw["nTrans"]
    nPbus::Int64 = raw["nPActBus"]
    nQbus::Int64 = raw["nQActBus"]
    Pbus = raw["PActBus"]
    Qbus = raw["QActBus"]

    if C == "N"
        CC = 1
        V_0 = Lininit["V"]
        θ_0 = Lininit["θ"]
        exc = [[]]
        Jacob = Lininit["NJacob"]
    elseif C == "CL"
        CC = raw["nConL"]
        V_0 = Lininit["V_c_Line"]
        θ_0 = Lininit["θ_c_Line"]
        exc = raw["ConL"]."#"
        Jacob = Lininit["CJacob"]
    else
        error("Unsupported C = $C")
    end

    # build some constants / sizes
    TR = ones(Float64, nline + ntrans)
    if ntrans != 0
        TR[nline+1:end] .= raw["Trans"]."ratio"
    end

    Names = ["From","To","g (pu)","b (pu)","g_sh (pu)","b_sh (pu)"]
    df::DataFrame = vcat(raw["Line"][:,Names], raw["Trans"][:,Names])
    g = df."g (pu)"
    b = df."b (pu)"
    fb = df."From"
    tb = df."To"

    # precomputed voltage/angle differences (shaped nTP x nLine x CC or similar)
    Δθ = @. θ_0[:, fb, :]
    Δθ .= θ_0[:, fb, :] .- θ_0[:, tb, :]   # Δθ[t, line, c] = θ_0[t,fb,c] - θ_0[t,tb,c]
    Vn = V_0[:, fb, :]
    Vm = V_0[:, tb, :]

    # note: line admittance scalar per line -> replicate into a 3D array for indexing convenience
    Line_admittance = @. g.^2 .+ b.^2
    # make a (nTP x nLine x CC) "repeated" view (cheap)
    Line_admittancemat = reshape(repeat(Line_admittance', outer=(nTP,1,CC)), (nTP, length(Line_admittance), CC))

    # Cartesian indices for assigning into sparse Pa arrays
    indexfb = CartesianIndex.(1:(nline + ntrans), fb)
    indextb = CartesianIndex.(1:(nline + ntrans), tb)

    # --- allocate Pa_IV / Pa_Iθ (one unique sparse matrix per (t,c)) ---
    pa_template = spzeros(Float64, nline + ntrans, nbus)
    Pa_IV = [copy(pa_template) for _ in 1:nTP, _ in 1:CC]
    Pa_Iθ = [copy(pa_template) for _ in 1:nTP, _ in 1:CC]

    # --- Preallocate SENSDict entries (unique matrix per (t,c)) ---
    SENSDict = Dict(
        "S_IP" => [zeros(Float64, nline + ntrans, nPbus) for _ in 1:nTP, _ in 1:CC],
        "S_IQ" => [zeros(Float64, nline + ntrans, nQbus) for _ in 1:nTP, _ in 1:CC],
        "S_VP" => [zeros(Float64, nPbus, nPbus) for _ in 1:nTP, _ in 1:CC],
        "S_VQ" => [zeros(Float64, nQbus, nQbus) for _ in 1:nTP, _ in 1:CC],
        "S_θP" => [zeros(Float64, nPbus, nPbus) for _ in 1:nTP, _ in 1:CC],
        "S_θQ" => [zeros(Float64, nQbus, nQbus) for _ in 1:nTP, _ in 1:CC]
    )

    # --- threaded compute ---
    Threads.@threads for t in 1:nTP
        # local aliases for speed (avoids global lookups)
        LineA_t = Line_admittancemat  # we index it with [t, :, c]
        Vn_t = Vn
        Vm_t = Vm
        Δθ_t = Δθ
        nbus_local = nbus

        # local temporaries used inside inner loop to reduce allocations
        # note: these are small dense vectors for a single line-index slice
        tmp_linevals = similar(Line_admittance, length(Line_admittance))  # reused for elementwise multiplications

        for c in 1:CC
            # --- compute Pa_IV and Pa_Iθ row values for time t, cond c ---
            # We'll compute the dense vector of values for all lines then assign into sparse matrix entries
            @inbounds begin
                # values for "from" side
                @. tmp_linevals = LineA_t[t, :, c] * ((2*Vn_t)[t, :, c] - 2*(Vm_t)[t, :, c]*cos( (Δθ_t)[t, :, c] ))
                for (idx, val) in enumerate(tmp_linevals)     # assign into sparse matrix using CartesianIndex
                    Pa_IV[t,c][CartesianIndex(idx, fb[idx])] = val
                end

                # values for "to" side (other term)
                @. tmp_linevals = LineA_t[t, :, c] * ((2*Vm_t)[t, :, c] - 2*(Vn_t)[t, :, c]*cos((Δθ_t)[t, :, c]))
                for (idx, val) in enumerate(tmp_linevals)
                    Pa_IV[t,c][CartesianIndex(idx, tb[idx])] = val
                end

                # term for angle part
                @. tmp_linevals = LineA_t[t, :, c] * (2 * (Vn_t)[t, :, c] * (Vm_t)[t, :, c] * sin((Δθ_t)[t, :, c]))
                for (idx, val) in enumerate(tmp_linevals)
                    Pa_Iθ[t,c][CartesianIndex(idx, fb[idx])] = val
                    Pa_Iθ[t,c][CartesianIndex(idx, tb[idx])] = -val
                end
            end

            # --- solve linear system to get SENStemp matrices ---
            # If Jacob[t,c] entries are views into a bigger structure, uncomment copy() to be safe:
            Jcb = Jacob[t,c]
            # Jcb = copy(Jacob[t,c])   # <--- uncomment only if Jacob entries are not independent
            fJcb = lu(Jcb)
            n = size(Jcb, 1)

            I_sub_Pbus = I(n)[:, Pbus]
            I_sub_Qbus = I(n)[:, nbus_local .+ Qbus]

            # SENStemp_* are (n x nPbus) and (n x nQbus)
            SENStemp_Pbus = fJcb \ I_sub_Pbus
            SENStemp_Qbus = fJcb \ I_sub_Qbus

            # --- compute S_IP and S_IQ via mul! directly into the preallocated unique matrices ---
            S_IP_mat = SENSDict["S_IP"][t,c]   # unique preallocated matrix
            S_IQ_mat = SENSDict["S_IQ"][t,c]

            # multiply into them (zero-copy, in-place)
            # S_IP = Pa_IV * SENStemp_Pbus[nbus+1:end, :]  +  Pa_Iθ * SENStemp_Pbus[1:nbus, :]
            @inbounds begin
                mul!(S_IP_mat, Pa_IV[t,c], view(SENStemp_Pbus, nbus_local+1:size(SENStemp_Pbus,1), :), 1, 0)
                mul!(S_IP_mat, Pa_Iθ[t,c], view(SENStemp_Pbus, 1:nbus_local, :), 1, 1)

                mul!(S_IQ_mat, Pa_IV[t,c], view(SENStemp_Qbus, nbus_local+1:size(SENStemp_Qbus,1), :), 1, 0)
                mul!(S_IQ_mat, Pa_Iθ[t,c], view(SENStemp_Qbus, 1:nbus_local, :), 1, 1)
            end

            # apply exclusions (mutates only the unique (t,c) matrices)
            @inbounds @. S_IP_mat[exc[c], :] = 0
            @inbounds @. S_IQ_mat[exc[c], :] = 0

            # store the SENStemp sub-blocks (these are small; keep as copies)
            SENSDict["S_VP"][t,c] = copy(view(SENStemp_Pbus, nbus_local+1:size(SENStemp_Pbus,1), :))
            SENSDict["S_VQ"][t,c] = copy(view(SENStemp_Qbus, nbus_local+1:size(SENStemp_Qbus,1), :))
            SENSDict["S_θP"][t,c] = copy(view(SENStemp_Pbus, 1:nbus_local, :))
            SENSDict["S_θQ"][t,c] = copy(view(SENStemp_Qbus, 1:nbus_local, :))
        end
    end

    return SENSDict
end


# function Cal_SENS_Cascade(Lininit::Dict{String,Any}, ctdall, nTP::Int64, raw::Dict{String, Any})
function Cal_SENS_Cascade(Jcbs,lineind,t,raw::Dict{String, Any}, Lininit::Dict{String,Any},reprow)
    # lineind=[j];t=i;
    # lineind=casreps[t,5][1,"LL"];Jcb=casreps[t,5][1,"Jcbs"]
    # --- unpack / init ---
    # reprow."busset"
    # raw["nBus"]
    nbus::Int64 = Int64(size(Jcbs,1)/2)
    nbus_org::Int64 =raw["nBus"]
    nline::Int64 = size(lineind,1)
    ntrans::Int64 = 0
    # nPbus::Int64 = raw["nPActBus"]
    # nQbus::Int64 = raw["nQActBus"]
    # Pbus = raw["PActBus"]
    # Qbus = raw["QActBus"]

    Pbus=findall(x-> x in raw["PActBus"] ,reprow."busset")
    Qbus=findall(x-> x in raw["QActBus"] ,reprow."busset")
    nPbus::Int64 = length(Pbus)
    nQbus::Int64 = length(Qbus)
    Pbusindex=findall(x-> x in reprow."busset" ,raw["PActBus"])
    Qbusindex=findall(x-> x in reprow."busset" ,raw["QActBus"])


    # CC = 1
    V_0 = zeros(Float64, nbus_org)
    θ_0 = zeros(Float64, nbus_org)
    V_0[reprow."busset"] .= reprow."Vmag"
    θ_0[reprow."busset"] .= reprow."Vang"
    TR = ones(Float64, nline + ntrans)
    if ntrans != 0
        TR[nline+1:end] .= raw["Trans"]."ratio"
    end

    Names = ["From","To","g (pu)","b (pu)","g_sh (pu)","b_sh (pu)"]
    # df::DataFrame = vcat(raw["Line"][:,Names], raw["Trans"][:,Names])
    df::DataFrame = raw["Line"][lineind,Names]
    g = df."g (pu)"
    b = df."b (pu)"
    fb = df."From"
    tb = df."To"

    # precomputed voltage/angle differences (shaped nTP x nLine x CC or similar)
    # Δθ = @. θ_0[t, fb, 1]
    # Δθ .= θ_0[t, fb, 1] .- θ_0[t, tb, 1]   # Δθ[t, line, c] = θ_0[t,fb,c] - θ_0[t,tb,c]
    # Vn = V_0[t, fb, 1]
    # Vm = V_0[t, tb, 1]

    Δθ = @. θ_0[fb]
    Δθ .= θ_0[fb] .- θ_0[tb]   # Δθ[t, line, c] = θ_0[t,fb,c] - θ_0[t,tb,c]
    Vn = V_0[fb]
    Vm = V_0[tb]

    # note: line admittance scalar per line -> replicate into a 3D array for indexing convenience
    Line_admittance = @. g.^2 .+ b.^2
    # --- allocate Pa_IV / Pa_Iθ (one unique sparse matrix per (t,c)) ---
    pa_template = spzeros(Float64, nline + ntrans, nbus_org)
    Pa_IV = copy(pa_template)
    Pa_Iθ = copy(pa_template)

    tmp_linevals = similar(Line_admittance, length(Line_admittance))  # reused for elementwise multiplications


    # --- compute Pa_IV and Pa_Iθ row values
    # We'll compute the dense vector of values for all lines then assign into sparse matrix entries
    @inbounds begin
        # values for "from" side
        @. tmp_linevals = Line_admittance * ((2*Vn) - 2*(Vm)*cos((Δθ)))
        for (idx, val) in enumerate(tmp_linevals)     # assign into sparse matrix using CartesianIndex
            Pa_IV[CartesianIndex(idx, fb[idx])] = val
        end

        # values for "to" side (other term)
        @. tmp_linevals = Line_admittance * ((2*Vm) - 2*(Vn)*cos((Δθ)))
        for (idx, val) in enumerate(tmp_linevals)
            Pa_IV[CartesianIndex(idx, tb[idx])] = val
        end

        # term for angle part
        @. tmp_linevals = Line_admittance * (2 * (Vn) * (Vm) * sin((Δθ)))
        for (idx, val) in enumerate(tmp_linevals)
            Pa_Iθ[CartesianIndex(idx, fb[idx])] = val
            Pa_Iθ[CartesianIndex(idx, tb[idx])] = -val
        end
    end

    fJcb = lu(Jcbs)
    n = size(Jcbs, 1)

    I_sub_Pbus = I(n)[:, Pbus]
    I_sub_Qbus = I(n)[:, nbus .+ Qbus]

    # SENStemp_* are (n x nPbus) and (n x nQbus)
    SENStemp_Pbus = fJcb \ I_sub_Pbus
    SENStemp_Qbus = fJcb \ I_sub_Qbus

    # --- compute S_IP and S_IQ via mul! directly into the preallocated unique matrices ---
    S_IP_mat = zeros(Float64, nline + ntrans, raw["nPActBus"])
    S_IQ_mat = zeros(Float64, nline + ntrans, raw["nQActBus"])

    S_IP_mat_slice = @view S_IP_mat[:, Pbusindex]
    S_IQ_mat_slice = @view S_IQ_mat[:, Qbusindex]
    Pa_IV_slice = @view Pa_IV[:, reprow."busset"]
    Pa_Iθ_slice = @view Pa_Iθ[:, reprow."busset"]
    @inbounds begin
        mul!(S_IP_mat_slice, Pa_IV_slice, view(SENStemp_Pbus, nbus+1:size(SENStemp_Pbus,1), :), 1, 0)
        mul!(S_IP_mat_slice, Pa_Iθ_slice, view(SENStemp_Pbus, 1:nbus, :), 1, 1)

        mul!(S_IQ_mat_slice, Pa_IV_slice, view(SENStemp_Qbus, nbus+1:size(SENStemp_Qbus,1), :), 1, 0)
        mul!(S_IQ_mat_slice, Pa_Iθ_slice, view(SENStemp_Qbus, 1:nbus, :), 1, 1)
    end
    

    # hcat(S_IP_mat,S_IQ_mat)
    return hcat(S_IP_mat,S_IQ_mat)
end
# function Cal_SENS_new(Lininit::Dict{String,Any}, ctdall, nTP::Int64, raw::Dict{String, Any}; C::String="N")
#     # Extract system parameters
#     nbus, nline, ntrans = raw["nBus"], raw["nLine"], raw["nTrans"]
#     nPbus, nQbus = raw["nPActBus"], raw["nQActBus"]
#     Pbus, Qbus = raw["PActBus"], raw["QActBus"]

#     # Determine contingency case
#     if C == "N"
#         CC = 1
#         V_0, θ_0 = Lininit["V"], Lininit["θ"]
#         exc = [[]]
#         SENS = Lininit["NSens"]
#     elseif C == "CL"
#         CC = raw["nConL"]
#         V_0, θ_0 = Lininit["V_c_Line"], Lininit["θ_c_Line"]
#         exc = raw["ConL"]["#"]
#         SENS = Lininit["CSens"]
#     end

#     # Transformer ratio
#     TR = ones(Float64, nline + ntrans)
#     if ntrans != 0
#         TR[nline+1:end] .= raw["Trans"]["ratio"]
#     end

#     # Sparse block allocations (avoid shared references)
#     Pa_IV  = [spzeros(Float64, nline + ntrans, nbus) for _ in 1:nTP, _ in 1:CC]
#     Pa_Iθ  = [spzeros(Float64, nline + ntrans, nbus) for _ in 1:nTP, _ in 1:CC]
#     pa_blockP  = [spzeros(Float64, nline + ntrans, nPbus) for _ in 1:nTP, _ in 1:CC]
#     pa_block2P = [spzeros(Float64, nPbus, nPbus) for _ in 1:nTP, _ in 1:CC]
#     pa_blockQ  = [spzeros(Float64, nline + ntrans, nQbus) for _ in 1:nTP, _ in 1:CC]
#     pa_block2Q = [spzeros(Float64, nQbus, nQbus) for _ in 1:nTP, _ in 1:CC]

#     # Combine blocks into dictionary
#     SENSDict = Dict(
#         "S_IP"  => pa_blockP,
#         "S_IQ"  => pa_blockQ,
#         "S_VP"  => pa_block2P,
#         "S_VQ"  => pa_block2Q,
#         "S_θP"  => pa_block2P,
#         "S_θQ"  => pa_block2Q
#     )

#     # Extract line parameters
#     Names = ["From", "To", "g (pu)", "b (pu)", "g_sh (pu)", "b_sh (pu)"]
#     df = vcat(raw["Line"][:, Names], raw["Trans"][:, Names])
#     g::Vector{Float64} = df[!, "g (pu)"]; b::Vector{Float64} = df[!, "b (pu)"];
#     fb::Vector{Int64} = df[!, "From"]; tb::Vector{Int64} = df[!, "To"];

#     # Compute line parameters
#     Δθ = @. θ_0[:, fb, :] - θ_0[:, tb, :]
#     Vn, Vm = V_0[:, fb, :], V_0[:, tb, :]
#     Line_admittance = @. g^2 + b^2
#     Line_admittancemat = (Line_admittance')[ones(Int8, nTP), :, ones(Int8, CC)]

#     indexfb = CartesianIndex.(1:nline+ntrans, fb)
#     indextb = CartesianIndex.(1:nline+ntrans, tb)

#     thresholdpercent = Lininit["FilterIndex"][1]

#     # Compute sensitivity matrices
#     for t = 1:nTP, c = 1:CC
#         Pa_IV[t, c][indexfb] = (@. Line_admittancemat * ((2 * Vn) - 2 * Vm * cos.(Δθ)))[t, :, c]
#         Pa_IV[t, c][indextb] = (@. Line_admittancemat * ((2 * Vm) - 2 * Vn * cos.(Δθ)))[t, :, c]
#         term = (@. Line_admittancemat * (2 * Vn * Vm * sin.(Δθ)))[t, :, c]
#         Pa_Iθ[t, c][indexfb] = term
#         Pa_Iθ[t, c][indextb] = -term

#         SENStemp = SENS[t, c]

#         SENSDict["S_IP"][t, c] = Pa_IV[t, c] * SENStemp[nbus+1:end, Pbus] + Pa_Iθ[t, c] * SENStemp[1:nbus, Pbus]
#         SENSDict["S_IQ"][t, c] = Pa_IV[t, c] * SENStemp[nbus+1:end, nbus .+ Qbus] + Pa_Iθ[t, c] * SENStemp[1:nbus, nbus .+ Qbus]

#         @. SENSDict["S_IP"][t, c][exc[c], :] = 0
#         @. SENSDict["S_IQ"][t, c][exc[c], :] = 0

#         SENSDict["S_VP"][t, c] = SENStemp[nbus+1:end, Pbus]
#         SENSDict["S_VQ"][t, c] = SENStemp[nbus+1:end, nbus .+ Qbus]
#         SENSDict["S_θP"][t, c] = SENStemp[1:nbus, Pbus]
#         SENSDict["S_θQ"][t, c] = SENStemp[1:nbus, nbus .+ Qbus]
#     end

#     return SENSDict
# end



# function Cal_SENS_Simple_bad(Lininit::Dict{String,Any},ctdall,nTP::Int64,raw::Dict{String, Any};C::String="N")
#     nbus::Int64=raw["nBus"];nline::Int64=raw["nLine"];ntrans::Int64=raw["nTrans"];nPbus::Int64=raw["nPActBus"];nQbus::Int64=raw["nQActBus"];
#     icondens=50;
#     Pbus=raw["PActBus"];Qbus=raw["QActBus"];
#     if C=="N"
#         CC=1;
#         V_0::Array{Float64}=Lininit["V"];θ_0::Array{Float64}=Lininit["θ"];
#         exc=[[]]
#         SENS::Array{Array{Float64}}=Lininit["NSens"];
#     elseif C=="CL"
#         CC::Int64=raw["nConL"];
#         V_0=Lininit["V_c_Line"];θ_0=Lininit["θ_c_Line"];#Union{Array{Float64, 3}, Matrix{Float64},Float64}
#         exc=raw["ConL"]."#";
#         SENS=Lininit["CSens"];#::Union{Matrix{Matrix{Float64}},Vector{Matrix{Float64}}}
#     end
#     TR=ones(Float64,nline+ntrans)
#     if ntrans!=0 TR[nline+1:end]= raw["Trans"]."ratio" end
#     block=spzeros(Float64,nline+ntrans,nbus)
#     pa_block=fill(block,nTP,CC)
#     block2=spzeros(Float64,nbus,nbus)
#     pa_block2=fill(block2,nTP,CC)

#     Pa_IV=copy(pa_block);Pa_Iθ=copy(pa_block);

#     Names=["From","To","g (pu)","b (pu)","g_sh (pu)","b_sh (pu)"]
#     df::DataFrame=vcat(raw["Line"][:,Names],raw["Trans"][:,Names])
#     # df=raw["Line"];
#     g::Vector{Float64}=df."g (pu)";b::Vector{Float64}=df."b (pu)";
#     fb::Vector{Int64}=df."From";tb::Vector{Int64}=df."To";
#     Δθ=@. θ_0[:,fb,:]-θ_0[:,tb,:];
#     Vn=V_0[:,fb,:];Vm=V_0[:,tb,:];
#     term=@. 2*Vn-2*Vm*cos(Δθ)
#     Line_admittance= @. g^2+b^2
#     Line_admittancemat=(Line_admittance')[ones(Int8,nTP),:,ones(Int8,CC)]
#     indexfb=CartesianIndex.(1:nline+ntrans, fb)
#     indextb=CartesianIndex.(1:nline+ntrans, tb)

#     blockP=spzeros(Float64,nline+ntrans,nPbus)
#     pa_blockP=fill(copy(blockP),nTP,CC)
#     block2P=spzeros(Float64,nPbus,nPbus)
#     pa_block2P=fill(copy(block2P),nTP,CC)

#     blockQ=spzeros(Float64,nline+ntrans,nQbus)
#     pa_blockQ=fill(copy(blockQ),nTP,CC)
#     block2Q=spzeros(Float64,nQbus,nQbus)
#     pa_block2Q=fill(copy(block2Q),nTP,CC)

#     SENSDict = Dict(
#         "S_IP" => copy(pa_blockP),
#         "S_IQ" => copy(pa_blockQ),
#         "S_VP" => copy(pa_block2P),
#         "S_VQ" => copy(pa_block2Q),
#         "S_θP" => copy(pa_block2P),
#         "S_θQ" => copy(pa_block2Q)
#     )
#     thresholdpercent = Lininit["FilterIndex"][1]
#     # ,c=1:CC
#     term_pre=@. Line_admittancemat*(2*Vn*Vm*sin(Δθ))
#     PA_IV_pre1=@. Line_admittancemat*((2*Vn)-2*Vm.*cos.(Δθ))
#     PA_IV_pre2=@. Line_admittancemat*((2*Vm)-2*Vn*cos(Δθ))

#     instances=[(t,c) for t=1:nTP, c=1:CC]
#     for (t,c) in instances

#             Pa_IV[t,c][indexfb]=(PA_IV_pre1)[t,:,c]
#             Pa_IV[t,c][indextb]=(PA_IV_pre2)[t,:,c]
#             term= (term_pre)[t,:,c]
#             Pa_Iθ[t,c][indexfb]=term
#             Pa_Iθ[t,c][indextb]=-term
#             SENStemp=SENS[t,c]
#             SENSDict["S_VP"][t,c]= SENStemp[nbus+1:end,Pbus]
#             SENSDict["S_VQ"][t,c]= SENStemp[nbus+1:end,nbus.+Qbus]
#             SENSDict["S_θP"][t,c]= SENStemp[1:nbus,Pbus]
#             SENSDict["S_θQ"][t,c]= SENStemp[1:nbus,nbus.+Qbus]
#             Pa_Iθ_slot=Pa_Iθ[t,c]
#             Pa_IV_slot=Pa_IV[t,c]
#             SENSDict["S_IP"][t,c]= (Pa_IV_slot*SENStemp[nbus+1:end,Pbus]+ Pa_Iθ_slot*SENStemp[1:nbus,Pbus]);
#             SENSDict["S_IQ"][t,c]= (Pa_IV_slot*SENStemp[nbus+1:end,nbus.+Qbus]+ Pa_Iθ_slot*SENStemp[1:nbus,nbus.+Qbus]);
#             @. SENSDict["S_IP"][t,c][exc[c],:]=0
#             @. SENSDict["S_IQ"][t,c][exc[c],:]=0
#     end
#     return SENSDict

# end

function Power_Inj(nTP,raw,ctdall,xij,xji,xijT,xjiT,LinMethod;C="N")
    # nTP,raw,xij,xji,xijT,xjiT=nTP,raw,Nqij_L,Nqji_L,Nqij_T,Nqji_T
    # xij,xji,xijT,xjiT=CLpij_L,CLpji_L,CLpij_T,CLpji_T
    # temp=Model()
    if C=="N"
        CC=1;
    elseif C=="CL"
        CC=raw["nConL"];
    end

    ctdlineF = @view ctdall."linef"[:]
    ctdlineT = @view ctdall."linet"[:]
    ctdtransF = @view ctdall."transf"[:]
    ctdtransT = @view ctdall."transt"[:]


    # @expression(temp, xinj[t=1:nTP, i=1:raw["nBus"],c=1:CC], sum(xij[t,n,c] for n in ctdlineF[i])
    #                                             + sum(xji[t,n,c] for n in ctdlineT[i])
    #                                             + sum(xijT[t,n,c] for n in ctdtransF[i])
    #                                             + sum(xjiT[t,n,c] for n in ctdtransT[i])
    #
    if LinMethod=="None"
        xinj=zeros(QuadExpr,nTP,raw["nBus"],CC)
    else
        xinj=zeros(AffExpr,nTP,raw["nBus"],CC)
    end
    # sum(xij[:,ctdlineF[i],:] ,dims=2)
    for i in 1:raw["nBus"]
        xinj[:,i,:]=sum(xij[:,ctdlineF[i],:] ,dims=2)+sum(xji[:,ctdlineT[i],:],dims=2)+sum(xijT[:,ctdtransF[i],:],dims=2)+sum(xjiT[:,ctdtransT[i],:],dims=2)
    end
    # @expression(temp, xinj[t=1:nTP, i=1:raw["nBus"],c=1:CC], sum(xij[t,n,c] for n in 1:raw["nLine"] if (raw["Line"]."From"[n]==i))
    #                                             + sum(xji[t,n,c] for n in 1:raw["nLine"] if (raw["Line"]."To"[n]==i))
    #                                             + sum(xijT[t,n,c] for n in 1:raw["nTrans"] if (raw["Trans"]."From"[n]==i))
    #                                             + sum(xjiT[t,n,c] for n in 1:raw["nTrans"] if (raw["Trans"]."To"[n]==i))
    #                                             )
    return xinj
end


function ConstraintPars(LinMethod,SCOPF,Vars,Parlist,nTP,raw,relaxconfac=0.0;States="Norm")
    if LinMethod=="None"
        SCOPF= Non_Linear_CP(SCOPF,Vars,Parlist,nTP,raw,States,relaxconfac)
    elseif LinMethod=="Yang"
        SCOPF=Linear_Yang_CP(SCOPF,Vars,Parlist,nTP,raw,States,relaxconfac)
    elseif LinMethod=="SOCP"
        SCOPF=Linear_SOCP_CP(SCOPF,Vars,Parlist,nTP,raw,States,relaxconfac)
    elseif LinMethod=="SENS"
        # SCOPF= Linear_SENS_CP_new(SCOPF,Vars,Parlist,nTP,raw,States,SENSDict,Lininit)
    end

    return SCOPF
end

function Define_Variable_Normal(SCOPF::Model,nTP::Int64,raw::Dict{String, Any},LinMethod::String,Lininit::Dict{String,Any})
    if LinMethod=="None"
        @variables(
            SCOPF,
            begin
                e[t = 1:nTP, i = 1:raw["nBus"]]
                f[t = 1:nTP, i = 1:raw["nBus"]]
            end
        )
        for  t=1:nTP, i=1:raw["nBus"]
            set_start_value(e[t,i], 1.00)
            set_start_value(f[t,i], 0.00)
        end
        Vars=e,f
    elseif LinMethod=="Yang"
        @variables(
            SCOPF,
            begin
                V2[t = 1:nTP, i = 1:raw["nBus"]]
                θ[t = 1:nTP, i = 1:raw["nBus"]]
            end
        )
        for  t=1:nTP, i=1:raw["nBus"]
            set_start_value(V2[t,i], Lininit["V"][t,i]^2)
            set_start_value(θ[t,i], Lininit["θ"][t,i])
        end
        Vars=V2,θ
    elseif LinMethod== "SOCP"
        @variables(
            SCOPF,
            begin
                V2[t = 1:nTP, i = 1:raw["nBus"]]
                A[t = 1:nTP, i = 1:raw["nLine"]]
                B[t = 1:nTP, i = 1:raw["nLine"]]
            end
        )
        for t=1:nTP, i=1:raw["nBus"]
            set_start_value(V2[t,i], Lininit["V"][t,i]^2)
        end
        A_0=zeros(nTP,raw["nLine"]);B_0=zeros(nTP,raw["nLine"])
        for t=1:nTP, i=1:raw["nLine"]
            A_0[t,i]=Lininit["V"][t,raw["Line"][i,"From"]]*Lininit["V"][t,raw["Line"][i,"To"]]*cos(Lininit["θ"][t,raw["Line"][i,"From"]]-Lininit["θ"][t,raw["Line"][i,"To"]])
            B_0[t,i]=Lininit["V"][t,raw["Line"][i,"From"]]*Lininit["V"][t,raw["Line"][i,"To"]]*sin(Lininit["θ"][t,raw["Line"][i,"From"]]-Lininit["θ"][t,raw["Line"][i,"To"]])
        end
        for t=1:nTP, i=1:raw["nLine"]
            set_start_value(A[t,i], A_0[t,i])
            set_start_value(B[t,i], B_0[t,i])
        end
        Vars=V2,A,B
    elseif LinMethod== "SENS"
        @variables(
            SCOPF,
            begin
                ΔPinjP[t = 1:nTP, i = 1:raw["nPActBus"]]>=0
                ΔPinjN[t = 1:nTP, i = 1:raw["nPActBus"]]>=0
                ΔQinjP[t = 1:nTP, i = 1:raw["nQActBus"]]>=0
                ΔQinjN[t = 1:nTP, i = 1:raw["nQActBus"]]>=0
            end
        )
        set_start_value.(ΔPinjP, 0)
        set_start_value.(ΔPinjN, 0)
        set_start_value.(ΔQinjP, 0)
        set_start_value.(ΔQinjN, 0)
        # fix.(ΔPinj, 0)
        # fix.(ΔQinj,0)

        # for  t=1:nTP, i in setdiff(1:raw["nBus"],union(raw["RES"]."Bus",raw["ESS"]."Bus",raw["Load"]."Bus_N",raw["Gen"]."Bus",raw["Shunt"]."Bus"))
        #     fix(ΔPinj[t,i],0)
        #     fix(ΔQinj[t,i],0)
        # end
        # fix(ΔQinj[1,4],0)
        Vars=ΔPinjP,ΔPinjN,ΔQinjP,ΔQinjN
    end
    return SCOPF,Vars
end

function Define_Variable_Contingency_L(SCOPF::Model,nTP::Int64,raw::Dict{String, Any},LinMethod::String,Lininit::Dict{String,Any})
    if LinMethod=="None"
        @variables(
            SCOPF,
            begin
                e_c_Line[t = 1:nTP, i = 1:raw["nBus"],cL = 1:raw["nConL"]]
                f_c_Line[t = 1:nTP, i = 1:raw["nBus"],cL = 1:raw["nConL"]]
            end
        )
        for  t=1:nTP, i=1:raw["nBus"],cL = 1:raw["nConL"]
            set_start_value(e_c_Line[t,i,cL], 1.00)
            set_start_value(f_c_Line[t,i,cL], 0.00)
        end
        Vars=e_c_Line,f_c_Line
    elseif LinMethod=="Yang"
        @variables(
            SCOPF,
            begin
                V2_c_Line[t = 1:nTP, i = 1:raw["nBus"],cL = 1:raw["nConL"]]
                θ_c_Line[t = 1:nTP, i = 1:raw["nBus"],cL = 1:raw["nConL"]]
            end
        )
        for  t=1:nTP, i=1:raw["nBus"],cL = 1:raw["nConL"]
            set_start_value(V2_c_Line[t,i,cL], Lininit["V_c_Line"][t,i,cL].^2)
            set_start_value(θ_c_Line[t,i,cL], Lininit["θ_c_Line"][t,i,cL])
        end
        Vars=V2_c_Line,θ_c_Line
    elseif LinMethod== "SOCP"
        @variables(
            SCOPF,
            begin
                V2_c_Line[t = 1:nTP, i = 1:raw["nBus"],cL = 1:raw["nConL"]]
                A_c_Line[t = 1:nTP, i = 1:raw["nLine"],cL = 1:raw["nConL"]]
                B_c_Line[t = 1:nTP, i = 1:raw["nLine"],cL = 1:raw["nConL"]]
            end
        )
        for t=1:nTP, i=1:raw["nBus"],cL = 1:raw["nConL"]
            set_start_value(V2_c_Line[t,i,cL], Lininit["V_c_Line"][t,i,cL]^2)
        end
        A_0=zeros(nTP,raw["nLine"],raw["nConL"]);B_0=zeros(nTP,raw["nLine"],raw["nConL"])
        for t=1:nTP, i=1:raw["nLine"] ,cL = 1:raw["nConL"]
            A_0[t,i,cL]=Lininit["V"][t,raw["Line"][i,"From"],cL]*Lininit["V"][t,raw["Line"][i,"To"],cL]*cos(Lininit["θ"][t,raw["Line"][i,"From"],cL]-Lininit["θ"][t,raw["Line"][i,"To"],cL])
            B_0[t,i,cL]=Lininit["V"][t,raw["Line"][i,"From"],cL]*Lininit["V"][t,raw["Line"][i,"To"],cL]*sin(Lininit["θ"][t,raw["Line"][i,"From"],cL]-Lininit["θ"][t,raw["Line"][i,"To"],cL])
        end
        for t=1:nTP, i=1:raw["nLine"],cL = 1:raw["nConL"]
            set_start_value(A_c_Line[t,i,cL], A_0[t,i,cL])
            set_start_value(B_c_Line[t,i,cL], B_0[t,i,cL])
        end
        Vars= V2_c_Line,A_c_Line,B_c_Line
    elseif LinMethod== "SENS"
        @variables(
            SCOPF,
            begin
                ΔPinjP_c_Line[t = 1:nTP, i = 1:raw["nPActBus"],cL = 1:raw["nConL"]]>=0
                ΔPinjN_c_Line[t = 1:nTP, i = 1:raw["nPActBus"],cL = 1:raw["nConL"]]>=0
                ΔQinjP_c_Line[t = 1:nTP, i = 1:raw["nQActBus"],cL = 1:raw["nConL"]]>=0
                ΔQinjN_c_Line[t = 1:nTP, i = 1:raw["nQActBus"],cL = 1:raw["nConL"]]>=0
            end
        )
            set_start_value.(ΔPinjP_c_Line, 0)
            set_start_value.(ΔPinjN_c_Line, 0)
            set_start_value.(ΔQinjP_c_Line, 0)
            set_start_value.(ΔQinjN_c_Line, 0)
            # fix.(ΔPinj_c_Line, 0)
            # fix.(ΔQinj_c_Line, 0)
        # for  t=1:nTP,cL = 1:raw["nConL"], i in setdiff(1:raw["nBus"],union(raw["RES"]."Bus",raw["ESS"]."Bus",raw["Load"]."Bus_N",raw["Gen"]."Bus",raw["Shunt"]."Bus"))
        #     fix(ΔPinj_c_Line[t,i,cL],0)
        #     fix(ΔQinj_c_Line[t,i,cL],0)
        # end
        Vars=ΔPinjP_c_Line,ΔPinjN_c_Line,ΔQinjP_c_Line,ΔQinjN_c_Line
    end
    return SCOPF,Vars
end



function Flow_Expr(LinMethod::String,Lininit::Dict{String,Any},nTP::Int64,raw::Dict{String, Any},Vars::Tuple, TorL::Char; Cont::DataFrame=DataFrame(),inv::Char='t',SENSDict::Dict{String, Array{Float64, 4}}=Dict{String, Array{Float64, 4}}())
    # LinMethod,Lininit,nTP,raw,Vars, TorL=LinMethod,Lininit,nTP,raw,Vars,'L';Cont=DataFrame();inv::Char='t';SENSDict=SensModel
    # Cont=raw["ConL"]
    # TorL='L'
    if TorL=='L' df=raw["Line"]; elseif TorL=='T' df=raw["Trans"] end
    nl::Int64=size(df,1)
    if nl==0 return (0,0,0,0); end
    cls::Int64=size(Cont,1)
    if cls==0
        cl=1;exc=[Vector{Int}() for _ in 1:cl];
        State="Norm"
    else
        cl=cls;exc::Vector{Int64}=Cont."#";
        State="Cont"
    end
    if TorL=='L'
        TR=ones(Float64,nl)
    elseif TorL=='T'
        exc=zeros(cl);TR=df."ratio"
    end
    # gsh=df."g_sh (pu)";bsh=df."b_sh (pu)";g=df."g (pu)";b=df."b (pu)";

    if LinMethod=="None"
        OutExs= Non_Linear_PB(nTP,nl,cl,exc,TR,df,Vars,TorL)
    elseif LinMethod=="Yang"
        OutExs= Linear_Yang_PB(Lininit,nTP,nl,cl,exc,TR,df,Vars,TorL,State)
    elseif LinMethod=="SOCP"
        OutExs=Linear_SOCP_PB(Lininit,nTP,nl,cl,exc,TR,df,Vars,TorL,State)
    elseif LinMethod=="SENS"
        OutExs=Linear_SENS_PB(Lininit,nTP,nl,cl,exc,Vars,TorL,State,SENSDict,raw,Cont)
    end

    return OutExs
end

function Flow_Expr_Simple(LinMethod::String,Lininit::Dict{String,Any},nTP::Int64,raw::Dict{String, Any},Vars::Tuple, TorL::Char; Cont::DataFrame=DataFrame(),inv::Char='t')
    # LinMethod,Lininit,nTP,raw,Vars, TorL=LinMethod,Lininit,nTP,raw,Vars,'L';Cont=DataFrame();inv::Char='t'
    # Cont=raw["ConL"]
    # TorL='L'
    if TorL=='L' df=raw["Line"]; elseif TorL=='T' df=raw["Trans"] end
    nl::Int64=size(df,1)
    cls::Int64=size(Cont,1)
    if cls==0
        cl=1;exc=[Vector{Int}() for _ in 1:cl];
        State="Norm"
    else
        cl=cls;
        exc=[Vector{Int}() for _ in 1:cl]
        if TorL=='L'
            for i in 1:cls
                exc[i]=findall(x-> x==Cont[i,"#"] , df."Line_N")
            end
        end
        State="Cont"
    end
    if nl==0 return (zeros(nTP,nl,cl),zeros(nTP,nl,cl),zeros(nTP,nl,cl),zeros(nTP,nl,cl)); end
    if TorL=='L'
        TR=ones(Float64,nl)
    elseif TorL=='T'
        exc=[Vector{Int}() for _ in 1:cl];TR=df."ratio"
    end
    # gsh=df."g_sh (pu)";bsh=df."b_sh (pu)";g=df."g (pu)";b=df."b (pu)";

    if LinMethod=="None"
        OutExs= Non_Linear_PB(nTP,nl,cl,exc,TR,df,Vars,TorL)
    elseif LinMethod=="Yang"
        OutExs= Linear_Yang_PB(Lininit,nTP,nl,cl,exc,TR,df,Vars,TorL,State)
    elseif LinMethod=="SOCP"
        OutExs=Linear_SOCP_PB(Lininit,nTP,nl,cl,exc,TR,df,Vars,TorL,State)
    elseif LinMethod=="SENS"
        OutExs=Linear_SENS_PB_Simple(Lininit,nTP,nl,cl,exc,Vars,TorL,State,raw,Cont)
    end

    return OutExs
end


function Non_Linear_PB(nTP,nl,cl,exc,TR,df::DataFrame,Vars,TorL,inv='t')
    # g=@view df[:,"g (pu)"];b=@view df[:,"b (pu)"];
    fb=@view df[:,"From"];tb=@view df[:,"To"]
    xe,xf=Vars
    invTR=(1 ./ TR)'[ones(Int8,nTP),:,ones(Int8,cl)]
    invTR2=(invTR.^2)
    gsh_2=(df[:,"g_sh (pu)"]/2)'[ones(Int8,nTP),:,ones(Int8,cl)]
    bsh_2=(df[:,"b_sh (pu)"]/2)'[ones(Int8,nTP),:,ones(Int8,cl)]
    g=(df[:,"g (pu)"])'[ones(Int8,nTP),:,ones(Int8,cl)]
    b=(df[:,"b (pu)"])'[ones(Int8,nTP),:,ones(Int8,cl)]

    temp=Model(); # Create a temporary model to remove the needless expression
    V2=@. (xe^2) + (xf^2)
    V2f=V2[:,fb,:]
    Vcos=@. (xe[:,fb,:] * xe[:,tb,:]) + (xf[:,fb,:] * xf[:,tb,:])
    Vsin=@. (xf[:,fb,:] * xe[:,tb,:]) - (xe[:,fb,:] * xf[:,tb,:])

        sh= @. invTR2*(gsh_2*V2f)
        sr1=@. invTR2*(g*V2f)
        sr2=@. invTR*( -g*Vcos - b*Vsin)
        # @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl], invTR2[n]*gsh_2[n]*V2f[t,n,c])
        # @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl], invTR2[n]*(g[n])*V2f[t,n,c])
        # @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl],  invTR[n]*(
        #                                     -(g[n])*Vcos[t,n,c]
        #                                     -(b[n])*Vsin[t,n,c]))
        @expression(temp, Npij, sh.+sr1.+sr2)
    # unregister.(temp,[:sh,:sr1,:sr2])

        sh=@. invTR2*(-bsh_2*V2f)
        sr1=@. invTR2*(-b*V2f)
        sr2=@. invTR*( b*Vcos - g*Vsin)
        # @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl], invTR2[n]*(-bsh_2[n])*V2f[t,n,c])
        # @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl], invTR2[n]*(-b[n])*V2f[t,n,c])
        # @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl],  invTR[n]*(
        #                                     (b[n])*Vcos[t,n,c]
        #                                     -(g[n])*Vsin[t,n,c]))
        @expression(temp, Nqij, sh.+sr1.+sr2)

        for c in 1:cl
            if exc[c]!=[]
                Npij[:,exc[c],c].=0
                Nqij[:,exc[c],c].=0
            end
        end
    if(inv=='t')
        # unregister.(temp,[:sh,:sr1,:sr2])
        V2t=V2[:,tb,:]

            sh=@. gsh_2*V2t
            sr1=@. g*V2t
            sr2=@. invTR*( -g*Vcos + b*Vsin)
            # @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl], gsh_2[n]*V2t[t,n,c])
            # @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl], (g[n])*V2t[t,n,c])
            # @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl], (invTR[n])*(
            #                                         -(g[n])*Vcos[t,n,c]
            #                                         +(b[n])*Vsin[t,n,c]))
            @expression(temp, Npji, sh.+sr1.+sr2)



        # unregister.(temp,[:sh,:sr1,:sr2])
            sh=@. -bsh_2*V2t
            sr1=@. -b*V2t
            sr2=@. invTR*( b*Vcos + g*Vsin)

            # @expression(temp,sh[t=1:nTP,n=1:nl,c=1:cl], (-bsh_2[n])*V2t[t,n,c])
            # @expression(temp,sr1[t=1:nTP,n=1:nl,c=1:cl], (-b[n])*V2t[t,n,c])
            # @expression(temp,sr2[t=1:nTP,n=1:nl,c=1:cl],  invTR[n]*(
            #                                         (b[n])*Vcos[t,n,c]
            #                                     +(g[n])*Vsin[t,n,c]))
            @expression(temp, Nqji, sh.+sr1.+sr2)

            for c in 1:cl
                if exc[c]!=[]
                    Npji[:,exc[c],c].=0
                    Nqji[:,exc[c],c].=0
                end
            end

    end
    if TorL == 'L'
        gb2= (g.^2+b.^2)
        ΔVline=((xe[:,fb,:] .- xe[:,tb,:]).^2) .+ ((xf[:,fb,:] .- xf[:,tb,:]).^2)
        NI2=@. gb2*ΔVline
        @expression(temp,NIflow2, NI2)# lineflow current limitation
        for c in 1:cl
            if exc[c]!=[]
                NIflow2[:,exc[c],c].=0
            end
        end
    end
    if TorL == 'L'
        return (Npij, Nqij, Npji, Nqji, NIflow2)
    elseif TorL == 'T'
        return (Npij, Nqij, Npji, Nqji,0)
    end
end
function compute_warmstart_line_params(Δθ0, V0f,V0t,g, b)

    s1 = cos.(Δθ0)
    s0 =@. (sin(Δθ0) - Δθ0 * cos(Δθ0))
    c1 = -sin.(Δθ0)
    c0 =@. (cos(Δθ0) + Δθ0 * sin(Δθ0))

    g_p_2 = @. (g * c0 + b * s0 + (g * c1 + b * s1) * Δθ0)/2
    b_p = @. ((g * c1 + b * s1) * V0f * V0t)
    b_q_2 = @. ((-g * s0 + b * c0) - (g * s1 - b * c1) * Δθ0)/2
    g_q = @. ((g * s1 - b * c1) * V0f * V0t)

    return (g_p_2,b_p,g_q,b_q_2)
end
function Linear_Yang_PB(Lininit::Dict{String,Any},nTP,nl,cl,exc,TR,df,Vars,TorL,State,inv='t')
    if State=="Norm"
        V_0=Lininit["V"];θ_0=Lininit["θ"];
    elseif State=="Cont"
        V_0=Lininit["V_c_Line"];θ_0=Lininit["θ_c_Line"];
    end
    warmstartB=Lininit["warmstart_yang"][1]

    fb=@view df[:,"From"];tb=@view df[:,"To"]
    xV2,xθ=Vars
    invTR=(1 ./ TR)'[ones(Int8,nTP),:,ones(Int8,cl)]
    invTR2=(invTR.^2)
    gsh_2=(df[:,"g_sh (pu)"]/2)'[ones(Int8,nTP),:,ones(Int8,cl)]
    bsh_2=(df[:,"b_sh (pu)"]/2)'[ones(Int8,nTP),:,ones(Int8,cl)]
    g=(df[:,"g (pu)"])'[ones(Int8,nTP),:,ones(Int8,cl)]
    b=(df[:,"b (pu)"])'[ones(Int8,nTP),:,ones(Int8,cl)]
    b_2=@. b/2
    g_2=@. g/2

    temp=Model(); # Create a temporary model to remove the needless expression
    V2f=xV2[:,fb,:]
    V2t=xV2[:,tb,:]
    V0f=V_0[:,fb,:]
    V0t=V_0[:,tb,:]
    ΔV0=@. V0f-V0t
    ΔV2=@. V2f-V2t
    Vsum0=@. V0f+V0t
    Vsum=@. V2f+V2t
    Δθ= @. xθ[:,fb,:] - xθ[:,tb,:]
    Δθ0= @. θ_0[:,fb,:] - θ_0[:,tb,:]

    if warmstartB==1 g_p_2,b_p,g_q,b_q_2=compute_warmstart_line_params(Δθ0, V0f,V0t,g, b) end

    sh= @. invTR2*(gsh_2*V2f)
    sr1=@. invTR2*(g*V2f)
    if warmstartB==1
        sr2=@. invTR*( -g_p_2*Vsum - b_p*Δθ)
    else
        sr2=@. invTR*( -g_2*Vsum - b*Δθ)
    end
    lossNpij_V=@. invTR*(g*(ΔV0/Vsum0)*ΔV2 -g_2*(ΔV0^2))
    loss_θ=@. invTR*( g*Δθ0*Δθ - g_2*(Δθ0^2))
    # lossij=@. lossNpij_V + loss_θ
    @expression(temp, Npij, sh .+ sr1 .+ sr2 .+ lossNpij_V .+ loss_θ)

    sh=@. invTR2*(-bsh_2*V2f)
    sr1=@. invTR2*(-b*V2f)
    if warmstartB == 1
        sr2=@. invTR*( b_q_2*Vsum - g_q*Δθ)
    else
        sr2=@. invTR*( b_2*Vsum - g*Δθ)
    end
    loss_V=@. invTR*( -b*(ΔV0/Vsum0)*ΔV2 +b_2*(ΔV0^2))
    loss_θ=@. invTR*( -b*Δθ0*Δθ + b_2*(Δθ0^2))
    # lossij=@. loss_V + loss_θ
    @expression(temp, Nqij, sh .+ sr1 .+ sr2 .+ loss_V .+ loss_θ)

    for c in 1:cl
        if exc[c]!=[]

            Npij[:,exc[c],c].=0
            Nqij[:,exc[c],c].=0
        end
    end

    if(inv=='t')

        sh= @. gsh_2*V2t
        sr1=@. g*V2t
        if warmstartB == 1
            sr2=@. invTR*( -g_p_2*Vsum + b_p*Δθ)
        else
            sr2=@. invTR*( -g_2*Vsum + b*Δθ)
        end
        lossNpji_V=@. invTR*(g*(ΔV0/Vsum0)*ΔV2 -g_2*(ΔV0^2))
        loss_θ=@. invTR*( g*Δθ0*Δθ - g_2*(Δθ0^2))
        # lossij=@. lossNpji_V + loss_θ
        @expression(temp, Npji, sh .+ sr1 .+ sr2 .+ lossNpji_V .+ loss_θ)

        sh=@. -bsh_2*V2t
        sr1=@. -b*V2t
        if warmstartB == 1
            sr2=@. invTR*( b_q_2*Vsum + g_q*Δθ)
        else
            sr2=@. invTR*( b_2*Vsum  + g*Δθ)
        end
        loss_V=@. invTR*(-b*(ΔV0/Vsum0)*ΔV2 + b_2*(ΔV0^2))
        loss_θ=@. invTR*( -b*Δθ0*Δθ + b_2*(Δθ0^2))
        # lossji=@. loss_V + loss_θ
        @expression(temp, Nqji, sh .+ sr1 .+ sr2 .+ loss_V .+ loss_θ)

        for c in 1:cl
            if exc[c]!=[]
                Npji[:,exc[c],c].=0
                Nqji[:,exc[c],c].=0
            end
        end

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
        @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : gsh_2[n]*(xV2[t,fb[n],c])) #Line flow expression
        @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (g[n])*(xV2[t,fb[n],c])) #Line flow expression
        @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 :  (TR[n])*(
                                    -(g[n])*(xA[t,n,c])
                                    -(b[n])*(xB[t,n,c]))) #Line flow expression

        @expression(temp, Npij, sh.+sr1.+sr2)
    temp=Model()
        @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-bsh_2[n])*(xV2[t,fb[n],c])) #Line flow expression
        @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (-b[n])*(xV2[t,fb[n],c])) #Line flow expression

        @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                                                            +(b[n])*(xA[t,n,c])
                                                                            -(g[n])*(xB[t,n,c]))) #Line flow expression
        @expression(temp, Nqij, sh .+ sr1 .+ sr2)

    if(inv=='t')
        temp=Model()
            @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*gsh_2[n]*(xV2[t,tb[n],c])) #Line flow expression
            @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(g[n])*(xV2[t,tb[n],c])) #Line flow expression

            @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                                                            -(g[n])*(xA[t,n,c])
                                                                            +(b[n])*(xB[t,n,c])))#Line flow expression
            @expression(temp, Npji, sh .+ sr1 .+ sr2)


        temp=Model()
            @expression(temp, sh[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(-bsh_2[n])*(xV2[t,tb[n],c])) #Line flow expression
            @expression(temp, sr1[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n]^2)*(-b[n])*(xV2[t,tb[n],c])) #Line flow expression

            @expression(temp, sr2[t=1:nTP,n=1:nl,c=1:cl],exc[c]==n ? 0 : (TR[n])*(
                                                                            +(g[n])*(xA[t,n,c])
                                                                            +(b[n])*(xB[t,n,c])))#Line flow expression

            @expression(temp, Nqji, sh .+ sr1 .+ sr2)

    end

    if TorL == 'L'
        @expression(temp,NIflow2[t=1:nTP,n=1:nl,c=1:cl], exc[c]==n ? 0 : (g[n]^2 + b[n]^2)
                                            *(xV2[t,fb[n],c]+xV2[t,tb[n],c]-2*xA[t,n,c]))# lineflow current limitation
    end

    if TorL == 'L'
        return (Npij, Nqij, Npji, Nqji, NIflow2)
    elseif TorL == 'T'
        return (Npij, Nqij, Npji, Nqji)
    end

end
function Linear_SENS_PB(Lininit::Dict{String,Any},nTP::Int64,nl,cl,exc,Vars,TorL,State,SENSDict::Dict{String, Array{Float64, 4}},raw::Dict{String, Any},Cont,inv='t')
    xΔPinj,xΔQinj=Vars
    if State=="Norm"
        V_0,θ_0=Lininit["V"],Lininit["θ"];
    elseif State=="Cont"
        V_0,θ_0=Lininit["V_c_Line"],Lininit["θ_c_Line"];
    end
    initialVar=(V_0.*cos.(θ_0), V_0.*sin.(θ_0));
    pij, qij, pji, qji,I2 =Flow_Expr("None",Dict{String,Any}(),nTP,raw,initialVar,TorL;Cont=Cont);
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

function Linear_SENS_PB_Simple(Lininit::Dict{String,Any},nTP::Int64,nl,cl,exc,Vars,TorL,State,raw::Dict{String, Any},Cont,inv='t')
    xΔPinjP,xΔPinjN,xΔQinjP,xΔQinjN=Vars
    xΔPinj=xΔPinjP-xΔPinjN
    xΔQinj=xΔQinjP-xΔQinjN
    if State=="Norm"
        V_0,θ_0=  reshape(Lininit["V"],size(Lininit["V"])...,1), reshape(Lininit["θ"],size(Lininit["θ"])...,1);
    elseif State=="Cont"
        V_0,θ_0=Lininit["V_c_Line"],Lininit["θ_c_Line"];
    end
    initialVar=(V_0.*cos.(θ_0), V_0.*sin.(θ_0));
    pij, qij, pji, qji,I2 =Flow_Expr_Simple("None",Dict{String,Any}(),nTP,raw,initialVar,TorL;Cont=Cont);

    if TorL=='L' n_1=0 elseif TorL=='T' n_1=raw["nLine"]  end

    Npij=zeros(AffExpr,nTP,nl,cl)
    Nqij=zeros(AffExpr,nTP,nl,cl)
    Npij.=pij
    Nqij.=qij
    if TorL == 'L'
        # Nploss=[AffExpr(0.0) for t=1:nTP,n=1:nl,c=1:cl]
        NIflow2=zeros(AffExpr,nTP,nl,cl)
        # NIflow2.+=I2
    end
    if (inv=='t')
        Npji=zeros(AffExpr,nTP,nl,cl)
        Nqji=zeros(AffExpr,nTP,nl,cl)
        Npji.=pji
        Nqji.=qji
    end
    for c=1:cl
        Npij[:,exc[c].- n_1,c].=0
        Nqij[:,exc[c].- n_1,c].=0
        Npji[:,exc[c].- n_1,c].=0
        Nqji[:,exc[c].- n_1,c].=0
        # if TorL == 'L'
        #     for t=1:nTP
        #         add_to_expression!.(NIflow2[t,:,c],SENSDict["S_IP"][t,c][(1:nl).+n_1,:]*xΔPinj[t,:,c])
        #         add_to_expression!.(NIflow2[t,:,c],SENSDict["S_IQ"][t,c][(1:nl).+n_1,:]*xΔQinj[t,:,c])
        #     end
        #     NIflow2[:,exc[c].- n_1,c].=0
        # end
    end
    if TorL == 'L'
        Parlist= NIflow2,I2,V_0,θ_0
        return (Npij, Nqij, Npji, Nqji, Parlist)
    elseif TorL == 'T'
        return (Npij, Nqij, Npji, Nqji)
    end
end


function Non_Linear_CP(SCOPF,Vars,Parlist,nTP,raw,States,relaxconfac=0.0)

    xe,xf=Vars
    NIflow2=Parlist
    F2maxtemp=@. raw["Line"]."RATE_A"^2
    slackbus=raw["Bus"][raw["Bus"]."Bus_Type".==3,"Bus_N"]
    V2mag=@. xe^2 + xf^2
    θmag=@. atan(xf,xe)
    if States=="Norm"
        cl=1
        V2maxtemp=@. raw["Bus"]."V_max (pu)"^2
        V2mintemp=@. raw["Bus"]."V_min (pu)"^2
        @expression(SCOPF,Line_flow,NIflow2)
        @constraint(SCOPF,NLine_flow2[t=1:nTP, n=1:raw["nLine"],c=1:cl],NIflow2[t,n,c]<= F2maxtemp[n])
        @constraint(SCOPF,NV2mag[t=1:nTP, n=1:raw["nBus"],c=1:cl], V2mintemp[n]<= V2mag[t,n,c]<=V2maxtemp[n])
        @expression(SCOPF,NV2magexp, V2mag)
        @expression(SCOPF,Nθexp, θmag)
        @constraint(SCOPF,NslackVangle[t=1:nTP, n=slackbus,c=1:cl], xf[t,n,c]==0)
    elseif States=="Cont"
        cl=raw["nConL"]
        V2maxtemp=@. (raw["Bus"]."V_max (pu)"*(1+relaxconfac))^2
        V2mintemp=@. (raw["Bus"]."V_min (pu)"*(1-relaxconfac))^2
        @expression(SCOPF,Line_flow_c_Line,NIflow2)
        @constraint(SCOPF,CLine_flow2[t=1:nTP, n=1:raw["nLine"],c=1:cl],NIflow2[t,n,c]<= F2maxtemp[n]*(1+relaxconfac))
        @constraint(SCOPF,CV2mag[t=1:nTP, n=1:raw["nBus"],c=1:cl], V2mintemp[n]<= V2mag[t,n,c]<=V2maxtemp[n])
        @expression(SCOPF,CV2magexp, V2mag)
        @expression(SCOPF,Cθexp, θmag)
        @constraint(SCOPF,CslackVangle[t=1:nTP, n=1:raw["nBus"],c=1:cl;raw["Bus"]."Bus_Type"[n]==3], xf[t,n,c]==0)
    end
    return SCOPF
end
function Linear_Yang_CP(SCOPF,Vars,Parlist,nTP,raw,States,relaxconfac=0.0)
    xV2,xθ=Vars
    Npij, Nqij,lossNpij_V,lossNpji_V=Parlist
    slackbus=raw["Bus"][raw["Bus"]."Bus_Type".==3,"Bus_N"]
    Fmaxtemp=raw["Line"]."RATE_A"
    if States=="Norm"
        cl=1
        pi_6=ones(nTP,raw["nLine"],cl)*tan(π/6)
        FmaxtempMAT=(Fmaxtemp)'[ones(Int8,nTP),:,ones(Int8,cl)]
        V2maxtemp=@. raw["Bus"]."V_max (pu)"^2
        V2mintemp=@. raw["Bus"]."V_min (pu)"^2
        @constraint(SCOPF,NslackVangle[t=1:nTP,c=1:cl], xθ[t,slackbus,c].==0)
        @constraint(SCOPF,NV2mag1[t=1:nTP, n=1:raw["nBus"],c=1:cl], V2mintemp[n] <=xV2[t,n,c])
        @constraint(SCOPF,NV2mag2[t=1:nTP, n=1:raw["nBus"],c=1:cl], xV2[t,n,c]<=V2maxtemp[n] )
        @expression(SCOPF,NV2magexp[t=1:nTP, n=1:raw["nBus"],c=1:cl], xV2[t,n,c])
        @expression(SCOPF,Nθexp[t=1:nTP, n=1:raw["nBus"],c=1:cl], xθ[t,n,c])
        @constraint(SCOPF, Nflow1, Npij .+ (pi_6.*Nqij) .<= FmaxtempMAT)
        @constraint(SCOPF, Nflow2, Npij .- (pi_6.*Nqij) .<= FmaxtempMAT)
        @constraint(SCOPF, Nflow3, Npij .+ (pi_6.*Nqij) .>= -FmaxtempMAT)
        @constraint(SCOPF, Nflow4, Npij .- (pi_6.*Nqij) .>= -FmaxtempMAT)
        if relaxpermit==false
            @constraint(SCOPF,Nloss1, lossNpij_V .>= 0)
            @constraint(SCOPF,Nloss2, lossNpji_V .>= 0)
        end
        Lflow=zeros(AffExpr, nTP, raw["nLine"], cl)
        Lflow=@. Npij^2+Nqij^2
        @expression(SCOPF,Line_flow,Lflow)
    elseif States=="Cont"
        cl=raw["nConL"]
        pi_6=ones(nTP,raw["nLine"],cl)*tan(π/6)
        FmaxtempMAT=(Fmaxtemp)'[ones(Int8,nTP),:,ones(Int8,cl)].*(1+relaxconfac)
        V2maxtemp=@. (raw["Bus"]."V_max (pu)"*(1+relaxconfac))^2
        V2mintemp=@. (raw["Bus"]."V_min (pu)"*(1-relaxconfac))^2
        @constraint(SCOPF,CslackVangle[t=1:nTP,c=1:cl], xθ[t,slackbus,c].==0)
        @constraint(SCOPF,CV2mag1[t=1:nTP, n=1:raw["nBus"],c=1:cl], V2mintemp[n] <=xV2[t,n,c])
        @constraint(SCOPF,CV2mag2[t=1:nTP, n=1:raw["nBus"],c=1:cl], xV2[t,n,c]<=V2maxtemp[n] )
        @expression(SCOPF,CV2magexp[t=1:nTP, n=1:raw["nBus"],c=1:cl], xV2[t,n,c])
        @expression(SCOPF,Cθexp[t=1:nTP, n=1:raw["nBus"],c=1:cl], xθ[t,n,c])
        @constraint(SCOPF, Cflow1, Npij .+ (pi_6.*Nqij) .<= FmaxtempMAT)
        @constraint(SCOPF, Cflow2, Npij .- (pi_6.*Nqij) .<= FmaxtempMAT)
        @constraint(SCOPF, Cflow3, Npij .+ (pi_6.*Nqij) .>= -FmaxtempMAT)
        @constraint(SCOPF, Cflow4, Npij .- (pi_6.*Nqij) .>= -FmaxtempMAT)
        if relaxpermit==false
            @constraint(SCOPF,Closs1, lossNpij_V.>=0)
            @constraint(SCOPF,Closs2, lossNpji_V.>=0)
        end
        Lflow=zeros(AffExpr, nTP, raw["nLine"], cl)
        Lflow=@. Npij^2+Nqij^2
        @expression(SCOPF,Line_flow_c_Line,Lflow)
    end
    return SCOPF
end
function Linear_SOCP_CP(SCOPF,Vars,Parlist,nTP,raw,State,relaxconfac=0.0)
    xV2,xA,xB=Vars
    NIflow2=Parlist
    if States=="Norm"
        cl=1
        @expression(SCOPF,Line_flow[t=1:nTP, n=1:raw["nLine"],c=1:cl],NIflow2[t,n,c])
        @constraint(SCOPF,NLine_flow2[t=1:nTP, n=1:raw["nLine"],c=1:cl],NIflow2[t,n,c]<= F2maxtemp[n])
        @constraint(SCOPF,NV2mag[t=1:nTP, n=1:raw["nBus"],c=1:cl], V2mintemp[n]<= xV2[t,n,c]<=V2maxtemp[n])
        @expression(SCOPF,NV2magexp[t=1:nTP, n=1:raw["nBus"],c=1:cl], xV2[t,n,c]^2)
        @expression(SCOPF,Ndiffθexp[t=1:nTP, n=1:raw["nLine"],c=1:cl], atan(xB[t,n,c],xA[t,n,c]))
        # @constraint(SCOPF,NslackVangle[t=1:nTP, n=1:raw["nBus"],c=1:cl;raw["Bus"]."Bus_Type"[n]==3], xf[t,n,c]==0)
    elseif States=="Cont"
        cl=raw["nConL"]
        @expression(SCOPF,Line_flow_c_Line[t=1:nTP, n=1:raw["nLine"],c=1:cl],NIflow2[t,n,c])
        @constraint(SCOPF,CLine_flow2[t=1:nTP, n=1:raw["nLine"],c=1:cl],NIflow2[t,n,c]<= F2maxtemp[n]*(1+relaxconfac))
        @constraint(SCOPF,CV2mag[t=1:nTP, n=1:raw["nBus"],c=1:cl], (raw["Bus"]."V_min (pu)"[n]*(1-relaxconfac))^2<= xV2[t,n,c]<=(raw["Bus"]."V_max (pu)"[n]*(1+relaxconfac))^2)
        @expression(SCOPF,CV2magexp[t=1:nTP, n=1:raw["nBus"],c=1:cl], xV2[t,n,c])
        @expression(SCOPF,Cdiffθexp[t=1:nTP, n=1:raw["nLine"],c=1:cl], atan(xB[t,n,c],xA[t,n,c]))
        # @constraint(SCOPF,CslackVangle[t=1:nTP, n=1:raw["nBus"],c=1:cl;raw["Bus"]."Bus_Type"[n]==3], xf[t,n,c]==0)
    end
    return SCOPF
end
function Linear_SENS_CP(SCOPF,Vars,Parlist,nTP,raw::Dict{String, Any},States,SENSDict::Dict{String,Any},Lininit::Dict{String,Any})

    xΔPinjP,xΔPinjN,xΔQinjP,xΔQinjN=Vars
    xΔPinj=xΔPinjP-xΔPinjN
    xΔQinj=xΔQinjP-xΔQinjN
    NIflow2 , I2_0 , V_0 , θ_0=Parlist

    F2maxtemp=@. raw["Line"]."RATE_A"^2
    slackbus=raw["Bus"][raw["Bus"]."Bus_Type".==3,"Bus_N"][1]

    if States=="Norm"
        cl=1
    # Technical constraints
    #---------------------------------------



        block=zeros(AffExpr,raw["nLine"])
        ΔI2=fill(block,nTP,cl);
        I2=zeros(AffExpr,nTP,raw["nLine"],cl);
        I2.=I2_0;
        for t=1:nTP,c=1:cl
            ΔI2[t,c]=SENSDict["S_IP"][t,c][1:raw["nLine"],:]*xΔPinj[t,:,c] + SENSDict["S_IQ"][t,c][1:raw["nLine"],:]*xΔQinj[t,:,c]
            I2[t,:,c]=I2_0[t,:,c]+ΔI2[t,c]
        end
        @expression(SCOPF,Line_flow,I2)
        # map_coefficients_inplace!.(x->round(x,digits=2),Line_flow)
        # for n=1:raw["nLine"]
        #     minthreshold=0.01*(F2maxtemp[n])/(TRfree+1)
        #     map_coefficients_inplace!.(x->ifelse.(abs(x) .> minthreshold, x, 0.0),Line_flow[:,n,:])
        # end
        # drop_zeros!.(Line_flow)

        if relaxpermit==true
            @variable(SCOPF, slk_Imax[t=1:nTP,n=1:raw["nLine"], c=1:cl]>=0)
            @variable(SCOPF, slk_Vmax[t=1:nTP,n=1:raw["nBus"], c=1:cl]>=0)
            @variable(SCOPF, slk_Vmin[t=1:nTP,n=1:raw["nBus"], c=1:cl]>=0)

            # @constraint(SCOPF,slk_Vminlimit, slk_Vmin.<=0.0507164)
            # @constraint(SCOPF,slk_Vmaxlimit, slk_Vmax.<=0.0507164)
            # @constraint(SCOPF,slk_Imaxlimit, slk_Imax.<=0.0507164)
        end
        # map_coefficients_inplace!.(x->ifelse.(maxthreshold .> x .> minthreshold, x, 0.0),Line_flow)
        # drop_zeros!.(Line_flow)
        Imaxlimit0=(F2maxtemp')[ones(Int8,nTP),:,ones(Int8,cl)]
        #        sum(abs.(SENSDict["S_VQ"][1][5,:]))+sum(abs.(SENSDict["S_VP"][1][5,:]))
        ovr_Sens=zeros(nTP,raw["nLine"],cl)
        ovr_Sens_P=zeros(nTP,raw["nLine"],cl)
        ovr_Sens_Q=zeros(nTP,raw["nLine"],cl)
        for i=1:raw["nLine"],t=1:nTP,c=1:cl
            ovr_Sens[t,i,c]=sum(abs.(SENSDict["S_IQ"][t,c][i,:]))*Lininit["TR"][2] +sum(abs.(SENSDict["S_IP"][t,c][i,:]))*Lininit["TR"][1]
            ovr_Sens_P[t,i,c]=sum(abs.(SENSDict["S_IP"][t,c][i,:]))
            ovr_Sens_Q[t,i,c]=sum(abs.(SENSDict["S_IQ"][t,c][i,:]))
        end

        # abs.(SENSDict["S_IP"][1,1,1]).>0.01


        # Lininit["TR"][1]=minimum([(1.5*round(maximum(max.((I2_0.-Imaxlimit0),0)./ovr_Sens_P)*20,digits=2)+0.01),2])

        Imaxindex1=findall(x->x==1,Imaxlimit0.-I2_0.<=1*ovr_Sens .&& ovr_Sens.>0.001*Lininit["TR"][1])
        Imaxindex2=findall(x->x==1,(I2_0./Imaxlimit0).>=1*ones(nTP,raw["nLine"],cl))
        Imaxindex=union([Imaxindex1;Imaxindex2])

        # setdiff!(Imaxindex,findall(x->x==0,NIflow2-I2_0))
        if relaxpermit==true
            @expression(SCOPF,Imaxlimit, Imaxlimit0.*(1 .+ slk_Imax))
        else
            @expression(SCOPF,Imaxlimit, Imaxlimit0)
        end
        linear_Imaxindex = [LinearIndices(Line_flow)[ind] for ind in Imaxindex]
        @constraint(SCOPF,NLine_flow2[ind in linear_Imaxindex], Line_flow[ind]<= Imaxlimit[ind])
        block=zeros(AffExpr,raw["nBus"])
        ΔV=fill(block,nTP,cl);Δθ=fill(block,nTP,cl);
        V=zeros(AffExpr,nTP,raw["nBus"],cl);θ=zeros(AffExpr,nTP,raw["nBus"],cl)
        V.=V_0;θ.=θ_0
        for t=1:nTP,c=1:cl
            ΔV[t,c]=SENSDict["S_VP"][t,c]*xΔPinj[t,:,c]+ SENSDict["S_VQ"][t,c]*xΔQinj[t,:,c]
            V[t,:,c]=V_0[t,:,c]+ΔV[t,c]
            Δθ[t,c]=SENSDict["S_θP"][t,c]*xΔPinj[t,:,c]+ SENSDict["S_θQ"][t,c]*xΔQinj[t,:,c]
            θ[t,:,c]=θ_0[t,:,c]+Δθ[t,c]
        end
        @expression(SCOPF,NV2magexp, V)
        # minthreshold=0.001*raw["Bus"]."V_min (pu)"[1]/(0.1*(TRfree+1))
        # map_coefficients_inplace!.(x->ifelse.(abs(x) .> minthreshold, x, 0.0),NV2magexp)
        # drop_zeros!.(NV2magexp)
        # map_coefficients_inplace!.(x->round(x,digits=4),NV2magexp)

        ovr_Sens_V=zeros(nTP,raw["nBus"],cl)
        ovr_Sens_V_P=zeros(nTP,raw["nBus"],cl)
        ovr_Sens_V_Q=zeros(nTP,raw["nBus"],cl)
        ovr_Sens_V10=zeros(nTP,raw["nBus"],cl)
        for i=1:raw["nBus"],t=1:nTP,c=1:cl
            ovr_Sens_V[t,i,c]=sum(abs.(SENSDict["S_VQ"][t,c][i,:])) *Lininit["TR"][2]+sum(abs.(SENSDict["S_VP"][t,c][i,:]))*Lininit["TR"][1]
            ovr_Sens_V_P[t,i,c]=sum(abs.(SENSDict["S_VP"][t,c][i,:]))
            ovr_Sens_V_Q[t,i,c]=sum(abs.(SENSDict["S_VQ"][t,c][i,:]))
            sorted=sort(nonzeros([abs.(SENSDict["S_VQ"][t,c][i,:])*Lininit["TR"][2];abs.(SENSDict["S_VP"][t,c][i,:])]*Lininit["TR"][1]))
            ovr_Sens_V10[t,i,c]=sum(sorted[end-min(length(sorted)-1,10):end])
        end
        # ovr_Sens_V=ovr_Sens_V_P.*Lininit["TR"][1]+ovr_Sens_V_Q.*Lininit["TR"][2]


        Vminlimit0=(raw["Bus"]."V_min (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
        Vmaxlimit0=(raw["Bus"]."V_max (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]


        #   # round(maximum(max.((I2_0.-Imaxlimit0),0)./ovr_Sens_P)*20,digits=2)+0.01
        #   findmax(max.((V_0.-Vmaxlimit0),0))
        #   ovr_Sens_V_Q[20,37,1]

        #   findmax(max.((V_0.-Vmaxlimit0),0)./ovr_Sens_V_Q)

        #   (V_0.-Vmaxlimit0)[20,27,1]
        #   # median(filter(x-> x>0,(max.((V_0.-Vmaxlimit0),0))))


        Vmaxindex1= findall(x->x==1,Vmaxlimit0.-V_0 .<=1*ovr_Sens_V10 .&& ovr_Sens_V.>0.01*Lininit["TR"][2])
        Vmaxindex2=findall(x->x==1,V_0.>=(raw["Bus"]."V_max (pu)".-0.00)'[ones(Int8,nTP),:,ones(Int8,cl)])
        Vmaxindex=union([Vmaxindex1;Vmaxindex2])
        Vminindex1= findall(x->x==1,V_0.- Vminlimit0.<=1*ovr_Sens_V10 .&& ovr_Sens_V.>0.01*Lininit["TR"][2])
        Vminindex2=findall(x->x==1,V_0.<=(raw["Bus"]."V_min (pu)".+0.00)'[ones(Int8,nTP),:,ones(Int8,cl)])
        Vminindex=union([Vminindex1;Vminindex2])
        # setdiff!(Vminindex,findall(x->x==0,V-V_0))
        # Vmaxindex=findall(x->x==1,V_0.>=(raw["Bus"]."V_max (pu)".-0.04)'[ones(Int8,nTP),:,ones(Int8,cl)])
        # setdiff!(Vmaxindex,findall(x->x==0,V-V_0))

        if relaxpermit==true
            @expression(SCOPF,Vminlimit, Vminlimit0.*(1 .- slk_Vmin))
            @expression(SCOPF,Vmaxlimit, Vmaxlimit0.*(1 .+ slk_Vmax))
        else
            @expression(SCOPF,Vminlimit, Vminlimit0)
            @expression(SCOPF,Vmaxlimit, Vmaxlimit0)
        end
        linear_Vminindex = [LinearIndices(NV2magexp)[ind] for ind in Vminindex]
        linear_Vmaxindex = [LinearIndices(NV2magexp)[ind] for ind in Vmaxindex]

        @constraint(SCOPF,NV2mag1[ind in linear_Vminindex], Vminlimit[ind]<= NV2magexp[ind])
        @constraint(SCOPF,NV2mag2[ind in linear_Vmaxindex], NV2magexp[ind]<=Vmaxlimit[ind])
        @expression(SCOPF,Nθexp, θ)
        @constraint(SCOPF,NslackVangle, Nθexp[:,slackbus,:].==0)
    # Trust region and modelling constraints
    #---------------------------------------
    # if relaxpermit==true
    #     @variable(SCOPF, slk_TR_P[t=1:nTP,n=1:raw["nPActBus"], c=1:cl]>=0)
    #     @variable(SCOPF, slk_TR_Q[t=1:nTP,n=1:raw["nPActBus"], c=1:cl]>=0)
    #     @variable(SCOPF, slk_SB_P[t=1:nTP,n=1:raw["nQActBus"], c=1:cl]>=0)
    #     @variable(SCOPF, slk_SB_Q[t=1:nTP,n=1:raw["nQActBus"], c=1:cl]>=0)
    # end
    if false #relaxpermit==true
        @variable(SCOPF, N_TRfreeP[t=1:nTP,i=1:raw["nPActBus"],c=1:cl]>=0,start=0)
        @variable(SCOPF, N_TRfreeQ[t=1:nTP,n=1:raw["nQActBus"],c=1:cl]>=0,start=0)
        @variable(SCOPF, N_SBfreeP[t=1:nTP,c=1:cl]>=0,start=0)
        @variable(SCOPF, N_SBfreeQ[t=1:nTP,c=1:cl]>=0,start=0)
    else
        @expression(SCOPF, N_TRfreeP, zeros(nTP,raw["nPActBus"],cl))
        @expression(SCOPF, N_TRfreeQ, zeros(nTP,raw["nQActBus"],cl))
        @expression(SCOPF, N_SBfreeP, zeros(nTP,cl))
        @expression(SCOPF, N_SBfreeQ, zeros(nTP,cl))
    end

        @constraint(SCOPF,NTrust_RegionP1,xΔPinj.<=Lininit["TR"][1].*(N_TRfreeP.+1))
        @constraint(SCOPF,NTrust_RegionP2,-Lininit["TR"][1].*(N_TRfreeP.+1).<=xΔPinj)
        @constraint(SCOPF,NTrust_RegionQ1,xΔQinj.<=Lininit["TR"][2].*(N_TRfreeQ.+1))
        @constraint(SCOPF,NTrust_RegionQ2,-Lininit["TR"][2].*(N_TRfreeQ.+1).<=xΔQinj)

        # @constraint(SCOPF,NTrust_RegionP1,xΔPinj.<=InitialTR .* (slk_TR_P .+ 1))
        # @constraint(SCOPF,NTrust_RegionP2,-InitialTR .* (slk_TR_P .+ 1).<=xΔPinj)
        # @constraint(SCOPF,NTrust_RegionQ1,xΔQinj.<=InitialTR .* (slk_TR_Q .+ 1))
        # @constraint(SCOPF,NTrust_RegionQ2,-InitialTR .* (slk_TR_Q .+ 1).<=xΔQinj)

        @constraint(SCOPF,Sensitivity_balance_P1[t=1:nTP,c=1:cl],sum(xΔPinj[t,:,c])<=Lininit["TRB"][1]*(N_SBfreeP[t,c]+1))
        @constraint(SCOPF,Sensitivity_balance_Q1[t=1:nTP,c=1:cl],sum(xΔQinj[t,:,c])<=Lininit["TRB"][1]*(N_SBfreeQ[t,c]+1))
        @constraint(SCOPF,Sensitivity_balance_P2[t=1:nTP,c=1:cl],-Lininit["TRB"][1]*(N_SBfreeP[t,c]+1)<=sum(xΔPinj[t,:,c]))
        @constraint(SCOPF,Sensitivity_balance_Q2[t=1:nTP,c=1:cl],-Lininit["TRB"][1]*(N_SBfreeQ[t,c]+1)<=sum(xΔQinj[t,:,c]))

        # @constraint(SCOPF, TRLimit1[t=1:nTP,c=1:cl], sum(xΔPinj[t,:,c].^2)<=36)
        # @constraint(SCOPF, TRLimit2[t=1:nTP,c=1:cl], sum(xΔQinj[t,:,c].^2)<=36)
        # fix.(xΔPinj,0;force=true)
        # fix.(xΔQinj,0;force=true)

    elseif States=="Cont"
        cl=raw["nConL"]
    # Technical constraints
    #---------------------------------------
        block=zeros(AffExpr,raw["nLine"])
        ΔI2=fill(block,nTP,cl);
        I2=zeros(AffExpr,nTP,raw["nLine"],cl);
        I2.=I2_0;
        for t=1:nTP,c=1:cl
            ΔI2[t,c]=SENSDict["S_IP"][t,c][1:raw["nLine"],:]*xΔPinj[t,:,c] + SENSDict["S_IQ"][t,c][1:raw["nLine"],:]*xΔQinj[t,:,c]
            I2[t,:,c]=I2_0[t,:,c]+ΔI2[t,c]
        end
        @expression(SCOPF,Line_flow_c_Line,I2)
        # map_coefficients_inplace!.(x->round(x,digits=4),Line_flow_c_Line)
        # drop_zeros!.(Line_flow_c_Line)
        if relaxpermit==true
            @variable(SCOPF, slk_Imax_c[t=1:nTP,n=1:raw["nLine"],c=1:cl]>=0)
            @variable(SCOPF, slk_Vmax_c[t=1:nTP,n=1:raw["nBus"],c=1:cl]>=0)
            @variable(SCOPF, slk_Vmin_c[t=1:nTP,n=1:raw["nBus"],c=1:cl]>=0)


            # @constraint(SCOPF,slk_Vmin_climit, slk_Vmin_c.<=0.185724)
            # @constraint(SCOPF,slk_Vmax_climit, slk_Vmax_c.<=0.185724)
            # @constraint(SCOPF,slk_Imax_climit, slk_Imax_c.<=0.185724)
        end
        # fix.(slk_Vmin_c,0; force=true)
        Imaxlimit0=(F2maxtemp')[ones(Int8,nTP),:,ones(Int8,cl)]
        ovr_Sens=zeros(nTP,raw["nLine"],cl)
        for i=1:raw["nLine"],t=1:nTP,c=1:cl
            ovr_Sens[t,i,c]=sum(abs.(SENSDict["S_IQ"][t,c][i,:]))*Lininit["TR"][2] +sum(abs.(SENSDict["S_IP"][t,c][i,:]))*Lininit["TR"][1]
        end

        Imaxindex1=findall(x->x==1,Imaxlimit0.-I2_0.<=1*ovr_Sens .&& ovr_Sens.>0.001*Lininit["TR"][1])
        Imaxindex2=findall(x->x==1,(I2_0./Imaxlimit0).>=1*ones(nTP,raw["nLine"],cl))
        Imaxindex=union([Imaxindex1;Imaxindex2])

        # Imaxlimit0=(F2maxtemp')[ones(Int8,nTP),:,ones(Int8,cl)]
        # Imaxindex=findall(x->x==1,(I2_0./Imaxlimit0).>=0.25*ones(nTP,raw["nLine"],cl))
        # setdiff!(Imaxindex,findall(x->x==0,NIflow2-I2_0))
        if relaxpermit==true
            @expression(SCOPF,CImaxlimit, Imaxlimit0.*(1 .+ slk_Imax_c))
        else
            @expression(SCOPF,CImaxlimit, Imaxlimit0)
        end
        linear_Imaxindex = [LinearIndices(Line_flow_c_Line)[ind] for ind in Imaxindex]
        @constraint(SCOPF,CLine_flow1[ind in linear_Imaxindex], Line_flow_c_Line[ind]<= CImaxlimit[ind])

        block=zeros(AffExpr,raw["nBus"])
        ΔV=fill(block,nTP,cl);Δθ=fill(block,nTP,cl);
        V=zeros(AffExpr,nTP,raw["nBus"],cl);θ=zeros(AffExpr,nTP,raw["nBus"],cl)
        V.=V_0;θ.=θ_0

        for t=1:nTP,c=1:cl
            ΔV[t,c]=SENSDict["S_VP"][t,c]*xΔPinj[t,:,c]+ SENSDict["S_VQ"][t,c]*xΔQinj[t,:,c]
            V[t,:,c]=V_0[t,:,c]+ΔV[t,c]
            Δθ[t,c]=SENSDict["S_θP"][t,c]*xΔPinj[t,:,c]+ SENSDict["S_θQ"][t,c]*xΔQinj[t,:,c]
            θ[t,:,c]=θ_0[t,:,c]+Δθ[t,c]
        end
        @expression(SCOPF,CV2magexp, V)
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
        ovr_Sens_V=zeros(nTP,raw["nBus"],cl)
        ovr_Sens_V10=zeros(nTP,raw["nBus"],cl)
        for i=1:raw["nBus"],t=1:nTP,c=1:cl
            ovr_Sens_V[t,i,c]=sum(abs.(SENSDict["S_VQ"][t,c][i,:])) *Lininit["TR"][2]+sum(abs.(SENSDict["S_VP"][t,c][i,:]))*Lininit["TR"][1]
            sorted=sort(nonzeros([abs.(SENSDict["S_VQ"][t,c][i,:])*Lininit["TR"][2];abs.(SENSDict["S_VP"][t,c][i,:])]*Lininit["TR"][1]))
            ovr_Sens_V10[t,i,c]=sum(sorted[end-min(length(sorted)-1,10):end])
        end
        Vminlimit0=(raw["Bus"]."V_min (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
        Vmaxlimit0=(raw["Bus"]."V_max (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]

        Vmaxindex1= findall(x->x==1,Vmaxlimit0.-V_0 .<=1*ovr_Sens_V10 .&& ovr_Sens_V.>0.01*Lininit["TR"][2])
        Vmaxindex2= findall(x->x==1,V_0.>=(raw["Bus"]."V_max (pu)".-0.00)'[ones(Int8,nTP),:,ones(Int8,cl)])
        Vmaxindex= union([Vmaxindex1;Vmaxindex2])
        Vminindex1= findall(x->x==1,V_0.- Vminlimit0.<=1*ovr_Sens_V10 .&& ovr_Sens_V.>0.01*Lininit["TR"][2])
        Vminindex2= findall(x->x==1,V_0.<=(raw["Bus"]."V_min (pu)".+0.00)'[ones(Int8,nTP),:,ones(Int8,cl)])
        Vminindex= union([Vminindex1;Vminindex2])

        if relaxpermit==true
            @expression(SCOPF,CVminlimit, Vminlimit0.*(1 .- slk_Vmin_c))
            @expression(SCOPF,CVmaxlimit, Vmaxlimit0.*(1 .+ slk_Vmax_c))
        else
            @expression(SCOPF,CVminlimit, Vminlimit0)
            @expression(SCOPF,CVmaxlimit, Vmaxlimit0)
        end
        linear_Vminindex = [LinearIndices(CV2magexp)[ind] for ind in Vminindex]
        linear_Vmaxindex = [LinearIndices(CV2magexp)[ind] for ind in Vmaxindex]
        @constraint(SCOPF,CV2mag1[ind in linear_Vminindex], CVminlimit[ind]<= CV2magexp[ind])
        @constraint(SCOPF,CV2mag2[ind in linear_Vmaxindex], CV2magexp[ind]<=CVmaxlimit[ind])
        @expression(SCOPF,Cθexp, θ)
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

        @constraint(SCOPF,CTrust_RegionP1,xΔPinj.<=Lininit["TR"][1].*(C_TRfreeP.+1));
        @constraint(SCOPF,CTrust_RegionP2,-Lininit["TR"][1].*(C_TRfreeP.+1).<=xΔPinj);
        @constraint(SCOPF,CTrust_RegionQ1,xΔQinj.<=Lininit["TR"][2].*(C_TRfreeQ.+1));
        @constraint(SCOPF,CTrust_RegionQ2,-Lininit["TR"][2].*(C_TRfreeQ.+1).<=xΔQinj);


        @constraint(SCOPF,CSensitivity_balance_P1[t=1:nTP,c=1:cl],sum(xΔPinj[t,:,c])<=Lininit["TRB"][1]*(C_SBfreeP[t,c]+1));
        @constraint(SCOPF,CSensitivity_balance_Q1[t=1:nTP,c=1:cl],sum(xΔQinj[t,:,c])<=Lininit["TRB"][1]*(C_SBfreeQ[t,c]+1));
        @constraint(SCOPF,CSensitivity_balance_P2[t=1:nTP,c=1:cl],-Lininit["TRB"][1]*(C_SBfreeP[t,c]+1)<=sum(xΔPinj[t,:,c]));
        @constraint(SCOPF,CSensitivity_balance_Q2[t=1:nTP,c=1:cl],-Lininit["TRB"][1]*(C_SBfreeQ[t,c]+1)<=sum(xΔQinj[t,:,c]));

        # fix.(xΔPinj,0;force=true)
        # fix.(xΔQinj,0;force=true)

        # @constraint(SCOPF, CTRLimit1[t=1:nTP,c=1:cl], sum(xΔPinj[t,:,c].^2)<=36)
        # @constraint(SCOPF, CTRLimit2[t=1:nTP,c=1:cl], sum(xΔQinj[t,:,c].^2)<=36)
    end
    return SCOPF
end

function Linear_SENS_CP_new(SCOPF,Vars,Parlist,nTP,raw::Dict{String, Any},States,SENSDict::Dict{String, Matrix{Matrix{Float64}}},Lininit::Dict{String,Any},indices_,acc,TR,Base_TR_N,inj_bound,Worst_Vio)
    # SCOPF,Vars,Parlist,nTP,raw,States,SENSDict,Lininit,indices_,acc,TR,Base_TR_N,inj_bound=LSCOPF,Vars_N,Parlist_N,nTP,raw,"Norm",SENSDict_N,Lininit,Norm_indices,acc,TR,Base_TR_N,inj_bound_N
    # SCOPF,Vars,Parlist,nTP,raw,States,SENSDict,Lininit,indices_,acc,TR,Base_TR_N,inj_bound=LSCOPF,Vars_C,Parlist_C,nTP,raw,"Cont",SENSDict_C,Lininit,Cont_indices,acc,TR,Base_TR_C,inj_bound_C
    #TR_act=TR_act_N
    # TR=TR_act_N
    # if W_vio[1]<=0.1 && W_vio[2]<=0.01
    #     W_vio./=5
    # end
    xΔPinjP,xΔPinjN,xΔQinjP,xΔQinjN=Vars
    xΔPinj=xΔPinjP-xΔPinjN
    xΔQinj=xΔQinjP-xΔQinjN

    xΔPinjm=xΔPinjP+xΔPinjN
    xΔQinjm=xΔQinjP+xΔQinjN

    NIflow2 , I2_0 , V_0 , θ_0=Parlist
    # ovr_Sens_I_P,ovr_Sens_I_Q,ovr_Sens_V_P,ovr_Sens_V_Q=Ovr_Sens
    F2maxtemp=@. (raw["Line"]."RATE_A")^2
    slackbus=raw["Bus"][raw["Bus"]."Bus_Type".==3,"Bus_N"][1]
    # GBus=raw["Gen"]."Bus"
    # GBus=zeros(Int,raw["nBus"])
    # GBus[raw["Gen"]."Bus"].=1

    # TR1=TR_act[1]
    # TR2=TR_act[2]

    TR1=TR[1]
    TR2=TR[2]


    set_upper_bound.(xΔPinjP, TR[1])
    set_upper_bound.(xΔPinjN, TR[1])
    set_upper_bound.(xΔQinjP, 1.5*TR[2])
    set_upper_bound.(xΔQinjN, 1.5*TR[2])

    if  States=="Norm" cl=1 elseif States=="Cont" cl=raw["nConL"] end

    # TR1=zeros(nTP,raw["nPActBus"],cl)
    # TR2=zeros(nTP,raw["nQActBus"],cl)
    # TR1.=TR[1]
    # TR2.=TR[2]
    # TRact_1=TR_act[1]
    # TRact_2=TR_act[2]
    inj_bound_P=inj_bound[1:2]
    inj_bound_Q=inj_bound[3:4]

    # inj_bound_P[2]
    Gen_ind=indexin(raw["Gen"]."Bus", raw["PActBus"])
    # TRact_1[:,Gen_ind,:]=TRact_1[:,Gen_ind,:].+Base_TR_N
    Imaxlimit0=(F2maxtemp')[ones(Int8,nTP),:,ones(Int8,cl)]
    Imaxlimit0_sqrt=sqrt.(Imaxlimit0)
    accI2=2*acc[1]*Imaxlimit0_sqrt

    if  States=="Norm"
        # @constraint(SCOPF,NTrust_RegionP1,xΔPinj .<= max.(1e-4,inj_bound_P[1]))
        # @constraint(SCOPF,NTrust_RegionP2,min.(-1e-4,inj_bound_P[2]).<= xΔPinj)
        # @constraint(SCOPF,NTrust_RegionQ1,xΔQinj .<= max.(1e-4,inj_bound_Q[1]))
        # @constraint(SCOPF,NTrust_RegionQ2,min.(-1e-4,inj_bound_Q[2]) .<= xΔQinj)
        @constraint(SCOPF,NTrust_RegionP1,xΔPinj .<= TR1)
        @constraint(SCOPF,NTrust_RegionP2,-TR1.<= xΔPinj)
        extraTR2 = similar(xΔQinj)
        extraTR2.=0.0
        shed_elem_Q=(SCOPF[:shed]-Lininit["Pshed"]).*reshape(raw["Load"]."PF" , 1, :)
        for t in 1:nTP, cL in 1:cl, i in raw["nQActBus"]
            load_idx =findall(x->x==raw["QActBus"][i],raw["Load"]."Bus_N")
            extraTR2[t, i, cL] = isempty(load_idx) ? 0.0 : shed_elem_Q[t, load_idx[1], cL]
        end      
          extraTR2.=0.0


        @constraint(SCOPF,NTrust_RegionQ1,xΔQinj .<= (TR2 .- extraTR2))
        @constraint(SCOPF,NTrust_RegionQ2,-(TR2 .- extraTR2) .<= xΔQinj)
        # @constraint(SCOPF,NTrust_RegionQ1,xΔQinj .<= TR2)
        # @constraint(SCOPF,NTrust_RegionQ2,-TR2 .<= xΔQinj)

        # @constraint(SCOPF,NMax_redispatchP[t=1:nTP,c=1:cl],sum(xΔPinjm[t,:,c])<=50)
        # @constraint(SCOPF,NMax_redispatchQ[t=1:nTP,c=1:cl],sum(xΔQinjm[t,:,c])<=8)
        @constraint(SCOPF,Sensitivity_balance_P1[t=1:nTP,c=1:cl],sum(xΔPinj[t,:,c])<=(TR[1]/10000))
        @constraint(SCOPF,Sensitivity_balance_Q1[t=1:nTP,c=1:cl],sum(xΔQinj[t,:,c])<=(TR[2]*raw["nBus"]*0.1))
        @constraint(SCOPF,Sensitivity_balance_P2[t=1:nTP,c=1:cl],-(TR[1]/10000)<=sum(xΔPinj[t,:,c]))
        @constraint(SCOPF,Sensitivity_balance_Q2[t=1:nTP,c=1:cl],-(TR[2]*raw["nBus"]*0.1)<=sum(xΔQinj[t,:,c]))
    elseif States=="Cont"
        @constraint(SCOPF,CTrust_RegionP1,xΔPinj .<= TR1)
        @constraint(SCOPF,CTrust_RegionP2,-TR1.<= xΔPinj)

        extraTR2 = similar(xΔQinj)
        extraTR2.=0.0
        shed_elem_Q=(SCOPF[:shed_c_Line]-Lininit["Pshed_c_Line"]).*reshape(raw["Load"]."PF" , 1, :)
        for t in 1:nTP, cL in 1:cl, i in 1:raw["nQActBus"]
            load_idx =findall(x->x==raw["QActBus"][i],raw["Load"]."Bus_N")
            extraTR2[t, i, cL] = isempty(load_idx) ? 0.0 : shed_elem_Q[t, load_idx[1], cL]
        end
                extraTR2.=0.0


        # shed_elem_Q[19,14,21]
        # extraTR2[19,5,21]

        @constraint(SCOPF,CTrust_RegionQ1,xΔQinj .<= (TR2 .- extraTR2)) #TR2)
        @constraint(SCOPF,CTrust_RegionQ2,-(TR2 .- extraTR2) .<= xΔQinj)

        # @constraint(SCOPF,CTrust_RegionQ1,xΔQinj .<= TR2)
        # @constraint(SCOPF,CTrust_RegionQ2,-TR2 .<= xΔQinj)
        # @constraint(SCOPF,CMax_redispatchP[t=1:nTP,c=1:cl],sum(xΔPinjm[t,:,c])<=50)
        # @constraint(SCOPF,CMax_redispatchQ[t=1:nTP,c=1:cl],sum(xΔQinjm[t,:,c])<=8)
        @constraint(SCOPF,CSensitivity_balance_P1[t=1:nTP,c=1:cl],sum(xΔPinj[t,:,c])<=(TR[1]/10000))
        @constraint(SCOPF,CSensitivity_balance_Q1[t=1:nTP,c=1:cl],sum(xΔQinj[t,:,c])<=(TR[2]*raw["nBus"]*0.1))
        @constraint(SCOPF,CSensitivity_balance_P2[t=1:nTP,c=1:cl],-(TR[1]/10000)<=sum(xΔPinj[t,:,c]))
        @constraint(SCOPF,CSensitivity_balance_Q2[t=1:nTP,c=1:cl],-(TR[2]*raw["nBus"]*0.1)<=sum(xΔQinj[t,:,c]))
    end
    # ovr_Sens_I=ovr_Sens_I_P*(TR1/10) + ovr_Sens_I_Q*(TR2/1.1)
    # ovr_Sens_V= ovr_Sens_V_P*TR1+ovr_Sens_V_Q*TR2
    # ovr_Sens_V= ovr_Sens_V_Q*(TR2/1.1) #+ ovr_Sens_V_P*TR1
    # ovr_Sens_V_G= ovr_Sens_V_Q*TR2+ ovr_Sens_V_P*TR1


    # minline=Int(floor(max(0.02*raw["nLine"],10)))+1

    # num0=sum(I2_0.>=Imaxlimit0)
    # num1=0
    # num=1000
    # TR_adjust=1
    # I2_1= I2_0 .+ (ovr_Sens_I_Q.*(TR2/1.1))
    # I2_2= I2_0 .- (ovr_Sens_I_Q.*(TR2/1.1))
    # while num>3*nTP*cl && num1==0
    #     index=I2_1+(ovr_Sens_I_P.*(TR1/(10*TR_adjust))) .>=Imaxlimit0
    #     index1=I2_2.-(ovr_Sens_I_P*(TR1/(10*TR_adjust))) .>=Imaxlimit0
    #     num1=sum(index1)
    #     num=sum(index)-num0
    #     TR_adjust += 0.1
    # end

    # .*(TR1/10)

    # minline=1
    # minline=3107
    # Imaxindex,elements=con_indices((Imaxlimit0.-I2_0),ovr_Sens_I, minline,nTP,cl)
    # CartesianIndex(22,46,30) in Imaxindex
    # LinearIndices(I2_0)[22,46,30]
    Imaxindex=indices_[1]
    linear_Imaxindex = LinearIndices(I2_0)[Imaxindex]
    Imaxlimit=Imaxlimit0[linear_Imaxindex]

    # findall(x->x==CartesianIndex(22,46,30),Imaxindex) --->6426
    # findall(x->x==CartesianIndex(22,46,30),Imaxindex) --->6426
    # findall(x-> x==42971,linear_Imaxindex)
    # linear_Imaxindex[75]
    # I2_0[Imaxindex[6426]]

    # CartesianIndex(22,46,30) in SCOPF[:slkDataFrame]["slkImax"]

    # findall(x->x==CartesianIndex(22,46,30),SCOPF[:slkDataFrame]["slkImax"])

    # LinearIndices(I2_0)[CartesianIndex(3,2265,2)]
    # IndexDictI=Dict(tuple.(linear_Imaxindex,tuple.(getindex.(Imaxindex,1),getindex.(Imaxindex,3))))
    # IndexDictI[42971]
    # ovr_Sens_I_P[1,2665,1]
    # (ovr_Sens_I_P*(TR1))[1,2665,1]
    # linear_Imaxindex
    # 6426 in linear_Imaxindexrelax
    linear_Imaxindexrelax=findall(I2_0[linear_Imaxindex].-Imaxlimit .> 1e-6)
    # 6426 in linear_Imaxindexrelax
    # unregister.(SCOPF,:slk_Imax)
    if relaxpermit==true
        if  States=="Norm"
            SCOPF[:slkDataFrame]["slkImax"]= filter(x->(I2_0.-Imaxlimit0)[x] .> 1e-6,Imaxindex)
            @variable(SCOPF,1>=slk_Imax[ind in 1:length(linear_Imaxindexrelax)]>=0)  #Worst_Vio[1]>=
            # @variable(SCOPF,1>=slk_Imax[ind in linear_Imaxindex]>=0)  #Worst_Vio[1]>=
            # @variable(SCOPF, max(SCOPF[:W_ramp],maximum(Worst_Vio))>=slk_Norm>=0)
            # @constraint(SCOPF, Norm_slack1,slk_Norm.>=slk_Imax)
            slkall=zeros(AffExpr,length(linear_Imaxindex))
            slkall[linear_Imaxindexrelax].=slk_Imax.*(I2_0[linear_Imaxindex][linear_Imaxindexrelax].-Imaxlimit[linear_Imaxindexrelax])

            Imaxlimit= Imaxlimit.+slkall
            # Imaxlimit = Imaxlimit.+slk_Imax.* max.(0,I2_0[linear_Imaxindex]-Imaxlimit)

            # collect(slk_Imax)
            # sum(I2_0[linear_Imaxindex].-Imaxlimit .>=0)
        elseif States=="Cont"
            SCOPF[:slkDataFrame]["slkImax_c"]= filter(x->(I2_0.-Imaxlimit0)[x] .> 1e-6,Imaxindex)
            @variable(SCOPF, 1>=slk_Imax_c[ind in 1:length(linear_Imaxindexrelax)]>=0) #Worst_Vio[1]>=
            # @variable(SCOPF, max(SCOPF[:W_ramp],maximum(Worst_Vio))>=slk_Cont>=0)
            # @constraint(SCOPF, Cont_slack1,slk_Cont.>=slk_Imax_c)
            # Imaxlimit .*= (1 .+ slk_Imax_c)
            # Imaxlimit = (1 .+ slk_Imax_c).*Imaxlimit
            slkall=zeros(AffExpr,length(linear_Imaxindex))
            slkall[linear_Imaxindexrelax].=slk_Imax_c.*(I2_0[linear_Imaxindex][linear_Imaxindexrelax].-Imaxlimit[linear_Imaxindexrelax])
            Imaxlimit= Imaxlimit.+slkall

            # Imaxlimit=Imaxlimit.+ slk_Imax_c.* max.(0,(I2_0[linear_Imaxindex].-Imaxlimit))
        end

    end
    # Imaxlimit[6426]
    # Imaxlimit[75]
    I2=AffExpr.(I2_0)
    elements = group(t -> (t[1], t[3]), getfield.(Imaxindex,1))
    # CartesianIndices(elements)
    # Threads.@threads
    # ind=collect(keys(elements))[1]
    # ind=(22,30)
    Threads.@threads for ind in  collect(keys(elements))
            lines=getfield.(elements[ind],2)
            t,c=ind
            # SENSDict["S_IP"][t,c][lines,:].*(TR1[t,:,c])'
            # Pcoef=filter_to_sparse(SENSDict["S_IP"][t,c][lines,:], x -> abs(x) .> (acc[1]/100),(TRact_1[t,:,c])' ./(2*Imaxlimit0_sqrt[t,lines,c]))
            signA=sign.(SENSDict["S_IP"][t,c][lines,:])
            Pos=((signA.+1)./2).*(inj_bound_P[1][t,:,c])'
            Neg=((signA.-1)./2).*(inj_bound_P[2][t,:,c])'
            PosC=((signA.+1)./2).*(inj_bound_P[2][t,:,c])'
            NegC=((signA.-1)./2).*(inj_bound_P[1][t,:,c])'
            TRact_1=Pos+Neg
            TRact_1_C=PosC+NegC

            signA=sign.(SENSDict["S_IQ"][t,c][lines,:])
            Pos=((signA.+1)./2).*(inj_bound_Q[1][t,:,c])'
            Neg=((signA.-1)./2).*(inj_bound_Q[2][t,:,c])'
            PosC=((signA.+1)./2).*(inj_bound_Q[2][t,:,c])'
            NegC=((signA.-1)./2).*(inj_bound_Q[1][t,:,c])'
            TRact_2=Pos+Neg
            TRact_2_C=PosC+NegC
            # Pcoef=filter_to_sparse_V_2(SENSDict["S_IP"][t,c][lines,:],filter_acc,TRact_1;neighborhood=80,multiplier2=TR1,neighborhood2=30)#
            # Qcoef=filter_to_sparse_V_2(SENSDict["S_IQ"][t,c][lines,:],filter_acc,TRact_2;neighborhood=30,multiplier2=TR2,neighborhood2=30)
            # Pcoef=filter_to_sparse_V_3(SENSDict["S_IP"][t,c][lines,:],accI2[t,lines,c],[TRact_1, TRact_1_C];neighborhood=[30,30])#
            # Qcoef=filter_to_sparse_V_3(SENSDict["S_IQ"][t,c][lines,:],accI2[t,lines,c],[TRact_2, TRact_2_C];neighborhood=[30,30])
            # Pcoef=filter_to_sparse_V_3(SENSDict["S_IP"][t,c][lines,:],0,[TRact_1, TRact_1_C];neighborhood=[90,30])#
            # Qcoef=filter_to_sparse_V_3(SENSDict["S_IQ"][t,c][lines,:],0,[TRact_2, TRact_2_C];neighborhood=[90,30])
            # Pcoef=filter_to_sparse_V_3(SENSDict["S_IP"][t,c][lines,:],0.0001,[TR1, TR1])#
            # Qcoef=filter_to_sparse_V_3(SENSDict["S_IQ"][t,c][lines,:],0.0001,[TR2, TR2])
            # filter_percent=0.05 / ((TR[1]+TR[2]>=0.5 || sum(acc)>1e-2) ? 1 : 10)
            # filter_percent=min(0.01+ 0.05*min(TR[1],8)/2, 0.01+ 0.05*min(TR[2],0.04)/0.005)
            # filter_percent=min(0.005+ 0.1*min(TR[1],2)/2, 0.005+ 0.1*min(TR[2],0.02)/0.02)
            # filter_percent=min(0.05*min(TR[1],8), 10*min(TR[2],0.04))

            # filter_percent=min(0.05*TR[1], 5*TR[2])
            # filter_percent=max(min(0.005*TR[1], 1*TR[2]),0.02)

            # filter_percent=max(0.01,min(0.02*TR[1], 10*TR[2]))

            # filter_percent=10* acc[1]
            filter_percent=0.001*TR[1]
            # filter_percent=4*min(0.05*min(TR[1],8), 10*min(TR[2],0.04))


            # Pcoef,Qcoef=filter_to_sparse_V_5([SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],filter_percent,[TRact_1, TRact_1_C,TRact_2, TRact_2_C])
            # Pcoef,Qcoef=filter_to_sparse_V_5([SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],filter_percent,[TRact_1, TRact_1_C,TRact_2, TRact_2_C])
            Pcoef,Qcoef=filter_to_sparse_V_13([SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],filter_percent,[TRact_1, TRact_1_C,TRact_2, TRact_2_C])

            # sum(Pcoef[5,:].*(value.(model[:ΔPinjP_c_Line]).-value.(model[:ΔPinjN_c_Line]))[22,:,30])
            # Qcoef[5,:]

            #     SENSDict["S_IP"][t,c][46,:]
            #     SENSDict_C["S_IP"][22,30][46,:]

            #     SENSDict["S_IQ"][t,c][46,:]
            #     SENSDict_C["S_IQ"][22,30][46,:]


            #     aa=SENSDict_C["S_IP"][22,30][46,:].*(value.(model[:ΔPinjP_c_Line]).-value.(model[:ΔPinjN_c_Line]))[22,:,30]
            #     # bb=SENSDict_C["S_IQ"][22,30][46,:].*(value.(model[:ΔQinjP_c_Line]).-value.(model[:ΔQinjN_c_Line]))[22,:,30]
            #     # sum(bb)
            #     # sum(aa)
            # error=1e-5
            # Pcoef,Qcoef=filter_to_sparse_V_12([SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],error,[TRact_1, TRact_1_C,TRact_2, TRact_2_C])

            # Pcoef=round.(SENSDict["S_IP"][t,c][lines,:],digits=5)
            # Qcoef=round.(SENSDict["S_IQ"][t,c][lines,:],digits=5)

            # Pcoef,Qcoef=filter_to_sparse_V_10([SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],filter_percent,[TRact_1, TRact_1_C,TRact_2, TRact_2_C])

            # maximum(abs.(Pcoef))
            # minimum(abs.(Pcoef))
            # maximum(abs.(Qcoef))
            # minimum(abs.(Qcoef))


            SCOPF[:SENS_filter_term_I]+=(nnz(Pcoef)+ nnz(Qcoef))
            SCOPF[:SENS_term_I]+= (length(SENSDict["S_IP"][t,c][lines,:])+ length(SENSDict["S_IQ"][t,c][lines,:]))

            # Pcoef=filter_to_sparse_V_4(SENSDict["S_IP"][t,c][lines,:],0.05,[TRact_1, TRact_1_C])
            # Qcoef=filter_to_sparse_V_4(SENSDict["S_IQ"][t,c][lines,:],0.05,[TRact_2, TRact_2_C])
            ΔI2=Pcoef*xΔPinj[t,:,c] + Qcoef*xΔQinj[t,:,c]
            I2[t,lines,c]+=ΔI2
    end
    # I2[3,2265,2]

    # A=SENSDict["S_IP"][t,c][lines,:]
    # multiplier=(TRact_1[t,:,c])' ./(2*Imaxlimit0_sqrt[t,lines,c])
    # devs=A.*multiplier
    # AAA=sparse(delete_elements!(A.*multiplier,50*acc[1]))
    # return (AAA.!=0).*A
    # I2[22,46,30]
    # I2[40774]
    # linear_Imaxindex[6426]

    if  States=="Norm" @expression(SCOPF,Line_flow,I2) elseif States=="Cont" @expression(SCOPF,Line_flow_c_Line,I2) end
    # Line_flow[Imaxindex]
    if  States=="Norm"
        @constraint(SCOPF,NLine_flow1, Line_flow[linear_Imaxindex].<= Imaxlimit)
    elseif States=="Cont"
        @constraint(SCOPF,CLine_flow1, Line_flow_c_Line[linear_Imaxindex].<= Imaxlimit)
    end

    # CLine_flow1[6426]

    # value(model[:CLine_flow1][6426])
    # value(model[:slk_Imax_c][751])
    # value.(model[:Line_flow_c_Line][40774])

    # CLine_flow1[75]
    # Imaxlimit[75]
    # Line_flow_c_Line[linear_Imaxindex[75]]
    # # Voltage constraints
    Vminlimit0=(raw["Bus"]."V_min (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
    Vmaxlimit0=(raw["Bus"]."V_max (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
    # Vminlimit0[:,Not(Gen_ind),:] .-= acc[2]
    # Vmaxlimit0[:,Not(Gen_ind),:] .+= acc[2]
    # minbus=Int(floor(max(0.02*raw["nBus"],10)))+1
    # minbus=1
    # Vmaxindex,Vmaxelements=con_indices((Vmaxlimit0.-V_0),ovr_Sens_V, minbus,nTP,cl,ovr_Sens_V_G,GBus)
    # Vminindex,Vminelements=con_indices((V_0.- Vminlimit0),ovr_Sens_V, minbus,nTP,cl,ovr_Sens_V_G,GBus)
    Vmaxindex,Vminindex= indices_[2], indices_[3]
    linear_Vminindex = LinearIndices(V_0)[Vminindex]
    linear_Vmaxindex = LinearIndices(V_0)[Vmaxindex]
    Vindex=union(Vmaxindex,Vminindex)

    Vmaxlimit=Vmaxlimit0[linear_Vmaxindex]
    Vminlimit=Vminlimit0[linear_Vminindex]

    # IndexDictVmax=Dict(tuple.(linear_Vmaxindex,tuple.(getindex.(Vmaxindex,1),getindex.(Vmaxindex,3))))
    # IndexDictVmin=Dict(tuple.(linear_Vminindex,tuple.(getindex.(Vminindex,1),getindex.(Vminindex,3))))
    linear_Vmaxindexrelax=findall(V_0[linear_Vmaxindex].-Vmaxlimit .> 1e-5)
    linear_Vminindexrelax=findall(Vminlimit .- V_0[linear_Vminindex] .> 1e-5)
    # setdiff!(Imaxindex,findall(x->x==0,NIflow2-I2_0))
    # unregister.(SCOPF,:slk_Vmax)
    # unregister.(SCOPF,:slk_Vmin)

    if relaxpermit==true
        if  States=="Norm"
            SCOPF[:slkDataFrame]["slkVmax"]= filter(x->(V_0.-Vmaxlimit0)[x] .> 1e-5,Vmaxindex)
            SCOPF[:slkDataFrame]["slkVmin"]= filter(x->(Vminlimit0.-V_0)[x] .> 1e-5,Vminindex)
            @variable(SCOPF, 1>=slk_Vmax[ind in 1:length(linear_Vmaxindexrelax)]>=0) #Worst_Vio[2]>=
            @variable(SCOPF, 1>=slk_Vmin[ind in 1:length(linear_Vminindexrelax)]>=0) #Worst_Vio[2]>=
            # @constraint(SCOPF, Norm_slack2, slk_Norm.>=slk_Vmax)
            # @constraint(SCOPF, Norm_slack3, slk_Norm.>=slk_Vmin)
            # Vmaxlimit = (1 .+ slk_Vmax).*Vmaxlimit
            # Vminlimit = (1 .- slk_Vmin).*Vminlimit
            slkall=zeros(AffExpr,length(linear_Vmaxindex))
            slkall[linear_Vmaxindexrelax] .= slk_Vmax.*max.(0.001,(V_0[linear_Vmaxindex][linear_Vmaxindexrelax].-Vmaxlimit[linear_Vmaxindexrelax]))
            Vmaxlimit= Vmaxlimit.+slkall

            slkall=zeros(AffExpr,length(linear_Vminindex))
            slkall[linear_Vminindexrelax] .= slk_Vmin.*max.(0.001,(Vminlimit[linear_Vminindexrelax].-V_0[linear_Vminindex][linear_Vminindexrelax]))
            Vminlimit=Vminlimit.-slkall

            # Vmaxlimit=Vmaxlimit.+ slk_Vmax .* max.(0,(V_0[linear_Vmaxindex].-Vmaxlimit))
            # Vminlimit=Vminlimit.- slk_Vmin .* max.(0,(Vminlimit .- V_0[linear_Vminindex]))
            # Vmaxlimit=Vmaxlimit.+slk_Vmax
            # Vminlimit=Vminlimit.-slk_Vmin
        elseif States=="Cont"
            SCOPF[:slkDataFrame]["slkVmax_c"]= filter(x->(V_0.-Vmaxlimit0)[x] .> 1e-5,Vmaxindex)
            SCOPF[:slkDataFrame]["slkVmin_c"]= filter(x->(Vminlimit0.-V_0)[x] .> 1e-5,Vminindex)
            @variable(SCOPF, 1>=slk_Vmax_c[ind in 1:length(linear_Vmaxindexrelax)]>=0) #Worst_Vio[2]>=
            @variable(SCOPF, 1>=slk_Vmin_c[ind in 1:length(linear_Vminindexrelax)]>=0) #Worst_Vio[2]>=
            # @constraint(SCOPF, Cont_slack2,slk_Cont.>=slk_Vmax_c)
            # @constraint(SCOPF, Cont_slack3,slk_Cont.>=slk_Vmin_c)
            # Vmaxlimit = (1 .+ slk_Vmax_c).*Vmaxlimit
            # Vminlimit = (1 .- slk_Vmin_c).*Vminlimit
            slkall=zeros(AffExpr,length(linear_Vmaxindex))
            slkall[linear_Vmaxindexrelax] .= slk_Vmax_c.*max.(0.001,(V_0[linear_Vmaxindex][linear_Vmaxindexrelax].-Vmaxlimit[linear_Vmaxindexrelax]))
            Vmaxlimit= Vmaxlimit.+slkall

            slkall=zeros(AffExpr,length(linear_Vminindex))
            slkall[linear_Vminindexrelax] .= slk_Vmin_c.*max.(0.001,(Vminlimit[linear_Vminindexrelax].-V_0[linear_Vminindex][linear_Vminindexrelax]))
            Vminlimit=Vminlimit.-slkall

            # Vmaxlimit=Vmaxlimit.+slk_Vmax_c .* max.(0,(V_0[linear_Vmaxindex].-Vmaxlimit))
            # Vminlimit=Vminlimit.-slk_Vmin_c .* max.(0,(Vminlimit.- V_0[linear_Vminindex]))
            # Vmaxlimit=Vmaxlimit.+slk_Vmax_c
            # Vminlimit=Vminlimit.-slk_Vmin_c
        end
    end

    V=AffExpr.(V_0)
    # V2=AffExpr.(V_0);θ2=AffExpr.(θ_0)
    #
    # union(getfield.(Vminindex,1),getfield.(Vmaxindex,1))
    Velements = group(t -> (t[1], t[3]),  union(getfield.(Vminindex,1),getfield.(Vmaxindex,1)))
    # (5,2470,1) in Velements[(5,1)]
    # key_=collect(keys(Velements))[2]
    # key_= (14,30)
    # CartesianIndex(14,60,30) in Vmaxindex
    # elems[1:15]
    Threads.@threads for key_ in collect(keys(Velements))
            elems=getfield.(Velements[key_],2)
            t,c=key_
            # Gelems=intersect(GBus,elems)
            signA=sign.(SENSDict["S_VP"][t,c][elems,:])
            Pos=((signA.+1)./2).*(inj_bound_P[1][t,:,c])'
            Neg=((signA.-1)./2).*(inj_bound_P[2][t,:,c])'

            PosC=((signA.+1)./2).*(inj_bound_P[2][t,:,c])'
            NegC=((signA.-1)./2).*(inj_bound_P[1][t,:,c])'
            TRact_1=Pos+Neg
            TRact_1_C=PosC+NegC

            # sum(abs.(Pos).>0.001 .&& abs.(Neg).>0.001)

            # TRact_1_C=PosC+NegC
            signA=sign.(SENSDict["S_VQ"][t,c][elems,:])
            Pos=((signA.+1)./2).*(inj_bound_Q[1][t,:,c])'
            Neg=((signA.-1)./2).*(inj_bound_Q[2][t,:,c])'
            PosC=((signA.+1)./2).*(inj_bound_Q[2][t,:,c])'
            NegC=((signA.-1)./2).*(inj_bound_Q[1][t,:,c])'
            TRact_2=Pos+Neg
            TRact_2_C=PosC+NegC

            # Pcoef_V=filter_to_sparse(SENSDict["S_VP"][t,c][elems,:], x -> abs(x) >  (acc[2]/5),(TRact_1[t,:,c])')
            # Pcoef_V=filter_to_sparse_V_2(SENSDict["S_VP"][t,c][elems,:],acc[2],TRact_1;neighborhood=50)
            # Pcoef_V=filter_to_sparse_V_3(SENSDict["S_VP"][t,c][elems,:],acc[2],[TRact_1,TRact_1_C];neighborhood=[30,30])
            # Pcoef_V=filter_to_sparse_V_3(SENSDict["S_VP"][t,c][elems,:],0,[TRact_1,TRact_1_C];neighborhood=[90,30])
            # Pcoef_V=filter_to_sparse_V_3(SENSDict["S_VP"][t,c][elems,:],0.0001,[TR1,TR1])
            # Qcoef_V=filter_to_sparse(SENSDict["S_VQ"][t,c][elems,:], x -> abs(x) >  (acc[2]/5),(TRact_2[t,:,c])')
            # Qcoef_V=filter_to_sparse_V_2(SENSDict["S_VQ"][t,c][elems,:],acc[2],TRact_2;neighborhood=50) #;neighborhood=10
            # Qcoef_V=filter_to_sparse_V_2(SENSDict["S_VQ"][t,c][elems,:],acc[2],TRact_2;neighborhood=30,multiplier2=TR2,neighborhood2=30) #;neighborhood=10
            # Qcoef_V=filter_to_sparse_V_3(SENSDict["S_VQ"][t,c][elems,:],acc[2],[TRact_2,TRact_2_C];neighborhood=[30,30]) #;neighborhood=10
            # Qcoef_V=filter_to_sparse_V_3(SENSDict["S_VQ"][t,c][elems,:],0,[TRact_2,TRact_2_C];neighborhood=[90,30]) #;neighborhood=10
            # Qcoef_V=filter_to_sparse_V_3(SENSDict["S_VQ"][t,c][elems,:],0.0001,[TR2,TR2]) #;neighborhood=10

            # filter_percent=max(0.001,(min(0.5*TR[2],0.001)))
            # # filter_percent=10*acc[2]
            filter_percent=0.0001*TR[2]

# evenly distribute
            # 0.0001
            # 0.04


            # filter_percent=min(0.01+ 0.1*min(TR[1],4)/2, 0.01+ 0.1*min(TR[2],0.04)/0.02)

            # Pcoef_V,Qcoef_V=filter_to_sparse_V_5([SENSDict["S_VP"][t,c][elems,:],SENSDict["S_VQ"][t,c][elems,:]],filter_percent,[TRact_1, TRact_1_C,TRact_2, TRact_2_C])

            Pcoef_V,Qcoef_V=filter_to_sparse_V_13([SENSDict["S_VP"][t,c][elems,:],SENSDict["S_VQ"][t,c][elems,:]],filter_percent,[TRact_1, TRact_1_C,TRact_2, TRact_2_C])
            # error=1e-10
            # Pcoef_V,Qcoef_V=filter_to_sparse_V_12([SENSDict["S_VP"][t,c][elems,:],SENSDict["S_VQ"][t,c][elems,:]],error,[TRact_1, TRact_1_C,TRact_2, TRact_2_C])

            # Pcoef_V=sparse(round.(SENSDict["S_VP"][t,c][elems,:],digits=5))
            # Pcoef_V=sparse(round.(SENSDict["S_VQ"][t,c][elems,:],digits=5))

            # aa=hcat(Pcoef_V,Qcoef_V)[11,:]
            # bb=hcat(SENSDict["S_VP"][t,c][elems,:],SENSDict["S_VQ"][t,c][elems,:])[11,:]

            # bb.-aa

            SCOPF[:SENS_filter_term_V]+=(nnz(Pcoef_V)+ nnz(Qcoef_V))
            SCOPF[:SENS_term_V]+= (length(SENSDict["S_VP"][t,c][elems,:])+ length(SENSDict["S_VQ"][t,c][elems,:]))

            # Pcoef_V=filter_to_sparse_V_4(SENSDict["S_VP"][t,c][elems,:],0.05,[TRact_1, TRact_1_C])
            # Qcoef_V=filter_to_sparse_V_4(SENSDict["S_VQ"][t,c][elems,:],0.05,[TRact_2, TRact_2_C]) #;neighborhood=10

            # maximum(abs.(Qcoef_V))
            # sum(SENSDict["S_VP"][t,c][elems,:].*TRact_1,dims=2)
            # sum(SENSDict["S_VQ"][t,c][elems,:].*TRact_2,dims=2)
            # findall(x->x==2470,elems)
            # Pcoef_V[124,:]
            # maximum(SENSDict["S_VQ"][t,c][elems,:])
            # Pcoef_V=SENSDict["S_VP"][t,c][elems,:]
            # Qcoef_V=SENSDict["S_VQ"][t,c][elems,:]
            # Pcoef_V[ Gelems,:]=SENSDict["S_VP"][t,c][Gelems,:]
            # Qcoef_V[ Gelems,:]=SENSDict["S_VQ"][t,c][Gelems,:]
            # findall(x->x==CartesianIndex(5,2470,1), Vminindex)

            ΔV=Pcoef_V*xΔPinj[t,:,c].+ Qcoef_V*xΔQinj[t,:,c]
            V[t,elems,c].+=ΔV
    end

    if  States=="Norm" @expression(SCOPF,NV2magexp, V) elseif States=="Cont" @expression(SCOPF,CV2magexp, V) end

    if  States=="Norm"
        @constraint(SCOPF,NV2mag1, Vminlimit.<= NV2magexp[linear_Vminindex])
        @constraint(SCOPF,NV2mag2, NV2magexp[linear_Vmaxindex].<=Vmaxlimit)
        # @expression(SCOPF,Nθexp, θ)
        # @constraint(SCOPF,NslackVangle, Nθexp[:,slackbus,:].==0)
    elseif States=="Cont"
        @constraint(SCOPF,CV2mag1, Vminlimit.<= CV2magexp[linear_Vminindex])
        @constraint(SCOPF,CV2mag2, CV2magexp[linear_Vmaxindex].<=Vmaxlimit)
        # @expression(SCOPF,Cθexp, θ)
        # @constraint(SCOPF,CslackVangle, Cθexp[:,slackbus,:].==0)
    end
    # model[:NV2magexp][5,2470,1]
    # NV2magexp[linear_Vminindex[211]]
    return SCOPF
end

function Sensitivity_constraints(Glob, LSCOPF,CascadeInfo,N_Var,C_Var,OPFon,raw,Lininit,ctdall)
    Td=Glob["Td"]
    CCT=Glob["CCT"]
    nTP=Glob["nTL"]
    
    SENSDict_All=Dict()
    Parlist_N=N_Var[14];Vars_N=N_Var[13]
    Parlist_C=C_Var[11];Vars_C=C_Var[10]

    if OPFon==0
        I2_0=cat(Parlist_N[2],Parlist_C[2],dims=3);V_0=cat(Parlist_N[3],Parlist_C[3],dims=3);θ_0=cat(Parlist_N[4],Parlist_C[4],dims=3)
        MaxPup_N,MinPup_N,MaxQup_N,MinQup_N=N_Var[15]
        MaxPup_C,MinPup_C,MaxQup_C,MinQup_C=C_Var[12]
        PBound=[cat(MaxPup_N[:,raw["PActBus"],:],MaxPup_C[:,raw["PActBus"],:];dims=3),cat(MinPup_N[:,raw["PActBus"],:],MinPup_C[:,raw["PActBus"],:];dims=3)]
        QBound=[cat(MaxQup_N[:,raw["QActBus"],:],MaxQup_C[:,raw["QActBus"],:];dims=3),cat(MinQup_N[:,raw["QActBus"],:],MinQup_C[:,raw["QActBus"],:];dims=3)]
        SENSDict_N=Cal_SENS_new2(Lininit,ctdall,nTP,raw;C="N");
        SENSDict_C=Cal_SENS_new2(Lininit,ctdall,nTP,raw;C="CL");
        # SENSDict_N_0==SENSDict_N
        # SENSDict_C_0==SENSDict_C
        for key in keys(SENSDict_N)
            SENSDict_All[key]=hcat(SENSDict_N[key],SENSDict_C[key])
        end
        cl=1+raw["nConL"]
    else
        # SENSDict_All=deepcopy(SENSDict_N)
        MaxPup_N,MinPup_N,MaxQup_N,MinQup_N=N_Var[15]
        PBound=[MaxPup_N[:,raw["PActBus"],:],MinPup_N[:,raw["PActBus"],:]]
        QBound=[MaxQup_N[:,raw["QActBus"],:],MinQup_N[:,raw["QActBus"],:]]
        SENSDict_All=Cal_SENS_new2(Lininit,ctdall,nTP,raw;C="N");
        # JDict_All=Lininit["NJacob"]
        I2_0=Parlist_N[2];V_0=Parlist_N[3];θ_0=Parlist_N[4]
        cl=1
    end

    # function find_ramp_rate_violations(Lininit,raw)
    if size(Lininit["PG"],1)>1
        PreT = vcat(Lininit["PG"][2:end,:], Lininit["PG"][end:end,:])  # Append row of zeros
        nextT=vcat( Lininit["PG"][1:1,:],Lininit["PG"][1:end-1,:])

        adiff=nextT-Lininit["PG"]
        bdiff=Lininit["PG"]-PreT

        Base_TR_N=round.(max.(abs.(adiff).- Td*(raw["Gen"].Ramp)',abs.(bdiff).- Td*(raw["Gen"].Ramp)',abs.(adiff+bdiff).-(2*Td*raw["Gen"].Ramp)',0),digits=4)
        if OPFon==0
            cdiff=Lininit["PG_c_Line"].-Lininit["PG"]

            Base_TR_C=round.(max.(abs.(cdiff).-CCT*(raw["Gen"].Ramp)',0),digits=4)
        else
            Base_TR_C=zeros(size(Lininit["PG"]))
        end
    else
        Base_TR_N=zeros(size(Lininit["PG"]))
        Base_TR_C=zeros(size(Lininit["PG"]))
    end
    #     return Base_TR_N,Base_TR_C
    # end

    # maximum(Base_TR_N)


    # plot(SENSDict_N["S_VP"][1][:,100])
    # plot!(SENSDict_N["S_VQ"][1][:,100])

    # SENSDict_N["S_IP"][1]
    # SENSDict_N["S_VQ"][1]

    # filter_to_sparse( SENSDict_N["S_VP"][1], x -> abs(x) > (0.001/TR[1]))
    # sort(abs.(SENSDict_N["S_VP"][1][1,:]))
    # findmin(PBound[1])
    # findmax(PBound[2])
    # findmin(QBound[1])
    # findmax(QBound[2])
    # raw["PActBus"][2]
    # ctdall[5,:]

    # aa=SENSDict_C["S_IP"][22,31][46,:].*(value.(model[:ΔPinjP_c_Line]).-value.(model[:ΔPinjN_c_Line]))[22,:,30]
    # bb=SENSDict_C["S_IQ"][22,31][46,:].*(value.(model[:ΔQinjP_c_Line]).-value.(model[:ΔQinjN_c_Line]))[22,:,30]
    # sum(bb)
    # sum(aa)
    # V_0[22,46,31]
    # I2_0[22,46,31]

    # @time TR,Indices,acc,Injlimit,Worst_Vio=Set_TR_new(Base_TR_N,Base_TR_C,Lininit,SENSDict_All,I2_0,V_0,raw,PBound,QBound,;cl=cl)
    # @time TR,Indices,acc,Injlimit,Worst_Vio=Set_TR_new2(Base_TR_N,Base_TR_C,SENSDict_All,I2_0,V_0,raw,PBound,QBound,;cl=cl)
    TR,Indices,acc,Injlimit,Worst_Vio=Set_TR_new2(nTP,Base_TR_N,Base_TR_C,SENSDict_All,I2_0,V_0,raw,PBound,QBound,CascadeInfo;cl=cl)



    # TR=[1,0.2]
    # acc=0.0
    # Indices[1]




    # using Cthulhu

    LSCOPF[:Num_Const_I]=nTP*raw["nLine"]* cl
    LSCOPF[:Num_Const_V]=nTP*raw["nBus"]*cl*2
    LSCOPF[:Num_filter_Const_I]=length(Indices[1])
    LSCOPF[:Num_filter_Const_V]=length(Indices[2])+length(Indices[3])

    # TR,Indices,acc,Injlimit,Worst_Vio=[TR_P, TR_Q], [ImaxIndex, VmaxIndex, VminIndex],[acc_I,acc_V],Injlimit,[Pworst_vio,Qworst_vio]
    LSCOPF[:acc]=acc
    LSCOPF[:W_VIO]=Worst_Vio
    # LSCOPF[:W_ramp]=max(maximum(Base_TR_N),maximum(Base_TR_C))


    # #
    # findmax(x->x>0.1,value.(model[:ΔPinjP]-model[:ΔPinjN]).-inj_bound_N[1])


    # value.(model[:ΔPinjP]-model[:ΔPinjN])[1,201]
    # inj_bound_N[1][1,201]
    # PBound[1][1,201,1]
    # raw["PActBus"][201]
    # ctdall[2123,:]
    # raw["Load"][1477,:]
    # OP_init["Pshed"][:,1477]
    # Indices[1]
    # CartesianIndex(18,43,2) in Indices[]
    # Indices[]
    # CartesianIndex(8,2740,7) in Indices[2]
    # Indices[3]
    # PQlimits=C_Var[12]
    # N_Var[15][1]
    # N_Var[15][2]
    # maxP,minP,maxQ,minQ= N_Var[15]
    # Pbus=raw["PActBus"],Qbus=raw["QActBus"]
    # P=maxP[:,Pbus];P=minP[:,Pbus];maxQ=maxQ[:,Qbus];minQ=minQ[:,Qbus]
    Lininit["TR"]=TR
    # Ovr_Sens_N=[Ovr_Sens[i][:,:,1] for i in eachindex(Ovr_Sens)]
    # C_LSCOPF=copy(LSCOPF)
    # L1=copy(LSCOPF)
    # LSCOPF=L1
    Norm_indices, Cont_indices=Seprate_Norm_Cont_indices(Indices)
    # Norm_indices[3]
    # CartesianIndex(1,2470,1) in Norm_indices[3]

    # TR_act_N=[TR_act[1][:,:,1],TR_act[2][:,:,1]]
    # TR_act_C=[TR_act[1][:,:,2:end],TR_act[2][:,:,2:end]]

    inj_bound_N=[Injlimit[1][1][:,:,1],Injlimit[1][2][:,:,1],Injlimit[2][1][:,:,1],Injlimit[2][2][:,:,1]]
    inj_bound_C=[Injlimit[1][1][:,:,2:end],Injlimit[1][2][:,:,2:end],Injlimit[2][1][:,:,2:end],Injlimit[2][2][:,:,2:end]]

    LSCOPF[:SENS_filter_term_I]=0;LSCOPF[:SENS_term_I]=0;LSCOPF[:SENS_filter_term_V]=0;LSCOPF[:SENS_term_V]=0
    LSCOPF[:slkDataFrame]=Dict("slkVmax"=>[], "slkVmin"=>[], "slkImax"=>[], "slkImax_c"=>[], "slkVmax_c"=>[], "slkVmin_c"=>[])

    # CartesianIndex(6,183,1) in Norm_indices[2]
    # CartesianIndex(6,2528,2) in Cont_indices[2]
    LSCOPF= Linear_SENS_CP_new(LSCOPF,Vars_N,Parlist_N,nTP,raw,"Norm",SENSDict_N,Lininit,Norm_indices,acc,TR,Base_TR_N,inj_bound_N,Worst_Vio)
    if OPFon==0
        # Ovr_Sens_C=[Ovr_Sens[i][:,:,2:raw["nConL"]+1] for i in eachindex(Ovr_Sens)]
        LSCOPF= Linear_SENS_CP_new(LSCOPF,Vars_C,Parlist_C,nTP,raw,"Cont",SENSDict_C,Lininit,Cont_indices,acc,TR,Base_TR_C,inj_bound_C,Worst_Vio)
    end
    LSCOPF[:modTR]=TR
    return LSCOPF
end


function calculate_I_sol(e,f,raw)
    # e=LSCOPF[:e]
    # f=LSCOPF[:f]
    fb=raw["Line"]."From"
    tb=raw["Line"]."To"
    g=raw["Line"]."g (pu)"
    b=raw["Line"]."b (pu)"
    gb2= (g.^2+b.^2)
    ΔVline=((e[:,fb,:] .- e[:,tb,:]).^2) .+ ((f[:,fb,:] .- f[:,tb,:]).^2)
    I_sol=@. gb2' *ΔVline
    return I_sol
end

