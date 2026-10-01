function SEVF(LNS, GT, RC, BT)
    return LNS+BT+1*(GT+RC)
end

function SEVF_casreps(casrep;stage=nothing)
    if isnothing(stage)
        LNS, GT, RC, BT = casrep[end, "LNS_T"], casrep[end, "GT"], casrep[end, "RC"], casrep[end, "BT"]
    else
        LNS, GT, RC, BT = casrep[stage, "LNS_T"], casrep[stage, "GT"], casrep[stage, "RC"], casrep[stage, "BT"]
    end
    return SEVF(LNS, GT, RC, BT)
end

function probdist(VMs; min=Glob["I_relax"], max=Glob["I_relaxmax"])
    mVM=maximum(VMs,init=0)
    if mVM<min
        return 0
    elseif mVM>max
        return 1
    else
        return (mVM^2-min^2)/(max^2-min^2)
    end
end

function probrisk(mVM; min=Glob["I_relax"], max=Glob["I_relaxmax"])
    if mVM<min
        return 0
    elseif mVM>max
        return 1
    else
        return (mVM^2-min^2)/(max^2-min^2)
    end
end

function Esev_cal_prob(casrep)
    VMs=casrep[1, "VM"]
    p=probrisk.(VMs)
    r= sum(p;init=0)
    return r
end

function Esev_cal(casrep)
    p=casrep[1, "prob"]
    r= p*casrep[end,"Sev"]+ (1-p) *casrep[1,"Sev"]
    return r
end

function Esev_cal_exact(casrep)
    si=size(casrep,1)
    prob=1 .-casrep[:, "prob"]
    for i in 1:si
        for j=1:i-1
            prob[i]*= casrep[j,"prob"]
        end
    end
    r= sum(prob.*casrep[:,"Sev"])
    return r
end



function Rcas_cal(casreps)
    p_k= [Glob["conprob_"]/size(casreps,2) for i in 1:size(casreps,2)]
    α_t= [1-t*0.01 for t in 1:size(casreps,1)]
    rvalues=Esev_cal.(casreps)
    rvalues_exact=Esev_cal_exact.(casreps)
    rvalues_prob=Esev_cal_prob.(casreps)
    expected_r= p_k'.*rvalues
    expected_r_exact= p_k'.*rvalues_exact
    Rcas_simple=sum(expected_r.*α_t)
    Rcas_exact= sum(expected_r_exact.*α_t)
    Rcas_prob= sum(rvalues_prob.*α_t)
    R_OP_max_simple= maximum(rvalues)
    R_OP_max_exact= maximum(rvalues_exact)
    R_OP_max_prob= maximum(rvalues_prob)
    return [Rcas_prob, Rcas_exact, Rcas_simple], [R_OP_max_prob, R_OP_max_exact, R_OP_max_simple]
end


function Rcas_cal_RS(RS)
    p_k= [Glob["conprob_"]/size(RS,2) for i in 1:size(RS,2)]
    α_t= [1-t*0.01 for t in 1:size(RS,1)]
    expected_r= p_k'.*RS
    Rcas=sum(expected_r.*α_t)
    return Rcas
end
# function Risk_OP_cal(casrep)

# end

