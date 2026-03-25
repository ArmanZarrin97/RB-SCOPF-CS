include("power_flow_IND.jl")
include("sub_function.jl")
include("middle_loop.jl")
include("island.jl")
include("CSC.jl")

    df_min = 0.02;df_max = 0.01;df_eps = 0.001;df_minmin = 0.06;df_maxmax = 0.06;rcf_max = 0.021;rcf_min = 0.021;eps_g = 0.01;eps_v = 0.01;
    epsilon_pq = 0.05;Vth = 0.9;sb = 100;epsilon_pq = 0.05;

    # Vth = 0.92
mutable struct Globcascade
    LdI::Vector{Float64}
    LdI_v::Vector{Float64}
    Ldsh_v::Vector{Float64}
    Ldsh_f::Vector{Float64}
    Gt::Vector{Float64}
    nBus::Int64
    print::Bool
    Opt_sens::DataFrame
end
function cascading_sim(Bs::DataFrame, Ls::DataFrame, Ts::DataFrame, Gs::DataFrame, Lds::DataFrame, Shus::DataFrame,RENs::DataFrame,BATs::DataFrame, Lout::Vector{Int64}; printcascade=false)
    # B, L, T, G, Ld, Shu,REN,BAT, Lout=CSysData["Bus"], CSysData["Line"], CSysData["Trans"], CSysData["Gen"], CSysData["Load"], CSysData["Shunt"],CSysData["RES"],CSysData["ESS"],Int64[]
    # Bs = deepcopy(B)
    # Ls = deepcopy(L)
    # Ts = deepcopy(T)
    # Gs = deepcopy(G)
    # Lds = deepcopy(Ld)
    # Shus = deepcopy(Shu)
    # RENs = deepcopy(REN)
    # BATs = deepcopy(BAT)

    # (Lds, Gs) = initiate_IND(Bs, Ls, Ts, Gs, Lds, Shus); # Please deactive this, and perform the initialization manually instead
    # # OP_Valid["V"][1, Ld[:, 1]
    # Globcascade.OptSens=DataFrame()
    LdI = zeros(size(Lds, 1));LdI_v = zeros(size(Lds, 1));
    Gt = zeros(size(Gs, 1));nBus = size(Bs, 1)
    Ldsh_v = zeros(size(Lds, 1))
    Ldsh_f = zeros(size(Lds, 1))
    Gcas= Globcascade(LdI,LdI_v,Ldsh_v,Ldsh_f,Gt,nBus,printcascade,DataFrame())
    (Rn0, H0, u0, B0, L0, T0, G0, Ld0, Shu0, new_L, LNS_v, LNS_f, LNS, tot_L, Line_vio, LNS_stage, JCB, casrep)=CSC(Bs, Ls, Ts, Gs, Lds, Shus, RENs, BATs, Lout,Gcas);
    # V, θ, p_inj, q_inj, convergence_indicator, Iij, percent, PGGG, QGGG, Jacob 



    # Rn0[1][1]["vmag_i"], Rn0[1][1]["vang_i"],Rn0[1][1]["p_i"],Rn0[1][1]["q_i"],Rn0[1][1]["cg_i"],Rn0[1][1]["Iij"],Rn0[1][1]["percentIij"],G0[:,"Pg (pu)"],G0[:,"Qg (pu)"],JCB[1]

    # println(collect(keys(Rn0[1][1])))

    # return (Rn0, H0, u0, B0, L0, T0, G0, Ld0, Shu0, new_L, LNS_v, LNS_f, LNS, tot_L, Line_vio, LNS_stage, JCB)
    return Rn0[1][1]["vmag_i"], Rn0[1][1]["vang_i"],Rn0[1][1]["p_i"],Rn0[1][1]["q_i"],Rn0[1][1]["cg_i"],Rn0[1][1]["Iij"],Rn0[1][1]["percentIij"],G0[:,"Pg (pu)"],G0[:,"Qg (pu)"],JCB[1],LNS_stage,casrep

end

# function cascading_sim(Bs, Ls, Ts, Gs, Lds, Shus,REN,BAT, Lout)
#    @time CSC(Bs, Ls, Ts, Gs, Lds, Shus, Lout);
# end

# @time cascading_sim(NSysData["Bus"], NSysData["Line"], NSysData["Trans"], NSysData["Gen"], NSysData["Load"], NSysData["Shunt"],NSysData["RES"],NSysData["ESS"],[23]);
