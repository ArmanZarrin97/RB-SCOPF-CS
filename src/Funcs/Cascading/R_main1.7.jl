using JuMP, OdsIO, MathOptInterface, Dates, LinearAlgebra, JLD2, XLSX, BenchmarkTools, GraphRecipes, Polynomials
using StatsBase, SparseArrays, Profile, Random, Distributions, Base.Threads


println()
printstyled("Hiiiii, This is Power System Cascading Outage Simulation Tool (PSCOST)"; color=:13, blink=true, bold=false)
println()
print("NewEngland39")
nam = "NewEngland39"
printstyled("
><><>< ----------------------------------------------------------------------------------------><><><"; color=:16, blink=true)

# B, L, T, G, Ld, Shu,REN,BAT, Lout=Matrix{Any}(NSysData["Bus"]), Matrix{Any}(NSysData["Line"]), Matrix{Any}(NSysData["Trans"]), Matrix{Any}(NSysData["Gen"]), Matrix{Any}(NSysData["Load"]), Matrix{Any}(NSysData["Shunt"]),Matrix{Any}(NSysData["RES"]),Matrix{Any}(NSysData["ESS"]),[]
B, L, T, G, Ld, Shu,REN,BAT, Lout=NSysData["Bus"], NSysData["Line"], NSysData["Trans"], NSysData["Gen"], NSysData["Load"], NSysData["Shunt"],NSysData["RES"],NSysData["ESS"],[]

# OP_Valid["V"][1, Ld[:, 1]]

df_min = 0.02;df_max = 0.01;df_eps = 0.001;df_minmin = 0.06;df_maxmax = 0.06;rcf_max = 0.021;rcf_min = 0.021;eps_g = 0.01;eps_v = 0.01;
epsilon_pq = 0.05;Vth = 0.9;LdI = zeros(size(Ld, 1));Ldsh_f = zeros(size(Ld, 1));Ldsh_v = zeros(size(Ld, 1));LdI_v = zeros(size(Ld, 1));
Gt = zeros(size(G, 1));nBus = size(B, 1);sb = 100;epsilon_pq = 0.05;

# # ----------------------------------------------------------------------------------------
# 3.--------------------Apply contingencies (initiating events/outages) --------------------
# # ----------------------------------------------------------------------------------------
Lout = [23]

#--------------------------------------------------------------------------------------------------
# # -----------------------------------------------------------------------------------------------
# # --------------------The proposed algorithm for cascading --------------------------------------
# # ----------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

include("power_flow_IND.jl")
include("sub_function.jl")
include("middle_loop.jl")
include("island.jl")
include("CSC.jl")

# ires = findall3(x -> x != 0, L[:, 8])
# L=L[ires, :]
# L[:, 9] .= collect(1:size(L,1))  

# ires = findall3(x -> x != 0, T[:, 10])
# T=T[ires, :]
# T[:, 11] .= collect(1:size(T,1))

Bs = deepcopy(B)
Ls = deepcopy(L)
Ts = deepcopy(T)
Gs = deepcopy(G)
Lds = deepcopy(Ld)
Shus = deepcopy(Shu)

(Lds, Gs) = initiate_IND(Bs, Ls, Ts, Gs, Lds, Shus); # Please deactive this, and perform the initialization manually instead
# OP_Valid["V"][1, Ld[:, 1]]

@time (Rn0, H0, u0, B0, L0, T0, G0, Ld0, Shu0, new_L, LNS_v, LNS_f, LNS, tot_L, Line_vio, LNS_stage, JCB) = CSC(Bs, Ls, Ts, Gs, Lds, Shus, Lout);

JCB[1][2,:]
JCB[2][2,:]

# gggg=bbbbbbbbbb
# @time (report, LNS_v, LNS_f, LNS, tot_L, LNS_stage, JCB) = run_CSC(Bs, Ls, Ts, Gs, Lds, Shus, Lout)


