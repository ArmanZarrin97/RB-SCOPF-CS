using JuMP, OdsIO, MathOptInterface, Dates, LinearAlgebra, JLD2, XLSX, BenchmarkTools, GraphRecipes, Polynomials
using StatsBase, SparseArrays, Profile, Random, Distributions, Base.Threads


# 1. -------------------Accessing current folder directory--------------------------
println(pwd());
cd(dirname(@__FILE__))
print(pwd());
printstyled("
><><>< ----------------------------------------------------------------------------------------><><><"; color=:16, blink=true)
println()
printstyled("Hiiiii, This is Power System Cascading Outage Simulation Tool (PSCOST)"; color=:13, blink=true, bold=false)
println()


print("NewEngland39")
nam = "NewEngland39"

# println()
# printstyled("PSCOST is triggered by outage of line (s) (please specify line (s) number): "; color=:12, blink=true)
# Lout = readline()
printstyled("
><><>< ----------------------------------------------------------------------------------------><><><"; color=:16, blink=true)
global filename = "input_data/$(nam).ods"
# # ----------------------------------------------------------------------------------------
# 2.--------------------Read Data from the Excel file---------------------------------------
# # ----------------------------------------------------------------------------------------
# sheetname = "Buses";
# fieldsB = ["Bus_N", "Bus_Type", "Area", "V_Nom (kV)", "V_max (pu)", "V_min (pu)", "V_init_re (pu)", "V_init_im (pu)"]; # Fields that have to be read from the file
# raw_dataB = ods_readall(filename; sheetsNames=[sheetname], innerType="Matrix")
# raw_dataB = raw_dataB[sheetname]   # Conversion from Dict to Array
# B = raw_dataB[2:end, 1:end]


# sheetname = "Lines";
# fieldsL = ["From", "To", "g (pu)", "b (pu)", "g_sh (pu)", "b_sh (pu)", "RATE_A", "br_status"]; # Fields that have to be read from the file
# raw_dataL = ods_readall(filename; sheetsNames=[sheetname], innerType="Matrix")
# raw_dataL = raw_dataL[sheetname]   # Conversion from Dict to Array
# L = raw_dataL[2:end, 1:end]


# sheetname = "Transformers";
# fieldsT = ["From", "To", "g", "b", "bsh", "ratio", "rmin", "rmax", "Snom", "BStatus0"];
# raw_dataT = ods_readall(filename; sheetsNames=[sheetname], innerType="Matrix")
# raw_dataT = raw_dataT[sheetname]   # Conversion from Dict to Array
# T = raw_dataT[2:end, 1:end]


# sheetname = "Gens";
# fieldsG = ["Bus", "Pg (pu)", "Qg (pu)", "Q_max (pu)", "Q_min (pu)", "V_set (pu)", "P_max (pu)", "P_min (pu)", "cost_a", "cost_b", "cost_c"];
# raw_dataG = ods_readall(filename; sheetsNames=[sheetname], innerType="Matrix")
# raw_dataG = raw_dataG[sheetname]   # Conversion from Dict to Array
# G = raw_dataG[2:end, 1:end]

# sheetname = "Loads";
# raw_dataLd = ods_readall(filename; sheetsNames=[sheetname], innerType="Matrix")
# raw_dataLd = raw_dataLd[sheetname]   # Conversion from Dict to Array
# Ld = raw_dataLd[2:end, 1:end]




# sheetname = "shunts";
# fieldsSh = ["Bus", "bsh0", "bshmin", "bshmax"]
# raw_dataSh = ods_readall(filename; sheetsNames=[sheetname], innerType="Matrix")
# raw_dataSh = raw_dataSh[sheetname]   # Conversion from Dict to Array
# Shu = raw_dataSh[2:end, 1:end]

B, L, T, G, Ld, Shu,REN,BAT, Lout=Matrix{Any}(NSysData["Bus"]), Matrix{Any}(NSysData["Line"]), Matrix{Any}(NSysData["Trans"]), Matrix{Any}(NSysData["Gen"]), Matrix{Any}(NSysData["Load"]), Matrix{Any}(NSysData["Shunt"]),Matrix{Any}(NSysData["RES"]),Matrix{Any}(NSysData["ESS"]),[]


ires = findall3(x -> x != 0, L[:, 8])
L=L[ires, :]
L[:, 9] .= collect(1:size(L,1))  

ires = findall3(x -> x != 0, T[:, 10])
T=T[ires, :]
T[:, 11] .= collect(1:size(T,1))





df_min = 0.02
df_max = 0.01
df_eps = 0.001

df_minmin = 0.06
df_maxmax = 0.06
rcf_max = 0.021
rcf_min = 0.021
eps_g = 0.01
eps_v = 0.01
epsilon_pq = 0.05
Vth = 0.9
LdI = zeros(size(Ld, 1))
Ldsh_f = zeros(size(Ld, 1))
Ldsh_v = zeros(size(Ld, 1))
LdI_v = zeros(size(Ld, 1))
Gt = zeros(size(G, 1))
nBus = size(B, 1)
sb = 100
epsilon_pq = 0.05


# # ----------------------------------------------------------------------------------------
# 3.--------------------Apply contingencies (initiating events/outages) --------------------
# # ----------------------------------------------------------------------------------------
Lout = [8]

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

  
ires = findall3(x -> x != 0, L[:, 8])
L=L[ires, :]
L[:, 9] .= collect(1:size(L,1))  

ires = findall3(x -> x != 0, T[:, 10])
T=T[ires, :]
T[:, 11] .= collect(1:size(T,1))

Bs = deepcopy(B)
Ls = deepcopy(L)
Ts = deepcopy(T)
Gs = deepcopy(G)
Lds = deepcopy(Ld)
Shus = deepcopy(Shu)

(Lds, Gs) = initiate_IND(Bs, Ls, Ts, Gs, Lds, Shus) # Please deactive this, and perform the initialization manually instead

@time (Rn0, H0, u0, B0, L0, T0, G0, Ld0, Shu0, new_L, LNS_v, LNS_f, LNS, tot_L, Line_vio, LNS_stage, JCB) = CSC(Bs, Ls, Ts, Gs, Lds, Shus, Lout);


# gggg=bbbbbbbbbb
# @time (report, LNS_v, LNS_f, LNS, tot_L, LNS_stage, JCB) = run_CSC(Bs, Ls, Ts, Gs, Lds, Shus, Lout)


