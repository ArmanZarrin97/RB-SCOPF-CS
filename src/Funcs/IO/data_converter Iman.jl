using PowerModels, OdsIO
using DataFrames

file_name = "case60nordic.m"
data=parse_file(file_name;import_all=true)
PM_data=read_Powermodels(data)
fieldDict,conDict=names_TESTIFY()
SysData= convert_TESTIFY(PM_data,conDict,fieldDict)
SysData=modifications!(SysData)

function modifications!(SysData)
    if file_name=="case60nordic.m"
        #others
        SysData["Gen"]."Pg (pu)"./=2
        SysData["Gen"]."P_max (pu)"./=2

        #contingencies
        SysData["ConL"]=DataFrame(
        cont=   [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33][1,:],
        From=   [1	2	3	25	5	5	6	6	7	10	26	26	26	26	23	23	27	28	31	14	14	15	15	16	16	16	36	37	37	17	20	21	12][1,:],
        To=     [3 24	24	4	7	9	8	9	8	11	23	27	28	12	28	12	29	31	29	36	20	16	36	36	17	18	37	19	21	18	21	22	13][1,:],
        index=  [1	3	5	7	9	11	13	15	16	18	20	21	22	23	24	25	26	29	31	40	41	42	43	44	45	46	47	49	51	52	53	54	56][1,:]
        )
        #storage
        ESSdata=
        [5	0	0	1	2200	400	400	0.95	0.95	1000	-50	70	0	0.01	0	0	1	5.00E-03	0	1	660	1100	0	4	0	0	0
        17	0	0	1	2200	400	400	0.95	0.95	1000	-50	70	0	0.01	0	0	1	5.00E-03	0	1	660	1100	0	4	0	0	0
        31	0	0	1	2200	400	400	0.95	0.95	1000	-50	70	0	0.01	0	0	1	5.00E-03	0	1	660	1100	0	4	0	0	0]
        propertynames(SysData["ESS"])
        append!(SysData["ESS"],DataFrame(ESSdata, propertynames(SysData["ESS"])))

        #Renewable
        RESdata=
        [1	38	7.2	1	0	0
        2	39	1.7	1	0	0
        3	40	6.3	1	0	0
        4	48	5.7	1	0	0
        5	51	6.3	1	0	0]
        propertynames(SysData["RES"])
        append!(SysData["RES"],DataFrame(RESdata, propertynames(SysData["RES"])))

        ##!!!### be careful about modelling ratio of lines

        #flexible Loads
        SysData["Load"]


    end
    return SysData
end

SysData=moussa_PF_modify!(SysData)

function moussa_PF_modify!(SysData)
    #generator PARAMETERS for PF
    coln=size(SysData["Gen"],2)
    SysData["Gen"][!,"slack"].=1
    SysData["Gen"][!,"Step"].=0.02
    SysData["Gen"][!,"Trip_priority"]=SysData["Gen"][!,"G_N"]
    SysData["Gen"][!,"gen_step"]=(SysData["Gen"][!,"P_max (pu)"]-SysData["Gen"][!,"P_min (pu)"])/5
    naddcol=size(SysData["Gen"],2)-coln
    SysData["Gen"]=SysData["Gen"][:,[(1:11)...,(end-naddcol+1:end)...,(12:coln)...]]
    return SysData
end


standardize_cost_terms!(data, order=2)

keys(data["gen"]["1"])

calc_thermal_limits!(data)

data["branch"]["1"]["transformer"]
# ref =build_ref(data)[:nw][0]

buses=[data["bus"][i]["bus_i"] for i in keys(data["bus"])]
type=[data["bus"][i]["bus_type"] for i in keys(data["bus"])]
area=[data["bus"][i]["zone"] for i in keys(data["bus"])]
v_nom=[data["bus"][i]["base_kv"] for i in keys(data["bus"])]
v_max=[data["bus"][i]["vmax"] for i in keys(data["bus"])]
v_min=[data["bus"][i]["vmin"] for i in keys(data["bus"])]
va=[data["bus"][i]["va"] for i in keys(data["bus"])]
vm=[data["bus"][i]["vm"] for i in keys(data["bus"])]
v_real=vm.*cos.(va)
v_imag=vm.*sin.(va)

fr_branch=[]
to_branch=[]
g_branch=[]
b_branch=[]
g_sh_branch=[]
b_sh_branch=[]
rate_a_branch=[]
br_stat=[]

fr_branch_tr=[]
to_branch_tr=[]
g_tr=[]
b_tr=[]
# g_sh_tr=[]
b_sh_tr=[]
ratio_tr=[]
shift_angle=[]
rmin_tr=[]
rmax_tr=[]
rate_a_tr=[]
br_stat_tr=[]

for i in keys(data["branch"])
    # if data["branch"][i]["br_status"]==1
    if !data["branch"][i]["transformer"]
        from_br=data["branch"][i]["f_bus"]
        to_br= data["branch"][i]["t_bus"]
        r=         data["branch"][i]["br_r"]
        x=         data["branch"][i]["br_x"]
        g=r./(r.^2+x.^2)
        b=-x./(r.^2+x.^2)
        g_sh=       2*data["branch"][i]["g_fr"]
        b_sh=       2*data["branch"][i]["b_fr"]
        rate_a=     data["branch"][i]["rate_a"]
        br_status=  data["branch"][i]["br_status"]
        push!(fr_branch,from_br)
        push!(to_branch,to_br)
        push!(g_branch,g)
        push!(b_branch,b)
        push!(g_sh_branch,g_sh)
        push!(b_sh_branch,b_sh)
        push!(rate_a_branch,rate_a)
        push!(br_stat,br_status)
    else
        from_br=data["branch"][i]["f_bus"]
        to_br= data["branch"][i]["t_bus"]
        r=         data["branch"][i]["br_r"]
        x=         data["branch"][i]["br_x"]
        g=r./(r.^2+x.^2)
        b=-x./(r.^2+x.^2)
        # g_sh=       2*data["branch"][i]["g_fr"]
        b_sh=       2*data["branch"][i]["b_fr"]
        ratio=data["branch"][i]["tap"]
        shift=data["branch"][i]["shift"]
        rmin=1
        rmax=1

        rate_a=     data["branch"][i]["rate_a"]
        br_status_tr=  data["branch"][i]["br_status"]

        push!(g_tr,g)
        push!(b_tr,b)
        # push!(g_sh_tr,g_sh)
        push!(b_sh_tr,b_sh)
        push!(ratio_tr,ratio)
        push!(shift_angle,shift)
        push!(rmin_tr,rmin)
        push!(rmax_tr,rmax)
        push!(rate_a_tr,rate_a)
        push!(br_stat_tr,br_status_tr)
        push!(fr_branch_tr,from_br)
        push!(to_branch_tr,to_br)
    end
# end
end
# tr=[fr_branch_tr to_branch_tr]
#
# tr_new=[]
# for i in 1:size(tr[:,1],1)
#     for j in 1:size(tr[:,1],1)
#         if tr[i,:]==tr[j,:] && i!=j
#             tr[j,:]=[0 0]
#         end
#
#     end
#     push!(tr_new,tr)
# end
# tr_final=[setdiff(tr_new[1][:,1], 0) setdiff(tr_new[1][:,2], 0)]
# fr_tr=tr_final[:,1]
# to_tr=tr_final[:,2]

load_bus_map=[data["load"][i]["load_bus"] for i in keys(data["load"])]
pd_load=     [data["load"][i]["pd"] for i in keys(data["load"])]
qd_load=     [data["load"][i]["qd"] for i in keys(data["load"])]
columns=[]
for i in keys(data["load"])
    push!(columns,1)
end

gen_bus_map=[data["gen"][i]["gen_bus"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
pg_gen=[data["gen"][i]["pg"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
qg_gen=[data["gen"][i]["qg"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
gen_stat=[data["gen"][i]["gen_status"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
pmin_gen=[data["gen"][i]["pmin"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
pmax_gen=[data["gen"][i]["pmax"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
qmin_gen=[data["gen"][i]["qmin"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
qmax_gen=[data["gen"][i]["qmax"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
vg_gen  =[data["gen"][i]["vg"] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]


# gen_bus_map=[data["gen"][i]["gen_bus"] for i in keys(data["gen"]) ]
# pg_gen=[data["gen"][i]["pg"] for i in keys(data["gen"]) ]
# qg_gen=[data["gen"][i]["qg"] for i in keys(data["gen"]) ]
# gen_stat=[data["gen"][i]["gen_status"] for i in keys(data["gen"]) ]
# pmin_gen=[data["gen"][i]["pmin"] for i in keys(data["gen"]) ]
# pmax_gen=[data["gen"][i]["pmax"] for i in keys(data["gen"]) ]
# qmin_gen=[data["gen"][i]["qmin"] for i in keys(data["gen"]) ]
# qmax_gen=[data["gen"][i]["qmax"] for i in keys(data["gen"]) ]
# vg_gen  =[data["gen"][i]["vg"] for i in keys(data["gen"]) ]
# %  CWS (Connection with Spain),
# %  FOG (Fossil Gas),
# %  FHC (Fossil Hard Coal),
# %  HWR (Hydro Water Reservoir),
# %  HPS (Hydro Pumped Storage),
# %  HRP (Hydro Run-of-river and poundage),
# %  SH1 (Small Hydro - P ≤ 10 MW),
# %  SH3 (Small Hydro - 10 MW < P ≤ 30 MW),
# %  PVP (Photovoltaic power plant),
# %  WON (Wind onshore),
# %  WOF (Wind offshore),
# %  MAR (Marine),
# %  OTH (Other thermal, such as geothermal, biomass, biogas, Municipal solid waste and CHP renewable and non-renewable)
# cost_b=[]
# for i in keys(data["gen"])
#     # if data["gen"][i]["gen_status"]==1
#     if data["gen"][i]["col_1"]=="CWS"
#         cost=1
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="FOG"
#         cost=2
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="FHC"
#         cost=3
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="HWR"
#         cost=4
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="HPS"
#         cost=5
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="HRP"
#         cost=6
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="SH1"
#         cost=7
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="SH3"
#         cost=8
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="PVP"
#         cost=9
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="WON"
#         cost=10
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="WOF"
#         cost=11
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="MAR"
#         cost=12
#         push!(cost_b,cost)
#     elseif data["gen"][i]["col_1"]=="OTH"
#         cost=13
#         push!(cost_b,cost)
#     end
# # end
# end
#




cost_a=[ data["gen"][i]["cost"][1]/(data["baseMVA"]^2) for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
cost_b=[ data["gen"][i]["cost"][2]/(data["baseMVA"]) for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]
cost_c=[ data["gen"][i]["cost"][3] for i in keys(data["gen"]) if data["gen"][i]["gen_status"]==1]

ods_write("input_data\\$nam.ods",Dict(("Buses",1,1)=>
["Bus_N" "Bus_Type"	"Area" "V_Nom (kV)"  "V_max (pu)" "V_min (pu)" "V_init_re (pu)" "V_init_im (pu)"]))

ods_write("input_data\\$nam.ods",Dict(("Lines",1,1)=>
["From" "To" "g (pu)" "b (pu)"  "g_sh (pu)" "b_sh (pu)" "RATE_A" "br_status"]))

ods_write("input_data\\$nam.ods",Dict(("Gens",1,1)=>
["Bus" "Pg (pu)" "Qg (pu)" "Q_max (pu)"  "Q_min (pu)" "V_set (pu)" "P_max (pu)" "P_min (pu)" "cost_a"   "cost_b"  "cost_c"
 ]))


ods_write("input_data\\$nam.ods",Dict(("Loads",1,1)=>
["Bus_N" "Pd (pu)" "Qd (pu)" ]))


ods_write("input_data\\$nam.ods",Dict(("Transformers",1,1)=>
["From" "To" "g" "b"  "bsh" "ratio" "rmin" "rmax"  "Snom" "BStatus0" ]))

ods_write("input_data\\$nam.ods",Dict(("shunts",1,1)=>["Bus"	"bsh0"	"bshmin"	"bshmax"]))
if !isempty(data["shunt"])
shnt_bus=[data["shunt"][i]["shunt_bus"] for i in keys(data["shunt"])]
shnt_val=[data["shunt"][i]["bs"] for i in keys(data["shunt"])]
shnt_min=[]
shnt_max=[]
for i in keys(data["shunt"])
    if data["shunt"][i]["bs"]>0
        push!(shnt_max, 2*data["shunt"][i]["bs"])
        push!(shnt_min, 0)
    elseif data["shunt"][i]["bs"]<0
        push!(shnt_max, 0)
        push!(shnt_min, 2*data["shunt"][i]["bs"])
    end
end

ods_write("input_data\\$nam.ods",Dict(("shunts",2,1)=>[shnt_bus shnt_val shnt_min shnt_max]))
end


ods_write("input_data\\$nam.ods",Dict(("Buses",2,1)=>[buses zeros(size(buses,1))]))
ods_write("input_data\\$nam.ods",Dict(("Buses",2,2)=>[type zeros(size(type,1))]))
ods_write("input_data\\$nam.ods",Dict(("Buses",2,3)=>[area zeros(size(area,1))]))
ods_write("input_data\\$nam.ods",Dict(("Buses",2,4)=>[v_nom zeros(size(v_nom,1))]))
ods_write("input_data\\$nam.ods",Dict(("Buses",2,5)=>[v_max zeros(size(v_max,1))]))
ods_write("input_data\\$nam.ods",Dict(("Buses",2,6)=>[v_min zeros(size(v_min,1))]))
ods_write("input_data\\$nam.ods",Dict(("Buses",2,7)=>[v_real v_imag]))
# ods_write("input_data\\$nam.ods",Dict(("Buses",2,8)=>[v_imag zeros(size(v_imag,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Buses",2,8)=>[v_imag nothing]))

# ods_write("input_data\\$nam.ods",Dict(("Buses",2,9)=>[nothing nothing]))

ods_write("input_data\\$nam.ods",Dict(("Lines",2,1)=>[fr_branch zeros(size(fr_branch,1))]))
ods_write("input_data\\$nam.ods",Dict(("Lines",2,2)=>[to_branch zeros(size(to_branch,1))]))
ods_write("input_data\\$nam.ods",Dict(("Lines",2,3)=>[g_branch zeros(size(g_branch,1))]))
ods_write("input_data\\$nam.ods",Dict(("Lines",2,4)=>[b_branch zeros(size(b_branch,1))]))
ods_write("input_data\\$nam.ods",Dict(("Lines",2,5)=>[g_sh_branch zeros(size(g_sh_branch,1))]))
ods_write("input_data\\$nam.ods",Dict(("Lines",2,6)=>[b_sh_branch zeros(size(b_sh_branch,1))]))
ods_write("input_data\\$nam.ods",Dict(("Lines",2,7)=>[rate_a_branch br_stat]))
# ods_write("input_data\\$nam.ods",Dict(("Lines",2,8)=>[br_stat zeros(size(br_stat,1))]))

# ods_write("input_data\\$nam.ods",Dict(("Transformers",2,1)=>)
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,1)=>[fr_branch_tr zeros(size(fr_branch_tr,1))]))
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,2)=>[to_branch_tr zeros(size(to_branch_tr,1))]))
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,3)=>[g_tr zeros(size(g_tr,1))]))
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,4)=>[b_tr zeros(size(b_tr,1))]))
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,5)=>[b_sh_tr zeros(size(b_sh_tr,1))]))
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,6)=>[ratio_tr zeros(size(ratio_tr,1))]))
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,7)=>[rmin_tr zeros(size(rmin_tr,1))]))
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,8)=>[rmax_tr zeros(size(rmax_tr,1))]))
ods_write("input_data\\$nam.ods",Dict(("Transformers",2,9)=>[rate_a_tr br_stat_tr]))
# ods_write("input_data\\$nam.ods",Dict(("Transformers",2,10)=>[br_stat_tr zeros(size(br_stat_tr,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Transformers",2,11)=>[shift_angle zeros(size(shift_angle,1))]))

ods_write("input_data\\$nam.ods",Dict(("Loads",2,1)=>[load_bus_map pd_load  qd_load]))
# ods_write("input_data\\$nam.ods",Dict(("Loads",2,2)=>[]))
# ods_write("input_data\\$nam.ods",Dict(("Loads",2,3)=>[qd_load  zeros(size(qd_load,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Loads",2,4)=>[columns  zeros(size(columns,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Loads",2,5)=>[columns  zeros(size(columns,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Loads",2,6)=>[columns  zeros(size(columns,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Loads",2,7)=>[columns  zeros(size(columns,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Loads",2,8)=>[columns  zeros(size(columns,1))]))


ods_write("input_data\\$nam.ods",Dict(("Gens",2,1)=>[gen_bus_map  zeros(size(gen_bus_map,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,2)=>[pg_gen  zeros(size(pg_gen,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,3)=>[qg_gen  zeros(size(qg_gen,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,4)=>[qmax_gen  zeros(size(qmax_gen,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,5)=>[qmin_gen  zeros(size(qmin_gen,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,6)=>[vg_gen  zeros(size(vg_gen,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Gens",2,7)=>[gen_stat  zeros(size(gen_stat,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,7)=>[pmax_gen  zeros(size(pmax_gen,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,8)=>[pmin_gen  zeros(size(pmin_gen,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,9)=>[cost_a  zeros(size(cost_a,1))]))
ods_write("input_data\\$nam.ods",Dict(("Gens",2,10)=>[cost_b  cost_c]))
# ods_write("input_data\\$nam.ods",Dict(("Gens",2,11)=>[cost_c  zeros(size(cost_c,1))]))
ods_write("input_data\\$nam.ods",Dict(("Base_MVA",1,1)=>["base_MVA" 0]))
ods_write("input_data\\$nam.ods",Dict(("Base_MVA",2,1)=>[data["baseMVA"] 0]))



#
# ods_write("input_data\\$(nam)_PROF.ods",Dict(("Dimension",1,1)=>["base_MVA" "NTP"	"NSC"	"RES_MAG"	"Load_MAG"	"Line_MAG"	"PF"	"load_correction"]))
#
# ods_write("input_data\\$(nam)_PROF.ods",Dict(("P_Profiles_Load",1,1)=>
# # ["Bus"	"t1"	"t2"	"t3"	"t4"	"t5"	"t6"	"t7"	"t8"	"t9"	"t10"	"t11"	"t12"	"t13"	"t14"	"t15"	"t16"	"t17"	"t18"	"t19"	"t20"	"t21"	"t22"	"t23"	"t24"]))
# ["Bus"	"t1"	]))
#  # ods_write("input_data\\$nam_PROF.ods",Dict(("Base_MVA",1,1)=>["base_MVA" 0]))
# # ods_write("input_data\\$nam.ods",Dict(("Base_MVA",2,1)=>[data["baseMVA"]  zeros(size(data["baseMVA"],1))]))
#
# ods_write("input_data\\$nam.ods",Dict(("P_Profiles_Load",1,1)=>
# # ["Bus"	"t1"	"t2"	"t3"	"t4"	"t5"	"t6"	"t7"	"t8"	"t9"	"t10"	"t11"	"t12"	"t13"	"t14"	"t15"	"t16"	"t17"	"t18"	"t19"	"t20"	"t21"	"t22"	"t23"	"t24"]))
# ["Bus"	"t1"	]))
#
# ods_write("input_data\\$nam.ods",Dict(("Q_Profiles_Load",1,1)=>
# # ["Bus"	"t1"	"t2"	"t3"	"t4"	"t5"	"t6"	"t7"	"t8"	"t9"	"t10"	"t11"	"t12"	"t13"	"t14"	"t15"	"t16"	"t17"	"t18"	"t19"	"t20"	"t21"	"t22"	"t23"	"t24"]))
# ["Bus"	"t1"	]))
#
# ods_write("input_data\\$nam.ods",Dict(("contingencies",1,1)=>["cont"	"From"	"To"]))
#
# ods_write("input_data\\$nam.ods",Dict(("Storage",1,1)=>["Bus"	"Ps"	"Qs"	"Energy (MWh)"	"E_rating (MWh)"	"Charge_rating (MW)"	"Discharge_rating (MW)"	"Charge_efficiency"	"Discharge_efficiency"	"Thermal_rating (MVA)"	"Qmin (MVAr)"	"Qmax (MVAr)"	"R"	"X"	"P_loss"	"Q_loss"	"Status"	"soc_initial"	"soc_min"	"soc_max"	"E_rating_min (MWh)"	"E_initial (MWh)"	"cost_a"	"cost_b"	"cost_c"
# ]))
#
# ods_write("input_data\\$nam.ods",Dict(("P_Profiles_Load",2,1)=>[load_bus_map zeros(size(load_bus_map,1))]))
# ods_write("input_data\\$nam.ods",Dict(("P_Profiles_Load",2,2)=>[pd_load  zeros(size(pd_load,1))]))
#
# ods_write("input_data\\$nam.ods",Dict(("Q_Profiles_Load",2,1)=>[load_bus_map zeros(size(load_bus_map,1))]))
# ods_write("input_data\\$nam.ods",Dict(("Q_Profiles_Load",2,2)=>[qd_load  zeros(size(qd_load,1))]))
