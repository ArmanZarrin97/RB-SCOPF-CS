using PowerModels, OdsIO
using DataFrames

function Dict_to_DataFrame(DataDict::Dict{String,Any},fieldDict_)
    # DataDict,fieldDict_=data[conDict[i]],fieldDict[i]
    # DataDict,fieldDict_=data[i],[]
    df=DataFrame(id=sort(parse.(Int,keys(DataDict))))
    if !isempty(DataDict)
        if !isempty(fieldDict_)
            for key in eachindex(fieldDict_[1])
                df[!,fieldDict_[1][key]]=[DataDict[i][fieldDict_[2][key]] for i in string.(df."id")]
            end
            df=df[!,2:end]
        else
            for key in keys(DataDict[string.(df."id")[1]])
                df[!,key]=[DataDict[i][key] for i in string.(df."id")]
            end
        end
    end
    return df
end
function read_Powermodels(data)
    SysData= Dict{String,Any}()
    for i in keys(data)
        # i=collect(keys(data))[1]
        if !isempty(data[i])
            if data[i] isa Dict
                SysData[i]=Dict_to_DataFrame(data[i],[])
            else
                SysData[i]=data[i]
            end
        else
            SysData[i]=DataFrame()
        end
    end
    return SysData
end

function names_TESTIFY()
    Sheets=["Bus","Branch","Gen","RES","Load","Shunt","ConL","ConG","ESS","DSO","BaseMVA"]
    dataname=["bus","branch","gen","","load","shunt","","","storage","","baseMVA"]

    fields=[["Bus_N" "Bus_Type"	"Area" "V_Nom (kV)"	"V_max (pu)" "V_min (pu)" "V_mag (pu)" "V_angle (pu)"],
            ["From" "To" "g (pu)" "b (pu)"	"g_sh (pu)" "b_sh (pu)" "ratio" "rmin" "rmax" "RATE_A" "br_status" "Trans" "B_N"],
            ["Bus" "Pg (pu)" "Qg (pu)" "Q_max (pu)" "Q_min (pu)" "V_set (pu)" "P_max (pu)" "P_min (pu)" "cost" "G_N" "status" "Ramp"], #"slack" "Step" "Trip_priority" "gen_step"
            ["RES_N" "Bus" "Cap" "Tech" "PG" "QG"],
            ["Bus_N" "Pd (pu)" "Qd (pu)"],
            ["Bus" "bsh0" "gsh0"],
            ["cont" "From" "To" "index"],
            ["cont" "Bus"],
            ["Bus" "Ps" "Qs" "Energy (MWh)" "E_rating (MWh)" "Charge_rating (MW)" "Discharge_rating (MW)" "Charge_efficiency" "Discharge_efficiency" "Thermal_rating (MVA)" "Qmin (MVAr)" "Qmax (MVAr)" "R" "X" "P_loss" "Q_loss" "br_status" "soc_initial" "soc_min" "soc_max" "E_rating_min (MWh)" "E_initial (MWh)" "cost_a" "cost_b" "cost_c" "PG" "QG"],
            ["PG" "QG" "P_min" "P_max" "Q_min" "Q_max" "P_P" "P_Q"],
            ["base_MVA" "nTL" "Td" "Inc_coef" "Dec_coef" "Eresponse" "Mode" "conprob" "Ptype" "init" "LSGCpermit" "nTS" "nP" "Lscost" "Rccost"]]

    fieldsname=[["bus_i" "bus_type"	"zone" "base_kv" "vmax" "vmin" "vm" "va"],
                ["f_bus" "t_bus" "br_g" "br_b"	"g_sh" "b_sh" "tap" "angmin" "angmax" "rate_a" "br_status" "transformer" "index"],
                ["gen_bus" "pg" "qg" "qmax" "qmin" "vg" "pmax" "pmin" "cost" "id" "gen_status" "ramp_30"],
                ["res_id" "bus" "cap" "tech" "pg" "qg"],
                ["load_bus" "pd" "qd"],
                ["shunt_bus" "bs" "gs"],
                ["from" "to" "cont" "#"],
                ["cont" "bus"],
                ["bus" "ps" "qs" "energy" "e_rating" "charge_rating" "discharge_rating" "charge_efficiency" "discharge_efficiency" "thermal_rating" "qmin" "qmax" "r" "x" "p_loss" "q_loss" "br_status" "soc_initial" "soc_min" "soc_max" "e_rating_min" "e_initial" "cost_a" "cost_b" "cost_c" "PG" "QG"],
                ["pg" "qg" "pmin" "pmax" "qmin" "qmax" "p_p" "p_q"],
                ["base_mva" "ntl" "td" "inc_coef" "dec_coef" "eresponse" "mode" "conprob" "ptype" "init" "lsgcpermit" "nts" "np" "lscost" "rccost"]]
    fieldDict=Dict(zip(Sheets,zip(fields,fieldsname)))
    conDict=Dict(zip(Sheets,dataname))
    return fieldDict,conDict
end

function convert_TESTIFY(PM_data,conDict,fieldDict)
    # Line and Bus are OK currently

    PM_data["branch"]."br_g"=round.(real(1 ./(PM_data["branch"]."br_r"+PM_data["branch"]."br_x"*im)),digits=9)
    PM_data["branch"]."br_b"=round.(imag(1 ./(PM_data["branch"]."br_r"+PM_data["branch"]."br_x"*im)),digits=9)
    PM_data["branch"]."b_sh"=round.(PM_data["branch"]."b_fr" .+ PM_data["branch"]."b_to",digits=9)
    PM_data["branch"]."g_sh"=round.(PM_data["branch"]."g_fr" .+ PM_data["branch"]."g_to",digits=9)

    # PM_data["gen"]."cost"=round.(PM_data["branch"]."g_fr" .+ PM_data["branch"]."g_to",digits=9)

    SysData= Dict{String,Any}()
    Sheets=collect(keys(conDict))
    for i in Sheets
        # i=Sheets[9]
        if !isempty(conDict[i])
            DataDict=PM_data[conDict[i]]
            if DataDict isa DataFrame #&& !isempty(DataDict)
                fieldDict_=fieldDict[i]
                if isempty(DataDict) df=DataFrame(id=[]) else df=DataFrame(id=DataDict."id") end
                for key in eachindex(fieldDict_[1])
                    if isempty(fieldDict_[2][key]) || !hasproperty(DataDict,Symbol(fieldDict_[2][key]))
                        df[!,fieldDict_[1][key]]= []
                    else
                        df[!,fieldDict_[1][key]]= DataDict[!,fieldDict_[2][key]]
                    end
                end
                df=df[!,2:end]
                SysData[i]=df
            else
                df=DataFrame(id=1)
                df[!,fieldDict[i][1][1]]=[PM_data[conDict[i]]]
                df=df[!,2:end]
                SysData[i]=df
            end
        else
            df=DataFrame(id=[])
            fieldDict_=fieldDict[i]
            for key in eachindex(fieldDict_[1])
                    df[!,fieldDict_[1][key]]= []
            end
            df=df[!,2:end]
            SysData[i]=df
        end
    end
    return SysData
end


# file_name = "case60nordic.m"
# data=parse_file(file_name;import_all=true)
# PM_data=read_Powermodels(data)
# fieldDict,conDict=names_TESTIFY()
# SysData= convert_TESTIFY(PM_data,conDict,fieldDict)


