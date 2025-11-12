# SYSTEM DATA
using OdsIO,DataFrames

#####################################################################
####################### Data Readers Function #######################
function Set_global_parameters(Test_Sys)

    global Path= "Test Systems/"*Test_Sys*"Data/"
    global SysData=SysData_Reader(Path)
    # ProfileCanvas.@profview SysData=SysData_Reader(Path)
    global Td=SysData["Base_MVA"]."Td"[1] ## Time duration in Day-ahead
    global Inc_coef=SysData["Base_MVA"]."Inc_coef"[1] ## FL Increase coefficient
    global Dec_coef=SysData["Base_MVA"]."Dec_coef"[1] ## FL Decrease coefficient
    global conprob_=SysData["Base_MVA"]."conprob"[1] ## each Contingency probability
    global SBase=SysData["Base_MVA"]."base_MVA"[1] ## Base MVA
    global Mode=SysData["Base_MVA"]."Mode"[1]   ## Mode of Secuirty: Preventive or Corrective
    global LSpermit=SysData["Base_MVA"]."LSGCpermit"[1] ## Permit of Load Shedding in Contingency cases
    global nTS=SysData["Base_MVA"]."nTS"[1] ## Number of Time periods in one hour in Real-time
    global nP=SysData["Base_MVA"]."nP"[1] ## size of look-ahead Horizon in Real-time in Hours
    global CCT=1/nTS; ## Time duration in Real-time in hours
    global ϵ=1e-3
    global v_relax_factor_Line_min = 0;
    global v_relax_factor_Line_max = 0;
        #################
            #modifications
            if Test_Sys=="Polish/"
                gout=findall(x-> (x."Q_max (pu)"-x."Q_min (pu)")<=0.01 ,eachrow(SysData["Gen"]))
                for i in eachrow(SysData["Gen"])
                    if i."P_min (pu)"!= i."P_max (pu)"
                        i."P_min (pu)"*=0.9
                        i."P_max (pu)"*=1.1
                    end
                    # if i."Q_max (pu)">=0
                    #     i."Q_max (pu)"*=1.5
                    # end
                    # if i."Q_min (pu)"<=0
                    #     i."Q_min (pu)"*=1.5
                    # end
                end
                SysData["Bus"][SysData["Gen"][gout,:Bus],:Bus_Type].=1
                SysData["Gen"]=SysData["Gen"][Not(gout),:]
                SysData["Gen"].G_N .= 1: size(SysData["Gen"],1)
                # SysData["Gen"].Ramp.*=3
                # SysData["ConL"]=SysData["ConL"][1:2,:]
            end


            if Test_Sys=="2737/"
                # SysData["Gen"].cost_b ./=100
                actline=findall(x->x==1,SysData["Line"]."br_status")
                SysData["Line"]=SysData["Line"][actline,:]
                SysData["Line"].Line_N .= 1: size(SysData["Line"],1)
                SysData["ConL"]."#"=[findall(x-> x.From==x1.From && x.To==x1.To,eachrow(SysData["Line"]))[1] for x1 in eachrow(SysData["ConL"])]
                lineindex=findall(x-> x>1.0, SysData["Line"].RATE_A)
                SysData["Line"].RATE_A .*=10
                # SysData["Line"]."g (pu)" .*=10
                # SysData["Line"]."b (pu)" .*=10
                SysData["Line"].RATE_A[lineindex] .*= (0.6/10)
                SysData["Line"]."g (pu)"[lineindex] .*= (0.6)
                SysData["Line"]."b (pu)"[lineindex] .*= (0.6)

                # admittance= @. complex(SysData["Line"]."g (pu)", SysData["Line"]."b (pu)")
                # Z= 1 ./ admittance
                # threshold = 1e-4
                # Z .= complex.(max.(real(Z), threshold), max.(imag(Z), threshold))
                # newadmittance= 1 ./Z
                # SysData["Line"]."g (pu)"=real(newadmittance)
                # SysData["Line"]."b (pu)"=imag(newadmittance)

                #
                # gout=findall(x-> x."Q_max (pu)"==x."Q_min (pu)" && x."P_max (pu)" == x."P_min (pu)",eachrow(SysData["Gen"]))
                gout=findall(x-> x."Q_max (pu)"==x."Q_min (pu)" ,eachrow(SysData["Gen"]))
                SysData["Bus"][SysData["Gen"][gout,:Bus],:Bus_Type].=1
                SysData["Gen"]=SysData["Gen"][Not(gout),:]
                SysData["Gen"].G_N .= 1: size(SysData["Gen"],1)
                for i in eachrow(SysData["Gen"])
                    if i."P_min (pu)"!= i."P_max (pu)"
                        i."P_min (pu)"/=2
                        i."P_max (pu)"*=2
                    end
                end
                SysData["Gen"].gen_step = (SysData["Gen"]."P_max (pu)"-SysData["Gen"]."P_min (pu)")/5
            end
        #################
    ################################### Profile Reader #################################
    True_Profile = True_Profile_Reader() # Read Real profiles
    ## new Profile generation ##
    Profile_Generator(new_profile,True_Profile)
    ## Read Profiles ##
    global Profile_Data= Profile_Reader(Path)
    ######################----------------------------------------#######################
    ######################----------------Initialization----------#######################
    global Init_Data=InitData_Reader(Path,["Zero"])
    ######################----------------------------------------#######################
    global conprob_=0.05;
    global TRfree=0.0
    global TRfreeC=0.0
    global relaxpermit=false
    global InitialTRB=1e-4
    global InitialTR=3.0
    # global Adaptive_TR=[0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 3.0]
    global Adaptive_TRB=[0.0001,0.1,0.5,1]
    # global Adaptive_TR=[3.0, 1.0, 0.5, 0.1, 0.05, 0.01, 0.005] #for 5bus
    global Adaptive_TR=[[1.0,0.001], [0.5,0.0005], [0.25,0.0002], [0.1,0.0001], [0.05,0.00002], [0.025,0.00001],[0.01,0.00001]] #for 5bus
    global filter_factor=[0.7,0.7,0.8,0.9,1,1,1]

    # # SysData["Load"]."Pd (pu)" .*= 1.2
    # # SysData["Load"]."Qd (pu)" .*= 1.2
    # # SysData["Line"].RATE_A.*=0.65
    # # SysData["Line"]."g (pu)" .*= 0.65
    # # SysData["Line"]."b (pu)" .*=0.65
    # # SysData["Line"]."b_sh (pu)" .*=0.65
    # # SysData["Line"].RATE_A[771:772]./=0.3
    # # SysData["Line"].RATE_A[1325] /= 0.3
    # # SysData["Line"].RATE_A[655] /=0.3
    # # SysData["Line"].RATE_A[922] /=0.3
    # # SysData["Line"].RATE_A[923] /=0.3
    # # SysData["Line"].RATE_A[2606] /=0.3
    # # SysData["Line"].RATE_A[754] /=0.3
    # # Mode="Preventive"
end
function SysData_Reader(Path::String)
    data=["Gen","RES", "Load", "Line", "Bus", "ConL", "ConG", "Trans", "Shunt","ESS","Base_MVA"]
    SysData=Dict{String,DataFrame}()
    for x in data
        try SysData[x]= ods_read(Path*"SysData",sheetName=x, retType= "DataFrame")
        catch
            error("SysData file does not exist.")
        end
    end
    SysData["FL"]=filter(x->x."ShFrOK"==1,SysData["Load"])[:,["Bus_N","D_N","ShFr","ShFrMin","ShFrMax","ShFrOK","ShCost","PF"]]
    return SysData
end

function InitData_Reader(Path,Stage)
    if "Zero" in Stage  Init_Data=init_zero(nTL) end
    initdata=["QG","PG","Pinc","Pdec","UE","Pcurt","Pcurtc","Pshed","Pch","Pdis","soc","Tcost","Vmag","Vang"]
    Init_Data=Dict()
    for y in Stage
        Init_Data[y]=Dict()
        if y=="Zero"  Init_Data[y]=init_zero(nTL)
        else
            for x in initdata
                try Init_Data[y][x]= ods_read(Path*y,sheetName=x, retType= "DataFrame")
                catch
                    error("Error in reading initialization file. solution: start with doing DA_OPF first. check if the init file has all fields") 
                end
            end
        end
    end
    return Init_Data
end

function init_zero(nTL)

    name=Path*"Zero"
    outlist=["PG", "QG","Pinc","Pdec","UE","Pcurt","Pcurtc","Pshed","Pch","Pdis","soc","Tcost","Vmag","Vang"];
    Out= Dict(
        "PG"=> DataFrame( [["t$i" for i in 1:nTL] zeros(nTL,size(SysData["Gen"],1))],Symbol.(["t";1:size(SysData["Gen"],1)])),
        "QG"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(SysData["Gen"],1))],Symbol.(["t";1:size(SysData["Gen"],1)])),
        "Pinc"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(filter(x->x."ShFrOK"==1,SysData["Load"]),1))],Symbol.(["t";1:size(filter(x->x."ShFrOK"==1,SysData["Load"]),1)])),
        "Pdec"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(filter(x->x."ShFrOK"==1,SysData["Load"]),1))],Symbol.(["t";1:size(filter(x->x."ShFrOK"==1,SysData["Load"]),1)])),
        "Pshed"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(SysData["Load"],1))],Symbol.(["t";1:size(SysData["Load"],1)])),
        "Pcurt"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(SysData["RES"],1))],Symbol.(["t";1:size(SysData["RES"],1)])),
        "Pcurtc"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(SysData["Gen"],1))],Symbol.(["t";1:size(SysData["Gen"],1)])),
        "Pch"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(filter(x->x."br_status"==1,SysData["ESS"]),1))],Symbol.(["t";1:size(filter(x->x."br_status"==1,SysData["ESS"]),1)])),
        "Pdis"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(filter(x->x."br_status"==1,SysData["ESS"]),1))],Symbol.(["t";1:size(filter(x->x."br_status"==1,SysData["ESS"]),1)])),
        "soc"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(filter(x->x."br_status"==1,SysData["ESS"]),1))],Symbol.(["t";1:size(filter(x->x."br_status"==1,SysData["ESS"]),1)])),
        "UE"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(filter(x->x."ShFrOK"==1,SysData["Load"]),1))],Symbol.(["t";1:size(filter(x->x."ShFrOK"==1,SysData["Load"]),1)])),
        "Tcost"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,1)],Symbol.(["t";1:1])),
        "Vmag"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(SysData["Bus"],1))],Symbol.(["t";1:size(SysData["Bus"],1)])),
        "Vang"=> DataFrame([["t$i" for i in 1:nTL] zeros(nTL,size(SysData["Bus"],1))],Symbol.(["t";1:size(SysData["Bus"],1)]))
        )
    ## Set initial SoC of battery. It can be changed from Excel file(SysData)
    if size(filter(x->x."br_status"==1,SysData["ESS"]),1)!=0 Out["soc"][:,2:end].= repeat(SysData["ESS"][:,"E_initial (MWh)"]'./SBase, nTL) end
    return Out
    # for i in outlist
    #     ods_write(name,Dict((i,1,1)=> DataFrame((t=["t$i" for i in 1:nTL]))))
    #     ods_write(name,Dict((i,1,2)=>Out[!,i][1]))
    # end
end

function Profile_Reader(Path)
    Profile_Data=Dict()
    Profile_Data["DA"]= load(Path*"DA_Profiles.jld2")["DA_Profiles"]
    Profile_Data["RT"]= load(Path*"RT_Profiles.jld2")["RT_Profiles"]
    return Profile_Data
end

function Profile_Generator(new_profile,True_Profile)
    if new_profile[1]
        DA_Profiles=DA_Profile_Generator(True_Profile,nTL;ferror=errortype,Data_=[("Load",ERORR["DA"][1]),("Wind",ERORR["DA"][2])])
    end
    if new_profile[2]
        RT_Profiles=RT_Profile_Generator(True_Profile;ferror=errortype,Data=[("Load",ERORR["RT"][1]),("Wind",ERORR["RT"][2])])
    elseif new_profile[2]=="same"
        RT_Profiles=RT_Profile_Generator(DA_Profiles;ferror=errortype,Data=[("Load",0),("Wind",0)])
    end
end

function Data_Reader(;Ptype="DA",init="DA_OPF",t_0=0)

    if Ptype=="DA" @assert init in ["DA_OPF","Zero"] && t_0==0 end
    if Ptype=="RT" @assert init in ["RT_OPF","DA_SCOPF"] && t_0>=1 end

    data=["Gen","RES", "Load", "Line", "Bus", "ConL", "ConG", "Trans", "Shunt","ESS"]
    profdata=["PL","QL","Wind"]
    initdata=["QG","PG","Pinc","Pdec","UE","Pcurt","Pcurtc","Pshed","Pch","Pdis","soc","Tcost","Vmag","Vang"]

    raw= Dict{String,Any }()
    for x in data
        if x in ["Line","ESS","Trans"]
            raw[x]=SysData[x][SysData[x]."br_status" .== 1,:]
        else
            raw[x]=SysData[x]
        end
        raw["n"*x]=size(raw[x],1)
    end

    raw["tratio"]=raw["Trans"]."ratio";
    raw["sbase"]=SysData["Base_MVA"]."base_MVA"[1]
    raw["LScost"]=SysData["Base_MVA"]."Lscost"[1]
    raw["RCcost"]=SysData["Base_MVA"]."Rccost"[1]
    #Felixable Load
    iFL  =  findall(x->x==1,raw["Load"][!,"ShFrOK"])
    raw["FL"]=DataFrame(Bus=raw["Load"][iFL,"Bus_N"],Cost=raw["Load"][iFL,"ShCost"],iFL=iFL)
    raw["nFL"]     =  size(iFL,1)

    initraw=Dict{String, DataFrame}()
    for x in initdata
        if Ptype=="DA" || init=="RT_OPF"
            initraw[x]=Init_Data[init][x][1:nTP,2:end]
        else
            initraw[x]=Init_Data[init][x][[(Int(floor((t_0+i-2)*TD))+1) for i=1:nTP if (Int(floor((t_0+i-2)*TD))+1) <25],2:end]
        end
    end

    profraw=Dict{String, DataFrame}()
    for x in profdata
        if Ptype=="DA"
            profraw[x]= Profile_Data[Ptype][x]
        else
            profraw[x]= Profile_Data[Ptype][t_0,x]
        end
    end

    # raw["Scen"]=profraw["Wind"][:,1]
    # raw=DataFrame(raw)[1,:]
    return raw,initraw,profraw
end

function Load_Data(nTP,raw,profraw,ctdall)
    Proffin=Dict{String,Array{Float64}}()
    #------------------Power Demand with load index
    Proffin["PL"]=(Matrix(profraw["PL"])*(raw["Load"]."Pd (pu)")')[1:nTP,:]
    Proffin["QL"]=(Matrix(profraw["QL"])*(raw["Load"]."Qd (pu)")')[1:nTP,:]

    # -----------------Power Demand in Buses
    ctdload = @view ctdall."load"[:]
    Proffin["PD"]=zeros(nTP,raw["nBus"])
    Proffin["QD"]=zeros(nTP,raw["nBus"])
    Proffin["PD"]=[sum(init=0,Proffin["PL"][t,j] for j in ctdload[i]) for t in 1:nTP,i in 1:raw["nBus"]]
    Proffin["QD"]=[sum(init=0,Proffin["QL"][t,j] for j in ctdload[i]) for t in 1:nTP,i in 1:raw["nBus"]]

    # Proffin["PD1"]=[sum(init=0,Proffin["PL"][t,j] for j in 1:raw["nLoad"] if raw["Load"][j,"Bus_N"]==i) for t in 1:nTP,i in 1:raw["nBus"]]
    # Proffin["QD1"]=[sum(init=0,Proffin["QL"][t,j] for j in 1:raw["nLoad"] if raw["Load"][j,"Bus_N"]==i) for t in 1:nTP,i in 1:raw["nBus"]]
    # # PG0=[sum(init=0,raw["P_Prof_G"][j,t] for j in 1:raw["nGen"] if raw["Gen"][j,"Bus"]==i) for t in 1:nTP,i in 1:raw["nBus"]]
    # QG0=[sum(init=0,raw["Q_Prof_G"][j,t] for j in 1:raw["nGen"] if raw["Gen"][j,"Bus"]==i) for t in 1:nTP,i in 1:raw["nBus"]]

    # PG0=[0 for t in 1:nTP,i in 1:raw["nBus"]]
    # QG0=[0 for t in 1:nTP,i in 1:raw["nBus"]]
    # Renewable Data
    PRES = zeros(Float64,(nTP,raw["nRES"]))
    for i=1:raw["nRES"]     PRES[:,i]=raw["RES"]."Cap"[i]*Vector(profraw["Wind"][:,1])[1:nTP]            end
    Proffin["Wind"]=PRES
    return Proffin
end
#####################################################################
####################### --------------------- #######################
#####################################################################
####################### Output write Function #######################

function ODSW(address,out)
    r=2;
    for i in names(out)
        if typeof(out[!,i])!=Vector{DataFrame}
            ods_write(address,Dict(("Main",r,1)=>[[i] out[!,i]]))
            r+=1;
        else
            ods_write(address,Dict((i,1,2)=>out[!,i][1]))
            if i!="C_report"
                ods_write(address,Dict((i,1,1)=>DataFrame((t=["t$i" for i in 1:nTP]))))
            end
        end
    end
end

function set_initial(name,out)
    outlist=["PG", "QG","Pinc","Pdec","UE","Pcurt","Pcurtc","Pshed","Pch","Pdis","soc","Tcost","Vmag","Vang"];
    for i in outlist
        ods_write(name,Dict((i,1,1)=> DataFrame((t=["t$i" for i in 1:nTP]))))
        ods_write(name,Dict((i,1,2)=>out[!,i][1]))
    end
end

function print_report()
    print("\n Normal operation: \n\n")
    print(Violation_rep()[1]["summary"],"\n")
    print("\n Worst Contingency case: \n\n")
    print(raw["ConL"][out.C_report[1].Con_N[1],:],"\n\n")
    print(Violation_rep( out.C_report[1].Con_N[1])[1]["summary"],"\n")
end

function Violation_rep(c=0)
    L_Bus=filter(x->x."Bus_Type"==1,raw["Bus"])."Bus_N"
    if c==0
        Vcon=SCOPF[:NNvconstraint]
        Fcon=SCOPF[:NNlflow]
        maxflow=sort(myfilter(x-> dual(x)<-0.01,Fcon),by=dual)
        E_maxflow=sort(myfindall(x-> dual(x)<-0.01,Fcon),by=x->dual(Fcon[x]))
    else
        Vcon=SCOPF[:CNvconstraint][:,:,c]
        Fcon=SCOPF[:CNlflow][:,:,c]
        maxflow=sort(filter(x-> dual(x)<-0.01,Fcon),by=dual)
        E_maxflow=sort(findall(x-> dual(x)<-0.01,Fcon),by=x->dual(Fcon[x]))
    end
    minV=sort(filter(x->dual(x)>0.01,Vcon[:,L_Bus]),by=dual,rev=true)
    maxV=sort(filter(x->dual(x)<-0.01,Vcon[:,L_Bus]),by=dual)

    E_minV=sort(findall(x->dual(x)>0.01,Vcon[:,L_Bus]),by=x->dual(Vcon[x]),rev=true)
    E_maxV=sort(findall(x->dual(x)<-0.01,Vcon[:,L_Bus]),by=x->dual(Vcon[x]))

    maxV_T=Set(i[1] for i in E_maxV)
    maxV_B=Set(i[2] for i in E_maxV)
    minV_T=Set(i[1] for i in E_minV)
    minV_B=Set(i[2] for i in E_minV)
    maxflow_T=Set(i[1] for i in E_maxflow)
    maxflow_L=Set(i[2] for i in E_maxflow)

    report=Dict()
    report["summary"]=DataFrame(

                type = ["max voltage violation","min voltage violation","line flow violation","total"],
                num_component= [length(maxV_B),length(minV_B),length(maxflow_L),length(maxV_B)+length(minV_B)+length(maxflow_L)],
                num_time=[ length(maxV_T),length(minV_T),length(maxflow_T),NaN],
                num_violation=[length(maxV),length(minV),length(maxflow),length(maxV)+length(minV)+length(maxflow)],
                )

    report["Elements"]=DataFrame(
                type = ["max voltage violation","min voltage violation","line flow violation"],
                bus = [tuple(maxV_B...),tuple(minV_B...),tuple(maxflow_L...)],
                time = [tuple(maxV_T...),tuple(minV_T...),tuple(maxflow_T...)],
                violation= Set.(map(x-> cartconvert.(x) ,[E_maxV,E_minV,E_maxflow])),
                )

    return report,(minV,maxV,maxflow,E_minV,E_maxV,E_maxflow,maxV_T,maxV_B,minV_T,minV_B,maxflow_T,maxflow_L)
end

function create_con_report(c,SCOPF,raw)
    if c!=0
        N_sheets=["cost","real_v","imag_v","mag_v","ang_v","Line_flow","PG","QG","Pdis","Pch","soc","Pinc","Pdec","Pcurt","Pcurtc","Pshed"]
        N_out=Get_ConResult(c,SCOPF,raw)
        for i in N_sheets
            ods_write("Result/"*Test_Sys*"CLR/"*string(c),Dict((i,1,2)=>N_out[!,i][1]))
            ods_write("Result/"*Test_Sys*"CLR/"*string(c),Dict((i,1,1)=>DataFrame((t=["t$i" for i in 1:nTP]))))
        end
    end
    reg=Violation_rep(c)[1]["Elements"]
    reg.violation=setdiff.(Violation_rep(c)[1]["Elements"].violation,Violation_rep(0)[1]["Elements"].violation)
    reg2=Violation_rep(c)[1]["Elements"]
    reg2.violation=setdiff.(Violation_rep(0)[1]["Elements"].violation,Violation_rep(c)[1]["Elements"].violation)
    out=DataFrame(
        violation= Violation_rep(c)[1]["Elements"],
        summary= Violation_rep(c)[1]["summary"],
        newV=reg[:,[1,4]],
        releavedV=reg2[:,[1,4]],
                )
    sheets=["summary","violation","newV","releavedV"]
    for i in sheets
        ods_write("Result/"*Test_Sys*"CLR/"*string(c),Dict((i,1,1)=>out[!,i][1]))
    end

    return out
end
#####################################################################
####################### --------------------- #######################
#####################################################################
####################### Miscellaneous Function #######################

function dual_check(model,key1="",key2="")
    g=all_constraints(model,include_variable_in_set_constraints=true)
    VC=filter(x->abs(dual(x))>=0.1 && occursin(key1,name(x)) && occursin(key2,name(x)),g)
    return VC
end

function myfindall(condition, array)
    buffer=[]
    for (i,j) in eachindex(array)
        if condition(array[i,j])
            push!(buffer,(i,j))
        end
    end
    return buffer
end
function myfilter(condition, array)
    buffer=[]
    for (i,j) in eachindex(array)
        if condition(array[i,j])
            push!(buffer,array[i,j])
        end
    end
    return buffer
end
function cartconvert(x)
    return (x[1],x[2])
end

#####################################################################
####################### --------------------- #######################
#####################################################################
####################### Profile maker Function #######################
function True_Profile_Reader()
    profdata=["PL","QL","Wind"];realprof=Dict();
    for x in profdata
        try
            realprof[x] = ods_read(Path * "Scenarios", sheetName = x, retType = "DataFrame")[:, 2:end][:, 1:1]
        catch
            error("Scenario file does not exist.") end
    end
    return realprof
end


function DA_Profile_Generator(True_Profile,nTS=24;ferror="Normal",Data_=[("Load",0.1),("Wind",0.1)])

    profdata=["PL","QL","Wind"];realprof=True_Profile;
    newprof=Dict();
    for (D,P) in Data_
        if ferror=="Normal"     noise=((P/4)*randn(nTS))
        elseif ferror=="Uniform"    noise=(P*(2rand(nTS).-1)) end

        if D=="Load"
            for v in ["PL","QL"]   newprof[v]=round.((noise.+1).*realprof[v][1:nTS,:],digits=4) end
        elseif D=="Wind"
            newprof[D]=round.((noise.+1).*realprof[D][1:nTS,:],digits=4)
        end
    end

    #write to excel
    for i in profdata
        ods_write(Path*"DA_Profiles",Dict((i,1,2)=>newprof[i]))
        ods_write(Path*"DA_Profiles",Dict((i,1,1)=>DataFrame((t=["t$j" for j in 1:nTS]))))
    end
    #write to jld2
    DA_Profiles=newprof
    @save Path*"DA_Profiles.jld2" DA_Profiles
    return newprof
end

function RT_Profile_Generator(True_Profile;MPC="linear",ferror="Normal",Data=[("Load",0.05),("Wind",0.05)])
    TD=CCT
    RT_Profiles=DataFrame(t=1:nTL*nTS,PL=empty(DataFrame()),QL=empty(DataFrame()), Wind=empty(DataFrame()))
    profdata=["PL","QL","Wind"];realprof=True_Profile;
    Load=["PL","QL"]


    for t0=1:nTL*nTS
        nTP=min(nTS*nP,nTL*nTS-t0+1)
        newprof=Dict()
        temp=Dict()
        t=Int.([floor((t0+i-2)*TD)+1 for i=1:nTP])

        for D in profdata temp[D]= realprof[D][t,:] end
        for (D,P) in Data
            if ferror=="Normal"
                noise=((P/4)*randn(nTP))
            elseif ferror=="Uniform"
                noise=(P*(2rand(nTP).-1))
            end

            if MPC=="nothing"
                noisempc=noise
            elseif MPC=="exp"
                noisempc= noise.*[exp(0.3*(t-12)) for t=1:nTP]
            elseif MPC=="linear"
                noisempc= noise.*[t/12 for t=1:nTP]
            end

            if D=="Load"
                for v in Load
                    newprof[v]=round.((noisempc.+1).*temp[v],digits=4)
                end
            else
                newprof[D]=round.((noisempc.+1).*temp[D],digits=4)
            end
        end

        for i in profdata
            RT_Profiles[t0,i]=newprof[i]
        end
    end
    @save Path*"RT_Profiles.jld2" RT_Profiles#---> saved in Path*"RT_Profiles
    return RT_Profiles
end
#####################################################################
####################### --------------------- #######################
#####################################################################
########################### Plot Function ###########################

function plot_cost_time(out)
    y=out."T_report"[1][:,1]
    tlim=size(y,1)
    x=range(1,tlim,length=tlim)
    plot(x,y,
        title="Expected cost per time",
        label="Expected cost",
        linewidth=3,
        ylabel="Euro",
        xlabel="Time (h)",
        ylims=(0,round(findmax(y)[1]*1.1)),
        xminorticks=5,
        )
    savefig("Figures/cost-time")
end

function plot_cost_con(out)
    y=out."C_report"[1][:,2]
    tlim=size(y,1)
    x=range(1,tlim,length=tlim)
    plot(x,y,
        title="Cost per contingencies",
        label="Cost",
        linewidth=3,
        ylabel="Euro",
        xlabel="Time (h)",
        ylims=(0,round(findmax(y)[1]*1.1)),
        xminorticks=5,
        )
        savefig("Figures/cost-con")
end

function Plot_profiles(profraw,newprof)
    Data=["PL","QL","Wind"]
    Td=1/6;
    for d in Data
        a=Int.(floor.(range(0,24-Td,step=Td))).+1
        if d=="Wind"
            plot(range(1,144,step=1),profraw[d][a,1])
            plot!(range(1,144,length=144),newprof[d][:,1])
        else
            plot(range(1,144,step=1),sum.(eachrow(profraw[d]))[a],ylim=[minimum(sum.(eachrow(profraw[d])))*0.95,maximum(sum.(eachrow(profraw[d])))*1.05])
            plot!(range(1,144,length=144),sum.(eachrow(newprof[d]))[a],ylim=[minimum(sum.(eachrow(newprof[d])))*0.95,maximum(sum.(eachrow(newprof[d])))*1.05])
        end
        savefig("Figures/"*d)
    end

end

function Plotline(y::Vector,yleg::Vector;ver=1,them=:default,xticks1=0:12:144,title="title",ylab="ylab",xlab="xlab",background_color = RGB(0.85, 0.87, 0.91),fontsize1=12)
    theme(them)
    plot()
    # markers = [:none,:circle, :square, :diamond, :star]
    markers = [:none,:none, :none, :none, :none]
    colors=[:blue, :red, :green, :purple]
    # lines=[:solid, :dash, :dot, :dashdot]
    lines=[:solid,:solid,:solid,:solid]
    if ver==1
        tlim=size(y[1],1)
        x=range(1,tlim,length=tlim)
        for i in eachindex(y)
            plot!(x, y[i], label=yleg[i],xticks=xticks1 ,line=lines[i], color=colors[i],marker=markers[i])
        end
        # ylims=ylims,
        # xminorticks=5,
    elseif ver==2
        tlim=size(y[1],1)
        x=range(1,tlim,length=tlim)
        for i in eachindex(y)
            plot!(x, y[i], label=yleg[i],xticks=xticks1 ,line=lines[i], color=colors[i],marker=markers[i])
        end
    end

    xlabel!(xlab,fontsize=fontsize1)
    ylabel!(ylab,fontsize=fontsize1)
    title!(title,fontsize=fontsize1)
    # plot!(bgcolor=background_color)
end
#####################################################################
####################### --------------------- #######################
#####################################################################
########################### .... Function ###########################
function ctdMaker(raw)
    ctdall=DataFrame(linef = [0 for _ in 1:raw["nBus"]])
    # Lines and Transformers
    A = spzeros(Int64, raw["nBus"], raw["nLine"])
    B = spzeros(Int64, raw["nBus"], raw["nLine"])
    C = spzeros(Int64, raw["nBus"], raw["nTrans"])
    D = spzeros(Int64, raw["nBus"], raw["nTrans"])
    BusLineF=@view raw["Line"][:,"From"]
    BusLineT=@view raw["Line"][:,"To"]
    BusTransF=@view raw["Trans"][:,"From"]
    BusTransT=@view raw["Trans"][:,"To"]
    for j in axes(A,2)
        A[BusLineF[j], j] = 1
        B[BusLineT[j], j] = 1
    end
    for j in axes(C,2)
        C[BusTransF[j], j] = 1
        D[BusTransT[j], j] = 1
    end
    ctdall[!,:linef] = find_ones(A)
    ctdall[!,:linet] = find_ones(B)
    ctdall[!,:transf] = find_ones(C)
    ctdall[!,:transt] = find_ones(D)
    # Loads
    A = spzeros(Int64, raw["nBus"], raw["nLoad"])
    BusLoad=@view raw["Load"][:,"Bus_N"]
    for j in axes(A,2)
        A[BusLoad[j], j] = 1
    end
    ctdall[!,:load] = find_ones(A)
    # Generators
    A=spzeros(Int64, raw["nBus"], raw["nGen"])
    BusGen=@view raw["Gen"][:,"Bus"]
    for j in axes(A,2)
        A[BusGen[j], j] = 1
    end
    ctdall[!,:gen] = find_ones(A)
    # Shunts
    A=spzeros(Int64, raw["nBus"], raw["nShunt"])
    BusShunt=@view raw["Shunt"][:,"Bus"]
    for j in axes(A,2)
        A[BusShunt[j], j] = 1
    end
    ctdall[!,:shunt] = find_ones(A)
    # ESS
    A=spzeros(Int64, raw["nBus"], raw["nESS"])
    BusESS=@view raw["ESS"][:,"Bus"]
    for j in axes(A,2)
        A[BusESS[j], j] = 1
    end
    ctdall[!,:ess] = find_ones(A)
    # RES
    A=spzeros(Int64, raw["nBus"], raw["nRES"])
    BusRES=@view raw["RES"][:,"Bus"]
    for j in axes(A,2)
        A[BusRES[j], j] = 1
    end
    ctdall[!,:res] = find_ones(A)
    # FL
    A=spzeros(Int64, raw["nBus"], raw["nFL"])
    BusFL=@view raw["FL"][:,"Bus"]
    for j in axes(A,2)
        A[BusFL[j], j] = 1
    end
    ctdall[!,:fl] = find_ones(A)
    ctdall[!,:nbus] = 1:raw["nBus"]
    return ctdall
end

function find_variable_in_model_constraints(model,variable)

    g=all_constraints(model,include_variable_in_set_constraints=true)
    for i in g
        if coefficient(constraint_object(i).func,variable)!=0
            println(i)
        end
    end
end

function Make_Iter_report(OP_NLP)
    # V_1,θ_1,PG_1,QG_1,cost1=Matrix(out."mag_v"[1]),Matrix(out."ang_v"[1]),Matrix(out."PG"[1]),Matrix(out."QG"[1]),out."OFValue_WP"[1]
    V_1,θ_1,PG_1,QG_1,cost1=OP_NLP["V"],OP_NLP["θ"],OP_NLP["PG"],OP_NLP["QG"],OP_NLP["obj_WP"]
    if Problem_Type=="SCOPF"
        V_c_Line_1,θ_c_Line_1,PG_c_Line_1,QG_c_Line_1= OP_NLP["V_c_Line"],OP_NLP["θ_c_Line"],OP_NLP["PG_c_Line"],OP_NLP["QG_c_Line"]
        # sqrt.(value.(model[:e_c_Line].^2 +model[:f_c_Line].^2)),atan.(value.(model[:f_c_Line]),value.(model[:e_c_Line])),value.(model[:PG_c_Line]),value.(model[:QG_c_Line])
    end
    if Report=="Florin"
        Head=["Initialization";"Method";"Time";"Res_Time";"Iter";"TR";"P_R";"Q_R";"OptGap";"Optimizer_Obj";"real_Obj";"Max_I_N";"Max_I_C";"Max_Vmax_N";"Max_Vmax_C";"Max_Vmin_N";"Max_Vmin_C";"Max_Ramp_N";"Max_Ramp_C";"Smax_N_N";"Vmin_N_N";"Vmax_N_N";"Ramp_N_N";"Smax_C_N";"Vmin_C_N";"Vmax_C_N";"Ramp_C_N";]
        Iter_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0 OP_NLP["time"] "-" "-" "-" "-" "-" cost1 cost1 "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"]
    elseif Report=="feasiblity"
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"TR";"FI";"OptGap";"acc";"Max_Inf_N";"Max_Inf_C";"Mean_Inf_N";"Mean_Inf_C";"Num_Inf_N";"Num_Inf_C";"Optimizer_Obj";"real_Obj";"Smax_Abs";"Vmin_Abs";"Vmax_Abs";"Ramp_Abs";"Smax_C_Abs";"Vmin_C_Abs";"Vmax_C_Abs";"Ramp_C_Abs";"Smax_N";"Vmin_N";"Vmax_N";"Ramp_N";"Smax_C";"Vmin_C";"Vmax_C";"Ramp_C";"Smax_N_M";"Vmin_N_M";"Vmax_N_M";"Ramp_N_M";"Smax_C_M";"Vmin_C_M";"Vmax_C_M";"Ramp_C_M";"Smax_N_N";"Vmin_N_N";"Vmax_N_N";"Ramp_N_N";"Smax_C_N";"Vmin_C_N";"Vmax_C_N";"Ramp_C_N";"Num_C_I"; "Num_C_V";"Num_S_I"; "Num_S_V"]
        Iter_Rep=DataFrame(repeat([[]],length(Head)), Head)
         #------------------- Print Optimal point and first iteration
         newrow=["-" "Non-Linear" 0 OP_NLP["time"] "-" "-" "-" "-" 0 "-" 0 0 0 0 0 0 cost1 "-"  "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"]
    elseif Report=="full"
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"TR";"FI";"OptGap";"Max_Inf_N";"Max_Inf_C";"Optimizer_Obj";"real_Obj";["PG$i" for i in 1:size(SysData["Gen"],1)];["QG$i" for i in 1:size(SysData["Gen"],1)];["V$i" for i in 1:size(SysData["Bus"],1)];["θ$i" for i in 1:size(SysData["Bus"],1)]]
        Iter_Rep=DataFrame(repeat([[]],length(Head)), Head)
        #------------------- Print Optimal point and first iteration
        newrow=["-" "Non-Linear" 0 OP_NLP["time"] 0 1 0 0 0 0 0 cost1 "-" PG_1 QG_1 V_1 θ_1]
    elseif Problem_Type=="SCOPF"
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"TR";"FI";"OptGap";"Max_Inf_N";"Max_Inf_C";"Optimizer_Obj";"real_Obj";"PG_RMSE(p.u)";"QG_RMSE(p.u)";"V_RMSE(p.u)";"θ_RMSE(Rad)";"PG_c_RMSE(p.u)";"QG_c_RMSE(p.u)";"V_c_RMSE(p.u)";"θ_c_RMSE(Rad)";]
        Iter_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0 OP_NLP["time"] 0 1 0 0 0 0 0 "-" OP_NLP["obj_WP"]  0 0 0 0 0 0 0 0]
    elseif Problem_Type=="OPF"
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"TR";"FI";"OptGap";"Max_Inf_N";"Max_Inf_C";"Optimizer_Obj";"real_Obj";"PG_RMSE(p.u)";"QG_RMSE(p.u)";"V_RMSE(p.u)";"θ_RMSE(Rad)";]
        Iter_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0 OP_NLP["time"] 0 1 0 0 0 0 0 "-"  OP_NLP["obj_WP"] 0 0 0 0]
    end
    push!(Iter_Rep,newrow[1,:])
    return Iter_Rep
end

function Add_Iter_report!(Runcode,Iter_Rep,OP,cost,OptGap,Vio_Sum_N,Vio_Sum_C,Vio_general_N,Vio_general_C,iter,sol_time,res_time,Filter,R_P=0,R_Q=0;Convergence=true)
    # Runcode,Iter_Rep,out,OP,cost,OptGap,Vio_Sum_N,Vio_Sum_C,Vio_general_N,Vio_general_C,iter,sol_time,res_time=Runcode,Iter_Rep,init_out,Best_solution,mincost,minOptGap,Best_vio_sum_N,Best_vio_sum_C,Best_vio_general_N,Best_vio_general_C,iter,0,0
    max_N= [0.0,0.0, 0.0, 0.0]
    max_C=[0.0,0.0, 0.0, 0.0]
    mean_N=[0.0,0.0, 0.0, 0.0]
    mean_C=[0.0,0.0, 0.0, 0.0]
    Num_N=[0,0,0,0]
    Num_C=[0,0,0,0]
    Abs_N=[0.0,0.0, 0.0, 0.0]
    Abs_C=[0.0,0.0, 0.0, 0.0]
    if iter !=0 TRB=OP["FilterIndex"] else TR=0;TRB=0 end
    TR="["*join(string.(OP["TR"]),",")*"]"
    ACC="["*join(string.(OP["acc"]),",")*"]"
    if Convergence==true
        max_N[1:size(Vio_general_N[1,1],1)] .= maximum(hcat([df[!, :max] for df in Vio_general_N]...),dims=2)
        Num_N[1:size(Vio_general_N[1,1],1)] .= sum(hcat([df[!, :N] for df in Vio_general_N]...),dims=2)
        Abs_N[1:size(Vio_general_N[1,1],1)] .= maximum(maximum.(hcat([df[!, :abs] for df in Vio_general_N]...),init=0),dims=2)
        mean_N[1:size(Vio_general_N[1,1],1)] .= sum(hcat([df."mean" .*df."N" for df in Vio_general_N]...),dims=2)./ Num_N[1:size(Vio_general_N[1,1],1)]


        if Problem_Type=="SCOPF"
            max_C[1:size(Vio_general_C[1,1,1],1)] = maximum(hcat([df."max" for df in Vio_general_C]...),dims=2)
            Abs_C[1:size(Vio_general_C[1,1,1],1)] = maximum(maximum.(hcat([df."abs" for df in Vio_general_C]...),init=0),dims=2)
            Num_C[1:size(Vio_general_C[1,1,1],1)] = sum(hcat([df."N" for df in Vio_general_C]...),dims=2)
            mean_C[1:size(Vio_general_C[1,1,1],1)] = sum(hcat([df."mean" .*df."N" for df in Vio_general_C]...),dims=2)./ Num_C[1:size(Vio_general_C[1,1,1],1)]
        end
    end

    if Problem_Type == "OPF"
        Vio_Sum_C=DataFrame(Max=0.0,abs=0.0,Mean=0.0,NVio=0)
    end

    if Report=="Florin"
        R_P_T="["*join(string.(R_P),",")*"]"
        R_Q_T="["*join(string.(R_Q),",")*"]"
        # Head=["Initialization";"Method";"Time";"Res_Time";"Iter";"TR";"OptGap";"Optimizer_Obj";"real_Obj";"Max_I_N";"Max_I_C";"Max_Vmax_N";"Max_Vmax_C";"Max_Vmin_N";"Max_Vmin_C";"Max_Ramp_N";"Max_Ramp_C";"Smax_N_N";"Vmin_N_N";"Vmax_N_N";"Ramp_N_N";"Smax_C_N";"Vmin_C_N";"Vmax_C_N";"Ramp_C_N";]
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  sol_time res_time iter TR R_P_T R_Q_T OptGap OP["obj_WP"] round(cost,digits=5) 100*Abs_N[3] 100*Abs_C[3] Abs_N[1] Abs_C[1] Abs_N[2] Abs_C[2] Abs_N[4] Abs_C[4] Num_N[3] Num_C[3] Num_N[1] Num_C[1] Num_N[2] Num_C[2] Num_N[4] Num_C[4]]
    elseif Report=="feasiblity"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time res_time iter TR TRB OptGap ACC Vio_Sum_N."Max"[1] Vio_Sum_C."Max"[1] Vio_Sum_N."Mean"[1] Vio_Sum_C."Mean"[1] Vio_Sum_N."NVio"[1] Vio_Sum_C."NVio"[1] OP["obj_WP"] round(cost,digits=5) Abs_N[3] Abs_N[2] Abs_N[1] Abs_N[4] Abs_C[3] Abs_C[2] Abs_C[1] Abs_C[4] max_N[3] max_N[2] max_N[1] max_N[4] max_C[3] max_C[2] max_C[1] max_C[4] mean_N[3] mean_N[2] mean_N[1] mean_N[4] mean_C[3] mean_C[2] mean_C[1] mean_C[4] Num_N[3] Num_N[2] Num_N[1] Num_N[4] Num_C[3] Num_C[2] Num_C[1] Num_C[4] Tuple(Filter[1]) Tuple(Filter[2]) Tuple(Filter[3]) Tuple(Filter[4])]
    elseif Report=="full"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time res_time iter TR TRB OptGap Vio_Sum_N."Max"[1] Vio_Sum_C."Max"[1] OP["obj_WP"] round(cost,digits=5) round.(OP["PG"],digits=5) round.(OP["QG"],digits=5) round.(OP["V"],digits=5) round.(OP["θ"],digits=5)]
    elseif Problem_Type=="SCOPF"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time res_time iter TR TRB OptGap Vio_Sum_N."Max"[1] Vio_Sum_C."Max"[1] OP["obj_WP"] round(cost,digits=5) sqrt(mean((round.(OP["PG"],digits=5).-PG_1).^2)) sqrt(mean((round.(OP["QG"],digits=5).-QG_1).^2)) sqrt(mean((round.(OP["V"],digits=5).-V_1).^2)) sqrt(mean((round.(OP["θ"],digits=5).-θ_1).^2)) sqrt(mean((round.(OP["PG_c_Line"],digits=5).-PG_c_Line_1).^2)) sqrt(mean((round.(OP["QG_c_Line"],digits=5).-QG_c_Line_1).^2)) sqrt(mean((round.(OP["V_c_Line"],digits=5).-V_c_Line_1).^2)) sqrt(mean((round.(OP["θ_c_Line"],digits=5).-θ_c_Line_1).^2))]
    elseif Problem_Type == "OPF"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time res_time iter TR TRB OptGap Vio_Sum_N."Max"[1] Vio_Sum_C."Max"[1] OP["obj_WP"] round(cost,digits=5) sqrt(mean((round.(OP["PG"],digits=5).-PG_1).^2)) sqrt(mean((round.(OP["QG"],digits=5).-QG_1).^2)) sqrt(mean((round.(OP["V"],digits=5).-V_1).^2)) sqrt(mean((round.(OP["θ"],digits=5).-θ_1).^2))]
    end
    push!(Iter_Rep,newrow[1,:])
    return Iter_Rep
end

# function Set_TR!(OP_init,Vio_general_C,Vio_general_N)

#     max_N= [0.0,0.0, 0.0, 0.0]
#     max_C=[0.0,0.0, 0.0, 0.0]
#     # if iter !=0 TRB=OP_init["FilterIndex"] else TR=0;TRB=0 end
#     # TR="["*join(string.(OP_init["TR"]),",")*"]"
#     # findmax(hcat([df[!, :max] for df in Vio_general_N]...)[1:2,:])
#     # findmax(hcat([df."max" for df in Vio_general_C]...)[1:2,:])
#     max_N[1:size(Vio_general_N[1,1],1)] .= maximum(hcat([df[!, :max] for df in Vio_general_N]...),dims=2)
#     if Problem_Type=="SCOPF"
#         max_C[1:size(Vio_general_C[1,1,1],1)] = maximum(hcat([df."max" for df in Vio_general_C]...),dims=2)
#     end

#     Vdev_N=maximum(max_N[1:2])
#     Vdev_C=maximum(max_C[1:2])

#     Pdev_N=maximum(max_N[3:4])
#     Pdev_C=maximum(max_C[3:4])

#     max_dev=[Vdev_N,Vdev_C,Pdev_N,Pdev_C]
#     OP_init["MaxDev"]=max_dev
#     return OP_init
# end

function find_total_redispatch!(Best_solution,OP_init0)
    maxnormalredispatchP=0
    maxconredispatchP=0
    Gres,ESSres,Curtres,Shedres,Incres= zeros(nTP,1),zeros(nTP,1),zeros(nTP,1),zeros(nTP,1),zeros(nTP,1)
    GresQ,ShedresQ= zeros(nTP,1),zeros(nTP,1)

    if !isempty(Best_solution["PG"])
        Gres=  sum(abs.(Best_solution["PG"]-OP_init0["PG"]),dims=2)
        GresQ= sum(abs.(Best_solution["QG"]-OP_init0["QG"]),dims=2)
    end
    if !isempty(Best_solution["Pch"])
        ESSres=sum(abs.(Best_solution["Pch"]-Best_solution["Pdis"]-OP_init0["Pch"] + OP_init0["Pdis"]),dims=2)
    end
    if !isempty(Best_solution["Pcurt"])
        Curtres=sum(abs.(Best_solution["Pcurt"]-OP_init0["Pcurt"]),dims=2)
    end
    if !isempty(Best_solution["Pshed"])
        Shedres=sum(abs.(Best_solution["Pshed"]-OP_init0["Pshed"]),dims=2)
        ShedresQ=sum(abs.((Best_solution["Pshed"]-OP_init0["Pshed"]).* reshape(SysData["Load"]."PF" , 1, :)),dims=2)
    end
    if !isempty(Best_solution["Pinc"])
        Incres=sum(abs.(Best_solution["Pinc"]-Best_solution["Pdec"]-OP_init0["Pinc"]+OP_init0["Pdec"]),dims=2)
    end
    maxnormalredispatchP=maximum(Gres+ESSres+Curtres+Shedres+Incres)
    maxnormalredispatchQ=maximum(GresQ+ShedresQ)

    if Problem_Type=="SCOPF"
        ncon=size(SysData["ConL"],1)
        GresC,ESSresC,CurtresC,ShedresC,IncresC= zeros(nTP,ncon),zeros(nTP,ncon),zeros(nTP,ncon),zeros(nTP,ncon),zeros(nTP,ncon)
        GresCQ,ShedresCQ= zeros(nTP,ncon),zeros(nTP,ncon)
        if !isempty(Best_solution["PG_c_Line"])
            GresC=sum(abs.(Best_solution["PG_c_Line"]-OP_init0["PG_c_Line"]),dims=2)[:,1,:]
            GresCQ=sum(abs.(Best_solution["QG_c_Line"]-OP_init0["QG_c_Line"]),dims=2)[:,1,:]
        end
        if !isempty(Best_solution["Pch_c_Line"])
            ESSresC=sum(abs.(Best_solution["Pch_c_Line"]-Best_solution["Pdis_c_Line"]-OP_init0["Pch_c_Line"] + OP_init0["Pdis_c_Line"]),dims=2)[:,1,:]
        end
        if !isempty(Best_solution["Pcurt_c_Line"])
            CurtresC=sum(abs.(Best_solution["Pcurt_c_Line"]-OP_init0["Pcurt_c_Line"]),dims=2)[:,1,:]
        end
        if !isempty(Best_solution["Pshed_c_Line"])
            ShedresC=sum(abs.(Best_solution["Pshed_c_Line"]-OP_init0["Pshed_c_Line"]),dims=2)[:,1,:]
            ShedresCQ=sum(abs.((Best_solution["Pshed_c_Line"]-OP_init0["Pshed_c_Line"]).* reshape(SysData["Load"]."PF" , 1, :)),dims=2)[:,1,:]
        end
        if !isempty(Best_solution["Pinc_c_Line"])
            IncresC=sum(abs.(Best_solution["Pinc_c_Line"]-Best_solution["Pdec_c_Line"]-OP_init0["Pinc_c_Line"]+OP_init0["Pdec_c_Line"]),dims=2)[:,1,:]
        end
        maxconredispatchP=maximum(GresC+ESSresC+CurtresC+ShedresC+IncresC)
        maxconredispatchQ=maximum(GresCQ+ShedresCQ)


    end
    return maxnormalredispatchP,maxconredispatchP,maxnormalredispatchQ,maxconredispatchQ
end

# function SET_TR2!(Lininit,SENSDict)
#     # devcapV_pos=minimum([quantile(sum(sparse(findnz(SENSDict["S_VQ"][i])...).*(SENSDict["S_VQ"][i] .> 0),dims=2)[:,1],0.9) for i=1:nTP])
#     # devcapV_neg=minimum([quantile(abs.(sum(sparse(findnz(SENSDict["S_VQ"][i])...).*(SENSDict["S_VQ"][i] .< 0),dims=2)[:,1]),0.9) for i=1:nTP])

#     devcapV_pos=minimum([sum(eachrow(filter(x-> x>=0,SENSDict["S_VQ"][i]))) for i=1:nTP])
#     devcapV_neg=minimum([abs.(sum(eachrow(filter(x-> x<=0,SENSDict["S_VQ"][i])))) for i=1:nTP])
#     SENSDict1=deepcopy(SENSDict)
#     A=SENSDict1["S_IP"][i]
#     sum(sparse(findnz(A)...).*(A .> 0),dims=2)

#     filter!(x-> x<0,eachrow(SENSDict1["S_IP"][1])[1])
#     devcapV_neg=minimum([minimum(abs.(sum(eachrow(filter(x-> x<=0,SENSDict["S_IP"][1]))))) for i=1:nTP])

#     maxvio_op=maximum(Lininit["MaxDev"][1:2])
#     Lininit["TR"][2]=maxvio_op/minimum([devcapV_pos,devcapV_neg])
# end

function sumJuMPi(A)
    answer=zero(AffExpr)
    if !isempty(A)
        for i in A
            add_to_expression!(answer,i)
        end
    end
    return answer
end

function filter_to_sparse(A::Matrix{T}, f::Function,multiplier=1) where T
    # f=  x -> abs(x) .> ( (acc[1]/10)./(TR1[t,:,c]))
    # A=SENSDict["S_IP"][t,c][lines,:]
    # A.>=( (acc[1]/10)./(TR1[t,:,c]))
    idx = findall(f, A.*multiplier)
    rows = [I[1] for I in idx]
    cols = [I[2] for I in idx]
    values = A[idx]
    return sparse(rows, cols, values, size(A)...)
end

function filter_to_sparse_V_2(A::Matrix{T}, error_,multiplier=1;neighborhood=size(A,2),multiplier2=0,neighborhood2=20) where T
    # A,error_,multiplier=SENSDict["S_VQ"][t,c][elems,:],filter_acc,TRact_2
    # A,error_,multiplier=SENSDict["S_IP"][t,c][lines,:],filter_acc,TRact_1./(2*Imaxlimit0_sqrt[t,lines,c]);neighborhood=80
    AAA=sparse(delete_elements!(abs.(A.*multiplier),error_,neighborhood))
    if multiplier2 != 0
        BBB=sparse(delete_elements!(abs.(A.*multiplier2),error_,neighborhood2))
    else
        BBB=AAA
    end
    return ((AAA.!=0) .| (BBB.!=0)).*A
end
function filter_to_sparse_V_3(A::Matrix{T}, error_,multiplier=[1,1];neighborhood=[size(A,2),size(A,2)]) where T
    # A,error_,multiplier=SENSDict["S_VQ"][t,c][elems,:],filter_acc,TRact_2
    # A,error_,multiplier=SENSDict["S_IP"][t,c][lines,:],filter_acc,TRact_1./(2*Imaxlimit0_sqrt[t,lines,c]);neighborhood=80
    AAA=sparse(delete_elements!(abs.(A.*multiplier[1]),error_,neighborhood[1]))
    CCC=sparse(delete_elements!(abs.(A.*multiplier[2]),error_,neighborhood[2]))
    # if multiplier2 != 0
    #     BBB=sparse(delete_elements!(abs.(A.*multiplier2),error_,neighborhood2))
    # else
    #     BBB=AAA
    # end
    return ((AAA.!=0) .| (CCC.!=0)).*A
end
function filter_to_sparse_V_4(A::Matrix{T}, error_,multiplier=[1,1];neighborhood=[size(A,2),size(A,2)]) where T
    # A,error_,multiplier=SENSDict["S_VQ"][t,c][elems,:],filter_acc,TRact_2
    # A,error_,multiplier=SENSDict["S_IP"][t,c][lines,:],filter_acc,TRact_1./(2*Imaxlimit0_sqrt[t,lines,c]);neighborhood=80
    AAA=sparse(delete_elements_pareto!(abs.(A.*multiplier[1]),error_,neighborhood[1]))
    CCC=sparse(delete_elements_pareto!(abs.(A.*multiplier[2]),error_,neighborhood[2]))
    # if multiplier2 != 0
    #     BBB=sparse(delete_elements!(abs.(A.*multiplier2),error_,neighborhood2))
    # else
    #     BBB=AAA
    # end
    return ((AAA.!=0) .| (CCC.!=0)).*A
end

function filter_to_sparse_V_5(A::Vector{Matrix{T}}, error_,multiplier=[1,1,1,1];neighborhood=[max(size(A[1],2),size(A[2],2)),max(size(A[1],2),size(A[2],2))]) where T
    # A,error_,multiplier=[SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],0.05,[TRact_1, TRact_1_C , TRact_2, TRact_2_C]

    # hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])
    # hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])
    AAA=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])),error_,neighborhood[1]))
    CCC=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])),error_,neighborhood[2]))
    coefs=((AAA.!=0) .| (CCC.!=0)).*hcat(A[1],A[2])
    return coefs[:,1:size(A[1],2)],coefs[:,size(A[1],2)+1:end]
end
function filter_to_sparse_V_8(A::Vector{Matrix{T}}, error_,multiplier=[1,1,1,1];neighborhood=[max(size(A[1],2),size(A[2],2)),max(size(A[1],2),size(A[2],2))]) where T
    # A,error_,multiplier=[SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],0.05,[TRact_1, TRact_1_C , TRact_2, TRact_2_C]

    # hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])
    # hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])
    AAA=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])),error_,neighborhood[1]))
    CCC=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])),error_,neighborhood[2]))
    DDD=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2])),error_,neighborhood[1]))
    # FFF=sparse(round.(hcat(multiplier[2],multiplier[4]),digits=3))
    
    coefs=((AAA.!=0) .| (CCC.!=0) .| (DDD.!=0)).*hcat(A[1],A[2])
    return coefs[:,1:size(A[1],2)],coefs[:,size(A[1],2)+1:end]
end
function filter_to_sparse_V_13(A::Vector{Matrix{T}}, error_,multiplier=[1,1,1,1];neighborhood=[max(size(A[1],2),size(A[2],2)),max(size(A[1],2),size(A[2],2))]) where T
    # A,error_,multiplier=[SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],0.05,[TRact_1, TRact_1_C , TRact_2, TRact_2_C]

    # hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])
    # hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])
    AAA=sparse(delete_elements_pareto2!(abs.(hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])),error_,neighborhood[1]))
    CCC=sparse(delete_elements_pareto2!(abs.(hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])),error_,neighborhood[2]))
    # DDD=sparse(delete_elements_pareto2!(abs.(hcat(A[1],A[2])),error_,neighborhood[1]))
    # FFF=sparse(round.(hcat(multiplier[2],multiplier[4]),digits=3))
    
    coefs=((AAA.!=0) .| (CCC.!=0)).*hcat(A[1],A[2])
    return coefs[:,1:size(A[1],2)],coefs[:,size(A[1],2)+1:end]
end

function filter_to_sparse_V_7(A::Vector{Matrix{T}}, error_,multiplier=[1,1,1,1];neighborhood=[max(size(A[1],2),size(A[2],2)),max(size(A[1],2),size(A[2],2))]) where T
    # A,error_,multiplier=[SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],0.05,[TRact_1, TRact_1_C , TRact_2, TRact_2_C]

    # hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])
    # hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])
    AAA1=sparse(delete_elements_pareto!(abs.(hcat(A[1]).* hcat(multiplier[1])),error_,neighborhood[1]))
    AAA2=sparse(delete_elements_pareto!(abs.(hcat(A[2]).* hcat(multiplier[3])),error_,neighborhood[1]))
    AAA3=hcat(AAA1,AAA2)
    CCC1=sparse(delete_elements_pareto!(abs.(hcat(A[1]).* hcat(multiplier[2])),error_,neighborhood[2]))
    CCC2=sparse(delete_elements_pareto!(abs.(hcat(A[2]).* hcat(multiplier[4])),error_,neighborhood[2]))
    CCC3=hcat(CCC1,CCC2)
    AAA=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])),error_,neighborhood[1]))
    CCC=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])),error_,neighborhood[2]))
    coefs=((AAA.!=0).| (AAA3.!=0) .| (CCC.!=0) .| (CCC3.!=0)).*hcat(A[1],A[2])
    return coefs[:,1:size(A[1],2)],coefs[:,size(A[1],2)+1:end]
end

function filter_to_sparse_V_10(A::Vector{Matrix{T}}, error_,multiplier=[1,1,1,1];neighborhood=[max(size(A[1],2),size(A[2],2)),max(size(A[1],2),size(A[2],2))]) where T
    # A,error_,multiplier=[SENSDict["S_IP"][t,c][lines,:],SENSDict["S_IQ"][t,c][lines,:]],0.05,[TRact_1, TRact_1_C , TRact_2, TRact_2_C]

    # hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])
    # hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])
    # AAA=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2]).* hcat(multiplier[1],multiplier[3])),error_,neighborhood[1]))
    # CCC=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2]).* hcat(multiplier[2],multiplier[4])),error_,neighborhood[2]))
    DDD=sparse(delete_elements_pareto!(abs.(hcat(A[1],A[2])),error_,neighborhood[1]))

    # minimum(hcat(A[1],A[2]).-((DDD.!=0)).*hcat(A[1],A[2]))
    # maximum(hcat(A[1],A[2]))

    
    coefs=((DDD.!=0)).*hcat(A[1],A[2])
    return coefs[:,1:size(A[1],2)],coefs[:,size(A[1],2)+1:end]
end


function delete_elements_pareto!(mat::Matrix{Float64}, threshold,neighborhood)
    # mat,threshold=abs.(A.*multiplier),error_
    m, n = size(mat)
    # neighborhood=10
    for row in 1:m
        sum_deleted = 0.0
        idx = sortperm(mat[row, :])  # Get sorted indices (ascending)
        if length(threshold)>1
            maxlevel=threshold[row]
        else
            maxlevel= threshold* sum(mat[row, :])
        end
        # idx_rev = sortperm(mat[row, :],rev=:true)
        # Vectorized loop replacement
        # plot(cumsum(mat[row, idx]))
        # (mat[row, idx_rev])
        # plot(cumsum(mat[row, idx_rev]).>=3*threshold)
        # plot(cumsum(mat[row, idx]) .<=threshold)
        # plot(cumsum(mat[row, idx])./threshold)
        # Ne_index=similar(mask1)
        # cumsum(mat[row, idx])[end-30:end]
        mask1 = cumsum(mat[row, idx]) .<= maxlevel  # Identify elements to delete
        mask2=similar(mask1)
        mask2.=1
        mask2[max(begin,end-3*neighborhood):end].=0
        mask = mask1.|mask2

        sum_deleted = sum(mat[row, idx[mask]])      # Compute total deleted sum
        mat[row, idx[mask]] .= 0.0                  # Set selected elements to 0
    end

    return mat
end

function delete_elements_pareto2!(mat::Matrix{Float64}, threshold,neighborhood)
    # mat,threshold=abs.(A.*multiplier),error_
    m, n = size(mat)
    # neighborhood=10
    for row in 1:m
        sum_deleted = 0.0
        idx = sortperm(mat[row, :])  # Get sorted indices (ascending)
        if length(threshold)>1
            maxlevel=threshold[row]
        else
            maxlevel= threshold
        end
        # idx_rev = sortperm(mat[row, :],rev=:true)
        # Vectorized loop replacement
        # plot(cumsum(mat[row, idx]))
        # (mat[row, idx_rev])
        # plot(cumsum(mat[row, idx_rev]).>=3*threshold)
        # plot(cumsum(mat[row, idx]) .<=threshold)
        # plot(cumsum(mat[row, idx])./threshold)
        # Ne_index=similar(mask1)
        # cumsum(mat[row, idx])[end-30:end]
        mask1 = cumsum(mat[row, idx]) .<= maxlevel  # Identify elements to delete
        mask2=similar(mask1)
        mask2.=1
        mask2[max(begin,end-3*neighborhood):end].=0
        mask = mask1.|mask2

        sum_deleted = sum(mat[row, idx[mask]])      # Compute total deleted sum
        mat[row, idx[mask]] .= 0.0                  # Set selected elements to 0
    end

    return mat
end

function delete_elements!(mat::Matrix{Float64}, threshold,neighborhood)
    # mat,threshold=abs.(A.*multiplier),error_
    m, n = size(mat)
    # neighborhood=10
    for row in 1:m
        sum_deleted = 0.0
        idx = sortperm(mat[row, :])  # Get sorted indices (ascending)
        if length(threshold)>1
            maxlevel=threshold[row]
        else
            maxlevel=threshold
        end
        # idx_rev = sortperm(mat[row, :],rev=:true)
        # Vectorized loop replacement
        # plot(cumsum(mat[row, idx]))
        # (mat[row, idx_rev])
        # plot(cumsum(mat[row, idx_rev]).>=3*threshold)
        # plot(cumsum(mat[row, idx]) .<=threshold)
        # plot(cumsum(mat[row, idx])./threshold)
        # Ne_index=similar(mask1)
        # cumsum(mat[row, idx])[end-30:end]
        mask1 = cumsum(mat[row, idx]) .<= maxlevel  # Identify elements to delete
        mask2=similar(mask1)
        mask2.=1
        mask2[max(begin,end-3*neighborhood):end].=0
        mask = mask1.|mask2

        sum_deleted = sum(mat[row, idx[mask]])      # Compute total deleted sum
        mat[row, idx[mask]] .= 0.0                  # Set selected elements to 0
    end

    return mat
end

function findTR_Q(V_0,Vmaxlimit0,Vminlimit0,S_VQ,QBound,raw,cl;Improve_Rate=5)

    z0=sum((V_0-Vmaxlimit0).>=0)+sum(Vminlimit0-V_0.>=0)
    zmax=round(size(Vmaxlimit0)[1]*size(Vmaxlimit0)[2]*size(Vmaxlimit0)[3]*Q_con_percent)+1
    # Vdown=(V_0-Vmaxlimit0)./ovr_Sens_VQ
    # Vup=(Vminlimit0-V_0)./ovr_Sens_VQ
    # TR=0.001
    # TR_max=0.05
    # TR_min=0
    # neighoborhood=50
    TR=0;TR_max=0.2;TR_min=0.0
    iter=0
    Itermax=10
    Vmaxinrisk=0
    Vmininrisk=0
    Vmaxviolated=0
    Vminviolated=0
    worst_vio=max(maximum(V_0-Vmaxlimit0),maximum(-(V_0-Vminlimit0)))
    # target_acc=(worst_vio/Improve_Rate)
    adj_acc=0
    acc=0
    best_acc=deepcopy(worst_vio)
    target_acc=0
    Vmaxlimit_new=deepcopy(Vmaxlimit0)
    # Vmaxlimit_new[:,:,2:end].+=target_acc
    # Vmaxlimit_new[:,Not(raw["Gen"].Bus),2:end].+=target_acc
    # Vminlimit_new=deepcopy(Vminlimit0)
    # Vminlimit_new[:,:,2:end].-=target_acc
    # Vminlimit_new[:,Not(raw["Gen"].Bus),2:end].-=target_acc
    # !(raw["Gen"].Bus)

    ovr_Sens_VQ_up = zeros(nTP, raw["nBus"], cl);ovr_Sens_VQ_dn = zeros(nTP, raw["nBus"], cl)
    ovr_Sens_VQ_up1 = zeros(nTP, raw["nBus"], cl);ovr_Sens_VQ_dn1 = zeros(nTP, raw["nBus"], cl)

    # SVQ_tc=zeros(raw["nQActBus"])
    # up_max = similar(SVQ_tc);dn_max = similar(SVQ_tc);up=similar(SVQ_tc);dn=similar(SVQ_tc)
    Qinj_up,Qinj_dn=0,0
    y=0
    z=0
    TR=0.01
    # target_acc0=worst_vio*adj_acc
    # Vmaxlimit_new0=Vmaxlimit0.+(target_acc0)
    # Vminlimit_new0=Vminlimit0.-(target_acc0)
    Vmaxlimit_new=Vmaxlimit0
    Vminlimit_new=Vminlimit0
    report_TR=DataFrame(worst_vio=0.0,acc=0.0,TR=0.0,z=z0,y=z0,z_acc=0.0,y_acc=0.0,plotdata=[0.0])
    deleteat!(report_TR,1)
    z_acc=0
    y_acc=0
    TR_best=0.0
    # while adj_acc<=0.9
    #     # target_acc=max(worst_vio*adj_acc,0.0001)
    #     target_acc=worst_vio*adj_acc
    #     Vmaxlimit_new=Vmaxlimit0.+(target_acc)
    #     Vminlimit_new=Vminlimit0.-(target_acc)
    #     Vmaxlimit_new[:,raw["Gen"].Bus,:].-=target_acc
    #     Vminlimit_new[:,raw["Gen"].Bus,:].+=target_acc
    #     iter=0
    #     TR_min=0
    #     iter=0
    #     if y!=0
    #         TR_max=0.5
    #     end
        while TR_max - TR_min > 1e-7  # Stop when the range is small enough
            # Update TR using bisection
            if iter==0
                nothing
            elseif iter>=Itermax
                TR=TR_best
            # elseif (z_acc/y_acc)>5 && y==0
            #     TR=min(TR_max,max(TR/(z_acc/y_acc),1.2 *TR_min))
            else
                TR = (TR_min + TR_max) / 2
            end

            # sum(QBound[1].<QBound[2])
            # sum(PBound[1].<PBound[2])
            Qinj_up=min.(TR,QBound[1])
            Qinj_dn=max.(-TR,QBound[2])
            for c in 1:cl, t in 1:nTP
                @views S_VQ_tc = S_VQ[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
                up=Qinj_up[t,:,c]
                dn=Qinj_dn[t,:,c]
                Threads.@threads for i in 1:raw["nBus"]
                # ThreadsX.foreach(1:raw["nLine"]) do i
                    SVQ_tc = S_VQ_tc[i, :]
                    # up_max = map(func_sens, SVQ_tc, Qinj_up, Qinj_dn)
                    # dn_max = map(func_sens, SVQ_tc, Qinj_up, Qinj_dn)
                    # @. up_max = SVQ_tc * ifelse(SVQ_tc > 0,up, dn)
                    # @. dn_max = SVQ_tc * ifelse(SVQ_tc > 0, dn, up)
                    up_max = SVQ_tc .* ifelse.(SVQ_tc .> 0, up, dn)
                    dn_max = SVQ_tc .* ifelse.(SVQ_tc .> 0, dn, up)

                    ovr_Sens_VQ_up[t, i, c] = sum(up_max) + 1e-6
                    ovr_Sens_VQ_dn[t, i, c] = sum(dn_max) + 1e-6

                    # ovr_Sens_VQ_up1[t, i, c] = sum(partialsort(up_max, rev=true, 1:neighoborhood)) + 1e-6
                    # ovr_Sens_VQ_dn1[t, i, c] = sum(partialsort(dn_max, 1:neighoborhood)) + 1e-6
                end
            end
            # Vmaxviolated= V_0.>=Vmaxlimit_new
            # Vminviolated= V_0.<=Vminlimit_new

            # vec(V_0+ovr_Sens_VQ_up.-Vmaxlimit_new)
            # vec(-(V_0+ovr_Sens_VQ_dn.-Vminlimit_new))
            # v = vcat( vec(V_0+ovr_Sens_VQ_up.-Vmaxlimit_new), vec(-(V_0+ovr_Sens_VQ_dn.-Vminlimit_new)))
            # plot!(sort(vcat( vec(V_0+ovr_Sens_VQ_up.-Vmaxlimit_new),vec(-(V_0+ovr_Sens_VQ_dn.-Vminlimit_new))),rev=true))
            z_acc=sort(vcat( vec(V_0+ovr_Sens_VQ_up/1 .-Vmaxlimit_new),
            vec(-(V_0+ovr_Sens_VQ_dn/1 .-Vminlimit_new))),rev=true)[min(end,Int(zmax))]
            # y_acc=max(maximum(V_0+(ovr_Sens_VQ_dn/Q_eff).-Vmaxlimit_new),maximum(-(V_0+(ovr_Sens_VQ_up/Q_eff).-Vminlimit_new)),0)
            y_acc=max(maximum(V_0+(ovr_Sens_VQ_dn/Q_eff).-Vmaxlimit_new),maximum(-(V_0+(ovr_Sens_VQ_up/Q_eff).-Vminlimit_new)))
            acc=max(z_acc,y_acc)
            Vmaxinrisk = V_0+ovr_Sens_VQ_up/1 .>= Vmaxlimit_new.+max(z_acc,-1e-4)
            Vmininrisk = V_0+ovr_Sens_VQ_dn/1 .<= Vminlimit_new.-max(z_acc,-1e-4)
            z= sum(Vmaxinrisk) + sum(Vmininrisk)
            y = sum(V_0+(ovr_Sens_VQ_dn/Q_eff) .> Vmaxlimit_new.+y_acc) + sum(V_0+(ovr_Sens_VQ_up/Q_eff) .<Vminlimit_new.-y_acc)
            push!(report_TR,[worst_vio,acc,TR,z,y,z_acc,y_acc,0.0])
            # plotdata=max.(0,sort(vcat( vec(V_0+ovr_Sens_VQ_up.-Vmaxlimit_new),vec(-(V_0+ovr_Sens_VQ_dn.-Vminlimit_new))),rev=true))
            # push!(report_TR,[worst_vio,acc,TR,z,y,z_acc,y_acc,plotdata],promote=true)
            if iter>=Itermax
                break
            end

            if acc<=best_acc
                best_acc=acc
                TR_best=TR
            end
            if z_acc>=y_acc|| z_acc>=0.04  #  || y_acc<=-0.01
                TR_max = TR  # Reduce TR
            else
                # if z>zmax
                #    break
                # end
                TR_min = TR  # Increase TR
            end
            iter+=1
        end
    #     Rem_vio=max(maximum(V_0+ovr_Sens_VQ_dn .-Vmaxlimit_new),maximum(-(V_0+ovr_Sens_VQ_up .-Vminlimit_new)))
    #     target_acc=worst_vio*adj_acc+max(Rem_vio,0)
    #     push!(report_TR,[target_acc,TR,z,y,adj_acc])
    #     if z>zmax
    #         adj_acc+=0.1
    #     else
    #         break
    #     end
    # end
    # (V_0+ovr_Sens_VQ_up/1)[21,3,5]
    VmaxIndex=findall(Vmaxinrisk)
    VminIndex=findall(Vmininrisk)
    display(report_TR)
    return TR,VmaxIndex,VminIndex,max(0.0,z_acc),[Qinj_up,Qinj_dn],worst_vio,report_TR
end
# knee_point(report_TR."worst_vio", report_TR."z")
# plot(report_TR."worst_vio", report_TR."z")

# knee_point((worst_vio.-report_TR."worst_vio")./worst_vio, report_TR."z")
# plot(worst_vio.-report_TR."worst_vio", report_TR."z")
# (worst_vio.-report_TR."worst_vio")./worst_vio
# ods_write("report_TR_Q", Dict(("data",1,1) => report_TR))

# Dict(report_TR)


# Function to estimate the knee point based on maximum curvature
function knee_point(x, y)
    #x=worst_vio.-report_TR."worst_vio"
    #y=report_TR."z"
    n = length(x)
    curvature = zeros(n)

    for i in 2:n-1
        x1, y1 = x[i-1], y[i-1]
        x2, y2 = x[i], y[i]
        x3, y3 = x[i+1], y[i+1]

        # Compute the curvature using the determinant method
        num = abs((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1))
        den = ((x2 - x1)^2 + (y2 - y1)^2)^1.5
        curvature[i] = num / (den + 1e-9)  # Adding a small term to avoid division by zero
    end

    # Find the index with maximum curvature
    knee_idx = argmax(curvature)
    return (x[knee_idx], y[knee_idx]), knee_idx
end

function findTR_P(Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,sens_IQ,PBound,raw,cl)

    Imaxlimit0_sqrt=sqrt.(Imaxlimit0)
    # z0=sum((I2_0+sens_IQ-Imaxlimit0).>=0)
    z0=0
    zmax=round(size(Imaxlimit0)[1]*size(Imaxlimit0)[2]*size(Imaxlimit0)[3]*P_con_percent)+1
    # Idown=(I2_0+sens_IQ-Imaxlimit0)./ovr_Sens_IP
    y=0
    z=0
    # TR=0.05
    # TR_max=0.5
    # TR_min=0
    # TR=0;TR_min=max(maximum(Base_TR_N),maximum(Base_TR_C))
    # TR_max=max(3*maximum(Base_TR_N),maximum(Base_TR_C),1);

    TR=0;TR_min=0.0
    TR_max=5;

    iter=0
    Iter_max=10
    atrisk=0
    violated=0
    worst_vio=maximum(I2_0-Imaxlimit0)
    worst_vio_sqrt=max.(maximum(sqrt.(I2_0)-Imaxlimit0_sqrt),maximum(Base_TR_N),maximum(Base_TR_C))
    best_acc=deepcopy(worst_vio_sqrt)
    acc=0
    # worst_vio_sqrt=maximum(sqrt.(I2_0)-Imaxlimit0_sqrt)

    # target_acc=(worst_vio/Improve_Rate)
    # if ramp_vio> worst_vio_sqrt
    #     Improve_Rate=1.5
    # end
    # target_acc_sqrt=(worst_vio_sqrt/Improve_Rate)
    adj_acc=0
    target_acc_sqrt=0
    Imaxlimit_new=deepcopy(Imaxlimit0)
    # neighoborhood=10

    # Imaxlimit_cache=((Imaxlimit0_sqrt.+target_acc_sqrt).^2)
    # Imaxlimit_new[:,:,2:end].=Imaxlimit_cache[:,:,2:end]

    ovr_Sens_IP_up = zeros(nTP, raw["nLine"], cl);ovr_Sens_IP_dn = zeros(nTP, raw["nLine"], cl)
    # SIP_tc=zeros(raw["nPActBus"])
    # up_max = similar(SIP_tc);dn_max = similar(SIP_tc);
    # up=similar(SIP_tc);dn=similar(SIP_tc)
    Pinj_up,Pinj_dn=0,0
    Gen_ind=indexin(raw["Gen"]."Bus", raw["PActBus"])
    # report_TR=DataFrame(worst_vio=0.0,TR=0.0,z=z0,y=z0,adj_acc=0.0)
    # deleteat!(report_TR,1)
    report_TR=DataFrame(worst_vio_sqrt=0.0,acc=0.0,TR=0.0,z=z0,y=z0,z_acc=0.0,y_acc=0.0,ramp_acc=0.0)
    deleteat!(report_TR,1)
    # TR=max(1,3*maximum(Base_TR_N),maximum(Base_TR_C))/2
    TR=1
    z_acc=0
    y_acc=0
    ramp_acc=0
    TR_best=0.0

    # max(maximum(Base_TR_N),maximum(Base_TR_C))
    # TR=5
    # while adj_acc<0.9
    #     target_acc_sqrt=worst_vio_sqrt*adj_acc
    #     Imaxlimit_new=(Imaxlimit0_sqrt.+target_acc_sqrt).^2
    #     # TR_max=0.5
    #     TR_min=min(max(maximum(Base_TR_N),maximum(Base_TR_C),0.01),TR_max-2e-4)
    #     iter=0
    # TR=0;TR_min=0
    # TR_max=0.001
        while TR_max - TR_min > 1e-7 # Stop when the range is small enough
            if iter==0
                nothing
            elseif iter>=Iter_max
                TR=TR_best
            elseif (z_acc)/max(y_acc,ramp_acc)>5 && y==0
                TR=max(TR/((z_acc)/max(y_acc)),1.2 *TR_min)
            else
                TR = (TR_min + TR_max) / 2
            end

            # TRmat=similar(PBound[1])
            # TRmat.=TR
            # TRmat[:,Gen_ind,:].+=Base_TR_N
            # Pinj_up=min.(max(0,TR - max(maximum(Base_TR_N),maximum(Base_TR_C))),PBound[1])
            Pinj_up=min.(TR,PBound[1])

            # Pinj_up[:,Gen_ind,:]=max.(Pinj_up[:,Gen_ind,:],2*Base_TR_N)
            # Pinj_dn=max.(min(0,-TR + max(maximum(Base_TR_N),maximum(Base_TR_C))),-TR,PBound[2])
            Pinj_dn=max.(-TR,PBound[2])

            # Pinj_dn[:,Gen_ind,:]=min.(Pinj_dn[:,Gen_ind,:],-2*Base_TR_N)
            for c in 1:cl, t in 1:nTP
                @views S_IP_tc = S_IP[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
                up=Pinj_up[t,:,c]
                dn=Pinj_dn[t,:,c]
                Threads.@threads for i in 1:raw["nLine"]
                # ThreadsX.foreach(1:raw["nLine"]) do i
                    SIP_tc = S_IP_tc[i, :]
                    # up_max = map(func_sens, SVQ_tc, Qinj_up, Qinj_dn)
                    # dn_max = map(func_sens, SVQ_tc, Qinj_up, Qinj_dn)
                    # @. up_max = SIP_tc .* ifelse(SIP_tc-1e-7 > 0, up, dn)
                    # @. dn_max = SIP_tc .* ifelse(SIP_tc-1e-7 > 0, dn, up)
                    up_max = SIP_tc .* ifelse.(SIP_tc .> 0, up, dn)
                    dn_max = SIP_tc .* ifelse.(SIP_tc .> 0, dn, up)

                    ovr_Sens_IP_up[t, i, c] = sum(up_max) + 1e-6
                    ovr_Sens_IP_dn[t, i, c] = sum(dn_max) + 1e-6
                    # ovr_Sens_IP_up[t, i, c] = sum(up_max) + 1e-6
                    # ovr_Sens_IP_dn[t, i, c] = sum(partialsort(dn_max, 1:neighoborhood)) + 1e-6

                end
            end

            # maximum(I2_0+ovr_Sens_IP_up.-Imaxlimit_new)
            # plot(sort(vec(sqrt.(I2_0+ovr_Sens_IP_up).-Imaxlimit0_sqrt),rev=true))
            # z_acc=sort(vec(sqrt.(max.(0,I2_0+ sens_IQ/5+ovr_Sens_IP_up/2)).-Imaxlimit0_sqrt),rev=true)[Int(zmax)]
            z_acc=sort(vec(sqrt.(max.(0,I2_0+ ovr_Sens_IP_up/1)).-Imaxlimit0_sqrt),rev=true)[min(end,Int(zmax))]

            # ramp_acc=max(maximum(Base_TR_N),maximum(Base_TR_C))-(TR)
            ramp_acc=0
            y_acc=max(maximum(sqrt.(max.(0,I2_0+ovr_Sens_IP_dn./P_eff)).-Imaxlimit0_sqrt),-0.1)
            acc=max(z_acc,y_acc)
            # Vmaxinrisk = V_0+ovr_Sens_VQ_up .>= Vmaxlimit_new.+z_acc
            # Vmininrisk = V_0+ovr_Sens_VQ_dn .<= Vminlimit_new.-z_acc
            # atrisk = I2_0+ovr_Sens_IP_up/2 .>= (Imaxlimit0_sqrt .+ max(z_acc,-1e-4)).^2
            atrisk = I2_0+ovr_Sens_IP_up/1 .>= (Imaxlimit0_sqrt .+ max(z_acc,-1e-4)).^2

            # violated= I2_0 .>= Imaxlimit_new
            # z= sum(atrisk)
            z = sum(atrisk)
            y = sum(I2_0+ovr_Sens_IP_dn./P_eff .>= (Imaxlimit0_sqrt .+ y_acc).^2)

            push!(report_TR,[worst_vio_sqrt,acc,TR,z,y,z_acc,y_acc,ramp_acc])
            if iter>=Iter_max
                break
            end

            if acc<=best_acc
                best_acc=acc
                TR_best=TR
            end
            if z_acc>=max(y_acc) || y_acc<=-0.1
                TR_max = TR  # Reduce TR
            else
                TR_min = TR  # Increase TR
            end
            iter+=1
        end
        # push!(report_TR,[target_acc_sqrt,TR,z,y,adj_acc])
    #     if z>zmax
    #         adj_acc+=0.1
    #     else
    #         break
    #     end
    # end
    # if y!=0
    #     target_acc_sqrt=maximum(sqrt.(max.(0,I2_0+(ovr_Sens_IP_dn)))-Imaxlimit0_sqrt)

    #     Imaxlimit_cache=((Imaxlimit0_sqrt.+(target_acc_sqrt/Improve_Rate)).^2)
    #     Imaxlimit_new[:,:,2:end].=Imaxlimit_cache[:,:,2:end]

    #     # Imaxlimit_new=(Imaxlimit0_sqrt.+target_acc_sqrt).^2
    #     atrisk= (I2_0+sens_IQ+ovr_Sens_IP_up .>= Imaxlimit_new)
    #     z = sum(atrisk)
    # end
    display(report_TR)
    ImaxIndex=findall(atrisk)
    return TR,ImaxIndex,max(0.0,z_acc),[Pinj_up,Pinj_dn],worst_vio,worst_vio_sqrt,report_TR
end
# knee_point(report_TR."worst_vio", report_TR."z")
# plot(report_TR."worst_vio", report_TR."z")

# knee_point((worst_vio.-report_TR."worst_vio")./worst_vio, report_TR."z")
# plot(worst_vio.-report_TR."worst_vio", report_TR."z")
# (worst_vio.-report_TR."worst_vio")./worst_vio
# ods_write("report_TR_P", Dict(("data",1,1) => report_TR))

function Seprate_Norm_Cont_indices(Indices)
    Norm_indices=similar(Indices)
    Cont_indices=similar(Indices)
    # set_= Indices[1]
    # ind=1
    for (ind,set_) in enumerate(Indices)
        Norm_indices[ind]=[x for x in set_ if x[3]==1]
        tuple_=[x for x in set_ if x[3]!=1]
        Cont_indices[ind]=map(x -> CartesianIndex(x[1], x[2], x[3] - 1), tuple_)
    end
    return Norm_indices, Cont_indices
end

function Worst_sens_IQ(Qinjlimit,S_IQ,raw,cl)
    # up_max = zeros(nTP, raw["nLine"], cl);
    up_max=zeros(raw["nQActBus"])
    Sens_IQ = zeros(nTP, raw["nLine"], cl)
    # SIQ_tc=zeros(raw["nQActBus"])
    for c in 1:cl, t in 1:nTP
        @views S_IQ_tc = S_IQ[t, c]
        up=Qinjlimit[1][t,:,c]
        dn=Qinjlimit[2][t,:,c]
        for i in 1:raw["nLine"]
            SIQ_tc = S_IQ_tc[i, :]
            @. up_max = ifelse(SIQ_tc > 0, SIQ_tc * up, SIQ_tc * dn)
            Sens_IQ[t, i, c] = sum(up_max) + 1e-6
        end
    end
    return Sens_IQ
end

function read_iter_rep(Problem_Type,Runcode)
        Result= ods_read(Pathplot*Problem_Type*string(Runcode),sheetName="General", retType= "DataFrame")
        return Result
end
function read_iter_rep_temp(Problem_Type,Runcode)
    # Runcode=rr
        Result= ods_read(Pathplot*Problem_Type*Runcode,sheetName="General", retType= "DataFrame")
        return Result
end

# Iter_Rep=read_iter_rep(Problem_Type,Runcode)
function plot_Iter_rep(Problem_Type,code;scale_="on")
    line1, line2, line5=0,0,0
    # code=code_
    datadict=Dict()
    for c in code[2]
        try
            raw=read_iter_rep(Problem_Type,(code[1],c))
            fi=findfirst(x->x=="end",raw[:,1])
            datadict[c]=coalesce.(raw[2:fi-1,:],0)
        catch
            continue
        end
    end
    titles=["Yang", "Hybrid", "Sensitivity"]

        # read_iter_rep(Problem_Type,code)
        # data=Iter_Rep[2:end-1,:]
        # time_method1=data[:,"Res_Time"]
        # max_V_violation_N=max.(data[:,"Vmax_Abs"],data[:,"Vmin_Abs"])
        # max_V_violation_C=max.(data[:,"Vmax_C_Abs"],data[:,"Vmin_C_Abs"])
        # Max_I_violation_N=data[:,"Smax_Abs"]
        # Max_I_violation_C=data[:,"Smax_C_Abs"]
        # opt_gap=data[:,"OptGap"]
    # Create the figure
    fig = Figure(resolution=1.5 .*(900, 600))

    # Prepare axis holders for each row
    axes_row1 = Axis[]  # Voltage Violation axes
    axes_row2 = Axis[]  # Current Violation axes
    axes_row3 = Axis[]  # Objective Gap axes
    axes_row4 = Axis[]  # ramp violation axes

    # c=5;i=2
    for (i, c) in enumerate(code[2])
        if i==1
            labley =["Voltage Violation (p.u.)","Current Violation (p.u.)","Ramp Violation (p.u.)","Objective Dev (%)"]
        else
            labley =["" "" "" ""]
        end
        if c in keys(datadict)
            data = datadict[c]
            time_method1 = data[:, "Res_Time"]
            max_V_violation_N = max.(data[:, "Vmax_Abs"], data[:, "Vmin_Abs"])
            max_V_violation_C = max.(data[:, "Vmax_C_Abs"], data[:, "Vmin_C_Abs"])
            Max_I_violation_N = data[:, "Smax_Abs"]
            Max_I_violation_C = data[:, "Smax_C_Abs"]
            Max_ramp_violation_N = data[:, "Ramp_Abs"]
            Max_ramp_violation_C = data[:, "Ramp_C_Abs"]
            opt_gap = data[:, "OptGap"]
        end

        # Primary axis
        ax = Axis(fig[1, i],
            ylabel = labley[1],
            title = titles[i],
            backgroundcolor = :white,
            xticklabelsize = 12, yticklabelsize = 12,
            xgridstyle = :dash, ygridstyle = :dash,
            xgridcolor = :gray, ygridcolor = :gray,
            ygridvisible = true)
        push!(axes_row1, ax)

        # Second row axis
        ax2 = Axis(fig[2, i],
            ylabel = labley[2],
            ygridvisible = true)
        push!(axes_row2, ax2)

        # Third row axis
        ax3 = Axis(fig[3, i],
            xlabel = "Time (seconds)",
            ylabel = labley[3],
            ygridvisible = true)
        push!(axes_row3, ax3)
        
        # Fourth row axis
        ax4 = Axis(fig[4, i],
            xlabel = "Time (seconds)",
            ylabel = labley[4],
            ygridvisible = true)
        push!(axes_row4, ax4)

        # Plot data
        if c in keys(datadict)
            for t in time_method1
                vlines!(ax, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
                vlines!(ax2, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
                vlines!(ax3, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
                vlines!(ax4, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
            end

            line1=scatterlines!(ax, time_method1, max_V_violation_N, label = "Max Violation (Method 1)", color = :red, linewidth = 2, markersize = 8, marker = :diamond)
            line2=scatterlines!(ax, time_method1, max_V_violation_C, label = "Max Violation (Method 1)", color = :orange, linewidth = 2, markersize = 8, linestyle = :dash, marker = :diamond)

            scatterlines!(ax2, time_method1, Max_I_violation_N, label = "Max Violation (Method 2)", color = :red, linewidth = 2, markersize = 8)
            scatterlines!(ax2, time_method1, Max_I_violation_C, label = "Max Violation (Method 2)", color = :orange, linewidth = 2, linestyle = :dash, markersize = 8)

            scatterlines!(ax3, time_method1, Max_ramp_violation_N, label = "Max Violation (Method 2)", color = :red, linewidth = 2, markersize = 8)
            scatterlines!(ax3, time_method1, Max_ramp_violation_C, label = "Max Violation (Method 2)", color = :orange, linewidth = 2, linestyle = :dash, markersize = 8)

            line5=scatterlines!(ax4, time_method1, opt_gap, label = "OptGap (Method 2)", color = :green, linewidth = 2, linestyle = :dash, markersize = 8)
        end
    end

    # # 🔗 Link axes for each row to use the same scale
    if scale_=="on"
        linkxaxes!(axes_row1...)  # Voltage Violation axes
        linkxaxes!(axes_row2...)  # Current Violation axes
        linkxaxes!(axes_row3...)  # Opt Gap axes
        linkxaxes!(axes_row4...)  # ramp violation axes
        linkyaxes!(axes_row1...)  # Voltage Violation axes
        linkyaxes!(axes_row2...)  # Current Violation axes
        linkyaxes!(axes_row3...)  # Opt Gap axes
        linkyaxes!(axes_row4...)  # ramp violation axes
    end
    fig
    Legend(fig[:, length(code[2]) +1],
    [line1, line2, line5],
    ["Normal Operation", "Contingency cases", "Obj Deviation"],
    patchsize = (35, 35), rowgap = 10)

    return fig
end
function plot_Iter_rep2(Problem_Type,code;scale_="on")
    line1, line2, line5=0,0,0
    # code=code_
    datadict=Dict()
    for c in code[2]
        try
            raw=read_iter_rep(Problem_Type,(code[1],c))
            fi=findfirst(x->x=="end",raw[:,1])
            datadict[c]=coalesce.(raw[2:fi-1,:],0)
        catch
            continue
        end
    end
    titles=["[17]", "LSB", "Hybrid LSB"]

        # read_iter_rep(Problem_Type,code)
        # data=Iter_Rep[2:end-1,:]
        # time_method1=data[:,"Res_Time"]
        # max_V_violation_N=max.(data[:,"Vmax_Abs"],data[:,"Vmin_Abs"])
        # max_V_violation_C=max.(data[:,"Vmax_C_Abs"],data[:,"Vmin_C_Abs"])
        # Max_I_violation_N=data[:,"Smax_Abs"]
        # Max_I_violation_C=data[:,"Smax_C_Abs"]
        # opt_gap=data[:,"OptGap"]
    # Create the figure
    fig = Figure(resolution=1 .*(900, 600))

    # Prepare axis holders for each row
    axes_row1 = Axis[]  # Voltage Violation axes
    axes_row2 = Axis[]  # Current Violation axes
    axes_row3 = Axis[]  # Objective Gap axes
    axes_row4 = Axis[]  # ramp violation axes

    # c=5;i=2
    for (i, c) in enumerate(code[2])
        if i==1
            labley =["Voltage Violation (p.u.)","Current Violation (p.u.)","Ramp Violation (p.u.)","Objective Dev (%)"]
        else
            labley =["" "" "" ""]
        end
        if c in keys(datadict)
            data = datadict[c]
            time_method1 = data[:, "Res_Time"]
            max_V_violation_N = max.(data[:, "Vmax_Abs"], data[:, "Vmin_Abs"])
            max_V_violation_C = max.(data[:, "Vmax_C_Abs"], data[:, "Vmin_C_Abs"])
            Max_I_violation_N = data[:, "Smax_Abs"]
            Max_I_violation_C = data[:, "Smax_C_Abs"]
            Max_ramp_violation_N = data[:, "Ramp_Abs"]
            Max_ramp_violation_C = data[:, "Ramp_C_Abs"]
            opt_gap = data[:, "OptGap"]
        end

        # Primary axis
        ax = Axis(fig[1, i],
            ylabel = labley[1],
            title = titles[i],
            # backgroundcolor = :white,
            # xticklabelsize = 12, yticklabelsize = 12,
            # xgridstyle = :dash, ygridstyle = :dash,
            # xgridcolor = :gray, ygridcolor = :gray,
            ygridvisible = true)
        push!(axes_row1, ax)

        # Second row axis
        ax2 = Axis(fig[2, i],
            xlabel = "Time (seconds)",
            ylabel = labley[2],
            ygridvisible = true)
        push!(axes_row2, ax2)

        # Third row axis
        # ax3 = Axis(fig[3, i],
        #     xlabel = "Time (seconds)",
        #     ylabel = labley[3],
        #     ygridvisible = true)
        # push!(axes_row3, ax3)
        
        # Fourth row axis
        # ax4 = Axis(fig[3, i],
        #     xlabel = "Time (seconds)",
        #     ylabel = labley[4],
        #     ygridvisible = true)
        # push!(axes_row4, ax4)

        # Plot data
        if c in keys(datadict)
            for t in time_method1
                vlines!(ax, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
                vlines!(ax2, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
                # vlines!(ax3, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
                # vlines!(ax4, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
            end

            line1=scatterlines!(ax, time_method1, max_V_violation_N, label = "Max Violation (Method 1)", color = :red, linewidth = 2, markersize = 8, marker = :diamond)
            line2=scatterlines!(ax, time_method1, max_V_violation_C, label = "Max Violation (Method 1)", color = :orange, linewidth = 2, markersize = 8, linestyle = :dash, marker = :diamond)

            scatterlines!(ax2, time_method1, Max_I_violation_N, label = "Max Violation (Method 2)", color = :red, linewidth = 2, markersize = 8)
            scatterlines!(ax2, time_method1, Max_I_violation_C, label = "Max Violation (Method 2)", color = :orange, linewidth = 2, linestyle = :dash, markersize = 8)

            # scatterlines!(ax3, time_method1, Max_ramp_violation_N, label = "Max Violation (Method 2)", color = :red, linewidth = 2, markersize = 8)
            # scatterlines!(ax3, time_method1, Max_ramp_violation_C, label = "Max Violation (Method 2)", color = :orange, linewidth = 2, linestyle = :dash, markersize = 8)

            # line5=scatterlines!(ax4, time_method1, opt_gap, label = "OptGap (Method 2)", color = :green, linewidth = 2, linestyle = :dash, markersize = 8)
        end
    end

    # # 🔗 Link axes for each row to use the same scale
    if scale_=="on"
        linkxaxes!(axes_row1...)  # Voltage Violation axes
        linkxaxes!(axes_row2...)  # Current Violation axes
        # linkxaxes!(axes_row3...)  # Opt Gap axes
        # linkxaxes!(axes_row4...)  # ramp violation axes
        linkyaxes!(axes_row1...)  # Voltage Violation axes
        linkyaxes!(axes_row2...)  # Current Violation axes
        # linkyaxes!(axes_row3...)  # Opt Gap axes
        # linkyaxes!(axes_row4...)  # ramp violation axes
    end
    fig
    # Legend(fig[:, length(code[2]) +1],
    # [line1, line2, line5],
    # ["Normal Operation", "Contingency cases", "Obj Deviation"],
    # patchsize = (35, 35), rowgap = 10)

    Legend(fig[3, :],
    [line1, line2],
    ["Normal Operation", "Contingency cases"],
    patchsize = (35, 35), rowgap = 10,orientation = :horizontal)

    return fig
end

function findIcoupled(TR_Q,I2_0,Imax_IP,ImaxIndex,Imaxlimit0,S_IQ,QBound,raw,cl,acc_I=-1e-4)

    Imaxlimit0_sqrt=sqrt.(Imaxlimit0)
    zmax=round(size(Imaxlimit0)[1]*size(Imaxlimit0)[2]*size(Imaxlimit0)[3]*P_con_percent_coupled)+1
    atrisk=0
    acc=0
    adj_acc=0
    target_acc_sqrt=0
    Imaxlimit_new=deepcopy(Imaxlimit0)
    ovr_Sens_IQ_up = zeros(nTP, raw["nLine"], cl);ovr_Sens_IQ_dn = zeros(nTP, raw["nLine"], cl)
    Qinj_up,Qinj_dn=0,0
    TR=1
    z_acc=0
    y_acc=0
    ramp_acc=0
    TR_best=0.0

    Qinj_up=min.(TR_Q,QBound[1])
    Qinj_dn=max.(-TR_Q,QBound[2])

    for c in 1:cl, t in 1:nTP
        @views S_IQ_tc = S_IQ[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
        up=Qinj_up[t,:,c]
        dn=Qinj_dn[t,:,c]
        Threads.@threads for i in 1:raw["nLine"]
        # ThreadsX.foreach(1:raw["nLine"]) do i
            SIQ_tc = S_IQ_tc[i, :]
            up_max = SIQ_tc .* ifelse.(SIQ_tc .> 0, up, dn)
            dn_max = SIQ_tc .* ifelse.(SIQ_tc .> 0, dn, up)

            ovr_Sens_IQ_up[t, i, c] = sum(up_max) + 1e-6
            ovr_Sens_IQ_dn[t, i, c] = sum(dn_max) + 1e-6
        end
    end

    I_max_IQ=sqrt.(max.(0,Imax_IP+ ovr_Sens_IQ_up/1))
    I_max_IQ[ImaxIndex].=0
    z_acc=sort(vec(I_max_IQ.-Imaxlimit0_sqrt),rev=true)[min(end,Int(zmax))]

    # z_acc=sort(vec(I_max_IQ.-Imaxlimit0_sqrt),rev=true)[min(end,Int(zmax))]

    # atrisk[14,2,5]
    # (Imax_IP+ ovr_Sens_IQ_up/1)[14,2,5]
    # acc=sort(vec(sqrt.(max.(0,I2_0+ ovr_Sens_IQ_up/1)).-Imaxlimit0_sqrt),rev=true)[min(end,Int(zmax))]

    # Plots.plot(sort(vec(sqrt.(max.(0,I2_0+ ovr_Sens_IQ_up/1)).-Imaxlimit0_sqrt),rev=true))
    # Vmaxinrisk = V_0+ovr_Sens_VQ_up .>= Vmaxlimit_new.+z_acc
    # Vmininrisk = V_0+ovr_Sens_VQ_dn .<= Vminlimit_new.-z_acc
    # atrisk = I2_0+ovr_Sens_IP_up/2 .>= (Imaxlimit0_sqrt .+ max(z_acc,-1e-4)).^2
    # atrisk = I2_0+ovr_Sens_IQ_up/1 .>= (Imaxlimit0_sqrt .+ max(z_acc,acc_I)).^2
    atrisk = I_max_IQ .>= (Imaxlimit0_sqrt .+ max(z_acc,acc_I))

    # ((Imaxlimit0_sqrt .+ max(z_acc,acc_I)).^2)[14,2,5]
    # violated= I2_0 .>= Imaxlimit_new
    # z= sum(atrisk)
    # z = sum(atrisk)
    ImaxIndex_Q= findall(atrisk)
    println("Number of lines in risk due to Q adjustment: ", length(ImaxIndex_Q))
    println("line acc due to Q adjustment: ", z_acc)
    return ImaxIndex_Q
end

function findVcoupled(TR_P,V_0,Vmaxlimit0,Vminlimit0,S_VP,PBound,raw,cl,acc_V=1e-4)

    zmax=round(size(Vmaxlimit0)[1]*size(Vmaxlimit0)[2]*size(Vmaxlimit0)[3]*Q_con_percent_coupled)+1

    Vmaxinrisk=0
    Vmininrisk=0
    acc=0
    ovr_Sens_VP_up = zeros(nTP, raw["nBus"], cl);ovr_Sens_VP_dn = zeros(nTP, raw["nBus"], cl)
    Pinj_up,Pinj_dn=0,0
    z_acc=0

    Pinj_up=min.(TR_P, PBound[1])
    Pinj_dn=max.(-TR_P, PBound[2])

    for c in 1:cl, t in 1:nTP
        @views S_VP_tc = S_VP[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
        up=Pinj_up[t,:,c]
        dn=Pinj_dn[t,:,c]
        for i in 1:raw["nBus"]
        # ThreadsX.foreach(1:raw["nLine"]) do i
            SVP_tc = S_VP_tc[i, :]
            up_max = SVP_tc .* ifelse.(SVP_tc .> 0, up, dn)
            dn_max = SVP_tc .* ifelse.(SVP_tc .> 0, dn, up)

            ovr_Sens_VP_up[t, i, c] = sum(up_max) + 1e-6
            ovr_Sens_VP_dn[t, i, c] = sum(dn_max) + 1e-6
        end
    end
    acc=sort(vcat( vec(V_0+ovr_Sens_VP_up/1 .-Vmaxlimit0),
    vec(-(V_0+ovr_Sens_VP_dn/1 .-Vminlimit0))),rev=true)[min(end,Int(zmax))]

    Vmaxinrisk = V_0+ovr_Sens_VP_up .>= Vmaxlimit0.+max(acc,acc_V,1e-4)
    Vmininrisk = V_0+ovr_Sens_VP_dn .<= Vminlimit0.-max(acc,acc_V,1e-4)
    # atrisk = I2_0+ovr_Sens_IP_up/2 .>= (Imaxlimit0_sqrt .+ max(z_acc,-1e-4)).^2
    # atrisk = I2_0+ovr_Sens_IQ_up/1 .>= (Imaxlimit0_sqrt .+ max(z_acc,-1e-4)).^2

    # violated= I2_0 .>= Imaxlimit_new
    # z= sum(atrisk)
    # z = sum(atrisk)
    Vmaxindex_P=findall(Vmaxinrisk)
    Vminindex_P=findall(Vmininrisk)
    println("Number of voltage in risk due to P adjustment: ", length(Vmaxindex_P)+length(Vminindex_P))
    println("voltage acc due to P adjustment: ", acc)
    return Vmaxindex_P,Vminindex_P
end

function findTR_P_coupled(Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,V_0,Vmaxlimit0,Vminlimit0,S_VP,PBound,raw,cl,acc_V0)

    # acc_V0=0.05
    # S_VP=keep_top_bottom!.(S_VP, 0.02)
    # S_VP=keep_top_abs!.(S_VP, 0.05)
    # S_VP=keep_K_percent!.(S_VP,0.5)
    # S_VP=keep_K_percent_twoside!.(S_VP,0.7)

    # keep_K_percent_twoside!

    # S_VP[1,1][:,1][51:60]
    # @profview begin
    Imaxlimit0_sqrt=sqrt.(Imaxlimit0)
    # z0=sum((I2_0+sens_IQ-Imaxlimit0).>=0)
    z0=0
    zmax=round(size(Imaxlimit0)[1]*size(Imaxlimit0)[2]*size(Imaxlimit0)[3]*P_con_percent)+1
    zmax_V=round(size(Vmaxlimit0)[1]*size(Vmaxlimit0)[2]*size(Vmaxlimit0)[3]*Q_con_percent_coupled)+1
    # Idown=(I2_0+sens_IQ-Imaxlimit0)./ovr_Sens_IP
    y=0
    z=0
    # TR=0.05
    # TR_max=0.5
    # TR_min=0
    # TR=0;TR_min=max(maximum(Base_TR_N),maximum(Base_TR_C))
    # TR_max=max(3*maximum(Base_TR_N),maximum(Base_TR_C),1);

    TR=0;TR_min=0.0
    TR_max=5;

    acc_V=min(acc_V0,max(maximum(V_0-Vmaxlimit0),maximum(-(V_0-Vminlimit0))))
    # acc_V=max(maximum(V_0-Vmaxlimit0),maximum(-(V_0-Vminlimit0)))


    iter=0
    Iter_max=10
    atrisk=0
    Vmaxinrisk=0;Vmininrisk=0;
    Vmax_VP=0;Vmin_VP=0;Imax_IP=0;
    violated=0
    worst_vio=maximum(I2_0-Imaxlimit0)
    worst_vio_sqrt=max.(maximum(sqrt.(I2_0)-Imaxlimit0_sqrt),maximum(Base_TR_N),maximum(Base_TR_C))
    best_acc=deepcopy(worst_vio_sqrt)
    acc=0
    z_acc_VP=0
    # worst_vio_sqrt=maximum(sqrt.(I2_0)-Imaxlimit0_sqrt)

    # target_acc=(worst_vio/Improve_Rate)
    # if ramp_vio> worst_vio_sqrt
    #     Improve_Rate=1.5
    # end
    # target_acc_sqrt=(worst_vio_sqrt/Improve_Rate)
    adj_acc=0
    target_acc_sqrt=0
    Imaxlimit_new=deepcopy(Imaxlimit0)
    genbus=raw["Gen"].Bus
    genbus_s = Set(genbus)  # O(1) lookup
    allbus_s= Set(1:raw["nBus"])
    # neighoborhood=10

    # Imaxlimit_cache=((Imaxlimit0_sqrt.+target_acc_sqrt).^2)
    # Imaxlimit_new[:,:,2:end].=Imaxlimit_cache[:,:,2:end]

    ovr_Sens_IP_up = zeros(nTP, raw["nLine"], cl);ovr_Sens_IP_dn = zeros(nTP, raw["nLine"], cl)
    ovr_Sens_VP_up = zeros(nTP, raw["nBus"], cl);ovr_Sens_VP_dn = zeros(nTP, raw["nBus"], cl)

    # ovr_Sens_VP_upmat = [zeros(raw["nBus"], raw["nPActBus"]) for i in 1:nTP, j in 1:cl]  
    # ovr_Sens_VP_dnmat = [zeros(raw["nBus"], raw["nPActBus"]) for i in 1:nTP, j in 1:cl]    
    # VP_up = [zeros(raw["nBus"], raw["nPActBus"]) for i in 1:nTP, j in 1:cl]  
    # VP_dn = [zeros(raw["nBus"], raw["nPActBus"]) for i in 1:nTP, j in 1:cl] 
    # upelem = [zeros(raw["nBus"], raw["nPActBus"]) for i in 1:nTP, j in 1:cl]  
    # dnelem = [zeros(raw["nBus"], raw["nPActBus"]) for i in 1:nTP, j in 1:cl]       

        # Dimensions
    nBus      = raw["nBus"]
    nPActBus  = raw["nPActBus"]

    # Allocate once, contiguous in memory
    ovr_Sens_VP_upmat = zeros(Float64, nBus, nPActBus, nTP, cl)
    ovr_Sens_VP_dnmat = zeros(Float64, nBus, nPActBus, nTP, cl)
    VP_up             = similar(ovr_Sens_VP_upmat)
    VP_dn             = similar(ovr_Sens_VP_upmat)
    upelem            = similar(ovr_Sens_VP_upmat)
    dnelem            = similar(ovr_Sens_VP_upmat)
    VP_up1      =falses(nBus, nPActBus, nTP, cl)
    VP_dn1      =similar(VP_up1)
    # VP_up3      =similar(VP_up1)
    # VP_dn3      =similar(VP_up1)
    # VP_up4      =similar(VP_up1)
    # VP_dn4      =similar(VP_up1)


    DV_up=V_0.- (acc_V .+ Vmaxlimit0)
    # DV_up=V_0.- (0 .+ Vmaxlimit0)

    # DV_up_4D=reshape(DV_up, size(DV_up,2),1,size(DV_up,1),size(DV_up,3))

    DV_dn=(V_0.- (-acc_V .+ Vminlimit0))
    # DV_dn=(V_0.- (-0 .+ Vminlimit0))

    # DV_dn_4D=reshape(DV_dn, size(DV_dn,2),1,size(DV_dn,1),size(DV_dn,3))
    DV_up_4D = zeros(Float64, nBus, 1, nTP, cl)
    DV_dn_4D = zeros(Float64, nBus, 1, nTP, cl)
    S_VP_4D = zeros(Float64, nBus, nPActBus, nTP, cl)
    Threads.@threads for t in 1:nTP
        for c in 1:cl
            @views S_VP_4D[:, :, t, c] .= S_VP[t, c]
            @views DV_up_4D[:, 1, t, c] .= DV_up[t, :, c]  # direct assignment
            @views DV_dn_4D[:, 1, t, c] .= DV_dn[t, :, c]  # direct assignment


        end
    end
    S_VP1=keep_K_percent_twoside_4D(S_VP_4D,0.99)


    # SIP_tc=zeros(raw["nPActBus"])
    # up_max = similar(SIP_tc);dn_max = similar(SIP_tc);
    # up=similar(SIP_tc);dn=similar(SIP_tc)
    Pinj_up,Pinj_dn=0,0
    Gen_ind=indexin(raw["Gen"]."Bus", raw["PActBus"])
    # report_TR=DataFrame(worst_vio=0.0,TR=0.0,z=z0,y=z0,adj_acc=0.0)
    # deleteat!(report_TR,1)
    report_TR=DataFrame(worst_vio_sqrt=0.0,acc=0.0,TR=0.0,z=z0,y=z0,z_acc=0.0,y_acc=0.0,ramp_acc=0.0,z_V=0,acc_VP=0.0)
    deleteat!(report_TR,1)
    # TR=max(1,3*maximum(Base_TR_N),maximum(Base_TR_C))/2
    TR=1
    z_acc=0
    y_acc=0
    ramp_acc=0
    TR_best=0.0

    # Report1=deepcopy(report_TR)
    # max(maximum(Base_TR_N),maximum(Base_TR_C))
    # TR=5
    # while adj_acc<0.9
    #     target_acc_sqrt=worst_vio_sqrt*adj_acc
    #     Imaxlimit_new=(Imaxlimit0_sqrt.+target_acc_sqrt).^2
    #     # TR_max=0.5
    #     TR_min=min(max(maximum(Base_TR_N),maximum(Base_TR_C),0.01),TR_max-2e-4)
    #     iter=0
    # TR=0;TR_min=0
    # TR_max=0.001
    # TR=1
    # Iter_max=0
        while TR_max - TR_min > 1e-7# Stop when the range is small enough
            if iter==0
                nothing
            elseif iter>=Iter_max
                TR=TR_best
            elseif (z_acc)/max(y_acc,ramp_acc)>5 && y==0
                TR=max(TR/((z_acc)/max(y_acc)),1.2 *TR_min)
            else
                TR = (TR_min + TR_max) / 2
            end

            VP_up1.=0
            VP_dn1.=0
            # TRmat=similar(PBound[1])
            # TRmat.=TR
            # TRmat[:,Gen_ind,:].+=Base_TR_N
            # Pinj_up=min.(max(0,TR - max(maximum(Base_TR_N),maximum(Base_TR_C))),PBound[1])
            Pinj_up=min.(TR,PBound[1])

            # Pinj_up[:,Gen_ind,:]=max.(Pinj_up[:,Gen_ind,:],2*Base_TR_N)
            # Pinj_dn=max.(min(0,-TR + max(maximum(Base_TR_N),maximum(Base_TR_C))),-TR,PBound[2])
            Pinj_dn=max.(-TR,PBound[2])

            # sum(Pinj_dn.>0.01)

            # Pinj_dn[:,Gen_ind,:]=min.(Pinj_dn[:,Gen_ind,:],-2*Base_TR_N)
            Threads.@threads for t in 1:nTP
                for c in 1:cl
                    @views S_IP_tc = S_IP[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
                    @views S_VP_tc = S_VP[t,c]
                    @views up = Pinj_up[t, :, c]  # 1 × nPActBus
                    @views dn = Pinj_dn[t, :, c]
                    
                    mask = S_IP_tc .> 0  # nLine × nPActBus
                    up_max = S_IP_tc .* (mask .* up' .+ .!mask .* dn')
                    dn_max = S_IP_tc .* (mask .* dn' .+ .!mask .* up')
                    
                    @views ovr_Sens_IP_up[t, :, c] .= sum(up_max, dims=2)[1:raw["nLine"]] .+ 1e-6
                    @views ovr_Sens_IP_dn[t, :, c] .= sum(dn_max, dims=2)[1:raw["nLine"]] .+ 1e-6

                    mask_V = S_VP_tc .> 0  # nLine × nPActBus
                    up_max_V = S_VP_tc .* (mask_V .* up' .+ .!mask_V .* dn')
                    dn_max_V = S_VP_tc .* (mask_V .* dn' .+ .!mask_V .* up')

                    @views ovr_Sens_VP_up[t, :, c] .= sum(up_max_V, dims=2) .+ 1e-6
                    @views ovr_Sens_VP_dn[t, :, c] .= sum(dn_max_V, dims=2) .+ 1e-6

                    @views ovr_Sens_VP_upmat[:, :, t, c] .= up_max_V
                    @views ovr_Sens_VP_dnmat[:, :, t, c] .= dn_max_V
                end
            end


            
            # keep_top!(ovr_Sens_VP_upmat[1,5].+( V_0[1,:,5].-Vmaxlimit0[1,:,5]),0.2).!=0 *

            # Threads.@threads for i in axes(ovr_Sens_VP_upmat,1)
            #     for j in axes(ovr_Sens_VP_upmat,2)
            #     VP_up[i,j]= ovr_Sens_VP_upmat[i,j] .+ (V_0[i,:,j] .- (acc_V.+Vmaxlimit0[i,:,j]))
            #     VP_dn[i,j]= -ovr_Sens_VP_dnmat[i,j] .- (V_0[i,:,j] .- (-acc_V.+Vminlimit0[i,:,j]))
            #     end
            # end

            @views VP_up .= ovr_Sens_VP_upmat .+ DV_up_4D
            @views VP_dn .= -ovr_Sens_VP_dnmat .- DV_dn_4D


            # VP_up= [ovr_Sens_VP_upmat[i,j] .+ (V_0[i,:,j] .- (acc_V.+Vmaxlimit0[i,:,j]))  for i in axes(ovr_Sens_VP_upmat,1), j in axes(ovr_Sens_VP_upmat,2)]
            # VP_dn= [-ovr_Sens_VP_dnmat[i,j] .- (V_0[i,:,j] .- (-acc_V.+Vminlimit0[i,:,j])) for i in axes(ovr_Sens_VP_dnmat,1), j in axes(ovr_Sens_VP_dnmat,2)]

            # VP_up3=deepcopy(VP_up)
            # VP_dn3=deepcopy(VP_dn)
            # VP_up4=deepcopy(VP_up)
            # VP_dn4=deepcopy(VP_dn)
            # VP_up3= [(V_0[i,:,j] .- (acc_V.+Vmaxlimit0[i,:,j])) ./ -ovr_Sens_VP_upmat[i,j]  for i in axes(ovr_Sens_VP_upmat,1), j in axes(ovr_Sens_VP_upmat,2)]
            # VP_dn3= [(V_0[i,:,j] .- (-acc_V.+Vminlimit0[i,:,j])) ./ (ovr_Sens_VP_dnmat[i,j]) for i in axes(ovr_Sens_VP_dnmat,1), j in axes(ovr_Sens_VP_dnmat,2)]

            # S_VP1=deepcopy(S_VP)
            # GC.gc()

            # S_VP1[1,1]
            # VP_up1= [ovr_Sens_VP_upmat[i,j].+1e-6 for i in axes(ovr_Sens_VP_upmat,1), j in axes(ovr_Sens_VP_upmat,2)]
            # VP_dn1= [-ovr_Sens_VP_dnmat[i,j].+1e-6 for i in axes(ovr_Sens_VP_dnmat,1), j in axes(ovr_Sens_VP_dnmat,2)]
       
            # keep_K_percent_new!.(VP_up,0.5)
            # keep_K_percent_new!.(VP_dn,0.5)
            @views VP_up1=keep_K_percent_new_4D!(VP_up1,VP_up,0.99)
            @views VP_dn1=keep_K_percent_new_4D!(VP_dn1,VP_dn,0.99)

            # @views VP_up1=keep_top_4D!(VP_up1,VP_up,3)
            # @views VP_dn1=keep_top_4D!(VP_dn1,VP_dn,3)
            # VP_up4=keep_top_idx.(VP_up,3)
            # VP_dn4=keep_top_idx.(VP_dn,3)

            # @views VP_up1=keep_top_4D!(VP_up1,VP_up,5)
            # @views VP_dn1=keep_top_4D!(VP_dn1,VP_dn,5)

            # VP_up4=keep_top_idx.(ovr_Sens_VP_upmat,1)
            # VP_dn4=keep_top_idx.(-ovr_Sens_VP_dnmat,1)

            # @views VP_up1 .=keep_top_4D!(VP_up1,ovr_Sens_VP_upmat,1)
            # @views VP_dn1 .=keep_top_4D!(VP_dn1,-ovr_Sens_VP_dnmat,1)

            # # # ovr_Sens_VP_upmat[:,:,1,1]

            # # # VP_up3=keep_top_row_gen_idx.(VP_up,Ref(genbus_s),1)
            # # # VP_dn3=keep_top_row_gen_idx.(VP_dn,Ref(genbus_s),1)
            # @views VP_up1 .=keep_top_row_gen_4D!(VP_up1,VP_up,genbus_s,1)
            # @views VP_dn1 .=keep_top_row_gen_4D!(VP_dn1,VP_dn,genbus_s,1)
            # @views VP_up1 .=keep_top_row_gen_4D!(VP_up1,VP_up,allbus_s,1)
            # @views VP_dn1 .=keep_top_row_gen_4D!(VP_dn1,VP_dn,allbus_s,1)


            # aa=(VP_up[20,15].!=0)

            # sum(aa,dims=1)

            # keep_top_bottom_row!.(VP_up1,0.01)
            # keep_top_bottom_row!.(VP_dn1,0.01)
            # keep_K_percent_new!.(VP_up1,0.9)
            # keep_K_percent_new!.(VP_dn1,0.9)
            # keep_K_percent_new!.(VP_up3,0.9)
            # keep_K_percent_new!.(VP_dn3,0.9)
            # keep_K_percent_new!.(VP_up1,0.75)
            # keep_K_percent_new!.(VP_dn1,0.75)



            # VP_up[i,j][57,21]
            # sum(ovr_Sens_VP_upmat[i,j][41,:])
            # findmax(ovr_Sens_VP_upmat[i,j][58,:])
            # ovr_Sens_VP_up[4,58,17]
            # V_0[11,41,16]
            # ovr_Sens_VP_up[7,58,13]

            # findmax(ovr_Sens_VP_upmat[i,j][:,21])

            # Plots.plot(ovr_Sens_VP_upmat[i,j][:,2])
            # Plots.plot(VP_up1[i,j][:,8])
            # findmax(VP_up[i,j][:,21])
            # Plots.plot(S_VP1[i,j][:,8])

            # # raw["PActBus"][21]
            # upelem[12,17][48,21]
            @. upelem=(VP_up1||S_VP1 ) *ovr_Sens_VP_upmat #||S_VP1 
            @. dnelem=(VP_dn1||S_VP1 ) *ovr_Sens_VP_dnmat #||S_VP1
            # ovr_Sens_VP_up.= reshape(sum(upelem,dims=2),nTP,raw["nBus"],cl)
            # ovr_Sens_VP_dn.= reshape(sum(dnelem,dims=2),nTP,raw["nBus"],cl)
            ovr_Sens_VP_up.= permutedims(sum(upelem,dims=2)[:,1,:,:], (2, 1, 3)) 
            ovr_Sens_VP_dn.= permutedims(sum(dnelem,dims=2)[:,1,:,:], (2, 1, 3)) 


            # typeof(VP_up3)

            # Threads.@threads for i in axes(ovr_Sens_VP_upmat,1) 
            #     for j in axes(ovr_Sens_VP_upmat,2)
            #         upelem[i,j]= ((VP_up4[i,j].!=0) .|| (VP_up3[i,j].!=0) .|| (VP_up1[i,j].!=0) .|| (S_VP1[i,j].!=0)) .* ovr_Sens_VP_upmat[i,j]
            #         # upelem[i,j][genbus,:]=ovr_Sens_VP_upmat[i,j][genbus,:]
            #         dnelem[i,j]= ((VP_dn4[i,j].!=0) .|| (VP_dn3[i,j].!=0) .|| (VP_dn1[i,j].!=0) .|| (S_VP1[i,j].!=0)) .* ovr_Sens_VP_dnmat[i,j]
            #         # dnelem[i,j][genbus,:]=ovr_Sens_VP_dnmat[i,j][genbus,:]
            #         upelemsum=sum(upelem[i,j],dims=2)
            #         dnelemsum=sum(dnelem[i,j],dims=2)
            #         ovr_Sens_VP_up[i, :, j].= upelemsum
            #         ovr_Sens_VP_dn[i, :, j].= dnelemsum
            #     end
            # end

            # upelem= [((VP_up3[i,j].!=0) .|| (VP_up[i,j].!=0) .|| (S_VP1[i,j].!=0)) .* ovr_Sens_VP_upmat[i,j] for i in axes(ovr_Sens_VP_upmat,1), j in axes(ovr_Sens_VP_upmat,2)]
            # dnelem= [((VP_dn3[i,j].!=0) .|| (VP_dn[i,j].!=0) .|| (S_VP1[i,j].!=0)) .* ovr_Sens_VP_dnmat[i,j] for i in axes(ovr_Sens_VP_upmat,1), j in axes(ovr_Sens_VP_upmat,2)]

            # upelemsum= sum.(upelem,dims=2)
            # dnelemsum= sum.(dnelem,dims=2)

            # i=1;j=7
            # ovr_Sens_VP_up[t, i, c]=sum(ovr_Sens_VP_upmat[t,c],dims=2)[i]
            #             b[1,1]

            # for t in axes(ovr_Sens_VP_upmat,1), c in axes(ovr_Sens_VP_upmat,2)
            #     ovr_Sens_VP_up[t, :, c].=upelemsum[t,c]
            #     ovr_Sens_VP_dn[t, :, c].=dnelemsum[t,c]
            # end
            # findall(Vmaxinrisk)
            # Vmaxinrisk[1,:,7]


            # V_0[1,:,1].-Vmaxlimit0[1,:,1]

            # keep_K_percent!(ovr_Sens_VP_upmat[1,1],0.7)

            vals_IP = vec(sqrt.(max.(0, I2_0 .+ ovr_Sens_IP_up)) .- Imaxlimit0_sqrt)
            k_IP = min(length(vals_IP), Int(zmax))
            z_acc = partialsort(vals_IP, 1:k_IP; rev=true)[end]

            vals_VP = vcat(vec(V_0 .+ ovr_Sens_VP_up .- Vmaxlimit0),
               vec(-(V_0 .+ ovr_Sens_VP_dn .- Vminlimit0)))
            k_VP = min(length(vals_VP), Int(zmax_V))
            z_acc_VP = partialsort(vals_VP, 1:k_VP; rev=true)[end]

            # Plots.plot(filter(x-> x>=0,sort(vcat( vec(V_0+ovr_Sens_VP_up/1 .-Vmaxlimit0),
            # vec(-(V_0+ovr_Sens_VP_dn/1 .-Vminlimit0))),rev=true)))

            Vmax_VP=V_0+ovr_Sens_VP_up
            Vmin_VP=V_0+ovr_Sens_VP_dn
            Vmaxinrisk = Vmax_VP .>= Vmaxlimit0.+max(z_acc_VP,1e-4)
            Vmininrisk = Vmin_VP .<= Vminlimit0.-max(z_acc_VP,1e-4)

            # ramp_acc=max(maximum(Base_TR_N),maximum(Base_TR_C))-(TR)
            ramp_acc=0
            y_acc=max(maximum(sqrt.(max.(0,I2_0+ovr_Sens_IP_dn./P_eff)).-Imaxlimit0_sqrt),-0.1)
            acc=max(z_acc,y_acc)
            # Vmaxinrisk = V_0+ovr_Sens_VQ_up .>= Vmaxlimit_new.+z_acc
            # Vmininrisk = V_0+ovr_Sens_VQ_dn .<= Vminlimit_new.-z_acc
            # atrisk = I2_0+ovr_Sens_IP_up/2 .>= (Imaxlimit0_sqrt .+ max(z_acc,-1e-4)).^2
            Imax_IP = I2_0+ovr_Sens_IP_up/1
            atrisk = Imax_IP .>= (Imaxlimit0_sqrt .+ max(z_acc,1e-4)).^2

            # violated= I2_0 .>= Imaxlimit_new
            # z= sum(atrisk)
            z = sum(atrisk)
            z_V=sum(Vmaxinrisk)+sum(Vmininrisk)
            y = sum(I2_0+ovr_Sens_IP_dn./P_eff .>= (Imaxlimit0_sqrt .+ y_acc).^2)

            push!(report_TR,[worst_vio_sqrt,acc,TR,z,y,z_acc,y_acc,ramp_acc,z_V,z_acc_VP])
            if iter>=Iter_max
                break
            end

            if acc<=best_acc && z_acc_VP <= acc_V
                best_acc=max(0,acc)
                TR_best=TR
            end
            if z_acc>=max(y_acc) || y_acc<=0 || z_acc_VP >= acc_V
                TR_max = TR  # Reduce TR
            else
                TR_min = TR  # Increase TR
            end
            iter+=1
        end
        # push!(report_TR,[target_acc_sqrt,TR,z,y,adj_acc])
    #     if z>zmax
    #         adj_acc+=0.1
    #     else
    #         break
    #     end
    # end
    # if y!=0
    #     target_acc_sqrt=maximum(sqrt.(max.(0,I2_0+(ovr_Sens_IP_dn)))-Imaxlimit0_sqrt)

    #     Imaxlimit_cache=((Imaxlimit0_sqrt.+(target_acc_sqrt/Improve_Rate)).^2)
    #     Imaxlimit_new[:,:,2:end].=Imaxlimit_cache[:,:,2:end]

    #     # Imaxlimit_new=(Imaxlimit0_sqrt.+target_acc_sqrt).^2
    #     atrisk= (I2_0+sens_IQ+ovr_Sens_IP_up .>= Imaxlimit_new)
    #     z = sum(atrisk)
    # end
    # TR=0.1
    display(report_TR)
    Vmaxindex_P = findall(Vmaxinrisk)
    Vminindex_P = findall(Vmininrisk)
    println("Number of voltage in risk due to P adjustment: ", length(Vmaxindex_P)+length(Vminindex_P))
    println("voltage acc left out due to P adjustment: ", z_acc_VP)
    ImaxIndex = findall(atrisk)
    # ImaxIndex = findall(.~isnan.(I2_0));    z_acc=0;

    return TR,ImaxIndex,Vmaxindex_P,Vminindex_P,max(0.0,z_acc),[Pinj_up,Pinj_dn],worst_vio,worst_vio_sqrt,report_TR,Vmax_VP,Vmin_VP,Imax_IP
end

function findTR_Q_coupled(V_0,Vmax_VP,Vmin_VP,VmaxIndex_P,VminIndex_P,Vmaxlimit0,Vminlimit0,S_VQ,QBound,raw,cl;Improve_Rate=5)

   
    # S_VQ=keep_top_bottom!.(S_VQ, 0.4)
    # S_VQ=keep_top_abs!.(S_VQ, 0.6)
    # S_VQ_org =deepcopy(S_VQ)
    # S_VQ=keep_top_abs!.(S_VQ, 0.25)
    # S_VQ=keep_K_percent_twoside!.(S_VQ,0.8)

    # S_VQ[1,1]
    # Plots.plot(S_VQ[1,1][:,32])
    z0=sum((V_0-Vmaxlimit0).>=0)+sum(Vminlimit0-V_0.>=0)
    zmax=round(size(Vmaxlimit0)[1]*size(Vmaxlimit0)[2]*size(Vmaxlimit0)[3]*Q_con_percent)+1
    # Vdown=(V_0-Vmaxlimit0)./ovr_Sens_VQ
    # Vup=(Vminlimit0-V_0)./ovr_Sens_VQ
    # TR=0.001
    # TR_max=0.05
    # TR_min=0
    # neighoborhood=50

    # Vmax_VP[VmaxIndex_P].=V_0[VmaxIndex_P]
    # Vmin_VP[VminIndex_P].=V_0[VminIndex_P]

    TR=0;TR_max=0.5;TR_min=0.0
    iter=0
    Itermax=10
    Vmaxinrisk=0
    Vmininrisk=0
    Vmaxviolated=0
    Vminviolated=0
    # worst_vio=max(maximum(Vmax_VP-Vmaxlimit0),maximum(-(Vmin_VP-Vminlimit0)))
    worst_vio=max(maximum(V_0-Vmaxlimit0),maximum(-(V_0-Vminlimit0)))

    # target_acc=(worst_vio/Improve_Rate)
    adj_acc=0
    acc=0
    best_acc=deepcopy(worst_vio)
    target_acc=0
    Vmaxlimit_new=deepcopy(Vmaxlimit0)
    # Vmaxlimit_new[:,:,2:end].+=target_acc
    # Vmaxlimit_new[:,Not(raw["Gen"].Bus),2:end].+=target_acc
    # Vminlimit_new=deepcopy(Vminlimit0)
    # Vminlimit_new[:,:,2:end].-=target_acc
    # Vminlimit_new[:,Not(raw["Gen"].Bus),2:end].-=target_acc
    # !(raw["Gen"].Bus)
    genbus=raw["Gen"].Bus
    genbus_s = Set(genbus)  # O(1) lookup
    allbus=1:raw["nBus"]
    allbus_s = Set(allbus)  # O(1) lookup


    nBus      = raw["nBus"]
    nQActBus  = raw["nQActBus"]
    ovr_Sens_VQ_up = zeros(nTP, raw["nBus"], cl);ovr_Sens_VQ_dn = zeros(nTP, raw["nBus"], cl)
    ovr_Sens_VQ_up_org = zeros(nTP, raw["nBus"], cl);ovr_Sens_VQ_dn_org = zeros(nTP, raw["nBus"], cl)

    ovr_Sens_VQ_upmat = zeros(Float64, nBus, nQActBus, nTP, cl)
    ovr_Sens_VQ_dnmat = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_up             = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_dn             = zeros(Float64, nBus, nQActBus, nTP, cl)
    upelem            = zeros(Float64, nBus, nQActBus, nTP, cl)
    dnelem            = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_up1      =falses(nBus, nQActBus, nTP, cl)
    VQ_dn1      =falses(nBus, nQActBus, nTP, cl)

    DV_up=Vmax_VP.- Vmaxlimit0
    DV_dn=Vmin_VP.- Vminlimit0

    DV_up_4D = zeros(Float64, nBus, 1, nTP, cl)
    DV_dn_4D = zeros(Float64, nBus, 1, nTP, cl)
    S_VQ_4D = zeros(Float64, nBus, nQActBus, nTP, cl)
    Threads.@threads for t in 1:nTP
        for c in 1:cl
            @views S_VQ_4D[:, :, t, c] .= S_VQ[t, c]
            @views DV_up_4D[:, 1, t, c] .= DV_up[t, :, c]  # direct assignment
            @views DV_dn_4D[:, 1, t, c] .= DV_dn[t, :, c]  # direct assignment

        end
    end
    S_VQ1=keep_K_percent_twoside_4D(S_VQ_4D,0.9)
    



    # ovr_Sens_VQ_upmat = [zeros(raw["nBus"], raw["nQActBus"]) for i in 1:nTP, j in 1:cl]  
    # ovr_Sens_VQ_dnmat = [zeros(raw["nBus"], raw["nQActBus"]) for i in 1:nTP, j in 1:cl]      
    # VQ_up = [zeros(raw["nBus"], raw["nQActBus"]) for i in 1:nTP, j in 1:cl]  
    # VQ_dn = [zeros(raw["nBus"], raw["nQActBus"]) for i in 1:nTP, j in 1:cl] 
    # upelem = [zeros(raw["nBus"], raw["nQActBus"]) for i in 1:nTP, j in 1:cl]  
    # dnelem = [zeros(raw["nBus"], raw["nQActBus"]) for i in 1:nTP, j in 1:cl]     

    # SVQ_tc=zeros(raw["nQActBus"])
    # up_max = similar(SVQ_tc);dn_max = similar(SVQ_tc);up=similar(SVQ_tc);dn=similar(SVQ_tc)
    Qinj_up,Qinj_dn=0,0
    y=0
    z=0
    TR=0.01
    # target_acc0=worst_vio*adj_acc
    # Vmaxlimit_new0=Vmaxlimit0.+(target_acc0)
    # Vminlimit_new0=Vminlimit0.-(target_acc0)
    Vmaxlimit_new=Vmaxlimit0
    Vminlimit_new=Vminlimit0
    report_TR= DataFrame(worst_vio=0.0,acc=0.0,TR=0.0,z=z0,y=z0,z_acc_Q=0.0,y_acc_Q=0.0,z_acc_P=0.0,y_acc_P=0.0,plotdata=[0.0])
    deleteat!(report_TR,1)
    z_acc_P=0;z_acc_Q=0
    y_acc_P=0;y_acc_Q=0

    TR_best=0.0
    # while adj_acc<=0.9
    #     # target_acc=max(worst_vio*adj_acc,0.0001)
    #     target_acc=worst_vio*adj_acc
    #     Vmaxlimit_new=Vmaxlimit0.+(target_acc)
    #     Vminlimit_new=Vminlimit0.-(target_acc)
    #     Vmaxlimit_new[:,raw["Gen"].Bus,:].-=target_acc
    #     Vminlimit_new[:,raw["Gen"].Bus,:].+=target_acc
    #     iter=0
    #     TR_min=0
    #     iter=0
    # TR=0
    #     if y!=0
    #         TR_max=0.5
    #     end
    # TR=0.03
    # Itermax=0
    while TR_max - TR_min > 1e-7  # Stop when the range is small enough
        VQ_up1.=0
        VQ_dn1.=0
        # Update TR using bisection
        if iter==0
            nothing
        elseif iter>=Itermax
            TR=TR_best
        # elseif (z_acc/y_acc)>5 && y==0
        #     TR=min(TR_max,max(TR/(z_acc/y_acc),1.2 *TR_min))
        else
            TR = (TR_min + TR_max) / 2
        end

        # sum(QBound[1].<QBound[2])
        # sum(PBound[1].<PBound[2])


        # findessential()

        # function findessential(V_0,Vmaxlimit0,Vminlimit0)
        #     Vmaxlimit0.+acc.-V_0
        # end
        Qinj_up=min.(TR,QBound[1])
        Qinj_dn=max.(-TR,QBound[2])
        Threads.@threads for t in 1:nTP
            for c in 1:cl
            @views S_VQ_tc = S_VQ[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
            @views up = Qinj_up[t, :, c]  # 1 × nPActBus
            @views dn = Qinj_dn[t, :, c]

            mask_V = S_VQ_tc .> 0  # nBus × nQActBus
            up_max_V = S_VQ_tc .* (mask_V .* up' .+ .!mask_V .* dn')
            dn_max_V = S_VQ_tc .* (mask_V .* dn' .+ .!mask_V .* up')

            @views ovr_Sens_VQ_up[t, :, c] .= sum(up_max_V, dims=2) .+ 1e-6
            @views ovr_Sens_VQ_dn[t, :, c] .= sum(dn_max_V, dims=2) .+ 1e-6
            @views ovr_Sens_VQ_up_org[t, :, c] .= sum(up_max_V, dims=2) .+ 1e-6
            @views ovr_Sens_VQ_dn_org[t, :, c] .= sum(dn_max_V, dims=2) .+ 1e-6

            @views ovr_Sens_VQ_upmat[:, :, t, c] .= up_max_V
            @views ovr_Sens_VQ_dnmat[:, :, t, c] .= dn_max_V
            end
        end


        @views VQ_up .= ovr_Sens_VQ_upmat .+ DV_up_4D
        @views VQ_dn .= -ovr_Sens_VQ_dnmat .- DV_dn_4D

        # size(ovr_Sens_VQ_upmat)


        # VQ_up[:,:,5,9]== (VQ_up6[5,9])
        # DV_up_4D[60,1,5,9]
        # DV_up[5,60,9]
        # sum(ovr_Sens_VQ_upmat[14,:,10,32])


        #  Threads.@threads for i in axes(ovr_Sens_VQ_dnmat,1)
        #     for j in axes(ovr_Sens_VQ_upmat,2)
        #     VQ_up[i,j]= ovr_Sens_VQ_upmat[i,j] .+ (Vmax_VP[i,:,j] .- Vmaxlimit0[i,:,j])
        #     VQ_dn[i,j]= -1*ovr_Sens_VQ_dnmat[i,j] .- (Vmin_VP[i,:,j] .- Vminlimit0[i,:,j])
        #     end
        # end
        # VQ_dn[13,31][17,:]
        # VQ_up= [1*ovr_Sens_VQ_upmat[i,j] .+ (Vmax_VP[i,:,j] .- Vmaxlimit0[i,:,j]) for i in axes(ovr_Sens_VQ_upmat,1), j in axes(ovr_Sens_VQ_upmat,2)]
        # VQ_dn= [-1*ovr_Sens_VQ_dnmat[i,j] .- (Vmin_VP[i,:,j] .- Vminlimit0[i,:,j]) for i in axes(ovr_Sens_VQ_dnmat,1), j in axes(ovr_Sens_VQ_dnmat,2)] 
        # VQ_up3=deepcopy(VQ_up)
        # VQ_dn3=deepcopy(VQ_dn)
        # VQ_up4=deepcopy(VQ_up)
        # VQ_dn4=deepcopy(VQ_dn)

        # VQ_up3= [(V_0[i,:,j] .- (Vmaxlimit0[i,:,j])) ./ -ovr_Sens_VQ_upmat[i,j]  for i in axes(ovr_Sens_VQ_upmat,1), j in axes(ovr_Sens_VQ_upmat,2)]
        # VQ_dn3= [(V_0[i,:,j] .- (Vminlimit0[i,:,j])) ./ (ovr_Sens_VQ_dnmat[i,j]) for i in axes(ovr_Sens_VQ_dnmat,1), j in axes(ovr_Sens_VQ_dnmat,2)]

        # VQ_up1= [ovr_Sens_VQ_upmat[i,j].+1e-6 for i in axes(ovr_Sens_VQ_upmat,1), j in axes(ovr_Sens_VQ_upmat,2)]
        # VQ_dn1= [-ovr_Sens_VQ_dnmat[i,j].+1e-6 for i in axes(ovr_Sens_VQ_dnmat,1), j in axes(ovr_Sens_VQ_dnmat,2)]
        # keep_K_percent_new!.(VQ_up,0.5)
        # keep_K_percent_new!.(VQ_dn,0.5)

        # S_VQ1=deepcopy(S_VQ)

        # VQ_up1=keep_K_percent_new_idx.(VQ_up,0.0)
        # VQ_dn1=keep_K_percent_new_idx.(VQ_dn,0.0)
        # # keep_top_bottom_row!.(VQ_up1,0.01)
        # # keep_top_bottom_row!.(VQ_dn1,0.01)

        # # VQ_up4=keep_top_idx.(VQ_up,10)
        # # VQ_dn4=keep_top_idx.(VQ_dn,10)

        # VQ_up4=keep_top_idx.(ovr_Sens_VQ_upmat,10)
        # VQ_dn4=keep_top_idx.(-ovr_Sens_VQ_dnmat,10)
        # VQ_up3=keep_top_row_gen_idx.(VQ_up,Ref(genbus_s),5)
        # VQ_dn3=keep_top_row_gen_idx.(VQ_dn,Ref(genbus_s),5)


        @views VQ_up1=keep_K_percent_new_4D!(VQ_up1,VQ_up,0.0)
        @views VQ_dn1=keep_K_percent_new_4D!(VQ_dn1,VQ_dn,0.0)

        @views VQ_up1=keep_top_4D!(VQ_up1,VQ_up,10)
        @views VQ_dn1=keep_top_4D!(VQ_dn1,VQ_dn,10)

        # VQ_up1[:,:,5,9]== (VQ_up5[5,9].!=0)


        # VQ_up1[59,:,5,9]
        # VQ_up5[5,9][59,:]
        # VQ_up[59,:,5,9]

        # VP_up4=keep_top_idx.(VP_up,3)
        # VP_dn4=keep_top_idx.(VP_dn,3)

        # VP_up4=keep_top_idx.(ovr_Sens_VP_upmat,1)
        # VP_dn4=keep_top_idx.(-ovr_Sens_VP_dnmat,1)

        # @views VQ_up1 =keep_top_4D!(VQ_up1,ovr_Sens_VQ_upmat,10)
        # @views VQ_dn1 =keep_top_4D!(VQ_dn1,-ovr_Sens_VQ_dnmat,10)
        # # @views VQ_up1 =keep_top_4D!(VQ_up1,S_VQ_4D,5)
        # # @views VQ_dn1 =keep_top_4D!(VQ_dn1,S_VQ_4D,5)
        # # @views VQ_up1 =keep_top_4D!(VQ_up1,-S_VQ_4D,5)
        # # @views VQ_dn1 =keep_top_4D!(VQ_dn1,-S_VQ_4D,5)

        # # VQ_up1[:,:,5,9]== (VQ_up4[5,9].!=0)




        # # ovr_Sens_VP_upmat[:,:,1,1]

        # # VP_up3=keep_top_row_gen_idx.(VP_up,Ref(genbus_s),1)
        # # VP_dn3=keep_top_row_gen_idx.(VP_dn,Ref(genbus_s),1)
        
        # @views VQ_up1 =keep_top_row_gen_4D!(VQ_up1,VQ_up,genbus_s,genprox)
        # @views VQ_dn1 =keep_top_row_gen_4D!(VQ_dn1,VQ_dn,genbus_s,genprox)
        # @views VQ_up1 =keep_top_row_gen_4D!(VQ_up1,VQ_up,allbus_s,10)
        # @views VQ_dn1 =keep_top_row_gen_4D!(VQ_dn1,VQ_dn,allbus_s,10)

        # VQ_up1[:,:,5,9]== (VQ_up3[5,9].!=0)

        # ceil(1.5)




        # keep_K_percent_new!.(VQ_up1,0.9)
        # keep_K_percent_new!.(VQ_dn1,0.9)
        # keep_K_percent_new!.(VQ_up3,0.8)
        # keep_K_percent_new!.(VQ_dn3,0.8)
        # # keep_K_percent_new!.(VQ_up1,0.05)
        # keep_K_percent_new!.(VQ_dn1,0.05)


        # Plots.plot(ovr_Sens_VQ_upmat[14,:,10,32])
        # ovr_Sens_VQ_up[10,14,32]
        # findmax(ovr_Sens_VQ_upmat[14,:,10,32])
        # raw["QActBus"][27:29]
        # Plots.plot(-S_VQ_4D[:,27,10,32])

        # S_VQ_4D[10:20,27,10,32]


        # (VQ_up[i,j][43,:].!=0)
        # findmax(upelem[18,2][43,:])

        # aa=(VQ_up[i,j].!=0) .|| (VQ_up3[i,j].!=0) .|| (S_VQ1[i,j].!=0)
        # # sum(aa,dims=1)

        # aa[46:60,:]

        # (((VQ_up[i,j].!=0) .|| (VQ_up1[i,j].!=0)))[43,:]
        # (((VQ_up[i,j].!=0) .|| (VQ_up1[i,j].!=0)) .* ovr_Sens_VQ_upmat[i,j])[43,:]

        @. upelem=(VQ_up1 ||S_VQ1) *ovr_Sens_VQ_upmat
        @. dnelem=(VQ_dn1 ||S_VQ1) *ovr_Sens_VQ_dnmat


        # sum(upelem,dims=2)[:,1,:,:]

        ovr_Sens_VQ_up.= permutedims(sum(upelem,dims=2)[:,1,:,:], (2, 1, 3)) 
        ovr_Sens_VQ_dn.= permutedims(sum(dnelem,dims=2)[:,1,:,:], (2, 1, 3)) 


        # aa=permutedims(sum(upelem,dims=2)[:,1,:,:], (2, 1, 3))  # (24, 60, 34)
        # aa[5,:,9]== ovr_Sens_VQ_up_cop[5,:,9]      
        # ovr_Sens_VQ_up==ovr_Sens_VQ_up_cop

        # sum(upelem,dims=2)[:,:,5,9]==ovr_Sens_VQ_up_cop[5,:,9]

        # Threads.@threads for i in axes(ovr_Sens_VQ_upmat,1) 
        #     for j in axes(ovr_Sens_VQ_upmat,2)
        #         upelem[i,j]= ((VQ_up4[i,j].!=0) .|| (VQ_up1[i,j].!=0) .|| (VQ_up3[i,j].!=0) .|| (S_VQ1[i,j].!=0)) .* ovr_Sens_VQ_upmat[i,j]
        #         # upelem[i,j][genbus,:]=ovr_Sens_VQ_upmat[i,j][genbus,:]
        #         dnelem[i,j]= ((VQ_dn4[i,j].!=0) .||(VQ_dn1[i,j].!=0) .|| (VQ_dn3[i,j].!=0) .|| (S_VQ1[i,j].!=0)) .* ovr_Sens_VQ_dnmat[i,j]
        #         # dnelem[i,j][genbus,:]=ovr_Sens_VQ_dnmat[i,j][genbus,:]
        #         upelemsum=sum(upelem[i,j],dims=2)
        #         dnelemsum=sum(dnelem[i,j],dims=2)
        #         ovr_Sens_VQ_up[i, :, j].= upelemsum
        #         ovr_Sens_VQ_dn[i, :, j].= dnelemsum
        #     end
        # end
        # upelem= [((VQ_up[i,j].!=0) .|| (VQ_up3[i,j].!=0) .|| (S_VQ1[i,j].!=0)) .* ovr_Sens_VQ_upmat[i,j] for i in axes(ovr_Sens_VQ_upmat,1), j in axes(ovr_Sens_VQ_upmat,2)]
        # upelemsum=sum.(upelem,dims=2)
        # dnelem= [((VQ_dn[i,j].!=0) .|| (VQ_dn3[i,j].!=0) .|| (S_VQ1[i,j].!=0)) .* ovr_Sens_VQ_dnmat[i,j] for i in axes(ovr_Sens_VQ_dnmat,1), j in axes(ovr_Sens_VQ_dnmat,2)]
        # dnelemsum=sum.(dnelem,dims=2)

        # i=1;j=7
        # ovr_Sens_VP_up[t, i, c]=sum(ovr_Sens_VP_upmat[t,c],dims=2)[i]

        
        #             b[1,1]

        # for t in axes(ovr_Sens_VQ_upmat,1), c in axes(ovr_Sens_VQ_upmat,2)
        #     ovr_Sens_VQ_up[t, :, c].= upelemsum[t,c]
        #     ovr_Sens_VQ_dn[t, :, c].= dnelemsum[t,c]
        # end
        # Vmaxviolated= V_0.>=Vmaxlimit_new
        # Vminviolated= V_0.<=Vminlimit_new

        # vec(V_0+ovr_Sens_VQ_up.-Vmaxlimit_new)
        # vec(-(V_0+ovr_Sens_VQ_dn.-Vminlimit_new))
        # v = vcat( vec(V_0+ovr_Sens_VQ_up.-Vmaxlimit_new), vec(-(V_0+ovr_Sens_VQ_dn.-Vminlimit_new)))
        # plot!(sort(vcat( vec(V_0+ovr_Sens_VQ_up.-Vmaxlimit_new),vec(-(V_0+ovr_Sens_VQ_dn.-Vminlimit_new))),rev=true))
        V_max_VPQ=(Vmax_VP+ovr_Sens_VQ_up/1 .-Vmaxlimit_new)
        V_max_VQ=(V_0+ovr_Sens_VQ_up/1 .-Vmaxlimit_new)
        V_max_VQ[VmaxIndex_P].=-10
        V_max_VPQ[VmaxIndex_P].=-10

        V_min_VPQ=-(Vmin_VP+ovr_Sens_VQ_dn/1 .-Vminlimit_new)
        V_min_VQ=-(V_0+ovr_Sens_VQ_dn/1 .-Vminlimit_new)
        V_min_VQ[VminIndex_P].=-10
        V_min_VPQ[VminIndex_P].=-10

        # V_max_VQ[15,14,21]

    # V_max_VPQ[11,41,16]
    # V_max_VQ[11,41,16]


        vals_VP = vcat( vec(V_max_VPQ),vec((V_min_VPQ)))
        k_VP = min(length(vals_VP), Int(zmax))
        z_acc_P = partialsort(vals_VP, 1:k_VP; rev=true)[end]

        # z_acc_P=sort(vcat( vec(V_max_VPQ),vec((V_min_VPQ))),rev=true)[min(end,Int(zmax))]

        # z_acc_Q=sort(vcat( vec(V_max_VQ),vec(V_min_VQ)),rev=true)[min(end,Int(zmax))]
        vals_VQ = vcat( vec(V_max_VQ),vec((V_min_VQ)))
        k_VQ = min(length(vals_VQ), Int(zmax))
        z_acc_Q = partialsort(vals_VQ, 1:k_VQ; rev=true)[end]

        # z_acc=sort(vcat( vec(Vmax_VP+ovr_Sens_VQ_up/1 .-Vmaxlimit_new),
        # vec(-(Vmin_VP+ovr_Sens_VQ_dn/1 .-Vminlimit_new))),rev=true)[min(end,Int(zmax))]
        # y_acc=max(maximum(V_0+(ovr_Sens_VQ_dn/Q_eff).-Vmaxlimit_new),maximum(-(V_0+(ovr_Sens_VQ_up/Q_eff).-Vminlimit_new)),0)
        y_acc_Q=max(maximum(V_0+(ovr_Sens_VQ_dn_org/Q_eff).-Vmaxlimit_new),maximum(-(V_0+(ovr_Sens_VQ_up_org/Q_eff).-Vminlimit_new)),0)

        y_acc_P=max(maximum(Vmax_VP+(ovr_Sens_VQ_dn_org/Q_eff).-Vmaxlimit_new),maximum(-(Vmin_VP+(ovr_Sens_VQ_up_org/Q_eff).-Vminlimit_new)),0)
        # y_acc_P=max(maximum(Vmax_VP+(ovr_Sens_VQ_dn/1).-Vmaxlimit_new),maximum(-(Vmin_VP+(ovr_Sens_VQ_up/1).-Vminlimit_new)),0)

        acc=max(-0.0,y_acc_Q)
        Vmaxinrisk_P = Vmax_VP+ovr_Sens_VQ_up/1 .>= Vmaxlimit_new.+max(z_acc_P,1e-4)
        Vmininrisk_P = Vmin_VP+ovr_Sens_VQ_dn/1 .<= Vminlimit_new.-max(z_acc_P,1e-4)

        Vmaxinrisk_Q = V_0+ovr_Sens_VQ_up/1 .>= Vmaxlimit_new.+max(z_acc_Q,1e-4)
        Vmininrisk_Q = V_0+ovr_Sens_VQ_dn/1 .<= Vminlimit_new.-max(z_acc_Q,1e-4)

        Vmaxinrisk= Vmaxinrisk_P #.|| Vmaxinrisk_Q
        Vmininrisk= Vmininrisk_P #.|| Vmininrisk_Q

        # 1 .|| 0
        z= sum(Vmaxinrisk) + sum(Vmininrisk)
        y = sum(Vmax_VP+(ovr_Sens_VQ_dn/Q_eff) .> Vmaxlimit_new.+y_acc_P) + sum(Vmin_VP+(ovr_Sens_VQ_up/Q_eff) .<Vminlimit_new.-y_acc_P)
        push!(report_TR,[worst_vio,acc,TR,z,y,z_acc_Q,y_acc_Q,z_acc_P,y_acc_P,0.0])
        # plotdata=[max.(0,sort(vcat( vec(V_0+ovr_Sens_VQ_up.-Vmaxlimit_new),vec(-(V_0+ovr_Sens_VQ_dn.-Vminlimit_new))),rev=true)),max.(0,sort(vcat( vec(V_0.+ovr_Sens_VQ_dn_org/Q_eff.-Vmaxlimit_new),vec(-(V_0.+ovr_Sens_VQ_up_org/Q_eff.-Vminlimit_new))),rev=true))]
        # push!(report_TR,[worst_vio,acc,TR,z,y,z_acc_Q,y_acc_Q,z_acc_P,y_acc_P,plotdata],promote=true)
        if iter>=Itermax
            break
        end

        if acc<=best_acc && z_acc_Q <= y_acc_Q && z_acc_P<=worst_vio  #&& max(z_acc_Q,z_acc_P) <= max(0.05,y_acc_Q)
            best_acc=acc
            TR_best=TR
        end
        if z_acc_P>=worst_vio  || ((y_acc_P<=max(0.001,maximum(Vmax_VP),maximum(-Vmin_VP))-0.001) && max(0,z_acc_Q)>=y_acc_Q) #|| z_acc_Q >= worst_vio #|| max(z_acc_P,z_acc_Q)>= max(0.05,y_acc_Q)
            TR_max = TR  # Reduce TR
        else
            # if z>zmax
            #    break
            # end
            TR_min = TR  # Increase TR
        end
        iter+=1
    end
    #     Rem_vio=max(maximum(V_0+ovr_Sens_VQ_dn .-Vmaxlimit_new),maximum(-(V_0+ovr_Sens_VQ_up .-Vminlimit_new)))
    #     target_acc=worst_vio*adj_acc+max(Rem_vio,0)
    #     push!(report_TR,[target_acc,TR,z,y,adj_acc])
    #     if z>zmax
    #         adj_acc+=0.1
    #     else
    #         break
    #     end
    # end
    # (V_0+ovr_Sens_VQ_up/1)[21,3,5]
    # TR=0.05
    VmaxIndex=findall(Vmaxinrisk)
    VminIndex=findall(Vmininrisk)

    # VmaxIndex=findall(.~isnan.(V_0))
    # VminIndex=findall(.~isnan.(V_0));z_acc_Q=0.0;
    display(report_TR)
    return TR,VmaxIndex,VminIndex,max(0.0,z_acc_Q),[Qinj_up,Qinj_dn],worst_vio,report_TR
end
function keep_top_row!(A::AbstractMatrix, k::Int=2)
    nrows, ncols = size(A)
    # k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    for j in 1:nrows
        row = view(A, j, :)

        # Find indices of k smallest and k largest elements
        # low_idx  = partialsortperm(row, 1:k; rev=false)
        high_idx = partialsortperm(row, 1:k; rev=true)

        # Build a boolean mask for the column
        mask = falses(ncols)
        # mask[low_idx] .= true
        mask[high_idx] .= true

        # Zero out in place
        @inbounds for i in 1:ncols
            if !mask[i]
                row[i] = 0
            end
        end
    end
    return A
end


function keep_top!(A::AbstractMatrix, k::Int64=5)
    nrows, ncols = size(A)
    # k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    for j in 1:ncols
        col = view(A, :, j)

        # Find indices of k smallest and k largest elements
        # low_idx  = partialsortperm(col, 1:k; rev=false)
        high_idx = partialsortperm(col, 1:k; rev=true)

        # Build a boolean mask for the column
        mask = falses(nrows)
        # mask[low_idx] .= true
        mask[high_idx] .= true

        # Zero out in place
        @inbounds for i in 1:nrows
            if !mask[i]
                col[i] = 0
            end
        end
    end
    return A
end


function keep_top_abs!(A::AbstractMatrix, frac::Float64=0.1)
    nrows, ncols = size(A)
    k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    for j in 1:ncols
        col = view(A, :, j)

        # Find indices of k smallest and k largest elements
        # low_idx  = partialsortperm(col, 1:k; rev=false)
        high_idx = partialsortperm(abs.(col), 1:k; rev=true)

        # Build a boolean mask for the column
        mask = falses(nrows)
        # mask[low_idx] .= true
        mask[high_idx] .= true

        # Zero out in place
        @inbounds for i in 1:nrows
            if !mask[i]
                col[i] = 0
            end
        end
    end
    return A
end

function keep_K_percent!(A::AbstractMatrix, percent::Float64=0.8)
    # A=S_VP[1,1]
    nrows, ncols = size(A)
    # k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    for j in 1:ncols
        col = view(A, :, j)

        # Find indices of k smallest and k largest elements
        # low_idx  = partialsortperm(col, 1:k; rev=false)
        # high_idx = partialsortperm(col, 1:k; rev=true)

        th=maximum(abs.(col))*percent
        idx=findall(x-> abs(x)>=th,col)


        # Build a boolean mask for the column
        mask = falses(nrows)
        mask[idx].=true

        # mask[low_idx] .= true
        # mask[high_idx] .= true

        # Zero out in place
        @inbounds for i in 1:nrows
            if !mask[i]
                col[i] = 0
            end
        end
    end
    return A
end

function keep_K_percent_new!(A::AbstractMatrix, percent::Float64=0.8)
    # A=S_VP[1,1]
    nrows, ncols = size(A)
    # k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    for j in 1:ncols
        col = view(A, :, j)

        # Find indices of k smallest and k largest elements
        # low_idx  = partialsortperm(col, 1:k; rev=false)
        # high_idx = partialsortperm(col, 1:k; rev=true)

        th=maximum(col)*percent
        idx=findall(x-> x>=th,col)


        # Build a boolean mask for the column
        mask = falses(nrows)
        mask[idx].=true

        # mask[low_idx] .= true
        # mask[high_idx] .= true

        # Zero out in place
        @inbounds for i in 1:nrows
            if !mask[i]
                col[i] = 0
            end
        end
    end
    return A
end


function keep_K_percent_twoside!(A::AbstractMatrix, percent::Float64=0.8)
    # A=S_VP[1,1]
    nrows, ncols = size(A)
    # k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    for j in 1:ncols
        col = view(A, :, j)

        # Find indices of k smallest and k largest elements
        # low_idx  = partialsortperm(col, 1:k; rev=false)
        # high_idx = partialsortperm(col, 1:k; rev=true)

        thmax=max(maximum(col)*percent,1e-4)
        thmin=min(minimum(col)*percent,-1e-4)
        high_idx=findall(x-> x>=thmax,col)
        low_idx=findall(x-> x<=thmin,col)


        # Build a boolean mask for the column
        mask = falses(nrows)
        # mask[idx].=true

        mask[low_idx] .= true
        mask[high_idx] .= true

        # Zero out in place
        @inbounds for i in 1:nrows
            if !mask[i]
                col[i] = 0
            end
        end
    end
    return A
end

function keep_K_percent_twoside_idx(A::AbstractMatrix, percent::Float64=0.8)
    # A=S_VP[1,1]
    B = zeros(Float64, size(A))
    nrows, ncols = size(B)
    # k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    Threads.@threads for j in 1:ncols
        col = view(A, :, j)
        # Find indices of k smallest and k largest elements
        # low_idx  = partialsortperm(col, 1:k; rev=false)
        # high_idx = partialsortperm(col, 1:k; rev=true)

        thmax=max(maximum(col)*percent,1e-4)
        thmin=min(minimum(col)*percent,-1e-4)
        high_idx=findall(x-> x>=thmax,col)
        low_idx=findall(x-> x<=thmin,col)
        B[low_idx,j] .= 1
        B[high_idx,j] .= 1
    end
    return B
end

function keep_top_row_gen!(A::AbstractMatrix,genbus_s, k::Int=5)
    # A=VP_up3[1,1] , k=5
    nrows, ncols = size(A)
    # k = min(k, ncols)
    genbus_set = collect(genbus_s)  # O(1) lookup

    Threads.@threads for j in 1:nrows
        row = view(A, j, :)
        if j in genbus_set
            keep = partialsortperm(row, 1:k; rev=true)  # thread-local allocation

            # zero everything else
            @inbounds for idx in 1:ncols
                if !(idx in keep)   # small O(k) lookup
                    row[idx] = 0
                end
            end
        else
            row .= 0
        end
    end

    return A
end
function keep_top_row_gen_idx(A::AbstractMatrix,genbus_s, k::Int=5)
    # A=VP_up3[1,1] , k=5
    B = zeros(Float64, size(A))
    nrows, ncols = size(B)
    # k = min(k, ncols)
    genbus_set = collect(genbus_s)  # O(1) lookup

    Threads.@threads for j in genbus_set
        row = view(A, j, :)
            keep = partialsortperm(row, 1:k; rev=true)  # thread-local allocation
            B[j,keep] .= 1
    end
    return B
end
function keep_top_row_gen_4D!(B::BitArray{4}, A::AbstractArray{T,4}, genbus_s, k::Int=5) where T
    nBus, nPActBus, nTP, cl = size(A)
    # B = falses(nBus, nPActBus, nTP, cl)  # BitArray initialized with false (0)

    genbus_set = collect(genbus_s)  # ensure array for iteration

    for t in 1:nTP, c in 1:cl
        @views slice = A[:, :, t, c]
        @views out   = B[:, :, t, c]
        # Thread over genbus_set to parallelize the row processing
        for j in genbus_set
            row = @view slice[j, :]
            keep = partialsortperm(row, 1:min(k, length(row)); rev=true)
            out[j, keep] .= 1
        end
    end

    return B
end
function keep_top_4D!(B::BitArray{4}, A::AbstractArray{T,4}, k::Int=5) where T
    # A=ovr_Sens_VP_upmat;k=5
    nBus, nPActBus, nTP, cl = size(A)
    # B = falses(nBus, nPActBus, nTP, cl)  # BitArray initialized with false (0)
    for t in 1:nTP, c in 1:cl
        @views slice = A[:, :, t, c]
        @views out   = B[:, :, t, c]
        # Thread over genbus_set to parallelize the row processing
        Threads.@threads for j in 1:nPActBus
            col = @view slice[:, j]
            keep = partialsortperm(col, 1:min(k, length(col)); rev=true)
            out[keep, j] .= 1
        end
    end

    return B
end



function keep_top_4D_V2!(
    B::BitArray{4},
    A::AbstractArray{T,4},
    k::Int = 5
) where {T<:Real}
    nBus, nPActBus, nTP, cl = size(A)
    k = min(k, nBus)

    # Parallelize across the outer dimensions (t, c)
    @threads for idx in 1:(nTP * cl)
        t = Int((idx - 1) % nTP) + 1
        c = Int((idx - 1) ÷ nTP) + 1

        @views slice = A[:, :, t, c]
        @views out   = B[:, :, t, c]

        # Reuse one temporary buffer per thread
        temp_idx = Vector{Int}(undef, k)
        temp_val = Vector{T}(undef, k)

        @inbounds for j in 1:nPActBus
            col = @view slice[:, j]

            # In-place top-k selection (no allocations)
            nk = partial_topk!(temp_idx, temp_val, col, k)
            @inbounds for p in 1:nk
                out[temp_idx[p], j] = true
            end
        end
    end

    return B
end

@inline function partial_topk!(
    idxbuf::Vector{Int},
    valbuf::Vector{T},
    v::AbstractVector{T},
    k::Int
) where {T<:Real}
    n = length(v)
    k = min(k, n)

    # Initialize buffers with first k elements
    @inbounds for i in 1:k
        idxbuf[i] = i
        valbuf[i] = v[i]
    end
    # Sort descending initially
    sortperm!(idxbuf, valbuf; rev=true)

    # Replace if larger element found
    @inbounds for i in (k+1):n
        x = v[i]
        if x > valbuf[end]
            valbuf[end] = x
            idxbuf[end] = i
            sortperm!(idxbuf, valbuf; rev=true)
        end
    end

    return k
end

function keep_top_idx(A::AbstractMatrix, k::Int64=5)
    nrows, ncols = size(A)
    B = zeros(Float64, size(A))
    # k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    Threads.@threads for j in 1:ncols
        col = view(A, :, j)

        # Find indices of k smallest and k largest elements
        # low_idx  = partialsortperm(col, 1:k; rev=false)
        high_idx = partialsortperm(col, 1:k; rev=true)

        # # Build a boolean mask for the column
        # mask = falses(nrows)
        # # mask[low_idx] .= true
        # mask[high_idx] .= true
        B[high_idx,j].=1
    end
    return B
end

function keep_K_percent_new_idx(A::AbstractMatrix, percent::Float64=0.8)
    # A=S_VP[1,1]
    nrows, ncols = size(A)
    B = zeros(Float64, size(A))
    # k = max(1, ceil(Int, frac * nrows))  # number of elements to keep

    Threads.@threads for j in 1:ncols
        col = view(A, :, j)
        th=maximum(col)*percent
        idx=findall(x-> x>=th,col)

        B[idx,j].=1.0

    end
    return B
end

function keep_K_percent_new_idx_V2!(B::Array{Float64,4}, A::Array{Float64,4}, percent::Float64=0.8)
    nrows, ncols, nTP, cl = size(A)
    fill!(B, 0.0)

    Threads.@threads for idx in 1:(nTP * cl * ncols)
        # unravel flat index → (j, t, c)
        j  = (idx - 1) % ncols + 1
        tc = div(idx - 1, ncols) + 1
        t  = (tc - 1) % nTP + 1
        c  = div(tc - 1, nTP) + 1

        col = @view A[:, j, t, c]
        th  = maximum(col) * percent
        @inbounds for i in 1:nrows
            if col[i] >= th
                B[i, j, t, c] = 1.0
            end
        end
    end

    return B
end

function keep_K_percent_twoside_4D(S::AbstractArray{T,4}, percent::Float64=0.8) where T
# S: (nBus, nPActBus, nTP, cl)
# percent: threshold fraction
# Returns a Float64 array of same shape: 1 where element is in top/bottom fraction, else 0

# Compute max/min per column (along rows, i.e., nBus)
thmax = maximum(S, dims=1) .* percent  # shape (1, nPActBus, nTP, cl)
thmin = minimum(S, dims=1) .* percent  # shape (1, nPActBus, nTP, cl)

# Create mask: 1 where >= thmax or <= thmin
B = @. (S >= thmax) | (S <= thmin)
return B
end

function keep_K_percent_new_4D!(B::BitArray{4}, S::AbstractArray{T,4}, percent::Float64=0.8) where T

thmax = maximum(S, dims=1) .* percent  # shape (1, nPActBus, nTP, cl)
@. B |= (S >= thmax) 
# return B
end

function keep_K_percent_new_4D_V2!(B::BitArray{4}, S::AbstractArray{T,4}, percent::Float64=0.8) where T
    @inbounds for j2 in axes(S,2), j3 in axes(S,3), j4 in axes(S,4)
        th = percent * maximum(@view S[:, j2, j3, j4])
        @inbounds @simd for j1 in axes(S,1)
            if S[j1, j2, j3, j4] >= th
                B[j1, j2, j3, j4] = true
            end
        end
    end
    return B
end

function findTR_Q_coupled2(V_0,Vmax_VP,Vmin_VP,VmaxIndex_P,VminIndex_P,Vmaxlimit0,Vminlimit0,S_VQ,QBound,raw,cl;Improve_Rate=5)
    
    z0=sum((V_0-Vmaxlimit0).>=0)+sum(Vminlimit0-V_0.>=0)
    zmax=round(size(Vmaxlimit0)[1]*size(Vmaxlimit0)[2]*size(Vmaxlimit0)[3]*Q_con_percent)+1
    TR=0;TR_max=0.5;TR_min=0.0
    iter=0
    Itermax=10
    Vmaxinrisk=falses(size(V_0));Vmininrisk=falses(size(V_0));
    Vmaxinrisk_P=falses(size(V_0));Vmininrisk_P=falses(size(V_0));
    # Vmaxinrisk_Q=falses(size(V_0));Vmaxinrisk_Q=falses(size(V_0));
    Vmaxviolated=0
    Vminviolated=0
    # worst_vio=max(maximum(Vmax_VP-Vmaxlimit0),maximum(-(Vmin_VP-Vminlimit0)))
    worst_vio=max(maximum(V_0-Vmaxlimit0),maximum(-(V_0-Vminlimit0)))

    # target_acc=(worst_vio/Improve_Rate)
    adj_acc=0
    acc=0
    best_acc=deepcopy(worst_vio)
    target_acc=0
    Vmaxlimit_new=deepcopy(Vmaxlimit0)
    genbus=raw["Gen"].Bus
    genbus_s = Set(genbus)  # O(1) lookup
    allbus=1:raw["nBus"]
    allbus_s = Set(allbus)  # O(1) lookup

    V0_length= length(V_0)
    nBus      = raw["nBus"]
    nQActBus  = raw["nQActBus"]
    ovr_Sens_VQ_up = zeros(nTP, raw["nBus"], cl);ovr_Sens_VQ_dn = zeros(nTP, raw["nBus"], cl)
    ovr_Sens_VQ_up_org = zeros(nTP, raw["nBus"], cl);ovr_Sens_VQ_dn_org = zeros(nTP, raw["nBus"], cl)

    ovr_Sens_VQ_upmat = zeros(Float64, nBus, nQActBus, nTP, cl)
    ovr_Sens_VQ_dnmat = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_up             = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_dn             = zeros(Float64, nBus, nQActBus, nTP, cl)
    # upelem            = zeros(Float64, nBus, nQActBus, nTP, cl)
    # dnelem            = zeros(Float64, nBus, nQActBus, nTP, cl)
    VQ_up1      =falses(nBus, nQActBus, nTP, cl)
    VQ_dn1      =falses(nBus, nQActBus, nTP, cl)
    V_max_VPQ   = similar(ovr_Sens_VQ_up);V_min_VPQ   = similar(ovr_Sens_VQ_dn)
    V_max_VQ   = similar(ovr_Sens_VQ_up);V_min_VQ   = similar(ovr_Sens_VQ_dn)
    vals_VP=zeros(Float64,2* V0_length);vals_VQ=zeros(Float64,2* V0_length)

    DV_up=Vmax_VP.- Vmaxlimit0
    DV_dn=Vmin_VP.- Vminlimit0

    DV_up_4D = zeros(Float64, nBus, 1, nTP, cl)
    DV_dn_4D = zeros(Float64, nBus, 1, nTP, cl)
    S_VQ_4D = zeros(Float64, nBus, nQActBus, nTP, cl)
    Threads.@threads for t in 1:nTP
        for c in 1:cl
            @views S_VQ_4D[:, :, t, c] .= S_VQ[t, c]
            @views DV_up_4D[:, 1, t, c] .= DV_up[t, :, c]  # direct assignment
            @views DV_dn_4D[:, 1, t, c] .= DV_dn[t, :, c]  # direct assignment
        end
    end
    S_VQ1=keep_K_percent_twoside_4D(S_VQ_4D,0.9)
    Qinj_up,Qinj_dn=0,0
    y=0
    z=0
    TR=0.01
    Vmaxlimit_new=Vmaxlimit0
    Vminlimit_new=Vminlimit0
    report_TR= DataFrame(worst_vio=0.0,acc=0.0,TR=0.0,z=z0,y=z0,z_acc_Q=0.0,y_acc_Q=0.0,z_acc_P=0.0,y_acc_P=0.0,plotdata=[0.0])
    deleteat!(report_TR,1)
    z_acc_P=0;z_acc_Q=0
    y_acc_P=0;y_acc_Q=0
    S_VQ_pos = max.(S_VQ_4D, 0.0)
    S_VQ_neg = min.(S_VQ_4D, 0.0)

    TR_best=0.0
    while TR_max - TR_min > 1e-7  # Stop when the range is small enough
        VQ_up1.=0
        VQ_dn1.=0
        # Update TR using bisection
        if iter==0
            nothing
        elseif iter>=Itermax
            TR=TR_best
        else
            TR = (TR_min + TR_max) / 2
        end
        Qinj_up=min.(TR,QBound[1])
        Qinj_dn=max.(-TR,QBound[2])
        Threads.@threads for t in 1:nTP
            for c in 1:cl
            # @views S_VQ_tc = S_VQ[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
            @views up = Qinj_up[t, :, c]  # 1 × nPActBus
            @views dn = Qinj_dn[t, :, c]

            @views @. ovr_Sens_VQ_upmat[:, :, t, c] = S_VQ_pos[:,:,t, c]* up' + S_VQ_neg[:,:,t, c]* dn'
            @views @. ovr_Sens_VQ_dnmat[:, :, t, c] = S_VQ_pos[:,:,t, c]* dn' + S_VQ_neg[:,:,t, c]* up'


            @views ovr_Sens_VQ_up[t, :, c] .= sum(ovr_Sens_VQ_upmat[:, :, t, c], dims=2)
            @views ovr_Sens_VQ_dn[t, :, c] .= sum(ovr_Sens_VQ_dnmat[:, :, t, c], dims=2)
            # @views ovr_Sens_VQ_upmat[:, :, t, c] .= up_max_V
            # @views ovr_Sens_VQ_dnmat[:, :, t, c] .= dn_max_V
            end
        end
        @views @. ovr_Sens_VQ_up_org = ovr_Sens_VQ_up + 1e-6
        @views @. ovr_Sens_VQ_dn_org = ovr_Sens_VQ_dn + 1e-6


        @views @. VQ_up = ovr_Sens_VQ_upmat + DV_up_4D
        @views @. VQ_dn = -ovr_Sens_VQ_dnmat - DV_dn_4D
        keep_K_percent_new_4D_V2!(VQ_up1,VQ_up,0.0)
        keep_K_percent_new_4D_V2!(VQ_dn1,VQ_dn,0.0)

        keep_top_4D_V2!(VQ_up1,VQ_up,10)
        keep_top_4D_V2!(VQ_dn1,VQ_dn,10)
        @. VQ_up1 |=S_VQ1
        @. VQ_dn1 |=S_VQ1
        @. ovr_Sens_VQ_upmat *= VQ_up1 
        @. ovr_Sens_VQ_dnmat *= VQ_dn1

        @views ovr_Sens_VQ_up.= permutedims(sum(ovr_Sens_VQ_upmat,dims=2)[:,1,:,:], (2, 1, 3)) 
        @views ovr_Sens_VQ_dn.= permutedims(sum(ovr_Sens_VQ_dnmat,dims=2)[:,1,:,:], (2, 1, 3)) 
        @. V_max_VPQ=(Vmax_VP+ovr_Sens_VQ_up/1 .-Vmaxlimit_new)
        @. V_max_VQ=(V_0+ovr_Sens_VQ_up/1 .-Vmaxlimit_new)
        @. V_max_VQ[VmaxIndex_P].=-10
        @. V_max_VPQ[VmaxIndex_P].=-10

        @. V_min_VPQ=-(Vmin_VP+ovr_Sens_VQ_dn/1 .-Vminlimit_new)
        @. V_min_VQ=-(V_0+ovr_Sens_VQ_dn/1 .-Vminlimit_new)
        @. V_min_VQ[VminIndex_P].=-10
        @. V_min_VPQ[VminIndex_P].=-10
        @views @. vals_VP[1:V0_length] = V_max_VPQ[:]
        @views @. vals_VP[1+V0_length:end] = V_min_VPQ[:]
        k_VP = min(length(vals_VP), Int(zmax))
        z_acc_P = partialsort(vals_VP, 1:k_VP; rev=true)[end]
        vals_VQ = vcat( vec(V_max_VQ),vec((V_min_VQ)))

        @views @. vals_VQ[1:V0_length] = V_max_VQ[:]
        @views @. vals_VQ[1+V0_length:end] = V_min_VQ[:]

        k_VQ = min(length(vals_VQ), Int(zmax))
        z_acc_Q = partialsort(vals_VQ, 1:k_VQ; rev=true)[end]
        y_acc_Q=max(maximum(V_0+(ovr_Sens_VQ_dn_org/Q_eff).-Vmaxlimit_new),maximum(-(V_0+(ovr_Sens_VQ_up_org/Q_eff).-Vminlimit_new)),0)

        y_acc_P=max(maximum(Vmax_VP+(ovr_Sens_VQ_dn_org/Q_eff).-Vmaxlimit_new),maximum(-(Vmin_VP+(ovr_Sens_VQ_up_org/Q_eff).-Vminlimit_new)),0)

        acc=max(-0.0,y_acc_Q)
        @. Vmaxinrisk_P = Vmax_VP+ovr_Sens_VQ_up/1 >= Vmaxlimit_new+max(z_acc_P,1e-4)
        @. Vmininrisk_P = Vmin_VP+ovr_Sens_VQ_dn/1 <= Vminlimit_new-max(z_acc_P,1e-4)

        # @. Vmaxinrisk_Q = V_0+ovr_Sens_VQ_up/1 .>= Vmaxlimit_new.+max(z_acc_Q,1e-4)
        # @. Vmininrisk_Q = V_0+ovr_Sens_VQ_dn/1 .<= Vminlimit_new.-max(z_acc_Q,1e-4)

        @. Vmaxinrisk= Vmaxinrisk_P #.|| Vmaxinrisk_Q
        @. Vmininrisk= Vmininrisk_P #.|| Vmininrisk_Q
        z= sum(Vmaxinrisk) + sum(Vmininrisk)
        y = sum(Vmax_VP+(ovr_Sens_VQ_dn/Q_eff) .> Vmaxlimit_new.+y_acc_P) + sum(Vmin_VP+(ovr_Sens_VQ_up/Q_eff) .<Vminlimit_new.-y_acc_P)
        push!(report_TR,[worst_vio,acc,TR,z,y,z_acc_Q,y_acc_Q,z_acc_P,y_acc_P,0.0])
        if iter>=Itermax
            break
        end

        if acc<=best_acc && z_acc_Q <= y_acc_Q && z_acc_P<=worst_vio  #&& max(z_acc_Q,z_acc_P) <= max(0.05,y_acc_Q)
            best_acc=acc
            TR_best=TR
        end
        if z_acc_P>=worst_vio  || ((y_acc_P<=max(0.001,maximum(Vmax_VP),maximum(-Vmin_VP))-0.001) && max(0,z_acc_Q)>=y_acc_Q) #|| z_acc_Q >= worst_vio #|| max(z_acc_P,z_acc_Q)>= max(0.05,y_acc_Q)
            TR_max = TR  # Reduce TR
        else
            TR_min = TR  # Increase TR
        end
        iter+=1
    end
    VmaxIndex=findall(Vmaxinrisk)
    VminIndex=findall(Vmininrisk)
    display(report_TR)
    return TR,VmaxIndex,VminIndex,max(0.0,z_acc_Q),[Qinj_up,Qinj_dn],worst_vio,report_TR
end

function findTR_P_coupled2(Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,V_0,Vmaxlimit0,Vminlimit0,S_VP,PBound,raw,cl,acc_V0)

    # acc_V0=0.05
    Imaxlimit0_sqrt=sqrt.(Imaxlimit0)
    z0=0
    zmax=round(size(Imaxlimit0)[1]*size(Imaxlimit0)[2]*size(Imaxlimit0)[3]*P_con_percent)+1
    zmax_V=round(size(Vmaxlimit0)[1]*size(Vmaxlimit0)[2]*size(Vmaxlimit0)[3]*Q_con_percent_coupled)+1
    y=0
    z=0
    TR=0;TR_min=0.0
    TR_max=5;
    V0_length= length(V_0)
    acc_V=min(acc_V0,max(maximum(V_0-Vmaxlimit0),maximum(-(V_0-Vminlimit0))))
    iter=0
    Iter_max=10
    # atrisk=0
    atrisk= falses(size(I2_0));
    # Vmaxinrisk=0;Vmininrisk=0;
    Vmaxinrisk=falses(size(V_0));Vmininrisk=falses(size(V_0));
    # Vmax_VP=0;Vmin_VP=0;Imax_IP=0;
    Vmax_VP=zeros(size(V_0));Vmin_VP=zeros(size(V_0)); Imax_IP=zeros(size(I2_0));
    violated=0
    worst_vio=maximum(I2_0-Imaxlimit0)
    worst_vio_sqrt=max.(maximum(sqrt.(I2_0)-Imaxlimit0_sqrt),maximum(Base_TR_N),maximum(Base_TR_C))
    best_acc=deepcopy(worst_vio_sqrt)
    acc=0
    z_acc_VP=0
    adj_acc=0
    target_acc_sqrt=0
    Imaxlimit_new=deepcopy(Imaxlimit0)
    genbus=raw["Gen"].Bus
    genbus_s = Set(genbus)  # O(1) lookup
    allbus_s= Set(1:raw["nBus"])
    ovr_Sens_IP_up = zeros(nTP, raw["nLine"], cl);ovr_Sens_IP_dn = zeros(nTP, raw["nLine"], cl)
    ovr_Sens_VP_up = zeros(nTP, raw["nBus"], cl);ovr_Sens_VP_dn = zeros(nTP, raw["nBus"], cl)
    nBus      = raw["nBus"]
    nPActBus  = raw["nPActBus"]

    # Allocate once, contiguous in memory
    ovr_Sens_VP_upmat = zeros(Float64, nBus, nPActBus, nTP, cl)
    ovr_Sens_VP_dnmat = zeros(Float64, nBus, nPActBus, nTP, cl)
    ovr_Sens_IP_upmat = zeros(Float64, raw["nLine"]+raw["nTrans"], nPActBus, nTP, cl)
    ovr_Sens_IP_dnmat = zeros(Float64, raw["nLine"]+raw["nTrans"], nPActBus, nTP, cl)
    VP_up             = similar(ovr_Sens_VP_upmat)
    VP_dn             = similar(ovr_Sens_VP_upmat)
    # upelem            = similar(ovr_Sens_VP_upmat)
    # dnelem            = similar(ovr_Sens_VP_upmat)
    VP_up1      =falses(nBus, nPActBus, nTP, cl)
    VP_dn1      =similar(VP_up1)
    vals_IP= zeros(Float64, length(I2_0))
    vals_VP= zeros(Float64, 2* length(V_0))

    DV_up = V_0.- (acc_V .+ Vmaxlimit0)
    DV_dn = V_0.- (-acc_V .+ Vminlimit0)
    DV_up_4D = zeros(Float64, nBus, 1, nTP, cl)
    DV_dn_4D = zeros(Float64, nBus, 1, nTP, cl)
    S_VP_4D = zeros(Float64, nBus, nPActBus, nTP, cl)
    S_IP_4D = zeros(Float64, raw["nLine"]+raw["nTrans"], nPActBus, nTP, cl)
    Threads.@threads for t in 1:nTP
        for c in 1:cl
            @views S_VP_4D[:, :, t, c] .= S_VP[t, c]
            @views S_IP_4D[:, :, t, c] .= S_IP[t, c]
            @views DV_up_4D[:, 1, t, c] .= DV_up[t, :, c]  # direct assignment
            @views DV_dn_4D[:, 1, t, c] .= DV_dn[t, :, c]  # direct assignment
        end
    end
    S_VP1=keep_K_percent_twoside_4D(S_VP_4D,0.99)
    Pinj_up,Pinj_dn=0,0
    Gen_ind=indexin(raw["Gen"]."Bus", raw["PActBus"])
    report_TR=DataFrame(worst_vio_sqrt=0.0,acc=0.0,TR=0.0,z=z0,y=z0,z_acc=0.0,y_acc=0.0,ramp_acc=0.0,z_V=0,acc_VP=0.0)
    deleteat!(report_TR,1)
    TR=1
    z_acc=0
    y_acc=0
    ramp_acc=0
    TR_best=0.0
    S_VP_pos = max.(S_VP_4D, 0.0)
    S_VP_neg = min.(S_VP_4D, 0.0)
    S_IP_pos = max.(S_IP_4D, 0.0)
    S_IP_neg = min.(S_IP_4D, 0.0)
        while TR_max - TR_min > 1e-7# Stop when the range is small enough
            if iter==0
                nothing
            elseif iter>=Iter_max
                TR=TR_best
            elseif (z_acc)/max(y_acc,ramp_acc)>5 && y==0
                TR=max(TR/((z_acc)/max(y_acc)),1.2 *TR_min)
            else
                TR = (TR_min + TR_max) / 2
            end

            VP_up1.=0
            VP_dn1.=0
            Pinj_up=min.(TR,PBound[1])
            Pinj_dn=max.(-TR,PBound[2])
            Threads.@threads for t in 1:nTP
                for c in 1:cl
                    # @views S_IP_tc = S_IP[t, c]  # Extract the 2D matrix (nBus × some_dim) only once
                    # @views S_VP_tc = S_VP[t,c]
                    @views up = Pinj_up[t, :, c]  # 1 × nPActBus
                    @views dn = Pinj_dn[t, :, c]
                    
                    @views @. ovr_Sens_IP_upmat[:, :, t, c] = S_IP_pos[:,:,t, c]* up' + S_IP_neg[:,:,t, c]* dn'
                    @views @. ovr_Sens_IP_dnmat[:, :, t, c] = S_IP_pos[:,:,t, c]* dn' + S_IP_neg[:,:,t, c]* up'
                    @views ovr_Sens_IP_up[t, :, c] .= sum(ovr_Sens_IP_upmat[:, :, t, c], dims=2)[1:raw["nLine"]] .+ 1e-6
                    @views ovr_Sens_IP_dn[t, :, c] .= sum(ovr_Sens_IP_dnmat[:, :, t, c], dims=2)[1:raw["nLine"]] .+ 1e-6

                    @views @. ovr_Sens_VP_upmat[:, :, t, c] .= S_VP_pos[:,:,t, c].* up' .+ S_VP_neg[:,:,t, c].* dn'
                    @views @. ovr_Sens_VP_dnmat[:, :, t, c] .= S_VP_pos[:,:,t, c].* dn' .+ S_VP_neg[:,:,t, c].* up'
                    @views ovr_Sens_VP_up[t, :, c] .= sum(ovr_Sens_VP_upmat[:, :, t, c], dims=2) .+ 1e-6
                    @views ovr_Sens_VP_dn[t, :, c] .= sum(ovr_Sens_VP_dnmat[:, :, t, c], dims=2) .+ 1e-6
                end
            end

            @views @. VP_up = ovr_Sens_VP_upmat + DV_up_4D
            @views @. VP_dn = -ovr_Sens_VP_dnmat - DV_dn_4D

            keep_K_percent_new_4D_V2!(VP_up1,VP_up,0.99)
            keep_K_percent_new_4D_V2!(VP_dn1,VP_dn,0.99)

            VP_up1 .|= S_VP1
            VP_dn1 .|= S_VP1
            @. ovr_Sens_VP_upmat *=VP_up1 #||S_VP1 
            @. ovr_Sens_VP_dnmat *=VP_dn1 #||S_VP1

            @views ovr_Sens_VP_up .= permutedims(sum(ovr_Sens_VP_upmat,dims=2)[:,1,:,:], (2, 1, 3)) 
            @views ovr_Sens_VP_dn .= permutedims(sum(ovr_Sens_VP_dnmat,dims=2)[:,1,:,:], (2, 1, 3)) 

            @views @. vals_IP = sqrt(max(0.0, I2_0[:] + ovr_Sens_IP_up[:])) - Imaxlimit0_sqrt[:]
            k_IP = min(length(vals_IP), Int(zmax))
            z_acc = partialsort(vals_IP, 1:k_IP; rev=true)[end]

            # @views vals_VP .= vcat(vec(V_0 .+ ovr_Sens_VP_up .- Vmaxlimit0),
            #    vec(-(V_0 .+ ovr_Sens_VP_dn .- Vminlimit0)))
            
             @views @. vals_VP[1:V0_length] = (V_0 + ovr_Sens_VP_up - Vmaxlimit0)[:]
             @views @. vals_VP[V0_length+1:end] = (-(V_0 .+ ovr_Sens_VP_dn .- Vminlimit0))[:]
        
            k_VP = min(length(vals_VP), Int(zmax_V))
            z_acc_VP = partialsort(vals_VP, 1:k_VP; rev=true)[end]

            @views Vmax_VP .= V_0 .+ ovr_Sens_VP_up
            @views Vmin_VP .= V_0 .+ ovr_Sens_VP_dn
            Vmaxinrisk .= Vmax_VP .>= Vmaxlimit0.+max(z_acc_VP,1e-4)
            Vmininrisk .= Vmin_VP .<= Vminlimit0.-max(z_acc_VP,1e-4)

        

            ramp_acc=0
            y_acc=max(maximum(sqrt.(max.(0,I2_0+ovr_Sens_IP_dn./P_eff)).-Imaxlimit0_sqrt),-0.1)
            acc=max(z_acc,y_acc)
            @. Imax_IP = I2_0 + ovr_Sens_IP_up/1
            @. atrisk = Imax_IP >= (Imaxlimit0_sqrt + max(z_acc,1e-4))^2

            z = sum(atrisk)
            z_V = sum(Vmaxinrisk)+sum(Vmininrisk)
            y = sum(@. (I2_0.+ovr_Sens_IP_dn./P_eff .>= (Imaxlimit0_sqrt .+ y_acc).^2))

            push!(report_TR,[worst_vio_sqrt,acc,TR,z,y,z_acc,y_acc,ramp_acc,z_V,z_acc_VP])
            if iter>=Iter_max
                break
            end

            if acc<=best_acc && z_acc_VP <= acc_V
                best_acc=max(0,acc)
                TR_best=TR
            end
            if z_acc>=max(y_acc) || y_acc<=0 || z_acc_VP >= acc_V
                TR_max = TR  # Reduce TR
            else
                TR_min = TR  # Increase TR
            end
            iter+=1
        end
    display(report_TR)
    Vmaxindex_P = findall(Vmaxinrisk)
    Vminindex_P = findall(Vmininrisk)
    println("Number of voltage in risk due to P adjustment: ", length(Vmaxindex_P)+length(Vminindex_P))
    println("voltage acc left out due to P adjustment: ", z_acc_VP)
    ImaxIndex = findall(atrisk)
    # ImaxIndex = findall(.~isnan.(I2_0));    z_acc=0;
    return TR,ImaxIndex,Vmaxindex_P,Vminindex_P,max(0.0,z_acc),[Pinj_up,Pinj_dn],worst_vio,worst_vio_sqrt,report_TR,Vmax_VP,Vmin_VP,Imax_IP
end
