function DA_OPF(LinMethod="None",OP_init=Dict{String,Any}(),init__="Zero") #Array{Union{Float64,Array{Float64}}}
    # LinMethod="None";OP_init=Dict{String,Any}()
    Ptype_="DA";init_=init__;
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
    if loadbusinjection=="true"
        raw["PActBus"]= filter(x-> !isempty(x."load")||!isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
        raw["QActBus"]= filter(x-> !isempty(x."load")||!isempty(x."shunt") || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
    else
        raw["PActBus"]= filter(x-> !isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
        raw["QActBus"]= filter(x->  !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus] #!isempty(x."shunt") ||
    end
    raw["nPActBus"]=size(raw["PActBus"])[1];raw["nQActBus"]=size(raw["QActBus"])[1]
    raw["Bus"][raw["PActBus"],"PActive"].=1;raw["Bus"][raw["QActBus"],"QActive"].=1

    glob_input=(nTP,TD,raw,initraw,initraw,Mode,Inc_coef,Dec_coef,CCT,LSpermit,ctdall)

    Proffin=Load_Data(nTP,raw,profraw,ctdall) # Prepare Load and Shunt data and generation setpoints
    SCOPF,N_Var,C_Var=try
       L_MPSCOPF(glob_input,Proffin;fixed_SOC_=1,conprob=conprob_,OPFon=1,LinMethod=LinMethod,Lininit=OP_init);
        # fixed_SOC_=1;conprob=conprob_;OPFon=1;LinMethod=LinMethod;Lininit=OP_init
    catch e
        @warn "Error in MPSCOPF function: $e"
        global debug_hint=glob_input,Proffin
        rethrow(e)
    end

    OPFout,Initials=try  Get_Result(SCOPF,raw,initraw;OPFon=1)
    catch e
        @warn "Error in Get_Result function: $e"
        global debug_hint= SCOPF,N_Var,C_Var
        rethrow(e)
    end
    if printlevel[2] ==1
        try
            ODSW(Path*"/../Result/DAOPF",OPFout)
            set_initial(Path*"DA_OPF",Initials)
        catch e
            @warn "Error in ODSW or set_initial function: $e"
            global debug_hint = SCOPF,OPFout,Initials
            rethrow(e)
        end
        ods_write(Path*"/../Result/Overview",Dict(("Obj",1,1)=>[["DA_OPF"] value(SCOPF[:real_objective])]))
    end
    return OPFout,SCOPF
end

function DA_SCOPF(Glob, casreps=nothing,h_cas=DataFrame();LinMethod="None",OP_init=Dict{String,Any}(),init__="Zero",relaxconfac=0.0) #Array{Union{Float64,Array{Float64}}}

    #init__="Zero";LinMethod="None";OP_init=Dict{String,Array{Union{Float64,Array{Float64}}}}()
    Ptype_="DA";init_=init__;
    nTP=Glob["nTL"];
    TD=Glob["Td"];
    raw,initraw,profraw= try Data_Reader(Glob;Ptype=Ptype_ , init=init_)
    catch e
        @warn "Problem with running Data_Reader function"
        rethrow(e)
    end
    # if !isempty(OP_init)
    #     for key in collect(intersect(keys(initraw),keys(OP_init)))
    #         initraw[key]=convert.(Float64,initraw[key])
    #         initraw[key].=OP_init[key]
    #     end
    # end
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
    if Glob["loadbusinjection"]=="true"
        raw["PActBus"]= filter(x-> !isempty(x."load")||!isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
        raw["QActBus"]= filter(x-> !isempty(x."load")||!isempty(x."shunt") || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
    else
        raw["PActBus"]= filter(x-> !isempty(intersect(x."load",Activeloadbus))||!isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
        raw["QActBus"]= filter(x-> !isempty(intersect(x."load",raw["PActBus"])) || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus] #||!isempty(x."shunt")
        # raw["PActBus"]= filter(x-> !isempty(x."ess")||!isempty(x."res")||!isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
        # raw["QActBus"]= filter(x-> !isempty(x."shunt") || !isempty(x."gen")||!isempty(x."fl"),ctdall)[!,:nbus]
    end
    raw["nPActBus"]=size(raw["PActBus"])[1];raw["nQActBus"]=size(raw["QActBus"])[1]
    raw["Bus"][raw["PActBus"],"PActive"].=1;raw["Bus"][raw["QActBus"],"QActive"].=1

    glob_input=(nTP,TD,raw,initraw,initraw,Glob["Mode"],Glob["Inc_coef"],Glob["Dec_coef"],Glob["CCT"],Glob["LSpermit"],ctdall)

    Proffin=Load_Data(nTP,raw,profraw,ctdall) # Prepare Load and Shunt data and generation setpoints
    SCOPF,N_Var,C_Var=try
        L_MPSCOPF(Glob,glob_input,Proffin,casreps,h_cas;conprob=Glob["conprob_"],LinMethod=LinMethod,Lininit=OP_init,relaxconfac=relaxconfac);
        # conprob=Glob["conprob_"];LinMethod=LinMethod;Lininit=OP_init
    catch e
        @warn "Error in MPSCOPF function: $e"
        debug_hint=glob_input,Proffin
        rethrow(e)
    end

    OPFout,Initials=try  Get_Result(SCOPF,raw,initraw)
    catch e
        @warn "Error in Get_Result function: $e"
        debug_hint= SCOPF,N_Var,C_Var
        rethrow(e)
    end

if printlevel[2] ==1
    try
        ODSW(Glob["Path"]*"/../Result/DASCOPF",OPFout)
        set_initial(Glob["Path"]*"DA_SCOPF",Initials)
    catch e
        @warn "Error in ODSW or set_initial function: $e"
        debug_hint = SCOPF,OPFout,Initials
        rethrow(e)
    end
    ods_write(Glob["Path"]*"/../Result/Overview",Dict(("Obj",2,1)=>[["DA_SCOPF"] value(SCOPF[:real_objective])]))
end
    return OPFout,SCOPF
end

function Post_Contingency_SP_OPF(raw_temp,CSP,Prof,t)
    raw_temp["Gen"]."Ramp".=SysData["Gen"]."Ramp"; #Inverse #Mousa selections: ask him for details
    CSP_temp=Dict(key => CSP[key][t:t,:] for key in keys(CSP))
    local Prof_temp=Dict(key => Prof[key][t:t,:] for key in keys(Prof))
    local glob_input=(1,TD,raw_temp,CSP_temp,CSP_temp,Mode,Inc_coef,Dec_coef,CCT,LSpermit)
    result,-,SCOPF=MPSCOPF(glob_input,Prof_temp;fixed_SOC_=0,conprob=0,OPFon=1,Winit=1)

    SCOPF,N_Var,C_Var=try
        MPSCOPF(glob_input,Prof_temp;fixed_SOC_=0,conprob=0,OPFon=1,Winit=1)
    catch e
        @warn "Error in MPSCOPF function: $e"
        global debug_hint=glob_input,Prof_temp
        rethrow(e)
    end

    result,- =try  Get_Result(SCOPF,raw_temp,CSP_temp;OPFon=1)
    catch e
        @warn "Error in Get_Result function: $e"
        global debug_hint= SCOPF,N_Var,C_Var
        rethrow(e)
    end
    return result,primal_status(SCOPF)
end

function L_RT_OPF!(initial_point,Initials,t0,LinMethod,Lininit_0)
    Ptype_="RT";init_="DA_SCOPF";
    raw,initraw,profraw= Data_Reader(;Ptype=Ptype_ , init=init_,t_0=t0)
    winit=0

    if(initial_point==1)
        winit=1
        set_keys=collect(keys(initraw))
        initial=Dict()
        for x in set_keys
            initial[x]=Initials[:,x][1]
        end
    else
        initial=initraw
        initial["UE"].=0
        for i=1:raw["nESS"] initial["soc"][:,i].=raw["ESS"]."E_initial (MWh)"[i]/raw["sbase"] end
    end

    Ncost0,Vio_Sum_N0,OP_init0_NV0,Vio_general_N0=PF_Valid_N(Lininit_0)
    if Problem_Type=="SCOPF"
        AVG_C_cost0,Vio_Sum_C0,Lininit,Vio_general_C0=PF_Valid_C(OP_init0_NV0)
    elseif Problem_Type=="OPF"
        Lininit=OP_init0_NV0
        Vio_Sum_C0=DataFrame(Nins=0,NVio=0,Mean=0,Max=0,Type="")
        AVG_C_cost0=Ncost0
        min_infmax_C=0
    end

    ctdall=ctdMaker(raw)
    glob_input=(nTP,TD,raw,initraw,initraw,Mode,Inc_coef,Dec_coef,CCT,LSpermit,ctdall)
    Proffin=Load_Data(nTP,raw,profraw,ctdall) # Prepare Load and Shunt data and generation setpoints
    OPF,N_Var,C_Var=try
        L_MPSCOPF(glob_input,Proffin;fixed_SOC_=2,conprob=0,OPFon=1,Winit=winit,t0=t0,WP=0,LinMethod=LinMethod,Lininit=Lininit);
        #fixed_SOC_=2;conprob=0;OPFon=1;Winit=winit;t0=t0;WP=0
        # MPSCOPF(glob_input,Proffin;fixed_SOC_=2,conprob=0,OPFon=1,Winit=winit,t0=t0,WP=0);
    catch e
        @warn "Error in MPSCOPF function: $e"
        global debug_hint=glob_input,Proffin
        rethrow(e)
    end

    OPFout,Initials=try  Get_Result(OPF,raw,initraw;OPFon=1)
    catch e
        @warn "Error in Get_Result function: $e"
        global debug_hint= OPF,N_Var,C_Var
        rethrow(e)
    end

    if printlevel!=0
        try
            ODSW(Path*"/../Result/RTOPF",OPFout)
            set_initial(Path*"RTOPF",Initials)
        catch e
            @warn "Error in ODSW or set_initial function: $e"
            global debug_hint = SCOPF,OPFout,Initials
            rethrow(e)
        end
        ods_write(Path*"/../Result/Overview",Dict(("Obj",3,1)=>[["RT_OPF"] value(OPF[:real_objective])]))
    end
    return OPFout,Initials, solve_time(OPF),OPF
end


