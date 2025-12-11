function Initialize(type_,out,init_model,DA_out,DA_model)
    #type_,out=Linearization_INIT[Runcode[1]],init_out
    V_1,θ_1,PG_1,QG_1,cost1=sqrt.(value.(init_model[:f]).^2 + value.(init_model[:e]).^2) ,atan.(value.(init_model[:f]),value.(init_model[:e])),value.(init_model[:PG]),value.(init_model[:QG]),out."OFValue_WP"[1]
    shed_1,curt_1,curtc_1=value.(init_model[:shed]),value.(init_model[:curt]),value.(init_model[:curtc])
    Pinc_1,Pdec_1,Pch_1,Pdis_1=value.(init_model[:P_inc]),value.(init_model[:P_dec]),value.(init_model[:P_ch]),value.(init_model[:P_dis])
    UE_1,soc_1=value.(init_model[:UE]),value.(init_model[:soc])
    if Problem_Type in ["SCOPF", "RBOPF"]
        V_c_Line_1,θ_c_Line_1,PG_c_Line_1,QG_c_Line_1=sqrt.(value.(init_model[:e_c_Line].^2 +init_model[:f_c_Line].^2)),atan.(value.(init_model[:f_c_Line]),value.(init_model[:e_c_Line])),value.(init_model[:PG_c_Line]),value.(init_model[:QG_c_Line])
        shed_c_Line_1,curt_c_Line_1,curtc_c_Line_1=value.(init_model[:shed_c_Line]),value.(init_model[:curt_c_Line]),value.(init_model[:curtc_c_Line])
        Pinc_c_Line_1,Pdec_c_Line_1,Pch_c_Line_1,Pdis_c_Line_1=value.(init_model[:P_inc_c_Line]),value.(init_model[:P_dec_c_Line]),value.(init_model[:P_ch_c_Line]),value.(init_model[:P_dis_c_Line])
    end
    if type_=="Flat"
        V_0=ones(nTL,size(SysData["Bus"],1));V_c_Line_0=ones(nTL,size(SysData["Bus"],1),size(SysData["ConL"],1));
        θ_0=zeros(nTL,size(SysData["Bus"],1));θ_c_Line_0=zeros(nTL,size(SysData["Bus"],1),size(SysData["ConL"],1));
        PG_0=zeros(nTL,size(SysData["Gen"],1));PG_c_Line_0=zeros(nTL,size(SysData["Gen"],1),size(SysData["ConL"],1));
        QG_0=zeros(nTL,size(SysData["Gen"],1));QG_c_Line_0=zeros(nTL,size(SysData["Gen"],1),size(SysData["ConL"],1));
        shed_c_Line_0=zeros(nTL,size(SysData["Load"],1),size(SysData["ConL"],1));curt_c_Line_0=zeros(nTL,size(SysData["RES"],1),size(SysData["ConL"],1));curtc_c_Line_0=zeros(nTL,size(SysData["Gen"],1),size(SysData["ConL"],1));
        Pinc_c_Line_0=zeros(nTL,size(SysData["FL"],1),size(SysData["ConL"],1));Pdec_c_Line_0=zeros(nTL,size(SysData["FL"],1),size(SysData["ConL"],1));Pch_c_Line_0=zeros(nTL,size(SysData["ESS"],1),size(SysData["ConL"],1));Pdis_c_Line_0=zeros(nTL,size(SysData["ESS"],1),size(SysData["ConL"],1));
        shed_0=zeros(nTL,size(SysData["Load"],1));curt_0=zeros(nTL,size(SysData["RES"],1));curtc_0=zeros(nTL,size(SysData["Gen"],1));
        Pinc_0=zeros(nTL,size(SysData["FL"],1));Pdec_0=zeros(nTL,size(SysData["FL"],1));Pch_0=zeros(nTL,size(SysData["ESS"],1));Pdis_0=zeros(nTL,size(SysData["ESS"],1));
        UE_0=zeros(nTL,size(SysData["FL"],1));soc_0=zeros(nTL,size(SysData["ESS"],1));
        for t in 1:nTL, cl in SysData["ConL"]."cont"
            PG_0[t,:].=(SysData["Gen"]."P_max (pu)"+SysData["Gen"]."P_min (pu)")./2
            QG_0[t,:].=(SysData["Gen"]."Q_max (pu)"+SysData["Gen"]."Q_min (pu)")./2
            if Problem_Type in ["SCOPF", "RBOPF"]
                PG_c_Line_0[t,:,cl].=(SysData["Gen"]."P_max (pu)"+SysData["Gen"]."P_min (pu)")./2
                QG_c_Line_0[t,:,cl].=(SysData["Gen"]."Q_max (pu)"+SysData["Gen"]."Q_min (pu)")./2
            end
        end
        obj=0
        obj_WP=0
        sol_time=0
    elseif type_=="Exact"
        V_0,θ_0,PG_0,QG_0,cost=V_1,θ_1,PG_1,QG_1,cost1
        shed_0,curt_0,curtc_0=shed_1,curt_1,curtc_1
        Pinc_0,Pdec_0,Pch_0,Pdis_0=Pinc_1,Pdec_1,Pch_1,Pdis_1
        UE_0,soc_0=UE_1,soc_1
        if Problem_Type in ["SCOPF", "RBOPF"]
            V_c_Line_0,θ_c_Line_0,PG_c_Line_0,QG_c_Line_0=V_c_Line_1,θ_c_Line_1,PG_c_Line_1,QG_c_Line_1 
            shed_c_Line_0,curt_c_Line_0,curtc_c_Line_0=shed_c_Line_1,curt_c_Line_1,curtc_c_Line_1
            Pinc_c_Line_0,Pdec_c_Line_0,Pch_c_Line_0,Pdis_c_Line_0=Pinc_c_Line_1,Pdec_c_Line_1,Pch_c_Line_1,Pdis_c_Line_1
        end
        obj=value(init_model[:real_objective])
        obj_WP=out."OFValue_WP"[1]
        sol_time=solve_time(init_model)
    elseif type_[1:2]=="DA"

        if Problem_Type in ["SCOPF", "RBOPF"]
            V_c_Line_0,θ_c_Line_0,PG_c_Line_0,QG_c_Line_0=sqrt.(value.(DA_model[type_][:e_c_Line].^2 +DA_model[type_][:f_c_Line].^2)),atan.(value.(DA_model[type_][:f_c_Line]),value.(DA_model[type_][:e_c_Line])),value.(DA_model[type_][:PG_c_Line]),value.(DA_model[type_][:QG_c_Line])
            shed_c_Line_0,curt_c_Line_0,curtc_c_Line_0=value.(DA_model[type_][:shed_c_Line]),value.(DA_model[type_][:curt_c_Line]),value.(DA_model[type_][:curtc_c_Line])
            Pinc_c_Line_0,Pdec_c_Line_0,Pch_c_Line_0,Pdis_c_Line_0=value.(DA_model[type_][:P_inc_c_Line]),value.(DA_model[type_][:P_dec_c_Line]),value.(DA_model[type_][:P_ch_c_Line]),value.(DA_model[type_][:P_dis_c_Line])
        end
        V_0,θ_0,PG_0,QG_0=sqrt.(value.(DA_model[type_][:e].^2 +DA_model[type_][:f].^2)),atan.(value.(DA_model[type_][:f]),value.(DA_model[type_][:e])),value.(DA_model[type_][:PG]),value.(DA_model[type_][:QG]),DA_out[type_]."OFValue_WP"[1]
        shed_0,curt_0,curtc_0=value.(DA_model[type_][:shed]),value.(DA_model[type_][:curt]),value.(DA_model[type_][:curtc])
        Pinc_0,Pdec_0,Pch_0,Pdis_0=value.(DA_model[type_][:P_inc]),value.(DA_model[type_][:P_dec]),value.(DA_model[type_][:P_ch]),value.(DA_model[type_][:P_dis])
        UE_0,soc_0=value.(DA_model[type_][:UE]),value.(DA_model[type_][:soc])
        obj=value(DA_model[type_][:real_objective])
        obj_WP=DA_out[type_]."OFValue_WP"[1]
        sol_time=solve_time(DA_model[type_])
        # ods_write(Path*"../Result/LIN_Sum_"*Problem_Type,Dict(("General",20,1)=>DataFrame(["-" "Non-Linear" 0 1 0 0 cost PG_0 QG_0 V_0 θ_0],   ["Initialization";"Method";"Runcode";"Iter";"OptGap";"MinInfIns";"Obj";["PG$i" for i in 1:size(SysData["Gen"],1)];["QG$i" for i in 1:size(SysData["Gen"],1)];["V$i" for i in 1:size(SysData["Bus"],1)];["θ$i" for i in 1:size(SysData["Bus"],1)]])))
    end 
    Lin_init=Dict{String,Any}() #Union{Float64,Array{Union{Float64,Array{Float64}}}}
    Lin_init["V"]=V_0;Lin_init["θ"]=θ_0;Lin_init["PG"]=PG_0;Lin_init["QG"]=QG_0;
    Lin_init["Pshed"]=shed_0;Lin_init["Pcurt"]=curt_0;Lin_init["Pcurtc"]=curtc_0;
    Lin_init["Pinc"]=Pinc_0;Lin_init["Pdec"]=Pdec_0;Lin_init["Pch"]=Pch_0;Lin_init["Pdis"]=Pdis_0;
    Lin_init["UE"]=UE_0;Lin_init["soc"]=soc_0;
    if Problem_Type in ["SCOPF", "RBOPF"]
        Lin_init["V_c_Line"]=V_c_Line_0;Lin_init["θ_c_Line"]=θ_c_Line_0;Lin_init["PG_c_Line"]=PG_c_Line_0;Lin_init["QG_c_Line"]=QG_c_Line_0;
        Lin_init["Pshed_c_Line"]=shed_c_Line_0;Lin_init["Pcurt_c_Line"]=curt_c_Line_0;Lin_init["Pcurtc_c_Line"]=curtc_c_Line_0;
        Lin_init["Pinc_c_Line"]=Pinc_c_Line_0;Lin_init["Pdec_c_Line"]=Pdec_c_Line_0;Lin_init["Pch_c_Line"]=Pch_c_Line_0;Lin_init["Pdis_c_Line"]=Pdis_c_Line_0;
    end
    Lin_init["TR"]=Adaptive_TR[1]
    Lin_init["TRB"]=[Adaptive_TRB[1]]
    Lin_init["FilterIndex"]=[filter_factor[1]]
    Lin_init["warmstart_yang"]=[0.0]
    Lin_init["obj"]=obj
    Lin_init["obj_WP"]=obj_WP
    Lin_init["time"]=sol_time
    return Lin_init
end

function Optimizer(LinMethod,OP_init,casreps=nothing)
    # LinMethod=Linearization_TYPE[Runcode[2]];casreps=nothing
    OP=Dict{String,Any}()
    if Problem_Type=="OPF"
        out,model=DA_OPF(LinMethod=LinMethod,OP_init=OP_init)
    elseif Problem_Type in ["SCOPF", "RBOPF"]
        if Problem_Type == "RBOPF" relaxconfac=0.2 else relaxconfac=0.0 end
        out,model= DA_SCOPF(casreps;LinMethod=LinMethod,OP_init=OP_init,relaxconfac=relaxconfac)
        OP["V_c_Line"]=sqrt.(value.(model[:CV2magexp])) ;OP["θ_c_Line"]=value.(model[:Cθexp]);OP["PG_c_Line"]=value.(model[:PG_c_Line]);OP["QG_c_Line"]=value.(model[:QG_c_Line])
        OP["Pshed_c_Line"]=value.(model[:shed_c_Line]);OP["Pcurt_c_Line"]=value.(model[:curt_c_Line]);OP["Pcurtc_c_Line"]=value.(model[:curtc_c_Line])
        OP["Pinc_c_Line"]=value.(model[:P_inc_c_Line]);OP["Pdec_c_Line"]=value.(model[:P_dec_c_Line]);OP["Pch_c_Line"]=value.(model[:P_ch_c_Line]);OP["Pdis_c_Line"]=value.(model[:P_dis_c_Line])
        OP["soc_c_Line"]=value.(model[:soc_c_Line]);
        # OP["UE_c_Line"]=value.(model[:UE_c_Line])
        if LinMethod=="SENS" OP["V_c_Line"]=value.(model[:CV2magexp]); end
    end
    OP["V"]=Matrix(out[1,"mag_v"]);OP["θ"]=Matrix(out[1,"ang_v"]);OP["PG"]=Matrix(out[1,"PG"]);OP["QG"]=Matrix(out[1,"QG"])
    OP["Pshed"]=Matrix(out[1,"Pshed"]);OP["Pcurt"]=Matrix(out[1,"Pcurt"]);OP["Pcurtc"]=Matrix(out[1,"Pcurtc"])
    OP["Pinc"]=Matrix(out[1,"Pinc"]);OP["Pdec"]=Matrix(out[1,"Pdec"]);OP["Pch"]=Matrix(out[1,"Pch"]);OP["Pdis"]=Matrix(out[1,"Pdis"])
    OP["soc"]=Matrix(out[1,"soc"]);
    OP["TR"]=OP_init["TR"];OP["TRB"]=OP_init["TRB"];OP["FilterIndex"]=OP_init["FilterIndex"]
    OP["acc"]=model[:acc]
    OP["obj"]= value(model[:real_objective])
    OP["obj_WP"]=out."OFValue_WP"[1]
    OP["time"]=solve_time(model)
    #OP["UE"]=Matrix(out[1,"UE"])
    return OP,out,model
end

function fix_FL!(OP)
  #  OP=OP_init
    Prof=Profile_Data[Ptype_]
    for i in eachindex(OP["Pinc"]) 
        if OP["Pinc"][i]>(Inc_coef*(Matrix(Prof["PL"])*(SysData["Load"]."Pd (pu)")'))[:,SysData["FL"]."D_N"][i]
            OP["Pinc"][i]=(Inc_coef*(Matrix(Prof["PL"])*(SysData["Load"]."Pd (pu)")'))[:,SysData["FL"]."D_N"][i]
        end
        if OP["Pdec"][i]>(Dec_coef*(Matrix(Prof["PL"])*(SysData["Load"]."Pd (pu)")'))[:,SysData["FL"]."D_N"][i]
            OP["Pdec"][i]=(Dec_coef*(Matrix(Prof["PL"])*(SysData["Load"]."Pd (pu)")'))[:,SysData["FL"]."D_N"][i]
        end
    end
    if Problem_Type in ["SCOPF", "RBOPF"]    
        si=size(OP["Pinc_c_Line"])
        for t=1:si[1],i=1:si[2],c=1:si[3]
            if OP["Pinc_c_Line"][t,i,c]>(Inc_coef*(Matrix(Prof["PL"])*(SysData["Load"]."Pd (pu)")'))[:,SysData["FL"]."D_N"][t,i]
                OP["Pinc_c_Line"][t,i,c]=(Inc_coef*(Matrix(Prof["PL"])*(SysData["Load"]."Pd (pu)")'))[:,SysData["FL"]."D_N"][t,i]
            end
            if OP["Pdec_c_Line"][t,i,c]>(Dec_coef*(Matrix(Prof["PL"])*(SysData["Load"]."Pd (pu)")'))[:,SysData["FL"]."D_N"][t,i]
                OP["Pdec_c_Line"][t,i,c]=(Dec_coef*(Matrix(Prof["PL"])*(SysData["Load"]."Pd (pu)")'))[:,SysData["FL"]."D_N"][t,i]
            end
        end
    end
    return OP
end

function Prepare_Dataset(OP::Dict{String,Any},TSysData::Dict{String, DataFrame},t;i=0,C="N",I_relax=1) #Array{Union{Float64,Array{Float64}}}
    # OP=OP_NValid
    # OP=OP_Valid
    NSysData=deepcopy(TSysData)
    T=select(NSysData["Trans"],Not("g_sh (pu)"));
    T."g_sh (pu)"=NSysData["Trans"]."g_sh (pu)" #bringing shunt conductance to last column
    NSysData["Trans"]=T
    NSysData["Line"]=NSysData["Line"][NSysData["Line"]."br_status".==1,:]
    NSysData["ESS"]=NSysData["ESS"][NSysData["ESS"]."br_status".==1,:]
    NSysData["Trans"]=NSysData["Trans"][NSysData["Trans"]."br_status".==1,:]
    Prof=Profile_Data[Ptype_]

    NSysData["RES"]."PG"=convert.(Float64,NSysData["RES"]."PG")
    NSysData["ESS"]."PG"=convert.(Float64,NSysData["ESS"]."PG")
    if isempty(Prof["Wind"]) NSysData["RES"]."PG".=0.0 else NSysData["RES"]."PG".= (Matrix(Prof["Wind"])*(NSysData["RES"]."Cap")')[t,:] end;
    if isempty(Prof["PL"]) NSysData["Load"]."Pd (pu)".=0.0 else NSysData["Load"]."Pd (pu)".= (Matrix(Prof["PL"])*(NSysData["Load"]."Pd (pu)")')[t,:] end;
    if isempty(Prof["QL"]) NSysData["Load"]."Qd (pu)".=0.0 else NSysData["Load"]."Qd (pu)".= (Matrix(Prof["QL"])*(NSysData["Load"]."Qd (pu)")')[t,:] end;

    if C=="N"
        NSysData["Gen"]."Pg (pu)".=     OP["PG"][t,:];NSysData["Gen"]."Qg (pu)".=OP["QG"][t,:];
        NSysData["Gen"]."V_set (pu)".=  OP["V"][t,:][SysData["Gen"]."Bus"];
        ### to respect ramp rate
        # ramp=NSysData["Gen"]."Ramp".*Td
        # if t!=1
        #     rampmax=OP["PG"][t-1,:]+ramp
        #     rampmin=OP["PG"][t-1,:]-ramp
        #     NSysData["Gen"]."P_max (pu)"=min.(NSysData["Gen"]."P_max (pu)",rampmax)
        #     NSysData["Gen"]."P_min (pu)"=max.(NSysData["Gen"]."P_min (pu)",rampmin)
        # end

        ###

        if isempty(OP["Pshed"]) Shed=   OP["Pshed"] else Shed=OP["Pshed"][t,:] end;
        if isempty(OP["Pcurt"]) Curt=   OP["Pcurt"] else Curt=OP["Pcurt"][t,:] end;
        if isempty(OP["Pcurtc"]) Curtc= OP["Pcurtc"] else Curtc=OP["Pcurtc"][t,:] end

        # if isempty(OP["Pshed"]) Shed=   OP["Pshed"] else OP["Pshed"][t,:]=min.(OP["Pshed"][t,:],NSysData["Load"]."Pd (pu)");Shed=OP["Pshed"][t,:] end;
        # if isempty(OP["Pcurt"]) Curt=   OP["Pcurt"] else OP["Pcurt"][t,:]=min.(OP["Pcurt"][t,:],NSysData["RES"]."PG");Curt=OP["Pcurt"][t,:] end;
        # if isempty(OP["Pcurtc"]) Curtc= OP["Pcurtc"] else OP["Pcurtc"][t,:]=min.(OP["Pcurtc"][t,:],NSysData["Gen"]."Pg (pu)");Curtc=OP["Pcurtc"][t,:] end

        if isempty(OP["Pinc"])  Pinc=   OP["Pinc"]  else Pinc=OP["Pinc"][t,:] end;
        if isempty(OP["Pdec"])  Pdec=   OP["Pdec"]  else Pdec=OP["Pdec"][t,:] end;
        if isempty(OP["Pdis"])  Pdis=   OP["Pdis"]  else Pdis=OP["Pdis"][t,:] end;
        if isempty(OP["Pch"])   Pch=    OP["Pch"]   else Pch=OP["Pch"][t,:] end;
    elseif C=="ConL"
        # NSysData["Line"]=filter(x-> x["Line_N"]!=TSysData["ConL"][i,"#"] , TSysData["Line"])
        line_outage= findall(x-> x==NSysData["ConL"][i,"#"] , NSysData["Line"]."Line_N")
        NSysData["Line"]=NSysData["Line"][Not(line_outage),:]

        NSysData["Gen"]."Pg (pu)".=     OP["PG_c_Line"][t,:,i];NSysData["Gen"]."Qg (pu)".=OP["QG_c_Line"][t,:,i];
        NSysData["Gen"]."V_set (pu)".=  OP["V_c_Line"][t,:,i][SysData["Gen"]."Bus"];

        # ### to respect ramp rate
        # ramp=NSysData["Gen"]."Ramp".*CCT
        # rampmax=OP["PG"][t,:]+ramp
        # rampmin=OP["PG"][t,:]-ramp

        # NSysData["Gen"]."P_max (pu)"=min.(NSysData["Gen"]."P_max (pu)",rampmax)
        # NSysData["Gen"]."P_min (pu)"=max.(NSysData["Gen"]."P_min (pu)",rampmin)
        ###


        if isempty(OP["Pshed_c_Line"]) Shed=   OP["Pshed_c_Line"] else Shed=OP["Pshed_c_Line"][t,:,i] end;
        if isempty(OP["Pcurt_c_Line"]) Curt=   OP["Pcurt_c_Line"] else Curt=OP["Pcurt_c_Line"][t,:,i] end;
        if isempty(OP["Pcurtc_c_Line"]) Curtc=  OP["Pcurtc_c_Line"] else Curtc=OP["Pcurtc_c_Line"][t,:,i] end;
        # if isempty(OP["Pshed_c_Line"]) Shed=   OP["Pshed_c_Line"] else OP["Pshed_c_Line"][t,:,i]=min.(OP["Pshed_c_Line"][t,:,i],NSysData["Load"]."Pd (pu)");Shed=OP["Pshed_c_Line"][t,:,i] end;
        # if isempty(OP["Pcurt_c_Line"]) Curt=   OP["Pcurt_c_Line"] else OP["Pcurt_c_Line"][t,:,i]=min.(OP["Pcurt_c_Line"][t,:,i],NSysData["RES"]."PG");Curt=OP["Pcurt_c_Line"][t,:,i] end;
        # if isempty(OP["Pcurtc_c_Line"]) Curtc=  OP["Pcurtc_c_Line"] else OP["Pcurtc_c_Line"][t,:,i]=min.(OP["Pcurtc_c_Line"][t,:,i],NSysData["Gen"]."Pg (pu)");Curtc=OP["Pcurtc_c_Line"][t,:,i] end;
        if isempty(OP["Pinc_c_Line"])  Pinc=   OP["Pinc_c_Line"]  else Pinc=OP["Pinc_c_Line"][t,:,i] end;
        if isempty(OP["Pdec_c_Line"])  Pdec=   OP["Pdec_c_Line"]  else Pdec=OP["Pdec_c_Line"][t,:,i] end;
        if isempty(OP["Pdis_c_Line"])  Pdis=   OP["Pdis_c_Line"]  else Pdis=OP["Pdis_c_Line"][t,:,i] end;
        if isempty(OP["Pch_c_Line"])   Pch=    OP["Pch_c_Line"]   else Pch=OP["Pch_c_Line"][t,:,i] end;
    elseif C=="CASCADE"
        line_outage= findall(x-> x==NSysData["ConL"][i,"#"] , NSysData["Line"]."Line_N")
        NSysData["Line"]=NSysData["Line"][Not(line_outage),:]

        NSysData["Line"]."RATE_A"*=I_relax

        NSysData["Gen"]."Pg (pu)".=     OP["PG"][t,:];NSysData["Gen"]."Qg (pu)".=OP["QG"][t,:];
        NSysData["Gen"]."V_set (pu)".=  OP["V"][t,:][SysData["Gen"]."Bus"];

        if isempty(OP["Pshed"]) Shed=   OP["Pshed"] else Shed=OP["Pshed"][t,:] end;
        if isempty(OP["Pcurt"]) Curt=   OP["Pcurt"] else Curt=OP["Pcurt"][t,:] end;
        if isempty(OP["Pcurtc"]) Curtc= OP["Pcurtc"] else Curtc=OP["Pcurtc"][t,:] end
        if isempty(OP["Pdis"])  Pdis=   OP["Pdis"]  else Pdis=OP["Pdis"][t,:] end;
        if isempty(OP["Pch"])   Pch=    OP["Pch"]   else Pch=OP["Pch"][t,:] end;
        if isempty(OP["Pinc"])  Pinc=   OP["Pinc"]  else Pinc=OP["Pinc"][t,:] end;
        if isempty(OP["Pdec"])  Pdec=   OP["Pdec"]  else Pdec=OP["Pdec"][t,:] end;
    end

    if !isempty(Pinc) NSysData["Load"]."Pd (pu)"[SysData["FL"]."D_N"].+= ((Pinc)+ (Pdec*-1)) end;
    if !isempty(Pdis) NSysData["ESS"]."PG".= Pdis+ (-1* Pch) end;
    if !isempty(Curt) NSysData["RES"]."PG".+= (-1*Curt); end
    if !isempty(Shed) NSysData["Load"]."Pd (pu)".+= -1*Shed; end
    if !isempty(Shed) NSysData["Load"]."Qd (pu)".+= -1*(Shed.*NSysData["Load"]."PF"); end
    if !isempty(Curtc) NSysData["Gen"]."Pg (pu)".-= Curtc; end

    return NSysData
end
function pre_prepration!(OP,NSysData;nCon=0,nT=0,C="N")
    Prof=Profile_Data[Ptype_]
    if !isempty(Prof["Wind"])  Ren_T= (Matrix(Prof["Wind"])*(NSysData["RES"]."Cap")') end;
    if !isempty(Prof["PL"])  Ld_T= (Matrix(Prof["PL"])*(NSysData["Load"]."Pd (pu)")') end;
    for t=1:nT
        Ren=Ren_T[t,:]
        Ld=Ld_T[t,:]
        if C=="N"
            PG= OP["PG"][t,:];
            if !isempty(OP["Pshed"]) OP["Pshed"][t,:]=min.(OP["Pshed"][t,:],Ld) end;
            if !isempty(OP["Pcurt"]) OP["Pcurt"][t,:]=min.(OP["Pcurt"][t,:],Ren) end;
            if !isempty(OP["Pcurtc"]) OP["Pcurtc"][t,:]=min.(OP["Pcurtc"][t,:],PG) end
        elseif C=="ConL"
            for i=1:nCon
                PG= OP["PG_c_Line"][t,:,i];
                if !isempty(OP["Pshed_c_Line"]) OP["Pshed_c_Line"][t,:,i]=min.(OP["Pshed_c_Line"][t,:,i],Ld) end;
                if !isempty(OP["Pcurt_c_Line"]) OP["Pcurt_c_Line"][t,:,i]=min.(OP["Pcurt_c_Line"][t,:,i],Ren) end;
                if !isempty(OP["Pcurtc_c_Line"]) OP["Pcurtc_c_Line"][t,:,i]=min.(OP["Pcurtc_c_Line"][t,:,i],PG) end;
            end
        end
    end
    return OP
end

function SLA(Runcode,OP_init0,init_model,init_out;BItermax=4)
    #BItermax=7;
    iter=0;
    BIter=0;
    sol_time=0;
    res_time=0;
    Best_Iter=0;
    Best_sol_time=0;
    Best_res_time=0;
    #----------------------------Check initial operating point -------------------------------------------
    Ncost0,Vio_Sum_N0,Vio_general_N0,AVG_C_cost0,Vio_Sum_C0,OP_init0V,Vio_general_C0,min_infmax_C=0,0,0,0,0,0,0,0
    try
        Ncost0,Vio_Sum_N0,OP_init0_NV0,Vio_general_N0=PF_Valid_N(OP_init0,Runcode)
            # @descend PF_Valid_N(OP_init0)
            # ProfileView.@profview PF_Valid_N(OP_init0)
        # findall(x-> x>=0.01, OP_init0_NV0["QG"]-OP_init0["QG"])
        if Problem_Type in ["SCOPF", "RBOPF"]
            AVG_C_cost0,Vio_Sum_C0,OP_init0V,Vio_general_C0=PF_Valid_C(OP_init0_NV0,Runcode)
        elseif Problem_Type=="OPF"
            OP_init0V=OP_init0_NV0
            Vio_Sum_C0=DataFrame(Nins=0,NVio=0,Mean=0,Max=0,Type="")
            AVG_C_cost0=Ncost0
            min_infmax_C=0
        end
    catch e
        @warn "PF not converged"
        if Linearization_TYPE[Runcode[2]] in ["SENS","Yang"]
            @info "Switch to DC for one iteration"
            while true
                try OP_init0,out,model=Optimizer(Linearization_TYPE[1],OP_init0)
                catch e
                    if relaxpermit==true
                        rethrow(e)
                    else
                        relaxpermit=true
                        continue
                    end
                end
            end
            relaxpermit=false

            # OP_init0=Initialize("Exact",out,model,[],[])
            Ncost0,Vio_Sum_N0,OP_init0_NV0,Vio_general_N0=PF_Valid_N(OP_init0,Runcode)
            if Problem_Type in ["SCOPF", "RBOPF"]
                AVG_C_cost0,Vio_Sum_C0,OP_init0V,Vio_general_C0=PF_Valid_C(OP_init0_NV0,Runcode)
            elseif Problem_Type=="OPF"
                OP_init0V=OP_init0_NV0
                Vio_Sum_C0=DataFrame(Nins=0,NVio=0,Mean=0,Max=0,Type="")
                AVG_C_cost0=Ncost0
                min_infmax_C=0
            end
        else
            rethrow(e)
        end
    end
    Best_solution= OP_init0V
    Best_vio_general_N=Vio_general_N0; Best_vio_sum_N=Vio_Sum_N0

    cost0=Ncost0.+(AVG_C_cost0.-Ncost0)*conprob_
    mincost=sum(cost0)
    minOptGap= (sum(cost0)-value(init_model[:real_objective]))/(value(init_model[:real_objective]))
    min_infmax_N=maximum(Vio_Sum_N0."Max")
    if Problem_Type in ["SCOPF", "RBOPF"] Best_vio_general_C=Vio_general_C0;Best_vio_sum_C=Vio_Sum_C0;min_infmax_C=maximum(Vio_Sum_C0."Max")
    else min_infmax_C=0;Best_vio_general_C=0;Best_vio_sum_C=0 end


    OP_init=deepcopy(OP_init0V)
    OP=0;
    V_1,θ_1,PG_1,QG_1,cost1=Matrix(init_out."mag_v"[1]),Matrix(init_out."ang_v"[1]),Matrix(init_out."PG"[1]),Matrix(init_out."QG"[1]),init_out."OFValue_WP"[1]
    if Problem_Type in ["SCOPF", "RBOPF"]
        V_c_Line_1,θ_c_Line_1,PG_c_Line_1,QG_c_Line_1=sqrt.(value.(init_model[:e_c_Line].^2 +init_model[:f_c_Line].^2)),atan.(value.(init_model[:f_c_Line]),value.(init_model[:e_c_Line])),value.(init_model[:PG_c_Line]),value.(init_model[:QG_c_Line])
    end
    if FullReport==1
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"TR";"TRB";"OptGap";"Max_Inf_N";"Max_Inf_C";"Optimizer_Obj";"real_Obj";["PG$i" for i in 1:size(SysData["Gen"],1)];["QG$i" for i in 1:size(SysData["Gen"],1)];["V$i" for i in 1:size(SysData["Bus"],1)];["θ$i" for i in 1:size(SysData["Bus"],1)]]
        Iter_Rep=DataFrame(repeat([[]],length(Head)), Head)
        #------------------- Print Optimal point and first iteration
        newrow=["-" "Non-Linear" 0 solve_time(init_model) 0 1 0 0 0 0 0 cost1 "-" PG_1 QG_1 V_1 θ_1]
    elseif Problem_Type in ["SCOPF", "RBOPF"]
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"TR";"TRB";"OptGap";"Max_Inf_N";"Max_Inf_C";"Optimizer_Obj";"real_Obj";"PG_RMSE(p.u)";"QG_RMSE(p.u)";"V_RMSE(p.u)";"θ_RMSE(Rad)";"PG_c_RMSE(p.u)";"QG_c_RMSE(p.u)";"V_c_RMSE(p.u)";"θ_c_RMSE(Rad)";]
        Iter_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0 solve_time(init_model) 0 1 0 0 0 0 0 "-" init_out."OFValue_WP"[1]  0 0 0 0 0 0 0 0]
    elseif Problem_Type=="OPF"
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"TR";"TRB";"OptGap";"Max_Inf_N";"Max_Inf_C";"Optimizer_Obj";"real_Obj";"PG_RMSE(p.u)";"QG_RMSE(p.u)";"V_RMSE(p.u)";"θ_RMSE(Rad)";]
        Iter_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0 solve_time(init_model) 0 1 0 0 0 0 0  "-"  init_out."OFValue_WP"[1] 0 0 0 0]
    end
    push!(Iter_Rep,newrow[1,:])
    if FullReport==1
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) 0 0 iter 0 0 minOptGap min_infmax_N min_infmax_C init_out."OFValue_WP"[1] round(mincost,digits=5) round.(Best_solution["PG"],digits=5) round.(Best_solution["QG"],digits=5) round.(Best_solution["V"],digits=5) round.(Best_solution["θ"],digits=5)]
    elseif Problem_Type in ["SCOPF", "RBOPF"]
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) 0 0 iter 0 0 minOptGap min_infmax_N min_infmax_C init_out."OFValue_WP"[1] round(mincost,digits=5) sqrt(mean((round.(Best_solution["PG"],digits=5).-PG_1).^2)) sqrt(mean((round.(Best_solution["QG"],digits=5).-QG_1).^2)) sqrt(mean((round.(Best_solution["V"],digits=5).-V_1).^2)) sqrt(mean((round.(Best_solution["θ"],digits=5).-θ_1).^2)) sqrt(mean((round.(Best_solution["PG_c_Line"],digits=5).-PG_c_Line_1).^2)) sqrt(mean((round.(Best_solution["QG_c_Line"],digits=5).-QG_c_Line_1).^2)) sqrt(mean((round.(Best_solution["V_c_Line"],digits=5).-V_c_Line_1).^2)) sqrt(mean((round.(Best_solution["θ_c_Line"],digits=5).-θ_c_Line_1).^2))]
    elseif Problem_Type == "OPF"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) 0 0 iter 0 0 minOptGap min_infmax_N min_infmax_C init_out."OFValue_WP"[1] round(mincost,digits=5) sqrt(mean((round.(Best_solution["PG"],digits=5).-PG_1).^2)) sqrt(mean((round.(Best_solution["QG"],digits=5).-QG_1).^2)) sqrt(mean((round.(Best_solution["V"],digits=5).-V_1).^2)) sqrt(mean((round.(Best_solution["θ"],digits=5).-θ_1).^2))]
    end
    push!(Iter_Rep,newrow[1,:])

    if Runcode[1]==2
        Best_vio_general_C=0; Best_vio_sum_C=0
        mincost=1e100
        minOptGap=100
        min_infmax_N=100
        if Problem_Type in ["SCOPF", "RBOPF"] min_infmax_C=1e100 else min_infmax_C=0 end
    end
    #--------------------------------------------------
    # global DP=OP_init0_NV0["Pinj"].-OP_init["Pinj"]
    # global DQ=OP_init0_NV0["Qinj"].-OP_init["Qinj"]


    #New
    Vio_general_C=Vio_general_C0
    Vio_general_N=Vio_general_N0

    while iter <30
        iter=iter+1
        OP,out,model=0,0,0
        # OP,out,model=Optimizer(Linearization_TYPE[3],OP_init)
        # model[:NEactive_power_balance][1,57]
        # show_constraints_summary( ,model)
        # Violations=(Vio_general_N,Vio_general_C)
        # OP,out,model=Optimizer(Linearization_TYPE[Runcode[2]],OP_init,Violations)
        while true
            # OP,out,model=Optimizer(Linearization_TYPE[Runcode[2]],OP_init)
                # OP,out,model=Optimizer(Linearization_TYPE[3],OP_init)
            try OP,out,model=Optimizer(Linearization_TYPE[Runcode[2]],OP_init)
            catch e
                # @warn "Expanding trust region"
                # # global TRfree=(1+TRfree)*1.1-1
                # global TRfree+=1
                # # if TRfree>10
                # #     @error "Too expanded trust region"
                # #     rethrow(e)
                # # end
                if relaxpermit==true
                    rethrow(e)
                    break
                end
                @warn "There is no feasible solution in the trust region"
                @info "Relaxing the technical constraints!"
                global relaxpermit=true
                BIter=0
                continue
            end
            if relaxpermit==true res_time+=solve_time(model) else sol_time+=solve_time(model) end
            if termination_status(model)==MOI.OPTIMAL
                if  relaxpermit==true && Linearization_TYPE[Runcode[2]]=="SENS"
                    if (Problem_Type=="OPF" && value.(sum(model[:slk_Imax])+sum(model[:slk_Vmax])+sum(model[:slk_Vmin])+sum(model[:slk_NNramp1])+sum(model[:slk_NNramp2]))<1e-3) ||(Problem_Type in ["SCOPF", "RBOPF"] && value.(sum(model[:slk_Imax])+sum(model[:slk_Vmax])+sum(model[:slk_Vmin])+sum(model[:slk_NNramp1])+sum(model[:slk_NNramp2])+sum(model[:slk_Imax_c])+ sum(model[:slk_Vmax_c])+sum(model[:slk_Vmin_c])+sum(model[:slk_CNramp1])+sum(model[:slk_CNramp2]))<1e-3) #
                    #value.(sum(model[:slk_Imax])+sum(model[:slk_Vmax])+sum(model[:slk_Vmin])+sum(model[:slk_NNramp1])+sum(model[:slk_NNramp2]))
                    #value.(sum(model[:slk_Imax_c])+ sum(model[:slk_Vmax_c])+sum(model[:slk_Vmin_c])+sum(model[:slk_CNramp1])+sum(model[:slk_CNramp2]))
                    #value.(sum(model[:slk_Imax_c]))+ value.(sum(model[:slk_Vmax_c]))+ value.(sum(model[:slk_Vmin_c]))+ value.(sum(model[:slk_CNramp1]))+ value.(sum(model[:slk_CNramp2]))
                        @info "Restoration is compeleted!"
                        global relaxpermit=false;BIter=0
                    end
                end
                break
            else @warn "Optimizer is not converged"; @error "termination status: $(termination_status(model))"; end
        end


        # A= OP_init0["PG"];B= OP_init["PG"];C= OP["PG"];
        # A= OP_init0["QG"] ;B= OP_init["QG"];C= OP["QG"]


        Ncost,Vio_Sum_N,OP_NValid,Vio_general_N,convergence_indicator=PF_Valid_N(OP,Runcode)
        # if convergence_indicator !=1 @error("PF is not converging anymore iteration $iter"); break; end
        if Problem_Type=="OPF"
            OP_Valid=OP_NValid
            AVG_C_cost=Ncost
            Vio_general_C=0
            Vio_Sum_C=DataFrame(Nins=0,NVio=0,Mean=0,Max=0,Type="")
        elseif Problem_Type in ["SCOPF", "RBOPF"]
            AVG_C_cost,Vio_Sum_C,OP_Valid,Vio_general_C=PF_Valid_C(OP_NValid,Runcode)
        end


        cost=Ncost.+(AVG_C_cost.-Ncost)*conprob_
        OptGap= (sum(cost)-value(init_model[:real_objective]))/(value(init_model[:real_objective]))
        Maxvio_C=maximum(Vio_Sum_C."Max")
        Maxvio_N=maximum(Vio_Sum_N."Max")
        #print out
        if FullReport==1
            newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time res_time iter OP_init["TR"] OP_init["TRB"] OptGap Maxvio_N Maxvio_C out."OFValue_WP"[1] round(sum(cost),digits=5) round.(OP_Valid["PG"],digits=5) round.(OP_Valid["QG"],digits=5) round.(OP_Valid["V"],digits=5) round.(OP_Valid["θ"],digits=5)]
        elseif Problem_Type in ["SCOPF", "RBOPF"]
            newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time res_time iter OP_init["TR"] OP_init["TRB"] OptGap Maxvio_N Maxvio_C out."OFValue_WP"[1] round(sum(cost),digits=5) sqrt(mean((round.(OP_Valid["PG"],digits=5).-PG_1).^2)) sqrt(mean((round.(OP_Valid["QG"],digits=5).-QG_1).^2)) sqrt(mean((round.(OP_Valid["V"],digits=5).-V_1).^2)) sqrt(mean((round.(OP_Valid["θ"],digits=5).-θ_1).^2)) sqrt(mean((round.(OP_Valid["PG_c_Line"],digits=5).-PG_c_Line_1).^2)) sqrt(mean((round.(OP_Valid["QG_c_Line"],digits=5).-QG_c_Line_1).^2)) sqrt(mean((round.(OP_Valid["V_c_Line"],digits=5).-V_c_Line_1).^2)) sqrt(mean((round.(OP_Valid["θ_c_Line"],digits=5).-θ_c_Line_1).^2))]
        elseif Problem_Type == "OPF"
            newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time res_time iter OP_init["TR"] OP_init["TRB"] OptGap Maxvio_N Maxvio_C out."OFValue_WP"[1] round(sum(cost),digits=5) sqrt(mean((round.(OP_Valid["PG"],digits=5).-PG_1).^2)) sqrt(mean((round.(OP_Valid["QG"],digits=5).-QG_1).^2)) sqrt(mean((round.(OP_Valid["V"],digits=5).-V_1).^2)) sqrt(mean((round.(OP_Valid["θ"],digits=5).-θ_1).^2))]
        end
        push!(Iter_Rep,newrow[1,:])
        #########
        if  (Maxvio_N <= max(min_infmax_N,1e-4) && ((1e-4<= min_infmax_N-Maxvio_N ) || (Maxvio_C <= max(min_infmax_C,1e-4) &&((1e-4<= min_infmax_C-Maxvio_C )||((1e-4<minOptGap-OptGap)&&(OptGap>10e-6))))))
            for k in collect(keys(OP_Valid))
                OP_init[k]=OP_Valid[k]
            end
        else
            BIter+=1
            if Linearization_TYPE[Runcode[2]]=="SENS"
                if relaxpermit==false
                    OP_init["TRB"]=[Adaptive_TRB[1]]
                    if BIter>=size(Adaptive_TR,1)
                        @error "No more iterations needed"
                        break;
                    end
                    OP_init["TR"]=[Adaptive_TR[BIter+1]]
                else
                    if BIter>=size(Adaptive_TRB,1)
                        @error "No feasible solution found"
                        break;
                    end
                    OP_init["TRB"]=[Adaptive_TRB[BIter+1]]
                end
            else
                for k in collect(keys(OP_Valid))
                    OP_init[k]=OP_Valid[k]
                end
                if BIter>=BItermax
                    @error "No more iterations needed"
                    break;
                end
            end
        end
        if Maxvio_N <= max(min_infmax_N,1e-4) && ((1e-6<= min_infmax_N-Maxvio_N ) || (Maxvio_C <= max(min_infmax_C,1e-4) &&((1e-6<= min_infmax_C-Maxvio_C )||(1e-4<mincost-sum(cost)))))
            # (1e-3<= min_infmax_N-Maxvio_N ) || (1e-3<= min_infmax_C-Maxvio_C )||(1e-1<mincost-sum(cost))
            mincost=sum(cost)
            minOptGap=OptGap
            min_infmax_C=Maxvio_C
            min_infmax_N=Maxvio_N
            Best_solution= OP_Valid
            Best_vio_general_C=Vio_general_C
            Best_vio_sum_C=Vio_Sum_C
            Best_vio_general_N=Vio_general_N
            Best_vio_sum_N=Vio_Sum_N
            Best_Iter=iter
            Best_sol_time=sol_time
            Best_res_time=res_time
            # global TRfree= TRfree/2
        end
        # Iter_Rep,BIter#,value.(sum(model[:slk_CNramp1])+sum(model[:slk_CNramp2]))
    end
    push!(Iter_Rep,repeat(["end"],size(Iter_Rep,2)))
    ods_write(Path*"../Result/Iter_Report/"*Problem_Type*string(Runcode),Dict(("General",1,1)=>Iter_Rep))
    return mincost,Best_solution,Best_Iter,min_infmax_N,min_infmax_C,minOptGap,Best_vio_general_N,Best_vio_sum_N,Best_vio_general_C,Best_vio_sum_C,Best_sol_time,Best_res_time
end

function SLA_clean(Runcode,OP_init0,OP_NLP;BItermax=4)
    #BItermax=1;
    iter=0;BIter=0;sol_time=0;res_time=0;Best_Iter=0;Best_sol_time=0;Best_res_time=0;LNSs=0;
    #----------------------------Check initial operating point -------------------------------------------
    Ncost0,Vio_Sum_N0,Vio_general_N0,AVG_C_cost0,Vio_Sum_C0,OP_init0V,Vio_general_C0,min_infmax_C,casreps0,casreps_NLP=0,0,0,0,0,0,0,0,0,0
    try
        Ncost0,Vio_Sum_N0,Vio_general_N0,AVG_C_cost0,Vio_Sum_C0,OP_init0V,Vio_general_C0,casreps0=AC_PFCheck(OP_init0,Runcode);

        a0,a1,a2,a3,a4,a5,a6,casreps_NLP=AC_PFCheck(OP_NLP,Runcode);


        # LNSs=[casreps[t,i][end,"LNS_T"] for t in 1:nTL, i in 1:size(SysData["ConL"],1)]
        # findall(x-> size(x,1)>1, casreps)
        # casreps[1,5]
        # Vio_general_C0[1,1]
        GC.gc()
    catch e
        @warn "PF not converged"
        if Linearization_TYPE[Runcode[2]] in ["SENS","Yang"]
            @info "Switch to DC for one iteration"
            while true
                try OP_init0,out,model=Optimizer(Linearization_TYPE[1],OP_init0)
                catch e
                    if relaxpermit==true rethrow(e) else relaxpermit=true;continue end
                end
                break;
            end
            relaxpermit=false
            Ncost0,Vio_Sum_N0,Vio_general_N0,AVG_C_cost0,Vio_Sum_C0,OP_init0V,Vio_general_C0, casreps=AC_PFCheck(OP_init0,Runcode)
        else rethrow(e) end
    end
    Best_solution= OP_init0V;Best_vio_general_N=Vio_general_N0; Best_vio_sum_N=Vio_Sum_N0


    cost0=Ncost0.+(AVG_C_cost0.-Ncost0)*conprob_
    #sum(Ncost0)
    #sum((AVG_C_cost0.-Ncost0))
    mincost=sum(cost0);minOptGap= 100*((sum(cost0)-OP_NLP["obj"])/(OP_NLP["obj"]))
    min_infmax_N=maximum(Vio_Sum_N0."Max")
    if Problem_Type in ["SCOPF", "RBOPF"] Best_vio_general_C=Vio_general_C0;Best_vio_sum_C=Vio_Sum_C0;min_infmax_C=maximum(Vio_Sum_C0."Max")
    else min_infmax_C=0;Best_vio_general_C=0;Best_vio_sum_C=0 end

    OP_init=deepcopy(OP_init0V);OP=0;
    casreps=deepcopy(casreps0)
    LNS_cas0,LNSs0=LNS_cal(casreps0)
    LNS_cas_NLP,LNSs_NLP=LNS_cal(casreps_NLP)
 
    Iter_Rep=Make_Iter_report(OP_NLP,LNS_cas_NLP)
    # Vio_general_N=Best_vio_general_N;Vio_general_C=Best_vio_general_C;OP=Best_solution
    # Best_solution["acc"]=[0.0 0.0]
    Iter_Rep=Add_Iter_report!(Runcode,Iter_Rep,Best_solution,mincost,minOptGap,Best_vio_sum_N,Best_vio_sum_C,Best_vio_general_N,Best_vio_general_C,iter,LNS_cas0,0,0,[[0 0] [0 0] [0 0] [0 0]])
    # if Runcode[1]==2
    #     Best_vio_general_C=0; Best_vio_sum_C=0
    #     mincost=1e100;minOptGap=10000;min_infmax_N=100
    #     if Problem_Type in ["SCOPF", "RBOPF"] min_infmax_C=1e100 else min_infmax_C=0 end
    # end
    #New
    Vio_general_C=Vio_general_C0;Vio_general_N=Vio_general_N0

    # Vio_general_C0[6,2]
    # sum(sum.([sum.(Vio_general_N[t,1]."percent") for t=1:nTP]))
    # Normal_infeasibility=0;
    # Contingency_infeasibility=0
    OP_Valid=Dict{String, Any}();cost=0;OptGap=0;Maxvio_C=0;Maxvio_N=0;Vio_Sum_C=0;Vio_Sum_N=0
    R_P=0;R_Q=0
    OP_init["TRB"][1]=maximum(sum(OP_init["PG"],dims=2))*0.0001

    OP,out,model=0,0,0
    LNS_cas=LNS_cas0

    while iter <10 && res_time<=600
        iter=iter+1;
        while true
            try
                if Runcode[2]!=5
                    if Cascade_Opt
                        OP,out,model=Optimizer(Linearization_TYPE[Runcode[2]],deepcopy(OP_init),casreps)
                    else
                        OP,out,model=Optimizer(Linearization_TYPE[Runcode[2]],deepcopy(OP_init))
                    end
                    # value.(model[:slack_Cascade_LP])
                    # maximum(value.(model[:slack_Cascade_LP]))

                elseif iter<2
                    OP,out,model=Optimizer(Linearization_TYPE[1],deepcopy(OP_init))
                else
                    if Cascade_Opt
                        OP,out,model=Optimizer(Linearization_TYPE[3],deepcopy(OP_init),casreps)
                    else
                        OP,out,model=Optimizer(Linearization_TYPE[3],deepcopy(OP_init))
                    end
                end

                # println("slk_normal:",value.(model[:slk_Norm]))
                # println("slk_cont:",value.(model[:slk_Cont]))

                # maximum(value.(model[:ΔPinjP]))
                # maximum(value.(model[:ΔPinjN]))
                # maximum(value.(model[:ΔQinjP]))
                # maximum(value.(model[:ΔQinjN]))
                # OP["TR"]
                OP["Pshed"]=max.(0,OP["Pshed"])
                OP["Pcurt"]=max.(0,OP["Pcurt"])
                OP["Pcurtc"]=max.(0,OP["Pcurtc"])
                OP["Pshed_c_Line"]=max.(0,OP["Pshed_c_Line"])
                OP["Pcurt_c_Line"]=max.(0,OP["Pcurt_c_Line"])
                OP["Pcurtc_c_Line"]=max.(0,OP["Pcurtc_c_Line"])

                # minimum(OP["Pshed"])
                # findmax(OP["V_c_Line"])
                # OP["V_c_Line"][7,20,31]

                # model[:ΔPinjP_c_Line][14,:,9]
                # OP_init["PG_c_Line"][14,11,9]
                # OP_init["QG_c_Line"][14,11,9]
                # OP["PG_c_Line"][14,11,9]
                # OP["QG_c_Line"][14,11,9]

                # findmin(OP["V_c_Line"])
                # value.(model[:CV2magexporg][10,14,21])
                # OP_init["V_c_Line"][8,43,12]

                # OP_Valid["V_c_Line"][7,50,31]
                # OP_init["V_c_Line"][7,50,31]
                # OP_NLP["V_c_Line"][7,50,31]
                # OP_Valid["V_c_Line"][7,20,31]
                # OP_init["V_c_Line"][7,20,31]

                # OP["V_c_Line"][15,14,20]



                if Linearization_TYPE[Runcode[2]]=="SENS"|| (Linearization_TYPE[Runcode[2]]=="Hybrid" && iter>1)
                    print( DataFrame(
                                obj=value(objective_function(model)),
                                real_obj=value.(model[:real_objective]),
                                TR_relax=value.(model[:cost_TR_relaxN]),
                                TR_relaxC=value.(model[:cost_TR_relaxC]),
                                TR_PenN=value.(model[:cost_TR_PenN]),
                                TR_PenC=value.(model[:cost_TR_PenC])
                                ))
                    R_P=[maximum(sum(value.(model[:ΔPinjP]+model[:ΔPinjN]),dims=2)),minimum(sum(value.(model[:ΔPinjP]+model[:ΔPinjN]),dims=2))]
                    R_Q=[maximum(sum(value.(model[:ΔQinjP]+model[:ΔQinjN]),dims=2)),minimum(sum(value.(model[:ΔQinjP]+model[:ΔQinjN]),dims=2))]
                    # sum(value.(model[:ΔPinjP]-model[:ΔPinjN]),dims=2)
                    # sum(value.(model[:ΔQinjP]-model[:ΔQinjN]),dims=2)
                else
                    R_P,R_Q=[0,0],[0,0]
                end
                println("R_P:",R_P)
                println("R_Q:",R_Q)
            catch e
                if relaxpermit==true 
                    # rethrow(e);
                    # Iter_Rep=Add_Iter_report!(Runcode,Iter_Rep,OP_Valid,sum(cost),OptGap,Vio_Sum_N,Vio_Sum_C,Vio_general_N,Vio_general_C,iter,sol_time,res_time;Convergence=false)
                    break 
                end
                @warn "There is no feasible solution in the trust region";@info "Relaxing the technical constraints!"
                global relaxpermit=true;continue
            end
            if relaxpermit==true res_time+=solve_time(model) else sol_time+=solve_time(model) end
            if termination_status(model)==MOI.OPTIMAL
                if  relaxpermit==true && (Linearization_TYPE[Runcode[2]]=="SENS"|| (Linearization_TYPE[Runcode[2]]=="Hybrid" && iter>1))
                    Normal_infeasibility=value.(sumJuMPi(model[:slk_Imax])+sumJuMPi(model[:slk_Vmax])+sumJuMPi(model[:slk_Vmin])+sumJuMPi(model[:slk_NNramp]))
                    println(Normal_infeasibility)
                    # Norm_Wrmap=value.(sumJuMPi(model[:slk_Imax])+sumJuMPi(model[:slk_Vmax])+sumJuMPi(model[:slk_Vmin]))
                    # Cont_Wramp=value.(sumJuMPi(model[:slk_Imax_c])+ sumJuMPi(model[:slk_Vmax_c])+sumJuMPi(model[:slk_Vmin_c]))

                    if Problem_Type in ["SCOPF", "RBOPF"] Contingency_infeasibility=value.(sumJuMPi(model[:slk_Imax_c])+ sumJuMPi(model[:slk_Vmax_c])+sumJuMPi(model[:slk_Vmin_c])+sumJuMPi(model[:slk_CNramp]))

                        println("slk_normal_P:",value.(model[:slk_Norm]))
                        println("slk_cont_P:",value.(model[:slk_Cont]))
                        println("slk_normal_Q:",value.(model[:slk_Norm1]))
                        println("slk_cont_Q:",value.(model[:slk_Cont1]))
                        if Cascade_Opt
                            println("slk_cascade_LP:",maximum(value.(model[:slack_Cascade_LP]);init=0.0))
                        end

               
                    # value.(model[:slk_Vmin])
                    # println("slk_normal:",value.(model[:slk_Norm]))
                    # println("slk_cont:",value.(model[:slk_Cont]))
                    # model[:slkDataFrame]["slkVmax_c"][557]
                    # CartesianIndex(10,14,21) in model[:slkDataFrame]["slkVmax_c"]
                    # value.(model[:slk_Imax_c])
                    # SysData["Bus"][14,5:6].=[1.1,0.9]
                    # model[:CV2magexporg][7,50,31]
                    # value.(model[:slk_Cont])
                    # findmax(value.(model[:slk_Vmax_c]))
                    # value.(model[:slk_Vmin_c])[63]
                    # model[:I_c_Line][24,3,6]
                    # model[:Line_flow_c_Line][24,3,6]


                    # xxx=sort(value.(model[:slk_Vmax_c]))
                    # Plots.plot(xxx)
                    # xxx=findall(x-> x>0.5, value.(model[:slk_Vmax_c]))
                    # model[:W_VIO]
                    # OP["acc"]

                    else Contingency_infeasibility=0 end
                    println(Contingency_infeasibility)
                    # if Problem_Type in ["SCOPF", "RBOPF"] t_infes=Normal_infeasibility+Contingency_infeasibility else t_infes=Normal_infeasibility end
                    # if t_infes<1e-6
                    # @info "Restoration is compeleted!";global relaxpermit=true;
                    # end

                    # if max(value.(model[:slk_Norm]),value.(model[:slk_Cont]))<=0.001
                    #     freeV=0
                    # end
                end
                break
            else @warn "Optimizer is not converged"; @error "termination status: $(termination_status(model))"; end
        end
        try
            Ncost,Vio_Sum_N,Vio_general_N,AVG_C_cost,Vio_Sum_C,OP_Valid,Vio_general_C,casreps=AC_PFCheck(OP,Runcode)
            LNS_cas,LNSs=LNS_cal(casreps)
            # LNSs=[casreps[t,i][end,"LNS_T"] for t in 1:nTL, i in 1:size(SysData["ConL"],1)]
            # findall(x-> size(x,1)>1, casreps)
            # OP["Pshed"][:,1872]
            # Vio_general_C[17,10]
            # Vio_general_C[6,26]."abs"
            # Vio_general_N[14,1]
            # Vio_general_C[19,10]

            GC.gc()
            cost= Ncost.+(AVG_C_cost.-Ncost)*conprob_
            #sum(Ncost)
            OptGap= 100*((sum(cost)-OP_NLP["obj"])/(OP_NLP["obj"]))
            Maxvio_C=maximum(Vio_Sum_C."Max");Maxvio_N=maximum(Vio_Sum_N."Max")
            OP["TR"]=model[:modTR]
            # OP["acc"]
            if Linearization_TYPE[Runcode[2]]=="SENS"
                Filter=[[model[:Num_filter_Const_I] model[:Num_Const_I]],
                        [model[:Num_filter_Const_V] model[:Num_Const_V]],
                        [model[:SENS_filter_term_I] model[:SENS_term_I]] ,
                        [model[:SENS_filter_term_V] model[:SENS_term_V]]]
            else
                Filter=[[0 0] [0 0] [0 0] [0 0]]
            end

                    # Tuple(Filter[1])
            Iter_Rep=Add_Iter_report!(Runcode,Iter_Rep,OP_Valid,sum(cost),OptGap,Vio_Sum_N,Vio_Sum_C,Vio_general_N,Vio_general_C,iter,LNS_cas,sol_time,res_time,Filter,R_P,R_Q)
            display(Iter_Rep)
            display(Iter_Rep[:,15:end])

            # findmax(OP["V_c_Line"].-(SysData["Bus"]."V_max (pu)")')
            # OP["V_c_Line"][1,43,10]-(SysData["Bus"]."V_max (pu)")[43]
            # value.(model[:CV2magexporg][18,41,1])
            # model[:I_c_Line][22,46,30]
            # value.(model[:Line_flow_c_Line][7,30,18])
            # model[:Line_flow][14,1]
            # model[:I_cal_c_Line][7,30,18]
            # print(model[:Line_flow_c_Line][22,46,30])
                            # aa= calculate_I_sol( model[:e_c_Line],model[:f_c_Line],raw)
                            # aa[22,46,30]
                            # model
                            # OP_Valid["V_c_Line"][22,:,30]
                            # OP["V_c_Line"][22,:,30]

            # raw["QActBus"]
            # SysData["Line"]."RATE_A"[57].^2
            # maximum(-OP["V"].+(SysData["Bus"]."V_min (pu)")')
        catch e
            # rethrow(e)
            @error "PF not converged"
            cost=0
            OptGap= 100
            Maxvio_C=100;Maxvio_N=100

            OP_Valid=OP
            Iter_Rep = Add_Iter_report!(Runcode,Iter_Rep,OP_Valid,sum(cost),OptGap,Vio_Sum_N,Vio_Sum_C,Vio_general_N,Vio_general_C,iter,LNS_cas,sol_time,res_time;Convergence=false)
        end
               #########
        if true || (max(Maxvio_N,Maxvio_C) <= max(min_infmax_N,min_infmax_C,1e-4) && ((1e-4<= max(min_infmax_N,min_infmax_C)-max(Maxvio_N,Maxvio_C)) ||((1e-2<minOptGap-OptGap)&&(OptGap>10e-4))))
            for k in collect(keys(OP_Valid))
                OP_init[k]=OP_Valid[k]
            end 
        else
            BIter+=1
            if Linearization_TYPE[Runcode[2]]=="SENS"
                # if BIter>=2 #size(Adaptive_TR,1)
                    @error "No more iterations needed"
                    break;
                # end
                # OP_init["TR"]=Adaptive_TR[BIter+1]
                # OP_init["FilterIndex"]=[filter_factor[BIter+1]]
            else
                for k in collect(keys(OP_Valid))
                    OP_init[k]=OP_Valid[k]
                end
                if BIter>=BItermax || OptGap==100
                    @error "No more iterations needed"
                    break;
                end
            end
        end
        if (max(Maxvio_N,Maxvio_C) <= max(min_infmax_N,min_infmax_C,1e-4) && ((1e-6<= max(min_infmax_N,min_infmax_C)-max(Maxvio_N,Maxvio_C)) ||((1e-4<minOptGap-OptGap)&&(OptGap>10e-6))))
            # (1e-3<= min_infmax_N-Maxvio_N ) || (1e-3<= min_infmax_C-Maxvio_C )||(1e-1<mincost-sum(cost))
            mincost=sum(cost);minOptGap=OptGap
            min_infmax_C=Maxvio_C;min_infmax_N=Maxvio_N
            Best_solution= OP_Valid
            Best_vio_general_C=Vio_general_C;Best_vio_sum_C=Vio_Sum_C
            Best_vio_general_N=Vio_general_N;Best_vio_sum_N=Vio_Sum_N
            Best_Iter=iter
            Best_sol_time=sol_time;Best_res_time=res_time
        end
        if max(min_infmax_N,min_infmax_C,LNS_cas)<1e-4 || ( Cascade_Opt && max(min_infmax_N,min_infmax_C)<1e-4) 
            @info "Violations are removed"
            break
        end
    end
    push!(Iter_Rep,repeat(["end"],size(Iter_Rep,2)))
    # Iter_Rep2=deepcopy(Iter_Rep)
    # Iter_Rep[1:7,:]
    ods_write(Path*"../Result/Iter_Report/"*Problem_Type*string(Runcode)*"_0",Dict(("General",1,1)=>Iter_Rep))
    # ods_write(Path*"../Result/Iter_Report/"*Problem_Type*string(Runcode),Dict(("General",1,1)=>Iter_Rep))
    return mincost,Best_solution,Best_Iter,min_infmax_N,min_infmax_C,minOptGap,Best_vio_general_N,Best_vio_sum_N,Best_vio_general_C,Best_vio_sum_C,Best_sol_time,Best_res_time
end

function SLP_Report_Prepare(init_model,init_out)
    sol_time1=solve_time(init_model)
    V_1,θ_1,PG_1,QG_1,cost1=Matrix(init_out."mag_v"[1]),Matrix(init_out."ang_v"[1]),Matrix(init_out."PG"[1]),Matrix(init_out."QG"[1]),init_out."OFValue_WP"[1]
    if Problem_Type in ["SCOPF", "RBOPF"]
        V_c_Line_1,θ_c_Line_1,PG_c_Line_1,QG_c_Line_1=sqrt.(value.(init_model[:e_c_Line].^2 +init_model[:f_c_Line].^2)),atan.(value.(init_model[:f_c_Line]),value.(init_model[:e_c_Line])),value.(init_model[:PG_c_Line]),value.(init_model[:QG_c_Line])
    end
    if Report =="feasiblity" || Report=="Florin"
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"Max_N_redisP";"Max_N_redisQ";"Max_C_redisP";"Max_C_redisQ";"OptGap";"Max_Inf_N";"Max_Inf_C";"Obj";"Smax_N";"Vmax_N";"Vmin_N";"Ramp_N";"Smax_C";"Vmax_C";"Vmin_C";"Ramp_C"]
        SLP_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0 sol_time1 0 1 "-" "-" "-" "-" 0 "-" "-" cost1 "-" "-" "-" "-" "-" "-" "-" "-" ]
    elseif Report =="full"
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"OptGap";"Max_Inf_N";"Max_Inf_C";"Obj";["PG$i" for i in 1:size(SysData["Gen"],1)];["QG$i" for i in 1:size(SysData["Gen"],1)];["V$i" for i in 1:size(SysData["Bus"],1)];["θ$i" for i in 1:size(SysData["Bus"],1)]]
        SLP_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0  sol_time1 0 1 0 0 0 cost1 PG_1 QG_1 V_1 θ_1]
    elseif Problem_Type in ["SCOPF", "RBOPF"]
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"OptGap";"Max_Inf_N";"Max_Inf_C";"Obj";"PG_Error";"QG_Error";"V_Error";"θ_Error";"PG_c_Error";"QG_c_Error";"V_c_Error";"θ_c_Error"]
        SLP_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0 sol_time1 0 1 0 0 0 init_out."OFValue_WP"[1] 0 0 0 0 0 0 0 0]
    elseif Problem_Type=="OPF"
        Head=["Initialization";"Method";"Runcode";"Time";"Res_Time";"Iter";"OptGap";"Max_Inf_N";"Max_Inf_C";"Obj";"PG_Error";"QG_Error";"V_Error";"θ_Error"]
        SLP_Rep=DataFrame(repeat([[]],length(Head)), Head)
        newrow=["-" "Non-Linear" 0 sol_time1 0 1 0 0 0 init_out."OFValue_WP"[1] 0 0 0 0]
    end
    push!(SLP_Rep,newrow[1,:])
    return SLP_Rep
end

function Add_SLP_report!(SLP_Rep,Runcode,Best_solution,mincost,minOptGap,min_infmax_N,min_infmax_C,Best_vio_general_N,Best_vio_general_C,iter,sol_time,restoration_time,Max_N_redis,Max_C_redis,Max_N_redisQ,Max_C_redisQ)

    max_N = maximum(hcat([df[!, :max] for df in Best_vio_general_N]...),dims=2)
    if Problem_Type in ["SCOPF", "RBOPF"]
        max_C = maximum(hcat([df[!, :max] for df in Best_vio_general_C]...),dims=2)
    else
        max_C=[0.0,0.0, 0.0, 0.0]
    end

    if Report=="feasiblity" || Report=="Florin"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time restoration_time iter Max_N_redis Max_N_redisQ Max_C_redis  Max_C_redisQ minOptGap min_infmax_N min_infmax_C round(mincost,digits=5) max_N[3] max_N[2] max_N[1] max_N[4] max_C[3] max_C[2] max_C[1] max_C[4]]
    elseif Report=="full"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time restoration_time iter minOptGap min_infmax_N min_infmax_C round(mincost,digits=5) round.(Best_solution["PG"],digits=5) round.(Best_solution["QG"],digits=5) round.(Best_solution["V"],digits=5) round.(Best_solution["θ"],digits=5)]
    elseif Problem_Type == "SCOPF"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time restoration_time iter minOptGap min_infmax_N min_infmax_C round(mincost,digits=5) sqrt(sum((round.(Best_solution["PG"],digits=5).-PG_1).^2)) sqrt(sum((round.(Best_solution["QG"],digits=5).-QG_1).^2)) sqrt(sum((round.(Best_solution["V"],digits=5).-V_1).^2)) sqrt(sum((round.(Best_solution["θ"],digits=5).-θ_1).^2)) sqrt(sum((round.(Best_solution["PG_c_Line"],digits=5).-PG_c_Line_1).^2)) sqrt(sum((round.(Best_solution["QG_c_Line"],digits=5).-QG_c_Line_1).^2)) sqrt(sum((round.(Best_solution["V_c_Line"],digits=5).-V_c_Line_1).^2)) sqrt(sum((round.(Best_solution["θ_c_Line"],digits=5).-θ_c_Line_1).^2))]
    elseif Problem_Type =="OPF"
        newrow=[Linearization_INIT[Runcode[1]] Linearization_TYPE[Runcode[2]]  Tuple(Runcode) sol_time restoration_time iter minOptGap min_infmax_N min_infmax_C round(mincost,digits=5) sqrt(sum((round.(Best_solution["PG"],digits=5).-PG_1).^2)) sqrt(sum((round.(Best_solution["QG"],digits=5).-QG_1).^2)) sqrt(sum((round.(Best_solution["V"],digits=5).-V_1).^2)) sqrt(sum((round.(Best_solution["θ"],digits=5).-θ_1).^2))]
    end
    push!(SLP_Rep,newrow[1,:],promote=true)

    return SLP_Rep
end


function Initial_Set_points()
    #exact non-linear solution
    if Problem_Type=="OPF"
        init_out,init_model=DA_OPF()
    elseif Problem_Type=="SCOPF"
        init_out,init_model=DA_SCOPF()
    elseif Problem_Type=="RBOPF"
        init_out,init_model=DA_SCOPF(;relaxconfac=0.2)
    end

    #DA solution with error
    DA_out,DA_model=Dict(),Dict()
    for i in sort(filter(x-> x[1:2]=="DA",collect(values(Linearization_INIT))))
        global ERORR["DA"]=[0.01 0.01]*parse(Int, i[3:end])
        Set_global_parameters(Test_Sys)
        if Problem_Type=="OPF"
            DA_out[i],DA_model[i]=DA_OPF()
        elseif Problem_Type=="SCOPF"
            DA_out[i],DA_model[i]=DA_SCOPF()
        elseif Problem_Type=="RBOPF"
            DA_out[i],DA_model[i]=DA_SCOPF(;relaxconfac=0.2)
        end
    end
    global ERORR["DA"]=[0.00 0.00]
    Set_global_parameters(Test_Sys)
    return DA_out,DA_model,init_out,init_model
end
function Partials(Jacob,V,isl,nBus)
    # Jac1=deepcopy(Jacob)
    Jac1=Jacob
    @. Jac1[:,nBus+1:end] /= V'
    Jcbfactorized1 = lu(Jac1, check=false)
    Sensivities=inv(Jcbfactorized1)
    Sensivities[isl,:].=0
    return Sensivities
end

function Set_TR(SENSDict_All,I2_0,V_0,raw;cl=1)
    F2maxtemp=@. raw["Line"]."RATE_A"^2
    Imaxlimit0=(F2maxtemp')[ones(Int8,nTP),:,ones(Int8,cl)]
    Vminlimit0=(raw["Bus"]."V_min (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
    Vmaxlimit0=(raw["Bus"]."V_max (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]

    ovr_Sens_IP=zeros(nTP,raw["nLine"],cl)
    ovr_Sens_IQ=zeros(nTP,raw["nLine"],cl)
    ovr_Sens_IP_worst=zeros(nTP,raw["nLine"],cl)
    ovr_Sens_IQ_worst=zeros(nTP,raw["nLine"],cl)

    # ovr_Sens_IP_N=zeros(nTP,raw["nLine"],cl)
    # ovr_Sens_IQ_P=zeros(nTP,raw["nLine"],cl)
    # ovr_Sens_IP_P=zeros(nTP,raw["nLine"],cl)
    # ovr_Sens_IQ_N=zeros(nTP,raw["nLine"],cl)
    # filter(x->x>0,SENSDict_All["S_IP"][t,c][i,:])
    for i=1:raw["nLine"],t=1:nTP,c=1:cl
        # # ovr_Sens[t,i,c]=sum(abs.(SENSDict["S_IQ"][t,c][i,:]))*Lininit["TR"][2] +sum(abs.(SENSDict["S_IP"][t,c][i,:]))*Lininit["TR"][1]
        ovr_Sens_IP[t,i,c]=sum(abs.(SENSDict_All["S_IP"][t,c][i,:]))+0.000001
        ovr_Sens_IQ[t,i,c]=sum(abs.(SENSDict_All["S_IQ"][t,c][i,:]))+0.000001
        # ovr_Sens_IP[t,i,c]=minimum([sum(filter(x->x>0,SENSDict_All["S_IP"][t,c][i,:])),-sum(filter(x->x<0,SENSDict_All["S_IP"][t,c][i,:]))])+0.000001
        # ovr_Sens_IQ[t,i,c]=minimum([sum(filter(x->x>0,SENSDict_All["S_IQ"][t,c][i,:])),-sum(filter(x->x<0,SENSDict_All["S_IQ"][t,c][i,:]))])+0.000001
        # ovr_Sens_IP_worst[t,i,c]=abs(sum(filter(x->x>0,SENSDict_All["S_IP"][t,c][i,:]))+sum(filter(x->x<0,SENSDict_All["S_IP"][t,c][i,:])))+0.000001
        # ovr_Sens_IQ_worst[t,i,c]=abs(sum(filter(x->x>0,SENSDict_All["S_IQ"][t,c][i,:]))+sum(filter(x->x<0,SENSDict_All["S_IQ"][t,c][i,:])))+0.000001
        # # ovr_Sens_IQ[t,i,c]=abs(sum(filter(x->x>0,SENSDict_All["S_IQ"][t,c][i,:]))+sum(filter(x->x<0,SENSDict_All["S_IQ"][t,c][i,:])))+0.000001
        # ovr_Sens_IP_P[t,i,c]=sum(filter(x->x>0,SENSDict_All["S_IP"][t,c][i,:]))+0.000001
        # ovr_Sens_IP_N[t,i,c]=sum(filter(x->x<0,SENSDict_All["S_IP"][t,c][i,:]))-0.000001
        # ovr_Sens_IQ_P[t,i,c]=sum(filter(x->x>0,SENSDict_All["S_IQ"][t,c][i,:]))+0.000001
        # ovr_Sens_IQ_N[t,i,c]=sum(filter(x->x<0,SENSDict_All["S_IQ"][t,c][i,:]))-0.000001
    end
    TR_P=minimum([(round(10*maximum(max.((I2_0.-Imaxlimit0),0)./ovr_Sens_IP),digits=6)+0.000001),2])

    ovr_Sens_VP=zeros(nTP,raw["nBus"],cl)
    ovr_Sens_VQ=zeros(nTP,raw["nBus"],cl)
    for i=1:raw["nBus"],t=1:nTP,c=1:cl
        ovr_Sens_VP[t,i,c]=sum(abs.(SENSDict_All["S_VP"][t,c][i,:]))+0.000001
        ovr_Sens_VQ[t,i,c]=sum(abs.(SENSDict_All["S_VQ"][t,c][i,:]))+0.000001
        # ovr_Sens_VP[t,i,c]=abs(sum(filter(x->x>0,SENSDict_All["S_VP"][t,c][i,:]))+sum(filter(x->x<0,SENSDict_All["S_VP"][t,c][i,:])))+0.000001
        # ovr_Sens_VQ[t,i,c]=abs(sum(filter(x->x>0,SENSDict_All["S_VQ"][t,c][i,:]))+sum(filter(x->x<0,SENSDict_All["S_VQ"][t,c][i,:])))+0.000001
        # # ovr_Sens_VQ[t,i,c]=abs(SENSDict_All["S_VQ"][t,c][i,:])+0.000001
    end
    TR_Q=minimum([(round(1.1*maximum(max.((V_0.-Vmaxlimit0),(Vminlimit0.-V_0),0)./ovr_Sens_VQ),digits=6)+0.000001),2])
    return [TR_P,TR_Q],[ovr_Sens_IP,ovr_Sens_IQ,ovr_Sens_VP,ovr_Sens_VQ]
end

function Set_TR_new(Base_TR_N,Base_TR_C,Lininit,SENSDict_All, I2_0, V_0, raw,PBound,QBound; cl=1,Improve_Rate=10)
    # Precompute limits
    F2maxtemp = (raw["Line"]."RATE_A") .^ 2
    Imaxlimit0=(F2maxtemp')[ones(Int8,nTP),:,ones(Int8,cl)]
    Vminlimit0=(raw["Bus"]."V_min (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
    Vmaxlimit0=(raw["Bus"]."V_max (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]

    # I2_0[3,2265,3]
    # Imaxlimit0[3,2265,3]
    # Allocate output arrays
    # ovr_Sens_IP = zeros(nTP, raw["nLine"], cl)
    # ovr_Sens_IQ = zeros(nTP, raw["nLine"], cl)
    # ovr_Sens_VP = zeros(nTP, raw["nBus"], cl)
    # ovr_Sens_VQ = zeros(nTP, raw["nBus"], cl)

    # plot(reshape(ovr_Sens_VQ,1,:)[1,:],title="Sensitivity of V",label="to Q")
    # plot!(reshape(ovr_Sens_VP,1,:)[1,:],linestyle=:dash,label="to P",linewidth=0.05,linealpha=0.5)
    # savefig("Sensitivity_V.png")
    # ovr_Sens_IP[:,:,1]
    # IP_1=filter_to_sparse(ovr_Sens_IP[:,:,1]./10e4, x-> x>=0.0001)
    # findall(x-> x>=0.1, ovr_Sens_IP)
    # # plot(reshape(ovr_Sens_IP./10e4,1,:)[1,:],title="Sensitivity of I",label="to P",linewidth=0.1,yscale=:log10)
    # plot!(reshape(ovr_Sens_IQ./10e4,1,:)[1,:],linestyle=:dash,label="to Q",linewidth=0.05,linealpha=1,yscale=:log10)
    # savefig("Sensitivity_I.png")
    # plot(reshape(ovr_Sens_IP,1,:)[1,:])
    # plot!(reshape(ovr_Sens_IQ,1,:)[1,:])
    # plot(reshape(IP_1,1,:)[1,:],title="Sensitivity of I",label="to P",linewidth=0.1,yscale=:log2)


    # IP=Array{Vector{Int}}(undef, nTP,raw["nLine"], cl)

    # Extract sensitivity data
    S_IP = SENSDict_All["S_IP"]
    S_IQ = SENSDict_All["S_IQ"]
    S_VP = SENSDict_All["S_VP"]
    S_VQ = SENSDict_All["S_VQ"]



    # plot(S_VQ[1,1][17,:])
    # findall(x-> abs(x)>=0.02,S_VQ[1,1][17,:])

    # mean(abs.(S_VQ[1,1][17,:]))

    # raw["QActBus"][findall(x-> abs(x)>=0.02,S_VQ[1,1][17,:])]
    # Compute ovr_Sens_VP and ovr_Sens_VQ (bus sensitivities)
    # for c in 1:cl, t in 1:nTP, i in 1:raw["nBus"]
    #     @views SVP_tc = S_VP[t, c][i, :]
    #     @views SVQ_tc = S_VQ[t, c][i, :]

    #     ovr_Sens_VP[t, i, c] = sum(abs.(SVP_tc)) + 1e-6
    #     ovr_Sens_VQ[t, i, c] = sum(abs.(SVQ_tc)) + 1e-6
    # end

    # Compute TR_Q
    # TR_Q = minimum([round(1.1 * maximum(max.((V_0 .- Vmaxlimit0), (Vminlimit0 .- V_0), 0) ./ ovr_Sens_VQ), digits=6) + 1e-6, 2])
    # TR_Q,VmaxIndex,VminIndex,acc_V,Qinjlimit,Qworst_vio,report_TR_Q=findTR_Q(V_0,Vmaxlimit0,Vminlimit0,S_VQ,QBound,raw,cl)


    # CartesianIndex(5,2470,1) in VminIndex
    # if max(maximum(Base_TR_N),maximum(Base_TR_C))>0.5
    #     println("Ramp rate limit violations are excessive")
    #     # TR_Q=0;
    #     empty!(VmaxIndex)
    #     empty!(VminIndex)
    #     # acc_V=1
    #     # Qinjlimit[1].=0
    #     # Qinjlimit[2].=0
    # end
    # CartesianIndex(6,183,1) in VmaxIndex
    # CartesianIndex(6,2528,3) in VmaxIndex
    # CartesianIndex(8,2740,7) in VmaxIndex
    # V_0[8,2740,7]
    # Vmaxlimit0[8,2740,7]

    # TR_Q= max(TR_Q,2,Lininit["TR"][2]*0.9)
    # Compute ovr_Sens_IP and ovr_Sens_IQ (line sensitivities)
    # for c in 1:cl, t in 1:nTP, i in 1:raw["nLine"]
    #     @views SIP_tc = S_IP[t, c][i, :]

    #     @views SIQ_tc = S_IQ[t, c][i, :]

    #     ovr_Sens_IP[t, i, c] = sum(abs.(SIP_tc)) + 1e-6
    #     ovr_Sens_IQ[t, i, c] = sum(abs.(SIQ_tc)) + 1e-6
    # end

    # sens_IQ=Worst_sens_IQ(Qinjlimit,S_IQ,raw,cl)

    # # ramp = abs.(Lininit["PG"][2:end, :] - Lininit["PG"][1:end-1, :])
    # # ramp_c=abs.(Lininit["PG_c_Line"] .- Lininit["PG"])
    # # ramp_vio=max(maximum(ramp .- (raw["Gen"]."Ramp")'),maximum(ramp_c.- (raw["Gen"]."Ramp")'))
    # TR_P,ImaxIndex,acc_I,Pinjlimit,P2worst_vio,Pworst_vio,report_TR_P=findTR_P(Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,sens_IQ,PBound,raw,cl)

    # ImaxIndex_Q=findIcoupled(TR_Q,I2_0,Imaxlimit0,S_IQ,QBound,raw,cl,acc_I)
    # VmaxIndex_P,VminIndex_P=findVcoupled(TR_P,V_0,Vmaxlimit0,Vminlimit0,S_VP,PBound,raw,cl,acc_V)

 
    
        # TR_P,ImaxIndex,VmaxIndex_P,VminIndex_P,acc_I,Pinjlimit,P2worst_vio,Pworst_vio,report_TR_P,Vmax_VP,Vmin_VP=findTR_P_coupled(Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,V_0,Vmaxlimit0,Vminlimit0,S_VP,PBound,raw,cl,acc_V)
    TR_P,ImaxIndex,VmaxIndex_P,VminIndex_P,acc_I,Pinjlimit,P2worst_vio,Pworst_vio,report_TR_P,Vmax_VP,Vmin_VP,Imax_IP = findTR_P_coupled( Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,V_0,Vmaxlimit0,Vminlimit0,deepcopy(S_VP),PBound,raw,cl,0.05)
    GC.gc()
    # TR_P=0.5
    # VmaxIndex_P[15,14,20]
    # Vmax_VP[17,32,20] 
    # V_0[17,32,]    
    # Vmin_VP
    # S_VP[18,2][41,1:10]

    # Vmax_VP[18,41,2]
    # sum(V_0.>Vmax_VP)
    # sum(V_0.+0.00001.<Vmin_VP)

    TR_Q,VmaxIndex,VminIndex,acc_V,Qinjlimit,Qworst_vio,report_TR_Q= findTR_Q_coupled(V_0,Vmax_VP,Vmin_VP,VmaxIndex_P,VminIndex_P,Vmaxlimit0,Vminlimit0,deepcopy(S_VQ),QBound,raw,cl)
    GC.gc()
    # CartesianIndex(1,41,14) in VmaxIndex
    # TR_Q,VmaxIndex,VminIndex,acc_V,Qinjlimit,Qworst_vio,report_TR_Q=findTR_Q_coupled(V_0,V_0,V_0,VmaxIndex_P,VminIndex_P,Vmaxlimit0,Vminlimit0,deepcopy(S_VQ),QBound,raw,cl)

    # ImaxIndex_Q=findIcoupled(TR_Q,I2_0,Imax_IP,ImaxIndex,Imaxlimit0,S_IQ,QBound,raw,cl,acc_I)

    # Imax_IP[14,2,]
    # raw["Line"]

    #  BSON.@save Path*"../Result/TR_figure_result.bson" report_TR_Q


    # ImaxIndext=union(ImaxIndex,ImaxIndex_Q)
    VmaxIndext=union(VmaxIndex,VmaxIndex_P)
    VminIndext=union(VminIndex,VminIndex_P)

    # TR_P,ImaxIndex,acc_I,Pinjlimit,P2worst_vio,Pworst_vio,report_TR_P=findTR_P_coupled(Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,V_0,Vmaxlimit0,Vminlimit0,S_VP,PBound,raw,cl)

 

    ImaxIndext=ImaxIndex
    # VmaxIndext=VmaxIndex
    # VminIndext=VminIndex

    # TR_act=[ max.(abs.(Pinjlimit[1]),abs.(Pinjlimit[2])) ,max.(abs.(Qinjlimit[1]),abs.(Qinjlimit[2]))]

    Injlimit=[Pinjlimit,Qinjlimit]
    # TR_P = minimum([round(2.5 * maximum(max.((I2_0 .- Imaxlimit0), 0) ./ ovr_Sens_IP), digits=6) + 1e-6, 2])

    return [TR_P, TR_Q], [ImaxIndext, VmaxIndext, VminIndext],[acc_I,acc_V],Injlimit,[P2worst_vio,Qworst_vio,Pworst_vio]

end

function Set_TR_new2(Base_TR_N,Base_TR_C,SENSDict_All, I2_0, V_0, raw,PBound,QBound,CascadeInfo; cl=1,Improve_Rate=10)
    # Precompute limits
    F2maxtemp = (raw["Line"]."RATE_A") .^ 2
    Imaxlimit0=(F2maxtemp')[ones(Int8,nTP),:,ones(Int8,cl)]
    Vminlimit0=(raw["Bus"]."V_min (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]
    Vmaxlimit0=(raw["Bus"]."V_max (pu)")'[ones(Int8,nTP),:,ones(Int8,cl)]

    # Extract sensitivity data
    S_IP = SENSDict_All["S_IP"]
    S_IQ = SENSDict_All["S_IQ"]
    S_VP = SENSDict_All["S_VP"]
    S_VQ = SENSDict_All["S_VQ"]


    TR_P,ImaxIndex,VmaxIndex_P,VminIndex_P,acc_I,Pinjlimit,P2worst_vio,Pworst_vio,report_TR_P,Vmax_VP,Vmin_VP,Imax_IP = findTR_P_coupled2( CascadeInfo, Base_TR_N,Base_TR_C,I2_0,Imaxlimit0,S_IP,V_0,Vmaxlimit0,Vminlimit0,deepcopy(S_VP),PBound,raw,cl,0.05)
    TR_Q,VmaxIndex,VminIndex,acc_V,Qinjlimit,Qworst_vio,report_TR_Q= findTR_Q_coupled2(V_0,Vmax_VP,Vmin_VP,VmaxIndex_P,VminIndex_P,Vmaxlimit0,Vminlimit0,deepcopy(S_VQ),QBound,raw,cl)
    # GC.gc()
    VmaxIndext=union(VmaxIndex,VmaxIndex_P)
    VminIndext=union(VminIndex,VminIndex_P)
    ImaxIndext=ImaxIndex
 
    Injlimit=[Pinjlimit,Qinjlimit]
    return [TR_P, TR_Q], [ImaxIndext, VmaxIndext, VminIndext],[acc_I,acc_V],Injlimit,[P2worst_vio,Qworst_vio,Pworst_vio]
end

function find_state_variables(LSCOPF,Lininit,raw,Mode="Norm")
    isl=findall(x-> x==3,raw["Bus"]."Bus_Type")[1]
    if Mode=="Norm"
        Jacob=Lininit["NJacob"]
        V=Lininit["V"]
        θ=Lininit["θ"]
        pinj=value.(LSCOPF[:ΔPinjP]-LSCOPF[:ΔPinjN])
        qinj=value.(LSCOPF[:ΔQinjP]-LSCOPF[:ΔQinjN])
        cl=1
    elseif Mode=="Cont"
        Jacob=Lininit["CJacob"]
        V=Lininit["V_c_Line"]
        θ=Lininit["θ_c_Line"]
        pinj=value.(LSCOPF[:ΔPinjP_c_Line]-LSCOPF[:ΔPinjN_c_Line])
        qinj=value.(LSCOPF[:ΔQinjP_c_Line]-LSCOPF[:ΔQinjN_c_Line])
        cl=raw["nConL"]
    end
    for t in 1:nTP
        for c in 1:cl
            injcompl=zeros(2*raw["nBus"])
            injcompl[raw["PActBus"]].=pinj[t,:,c]
            injcompl[raw["QActBus"].+raw["nBus"]].=qinj[t,:,c]
            Jcb=lu(Jacob[t,c])
            Δx=Jcb \ (injcompl)
            Δx[isl]=0
            θ[t,:,c]=θ[t,:,c]+ Δx[1:raw["nBus"]]
            V[t,:,c]=V[t,:,c]+ Δx[raw["nBus"]+1:end]
        end
    end
    return θ,V
    # Plots.plot(Δx[raw["nBus"]+1:end])
    # Δx[raw["nBus"]+15]
    # V[t,14,c]
    

end

# SENSDict_N["S_VQ"][1,1]
# Jcbinvers=Jcb\I

# Jcbinvers[raw["nBus"]+1:end,raw["QActBus"].+raw["nBus"]]

function con_indices(parameter,sens, maxnum,nTP,cl,sens_G=0,Gbus=0)
    #parameter,sens, maxnum,nTP,cl= (Imaxlimit0.-I2_0),ovr_Sens_I,minline,nTP,cl
    # parameter,sens, maxnum,nTP,cl,sens_G,Gbus=(Vmaxlimit0.-V_0),ovr_Sens_V, minbus,nTP,cl,ovr_Sens_V_G,GBus
    Gbus=0
    vio=parameter.-sens
    if Gbus!=0
        vioG= parameter.-sens_G
    end
    index_out=Vector{CartesianIndex{3}}()
    elements = Array{Vector{Int}}(undef, nTP, cl)  # 3×3 matrix where each element is a vector
    for t=1:nTP
        for c=1:cl
            poten=partialsortperm(vio[t,:,c], 1:min(sum(vio.<0),maxnum))
            viola=findall(parameter[t,:,c].<=0.01)
            violaG= (Gbus==0) ? [] : findall(vioG[t,:,c].<=0 .&& Bool.(Gbus))
            indices_=union(poten,viola,violaG)
            elements[t,c]=indices_
            carind=CartesianIndex.([(t,ind,c) for ind in indices_])
            append!(index_out,carind)
        end
    end
    return index_out,elements
end