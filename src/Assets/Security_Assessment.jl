function Contingency_Screen!(initial_point,Set_points,t0)

    global ϵ=1e-4
    Infes_Conv_rep=[[],[],[]]
    ##########-----------------Read Data-----------------##########
    Ptype_="RT";init_="DA_SCOPF";# CSP=Current Set Points: init_="DA_SCOPF" ==> DA_SCOPF Set points
    raw,CSP,profraw= try Data_Reader(;Ptype=Ptype_ , init=init_,t_0=t0)
    catch e
        @warn "Problem with running Data_Reader function"
        rethrow(e)
    end
    Prof=Load_Data(nTP,raw,profraw) # Prepare Load and Shunt data and generation setpoints
    #####------Set Initial point based on measurements-------######
    if(initial_point==1)
        set_keys=collect(keys(CSP))
        CSP=Dict()
        for x in set_keys
            CSP[x]=Set_points[:,x][1]
        end
    end
    ##########-----------------Intialization by zero-----------------##########

    vio = DataFrame([Symbol("Con$i") => [DataFrame() for _ in 1:nTP] for i in 1:size(raw["ConL"],1)])
    Stage2out = DataFrame([Symbol("Con$i") => [DataFrame() for _ in 1:nTP] for i in 1:size(raw["ConL"],1)])
    SA=DataFrame([Symbol("Con$i") => [DataFrame(Pre_vio=DataFrame(),OPF=DataFrame(),Post_vio=DataFrame()) for _ in 1:nTP] for i in 1:size(raw["ConL"],1)])
    ##########-----------------Main Part-----------------##########

    for i in ["N",eachrow(raw["ConL"])...]
        # L_n=L
        for t=1:nTP
            # Apply contingency
            raw_temp=deepcopy(raw)
            if i!="N"
                raw_temp["Line"]= filter(row->row."Line_N"!=i."#",raw["Line"])
                raw_temp["nLine"]-=1
            end
            # Power Flow
            ov,uv,lo,lo_index,vmag_i,vang_i,ind_converge,Iij,percent=PF(raw_temp,CSP,Prof,t)
            # Check Violations
            if i=="N" # Normal Operation Violations
                if size(ov)[1]+size(uv)[1]+size(lo)[1]!=0
                    Infes_Conv_rep[3]=[(t,[size(ov)[1] size(uv)[1] size(lo)[1]]),percent]
                elseif sum((Matrix(CSP["Vmag"])[t,:]'.-vmag_i').<ϵ)!=60
                    Infes_Conv_rep[3]=[(t,maximum((Matrix(CSP["Vmag"])[t,:]'.-vmag_i')))]
                elseif sum((Matrix(CSP["Vang"])[t,:]'.-vang_i').<ϵ)==60
                    Infes_Conv_rep[3]=[(t,maximum((Matrix(CSP["Vang"])[t,:]'.-vang_i')))]
                end
                if ERORR["RT"]==[0 0] #To Check Power Flow Performance
                    @assert size(ov)[1]+size(uv)[1]+size(lo)[1]==0
                    @assert sum((Matrix(CSP["Vmag"])[t,:]'.-vmag_i').<ϵ)==60
                    @assert sum((Matrix(CSP["Vang"])[t,:]'.-vang_i').<ϵ)==60
                end
            else
                SA[t,i[1]][1,"Pre_vio"]=DataFrame(Number=[size(ov)[1]+size(uv)[1]+size(lo)[1]],Over_Voltage=DataFrame([(ov[j],vmag_i[ov[j]]) for j in eachindex(ov)]),Under_Voltage=DataFrame([(uv[j],vmag_i[uv[j]]) for j in eachindex(uv)]),Line_Over=DataFrame([(lo[j],Iij[lo_index[j]],percent[lo_index[j]]) for j in eachindex(lo_index)]),ind=ind_converge)
                flag=1
                result=Post_Contingency_SP_OPF(raw_temp,CSP,Prof,t)
                SA[t,i[1]][1,"OPF"]=result[1]
                SA[t,i[1]][1,"Post_vio"]=DataFrame(Shed= sum(result[1]."Pshed"[1][1,:]), Curt=sum(result[1]."Pcurt"[1][1,:]),cost=result[1]."OFValue")
                SA[t,i[1]][1,"Post_vio"]."Problem".=sum(result[1]."Pshed"[1][1,:])<1e-4 ? 0 : 1
                if ind_converge!=1
                    SA[t,i[1]][1,"Post_vio"]."Problem".=1
                    # newline="\n t= $t \t contingency= $(i[1])"
                    # Infes_Conv_rep[1]*=newline
                    push!(Infes_Conv_rep[1],(t,i[1]))
                end

                if (Int(result[2])==3)
                    # newline="\n t= $t \t contingency= $(i[1])"
                    # Infes_Conv_rep[2]*=newline
                    push!(Infes_Conv_rep[2],(t,i[1]))
                    SA[t,i[1]][1,"Post_vio"]."Problem".=1
                end
            end

        end
    end

    # findall((x,y)-> vio[1,:][x]."ind"==0,1:nCon)

    # vioindex=findall(x-> ~isempty(x[1,4]) ,Matrix(vio))

    # Set([a[i][1] for i in eachindex(vioindex)])
    # Set([a[i][2] for i in eachindex(vioindex)])

    # global Test_Sys="Nordic/WD/"
    # global Path= "Data/"*Test_Sys
    # init="RT_OPF"
    # initdata=["QG","PG","Pinc","Pdec","Pcurt","Pshed","Pch","Pdis","soc","Tcost"]
    # initraw=Dict()
    # for x in initdata
    #     initraw[x]=ods_read(Path*init,sheetName=x, retType= "DataFrame")
    # end
    
    
    sum(SA[t,i][1,"Pre_vio"]."Number" for t in 1:nTP,i in 1:raw["nConL"])
    sum(SA[t,i][1,"Post_vio"]."Problem" for t in 1:nTP,i in 1:raw["nConL"])
    
    # findall(i->[1,i][1,"Post_vio"]."Problem"==1,SA)
    
    tuples = [(t, i) for t in 1:nTP for i in 1:raw["nConL"] if SA[t,i][1,"Post_vio"]."Problem"[1] == 1]
    
    
    SA[1,14][1,"Post_vio"]
    
    # global ttmp=()
    # global itmp=()
    
    # for (t,i) in tuples
    #     ttmp=(ttmp...,t)
    #     itmp=(itmp...,i)
    # end
    Cost_Post_t1=sum(max(0,SA[t,i[1]][1,"Post_vio"]."cost"[1]) for t in 1:1,i in 1:raw["nConL"] )/raw["nConL"]
    
    ttmp = [t for (t, i) in tuples]
    itmp = [i for (t, i) in tuples]
    itmp1 = [i for (t, i) in tuples if t==1]
    
    cr_t=Set(ttmp)
    cr_con=Set(itmp)
    cr_con1=Set(itmp1)
    
    return cr_con,cr_con1,Cost_Post_t1,Infes_Conv_rep
    end