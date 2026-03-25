
# export Get_Result,Get_ConResult

function Get_Result(SCOPF,raw,initraw;OPFon=0)

    if (OPFon==0 && raw["nConL"]!=0)
        C_rep=DataFrame(Con_N=sort(findall(x -> x>-1e10,value.(SCOPF[:cost_Cor_C])),rev=true))
        C_rep."index"=raw["ConL"]."cont"[C_rep.Con_N,]
        C_rep."cost"=round.(value.(SCOPF[:cost_Cor_C]),digits=4)[C_rep.Con_N]
        C_rep."cost_CS"=round.(value.(SCOPF[:cost_CS_C]),digits=4)[C_rep.Con_N]
        C_rep."cost_R"=round.(sum(value.(SCOPF[:cost_Flex_C]),dims=1),digits=4)[C_rep.Con_N]
        C_rep."cost_Gen"=round.(sum(value.(SCOPF[:cost_Gen_C]),dims=(1,2)),digits=4)[C_rep.Con_N]
        # C_rep."cost_ESS"=sum(value.(SCOPF[:cost_ESS_C]),dims=(1,2))[C_rep.Con_N]
        C_rep."cost_FL"=round.(sum(value.(SCOPF[:cost_FL_C]),dims=(1,2);init=0.0),digits=4)[C_rep.Con_N]
        C_rep."PG(pu)"=reshape(round.(sum(abs.(value.(SCOPF[:PG_c_Line])),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep."QG(pu)"=reshape(round.(sum(abs.(value.(SCOPF[:QG_c_Line])),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep."Pinc(pu)"=reshape(round.(sum(value.(SCOPF[:P_inc_c_Line]),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep."Pdec(pu)"=reshape(round.(sum(value.(SCOPF[:P_dec_c_Line]),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep."Pch(pu)"=reshape(round.(sum(value.(SCOPF[:P_ch_c_Line]),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep."Pdis(pu)"=reshape(round.(sum(value.(SCOPF[:P_dis_c_Line]),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep."shed(pu)"=reshape(round.(sum(value.(SCOPF[:shed_c_Line]),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep."curt(pu)"=reshape(round.(sum(value.(SCOPF[:curt_c_Line]),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep."curtc(pu)"=reshape(round.(sum(value.(SCOPF[:curtc_c_Line]),dims=(1,2)),digits=4),(raw["nConL"],))[C_rep.Con_N]
        C_rep=sort(C_rep, :cost, rev=true)
    else C_rep=DataFrame() end

    T_rep=DataFrame(Cost=round.(value.(SCOPF[:cost_T]),digits=4),
                    Cost_Pre=round.(value.(SCOPF[:cost_Pre_T]),digits=4),
                    Cost_Post=round.(value.(SCOPF[:cost_Post_T]),digits=4),
                    )

    MPOutput= DataFrame(    OFValue=round.(JuMP.objective_value(SCOPF),digits=4),
                            OFValue_WP=round.(sum(value.(SCOPF[:cost_T])),digits=4),
                            Penalty=round.(sum(value.(SCOPF[:cost_Pen])),digits=4),
                            solve_time=round.(JuMP.solve_time(SCOPF),digits=4),
                            Pre_Cost=round.(sum(value.(SCOPF[:cost_Pre_T])),digits=4),
                            Post_Cost=round.(sum(value.(SCOPF[:cost_Post_T])),digits=4),
                            Pre_Cost_CS=round.(sum(value.(SCOPF[:cost_Shed]+SCOPF[:cost_Curt]+SCOPF[:cost_Curtc])),digits=4),
                            Pre_Cost_R=round.(sum(value.(SCOPF[:cost_Flex])),digits=4),
                            real_v=DataFrame(round.(value.(SCOPF[:e]),digits=4),Symbol.(1:raw["nBus"])),
                            imag_v=DataFrame(round.(value.(SCOPF[:f]),digits=4),Symbol.(1:raw["nBus"])),
                            mag_v=DataFrame(round.(sqrt.(abs.((value.(SCOPF[:e]).^2) + (value.(SCOPF[:f]).^2))),digits=4),Symbol.(1:raw["nBus"])),
                            ang_v=DataFrame(round.(atan.(value.(SCOPF[:f]),value.(SCOPF[:e])),digits=4),Symbol.(1:raw["nBus"])),
                            Line_flow=DataFrame(round.(sqrt.(abs.(value.(SCOPF[:Line_flow][:,:,1]))),digits=4),Symbol.(1:raw["nLine"])),
                            PG=DataFrame(round.(value.(SCOPF[:PG]),digits=4),Symbol.(1:raw["nGen"])),
                            QG=DataFrame(round.(value.(SCOPF[:QG]),digits=4),Symbol.(1:raw["nGen"])),
                            Pdis=DataFrame(round.(value.(SCOPF[:P_dis]),digits=4),Symbol.(1:raw["nESS"])),
                            Pch=DataFrame(round.(value.(SCOPF[:P_ch]),digits=4),Symbol.(1:raw["nESS"])),
                            soc=DataFrame(round.(value.(SCOPF[:soc]),digits=4),Symbol.(1:raw["nESS"])),
                            UE=DataFrame(round.(value.(SCOPF[:UE]),digits=4),Symbol.(1:raw["nFL"])),
                            Pinc=DataFrame(round.(value.(SCOPF[:P_inc]),digits=4),Symbol.(1:raw["nFL"])),
                            Pdec=DataFrame(round.(value.(SCOPF[:P_dec]),digits=4),Symbol.(1:raw["nFL"])),
                            Pcurt=DataFrame(round.(value.(SCOPF[:curt]),digits=4),Symbol.(1:raw["nRES"])),
                            Pcurtc=DataFrame(round.(value.(SCOPF[:curtc]),digits=4),Symbol.(1:raw["nGen"])),
                            Pshed=DataFrame(round.(value.(SCOPF[:shed]),digits=4),Symbol.(1:raw["nLoad"])),
                            Shunt=DataFrame(round.((raw["Shunt"]."bsh0")',digits=4),Symbol.(1:raw["nShunt"])),
                            C_report=C_rep,
                            T_report=T_rep,
    )
    Initials= DataFrame(
                            PG= MPOutput."PG"[1].+initraw["PG"],
                            QG= MPOutput."QG"[1].+initraw["QG"],
                            Pinc= MPOutput."Pinc"[1].+initraw["Pinc"],
                            Pdec= MPOutput."Pdec"[1].+initraw["Pdec"],
                            Pshed= MPOutput."Pshed"[1].+initraw["Pshed"],
                            Pcurt= MPOutput."Pcurt"[1].+initraw["Pcurt"],
                            Pcurtc= MPOutput."Pcurtc"[1].+initraw["Pcurtc"],
                            Pch= MPOutput."Pch"[1].+initraw["Pch"],
                            Pdis= MPOutput."Pdis"[1].+initraw["Pdis"],
                            soc= MPOutput."soc"[1],
                            UE= MPOutput."UE"[1],
                            Tcost= DataFrame(cost=value.(SCOPF[:cost_Pre_T])),
                            Vmag=MPOutput."mag_v"[1],
                            Vang=MPOutput."ang_v"[1]
                            )
    return MPOutput,Initials
end

function Get_ConResult(cnum,SCOPF,raw)

    Conout= DataFrame(
                            cost=DataFrame(cost=round.(value.(SCOPF[:cost_Shed_C]+SCOPF[:cost_Curt_C]+SCOPF[:cost_Curtc_C]+SCOPF[:cost_Flex_C]),digits=4)[:,cnum],
                            cost_shed=round.(value.(SCOPF[:cost_Shed_C]),digits=4)[:,cnum],
                            cost_curt=round.(value.(SCOPF[:cost_Curt_C]),digits=4)[:,cnum],
                            cost_curtc=round.(value.(SCOPF[:cost_Curtc_C]),digits=4)[:,cnum],
                            cost_flex=round.(value.(SCOPF[:cost_Flex_C]),digits=4)[:,cnum]                            ),
                            real_v=DataFrame(round.(value.(SCOPF[:e_c_Line])[:,:,cnum],digits=4),Symbol.(1:raw["nBus"])),
                            imag_v=DataFrame(round.(value.(SCOPF[:f_c_Line])[:,:,cnum],digits=4),Symbol.(1:raw["nBus"])),
                            mag_v=DataFrame(round.(sqrt.(abs.((value.(SCOPF[:e_c_Line])[:,:,cnum].^2) + (value.(SCOPF[:f_c_Line])[:,:,cnum].^2))),digits=4),Symbol.(1:raw["nBus"])),
                            ang_v=DataFrame(round.(atan.(value.(SCOPF[:f_c_Line])[:,:,cnum],value.(SCOPF[:e_c_Line])[:,:,cnum]),digits=4),Symbol.(1:raw["nBus"])),
                            Line_flow=DataFrame(round.(sqrt.(abs.(value.(SCOPF[:Line_flow_c_Line])))[:,:,cnum],digits=4),Symbol.(1:raw["nLine"])),
                            PG=DataFrame(round.(value.(SCOPF[:PG_c_Line].-SCOPF[:PG])[:,:,cnum],digits=4),Symbol.(1:raw["nGen"])),
                            QG=DataFrame(round.(value.(SCOPF[:QG_c_Line].-SCOPF[:QG])[:,:,cnum],digits=4),Symbol.(1:raw["nGen"])),
                            Pdis=DataFrame(round.(value.(SCOPF[:P_dis_c_Line].-SCOPF[:P_dis])[:,:,cnum],digits=4),Symbol.(1:raw["nESS"])),
                            Pch=DataFrame(round.(value.(SCOPF[:P_ch_c_Line].-SCOPF[:P_ch])[:,:,cnum],digits=4),Symbol.(1:raw["nESS"])),
                            soc=DataFrame(round.(value.(SCOPF[:soc_c_Line])[:,:,cnum],digits=4),Symbol.(1:raw["nESS"])),
                            Pinc=DataFrame(round.(value.(SCOPF[:P_inc_c_Line].-SCOPF[:P_inc])[:,:,cnum],digits=4),Symbol.(1:raw["nFL"])),
                            Pdec=DataFrame(round.(value.(SCOPF[:P_dec_c_Line].-SCOPF[:P_dec])[:,:,cnum],digits=4),Symbol.(1:raw["nFL"])),
                            Pcurt=DataFrame(round.(value.(SCOPF[:curt_c_Line].-SCOPF[:curt])[:,:,cnum],digits=4),Symbol.(1:raw["nRES"])),
                            Pcurtc=DataFrame(round.(value.(SCOPF[:curtc_c_Line].-SCOPF[:curtc])[:,:,cnum],digits=4),Symbol.(1:raw["nGen"])),
                            Pshed=DataFrame(round.(value.(SCOPF[:shed_c_Line].-SCOPF[:shed])[:,:,cnum],digits=4),Symbol.(1:raw["nLoad"])),
    )
    return Conout
end