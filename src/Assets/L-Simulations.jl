function MPC_OPF_Simulation(LinMethod="None",Lininit_0=Dict{String, Array{Union{Matrix{Float64},Float64}}}())
    Initials=0
    for t=1:nTL*nTS
        global nTP=min(nTS*nP,nTL*nTS-t+1)
        t0=t
        print(t0)
        if (t0==1) initial_point=0 else initial_point=1  end


        # cr_con,cr_con1,Cost_Post_t1=Contingency_Screen!(initial_point,Initials,t0)
        # OPFres[t,"Before"]=collect(cr_con1)

        OPFout,Initials,timesolve,OPFmodel=L_RT_OPF!(initial_point,Initials,t0,LinMethod,Lininit_0)


        RT_result."case1"[1]+=OPFout."T_report"[1][1,1]
        # RT_result."State"[1]."case1"[1]=Initials
        OPFres[t,"cost"]=OPFout."T_report"[1][1,1]
        OPFres[t,"State"]=Initials
        OPFres[t,"time"]=timesolve
        initial_point=1;nTP=1;
        cr_con,cr_con1,Cost_Post_t1,Infes_Conv_rep=Contingency_Screen!(initial_point,Initials,t0)
        OPFres[t,"EConcost"]=Cost_Post_t1*0.05
        OPFres[t,"After"]=collect(cr_con1)
        OPFres[t,"InfRep"]= Infes_Conv_rep
    end
end