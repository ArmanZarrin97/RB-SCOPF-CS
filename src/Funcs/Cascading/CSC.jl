# In the name of Allah, the most beneficint, the most merciful

function CSC(B::DataFrame, L::DataFrame, T::DataFrame, G::DataFrame, Ld::DataFrame, Shu::DataFrame, REN::DataFrame, BAT::DataFrame, Lout::Vector{Int64},Gcas)
    # B, L, T, G, Ld, Shu, REN, BAT, Lout=deepcopy(Bs), deepcopy(Ls), deepcopy(Ts), deepcopy(Gs), deepcopy(Lds), deepcopy(Shus), deepcopy(RENs), deepcopy(BATs), Int64[]
    casrep=DataFrame(Stage=[], LL=[],VM=[], prob=[], LNS_V=[], LNS_F=[], LNS_T=[], GT=[], RC=[], BT=[], Sev=[], Jcbs=[],busset=[],LNSset=[], Vmag=[], Vang=[])
    comcas_i=DataFrame(Stage=[], LL=[],VM=[], Jcbs=[],busset=[], Vmag=[], Vang=[])
    Rn = []
    Lvio = []
    Lvior=[]
    convergence=0
    Lvio_stage = []
    Lvior_stage = []
    Lvio_tot = []
    H = 1
    u = []
    Li = deepcopy(L)
    JCB = []
    Line_vio=[]
    Lres = setdiff(collect(1:size(L, 1)), Lout)
    L = L[Lres, :]
    LNS_stage =[]

    while true
        # for i in 1
        if Gcas.print
            println()
            printstyled("
    ------------------------|O|.|O|.|O|  CASCADE-$(H))  |O|.|O|.|O|-----------------------------------------"; color=:20, blink=true)
            println()
        end


        (r, n_i, u, B, L, T, G, Ld, Shu, Jcb) = total_i(B, L, T, G, Ld, Shu, REN, BAT,Gcas);

        for i = 1:n_i
            main_name_L = r[i]["Line_i"]
            Lvio_stage = r[i]["vio_i"]
            Lvior_stage = r[i]["vior_i"]
            Lvio = main_name_L[r[i]["vio_i"]]
            Lvior= main_name_L[r[i]["vior_i"]]
            convergence= r[i]["cg_i"]
        end
        push!(Rn, r)
        push!(JCB, Jcb)
        push!(Line_vio, Lvio)

        LNSset=Gcas.Ldsh_f.+Gcas.Ldsh_v
        LNS_f = sum(Gcas.Ldsh_f)
        LNS_v = sum(Gcas.Ldsh_v)
        LNS = LNS_f + LNS_v
        GT= sum(Gcas.Gt,init=0.0)
        BT= sum(Gcas.BT_L,init=0.0)+sum(Gcas.BT_G,init=0.0)
        RC= sum(Gcas.RC,init=0.0)
        Lvio_tot = union(Lvio_tot, Lvio, Lout)
        Lvio_tot = sort(Lvio_tot)
        notLres_tot=findall( in(Set(Lvio_tot)),Li."Line_N")
        # notLres=findall( in(Set(Lvio)),L."Line_N")

        push!(LNS_stage, LNS)
        # push!(casrep,[H-1 [Lvio] [if n_i==0; [] else r[n_i]["percentIij"][Lvio_stage] end] LNS_v LNS_f LNS [Jcb] [if n_i==0; [] else r[n_i]["Bus_i"] end] [LNSset] [if n_i==0; [] else r[n_i]["vmag_i"] end] [if n_i==0; [] else r[n_i]["vang_i"] end]])
        if n_i==0 || convergence==0
            LNS= sum(Ld[:,2])
            GT= sum(G[:,2])
            RC= sum(REN[:,2])
            BT= sum(BAT[:,2])
            Sev=SEVF(LNS, GT, RC, BT)
            push!(casrep,[H-1 [Lvio] [[]] 1 LNS_v LNS_f LNS GT RC BT Sev [Jcb] [[]] [LNSset] [[]] [[]]])
            # push!(comcas_i,[H-1 [Lvior] [[]][Jcb] [[]] [[]] [[]]])
            if Gcas.print
                printstyled("No Powerflow solution"; color=:10, blink=true)
            end
            break
        else
            Sev=SEVF(LNS, GT, RC, BT)
            push!(casrep,[H-1 [Lvio] [r[n_i]["percentIij"][Lvio_stage]] 1 LNS_v LNS_f LNS GT RC BT Sev [Jcb] [r[n_i]["Bus_i"]] [LNSset] [r[n_i]["vmag_i"]] [r[n_i]["vang_i"]]])
            # push!(casrep,[H-1 [Lvio] [r[n_i]["percentIij"][Lvio_stage]] LNS_v LNS_f LNS [Jcb] [r[n_i]["Bus_i"]] [LNSset] [r[n_i]["vmag_i"]] [r[n_i]["vang_i"]]])
            push!(comcas_i,[H-1 [Lvior] [r[n_i]["percentIij"][Lvior_stage]] [Jcb] [r[n_i]["Bus_i"]] [r[n_i]["vmag_i"]] [r[n_i]["vang_i"]]])

        end
        # push!(casrep,[H-1 [Lvio] [r[n_i]["percentIij"][Lvio_stage]] LNS_v LNS_f LNS [Jcb] [r[n_i]["Bus_i"]] [LNSset] [r[n_i]["vmag_i"]] [r[n_i]["vang_i"]]])

        
        line_viol = length(Lvio)

        if line_viol != 0 && n_i !=0
            if Gcas.print 
                printstyled("NEW line violation is observed!"; color=:1, blink=true)
            end
            # Lvio_tot = union(Lvio_tot, Lvio, Lout)
            # Lvio_tot = sort(Lvio_tot)
            # # Lres = setdiff(collect(1:size(Li, 1)), Lvio_tot)
            # # L = Li[Lres, :]
            # notLres=findall( in(Set(Lvio_tot)),Li."Line_N")
            L = Li[Not(notLres_tot), :]
            if Gcas.print
                println()
                printstyled("New   lines outaged: $Lvio"; color=:1, blink=true)
                println()
                printstyled("Total lines outaged: $Lvio_tot"; color=:1, blink=true)
            end
        elseif line_viol == 0 || size(L, 1) == 0
            if Gcas.print
                printstyled("NO line violation is observed"; color=:10, blink=true)
            end
            break
        else 
            if Gcas.print
                printstyled("No Powerflow solution"; color=:10, blink=true)
            end
            break
        end
 
        H += 1
    end

    new_L = length(Lvio)
    tot_L = length(Lvio_tot)

    if Gcas.print
        printstyled("
    ><><>< ----------------------------------------------------------------------------------------><><><"; color=:16, blink=true)
        println()
        printstyled("Simulation is done for a single-contingency"; color=:12, blink=true)
        println()
    end


    Ldsh = Gcas.Ldsh_f + Gcas.Ldsh_v
    LNS_f = sum(Gcas.Ldsh_f) * sb
    LNS_v = sum(Gcas.Ldsh_v) * sb
    LNS = LNS_f + LNS_v


    Gcas.Ldsh_v[1:size(Ld, 1)] = zeros(size(Ld, 1))
    Gcas.Ldsh_f[1:size(Ld, 1)] = zeros(size(Ld, 1))
    # Gcas.Gt[1:size(G, 1)] = zeros(size(G, 1))
    # Gcas.RC[1:size(REN, 1)] = zeros(size(REN, 1))
    # Gcas.BT[1:size(BAT, 1)] = zeros(size(BAT, 1))

    return Rn, H, u, B, L, T, G, Ld, Shu, new_L, LNS_v, LNS_f, LNS, tot_L, Line_vio, LNS_stage, JCB, casrep, comcas_i
end
