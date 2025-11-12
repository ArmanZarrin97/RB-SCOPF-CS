# In the name of Allah, the most beneficint, the most merciful

function CSC(B::DataFrame, L::DataFrame, T::DataFrame, G::DataFrame, Ld::DataFrame, Shu::DataFrame, REN::DataFrame, BAT::DataFrame, Lout::Vector{Int64},Gcas)
    # B, L, T, G, Ld, Shu, Lout=Bs, Ls, Ts, Gs, Lds, Shus, Lout
    casrep=DataFrame(Stage=[], LL=[], LNS_V=[], LNS_F=[], LNS_T=[], Jcbs=[])
    Rn = []
    Lvio = []
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
            Lvio = main_name_L[r[i]["vio_i"]]
        end
        push!(Rn, r)
        push!(JCB, Jcb)
        push!(Line_vio, Lvio)


        LNS_f = sum(Gcas.Ldsh_f) 
        LNS_v = sum(Gcas.Ldsh_v) 
        LNS = LNS_f + LNS_v
        Lvio_tot = union(Lvio_tot, Lvio, Lout)
        Lvio_tot = sort(Lvio_tot)
        notLres=findall( in(Set(Lvio_tot)),Li."Line_N")

        push!(LNS_stage, LNS)
        push!(casrep,[H-1 [notLres] LNS_v LNS_f LNS [Jcb]])

        line_viol = length(Lvio)

        if line_viol != 0
            if Gcas.print 
                printstyled("NEW line violation is observed!"; color=:1, blink=true)
            end
            # Lvio_tot = union(Lvio_tot, Lvio, Lout)
            # Lvio_tot = sort(Lvio_tot)
            # # Lres = setdiff(collect(1:size(Li, 1)), Lvio_tot)
            # # L = Li[Lres, :]
            # notLres=findall( in(Set(Lvio_tot)),Li."Line_N")
            L = Li[Not(notLres), :]
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

    return Rn, H, u, B, L, T, G, Ld, Shu, new_L, LNS_v, LNS_f, LNS, tot_L, Line_vio, LNS_stage, JCB, casrep
end
