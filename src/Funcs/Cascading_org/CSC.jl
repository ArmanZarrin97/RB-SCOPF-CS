# In the name of Allah, the most beneficint, the most merciful

function CSC(B, L, T, G, Ld, Shu, Lout)
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
        println()
        printstyled("
------------------------|O|.|O|.|O|  CASCADE-$(H))  |O|.|O|.|O|-----------------------------------------"; color=:20, blink=true)
        println()


        (r, n_i, u, B, L, T, G, Ld, Shu, Jcb) = total_i(B, L, T, G, Ld, Shu)

        for i = 1:n_i
            main_name_L = r[i]["Line_i"]
            Lvio = main_name_L[r[i]["vio_i"]]
        end
        push!(Rn, r)
        push!(JCB, Jcb)
        push!(Line_vio, Lvio)


        LNS_f = sum(Ldsh_f) 
        LNS_v = sum(Ldsh_v) 
        LNS = LNS_f + LNS_v
        push!(LNS_stage, LNS)

        line_viol = length(Lvio)

        if line_viol != 0
            printstyled("NEW line violation is observed!"; color=:1, blink=true)
            Lvio_tot = union(Lvio_tot, Lvio, Lout)
            Lvio_tot = sort(Lvio_tot)
            Lres = setdiff(collect(1:size(Li, 1)), Lvio_tot)
            L = Li[Lres, :]
            println()
            printstyled("New   lines outaged: $Lvio"; color=:1, blink=true)
            println()
            printstyled("Total lines outaged: $Lvio_tot"; color=:1, blink=true)
        elseif line_viol == 0 || size(L, 1) == 0
            printstyled("NO line violation is observed"; color=:10, blink=true)
            break
        end

        H += 1
    end

    new_L = length(Lvio)
    tot_L = length(Lvio_tot)

    printstyled("
><><>< ----------------------------------------------------------------------------------------><><><"; color=:16, blink=true)
    println()
    printstyled("Simulation is done for a single-contingency"; color=:12, blink=true)
    println()


    Ldsh = Ldsh_f + Ldsh_v
    LNS_f = sum(Ldsh_f) * sb
    LNS_v = sum(Ldsh_v) * sb
    LNS = LNS_f + LNS_v


    Ldsh_v[1:size(Lds, 1)] = zeros(size(Lds, 1))
    Ldsh_f[1:size(Lds, 1)] = zeros(size(Lds, 1))

    return Rn, H, u, B, L, T, G, Ld, Shu, new_L, LNS_v, LNS_f, LNS, tot_L, Line_vio, LNS_stage, JCB
end
