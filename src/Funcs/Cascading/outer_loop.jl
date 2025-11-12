# In the name of Allah, the most beneficint, the most merciful
function outer(B, L, T, G, Ld, Shu, Lcnd, Tcnd, Gcnd, Ldcnd)

    n_ini = length(Lcnd) + length(Tcnd) + length(Gcnd) + length(Ldcnd)

    nL = length(Lcnd)
    nT = length(Tcnd)
    nG = length(Gcnd)
    nLd = length(Ldcnd)

    r_result = []

    Ltemp = deepcopy(L)
    Ttemp = deepcopy(T)
    Gtemp = deepcopy(G)
    Ldtemp = deepcopy(Ld)

    io = 1

    n_ini=1
    while io <= n_ini
        #-------------------------------------------------------------------------------------
        # [-][-][-][-][-][-][-][-][-][-][-][-]  Select an initial outage and run for it [-][-]
        #-------------------------------------------------------------------------------------
        # if 1 <= io && io <= nL
        #     # Lres = setdiff(collect(1:size(L, 1)), Lcnd[io])
        #     # Ltemp = L[Lres, :]

        #     rt = CSC(B, Ltemp, Ttemp, Gtemp, Ldtemp, Shu)

        #     push!(r_result, rt)
        #     Ltemp = deepcopy(L)

        # elseif (nL + 1) <= io && io <= (nL + nT)
        #     Tres = setdiff(collect(1:size(T, 1)), Tcnd[io-nL])
        #     Ttemp = T[Tres, :]

        #     rt = CSC(B, Ltemp, Ttemp, Gtemp, Ldtemp, Shu)
        #     push!(r_result, rt)
        #     Ttemp = deepcopy(T)

        # elseif (nL + nT + 1) <= io && io <= (nL + nT + nG)

        #     # !!! dont forget the copuling between some gens and trsns!!!!!
        #     # !!! dont forget the copuling between some gens and trsns!!!!!
        #     # !!! dont forget the copuling between some gens and trsns!!!!!
        #     # !!! dont forget the copuling between some gens and trsns!!!!!
        #     # !!! dont forget the copuling between some gens and trsns!!!!!

        #     Gres = setdiff(collect(1:size(G, 1)), Gcnd[io-nL-nT])
        #     Gtemp = G[Gres, :]

        #     rt = CSC(B, Ltemp, Ttemp, Gtemp, Ldtemp, Shu)
        #     push!(r_result, rt)
        #     Gtemp = deepcopy(G)

        # elseif (nL + nT + nG + 1) <= io && io <= n_ini
        #     Ldres = setdiff(collect(1:size(Ld, 1)), Ldcnd[io-nL-nT-nG])
        #     Ldtemp = Ld[Ldres, :]

        #     rt = CSC(B, Ltemp, Ttemp, Gtemp, Ldtemp, Shu)
        #     push!(r_result, rt)
        #     Ldtemp = deepcopy(Ld)

        # end

        #-------------------------------------------------------------------------------------
        # [-][-][-][-][-][-][-][-][-][-][-][-] [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
        #-------------------------------------------------------------------------------------






        (rn, B, Lo, To, Go, Ldo, Shu) = CSC(B, Ltemp, Ttemp, Gtemp, Ldtemp, Shu)
        push!(r_result, rn)

        io += 1
        
        print("io is: ")
        println(io) 
    end




    return r_result, io

end



# Lcnd = [46, 52]
# Lcnd = collect(1:size(L, 1))
# Tcnd = [24, 25]
# Gcnd = []
# Ldcnd = []
# (r_result, io) = outer(B, L, T, G, Ld, Shu, Lcnd, Tcnd, Gcnd, Ldcnd)
# println(" Success!")