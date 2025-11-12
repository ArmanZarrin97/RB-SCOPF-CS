# In the name of Allah, the most beneficint, the most merciful

function total_i(B, L, T, G, Ld, Shu, REN, BAT, Gcas)

    #><><><><><><><><><><><><><><><><><><><><><><><><>

    # First: detect islands <><><><><><><><><><><><><>

    # ><><><><><><><><><><><><><><><><><><><><><><><><

    nBus = size(B, 1)
    n_i = 0
    Jcb = []

    A = spzeros(Int64, nBus, nBus)
    A = A_net(nBus, L, T)
    A.nzval .= 1  # Set all nonzero elements to 1 efficiently

    ctd = find_ones(A)
    y_bus = Y_bus2(nBus, L, T, ctd)

    buses =  Set(B[:, 1])
    lines =  Tuple{Int, Int}.(eachrow(L[:, 1:2])) 
    trans =  Tuple{Int, Int}.(eachrow(T[:, 1:2])) 
    lines =  vcat(lines, trans)
    
    u = find_islands(buses, lines)


    #><><><><><><><><><><><><><><><><><><><><><><><><>

    # Second: for various islands, for each (csc_i) is
    # is performed as follows: ><><><><><><><><><><><>
    # - re-build its corresponding network, ><><><><><
    # - solve it, ><><><><><><><><><><><><><><><><><><
    # - feed back the results ><><><><><><><><><><><><

    # ><><><><><><><><><><><><><><><><><><><><><><><><

    r = []
    blkout = 0

    # __r = Dict(
    #     "S0" => [],
    #     "vio_i" => [],
    #     "Jcb" => [],
    #     "LNS_f" => [],
    #     "LNS_v" => [],
    #     "buses" => [],
    # )
    for i in 1:length(u)

        bi = u[i]
        bi = sort(bi)
        # push!(__r["buses"], bi)
        (G, Ld, flag, Gcas) = csc_i(bi, G, Ld, Gcas);
        if Gcas.print
            printstyled("- island#$i, flag=$flag, {$bi}:"; color=:blue, blink=true)
            println()
            printstyled("  $(length(bi))-bus"; color=:blue, blink=true)
            println()
            if flag==4 || flag==5
                println("blackout of island, due to extreme frequency issue")
                # push!(__r["Jcb"], NaN)
            end
        end
        if flag == 0
           (r_i, B, L, T, G, Ld, Shu, blkout, Jcb, Gcas) = csc_ii(bi, B, L, T, G, Ld, Shu, REN, BAT, y_bus, A, Gcas);
            # push!(__r["Jcb"], Jcb);push!(__r["S0"], r_i["sij_i"]);push!(__r["vio_i"], r_i["vio_i"])
            # push!(__r["LNS_f"], r_i["LNS_f"]);push!(__r["LNS_v"], r_i["LNS_v"])
            if blkout == 0
                n_i += 1
                push!(r, r_i)
            else
                if Gcas.print
                    printstyled("- island#$i, flag=1, {$bi}:"; color=:blue, blink=true)
                    printstyled("- island#$i, flag=1, {$bi}:"; color=:blue, blink=true)
                    println()
                end
            end
        end
    end

    # ><><><><><><><><><><><><><><><><><><><><><><><><
    # ><><><><><><><><><><><><><><><><><><><><><><><><


    return r, n_i, u, B, L, T, G, Ld, Shu, Jcb

end

