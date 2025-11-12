# In the name of Allah, the most beneficint, the most merciful

function total_i(B, L, T, G, Ld, Shu)

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
    y_bus = Y_bus1(nBus, L, T, ctd)

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

    for i in 1:length(u)

        bi = u[i]
        bi = sort(bi)
        (G, Ld, flag) = csc_i(bi, G, Ld)
        printstyled("- island#$i, flag=$flag, {$bi}:"; color=:blue, blink=true)
        println()
        printstyled("  $(length(bi))-bus"; color=:blue, blink=true)
        println()
        if flag==4 || flag==5
            println("blackout of island, due to extreme frequency issue")
        end
        if flag == 0
           (r_i, B, L, T, G, Ld, Shu, blkout, Jcb) = csc_ii(bi, B, L, T, G, Ld, Shu, y_bus, A)
            if blkout == 0
                n_i += 1
                push!(r, r_i)
            else
                printstyled("- island#$i, flag=1, {$bi}:"; color=:blue, blink=true)
                printstyled("- island#$i, flag=1, {$bi}:"; color=:blue, blink=true)
                println()
            end
        end
    end

    # ><><><><><><><><><><><><><><><><><><><><><><><><
    # ><><><><><><><><><><><><><><><><><><><><><><><><


    return r, n_i, u, B, L, T, G, Ld, Shu, Jcb

end

