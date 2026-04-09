# In the name of Allah, the most beneficint, the most merciful



function csc_ii(bi, B, L, T, G, Ld, Shu, REN, BAT, y_bus, A, Gcas)
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    nBus = size(B, 1)
    nLine = size(L, 1)
    Bus_i = bi
    nBus_i = length(bi)
    vmag_i = ones(nBus_i)
    vang_i = zeros(nBus_i)
    vang_i = zeros(nBus_i)
    p_i = zeros(nBus_i)
    q_i = zeros(nBus_i)
    cg_i = zeros(1)
    pim_i = zeros(1)
    itr_i = zeros(1)
    r_i = []
    blkout = 0
    #----------------------------------------------------------------------------------------------------------
    #                               STEP 1: re-build island                 
    #              (what is the specific <B, L, T, G, Ld, Shu> for the supposed island)      
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    # (B_i, L_i, T_i, G_i, Ld_i, Shu_i, y_bus_i, A_i) = build_i(bi, B, L, T, G, Ld, Shu, y_bus, A)

    (B_i, L_i, T_i, G_i, Ld_i, Shu_i, BAT_i, REN_i, y_bus_i, A_i) = build_is(bi, B, L, T, G, Ld, Shu, BAT, REN, y_bus, A)


    Line_i = L_i[:, 9]
    Tr_i = T_i[:, 11]
    gi = G_i[:, 16]
    di = Ld_i[:, 14]
    reni= REN_i[:, 1]
    bati= BAT_i[:, 28]

    ctd = find_ones(A_i)

    #----------------------------------------------------------------------------------------------------------
    #                               STEP 2: perform UFLS and pre-UVLS                    
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------

    # (G_i, Ld_i, blkout) = pre_uvls(y_bus_i, B_i, L_i, G_i, Ld_i, Shu_i, ctd)
    (G_i, Ld_i, df_i, blkout, Gcas, BAT_i, REN_i) = floop3(G_i, Ld_i,Gcas, BAT_i, REN_i);




    if blkout == 0
        #----------------------------------------------------------------------------------------------------------
        #                               STEP 3: perform UVLS             
        #----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------
        (B_i, L_i, T_i, G_i, Ld_i, Shu_i, REN_i, BAT_i, vmag_i, vang_i, delta_mismatch_inj, mmm_i, p_i, q_i, p_V, q_V, cg_i, pim_i, df_i, fm_i, sm_i, Gs, Bs, itr_i, sij_i, vio_i, Iij, percentIij, Jcb, Gcas) = assessV4(B_i, L_i, T_i, G_i, Ld_i, Shu_i, REN_i, BAT_i, y_bus_i, ctd, Gcas)
      


        #----------------------------------------------------------------------------------------------------------
        #                               STEP 4: return island values            
        #----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------
        order=0
        if Gcas.print
            println(" after uvls,   the frequency will be changed :       $df_i pu")
        end
        r_i = Dict("Bus_i" => Bus_i, "Line_i" => Line_i, "Tr_i" => Tr_i, "vmag_i" => vmag_i, "vang_i" => vang_i, "p_i" => p_i, "q_i" => q_i,
            "cg_i" => cg_i, "pim_i" => pim_i, "itr_i" => itr_i, "sij_i" => sij_i, "vio_i" => vio_i, "df_i" => df_i, "Iij"=>Iij, "percentIij"=>percentIij, "cg_i"=>cg_i)

    end


    G[gi, 2] = G_i[:, 2]
    Ld[di, 2] = Ld_i[:, 2]
    Ld[di, 3] = Ld_i[:, 3]
    REN[reni, "PG"] = REN_i[:, "PG"]
    REN[reni, "QG"] = REN_i[:, "QG"]
    BAT[bati, "Ps"] = BAT_i[:, "Ps"]
    BAT[bati, "Qs"] = BAT_i[:, "Qs"]

    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------


    return r_i, B, L, T, G, Ld, Shu, REN, BAT, blkout, Jcb, Gcas
end


#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#--------------------------------------DO SPECIAL CASES LIKE THIS:-----------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------

function csc_i(bi, G, Ld, BAT, REN, Gcas)

    flag = 0
    bi_set = Set(bi)
    gi = findall(x -> x in bi_set, G[:, 1])  # which rows of G?
    di = findall(x -> x in bi_set, Ld[:, 1]) # which rows of Ld?
    bati = findall(x -> x in bi_set, BAT[:, 1]) # which rows of BAT?
    reni = findall(x -> x in bi_set, REN[:, 2]) # which rows of REN?
    G_i = G[gi, :]                         # one bus consisting of g&d
    Ld_i = Ld[di, :]
    Bat_i = BAT[bati, :]
    Ren_i= REN[reni, :]
    dmain=Ld[di, 14]
    gmain=G[gi, 16]

    if length(gi) == 0 || length(di) == 0
        BAT[bati, "Ps"] .= zeros(length(bati), 1)
        BAT[bati, "Qs"] .= zeros(length(bati), 1)
        REN[reni, "PG"] .= zeros(length(reni), 1)
        REN[reni, "QG"] .= zeros(length(reni), 1)
    end
        

    if length(gi) == 0 && length(di) == 0      # no g/ no d
        flag = 1
    elseif length(gi) == 0 && length(di) != 0  # no g
        Gcas.LdI[dmain] = 4 * ones(length(di), 1)
        Gcas.Ldsh_f[dmain] =Gcas.Ldsh_f[dmain]+ Ld[di, 2]
        Ld[di, 2] .= zeros(length(di), 1)
        # Ld[[3,4,5],2] = zeros(3, 1)
        Ld[di, 3] .= zeros(length(di), 1)
        flag = 2
    elseif length(gi) != 0 && length(di) == 0  # no d
        Gcas.Gt[gmain] = Gcas.Gt[gmain] + G[gi, 2] 
        G[gi, 2] .= zeros(length(gi), 1)
        flag = 3
    elseif sum(G_i[:, 2])+sum(Ren_i[:, "PG"])+sum(filter(x -> x > 0, Bat_i[:, "Ps"]);init=0) == 0 && sum(Ld_i[:, 2])-sum(filter(x -> x < 0, Bat_i[:, "Ps"]);init=0) == 0   # every thing off
        flag = 4
    
    elseif fr2(G_i, Ld_i, BAT, REN) > df_maxmax || fr2(G_i, Ld_i, BAT, REN) < -df_minmin
        Gcas.LdI[dmain] = 4 * ones(length(di), 1)
        Gcas.Ldsh_f[dmain] =Gcas.Ldsh_f[dmain]+ Ld[di, 2]
        Gcas.Gt[gmain] = Gcas.Gt[gmain] + G[gi, 2] 
        Ld[di, 2] .= zeros(length(di), 1)
        Ld[di, 3] .= zeros(length(di), 1)
        G[gi, 2] .= zeros(length(gi), 1)
        BAT[bati, "Ps"] .= zeros(length(bati), 1)
        BAT[bati, "Qs"] .= zeros(length(bati), 1)
        REN[reni, "PG"] .= zeros(length(reni), 1)
        REN[reni, "QG"] .= zeros(length(reni), 1)
        flag = 5
    
    elseif length(bi) == 1 && G[gi, 1] == Ld[di, 1]
        G_i = G[gi, :]                         # one bus consisting of g&d
        Ld_i = Ld[di, :]
        (G_i, Ld_i, df_i, blkout, Gcas) = floop3(G, Ld, Gcas, BAT, REN)
        G[gi, 2] = G_i[:, 2]
        Ld[di, 2] = Ld_i[:, 2]
        Ld[di, 3] = Ld_i[:, 3]
        flag = 6
    end

    return G, Ld, BAT, REN, flag, Gcas
end

#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------


