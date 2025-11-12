methods_=[1,3,5]
StartP=6
global Pathplot=Path*"../Result/Iter_Report/toplot/" 
code=[StartP,methods_]
scale_="on"

line1, line2, line5=0,0,0
# code=code_
datadict=Dict()
for c in code[2]
    try
        raw=read_iter_rep(Problem_Type,(code[1],c))
        fi=findfirst(x->x=="end",raw[:,1])
        datadict[c]=coalesce.(raw[2:fi-1,:],0)
    catch
        continue
    end
end

row_titles = ["Voltage Violation (p.u.)", "Line Flow Violation (MW)"]
methods = ["[17]", "LSB", "Hybrid LSB"]
subplot_labels = ["(a)" "(b)" "(c)"; "(d)" "(e)" "(f)"]  # one per subplot
inch = 96
pt = 4/3
cm = inch / 2.54
fig = Figure(resolution=(900, 600), fontsize = 13pt)
# fig = Figure(size = (9inch, 7inch), fontsize = 12pt)

# Prepare axis holders for each row
axes_row1 = Axis[]  # Voltage Violation axes
axes_row2 = Axis[]  # Current Violation axes
# axes_row3 = Axis[]  # Objective Gap axes
# axes_row4 = Axis[]  # ramp violation axes

# c=5;i=2
axes = Matrix{Axis}(undef, 2, 3)
for (i, c) in enumerate(code[2])
    if c in keys(datadict)
        data = datadict[c]
        time_method1 = data[:, "Res_Time"]
        max_V_violation_N = max.(data[:, "Vmax_Abs"], data[:, "Vmin_Abs"])
        max_V_violation_C = max.(data[:, "Vmax_C_Abs"], data[:, "Vmin_C_Abs"])
        Max_I_violation_N = 100* data[:, "Smax_Abs"]
        Max_I_violation_C = 100* data[:, "Smax_C_Abs"]
        Max_ramp_violation_N = data[:, "Ramp_Abs"]
        Max_ramp_violation_C = data[:, "Ramp_C_Abs"]
        opt_gap = data[:, "OptGap"]
    end

    ax = Axis(fig[1, i];
    xlabel = "",
    ylabel = "",
    title = methods[i],
    titlealign = :center,
    titlesize = 18,
    titlefont = :bold,
    # xlabelsize = 12,
    # xticklabelsize = 11,
    # yticklabelsize = 11,
    xgridvisible = true,
    ygridvisible = true,
    xgridstyle = :dash,
    ygridstyle = :dash,
    backgroundcolor = (:white, 0.0),
    subtitle = subplot_labels[1,i]
    # subtitlesize = 12
    )   
    axes[1, i] = ax

    ax = Axis(fig[2, i];
    xlabel = "Time (seconds)",
    ylabel = "",
    title = "",
    titlealign = :center,
    titlesize = 18,
    titlefont = :bold,
    # xlabelsize = 12,
    # xticklabelsize = 11,
    # yticklabelsize = 11,
    xgridvisible = true,
    ygridvisible = true,
    xgridstyle = :dash,
    ygridstyle = :dash,
    backgroundcolor = (:white, 0.0),
    subtitle = subplot_labels[2,i]
    )   
    # text!(ax,"(a)",position=(xmax,ymax), align = (:right, :top), fontsize=15)  
    axes[2, i] = ax
    
    if c in keys(datadict)
        for t in time_method1
            vlines!(ax, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
            vlines!(ax2, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
            # vlines!(ax3, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
            # vlines!(ax4, t, color = :gray, linestyle = :dash, linewidth = 1, alpha = 0.3)
        end

        line1=scatterlines!(axes[1, i], time_method1, max_V_violation_N, label = "Max Violation (Method 1)", color = :red, linewidth = 2, markersize = 8, marker = :diamond)
        line2=scatterlines!(axes[1, i], time_method1, max_V_violation_C, label = "Max Violation (Method 1)", color = :orange, linewidth = 2, markersize = 8, linestyle = :dash, marker = :diamond)

        scatterlines!(axes[2, i], time_method1, Max_I_violation_N, label = "Max Violation (Method 2)", color = :red, linewidth = 2, markersize = 8)
        scatterlines!(axes[2, i], time_method1, Max_I_violation_C, label = "Max Violation (Method 2)", color = :orange, linewidth = 2, linestyle = :dash, markersize = 8)

        # scatterlines!(ax3, time_method1, Max_ramp_violation_N, label = "Max Violation (Method 2)", color = :red, linewidth = 2, markersize = 8)
        # scatterlines!(ax3, time_method1, Max_ramp_violation_C, label = "Max Violation (Method 2)", color = :orange, linewidth = 2, linestyle = :dash, markersize = 8)

        # line5=scatterlines!(ax4, time_method1, opt_gap, label = "OptGap (Method 2)", color = :green, linewidth = 2, linestyle = :dash, markersize = 8)
    end
end
if scale_=="on"
    # linkxaxes!(axes_row1...)  # Voltage Violation axes
    # linkxaxes!(axes_row2...)  # Current Violation axes
    # linkxaxes!(axes_row3...)  # Opt Gap axes
    # linkxaxes!(axes_row4...)  # ramp violation axes
    linkyaxes!(axes[1,:])  # Voltage Violation axes
    linkyaxes!(axes[2,:])  # Current Violation axes
    # linkyaxes!(axes_row3...)  # Opt Gap axes
    # linkyaxes!(axes_row4...)  # ramp violation axes
end

Label(fig[1, 0], row_titles[1];
fontsize = 18,
font = "sans-bold",
rotation = pi/2,        # vertical
halign = :right,
valign = :center,
tellheight = false,
padding = (0cm, -0.4cm, 0cm, 0cm))

Label(fig[2, 0], row_titles[2];
fontsize = 18,
font = "sans-bold",
rotation = pi/2,        # vertical
halign = :right,
valign = :center,
tellheight = false,
padding = (0cm, -0.4cm, 0cm, 0cm))

Legend(fig[3, :],
[line1, line2],
["Normal Operation", "Contingency cases"],
patchsize = (20, 20), rowgap = 10,orientation = :horizontal)

colgap!(fig.layout, 20)
rowgap!(fig.layout, 20)
fig

