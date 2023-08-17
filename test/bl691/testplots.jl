using JPerplex
using CairoMakie
fig = Figure(font = "Trebuchet",fontsize = 18)
ax = Axis(fig[1,1],
    xlabelsize = 28,
    ylabelsize = 28,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xticksmirrored = true,
    yticksmirrored = true,
    xtickalign = 1,
    ytickalign = 1,
    xminortickalign = 1,
    yminortickalign = 1,
    xticksize = 10,
    xtickwidth = 2,
    yticksize = 10,
    ytickwidth = 2,
    xminorticksize = 5,
    xminortickwidth = 2,
    yminorticksize = 5,
    yminortickwidth = 2,
    xgridvisible = false,
    ygridvisible = false
)

pseudo = getPseudosection("bl691")
x = getX.(pseudo.assemblages)
y = getY.(pseudo.assemblages)
xAx = pseudo.xAx
yAx = pseudo.yAx
if occursin("T (K)", pseudo.xAx) && tempInC
    x = x.-273.15
    xAx = "T (°C)"
end
if occursin("T (K)",pseudo.yAx) && tempInC
    y = y.-273.15
    yAx = "T (°C)"
end
key = 6
contourf!(x,y,filterGrid(pseudo,key),levels =0:0.5:1)
ax.xlabel = xAx
ax.ylabel = yAx

save("testfig"*string(key)*".svg",fig)