using JPerpleX
using CairoMakie
using DataFrames
# using Makie.Colors
# using Statistics


include("../../PlotDefaults.jl")
set_theme!(myTheme)
fig = Figure()
ax = Axis(fig[1,1])
ax.xlabel = pseudo.xaxis
ax.ylabel = pseudo.yaxis
xmin = minimum(x.(pseudo.assemblages))
ymin = minimum(y.(pseudo.assemblages))
xmax = maximum(x.(pseudo.assemblages))
ymax = maximum(y.(pseudo.assemblages))
ax.limits = (xmin,xmax,ymin,ymax)
pseudo = get_pseudosection("bl691",iscelsius=true, iskbar = true)

pseudosection!(ax,pseudo)

save("bl691.svg",fig)

output_assemblages("bl691_assemblages",pseudo)
# empty!(ax) #Empty axis but maintains extent, formatting, etc
contourdata = read_werami_output("bl691_1.tab",tempInC = true, pInKBar = true)

contour!(contourdata[!,1],contourdata[!,2],contourdata[!,3],labels=true)
save("bl691_melt_overlay.svg",fig)


