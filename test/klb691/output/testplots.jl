using JPerpleX
using CairoMakie
using DataFrames
# using Makie.Colors
# using Statistics

fig = Figure(font = "B612",fontsize = 18)
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
    ygridvisible = false,
    aspect = 1.0
)

pseudo = getPseudosection("klb691",tempInC=true, pInKBar = true)
# x = getX.(pseudo.assemblages)
# y = getY.(pseudo.assemblages)
# xAx = pseudo.xAx
# yAx = pseudo.yAx



# function minMaxPhaseVar(pseudo::PerplexGrid)
#     uniqueAsms = listUniqueAssemblages(pseudo.assemblages)
#     maxVar = 0
#     minVar = 100
#     for asm in uniqueAsms
#         if length(asm.phases) >maxVar
#             maxVar = length(asm.phases)
#         end
#         if length(asm.phases) <minVar
#             minVar = length(asm.phases)
#         end

#     end

#     return minVar, maxVar
# end
# uniqueAsms = listUniqueAssemblages(pseudo.assemblages)
# for i in range(1,lastindex(uniqueAsms))
    
#     minV, maxV = minMaxPhaseVar(pseudo)
#     colorVal = 1-(1-(length(uniqueAsms[i].phases)-minV)/(maxV-minV))*0.8
#     contourf!(ax,x,y,filterKeyArray(pseudo.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,Colors.HSV(0,0,colorVal)])
    
# end

# for i in range(1,lastindex(uniqueAsms))
#     contour!(ax,x,y,filterKeyArray(pseudo.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
#     iGrid = filterGrid(pseudo,i)
#     scatter!(ax,mean(getX.(iGrid)),mean(getY.(iGrid)), marker = :circle, strokecolor = :black,strokewidth = 1,color = :transparent)
#     text!(ax,mean(getX.(iGrid)),mean(getY.(iGrid)),text = string(i))
# end

# # contourf!(x,y,filterGrid(pseudo,4),levels =-0.5:1:1.5,colormap = [:transparent,:green])
# # contourf!(x,y,filterGrid(pseudo,5),levels =-0.5:1:1.5,colormap = [:transparent,:orange])
# # contourf!(x,y,filterGrid(pseudo,6),levels =-0.5:1:1.5,colormap = [:transparent,:red])
# # contourf!(x,y,filterGrid(pseudo,7),levels =-0.5:1:1.5,colormap = [:transparent,:blue])
# # contourf!(x,y,filterGrid(pseudo,8),levels =-0.5:1:1.5,colormap = [:transparent,:purple])
# # contourf!(x,y,filterGrid(pseudo,9),levels =-0.5:1:1.5,colormap = [:transparent,Colors.HSV(0,0,1)])
# # contour!(x,y,filterGrid(pseudo,4),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
# # contour!(x,y,filterGrid(pseudo,5),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
# # contour!(x,y,filterGrid(pseudo,6),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
# # contour!(x,y,filterGrid(pseudo,7),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
# # contour!(x,y,filterGrid(pseudo,8),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
# # contour!(x,y,filterGrid(pseudo,9),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
# ax.xlabel = xAx
# ax.ylabel = yAx
# xMin = minimum(getX.(pseudo.assemblages))
# yMin = minimum(getY.(pseudo.assemblages))
# xMax = maximum(getX.(pseudo.assemblages))
# yMax = maximum(getY.(pseudo.assemblages))
# ax.limits = (xMin,xMax,yMin,yMax)

plotPseudosection!(ax,pseudo)

save("klb691.svg",fig)

outputAssemblages("klb691_assemblages",pseudo)
# empty!(ax) #Empty axis but maintains extent, formatting, etc


