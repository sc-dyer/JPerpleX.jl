# myColours = ["#FF8CA6","#6D8FFF","#F9C52A","#A9EFA5", "#E6A7FB","#000000","#FF0000","#008000","#0000FF","#4B0082","#FF8C00"]
# myTheme = Theme(
#     fontsize = 20,
#     palette = (;color =myColours),
#     Axis = ( 
#         xlabelsize = 28,
#         ylabelsize = 28,
#         xminorticksvisible = true,
#         yminorticksvisible = true,
#         xticksmirrored = true,
#         yticksmirrored = true,
#         xtickalign = 1,
#         ytickalign = 1,
#         xminortickalign = 1,
#         yminortickalign = 1,
#         xticksize = 10,
#         xtickwidth = 2,
#         yticksize = 10,
#         ytickwidth = 2,
#         xminorticksize = 5,
#         xminortickwidth = 2,
#         yminorticksize = 5,
#         yminortickwidth = 2,
#         xgridvisible = false,
#         ygridvisible = false
#         )
# )
myColours = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]
myTheme = Theme(
    fontsize = 20,
    palette = (;color =myColours,marker = [:circle,:rect,:diamond,:cross, :xcross,:utriangle,:star5,:star8]),
    Poly = (;cycle = [:color]);
    Scatter = (;cycle = [:color,:marker],markersize = 15),
    Axis = ( 
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
)
