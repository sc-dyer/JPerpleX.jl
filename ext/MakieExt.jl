module MakieExt

using
    JPerpleX,
    Makie,
    Makie.Colors,
    ColorSchemes,
    Statistics,
    DocStringExtensions
    
import JPerpleX: pseudosection, pseudosection!

"""
$(TYPEDSIGNATURES)

Plot recipe to make a pseudosection from a 'PerplexGrid'. These will be plotted such that they 
can still be edited with a vector graphics editor. Fields are labeled by assemblage key, a list of unique 
assemblages can be retrieved with 'unique_assemblages' or they can be printed to a text file with 'output_assemblages'. 
"""
@recipe(Pseudosection) do scene
    Attributes(
        colormap = [Colors.HSV(0,0,0.2),Colors.HSV(0,0,1)],
        linewidth = 2,
        linecolor = :black,
        marker = :circle,
        markersize = 5,
        markerfill = :transparent,
        markerstroke = :black,
        markerstrokewidth = 1,
        labelfontsize = 12

    )
end


function Makie.plot!(pseudo::Pseudosection)
    
  

    pgrid = pseudo[1][]
    xs = x.(pgrid.assemblages)
    ys = y.(pgrid.assemblages)

    #Here we figure out the assemblage variance used in determining the field colour
    present_assemblages = unique_assemblages(pgrid.assemblages)
    max_phasevariance = 0
    min_phasevariance = 100
    for asm in present_assemblages
        if length(asm.phases) >max_phasevariance
            max_phasevariance = length(asm.phases)
        end
        if length(asm.phases) <min_phasevariance
            min_phasevariance = length(asm.phases)
        end
    end
    
    # colorramp= nothing
    # try
    #     colorramp = ColorSchemes.colorschemes[pseudo.colormap]
    # catch e
    #     colorramp = ColorScheme(pseudo.colormap)
    # end
    colorramp =  ColorScheme(pseudo.colormap[])
    #This part plots a filled contour around each unique assemblage
    #Ranges of the levels and colormaps are based on trial and error, do not change them
    #They should be roughly halfway beetween grid points in a smoothed line
    for i in range(1,lastindex(present_assemblages))
        colorvalue = 1-(length(present_assemblages[i].phases)-min_phasevariance)/(max_phasevariance-min_phasevariance)
        contourf!(pseudo,xs,ys,filterkeys(pgrid.assemblages,i),levels =-0.5:1:1.5,
        colormap = [:transparent,get(colorramp,colorvalue)])
    end

    #For added flair, we add outlines and labels to each assemblages
    for i in range(1,lastindex(present_assemblages))
        igrid = filtergrid(pgrid,i)
        scatter!(pseudo,mean(x.(igrid)),mean(y.(igrid)), marker = pseudo.marker, strokecolor = pseudo.markerstroke,strokewidth=pseudo.markerstrokewidth,
                    color = pseudo.markerfill,markersize = pseudo.markersize)
        text!(pseudo,mean(x.(igrid)),mean(y.(igrid)),text = string(i),fontsize = pseudo.labelfontsize)
    end
    #Making the contours seperate so they can be selected easily in post-processing
    for i in range(1,lastindex(present_assemblages))
        # println(string(uniqueAsms[i].key)*" = "*join(uniqueAsms[i].phases," "))
        contour!(pseudo,xs,ys,filterkeys(pgrid.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,pseudo.linecolor[],pseudo.linecolor[]],linewidth=pseudo.linewidth)
    end
    # pseudo.xlabel = pgrid.xaxis
    # pseudo.ylabel = pgrid.yaxis
    # xmin = minimum(x.(pgrid.assemblages))
    # ymin = minimum(y.(pgrid.assemblages))
    # xmax = maximum(x.(pgrid.assemblages))
    # ymax = maximum(y.(pgrid.assemblages))
    # pseudo.limits = (xmin,xmax,ymin,ymax)
end

end