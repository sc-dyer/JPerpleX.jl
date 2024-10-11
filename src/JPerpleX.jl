
"""
This is the main module for a more "advanced" Perple_X julia wrapper. The advantages of this 
library over a more "basic" wrapper is that it calls the compiled fortran functions in perplexwrap.f 
directly and *shouldn't* require repeated file i/o or command piping.

# Exports
$(EXPORTS)
"""
module JPerpleX
export 
    init_meemum,
    minimizepoint,
    Assemblage,
    PerplexGrid,
    get_pseudosection,
    unique_assemblages,
    x,
    y,
    key,
    filtergrid,
    filterkeys,
    pseudosection,
    pseudosection!,
    output_assemblages,
    read_werami_output
using
    DocStringExtensions,
    Reexport,
    Statistics,
    DataFrames

@reexport using PetroBase
# Write your package code here.

#Constants from PerpleX, these MUST match the equivalent values in perplex_parameters.h
const K0 = 25
const K1 = 3000000
const K2 = 100000
const K3 = 2000
const K5 = 14
const K8 = K5+2
const I8 = 28
const L7 = 2048
const H9 = 30
const MAX_COMPONENT_NAME_L = 5
const MAX_PHASE_NAME_L = 14
const PURE_PHASE_NAME_L = 8
const SOL_PHASE_ABBREV_L = 6
const VARIABLE_NAME_L = 8

is_init::Bool = false
# const almostgreys = loadcolorscheme(:almostgreys,[Colors.HSV(0,0,0.2),Colors.HSV(0,0,1)])
"""
$(TYPEDSIGNATURES)

Calls the 'initMeemum' subroutine from 'perplexwrap.f' and initialize meemum using the local 'datFile' for the model parameters. 
This will return an array of 'Component' variables.
"""
function init_meemum(datfile)
    #function wrapper for the initmeemum subroutine in perplexwrap.f
    #Returns list of components with their composition
    filename = datfile
    
    #Any anticipated output from Fortran has to be put into an array of the same shape
    compositions = fill(0.0,3,K5)
    componentmass = fill(0.0,K0)
    componentnames = rpad("",K5*MAX_COMPONENT_NAME_L)#The array of strings fed from fortran is provided as just a string that will need to be parsed
    ccall((:__perplexwrap_MOD_initmeemum,joinpath(@__DIR__,"perplexwrap.so")),
        Cvoid,(Cstring,Ref{Int32}, Cstring, Ref{Float64},Ref{Float64}), 
        filename, sizeof(filename),componentnames, compositions, componentmass)
    
    #Iterate through the composition matrix and parse out component names into an array of components
    #Foreseeable issue, if a composition is input with a value of 0, this will break
   
    components = Array{Component}([])
    for row in eachrow(compositions)
        if row[1] > 0
            thiscol = Array{Component}([])
            for i in 1:lastindex(row)
                
                if row[i] > 0
                    #Parse the names assuming constant length of each component name
                    componentname = String(rstrip(componentnames[(i-1)*MAX_COMPONENT_NAME_L+1:i*MAX_COMPONENT_NAME_L]))
                    thiscomponent = Component(componentname,componentmass[i],row[i])
                    push!(thiscol,thiscomponent)
                    
                end
            end
            if length(components) == 0
                components = thiscol
            else
                components = hcat(components,thiscol)
            end
            
        end
       
    end
    global is_init = true
    return components
    
end


"""
$(TYPEDSIGNATURES)

This is function runs the 'minimizePoint' function in 'perplexwrap.f ' 
for the provided composition ('comps') at the given pressure ('pres') and temperature ('temp')  in bars and °C. 
This will return a PetroSystem.
"""
function minimizepoint(comps,temperature,pressure; suppresswarn= false, X = NaN, μ1 = NaN, μ2 = NaN)
    if !is_init
        throw(ErrorException("You must run init_meemum() before minimizepoint() can be used"))
    else
        #WARNING!!!!!!! DO NOT CHANGE ANYTHING BELOW THIS COMMENT IF YOU DO NOT KNOW WHAT YOU ARE DOING
        #INPUT variables
        temperature += 273 #Convert to K
        componentnames = name.(comps)
        componentstring = ""
        #Convert the array of names into a format readable by fortran
        for name in componentnames
            componentstring *= rpad(name,MAX_COMPONENT_NAME_L)
        end
        componentstring = rpad(componentstring,K5*MAX_COMPONENT_NAME_L)
        systemcomposition= vcat(concentration.(comps),fill(0.0,K5-length(comps)))
        
        #OUTPUT variables
        #Any anticipated output from Fortran has to be put into an array of the same shape
        
        chempotentials = fill(0.0,K8)
        phasenames = rpad("",K5*MAX_PHASE_NAME_L)
        phaseproperties = fill(0.0,I8,K5)
        phasecompositions = fill(0.0,K0,K5)
        systemproperties = fill(0.0,I8)
        
        
        ccall((:__perplexwrap_MOD_minimizepoint,joinpath(@__DIR__,"perplexwrap.so")),
            Cvoid,(Cstring,Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Bool},Ref{Float64},Cstring,Ref{Float64},Ref{Float64},Ref{Float64}),
            componentstring,systemcomposition,pressure,temperature,X,μ1,μ2,suppresswarn,chempotentials,phasenames,phaseproperties,phasecompositions,
            systemproperties)


        #WARNING!!!!!!! DO NOT CHANGE ANYTHING ABOVE THIS COMMENT IF YOU DO NOT KNOW WHAT YOU ARE DOING

        
        #Make a new array of components identical to the input, but modify the chemical potential value
        newcomponents = Array{Component}([])
        for i in 1:lastindex(comps)
        push!(newcomponents,Component(comps[i],μ=chempotentials[i]))
        end


        #Make an array of Phase objects with the properties calculated from PerpleX
        #Iterate through the phaseProperties and connect it with the phase name and composition
        phasearray = Array{Phase}([])
        for i in 1:K5
            #Check to make sure a phase exists at this index by looking at molar volume
            if phaseproperties[1,i] > 0
                iphasecomposition = Array{Component}([])
                iphasename = rstrip(phasenames[(i-1)*MAX_PHASE_NAME_L+1:i*MAX_PHASE_NAME_L])
                for j in 1:lastindex(comps)
                    push!(iphasecomposition,Component(comps[j],phasecompositions[j,i]))
                end
                #Assign appropriate values from the array, based on comments in perplex code
                iphase = Phase(name = iphasename, 
                                composition = iphasecomposition, 
                                mol = phaseproperties[16,i], 
                                vol = phaseproperties[16,i]*phaseproperties[1,i], 
                                mass = phaseproperties[16,i]*phaseproperties[17,i],
                                ρ =phaseproperties[10,i],
                                molarmass = phaseproperties[17,i],
                                G = phaseproperties[11,i],
                                H =phaseproperties[2,i],
                                S = phaseproperties[15,i],
                                Cp = phaseproperties[12,i],
                                Vmol = phaseproperties[1,i],
                                Cp_Cv = phaseproperties[28,i],
                                α = phaseproperties[13,i],
                                β = phaseproperties[14,i] )
                push!(phasearray,iphase)
            end
        end

        system = PetroSystem(composition = newcomponents,
                            phases = phasearray,
                            mol = systemproperties[16],
                            vol = systemproperties[16]*systemproperties[1],
                            mass = systemproperties[16]*systemproperties[17],
                            ρ = systemproperties[10],
                            molarmass = systemproperties[17],
                            G = systemproperties[11],
                            H = systemproperties[2],
                            S = systemproperties[15],
                            Cp = systemproperties[12],
                            Vmol = systemproperties[1],
                            Cp_Cv = systemproperties[28],
                            α = systemproperties[13],
                            β = systemproperties[14])
        #return compoNames,cPotentials,phaseNames,phaseProps,phaseComps,sysProps
        return system
    end
end

"""
$(SIGNATURES)
This is a simple type defined to keep track of the phases present at each x and y coordinate in a 'PerplexGrid'. 
The variables of the x and y axis will be defined in the 'PerplexGrid'
$(TYPEDFIELDS)
"""
struct Assemblage
    "List of phase names"
    phases::Array{String}
    "X coordinate"
    x::Float64
    "Y coordinate"
    y::Float64
    "Integer key that corresponds to the assemblage"
    key::Int64
end

"""
$(TYPEDSIGNATURES)

Returns the 'x' of 'asm', useful for broadcasting
"""
function x(asm::Assemblage)
    return asm.x
end

"""
$(TYPEDSIGNATURES)

Returns the 'y' of 'asm', useful for broadcasting
"""
function y(asm::Assemblage)
    return asm.y
end

"""
$(TYPEDSIGNATURES)

Returns the 'key' of 'asm', useful for broadcasting
"""
function key(asm::Assemblage)
    return asm.key
end

"""
$(TYPEDSIGNATURES)

Returns an array of 'Assemblage' variables constructed from 'asms' where each item in the array has a unique key
"""
function unique_assemblages(asms)

    present_assemblages = Array{Assemblage}([])
    

    for asm in asms
        ispresent = false
        for uasm in present_assemblages
            if asm.key == uasm.key
                ispresent = true
            end
        end

        if !ispresent
            push!(present_assemblages,asm)#x and y values are arbitrary here
        end

    end

    return sort(present_assemblages,by = key)

end

"""
$(SIGNATURES)
This describes the computed grid output by vertex (currently only the assemblage). Where each element of 
'assemblages' is a point on the grid.
$(TYPEDFIELDS)
"""
struct PerplexGrid
    "Entire collection of points from the computation"
    assemblages::Array{Assemblage}
    "Variable name of the x-axis"
    xaxis::String
    "Variable name of the y-axis"
    yaxis::String
end


"""
$(TYPEDSIGNATURES)

This calls the 'pseudosection' function in 'perplexwrap.f. This will read the output of a 
vertex calculation of the given 'datFile' and provide a 'PerplexGrid'. If 'iscelsius' is set to 'true', 
temperature variables will be converted to °C. If 'iskbar' is set to 'true', pressure variables 
will be converted to kBar.
"""
function get_pseudosection(datfile; iscelsius = false, iskbar = false)

    #Might want to convert return type to a DataFrame for easier use?

    #WARNING!!!!!!! DO NOT CHANGE ANYTHING BELOW THIS COMMENT IF YOU DO NOT KNOW WHAT YOU ARE DOING
    #Start by initializing all the necessary variables
    #INPUT VARIABLES
    filename = datfile
    #OUTPUT VARIABLES
    grid = fill(Int32(0),L7,L7)
    grid2assemblage = fill(Int32(0),K2)
    assemblage2phase = fill(Int32(0),K5,K3)
    #String arrays are just very long strings
    purephases = rpad("",K1*PURE_PHASE_NAME_L)
    solutionphases = rpad("",H9*SOL_PHASE_ABBREV_L)
    numphases = fill(Int32(0),K3)
    #Apparently floats have to be passed as an array, no idea why
    xmin = [0.0]
    xmax = [0.0]
    xincrement = [0.0]
    ymin = [0.0]
    ymax = [0.0]
    yincrement = [0.0]
    xvariable = rpad("",VARIABLE_NAME_L)
    yvariable= rpad("",VARIABLE_NAME_L)

    ccall((:__perplexwrap_MOD_pseudosection,joinpath(@__DIR__,"perplexwrap.so")),
        Cvoid,(Cstring,Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},Cstring,Cstring,Ref{Int32},Ref{Float64},Ref{Float64},
        Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Cstring,Cstring),
        filename,sizeof(filename),grid,grid2assemblage,assemblage2phase,purephases,solutionphases,numphases,
        xmin,xmax,ymin,ymax,xincrement,yincrement,xvariable,yvariable)

    #WARNING!!!!!!! DO NOT CHANGE ANYTHING ABOVE THIS COMMENT IF YOU DO NOT KNOW WHAT YOU ARE DOING
    #Converting to just floats because i dont wanna type the extra "[1]"
    xmin = xmin[1]
    xmax = xmax[1]
    xincrement = xincrement[1]
    ymin = ymin[1]
    ymax = ymax[1]
    yincrement = yincrement[1]

    xvariable = rstrip(xvariable)
    yvariable = rstrip(yvariable)
    
    #This is the K to °C conversion
    if occursin("T(K)", xvariable) && iscelsius
        xmin = xmin-273.15
        xmax = xmax-273.15
        xvariable = "T(°C)"
    end

    if occursin("T(K)",yvariable) && iscelsius
        ymin = ymin-273.15
        ymax = ymax-273.15
        yvariable  = "T(°C)"
    end

    #bar to kbar conversion
    if occursin("P(bar)", xvariable) && iskbar
        xmin = xmin/1000
        xmax = xmax/1000
        xincrement = xincrement/1000
        xvariable  = "P(kBar)"
    end

    if occursin("P(bar)",yvariable) && iskbar
        ymin = ymin/1000
        ymax = ymax/1000
        yincrement = yincrement/1000
        yvariable  = "P(kBar)"
    end
    #Convert purePhases to list of strings
    phasename = rstrip(purephases[1:PURE_PHASE_NAME_L])
    purephase_array = Array{String}([])
    count = 1
    while length(phasename) > 0 && count*PURE_PHASE_NAME_L < length(purephases)
        push!(purephase_array,phasename)
        firstindex= count*PURE_PHASE_NAME_L+1
        lastindex = (count+1)*PURE_PHASE_NAME_L
        phasename = rstrip(purephases[firstindex:lastindex])
        count += 1
    end

    #Convert solPhases to list of strings
    phasename = rstrip(solutionphases[1:SOL_PHASE_ABBREV_L-1])
    solutionphase_array = Array{String}([])
    count = 1
    while length(phasename) > 0 && count*SOL_PHASE_ABBREV_L < length(solutionphases)
        push!(solutionphase_array,phasename)
        firstindex = count*SOL_PHASE_ABBREV_L+1
        lastindex = (count+1)*SOL_PHASE_ABBREV_L
        phasename = rstrip(solutionphases[firstindex:lastindex])
        count += 1
    end


    gridded_assemblage = Array{Assemblage}([])
    #Decoding the grid and assigning assemblages
    for ix in CartesianIndices(grid)#Best way to iterate through a matrix apparently?
        #ix[1] is the x-axis, ix[2] is the y-axis and these are integer indices in grid

        if grid[ix[1],ix[2]] > 0
            
            #This is the conversion from grid points to x-y points
            x = xmin + xincrement*(ix[1]-1)
            y = ymin + yincrement*(ix[2]-1)
            assemblage_key = grid2assemblage[grid[ix[1],ix[2]]]#This tells us what the assemblage is
            phaselist_index = assemblage2phase[:,assemblage_key]#This gives us the list of phases in this assemblage
          
            phaselist = Array{String}([])
            
            for i in range(1,numphases[assemblage_key])
                #If the index is positive its a solution, if not its a pure phase
                if phaselist_index[i]> 0
                    push!(phaselist,solutionphase_array[phaselist_index[i]])
                elseif phaselist_index[i] < 0
                    push!(phaselist,purephase_array[-phaselist_index[i]])
                end
            end

            assemblage = Assemblage(phaselist,x,y,assemblage_key)
            push!(gridded_assemblage,assemblage)

        end

    end
    
    return PerplexGrid(gridded_assemblage,rstrip(xvariable),rstrip(yvariable))
end

"""
$(TYPEDSIGNATURES)

This is a function used for plotting. It will return a list of integers of the same length as 'asms'. 
Where each index corresponds to the same asm. These will have a value of 1 if 'asm[i]' is the same as 
the 'filterKey' and 0 if not.
"""
function filterkeys(asms,filterkey)
#Returns an array of 1s and 0s where 1 is in the index in asms that has a key that matches filterKey
#used for plotting
    keys = key.(asms)
    
    for i in range(1,lastindex(keys))

        if keys[i] != filterkey
            keys[i] = 0
        else
            keys[i] = 1
        end

    end

    return keys

end

"""
$(TYPEDSIGNATURES)

This will take 'pGrid' and return an array of all the 'Assemblage' variables in 'pGrid.assemblages' that match 
the 'filterKey'.
"""
function filtergrid(pgrid,filterkey)
#Returns a list of assemblages with a key matching filterKey
    assemblages = Array{Assemblage}([])

    for assem in pgrid.assemblages
        if assem.key == filterkey
            push!(assemblages,assem)
        end
    end

    return assemblages

end


# function plotPseudosection!(ax::Axis,pseudo::PerplexGrid)
   
#     x = getX.(pseudo.assemblages)
#     y = getY.(pseudo.assemblages)

#     #Here we figure out the assemblage variance used in determining the field colour
#     uniqueAsms = listUniqueAssemblages(pseudo.assemblages)
#     maxPhaseVar = 0
#     minPhaseVar = 100
#     for asm in uniqueAsms
#         if length(asm.phases) >maxPhaseVar
#             maxPhaseVar = length(asm.phases)
#         end
#         if length(asm.phases) <minPhaseVar
#             minPhaseVar = length(asm.phases)
#         end
#     end

#     #This part plots a filled contour around each unique assemblage
#     #Ranges of the levels and colormaps are based on trial and error, do not change them
#     #They should be roughly halfway beetween grid points in a smoothed line
#     for i in range(1,lastindex(uniqueAsms))
#         colorVal = 1-(1-(length(uniqueAsms[i].phases)-minPhaseVar)/(maxPhaseVar-minPhaseVar))*0.8
#         contourf!(ax,x,y,filterKeyArray(pseudo.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,Colors.HSV(0,0,colorVal)])
#     end

#     #For added flair, we add outlines and labels to each assemblages
#     for i in range(1,lastindex(uniqueAsms))
#         iGrid = filterGrid(pseudo,i)
#         scatter!(ax,mean(getX.(iGrid)),mean(getY.(iGrid)), marker = :circle, strokecolor = :black,strokewidth = 1,color = :transparent,markersize = 5)
#         text!(ax,mean(getX.(iGrid)),mean(getY.(iGrid)),text = string(i))
#     end
#     #Making the contours seperate so they can be selected easily in post-processing
#     for i in range(1,lastindex(uniqueAsms))
#         # println(string(uniqueAsms[i].key)*" = "*join(uniqueAsms[i].phases," "))
#         contour!(ax,x,y,filterKeyArray(pseudo.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
#     end

#     ax.xlabel = pseudo.xAx
#     ax.ylabel = pseudo.yAx
#     xMin = minimum(getX.(pseudo.assemblages))
#     yMin = minimum(getY.(pseudo.assemblages))
#     xMax = maximum(getX.(pseudo.assemblages))
#     yMax = maximum(getY.(pseudo.assemblages))
#     ax.limits = (xMin,xMax,yMin,yMax)
# end


"""
$(TYPEDSIGNATURES)

Creates a text file called 'filename.txt' with a list of each assemblage in 'pseudo'.
"""
function output_assemblages(filename,pseudo)
    present_assemblages = unique_assemblages(pseudo.assemblages)
    
    writefile = open(filename*".txt","w")

    for asm in present_assemblages
        
        writeline = string(asm.key)*" = "*join(asm.phases," ")*"\n"
        write(writefile, writeline)
    end

    close(writefile)
end

"""
$(TYPEDSIGNATURES)

Reads the output file of a werami run and returns a 'DataFrame' that can be used to create contour plots. If 'iscelsius' is set to 'true', 
temperature variables will be converted to °C. If 'iskbar' is set to 'true', pressure variables will be converted to kBar.
"""
function read_werami_output(filename; iscelsius = false, iskbar = false)

    readfile = open(filename,"r")
    lines = readlines(readfile)

    #Assumes a specific file structure with variable headings at line 13
    headings = split(lines[13])
    #Set up arrays to read into
    xs = Array{Float64}([])
    ys = Array{Float64}([])
    variables = Array{Array{Float64}}([])
    for i in range(3,lastindex(headings))
        push!(variables,[])
    end

    #Parsing out the data
    for i in range(14,lastindex(lines))
        vals= split(lines[i])
        push!(xs,parse(Float64,vals[1]))
        push!(ys,parse(Float64,vals[2]))

        for j in range(3,lastindex(vals))
            val = parse(Float64,vals[j])
            if isnan(val)
                val = 0.0 #Problems when plotting with NaN with contours
            end
            push!(variables[j-2],val)
        end
    end

    #Put it all in a DataFrame
    df = DataFrame(headings[1]=>xs,headings[2]=>ys)
    for i in range(3,lastindex(headings))
        
        df[:,Symbol(headings[i])] = variables[i-2]
    end

    if iscelsius && hasproperty(df,"T(K)")
        rename!(df,Symbol("T(K)") => "T(°C)")
        df[!,Symbol("T(°C)")] .-= 273.15
    end

    if iskbar && hasproperty(df,"P(bar)")
        rename!(df,Symbol("P(bar)") => "P(kbar)")
        df[!,Symbol("P(kbar)")] ./= 1000
    end
    
    return df
end


#Generic functions for plotting extensions

function pseudosection end
function pseudosection! end
end
