
"""
This is the main module for a more "advanced" Perple_X julia wrapper. The advantages of this 
library over a more "basic" wrapper is that it calls the compiled fortran functions in perplexwrap.f 
directly and *shouldn't* require repeated file i/o or command piping.

# Exports
$(EXPORTS)
"""
module JPerpleX
export 
    initMeemum,
    minimizePoint,
    Assemblage,
    PerplexGrid,
    getPseudosection,
    listUniqueAssemblages,
    getX,
    getY,
    getKey,
    filterGrid,
    filterKeyArray,
    plotPseudosection!,
    outputAssemblages,
    readWeramiOutput
using
    DocStringExtensions,
    Reexport,
    CairoMakie,
    Statistics,
    Makie.Colors,
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
const maxCompNameL = 5
const maxPhaseNameL = 14
const purePhaseNameL = 8
const solPhaseAbbrL = 6
const varNameL = 8

"""
$(TYPEDSIGNATURES)

Calls the 'initMeemum' subroutine from 'perplexwrap.f' and initialize meemum using the local 'datFile' for the model parameters. 
This will return an array of 'Component' variables.
"""
function initMeemum(datFile::String)
    #function wrapper for the initmeemum subroutine in perplexwrap.f
    #Returns list of components with their composition
    fileName = datFile
    
    #Any anticipated output from Fortran has to be put into an array of the same shape
    compositions = fill(0.0,3,K5)
    componentMass = fill(0.0,K0)
    componentNames = rpad("",K5*maxCompNameL)#The array of strings fed from fortran is provided as just a string that will need to be parsed
    ccall((:__perplexwrap_MOD_initmeemum,joinpath(@__DIR__,"perplexwrap.so")),
        Cvoid,(Cstring,Ref{Int32}, Cstring, Ref{Float64},Ref{Float64}), 
        fileName, sizeof(fileName),componentNames, compositions, componentMass)
    
    #Iterate through the composition matrix and parse out component names into an array of components
    #Foreseeable issue, if a composition is input with a value of 0, this will break
   
    components = Array{Component}([])
    for row in eachrow(compositions)
        if row[1] > 0
            thisCol = Array{Component}([])
            for i in 1:lastindex(row)
                
                if row[i] > 0
                    #Parse the names assuming constant length of each component name
                    compName = String(rstrip(componentNames[(i-1)*maxCompNameL+1:i*maxCompNameL]))
                    thisComponent = Component(compName,componentMass[i],row[i])
                    push!(thisCol,thisComponent)
                    
                end
            end
            if length(components) == 0
                components = thisCol
            else
                components = hcat(components,thisRow)
            end
            
        end
       
    end

    return components
    
end


"""
$(TYPEDSIGNATURES)

This is function runs the 'minimizePoint' function in 'perplexwrap.f ' 
for the provided composition ('comps') at the given pressure ('pres') and temperature ('temp')  in bars and °C. 
This will return a PetroSystem.
"""
function minimizePoint(comps::Array{Component},temp::Real,pres::Real; suppressWarn::Bool = false)

    #WARNING!!!!!!! DO NOT CHANGE ANYTHING BELOW THIS COMMENT IF YOU DO NOT KNOW WHAT YOU ARE DOING
    #INPUT variables
    temp += 273 #Convert to K
    compoNames = name.(comps)
    compoString = ""
    #Convert the array of names into a format readable by fortran
    for name in compoNames
        compoString *= rpad(name,maxCompNameL)
    end
    compoString = rpad(compoString,K5*maxCompNameL)
    sysCompo = vcat(conc.(comps),fill(0.0,K5-length(comps)))
    
    #OUTPUT variables
    #Any anticipated output from Fortran has to be put into an array of the same shape
    
    cPotentials = fill(0.0,K8)
    phaseNames = rpad("",K5*maxPhaseNameL)
    phaseProps = fill(0.0,I8,K5)
    phaseComps = fill(0.0,K0,K5)
    sysProps = fill(0.0,I8)
    
    ccall((:__perplexwrap_MOD_minimizepoint,joinpath(@__DIR__,"perplexwrap.so")),
        Cvoid,(Cstring,Ref{Float64},Ref{Float64},Ref{Float64},Ref{Bool},Ref{Float64},Cstring,Ref{Float64},Ref{Float64},Ref{Float64}),
        compoString,sysCompo,pres,temp,suppressWarn,cPotentials,phaseNames,phaseProps,phaseComps,sysProps)


    #WARNING!!!!!!! DO NOT CHANGE ANYTHING ABOVE THIS COMMENT IF YOU DO NOT KNOW WHAT YOU ARE DOING

    
    #Make a new array of components identical to the input, but modify the chemical potential value
    newComps = Array{Component}([])
    for i in 1:lastindex(comps)
       push!(newComps,Component(comps[i],μ=cPotentials[i]))
    end


    #Make an array of Phase objects with the properties calculated from PerpleX
    #Iterate through the phaseProperties and connect it with the phase name and composition
    phaseArray = Array{Phase}([])
    for i in 1:K5
        #Check to make sure a phase exists at this index by looking at molar volume
        if phaseProps[1,i] > 0
            phaseCompo = Array{Component}([])
            phaseName = rstrip(phaseNames[(i-1)*maxPhaseNameL+1:i*maxPhaseNameL])
            for j in 1:lastindex(comps)
                push!(phaseCompo,Component(comps[j],phaseComps[j,i]))
            end
            #Assign appropriate values from the array, based on comments in perplex code
            myPhase = Phase(name = phaseName, 
                            compo = phaseCompo, 
                            mol = phaseProps[16,i], 
                            vol = phaseProps[16,i]*phaseProps[1,i], 
                            mass = phaseProps[16,i]*phaseProps[17,i],
                            ρ =phaseProps[10,i],
                            mMass = phaseProps[17,i],
                            G = phaseProps[11,i],
                            H =phaseProps[2,i],
                            S = phaseProps[15,i],
                            Cp = phaseProps[12,i],
                            Vmol = phaseProps[1,i],
                            Cp_Cv = phaseProps[28,i],
                            α = phaseProps[13,i],
                            β = phaseProps[14,i] )
            push!(phaseArray,myPhase)
        end
    end

    system = PetroSystem(compo = newComps,
                        phases = phaseArray,
                        mol = sysProps[16],
                        vol = sysProps[16]*sysProps[1],
                        mass = sysProps[16]*sysProps[17],
                        ρ = sysProps[10],
                        mMass = sysProps[17],
                        G = sysProps[11],
                        H = sysProps[2],
                        S = sysProps[15],
                        Cp = sysProps[12],
                        Vmol = sysProps[1],
                        Cp_Cv = sysProps[28],
                        α = sysProps[13],
                        β = sysProps[14])
    #return compoNames,cPotentials,phaseNames,phaseProps,phaseComps,sysProps
    return system
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
function getX(asm::Assemblage)
    return asm.x
end

"""
$(TYPEDSIGNATURES)

Returns the 'y' of 'asm', useful for broadcasting
"""
function getY(asm::Assemblage)
    return asm.y
end

"""
$(TYPEDSIGNATURES)

Returns the 'key' of 'asm', useful for broadcasting
"""
function getKey(asm::Assemblage)
    return asm.key
end

"""
$(TYPEDSIGNATURES)

Returns an array of 'Assemblage' variables constructed from 'asms' where each item in the array has a unique key
"""
function listUniqueAssemblages(asms::Array{Assemblage})

    uniqueAsms = Array{Assemblage}([])
    

    for asm in asms
        isPresent = false
        for uAsm in uniqueAsms
            if asm.key == uAsm.key
                isPresent = true
            end
        end

        if !isPresent
            push!(uniqueAsms,asm)#x and y values are arbitrary here
        end

    end

    return sort(uniqueAsms,by = getKey)

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
    xAx::String
    "Variable name of the y-axis"
    yAx::String
end


"""
$(TYPEDSIGNATURES)

This calls the 'pseudosection' function in 'perplexwrap.f. This will read the output of a 
vertex calculation of the given 'datFile' and provide a 'PerplexGrid'. If 'tempInC' is set to 'true', 
temperature variables will be converted to °C. If 'pInkBar' is set to 'true', pressure variables 
will be converted to kBar.
"""
function getPseudosection(datFile::String; tempInC::Bool = false, pInKBar::Bool = false)

    #Might want to convert return type to a DataFrame for easier use?

    #WARNING!!!!!!! DO NOT CHANGE ANYTHING BELOW THIS COMMENT IF YOU DO NOT KNOW WHAT YOU ARE DOING
    #Start by initializing all the necessary variables
    #INPUT VARIABLES
    fileName = datFile
    #OUTPUT VARIABLES
    grid = fill(Int32(0),L7,L7)
    gridToAssem = fill(Int32(0),K2)
    assemToPhase = fill(Int32(0),K5,K3)
    #String arrays are just very long strings
    purePhases = rpad("",K1*purePhaseNameL)
    solPhases = rpad("",H9*solPhaseAbbrL)
    numPhases = fill(Int32(0),K3)
    #Apparently floats have to be passed as an array, no idea why
    xMin = [0.0]
    xMax = [0.0]
    xInc = [0.0]
    yMin = [0.0]
    yMax = [0.0]
    yInc = [0.0]
    xVarName = rpad("",varNameL)
    yVarName = rpad("",varNameL)

    ccall((:__perplexwrap_MOD_pseudosection,joinpath(@__DIR__,"perplexwrap.so")),
        Cvoid,(Cstring,Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},Cstring,Cstring,Ref{Int32},Ref{Float64},Ref{Float64},
        Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Cstring,Cstring),
        fileName,sizeof(fileName),grid,gridToAssem,assemToPhase,purePhases,solPhases,numPhases,
        xMin,xMax,yMin,yMax,xInc,yInc,xVarName,yVarName)

    #WARNING!!!!!!! DO NOT CHANGE ANYTHING ABOVE THIS COMMENT IF YOU DO NOT KNOW WHAT YOU ARE DOING
    #Converting to just floats because i dont wanna type the extra "[1]"
    xMin = xMin[1]
    xMax = xMax[1]
    xInc = xInc[1]
    yMin = yMin[1]
    yMax = yMax[1]
    yInc = yInc[1]

    xVarName = rstrip(xVarName)
    yVarName = rstrip(yVarName)
    
    #This is the K to °C conversion
    if occursin("T(K)", xVarName) && tempInC
        xMin = xMin-273.15
        xMax = xMax-273.15
        xVarName  = "T(°C)"
    end

    if occursin("T(K)",yVarName) && tempInC
        yMin = yMin-273.15
        yMax = yMax-273.15
        yVarName  = "T(°C)"
    end

    #bar to kbar conversion
    if occursin("P(bar)", xVarName) && pInKBar
        xMin = xMin/1000
        xMax = xMax/1000
        xInc = xInc/1000
        xVarName  = "P(kBar)"
    end

    if occursin("P(bar)",yVarName) && pInKBar
        yMin = yMin/1000
        yMax = yMax/1000
        yInc = yInc/1000
        yVarName  = "P(kBar)"
    end
    #Convert purePhases to list of strings
    phaseName = rstrip(purePhases[1:purePhaseNameL])
    purePhaseArr = Array{String}([])
    count = 1
    while length(phaseName) > 0 && count*purePhaseNameL < length(purePhases)
        push!(purePhaseArr,phaseName)
        firstIndex = count*purePhaseNameL+1
        lastIndex = (count+1)*purePhaseNameL
        phaseName = rstrip(purePhases[firstIndex:lastIndex])
        count += 1
    end

    #Convert solPhases to list of strings
    phaseName = rstrip(solPhases[1:solPhaseAbbrL-1])
    solPhaseArr = Array{String}([])
    count = 1
    while length(phaseName) > 0 && count*solPhaseAbbrL < length(solPhases)
        push!(solPhaseArr,phaseName)
        firstIndex = count*solPhaseAbbrL+1
        lastIndex = (count+1)*solPhaseAbbrL
        phaseName = rstrip(solPhases[firstIndex:lastIndex])
        count += 1
    end


    griddedAssemblage = Array{Assemblage}([])
    #Decoding the grid and assigning assemblages
    for ix in CartesianIndices(grid)#Best way to iterate through a matrix apparently?
        #ix[1] is the x-axis, ix[2] is the y-axis and these are integer indices in grid

        if grid[ix[1],ix[2]] > 0
            
            #This is the conversion from grid points to x-y points
            x = xMin + xInc*(ix[1]-1)
            y = yMin + yInc*(ix[2]-1)
            asmKey = gridToAssem[grid[ix[1],ix[2]]]#This tells us what the assemblage is
            phaseListIndex = assemToPhase[:,asmKey]#This gives us the list of phases in this assemblage
          
            phaseList = Array{String}([])
            
            for i in range(1,numPhases[asmKey])
                #If the index is positive its a solution, if not its a pure phase
                if phaseListIndex[i]> 0
                    push!(phaseList,solPhaseArr[phaseListIndex[i]])
                elseif phaseListIndex[i] < 0
                    push!(phaseList,purePhaseArr[-phaseListIndex[i]])
                end
            end

            asm = Assemblage(phaseList,x,y,asmKey)
            push!(griddedAssemblage,asm)

        end

    end
    
    return PerplexGrid(griddedAssemblage,rstrip(xVarName),rstrip(yVarName))
end

"""
$(TYPEDSIGNATURES)

This is a function used for plotting. It will return a list of integers of the same length as 'asms'. 
Where each index corresponds to the same asm. These will have a value of 1 if 'asm[i]' is the same as 
the 'filterKey' and 0 if not.
"""
function filterKeyArray(asms::Array{Assemblage},filterKey::Integer)
#Returns an array of 1s and 0s where 1 is in the index in asms that has a key that matches filterKey
#used for plotting
    keys = getKey.(asms)
    
    for i in range(1,lastindex(keys))

        if keys[i] != filterKey
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
function filterGrid(pGrid::PerplexGrid,filterKey::Integer)
#Returns a list of assemblages with a key matching filterKey
    asms = Array{Assemblage}([])

    for assem in pGrid.assemblages
        if assem.key == filterKey
            push!(asms,assem)
        end
    end

    return asms

end

"""
$(TYPEDSIGNATURES)

This will take 'ax' and plot the pseudosection defined by 'pseudo'. These will be plotted such that they 
can still be edited with a vector graphics editor. Fields will be labeled by assemblage key, a list of unique 
assemblages can be retrieved with 'getUniqueAssemblages' or they can be printed to a text file with 'outputAssemblages'. 
The 'Axis' type is defined in the CairoMakie package, but it might work with Plots? The user should be able to 
modify the formatting of 'ax' before and after calling this function.
"""
function plotPseudosection!(ax::Axis,pseudo::PerplexGrid)
   
    x = getX.(pseudo.assemblages)
    y = getY.(pseudo.assemblages)

    #Here we figure out the assemblage variance used in determining the field colour
    uniqueAsms = listUniqueAssemblages(pseudo.assemblages)
    maxPhaseVar = 0
    minPhaseVar = 100
    for asm in uniqueAsms
        if length(asm.phases) >maxPhaseVar
            maxPhaseVar = length(asm.phases)
        end
        if length(asm.phases) <minPhaseVar
            minPhaseVar = length(asm.phases)
        end
    end

    #This part plots a filled contour around each unique assemblage
    #Ranges of the levels and colormaps are based on trial and error, do not change them
    #They should be roughly halfway beetween grid points in a smoothed line
    for i in range(1,lastindex(uniqueAsms))
        colorVal = 1-(1-(length(uniqueAsms[i].phases)-minPhaseVar)/(maxPhaseVar-minPhaseVar))*0.8
        contourf!(ax,x,y,filterKeyArray(pseudo.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,Colors.HSV(0,0,colorVal)])
    end

    #For added flair, we add outlines and labels to each assemblage
    
    for i in range(1,lastindex(uniqueAsms))
        # println(string(uniqueAsms[i].key)*" = "*join(uniqueAsms[i].phases," "))
        contour!(ax,x,y,filterKeyArray(pseudo.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
        iGrid = filterGrid(pseudo,i)
        scatter!(ax,mean(getX.(iGrid)),mean(getY.(iGrid)), marker = :circle, strokecolor = :black,strokewidth = 1,color = :transparent,markersize = 5)
        text!(ax,mean(getX.(iGrid)),mean(getY.(iGrid)),text = string(i))
    end

    ax.xlabel = pseudo.xAx
    ax.ylabel = pseudo.yAx
    xMin = minimum(getX.(pseudo.assemblages))
    yMin = minimum(getY.(pseudo.assemblages))
    xMax = maximum(getX.(pseudo.assemblages))
    yMax = maximum(getY.(pseudo.assemblages))
    ax.limits = (xMin,xMax,yMin,yMax)
end

"""
$(TYPEDSIGNATURES)

Creates a text file called 'fileName.txt' with a list of each assemblage in 'pseudo'.
"""
function outputAssemblages(fileName::String,pseudo::PerplexGrid)
    uniqueAsms = listUniqueAssemblages(pseudo.assemblages)
    
    writeFile = open(fileName*".txt","w")

    for asm in uniqueAsms
        
        wLine = string(asm.key)*" = "*join(asm.phases," ")*"\n"
        write(writeFile, wLine)
    end

    close(writeFile)
end

"""
$(TYPEDSIGNATURES)

Reads the output file of a werami run and returns a 'DataFrame' that can be used to create contour plots. If 'tempInC' is set to 'true', 
temperature variables will be converted to °C. If 'pInkBar' is set to 'true', pressure variables will be converted to kBar.
"""
function readWeramiOutput(fileName::String; tempInC::Bool = false, pInKBar::Bool = false)

    readFile = open(fileName,"r")
    lines = readlines(readFile)

    #Assumes a specific file structure with variable headings at line 13
    headings = split(lines[13])
    #Set up arrays to read into
    xVals = Array{Float64}([])
    yVals = Array{Float64}([])
    varVals = Array{Array{Float64}}([])
    for i in range(3,lastindex(headings))
        push!(varVals,[])
    end

    #Parsing out the data
    for i in range(14,lastindex(lines))
        vals = split(lines[i])
        push!(xVals,parse(Float64,vals[1]))
        push!(yVals,parse(Float64,vals[2]))

        for j in range(3,lastindex(vals))
            val = parse(Float64,vals[j])
            if isnan(val)
                val = 0.0 #Problems when plotting with NaN with contours
            end
            push!(varVals[j-2],val)
        end
    end

    #Put it all in a DataFrame
    df = DataFrame(headings[1]=>xVals,headings[2]=>yVals)
    for i in range(3,lastindex(headings))
        
        df[:,Symbol(headings[i])] = varVals[i-2]
    end

    if tempInC && hasproperty(df,"T(K)")
        rename!(df,Symbol("T(K)") => "T(°C)")
        df[!,Symbol("T(°C)")] .-= 273.15
    end

    if pInKBar && hasproperty(df,"P(bar)")
        rename!(df,Symbol("P(bar)") => "P(kbar)")
        df[!,Symbol("P(kbar)")] ./= 1000
    end
    
    return df
end

end
