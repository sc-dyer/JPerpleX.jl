
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
    plotPseudosection!
using
    DocStringExtensions,
    Reexport,
    CairoMakie,
    Statistics,
    Makie.Colors

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
    initMeemum(datFile)
Calls the initMeemum subroutine from perplexwrap.f and initialize meemum using the local datFile for the model parameters.
This will return a matrix of Components providing the bulk rock composition (in mol), component name, and its molar mass.
Each column in the matrix represents one set of compositions provided in the .dat file.
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
    minimizePoint(comps, pres, temp)
This is function runs meemum for the provided composition (comps) at the given pres and temp in bars and degree C. 
    This will return a list of stable phases and their properties as well as the bulk system properties
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


struct Assemblage
    phases::Array{String}
    x::Float64
    y::Float64
    key::Int64
end

function getX(asm::Assemblage)
    return asm.x
end

function getY(asm::Assemblage)
    return asm.y
end

function getKey(asm::Assemblage)
    return asm.key
end

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
struct PerplexGrid
    assemblages::Array{Assemblage}
    xAx::String
    yAx::String
end

function getPseudosection(datFile::String; tempInC::Bool = false, pInKBar = false)

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
        Cvoid,(Cstring,Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},Cstring,Cstring,Ref{Float64},Ref{Float64},
        Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Cstring,Cstring),
        fileName,sizeof(fileName),grid,gridToAssem,assemToPhase,purePhases,solPhases,
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
    
    if occursin("T(K)", xVarName) && tempInC
        xMin = xMin-273.15
        xMax = xMax-273.15
        xVarName  = "T (°C)"
    end

    if occursin("T(K)",yVarName) && tempInC
        yMin = yMin-273.15
        yMax = yMax-273.15
        yVarName  = "T (°C)"
    end

    if occursin("P(bar)", xVarName) && pInKBar
        xMin = xMin/1000
        xMax = xMax/1000
        xInc = xInc/1000
        xVarName  = "P (kBar)"
    end

    if occursin("P(bar)",yVarName) && pInKBar
        yMin = yMin/1000
        yMax = yMax/1000
        yInc = yInc/1000
        yVarName  = "P (kBar)"
    end
    #Convert purePhases to list of strings
    phaseName = rstrip(purePhases[1:purePhaseNameL])
    # println(purePhases[1:100])
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
    phaseName = phaseName = rstrip(solPhases[1:solPhaseAbbrL-1])
    solPhaseArr = Array{String}([])
    count = 1
    while length(phaseName) > 0 && count*solPhaseAbbrL < length(solPhases)
        push!(solPhaseArr,phaseName)
        firstIndex = count*solPhaseAbbrL+1
        lastIndex = (count+1)*solPhaseAbbrL
        phaseName = rstrip(solPhases[firstIndex:lastIndex])
        count += 1
    end
    # println(purePhaseArr[60:150])
    griddedAssemblage = Array{Assemblage}([])
    #Decoding the grid and assigning assemblages
    for ix in CartesianIndices(grid)

        if grid[ix[1],ix[2]] > 0
            
            #Might be reversed, got to try and find out I guess
            x = xMin + xInc*(ix[1]-1)
            y = yMin + yInc*(ix[2]-1)
            asmKey = gridToAssem[grid[ix[1],ix[2]]]#This tells us what the assemblage is
            phaseListIndex = assemToPhase[:,asmKey]#This gives us the list of phases in this assemblage
            phaseList = Array{String}([])
            
            for index in phaseListIndex
                #If the index is positive its a solution, if not its a pure phase
                if index > 0
                    push!(phaseList,solPhaseArr[index])
                elseif index < 0
                    push!(phaseList,purePhaseArr[-index])
                end
            end

            asm = Assemblage(phaseList,x,y,asmKey)
            push!(griddedAssemblage,asm)

        end

    end
    
    return PerplexGrid(griddedAssemblage,rstrip(xVarName),rstrip(yVarName))
end

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


function plotPseudosection!(ax::Axis,pseudo::PerplexGrid)
    #Will modify the axis provided to make a getPseudosection
    #User will define design elements
    
    #For now do heatmap
    x = getX.(pseudo.assemblages)
    y = getY.(pseudo.assemblages)
    xAx = pseudo.xAx
    yAx = pseudo.yAx

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


    for i in range(1,lastindex(uniqueAsms))
        
        colorVal = 1-(1-(length(uniqueAsms[i].phases)-minPhaseVar)/(maxPhaseV-minPhaseV))*0.8
        contourf!(ax,x,y,filterKeyArray(pseudo.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,Colors.HSV(0,0,colorVal)])
        
    end

    for i in range(1,lastindex(uniqueAsms))
        contour!(ax,x,y,filterKeyArray(pseudo.assemblages,i),levels =-0.5:1:1.5,colormap = [:transparent,:black,:black],linewidth=2)
        iGrid = filterGrid(pseudo,i)
        scatter!(ax,mean(getX.(iGrid)),mean(getY.(iGrid)), marker = :circle, strokecolor = :black,strokewidth = 1,color = :transparent)
        text!(ax,mean(getX.(iGrid)),mean(getY.(iGrid)),text = string(i))
    end

    ax.xlabel = xAx
    ax.ylabel = yAx
    xMin = minimum(getX.(pseudo.assemblages))
    yMin = minimum(getY.(pseudo.assemblages))
    xMax = maximum(getX.(pseudo.assemblages))
    yMax = maximum(getY.(pseudo.assemblages))
    ax.limits = (xMin,xMax,yMin,yMax)
end

end
