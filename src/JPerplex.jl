
"""
This is the main module for a more "advanced" Perple_X julia wrapper. The advantages of this
library over a more "basic" wrapper is that it calls the compiled fortran functions in perplexwrap.f
directly and *shouldn't* require repeated file i/o or command piping.

# Exports
$(EXPORTS)
"""
module JPerplex
export 
    initMeemum,
    minimizePoint,
    Assemblage,
    PerplexGrid,
    getPseudosection
using
    DocStringExtensions,
    Reexport

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
                    compName = rstrip(componentNames[(i-1)*maxCompNameL+1:i*maxCompNameL])
                    thisComponent = Component(compName,componentMass[i],row[i],0)
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
function minimizePoint(comps::Array{Component},pres::Real,temp::Real; suppressWarn::Bool = false)

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
    sysCompo = vcat(mol.(comps),fill(0.0,K5-length(comps)))
    
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

    system = SystemProperties(compo = newComps,
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
    return phaseArray, system
end


struct Assemblage
    phases::Array{String}
    x::Float64
    y::Float64
    key::Int64
end

struct PerplexGrid
    assemblages::Array{Assemblage}
    xAx::String
    yAx::String
end

function getPseudosection(datFile::String)

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
    println(purePhaseArr[60:150])
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

end
