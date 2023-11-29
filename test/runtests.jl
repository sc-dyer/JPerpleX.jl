#Tests adapted from main Perple_X repo
using JPerpleX
using CairoMakie
using Test

@testset "JPerplex.jl" begin
    # Write your tests here.
    #T
    @testset "Minimization test" begin

        readCompo = initMeemum("bl691_MEEMUM/bl691_MEEMUM")
        Na2O = Component("Na2O",61.9790,2.5416)
        MgO = Component("MgO",40.3040,8.110)
        Al2O3 = Component("Al2O3",101.9610,9.2646)
        K2O = Component("K2O",94.1960,0.1054)
        CaO = Component("CaO",56.0770,10.146)
        TiO2 = Component("TiO2",79.8660,1.346)
        FeO = Component("FeO",71.8440,10.1382)
        O2 = Component("O2",31.9990,0.488 )
        H2O = Component("H2O",18.0150,3.4162)
        SiO2 = Component("SiO2",60.0840,53.9559)

        testCompo = [Na2O,MgO,Al2O3,K2O,CaO,TiO2,FeO,O2,H2O,SiO2]
        for elem in testCompo
            @test elem ≈ readCompo[findChem(readCompo,elem)]
        end
        println("Init test done")
        T = 1050-273
        P = 4000
        calcSys = minimizePoint(readCompo, T, P)

        # Check a few variables to compare with expected output
        # println(calcSys.phases)
        @test calcSys.mMass ≈ 6162.817774
        @test calcSys.G ≈ -85883543
        @test calcSys.S ≈  12574.16208

        @test calcSys.phases[1].name == "cAmph(G)"
        @test calcSys.phases[1].mMass ≈  899.28334
        @test calcSys.phases[1].G ≈ -1.18096964e7
        
    end

    @testset "Pseudosection test" begin


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

        
        pseudo = getPseudosection("bl691/bl691",tempInC=true, pInKBar = true)

        uniqueAsm = listUniqueAssemblages(pseudo.assemblages)
        @test length(uniqueAsm) == 48
        @test uniqueAsm[8].key == 8
        plotPseudosection!(ax,pseudo)
        save("bl691/bl691.svg",fig)

        @test filesize("bl691/bl691.svg") ≈ filesize("bl691/output/bl691.svg")
        outputAssemblages("bl691/bl691_assemblages",pseudo)
        @test filesize("bl691/bl691_assemblages.txt") ≈ filesize("bl691/output/bl691_assemblages.txt")

        pseudo2 = getPseudosection("klb691/klb691",tempInC=true,pInKBar=true)
        empty!(ax)
        plotPseudosection!(ax,pseudo2)
        save("klb691/klb691.svg",fig)
        @test filesize("klb691/klb691.svg") ≈ filesize("klb691/output/klb691.svg")
        
    end
end
