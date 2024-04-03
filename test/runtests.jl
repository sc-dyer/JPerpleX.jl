#Tests adapted from main Perple_X repo
using JPerpleX
using CairoMakie
using Test

@testset "JPerplex.jl" begin
    # Write your tests here.
    #T
    @testset "Minimization test" begin

        readcompo = init_meemum("bl691_MEEMUM/bl691_MEEMUM")
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

        testcompo = [Na2O,MgO,Al2O3,K2O,CaO,TiO2,FeO,O2,H2O,SiO2]
        for elem in testcompo
            @test elem ≈ readcompo[findchemical(readcompo,elem)]
        end
        println("Init test done")
        T = 1050-273
        P = 4000
        calc_sys = minimizepoint(readcompo, T, P)

        # Check a few variables to compare with expected output
        # println(calcSys.phases)
        @test calc_sys.molarmass ≈ 6162.817774
        @test calc_sys.G ≈ -85883543
        @test calc_sys.S ≈  12574.16208

        @test calc_sys.phases[1].name == "cAmph(G)"
        @test calc_sys.phases[1].molarmass ≈  899.28334
        @test calc_sys.phases[1].G ≈ -1.18096964e7
        
    end

    @testset "Pseudosection test" begin

        include("PlotDefaults.jl")
       
        set_theme!(myTheme)
        fig = Figure()
        ax = Axis(fig[1,1])
        pseudo = get_pseudosection("bl691/bl691",iscelsius=true, iskbar = true)

        present_assemblages= unique_assemblages(pseudo.assemblages)
        @test length(present_assemblages) == 48
        @test present_assemblages[8].key == 8
        pseudosection!(ax,pseudo)
        ax.xlabel = pseudo.xaxis
        ax.ylabel = pseudo.yaxis
        xmin = minimum(x.(pseudo.assemblages))
        ymin = minimum(y.(pseudo.assemblages))
        xmax = maximum(x.(pseudo.assemblages))
        ymax = maximum(y.(pseudo.assemblages))
        ax.limits = (xmin,xmax,ymin,ymax)
        save("bl691/bl691.svg",fig)

        @test filesize("bl691/bl691.svg") ≈ filesize("bl691/output/bl691.svg")
        output_assemblages("bl691/bl691_assemblages",pseudo)
        @test filesize("bl691/bl691_assemblages.txt") ≈ filesize("bl691/output/bl691_assemblages.txt")

        pseudo2 = get_pseudosection("klb691/klb691",iscelsius=true,iskbar=true)
        empty!(ax)
        ax.xlabel = pseudo2.xaxis
        ax.ylabel = pseudo2.yaxis
        xmin = minimum(x.(pseudo2.assemblages))
        ymin = minimum(y.(pseudo2.assemblages))
        xmax = maximum(x.(pseudo2.assemblages))
        ymax = maximum(y.(pseudo2.assemblages))
        ax.limits = (xmin,xmax,ymin,ymax)
        pseudosection!(ax,pseudo2)
        save("klb691/klb691.svg",fig)
        @test filesize("klb691/klb691.svg") ≈ filesize("klb691/output/klb691.svg")
        
    end
end
