#Tests adapted from main Perple_X repo
using JPerpleX
using CairoMakie
using Test

@testset "JPerplex.jl" begin
    # Write your tests here.
    #T
    @testset "Minimization test" begin

        bl691 = init_meemum("bl691_MEEMUM/bl691_MEEMUM")
        readcompo = getcompo(bl691)
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
        calc_sys = minimizepoint(bl691, T, P)

        # Check a few variables to compare with expected output
        # println(calcSys.phases)
        @test calc_sys.molarmass ≈ 247.73907
        @test calc_sys.G ≈ -85883543
        @test calc_sys.S ≈  12574.16208

        @test calc_sys.phases[1].name == "cAmph(G)"
        @test calc_sys.phases[1].molarmass ≈  899.28334
        @test calc_sys.phases[1].G ≈ -1.18096964e7

        close_meemum!(bl691)

        h2o_23SD20A = init_meemum("23SD20A/23SD20A_uH2O")
        readcompo = getcompo(h2o_23SD20A)
        SiO2 = Component("SiO2",60.0840,73.9801)
        Al2O3 = Component("Al2O3",101.9610,8.89559)
        FeO = Component("FeO",71.8440,3.67942)
        MgO = Component("MgO",40.3040,1.06316)
        CaO = Component("CaO",56.0770,2.44288)
        Na2O = Component("Na2O",61.9790,2.98545)
        K2O = Component("K2O",94.1960,4.13531)
        TiO2 = Component("TiO2",79.8660,0.537338)
        O2 = Component("O2",31.9990,0.422828)
        testcompo = [SiO2, Al2O3, FeO, MgO, CaO, Na2O, K2O, TiO2, O2]

        for elem in testcompo
            @test elem ≈ readcompo[findchemical(readcompo,elem)]
        end
        @test findchemical(readcompo,"H2O") == 0
        println("Init test 2 done")

        T = 850
        P = 9000
        μ = -320000
        calc_sys = minimizepoint(h2o_23SD20A, T, P, μ1 = μ)

        h2o_test = Component("H2O",18.0150,39.1751259,μ = -320000)
        # Check a few variables to compare with expected output
        # println(calcSys.phases)
        syscompo = calc_sys.composition
        @test h2o_test ≈ syscompo[findchemical(syscompo,h2o_test)]
        for elem in testcompo
            @test elem ≈ syscompo[findchemical(syscompo,elem)]
        end
        
        @test round(calc_sys.molarmass,digits=2) ≈ 110.81
        @test round(calc_sys.G,sigdigits = 4)  ≈ -9.281e7
        @test round(calc_sys.S,digits=2) ≈  20590.94

        @test calc_sys.phases[1].name == "Cpx"
        @test round(calc_sys.phases[1].molarmass,digits=3) ≈  226.007
        @test round(calc_sys.phases[1].G,sigdigits=5) ≈ -3.0754e6
        close_meemum!(h2o_23SD20A)

        #Testing meemum on files used in Migmatites.jl tests
        hostlib = init_meemum("23SD20A_melt-test1/Host")
        host = minimizepoint(hostlib,800,9000,μ1 = -316240)
    
        hosth2o = getchemical(host.composition,"H2O")
    
        @test round(hosth2o.mol,digits=3) ≈ 0.414
        @test round(hosth2o.μ,sigdigits = 6) ≈ -316240
        @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15
        close_meemum!(hostlib)

        sourcelib = init_meemum("23SD20A_melt-test1/MeltSource")
        source = minimizepoint(sourcelib,875,10000)
        sourcemelt = getphase(source,"melt")[1]
        sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
        @test round(sourcemelth2o.mol,digits=5) ≈ 0.40999
        @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 6.89
        
        melt = minimizepoint(sourcelib,800,9000,composition=sourcemelt.composition.*100)
        melth2o = getchemical(melt.composition,"H2O")
        @test round(melth2o.μ,sigdigits=6) ≈ -316239
        @test round(melth2o.mol,digits=3) ≈ 40.999

        close_meemum!(sourcelib)

        meltlib = init_meemum("23SD20A_melt-test1/Melt")
        melt = minimizepoint(meltlib,800,9000,composition=sourcemelt.composition)
        melth2o = getchemical(melt.composition,"H2O")
        @test round(melth2o.μ,sigdigits=6) ≈ -316240
        @test round(melth2o.mol,digits=5) ≈ 0.40999
        close_meemum!(meltlib)

        hostlib = init_meemum("23SD20A_melt-test1/Host")
        host = minimizepoint(hostlib,800,9000,μ1 = -316240)
    
        hosth2o = getchemical(host.composition,"H2O")
    
        @test round(hosth2o.mol,digits=3) ≈ 0.414
        @test round(hosth2o.μ,sigdigits = 6) ≈ -316240
        @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15
        close_meemum!(hostlib)

        
    end

    @testset "Pseudosection test" begin

        include("PlotDefaults.jl")
       
        set_theme!(myTheme)
        fig = Figure(size = (600,450))
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

    @testset "Werami test" begin
        function feldspar_conditions(phase)
            if lowercase(phase.name) == "fsp"
                k2o = getchemical(phase.composition,"K2O")
                na2o = getchemical(phase.composition,"Na2O")
                cao = getchemical(phase.composition,"CaO")
            
                if mol(k2o*2) >0.1
                    return changename(phase,"Afs")
                else
                    return changename(phase,"Pl")
                end
            else
                return phase
            end
        end

        wdf, xvar, yvar = read_werami_output("20SD06_R2c_1_trimmed.phm", iscelsius = true, iskbar = true)

        petrodf = werami_to_petrosys(wdf,xvar, yvar, phasefunc = [feldspar_conditions])

        petrodfrows = eachrow(petrodf)

        testrow = petrodfrows[1]
        @test testrow[1] == 723 - 273.15
        @test testrow[2] == 2000.0/1000.0
        @test length(petrodfrows) == 180
        psys = testrow[3]
        @test length(psys.phases) == 9
        
        @test round(get_volprop(psys,"Cpx")*100,digits=2) ≈ 3.06 

        
    end

    
    @testset "ModeBox test" begin

        #Parameter functions for phase imports

        function feldspar_conditions(phase)
            if lowercase(phase.name) == "fsp"
                k2o = getchemical(phase.composition,"K2O")
                na2o = getchemical(phase.composition,"Na2O")
                cao = getchemical(phase.composition,"CaO")
            
                if mol(k2o*2) >0.1
                    return changename(phase,"Afs")
                else
                    return changename(phase,"Pl")
                end
            else
                return phase
            end
        end

        sourcelib = init_meemum("23SD20A_melt-test1/MeltSource")
        trange = range(800,900,length=10)
        sources = [minimizepoint(sourcelib,t,10000,phasefunc = [feldspar_conditions]) for t in trange]

        fig = Figure(size = (600,450))
        ax = Axis(fig[1,1])
    
        modebox!(ax,trange,sources)
        fig[1,2] = Legend(fig,ax)
        save("23SD20A_melt-test1/Sources.svg",fig)
    end
end
