include("general_plots.jl")
include("model.jl")

using Images
savepath = "../PLOTS/"
global ds = dsd["normal"]
global dsa = dsd["accum"]
global sims = dsa["sim"][:]

function plot_abstract(save=false)
        dsa = dsd["accum"]
        qi = dsa["qi"][:]
        qni = dsa["qni"][:]
        ri = calc_radius_ice(qi,qni) .* 1000
        rhi = dsa["RHi"][:]
        inps = dsa["INP"][:]
        caps = dsa["CAP"][:]
        
        msk3 = (inps .>= 1) .& (caps .>= 0.8)
        msk4 = (inps .< 0.5) .& (caps .< 0.4)

        rbins = (1e-3:0.01:0.4)
        qnibins = 10 .^ (2:0.1:9)
        rhibins = 80:1:135

        bp1 = blankp_clean()
        plot!(bp1,[],lw=1.5,c=:red,label="INP scaling and CAP large")
        plot!(bp1,[],lw=1.5,c=:cyan,label="INP scaling and CAP small")
        plot!(bp1,legend=:top,legend_columns=2)
        
        #p4 = plot(qi[msk3] .* 1000,st=:stephist,
        #          bins=qibins,c=:red,lw=1.5,
        #          xlabel=L"$q_i\quad[\mathrm{g/kg}]$",ylabel="counts")
        #plot!(p4,qi[msk4] .* 1000,st=:stephist,bins=qibins,c=:cyan,lw=1.5,label="CAP < 0.4")
        #plot!(p4,title="d)",title_location=:left,legend=:false,
        #      xscale=:log10,xticks=[1e-4,1e-3,1e-2,1e-1,1e0])

        p1 = plot(qni[msk3],st=:stephist, bins=qnibins, xscale=:log10, label=L"CAP $\geq$ 0.8",
                  c=:red,lw=1.5,
                  xlabel=L"$N_i\quad[\mathrm{1/kg}]$",ylabel="counts")
        plot!(p1,qni[msk4],st=:stephist, bins=qnibins,c=:cyan,lw=1.5)
        plot!(p1,title="Ice number @ end of ascent",title_location=:left,legend=:false,grid=false)
        
        p2 = plot(ri[msk3],st=:stephist,bins=rbins,c=:red,lw=1.5,
                  xlabel=L"$r_i\quad[\mathrm{mm}]$",ylabel="counts")
        plot!(p2,ri[msk4],st=:stephist,bins=rbins,c=:cyan,lw=1.5)
        plot!(p2,legend_columns=1,background_color_legend=nothing,
              title="Ice radius @ end of ascent",title_location=:left,legend=:false,
              grid=false)
        
        p3 = plot(rhi[msk3],st=:stephist,bins=rhibins,c=:red,lw=1.5,
                  xlabel=L"$RH_i\quad[\%]$",ylabel="counts")
        plot!(p3,rhi[msk4],st=:stephist,bins=rhibins,c=:cyan,lw=1.5)
        plot!(p3,legend_columns=1,background_color_legend=nothing,
              title="RHi @ end of ascent",title_location=:left,legend=:false,
              grid=false)

        lay = @layout [o{0.07h}
                       d e f]

        p_out = plot(bp1,p1,p2,p3,layout=lay,
                     titlefontsize=9,size=(800,350),dpi=:300,
                     bottom_margin=5mm,left_margin=5mm)

	if save
		savefig(p_out,savepath*"Abstract.png")
	else
		return p_out
	end
end

#----------------------

function get_average_hist(vr::String,bins)
        yvals = zeros((length(bins) .* 2) -1)
        xvals = zeros(length(yvals))
        miny = zeros(length(yvals))
        maxy = zeros(length(yvals))

        for i in 1:70
            p0 = plot(dsa[vr][sims .== i],st=:stephist,bins=bins)
            if i == 1
                xvals[:] = p0[1][1][:x][1:length(yvals)]
                miny = p0[1][1][:y][1:length(yvals)]
                maxy = p0[1][1][:y][1:length(yvals)]
            end
            ytemp = p0[1][1][:y][1:length(yvals)]
            yvals = yvals .+ ytemp
            for j in 1:length(yvals)
                if miny[j] > ytemp[j]
                    miny[j] = ytemp[j]
                end
                if maxy[j] < ytemp[j]
                    maxy[j] = ytemp[j]
                end
            end
        end
        return (xvals,yvals ./ 70,miny,maxy)
end

function plot_qv_end(save=false)
        sim1 = 65 - 19
        sim2 = 58 - 19
        qvmax = dsa["qv"][:][sims .== sim1] .* 1000
        qvmin = dsa["qv"][:][sims .== sim2] .* 1000
        bins = 0:0.01:1
        xvals,yvals,ymin,ymax = get_average_hist("qv",bins ./ 1000)

        p1 = plot(qvmax,st=:stephist,label="largest mean",c=:blue,lw=2,xlims=(0,1),
                  bins=bins)
	plot!(p1,qvmin,st=:stephist,label="smallest mean",c=:coral,lw=2,bins=bins)
        plot!(p1,xvals .* 1000,yvals,label="average",c=:grey,lw=1.5,
              fillrange=(ymin,ymax),fillalpha=:0.5)
	plot!(xlabel="qv at end of ascent [g/kg]",ylabel="Counts")
	return p1
end

function plot_t_end(save=false)
        sim1 = 80 - 19
        sim2 = 81 - 19
        tmax = dsa["t"][:][sims .== sim1] .- 273.15
        tmin = dsa["t"][:][sims .== sim2] .- 273.15
        bins = -70:0.5:-20
        xvals,yvals,ymin,ymax = get_average_hist("t",bins .+ 273.15)
        
        p1 = plot(tmax,st=:stephist,label="largest mean",c=:blue,lw=1.5,bins=bins)
        plot!(p1,tmin,st=:stephist,label="smallest mean",c=:coral,lw=1.5,bins=bins)
        plot!(p1,xvals .- 273.15,yvals,label="average",c=:grey,lw=1.5,
             fillrange=(ymin,ymax),fillalpha=:0.5)
	plot!(xlabel="T at end of ascent [°C]",ylabel="Counts")
	return p1
end

function plot_p_end(save=false)
        sim1 = 65 - 19
        sim2 = 48 - 19
        pmax = dsa["p"][:][sims .== sim1]
        pmin = dsa["p"][:][sims .== sim2]
        bins = 180:3:420
        xvals,yvals,ymin,ymax = get_average_hist("p",bins)

	p1 = plot(pmax,st=:stephist,label="largest mean",c=:blue,lw=1.5,bins=bins)
	plot!(p1,pmin,st=:stephist,label="smallest mean",c=:coral,lw=1.5,bins=bins)
        plot!(p1,xvals,yvals,label="average",c=:grey,lw=1.5,
             fillrange=(ymin,ymax),fillalpha=:0.5)
	plot!(xlabel="p at end of ascent [hPa]",ylabel="Counts")
	return p1
end

function plot_qi_end(save=false)
        sim1 = 33
        sim2 = 46
        qimax = dsa["qi"][:][sims .== sim1]
        qimin = dsa["qi"][:][sims .== sim2]
        #bins = 0:0.005:0.4
        bins = 10 .^ (-3.5:0.05:0)
        xvals,yvals,ymin,ymax = get_average_hist("qi",bins ./ 1000)

	p1 = plot(qimax .* 1000,st=:stephist,label="largest mean",c=:blue,lw=1.5,bins=bins,xscale=:log10)
	plot!(p1,qimin .* 1000,st=:stephist,label="smallest mean",c=:coral,lw=1.5,bins=bins)
        plot!(p1,xvals .* 1000,yvals,label="average",c=:grey,lw=1.5,
             fillrange=(ymin,ymax),fillalpha=:0.5)
        plot!(xlabel="qi at end of ascent [g/kg]",ylabel="Counts",xticks=[1e-3,1e-2,1e-1,1e0])
	return p1
end

function plot_RHi_end(save=false)
        i1 = 70 + 19
        i2 = 61 + 19
        sim1 = 70
        sim2 = 61
        RHimax = dsa["RHi"][:][sims .== sim1]
        RHimin = dsa["RHi"][:][sims .== sim2]
        bins = 85:0.75:140
        xvals,yvals,ymin,ymax = get_average_hist("RHi",bins)

        p1 = plot(RHimax,st=:stephist,label="largest mean",c=:blue,lw=1.5,bins=bins)
        plot!(p1,RHimin,st=:stephist,label="smalles mean",c=:coral,lw=1.5,bins=bins)
        plot!(p1,xlabel="RHi at end of ascent [%]",ylabel="Counts",title="a)",
              title_location=:left,xlims=(85,140))
        plot!(p1,xvals,yvals,label="average",c=:grey,lw=1.5,
             fillrange=(ymin,ymax),fillalpha=:0.5)
        return p1
end

function plot_lat_end(save=false)
        i1 = 70 + 19
        i2 = 29 + 19
        sim1 = 70
        sim2 = 29
        latmax = dsa["lat"][:][sims .== sim1]
        latmin = dsa["lat"][:][sims .== sim2]
        bins = 30:0.5:85
        xvals,yvals,ymin,ymax = get_average_hist("lat",bins)

        p1 = plot(latmax,st=:stephist,label="largest mean",c=:blue,lw=1.5,bins=bins)
        plot!(p1,latmin,st=:stephist,label="smalles mean",c=:coral,lw=1.5,bins=bins)
        plot!(p1,xlabel="lat at end of ascent [°N]",ylabel="Counts",title="a)",
              title_location=:left,xlims=(30,85))
        plot!(p1,xvals,yvals,label="average",c=:grey,lw=1.5,
             fillrange=(ymin,ymax),fillalpha=:0.5)
        return p1
end

function plot_lon_end(save=false)
        i1 = 48 + 19
        i2 = 36 + 19
        sim1 = 48
        sim2 = 36
        lonmax = dsa["lon"][:][sims .== sim1]
        lonmin = dsa["lon"][:][sims .== sim2]
        bins = -65:1:45
        xvals,yvals,ymin,ymax = get_average_hist("lon",bins)

        p1 = plot(lonmax,st=:stephist,label="largest mean",c=:blue,lw=1.5,bins=bins)
        plot!(p1,lonmin,st=:stephist,label="smalles mean",c=:coral,lw=1.5,bins=bins)
        plot!(p1,xlabel="lon at end of ascent [°E]",ylabel="Counts",title="c)",
              title_location=:left,xlims=(-65,45))
        plot!(p1,xvals,yvals,label="average",c=:grey,lw=1.5,
             fillrange=(ymin,ymax),fillalpha=:0.5)
        plot!(xticks=collect(-60:10:40))
        return p1
end

function plot_geographical(save=false)
        
        bp1 = blankp_clean()
	plot!(bp1,[],c=:blue,lw=1.5,label="largest mean",legend=:top)
	plot!(bp1,[],c=:coral,lw=1.5,label="smallest mean",legend_columns=4)
        plot!(bp1,[],c=:grey,lw=1.5,label="average")
        plot!(bp1,[],fillrange=([NaN],[NaN]),lw=0,label="min/max range",c=:grey,fillalpha=:0.5)

        p1 = plot_lon_end()
        plot!(p1,legend=false,title="a)",bottom_margin=3mm)

        p2 = plot_lat_end()
        plot!(p2,legend=false,title="b)",bottom_margin=3mm)

        bg = load("../PLOTS/lonlat.png")
        p3 = plot(bg,xlims=(0,2400),ylims=(0,1050),xaxis=false,yaxis=false,xticks=false,yticks=false,
                  title="c)",title_location=:left,top_margin=-1mm)
        lay = @layout [o{0.06h}
                       a b
                       c{0.55h}]
        p_out = plot(bp1,p1,p2,p3,layout=lay,size=(1000,700),dpi=:400,right_margin=5mm,
                     left_margin=5mm,titlefontsize=11,legendfontsize=10)
        if save
            savefig(p_out,savepath*"geographical.png")
        else
            return p_out
        end
end


function plot_histograms(save=false)
        ds = dsd["normal"]
	
        p1 = plot_t_end()
        plot!(p1,legend=false,title="a)",title_location=:left,top_margin=-2mm,bottom_margin=-2mm)
        p2 = plot_spread_vr("t","normal",1)
        plot!(p2,title="b)",title_location=:left,ylabel="T [°C]",top_margin=-2mm,bottom_margin=-2mm)
        
        p3 = plot_p_end()
	plot!(p3,legend=false,title="c)",title_location=:left,top_margin=-2mm,bottom_margin=-2mm)
        p4 = plot_spread_vr("p","normal",1)
        plot!(p4,title="d)",title_location=:left,ylabel="p [hPa]",top_margin=-2mm,bottom_margin=-2mm)
        
        p5 = plot_qv_end()
	plot!(p5,legend=false,title="e)",title_location=:left,top_margin=-2mm,bottom_margin=-2mm)
        p6 = plot_spread_vr("qv","normal",1000)
        plot!(p6,title="f)",title_location=:left,ylabel="qv [g/kg]",top_margin=-2mm,bottom_margin=-2mm)
       
        p7 = plot_qi_end()
        plot!(p7,legend=false,title="g)",title_location=:left,top_margin=-2mm,bottom_margin=-2mm)
        p8 = plot_spread_vr("qi","normal",1000)
        plot!(p8,title="h)",title_location=:left,ylabel="qi [g/kg]",top_margin=-2mm,bottom_margin=-2mm)
        
        p9 = plot_RHi_end()
	plot!(p9,legend=false,title="i)",title_location=:left,top_margin=-2mm,bottom_margin=-2mm)
        p10 = plot_spread_vr("RHi","normal",1)
        plot!(p10,title="j)",title_location=:left,ylabel="RHi [%]",top_margin=-2mm,bottom_margin=-2mm)
       
        bp1 = blankp_clean()
	plot!(bp1,[],c=:blue,lw=1.5,label="largest mean",legend=:top)
	plot!(bp1,[],c=:coral,lw=1.5,label="smallest mean",legend_columns=2)
        plot!(bp1,[],c=:grey,lw=1.5,label="average")
        plot!(bp1,[],fillrange=([NaN],[NaN]),lw=0,label="min/max range",c=:grey,fillalpha=:0.5)

        bp2 = blankp_clean()
	plot!(bp2,[],c=:black,lw=1.5,label="mean",legend=:top)
	plot!(bp2,[],c=:black,lw=1.0,label="median",legend=:top,ls=:dash)
        plot!(bp2,[],fillrange=([NaN],[NaN]),c=:black,
                                fillalpha=0.4,lw=0,label="25th-75th")
        plot!(bp2,[],fillrange=([NaN],[NaN]),c=:black,
                                fillalpha=0.2,lw=0,label="5th-95th",legend_columns=2)
	
        lay = @layout [o{0.05h} u{0.05h}
		        a b
			c d
                        e f
                        g h
                        i j]
	
        p_out = plot(bp1,bp2,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,layout=lay,dpi=:350,size=(800,800),
                     labelfontsize=:9,titlefontsize=:9,right_margin=5mm,left_margin=5mm)
 	
        if save
		savefig(p_out,savepath*"histograms.png")
	else
		return p_out
	end
end
#----------------------


#----------------------
function plot_qv_T_correlations(save=false,typ="normal")
	p1a = plot_cor_vr(["qv"],1,typ)
	plot!(p1a,colorbar=false,titlefontsize=12,labelfontsize=10,
	       xticks=false,yticks=([1],["mean qv"]),
	       title="c)",title_location=:left,top_margin=2mm)
	p1b = plot_cor_vr(["t"],1,typ)
	plot!(p1b,colorbar=false,titlefontsize=1,labelfontsize=10,title="",
	       xticks=false,yticks=([1],["mean T"]))
	p1c = plot_cor_vr(["qv"],7,typ)
	plot!(p1c,colorbar=false,titlefontsize=1,labelfontsize=10,title="",
	       xticks=false,yticks=([1],["95th qv"]))
	p1d = plot_cor_vr(["t"],7,typ)
	plot!(p1d,colorbar=false,titlefontsize=1,labelfontsize=10,title="",
	       yticks=([1],["95th T"]))


	cbar = spear_cbar()

	lay1 = @layout [a
			b
			c
			d
			e{0.15h}]

	p_out1 = plot(p1a,p1b,p1c,p1d,cbar,layout=lay1)

	ind_maxt = 65
	ind_mint = 25

        ds = dsd[typ]
	p_out2 = scatter(ds["t"][end,:] .- 273.15,ds["qv"][end,:] .* 1000,c=:black,
					 markerstrokewidth=0,markersize=5,label=false,
					 ylabel="qv 95th percentile [g/kg]",
					 xlabel="T 95th percentile [°C]",
					 title="a)",title_location=:left,
					 titlefontsize=12)
	
	scatter!(p_out2,[ds["t"][end,ind_maxt] - 273.15],[ds["qv"][end,ind_maxt] * 1000],
		  c=:cyan,label="max qv",markerstrokewidth=0)	
	scatter!(p_out2,[ds["t"][end,ind_mint] - 273.15],[ds["qv"][end,ind_mint] * 1000],
		  c=:coral,label="min qv",markerstrokewidth=0)	
	
	tmax = load_vr_end("t",ind_maxt)
	qvmax = load_vr_end("qv",ind_maxt)
	pmax = load_vr_end("p",ind_maxt)
	qvs_max = qv_sat_i.(tmax,pmax .* 100) .* 1000

	tmin = load_vr_end("t",ind_mint)
	qvmin = load_vr_end("qv",ind_mint)
	pmin = load_vr_end("p",ind_mint)
	qvs_min = qv_sat_i.(tmin,pmin .* 100) .* 1000	

	p_out3 = scatter(qvmax .* 1000,qvs_max,c=:cyan,
                                         markerstrokewidth=0,markersize=1,label="max qv",
                                         xlabel="qv [g/kg]",
                                         ylabel="qv_sat(T,p) [g/kg]",
					 title="b)",title_location=:left,
					 titlefontsize=12)

        scatter!(p_out3,qvmin .* 1000,qvs_min,c=:coral,
			      markerstrokewidth=0,markersize=1,label="min qv")


	lay = @layout [a b
			 c{0.35h}]
	p_out = plot(p_out2,p_out3,p_out1,layout=lay,size=(800,700),dpi=:300)
 	if save
		savefig(p_out,savepath*"Fig2.png")
	else
		return p_out
	end
end
#----------------------

#----------------------
function plot_CCN_figure(save=false,typ="normal")
	p1 = scat_vr_param_shaded("CCN","p_gl50",1,1,"hPa",typ,false,false)
	plot!(p1,title="a)",title_location=:left,titlefontsize=10,xlabel="CCN scaling factor",
              ylabel="p(50% gl.) [hPa]")

	p2 = scat_vr_param_shaded("CCN","T_gl50",1,1,"°C",typ,false,false)
	plot!(p2,title="b)",title_location=:left,titlefontsize=10,xlabel="CCN scaling factor",
              ylabel="T(50% gl.) [°C]")

	p3 = scat_vr_param_shaded("CCN","max_SLW",1,1000,"g/kg",typ,false,false)
	plot!(p3,ylabel="max SLW [g/kg]", title="c)",title_location=:left,titlefontsize=10,
              xlabel="CCN scaling factor")
	if typ == "glac"
            typ2 = "normal"
        elseif typ == "conv"
            typ2 = "conv"
        elseif typ == "normal"
            typ2 = "normal"
        end

	p4 = scat_vr_param_shaded("CCN","maxnH",1,1e-6,"1e6/kg",typ2,false,false)
	plot!(p4,ylabel="max Nhy [1e6/kg]", title="d)",title_location=:left,titlefontsize=10,
             xlabel="CCN scaling factor")

	p5 = plot_cor_param_vr(["CCN"],
			       ["max_qc","max_qs","max_qi","max_qg","max_qr"],
			       1,typ2)
        plot!(p5,title="e)",titlefontsize=10,title_location=:left,colorbar=false)#,
             #xticks=["max_qc","max_qs","max_qi","max_qg","max_qr"])
	#p6 = plot_cor_param_vr(["CAP"],
	#		       ["max_qc","max_qs","max_qi","max_qg","max_qr"],
	#		       1,typ2)
	#plot!(p6,colorbar=false,top_margin=-2mm)

	cbar = spear_cbar()
	blankp = blankp_5_95()
	
	lay = @layout [a b
                        c d
                        e{0.04h}
                        g{0.04h}]

	p_out = plot(p1,p2,p3,p4,p5,cbar,layout=lay,size=(800,650),labelfontsize=9,dpi=:350)
	if save
		savefig(p_out,savepath*"Fig3.png")
	else
		return p_out
	end
end
#----------------------

#----------------------

function plot_histograms_CAP(msk,save=false)
        dsa = dsd["accum"]
        qi = dsa["qi"][:]
        qni = dsa["qni"][:]
        ri = calc_radius_ice(qi,qni) .* 1000
        caps = dsa["CAP"][:]
        
        msk1 = (caps .>= 0.8) .& msk
        msk2 = (caps .< 0.4) .& msk

        rbins = (1e-3:0.01:0.4)
        qnibins = 10 .^ (2:0.1:9)
        taubins = 10 .^ (1:0.1:6)
        qibins = 10 .^ (-4:0.05:0.5)

        bp1 = blankp_clean()
        plot!(bp1,[],lw=1.5,c=:red,label=L"CAP $\geq$ 0.8")
        plot!(bp1,[],lw=1.5,c=:cyan,label=L"CAP $<$ 0.4")
        plot!(bp1,legend=:top,legend_columns=2)

	p1 = scat_vr_param_shaded("CAP","qi",2,1000,"g/kg","normal",false,false,:log10,(6e-4,3e-1))
	plot!(p1,title="a)",title_location=:left,ylabel=L"q_\mathrm{i}\quad [\mathrm{g/kg}]",
	       xlabel=L"\mathrm{CAP}",)	
        
        p2 = scat_vr_param_shaded("CAP","qni",2,1,"1/kg","normal",false,false,:log10,(1e3,5e7),true)
	plot!(p2,title="b)",title_location=:left,ylabel=L"N_\mathrm{i}\quad [1/\mathrm{kg}]",
               xlabel=L"\mathrm{CAP}",yticks=[1e3,1e4,1e5,1e6,1e7])
        
        p3 = scat_vr_param_shaded("CAP","r",2,1000,"mm","normal",false,false)
	plot!(p3,title="c)",title_location=:left,ylabel=L"r_\mathrm{i}\quad [\mathrm{mm}]",
	       xlabel=L"\mathrm{CAP}")

        p4 = plot(qi[msk1] .* 1000,st=:stephist,
                  bins=qibins,c=:red,lw=1.5,label=L"CAP $\geq$ 0.8",
                  xlabel=L"$q_i\quad[\mathrm{g/kg}]$",ylabel="counts")
        plot!(p4,qi[msk2] .* 1000,st=:stephist,bins=qibins,c=:cyan,lw=1.5,label="CAP < 0.4")
        plot!(p4,title="d)",title_location=:left,legend=:false,
              xscale=:log10,xticks=[1e-4,1e-3,1e-2,1e-1,1e0])

        p5 = plot(qni[msk1],st=:stephist, bins=qnibins, xscale=:log10, label=L"CAP $\geq$ 0.8",
                  c=:red,lw=1.5,
                  xlabel=L"$N_i\quad[\mathrm{1/kg}]$",ylabel="counts")
        plot!(p5,qni[msk2],st=:stephist, bins=qnibins, label="CAP < 0.4",c=:cyan,lw=1.5)
        plot!(p5,title="e)",title_location=:left,legend=:false)
        
        p6 = plot(ri[msk1],st=:stephist,bins=rbins,label=L"CAP $\geq$ 0.8",c=:red,lw=1.5,
                  xlabel=L"$r_i\quad[\mathrm{mm}]$",ylabel="counts")
        plot!(p6,ri[msk2],st=:stephist,bins=rbins,label="CAP < 0.4",c=:cyan,lw=1.5)
        plot!(p6,legend_columns=1,background_color_legend=nothing,
              title="f)",title_location=:left,legend=:false)


        lay = @layout [a b c
                       o{0.05h}
                       d e f]

        p_out = plot(p1,p2,p3,bp1,p4,p5,p6,layout=lay,
                     titlefontsize=10,size=(800,600),dpi=:300)

	if save
		savefig(p_out,savepath*"Fig4.png")
	else
		return p_out
	end
end
#----------------------


#----------------------
function plot_histograms_INP(msk,save=false)
        dsa = dsd["accum"]
        qi = dsa["qi"][:]
        qni = dsa["qni"][:]
        ri = calc_radius_ice(qi,qni) .* 1000
        inps = dsa["INP"][:]
        
        msk3 = (inps .>= 1) .& msk
        msk4 = (inps .< 0.5) .& msk

        rbins = (1e-3:0.01:0.4)
        qnibins = 10 .^ (2:0.1:9)
        taubins = 10 .^ (1:0.1:6)
        qibins = 10 .^ (-4:0.05:0.5)

        bp1 = blankp_clean()
        plot!(bp1,[],lw=1.5,c=:red,label=L"INP scaling $\geq$ 1")
        plot!(bp1,[],lw=1.5,c=:cyan,label=L"INP scaling $<$ 0.5")
        plot!(bp1,legend=:top,legend_columns=2)
        
        p1 = scat_vr_param_shaded("INP","qi",2,1000,"g/kg","normal",false,false,:log10,(6e-4,3e-1))
	plot!(p1,title="a)",title_location=:left,ylabel=L"q_\mathrm{i}\quad [\mathrm{g/kg}]",
	       xlabel=L"\mathrm{INP\,\,scaling}")
        plot!(p1,xscale=:log10,xlims=(0.01,20),xticks=[1e-2,1e-1,1,1e1])
        
        p2 = scat_vr_param_shaded("INP","qni",2,1,"1/kg","normal",false,false,:log10,(1e3,5e7),true)
	plot!(p2,title="b)",title_location=:left,ylabel=L"N_\mathrm{i}\quad [1/\mathrm{kg}]",
               xlabel=L"\mathrm{INP\,\,scaling}",yticks=[1e3,1e4,1e5,1e6,1e7])
        plot!(p2,xscale=:log10,xlims=(0.01,20),xticks=[1e-2,1e-1,1,1e1])
	
        p3 = scat_vr_param_shaded("INP","r",2,1000,"mm","normal",false,false)
	plot!(p3,title="c)",title_location=:left,ylabel=L"r_\mathrm{i}\quad [\mathrm{mm}]",
	       xlabel=L"\mathrm{INP\,\,scaling}")
        plot!(p3,xscale=:log10,xlims=(0.01,20),xticks=[1e-2,1e-1,1,1e1])

        p4 = plot(qi[msk3] .* 1000,st=:stephist,
                  bins=qibins,c=:red,lw=1.5,
                  xlabel=L"$q_i\quad[\mathrm{g/kg}]$",ylabel="counts")
        plot!(p4,qi[msk4] .* 1000,st=:stephist,bins=qibins,c=:cyan,lw=1.5,label="CAP < 0.4")
        plot!(p4,title="d)",title_location=:left,legend=:false,
              xscale=:log10,xticks=[1e-4,1e-3,1e-2,1e-1,1e0])

        p5 = plot(qni[msk3],st=:stephist, bins=qnibins, xscale=:log10, label=L"CAP $\geq$ 0.8",
                  c=:red,lw=1.5,
                  xlabel=L"$N_i\quad[\mathrm{1/kg}]$",ylabel="counts")
        plot!(p5,qni[msk4],st=:stephist, bins=qnibins,c=:cyan,lw=1.5)
        plot!(p5,title="e)",title_location=:left,legend=:false)
        
        p6 = plot(ri[msk3],st=:stephist,bins=rbins,c=:red,lw=1.5,
                  xlabel=L"$r_i\quad[\mathrm{mm}]$",ylabel="counts")
        plot!(p6,ri[msk4],st=:stephist,bins=rbins,c=:cyan,lw=1.5)
        plot!(p6,legend_columns=1,background_color_legend=nothing,
              title="f)",title_location=:left,legend=:false)


        lay = @layout [a b c
                       o{0.05h}
                       d e f]

        p_out = plot(p1,p2,p3,bp1,p4,p5,p6,layout=lay,
                     titlefontsize=10,size=(800,600),dpi=:300)

	if save
		savefig(p_out,savepath*"Fig5.png")
	else
		return p_out
	end
end

#not used
function plot_only_histograms_INP_CAP(msk,save=false)
        dsa = dsd["accum"]
        qi = dsa["qi"][:]
        qni = dsa["qni"][:]
        ri = calc_radius_ice(qi,qni) .* 1000
        inps = dsa["INP"][:]
        caps = dsa["CAP"][:]
        
        msk1 = (caps .>= 0.8) .& msk
        msk2 = (caps .< 0.4) .& msk 

        msk3 = (inps .>= 1) .& msk
        msk4 = (inps .< 0.5) .& msk

        rbins = (1e-3:0.01:0.4)
        qnibins = 10 .^ (2:0.1:9)
        taubins = 10 .^ (1:0.1:6)
        qibins = 10 .^ (-4:0.05:0.5)

        bp1 = blankp_clean()
        plot!(bp1,[],lw=1.5,c=:red,label=L"CAP $\geq$ 0.8")
        plot!(bp1,[],lw=1.5,c=:cyan,label=L"CAP $<$ 0.4")
        plot!(bp1,legend=:top,legend_columns=2)

        p1 = plot(qi[msk1] .* 1000,st=:stephist,
                  bins=qibins,c=:red,lw=1.5,label=L"CAP $\geq$ 0.8",
                  xlabel=L"$q_i\quad[\mathrm{g/kg}]$",ylabel="counts")
        plot!(p1,qi[msk2] .* 1000,st=:stephist,bins=qibins,c=:cyan,lw=1.5,label="CAP < 0.4")
        plot!(p1,title="a)",title_location=:left,legend=:false,
              xscale=:log10,xticks=[1e-4,1e-3,1e-2,1e-1,1e0])

        p2 = plot(qni[msk1],st=:stephist, bins=qnibins, xscale=:log10, label=L"CAP $\geq$ 0.8",
                  c=:red,lw=1.5,
                  xlabel=L"$N_i\quad[\mathrm{1/kg}]$",ylabel="counts")
        plot!(p2,qni[msk2],st=:stephist, bins=qnibins, label="CAP < 0.4",c=:cyan,lw=1.5)
        plot!(p2,title="b)",title_location=:left,legend=:false)
        
        p3 = plot(ri[msk1],st=:stephist,bins=rbins,label=L"CAP $\geq$ 0.8",c=:red,lw=1.5,
                  xlabel=L"$r_i\quad[\mathrm{mm}]$",ylabel="counts")
        plot!(p3,ri[msk2],st=:stephist,bins=rbins,label="CAP < 0.4",c=:cyan,lw=1.5)
        plot!(p3,legend_columns=1,background_color_legend=nothing,
              title="c)",title_location=:left,legend=:false)

        bp2 = blankp_clean()
        plot!(bp2,[],lw=1.5,c=:red,label=L"INP scaling $\geq$ 1")
        plot!(bp2,[],lw=1.5,c=:cyan,label=L"INP scaling $<$ 0.5")
        plot!(bp2,legend=:top,legend_columns=2)
        
        p4 = plot(qi[msk3] .* 1000,st=:stephist,
                  bins=qibins,c=:red,lw=1.5,
                  xlabel=L"$q_i\quad[\mathrm{g/kg}]$",ylabel="counts")
        plot!(p4,qi[msk4] .* 1000,st=:stephist,bins=qibins,c=:cyan,lw=1.5,label="CAP < 0.4")
        plot!(p4,title="d)",title_location=:left,legend=:false,
              xscale=:log10,xticks=[1e-4,1e-3,1e-2,1e-1,1e0])

        p5 = plot(qni[msk3],st=:stephist, bins=qnibins, xscale=:log10, label=L"CAP $\geq$ 0.8",
                  c=:red,lw=1.5,
                  xlabel=L"$N_i\quad[\mathrm{1/kg}]$",ylabel="counts")
        plot!(p5,qni[msk4],st=:stephist, bins=qnibins,c=:cyan,lw=1.5)
        plot!(p5,title="e)",title_location=:left,legend=:false)
        
        p6 = plot(ri[msk3],st=:stephist,bins=rbins,c=:red,lw=1.5,
                  xlabel=L"$r_i\quad[\mathrm{mm}]$",ylabel="counts")
        plot!(p6,ri[msk4],st=:stephist,bins=rbins,c=:cyan,lw=1.5)
        plot!(p6,legend_columns=1,background_color_legend=nothing,
              title="f)",title_location=:left,legend=:false)


        lay = @layout [o{0.05h}
                       a b c
                       o{0.05h}
                       d e f]

        p_out = plot(bp1,p1,p2,p3,bp2,p4,p5,p6,layout=lay,
                     titlefontsize=10,size=(800,600),dpi=:300)

	if save
		savefig(p_out,savepath*"SI_1.png")
	else
		return p_out
	end
end
#----------------------


#----------------------
#not used
function plot_qc_Nc_vs_INP(save=false)
		
	p1 = scatter(dp["INP"],ds["qc"][1,:] .* 1000,xlabel="INP scaling",ylabel="mean(qc) [g/kg]",
					     label=false,c=:grey,
					     left_margin=2mm,bottom_margin=2mm,
					     title="a)",title_location=:left,
                                             xscale=:log10,xticks=[1e-2,1e-1,1e0,1e1],xlims=(0.01,20))
	p2 = scatter(dp["INP"],ds["qnc"][1,:] .* 1e-6,xlabel="INP scaling",ylabel="mean(Nc) [1e6/kg]",
					     label=false,c=:grey,
					     left_margin=2mm,bottom_margin=2mm,
					     title="b)",title_location=:left,
                                             xscale=:log10,xticks=[1e-2,1e-1,1e0,1e1],xlims=(0.01,20))
	
	lay = @layout [a b]
	
	p_out = plot(p1,p2,layout=lay,size=(800,400),titlefontsize=10)
	if save
		savefig(p_out,savepath*"Fig6.png")
	else
		return p_out
	end
end
#----------------------


#----------------------
function plot_RHi_2d_hists(save=false,typ="normal")
        ds = dsd[typ]
	X = Matrix([dp["CAP"] dp["INP"] dp["CCN"]])
	Y = ds["RHi"][1,:]

	clims = (104,107)
        p1 = pdp_2d_heatmap(X,Y,1,2,"CAP","INP","\n RHi")
	plot!(p1,title="a)",title_location=:left,clims=clims,colorbar=:false,right_margin=-4mm,
              yscale=:log10,ylims=(0.01,20),yticks=[1e-2,1e-1,1,10],ylabel="INP scaling",left_margin=3mm)
	p2 = pdp_2d_heatmap(X,Y,1,3,"CAP","CCN scaling","\n RHi")
	plot!(p2,title="b)",title_location=:left,clims=clims,colorbar=:false,right_margin=-4mm)
	p3 = pdp_2d_heatmap(X,Y,2,3,"INP scaling","CCN scaling","\n RHi")
	plot!(p3,title="c)",title_location=:left,clims=clims,colorbar=:false,
              xscale=:log10,xlims=(0.01,20),xticks=[1e-2,1e-1,1,10],top_margin=0mm)	

	cx = (104:0.01:107)
	cx2 = 104:0.5:107
	cbar = heatmap(cx,cx,cx' .* ones(length(cx),1), legend=:none, yticks=:none, xticks=(cx2, string.(cx2)))
        plot!(cbar,xlabel="predicted mean(RHi) [%]",bottom_margin=5mm)
	

	lay = @layout [a b c
                       d{0.05h}]
	
	p_out = plot(p1,p2,p3,cbar,layout=lay,dpi=:300,size=(800,330),
		      labelfontsize=9,titlefontsize=10)
	if save
		savefig(p_out,savepath*"Fig8.png")
	else
		return p_out
	end
end
#----------------------

#----------------------
function plot_CAP_INP_CCN_scatters_RHi_tau(save=false,typ="normal")
        ds = dsd[typ]
        dsa = dsd["accum"]

        msk1 = dp["INP"] .>= 1
        msk2 = dp["INP"] .< 1
        
        p_out1 = scat_vr_param_shaded_2msks("CAP","RHi",2,1,msk1,msk2,
					    L"$\mathrm{INP} \geq 1$",L"$\mathrm{INP} < 1$","%",
                                            typ,false)
	plot!(p_out1,title="a)",title_location=:left,titlefontsize=:10,xlabel="CAP")
	plot!(p_out1,legend=:top,bottom_margin=-2mm,legend_columns=2)

	msk3 = dp["CAP"] .>= 0.5
	msk4 = dp["CAP"] .< 0.5
	p_out2 = scat_vr_param_shaded_2msks("INP","RHi",2,1,msk3,msk4,
					    L"$\mathrm{CAP} \geq 0.5$",L"$\mathrm{CAP} < 0.5$",
					    "%",typ,false)
	plot!(p_out2,title="b)",title_location=:left,titlefontsize=:10,xlabel="INP scaling")
        plot!(p_out2,legend=:top,xscale=:log10,xlims=(0.01,20),xticks=[1e-2,1e-1,1,1e1],
              bottom_margin=-2mm,legend_columns=2)


	p_out3 = scat_vr_param_shaded("CCN","RHi",2,1,"%",typ,false,false)
	plot!(p_out3,title="c)",title_location=:left,titlefontsize=:10,xlabel="CCN scaling")
	plot!(p_out3,legend=:false,background_color_legend=nothing,bottom_margin=-2mm,legend_columns=2)
        

        p_out4 = scat_vr_param_shaded_2msks("CAP","tau",2,1,msk1,msk2,
                                            L"\mathrm{INP} \geq 1",
                                            L"\mathrm{INP} < 1","s","normal",false)
        plot!(p_out4,yscale=:log10,ylims=(10,1e5),title="d)",title_location=:left,titlefontsize=:10,
             legend=:top,xlabel="CAP",bottom_margin=0mm,legend_columns=2,
             ylabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$")
        p_out5 = scat_vr_param_shaded_2msks("INP","tau",2,1,msk3,msk4,
                                            L"$\mathrm{CAP} \geq 0.5$",
                                            L"$\mathrm{CAP} < 0.5$","s","normal",false)
        plot!(p_out5,yscale=:log10,ylims=(10,1e5),title="e)",title_location=:left,titlefontsize=:10)
        plot!(p_out5,legend=:top,xscale=:log10,xlims=(0.01,20),xticks=[1e-2,1e-1,1,1e1],
              xlabel="INP scaling",bottom_margin=-0mm,legend_columns=2,
             ylabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$")
        p_out6 = scat_vr_param_shaded("CCN","tau",2,1,"s","normal",false,false)
        plot!(p_out6,yscale=:log10,ylims=(10,1e5),title="f)",title_location=:left,titlefontsize=:10,
              legend=:false,xlabel="CCN scaling",bottom_margin=0mm,legend_columns=2,
             ylabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$")


	lay = @layout [a b c
                       d e f]
	p_out = plot(p_out1,p_out2,p_out3,p_out4,p_out5,p_out6,layout=lay,size=(800,600),dpi=:300)
	if save
		savefig(p_out,savepath*"Fig9.png")
	else
		return p_out
	end
end
#----------------------


#----------------------
function plot_RHi_tau_CAPINP(msk=true,save=false)
        dsa = dsd["accum"]

        msk1 = (dsa["CAP"][:] .>= 0.8) .& (dsa["INP"][:] .>= 1) .& msk
        msk2 = (dsa["CAP"][:] .< 0.4) .& (dsa["INP"][:] .< 0.5) .& msk

        RHbins = 80:0.5:135
        taubins = 10 .^ (1:0.1:6)
        bp = blankp_clean()
        plot!(bp,[],lw=1.5,c=:red,label="CAP > 0.8 and INP > 1")
        plot!(bp,[],lw=1.5,c=:cyan,label="CAP < 0.4 and INP < 0.5")
        plot!(bp,legend=:top)
        
        p1 = plot(dsa["RHi"][msk1],st=:stephist,c=:red,lw=1.5,bins=RHbins,
                      normalize=:true,label=false,ylabel="normalized counts",
                      xlabel=L"$\mathrm{RH}_i\quad[\%]$")
        plot!(p1,dsa["RHi"][msk2],st=:stephist,c=:cyan,lw=1.5,bins=RHbins,
              normalize=:true,label=false,title="a)",title_location=:left)

        p2 = plot(dsa["tau"][msk1],st=:stephist,c=:red,lw=1.5,bins=taubins,
                      normalize=:false,label=false,ylabel="counts",
                      xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",
                      xscale=:log10)
        plot!(p2,dsa["tau"][msk2],st=:stephist,c=:cyan,lw=1.5,bins=taubins,
              normalize=:false,label=false,xticks=[10,1e2,1e3,1e4,1e5,1e6],
              title="b)",title_location=:left)

        lay = @layout [o{0.1h}
                       a b]
        p_out = plot(bp,p1,p2,layout=lay,dpi=:300,titlefontsize=10,size=(800,400),
                    bottom_margin=5mm,left_margin=5mm,right_margin=5mm)

 	if save
		savefig(p_out,savepath*"Fig10.png")
	else
		return p_out
	end
end
#----------------------


#----------------------

function plot_tau_RHi_hist_all(save=false)
        dsa = dsd["accum"] 
        tau_all = dsa["tau"][:]
        RHi_all = dsa["RHi"][:]
        p_out = plot_hist2d_2ars_general(tau_all,RHi_all,10 .^ (1:0.1:5),80:1:140,"tau","RHi")
        plot!(p_out,xscale=:log10,legend=false)
        plot!(xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",ylabel=L"$\mathrm{RH}_i\quad[\%]$",
             dpi=:350)
 	if save
		savefig(p_out,savepath*"Fig11.png")
	else
		return p_out
	end
end
#----------------------

#----------------------
function plot_glac_pres_temp_hist2d(save=false)
        dsa = dsd["accum"]
        tbins = collect(220:1:275) .- 273.15
        pbins = 250:5:600
        t6bins = 0:1:48
        taulabel = L"$\tau_{600}\quad [\mathrm{h}]$"
        plabel = L"$\mathrm{Glaciation\,\,pressure}\quad[\mathrm{hPa}]$"
        tlabel = L"$\mathrm{Glaciation\,\,temperature}\quad[\mathrm{°C}]$"

        p1 = plot_hist2d_2ars_general(dsa["t600"][:],dsa["p_gl"][:],t6bins,pbins,"","")
        plot!(p1,xlabel=taulabel,ylabel=plabel,legend=false,title="a)",
             title_location=:left,colorbar=false,left_margin=5mm,bottom_margin=5mm)
        p2 = plot_hist2d_2ars_general(dsa["t600"][:],dsa["t_gl"][:] .- 273.15,t6bins,tbins,"","")
        plot!(p2,xlabel=taulabel,ylabel=tlabel,legend=false,title="b)",
              title_location=:left,bottom_margin=5mm,lw=0.5)

        lay = @layout [a{0.4w} b{0.6w}]
        p_out = plot(p1,p2,layout=lay,dpi=:300,size=(1000,400))
 	if save
		savefig(p_out,savepath*"Fig12.png")
	else
		return p_out
	end
end
#----------------------

#----------------------
function plot_SAT_effects(save=false)

        p1 = scat_vr_param_shaded("SATAD","p",1,1,"1","conv",false,false)
        plot!(p1,xlabel="SAT perturbation",ylabel="pressure [hPa]",
              title="a)",title_location=:left,left_margin=3mm,bottom_margin=5mm)
        p2 = scat_vr_param_shaded("SATAD","t",1,1,"1","conv",false,false)
        plot!(p2,xlabel="SAT perturbation",ylabel="temperature [°C]",
              title="b)",title_location=:left,bottom_margin=5mm)
        p3 = scat_vr_param_shaded("SATAD","qv",2,1000,"1","conv",false,false)
        plot!(p3,xlabel="SAT perturbation",ylabel="qv [g/kg]",
              title="c)",title_location=:left,bottom_margin=5mm)

        lay = @layout [a b c]
        p_out = plot(p1,p2,p3,layout=lay,dpi=:300,size=(850,400),titlefontsize=10)
        if save
		savefig(p_out,savepath*"Fig14.png")
	else
		return p_out
	end
end
#----------------------

function plot_histograms_CCN_conv()
    dsa = dsd["accum"]
    t6 = dsa["t600"][:]
    mskc = t6 .< 10

    qi = dsa["qi"][:]
    qni = dsa["qni"][:]
    ri = calc_radius_ice(qi,qni) .* 1000
    ccns = dsa["CCN"][:]
    
    msk1 = (ccns .> 15) .& mskc
    msk2 = (ccns .< 5) .& mskc

    rbins = (1e-3:0.01:0.4)
    qibins = 10 .^ (-4:0.07:0.5)
    
    RHi = dsa["RHi"][:]
    tau = dsa["tau"][:]
    RHbins = 90:0.3:135
    taubins = 10 .^ (0:0.07:6)
    qnibins = 10 .^ (2:0.15:9.5);
    
    ccns = dsa["CCN"][:]
    caps = dsa["CAP"][:]
    
    tauticks = [1e0,1e1,1e2,1e3,1e4,1e5,1e6]
    ni = dsa["qni"][:]

    bp1 = blankp_clean()
    plot!(bp1,[],lw=1.5,c=:red,label=L"CCN $>$ 15")
    plot!(bp1,[],lw=1.5,c=:cyan,label=L"CCN $<$ 5")
    plot!(bp1,legend=:top,legend_columns=2)
        
    p1 = plot(qi[msk1] .* 1000,st=:stephist,
              bins=qibins,c=:red,lw=1.2,label=L"CAP $\geq$ 0.8",
              xlabel=L"$q_i\quad[\mathrm{g/kg}]$",ylabel="counts")
    plot!(p1,qi[msk2] .* 1000,st=:stephist,bins=qibins,c=:cyan,lw=1.2,label="CAP < 0.4")
    plot!(p1,title="a)",title_location=:left,legend=:false,
          xscale=:log10,xticks=[1e-4,1e-3,1e-2,1e-1,1e0])

    p2 = plot(qni[msk1],st=:stephist, bins=qnibins, xscale=:log10, label=L"CAP $\geq$ 0.8",
              c=:red,lw=1.2,
              xlabel=L"$N_i\quad[\mathrm{1/kg}]$",ylabel="counts")
    plot!(p2,qni[msk2],st=:stephist, bins=qnibins, label="CAP < 0.4",c=:cyan,lw=1.2)
    plot!(p2,title="b)",title_location=:left,legend=:false,
          xticks=[1e3,1e5,1e7,1e9])
    
    p3 = plot(ri[msk1],st=:stephist,bins=rbins,label=L"CAP $\geq$ 0.8",c=:red,lw=1.2,
              xlabel=L"$r_i\quad[\mathrm{mm}]$",ylabel="counts")
    plot!(p3,ri[msk2],st=:stephist,bins=rbins,label="CAP < 0.4",c=:cyan,lw=1.2)
    plot!(p3,legend_columns=1,background_color_legend=nothing,
          title="c)",title_location=:left,legend=:false)


    p4 = plot(tau[msk1],st=:stephist,bins=taubins,
               c=:red,lw=1.2,xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",
               ylabel="counts",label=false,xscale=:log10,xticks=tauticks)
    plot!(tau[msk2],lw=1.2,c=:cyan,bins=taubins,label=false,st=:stephist,
          title="d)",title_location=:left)
 
    p5 = plot(RHi[msk1],lw=1.2,c=:red,xlabel=L"$RH_i\quad [\%]$",bins=RHbins,
               st=:stephist,label=false)
    plot!(RHi[msk2],lw=1.2,c=:cyan,st=:stephist,bins=RHbins,label=false,
          title="e)",title_location=:left)

    lay = @layout [o{0.07h}
                   a b c
                   d e]

    p_out = plot(bp1,p1,p2,p3,p4,p5,layout=lay,dpi=:300,size=(800,550),titlefontsize=9,
                bottom_margin=5mm)
    return p_out 
end

#----------------------
function plot_histograms_5h(save=false)
        dsah = dsd["accum_hours"]
        dsa = dsd["accum"]
        qni = dsah["qni"][:]
        RHi = dsah["RHi"][:]
        dpres = dsa["p"][:] .- dsah["p"][:]
        mskdp = dpres .> 0
        inps = dsa["INP"][:]
        caps = dsa["CAP"][:]

        qnibins = 10 .^ (1:0.075:7.8)
        RHibins = 40:1:140
        xticks = [1e1,1e3,1e5,1e7]

        p1 = plot(qni[mskdp .& (inps .< 0.5)],st=:stephist,bins=qnibins,xscale=:log10,
                  c=:cyan,lw=1,label="INP sc. < 0.5",xlabel=L"$N_i\quad[\mathrm{kg}^{-1}]$",
                  ylabel="counts",title="a)",title_location=:left,xticks=xticks)
        plot!(p1,qni[mskdp .& (inps .> 1)],st=:stephist,bins=qnibins,xscale=:log10,
              c=:red,lw=1.5,label="INP sc. > 1",legend=:topright,titlefontsize=10,bottom_margin=5mm)
        
        p3 = plot(RHi[mskdp .& (inps .< 0.5)],st=:stephist,bins=RHibins,
                  c=:cyan,lw=1,label="INP sc. < 0.5",xlabel=L"$RH_i\quad[\%]$",
                  ylabel="counts",title="c)",title_location=:left)
        plot!(p3,RHi[mskdp .& (inps .> 1)],st=:stephist,bins=RHibins,
              c=:red,lw=1.5,label="INP sc. > 1",legend=:topleft,titlefontsize=10,bottom_margin=5mm)

        p2 = plot(qni[mskdp .& (caps .< 0.4)],st=:stephist,bins=qnibins,xscale=:log10,
                  c=:cyan,lw=1,label="CAP < 0.4",xlabel=L"$N_i\quad[\mathrm{kg}^{-1}]$",
                  ylabel="counts",title="b)",title_location=:left,xticks=xticks)
        plot!(p2,qni[mskdp .& (caps .> 0.8)],st=:stephist,bins=qnibins,xscale=:log10,
              c=:red,lw=1.5,label="CAP > 0.8",legend=:topright,titlefontsize=10,bottom_margin=5mm)
        
        p4 = plot(RHi[mskdp .& (caps .< 0.4)],st=:stephist,bins=RHibins,
                  c=:cyan,lw=1,label="CAP < 0.4",xlabel=L"$RH_i\quad[\%]$",
                  ylabel="counts",title="d)",title_location=:left)
        plot!(p4,RHi[mskdp .& (caps .> 0.8)],st=:stephist,bins=RHibins,
              c=:red,lw=1.5,label="CAP > 0.8",legend=:topleft,titlefontsize=10,bottom_margin=5mm)


        lay = @layout [a b
                       c d]
        p_out = plot(p1,p2,p3,p4,layout=lay,dpi=:300,size=(700,600))
        if save
		savefig(p_out,savepath*"Fig15.png")
	else
		return p_out
	end
end
#----------------------

#---------------SUPPLEMENTARY-----------------------
#
function plot_tau_RHi_hist_t600()
        dsa = dsd["accum"] 
        tau_all = dsa["tau"][:]
        RHi_all = dsa["RHi"][:]
        t6 = dsa["t600"][:]
        mskc = t6 .< 10

        p_out = plot_hist2d_2ars_general(tau_all[mskc],RHi_all[mskc],
                                         10 .^ (1:0.1:5),80:1:140,"tau","RHi")
        plot!(p_out,xscale=:log10,legend=false)
        plot!(xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",ylabel=L"$\mathrm{RH}_i\quad[\%]$")
        return p_out
end
#----------------------



function plot_max_qX_CAP_CCN(save=false)
        p1 = scat_vr_param_shaded("CAP","max_qi",2,1000,"1","normal",false,false)
        plot!(p1,title="a)",title_location=:left,ylabel="max(qi) [g/kg]",ylims=(0.1,0.5),xlabel="CAP",
              bottom_margin=5mm,left_margin=5mm)
    
        p2 = scat_vr_param_shaded("CAP","max_qg",2,1000,"1","normal",false,false)
        plot!(p2,title="b)",title_location=:left,ylabel="max(qg) [g/kg]",ylims=(0,1),xlabel="CAP",
             bottom_margin=5mm)

        p3 = scat_vr_param_shaded("CCN","max_qs",2,1000,"1","normal",false,false)
        plot!(p3,title="c)",title_location=:left,ylabel="max(qs) [g/kg]",ylims=(0,1),
              xlabel="CCN scaling",bottom_margin=5mm)
        lay = @layout [a b c]
        p_out = plot(p1,p2,p3,layout=lay,dpi=:300,size=(1000,400),titlefontsize=10)
        return p_out
end


function plot_spearman_fast_vs_slow(save=false)
       
        yt = ["mean Qv","mean T","mean p","mean Qi","mean Ni","mean RHi"]
        xt = ["SST","CCN","INP","CAP","SAT"]
        p1 = plot_cor_vr(["qv","t","p","qi","qni","RHi"],1,"conv")
        plot!(p1,title="a) fast trajectories",title_location=:left,colorbar=false,titlefontsize=10,
              yticks=(1:6,yt),xticks=(1:5,xt))

        p2 = plot_cor_vr(["qv","t","p","qi","qni","RHi"],1,"normal")
        plot!(p2,title="b) all trajectories",title_location=:left,colorbar=false,titlefontsize=10,
              yticks=(1:6,yt),xticks=(1:5,xt))

        cbar = spear_cbar()
        plot!(cbar,bottom_margin=2mm)

        lay = @layout [a b
                       c{0.05h}]

        p_out = plot(p1,p2,cbar,layout=lay,size=(800,400),dpi=:300)
        if save
                savefig(p_out,savepath*"spearmans_fast_slow.png")
        else
                return p_out
        end
end



function plot_scat_vr_after_hours(save=false)
	dsh = dsd["hours"]	

	p1 = scat_vr_param_shaded("CAP","RHi",1,1,"1","hours",false,false)
	plot!(p1,legend=:top,ylabel="RHi(5h) [%]",title="a)",title_location=:left)
	p2 = scat_vr_param_shaded("INP","RHi",1,1,"1","hours",false,false)
	plot!(p2,legend=:top,ylabel="RHi(5h) [%]",title="b)",title_location=:left,
             xscale=:log10)

	p3 = scat_vr_param_shaded("CAP","qi",2,1e6,"1","hours",false,false)
	plot!(p3,legend=:top,ylabel="qi(5h) [microg/kg]",title="c)",title_location=:left,
              yscale=:log10,ylims=(0.0001,25),yticks=[1e-3,1e-2,1e-1,1e0,1e1])
	p4 = scat_vr_param_shaded("INP","qi",2,1e6,"1","hours",false,false)
	plot!(p4,legend=:top,ylabel="qi(5h) [microg/kg]",title="d)",title_location=:left,
              yscale=:log10,ylims=(0.0001,25),yticks=[1e-3,1e-2,1e-1,1e0,1e1],
              xscale=:log10)
	
	p5 = scat_vr_param_shaded("CAP","qni",2,1,"1","hours",false,false)
	plot!(p5,legend=:top,ylabel="Ni(5h) [1/kg]",title="e)",title_location=:left,
              yscale=:log10,ylims=(0.01,1e5))
	p6 = scat_vr_param_shaded("INP","qni",2,1,"1","hours",false,false)
	plot!(p6,legend=:top,ylabel="Ni(5h) [1/kg]",title="f)",title_location=:left,
              xscale=:log10,yscale=:log10,ylims=(0.01,1e5))


	lay = @layout [a b
			 c d
			 e f]
	p_out = plot(p1,p2,p3,p4,p5,p6,layout=lay,size=(800,800),titlefontsize=10)

        return p_out
end


function plot_tau_RHi_mean()

        xbins = 10 .^ (1:0.1:5)
        ybins = 80:1:140
        ylims=(95,130)

        dsa = dsd["accum"]
        tau_all = dsa["tau"][:]
        RHi_all = dsa["RHi"][:]
        CAP_all = dsa["CAP"][:]
        INP_all = dsa["INP"][:]
        CCN_all = dsa["CCN"][:]
        qc = dsa["qc"][:]

        mskqc = qc .== 0
        #---masking CAP---
        msk1 = (CAP_all .>= 0.8) .& mskqc
        msk2 = (CAP_all .< 0.4) .& mskqc

        RHi1 = RHi_all[msk1]
        tau1 = tau_all[msk1]
        f10_1 = anyar_prcntl(RHi1,tau1,xbins,25)
        f90_1 = anyar_prcntl(RHi1,tau1,xbins,75)
        
        RHi2 = RHi_all[msk2]
        tau2 = tau_all[msk2]
        f10_2 = anyar_prcntl(RHi2,tau2,xbins,25)
        f90_2 = anyar_prcntl(RHi2,tau2,xbins,75)

        #---masking INP---
        msk3 = (INP_all .>= 1) .& mskqc
        msk4 = (INP_all .< 1) .& mskqc        

        RHi3 = RHi_all[msk3]
        tau3 = tau_all[msk3]
        f10_3 = anyar_prcntl(RHi3,tau3,xbins,25)
        f90_3 = anyar_prcntl(RHi3,tau3,xbins,75)
        
        RHi4 = RHi_all[msk4]
        tau4 = tau_all[msk4]
        f10_4 = anyar_prcntl(RHi4,tau4,xbins,25)
        f90_4 = anyar_prcntl(RHi4,tau4,xbins,75)
        
        bp1 = blankp_clean()
        plot!(bp1,[],lw=1.5,c=:red,label=L"CAP $\geq$ 0.8")
        plot!(bp1,[],lw=1.5,c=:cyan,label=L"CAP $<$ 0.4")
        plot!(bp1,legend=:top,legend_columns=2)

        bp2 = blankp_clean()
        plot!(bp2,[],lw=1.5,c=:red,label=L"INP $\geq$ 1")
        plot!(bp2,[],lw=1.5,c=:cyan,label=L"INP $<$ 1")
        plot!(bp2,legend=:top,legend_columns=2)

        p1 = plot(xbins[2:end],bin_two_ars(RHi1,tau1,xbins,nm.mean),
                     fillrange=(f10_1,f90_1),fillalpha=0.2,
                     lw=2,c=:red,label=false,xscale=:log10,xticks=[10,1e2,1e3,1e4,1e5],
                     title="a)",title_location=:left,
                    bottom_margin=5mm,left_margin=2mm)
        plot!(p1,xbins[2:end],bin_two_ars(RHi_all[msk1],tau_all[msk1],xbins,nm.median),
                     lw=2,c=:red,label=false,ls=:dash)
        plot!(p1,xbins[2:end],bin_two_ars(RHi_all[msk2],tau_all[msk2],xbins,nm.mean),
                     fillrange=(f10_2,f90_2),fillalpha=0.3,
                     lw=2,c=:cyan,label=false)
        plot!(p1,xbins[2:end],bin_two_ars(RHi_all[msk2],tau_all[msk2],xbins,nm.median),
                     lw=2,c=:cyan,label=false,ls=:dash)
        plot!(p1,xlabel=L"$\tau_\mathrm{sat,ice}\quad [\mathrm{s}]$",
              ylabel=L"$\mathrm{RH}_i\quad [\%]$",ylims=ylims,
              xlims=(xbins[1],xbins[end]))
       
        p2 = plot(xbins[2:end],bin_two_ars(RHi_all[msk3],tau_all[msk3],xbins,nm.mean),
                     fillrange=(f10_3,f90_3),fillalpha=0.2,
                     lw=2,c=:red,label=false,xscale=:log10,xticks=[10,1e2,1e3,1e4,1e5],
                     title="b)",title_location=:left,bottom_margin=5mm,left_margin=2mm)
        plot!(p2,xbins[2:end],bin_two_ars(RHi_all[msk3],tau_all[msk3],xbins,nm.median),
                     lw=2,c=:red,label=false,ls=:dash)
        plot!(p2,xbins[2:end],bin_two_ars(RHi_all[msk4],tau_all[msk4],xbins,nm.mean),
                     fillrange=(f10_4,f90_4),fillalpha=0.2,
                     lw=2,c=:cyan,label=false)
        plot!(p2,xbins[2:end],bin_two_ars(RHi_all[msk4],tau_all[msk4],xbins,nm.median),
                     lw=2,c=:cyan,label=false,ls=:dash)
        plot!(p2,xlabel=L"$\tau_\mathrm{sat,ice}\quad [\mathrm{s}]$",
              ylabel=L"$\mathrm{RH}_i\quad [\%]$",ylims=ylims,
             xlims=(xbins[1],xbins[end]))
        
        lay = @layout [a{0.1h} b{0.1h}
                       c d]
        p_out = plot(bp1,bp2,p1,p2,layout=lay,
                     dpi=:300,titlefontsize=10,size=(800,400))
        return p_out
end


function plot_RHi_hist_2msks()
        dsa = dsd["accum"]
        caps = dsa["CAP"][:]
        inps = dsa["INP"][:]
        RHi = dsa["RHi"][:]
        tau = dsa["tau_gl"][:]

        msk1 = (caps .< 0.4) .& (inps .< 1);
        msk2 = (caps .>= 0.8) .& (inps .>= 1);
        taubins = 10 .^ (1:0.1:6)

        bp1 = blankp_clean()
        plot!(bp1,[],lw=1.5,c=:red,label="CAP > 0.8 and INP > 1")
        plot!(bp1,[],lw=1.5,c=:cyan,label="CAP < 0.4 and INP < 1")
        plot!(bp1,legend=:top,legend_columns=2)

        p1 = plot(tau[msk1],st=:stephist,bins=taubins,xscale=:log10,c=:cyan,
                  xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",lw=1.5,
                  ylabel="counts",title="a)",title_location=:left,label=false)
        plot!(p1,tau[msk2],st=:stephist,bins=taubins,label=false,c=:red,lw=1.5,
              xticks=[10,1e2,1e3,1e4,1e5,1e6])

        p2 = plot(sort(tau[msk1]),(1:sum(msk1)) ./ sum(msk1),c=:cyan,label="CAP < 0.4 and INP < 1")
        plot!(p2,sort(tau[msk2]),(1:sum(msk2)) ./ sum(msk2),c=:red,label="CAP > 0.4 and INP > 1")
        plot!(p2,xlims=(taubins[1],taubins[end]),xscale=:log10,title="b)",title_location=:left,legend=false,
              xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$")


        p3 = plot(RHi[msk1],st=:stephist,normalize=false,label="CAP < 0.4 and INP < 1",
             xlims=(85,140),c=:cyan,lw=1.5,bins=(90:0.5:140),
             xlabel=L"$\mathrm{RH}_i\quad[\%]$",ylabel="counts")
        plot!(p3,RHi[msk2],st=:stephist,normalize=false,
              label="CAP > 0.4 and INP > 1",xlims=(85,140),c=:red,lw=1.5,bins=(90:0.5:140))
        plot!(p3,title="c)",title_location=:left,legend=false)

        p4 = plot(sort(RHi[msk1]),(1:sum(msk1)) ./ sum(msk1),c=:cyan,label="CAP < 0.4 and INP < 1")
        plot!(p4,sort(RHi[msk2]),(1:sum(msk2)) ./ sum(msk2),c=:red,label="CAP > 0.4 and INP > 1")
        plot!(p4,xlims=(85,140),title="d)",title_location=:left,legend=false,
              xlabel=L"$\mathrm{RH}_i\quad[\%]$")

        lay = @layout [o{0.1h}
                       a b
                       c d]
        p_out = plot(bp1,p1,p2,p3,p4,layout=lay,dpi=:300,size=(800,600),titlefontsize=10)
        return p_out
end

function plot_t_p_conv()
       dsa = dsd["accum"]
       dsc = dsd["conv"]
       t6 = dsa["t600"][:]
       sim = dsa["sim"][:]
       mskc = t6 .< 10
       msks = t6 .> 10

       pbins = 170:6:410
       tbins = 200:2:255
       #qvbin = 
       
       ix_pminmean = argmin(dsc["p"][1,:])
       ix_pminmed = argmin(dsc["p"][2,:])
       ix_pmaxmean = argmax(dsc["p"][1,:])
       ix_pmaxmed = argmax(dsc["p"][2,:])
       
       pminmean = dsa["p"][(sim .== ix_pminmean) .& mskc]
       pminmed = dsa["p"][(sim .== ix_pminmed) .& mskc]
       pmaxmean = dsa["p"][(sim .== ix_pmaxmean) .& mskc]
       pmaxmed = dsa["p"][(sim .== ix_pmaxmed) .& mskc]

       ix_tminmean = argmin(dsc["t"][1,:])
       ix_tminmed = argmin(dsc["t"][2,:])
       ix_tmaxmean = argmax(dsc["t"][1,:])
       ix_tmaxmed = argmax(dsc["t"][2,:])
       
       tminmean = dsa["t"][(sim .== ix_tminmean) .& mskc]
       tminmed = dsa["t"][(sim .== ix_tminmed) .& mskc]
       tmaxmean = dsa["t"][(sim .== ix_tmaxmean) .& mskc]
       tmaxmed = dsa["t"][(sim .== ix_tmaxmed) .& mskc]
       
       p1 = plot(pminmean,st=:stephist,normalize=true,bins=pbins,label="min mean")
       plot!(p1,pmaxmed,st=:stephist,normalize=true,bins=pbins,label="max median")
       plot!(p1,title="a)",title_location=:left,
             xlabel="pressure [hPa]",ylabel="normalized counts",
            legend_columns=2,legend=:outertop)

       p2 = plot(tminmean,st=:stephist,normalize=true,bins=tbins,label="min mean")
       plot!(p2,tmaxmean,st=:stephist,normalize=true,bins=tbins,label="max mean")
       plot!(p2,title="b)",title_location=:left,
             xlabel="temperature [K]",ylabel="normalized counts",
            legend_columns=2,legend=:outertop)
     
        print(size(tmaxmean),size(tmaxmed))

       lay = @layout [a b]
       p_out = plot(p1,p2,layout=lay,dpi=:300)
       return p_out

end


function plot_RHi_tau_CCN_conv()
    dsa = dsd["accum"]
    RHi = dsa["RHi"][:]
    tau = dsa["tau"][:]
    RHbins = 90:0.5:135
    taubins = 10 .^ (0.5:0.1:6)
    qnibins = 10 .^ (2:0.05:9);
    t6 = dsa["t600"][:]
    mskc = t6 .< 10
    ccns = dsa["CCN"][:]
    caps = dsa["CAP"][:]
    tauticks = [1e1,1e3,1e5]
    ni = dsa["qni"][:]


    p1 = scat_vr_param_shaded("CCN","tau",2,1,"s","conv",false,false)
    plot!(ylabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",xlabel="CCN scaling",
          title="a)",title_location=:left,yscale=:log10,yticks=[1e1,1e2,1e3,1e4])
    vline!([5],label=false,lw=1.5,c=:cyan)
    vline!([15],label=false,lw=1.5,c=:red)
    p1b = plot(tau[mskc .& (ccns .> 15)],st=:stephist,bins=taubins,
               c=:red,lw=1.5,xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",
               ylabel="counts",label=false,xscale=:log10,xticks=tauticks)
    plot!(tau[mskc .& (ccns .< 5)],lw=1.5,c=:cyan,bins=taubins,label=false,st=:stephist,
          title="e)",title_location=:left)

    p2 = scat_vr_param_shaded("CCN","qni",2,1,"s","conv",false,false)
    plot!(ylabel=L"$N_i\quad[\mathrm{kg}^{-1}]$",xlabel="CAP",
          title="b)",title_location=:left,yscale=:log10,label=false)

    p2b = plot(ni[mskc .& (ccns .> 5)],st=:stephist,bins=qnibins,
               xscale=:log10,c=:red,lw=1.5,xlabel=L"$N_i\quad[1/\mathrm{kg}]$",
               ylabel="counts",label=false)
    plot!(ni[mskc .& (ccns .< 15)],lw=1.5,c=:cyan,bins=qnibins,st=:stephist,label=false,
          title="f)",title_location=:left)
   
    p3 = scat_vr_param_shaded("CCN","RHi",1,1,"1","conv",false,false)
    plot!(ylabel="RHi [%]",xlabel="CCN scaling",title="c)",title_location=:left)
    vline!([5],label=false,lw=1.5,c=:cyan)
    vline!([15],label=false,lw=1.5,c=:red)
    
    p3b = plot(RHi[mskc .& (ccns .> 15)],lw=1.5,c=:red,xlabel="RHi [%]",bins=RHbins,
               st=:stephist,label=false)
    plot!(RHi[mskc .& (ccns .< 5)],lw=1.5,c=:cyan,st=:stephist,bins=RHbins,label=false,
          title="g)",title_location=:left)

    lay = @layout [a b c
                   d e f]

    p_out = plot(p2,p1,p3,p2b,p1b,p3b,layout=lay,dpi=:300,size=(800,600))
    return p_out 
end


function plot_histograms_CCN_conv_1(save=false)
        dsa = dsd["accum"]
        t6 = dsa["t600"][:]
        msk = t6 .< 10

        qi = dsa["qi"][:]
        qni = dsa["qni"][:]
        ri = calc_radius_ice(qi,qni) .* 1000
        ccns = dsa["CCN"][:]
        
        msk1 = (ccns .> 15) .& msk
        msk2 = (ccns .< 5) .& msk

        rbins = (1e-3:0.01:0.4)
        qnibins = 10 .^ (2:0.1:9)
        taubins = 10 .^ (1:0.1:6)
        qibins = 10 .^ (-4:0.05:0.5)

        bp1 = blankp_clean()
        plot!(bp1,[],lw=1.5,c=:red,label=L"CCN $>$ 15")
        plot!(bp1,[],lw=1.5,c=:cyan,label=L"CAP $<$ 5")
        plot!(bp1,legend=:top,legend_columns=2)

	p1 = scat_vr_param_shaded("CCN","qi",2,1000,"g/kg","normal",false,false,:log10,(6e-4,3e-1))
	plot!(p1,title="a)",title_location=:left,ylabel=L"q_\mathrm{i}\quad [\mathrm{g/kg}]",
	       xlabel=L"\mathrm{CCN}",)	
        
        p2 = scat_vr_param_shaded("CCN","qni",2,1,"1/kg","normal",false,false,:log10,(1e3,5e7),true)
	plot!(p2,title="b)",title_location=:left,ylabel=L"N_\mathrm{i}\quad [1/\mathrm{kg}]",
               xlabel=L"\mathrm{CCN}",yticks=[1e3,1e4,1e5,1e6,1e7])
        
        p3 = scat_vr_param_shaded("CCN","r",2,1000,"mm","normal",false,false)
	plot!(p3,title="c)",title_location=:left,ylabel=L"r_\mathrm{i}\quad [\mathrm{mm}]",
	       xlabel=L"\mathrm{CCN}")

        p4 = plot(qi[msk1] .* 1000,st=:stephist,
                  bins=qibins,c=:red,lw=1.5,label=L"CAP $\geq$ 0.8",
                  xlabel=L"$q_i\quad[\mathrm{g/kg}]$",ylabel="counts")
        plot!(p4,qi[msk2] .* 1000,st=:stephist,bins=qibins,c=:cyan,lw=1.5,label="CAP < 0.4")
        plot!(p4,title="d)",title_location=:left,legend=:false,
              xscale=:log10,xticks=[1e-4,1e-3,1e-2,1e-1,1e0])

        p5 = plot(qni[msk1],st=:stephist, bins=qnibins, xscale=:log10, label=L"CAP $\geq$ 0.8",
                  c=:red,lw=1.5,
                  xlabel=L"$N_i\quad[\mathrm{1/kg}]$",ylabel="counts")
        plot!(p5,qni[msk2],st=:stephist, bins=qnibins, label="CAP < 0.4",c=:cyan,lw=1.5)
        plot!(p5,title="e)",title_location=:left,legend=:false)
        
        p6 = plot(ri[msk1],st=:stephist,bins=rbins,label=L"CAP $\geq$ 0.8",c=:red,lw=1.5,
                  xlabel=L"$r_i\quad[\mathrm{mm}]$",ylabel="counts")
        plot!(p6,ri[msk2],st=:stephist,bins=rbins,label="CAP < 0.4",c=:cyan,lw=1.5)
        plot!(p6,legend_columns=1,background_color_legend=nothing,
              title="f)",title_location=:left,legend=:false)


        lay = @layout [a b c
                       o{0.05h}
                       d e f]

        p_out = plot(p1,p2,p3,bp1,p4,p5,p6,layout=lay,
                     titlefontsize=10,size=(800,600),dpi=:300)

	if save
		savefig(p_out,savepath*"Fig4.png")
	else
		return p_out
	end
end


function plot_histograms_CCN_conv_2()
    dsa = dsd["accum"]
    
    RHi = dsa["RHi"][:]
    tau = dsa["tau"][:]
    RHbins = 90:0.5:135
    taubins = 10 .^ (0.5:0.1:6)
    qnibins = 10 .^ (2:0.05:9);

    t6 = dsa["t600"][:]
    mskc = t6 .< 10
    
    ccns = dsa["CCN"][:]
    caps = dsa["CAP"][:]
    
    tauticks = [1e1,1e3,1e5]
    ni = dsa["qni"][:]

    bp1 = blankp_clean()
    plot!(bp1,[],lw=1.5,c=:red,label=L"CCN $>$ 15")
    plot!(bp1,[],lw=1.5,c=:cyan,label=L"CAP $<$ 5")
    plot!(bp1,legend=:top,legend_columns=2)

    p1 = scat_vr_param_shaded("CCN","tau",2,1,"s","conv",false,false)
    plot!(ylabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",xlabel="CCN scaling",
          title="a)",title_location=:left,yscale=:log10,yticks=[1e1,1e2,1e3,1e4])
    p1b = plot(tau[mskc .& (ccns .> 15)],st=:stephist,bins=taubins,
               c=:red,lw=1.5,xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",
               ylabel="counts",label=false,xscale=:log10,xticks=tauticks)
    plot!(tau[mskc .& (ccns .< 5)],lw=1.5,c=:cyan,bins=taubins,label=false,st=:stephist,
          title="c)",title_location=:left)
 
    p2 = scat_vr_param_shaded("CCN","RHi",1,1,"1","conv",false,false)
    plot!(ylabel="RHi [%]",xlabel="CCN scaling",title="b)",title_location=:left)
    
    p2b = plot(RHi[mskc .& (ccns .> 15)],lw=1.5,c=:red,xlabel="RHi [%]",bins=RHbins,
               st=:stephist,label=false)
    plot!(RHi[mskc .& (ccns .< 5)],lw=1.5,c=:cyan,st=:stephist,bins=RHbins,label=false,
          title="d)",title_location=:left)

    lay = @layout [a b
                   o{0.08h}
                   d e]

    p_out = plot(p1,p2,bp1,p1b,p2b,layout=lay,dpi=:300,size=(600,500),titlefontsize=9)
    return p_out 
end


function plot_dp_RHi_after_asc()
        dsa = dsd["accum"]
        dsah = dsd["accum_hours"]

        RHi = dsah["RHi"][:]
        dp = dsa["p"][:] .- dsah["p"][:]

        p_out = plot_hist2d_2ars_general(dp,RHi,-150:6:150,10:1.5:135,"","")
        plot!(xlabel=L"$p(t_\mathrm{WCB})-p(t_\mathrm{WCB}+ 5\,\mathrm{h})\quad[\mathrm{hPa}]$",
              ylabel=L"$\mathrm{RH}_i(t_\mathrm{WCB} + 5\,\mathrm{h})\quad[\%]$",
              legend=false,labelfontsize=10)
        return p_out

end

function plot_qc_35deg_CAP()
        p1 = plot([],label=false,xlabel="CAP",ylabel="Number of trajs",
                 title="Number of trajs with qc > 0 at -35°C")

        for i in 1:70
            print("$(i), ")
            mem = i + 19
            T = load_vr("t",mem)
            qc = load_vr("qc",mem)

            inds = [findfirst(T[k,:] .< (-35 + 273.15)) for k in 1:size(T)[1]]
            counter = 0
            for j in 1:length(inds)
                if inds[j] != nothing
                    if qc[j,inds[j]] > 0
                        counter = counter + 1
                    end
                end
            end
            scatter!(p1,[dp["CAP"][i]],[counter],label=false,c=:blue,markerstrokewidth=0)
        end
        return p1
end

function plot_qc_35deg_INP()
        p1 = plot([],label=false,xlabel="INP",ylabel="Number of trajs",
                 title="Number of trajs with qc > 0 at -35°C")

        for i in 1:70
            print("$(i), ")
            mem = i + 19
            T = load_vr("t",mem)
            qc = load_vr("qc",mem)

            inds = [findfirst(T[k,:] .< (-35 + 273.15)) for k in 1:size(T)[1]]
            counter = 0
            for j in 1:length(inds)
                if inds[j] != nothing
                    if qc[j,inds[j]] > 0
                        counter = counter + 1
                    end
                end
            end
            scatter!(p1,[dp["INP"][i]],[counter],label=false,c=:blue,markerstrokewidth=0)
        end
        plot!(p1,xscale=:log10)
        return p1
end

function plot_qnc_35deg_INP()
        p1 = plot([],label=false,xlabel="INP",ylabel="Number of trajs",
                 title="Nc average at -35°C")

        for i in 1:1:70
            print("$(i), ")
            mem = i + 19
            T = load_vr("t",mem)
            qc = load_vr("qc",mem)
            Nc = load_vr("qnc",mem)
            inds = [findfirst(T[k,:] .< (-35 + 273.15)) for k in 1:size(T)[1]]
            counter = 0
            for j in 1:length(inds)
                if inds[j] != nothing
                    if qc[j,inds[j]] > 0
                        counter = counter + Nc[j,inds[j]]
                    end
                end
            end
            scatter!(p1,[dp["INP"][i]],[counter],label=false,c=:blue,markerstrokewidth=0)
        end
        plot!(p1,xscale=:log10,yscale=:log10)
        return p1
end

function plot_qnc_35deg_CAP()
        p1 = plot([],label=false,xlabel="CAP",ylabel="Number of trajs",
                 title="Nc average at -35°C")

        for i in 1:1:70
            print("$(i), ")
            mem = i + 19
            T = load_vr("t",mem)
            qc = load_vr("qc",mem)
            Nc = load_vr("qnc",mem)
            inds = [findfirst(T[k,:] .< (-35 + 273.15)) for k in 1:size(T)[1]]
            counter = 0
            for j in 1:length(inds)
                if inds[j] != nothing
                    if qc[j,inds[j]] > 0
                        counter = counter + Nc[j,inds[j]]
                    end
                end
            end
            scatter!(p1,[dp["CAP"][i]],[counter],label=false,c=:blue,markerstrokewidth=0)
        end
        plot!(p1,yscale=:log10)
        return p1
end

function plot_histograms_SST(msk,save=false)
        dsa = dsd["accum"]
        qi = dsa["qi"][:]
        qni = dsa["qni"][:]
        ri = calc_radius_ice(qi,qni) .* 1000
        ssts = dsa["SST"][:]
        
        msk3 = (ssts .>= 1) .& msk
        msk4 = (ssts .< -1) .& msk

        rbins = (1e-3:0.01:0.4)
        qnibins = 10 .^ (2:0.1:9)
        taubins = 10 .^ (1:0.1:6)
        qibins = 10 .^ (-4:0.05:0.5)

        bp1 = blankp_clean()
        plot!(bp1,[],lw=1.5,c=:red,label=L"SST pert. $\geq$ 1")
        plot!(bp1,[],lw=1.5,c=:cyan,label=L"SST pert. $<$ -1")
        plot!(bp1,legend=:top,legend_columns=2)
        
        p1 = scat_vr_param_shaded("SST","qi",2,1000,"g/kg","normal",false,false,:log10,(6e-4,3e-1))
	plot!(p1,title="a)",title_location=:left,ylabel=L"q_\mathrm{i}\quad [\mathrm{g/kg}]",
	       xlabel=L"\mathrm{SST\,\,pert.}")
        plot!(p1,xscale=:identity,xlims=(-2,2),xticks=collect(-2:1:2))
        
        p2 = scat_vr_param_shaded("SST","qni",2,1,"1/kg","normal",false,false,:log10,(1e3,5e7),true)
	plot!(p2,title="b)",title_location=:left,ylabel=L"N_\mathrm{i}\quad [1/\mathrm{kg}]",
               xlabel=L"\mathrm{SST\,\,pert.}",yticks=[1e3,1e4,1e5,1e6,1e7])
        plot!(p2,xscale=:identity,xlims=(-2,2),xticks=collect(-2:1:2))
	
        p3 = scat_vr_param_shaded("SST","r",2,1000,"mm","normal",false,false)
	plot!(p3,title="c)",title_location=:left,ylabel=L"r_\mathrm{i}\quad [\mathrm{mm}]",
	       xlabel=L"\mathrm{SST\,\,pert.}")
        plot!(p3,xscale=:identity,xlims=(-2,2),xticks=collect(-2:1:2))

        p4 = plot(qi[msk3] .* 1000,st=:stephist,
                  bins=qibins,c=:red,lw=1.5,
                  xlabel=L"$q_i\quad[\mathrm{g/kg}]$",ylabel="counts")
        plot!(p4,qi[msk4] .* 1000,st=:stephist,bins=qibins,c=:cyan,lw=1.5,label="CAP < 0.4")
        plot!(p4,title="d)",title_location=:left,legend=:false,
              xscale=:log10,xticks=[1e-4,1e-3,1e-2,1e-1,1e0])

        p5 = plot(qni[msk3],st=:stephist, bins=qnibins, xscale=:log10, label=L"CAP $\geq$ 0.8",
                  c=:red,lw=1.5,
                  xlabel=L"$N_i\quad[\mathrm{1/kg}]$",ylabel="counts")
        plot!(p5,qni[msk4],st=:stephist, bins=qnibins,c=:cyan,lw=1.5)
        plot!(p5,title="e)",title_location=:left,legend=:false)
        
        p6 = plot(ri[msk3],st=:stephist,bins=rbins,c=:red,lw=1.5,
                  xlabel=L"$r_i\quad[\mathrm{mm}]$",ylabel="counts")
        plot!(p6,ri[msk4],st=:stephist,bins=rbins,c=:cyan,lw=1.5)
        plot!(p6,legend_columns=1,background_color_legend=nothing,
              title="f)",title_location=:left,legend=:false)


        lay = @layout [a b c
                       o{0.05h}
                       d e f]

        p_out = plot(p1,p2,p3,bp1,p4,p5,p6,layout=lay,
                     titlefontsize=10,size=(800,600),dpi=:300)

	if save
		savefig(p_out,savepath*"Fig5.png")
	else
		return p_out
	end
end


#------------------------------------------------------
#
function plot_SI1()
        p1 = scat_vr_param_shaded("CAP","T_gl99",1,1,"1","normal",false,false)
        plot!(p1,title="a)",title_location=:left,xlabel="CAP",ylabel="T @ full glaciation [°C]")
        p2 = scat_vr_param_shaded("CAP","p_gl99",1,1,"1","normal",false,false)
        plot!(p2,title="b)",title_location=:left,xlabel="CAP",ylabel="p @ full glaciation [hPa]")
        p3 = scat_vr_param_shaded("INP","T_gl99",1,1,"1","normal",false,false)
        plot!(p3,title="c)",title_location=:left,xlabel="INP scaling factor",xscale=:log10,
              ylabel="T @ full glaciation [°C]",xticks=[1e-2,1e-1,1e0,1e1,1e2])
        p4 = scat_vr_param_shaded("INP","p_gl99",1,1,"1","normal",false,false)
        plot!(p4,title="d)",title_location=:left,xlabel="INP scaling factor",xscale=:log10,
              ylabel="p @ full glaciation [hPa]",xticks=[1e-2,1e-1,1e0,1e1,1e2])

        lay = @layout [a b
                       c d]
        p_out = plot(p1,p2,p3,p4,layout=lay,dpi=:300,size=(700,600),titlefontsize=10,labelfontsize=10)
        return p_out
end

function plot_SI2()
    p1 = plot_peak_Ni_param("CAP")
    plot!(p1,title="a)",title_location=:left)
    p2 = plot_peak_Ni_param("INP")
    plot!(p2,title="b)",title_location=:left,xscale=:log10,xticks=[1e-2,1e-1,1e0,1e1,1e2])

    lay = @layout [a b]
    p_out = plot(p1,p2,layout=lay,dpi=:300,size=(600,400),bottom_margin=2mm)
    return p_out
end

function plot_peak_Ni_param(param::String)
        dsa = Dataset("../DATA/PPE_vars_accum.nc");
        Ni = dsa["qni"][:]
        sims = dsa["sim"][:]
        qnibins = 10 .^ (0:0.1:9)
        peak_Ni = zeros(70)
        for i in 1:70
            msk = sims .== i
            p1 = plot(Ni[msk],st=:stephist,bins=qnibins)
            peak_Ni[i] = p1[1][1][:x][argmax(p1[1][1][:y])]
        end
        closeall()
        p_out = scatter(dp[param],peak_Ni,xlabel=param,ylabel="peak Ni",yscale=:log10,
                       label=false)
        scatter!(p_out,[dp[param][33]],[peak_Ni[33]],label=false)
        return p_out
end

#----------------------
#RF model RHi prediction
function plot_SI3(save=false,typ="normal")
        ds = dsd[typ]
        X = Matrix([dp["CAP"] dp["INP"] dp["CCN"]])
        p1 = scatter_forest(X, ["CAP","INP","CCN"], ds["RHi"][1,:], "RHi [%]")
	
	p_tit = blankp_clean()
	annotate!(p_tit,0.5,0.5,L"$\mathrm{RF}^2_\mathrm{mean(RHi)}$",annotation_fontsize=15)

	lay = @layout [a{0.1h}
			b]
	
	p_out = plot(p_tit,p1,layout=lay,dpi=:300)
        
	if save
                savefig(p_out,savepath*"SI3.png")
        else
                return p_out
        end
end

#----------------------
#RHi vs RHi_glac
function plot_SI4()
        dsa = Dataset("../DATA/PPE_vars_accum.nc");
        p_out = histogram2d(dsa["RHi"][:],dsa["RHi_gl"][:],bins=(60:1:140,60:1:140),
                            xlabel="RHi @ end of ascent [%]",
                            ylabel="RHi @ full glaciation [%]",
                            colorbartitle="\n counts",right_margin=8mm)

        return p_out
end

#----------------------
#tau changing for CCN histogram
function plot_SI5()
        dsa = Dataset("../DATA/PPE_vars_accum.nc");
        ccns = dsa["CCN"][:]
        msk1 = ccns .< 5
        msk2 = ccns .> 15
        taubins = 10 .^ (0.5:0.1:6)
        tau = dsa["tau"][:]
        xticks = [1e1,1e2,1e3,1e4,1e5,1e6]
        p1 = plot(tau[msk1],st=:stephist,bins=taubins,xscale=:log10,label="CCN < 5",
                  c=:cyan,lw=1.5,xlabel=L"$\tau_\mathrm{sat,ice}\quad[\mathrm{s}]$",ylabel="counts",
                  xticks=xticks)
        plot!(p1,tau[msk2],st=:stephist,bins=taubins,label="CCN > 15",c=:red,lw=1.5)
        return p1
end


#----------------------
#Ni distribution 5h after the ascent for different RHi
function plot_SI6()
        dsah = Dataset("../DATA/PPE_vars_accum_hours.nc")
        RHi = dsah["RHi"][:]
        qni = dsah["qni"][:]
        qnibins =  10 .^ (0:0.05:8)

        msk1 = (RHi .> 95) .& (RHi .< 105)
        msk2 = (RHi .> 105)
        msk3 = (RHi .< 95)

        p_out = plot(qni[msk1],st=:stephist,bins=qnibins,xscale=:log10,
                     label=L"$RH_i(5\mathrm{h}) \in (95,105)$",
                     xlabel=L"$N_i\quad [1/\mathrm{kg}]$",c=:black,lw=1.5,
                     ylabel="counts",legendfontsize=10,dpi=:300,
                     xticks=[1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8])
        plot!(p_out,qni[msk2],st=:stephist,bins=qnibins,label=L"$RH_i(5\mathrm{h}) > 105$",c=:red,lw=1.5)
        plot!(p_out,qni[msk3],st=:stephist,bins=qnibins,label=L"$RH_i(5\mathrm{h}) < 95$",c=:blue,lw=1.5)
        
        return p_out
end


function plot_SI7()
        dsa = dsd["accum"]
        tgl = dsa["t_gl"][:] .- 273.15
        pgl = dsa["p_gl"][:]
        ssts = dsa["SST"][:]

        p1 = plot(tgl[ssts .< 0],st=:stephist,bins=-60:0.5:20,
                  normalize=true,xlabel="Temperature @ glaciation [°C]",
                  ylabel="normalized counts",label="SST < 0",
                  c=:cyan,lw=1,title="a)",title_location=:left)
        plot!(p1,tgl[ssts .>= 0],st=:stephist,bin=-60:0.5:20,
              normalize=true,c=:red,lw=1,label="SST >= 0")

        p2 = plot(pgl[ssts .< 0],st=:stephist,bins=200:3:700,
                  normalize=true,xlabel="Pressure @ glaciation [hPa]",
                  ylabel="normalized counts",label="SST < 0",
                  c=:cyan,lw=1,title="b)",title_location=:left)
        plot!(p2,pgl[ssts .>= 0],st=:stephist,normalize=true,bins=200:3:700,
              c=:red,lw=1,label="SST >= 0")

        lay = @layout [a b]
        p_out = plot(p1,p2,layout=lay,dpi=:300,size=(800,400))
        return p_out
end







