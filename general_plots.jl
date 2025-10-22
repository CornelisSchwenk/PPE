include("calcs.jl")
savepath = "../PPE_code/PLOTS/"

taulabel = L"\tau_{600}\quad [\mathrm{h}]"

labeldict=Dict()
labeldict[1] = "mean"
labeldict[2] = "median"
labeldict[3] = "std"
labeldict[4] = "5th"
labeldict[5] = "10th"
labeldict[6] = "90th"
labeldict[7] = "95th"

#------------------------------------------
#2d histograms
#------------------------------------------

#-----"normalize 2dim histogram"---
function normalize_hist2d(hist,dim=1)
        zh = hist[1][1][:z]
        zh2 = replace(reshape(collect(zh),size(zh)),NaN => 0.0)
        norm = sum(zh2,dims=dim)
        normed = replace(zh2 ./ norm,0.0 => NaN)
        hist[1][1][:z] = Surface{Matrix{Float64}}(normed)
        hist[1][1][:clims_calculated] = (0.0,1.0)
end

#-----"general 2dim histogram any two arrays"---
function plot_hist2d_2ars_general(ar1,ar2,xbins,ybins,xlabel,ylabel,loc=:top,pl_mm=true)
        h_ar = histogram2d((ar1,ar2),clims=(0,0.1),xlims=(xbins[1],xbins[end]),
                           ylims=(ybins[1],ybins[end]),xlabel=xlabel,ylabel=ylabel,
                           bins=(xbins,ybins),labelfontsize=14,
                           colorbartitle="\n Normalised number of trajectories",
                           right_margin=5mm,legend=loc)
        normalize_hist2d(h_ar)
        if pl_mm == true
                plot!(xbins[2:end],bin_two_ars(ar2,ar1,xbins,nm.mean),
                       label="mean",lw=:2,c=:cyan)
                plot!(xbins[2:end],bin_two_ars(ar2,ar1,xbins,nm.median),
                       label="median",lw=:2,c=:cyan,ls=:dash)
        end
        return h_ar
end

function plot_mean_median_2ars(ar1,ar2,xbins,ybins,xlabel,ylabel)
        p_out = plot(xbins[2:end],bin_two_ars(ar2,ar1,xbins,nm.mean),lw=2,c=:cyan,label=false)
        plot!(xbins[2:end],bin_two_ars(ar2,ar1,xbins,nm.median),lw=2,c=:cyan,ls=:dash)
        return p_out
end

#-------------------------------------------
#"Functions for plotting random forest models and their PDPs"
#-------------------------------------------

function partial_dependence_2d(X, Y, feature1_idx, feature2_idx, grid_points=50)
	#"Function for getting 2-dimensional partial dependence"
        #"X must be matrix of all variables for training! Not only for plotting"
	model = fr_model(X,Y)
        feature1_vals = range(minimum(X[:, feature1_idx]), 
			       stop=maximum(X[:, feature1_idx]), length=grid_points)
        feature2_vals = range(minimum(X[:, feature2_idx]), 
			       stop=maximum(X[:, feature2_idx]), length=grid_points)

        pdp_vals = zeros(grid_points, grid_points)

        for i in 1:grid_points
                for j in 1:grid_points
                        X_copy = copy(X)
                        X_copy[:, feature1_idx] .= feature1_vals[i]
                        X_copy[:, feature2_idx] .= feature2_vals[j]
                        pdp_vals[i, j] = mean(DecisionTree.predict(model, X_copy))
                end
        end

    return feature1_vals, feature2_vals, pdp_vals
end

#"Function for plotting 2-dimension partial dependence as heatmap"
function pdp_2d_heatmap(X,Y,feature1_idx, feature2_idx,xlabel="x",ylabel="y",ctit="Y")
	x_vals, y_vals, pdp_2d = partial_dependence_2d(X,Y, feature1_idx, feature2_idx)
	p_out = heatmap(x_vals, y_vals, pdp_2d', xlabel=xlabel,ylabel=ylabel,
			colorbar_title=ctit)
	return p_out
end

#"Function for plotting individual conditional expectation plot (not used)"
function ice_plot(X, Y, feature_idx, grid_points=50)
        model = fr_model(X,Y)

        feature_vals = range(minimum(X[:, feature_idx]), 
			      stop=maximum(X[:, feature_idx]), length=grid_points)
        ice_vals = zeros(size(X, 1), grid_points)

        for i in 1:size(X, 1)
                X_copy = copy(X)

                for (j, val) in enumerate(feature_vals)
                        X_copy[i, feature_idx] = val
                        ice_vals[i, j] = DecisionTree.predict(model, X_copy)[i]
                end
        end

        plot(feature_vals, ice_vals', alpha=0.3, 
			   legend=false, title="ICE Plot for Feature $(feature_idx)",
             xlabel="Feature $(feature_idx)", ylabel="Predicted Value", color=:blue)

        avg_ice = mean(ice_vals, dims=1)
        plot!(feature_vals, avg_ice[:], linewidth=3, color=:red, label="Average ICE (PDP)")
end

function plot_PDP_shaded(model, X, Y, feature_idx, xlabel, ylabel, grid_points=50)
        feature_vals = X[:, feature_idx]
        feature_fine = range(minimum(feature_vals), stop=maximum(feature_vals), length=grid_points)

        pdp_vals = zeros(grid_points)
        pdp_95 = zeros(grid_points)
        pdp_05 = zeros(grid_points)

        predicted_vals = DecisionTree.predict(model,X)
        Yp = predicted_vals
        for (i, val) in enumerate(feature_fine)
                X_copy = copy(X)
                X_copy[:, feature_idx] .= val
                out = DecisionTree.predict(model, X_copy)
                pdp_vals[i] = mean(out)
                pdp_95[i] = prcntl_1d(out,95)
                pdp_05[i] = prcntl_1d(out,5)
        end

        p1 = plot(feature_fine, pdp_vals,xlabel=xlabel,ylabel="Pred. $(ylabel)",label=false,lw=1,
                           c=:black,
                           ylims=(minimum(Y),maximum(Y)),
                           fillrange=(pdp_95,pdp_05),fillalpha=0.3)

        return p1
end

function scatter_forest(X,x_list,Y,vr::String,trees=100,depth=3,mss=10)
        fr = fr_model(X,Y,trees,depth,mss)
        y_pred = DecisionTree.predict(fr, X)
        y_true = Y

        ss_total = sum((y_true .- mean(y_true)) .^ 2)
        ss_residual = sum((y_true - y_pred) .^ 2)
        r_squared = round(1 - ss_residual / ss_total,digits=3)
        println("R^2 = $(r_squared)")
        RMSE = rmsd(y_true,y_pred)
        NRMSE = round(RMSE / std(y_true),digits=3)
        print("NRMSE = $(NRMSE)")
        p1 = scatter(y_true,y_pred,
                      xlabel=vr,ylabel="Pred. $(vr)",label=false,c=:cyan)

        plot!(p1,Y,Y,label=false,lw=1,c=:black,title="a)",title_location=:left,
               titlefontsize=9,labelfontsize=8)

        annotate!(p1,mean(Y), maximum(Y),
                              L"$\mathrm{NRMSE}$"*" = $(NRMSE) \n "*L"$R^2$"*" = $(r_squared)",
                              annotationfontsize=:8)
        pdps = Dict()
        alph = ["b","c","d","e","f"]
        for i in 1:length(x_list)
                pdps[i] = plot_PDP_shaded(fr,X,Y,i,x_list[i],vr)
                plot!(pdps[i],title="$(alph[i]))",title_location=:left,titlefontsize=9,
                       labelfontsize=7)
                if x_list[i] == "SAT"
                        plot!(pdps[i],xrotation=45)
                end
        end

        if length(x_list) == 1
                lay = @layout [a b]
                p_out = plot(p1,pdps[1],layout=lay)
        elseif length(x_list) == 2
                lay = @layout [a [b
                                   c]]
                p_out = plot(p1,pdps[1],pdps[2],layout=lay)
        elseif length(x_list) == 3
                lay = @layout [a [b c
                                    d]]
                p_out = plot(p1,pdps[1],pdps[2],pdps[3],layout=lay)
        elseif length(x_list) == 4
                lay = @layout [a{0.35w} [b c
                                  d e]]
                p_out = plot(p1,pdps[1],pdps[2],pdps[3],pdps[4],layout=lay)
        end
        plot!(p_out,tickfontsize=7)
        return p_out
end

#-------------------------------------------------------
#"General plotting functions to be used in plots.jl"
#-------------------------------------------------------

#Function for plotting the spread in values across simulations
function plot_spread_vr(vr::String,nm="normal",factor=1)
        ds = dsd[nm]
        vr_ar = ds[vr][1,:] .* factor
        sp = sortperm(vr_ar)
        vr_med = ds[vr][2,:][sp] .* factor

        frange1 = (ds[vr][4,:][sp] .* factor,ds[vr][end,:][sp] .* factor)
        frange2 = (ds[vr][5,:][sp] .* factor,ds[vr][6,:][sp] .* factor)
        if vr == "t"
            vr_ar = vr_ar .- 273.15
            vr_med = vr_med .- 273.15
            frange1 = (ds[vr][4,:][sp] .- 273.15, ds[vr][end,:][sp] .- 273.15)
            frange2 = (ds[vr][5,:][sp] .- 273.15, ds[vr][6,:][sp] .-273.15)
        end
        p_out = plot(1:70,vr_ar[sp],fillrange=frange1,
                     color=:black,fillalpha=0.2,label=false,
                    ylabel=vr,xlabel="Simulation #",
                    xticks=(1:5:70,collect(1:70)[sp][1:5:end]))
        
        plot!(p_out,1:70, vr_med,
              fillrange=frange2,
              color=:black,fillalpha=0.3,label=false,
              ls=:dash, lw=1)
        return p_out
end


#"Function for plotting spearman correlation of variables in vr_list with all parameters"
function plot_cor_vr(vr_list,ind::Int,nm::String)
	ds = dsd[nm]
	a = zeros(length(vr_list),5)
	for i in 1:length(vr_list)
		if nm != "hours"
			a[i,:] = [corspearman(ds[vr_list[i]][ind,:],dp[p]) for p in params_list]
		elseif nm == "hours"
			a[i,:] = [corspearman(ds[vr_list[i]][ind,:,end],dp[p]) for p in params_list]
		end
	end
	addon=", all trajs"
	if nm == "masked"
		addon=", t600 < 10h"
	elseif nm =="hours"
		addon=", after 10h"
	end
	p_out = heatmap(a,clims=(-1,1),c=:vik,title=labeldict[ind]*addon,
			xticks=(1:5,params_list),yticks=(1:length(vr_list),vr_list),
			colorbartitle="\n corspearman",right_margin=5mm)
	return p_out
	GC.gc()
end

#"plot of spearman cors between variables in vr_list and parameters in param_list" 
function plot_cor_param_vr(param_list,vr_list,ind::Int,nm::String)
	ds = dsd[nm]
	a = zeros(length(vr_list),length(param_list))

	for i in 1:length(vr_list)
		if nm != "hours"
			a[i,:] = [corspearman(ds[vr_list[i]][ind,:],dp[p]) 
				   for p in param_list]
		elseif nm == "hours"
			a[i,:] = [corspearman(ds[vr_list[i]][ind,:,end],dp[p]) 
				   for p in param_list]
		end
	end
	addon=", all trajs"
	if nm == "masked"
		addon=", t600 < 10h"
	elseif nm =="hours"
		addon=", after 10h"
	end
	a = a'
	p_out = heatmap(a,clims=(-1,1),c=:vik,
			yticks=(1:length(param_list),param_list),xticks=(1:length(vr_list),vr_list),
			colorbartitle="\n corspearman",right_margin=5mm)
	return p_out
	GC.gc()
end


#"scatter plot of variable as a function of parameter pert. with shaded area"
function scat_vr_param_shaded(param::String,vr::String,ind::Int,factor,
			       units="g/kg",nm="normal",msk=false,blankp=true,
                               yscale=:identity,ylims=:auto,fitx=false)
	#nm must be "normal", "masked" or "hours"
	ds = dsd[nm]
	ar1 = dp[param]
        offset = 0
        if vr in ["t","temp","T_gl50","T_gl99"]
            offset = -273.15
        end

	if nm != "hours"
                ar2 = (ds[vr][ind,:] .* factor) .+ offset
                f95 = (ds[vr][7,:] .* factor) .+ offset
		f90 = (ds[vr][6,:] .* factor) .+ offset
		f10 = (ds[vr][5,:] .* factor) .+ offset
		f05 = (ds[vr][4,:] .* factor) .+ offset
	elseif nm =="hours"
                ar2 = (ds[vr][ind,:,end] .* factor) .+ offset
                f95 = (ds[vr][7,:,end] .* factor) .+ offset
		f90 = (ds[vr][6,:,end] .* factor) .+ offset
		f10 = (ds[vr][5,:,end] .* factor) .+ offset
		f05 = (ds[vr][4,:,end] .* factor) .+ offset
	end
	
	if msk != false
		ar1 = ar1[msk]
		ar2 = ar2[msk]
		f95 = f95[msk]
		f90 = f90[msk]
		f10 = f10[msk]
		f05 = f05[msk]
	end

	data = DataFrame(X=ar1, Y=ar2)

	sp = sortperm(data.X)

	if blankp == true
		bp = blankp_5_95()
	end

	p_out1 = plot(data.X[sp],f95[sp],fillrange=(f05[sp],f95[sp]),
		      alpha=0.5,c=:grey,label=false,lw=0,yscale=yscale,ylims=ylims)

	plot!(p_out1,data.X[sp],f90[sp],fillrange=(f10[sp],f90[sp]),
	       alpha=0.8,c=:grey,label=false,lw=0,yscale=yscale,ylims=ylims)

	scatter!(p_out1,data.X,data.Y,label=false,ylabel="$(vr) [$(units)]",
		  xlabel="$(param) pert.",c=:grey,markerstrokewidth=0.5,markersize=3,
                 yscale=yscale,ylims=ylims)

        modd = Dict()
        if fitx==false
            ols = lm(@formula(Y ~ X), data)
            model1(x) = coef(ols)[1] + coef(ols)[2] * x
            modd[false] = model1
            R2 = round(r2(ols),digits=3)
            c = round(cor(data.X,data.Y),digits=3)
            cs = round(corspearman(data.X,data.Y),digits=3)
        elseif fitx==true
            ols = lm(@formula(Y ~ 1/X),data)
            model2(x) = coef(ols)[1] + coef(ols)[2] / x
            modd[true] = model2
            R2 = round(r2(ols),digits=3)
            c = round(cor(data.X,data.Y),digits=3)
            cs = round(corspearman(data.X,data.Y),digits=3)
        end
        
        #model = modd[fitx]
        #plot!(p_out1,sort(data.X),model.(sort(data.X)),c=:black,label=false,lw=2,
	#       bottom_margin=5mm,left_margin=3mm,
	#       right_margin=3mm,yscale=yscale,ylims=ylims)

	if blankp == true
		lay = @layout [a{0.05h}
				b]		
		p_out = plot(bp,p_out1,layout=lay)
	else 
		p_out = p_out1
	end
	return p_out
end

#"scatter plot of variable as a function of parameter pert. with shaded area"
#"and masked according to msk1 and msk2"
function scat_vr_param_shaded_2msks(param::String,vr::String,ind::Int,factor,msk1,msk2,
				     label1,label2,units="g/kg",nm="normal",blankp=true)
	ds = dsd[nm]
	ar1 = dp[param]
	
	if nm != "hours"
		ar2 = ds[vr][ind,:] .* factor
		f95 = ds[vr][7,:] .* factor
		f90 = ds[vr][6,:] .* factor
		f10 = ds[vr][5,:] .* factor
		f05 = ds[vr][4,:] .* factor
	elseif nm =="hours"
		ar2 = ds[vr][ind,:,end] .* factor
		f95 = ds[vr][7,:,end] .* factor
		f90 = ds[vr][6,:,end] .* factor
		f10 = ds[vr][5,:,end] .* factor
		f05 = ds[vr][4,:,end] .* factor
	end
	
	ar1_1 = ar1[msk1]
	ar2_1 = ar2[msk1]

	ar1_2 = ar1[msk2]
	ar2_2 = ar2[msk2]

	data = DataFrame(X=ar1,Y=ar2)
	data_1 = DataFrame(X=ar1_1, Y=ar2_1)
	data_2 = DataFrame(X=ar1_2, Y=ar2_2)
       
        data_fit1 = data_1
        data_fit2 = data_2

        if param == "INP"
            data_fit1 = DataFrame(X=log10.(ar1_1),Y=ar2_1)
            data_fit2 = DataFrame(X=log10.(ar1_2),Y=ar2_2)
        end

	sp = sortperm(data.X)
	sp_1 = sortperm(data_1.X)
	sp_2 = sortperm(data_2.X)

	if blankp == true
		bp = blankp_5_95()
	end

	p_out1 = plot(data.X[sp],f95[sp],fillrange=(f05[sp],f95[sp]),
		      alpha=0.4,c=:grey,label=false,lw=0)

	plot!(p_out1,data.X[sp],f90[sp],fillrange=(f10[sp],f90[sp]),
	       alpha=0.8,c=:grey,label=false,lw=0)

	scatter!(p_out1,data_1.X,data_1.Y,label=label1,ylabel="$(vr) [$(units)]",
		  xlabel="$(param) pert.",c=:red,markerstrokewidth=0.5,markersize=3,
		  background_color_legend=nothing)
	scatter!(p_out1,data_2.X,data_2.Y,label=label2,c=:cyan,markerstrokewidth=0.5,markersize=3)

	ols_1 = lm(@formula(Y ~ X), data_fit1)
	model_1(x) = coef(ols_1)[1] + coef(ols_1)[2] * x
	R2_1 = round(r2(ols_1),digits=3)
	c_1 = round(cor(data_fit1.X,data_fit1.Y),digits=3)
	cs_1 = round(corspearman(data_fit1.X,data_fit1.Y),digits=3)
       
        println("----")
        println("$(param), $(vr)")
        println("spearman for red: $(cs_1)") 

        if param != "INP"
            plot!(p_out1,sort(data_1.X),model_1.(sort(data_1.X)),c=:red,label=false,lw=1,
                   bottom_margin=5mm,left_margin=3mm,
                   right_margin=3mm)
        elseif param == "INP"
            plot!(p_out1,sort(data_1.X),model_1.(sort(data_fit1.X)),c=:red,
                  label=false,lw=1,bottom_margin=5mm,left_margin=3mm,right_margin=3mm)
        end

	ols_2 = lm(@formula(Y ~ X), data_fit2)
	model_2(x) = coef(ols_2)[1] + coef(ols_2)[2] * x
	R2_2 = round(r2(ols_2),digits=3)
	c_2 = round(cor(data_fit2.X,data_fit2.Y),digits=3)
	cs_2 = round(corspearman(data_fit2.X,data_fit2.Y),digits=3)
        
        println("spearman for cyan: $(cs_2)") 
        println("----")
        if param != "INP"
            plot!(p_out1,sort(data_2.X),model_2.(sort(data_2.X)),c=:cyan,label=false,lw=1,
                   bottom_margin=5mm,left_margin=3mm,
                   right_margin=3mm)
        elseif param == "INP"
             plot!(p_out1,sort(data_2.X),model_2.(sort(data_fit2.X)),c=:cyan,label=false,lw=1,
                   bottom_margin=5mm,left_margin=3mm,right_margin=3mm)
        end

	if blankp == true
		lay = @layout [a{0.05h}
				b]		
		p_out = plot(bp,p_out1,layout=lay)
	else 
		p_out = p_out1
	end
	return p_out
end

function scat_vr_param_shaded_2msks_diff(param::String,vr::String,ind::Int,factor,msk1,msk2,
				     label1,label2,units="g/kg",nm="normal",blankp=true)
	ds = dsd[nm]
	ar1 = dp[param]
	
	if nm != "hours"
		ar2 = ds[vr][ind,:] .* factor
		f95 = ds[vr][7,:] .* factor
		f90 = ds[vr][6,:] .* factor
		f10 = ds[vr][5,:] .* factor
		f05 = ds[vr][4,:] .* factor
	elseif nm =="hours"
		ar2 = ds[vr][ind,:,end] .* factor
		f95 = ds[vr][7,:,end] .* factor
		f90 = ds[vr][6,:,end] .* factor
		f10 = ds[vr][5,:,end] .* factor
		f05 = ds[vr][4,:,end] .* factor
	end
	
	ar1_1 = ar1[msk1]
	ar2_1 = ar2[msk1]

	ar1_2 = ar1[msk2]
	ar2_2 = ar2[msk2]

	data = DataFrame(X=ar1,Y=ar2)
	data_1 = DataFrame(X=ar1_1, Y=ar2_1)
	data_2 = DataFrame(X=ar1_2, Y=ar2_2)
	
	sp = sortperm(data.X)
	sp_1 = sortperm(data_1.X)
	sp_2 = sortperm(data_2.X)

	if blankp == true
		bp = blankp_5_95()
	end

        p_out1 = plot(ar1_1[sp_1],f95[msk1][sp_1],
                      fillrange=(f05[msk1][sp_1],f95[msk1][sp_1]),
		      alpha=0.3,c=:red,label=false,lw=0)
        plot!(p_out1,ar1_2[sp_2],f95[msk2][sp_2],
                      fillrange=(f05[msk2][sp_2],f95[msk2][sp_2]),
		      alpha=0.2,c=:cyan,label=false,lw=0)

        plot!(p_out1,ar1_1[sp_1],f90[msk1][sp_1],fillrange=(f10[msk1][sp_1],f90[msk1][sp_1]),
	       alpha=0.4,c=:red,label=false,lw=0)
        plot!(p_out1,ar1_2[sp_2],f90[msk2][sp_2],fillrange=(f10[msk2][sp_2],f90[msk2][sp_2]),
	       alpha=0.4,c=:cyan,label=false,lw=0)

	scatter!(p_out1,data_1.X,data_1.Y,label=label1,ylabel="$(vr) [$(units)]",
		  xlabel="$(param) pert.",c=:red,markerstrokewidth=0.5,markersize=3,
		  background_color_legend=nothing)


	scatter!(p_out1,data_2.X,data_2.Y,label=label2,c=:cyan,markerstrokewidth=0.5,markersize=3)

	ols_1 = lm(@formula(Y ~ X), data_1)
	model_1(x) = coef(ols_1)[1] + coef(ols_1)[2] * x
	R2_1 = round(r2(ols_1),digits=3)
	c_1 = round(cor(data_1.X,data_1.Y),digits=3)
	cs_1 = round(corspearman(data_1.X,data_1.Y),digits=3)
	
        plot!(p_out1,sort(data_1.X),model_1.(sort(data_1.X)),c=:red,label=false,lw=1,
	       bottom_margin=5mm,left_margin=3mm,
	       right_margin=3mm)

	ols_2 = lm(@formula(Y ~ X), data_2)
	model_2(x) = coef(ols_2)[1] + coef(ols_2)[2] * x
	R2_2 = round(r2(ols_2),digits=3)
	c_2 = round(cor(data_2.X,data_2.Y),digits=3)
	cs_2 = round(corspearman(data_2.X,data_2.Y),digits=3)
	
        plot!(p_out1,sort(data_2.X),model_2.(sort(data_2.X)),c=:cyan,label=false,lw=1,
	       bottom_margin=5mm,left_margin=3mm,
	       right_margin=3mm)

	if blankp == true
		lay = @layout [a{0.05h}
				b]		
		p_out = plot(bp,p_out1,layout=lay)
	else 
		p_out = p_out1
	end
	return p_out
end

#"plot the legend for shaded areas as stand-alone plot"
function blankp_5_95()
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
                               yticks=false,xlim=(0,1),ylim=(0,1))

	plot!(blankp,[],c=:grey,fillrange=([NaN],[NaN]),fillalpha=0.8,label="10-90th",lw=0,
	       legend=:top,legend_columns=2)
	plot!(blankp,[],c=:grey,fillrange=([NaN],[NaN]),fillalpha=0.4,label="5-95th",lw=0)
	return blankp
end

#"plot a blank plot"
function blankp_clean()
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
                               yticks=false,xlim=(0,1),ylim=(0,1))
	return blankp
end

#"Function for creating a standalone colorbar for the spearman correlation coefficient"
function spear_cbar()
        cbar = heatmap((-1:0.02:1)' .* ones(101,1), legend=:none,labelfontsize=8,
                                    yticks=:none,xticks=(1:10:101, string.(-1:0.2:1)),
                                    xlabel="Spearman correlation",c=:vik)
        return cbar
end
