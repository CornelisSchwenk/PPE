#---------------------------------------------
""" This script imports the main modules and function from other scripts to enable
the calculations and analyses preformed for this Paper.
"""
#---------------------------------------------

include("use_config.jl")
include("thermodynamic_functions.jl")
include("statistical_functions.jl")

using NPZ
using DelimitedFiles
using GLM
using NCDatasets
using Plots
using Plots.PlotMeasures
using Statistics
using DataStructures
using NaNMath; nm=NaNMath
using LaTeXStrings
using ColorSchemeTools
using DataFrames
using StatsBase

#path to data
data_path = "../DATA/"

#hydrometeor characters
moist = ["v","c","r","i","s","g"]

#parameter values
params = readdlm("../DATA/params.txt",',', Float64, '\n',skipstart=1)
params_list = ["SST", "CCN", "INP", "CAP", "SATAD"]
dp = Dict()
for i in 1:length(params_list)
	dp[params_list[i]] = params[:,i]
end

#Load netCDF files with relevant variables.
#Files are only loaded if they exist. If a file is written, then include this script again. 
#All files are written in DO_FULL_ANALYSIS.jl
dsd = Dict()
if isfile(data_path*"PPE_vars_all.nc")
    dsd["normal"] = Dataset(data_path*"PPE_vars_all.nc")
end
if isfile(data_path*"PPE_vars_end_hours.nc")
    dsd["hours"] = Dataset(data_path*"PPE_vars_end_hours.nc")
end
if isfile(data_path*"PPE_vars_accum.nc")
    dsd["accum"] = Dataset(data_path*"PPE_vars_accum.nc")
end
if isfile(data_path*"PPE_vars_accum_hours.nc")
    dsd["accum_hours"] = Dataset(data_path*"PPE_vars_accum_hours.nc")
end
if isfile(data_path*"PPE_vars_t600.nc")
    dsd["conv"] = Dataset(data_path*"PPE_vars_t600.nc")
end

#Function to load variable for a PPE member
function load_vr(vr::String,mem::Int,msk=false)
	ar = npzread(data_path*"sens_exp_ID0$(mem)/WCB_ar_$(vr).npy")
	upto = findfirst(vec(sum((!isnan).(ar),dims=2)) .== 0) -1
        if msk != false
            return ar[collect(1:upto)[msk],:]
        else
	    return ar[1:upto,:]
        end
end

function find_tau_x00(p,x)
        #"This function takes a 1 dim pressure array in hPa and returns ascent time"
        #"p must be in hPa and a 1 dim array"
        DT = Inf
        t = []
        for i in 1:(length(p)-1)
                for j in i:length(p)
                       dp = p[i] - p[j]
                       dt = j-i
                       if dp > x
                               if dt < DT
                                       push!(t,dt)
                                       DT=dt
                               end
                               break
                       end
               end
        end
        return t[end]
end

#function to calculate all tau_600 values for the WCB trajectories of a PPE member
function calc_t600(mem::Int)
	pres = load_vr("p",mem)
	return [find_tau_x00(pres[i,:],600) for i in 1:size(pres)[1]]
end

#function to determine end of ascent
function find_stop(p_ar)
	p = p_ar[61:end]
	st = sum((!isnan).(p))
	t6 = find_tau_x00(p,600)
	if t6+1 == st
		return st+60
	end
	for i in t6:st
		p_diff_min = 8
		if i == st
			return st+60
		elseif p[i] - p[i+1] < p_diff_min
			if i+60 < 61 + t6
				return 61+t6
			else
				return i+60
			end
		end
	end
end

#function to determine the beginning of ascent
function find_start(p_ar)
	p = p_ar[1:61]
	st = sum((isnan).(p))
	if st == 60
		return 61
	end
	for i in 61:-1:st
		p_diff_min = 8
		if i == st+1
			return st+1	
		elseif p[i-1] - p[i] < p_diff_min
			return i
		end
	end
end

#function to get indices for end of ascent plus ts amount of timesteps
function get_tWCB_plus_ts(mem::Int,ts::Int,msk=false)
	pres = load_vr("p",mem,msk)
	stops = [find_stop(pres[i,:]) for i in 1:size(pres)[1]]
	firstnan = [findfirst((!isnan).(pres[i,:])) for i in 1:size(pres)[1]]
	lastnan = [findfirst((isnan).(pres[i,firstnan[i]:end])) for i in 1:size(pres)[1]]
	max_stop = firstnan .+ lastnan .- 2
	inds = [CartesianIndex(i,stops[i]+ts) for i in 1:size(pres)[1]]
	msk2 = stops .+ ts .< max_stop
	return inds,msk2
end

#function to load variable at the end of the tau_WCB ascent
function load_vr_end(vr::String,mem::Int,msk=false)
	pres = load_vr("p",mem,msk)
	stops = [find_stop(pres[i,:]) for i in 1:size(pres)[1]]	
	stop_inds = [CartesianIndex(i,stops[i]) for i in 1:size(pres)[1]]
	ar = load_vr(vr,mem,msk)
	GC.gc()
	return ar[stop_inds]
end

#function to load variable at the start of the tau_WCB ascent
function load_vr_start(vr::String,mem::Int,msk=false)
	pres = load_vr("p",mem,msk)
	starts = [find_start(pres[i,:]) for i in 1:size(pres)[1]]	
	start_inds = [CartesianIndex(i,starts[i]) for i in 1:size(pres)[1]]
	ar = load_vr(vr,mem,msk)
	GC.gc()
	return ar[start_inds]
end
		
#function to load variable ts timesteps after the end of the tau_WCB ascent
function load_vr_ts(vr::String,mem::Int,ts::Int,msk=false)
	v = load_vr(vr,mem,msk)
	inds,msk2=get_tWCB_plus_ts(mem,ts,msk)
	return v[inds[msk2]]
end

#function to do the same but keep NaN values if no timesteps after ts
function load_vr_ts_withnan(vr::String,mem::Int,ts::Int,msk=false)
	v = load_vr(vr,mem,msk)
        len = size(v)[1]
        out = [NaN for i in 1:len]
        loopinds = collect(1:len)

	inds,msk2=get_tWCB_plus_ts(mem,ts,msk)
        for k in loopinds[msk2]
            out[k] = v[inds[k]]
        end
        return out
end

#function to calculate the total water content for all trajectories in a PPE member
function calc_qtot_all(mem::Int,msk=false)
	qtot = 0
	for hy in moist
		qtot = qtot .+ load_vr("q$(hy)",mem,msk)
	end
	return qtot
end

#function to calculate the total water content at the end of the ascent for a PPE member
function calc_qtot_end(mem::Int,msk=false)
	qtot = 0
	for hy in moist
		qtot = qtot .+ load_vr_end("q$(hy)",mem,msk)
	end
	return qtot
end

#function to calculate the total water content at the start of the ascent for a PPE member
function calc_qtot_start(mem::Int,msk=false)
	qtot = 0
	for hy in moist
		qtot = qtot .+ load_vr_start("q$(hy)",mem,msk)
	end
	return qtot
end

#function to calculate the total hydrometeor content for all trajectories in a PPE member
function calc_qhyd_all(mem::Int,msk=false)
	qhyd = 0
	for hy in moist[2:end]
		qhyd = qhyd .+ load_vr("q$(hy)",mem,msk)
	end
	return qhyd
end

#function to calculate total hydrometeor number mixing ratio for all trajectories in a PPE member
function calc_qnhyd_all(mem::Int,msk=false)
	qnhyd = 0
	for hy in ["c","i","s"]
		qnhyd = qnhyd .+ load_vr("qn$(hy)",mem,msk)
	end
	return qnhyd
end

#function to calculate total hydrometeor content at the end of the ascent for a PPE member
function calc_qhyd_end(mem::Int,msk=false)
	qhyd = 0
	for hy in moist[2:end]
		qhyd = qhyd .+ load_vr_end("q$(hy)",mem,msk)
		GC.gc()
	end
	return qhyd
end

#function to calculate total hydrometeor content at the start of the ascent for a PPE member
function calc_qhyd_start(mem::Int,msk=false)
	qhyd = 0
	for hy in moist[2:end]
		qhyd = qhyd .+ load_vr_start("q$(hy)",mem,msk)
	end
	return qhyd
end

#function to calculate the maximum hydrometeor content encountered during the ascent
#for each trajectory in a PPE member
function calc_max_hyd(mem::Int,msk=false)
	H = calc_qhyd_all(mem,msk)
	return [nm.maximum(H[i,:]) for i in 1:size(H)[1]]
end

#function to calculate the maximum hydrometeor number mixing ratio encountered during the ascent
#for each trajectory in a PPE member
function calc_max_nhyd(mem::Int,msk=false)
	H = calc_qnhyd_all(mem,msk)
	return [nm.maximum(H[i,:]) for i in 1:size(H)[1]]
end

#function to calculate the maximum value of a variable during the ascent for each trajectory
#in a PPE member
function calc_max_var(vr::String,mem::Int,msk=false)
	vr_ar = load_vr(vr,mem,msk)
	return [nm.maximum(vr_ar[i,:]) for i in 1:size(vr_ar)[1]]
end

#function to calculate the drying ratio for all trajectories in a PPE member
function calc_DR(mem::Int,msk=false)
	qtot0 = calc_qtot_start(mem,msk)
	qtot1 = calc_qtot_end(mem,msk)
	DR = (qtot0 .- qtot1) ./ qtot0
	GC.gc()
	return DR
end

#function to calculate the precipitation efficiency for all trajectories in a PPE member
function calc_PE(mem::Int,msk=false)
	max_H = calc_max_hyd(mem,msk)
	H_end = calc_qhyd_end(mem,msk)
	GC.gc()
	return 1 .- (H_end ./ max_H)
end

#get cloud number mixing ratio at cloud base
function get_qnc_cloud_base(mem)
    qnc = load_vr("qnc",mem)
    cb = zeros(size(qnc)[1])
    for i in 1:length(cb)
        cb[i] = qnc[i,findfirst(qnc[i,:] .> 0)]
    end
    return cb
end
