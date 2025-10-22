include("calcs.jl")
using Base.Threads

pth = "../DATA/"
ds_out_name = "PPE_vars_end_hours.nc"
ds_out = Dataset(pth*ds_out_name,"c")

defDim(ds_out,"mems",70) #ensemble members
defDim(ds_out,"stats",7) #mean, median, std, 5th, 10th, 90th 95th
defDim(ds_out,"hours",2) #5 and 10 hours after ascent

ds_out.attrib["title"] = "Statstics per ensemble member in the hours AFTER the ascent."
ds_out.attrib["description"] = "For physical variables the first dimension is ensemble member from 20 to 89, second dimension is statistics in this order: mean, median, std, 5th, 25th, 75th and 95th percentile. Third dimension is hours after the ascent as 2, 4, 6, 8 and 10 hours."

v_end = Dict()
v_start = Dict()

qv_out = defVar(ds_out,"qv",Float64,("stats","mems","hours"))
qv_out[:] = zeros(7,70,2)

p_out = defVar(ds_out,"p",Float64,("stats","mems","hours"))
p_out[:] = zeros(7,70,2)

T_out = defVar(ds_out,"t",Float64,("stats","mems","hours"))
T_out[:] = zeros(7,70,2)

RHi_out = defVar(ds_out,"RHi",Float64,("stats","mems","hours"))
RHi_out[:] = zeros(7,70,2)

RHw_out = defVar(ds_out,"RHw",Float64,("stats","mems","hours"))
RHw_out[:] = zeros(7,70,2)

alt_out = defVar(ds_out,"alt",Float64,("stats","mems","hours"))
alt_out[:] = zeros(7,70,2)

pv_out = defVar(ds_out,"pv",Float64,("stats","mems","hours"))
pv_out[:] = zeros(7,70,2)

qi_out = defVar(ds_out,"qi",Float64,("stats","mems","hours"))
qi_out[:] = zeros(7,70,2)

qni_out = defVar(ds_out,"qni",Float64,("stats","mems","hours"))
qni_out[:] = zeros(7,70,2)


global lock_obj = ReentrantLock()
global lock_obj2 = ReentrantLock()

#This one is different than the others bc it accounts for 3dim
function process_array3d(input_array,output_array,i,j)
        local_result = zeros(7)

        sorted_array = sort(input_array)

        local_result[1] = mean(sorted_array)
        local_result[2] = median(sorted_array)
        local_result[3] = std(sorted_array)
        local_result[4] = sorted_array[Int(round(length(sorted_array) * 5 / 100, digits=0))]
        local_result[5] = sorted_array[Int(round(length(sorted_array) * 25 / 100, digits=0))]
        local_result[6] = sorted_array[Int(round(length(sorted_array) * 75 / 100, digits=0))]
        local_result[7] = sorted_array[Int(round(length(sorted_array) * 95 / 100, digits=0))]

        #lock(lock_obj)
        output_array[:,i,j] .= local_result
        #unlock(lock_obj)
end

#You CAN mutlithread this code, but I make no promises. If you want to, 
#then uncomment below and comment further below. Also uncomment the locks 
#in the loop and in the function above.

#@threads for i in 1:70
for i in 1:70
	m = i + 19
	print("$(m): ")	
        for (j,ts) in enumerate([5,10])
		print("$(ts), ")
		Tmu = load_vr_ts("t",m,ts)
		pmu = load_vr_ts("p",m,ts)
		qvmu = load_vr_ts("qv",m,ts)
		RHim = RH_i.(Tmu,qvmu,pmu .* 100)
		RHwm = RH_w.(Tmu,qvmu,pmu .* 100)
                altm = load_vr_ts("alt",m,ts)
                pvm = load_vr_ts("pv",m,ts)
                qim = load_vr_ts("qi",m,ts)
                qnim = load_vr_ts("qni",m,ts)


                process_array3d(Tmu, T_out, i,j)
                process_array3d(pmu, p_out, i,j)
                process_array3d(qvmu, qv_out, i,j)
                process_array3d(RHim, RHi_out, i,j)
                process_array3d(RHwm, RHw_out, i,j)
                process_array3d(altm, alt_out, i,j)
                process_array3d(pvm, pv_out, i,j)
                process_array3d(qim, qi_out, i,j)
                process_array3d(qnim, qni_out, i,j)
        end
end
println(" ")
println("done")
close(ds_out)
