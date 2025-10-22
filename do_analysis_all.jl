include("calcs.jl")
using Base.Threads

println("---USING NUMBER OF THREADS---")
println(Threads.nthreads())
println("-----------------------------")

pth = "../DATA/"
ds_out_name = "PPE_vars_all.nc"
ds_out = Dataset(pth*ds_out_name,"c")

defDim(ds_out,"mems",70) #ensemble members
defDim(ds_out,"stats",7) #mean, median, std, 5th, 10th, 90th 95th

ds_out.attrib["title"] = "Statstics per ensemble member FOR ALL TRAJECTORIES."
ds_out.attrib["description"] = "For physical variables the first dimension is ensemble member from 20 to 89, second dimension is statistics in this order: mean, median, std, 5th, 25th, 75th and 95th percentile. For ntrajs its just one entry per ensemble member"

#---------------------------------------------------------
#       DEFINE VARIABLES FOR NETCDF FILE
#---------------------------------------------------------

#--------------------------------------------
# Variables at the start and end of the ascent
#--------------------------------------------

alt_start_out = defVar(ds_out,"alt_start",Float64,("stats","mems"))
alt_start_out[:] = zeros(7,70)
pv_start_out = defVar(ds_out,"pv_start",Float64,("stats","mems"))
pv_start_out[:] = zeros(7,70)
qv_start_out = defVar(ds_out,"qv_start",Float64,("stats","mems"))
qv_start_out[:] = zeros(7,70)
qc_start_out = defVar(ds_out,"qc_start",Float64,("stats","mems"))
qc_start_out[:] = zeros(7,70)
qr_start_out = defVar(ds_out,"qr_start",Float64,("stats","mems"))
qr_start_out[:] = zeros(7,70)
qi_start_out = defVar(ds_out,"qi_start",Float64,("stats","mems"))
qi_start_out[:] = zeros(7,70)
qs_start_out = defVar(ds_out,"qs_start",Float64,("stats","mems"))
qs_start_out[:] = zeros(7,70)
qg_start_out = defVar(ds_out,"qg_start",Float64,("stats","mems"))
qg_start_out[:] = zeros(7,70)
qnc_start_out = defVar(ds_out,"qnc_start",Float64,("stats","mems"))
qnc_start_out[:] = zeros(7,70)
qni_start_out = defVar(ds_out,"qni_start",Float64,("stats","mems"))
qni_start_out[:] = zeros(7,70)
qns_start_out = defVar(ds_out,"qns_start",Float64,("stats","mems"))
qns_start_out[:] = zeros(7,70)

p_start_out = defVar(ds_out,"p_start",Float64,("stats","mems"))
p_start_out[:] = zeros(7,70)
T_start_out = defVar(ds_out,"T_start",Float64,("stats","mems"))
T_start_out[:] = zeros(7,70)
Theta_start_out = defVar(ds_out,"Theta_start",Float64,("stats","mems"))
Theta_start_out[:] = zeros(7,70)
Theta_end_out = defVar(ds_out,"Theta_end",Float64,("stats","mems"))
Theta_end_out[:] = zeros(7,70)
d_Theta_out = defVar(ds_out,"d_Theta",Float64,("stats","mems"))
d_Theta_out[:] = zeros(7,70)
d_alt_out = defVar(ds_out,"d_alt",Float64,("stats","mems"))
d_alt_out[:] = zeros(7,70)

#Omit qv here, done further down
alt_end_out = defVar(ds_out,"alt",Float64,("stats","mems"))
alt_end_out[:] = zeros(7,70)
pv_end_out = defVar(ds_out,"pv",Float64,("stats","mems"))
pv_end_out[:] = zeros(7,70)
qc_end_out = defVar(ds_out,"qc",Float64,("stats","mems"))
qc_end_out[:] = zeros(7,70)
qr_end_out = defVar(ds_out,"qr",Float64,("stats","mems"))
qr_end_out[:] = zeros(7,70)
qi_end_out = defVar(ds_out,"qi",Float64,("stats","mems"))
qi_end_out[:] = zeros(7,70)
qs_end_out = defVar(ds_out,"qs",Float64,("stats","mems"))
qs_end_out[:] = zeros(7,70)
qg_end_out = defVar(ds_out,"qg",Float64,("stats","mems"))
qg_end_out[:] = zeros(7,70)
qnc_end_out = defVar(ds_out,"qnc",Float64,("stats","mems"))
qnc_end_out[:] = zeros(7,70)
qni_end_out = defVar(ds_out,"qni",Float64,("stats","mems"))
qni_end_out[:] = zeros(7,70)
qns_end_out = defVar(ds_out,"qns",Float64,("stats","mems"))
qns_end_out[:] = zeros(7,70)

#--------------------------------------------
# Number of trajectories in sim
#--------------------------------------------
ntrajs_out = defVar(ds_out,"ntrajs",Int,("mems",))
ntrajs_out[:] = Int.(zeros(70))

#--------------------------------------------
#       Variables at 95% glaciation
#--------------------------------------------
pgl05_out = defVar(ds_out,"p_gl05",Float64,("stats","mems"))
pgl05_out[:] = zeros(7,70)

tgl05_out = defVar(ds_out,"T_gl05",Float64,("stats","mems"))
tgl05_out[:] = zeros(7,70)

#--------------------------------------------
#       Variables at 50% glaciation
#--------------------------------------------
pgl50_out = defVar(ds_out,"p_gl50",Float64,("stats","mems"))
pgl50_out[:] = zeros(7,70)

tgl50_out = defVar(ds_out,"T_gl50",Float64,("stats","mems"))
tgl50_out[:] = zeros(7,70)

#--------------------------------------------
# Maximum variables during ascent (per trajectory)
#--------------------------------------------
MH_out =  defVar(ds_out,"maxH",Float64,("stats","mems"))
MH_out[:] = zeros(7,70)

MnH_out =  defVar(ds_out,"maxnH",Float64,("stats","mems"))
MnH_out[:] = zeros(7,70)

max_qc_out = defVar(ds_out,"max_qc",Float64,("stats","mems"))
max_qc_out[:] = zeros(7,70)

max_qr_out = defVar(ds_out,"max_qr",Float64,("stats","mems"))
max_qr_out[:] = zeros(7,70)

max_qi_out = defVar(ds_out,"max_qi",Float64,("stats","mems"))
max_qi_out[:] = zeros(7,70)

max_qs_out = defVar(ds_out,"max_qs",Float64,("stats","mems"))
max_qs_out[:] = zeros(7,70)

max_qg_out = defVar(ds_out,"max_qg",Float64,("stats","mems"))
max_qg_out[:] = zeros(7,70)

max_qnc_out = defVar(ds_out,"max_qnc",Float64,("stats","mems"))
max_qnc_out[:] = zeros(7,70)

max_qni_out = defVar(ds_out,"max_qni",Float64,("stats","mems"))
max_qni_out[:] = zeros(7,70)

max_qns_out = defVar(ds_out,"max_qns",Float64,("stats","mems"))
max_qns_out[:] = zeros(7,70)

max_SLW_out = defVar(ds_out,"max_SLW",Float64,("stats","mems"))
max_SLW_out[:] = zeros(7,70)

max_qfr_out = defVar(ds_out,"max_qfr",Float64,("stats","mems"))
max_qfr_out[:] = zeros(7,70)

#--------------------------------------------
# Radius and relaxation timescale
#--------------------------------------------
rad_out = defVar(ds_out,"r",Float64,("stats","mems"))
rad_out[:] = zeros(7,70)

Ncr_out = defVar(ds_out,"Ncr",Float64,("stats","mems"))
Ncr_out[:] = zeros(7,70)

tau_out = defVar(ds_out,"tau",Float64,("stats","mems"))
tau_out[:] = zeros(7,70)

#--------------------------------------------
#       Variables at 99.9% glaciation
#--------------------------------------------
pgl99_out = defVar(ds_out,"p_gl99",Float64,("stats","mems"))
pgl99_out[:] = zeros(7,70)

tgl99_out = defVar(ds_out,"T_gl99",Float64,("stats","mems"))
tgl99_out[:] = zeros(7,70)

qfr99_out = defVar(ds_out,"qfr99",Float64,("stats","mems"))
qfr99_out[:] = zeros(7,70)

qi99_out = defVar(ds_out,"qi99",Float64,("stats","mems"))
qi99_out[:] = zeros(7,70)

qni99_out = defVar(ds_out,"qni99",Float64,("stats","mems"))
qni99_out[:] = zeros(7,70)

ri99_out = defVar(ds_out,"ri99",Float64,("stats","mems"))
ri99_out[:] = zeros(7,70)

qs99_out = defVar(ds_out,"qs99",Float64,("stats","mems"))
qs99_out[:] = zeros(7,70)

qns99_out = defVar(ds_out,"qns99",Float64,("stats","mems"))
qns99_out[:] = zeros(7,70)

rs99_out = defVar(ds_out,"rs99",Float64,("stats","mems"))
rs99_out[:] = zeros(7,70)

qg99_out = defVar(ds_out,"qg99",Float64,("stats","mems"))
qg99_out[:] = zeros(7,70)

Ncr99_out = defVar(ds_out,"Ncr99",Float64,("stats","mems"))
Ncr99_out[:] = zeros(7,70)

qv99_out = defVar(ds_out,"qv99",Float64,("stats","mems"))
qv99_out[:] = zeros(7,70)

tau99_out = defVar(ds_out,"tau99",Float64,("stats","mems"))
tau99_out[:] = zeros(7,70)

RHi99_out = defVar(ds_out,"RHi99",Float64,("stats","mems"))
RHi99_out[:] = zeros(7,70)
#--------------------------------------------
# tau_600
#--------------------------------------------
t6_out = defVar(ds_out,"t600",Float64,("stats","mems"))
t6_out[:] = zeros(7,70)

#--------------------------------------------
# PE, DR, qv, p, T, RHi, RHw at the end of ascent
#--------------------------------------------
PE_out = defVar(ds_out,"PE",Float64,("stats","mems"))
PE_out[:] = zeros(7,70)

DR_out = defVar(ds_out,"DR",Float64,("stats","mems"))
DR_out[:] = zeros(7,70)

qv_out = defVar(ds_out,"qv",Float64,("stats","mems"))
qv_out[:] = zeros(7,70)

p_out = defVar(ds_out,"p",Float64,("stats","mems"))
p_out[:] = zeros(7,70)

T_out = defVar(ds_out,"t",Float64,("stats","mems"))
T_out[:] = zeros(7,70)

RHi_out = defVar(ds_out,"RHi",Float64,("stats","mems"))
RHi_out[:] = zeros(7,70)

RHw_out = defVar(ds_out,"RHw",Float64,("stats","mems"))
RHw_out[:] = zeros(7,70)

global lock_obj = ReentrantLock()
global lock_obj2 = ReentrantLock()

function process_array(input_array,output_array,i)
        local_result = zeros(7)
        
        sorted_array = sort(input_array)
        
        local_result[1] = nm.mean(sorted_array)
        local_result[2] = nm.median(sorted_array)
        local_result[3] = std(sorted_array)
        local_result[4] = sorted_array[Int(round(length(sorted_array) * 5 / 100, digits=0))]
        local_result[5] = sorted_array[Int(round(length(sorted_array) * 25 / 100, digits=0))]
        local_result[6] = sorted_array[Int(round(length(sorted_array) * 75 / 100, digits=0))]
        local_result[7] = sorted_array[Int(round(length(sorted_array) * 95 / 100, digits=0))]
        
        #lock(lock_obj)
        output_array[:,i] .= local_result
        #unlock(lock_obj)
end

#You CAN use threads if you want. If so, uncomment below and comment further below. Also uncomment the locks in the function above.
#@threads for i in 1:70

for i in 1:70
	m = i + 19
	print("$(m), ")

        t6 = Float64.(calc_t600(m))
        msk = t6 .< 100 #dummy mask, we are taking all trajs. Lazy coding. 

        lock(lock_obj2)
        ntrajs_out[i] = length(t6[msk])
        unlock(lock_obj2)
        Tmu = load_vr_end("t",m,msk)
        pmu = load_vr_end("p",m,msk)
        qvmu = load_vr_end("qv",m,msk)
 
	RHim = sort(RH_i.(Tmu,qvmu,pmu .* 100))
	RHwm = sort(RH_w.(Tmu,qvmu,pmu .* 100))
        
        PEm = sort(calc_PE(m,msk))
        DRm = sort(calc_DR(m,msk))
        MaxHm = sort(calc_max_hyd(m,msk))
        MaxnHm = sort(calc_max_nhyd(m,msk))
		
	process_array(t6, t6_out, i)
	process_array(Tmu, T_out, i)
	process_array(pmu, p_out, i)
	process_array(qvmu, qv_out, i)
	process_array(RHim, RHi_out, i)
	process_array(RHwm, RHw_out, i)
	process_array(PEm, PE_out, i)
	process_array(DRm, DR_out, i)
	process_array(MaxHm, MH_out, i)
	process_array(MaxnHm, MnH_out, i)
        
        alt_end = load_vr_end("alt",m,msk)
        pv_end  = load_vr_end("pv",m,msk)
        qc_end  = load_vr_end("qc",m,msk)
        qr_end  = load_vr_end("qr",m,msk)
        qi_end  = load_vr_end("qi",m,msk)
        qs_end  = load_vr_end("qs",m,msk)
        qg_end  = load_vr_end("qg",m,msk)
        qnc_end = load_vr_end("qnc",m,msk)
        qni_end = load_vr_end("qni",m,msk)
        qns_end = load_vr_end("qns",m,msk)

        process_array(alt_end, alt_end_out, i)
        process_array(pv_end, pv_end_out, i)
        process_array(qc_end, qc_end_out, i)
        process_array(qr_end, qr_end_out, i)
        process_array(qi_end, qi_end_out, i)
        process_array(qs_end, qs_end_out, i)
        process_array(qg_end, qg_end_out, i)
        process_array(qnc_end, qnc_end_out, i)
        process_array(qni_end, qni_end_out, i)
        process_array(qns_end, qns_end_out, i)
        
        alt_start = load_vr_start("alt",m,msk)
        pv_start  = load_vr_start("pv",m,msk)
        qv_start  = load_vr_start("qv",m,msk)
        qc_start  = load_vr_start("qc",m,msk)
        qr_start  = load_vr_start("qr",m,msk)
        qi_start  = load_vr_start("qi",m,msk)
        qs_start  = load_vr_start("qs",m,msk)
        qg_start  = load_vr_start("qg",m,msk)
        qnc_start = load_vr_start("qnc",m,msk)
        qni_start = load_vr_start("qni",m,msk)
        qns_start = load_vr_start("qns",m,msk)
        T_start   = load_vr_start("t",m,msk)
        p_start   = load_vr_start("p",m,msk)
        
        process_array(alt_start, alt_start_out, i)
        process_array(pv_start, pv_start_out, i)
        process_array(qv_start, qv_start_out, i)
        process_array(qc_start, qc_start_out, i)
        process_array(qr_start, qr_start_out, i)
        process_array(qi_start, qi_start_out, i)
        process_array(qs_start, qs_start_out, i)
        process_array(qg_start, qg_start_out, i)
        process_array(qnc_start, qnc_start_out, i)
        process_array(qni_start, qni_start_out, i)
        process_array(qns_start, qns_start_out, i)
        process_array(T_start, T_start_out, i)
        process_array(p_start, p_start_out, i)

        theta_end = Theta.(Tmu,pmu .* 100)
        theta_start = Theta.(T_start,p_start .* 100)
        dtheta = theta_end .- theta_start
        dalt = alt_end .- alt_start
        
        process_array(theta_end, Theta_end_out, i)
        process_array(theta_start, Theta_start_out, i)
        process_array(dtheta, d_Theta_out, i)
        process_array(dalt, d_alt_out, i)

        max_qc = sort(calc_max_var("qc",m,msk))
        max_qr = sort(calc_max_var("qr",m,msk))
        max_qi = sort(calc_max_var("qi",m,msk))
        max_qs = sort(calc_max_var("qs",m,msk))
        max_qg = sort(calc_max_var("qg",m,msk))
        
        process_array(max_qc,max_qc_out,i)
        process_array(max_qr,max_qr_out,i)
        process_array(max_qi,max_qi_out,i)
        process_array(max_qs,max_qs_out,i)
        process_array(max_qg,max_qg_out,i)

        max_qnc = sort(calc_max_var("qnc",m,msk))
        max_qni = sort(calc_max_var("qni",m,msk))
        max_qns = sort(calc_max_var("qns",m,msk))

        process_array(max_qnc,max_qnc_out,i)
        process_array(max_qni,max_qni_out,i)
        process_array(max_qns,max_qns_out,i)

        #-----------------------------------------

        p = load_vr("p",m,msk)
        T = load_vr("t",m,msk)
        qv = load_vr("qv",m,msk)
        qc = load_vr("qc",m,msk)
        qr = load_vr("qr",m,msk)
        qi = load_vr("qi",m,msk)
        qni = load_vr("qni",m,msk)
        qs = load_vr("qs",m,msk)
        qns = load_vr("qns",m,msk)
        qg = load_vr("qg",m,msk)

        qliq = qc .+ qr
        qfr = qi .+ qs .+ qg
        liq_frac = qliq ./ (qliq .+ qfr)

        qi_end = load_vr_end("qi",m,msk)
        qni_end = load_vr_end("qni",m,msk)

        tau = tau_rel(pmu,Tmu,qni_end,qi_end,dp["CAP"][i])
        
        r = calc_radius_ice(qi_end,qni_end)
        r_all = calc_radius_ice(qi,qni)
        Ncr = qni_end .* r

        SLW = qliq
        SLW[T .> 273.15] .= NaN

        r = r[(!).(isnan.(r))]
        Ncr = Ncr[(!).(isnan.(Ncr))]
        tau = tau[(!).(isnan.(tau))]
        max_SLW = [nm.maximum(SLW[k,:]) for k in 1:size(SLW)[1]]
        max_qfr = [nm.maximum(qfr[k,:]) for k in 1:size(qfr)[1]]

        process_array(r,rad_out,i)
        process_array(Ncr,Ncr_out,i)
        process_array(tau,tau_out,i)
        process_array(max_SLW,max_SLW_out,i)
        process_array(max_qfr,max_qfr_out,i)

        #-----------------------------------------

        sz = size(liq_frac)[1]
        p_gl05 = zeros(sz)
        T_gl05 = zeros(sz)
        p_gl50 = zeros(sz)
        T_gl50 = zeros(sz)
        p_gl99 = zeros(sz)
        T_gl99 = zeros(sz)
        p_gl99 = zeros(sz)
        T_gl99 = zeros(sz)
        qfr99 = zeros(sz)
        qi99 = zeros(sz)
        qni99 = zeros(sz)
        ri99 = zeros(sz)
        qs99 = zeros(sz)
        qns99 = zeros(sz)
        rs99 = zeros(sz)
        qg99 = zeros(sz)
        Ncr99 = zeros(sz)
        qv99 = zeros(sz)

        for j in 1:length(p_gl05)
                strt = find_start(p[j,:])
                if nm.minimum(liq_frac[j,strt:end]) > 0.05
                        p_gl05[j] = NaN
                        T_gl05[j] = NaN
                else
                        ind = findfirst(liq_frac[j,strt:end] .< 0.05) + strt - 1
                        p_gl05[j] = p[j,ind]
                        T_gl05[j] = T[j,ind]
                end
        end

        for j in 1:length(p_gl50)
                strt = find_start(p[j,:])
                if nm.minimum(liq_frac[j,strt:end]) > 0.5
                        p_gl50[j] = NaN
                        T_gl50[j] = NaN
                else
                        ind = findfirst(liq_frac[j,strt:end] .< 0.5) + strt - 1
                        p_gl50[j] = p[j,ind]
                        T_gl50[j] = T[j,ind]
                end
        end

        for j in 1:length(p_gl99)
                stop = find_stop(p[j,:])
                if (qc_end[j] != 0) | (nm.maximum(qc[j,:]) .== 0)
                        p_gl99[j] = NaN
                        T_gl99[j] = NaN
                        qfr99[j] = NaN
                        qi99[j] = NaN
                        qni99[j] = NaN
                        ri99[j] = NaN
                        qs99[j] = NaN
                        qns99[j] = NaN
                        rs99[j] = NaN
                        qg99[j] = NaN
                        Ncr99[j] = NaN
                        qv99[j] = NaN

                elseif qc_end[j] == 0
                        ind = stop + 2 - findfirst(reverse(qc[j,1:stop]) .> 0)
                        p_gl99[j] = p[j,ind]
                        T_gl99[j] = T[j,ind]
                        qfr99[j] = qfr[j,ind]
                        qi99[j] = qi[j,ind]
                        qni99[j] = qni[j,ind]
                        ri99[j] = r_all[j,ind]
                        qs99[j] = qs[j,ind]
                        qns99[j] = qns[j,ind]
                        rs99[j] = qs99[j] / qns99[j]
                        qg99[j] = qg[j,ind]
                        Ncr99[j] = qni[j,ind] * r_all[j,ind]
                        qv99[j] = qv[j,ind]
                end
        end

        T_gl05 = T_gl05[(!).(isnan.(T_gl05))]
        p_gl05 = p_gl05[(!).(isnan.(p_gl05))]
        T_gl50 = T_gl50[(!).(isnan.(T_gl50))]
        p_gl50 = p_gl50[(!).(isnan.(p_gl50))]
        T_gl99 = T_gl99[(!).(isnan.(T_gl99))]
        p_gl99 = p_gl99[(!).(isnan.(p_gl99))]

        process_array(T_gl05,tgl05_out,i)
        process_array(p_gl05,pgl05_out,i)
        process_array(T_gl50,tgl50_out,i)
        process_array(p_gl50,pgl50_out,i)
        process_array(T_gl99,tgl99_out,i)
        process_array(p_gl99,pgl99_out,i)

        qfr99 = qfr99[(!).(isnan.(qfr99))]
        qi99 = qi99[(!).(isnan.(qi99))]
        qni99 = qni99[(!).(isnan.(qni99))]
        ri99 = ri99[(!).(isnan.(ri99))]
        qs99 = qs99[(!).(isnan.(qs99))]
        qns99 = qns99[(!).(isnan.(qns99))]
        rs99 = rs99[(!).(isnan.(rs99))]
        qg99 = qg99[(!).(isnan.(qg99))]
        qv99 = qv99[(!).(isnan.(qv99))]
        Ncr99 = Ncr99[(!).(isnan.(Ncr99))]
        tau99 = tau_rel(p_gl99,T_gl99,qni99,qi99,dp["CAP"][i])
        RHi99 = RH_i.(T_gl99,qv99,p_gl99 .* 100)

        process_array(qfr99,qfr99_out,i)
        process_array(qi99,qi99_out,i)
        process_array(qni99,qni99_out,i)
        process_array(ri99,ri99_out,i)
        process_array(qs99,qs99_out,i)
        process_array(qns99,qns99_out,i)
        process_array(rs99,rs99_out,i)
        process_array(qg99,qg99_out,i)
        process_array(Ncr99,Ncr99_out,i)
        process_array(qv99,qv99_out,i)
        process_array(tau99,tau99_out,i)
        process_array(RHi99,RHi99_out,i)
end

println(" ")
println("done")
close(ds_out)
