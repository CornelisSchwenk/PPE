include("calcs.jl")
using Base.Threads

pth = "../DATA/"
ds_out_name = "PPE_vars_accum.nc"
ds_out = Dataset(pth*ds_out_name,"c")

num_trajs = 2799307

defDim(ds_out,"trajs",num_trajs)

ds_out.attrib["title"] = "Accumulated variables at end of ascent for all simulations"

ds = dsd["normal"]

p_out = defVar(ds_out,"p",Float64,("trajs",))
p_out[:] = zeros(num_trajs)

t6_out = defVar(ds_out,"t600",Float64,("trajs",))
t6_out[:] = zeros(num_trajs)

p_gl_out = defVar(ds_out,"p_gl",Float64,("trajs",))
p_gl_out[:] = zeros(num_trajs)

lon_out = defVar(ds_out,"lon",Float64,("trajs",))
lon_out[:] = zeros(num_trajs)

lat_out = defVar(ds_out,"lat",Float64,("trajs",))
lat_out[:] = zeros(num_trajs)

qidep_out = defVar(ds_out,"qidep",Float64,("trajs",))
qidep_out[:] = zeros(num_trajs)

t_out = defVar(ds_out,"t",Float64,("trajs",))
t_out[:] = zeros(num_trajs)

t_gl_out = defVar(ds_out,"t_gl",Float64,("trajs",))
t_gl_out[:] = zeros(num_trajs)

qv_out = defVar(ds_out,"qv",Float64,("trajs",))
qv_out[:] = zeros(num_trajs)

qv_gl_out = defVar(ds_out,"qv_gl",Float64,("trajs",))
qv_gl_out[:] = zeros(num_trajs)

qc_out = defVar(ds_out,"qc",Float64,("trajs",))
qc_out[:] = zeros(num_trajs)

qc_gl_out = defVar(ds_out,"qc_gl",Float64,("trajs",))
qc_gl_out[:] = zeros(num_trajs)

qi_out = defVar(ds_out,"qi",Float64,("trajs",))
qi_out[:] = zeros(num_trajs)

qi_gl_out = defVar(ds_out,"qi_gl",Float64,("trajs",))
qi_gl_out[:] = zeros(num_trajs)

qni_out = defVar(ds_out,"qni",Float64,("trajs",))
qni_out[:] = zeros(num_trajs)

qni_gl_out = defVar(ds_out,"qni_gl",Float64,("trajs",))
qni_gl_out[:] = zeros(num_trajs)

tau_out = defVar(ds_out,"tau",Float64,("trajs",))
tau_out[:] = zeros(num_trajs)

tau_gl_out = defVar(ds_out,"tau_gl",Float64,("trajs",))
tau_gl_out[:] = zeros(num_trajs)

RHi_out = defVar(ds_out,"RHi",Float64,("trajs",))
RHi_out[:] = zeros(num_trajs)

RHi_gl_out = defVar(ds_out,"RHi_gl",Float64,("trajs",))
RHi_gl_out[:] = zeros(num_trajs)

sim_num = defVar(ds_out,"sim",Float64,("trajs",))
sim_num[:] = zeros(num_trajs)

CAP_out = defVar(ds_out,"CAP",Float64,("trajs",))
CAP_out[:] = zeros(num_trajs)
INP_out = defVar(ds_out,"INP",Float64,("trajs",))
INP_out[:] = zeros(num_trajs)
SST_out = defVar(ds_out,"SST",Float64,("trajs",))
SST_out[:] = zeros(num_trajs)
CCN_out = defVar(ds_out,"CCN",Float64,("trajs",))
CCN_out[:] = zeros(num_trajs)
SAT_out = defVar(ds_out,"SAT",Float64,("trajs",))
SAT_out[:] = zeros(num_trajs)

CAPS = [dp["CAP"][i] for i in 1:70]
INPS = [dp["INP"][i] for i in 1:70]
SATS = [dp["SATAD"][i] for i in 1:70]
CCNS = [dp["CCN"][i] for i in 1:70]
SSTS = [dp["SST"][i] for i in 1:70]

fromto = Int.(zeros(70,2));
fromto[1,1] = 1
fromto[1,2] = ds["ntrajs"][1]

for i in 2:70
    fromto[i,1] = fromto[i-1,2]+1
    fromto[i,2] = fromto[i,1] -1 + ds["ntrajs"][i]
end

ntrajs = ds["ntrajs"][:]

global lock_obj = ReentrantLock()


#You CAN multithread this, but I make no promises
#@threads for i in 1:70

for i in 1:70
    m = i + 19
    print("$(m), ")
    
    ix1 = fromto[i,1]
    ix2 = fromto[i,2]

    t6 = calc_t600(m)
    lat = load_vr_end("lat",m)
    lon = load_vr_end("lon",m)
    qidep = load_vr_end("qidep",m)
    t = load_vr_end("t",m)
    p = load_vr_end("p",m)
    qv = load_vr_end("qv",m)
    qc = load_vr_end("qc",m)
    qni = load_vr_end("qni",m)
    qi = load_vr_end("qi",m)
    RHi = RH_i.(t,qv,p .* 100)
    tau = tau_rel(p,t,qni,qi,CAPS[i])

    #at glaciation
    qc_all = load_vr("qc",m)
    qv_all = load_vr("qv",m)
    qi_all = load_vr("qi",m)
    qni_all = load_vr("qni",m)
    p_all = load_vr("p",m)
    t_all = load_vr("t",m)
    
    temp_tgl   = zeros(length(t))
    temp_pgl   = zeros(length(t))
    temp_qvgl  = zeros(length(t))
    temp_qigl  = zeros(length(t))
    temp_qnigl = zeros(length(t))
    temp_qcgl  = zeros(length(t))
    temp_RHigl = zeros(length(t))
    temp_taugl = zeros(length(t))
    temp_qidepgl = zeros(length(lon))

    for j in 1:length(t6)
        stop = find_stop(p_all[j,:])
        if (qc[j] != 0) | (nm.maximum(qc_all[j,:]) .== 0) #if the parcel doesn't glaciate or if it never had any qc
            temp_tgl[j]   = NaN
            temp_pgl[j]   = NaN
            temp_qvgl[j]  = NaN
            temp_qigl[j]  = NaN
            temp_qnigl[j] = NaN
            temp_qcgl[j]  = NaN
            temp_RHigl[j] = NaN
            temp_taugl[j] = NaN
        elseif qc[j] == 0 #here I am taking the first index where qc is zero
            ix0 = stop + 2 - findfirst(reverse(qc_all[j,1:stop]) .> 0) 
   
            temp_tgl[j]   = t_all[j,ix0]
            temp_pgl[j]   = p_all[j,ix0]
            temp_qvgl[j]  = qv_all[j,ix0]
            temp_qigl[j]  = qi_all[j,ix0]
            temp_qnigl[j] = qni_all[j,ix0]
            temp_qcgl[j]  = qc_all[j,ix0]
            temp_RHigl[j] = RH_i(t_all[j,ix0],qv_all[j,ix0],p_all[j,ix0]*100)
            temp_taugl[j] = tau_rel_single(p_all[j,ix0],t_all[j,ix0],qni_all[j,ix0],qi_all[j,ix0],CAPS[i])
        end
    end

    #----lock----
    #This is not reccommended but you can use these locks when multithreading
    #lock(lock_obj)

    t_gl_out[ix1:ix2]   = temp_tgl  
    p_gl_out[ix1:ix2]   = temp_pgl  
    qv_gl_out[ix1:ix2]  = temp_qvgl 
    qi_gl_out[ix1:ix2]  = temp_qigl 
    qni_gl_out[ix1:ix2] = temp_qnigl
    qc_gl_out[ix1:ix2]  = temp_qcgl 
    RHi_gl_out[ix1:ix2] = temp_RHigl
    tau_gl_out[ix1:ix2] = temp_taugl
    
    lon_out[ix1:ix2] = lon
    lat_out[ix1:ix2] = lat
    qidep_out[ix1:ix2] = qidep
    t_out[ix1:ix2] = t
    p_out[ix1:ix2] = p
    qv_out[ix1:ix2] = qv
    qc_out[ix1:ix2] = qc
    qni_out[ix1:ix2] = qni
    qi_out[ix1:ix2] = qi
    RHi_out[ix1:ix2] = RHi
    tau_out[ix1:ix2] = tau
    t6_out[ix1:ix2] = t6

    CAP_out[ix1:ix2] = [CAPS[i] for k in 1:ntrajs[i]]
    INP_out[ix1:ix2] = [INPS[i] for k in 1:ntrajs[i]]
    SST_out[ix1:ix2] = [SSTS[i] for k in 1:ntrajs[i]]
    CCN_out[ix1:ix2] = [CCNS[i] for k in 1:ntrajs[i]]
    SAT_out[ix1:ix2] = [SATS[i] for k in 1:ntrajs[i]]
    sim_num[ix1:ix2] = [i for k in 1:ntrajs[i]]
    #unlock(lock_obj)
    #----unlock----
end
println("")
println("done")
close(ds_out)
