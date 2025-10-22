include("calcs.jl")
using Base.Threads

pth = "../DATA/"
ds_out_name = "PPE_vars_accum_hours.nc"
ds_out = Dataset(pth*ds_out_name,"c")

num_trajs = 2799307

defDim(ds_out,"trajs",num_trajs)

ds_out.attrib["title"] = "Accumulated variables 5 hours after the ascent for all simulations"

ds = dsd["normal"]

lon_out = defVar(ds_out,"lon",Float64,("trajs",))
lon_out[:] = zeros(num_trajs)

lat_out = defVar(ds_out,"lat",Float64,("trajs",))
lat_out[:] = zeros(num_trajs)

p_out = defVar(ds_out,"p",Float64,("trajs",))
p_out[:] = zeros(num_trajs)

t6_out = defVar(ds_out,"t600",Float64,("trajs",))
t6_out[:] = zeros(num_trajs)

t_out = defVar(ds_out,"t",Float64,("trajs",))
t_out[:] = zeros(num_trajs)

qv_out = defVar(ds_out,"qv",Float64,("trajs",))
qv_out[:] = zeros(num_trajs)

qc_out = defVar(ds_out,"qc",Float64,("trajs",))
qc_out[:] = zeros(num_trajs)

qi_out = defVar(ds_out,"qi",Float64,("trajs",))
qi_out[:] = zeros(num_trajs)

qni_out = defVar(ds_out,"qni",Float64,("trajs",))
qni_out[:] = zeros(num_trajs)

tau_out = defVar(ds_out,"tau",Float64,("trajs",))
tau_out[:] = zeros(num_trajs)

RHi_out = defVar(ds_out,"RHi",Float64,("trajs",))
RHi_out[:] = zeros(num_trajs)

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

#You CAN multithread this code, but I make no promises. If you want to,
#then uncomment below and comment further below

#@threads for i in 1:70
for i in 1:70
    m = i + 19
    print("$(m), ")
    
    ix1 = fromto[i,1]
    ix2 = fromto[i,2]
    ts = 10

    lon = load_vr_ts_withnan("lon",m,ts)
    lat = load_vr_ts_withnan("lat",m,ts)

    t6 = calc_t600(m)
    t = load_vr_ts_withnan("t",m,ts)
    p = load_vr_ts_withnan("p",m,ts)
    qv = load_vr_ts_withnan("qv",m,ts)
    qc = load_vr_ts_withnan("qc",m,ts)
    qni = load_vr_ts_withnan("qni",m,ts)
    qi = load_vr_ts_withnan("qi",m,ts)
    RHi = RH_i.(t,qv,p .* 100)
    tau = tau_rel(p,t,qni,qi,CAPS[i])

    #----lock----
    #uncomment this when multithreading. However, if there is an error,
    #the code will get stuck here. 
    
    #lock(lock_obj)

    lon_out[ix1:ix2] = lon
    lat_out[ix1:ix2] = lat
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
