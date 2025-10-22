using NaNMath; nm=NaNMath

#---------------------------------------------------------------------------
#	"Skript for the computation of thermodynamic variables etc."
#	"Created by Cornelis Schwenk 19.10.2023"
#	"c.schwenk@uni-mainz.de"
#---------------------------------------------------------------------------

function e_sat_i(T)
        #------------------------------
        #"This function calculates the saturation vapor pressure over ICE"
	#"T must be in Kelvin, result is given in Pa"
        #------------------------------
        if (T > 273.15) | (T < 110)
                return NaN
        else
                return exp(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)
	end
end

function e_sat_w(T)
        #------------------------------
        #"This function calculates the saturation vapor pressure over WATER"
	#"T must be in Kelvin, result is given in Pa"
        #------------------------------
        if (T > 332) | (T < 123)
                error("T out of bounds. T in (123,332)K")
        end
        return exp(54.842763 - (6763.22 / T) - 4.210*log(T) + 0.000367*T +
                   tanh(0.0415 * (T - 218.8))*(53.878 - (1331.22 / T) -
                                               9.44523*log(T) + 0.014025*T))
end

function qv_sat_w(T,p)
        #------------------------------
        #"This function calculates the saturation qv(T,p) over WATER"
	#"T must be in Kelvin, p in Pa and result is given in kg/kg"
        #------------------------------
        return e_sat_w(T)*0.622/p
end

function qv_sat_i(T,p)
        #------------------------------
        #"This function calculates the saturation qv(T,p) over ICE"
	#"T must be in Kelvin, p in Pa and result is given in kg/kg"
        #------------------------------
        return e_sat_i(T)*0.622/p
end

function RH_w(T,qv,p)
        #------------------------------
        #"This function calculates RH over WATER"
	#"T must be in Kelvin, p in Pa and qv in kg/kg"
        #"Result is given in %"
	#------------------------------
        e = (qv .* p) ./ 0.622
        return 100*e ./ e_sat_w.(T)
end

function RH_i(T,qv,p)
        #------------------------------
        #"This function calculates RH over ICE"
	#"T must be in Kelvin, p in Pa and qv in kg/kg"
        #"Result is given in %"
	#------------------------------
        e = (qv .* p) ./ 0.622
        return 100*e ./ e_sat_i.(T)
end

function Theta(T,p)
        #------------------------------
        #"This function calculates Potential Temperature"
	#"T must be in Kelvin, p in Pa"
        #"Result is given in Kelvin"
	#------------------------------
	return T .* ((100000 ./ p) .^ 0.286)
end

function Dv(T)
        #------------------------------
        #"This function calculates the self-diffusion coefficient of supercooled water"
	#"T must be in Kelvin"
        #"Result is given in m^2/s"
	#"Taken from https://doi.org/10.1103/PhysRevLett.76.2730"
	#------------------------------
	D_v = 13.93 .* ((T ./ 198.7) .- 1) .^ 2.73
	return (D_v .* 1e-9)
end

#---- Cloud parameters-----
global a_geo_ice = 0.835
global b_geo_ice = 0.390

global R_air = 287.05287 #"Specific gas constant for dry air"
function calc_radius_ice(q,qn)
	#-----------------------------
	#"This function calculates the radius of ice according to ICON"
	#"Look at Seifert and Beheng"
	#"radius is given in meters"
	#"qni in 1/kg and qi in kg/kg"
	#-----------------------------	
        m = q ./ qn
        r = 0.5 .* a_geo_ice .* (m .^ b_geo_ice)
	r[qn .< 100] .= NaN
	return r
end

function calc_radius_ice_single(q,qn)
	#-----------------------------
	#"This function calculates the radius of ice according to ICON"
	#"Look at Seifert and Beheng"
	#"radius is given in meters"
	#"qni in 1/kg and qi in kg/kg"
	#-----------------------------	
        m = q ./ qn
        r = 0.5 .* a_geo_ice .* (m .^ b_geo_ice)
	if qn < 100
            return NaN
        else
            return r
        end
end

function rho_air(p,T)
        #------------------------------
	#"This function calculates the density of air"
        #"p in pascal and T in Kelvin"
	#------------------------------	
	return p ./ (R_air .* T)
end

function tau_rel(p,T,qni,qi,CAP=0.5)
	#------------------------------
	#"This function calculates the relaxation timescale for"
	#"supersaturation over ice for an air parcel"
	#"T in kelvin, p in hPa, qni in 1/kg and qi in kg/kg"
	#------------------------------
        #D_v = 2.5e-5 .* ( T ./ 273.15 ) .^ 1.81
        D_v = 3e-5
        Ni = copy(qni)
	if minimum(Ni) < 100
		Ni[Ni .< 100] .= NaN
	end
	r = calc_radius_ice(qi,qni)
	rho = rho_air(p .* 100,T)
	return 1 ./ (4 .* pi .* Ni .* D_v .* CAP .* r .* rho)
end

function tau_rel_single(p,T,qni,qi,CAP=0.5)
	#------------------------------
	#"This function calculates the relaxation timescale for"
	#"supersaturation over ice for an air parcel"
	#"T in kelvin, p in hPa, qni in 1/kg and qi in kg/kg"
	#------------------------------
        #D_v = 2.5e-5 .* ( T ./ 273.15 ) .^ 1.81
        D_v = 3e-5
	Ni = qni
	if Ni < 100
		return NaN
	end
	r = calc_radius_ice_single(qi,qni)
	rho = rho_air(p .* 100,T)
	return 1 ./ (4 .* pi .* Ni .* D_v .* CAP .* r .* rho)
end
