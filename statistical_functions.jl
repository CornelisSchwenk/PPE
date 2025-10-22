function prcntl(ar,n::Int)
        #------------------------------
        #"This function gives you the n-th percentile"
        #"of a 2-dim array along the first dimension"
        #"This can handle NaNs"
        #------------------------------

        ar = sort(ar,dims=1)
        out = zeros(size(ar)[2])
        for j in 1:size(ar)[2]
                a = ar[:,j][(!isnan).(ar[:,j])]
                if length(a) == 0
                        out[j] = NaN
                else
                        i = Int(round(size(a)[1]*n/100,digits=0))
                        if i == 0
                                out[j] = NaN
                        else
                                out[j] = a[i]
                        end
                end
        end
        return out
end

function bin_two_ars(ar1,ar2,bins,func)
        #"This funciton bins func(ar1) into bins of ar2"
        #"used to calculate median and mean of ar1 as function of ar2"
        #"example: y_mean = t6_binning(ydata,xdata,xbins,mean)"
        v = []
        for t in 1:(length(bins)-1)
                if length(ar1[bins[t] .< ar2 .<= bins[t+1]]) == 0
                        push!(v,NaN)
                else
                        push!(v,func(ar1[bins[t] .< ar2 .<= bins[t+1]]))
                end
        end
        return v
end

function t6_binning(ar,t6,func)
        #"Bin func(values) into tau_600 bins."
        #"Used to calculate median and mean as function of t6."
        #"Example: CR_mean = t6_binning(CR,tau_600,mean)."
        #"Uses bin_two_ars function from statistical_functions.jl"
        bins = 0:48
        return bin_two_ars(ar,t6,bins,func)
end

function t6_binning_prcntl(ar,t6,prc)
	bins=0:48
	v = []
	for t in 1:(length(bins)-1)
		if length(ar[bins[t] .< t6 .<= bins[t+1]]) == 0
			push!(v,NaN)
		else
			push!(v,prcntl_1d(ar[bins[t] .< t6 .<= bins[t+1]],prc))
		end
	end
	return Float64.(v)
end

function prcntl_1d(ar,n::Int)
        #------------------------------
        #"This function gives you the n-th percentile"
        #"of a 1-dim array"
        #"This can handle NaNs"
        #------------------------------
        ar = sort(ar)
        a = ar[(!isnan).(ar)]
        if length(a) == 0
                return NaN
        else
                i = Int(round(length(a)*n/100,digits=0))
                return a[i]
        end
end

function anyar_prcntl(ar1,ar2,bins,prc)
        #------------------------------
        #"This function gives you the n-th percentile"
        #"of a 1-dim array ar1 binned into bins of ar2."
        #"This can handle NaNs"
        #------------------------------
        v = Float64[]
        for i in 1:(length(bins)-1)
                if length(ar1[bins[i] .< ar2 .<= bins[i+1]]) == 0
                        push!(v,NaN)
                else
                        s_ar = sort(ar1[bins[i] .< ar2 .< bins[i+1]])
                        push!(v,prcntl_1d(s_ar,prc))
                end
        end
        return v
end

