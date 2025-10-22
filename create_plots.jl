include("plotting_functions.jl")

function plot_Fig2(save=false)
    p_out = plot_geographical(false)
    if save
        savefig(p_out,savepath*"Fig2.png")
    else
        return p_out
    end
end

function plot_Fig3(save=false)
    p_out = plot_histograms(false)
    if save
        savefig(p_out,savepath*"Fig3.png")
    else
        return p_out
    end
end

function plot_Fig4(save=false)
    msk = true
    p_out = plot_histograms_CAP(msk,false)
    if save
        savefig(p_out,savepath*"Fig4.png")
    else
        return p_out
    end
end

function plot_Fig5(save=false)
    msk = true
    p_out = plot_histograms_INP(msk,false)
    if save
        savefig(p_out,savepath*"Fig5.png")
    else
        return p_out
    end
end

function plot_Fig6(save=false)
    p_out = plot_CCN_figure(false)
    if save
        savefig(p_out,savepath*"Fig6.png")
    else
        return p_out
    end
end

function plot_Fig7(save=false)
    p_out = plot_RHi_2d_hists(false)
    if save
        savefig(p_out,savepath*"Fig7.png")
    else
        return p_out
    end
end

function plot_Fig8(save=false)
    p_out = plot_CAP_INP_CCN_scatters_RHi_tau(false)
    if save
        savefig(p_out,savepath*"Fig8.png")
    else
        return p_out
    end
end

function plot_Fig9(save=false)
    p_out = plot_RHi_tau_CAPINP(true,false)
    if save
        savefig(p_out,savepath*"Fig9.png")
    else
        return p_out
    end
end

function plot_Fig10(save=false)
    p_out = plot_tau_RHi_hist_all(false)
    if save
        savefig(p_out,savepath*"Fig10.png")
    else
        return p_out
    end
end

function plot_Fig11(save=false)
    p_out = plot_glac_pres_temp_hist2d(false)
    if save
        savefig(p_out,savepath*"Fig11.png")
    else
        return p_out
    end
end

function plot_Fig12(save=false)
    p_out = plot_only_histograms_INP_CAP(dsa["t600"][:] .< 10,false)
    if save
        savefig(p_out,savepath*"Fig12.png")
    else
        return p_out
    end
end

function plot_Fig13(save=false)
        dsa = dsd["accum"]
        msk = dsa["t600"][:] .< 10
        p_out = plot_RHi_tau_CAPINP(msk,false)
        plot!(p_out,legendfontsize=10,labelfontsize=13,tickfontsize=10)
        if save
                savefig(p_out,savepath*"Fig13.png")
        else
                return p_out
        end
end

function plot_Fig14(save=false)
    p_out = plot_SAT_effects(false)
    if save
        savefig(p_out,savepath*"Fig14.png")
    else
        return p_out
    end
end

function plot_Fig15(save=false)
    p_out = plot_histograms_CCN_conv()
    if save
        savefig(p_out,savepath*"Fig15.png")
    else
        return p_out
    end
end

function plot_Fig16(save=false)
    p_out = plot_histograms_5h(false)
    if save
        savefig(p_out,savepath*"Fig16.png")
    else
        return p_out
    end
end

function plot_FigA1(save=false)
    p_out = plot_qv_T_correlations(false)
    if save
        savefig(p_out,savepath*"FigA1.png")
    else
        return p_out
    end
end

function plot_Supplement_1(save=false)
    p_out = plot_SI1()
    if save
        savefig(p_out,savepath*"SI1.png")
    else
        return p_out
    end
end

function plot_Supplement_2(save=false)
    p_out = plot_SI2()
    if save
        savefig(p_out,savepath*"SI2.png")
    else
        return p_out
    end
end


function plot_Supplement_3(save=false)
    p_out = plot_SI3()
    if save
        savefig(p_out,savepath*"SI3.png")
    else
        return p_out
    end
end


function plot_Supplement_4(save=false)
    p_out = plot_SI4()
    if save
        savefig(p_out,savepath*"SI4.png")
    else
        return p_out
    end
end


function plot_Supplement_5(save=false)
    p_out = plot_SI5()
    if save
        savefig(p_out,savepath*"SI5.png")
    else
        return p_out
    end
end


function plot_Supplement_6(save=false)
    p_out = plot_SI6()
    if save
        savefig(p_out,savepath*"SI6.png")
    else
        return p_out
    end
end

function plot_Supplement_7(save=false)
    p_out = plot_SI7()
    if save
        savefig(p_out,savepath*"SI7.png")
    else
        return p_out
    end
end
