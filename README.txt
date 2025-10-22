
This repository contains the code, data, and instructions to reproduce 
the results from our 2025 publication:  
 	
	"Microphysical parameter choices modulate ice content and 
	relative humidity in the outflow of a warm conveyor belt"

	doi for the pre-print: 10.5194/egusphere-2025-1816

The code (without data) will also be available at https://github.com/CornelisSchwenk/

---------------------------------------
---------------------------------------

Instructions:

	- Unzip the archive
	- Navigate to the "CODE/" directory
	- Make sure you have installed Julia
	- Install required packages by running the command

		julia --project=. install_packages.jl

	- Open the Julia REPL by running the command 

		julia --project=. 

	- Then create the plots from the paper by running these commands 

		>> include("create_plots.jl")
		>> plot_Fig2(true)
		>> plot_Fig3(true)
		>> and so on with Fig4, ... Fig16, FigA1, SI1, SI2, ..., SI7.

	- Fig1 is not included here because you need access to ERA5 data			
---------------------------------------
---------------------------------------

Notes on scripts and other files: 

	- general_plots.jl --> building block functions used often for plotting.

	- plotting_functions.jl --> defines many plotting functions.

	- create_plots.jl --> defines only those plots used in the paper. 
				Uses functions from previous two scripts.

	- do_analysis_....jl --> does additional analysis on trajectory data to create
				additional netCDF files with important filtered data.

	- model.jl --> defines the functions around the RF models used in the Paper.

	- calcs.jl --> defines several functions used for calculations. Defines constants,
			imports necessary modules.  

	- statistical_functions.jl --> defines some functions for statistical measures.

	- thermodynamic_functions.jl --> defines some functions for thermodynamic 
						calculations. 
	
	- use_config.jl --> loads useful packages when working with Julia from the REPL. 

	- cartopy_background.py --> only python script just for creating background
				map with borders for plotting. 

	- slurmfiles/ --> directory for output files when running scripts on cluster 
				using slurm. 

---------------------------------------
---------------------------------------

Notes on DATA:

The DATA/ directory contains not only the trajectory data, but also processed data
in the form of the files:

	- PPE_vars_all.nc
	- PPE_vars_end_hours.nc
	- PPE_vars_accum.nc
	- PPE_vars_accum_hours.nc
	- PPE_vars_t600.nc

These contain multiple variables such as the pressure of glaciation, other summary
statistics and so on which are used in the analysis. These files are created by 
the following scripts:

	- do_analysis_all.jl --> creates PPE_vars_all.nc
	- do_analysis_after_hours.jl --> creates PPE_vars_end_hours.nc
	- do_analysis_accum.jl --> creates PPE_vars_accum.nc
	- do_analysis_accum_hours.jl --> creates PPE_vars_accum_hours.nc
	- do_analysis_t600.jl --> creates PPE_vars_t600.nc

These can also be run using

	julia --project=. do_analysis_all.jl

but might take a while and you need sufficient RAM. You can also use multithreading:

	julia --project=. --threads X do_analysis_all.jl

with X the number of threads. However, you have to uncomment some lines and I can make
no promises that there won't be a bus error! 


The file "params.txt" contains information on the perturbed parameters.

Any additional data used in the scope of the study by Oertel 2025 https://doi.org/10.1002/qj.4986 which is relevant for this study, can be found under the doi: 10.35097/ecgs4f56mp3ymjmt


