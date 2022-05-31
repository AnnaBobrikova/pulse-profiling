# pulse-profiling

So, there are ~100 FITS files, all are needed for the interpolations. 
For each Te, we have 3 FITS tables: 1 for small tau_T (0.5-2.0), 1 for medium tau_T (2.1-3.0), 1 for large tau_T (3.1-3.5). As you can see, the bigger the tau_T, the longer it takes to calculate 1 table. 

Typical file that creates 1 set of tables is test.jl. T_e are coming from the job.sh, which submits the job to the cluster. 

Reading-and-interpolating routines are:
- on julia - interpolations_for_grid.jl
- on python, commented - interpolations_grid.py which is run from driver_for_interpolation.py

Pulse-profiling routine:
polpulse.py run from driver.py - Tuomo&Vlad's code in which I replaced the SimpleThomson calculations with reading-and-interpolating routine from interpolations_grid.py

Plotting:
poldeg_plot.py, beaming_plot.py - unchanged Tuomo's codes
plot_stokes_model.py - barely changed Tuomo's code

Bin-files are not really needed, these were simply created while fixing the bugs in pulse profiling. RN polpulse.py uses Compx_150.bin, I'll fix it soon, surely this vector can be calculated directy within the code. 

Reminder for myself:
- add te,tbb,tauT to calling compf in driver.py
- create a loop over them in driver.py if needed for MC simulations
- replace reading Compx_150.bin in polpulse.py
