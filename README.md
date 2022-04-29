# pulse-profiling

So, there are ~100 FITS files, all are needed for the interpolations. If you need an inside of how are the tables organized (so you would be able to get exact calculation results for chosen parameters), let me know

Python files:

driver_for_pics.py - running code, adapted for 3 additional variables (Te, tbb, tauT), otherwise unchanged

interpolations_grid.jl - the code that does interpolations for a grid of Te, Tbb, tauT. Results go to the specific FITS file. This code is not used anymore

plot_stokes_model.py - plots Figs 3-4 from https://www.aanda.org/articles/aa/abs/2021/02/aa39470-20/aa39470-20.html 

polpulse_for_pics.py - your code with reading through FITS tables and linearly interpolating over them. 
