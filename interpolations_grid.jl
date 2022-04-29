using FITSIO
using Interpolations

#introducing basic ranges and variables
t__e = range(40.0, step=4.0, stop=200.0)
t__bb = range(0.00015, step=0.00015, stop=0.003) 
tau__T = range(0.5, step=0.1, stop=3.5) 
NEnergy = 150
NZenith = 9

#this arrays contain the values of the parameters for which we will do the interpolations
Task_Te = [42.0, 122.0, 190.0]
Task_Tbb = [0.00025, 0.00158, 0.0028]
Task_tau = [0.55, 1.95, 3.25]

#we're collecting all the values from all the FITS files into these huge 5D tensors:
I_mighty = zeros(length(t__e), length(t__bb), length(tau__T), NEnergy, NZenith)
Q_mighty = zeros(length(t__e), length(t__bb), length(tau__T), NEnergy, NZenith)

#three loops in which we collect the data from FITS files, nothing interesting here. there are 3 loops, because a calculation of each of these FITS files take a week, and that's the maximum time my tasks can have on titan
p=1
for i in t__e
	f = FITS("CompSlab_$i.fits")
    println("with Compslab_$i.fits still works")
    for ii in 1:length(t__bb)
        for iii in 1:16
            data = read(f[(iii-1)*length(t__bb)+ii+1], "I")
            data2 = read(f[(iii-1)*length(t__bb)+ii+1], "Q")
            for kk in 1:NEnergy
                for kkk in 1:NZenith
                    I_mighty[p, ii, iii, kk, kkk] = data[kk, kkk+1]
                    Q_mighty[p, ii, iii, kk, kkk] = data2[kk, kkk+1]
                end
            end
        end
    end
    global p +=1
end

p=1
for i in t__e
	f = FITS("CompSlab_bigTe_$i.fits")
    println("with Compslab_bigTe_$i.fits still works")
    for ii in 1:length(t__bb)
        for iii in 1:10
            data = read(f[(iii-1)*length(t__bb)+ii+1], "I")
            data2 = read(f[(iii-1)*length(t__bb)+ii+1], "Q")
            for kk in 1:NEnergy
                for kkk in 1:NZenith
                    I_mighty[p, ii, iii+16, kk, kkk] = data[kk, kkk+1]
                    Q_mighty[p, ii, iii+16, kk, kkk] = data2[kk, kkk+1]
                end
            end
        end
    end
    global p +=1
end 

p=1
for i in t__e
	f = FITS("CompSlab_hugeT_$i.fits")
    println("with Compslab_hugeT_$i.fits still works")
    for ii in 1:length(t__bb)
        for iii in 1:5
            data = read(f[(iii-1)*length(t__bb)+ii+1], "I")
            data2 = read(f[(iii-1)*length(t__bb)+ii+1], "Q")
            for kk in 1:NEnergy
                for kkk in 1:NZenith
                    I_mighty[p, ii, iii+26, kk, kkk] = data[kk, kkk+1]
                    Q_mighty[p, ii, iii+26, kk, kkk] = data2[kk, kkk+1]
                end
            end
        end
    end
    global p +=1
end 

#maybe this all can be done with less amount of matrices, but I'm just glad it works this way
I_3D = zeros(length(t__e), length(t__bb), length(tau__T))
P_3D = zeros(length(t__e), length(t__bb), length(tau__T)) #these two will contain values for I and A for all Te, Tbb and TauT, but E and mu are fixed. 
I_2D = zeros(NEnergy, NZenith)
P_2D = zeros(NEnergy, NZenith) #these two will contain the interpolation results of log(I) and PD (Q/I) for specific values of Te, Tbb and TauT (given by the Task, one combination at a time) for all Es and mus.
real_I_2D = zeros(NEnergy, NZenith) 
Q_2D = zeros(NEnergy, NZenith) #with these two we're going back to values of I and Q parameters. 

ff = FITS("CompSlab_Int3D_testgrid.fits", "r+")

for te in Task_Te
    x = 1+ (te - 40.0)*(length(t__e)-1)/(200.0-40.0)
    for tau in Task_tau
        z = 1+ (tau - 0.5)*(length(tau__T)-1)/(3.5-0.5)
        for tb in Task_Tbb               
            y = 1+ (tb - 0.00015)*(length(t__bb)-1)/(0.003-0.00015) #cicles over the grid of Task values, so the interpolation can be done more than for 1 set of values
            for kk in 1:NEnergy
                for kkk in 1:NZenith
                    for i in 1:length(t__e)
                        for ii in 1:length(t__bb)
                            for iii in 1:length(tau__T)
                                I_3D[i, ii, iii]=log(I_mighty[i, ii, iii, kk, kkk])
                                P_3D[i, ii, iii]=Q_mighty[i, ii, iii, kk, kkk]/I_mighty[i, ii, iii, kk, kkk] #transforming the I and Q into logI and PD for better interpolation
                            end
                        end
                    end
                    itp_3D = interpolate(I_3D, BSpline(Linear()))
                    itp_Q_3D = interpolate(P_3D, BSpline(Linear())) #interpolations are done for fixed E and Mu
                    I_2D[kk, kkk] = itp_3D(x, y, z) 
                    P_2D[kk, kkk] = itp_Q_3D(x, y, z) #exact numbers for specific Te, Tbb and Tau are written into 2D martices in the places reserved for this combination of fixed E and Mu
                    real_I_2D[kk,kkk]=exp(I_2D[kk,kkk])
                    Q_2D[kk,kkk] = P_2D[kk,kkk]*real_I_2D[kk,kkk]     #and we're going back to I and Q               
                end
            end
            data = Dict("I"=>real_I_2D[:,:], "Q"=>Q_2D[:,:], "param"=>[te, tb, tau]) 
            write(ff, data)
        end
    end
end

