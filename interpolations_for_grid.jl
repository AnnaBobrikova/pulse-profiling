using FITSIO
using Interpolations

t__e = range(40.0, step=4.0, stop=200.0)
t__bb = range(0.00015, step=0.00015, stop=0.003) 
tau__T = range(0.5, step=0.1, stop=3.5) 
NEnergy = 150
NZenith = 9

#Task_Te = [42.0, 122.0, 190.0]
#Task_Tbb = [0.00025, 0.00158, 0.0028]
#Task_tau = [0.55, 1.95, 3.25]

Task_Te = [97.0]
Task_Tbb = [0.002]
Task_tau = [1.0]

#=x = 1+ (Task_Te - 40.0)*(length(t__e)-1)/(200.0-40.0)
y = 1+ (Task_Tbb - 0.00015)*(length(t__bb)-1)/(0.003-0.00015)
z = 1+ (Task_tau - 0.5)*(length(tau__T)-1)/(3.5-0.5) =#

I_mighty = zeros(length(t__e), length(t__bb), length(tau__T), NEnergy, NZenith)
Q_mighty = zeros(length(t__e), length(t__bb), length(tau__T), NEnergy, NZenith)
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

I_3D = zeros(length(t__e), length(t__bb), length(tau__T))
I_2D = zeros(NEnergy, NZenith)
P_3D = zeros(length(t__e), length(t__bb), length(tau__T))
P_2D = zeros(NEnergy, NZenith)
#I_functions = zeros(NEnergy, NZenith)
real_I_2D = zeros(NEnergy, NZenith)
Q_2D = zeros(NEnergy, NZenith)

#ff = FITS("CompSlab_Int3D_testgrid.fits", "r+")

for te in Task_Te
    x = 1+ (te - 40.0)*(length(t__e)-1)/(200.0-40.0)
    for tau in Task_tau
        z = 1+ (tau - 0.5)*(length(tau__T)-1)/(3.5-0.5)
        for tb in Task_Tbb               
            y = 1+ (tb - 0.00015)*(length(t__bb)-1)/(0.003-0.00015)
            for kk in 1:NEnergy
                for kkk in 1:NZenith
                    for i in 1:length(t__e)
                        for ii in 1:length(t__bb)
                            for iii in 1:length(tau__T)
                                I_3D[i, ii, iii]=log(I_mighty[i, ii, iii, kk, kkk])
                                P_3D[i, ii, iii]=Q_mighty[i, ii, iii, kk, kkk]/I_mighty[i, ii, iii, kk, kkk]
                            end
                        end
                    end
                    itp_3D = interpolate(I_3D, BSpline(Linear()))
                    itp_Q_3D = interpolate(P_3D, BSpline(Linear()))
                    I_2D[kk, kkk] = itp_3D(x, y, z)
                    P_2D[kk, kkk] = itp_Q_3D(x, y, z)
                    real_I_2D[kk,kkk]=exp(I_2D[kk,kkk])
                    Q_2D[kk,kkk] = P_2D[kk,kkk]*real_I_2D[kk,kkk]                    
                end
            end
            #data = Dict("I"=>real_I_2D[:,:], "Q"=>Q_2D[:,:], "param"=>[te, tb, tau])
            #write(ff, data)
            f = open("I_int_jul.txt", "w")
            for i in real_I_2D
                println(f,i)
            end
            close(f)

            f = open("Q_int_jul.txt", "w")
            for i in Q_2D
                println(f,i)
            end
            close(f)
            
	    println("one more grid point")
        end
    end
end
I = zeros(NEnergy, NZenith, 2)
I[:,:,1]=real_I_2D[:,:]
I[:,:,2]=Q_2D[:,:]

using PyCall 
#numpy = pyimport("numpy")
#Intensity = numpy.array(I)
Intensity = PyObject(I)
#muP = PyObject(mu)
#xP = PyObject(x)

outI = open("CompI_jl.bin","w")
#outx = open("Compx_jl.bin","w")
#outm = open("Compm_jl.bin","w")
Intensity.tofile(outI,format="%e")
#xP.tofile(outx,format="%e")
#muP.tofile(outm,format="%e")