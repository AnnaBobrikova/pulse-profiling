using ArgParse
using FITSIO

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--te"
            help = "give T_e"
            arg_type = Float64
            default = 0.08
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
#set values
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
    eval(Meta.parse("$(arg)=$(val)"))
end 

T__e = te/1000.0 

include("./cs_fast.jl")

using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab


t__bb = range(0.00015, step=0.00015, stop=0.003) 
tau__T = range(0.5, step=0.1, stop=2.0) 

IsothermalComptonAtmosphere.init_x(150) 
IsothermalComptonAtmosphere.init_μ(9) 
IsothermalComptonAtmosphere.init_Θ(T__e) # theta=kTe/mec2, as we need to set this one before init_atm
init_atmosphere()

f = FITS("CompSlab_$te.fits", "r+")

for ii in tau__T
    for iii in t__bb
        layers = convert(Int, (ii*50) ÷ 1)
        number = convert(Int, (ii*9 + 5) ÷ 1)
        IsothermalComptonAtmosphere.init_τ(layers, ii) #here come tau(nubmer of optical depth levels) and then tau_t, Thomson optical depth of thermalization
        IsothermalComptonAtmosphere.init_Θ(T__e, iii) # theta=kTe/mec2 and t=kTbb/mec2
        IsothermalComptonAtmosphere.set_ScatterNum(number) #20? orders of scattering icluded, 
        I = compute_slab()
        param = [T__e, iii, ii, layers, number]
        data = Dict("I"=>I[:,9:end,1], "Q"=>I[:,9:end,2], "param"=>param[:])
        write(f,data)
    end
end
