#!/home/david/julia-1.6.1/bin/julia
#import GZip
if !@isdefined(CommandLine_loaded)
  include("CommandLine.jl")
end
if !@isdefined(util_loaded)
  include("util.jl")
end
if !@isdefined(sort_loaded)
  include("sort.jl")
end

using Printf
using DelimitedFiles
using Random

mutable struct Parameters
  lambda::Float64
  omega::Array{Float64}
  decay_params::Dict{String,Vector{Float64}}
  #  rho::Array{Float64}
#  sigma::Array{Float64}
end

mutable struct HawkesPoint
  mark::String # the no. of the process that generated this point
  time::Float64
end
function Base.println(p::HawkesPoint)
  println("mark: $(p.mark) time: $(p.time)")
end

function Base.:<(x::HawkesPoint,y::HawkesPoint)
  return x.time < y.time
end

mutable struct Process
  process_no::Int64
  parent::Int64
  generation::Int64
  start_time::Float64
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "seed" => 12345,
     "in_file" => "hawkes_test_data.txt",
    "ndata" => 0,
    "out_file"=>"",
    "rho_0" => 1, 
    "sigma_0" => 2, # initial child process rate
    "lambda_0" => 1,
  )
  cl = get_vals(defaults,cmd_line) # replace defaults with command line values if they are specified
  println("parameters: $defaults")
  seed = defaults["seed"]
  in_file = defaults["in_file"]
  ndata = defaults["ndata"]
  out_file = defaults["out_file"]
  rho_0 = defaults["rho_0"]
  sigma_0 = defaults["sigma_0"]
  lambda_0 = defaults["lambda_0"]

  #Now read the data and initialize the parameters
  
  stream = tryopen(in_file) # this is from util.jl
  all_lines = readlines(stream)
  if(ndata == 0)
    ndata = length(all_lines)
  end
  nmarks = 0
  data = HawkesPoint[]
  println("$ndata data points:\n") 
  p0 = Parameters(lambda_0,[],Dict())
  for i in 1:ndata
    field = map(string,split(all_lines[i])) # split the ith line
    if !haskey(p0.decay_params,field[1])
      nmarks += 1
      p0.decay_params[field[1]] = [sigma_0,rho_0]
    end
    push!(data, HawkesPoint(field[1],myparse(Float64,field[2])))
    println(data[i])
  end
  println("\n$nmarks decay_params:\n",p0.decay_params)

  omega1 = Matrix{Float64}(undef,ndata,ndata)
#  k_hat = Matrix{Float64}(undef,ndata,ndata)
  rng = MersenneTwister(seed)
  sum = 0
  for i in 1:ndata
    push!(p0.omega,rand(rng))
    sum += p0.omega[i]
  end
  p0.omega ./= sum

  #OK, here we go
  k_hat = fill(0.0,(ndata,ndata))
  for i in 1:ndata
    t_i = data[i].time
    sum = 0
    for j in 1:i-1
      sigma = p0.decay_params[data[j].mark][1]
      rho = p0.decay_params[data[j].mark][2]
      if j == 1
        t_j = 0
        k_hat_0 = lambda*t_i
      else
        t_j = (j == 1 ? 0 : data[j].time)
        k_hat_0 = p0.sigma*(exp(-p0.rho*t_j) - exp(-p0.rho*t_i))/p0.rho
      end
      omega1[i,j] = p0.omega[j]*sqrt(k_hat_0)/(t_i-t_j)
      sum += omega1[i,j]
      if i > j
        k_hat[i,j] = omega1[i,j] + k_hat[i-1,j]
      end
    end
    for j in 1:i-1
      omega1[i,j] /= sum
    end
    # OK, we have the posterior state probabilities w.r.t. the prior parameters (p0)
  end
  
end #main


# execution begins here
if occursin("hawkes.jl",PROGRAM_FILE)
  #was hawkes.jl called by a command?
  println("calling main with ARGS = $ARGS")
  main(ARGS) 
else
  if isinteractive()
    print("enter command line: ")
    cmd = readline()
    main(map(string,split(cmd)))
  end
end 


  #  while true
  #     t -= log(1-rand(rng))/lambda_0;
  #     i += 1
  #     if t > total_time;break;end
  #     push!(data,HawkesPoint(next_process,parent,false,t))
  #     println("$i: $t")
  #   end
  #   next_process += 1
  # while i > 1
  #   i -= 1
  #   println(data[i])
  # end
 
