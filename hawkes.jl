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

# mutable struct Process
#   process_no::Int64
#   parent::Int64
#   generation::Int64
#   start_time::Float64
# end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "seed" => 12345,
     "in_file" => "hawkes_test_data.txt",
    "ndata" => 0,
    "out_file"=>"",
    "rho_0" => 1, 
    "sigma_0" => 2, # initial child process rate
    "lambda_0" => .1,
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
  processes = Dict{String},Vector{Int64}}
  println("$ndata data points:\n") 
  p0 = Parameters(lambda_0,[],Dict())
  for i in 1:ndata
    field = map(string,split(all_lines[i])) # split the ith line
    if !haskey(p0.decay_params,field[1])
      nmarks += 1
      p0.decay_params[field[1]] = [sigma_0,rho_0]
    end
    push!(data, HawkesPoint(field[1],myparse(Float64,field[2])))
    push!(processes[field[1]],i)
    println(data[i])
  end
  println("\n$nmarks decay_params:\n",p0.decay_params)

  omega1 = Matrix{Float64}(undef,ndata,ndata)
#  k_hat = Matrix{Float64}(undef,ndata,ndata)
  rng = MersenneTwister(seed)
  sum = p0.omega[1] = lambda_0
  for i in 2:ndata
    push!(p0.omega,rand(rng))
    sum += p0.omega[i]
  end
  p0.omega ./= sum

  #OK, here we go
  k_hat = fill(0.0,(ndata,ndata)) // posterior
  k_hat_0 = fill(0.0,ndata,ndata)) // prior
  omega1[1,1] = k_hat[1,1] = k_hat_0[1,1] = 1.0
  for i in 2:ndata
    t_i = data[i].time
    sum = 0
    for j in 1:i-1
      sigma = p0.decay_params[data[j].mark][1]
      rho = p0.decay_params[data[j].mark][2]
      if j == 1
        k_hat_0[i,1] = p0.lambda*t_i
      else
        t_ij = (j == 1 ? t_i : t_i - data[j].time)
        k_hat_0[i,j] = sigma*(1 - exp(-rho*t_ij))/rho
      end
      omega1[i,j] = p0.omega[j]*sqrt(k_hat_0[i,j])/t_ij #uses Stirling's approximation to Gamma
      sum += omega1[i,j]
    end
    for j in 1:i-1
      omega1[i,j] /= sum
      if i > j
        k_hat[i,j] = omega1[i,j] + k_hat[i-1,j]
      end
    end
    # OK, we have the posterior probability of state j at time i 
    # and the expected number of arrivals at time i, both w.r.t. the prior parameters (p0)
  end # for i

  #Now we re-estimate the parameters
  for j in 1:ndata
    p0.omega[j] = sum[omega1[:,j])/ndata
  end
  p0.lambda = p0.omega[1]

  #finally, we use non-linear least-squares to re-estimate sigma_m and rho_m
  for mark in p0.decay_params
    sigma = p0.decay_params[mark][1]
    rho = p0.decay_params[mark][2]
    A = [0 0]
    b = [0]
    for j in processes[mark]
      t_j = data[j].time
      for i in j+1:ndata
        t_ij = data[i].time - t_j
        A = vcat([k_hat_0[i,j]/sigma (t*(sigma/rho -k_hat[i,j]) - k_hat_0[i,j]/rho)])
        push!(b,k_hat[ij] - k_hat_0[i,j]) # residual
      end
    end
    delta = inverse(transpose(A)*A)*transpose(A)*b
    sigma += delta[1]
    rho += delta[2]
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
 
