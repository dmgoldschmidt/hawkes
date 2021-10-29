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

mutable struct HawkesPoint
  mark::String # the no. of the process that generated this point
  time::Float64
end
function Base.println(p::HawkesPoint)
  println("mark: $(p.process) time: $(p.time)")
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
    "total_time" => 10,
    "lambda_0" => 1, # base process rate
    "sigma_0" => 2, # initial child process rate (see below)
    "child_half_life" => .1, #half-life of child processes 
  )
  cl = get_vals(defaults,cmd_line) # replace defaults with command line values if they are specified
  println("parameters: $defaults")
  seed = defaults["seed"]
  in_file = defaults["in_file"]
  ndata = defaults["ndata"]
  out_file = defaults["out_file"]
  total_time = defaults["total_time"]
  lambda_0 = defaults["lambda_0"]
  sigma_0 = defaults["sigma_0"]
  child_half_life = defaults["child_half_life"]
  
  stream = tryopen(in_file) # this is from util.jl
  all_lines = readlines(stream)
  if(ndata == 0)
    ndata = length(all_lines)
  end
  decay_params = Dict{String,Vector{Float64}()
  nmarks = 0
  data = HawkesPoint[]
  println("$ndata data points:\n") 
  for i in 1:ndata
    field = map(string,split(all_lines[i])) # split the ith line
#    println("field: ",field)
    if !haskey(decay_params,field[1])
      nmarks += 1
      marks[field[1]] = [sigma_0,rho_0]
    end
    push!(data, HawkesPoint(field[1],myparse(Float64,field[2])))
    println(data[i])
  end
  println("\n$nmarks marks:\n",marks)

  omega = Vector{Float64}(undef,ndata)
  sigma = Vector{Float64}(undef,nmarks)
  rho = Vector{Float64}(undef,nmarks)
  omega1 = Matrix{Float64}(undef,ndata,ndata) 
  for i in 1:nmarks
    sigma[i] = sigma_0
    rho[i] = rho_0

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
 
