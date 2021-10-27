#!/home/david/julia-1.6.1/bin/julia
#import GZip
if !@isdefined(CommandLine_loaded)
  include("CommandLine.jl")
end
# if !@isdefined(sort_loaded)
#   include("sort.jl")
# end
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
  process::Int64
  time::Float64
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
    "file_dir" => "",
    "out_file"=>"hawkes_test_data.txt",
    "total_time" => 10,
    "lambda_0" => 1, # base process rate
    "sigma_0" => 2, # initial child process rate (see below)
    "child_half_life" => .1, #half-life of child processes 
  )
  cl = get_vals(defaults,cmd_line) # replace defaults with command line values if they are specified
  println("parameters: $defaults")
  seed = defaults["seed"]
  out_file = defaults["out_file"]
  total_time = defaults["total_time"]
  lambda_0 = defaults["lambda_0"]
  sigma_0 = defaults["sigma_0"]
  child_half_life = defaults["child_half_life"]

  data = HawkesPoint[]
  processes = Process[]
  next_process::Int64 = 0 #indexes the processes array
  next_index::Int64 = 0 #indexes the data array
  rng = MersenneTwister(seed)
  t = 0
  i = 0
  
  parent_process = Process(0,-1,-1,0) # dummy parent data for base process
  data_point = HawkesPoint(0,0)
  rho = 0 #fake values for base process
  sigma = lambda_0
  println("dummy parent for base process: $parent_process")
  
  while next_index <= length(data)
    if(next_index > 0) # skip the startup for the base process
      data_point = data[next_index]
      parent_process = processes[data_point.process]
      rho = log(2)/child_half_life # decay rate of child lambda
      sigma = sigma_0
    end
     # OK, we have a new parent. Set up the child process (this will be the base on the first time thru
    child_process = Process(length(processes)+1,parent_process.process_no,parent_process.generation+1,data_point.time)
    push!(processes,child_process) # save the child process
    t0 = child_process.start_time
    t = t0
    while true # generate the next point
      lambda = sigma*exp(-rho*(t-t0))/(2^(child_process.generation))
      y = rand(rng)
      t -= log(1-y)/lambda;
      if t > total_time;break;end
      push!(data,HawkesPoint(child_process.process_no,t))
    end # while true
    next_index += 1
  end # while(next_index <= length(data)
  heapsort(data)
  #  println(out_file,"data:\n")
  open(out_file, "w") do io
    for p in data
      write(io, "$(p.process)  $(p.time)\n")
    end
  end
  

  println("processes:\n")
  for p in processes
    println("$(p.process_no): parent: $(p.parent) generation: $(p.generation)  start_time: $(p.start_time)")
    end
end #main


# execution begins here
if occursin("simulate.jl",PROGRAM_FILE)
  #was simulate.jl called by a command?
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
 
