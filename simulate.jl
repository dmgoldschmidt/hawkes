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

using Printf
using DelimitedFiles
using Random

mutable struct HawkesPoint
  process::Int64
  time::Float64
end

mutable struct Process
  parent::Int64
  generation::Int64
  childen::Array{Int64}
  start_time::Float64
end

function has_children(process::Process)::Bool
  return length(process.children) > 0
end


function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "seed" => 12345,
    "file_dir" => "",
    "out_file"=>"hawkes_test_data.txt",
    "total_time" => 10,
    "k_hat_base" => 10, #expected total no. of base process events, so base rate is 10/10 = 1
    "k_hat_child" => 2, # expected total number of child process events
    "child_half_life" => .1, #half-life of child processes 
  )
  cl = get_vals(defaults,cmd_line) # replace defaults with command line values if they are specified
  println("parameters: $defaults")
  seed = defaults["seed"]
  out_file = defaults["out_file"]
  total_time = defaults["total_time"]
  k_hat_base = defaults["k_hat_base"]
  k_hat_child = defaults["k_hat_child"]
  child_half_file = defaults["child_half_life"]

  data = HawkesPoint[]
  processes = Process[]
  next_process::Int64 = 0 #indexes the processes array
  next_parent::Int64 = 1 #indexes the data array
  rng = MersenneTwister(seed)
  t = 0
  i = 0
  lambda_0 = k_hat_base/total_time;
  while next_parent < length(data)
    parent_point = data[next_parent]
    parent_process = processes[parent_point.process]
    if has_children(parent_process)
      next_parent += 1
      continue
    end
    # OK, we have a new parent. Set up the child process
    this_process = Process(next_parent,parent_process.generation+1,[],parent_point.time)
    push!(processes,this_process) # save the new process
    push!(parent_process.children,length(processes)) # add the index of the child process to the parent process  
    
    while true # generate the next point
      lambda = 2^(-this_process.generation)*sigma*exp(-rho*(t - this_process.start_time))
      t0 -= log(1-rand(rng))/lambda;
      if t0 > total_time;break;end
      t = t0
      push!(data,HawkesPoint(this_process,t))
      println("$t")
    end # while true
    println("\ndata for process $this_process")
    for d in data
      println(d)
    end # for
  end #while next_parent < length(data)


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
 
