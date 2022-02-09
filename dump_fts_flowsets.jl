#!/home/david/julia-1.6.5/bin/julia
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
using LinearAlgebra
using SpecialFunctions
using JLD2

function pretty_print(mat::Matrix{Float64}, digits = 5)
  (nrows,ncols) = size(mat)
  rnd = my_round(digits)
  s = ""
  for i in 1:nrows
    s = s*"\n"*string(map(rnd,mat[i,:]))
  end
  return s
end

mutable struct HawkesPoint
  time::Float64
  mark::String # the webip
end
function Base.println(p::HawkesPoint)
  println("mark: $(p.mark) time: $(p.time)")
end

function Base.:<(x::HawkesPoint,y::HawkesPoint)
  return x.time < y.time
end

mutable struct Flowset
  enip::String
  active::Bool
  start_time::Float64
  data::Array{HawkesPoint}
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "in_file" => "fts_flowsets.jld2",
  )
  cl = get_vals(defaults,cmd_line) # update defaults with command line values if they are specified
  #  println("parameters: $defaults")
  for (key,val) in defaults
    println("$key: $val")
  end
  # update defaults (if they appeared on the command line)
  in_file = defaults["in_file"]

  fts_flowsets = Dict{String,Flowset}()
  @load in_file fts_flowsets
  rnd = my_round(3) 
  for enip in keys(fts_flowsets)
    println("$(enip): $(length(fts_flowsets[enip].data)) events")
  end
  end #main

# execution begins here

if occursin("dump_fts_flowsets.jl",PROGRAM_FILE)
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
 
  # stream = tryopen(in_file) # this is from util.jl
  # all_lines = readlines(stream)
  # if(ndata == 0)
  #   ndata = length(all_lines)
  # end
  # nmarks = 0
  # data = HawkesPoint[]
  # children = Dict{String,Vector{Int64}}()
  # t = fill(0.0,ndata)
  # for i in 1:ndata
  #   field = map(string,split(all_lines[i])) # split the ith line into strings
  #   push!(data, HawkesPoint(field[1],myparse(Float64,field[2])))
  #   if !haskey(children,field[1])
  #     children[field[1]] = []
  #   end
  #   push!(children[field[1]],i) #the child process is associated with the mark of the parent
  #   #    println(children)
  #   # println(data[i])
  # end
