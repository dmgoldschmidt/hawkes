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
  mark::Vector{String} # the webip
end
function Base.println(p::HawkesPoint)
  println("mark: $(p.mark) time: $(p.time)")
end

function Base.:<(x::HawkesPoint,y::HawkesPoint)
  return x.time < y.time
end

mutable struct TimeSeries
  enip::String
  active::Bool
  start_time::Float64
  data::Vector{HawkesPoint}
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "in_file" => "fts_time_series.jld2",
    "enip_no" => 0,
    "verbose" => false,
  )
  cl = get_vals(defaults,cmd_line) # update defaults with command line values if they are specified
  #  println("parameters: $defaults")
  for (key,val) in defaults
    println("$key: $val")
  end
  # update defaults (if they appeared on the command line)
  in_file = defaults["in_file"]
  enip_no = defaults["enip_no"]
  verbose = defaults["verbose"]
  
  fts_time_series = Dict{String,TimeSeries}()
  @load in_file fts_time_series #load the data dictionary (keyed on enip)
  rnd = my_round(3)
  n = 0
  for enip in keys(fts_time_series)
    n += 1
    if enip_no >= 1 && n != enip_no; continue; end
    data = fts_time_series[enip].data
    println("$(enip): $(length(data)) events")
    if verbose
      for event in data
        println("$(event.time): $(event.mark)")
      end
    end # if verbose
  end #for enip in keys
end #main
  

# execution begins here

if occursin("dump_fts_time_series.jl",PROGRAM_FILE)
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

