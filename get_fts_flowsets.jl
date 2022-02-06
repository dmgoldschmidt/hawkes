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
using FileIO

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
  mark::String # the no. of the process that generated this point
  time::Float64
end
function Base.println(p::HawkesPoint)
  println("mark: $(p.mark) time: $(p.time)")
end

function Base.:<(x::HawkesPoint,y::HawkesPoint)
  return x.time < y.time
end

struct Flowset
  enip::String
  start_time::Float64
  data::Array{HawkesPoint}
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "in_file" => "wsa.raw.1M.txt",
    "out_file"=> "fts_flowsets.jld2",
    "rare_file" => "rare_webips.jld2",
    "max_flowsets" => Int64(50),
    "duration" => Float64(240.0),
  )
  cl = get_vals(defaults,cmd_line) # update defaults with command line values if they are specified
  #  println("parameters: $defaults")
  for (key,val) in defaults
    println("$key: $val")
  end
  # update defaults (if they appeared on the command line)
  seed = defaults["seed"]
  in_file = defaults["in_file"]
  rare_file = defaults["rare_file"]
  max_flowsets = defaults["max_flowsets"]
  duration = defaults["duration"]

  # now get rare_webips 
  rare_webips = Dict{String,Int64}
  if !isfile("rare_webips.jld2") #dictionary does not exist
    rare_stream = tryopen("rare_wbips.txt")
    for line in eachline(rare_stream)
      l = tuple(split(line)...)
      name = l[1]
      level = tryparse(Int64,l[2])
      rare_webips[name] = level
    end
    close(rare_stream)
    @save "rare_webips.jld2" rare_webips
   else #load dictionary from disk
    @load "rare_webips.jld2" rare_webips
  end
  
  #now read raw wsa data, find a rare_webip, and write 4 minutes of data
  fts_flowsets = Dict{String,Flowset}
  wsa_stream = tryopen(infile)
  readline(wsa_stream);readline(wsa_stream) #skip two header lines
  nflowsets = 0
  for nflowsets <= max_flowsets
    if !eof(wsa_stream)
      raw_line = readline(wsa_stream)
    else
      println(std_err,"EOF on $infile after $nflowsets flowsets read")
      exit(1)
    end
    line = map(String,split(raw_line,"|"))
    time = myparse(Float64,line[1])
    wbip = line[6]
    enip = line[4]
    if !haskey(fts_flowsets, enip)
      if wbip in rare_webips # start a new Flowset
        data = HawkesPoint[]
        fts_flowsets[enip] = Flowset(enip,time,data)
        nflowsets += 1
      else
        continue #nothing to do for this record
      end
    end
    #OK, we have a flowset for this enip
    flowset = fts_flowsets[enip]
    if time < flowset.start_time + duration # we're still recording
       push!(flowset.data, HawkesPoint(time,wbip))
    end
  end
  @save outfile fts_flowsets
end #main

# execution begins here

if occursin("get_fts_flowsets.jl",PROGRAM_FILE)
  #was get_fts_flowsets.jl called by a command?
  println("calling main with ARGS = $ARGS")
  main(ARGS) 
else
  if isinteractive()
    print("enter command line: ")
    cmd = readline()
    main(map(string,split(cmd)))
  end
end 
