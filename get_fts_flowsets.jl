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

struct HawkesPoint
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
    "in_file" => "wsa.raw.1M.txt",
    "out_file"=> "fts_flowsets.jld2",
    "rare_file" => "rare_webips.jld2",
    "max_flowsets" => Int64(5),
    "duration" => Float64(240.0),
  )
  cl = get_vals(defaults,cmd_line) # update defaults with command line values if they are specified
  #  println("parameters: $defaults")
  for (key,val) in defaults
    println("$key: $val")
  end
  # update defaults (if they appeared on the command line)
  in_file = defaults["in_file"]
  out_file = defaults["out_file"]
  rare_file = defaults["rare_file"]
  max_flowsets = defaults["max_flowsets"]
  duration = defaults["duration"]

  # now get rare_webips 
  rare_webips = Dict{String,Int64}()
  if !isfile("rare_webips.jld2") #dictionary does not exist
    rare_stream = tryopen("rare_wbips.txt")
    for line in eachline(rare_stream)
      l = tuple(split(line)...)
      name = String(l[1])
      level = tryparse(Int64,l[2])
      rare_webips[name] = level
    end
    close(rare_stream)
    @save "rare_webips.jld2" rare_webips
   else #load dictionary from disk
    @load "rare_webips.jld2" rare_webips
  end
  
  #now read raw wsa data, find a rare_webip, and write 4 minutes of data
  fts_flowsets = Dict{String,Flowset}()
  wsa_stream = tryopen(in_file)
  readline(wsa_stream);readline(wsa_stream) #skip two header lines
  nflowsets = 0
  nstarts = 0
  while nflowsets < max_flowsets
    if !eof(wsa_stream)
      raw_line = readline(wsa_stream)
    else
      println(std_err,"EOF on $infile after $nflowsets flowsets read")
      exit(1)
    end
    fields = map(String,split(raw_line,"|"))
    time = myparse(Float64,fields[1])
    wbip = fields[6]
    enip = fields[4]
#    println("read $enip at $time, wbip = $wbip")
    if  !(enip in keys(fts_flowsets))
      if wbip in keys(rare_webips) && nstarts < max_flowsets
        # start a new Flowset
        data = HawkesPoint[]
        fts_flowsets[enip] = Flowset(enip,true,time,data)
        println("found trigger $enip at $time")
        nstarts += 1
      else
        continue # we're not tracking this enip and the webip is common
      end
      # we started a new flowset for this enip
    end 
    #OK we have a flowset for this enip
    if time < fts_flowsets[enip].start_time + duration # and it hasn't expired
      push!(fts_flowsets[enip].data, HawkesPoint(time,wbip))
    elseif fts_flowsets[enip].active == true
      nflowsets += 1
      fts_flowsets[enip].active = false
      println("flowset $enip completed at time $time with $(length(fts_flowsets[enip].data)) Hawkes points")
    end #if
  end #read loop
  @save out_file fts_flowsets
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
