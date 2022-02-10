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
  mark::Array{String} # the webip
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
  data::HawkesPoint
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "in_file" => "wsa.raw.1M.txt",
    "out_file"=> "fts_time_series.jld2",
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
  fts_time_series = Dict{String,TimeSeries}()
  wsa_stream = tryopen(in_file)
  readline(wsa_stream);readline(wsa_stream) #skip two header lines
  nflowsets = 0
  nstarts = 0
  rnd = my_round(5)
  while max_flowsets > 0 ? nflowsets < max_flowsets : true
    if !eof(wsa_stream)
      raw_line = readline(wsa_stream)
    else
      println(stderr,"EOF on $in_file after $nflowsets flowsets read")
      break
    end
    fields = map(String,split(raw_line,"|"))
    time = [myparse(Float64,fields[1])]
    wbip = fields[6]
    enip = fields[4]
#    println("read $enip at $time, wbip = $wbip")
    if  !(enip in keys(fts_time_series))
      if wbip in keys(rare_webips) && (max_flowsets > 0 ? nstarts < max_flowsets : true)
        # start a new TimeSeries
        data = HawkesPoint[]
        fts_time_series[enip] = TimeSeries(enip,true,time,data)
        println("found trigger $enip at $time")
        nstarts += 1
      else
        continue # we're not tracking this enip and the webip is common
      end
      # we started a new TimeSeries for this enip
    end 
    #OK we have a TimeSeries for this enip
    time_series = fts_time_series[enip]
    if time < time_series.start_time + duration # and it hasn't expired
      n = length(time_series.data)
      if n > 1 && time_series.data[n].time == time_series.data[n-1].time
        # don't make a new event with the same time.  Just add the wbip to the existing mark
         push!(time_series.data[n].mark,wbip)
      else # make a new event
        push!(fts_time_series[enip].data, HawkesPoint(time,[wbip]))
      end 
    elseif fts_time_series[enip].active == true
      nflowsets += 1
      fts_time_series[enip].active = false
      delta_t = rnd(time - fts_time_series[enip].start_time)
      println("flowset $enip completed after $delta_t seconds with $(length(fts_time_series[enip].data)) Hawkes points")
    end #if
  end #read loop
  
  @save out_file fts_time_series
end #main

# execution begins here

if occursin("get_fts_time_series.jl",PROGRAM_FILE)
  #was get_fts_time_series.jl called by a command?
  println("calling main with ARGS = $ARGS")
  main(ARGS) 
else
  if isinteractive()
    print("enter command line: ")
    cmd = readline()
    main(map(string,split(cmd)))
  end
end 
