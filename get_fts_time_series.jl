#!/usr/local/bin/julia
import GZip
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
  trigger::String
  active::Bool
  start_time::Float64
  events::Vector{HawkesPoint}
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "in_file" => "../wsa/wsa.2013-07-01.gz",
    "rare_file_raw" => "../rare_wsa_complete.txt",
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
  rare_file_raw = defaults["rare_file_raw"]
  out_file = defaults["out_file"]   
  rare_file = defaults["rare_file"]
  max_flowsets = defaults["max_flowsets"]
  duration = defaults["duration"]

  # now get rare_webips 
  rare_webips = Dict{String,String}()
  if !isfile("rare_webips.jld2") #dictionary does not exist
    rare_stream = tryopen(rare_file_raw)
    for line in eachline(rare_stream)
      l = map(String,split(line,"|"))
      rare_webips[l[1]] = l[2]
    end
    close(rare_stream)
    @save "rare_webips.jld2" rare_webips
   else #load dictionary from disk
    @load "rare_webips.jld2" rare_webips
  end
  
  #now read raw wsa data, find a rare_webip, and write 4 minutes of data
  fts_time_series = Dict{String,TimeSeries}()
  wsa_stream = occursin(".gz",in_file) ? GZip.open(in_file) : open(in_file)
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
    x = myparse(Float64,fields[1])/1.0e9
    time = round(((((x*1000)%1)*100)%1)*10000,sigdigits = 7)
    wbip = fields[6]
    enip = fields[4]
#    println("read $enip at $time, wbip = $wbip")
    if  !(enip in keys(fts_time_series)) #have we seen this enip before?
      if wbip in keys(rare_webips) && (max_flowsets > 0 ? nstarts < max_flowsets : true)
        # start a new TimeSeries
        events = HawkesPoint[]
        fts_time_series[enip] = TimeSeries(enip,wbip,true,time,events)
        println("found trigger $wbip for  $enip at $time")
        nstarts += 1
      else
        continue # we're not tracking this enip and the webip is common
      end
      # we started a new TimeSeries for this enip (no events yet)
    end 
    #OK we have a TimeSeries for this enip
    time_series = fts_time_series[enip]
    if time < time_series.start_time + duration # and it hasn't expired
      n = length(time_series.events)
      if n > 1 && time == time_series.events[n].time
        # don't make a new event with the same time.  Just add the wbip to the existing mark
         push!(time_series.events[n].mark,wbip)
      else # make a new event and add it to the time series
        push!(fts_time_series[enip].events, HawkesPoint(time,[wbip]))
      end #if n>1 
    elseif fts_time_series[enip].active == true
      #time has expired.  Stop adding events. 
      nflowsets += 1
      fts_time_series[enip].active = false
      delta_t = rnd(time - fts_time_series[enip].start_time)
      println("flowset $enip completed after $delta_t seconds with $(length(fts_time_series[enip].events)) Hawkes points")
    end #if time < start_time
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
