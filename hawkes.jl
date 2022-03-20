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
#using DelimitedFiles
#using Random
#using LinearAlgebra
using SpecialFunctions
using JLD2
using Plots

function pretty_print(mat::Matrix{Float64}, digits = 5)
  (nrows,ncols) = size(mat)
  rnd = my_round(digits)
  s = ""
  for i in 1:nrows
    s = s*"\n"*string(map(rnd,mat[i,:]))
  end
  return s
end

mutable struct Parameters
  nevents::Int64
  nstates::Int64
  lambda::Float64
  rho::Float64
  sigma::Array{Float64}
  omega::Array{Float64}
end

function Base.println(p::Parameters)
  rnd = my_round(3)
  println("Parameters:\nlambda: $(rnd(p.lambda)) rho: $(rnd(p.rho)) sigma: $(rnd(p.sigma)) \nomega: $(pretty_print(p.omega))")
end

mutable struct HawkesPoint
  time::Float64
  mark::Vector{String} # the webips seen at this time
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
  events::Array{HawkesPoint}
end

function Omegas(omega1::Matrix{Float64},p::Parameters,t::Vector{Float64})
  # input old p.omega and t, output omega1, new p.omega
  omega = true
  (nrows,ncols) = size(omega1)
  score = 0.0
  omega1 .= 0.0
  omega1[1,1] = 1.0
  for i in 2:nrows
    for j in 1:min(i,ncols)
      if j == 1 # base process
        t_ij = t[i]
        khat = p.lambda*t[i]
      else
        t_ij = t[i] - t[i-j+1]
        khat = p.sigma[i-j+1]/p.rho*(1-exp(-p.rho*t_ij))
      end
      if khat <= 0 
        println(stderr,"at i = $i, j = $j: sigma = $(p.sigma[i-j+1]), rho = $(p.rho),khat = $khat")
        exit(1)
      end
      omega1[i,j] = p.omega[j]*exp(khat*log(khat) - khat - log(t_ij) - loggamma(khat))
      if isnan(omega1[i,j])
        println(stderr,"NaN at i = $i, j= $j, khat = $khat, p.lambda = $(p.lambda),  t[i] = $(t[i])")
        exit(1)
      end
#      println("omega1[$i,$j]: ",omega1[i,j])
    end
    rsum = sum(omega1[i,:])
    omega1[i,:] ./= rsum
    score += log(rsum)
    if rsum <= 0
      println("score == NaN at ($i,$j), khat = $khat")
      exit(1)
    end
  end
  if omega
    sums = sum(omega1,dims=1) # get column sums
    for j in 1:ncols; p.omega[j] = sums[j]/nrows; end
  end
#  println("column sums: $(p.omega)")
  return score
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "in_file" => "fts_time_series.jld2",
    "out_file" => "", # default is no output
    "plot_enip" => -1, # -1:  no plots, 0: plot all series.  n>0:  plot the n^th series read
    "nstates" => 10,
    "nenips" => 0, # 0 gets all time series in in_file. n>0 gets first n.
    "half_life" => .01, # child process intensity will decay to 2^{-10}*initial value after 10% of the interval  
    "sigma_0" => 2, # initial child process rate
    "lambda_0" => 1.0, # base process generates all events on average
    "t_0" => 0.0,
    "max_iters" => 10,
    "eps" => 1.0e-3,
  )
  cl = get_vals(defaults,cmd_line) # update defaults with command line values if they are specified
  #  println("parameters: $defaults")
  verbose = false
  if "v" in cl.option
    verbose = true
  end
  for (key,val) in defaults
    println("$key: $val")
  end
  # update defaults (if they appeared on the command line)
  in_file = defaults["in_file"] # where to read the .jld2 file
  out_file = defaults["out_file"]  # where to plot sigma vs time"
  plot_enip = defaults["plot_enip"] # which series to plot 
  nenips = defaults["nenips"] # how many series to process
  nstates = defaults["nstates"] # how many states to use
  half_life = defaults["half_life"] # decay constant for child processes
  sigma_0 = defaults["sigma_0"]
  lambda_0 = defaults["lambda_0"]
  t_0 = defaults["t_0"]
  max_iters = defaults["max_iters"]
  eps = defaults["eps"]
  if length(out_file) > 0
    out_stream = tryopen(out_file,"w")
  end
  
  #Now read the data and initialize the parameters
  fts_time_series = Dict{String, TimeSeries}()
  @load in_file fts_time_series
  nenips = nseries = length(keys(fts_time_series))
  avg_delta = 0.0

  for enip in keys(fts_time_series)
    if nenips == 0; break;end
    nenips -= 1
    events = fts_time_series[enip].events
    t_0 = fts_time_series[enip].start_time
    nevents = length(events)
    println("processing $enip with $nevents events")
    if t_0 == events[1].time # reduce t_0 by one tick
      t_0 -= (events[nevents].time - events[1].time)/nevents
    end
    t = fill(0.0,nevents)
    tot_time = events[nevents].time - t_0
    
    avg_delta_t = tot_time/nevents
    for i in 1:nevents
      t[i] = (events[i].time - t_0)/tot_time # normalize arrival times to [0,1]
      if t[i] < 0 || (i > 1 && t[i] < t[i-1])
        println(stderr, "negative arrival time or time(s) out of sequence at t[$i].  Bailing out.")
        exit(1)
      end
    end
    rnd = my_round(5)
    #    println("normalized arrival times:\n$(map(rnd,t))")
    lambda = lambda_0/avg_delta_t  # base process generates lambda_0*nevents per unit time 
    omega = fill(1.0/nstates,nstates)
#    rho = Vector{Float64}(undef,nevents)
    sigma = fill(sigma_0,nevents-1) 
    rho = log(2)/half_life #-log(.5)/avg_delta_t
    if verbose; println("rho: $rho");end
    params = Parameters(nevents,nstates,lambda,rho,sigma,omega)
    # println("omega: $(map(rnd,omega))")
    omega1 = Matrix{Float64}(undef,nevents,nstates)
    
    last_score = 0.0
    for niters in 1:max_iters # begin EM iteration ***********************
      score = Omegas(omega1,params,t) #= compute posterior probability matrix omega1
      and state probability vector params.omega =# 
      rnd = my_round(3)
      #    println(pretty_print(map(rnd,omega1)))
      if verbose; println("log likelihood at iteration $niters: $score");end
      if abs((score - last_score)/score) < eps
        println("relative score = $(abs(score-last_score)/score) < $eps.  Exiting")
        break
      end
      last_score = score
      if verbose; println("omega: $(map(rnd,params.omega))"); end
      
      if niters < max_iters # re-estimate params
        if verbose; println("Begin iteration $niters");end
        # now recompute sigmas for all child processes
        avg_sigma = 0.0
        for i in 1:nevents-1
          khat = 0.0
          m = min(nstates-1,nevents-i)
          for j in 1:m  # compute expected total no. of arrivals for child process i 
            khat += omega1[i+j, j+1]
            if omega1[i+j, j+1] <= 0
              println(stderr,"at i = $i, j = $j: omega1[$(i+j),$(j+1)] = $(omega1[i+j,j+1])")
              exit(1)
            end #if
          end #for j in 1:m
          params.sigma[i] = khat*rho #/(1-exp(-rho*(t[i+m]-t[i])))
          if khat > nstates; println(stderr,"khat = $khat at event $i");end
          avg_sigma += params.sigma[i]
        end #for i in 1:nevents-1
        avg_sigma /= nevents-1
      end # if niters < max_iters
    end # for niters < max_iters (end EM iterations)
    sum_sig = sum(params.sigma)
    len_sig = length(params.sigma)
    if verbose; println("sum sigma = $sum_sig, len sigma = $len_sig, mean sigma/rho = $(sum_sig/(len_sig*params.rho))");end
    if length(out_file) > 0
      println(out_stream, "\n*** series $(nseries-nenips):  $nevents events for enip $enip with trigger  $(fts_time_series[enip].trigger) at $(fts_time_series[enip].start_time)")
    end
    x_vals = Vector{Float64}(undef,nevents-1)
    y_vals = Vector{Float64}(undef,nevents-1)
    for i in 1:nevents-1
#     x_vals = fill(0.0,nevents-1)
#     y_vals = fill(0.0,nevents-1)
      max_p = 0.0
      state = 0
      for j in 1:nstates
        if omega1[i,j] > max_p; max_p = omega1[i,j]; state = j; end
      end
      x_vals[i] = events[i].time
      y_vals[i] = params.sigma[i]
      if length(out_file) > 0
        @printf(out_stream,"%.3f %.3f  %d (%.3f) %s ",x_vals[i], y_vals[i],state-1,max_p, events[i].mark)
        println(out_stream,rnd.(omega1[i,:]))
      end
    end #for i in 1:nevents-1
    
    if plot_enip == 0 || plot_enip == nseries - nenips
      opt::String = ""
      plot(x_vals,y_vals,show = true)
      print(stderr,"enter Q to quit, filename to save (.pdf will be appended), or return to continue:  ")
      opt = readline()
      if length(opt) > 0
        if opt[1] == 'Q'; exit(0);end
        if !occursin(".",opt); opt = opt*".pdf";end
        savefig(opt)
        println(stderr, "plot was saved in ",opt)
      end #if length(opt) > 0
    end #if plot_enip
  end #for enip in keys(fts_time_series)
  if out_stream != ""; close(out_stream); end
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

         # if i == 21
          #   println(stderr,"line 297: khat = $khat, sigma = $(params.sigma[i]), rho = $(rho[i]), t[$(i+m)] = $(t[i+m]), t[$i] = $(t[i])")
          # end
 

# function rho_dump(p::Parameters, omega1::Matrix{Float64}, t::Vector{Float64},
#                   rho_min::Float64, rho_max::Float64, delta::Float64)
#   rho_save = p.rho
#   sigma_save = p.sigma
#   p.rho = rho_min
#   while(p.rho <= rho_max)
#     dq,q = dq_drho(p,omega1,t)
#     println("rho: $(p.rho) q: $q dq:$dq")
#     p.rho += delta
#   end
#   p.rho = rho_save
#   p.sigma = sigma_save
# end

# function dsigma(p::Parameters,t::Vector{Float64})
#   S0 = 0.0
#   dS0 = 0.0
#   for j in 1:min(p.nevents,p.nstates)
#     t_j = t[p.nevents] - t[p.nevents-j+1]
#     e_j = exp(-p.rho*t_j)
#     S0 += 1-e_j
#     dS0 += t_j*e_j;
#   end
#   return ((p.nevents - p.lambda) - p.sigma*dS0)/S0
# end

# function sigma(p::Parameters,t::Vector{Float64})
#   sum = 0.0
#   for j in 1:min(p.nevents,p.nstates)
#     t_j = t[p.nevents] - t[p.nevents-j+1]
#     sum += 1-exp(-p.rho*t_j)
#   end
#   p.sigma = (p.nevents - p.lambda)*p.rho/sum
# end

# function dq_drho(p::Parameters, omega1::Matrix{Float64}, t::Vector{Float64})
#   # compute dq/drho and (with little extra effort) q(rho)
#   sigma(p,t)       # compute sigma(rho)
#   ds = dsigma(p,t) # and sigma'(rho)
#   q = dq = 0.0
    
#   for i = 2:p.nevents
#     for j = 2:min(i,p.nstates)
#       t_ij = t[i] - t[i-j+1]
#       e_ij = exp(-p.rho*t_ij)
#       khat = p.sigma/p.rho*(1-e_ij)
#       dk = (p.rho*ds - p.sigma)/(p.rho*p.rho)*(1-e_ij) + p.sigma/p.rho*t_ij*e_ij
#       q += omega1[i,j]*(khat*log(khat) - khat - loggamma(khat))
#       dq += omega1[i,j]*(log(khat) - digamma(khat))*dk
#  #     println("dq_drho($i,$j): omega1: $(omega1[i,j]) t_ij: $t_ij e_ij: $e_ij khat: $khat dk: $dk dq: $dq")
#     end
#   end
#   return (-dq,-q)
# end

# function update_params(p::Parameters, omega1::Matrix{Float64}, t::Vector{Float64}, rho_min::Float64 = .001,
#                        rho_max::Float64 = 100.0, eps::Float64 = 1.0e-5)
#   p.lambda = .1*p.nevents #p.nevents*p.omega[1] #  this is just \sum_{i=1}^nevents omega1[i,1]
#   rho0 = p.rho = rho_min
#   delta = 1.0
#   (dq0,q0) = dq_drho(p,omega1,t) # get initial values
#   if(dq0 < 0)
#     rho0 = p.rho = rho_max
#     delta = -1.0
#   end
#   dq = dq0

#   while rho_min <= p.rho <= rho_max && dq*dq0 > 0 # dq has not changed sign yet
#     rho0 = p.rho
#     dq0 = dq
#     (dq,q) = dq_drho(p,omega1,t) # this also recomputes sigma(rho)
#     println("$(p.rho)  $q $dq")
#     p.rho += delta # new value
#   end
#   if(p.rho == 0);  p.rho = rho_min; end
#   println("rho scan exit: dq($(p.rho)) = $dq")
# #  exit(0)

#   if dq*dq0 < 0 # dq changed sign between rho0 and p.rho = rho0 + delta. Find q_max by binary search
#     rho1 = p.rho 
#     niters = 0
#     while niters < 10 && (abs(rho0 - rho1) > eps || abs(dq) > eps) 
#       p.rho = (rho0+rho1)/2
#       println("p.rho: $(p.rho) rho0: $(rho0) rho1: $(rho1)")
#       (dq,q) = dq_drho(p,omega1,t) # compute dq (and sigma) at the mid-point
#       if dq == 0; break; end
#       if dq*dq0 > 0
#         rho1 = p.rho # dq(p.rho) and dq(rho1) have the same sign, so move rho1 to the midpoint
#       else
#         rho0 = p.rho # dq(p.rho) and dq(rho0) have the same sign, so move rho0 to the midpoint
#       end
#       niters += 1
#     end
#     return true
#   else # if dq never changed sign, then rho \approx rho_min (if dq < 0) or rho_max (if dq > 0)
#     println("No sign change found.  dq = $dq, rho = $(p.rho)")
#   end
# #  println("rho: $(p.rho) sigma: $(p.sigma) q: $q dq: $(dq)")
# #  exit(0)
# end

