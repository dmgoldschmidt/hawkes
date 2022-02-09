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

mutable struct Parameters
  ndata::Int64
  nstates::Int64
  lambda::Float64
  rho::Array{Float64}
  sigma::Array{Float64}
  omega::Array{Float64}
end

function Base.println(p::Parameters)
  rnd = my_round(3)
  println("Parameters:\nlambda: $(rnd(p.lambda)) rho: $(rnd(p.rho)) sigma: $(rnd(p.sigma)) \nomega: $(pretty_print(p.omega))")
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

function rho_dump(p::Parameters, omega1::Matrix{Float64}, t::Vector{Float64},
                  rho_min::Float64, rho_max::Float64, delta::Float64)
  rho_save = p.rho
  sigma_save = p.sigma
  p.rho = rho_min
  while(p.rho <= rho_max)
    dq,q = dq_drho(p,omega1,t)
    println("rho: $(p.rho) q: $q dq:$dq")
    p.rho += delta
  end
  p.rho = rho_save
  p.sigma = sigma_save
end

function dsigma(p::Parameters,t::Vector{Float64})
  S0 = 0.0
  dS0 = 0.0
  for j in 1:min(p.ndata,p.nstates)
    t_j = t[p.ndata] - t[p.ndata-j+1]
    e_j = exp(-p.rho*t_j)
    S0 += 1-e_j
    dS0 += t_j*e_j;
  end
  return ((p.ndata - p.lambda) - p.sigma*dS0)/S0
end

function sigma(p::Parameters,t::Vector{Float64})
  sum = 0.0
  for j in 1:min(p.ndata,p.nstates)
    t_j = t[p.ndata] - t[p.ndata-j+1]
    sum += 1-exp(-p.rho*t_j)
  end
  p.sigma = (p.ndata - p.lambda)*p.rho/sum
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
        khat = p.sigma[i-j+1]/p.rho[i-j+1]*(1-exp(-p.rho[i-j+1]*t_ij))
      end
      omega1[i,j] = p.omega[j]*exp(khat*log(khat) - khat - log(t_ij) - loggamma(khat))
#      println("omega1[$i,$j]: ",omega1[i,j])
    end
    rsum = sum(omega1[i,:])
    omega1[i,:] ./= rsum
    score += log(rsum)
    if score == NaN
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

function dq_drho(p::Parameters, omega1::Matrix{Float64}, t::Vector{Float64})
  # compute dq/drho and (with little extra effort) q(rho)
  sigma(p,t)       # compute sigma(rho)
  ds = dsigma(p,t) # and sigma'(rho)
  q = dq = 0.0
    
  for i = 2:p.ndata
    for j = 2:min(i,p.nstates)
      t_ij = t[i] - t[i-j+1]
      e_ij = exp(-p.rho*t_ij)
      khat = p.sigma/p.rho*(1-e_ij)
      dk = (p.rho*ds - p.sigma)/(p.rho*p.rho)*(1-e_ij) + p.sigma/p.rho*t_ij*e_ij
      q += omega1[i,j]*(khat*log(khat) - khat - loggamma(khat))
      dq += omega1[i,j]*(log(khat) - digamma(khat))*dk
 #     println("dq_drho($i,$j): omega1: $(omega1[i,j]) t_ij: $t_ij e_ij: $e_ij khat: $khat dk: $dk dq: $dq")
    end
  end
  return (-dq,-q)
end

function update_params(p::Parameters, omega1::Matrix{Float64}, t::Vector{Float64}, rho_min::Float64 = .001,
                       rho_max::Float64 = 100.0, eps::Float64 = 1.0e-5)
  p.lambda = .1*p.ndata #p.ndata*p.omega[1] #  this is just \sum_{i=1}^ndata omega1[i,1]
  rho0 = p.rho = rho_min
  delta = 1.0
  (dq0,q0) = dq_drho(p,omega1,t) # get initial values
  if(dq0 < 0)
    rho0 = p.rho = rho_max
    delta = -1.0
  end
  dq = dq0

  while rho_min <= p.rho <= rho_max && dq*dq0 > 0 # dq has not changed sign yet
    rho0 = p.rho
    dq0 = dq
    (dq,q) = dq_drho(p,omega1,t) # this also recomputes sigma(rho)
    println("$(p.rho)  $q $dq")
    p.rho += delta # new value
  end
  if(p.rho == 0);  p.rho = rho_min; end
  println("rho scan exit: dq($(p.rho)) = $dq")
#  exit(0)

  if dq*dq0 < 0 # dq changed sign between rho0 and p.rho = rho0 + delta. Find q_max by binary search
    rho1 = p.rho 
    niters = 0
    while niters < 10 && (abs(rho0 - rho1) > eps || abs(dq) > eps) 
      p.rho = (rho0+rho1)/2
      println("p.rho: $(p.rho) rho0: $(rho0) rho1: $(rho1)")
      (dq,q) = dq_drho(p,omega1,t) # compute dq (and sigma) at the mid-point
      if dq == 0; break; end
      if dq*dq0 > 0
        rho1 = p.rho # dq(p.rho) and dq(rho1) have the same sign, so move rho1 to the midpoint
      else
        rho0 = p.rho # dq(p.rho) and dq(rho0) have the same sign, so move rho0 to the midpoint
      end
      niters += 1
    end
    return true
  else # if dq never changed sign, then rho \approx rho_min (if dq < 0) or rho_max (if dq > 0)
    println("No sign change found.  dq = $dq, rho = $(p.rho)")
  end
#  println("rho: $(p.rho) sigma: $(p.sigma) q: $q dq: $(dq)")
#  exit(0)
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "seed" => 12345,
    "in_file" => "fts_flowsets.jld2",
    "max_data" => 50,
    "nstates" => 20,
    "nenips" => 0, # this gets all flowsets in in_file
    "out_file"=>"",
    "rho_0" => 1, 
    "sigma_0" => 2, # initial child process rate
    "lambda_0" => 10,
    "t_0" => 0.0,
    "max_iters" => 10,
    "eps" => 1.0e-3
  )
  cl = get_vals(defaults,cmd_line) # update defaults with command line values if they are specified
  #  println("parameters: $defaults")
  for (key,val) in defaults
    println("$key: $val")
  end
  # update defaults (if they appeared on the command line)
  seed = defaults["seed"]
  in_file = defaults["in_file"]
  nenips = defaults["nenips"]
  nstates = defaults["nstates"]
  out_file = defaults["out_file"]
  rho_0 = defaults["rho_0"]
  sigma_0 = defaults["sigma_0"]
  lambda_0 = defaults["lambda_0"]
  t_0 = defaults["t_0"]
  max_iters = defaults["max_iters"]
  eps = defaults["eps"]

  #Now read the data and initialize the parameters
  fts_flowsets = Dict{String,Flowset}()
  @load in_file fts_flowsets

  for enip in keys(fts_flowsets)
    data = fts_flowsets[enip].data
    t_0 = fts_flowsets[enip].start_time   
    ndata = length(data)

    t = fill(0.0,ndata)
    tot_time = data[ndata].time - t_0
    for i in 1:ndata
      t[i] = (data[i].time - t_0)/tot_time # normalize arrival times to [0,1]
      if t[i] < 0 || (i > 1 && t[i] <= t[i-1])
        println(stderr, "negative arrival time or time(s) out of sequence at t[$i].  Bailing out.")
        exit(1)
      end
    end
    
    # rnd = my_round(3)
    # println("normalized arrival times:\n$(map(rnd,t))")
    omega = fill(1.0/nstates,nstates)
    rho = Vector{Float64}(undef,ndata)
    sigma = fill(sigma_0,ndata)
    for i in 1:ndata-1
      j = min(i+nstates,ndata)
      rho[i] = 10*log(2)/(t[j]-t[i])        # sigma[i]*e^{-rho[i]*(t[j]-t[i])} = 2^{-10}sigma[i]
    end
    params = Parameters(ndata,nstates,lambda_0,rho,sigma,omega)
    println("omega: $(map(rnd,omega))")
    omega1 = Matrix{Float64}(undef,ndata,nstates)

    last_score = 0.0
    for niters in 1:max_iters # begin EM iteration ***********************
      score = Omegas(omega1,params,t) #= compute posterior probability matrix omega1
      and state probability vector params.omega =# 
      rnd = my_round(3)
      #    println(pretty_print(map(rnd,omega1)))
      println("log likelihood at iteration $niters: $score")
      if abs((score - last_score)/score) < eps
        println("relative score = $(abs(score-last_score)/score) < $eps.  Exiting")
        break
      end
      last_score = score
      
      if niters < max_iters # re-estimate params
        println("Begin iteration $niters")
        # now recompute sigmas for all complete child processes
        for i in 1:ndata-nstates-1
          khat = 0.0
          for j in 1:min(nstates-1,ndata - i) # compute expected total no. of arrivals for process i 
            khat += omega1[i+j, j+1]
          end
          params.sigma[i] = khat*rho[i]/(1-2^(-10))
        end
        params.lambda = params.omega[1]
        println("omega: $(map(rnd,params.omega))")
      end #re-estimation
    end # EM iteration
    println("updated parameters for $enip: lambda: $(rnd(params.lambda)) \nrho: $(rnd.(params.rho)) \nsigma: $(rnd.(params.sigma))")
    nenips -= 1
    if nenips == 0;break;end
  end #for enip in keys(fts_flowsets)
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
