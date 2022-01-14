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
  N::Int64
  lambda::Float64
  rho::Float64
  sigma::Float64
  omega::Array{Float64}
end

function Base.println(p::Parameters)
  rnd = my_round(3)
  println("Parameters:\nlambda: $(rnd(p.lambda)) rho: $(rnd(p.rho)) sigma: $(rnd(p.sigma)) \nomega: $(pretty_print(p.omega))")
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

function dsigma(p::Parameters,t::Vector{Float64})
  S0 = 0.0
  dS0 = 0.0
  for j in 1:p.N-1
    t_Nj = t[p.N] - t[j]
    e_Nj = exp(-p.rho*t_Nj)
    S0 += 1-e_Nj
    dS0 += t_Nj*e_Nj;
  end
  return ((p.N - p.lambda) - p.sigma*dS0)/S0
end

function sigma(p::Parameters,t::Vector{Float64})
  sum = 0.0
  for j in 1:p.N-1
    sum += 1-exp(-p.rho*(t[p.N]-t[j]))
  end
  p.sigma = (p.N - p.lambda)*p.rho/sum
end

function Omegas(omega1::Matrix{Float64},p::Parameters,t::Vector{Float64}) # input old p.omega and t, output omega1, new p.omega
  (nrows,ncols) = size(omega1)
  score = 0.0
  omega1 .= 0.0
  omega1[1,1] = 1.0
  for i in 2:nrows
    for j in 1:i
      if j == 1 # base process
        t_ij = t[i]
        khat = p.lambda*t[i]
      else
        t_ij = t[i] - t[j-1]
        khat = p.sigma/p.rho*(1-exp(-p.rho*t_ij))
      end
      omega1[i,j] = p.omega[j]*exp(khat*log(khat) - khat - log(t_ij) - loggamma(khat))
    end
    rsum = sum(omega1[i,:])
    omega1[i,:] ./= rsum
    score += log(rsum)
  end
  sums = sum(omega1,dims=1) # get column sums
  for j in 1:ncols; p.omega[j] = sums[j]/nrows; end
#  println("column sums: $(p.omega)")
  return score
end

function dq_comp(p::Parameters, omega1::Matrix{Float64}, t::Vector{Float64})
  sigma(p,t)      #compute sigma(rho)
  ds = dsigma(p,t)# and sigma'(rho)
  q = dq = 0
  for i = 1:p.N
    for j = 2:i-1
      t_ij = t[i] - t[j]
      e_ij = exp(-p.rho*t_ij)
      khat = p.sigma/p.rho*(1-e_ij)
      dk = (p.rho*ds - p.sigma)/(p.rho*p.rho)*(1-e_ij) + p.sigma/p.rho*t_ij*e_ij
      q += omega1[i,j]*log(khat)/2
      dq += omega1[i,j]/(2*khat)*dk
    end
  end
  return (dq,q)
end

function update_params(p::Parameters, omega1::Matrix{Float64}, t::Vector{Float64}, rho_max::Int64 = 100, eps::Float64 = 1.0e-5)
  p.lambda = p.N*p.omega[1] #  = \sum_{i=1}^N omega1[i,1]
  (dq,q) = dq_comp(p,omega1,t)
  rho1 = p.rho
  while p.rho <= rho_max && dq > 0
    rho1 = p.rho # old value
    p.rho += 10 # new value
    (dq,q) = dq_comp(p,omega1,t)
  end
  if dq < 0 # either dq went from + to - or dq was < 0 initially. If we hit rho_max with dq >= 0, we exit
    rho = p.rho # dq(rho) < 0 and dq(rho1) > 0 
    while abs(rho - rho1) > eps 
      p.rho = (rho+rho1)/2
      (dq,q) = dq_comp(p,omega1,t) # compute dq at the mid-point
      if dq == 0; break; end
      if dq < 0
        rho = p.rho
      else
        rho1 = p.rho
      end
    end
  end
#  println("rho: $(p.rho) sigma: $(p.sigma) q: $q dq: $(dq)")
#  exit(0)
end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "seed" => 12345,
    "in_file" => "hawkes_test_data.txt",
    "ndata" => 10,
    "out_file"=>"",
    "rho_0" => 1, 
    "sigma_0" => 2, # initial child process rate
    "lambda_0" => .1,
    "t_0" => 0.0,
    "max_iters" => 5,
  )
  cl = get_vals(defaults,cmd_line) # update defaults with command line values if they are specified
#  println("parameters: $defaults")
  for (key,val) in defaults
    println("$key: $val")
  end
  # update defaults (if they appeared on the command line)
  seed = defaults["seed"]
  in_file = defaults["in_file"]
  ndata = defaults["ndata"]
  out_file = defaults["out_file"]
  rho_0 = defaults["rho_0"]
  sigma_0 = defaults["sigma_0"]
  lambda_0 = defaults["lambda_0"]
  t_0 = defaults["t_0"]
  max_iters = defaults["max_iters"]

  #Now read the data and initialize the parameters
  stream = tryopen(in_file) # this is from util.jl
  all_lines = readlines(stream)
  if(ndata == 0)
    ndata = length(all_lines)
  end
  nmarks = 0
  data = HawkesPoint[]
  children = Dict{String,Vector{Int64}}()
  t = fill(0.0,ndata)
  for i in 1:ndata
    field = map(string,split(all_lines[i])) # split the ith line into strings
    push!(data, HawkesPoint(field[1],myparse(Float64,field[2])))
    if !haskey(children,field[1])
      children[field[1]] = []
    end
    push!(children[field[1]],i) #the child process is associated with the mark of the parent
#    println(children)
    println(data[i])
  end
  tot_time = data[ndata].time - t_0
  for i in 1:ndata
    t[i] = (data[i].time - t_0)/tot_time # normalize arrival times to [0,1]
    if t[i] < 0 || (i > 1 && t[i] <= t[i-1])
      println(stderr, "negative arrival time or time(s) out of sequence at t[$i].  Bailing out.")
      exit(1)
    end
  end
  rnd = my_round(3)
  println("normalized arrival times:\n$(map(rnd,t))")
  omega = [Float64(ndata+1-j) for j in 1:ndata]
  omega ./= sum(omega)
  params = Parameters(ndata,lambda_0,rho_0,sigma_0,omega)
  for (key,val) in children
    println("mark $(key): children: $(transpose(val))")
  end
  println("omega: $(map(rnd,omega))")

  omega1 = Matrix{Float64}(undef,ndata,ndata)
  omega = Vector{Float64}(undef,ndata)
  
  #OK, begin EM iteration
  for niters in 1:max_iters
    score = Omegas(omega1,params,t) #= compute posterior probability matrix omega1
                                   and state probability vector params.omega =# 
    rnd = my_round(3)
#    println(pretty_print(map(rnd,omega1)))
    println("log likelihood at iteration $niters: $score")

    if niters < max_iters # re-estimate params
      println("Begin iteration $niters")
      params.lambda = params.omega[1]
      update_params(params,omega1,t) # update lambda, rho, and sigma(rho)
      println("updated parameters: lambda: $(rnd(params.lambda)) rho: $(rnd(params.rho)) sigma: $(rnd(params.sigma))")
#      rnd = my_round(5)
      println("omega: $(map(rnd,params.omega))")
    end
  end
end

#   niters = 0
#   while niters <= max_iters
#     k_hat = fill(0.0,(ndata,ndata)) # posterior
#     k_hat_0 = fill(0.0,(ndata,ndata)) # prior
#     omega1[1,1] = k_hat[1,1] = k_hat_0[1,1] = 1.0
#     log_likelihood = 0.
#     for i in 2:ndata
#       t_i = data[i].time
#       row_sum = 0
#       for j in 1:i-1
#         sigma = p0.decay_params[data[j].mark][1]
#         rho = p0.decay_params[data[j].mark][2]
#         if j == 1
#           k_hat_0[i,1] = p0.lambda*t_i
#           t_ij = t_i
#         else
#           t_ij = (j == 1 ? t_i : t_i - data[j].time)
#           k_hat_0[i,j] = sigma*(1 - exp(-rho*t_ij))/rho
#           if k_hat_0[i,j] <= 0 || isnan(k_hat_0[i,j])
#             println("k_hat_0[$i,$j] = $(k_hat_0[i,j]). rho = $rho, t_ij = $t_ij, sigma = $sigma")
#             exit(0)
#           end
#         end
#         omega1[i,j] = p0.omega[j]*sqrt(k_hat_0[i,j])/t_ij #uses Stirling's approximation to Gamma
#         row_sum += omega1[i,j]
#       end #for j
      
#       if row_sum > 0
#         log_likelihood += log(row_sum) # accumulate the posterior log_liklihood of the data
#       else
#         println("row_sum = $row_sum at i = $i")
#         for j in 1:i-1
#           print("$(omega1[i,j]) ")
#         end
#         print("\n")
#         println("k_hat_0: $(k_hat_0[i,i-1]) t_i:$t_i omega[i-1]: $(p0.omega[i-1])")
#         exit(0)
#       end
#       for j in 1:i-1
#         omega1[i,j] /= row_sum
#         if i > j
#           k_hat[i,j] = omega1[i,j] + k_hat[i-1,j]
#         end
#       end
#       # OK, we have the posterior probability of state j at time i 
#       # and the expected number of arrivals at time i, both w.r.t. the prior parameters 
#     end # for i
#     println("log_likelihood at iteration $niters: $log_likelihood")

#     if niters < max_iters
#       println("begin iteration $niters ")
#       #Now we re-estimate the parameters
#       for j in 1:ndata
#         p0.omega[j] = sum(omega1[j+1:ndata,j])/ndata # posterior state parameters
#       end
#       # print("omega: ")
#       # for j in 1:5
#       #   print("$(p0.omega[j]) ")
#       # end
#       # print("\n")
#       # exit(0)
      
#       #re-estimate base process rate
#       p0.lambda = 0
#       for i in 1:ndata
#         p0.lambda += omega1[i,1]/data[i].time # posterior estimate up to time t_i
#       end
#       p0.lambda /= ndata # average rate

#       #finally, we use non-linear least-squares to re-estimate sigma_m and rho_m for every mark m
#       A = Matrix{Float64}(undef,10000,2)
#       b = Vector{Float64}(undef,10000)
#       for (mark,val) in p0.decay_params
# #        println("\nprocess[$mark] has the following children: $(children[mark])")
#         sigma = val[1]
#         rho = val[2]
#         k = 1 #restart A and b 
#         for j in children[mark] # toss in all the children of this mark (they all use the same sigma and rho)
#           t_j = data[j].time
#           for i in (j+1):ndata # get an approximate linear equation on delta_rho & delta_sigma at every subsequent time
#             t_ij = data[i].time - t_j
#             dk_dsigma = k_hat_0[i,j]/sigma
#             dk_drho = t_ij*(sigma/rho - k_hat_0[i,j]) - k_hat_0[i,j]/rho
#             A[k,1] = dk_dsigma;A[k,2] = dk_drho
#             b[k] = k_hat[i,j] - k_hat_0[i,j] # residual
#             k += 1
#           end
#         end
#         A0 = A[1:k,:] #OK, we got k equations
#         b0 = b[1:k]
#         normal_mat = transpose(A0)*A0
#  #       println("got $k equations. det(normal_mat) = $(det(normal_mat))")
#         delta = inv(normal_mat)*(transpose(A0)*b0)
#  #       println("delta_sigma: $(delta[1]) delta_rho: $(delta[2])")
#         f = .1*sigma/abs(delta[1])
#         if f < 1; delta[1] *= f; end
#         sigma += delta[1]
#         f = .1*rho/abs(delta[2])
#         if f < 1; delta[2] *= f; end
#         rho += delta[2]
# #        println("sigma: $sigma  rho: $rho")
#         ls_err = sqrt(dot(b0,b0)/k)
#         val[1] = sigma
#         val[2] = rho
# #        println("sigma/rho restimation for process $mark had rms error = $ls_err")
#       end #for mark
#     end # if niters < max_iters
#     println("re-estimation ended for iteration $niters")
#     niters += 1
#   end #while niters <= max_iters
#   println(p0)
#   # rnd = my_round(3)
#   # for i in 1:ndata
#   #   print("time $i: ")
#   #   for j in 1:i-1
#   #     print(map(rnd,omega1[i,j])," ")
#   #   end
#   #   print("\n")
#   # end            
# end #main



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
 
