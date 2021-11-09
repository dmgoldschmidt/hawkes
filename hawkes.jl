#!/home/david/julia-1.6.1/bin/julia
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

mutable struct Parameters
  lambda::Float64
  omega::Array{Float64}
  decay_params::Dict{String,Vector{Float64}}
end

function Base.println(p::Parameters)
  println("Parameters:\nlambda: $(p.lambda)\ndecay parameters:")
  for (mark,val) in p.decay_params
    println("mark: $mark sigma: $(val[1]) rho: $(val[2])")
  end
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

# mutable struct Process
#   process_no::Int64
#   parent::Int64
#   generation::Int64
#   start_time::Float64
# end

function main(cmd_line = ARGS)    
  defaults = Dict{String,Any}(
    "seed" => 12345,
    "in_file" => "hawkes_test_data.txt",
    "ndata" => 0,
    "out_file"=>"",
    "rho_0" => 1, 
    "sigma_0" => 2, # initial child process rate
    "lambda_0" => .1,
    "max_iters" => 3,
  )
  cl = get_vals(defaults,cmd_line) # replace defaults with command line values if they are specified
#  println("parameters: $defaults")
  for (key,val) in defaults
    println("$key: $val")
  end
  # reset defaults (if they appeared on the command line)
  seed = defaults["seed"]
  in_file = defaults["in_file"]
  ndata = defaults["ndata"]
  out_file = defaults["out_file"]
  rho_0 = defaults["rho_0"]
  sigma_0 = defaults["sigma_0"]
  lambda_0 = defaults["lambda_0"]
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
  println("$ndata data points:\n") 
  p0 = Parameters(lambda_0,[],Dict())
  for i in 1:ndata
    field = map(string,split(all_lines[i])) # split the ith line into strings
    if !haskey(p0.decay_params,field[1])
      nmarks += 1 # new mark, set nominal decay params
      p0.decay_params[field[1]] = [sigma_0,rho_0]
    end
    push!(data, HawkesPoint(field[1],myparse(Float64,field[2])))
    if !haskey(children,field[1])
      children[field[1]] = []
    end
    push!(children[field[1]],i) #the child process is associated with the mark of the parent
#    println(children)
    println(data[i])
  end
  for (key,val) in children
    println("mark $(key): children: $(transpose(val))")
  end

  omega1 = Matrix{Float64}(undef,ndata,ndata)
  #  k_hat = Matrix{Float64}(undef,ndata,ndata)
  rng = MersenneTwister(seed)
  push!(p0.omega,lambda_0)
  sum1 = lambda_0
  for i in 2:ndata
    push!(p0.omega,rand(rng))
    sum1 += p0.omega[i]
  end
  p0.omega ./= sum1

  #OK, begin EM iteration
  niters = 0
  while niters <= max_iters
    k_hat = fill(0.0,(ndata,ndata)) # posterior
    k_hat_0 = fill(0.0,(ndata,ndata)) # prior
    omega1[1,1] = k_hat[1,1] = k_hat_0[1,1] = 1.0
    log_likelihood = 0.
    for i in 2:ndata
      t_i = data[i].time
      row_sum = 0
      for j in 1:i-1
        sigma = p0.decay_params[data[j].mark][1]
        rho = p0.decay_params[data[j].mark][2]
        if j == 1
          k_hat_0[i,1] = p0.lambda*t_i
          t_ij = t_i
        else
          t_ij = (j == 1 ? t_i : t_i - data[j].time)
          k_hat_0[i,j] = sigma*(1 - exp(-rho*t_ij))/rho
          if k_hat_0[i,j] <= 0 || isnan(k_hat_0[i,j])
            println("k_hat_0[$i,$j] = $(k_hat_0[i,j]). rho = $rho, t_ij = $t_ij, sigma = $sigma")
            exit(0)
          end
        end
        omega1[i,j] = p0.omega[j]*sqrt(k_hat_0[i,j])/t_ij #uses Stirling's approximation to Gamma
        row_sum += omega1[i,j]
      end #for j
      
      if row_sum > 0
        log_likelihood += log(row_sum) # accumulate the posterior log_liklihood of the data
      else
        println("row_sum = $row_sum at i = $i")
        for j in 1:i-1
          print("$(omega1[i,j]) ")
        end
        print("\n")
        println("k_hat_0: $(k_hat_0[i,i-1]) t_i:$t_i omega[i-1]: $(p0.omega[i-1])")
        exit(0)
      end
      for j in 1:i-1
        omega1[i,j] /= row_sum
        if i > j
          k_hat[i,j] = omega1[i,j] + k_hat[i-1,j]
        end
      end
      # OK, we have the posterior probability of state j at time i 
      # and the expected number of arrivals at time i, both w.r.t. the prior parameters 
    end # for i
    println("log_likelihood at iteration $niters: $log_likelihood")

    if niters < max_iters
      println("begin iteration $niters ")
      #Now we re-estimate the parameters
      for j in 1:ndata
        p0.omega[j] = sum(omega1[j+1:ndata,j])/ndata # posterior state parameters
      end
      # print("omega: ")
      # for j in 1:5
      #   print("$(p0.omega[j]) ")
      # end
      # print("\n")
      # exit(0)
      
      #re-estimate base process rate
      p0.lambda = 0
      for i in 1:ndata
        p0.lambda += omega1[i,1]/data[i].time # posterior estimate up to time t_i
      end
      p0.lambda /= ndata # average rate

      #finally, we use non-linear least-squares to re-estimate sigma_m and rho_m for every mark m
      A = Matrix{Float64}(undef,10000,2)
      b = Vector{Float64}(undef,10000)
      for (mark,val) in p0.decay_params
#        println("\nprocess[$mark] has the following children: $(children[mark])")
        sigma = val[1]
        rho = val[2]
        k = 1 #restart A and b 
        for j in children[mark] # toss in all the children of this mark (they all use the same sigma and rho)
          t_j = data[j].time
          for i in (j+1):ndata # get an approximate linear equation on delta_rho & delta_sigma at every subsequent time
            t_ij = data[i].time - t_j
            dk_dsigma = k_hat_0[i,j]/sigma
            dk_drho = t_ij*(sigma/rho - k_hat_0[i,j]) - k_hat_0[i,j]/rho
            A[k,1] = dk_dsigma;A[k,2] = dk_drho
            b[k] = k_hat[i,j] - k_hat_0[i,j] # residual
            k += 1
          end
        end
        A0 = A[1:k,:] #OK, we got k equations
        b0 = b[1:k]
        normal_mat = transpose(A0)*A0
 #       println("got $k equations. det(normal_mat) = $(det(normal_mat))")
        delta = inv(normal_mat)*(transpose(A0)*b0)
 #       println("delta_sigma: $(delta[1]) delta_rho: $(delta[2])")
        f = .1*sigma/abs(delta[1])
        if f < 1; delta[1] *= f; end
        sigma += delta[1]
        f = .1*rho/abs(delta[2])
        if f < 1; delta[2] *= f; end
        rho += delta[2]
#        println("sigma: $sigma  rho: $rho")
        ls_err = sqrt(dot(b0,b0)/k)
        val[1] = sigma
        val[2] = rho
#        println("sigma/rho restimation for process $mark had rms error = $ls_err")
      end #for mark
    end # if niters < max_iters
    println("re-estimation ended for iteration $niters")
    niters += 1
  end #while niters <= max_iters
  println(p0)
  # rnd = my_round(3)
  # for i in 1:ndata
  #   print("time $i: ")
  #   for j in 1:i-1
  #     print(map(rnd,omega1[i,j])," ")
  #   end
  #   print("\n")
  # end            
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
 
