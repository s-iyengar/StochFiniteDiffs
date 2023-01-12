import Random

#Utility functions for Gillespie
function CDFsample(rng::Random.AbstractRNG, weights::Vector{Float64},n::Int, weightsum::Float64)
    t = rand(rng) * weightsum
    i = 1
    cuw = weights[1] #current cumulative sum
    while cuw < t && i < n #check if t in current part of weights, and you're in bounds
        i += 1
        @inbounds cuw += weights[i]
    end
    return i
end

function relative_error(x,y)
    abs(x-y)/y
end

#Direct Gillespie Simulation
#statevector: n-length vector of unsigned integers. Input is initial condition
#params: m-length list of numbers used in rate calculations
#Outcomes: array of size kxn (k reactions for system size n) with jumps in kth reaction (ex. [1 0;-1 0;0 1;0 -1 ])
#rate_calc: takes the state, params, and a k-d vector and fills vector with the rates of the system
#calcstor: any object you want to use to calculate things: ex. a struct with a vector to calculate the means
#updatestorage!: a function that has a method for your calcstor, calculating based on state, rates, params, and timestep
#rng: any generator, best is Random.MersenneTwister()
#numsteps: minimum number of times each reaction should occur before simulation stops
#maxsteps: maximum number of steps before stopping simulation
function direct_gillespie!(statevector::Vector{UInt},params::Vector{Float64},outcomes::Array{Int},rate_calc!::Function,
                            calcstor,updatestorage!::Function,rng::Random.AbstractRNG,numsteps::UInt,maxsteps::UInt)

    numreactions = size(outcomes)
    reactioncounters = fill(numsteps,numreactions)
    stoppingvec = reactioncounters .!= 0
    steps = 0
    rates = Vector{Float64}(undef,numreactions)
    zerorates = Vector{Bool}(undef,numreactions)

    while (steps <= maxsteps) & (any(stoppingvec))
        rate_calc!(statevector,params,rates)
        zerorates .= rates .== 0
        if all(zerorates)
            println("All rates trivially zero at state $(statevector) after $(steps) steps")
            return statevector,steps
        end
        totalrate = sum(rates)
        reaction = CDFsample(rng,rates,numreactions,totalrate)
        timestep = log(1/rand(rng))/totalrate
        delta = view(outcomes,reaction,:)
        statevector .+= delta
        updatestorage!(calcstor,statevector,rates,params,timestep)
        steps += 1
        reactioncounters[reaction] -= 1
        if reactioncounters[reaction] <= 0
            reactioncounters[reaction] = 0
            stoppingvec[reaction] = false
        end
    end
    return statevector,steps
end