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

function vectorreplace!(vec::Vector,index::UInt,newvalue)
    vec[index] = newvalue
    return vec
end

#To reduce allocations, you want to in place calc randvec and then take argmin
function exponential_sample!(rng::Random.AbstractRNG,weights::Vector{Float64},randvec::Vector{Float64})
    Random.rand!(rng,randvec)
    randvec .= -1 .* log.(randvec) ./ weights
    argmin(randvec)
end

function random_fill!(rng::Random.AbstractRNG,vec::Vector{Float64})
    for i in 1:length(a)
        vec[i] = rand(rng)
    end
    return nothing
end

#Direct Gillespie Simulation
#Outcomes: array of size kxn (k reactions for system size n) with jumps in kth reaction
function direct_gillespie(initialstate::Vector{Int64},params::Vector{Float64},numsteps::Int64,maxsteps::Int64,
                      outcomes::Array{Int64},calcstor_constructor,rng::Random.AbstractRNG)
    statevector = initialstate

    numreactions = size(outcomes)
    reactioncounters = fill(numsteps,numreactions)
    stoppingvec = reactioncounters .!= 0
    steps = 0
    rates = Vector{Float64}(undef,numreactions)
    zerorates = Vector{Bool}(undef,numreactions)

    calcstor = calcstor_constructor()
    while (steps <= maxsteps-1) & (any(stoppingvec))
        tt_protein_fb_rates!(statevector,params,rates)
        zerorates .= rates .== 0
        if all(zerorates)
            return "All rates trivially zero at state $(statevector) after $(steps) steps"
        end
        totalrate = sum(rates)
        reaction = CDFsample(rng,rates,numreactions,totalrate)
        timestep = log(1/rand(rng))/totalrate
        delta = view(outcomes,reaction,:)
        statevector .+= delta
        updatestorage!(calcstor,statevector,rates,timestep)
        steps += 1
        reactioncounters[reaction] -= 1
        if reactioncounters[reaction] <= 0
            reactioncounters[reaction] = 0
            stoppingvec[reaction] = false
        end
    end
    return calcstor
end