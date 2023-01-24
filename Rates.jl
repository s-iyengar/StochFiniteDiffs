function constantrate(param::Float64)
    return param
end

function linearrate(param::Float64,moleculelevel::Int64)
    return param*moleculelevel
end

function hillfunction(lambda::Float64,n::Float64,K::Float64,x::Int64)
    if n >= 0
        return lambda*(x^n)/(K^n + x^n)
    else
        N = abs(n)
        return lambda*(K^N)/(K^N+x^N)
    end
end

function transcription_translation_rates!(state::Vector{Int64},params::Vector{Float64},rates::Vector{Float64})
    rates[1] = constantrate(params[1])
    rates[2] = linearrate(params[2],state[1])
    rates[3] = linearrate(params[3],state[1])
    rates[4] = linearrate(params[4],state[2])
    return nothing
end

function cov_balance_comp(CiRmj,CjRmi,Dij,CiRpj,CjRpi)
    LHS = CiRmj+CjRmi
    RHS = CiRpj+CjRpi+Dij
    return relative_error(RHS,LHS)
end

function tt_protein_fb_rates!(state::Vector,params::Vector,rates::Vector)
    rates[1] = hillfunction(params[1],params[5],params[6],state[2])
    rates[2] = linearrate(params[2],state[1])
    rates[3] = linearrate(params[3],state[1])
    rates[4] = linearrate(params[4],state[2])
    return nothing
end