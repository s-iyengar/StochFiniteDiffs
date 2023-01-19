struct BalanceEquationData2D
    #For a 2D system, stores means, variances, covariance
    #Stores mean fluxes (plus and minus) for each state
    #stores covairances of fluxes with states
    weightsum::Vector{Float64}

    means::Vector{Float64}
    #Variances: Vx,Vy,Cxy
    statevariances::Array{Float64}
    #Rpx,Rmx,Rpy,Rmy
    meanfluxes::Vector{Float64}
    #StateFluxCovs
    statefluxcovs::Array{Float64}

    #Delta states: used when updating
    dmeans::Vector{Float64}
    #DeltaRates: used when updating
    dfluxes::Vector{Float64}
    #Statediff: used when updating
    diff::Vector{Float64}
    
    function BalanceEquationData2D()
        return new(zeros(1),zeros(2),zeros(2,2),zeros(4),zeros(2,4),zeros(2),zeros(4),zeros(2))
    end
end

function updatestorage!(storage::BalanceEquationData2D,statevector::Vector{Int64},rates::Vector{Int64},params::Vector{Float64},jumptime::Float64,simtime::Float64)
    rw = jumptime/simtime

    storage.dmeans .= statevector .- storage.means

    storage.dfluxes .= rates .- storage.meanfluxes

    storage.means .+= rw * storage.dmeans
    storage.meanfluxes .+= rw * storage.dfluxes

    storage.diff .= statevector .- storage.means

    storage.statevariances .+= weight * storage.dmeans .* storage.diff'

    #Rows are first untransposed vector, columns second transposed one
    #So here, two rows with four columns, with cols Rpx,Rmx,Rpy,Rmy
    storage.statefluxcovs .+= weight * storage.diff .* storage.dfluxes'
    return nothing
end

function balance_eq_check(statefluxcovs,D,i,j)
    #picks out CiRmj and CjRmi
    LHS = statefluxcovs[i,j+1]+statefluxcovs[j,i+1]
    RHS = statefluxcovs[i,j]+statefluxcovs[j,i+1]+D[i,j]
    return relative_error(RHS,LHS)
end