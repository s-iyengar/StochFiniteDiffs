#Online 2D state means, variances, covariance, and rate-state covariances
mutable struct BalanceEqData
    meanx::Float64
    meany::Float64
    weightsum::Float64
    VX::Float64
    VY::Float64
    CXY::Float64
    meanRpX::Float64
    meanRmX::Float64
    meanRpY::Float64
    meanRmY::Float64
    CXRmX::Float64
    CXRpX::Float64
    CYRmX::Float64
    CYRpX::Float64
    CXRmY::Float64
    CXRpY::Float64
    CYRmY::Float64
    CYRpY::Float64

    function BalanceEqData()
        return new(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    end
    
end
#Valid update function method.
function updatestorage!(storage::BalanceEqData,state::Vector{Int64},rates::Vector{Float64},params,weight)
    #States: [x,y]
    #Rates: [RpX,RmX,RpY,RmY]
    storage.weightsum += weight
    rw = weight/storage.weightsum

    x = state[1]
    y = state[2]
    RpX = rates[1]
    RmX = rates[2]
    RpY = rates[3]
    RmY = rates[4]


    dx = x - storage.meanx
    dy = y - storage.meany
    dRpX = RpX - storage.meanRpX
    dRmX = RmX - storage.meanRmX
    dRpY = RpY - storage.meanRpY
    dRmY = RmY - storage.meanRmY

    storage.meanx += (rw)*dx
    storage.meany += (rw)*dy
    storage.meanRpX += (rw)*dRpX
    storage.meanRmX += (rw)*dRmX
    storage.meanRpY += (rw)*dRpY
    storage.meanRmY += (rw)*dRmY

    xdiff = x - storage.meanx
    ydiff = y - storage.meany

    storage.VX += weight*dx*xdiff
    storage.VY += weight*dy*ydiff

    storage.CXY += weight*dx*ydiff

    storage.CXRmX += weight*dRmX*xdiff
    storage.CXRpX += weight*dRpX*xdiff
    storage.CYRmX += weight*dRmX*ydiff
    storage.CYRpX += weight*dRpX*ydiff
    storage.CXRmY += weight*dRmY*xdiff
    storage.CXRpY += weight*dRpY*xdiff
    storage.CYRmY += weight*dRmY*ydiff
    storage.CYRpY += weight*dRpY*ydiff
    return nothing
end
