#Online 2D state means, variances, covariance, and rate-state covariances
mutable struct BalanceEqData
    meanx::Float64
    meany::Float64
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

    meanF_y::Float64

    function BalanceEqData()
        return new(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    end
    
end
#Valid update function method.
function updatestorage!(storage::BalanceEqData,state::Vector{Int64},rates::Vector{Float64},params,weight,weightsum)
    #States: [x,y]
    #Rates: [RpX,RmX,RpY,RmY]
    rw = weight/weightsum

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
    dF_Y = hill_derivative(params[1],params[5],params[6],state[2]) - storage.meanF_y

    storage.meanx += (rw)*dx
    storage.meany += (rw)*dy
    storage.meanRpX += (rw)*dRpX
    storage.meanRmX += (rw)*dRmX
    storage.meanRpY += (rw)*dRpY
    storage.meanRmY += (rw)*dRmY

    storage.meanF_y += rw*dF_Y

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

function hill_derivative(lambda::Float64,n::Float64,K::Float64,x::Int64)::Float64
    if(x==0)
        return 0
    end
    return lambda*n*x^(n-1)*K^n/(K^n+x^n)^2
end

mutable struct Balance_with_Derivative
    balancedata::BalanceEqData
    meanF_y::Float64
    function Balance_with_Derivative()
        return new(BalanceEqData(),0,0)
    end
end

function updatestorage!(storage::Balance_with_Derivative,state::Vector{Int64},rates::Vector{Float64},params,weight,weightsum)
    updatestorage!(storage.balancedata,state,rates,params,weight,weightsum)
    #params: [lambda,beta_m,gamma,beta_p,n,K]
    dF_y = hill_derivative(params[1],params[5],params[6],state[2]) - storage.meanF_y
    storage.meanF_y += (weight/weightsum)*dF_y 
    return nothing
end

mutable struct BalanceEqData_nofb
    weightsum::Float64
    meanx::Float64
    meany::Float64
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


    function BalanceEqData_nofb()
        return new(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    end
    
end

function updatestorage!(storage::BalanceEqData_nofb,state::Vector{Int64},rates::Vector{Float64},params,weight,weightsum)
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

mutable struct BalanceEqData_Fx
    weightsum::Float64
    meanx::Float64
    meany::Float64
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

    meanF_x::Float64

    function BalanceEqData_Fx()
        return new(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    end
    
end
#Valid update function method.
function updatestorage!(storage::BalanceEqData_Fx,state::Vector{Int64},rates::Vector{Float64},params,weight,weightsum)
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
    dF_x = hill_derivative(params[1],params[5],params[6],state[1]) - storage.meanF_x

    storage.meanx += (rw)*dx
    storage.meany += (rw)*dy
    storage.meanRpX += (rw)*dRpX
    storage.meanRmX += (rw)*dRmX
    storage.meanRpY += (rw)*dRpY
    storage.meanRmY += (rw)*dRmY

    storage.meanF_x += rw*dF_x

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


mutable struct BalanceEqData_Fy
    weightsum::Float64
    meanx::Float64
    meany::Float64
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

    meanF_y::Float64

    function BalanceEqData_Fy()
        return new(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    end
    
end
#Valid update function method.
function updatestorage!(storage::BalanceEqData_Fy,state::Vector{Int64},rates::Vector{Float64},params,weight,weightsum)
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
    dF_y = hill_derivative(params[1],params[5],params[6],state[2]) - storage.meanF_y

    storage.meanx += (rw)*dx
    storage.meany += (rw)*dy
    storage.meanRpX += (rw)*dRpX
    storage.meanRmX += (rw)*dRmX
    storage.meanRpY += (rw)*dRpY
    storage.meanRmY += (rw)*dRmY

    storage.meanF_y += rw*dF_y

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