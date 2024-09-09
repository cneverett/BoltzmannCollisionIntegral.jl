"""
    PhaseSpaceFactorsSync1!(SMatrix,p1val,t1val,p2val,t2val)

Applies phase space volume element factors for 'SMatrix' terms in order to correctly apply 'SyncSymmetry' corrections. 
"""
function PhaseSpaceFactorsSync1!(SMatrix::Array{Float64,4},p1val::Vector{Float64},t1val::Vector{Float64},p2val::Vector{Float64},t2val::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(Smatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] *= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) #dp1dmu1 
        SMatrix[ll,kk,jj,ii] *= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) #dp2dmu2
    end

end # function

"""
    PhaseSpaceFactorsSync2!(SMatrix,p1val,t1val)

To follow 'PhaseSpaceFactorsSync1' and 'SyncSymmetry'. Correct phase spcae factors on 'SMatrix' for use in kinetic codes. 
"""
function PhaseSpaceFactorsSync2!(SMatrix::Array{Float64,4},p1val::Vector{Float64},t1val::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(Smatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] /= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) #dp1dmu1 
    end

end # function

"""
    SyncSymmetry!(SMatrix)

To follow 'PhaseSpaceFactorsSync1'. Synchrotron emission has a symmetry with respect to cos(theta) -> -cos(theta) for both initial particle and photon momenta.
"""
function SyncSymmetry!(SMatrix::Array{Float64,4})

    avgS = zeros(Float64,size(SMatrix))
    @. avgS = (SMatrix[:,end:-1:1,:,end:-1:1] + SMatrix)/2
    SMatrix .= avgS

end # function