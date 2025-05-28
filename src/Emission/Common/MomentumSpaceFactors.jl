"""
    PhaseSpaceFactorsSync1!(SMatrix,p1val,t1val,p2val,t2val)

Applies phase space volume element factors for 'SMatrix' terms in order to correctly apply 'SyncSymmetry' corrections. 
"""
function PhaseSpaceFactorsSync1!(SMatrix::Array{Float64,4},p1val::Vector{Float64},t1val::Vector{Float64},p2val::Vector{Float64},t2val::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(SMatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] *= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) #dp1dmu1 
        SMatrix[ll,kk,jj,ii] *= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) #dp2dmu2
    end

end # function

"""
    PhaseSpaceFactorsSync2!(SMatrix,p1val,t1val)

To follow 'PhaseSpaceFactorsSync1' and 'SyncSymmetry'. Correct phase space factors on 'SMatrix' for use in kinetic codes. 
"""
function PhaseSpaceFactorsSync2!(SMatrix::Array{Float64,4},p1val::Vector{Float64},t1val::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(SMatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] /= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) #dp1dmu1 
    end

end # function

function MomentumSpaceFactorsEmission!(LossMatrix1,GainMatrix2,GainMatrix3::Array{Float64,6},p3val::Vector{Float64},u3val::Vector{Float64},h3val::Vector{Float64})

    for h1 in axes(GainMatrix3,6), u1 in axes(GainMatrix3,5), p1 in axes(GainMatrix3,4), h3 in axes(GainMatrix3,3), u3 in axes(GainMatrix3,2), p3 in axes(GainMatrix3,1)
        GainMatrix3[p3,u3,h3,p1,u1,h1] *= (u3val[u3+1]-u3val[u3])*(p3val[p3+1]-p3val[p3])*(h3val[h3+1]-h3val[h3]) #dp3du3dh3
    end

    ## TO BE ADDED GAINMATRIX2 and LOSSMATRIX1

end

"""
    SyncSymmetry!(SMatrix)

To follow 'PhaseSpaceFactorsSync1'. Synchrotron emission has a symmetry with respect to cos(theta) -> -cos(theta) for both initial particle and photon momenta.
"""
function SymmetryEmission!(TMatrix1::Array{Float64,3},SMatrix2::Array{Float64,6},SMatrix3::Array{Float64,6})

    avgS = zeros(Float64,size(SMatrix2))
    @. avgS = (SMatrix2[:,end:-1:1,:,:,end:-1:1,:] + SMatrix2)/2
    SMatrix2 .= avgS

    avgS = zeros(Float64,size(SMatrix3))
    @. avgS = (SMatrix3[:,end:-1:1,:,:,end:-1:1,:] + SMatrix3)/2
    SMatrix3 .= avgS

end # function

function GainLossSymmetryEmission!(GainTotal2,GainTotal3,GainTally2,GainTally3,LossTotal1,LossTally1)

    GainTotal2Mirror = @view(GainTotal2[:,end:-1:1,:,:,end:-1:1,:])
    GainTotal3Mirror = @view(GainTotal3[:,end:-1:1,:,:,end:-1:1,:])
    LossTotal1Mirror = @view(LossTotal1[:,end:-1:1,:,:,end:-1:1,:])
    GainTally2Mirror = @view(GainTally2[:,end:-1:1,:,:,end:-1:1,:])
    GainTally3Mirror = @view(GainTally3[:,end:-1:1,:,:,end:-1:1,:])
    LossTally1Mirror = @view(LossTally1[:,end:-1:1,:,:,end:-1:1,:])

    @. GainTotal2 += GainTotal2Mirror
    @. GainTotal3 += GainTotal3Mirror
    @. LossTotal1 += LossTotal1Mirror
    @. GainTally2 += GainTally2Mirror
    @. GainTally3 += GainTally3Mirror
    @. LossTally1 += LossTally1Mirror


end # function