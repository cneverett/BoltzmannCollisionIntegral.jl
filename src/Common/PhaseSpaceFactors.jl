"""
    PhaseSpaceFactors1!(SMatrix,TMatrix,t3val,p1val,t1val,p2val,t2val,name1,name2)

Applies phase space volume element factors for 'SMatrix' and 'TMatrix' terms in order to correctly apply 'STSymmetry' corrections. 
"""
function PhaseSpaceFactors1!(SMatrix::Array{Float32,6},TMatrix::Array{Float32,4},t3val::Vector{Float32},p1val::Vector{Float32},t1val::Vector{Float32},p2val::Vector{Float32},t2val::Vector{Float32},name1::String,name2::String)

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays such that the symmetries can be applied.

    # Momentum space volume elements
    for ii in axes(SMatrix,6), jj in axes(SMatrix,5), kk in axes(SMatrix,4), ll in axes(SMatrix,3)
        for mm in axes(SMatrix,2), nn in 1:size(SMatrix,1)-1
            SMatrix[nn,mm,ll,kk,jj,ii] *= t3val[mm+1]-t3val[mm] # dmu3
            SMatrix[nn,mm,ll,kk,jj,ii] *= (t1val[kk]-t1val[kk+1])*(p1val[ll+1]-p1val[ll]) # dp1dmu1
            SMatrix[nn,mm,ll,kk,jj,ii] *= (t2val[ii]-t2val[ii+1])*(p2val[jj+1]-p2val[jj]) # dp2dmu2
            SMatrix[nn,mm,ll,kk,jj,ii] /= (1f0+Float32(name1==name2))
        end
        TMatrix[ll,kk,jj,ii] *= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) # dp2dmu2
        TMatrix[ll,kk,jj,ii] *= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) # dp1dmu1
    end

    # underflow bin size
    #= 
    for ii in 1:numt2, jj in 1:nump2, kk in 1:numt1, ll in 1:nump1, mm in 1:numt3
        SMatrix[1,mm,ll,kk,jj,ii] *= (t3val[mm+1]-t3val[mm]) # dmu3
        SMatrix[1,mm,ll,kk,jj,ii] *= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) # dp1dmu1
        SMatrix[1,mm,ll,kk,jj,ii] *= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) # dp2dmu2
        SMatrix[1,mm,ll,kk,jj,ii] /= (1f0+Float32(name1==name2))
    end
    =#

    # overflow bin size assumed to be up to 1*maximum p3val 
    for ii in axes(SMatrix,6), jj in axes(SMatrix,5), kk in axes(SMatrix,4), ll in axes(SMatrix,3), mm in axes(SMatrix,2)
        SMatrix[end,mm,ll,kk,jj,ii] *= (t3val[mm+1]-t3val[mm]) # dmu3
        SMatrix[end,mm,ll,kk,jj,ii] *= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) # dp1dmu1
        SMatrix[end,mm,ll,kk,jj,ii] *= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) # dp2dmu2
        SMatrix[end,mm,ll,kk,jj,ii] /= (1f0+Float32(name1==name2))
    end

    return nothing

end

"""
    PhaseSpaceFactors2!(SMatrix,TMatrix,t3val,p1val,t1val,p2val,t2val)

To follow 'PhaseSpaceFactors1' and 'STSymmetry'. Corrects phase space factors on 'SMatrix' and 'TMatrix' for use in kinetic codes.
Assumes f(x,p,μ)= constant
"""
function PhaseSpaceFactors2!(SMatrix::Array{Float32,6},TMatrix::Array{Float32,4},p3val::Vector{Float32},t3val::Vector{Float32},p1val::Vector{Float32},t1val::Vector{Float32})

    # Function that divides the S T elements by dp3dmu3 or equivilant to then be used in kinetic models

    # Momentum space volume elements
    for ii in axes(SMatrix,6), jj in axes(SMatrix,5), kk in axes(SMatrix,4), ll in axes(SMatrix,3)
        for mm in axes(SMatrix,2), nn in 1:size(SMatrix,1)-1
            if nn == 1 # must account for underflow values increasing bin size (perhaps not)
            SMatrix[nn,mm,ll,kk,jj,ii] /= (t3val[mm+1]-t3val[mm])*(p3val[nn+1]-p3val[nn])# dp3dmu3 
            else
            SMatrix[nn,mm,ll,kk,jj,ii] /= (t3val[mm+1]-t3val[mm])*(p3val[nn+1]-p3val[nn]) # dp3dmu3 
            end
        end
        TMatrix[ll,kk,jj,ii] /= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) # dp1dmu1      
    end

    # underflow bin size
    #= 
    for ii in 1:numt2, jj in 1:nump2, kk in 1:numt1, ll in 1:nump1, mm in 1:numt3
        SMatrix[1,mm,ll,kk,jj,ii] /= (t3val[mm+1]-t3val[mm])*(p3val[1]) # dp3dmu3
    end
    =#

    # overflow bin size assumed to be up to 1*maximum p3val 
    for ii in axes(SMatrix,6), jj in axes(SMatrix,5), kk in axes(SMatrix,4), ll in axes(SMatrix,3), mm in axes(SMatrix,2)
        SMatrix[end,mm,ll,kk,jj,ii] /= (t3val[mm+1]-t3val[mm])*(p3val[end]) # dp3dmu3
    end

    return nothing

end


"""
    STSymmetry!(SMatrix,TMatrix,t3val,p1val,t1val,p2val,t2val,mu1,mu2)

To follow 'PhaseSpaceFactors1'. Physical nature of binary interaction has certain symmetries. 'STSymmetry' uses these symmetries to improve MC sampling of 'SMatrix' and 'TMatrix'.
"""
function STSymmetry!(SMatrix::Array{Float32,6},TMatrix::Array{Float32,4},mu1::Float32,mu2::Float32)

    # The S and T matricies are symmetric in two ways. 
    # FIRST: they are ALWAYS symmetric with respect to θ->π-θ for all particle momentum states
    # SECOND: if the incident masses are equal (mu1==mu2) then S and T are symmetric to swapping the incident particles 

    avgT = 0f0
    avgS = zeros(Float32,size(SMatrix[:,:,1,1,1,1]))

    # SMatrix has the symmetry that if t1 t2 are mirrored in the t=pi/2 plane then t3 is also mirrored in pi/2 plane genreting a mirrored identical state
    for ii in axes(SMatrix,6), jj in axes(SMatrix,5), kk in axes(SMatrix,4), ll in axes(SMatrix,3)
        sizet1 = size(SMatrix)[4]
        sizet2 = size(SMatrix)[6]
        SView = @view(SMatrix[:,1:end,ll,kk,jj,ii])
        SViewAngleMirror = @view(SMatrix[:,end:-1:1,ll,sizet1-kk+1,jj,sizet2-ii+1])

        if (mu1==mu2) # Both first and second Symmetry true
            SView12Swap = @view(SMatrix[:,1:end,jj,ii,ll,kk])
            SViewAngleMirror12Swap = @view(SMatrix[:,end:-1:1,jj,sizet2-ii+1,ll,sizet1-kk+1])
            @. avgS = (SView + SViewAngleMirror + SView12Swap + SViewAngleMirror12Swap)/4
            @. SView = avgS
            @. SView12Swap = avgS
            @. SViewAngleMirror = avgS
            @. SViewAngleMirror12Swap = avgS
        else # only first symmetry true
            @. avgS = (SView + SViewAngleMirror)/2
            @. SView = avgS
            @. SViewAngleMirror = avgS
        end
        
    end 

    # TMatrix has the symmetry that if t1 t2 are mirrored in the t=pi/2 plane generating identical mirrored states
    for ii in axes(TMatrix,4), jj in axes(TMatrix,3), kk in axes(TMatrix,2), ll in axes(TMatrix,1)
        sizet1 = size(TMatrix)[2]
        sizet2 = size(TMatrix)[4]        

        if (mu1==mu2) # Both first and second Symmetry true
            avgT = (TMatrix[ll,kk,jj,ii] + TMatrix[jj,ii,ll,kk] + TMatrix[ll,sizet1-kk+1,jj,sizet2-ii+1] + TMatrix[jj,sizet2-ii+1,ll,sizet1-kk+1])/4
            TMatrix[ll,kk,jj,ii] = avgT
            TMatrix[jj,ii,ll,kk] = avgT
            TMatrix[ll,sizet1-kk+1,jj,sizet2-ii+1] = avgT
            TMatrix[jj,sizet2-ii+1,ll,sizet1-kk+1] = avgT
        else # only first symmetry true
            avgT = (TMatrix[ll,kk,jj,ii] + TMatrix[ll,sizet1-kk+1,jj,sizet2-ii+1])/2
            TMatrix[ll,kk,jj,ii] = avgT
            TMatrix[ll,sizet1-kk+1,jj,sizet2-ii+1] = avgT
        end
    end 
    
    return nothing

end

# ====== Currently Unused ================ #

    function SCorrection!(SMatrix::Array{Float32,6},TMatrix::Array{Float32,4},p3val::Vector{Float32},t3val::Vector{Float32},p1val::Vector{Float32},t1val::Vector{Float32},p2val::Vector{Float32},t2val::Vector{Float32})

        # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays

        SFull = zeros(size(SMatrix))
        TFull = zeros(size(TMatrix))

        # Momentum space volume elements
        for ii in 1:numt2
            for jj in 1:nump2
                for kk in 1:numt1
                    for ll in 1:nump1
                        for mm in 1:numt3
                            for nn in 1:nump3
                                SFull[nn+2,mm,ll,kk,jj,ii] = SMatrix[nn+2,mm,ll,kk,jj,ii] *(cospi(t3val[mm])-cospi(t3val[mm+1]))*(p3val[nn+1]-p3val[nn]) # d^2pvec3
                            end
                        end
                        TFull[ll,kk,jj,ii] = TMatrix[ll,kk,jj,ii] * (cospi(t1val[kk])-cospi(t1val[kk+1]))*(p1val[ll+1]-p1val[ll])
                    end
                end
            end
        end

        # underflow bin size 
        for ii in 1:numt2
            for jj in 1:nump2
                for kk in 1:numt1
                    for ll in 1:nump1
                        for mm in 1:numt3
                            SFull[1,mm,ll,kk,jj,ii] = SMatrix[1,mm,ll,kk,jj,ii] * (cospi(t3val[mm])-cospi(t3val[mm+1]))*(p3val[1]) # d^3pvec3
                        end
                    end
                end
            end
        end

        # overflow bin size assumed to be up to the next pval (if one exsisted)
        for ii in 1:numt2
            for jj in 1:nump2
                for kk in 1:numt1
                    for ll in 1:nump1
                        for mm in 1:numt3
                            SFull[2,mm,ll,kk,jj,ii] = SMatrix[2,mm,ll,kk,jj,ii] * (cospi(t3val[mm])-cospi(t3val[mm+1])) * (p3val[nump3+1]*10^(log10(p3val[nump3+1])-log10(p3val[nump3]))-p3val[nump3+1])  # d^3pvec3
                        end
                    end
                end
            end
        end

        testT = TFull
        testS = dropdims(sum(SFull,dims=(1,2)),dims=(1,2))
        SCor = zeros(size(testS))
        @. SCor = testT/testS

        #display(dropdims(sum(SCor,dims=(2,4)),dims=(2,4)))

        for ii in axes(SCor,4), jj in axes(SCor,3), kk in axes(SCor,2), ll in axes(SCor,1)
            SMatrix[:,:,ll,kk,jj,ii] *= SCor[ll,kk,jj,ii]
        end

    end

    function SCorrection2!(SMatrix::Array{Float32,6},TMatrix::Array{Float32,4})

        # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays

        SMatrixSum = dropdims(sum(SMatrix,dims=(1,2)),dims=(1,2))
        SCor = zeros(size(SMatrixSum))
        @. SCor = TMatrix/SMatrixSum

        #display(dropdims(sum(SCor,dims=(2,4)),dims=(2,4)))

        for ii in axes(SCor,4), jj in axes(SCor,3), kk in axes(SCor,2), ll in axes(SCor,1)
            SMatrix[:,:,ll,kk,jj,ii] *= SCor[ll,kk,jj,ii]
        end

    end

# ======================================== #