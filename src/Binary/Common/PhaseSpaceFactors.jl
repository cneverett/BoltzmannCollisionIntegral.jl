"""
    PhaseSpaceFactors1!(SMatrix3,SMatrix4,TMatrix,u3_bounds,u4_bounds,p1_bounds,u1_bounds,p2_bounds,u2_bounds,Indistinguishable_12)

Applies phase space volume element factors for 'SMatrix' and 'TMatrix' terms in order to correctly apply 'STSymmetry' corrections. 
"""
function PhaseSpaceFactors1!(SMatrix3::Array{Float64,6},SMatrix4::Array{Float64,6},TMatrix1::Array{Float64,4},u3_bounds::Vector{Float64},u4_bounds::Vector{Float64},p1_bounds::Vector{Float64},u1_bounds::Vector{Float64},p2_bounds::Vector{Float64},u2_bounds::Vector{Float64},Indistinguishable_12::Bool)

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays such that the symmetries can be applied.


    # Momentum space volume elements
    for ii in axes(SMatrix3,6), jj in axes(SMatrix3,5), kk in axes(SMatrix3,4), ll in axes(SMatrix3,3) # common axes
        for mm in axes(SMatrix3,2), nn in 1:size(SMatrix3,1)
            SMatrix3[nn,mm,ll,kk,jj,ii] *= u3_bounds[mm+1]-u3_bounds[mm] # dmu3
            SMatrix3[nn,mm,ll,kk,jj,ii] *= (u1_bounds[kk+1]-u1_bounds[kk])*(p1_bounds[ll+1]-p1_bounds[ll]) # dp1dmu1
            SMatrix3[nn,mm,ll,kk,jj,ii] *= (u2_bounds[ii+1]-u2_bounds[ii])*(p2_bounds[jj+1]-p2_bounds[jj]) # dp2dmu2
            SMatrix3[nn,mm,ll,kk,jj,ii] /= (1e0+Float64(Indistinguishable_12))
        end
        for mm in axes(SMatrix4,2), nn in 1:size(SMatrix4,1)
            SMatrix4[nn,mm,ll,kk,jj,ii] *= u4_bounds[mm+1]-u4_bounds[mm] # dmu4
            SMatrix4[nn,mm,ll,kk,jj,ii] *= (u1_bounds[kk+1]-u1_bounds[kk])*(p1_bounds[ll+1]-p1_bounds[ll]) # dp1dmu1
            SMatrix4[nn,mm,ll,kk,jj,ii] *= (u2_bounds[ii+1]-u2_bounds[ii])*(p2_bounds[jj+1]-p2_bounds[jj]) # dp2dmu2
            SMatrix4[nn,mm,ll,kk,jj,ii] /= (1e0+Float64(Indistinguishable_12))
        end
        TMatrix1[ll,kk,jj,ii] *= (u2_bounds[ii+1]-u2_bounds[ii])*(p2_bounds[jj+1]-p2_bounds[jj]) # dp2dmu2
        TMatrix1[ll,kk,jj,ii] *= (u1_bounds[kk+1]-u1_bounds[kk])*(p1_bounds[ll+1]-p1_bounds[ll]) # dp1dmu1
    end

    return nothing

end

"""
    PhaseSpaceFactors2!(SMatrix3,SMatrix4,TMatrix,p3_bounds,u3_bounds,p4_bounds,u4_bounds,p1_bounds,u1_bounds,p2_bounds,u2_bounds)

To follow 'PhaseSpaceFactors1' and 'STSymmetry'. Corrects phase space factors on 'SMatrix' and 'TMatrix' for use in kinetic codes.
Assumes f(x,p,μ)= constant
"""
function PhaseSpaceFactors2!(SMatrix3::Array{Float64,6},SMatrix4::Array{Float64,6},TMatrix1::Array{Float64,4},TMatrix2::Array{Float64,4},p3_bounds::Vector{Float64},u3_bounds::Vector{Float64},p4_bounds::Vector{Float64},u4_bounds::Vector{Float64},p1_bounds::Vector{Float64},u1_bounds::Vector{Float64},p2_bounds::Vector{Float64},u2_bounds::Vector{Float64})

    # Function that divides the S T elements by dp3dmu3 or equivalent to then be used in kinetic models

    # === SMatrix3 === #
    # Momentum space volume elements
    for ii in axes(SMatrix3,6), jj in axes(SMatrix3,5), kk in axes(SMatrix3,4), ll in axes(SMatrix3,3)
        for mm in axes(SMatrix3,2), nn in 1:size(SMatrix3,1)-1
            if nn == 1 # must account for underflow values increasing bin size 
                SMatrix3[nn,mm,ll,kk,jj,ii] /= (u3_bounds[mm+1]-u3_bounds[mm])*(p3_bounds[nn+1])# dp3dmu3 
            else
                SMatrix3[nn,mm,ll,kk,jj,ii] /= (u3_bounds[mm+1]-u3_bounds[mm])*(p3_bounds[nn+1]-p3_bounds[nn]) # dp3dmu3 
            end
        end
        TMatrix1[ll,kk,jj,ii] /= (u1_bounds[kk+1]-u1_bounds[kk])*(p1_bounds[ll+1]-p1_bounds[ll]) # dp1dmu1
        TMatrix2[jj,ii,ll,kk] /= (u2_bounds[ii+1]-u2_bounds[ii])*(p2_bounds[jj+1]-p2_bounds[jj]) # dp2dmu2           
    end

    # overflow bin size assumed to be up to 1*maximum p3_bounds 
    for ii in axes(SMatrix3,6), jj in axes(SMatrix3,5), kk in axes(SMatrix3,4), ll in axes(SMatrix3,3), mm in axes(SMatrix3,2)
        SMatrix3[end,mm,ll,kk,jj,ii] /= (u3_bounds[mm+1]-u3_bounds[mm])*(p3_bounds[end]) # dp3dmu3
    end

    # === SMatrix4 === #
    # Momentum space volume elements
    for ii in axes(SMatrix4,6), jj in axes(SMatrix4,5), kk in axes(SMatrix4,4), ll in axes(SMatrix4,3)
        for mm in axes(SMatrix4,2), nn in 1:size(SMatrix4,1)-1
            if nn == 1 # must account for underflow values increasing bin size 
                SMatrix4[nn,mm,ll,kk,jj,ii] /= (u4_bounds[mm+1]-u4_bounds[mm])*(p4_bounds[nn+1])# dp3dmu3 
            else
                SMatrix4[nn,mm,ll,kk,jj,ii] /= (u4_bounds[mm+1]-u4_bounds[mm])*(p4_bounds[nn+1]-p4_bounds[nn]) # dp3dmu3 
            end 
        end
    end

    # overflow bin size assumed to be up to 1*maximum p4_bounds 
    for ii in axes(SMatrix4,6), jj in axes(SMatrix4,5), kk in axes(SMatrix4,4), ll in axes(SMatrix4,3), mm in axes(SMatrix4,2)
        SMatrix4[end,mm,ll,kk,jj,ii] /= (u4_bounds[mm+1]-u4_bounds[mm])*(p4_bounds[end]) # dp3dmu3
    end

    return nothing

end


"""
    STSymmetry!(SMatrix3,SMatrix4,TMatrix,u3_bounds,mu1,mu2)

To follow 'PhaseSpaceFactors1'. Physical nature of binary interaction has certain symmetries. 'STSymmetry' uses these symmetries to improve MC sampling of 'SMatrix' and 'TMatrix'.
"""
function STSymmetry!(SMatrix3::Array{Float64,6},SMatrix4::Array{Float64,6},TMatrix1::Array{Float64,4},mu1::Float64,mu2::Float64)

    # The S and T matricies are symmetric in two ways. 
    # FIRST: they are ALWAYS symmetric with respect to θ->π-θ for all particle momentum states
    # SECOND: if the incident masses are equal (mu1==mu2) then S and T are symmetric to swapping the incident particles 

    avgT = 0e0
    avgS3 = zeros(Float64,size(SMatrix3[:,:,1,1,1,1]))

    # SMatrix has the symmetry that if t1 t2 are mirrored in the t=pi/2 plane then t3 is also mirrored in pi/2 plane genreting a mirrored identical state
    # === SMatrix3 === #
    for ii in axes(SMatrix3,6), jj in axes(SMatrix3,5), kk in axes(SMatrix3,4), ll in axes(SMatrix3,3)
        sizet1 = size(SMatrix3)[4]
        sizet2 = size(SMatrix3)[6]
        SView3 = @view(SMatrix3[:,1:end,ll,kk,jj,ii])
        SViewAngleMirror3 = @view(SMatrix3[:,end:-1:1,ll,sizet1-kk+1,jj,sizet2-ii+1])

        if (mu1==mu2) # Both first and second Symmetry true
            SView12Swap3 = @view(SMatrix3[:,1:end,jj,ii,ll,kk])
            SViewAngleMirror12Swap3 = @view(SMatrix3[:,end:-1:1,jj,sizet2-ii+1,ll,sizet1-kk+1])
            @. avgS3 = (SView3 + SViewAngleMirror3 + SView12Swap3 + SViewAngleMirror12Swap3)/4
            @. SView3 = avgS3
            @. SView12Swap3 = avgS3
            @. SViewAngleMirror3 = avgS3
            @. SViewAngleMirror12Swap3 = avgS3
        else # only first symmetry true
            @. avgS3 = (SView3 + SViewAngleMirror3)/2
            @. SView3 = avgS3
            @. SViewAngleMirror3 = avgS3
        end
        
    end 

    avgS4 = zeros(Float64,size(SMatrix4[:,:,1,1,1,1]))

    # === SMatrix4 === #
    for ii in axes(SMatrix4,6), jj in axes(SMatrix4,5), kk in axes(SMatrix4,4), ll in axes(SMatrix4,3)
        sizet1 = size(SMatrix4)[4]
        sizet2 = size(SMatrix4)[6]
        SView4 = @view(SMatrix4[:,1:end,ll,kk,jj,ii])
        SViewAngleMirror4 = @view(SMatrix4[:,end:-1:1,ll,sizet1-kk+1,jj,sizet2-ii+1])

        if (mu1==mu2) # Both first and second Symmetry true
            SView12Swap4 = @view(SMatrix4[:,1:end,jj,ii,ll,kk])
            SViewAngleMirror12Swap4 = @view(SMatrix4[:,end:-1:1,jj,sizet2-ii+1,ll,sizet1-kk+1])
            @. avgS4 = (SView4 + SViewAngleMirror4 + SView12Swap4 + SViewAngleMirror12Swap4)/4
            @. SView4 = avgS4
            @. SView12Swap4 = avgS4
            @. SViewAngleMirror4 = avgS4
            @. SViewAngleMirror12Swap4 = avgS4
        else # only first symmetry true
            @. avgS4 = (SView4 + SViewAngleMirror4)/2
            @. SView4 = avgS4
            @. SViewAngleMirror4 = avgS4
        end
             
    end

    # TMatrix has the symmetry that if t1 t2 are mirrored in the t=pi/2 plane generating identical mirrored states
    for ii in axes(TMatrix1,4), jj in axes(TMatrix1,3), kk in axes(TMatrix1,2), ll in axes(TMatrix1,1)
        sizet1 = size(TMatrix1)[2]
        sizet2 = size(TMatrix1)[4]        

        if (mu1==mu2) # Both first and second Symmetry true
            avgT = (TMatrix1[ll,kk,jj,ii] + TMatrix1[jj,ii,ll,kk] + TMatrix1[ll,sizet1-kk+1,jj,sizet2-ii+1] + TMatrix1[jj,sizet2-ii+1,ll,sizet1-kk+1])/4
            TMatrix1[ll,kk,jj,ii] = avgT
            TMatrix1[jj,ii,ll,kk] = avgT
            TMatrix1[ll,sizet1-kk+1,jj,sizet2-ii+1] = avgT
            TMatrix1[jj,sizet2-ii+1,ll,sizet1-kk+1] = avgT
        else # only first symmetry true
            avgT = (TMatrix1[ll,kk,jj,ii] + TMatrix1[ll,sizet1-kk+1,jj,sizet2-ii+1])/2
            TMatrix1[ll,kk,jj,ii] = avgT
            TMatrix1[ll,sizet1-kk+1,jj,sizet2-ii+1] = avgT
        end
    end 
    
    return nothing

end

# ====== Currently Unused ================ #

    function SCorrection!(SMatrix::Array{Float64,6},TMatrix::Array{Float64,4},p3_bounds::Vector{Float64},u3_bounds::Vector{Float64},p1_bounds::Vector{Float64},u1_bounds::Vector{Float64},p2_bounds::Vector{Float64},u2_bounds::Vector{Float64})

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
                                SFull[nn+2,mm,ll,kk,jj,ii] = SMatrix[nn+2,mm,ll,kk,jj,ii] *(cospi(u3_bounds[mm])-cospi(u3_bounds[mm+1]))*(p3_bounds[nn+1]-p3_bounds[nn]) # d^2pvec3
                            end
                        end
                        TFull[ll,kk,jj,ii] = TMatrix[ll,kk,jj,ii] * (cospi(u1_bounds[kk])-cospi(u1_bounds[kk+1]))*(p1_bounds[ll+1]-p1_bounds[ll])
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
                            SFull[1,mm,ll,kk,jj,ii] = SMatrix[1,mm,ll,kk,jj,ii] * (cospi(u3_bounds[mm])-cospi(u3_bounds[mm+1]))*(p3_bounds[1]) # d^3pvec3
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
                            SFull[2,mm,ll,kk,jj,ii] = SMatrix[2,mm,ll,kk,jj,ii] * (cospi(u3_bounds[mm])-cospi(u3_bounds[mm+1])) * (p3_bounds[nump3+1]*10^(log10(p3_bounds[nump3+1])-log10(p3_bounds[nump3]))-p3_bounds[nump3+1])  # d^3pvec3
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

    function SCorrection2!(SMatrix::Array{Float64,6},TMatrix::Array{Float64,4})

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