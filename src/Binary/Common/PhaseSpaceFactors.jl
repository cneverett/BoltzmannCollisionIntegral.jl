"""
    PhaseSpaceFactors1!(SMatrix3,SMatrix4,TMatrix,u3_bounds,u4_bounds,p1_bounds,u1_bounds,p2_bounds,u2_bounds,Indistinguishable_12)

Applies phase space volume element factors for 'SMatrix' and 'TMatrix' terms in order to correctly apply 'STSymmetry' corrections. 
"""
function PhaseSpaceFactors1!(SMatrix3::Array{Float64,9},SMatrix4::Array{Float64,9},TMatrix1::Array{Float64,6},u3_bounds::Vector{Float64},h3_bounds::Vector{Float64},u4_bounds::Vector{Float64},h4_bounds::Vector{Float64},p1_bounds::Vector{Float64},u1_bounds::Vector{Float64},h1_bounds::Vector{Float64},p2_bounds::Vector{Float64},u2_bounds::Vector{Float64},h2_bounds::Vector{Float64},Indistinguishable_12::Bool)

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays such that the symmetries can be applied.

    # Momentum space volume elements
    for h2 in axes(SMatrix3,9), u2 in axes(SMatrix3,8), p2 in axes(SMatrix3,7), h1 in axes(SMatrix3,6), u1 in axes(SMatrix3,5), p1 in axes(SMatrix3,4) # common axes
        for h3 in axes(SMatrix3,3), u3 in axes(SMatrix3,2), p3 in 1:size(SMatrix3,1)
            SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3]) # du3dh3
            SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u1_bounds[u1+1]-u1_bounds[u1])*(p1_bounds[p1+1]-p1_bounds[p1])*(h1_bounds[h1+1]-h1_bounds[h1]) # dp1du1dh1
            SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u2_bounds[u2+1]-u2_bounds[u2])*(p2_bounds[p2+1]-p2_bounds[p2])*(h2_bounds[h2+1]-h2_bounds[h2]) # dp2du2dh2
            SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (1e0+Float64(Indistinguishable_12))
        end
        for h4 in axes(SMatrix4,3), u4 in axes(SMatrix4,2), p4 in 1:size(SMatrix4,1)
            SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] *= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4]) # du4dh4
            SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] *= (u1_bounds[u1+1]-u1_bounds[u1])*(p1_bounds[p1+1]-p1_bounds[p1])*(h1_bounds[h1+1]-h1_bounds[h1]) # dp1du1dh1
            SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] *= (u2_bounds[u2+1]-u2_bounds[u2])*(p2_bounds[p2+1]-p2_bounds[p2])*(h2_bounds[h2+1]-h2_bounds[h2]) # dp2du2dh2
            SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] /= (1e0+Float64(Indistinguishable_12))
        end
        TMatrix1[p1,u1,h1,p2,u2,h2] *= (u2_bounds[u2+1]-u2_bounds[u2])*(p2_bounds[p2+1]-p2_bounds[p2])*(h2_bounds[h2+1]-h2_bounds[h2]) # dp2du2dh2
        TMatrix1[p1,u1,h1,p2,u2,h2] *= (u1_bounds[u1+1]-u1_bounds[u1])*(p1_bounds[p1+1]-p1_bounds[p1])*(h1_bounds[h1+1]-h1_bounds[h1]) # dp1du1dh1
    end

    return nothing

end

"""
    PhaseSpaceFactors2!(SMatrix3,SMatrix4,TMatrix,p3_bounds,u3_bounds,p4_bounds,u4_bounds,p1_bounds,u1_bounds,p2_bounds,u2_bounds)

To follow 'PhaseSpaceFactors1' and 'STSymmetry'. Corrects phase space factors on 'SMatrix' and 'TMatrix' for use in kinetic codes.
Assumes f(x,p,u,ϕ)= f(x,vec{p})/p^2=constant
"""
function PhaseSpaceFactors2!(SMatrix3::Array{Float64,9},SMatrix4::Array{Float64,9},TMatrix1::Array{Float64,6},TMatrix2::Array{Float64,6},p3_bounds::Vector{Float64},u3_bounds::Vector{Float64},h3_bounds::Vector{Float64},p4_bounds::Vector{Float64},u4_bounds::Vector{Float64},h4_bounds::Vector{Float64},p1_bounds::Vector{Float64},u1_bounds::Vector{Float64},h1_bounds::Vector{Float64},p2_bounds::Vector{Float64},u2_bounds::Vector{Float64},h2_bounds::Vector{Float64})

    # Function that divides the S T elements by dp3dmu3 or equivalent to then be used in kinetic models

    # === SMatrix3 === #
    # Momentum space volume elements
    for h2 in axes(SMatrix3,9), u2 in axes(SMatrix3,8), p2 in axes(SMatrix3,7), h1 in axes(SMatrix3,6), u1 in axes(SMatrix3,5), p1 in axes(SMatrix3,4) # common axes
        for h3 in axes(SMatrix3,3), u3 in axes(SMatrix3,2), p3 in axes(SMatrix3,1)
            if p3 == 1 # must account for underflow values increasing bin size 
                SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3])*(p3_bounds[p3+1])# dp3du3dh3
            elseif p3 == size(SMatrix3,1) # overflow bin size assumed to be 1*maximum p3_bounds 
                SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3])*(p3_bounds[end]) # dp3du3dh3
            else
                SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3])*(p3_bounds[p3+1]-p3_bounds[p3]) # dp3du3dh3 
            end
        end
        TMatrix1[p1,u1,h1,p2,u2,h2] /= (u1_bounds[u1+1]-u1_bounds[u1])*(h1_bounds[h1+1]-h1_bounds[h1])*(p1_bounds[p1+1]-p1_bounds[p1]) # dp1du1dh1
        TMatrix2[p2,u2,h2,p1,u1,h1] /= (u2_bounds[u2+1]-u2_bounds[u2])*(h2_bounds[h2+1]-h2_bounds[h2])*(p2_bounds[p2+1]-p2_bounds[p2]) # dp2du2dh2           
    end

    # === SMatrix4 === #
    # Momentum space volume elements
    for h2 in axes(SMatrix3,9), u2 in axes(SMatrix3,8), p2 in axes(SMatrix3,7), h1 in axes(SMatrix3,6), u1 in axes(SMatrix3,5), p1 in axes(SMatrix3,4) # common axes
        for h4 in axes(SMatrix4,3), u4 in axes(SMatrix4,2), p4 in axes(SMatrix4,1)
            if p4 == 1 # must account for underflow values increasing bin size 
                SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] /= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4])*(p4_bounds[p4+1])# dp4du4dh4 
            elseif p4 == size(SMatrix4,1) # overflow bin size assumed to be 1*maximum p4_bounds 
                SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] /= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4])*(p4_bounds[end]) # dp4du4dh4
            else
                SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] /= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4])*(p4_bounds[p4+1]-p4_bounds[p4]) # dp4du4dh4
            end 
        end
    end

    return nothing

end


"""
    STSymmetry!(SMatrix3,SMatrix4,TMatrix,u3_bounds,mu1,mu2)

To follow 'PhaseSpaceFactors1'. Physical nature of binary interaction has certain symmetries. 'STSymmetry' uses these symmetries to improve MC sampling of 'SMatrix' and 'TMatrix'.
"""
function STSymmetry!(SMatrix3::Array{Float64,9},SMatrix4::Array{Float64,9},TMatrix1::Array{Float64,6},mu1::Float64,mu2::Float64)

    # The S and T matrices are symmetric in two ways. 
    # FIRST: they are ALWAYS symmetric with respect to θ->π-θ for all particle momentum states
    # SECOND: if the incident masses are equal (mu1==mu2) then S and T are symmetric to swapping the incident particles 

    avgT = 0e0
    avgS3 = zeros(Float64,size(SMatrix3[:,:,1,1,1,1,1,1,1]))

    # SMatrix has the symmetry that if u1 u2 are mirrored in the u=0 (pi/2) plane then u3 is also mirrored in pi/2 plane generating a mirrored identical state
    # === SMatrix3 === #
    for h2 in axes(SMatrix3,9),u2 in axes(SMatrix3,8), p2 in axes(SMatrix3,7), h1 in axes(SMatrix3,6), u1 in axes(SMatrix3,5), p1 in axes(SMatrix3,4), h3 in axes(SMatrix3,3)

        sizeu1 = size(SMatrix3)[5]
        sizeu2 = size(SMatrix3)[8]
        SView3 = @view(SMatrix3[:,1:end,h3,p1,u1,h1,p2,u2,h2])
        SViewAngleMirror3 = @view(SMatrix3[:,end:-1:1,h3,p1,sizeu1-u1+1,h1,p2,sizeu2-u2+1,h2])

        if (mu1==mu2) # Both first and second Symmetry true
            SView12Swap3 = @view(SMatrix3[:,1:end,h3,p2,u2,h2,p1,u1,h1])
            SViewAngleMirror12Swap3 = @view(SMatrix3[:,end:-1:1,h3,p2,sizeu2-u2+1,h2,p1,sizeu1-u1+1,h1])
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

    avgS4 = zeros(Float64,size(SMatrix4[:,:,1,1,1,1,1,1,1]))

    # === SMatrix4 === #
    for h2 in axes(SMatrix4,9),u2 in axes(SMatrix4,8), p2 in axes(SMatrix4,7), h1 in axes(SMatrix4,6), u1 in axes(SMatrix4,5), p1 in axes(SMatrix4,4), h4 in axes(SMatrix4,3)
        sizeu1 = size(SMatrix4)[5]
        sizeu2 = size(SMatrix4)[8]
        SView4 = @view(SMatrix4[:,1:end,h4,p1,u1,h1,p2,u2,h2])
        SViewAngleMirror4 = @view(SMatrix4[:,end:-1:1,h4,p1,sizeu1-u1+1,h1,p2,sizeu2-u2+1,h2])

        if (mu1==mu2) # Both first and second Symmetry true
            SView12Swap4 = @view(SMatrix4[:,1:end,h4,p2,u2,h2,p1,u1,h1])
            SViewAngleMirror12Swap4 = @view(SMatrix4[:,end:-1:1,h4,p2,sizeu2-u2+1,h2,p1,sizeu1-u1+1,h1])
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

    # TMatrix has the symmetry that if u1 u2 are mirrored in the u=pi/2 plane generating identical mirrored states
    for h2 in axes(TMatrix1,6), u2 in axes(TMatrix1,5), p2 in axes(TMatrix1,4), h1 in axes(TMatrix1,3), u1 in axes(TMatrix1,2), p1 in axes(TMatrix1,1)
        sizeu1 = size(TMatrix1)[2]
        sizeu2 = size(TMatrix1)[5]        

        if (mu1==mu2) # Both first and second Symmetry true
            avgT = (TMatrix1[p1,u1,h1,p2,u2,h2] + TMatrix1[p2,u2,h2,p1,u1,h1] + TMatrix1[p1,sizeu1-u1+1,h1,p2,sizeu2-u2+1,h2] + TMatrix1[p2,sizeu2-u2+1,h2,p1,sizeu1-u1+1,h1])/4
            TMatrix1[p1,u1,h1,p2,u2,h2] = avgT
            TMatrix1[p2,u2,h2,p1,u1,h1] = avgT
            TMatrix1[p1,sizeu1-u1+1,h1,p2,sizeu2-u2+1,h2] = avgT
            TMatrix1[p2,sizeu2-u2+1,h2,p1,sizeu1-u1+1,h1] = avgT
        else # only first symmetry true
            avgT = (TMatrix1[p1,u1,h1,p2,u2,h2] + TMatrix1[p1,sizeu1-u1+1,h1,p2,sizeu2-u2+1,h2])/2
            TMatrix1[p1,u1,h1,p2,u2,h2] = avgT
            TMatrix1[p1,sizeu1-u1+1,h1,p2,sizeu2-u2+1,h2] = avgT
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