#= Script for running the ST integration and returning data arrays =#

    include("STMonteCarlo_Serial.jl")
    #include("..\\Common\\UsefulGridValueFunctions.jl")
    #include("..\\Common\\PhaseSpaceFactors.jl")
    using JLD2

function SpectraEvaluateSerial()

    # ========= Load/Create Files ========== #

        filePath = fileLocation*"\\"*fileName
        fileExist = isfile(filePath)

        if fileExist
            f = jldopen(filePath,"r+");
            SAtotal = f["STotal"];
            TAtotal = f["TTotal"];
            SAtally = f["STally"];
            TAtally = f["TTally"];
            #SMatrix = f["SMatrix"];
            #TMatrix = f["TMatrix"];
            close(f)
        else
            SAtotal = zeros(Float32,(nump3+1),numt3,nump1,numt1,nump2,numt2); 
            TAtotal = zeros(Float32,nump1,numt1,nump2,numt2);
            SAtally = zeros(UInt32,numt3,nump1,numt1,nump2,numt2);
            TAtally = zeros(UInt32,nump1,numt1,nump2,numt2);
        end

    # ====================================== #

    # ========= Pre-Allocate Arrays ======== #

        # pre-allocate arrays for momentum
        p3v = zeros(Float32,3,2); # two three array vector ((p3,t3,h1),(p3',t3',h1')) second corresponds to mirrored point in angle space
        p1v = zeros(Float32,3);
        p2v = zeros(Float32,3);

        # pre-allocate arrays for ST values
        ST = zeros(Float32,3); # [S,Sp,T]

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        STMonteCarloAxi_Serial!(SAtotal,TAtotal,SAtally,TAtally,p3v,p1v,p2v,ST)

    # ===================================== #

    # ===== Calculate S and T Matricies === #

        # preallocate
        SMatrix = zeros(Float32,(nump3+1),numt3,nump1,numt1,nump2,numt2);
        TMatrix = zeros(Float32,nump1,numt1,nump2,numt2);

        # divide element wise by tallys
        @inbounds for i in axes(SMatrix,1)
            @. @view(SMatrix[i,:,:,:,:,:]) = @view(SAtotal[i,:,:,:,:,:]) / SAtally
        end
        replace!(SMatrix,NaN=>0f0) # remove NaN caused by /0f0
        TMatrix = TAtotal ./ TAtally
        replace!(TMatrix,NaN=>0f0)

        # Angle / Momentum Ranges
        t3val = trange(t3l,t3u,numt3) # bounds of numt3 blocks
        t1val = trange(t1l,t1u,numt1)
        t2val = trange(t2l,t2u,numt2)
        p3val = prange(p3l,p3u,nump3)
        p1val = prange(p1l,p1u,nump1)
        p2val = prange(p2l,p2u,nump2)


        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix,TMatrix,t3val,p1val,t1val,p2val,t2val)    #applies phase space factors for symmetries
        STSymmetry!(SMatrix,TMatrix)                                        #initial states are symmetric -> apply symmetry of interaction to improve MC values
        PhaseSpaceFactors2!(SMatrix,TMatrix,p3val,t3val,p1val,t1val)    #corrects phase space factors for application in kinetic models
                                 
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 

    # ===================================== # 

    # ========== Save Arrays ============== #
        
        # have to delete data field and recreate cannot just update
        if fileExist    # i.e. not first time
            f = jldopen(filePath,"r+")
            Base.delete!(f,"STotal")
            Base.delete!(f,"TTotal") 
            Base.delete!(f,"STally")
            Base.delete!(f,"TTally")
            Base.delete!(f,"SMatrix")
            Base.delete!(f,"TMatrix")
            write(f,"STotal",SAtotal)
            write(f,"TTotal",TAtotal)
            write(f,"STally",SAtally)
            write(f,"TTally",TAtally)
            write(f,"SMatrix",SMatrix)
            write(f,"TMatrix",TMatrix)
        else    # create file
            f = jldopen(filePath,"w") # creates file
            write(f,"STotal",SAtotal)
            write(f,"TTotal",TAtotal)
            write(f,"STally",SAtally)
            write(f,"TTally",TAtally)
            write(f,"SMatrix",SMatrix)
            write(f,"TMatrix",TMatrix)
        end
        close(f)

        # --------- Saving Integration Parameters ------ #

        if fileExist==false # only on first time
            f = jldopen(filePath,"r+");
            write(f,"name1Data",eval(Symbol(name1*"Data")))
            write(f,"name2Data",eval(Symbol(name2*"Data")))
            write(f,"name3Data",eval(Symbol(name3*"Data")))
            write(f,"name4Data",eval(Symbol(name4*"Data")))
            close(f)
        end

        # ---------------------------------------------- #

    # ===================================== #

        return nothing

end #function

#====== Testing only =#
#SAtot = zeros(Float32,(nump3+2),numt3,nump1,numt1,nump2,numt2); 
#TAtot = zeros(Float32,nump1,numt1,nump2,numt2);
#AStal = zeros(UInt32,(nump3+2),numt3,nump1,numt1,nump2,numt2);
#ATtal = zeros(UInt32,nump1,numt1,nump2,numt2);

#=using BenchmarkTools

@benchmark STMonteCarloAxi!(SAtot,TAtot,Atal,p3v,p1v,p2v,ST)
@benchmark STMonteCarloAxi!(SAtot,TAtot,Atal,p3v,p1v,p2v,ST)

@btime STMonteCarloAxi!($SAtot,$TAtot,$AStal,$ATtal,$p3v,$p1v,$p2v,$ST)

@time STMonteCarloAxi!(SAtot,TAtot,Atal,p3v,p1v,p2v,ST)

STMonteCarloAxi!(SAtot,TAtot,Atal,p3v,p1v,p2v,ST) =#

#=====================#

#= 
function PhaseSpaceFactorstest!(SMatrix::Array{Float32,6},TMatrix::Array{Float32,4},p3val::Vector{Float32},t3val::Vector{Float32},p1val::Vector{Float32},t1val::Vector{Float32},p2val::Vector{Float32},t2val::Vector{Float32})

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays

    # Momentum space volume elements
    for ii in 1:numt2
        for jj in 1:nump2
            for kk in 1:numt1
                for ll in 1:nump1
                    for mm in 1:numt3
                        for nn in 1:nump3
                            SMatrix[nn+2,mm,ll,kk,jj,ii] *= (cospi(t3val[mm])-cospi(t3val[mm+1])) # 1/(p3val[nn+1]-p3val[nn]) # d^2pvec3
                            SMatrix[nn+2,mm,ll,kk,jj,ii] *= (cospi(t1val[kk])-cospi(t1val[kk+1]))*(p1val[ll+1]-p1val[ll])# d^3pvec3
                            SMatrix[nn+2,mm,ll,kk,jj,ii] *= (cospi(t2val[ii])-cospi(t2val[ii+1]))*(p2val[jj+1]-p2val[jj]) # d^3pvec4
                            SMatrix[nn+2,mm,ll,kk,jj,ii] *= (1f0+Float32(name3==name4))/(1f0+Float32(name1==name2))
                        end
                    end
                    TMatrix[ll,kk,jj,ii] *= (cospi(t2val[ii])-cospi(t2val[ii+1]))*(p2val[jj+1]-p2val[jj]) # d^3pvec4
                    TMatrix[ll,kk,jj,ii] *= (cospi(t1val[kk])-cospi(t1val[kk+1]))*(p1val[ll+1]-p1val[ll])
                    TMatrix[ll,kk,jj,ii] /= (1f0+Float32(name1==name2))
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
                        SMatrix[1,mm,ll,kk,jj,ii] *= (cospi(t3val[mm])-cospi(t3val[mm+1])) # 1/(p3val[1]) # d^3pvec3
                        SMatrix[1,mm,ll,kk,jj,ii] *= (cospi(t1val[kk])-cospi(t1val[kk+1]))*(p1val[ll+1]-p1val[ll]) # d^3pvec3
                        SMatrix[1,mm,ll,kk,jj,ii] *= (cospi(t2val[ii])-cospi(t2val[ii+1]))*(p2val[jj+1]-p2val[jj]) # d^3pvec4
                        SMatrix[1,mm,ll,kk,jj,ii] *= (1f0+Float32(name3==name4))/(1f0+Float32(name1==name2))
                    end
                end
            end
        end
    end

    # overflow bin size assumed to be up to 2*maximum p3val 
    for ii in 1:numt2
        for jj in 1:nump2
            for kk in 1:numt1
                for ll in 1:nump1
                    for mm in 1:numt3
                        SMatrix[2,mm,ll,kk,jj,ii] *= (cospi(t3val[mm])-cospi(t3val[mm+1])) # 1/(2*p3val[nump3+1]-p3val[nump3+1])  # d^3pvec3
                        SMatrix[2,mm,ll,kk,jj,ii] *= (cospi(t1val[kk])-cospi(t1val[kk+1]))*(p1val[ll+1]-p1val[ll])# d^3pvec3
                        SMatrix[2,mm,ll,kk,jj,ii] *= (cospi(t2val[ii])-cospi(t2val[ii+1]))*(p2val[jj+1]-p2val[jj]) # d^3pvec4
                        SMatrix[2,mm,ll,kk,jj,ii] *= (1f0+Float32(name3==name4))/(1f0+Float32(name1==name2))
                    end
                end
            end
        end
    end

end

function STCheck!(SMatrix::Array{Float32,6},TMatrix::Array{Float32,4},p3val::Vector{Float32},t3val::Vector{Float32},p1val::Vector{Float32},t1val::Vector{Float32},p2val::Vector{Float32},t2val::Vector{Float32})

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays

    # Momentum space volume elements
    for ii in 1:numt2
        for jj in 1:nump2
            for kk in 1:numt1
                for ll in 1:nump1
                    for mm in 1:numt3
                        for nn in 1:nump3
                            SMatrix[nn+2,mm,ll,kk,jj,ii] *= (cospi(t3val[mm])-cospi(t3val[mm+1]))*(p3val[nn+1]-p3val[nn]) # d^2pvec3
                        end
                    end
                    TMatrix[ll,kk,jj,ii] *= (cospi(t1val[kk])-cospi(t1val[kk+1]))*(p1val[ll+1]-p1val[ll])
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
                        SMatrix[1,mm,ll,kk,jj,ii] *= (cospi(t3val[mm])-cospi(t3val[mm+1]))*(p3val[1]) # d^3pvec3
                    end
                end
            end
        end
    end

    # overflow bin size assumed to be up to 2*maximum p3val 
    for ii in 1:numt2
        for jj in 1:nump2
            for kk in 1:numt1
                for ll in 1:nump1
                    for mm in 1:numt3
                        SMatrix[2,mm,ll,kk,jj,ii] *= (cospi(t3val[mm])-cospi(t3val[mm+1]))*(p3val[nump3+1]*10^(log10(p3val[nump3+1])-log10(p3val[nump3]))-p3val[nump3+1])  # d^3pvec3
                    end
                end
            end
        end
    end

end
 =#