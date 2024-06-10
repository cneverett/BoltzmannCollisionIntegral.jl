
# ================ Struct for Particle Data ================= #

    struct PD
        name::String        # three letter abbreviation

        # Momentum data
        pl::Float32         # lowest momentum value (log10 space)
        pu::Float32         # highest momentum value (log10 space)
        nump::Int64         # number of mometum divisions/grid cells
        
        # Angle Data
        tl::Float32         # lower angular bound (cos(theta) space)
        tu::Float32         # upper angular bound (cos(theta) space)
        numt::Int64         # number of angular divisions/grid cells

        # Physical Data
        mass::Float32          # particle mass in kg
        normmass::Float32         # reduced particle mass (wrt electron mass)
    end

# =========================================================== #

# ====================== Assign Data ======================== # 

    # Hard Sphere
    pl = -5f0; pu = 4f0; nump = 72;
    tl = -1f0; tu = 1f0; numt = 8;
    SphData = PD("Sph",pl,pu,nump,tl,tu,numt,1.672622f-27,1836.1528f0)

    # Electron
    EleData = PD("Ele",pl,pu,nump,tl,tu,numt,9.109383f-31,1f0)
    # Positron
    PosData = PD("Pos",pl,pu,nump,tl,tu,numt,9.109383f-31,1f0)
    # Photon 
    PhoData = PD("Pho",pl,pu,nump,tl,tu,numt,0f0,0f0)

    #ProData = PD("Pro",-5f0,4f0,36,[range(-5f0,4f0,36+1);],0f0,1f0,8,[0f0,1f0,8+1],1.672622f-27,1836.1528f0)

# =========================================================== # 

