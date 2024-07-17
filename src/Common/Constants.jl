#= 
This defines all the values to be taken as constant in the evaluation of the Monte Carlo Integration of the S and T Matricies.
Including particle properties and some domain bounds
=#

# Physical constants taken from https://physics.nist.gov/cuu/Constants/index.html

const c::Float32 = Float32(299792458);              # Speed of light [m s-1]
const σT::Float32 = Float32(6.6524587e-29);         # Thompson scattering cross section [m^2]
const RSph::Float32 = Float32(1e-15);               # Sphere radius used for hard sphere collisions [m]

# Particle Masses
const mEle::Float32 = Float32(9.10938356e-31);      # Mass of Electron [kg]
const mPos::Float32 = Float32(9.10938356e-31);      # Mass of Positron [kg]
const mPho::Float32 = Float32(0);                   # Mass of Photon [kg]
const mPro::Float32 = Float32(1.6726219e-27);       # Mass of Proton [kg]
const mSph::Float32 = Float32(1.6726219e-27);       # Mass of hard sphere [kg]

# Normalised Particle Masses (wrt electron mass)
const muEle::Float32 = Float32(1);                  # Reduced mass of Electron
const muPos::Float32 = Float32(1);                  # Reduced mass of Positron
const muPho::Float32 = Float32(0);                  # Reduced mass of Photon
const muPro::Float32 = Float32(1836.1528);          # Reduced mass of Proton
const muSph::Float32 = Float32(1836.1528);          # Reduced mass of hard sphere

# Domain bounds
const tl::Float32 = Float32(-1);                    # Lower bound for cos(theta)
const tu::Float32 = Float32(1);                     # Upper bound for cos(theta) 