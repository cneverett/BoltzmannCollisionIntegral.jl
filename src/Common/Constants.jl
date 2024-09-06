#= 
This defines all the values to be taken as constant in the evaluation of the Monte Carlo Integration of the S and T Matricies.
Including particle properties and some domain bounds
=#

# Physical constants taken from https://physics.nist.gov/cuu/Constants/index.html

const c::Float64 = Float64(299792458);              # Speed of light [m s-1]
const σT::Float64 = Float64(6.6524587e-29);         # Thompson scattering cross section [m^2]
const ħ::Float64 = Float64(1.0545718e-34);          # Reduced Planck constant [J s]
const μ0::Float64 = Float64(4*π*1e-7);              # Permeability of free space [kg m s-2 A-2]
const RSph::Float64 = Float64(1e-15);               # Sphere radius used for hard sphere collisions [m]
const q::Float64 = Float64(1.60217662e-19);        # Elementary charge [C]

# Particle Masses
const mEle::Float64 = Float64(9.10938356e-31);      # Mass of Electron [kg]
const mPos::Float64 = Float64(9.10938356e-31);      # Mass of Positron [kg]
const mPho::Float64 = Float64(0);                   # Mass of Photon [kg]
const mPro::Float64 = Float64(1.6726219e-27);       # Mass of Proton [kg]
const mSph::Float64 = Float64(1.6726219e-27);       # Mass of hard sphere [kg]

# Normalised Particle Masses (wrt electron mass)
const muEle::Float64 = Float64(1);                  # Reduced mass of Electron
const muPos::Float64 = Float64(1);                  # Reduced mass of Positron
const muPho::Float64 = Float64(0);                  # Reduced mass of Photon
const muPro::Float64 = Float64(1836.1528);          # Reduced mass of Proton
const muSph::Float64 = Float64(1836.1528);          # Reduced mass of hard sphere

# Normalised Particle Charges (wrt elementary charge)
const zEle::Float64 = Float64(-1);                  # Charge of Electron 
const zPos::Float64 = Float64(1);                   # Charge of Positron
const zPho::Float64 = Float64(0);                   # Charge of Photon
const zPro::Float64 = Float64(1);                   # Charge of Proton
const zSph::Float64 = Float64(0);                   # Charge of hard sphere

# Domain bounds
const tl::Float64 = Float64(-1);                    # Lower bound for cos(theta)
const tu::Float64 = Float64(1);                     # Upper bound for cos(theta) 