module Janus

# Janus.jl
# Simple implementation of Doorway-Window approach for time-dependent spectrscopy 

using DelimitedFiles
using LinearAlgebra: mul!


include("constants.jl")
export eV2Ha, Ha2eV, fs2au, au2fs, nm2eV, eV2nm, eV2invcm, invcm2eV

include("types.jl")
export SignalType, GSB, SE, ESA
export Envelope, Gaussian, Lorentzian, Sech
export Trajectory, LaserPulse, Dephasing

include("io.jl")
export load_trajectory
export delta_energy, n_states, n_steps

include("dw.jl")
export dipole_moment, spectral_envelope
export doorway, window, spectrum

end

