# Abstract types for dispatch

abstract type SignalType end
struct GSB <: SignalType end   # Ground State Bleach
struct SE  <: SignalType end   # Stimulated Emission
struct ESA <: SignalType end   # Excited State Absorption

abstract type Envelope end
struct Gaussian   <: Envelope end
struct Lorentzian <: Envelope end
struct Sech       <: Envelope end



"""
Trajectory data from a non-adiabatic (surface-hopping kind presumed) simulation.

Fields:
- `time::Vector{Float64}`: simulation times (fs)
- `active_state::Vector{Int}`: occupied electronic state (density?) at each step
- `energies::Matrix{Float64}`: adiabatic energies (n_states x n_steps), Hartree
- `oscillators::Matrix{Float64}`: oscillator strengths ((n_states-1) x n_steps)
"""
struct Trajectory
    time::Vector{Float64}
    active_state::Vector{Int}
    energies::Matrix{Float64}
    oscillators::Matrix{Float64}
end

n_states(t::Trajectory) = size(t.energies, 1)
n_steps(t::Trajectory)  = size(t.energies, 2)

"""ΔE_e = E_e - E_ground for all excited states at one frame."""
delta_energy(traj::Trajectory, frame::Int) =
    @views traj.energies[2:end, frame] .- traj.energies[1, frame]

"""ΔE matrix for all frames: (n_states-1) x n_steps."""
delta_energy(traj::Trajectory) =
    @views traj.energies[2:end, :] .- traj.energies[1:1, :]


"""
Laser pulse characteristics (envelope shape, duration, dispersion).
(Nb: central frequency passed elsewhere; this is just arb pulse!)

Fields:
- `envelope::E`: spectral envelope shape (Gaussian, Lorentzian, Sech, etc.)
- `τ::Float64`: pulse duration (fs)
- `D2`--`D5`: dispersion coefficients (fs^n), default 0
"""
@kwdef struct LaserPulse{E<:Envelope}
    envelope::E = Gaussian()
    τ::Float64
    D2::Float64 = 0.0
    D3::Float64 = 0.0
    D4::Float64 = 0.0
    D5::Float64 = 0.0
end

"""Dephasing rate for dispersed and 2D calculations."""
@kwdef struct Dephasing
    γ::Float64 = 0.01  # eV
end
# not used yet

# Nb: spectra envelope stuff needs to be down here, as depends on LaserPulse having been ready by Julia
"""
    spectral_envelope(Δω, τ_au, envelope, pulse)

Spectral field envelope E(Δω, τ) in the frequency domain.

Returns the field amplitude. Signal formulas use E² (proportional to intensity).
All envelopes include τ_au prefactor for normalisation. 
"""
function spectral_envelope(Δω, τ_au::Real, ::Gaussian, p::LaserPulse)
    @. τ_au * exp(
        -(Δω * τ_au)^2 / 4 - (
            _inv_fact2 * p.D2 * Δω^2 +
            _inv_fact3 * p.D3 * Δω^3 +
            _inv_fact4 * p.D4 * Δω^4 +
            _inv_fact5 * p.D5 * Δω^5
        )
    )
end

# Lorentzian
function spectral_envelope(Δω, τ_au::Real, ::Lorentzian, ::LaserPulse)
    @. τ_au * exp(-2 * abs(Δω) * τ_au)
end

# Sech
function spectral_envelope(Δω, τ_au::Real, ::Sech, ::LaserPulse)
    @. τ_au * sech(π * Δω * τ_au / 2)
end
