# magic broadcast with @. 
# factor 3//2 due to isotropic overage of dipole
dipole_moment_sq(f, ΔE) = @. 3//2 * f / ΔE

# DOORWAYS

# "Beware of Doors." - Neverwhere

# GSB 
function doorway(::GSB, traj::Trajectory, pump::LaserPulse, ω_pump::Real, frame::Int)
    ΔE  = delta_energy(traj, frame)
    μ²  = dipole_moment_sq(traj.oscillators[:, frame], ΔE)
    Δω  = @. ω_pump - ΔE
    τ_au = pump.τ * fs2au
    E   = spectral_envelope(Δω, τ_au, pump.envelope, pump)
    return sum(@. E^2 * μ²)   # returns a scalar, no mutation
end

# SE  (active_state is 1-indexed: 1=S₀, 2=S₁, … ; energies rows match directly)
function doorway(::SE, traj::Trajectory, pump::LaserPulse, ω_pump::Real, frame::Int)
    s   = traj.active_state[frame]
    s == 1 && return 0.0                       # no SE from ground state
    ΔE  = traj.energies[s, frame] - traj.energies[1, frame]
    μ²  = dipole_moment_sq(traj.oscillators[s - 1, frame], ΔE)
    Δω  = ω_pump - ΔE
    τ_au = pump.τ * fs2au
    E   = spectral_envelope(Δω, τ_au, pump.envelope, pump)
    return E^2 * μ²   # scalar, single active state
end

# WINDOWS

# GSB
function window(::GSB, traj::Trajectory, probe::LaserPulse, ω_probe::AbstractVector)
    ΔE  = delta_energy(traj)                              # (n_states-1) x n_steps
    μ²  = dipole_moment_sq(traj.oscillators, ΔE)          # (n_states-1) x n_steps
    n_t = n_steps(traj)
    n_ω = length(ω_probe)
    τ_au = probe.τ * fs2au

    W = zeros(n_t, n_ω)
    for (j, ω) in enumerate(ω_probe)
        Δω = @. ω .- ΔE                                   # (n_states-1) x n_steps
        E  = spectral_envelope(Δω, τ_au, probe.envelope, probe)
        @views W[:, j] .= vec(sum(@. E^2 * μ²; dims=1))  # sum over states
    end
    return W
end

# SE
function window(::SE, traj::Trajectory, probe::LaserPulse, ω_probe::AbstractVector)
    n_t = n_steps(traj)
    n_ω = length(ω_probe)
    τ_au = probe.τ * fs2au

    # Extract active-state ΔE and oscillator for each frame
    ΔE_active = Vector{Float64}(undef, n_t)
    f_active  = Vector{Float64}(undef, n_t)
    for t in 1:n_t
        s = traj.active_state[t]
        if s == 1                              # on ground state → no SE
            ΔE_active[t] = 0.0
            f_active[t]  = 0.0
        else
            ΔE_active[t] = traj.energies[s, t] - traj.energies[1, t]
            f_active[t]  = traj.oscillators[s - 1, t]
        end
    end

    μ² = dipole_moment_sq(f_active, ΔE_active)
    # Guard against Inf/NaN from zero ΔE (trajectory at ground state)
    μ² .= ifelse.(isfinite.(μ²), μ², 0.0)

    W = zeros(n_t, n_ω)
    for (j, ω) in enumerate(ω_probe)
        Δω = @. ω - ΔE_active
        E  = spectral_envelope(Δω, τ_au, probe.envelope, probe)
        @views W[:, j] .= @. E^2 * μ²
    end
    return W
end

# OK, generate the spectrum!

# Convenience function for same LaserPulse in both pump and probe
function spectrum(sig::SignalType, traj::Trajectory,
                  laser::LaserPulse, ω_pump::Real, ω_probe::AbstractVector)
    spectrum(sig, traj, laser, ω_pump, laser, ω_probe)
end

"""
Pump-probe spectrum via doorway-window formalism.

The doorway is the pump absorption probability at t=0 (frame 1) — a scalar.
The window evolves at each delay time — an (n_steps, n_probe) matrix.
Signal: S(T, ω_pr) = D(0) × W(T, ω_pr).
"""
function spectrum(sig::SignalType, traj::Trajectory,
                  pump::LaserPulse, ω_pump::Real,
                  probe::LaserPulse, ω_probe::AbstractVector)
    ω_pu_au = ω_pump * eV2Ha # uhm, should have really called this E_pump
    ω_pr_au = ω_probe .* eV2Ha
    D = doorway(sig, traj, pump, ω_pu_au, 1)   # scalar at pump time (frame 1)
    W = window(sig, traj, probe, ω_pr_au)
    return D .* W
end


