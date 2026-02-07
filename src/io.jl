# ZagHop ASCII format (?) currently we just have the examples with Wavemixings to test with

function load_trajectory(energy_file::AbstractString, oscillator_file::AbstractString)
    # Energy file: #header, then columns [time, cst, etot, epot, qe001, qe002, ...]
    raw = readdlm(energy_file; comments=true, comment_char='#')

    time         = @view raw[:, 1]
    active_state = Int.(@view raw[:, 2])

    # Columns 5:end are the adiabatic state energies (skip time, cst, etot, epot)
    # Transpose so rows = states, columns = time steps
    energies = permutedims(@view raw[:, 5:end])

    # Handle duplicate time steps (surface hops produce repeated entries)
    keep = _unique_timesteps(time, active_state)

    # Oscillator file: no header, columns are oscillator strengths per transition
    osc_raw = readdlm(oscillator_file; comments=true, comment_char='#')
    oscillators = permutedims(osc_raw)

    Trajectory(
        time[keep],
        active_state[keep],
        energies[:, keep],
        max.(oscillators[:, keep], 0.0),  # clamp negative oscillator strengths
    )
end

function _unique_timesteps(time, active_state)
    # When a surface hop occurs, ZagHop seems to write the same time step twice (old state
    # then new state). Keep only the last entry for each time.
    keep = trues(length(time))
    for i in 1:(length(time) - 1)
        if time[i] == time[i + 1]
            keep[i] = false
        end
    end
    return keep
end

