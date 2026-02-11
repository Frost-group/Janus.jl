using Janus
using Gnuplot

# Load trajectory data (same loader for GSB and SE!)
traj = load_trajectory("energy.dat", "oscill.dat") 
# ZagHop default filenames; so this is LZ SH, but who cares, use it for GSB!

# Laser pulse (5 fs Gaussian, no dispersion)
laser = LaserPulse(τ=5.0)

# Experimental parameters
ω_pump  = 3.0                              # eV
ω_probe = collect(range(1, 6; length=200)) # eV grid

# GSB spectrum
S_gsb = spectrum(GSB(), traj, laser, ω_pump, ω_probe)

@show S_gsb

@gp "set xlabel 'Time (steps)'" "set ylabel 'Energy (arb.)'" :-
@gp :- "set xrange [0:50]"
@gp :- "set yrange [1:200]"
@gp :- S_gsb' "with image title 'GSB'"
Gnuplot.save("gsb.png", term="pngcairo size 550,350 fontscale 0.8")

S_se = spectrum(SE(), traj, laser, ω_pump, ω_probe)

@show S_se

# Time-resolved 2D map
@gp "set xlabel 'Time (steps)'" "set ylabel 'Energy (arb.)'" :-
@gp :- "set xrange [0:50]"
@gp :- "set yrange [1:200]"
@gp :- S_se' "with image title 'SE'"
Gnuplot.save("se.png", term="pngcairo size 550,350 fontscale 0.8")

# different pump and probe lasers
pump_laser  = LaserPulse(τ=5.0, envelope=Gaussian())
probe_laser = LaserPulse(τ=10.0, envelope=Sech())
S_nd = spectrum(GSB(), traj_gs, pump_laser, ω_pump, probe_laser, ω_probe)

@gp "set xlabel 'Time (steps)'" "set ylabel 'Energy (arb.)'" :-
@gp :- "set xrange [0:50]"
@gp :- "set yrange [1:200]"
@gp :- S_nd' "with image title 'different pump and probe lasers'"
Gnuplot.save("gsb_weird_laser.png", term="pngcairo size 550,350 fontscale 0.8")

S_send = spectrum(SE(), traj_se, pump_laser, ω_pump, probe_laser, ω_probe)
@gp "set xlabel 'Time (steps)'" "set ylabel 'Energy (arb.)'" :-
@gp :- "set xrange [0:50]"
@gp :- "set yrange [1:200]"
@gp :- S_send' "with image title 'different pump and probe lasers'"
Gnuplot.save("se_weird_laser.png", term="pngcairo size 550,350 fontscale 0.8")

