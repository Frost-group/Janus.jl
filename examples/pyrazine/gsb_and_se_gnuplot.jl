using Janus
using Gnuplot

# Load trajectory data (same loader for GSB and SE!)
traj_gs = load_trajectory("gsb/gsb_ene.dat", "gsb/gsb_osc.dat")

# Laser pulse (5 fs Gaussian, no dispersion)
laser = LaserPulse(τ=5.0)

# Experimental parameters
ω_pump  = 3.0                             # eV
ω_probe = collect(range(1, 6; length=200)) # eV grid

# GSB spectrum
S_gsb = spectrum(GSB(), traj_gs, laser, ω_pump, ω_probe)

@show S_gsb

@gp "set xlabel 'Probe energy (eV)'" "set ylabel 'Signal (arb.)'" :-
@gp :- "set xrange [0:50]"
@gp :- "set yrange [1:200]"
@gp :- S_gsb' "with image title 'GSB'"
Gnuplot.save("gsb.png", term="pngcairo size 550,350 fontscale 0.8")

# SE -- same loader, different signal type
traj_se = load_trajectory("se/se_ene.dat", "se/se_osc.dat")
S_se = spectrum(SE(), traj_se, laser, ω_pump, ω_probe)

@show S_se

# Time-resolved 2D map
@gp "set xlabel 'Time step'" "set ylabel 'Probe index'" :-
@gp :- "set xrange [0:50]"
@gp :- "set yrange [1:200]"
@gp :- S_se' "with image title 'SE'"
Gnuplot.save("se.png", term="pngcairo size 550,350 fontscale 0.8")

# different pump and probe lasers
pump_laser  = LaserPulse(τ=5.0, envelope=Gaussian())
probe_laser = LaserPulse(τ=10.0, envelope=Sech())
S_nd = spectrum(GSB(), traj_gs, pump_laser, ω_pump, probe_laser, ω_probe)

@gp "set xlabel 'Time step'" "set ylabel 'Probe index'" :-
@gp :- "set xrange [0:50]"
@gp :- "set yrange [1:200]"
@gp :- S_nd' "with image title 'different pump and probe lasers'"
Gnuplot.save("gsb_weird_laser.png", term="pngcairo size 550,350 fontscale 0.8")
