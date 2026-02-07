const Ha2eV  = 27.2113862459      # Hartree -> eV
const eV2Ha  = 1.0 / Ha2eV        # eV -> Hartree
const fs2au  = 41.341374575751    # fs -> a.u. time
const au2fs  = 1.0 / fs2au        # a.u. time -> fs

# Careful! Energy is reciprocal of wavelength; these are dangerously named
const eV2nm = 1239.84
const nm2eV = 1.0 / eV2nm

const eV2invcm = 8065.5439373492107
const invcm2eV = 1.0 / eV2invcm

# Factorial constants for Taylor expansions; a little bit OTT perhaps
# I wonder if this actually gets any speed increase? Only used in Gaussian wavepackets? Are
# those called repeatedly for the Window?
const _inv_fact2 = 1 / factorial(2)
const _inv_fact3 = 1 / factorial(3)
const _inv_fact4 = 1 / factorial(4)
const _inv_fact5 = 1 / factorial(5)

