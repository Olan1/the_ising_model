# the_ising_model

## The 2D Ising Model:

<p>In statistical mechanics, the Ising model is used to model ferromagnetism. The
model analyses the magnetic dipole moments of atomic spins within a lattice.
These spins are represented by +1 for spin-up and -1 for spin-down. Each spin
within the lattice interacts with its neighbours, and influence their spin
direction. Parallel spin allignment reduces the total energy of the lattice
and is therefore more energetically favourable as the system tends to its
lowest energy state. Heat interferes with this tendency, allowing for higher
energy configurations to occur.</p>
<p>This programme will simulate the 2D Ising model. First a 2D numpy array will be
generated, representing a lattice. At each point will be a value of +/- 1,
representing spin-up and spin-down, respectively. The energy of each lattice
point will be calculated by calculating the sum of the products of the spins
with each of the points neighbours:</p>
                            E = −∑⟨i,j⟩ J.σi.σj
<p>Where E is the energy of the lattice point, i and j representing neighbouring
points, J is the spin-spin interaction, and σi and σj are the individual spins
of each lattice point.</p>
<p>Once the energy has been calculated, a random point in the lattice will be
chosen and its spin flipped. The energy of this new lattice will then be
calculated. The energy difference between both lattices will be calculated
using:</p>
                            Delta_E = E2 - E1
<p>Where Delta_E is the energy difference, E2 is the energy of the lattice with a
randomly flipped spin, and E1 is the energy of the original lattice.</p>
<p>For a temperature independent simulation, if Delta_E is less than zero, the new
lattice configuration is more energy favourable and will replace the original.
In a temperature dependent system, a temperature threshold function will be
used to determine if a lattice configuration which is less energetically
favourable will still replace the original lattice configuration:</p>
                            P = exp(-Delta_E/T)
</p>Where P is the probability distribution, Delat_E is the energy difference, and
T is the temperature.</p>
<p>If a randomly generated number r is less than P, the new configuration will
replace the original lattice configuration, regardless of Delta_E.
This series of steps will be run over a specified number of iterations at
temperature T for a randomly generated lattice. Once the lattice has reached
equilibrium (approximately) the lattice will be graphed using a colour map. A
graph of the lattice energy verses temperature will also be plotted. This will
be repeated over increasing temperatures. The results will be plotted and
animated using matplotlib.</p>
