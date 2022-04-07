# -*- coding: utf-8 -*-
"""
The 2D Ising Model:

In statistical mechanics, the Ising model is used to model ferromagnetism. The
model analyses the magnetic dipole moments of atomic spins within a lattice.
These spins are represented by +1 for spin-up and -1 for spin-down. Each spin
within the lattice interacts with its neighbours, and influence their spin
direction. Parallel spin allignment reduces the total energy of the lattice
and is therefore more energetically favourable as the system tends to its
lowest energy state. Heat interferes with this tendency, allowing for higher
energy configurations to occur.
This programme will simulate the 2D Ising model. First a 2D numpy array will be
generated, representing a lattice. At each point will be a value of +/- 1,
representing spin-up and spin-down, respectively. The energy of each lattice
point will be calculated by calculating the sum of the products of the spins
with each of the points neighbours:
                            E = −∑⟨i,j⟩ J.σi.σj
Where E is the energy of the lattice point, i and j representing neighbouring
points, J is the spin-spin interaction, and σi and σj are the individual spins
of each lattice point.
Once the energy has been calculated, a random point in the lattice will be
chosen and its spin flipped. The energy of this new lattice will then be
calculated. The energy difference between both lattices will be calculated
using:
                            Delta_E = E2 - E1
Where Delta_E is the energy difference, E2 is the energy of the lattice with a
randomly flipped spin, and E1 is the energy of the original lattice.
For a temperature independent simulation, if Delta_E is less than zero, the new
lattice configuration is more energy favourable and will replace the original.
In a temperature dependent system, a temperature threshold function will be
used to determine if a lattice configuration which is less energetically
favourable will still replace the original lattice configuration:
                            P = exp(-Delta_E/T)
Where P is the probability distribution, Delat_E is the energy difference, and
T is the temperature.
If a randomly generated number r is less than P, the new configuration will
replace the original lattice configuration, regardless of Delta_E.
This series of steps will be run over a specified number of iterations at
temperature T for a randomly generated lattice. Once the lattice has reached
equilibrium (approximately) the lattice will be graphed using a colour map. A
graph of the lattice energy verses temperature will also be plotted. This will
be repeated over increasing temperatures. The results will be plotted and
animated using matplotlib.
"""

# Import libraries/modules
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

#Variables
# Number of rows in lattice
n_rows = 10
# Number of columns in lattice
n_cols = 10
# Number of frames for animation
num_frames = 500
# Initialise list to store temperature values
T_vals = []
# Initialise list to store lattice energy values
Energy_vals = []
# Final temperature
Tf = 0.5

def create_lattice(n_rows, n_cols, flag=0):
    """
    Create a 2D array representing a 2D lattice
    
    ...
    
    Parameters
    ----------
    n_rows : TYPE - Integar
        DESCRIPTION - Number of rows in a 2D latice
    n_cols : TYPE - Integar
        DESCRIPTION - Number of columns in a 2D lattice
    flag : TYPE, optional - Choose lattice type
        DESCRIPTION. The default is 0.
                    0 -> Random +1 or -1 at each point
                    1 -> All spin up (+1)
                    2 -> All spin down (-1)
    Returns
    -------
    output_list : TYPE - 2D numpy array
        DESCRIPTION - Lattice with each value as +1 (spin up) or -1 (spin down)
    """
    # If flag is 1
    if flag == 1:
        # Create lattice with all points spin up (+1)
        return np.ones([n_rows, n_cols])
    # If flag is 2
    elif flag == 2:
        # Create lattice with all points spin down (-1)
        return -np.ones([n_rows, n_cols])
    # If flag is 0 (or not equal to 1 or 2)
    else:
        # Create lattice with each point randomly spin up (+1) or spin down (-1)
        return np.round(np.random.rand(n_rows, n_cols)) * 2 - 1
 
def calculate_random_coordinate(range_max, range_min = 0):
    """
    Calculate 2 random coordinates with a defined range
    
    ...
    
    Parameters
    ----------
    range_max : TYPE - Float
        DESCRIPTION - Final range value
    range_min : TYPE - Float, optional - Default is 0
        DESCRIPTION - Initial range value
    Returns
    -------
    TYPE - Float
        DESCRIPTION - A random coordinate between range_min and range_max
    """
    # Calculate and return random coordinate between specifified range
    return (range_max - range_min) * np.random.rand() - range_min

def flip_spin(lattice):
    '''
    Flip to spin of a random lattice point within the lattice
    
    ...
    
    Parameters
    ----------
    lattice : TYPE - 2D array
        DESCRIPTION - 2D array of lattice points
    Returns
    -------
    lattice : TYPE - 2D array
        DESCRIPTION - Lattice with the spin of a random lattice point flipped
    '''
    # Create new copy of lattice within function
    lattice = np.array(lattice)
    # Get the number of rows and columns in the lattice
    n_rows, n_cols = lattice.shape
    # Pick a random row from the lattice
    row = int(np.floor(calculate_random_coordinate(n_rows)))
    # Pick a random column from the lattice
    col = int(np.floor(calculate_random_coordinate(n_cols)))
    # Flip the spin of the lattice point at position rand_row, rand_col
    lattice[row, col] *= -1
    # Return updated lattice
    return lattice

def get_lattice_energy(lattice, J = -1):
    '''
    Calculate the total energy of the lattice
    
    ...
    
    Parameters
    ----------
    lattice : TYPE - 2D array
        DESCRIPTION - 2D array of lattice points
    J : TYPE - Integer
        DESCRIPTION - Default is -1
                    -1 -> -1 for ferromagnetic coupling
                    +1 -> +1 for antiferromagnetic coupling
    Returns
    -------
    energy : TYPE - Float
        DESCRIPTION - Eenergy of lattice point
    '''
    # Calculate the total energy of the lattice
    lattice_energy = sum(sum(                       # Sum the energy of all lattice points
                     J * (                          # Multiply each lattice shift by J
                          np.roll(lattice, 1, 0)    # Shift lattice down
                          + np.roll(lattice, -1, 0) # Shift lattice up
                          + np.roll(lattice, 1, 1)  # Shift lattice left
                          + np.roll(lattice, -1, 1) # Shift lattice right
                          ) * lattice               # Multiply each lattice shift by the original lattice
        ))
    # Return the energy per lattice point
    return lattice_energy/(2*lattice.size)

def test_spin_flip(lattice, T):
    '''
    Test if flipping the spin of a random lattice point decreases the total energy
    of the lattice and return the lowest energy lattice.
    
    ...
    
    Parameters
    ----------
    lattice : TYPE - 2D array
        DESCRIPTION - 2D array representing a 2D lattice of spins +/- 1
    T : TYPE - Float
        DESCRIPTION - Temperature
    Returns
    -------
    lattice : TYPE - 2D array
        DESCRIPTION - Return a new lattice if flipping the spin of the a random
                      point decreases the lattice energy, or if the random number
                      is less than the probability distribution, else return the
                      original lattice
    '''
    # Create new copy of lattice within function
    lattice = np.array(lattice)
    # Calcuate energy of the original lattice
    E1 = get_lattice_energy(lattice)
    # Flip the spin a random lattice point and save new lattice
    lattice_spin_flipped = flip_spin(lattice)
    # Calculate energy of the of the new lattice
    E2 = get_lattice_energy(lattice_spin_flipped)
    # Calculate the energy difference between the 2 lattices
    delta_E = E2 - E1
    # If the energy difference is less than zero
    if delta_E < 0:
        # Return the new lattice
        return lattice_spin_flipped
    # If the randomly generated number is less than the probability distribution
    elif np.random.rand() < np.exp(-delta_E/T):
        # Return the new lattice
        return lattice_spin_flipped
    # Else return the original lattice
    else:
        return lattice

def equilibrate(n_rows, n_cols, T):
    '''
    Run the test_spin_flip function for a calculated number of iterations to cause
    the lattice to reach a point of equilibrium.
    
    ...
    
    Parameters
    ----------
    n_rows : TYPE - Integer
        DESCRIPTION - Number of rows in lattice
    n_cols : TYPE - Integer
        DESCRIPTION - Number of columns in lattice
    T : TYPE - Float
        DESCRIPTION - Temperature
    Returns
    -------
    lattice : TYPE - 2D numpy array
        DESCRIPTION - Lattice with each value as +1 (spin up) or -1 (spin down)
    '''
    # Create lattice
    lattice = create_lattice(n_rows, n_cols, 0)
    # Determine number of iterations
    iterations = range(lattice.size * 15)
    # Loop through number of iterations
    for i in iterations:
        # Test if flipping spin at a random lattice point increases/decreases energy
        lattice = test_spin_flip(lattice, T)
    # Return equilibrated lattice
    return lattice

def animate(frame_number, frames, Tf, n_rows, n_cols):
    '''
    Animate the progression of the lattice as it moves towards equilibrium
    
    ...
    
    Parameters
    ----------
    frame_number : TYPE - Integer
        DESCRIPTION - Current frame
    frames : TYPE - Integer
        DESCRIPTION - Total number of frames
    Tf : TYPE - Float
        DESCRIPTION - Final temperature
    n_rows : TYPE - Integer
        DESCRIPTION - Number of rows in lattice
    n_cols : TYPE - Integer
        DESCRIPTION - Number of columns in lattice
    Returns
    -------
    None.
    '''
    # Calculate temperature
    T = frame_number/frames * Tf + 0.000001
    # Test if flipping the spin of random latice point increases/decreases lattice energy
    lattice = equilibrate(n_rows, n_cols, T)
    # Append current temperature to T_vals list
    T_vals.append(T)
    # Append the current lattice energy to Energy_vals list
    Energy_vals.append(get_lattice_energy(lattice))
    # Update data for regular plot (ax2)
    pl.set_data(T_vals, Energy_vals)
    # Map colours based on lattice 2D array
    im.set_array(lattice)


''' Plot data '''
# Create figure
fig = plt.figure()
# Figure title
plt.suptitle('The Ising Model')

# Create subplot 1
ax1 = fig.add_subplot(1, 2, 1)
# ax1 title
ax1.set_title('Lattice')
# Set ax1 x and y labels
ax1.set(xlabel='x', ylabel='y')
# Add gridlines to ax1
ax1.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2, animated=True)
# Place ticks at -0.5 to n_cols with increments of 1
ax1.set_xticks(np.arange(-0.5, n_cols, 1))
# Place ticks at -0.5 to n_rows with increments of 1
ax1.set_yticks(np.arange(-0.5, n_rows, 1))
# Remove tick labels from ax1 x and y axis
ax1.set_xticklabels([])
ax1.set_yticklabels([])

# Create subplot 2
ax2 = fig.add_subplot(1, 2, 2, xlim=(0, Tf), ylim=(-2, 1))
# Set ax2 title
ax2.set_title('Lattice Energy vs. Temperature')
# Set ax2 x and y labels
ax2.set(xlabel='Temperature', ylabel='Energy per Lattice Point')
# Add gridlines to ax2
ax2.grid(which='major', axis='both', linestyle='-', color='k', linewidth=0.25, animated=True)

# Create initial lattice
lattice = create_lattice(n_rows, n_cols, 0)
# Plot lattice as a colour chart
im = ax1.imshow(lattice)
# Get reference to the regular plot (ax2)
pl, = ax2.plot(T_vals, Energy_vals)

# Animate figure
anim = animation.FuncAnimation(fig,                 # Figure to animate
                               animate,             # Function for animation
                               frames = num_frames, # Number of frames
                               interval = 10,       # Time interval between frames(ms)
                               repeat = False,      # Do not repeat animation
                               fargs = (num_frames, Tf, n_rows, n_cols)) # Function arguments
