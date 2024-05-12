from ase.io import read, write
import numpy as np
import getopt, sys

print(sys.argv)

# get full passed arguments
argument_list = sys.argv[1:]

# define short and long list of options
smll_options = 'ht:n:'
long_options = ['help','traj=','number=']

try:
    arguments, options = getopt.getopt(argument_list, smll_options, long_options)
    
    for currArgument, currValue in arguments:
        print(currArgument,currValue    )
        if currArgument in ("-h","--help"):
            print('HELP!')
        elif currArgument in ("-t","--traj"):
            print('found base: ',currValue)
            base_flake = read(currValue)
            mol_name = currValue.split('.')[0]
        elif currArgument in ("-n","--number"):
            print('Produced flake radius: ', currValue)
            N = int(currValue)
except getopt.error as err:
    print(str(err))

# copy the flake to create the output base
output_flake = base_flake.copy()

# get unit cell parameters
base_cell = base_flake.get_cell()
base_cell_lengths = base_cell.lengths()
base_cell_anglesD = base_cell.angles()
base_cell_anglesR = np.deg2rad(base_cell_anglesD)

# assumption 1: the unit cell vector [0] is aligned with the cartesian x-coordinates
# assumption 2: the unit cell vectors [0,1] are aligned with the cartesian xy-plane
x_step = base_cell_lengths[0]
y_step = base_cell_lengths[1]*np.sin(np.deg2rad(60))

# assumption 3: I'm dealing with a hexagonal unit cell, therefore I can simply compute the repetition
# points based on equilateral triangle stacking. The objective is to form a hexagonal flake
c=0
for m in range(0,N):
#    print(m,c) # <--- used for debugging
    # because we're stacking triangles to form hexagons, each level up require one less horizontal repetition
    x_range = 2*N-m-1
    #I'll run the loop using integer repetitions of the step, thus we start on 1 and go up to range
    for n in range(0,x_range):
        # when m = 0, we're at the first line, no going up!
        if m != 0:
            x = n*x_step+c
            y = m*y_step
        else:
            x = n*x_step
            y = 0
#        print(n,[x,y,0]) # <--- used for debugging
        
        # this is a very ugly way to do the movements
        base_flake.translate([x,y,0])    # move up
        output_flake += base_flake       # add
        base_flake.translate([0,-2*y,0]) # move down
        output_flake += base_flake       # add
        base_flake.translate([-x,+y,0])  # move to origin

    # update the starting point correction
    c += x_step/2


# making it pretty!
output_flake.center()
output_flake.set_pbc([False, False, False])

# save file
write(f'flake_{mol_name}_{N}.traj',output_flake)
