import numpy as np
from pathlib import Path
import os

def generate_orthorhombic(ions_on_sides):
    """
    Ions_on_side=4 will generate points with x coordinates 0, 0.25, 0.5 and 0.75.
    Note that 1 is kind of another cell already. So, side/ions is the step size
    There will be prod(ions_on_side) points in total in the cell.
    Check separately whether you get charge NEUTRAL system!
    """

    step = 1.0 / np.array(ions_on_sides)

    pos = np.zeros((ions_on_sides[0] * ions_on_sides[1] * ions_on_sides[2], 3))
    # print("The total number of points in the cell is ", len(self.ions))

    row = 0
    for (i, j, k) in np.ndindex(ions_on_sides[0], ions_on_sides[1], ions_on_sides[2]):
        pos[row,] = np.array([i * step[0], j * step[1], k * step[2]])
        row = row + 1
    return pos


def cubic(ions_on_side):
    """
    generate points for the cubic grid
    """
    return generate_orthorhombic([ions_on_side, ions_on_side, ions_on_side])

class Phase:
    filedir = Path(__file__).resolve().parents[2] / 'data/'

    def __init__(self, phase_name):

        self.types = []

        self.radius = {}

        self.closest_distance = {}  # dict['O', 'O'] = 2.2 A

        self.charge = {}

        self.buck = {}

        self.name = phase_name

        self.location = self.filedir / phase_name

        self.garnet = False  # HACK:garnet case is treated separately

        with open(os.path.join(".", self.filedir / phase_name / 'radii.lib'), 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n')
                line = line.split()
                self.types.append(line[0])
                self.radius[line[0]] = float(line[1])
        print('Radii', self.radius)

        try:
            with open(os.path.join(".", self.filedir / phase_name / 'dist.lib'), 'r') as f:
                for line in f.readlines():
                    if line.startswith('#'):
                        continue
                    line = line.rstrip('\n')
                    line = line.split()
                    pair = (min(line[0], line[1]), max(line[0], line[1]))
                    self.closest_distance[pair] = float(line[2])
            print('Overriding Shannon radii with the following closest distances:', self.closest_distance)
        except IOError:
            print("No closest distances between the ions were provided")
            print("Relying on Shannon radii only")

        try:

            with open(os.path.join(".", self.filedir / phase_name / 'buck.lib'), 'r') as f:
                charge_lines = False
                buck_lines = False

                for line in f.readlines():
                    if line.startswith('#'):
                        continue
                    line = line.rstrip('\n')
                    if len(line) > 0:

                        if 'species' in line:
                            charge_lines = True
                            buck_lines = False
                            continue

                        if 'buck' in line:
                            buck_lines = True
                            charge_lines = False
                            continue

                        if charge_lines:
                            line = line.split()
                            self.charge[line[0]] = float(line[-1])

                        if buck_lines:
                            line = line.split()
                            pair = (min(line[0], line[2]), max(line[0], line[2]))
                            self.buck[pair] = {}
                            self.buck[pair]['par'] = list(map(float, line[4:7]))
                            self.buck[pair]['lo'] = float(line[7])
                            self.buck[pair]['hi'] = float(line[-1])
                            # print(line)

                    # print(len(line))
                    # line = line.split()
            print('Charges:', self.charge)
        except IOError:
            print("There is no buck file! I assume that we are dealing with the garnet problem.")

            # Handcrafted parameters, write a parser later on
            self.gar_param = {}
            self.gar_param[('Al', 'O')] = {'De': 0.361581, 'a0': 1.900442, 'r0': 2.164818, 'A': 0.9, 'lo': 0.0,
                                           'hi': 15.0}
            self.gar_param[('Ca', 'O')] = {'De': 0.030211, 'a0': 2.2413340, 'r0': 2.923245, 'A': 5.0, 'lo': 0.0,
                                           'hi': 15.0}
            self.gar_param[('O', 'O')] = {'De': 0.042395, 'a0': 1.379316, 'r0': 3.618701, 'A': 22.0, 'lo': 0.0,
                                          'hi': 15.0}
            self.gar_param[('O', 'Si')] = {'De': 0.340554, 'a0': 2.0067, 'r0': 2.1, 'A': 1.0, 'lo': 0.0, 'hi': 15.0}

            self.charge['Al'] = 1.8
            self.charge['Ca'] = 1.2
            self.charge['O'] = -1.2
            self.charge['Si'] = 2.4

            self.garnet = True

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

