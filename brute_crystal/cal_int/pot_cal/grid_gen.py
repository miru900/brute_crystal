import numpy as np

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
