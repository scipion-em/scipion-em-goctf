# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import numpy as np
from collections import OrderedDict

from pyworkflow.object import ObjectWrap
import pyworkflow.utils as pwutils
from pwem.constants import ALIGN_2D, ALIGN_3D, ALIGN_PROJ, ALIGN_NONE
from pwem.convert.transformations import translation_from_matrix
import pwem.emlib.metadata as md


CTF_DICT = OrderedDict([
       ("_defocusU", md.RLN_CTF_DEFOCUSU),
       ("_defocusV", md.RLN_CTF_DEFOCUSV),
       ("_defocusAngle", md.RLN_CTF_DEFOCUS_ANGLE)
       ])


class CoordinatesWriter:
    """ Simple class to write coordinates and CTF into a star file. """
    HEADER = """
data_

loop_
_rlnCoordinateX #1
_rlnCoordinateY #2
_rlnDefocusU #3
_rlnDefocusV #4
_rlnDefocusAngle #5
"""

    def __init__(self, filename):
        """ Filename where to write the coordinates. """
        pwutils.makePath(os.path.dirname(filename))  # Ensure path exists
        self._f = open(filename, 'w')
        self._f.write(self.HEADER)

    def writeRow(self, x, y, defU, defV, defAng):
        self._f.write(f"{x:.2f} {y:.2f} {defU:.2f} {defV:.2f} {defAng:.2f}\n")

    def close(self):
        self._f.close()


def rowToCtfModel(ctfRow, ctfModel):
    """ Create a CTFModel from a row of a meta """
    if ctfRow.containsAll(CTF_DICT):
        for attr, label in CTF_DICT.items():
            value = ctfRow.getValue(label)
            if not hasattr(ctfModel, attr):
                setattr(ctfModel, attr, ObjectWrap(value))
            else:
                getattr(ctfModel, attr).set(value)

        ctfModel.standardize()
    else:
        ctfModel = None

    return ctfModel


def getShifts(transform, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
    if alignType == ALIGN_NONE:
        return None

    inverseTransform = alignType == ALIGN_PROJ
    # only flip is meaningful if 2D case
    # in that case the 2x2 determinant is negative
    flip = False
    matrix = transform.getMatrix()
    if alignType == ALIGN_2D:
        # get 2x2 matrix and check if negative
        flip = bool(np.linalg.det(matrix[0:2, 0:2]) < 0)
        if flip:
            matrix[0, :2] *= -1.  # invert only the first two columns keep x
            matrix[2, 2] = 1.  # set 3D rot

    elif alignType == ALIGN_3D:
        flip = bool(np.linalg.det(matrix[0:3, 0:3]) < 0)
        if flip:
            matrix[0, :4] *= -1.  # now, invert first line including x
            matrix[3, 3] = 1.  # set 3D rot

    shifts = geometryFromMatrix(matrix, inverseTransform)

    return shifts


def geometryFromMatrix(matrix, inverseTransform):
    if inverseTransform:
        matrix = np.linalg.inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    return shifts
