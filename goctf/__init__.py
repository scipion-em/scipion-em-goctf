# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

import pwem

from .constants import *


__version__ = '3.0.0b1'
_references = ['Su2019']


class Plugin(pwem.Plugin):
    _homeVar = GOCTF_HOME
    _pathVars = [GOCTF_HOME]
    _supportedVersions = [V1_2_0]
    _url = "https://sites.google.com/umich.edu/min-su-ph-d/home"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(GOCTF_HOME, 'goctf-%s' % V1_2_0)

    @classmethod
    def getProgram(cls):
        """ Return the program binary that will be used. """
        return cls.getHome('bin', 'goctf')

    @classmethod
    def defineBinaries(cls, env):
        for v in cls._supportedVersions:
            env.addPackage('goctf', version=v,
                           tar='goctf_v%s.tgz' % v,
                           default=v == V1_2_0)
