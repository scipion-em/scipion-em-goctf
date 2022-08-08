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
from collections import OrderedDict
from enum import Enum

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA, SCIPION_DEBUG_NOCLEAN
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem import emlib
import pwem.emlib.metadata as md
from pwem.objects import SetOfParticles
from pwem.protocols import EMProtocol, ProtParticles

from .. import Plugin
from ..convert import CoordinatesWriter, rowToCtfModel, getShifts


class outputs(Enum):
    outputParticles = SetOfParticles


class ProtGoCTF(ProtParticles):
    """ Geometrically-Optimized CTF determination for single particle cryo-EM

    To find more information about goCTF visit:
    https://www.lsi.umich.edu/science/centers-technologies/cryo-electron-microscopy/research/goctf
    """
    _label = 'ctf refinement'
    _devStatus = BETA
    _possibleOutputs = outputs

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      important=True,
                      pointerCondition='hasCTF',
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles for local CTF refinement.')
        form.addParam('applyShifts', params.BooleanParam, default=False,
                      label='Apply particle shifts?',
                      help='Apply particle shifts from 2D alignment to '
                           'recalculate new coordinates. This can be useful '
                           'for re-centering particle coordinates.')
        form.addParam('inputMicrographs', params.PointerParam, important=True,
                      label='Input micrographs',
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrographs related to input particles.')

        form.addSection(label='Params')
        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      help='Set to 1 for no downsampling. Non-integer '
                           'downsample factors are possible. This downsampling '
                           'is only used for estimating the CTF and it does not '
                           'affect any further calculation. Ideally the estimation '
                           'of the CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) and '
                           'not occupying the whole power spectrum (since this '
                           'downsampling might entail aliasing).')

        form.addParam('windowSize', params.IntParam, default=512,
                      label='FFT box size (px)',
                      help='The dimensions (in pixels) of the amplitude '
                           'spectrum goCTF will compute. Smaller box '
                           'sizes make the fitting process significantly '
                           'faster, but sometimes at the expense of '
                           'fitting accuracy. If you see warnings '
                           'regarding CTF aliasing, consider '
                           'increasing this parameter.')

        group = form.addGroup('Search limits')
        line = group.addLine('Resolution (A)',
                             help='The CTF model will be fit to regions '
                                  'of the amplitude spectrum corresponding '
                                  'to this range of resolution.')
        line.addParam('lowRes', params.FloatParam, default=30., label='Min')
        line.addParam('highRes', params.FloatParam, default=5., label='Max')

        line = group.addLine('Defocus search range (A)',
                             help='Select _minimum_ and _maximum_ values for '
                                  'defocus search range (in A). Underfocus'
                                  ' is represented by a positive number.')
        line.addParam('minDefocus', params.FloatParam, default=5000.,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=50000.,
                      label='Max')
        group.addParam('stepDefocus', params.FloatParam, default=500.,
                       label='Defocus step (A)',
                       help='Step size for the defocus search.')

        group.addParam('doRefine', params.BooleanParam, default=True,
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Run per-particle refinement?")

        form.addParallelSection(threads=2, mpi=1)

    # -------------------------- STEPS functions -------------------------------
    def _createMicDict(self):
        """ Create a dictionary with all micrographs that
         are both in the input micrographs set and there
         are particles belonging to it.
         micName will be the key to that dict.
         """
        inputParticles = self.inputParticles.get()
        firstCoord = inputParticles.getFirstItem().getCoordinate()
        self.hasMicName = firstCoord.getMicName() is not None
        inputMicDict = {mic.getMicName(): mic.clone()
                        for mic in self._getMicrographs()}
        # Check now which if these mics have particles belonging
        self.micDict = OrderedDict()
        # match the mic from coord with micDict
        lastMicId = None
        # TODO: If this loop is too expensive for very large input datasets,
        # we could consider using the aggregate functions in the mapper
        for particle in inputParticles.iterItems(orderBy='_micId'):
            micId = particle.getMicId()
            if micId != lastMicId:  # Do no repeat check when this is the same mic
                micName = particle.getCoordinate().getMicName()
                if micName in inputMicDict:
                    self.micDict[micName] = inputMicDict[micName]
                lastMicId = micId

    def _insertAllSteps(self):
        self._createMicDict()
        self._defineArgs()

        convIdDeps = [self._insertFunctionStep('convertInputStep')]
        refineDeps = []

        for micName, mic in self.micDict.items():
            stepId = self._insertFunctionStep('refineCtfStep', mic.getFileName(),
                                              prerequisites=convIdDeps)
            refineDeps.append(stepId)

        self._insertFunctionStep('createOutputStep', prerequisites=refineDeps)

    def _iterParticlesMic(self, newMicCallback):
        """ Iterate through particles sorting by micId and only for
        those that are present in the input set of micrographs. """
        inputParts = self.inputParticles.get()
        lastMicId = None

        for particle in inputParts.iterItems(orderBy=['_micId', 'id']):
            coord = particle.getCoordinate()
            micId = particle.getMicId()
            micName = coord.getMicName()

            if micId != lastMicId:  # Do no repeat check when this is the same mic
                mic = self.micDict.get(micName, None)
                if mic is None:
                    self.warning(f"Skipping all particles from micrograph, "
                                 f"key {micName} not found")
                else:
                    newMicCallback(mic)  # Notify about a new micrograph found
                lastMicId = micId

            if mic is not None:
                yield particle

    def convertInputStep(self):
        inputParts = self.inputParticles.get()
        alignType = inputParts.getAlignment()
        inputMics = self._getMicrographs()

        scale = inputParts.getSamplingRate() / inputMics.getSamplingRate() / self.ctfDownFactor.get()
        doScale = abs(scale - 1.0 > 0.00001)
        if doScale:
            self.info(f"Scaling coordinates by a factor {scale:0.2f}")

        self._lastWriter = None
        coordDir = self._getTmpPath()

        def _newMic(mic):
            if self._lastWriter:
                self._lastWriter.close()
            micBase = pwutils.removeBaseExt(mic.getFileName())
            posFn = os.path.join(coordDir, micBase, micBase + '_go.star')
            self._lastWriter = CoordinatesWriter(posFn)

        for particle in self._iterParticlesMic(newMicCallback=_newMic):
            coord = particle.getCoordinate()
            x, y = coord.getPosition()
            if self.applyShifts:
                shifts = getShifts(particle.getTransform(), alignType)
                x, y = x - int(shifts[0]), y - int(shifts[1])
            if doScale:
                x, y = x * scale, y * scale
            ctf = particle.getCTF()
            self._lastWriter.writeRow(x, y, ctf.getDefocusU(),
                                      ctf.getDefocusV(), ctf.getDefocusAngle())

        if self._lastWriter:
            self._lastWriter.close()  # Close file writing for last mic

    def refineCtfStep(self, micFn):
        if not os.path.exists(micFn):
            raise FileNotFoundError("Missing input micrograph: %s" % micFn)

        downFactor = self.ctfDownFactor.get()
        # We convert the input micrograph on demand if not in .mrc
        micPath = self._getTmpPath(pwutils.removeBaseExt(micFn))
        ih = emlib.image.ImageHandler()
        micFnMrc = os.path.join(micPath, pwutils.replaceBaseExt(micFn, 'mrc'))

        if downFactor != 1:
            ih.scaleFourier(micFn, micFnMrc, downFactor)
        else:
            ih.convert(micFn, micFnMrc, emlib.DT_FLOAT)

        # Run goCTF
        try:
            self._params.update({
                'micFn': os.path.basename(micFnMrc),
                'goctfOut': self._getOutputPath(micFn, ext="_ctf.log"),
                'goctfPSD': self._getOutputPath(micFn, ext="_ctf.mrc")
            })
            args = self._args % self._params
            self.runJob(Plugin.getProgram(), args,
                        env=Plugin.getEnviron(),
                        cwd=micPath)

            # Let's clean the temporary mrc micrograph
            if not pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN):
                pwutils.cleanPath(micFnMrc)

        except:
            self.error(f"ERROR: goCTF has failed on {micFnMrc}")
            import traceback
            traceback.print_exc()

    def createOutputStep(self):
        inputParts = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        partSet.copyInfo(inputParts)
        self._rowList = None
        self._rowCounter = 0

        def _newMic(mic):
            micFn = mic.getFileName()
            micPath = self._getTmpPath(pwutils.removeBaseExt(micFn))
            ctfFn = os.path.join(micPath,
                                 self._getOutputPath(micFn, ext="_goCTF.star"))
            self._rowCounter = 0
            if os.path.exists(ctfFn):
                self._rowList = [row.clone() for row in md.iterRows(ctfFn)]
            else:
                self._rowList = None

        for particle in self._iterParticlesMic(newMicCallback=_newMic):
            if self._rowList is None:  # Ignore particles if no CTF
                continue
            newPart = particle.clone()
            row = self._rowList[self._rowCounter]
            self._rowCounter += 1
            rowToCtfModel(row, newPart.getCTF())
            partSet.append(newPart)

        self._defineOutputs(**{outputs.outputParticles.name: partSet})
        self._defineTransformRelation(self.inputParticles, partSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        summary = []

        if not hasattr(self, 'outputParticles'):
            summary.append("Output is not ready yet.")
        else:
            summary.append("CTF refinement of %d particles."
                           % self.inputParticles.get().getSize())

        return summary

    def _methods(self):
        if self.inputParticles.get() is None:
            return ['Input particles not available yet.']
        methods = "We refined the CTF of "
        methods += self.getObjectTag('inputParticles')
        methods += " using goCTF [Su2019]. "

        if self.hasAttribute('outputParticles'):
            methods += 'Output particles: %s' % self.getObjectTag('outputParticles')

        return [methods]

    # -------------------------- UTILS functions -------------------------------
    def _defineArgs(self):
        self.inputMics = self._getMicrographs()
        acq = self.inputMics.getAcquisition()

        self._params = {'voltage': acq.getVoltage(),
                        'sphericalAberration': acq.getSphericalAberration(),
                        'ampContrast': acq.getAmplitudeContrast(),
                        'samplingRate': self.inputMics.getSamplingRate() * self.ctfDownFactor.get(),
                        'windowSize': self.windowSize.get(),
                        'lowRes': self.lowRes.get(),
                        'highRes': self.highRes.get(),
                        'minDefocus': self.minDefocus.get(),
                        'maxDefocus': self.maxDefocus.get(),
                        'step_focus': self.stepDefocus.get(),
                        'doRefine': "yes" if self.doRefine else "no"
                        }

        self._args = """   << eof > %(goctfOut)s
%(micFn)s
%(goctfPSD)s
%(samplingRate)f
%(voltage)f
%(sphericalAberration)f
%(ampContrast)f
%(windowSize)d
%(lowRes)f
%(highRes)f
%(minDefocus)f
%(maxDefocus)f
%(step_focus)f
%(doRefine)s
eof\n
"""

    def _getOutputPath(self, micFn, ext):
        return pwutils.removeBaseExt(micFn) + ext

    def _getMicrographs(self):
        return self.inputMicrographs.get()
