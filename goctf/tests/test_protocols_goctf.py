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

from pwem.protocols import ProtImportMicrographs, ProtImportParticles
from pyworkflow.utils import magentaStr
from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from goctf.protocols import ProtGoCTF


class TestGoCTFBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micFn = cls.dataset.getFile('allMics')
        cls.partFn1 = cls.dataset.getFile('particles2')
        cls.partFn2 = cls.dataset.getFile('particles3')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, sphericalAberration):
        cls.protImport = ProtImportMicrographs(
            objLabel='import mics', filesPath=pattern,
            samplingRate=samplingRate,
            voltage=voltage, sphericalAberration=sphericalAberration)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        cls.assertIsNotNone(cls.protImport.outputMicrographs,
                            "SetOfMicrographs has not been produced.")
        return cls.protImport

    @classmethod
    def runImportParticles(cls, pattern, label, samplingRate, voltage):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         objLabel=label,
                                         sqliteFile=pattern,
                                         samplingRate=samplingRate,
                                         voltage=voltage,
                                         importFrom=ProtImportParticles.IMPORT_FROM_SCIPION)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        cls.assertIsNotNone(cls.protImport.outputParticles,
                            "SetOfParticles has not been produced.")
        return cls.protImport

    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern,
                                       samplingRate=1.237,
                                       voltage=300,
                                       sphericalAberration=2)

    @classmethod
    def runImportParticlesBPV(cls, pattern, label):
        """ Run an Import particles protocol. """
        return cls.runImportParticles(pattern,
                                      label=label,
                                      samplingRate=4.95,
                                      voltage=300)


class TestGoCTF(TestGoCTFBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestGoCTF.setData()
        print(magentaStr("\n==> Importing data - particles (no alignment):"))
        cls.protImport1 = cls.runImportParticlesBPV(cls.partFn1,
                                                    label='import particles (no alignment)')
        print(magentaStr("\n==> Importing data - particles (with alignment):"))
        cls.protImport2 = cls.runImportParticlesBPV(cls.partFn2,
                                                    label='import particles (with alignment)')
        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImportMics = cls.runImportMicrographBPV(cls.micFn)

    def testRunGoctf1(self):
        print(magentaStr("\n==> Testing goctf (1/2):"))
        protCTF = ProtGoCTF(objLabel='goctf (1/2)')
        protCTF.inputParticles.set(self.protImport1.outputParticles)
        protCTF.inputMicrographs.set(self.protImportMics.outputMicrographs)

        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputParticles,
                             "SetOfParticles has not been produced.")
        self.assertEqual(protCTF.inputParticles.get().getSize(),
                         protCTF.outputParticles.getSize())

    def testRunGoctf2(self):
        print(magentaStr("\n==> Testing goctf (2/2):"))
        protCTF = ProtGoCTF(objLabel='goctf (2/2)')
        protCTF.inputParticles.set(self.protImport2.outputParticles)
        protCTF.applyShifts.set(True)
        protCTF.inputMicrographs.set(self.protImportMics.outputMicrographs)

        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputParticles,
                             "SetOfParticles has not been produced.")
        self.assertEqual(protCTF.inputParticles.get().getSize(),
                         protCTF.outputParticles.getSize())
