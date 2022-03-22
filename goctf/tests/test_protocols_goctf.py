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
    def setData(cls, dataProject='goctf'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.partsFn = cls.dataset.getFile('particles')

    @classmethod
    def runImportMicrographs(cls):
        cls.protImport = ProtImportMicrographs(
            objLabel='import mics',
            filesPath=cls.micsFn,
            samplingRate=1.31)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        cls.assertIsNotNone(cls.protImport.outputMicrographs,
                            "SetOfMicrographs has not been produced.")
        return cls.protImport

    @classmethod
    def runImportParticles(cls):
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         objLabel='import parts',
                                         sqliteFile=cls.partsFn,
                                         samplingRate=1.31,
                                         importFrom=ProtImportParticles.IMPORT_FROM_SCIPION)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        cls.assertIsNotNone(cls.protImport.outputParticles,
                            "SetOfParticles has not been produced.")
        return cls.protImport


class TestGoCTF(TestGoCTFBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestGoCTF.setData()
        print(magentaStr("\n==> Importing data - particles:"))
        cls.protImportParts = cls.runImportParticles()
        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImportMics = cls.runImportMicrographs()

    def testRunGoctf1(self):
        print(magentaStr("\n==> Testing goctf:"))
        protCTF = ProtGoCTF()
        protCTF.inputParticles.set(self.protImportParts.outputParticles)
        protCTF.inputMicrographs.set(self.protImportMics.outputMicrographs)

        self.proj.launchProtocol(protCTF, wait=True)
        self.assertIsNotNone(protCTF.outputParticles,
                             "SetOfParticles has not been produced.")
        self.assertEqual(protCTF.inputParticles.get().getSize(),
                         protCTF.outputParticles.getSize())
