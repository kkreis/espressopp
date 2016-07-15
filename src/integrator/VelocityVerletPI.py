#  Copyright (C) 2012,2013,2014,2015
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


"""
******************************************
**espressopp.integrator.VelocityVerletPI**
******************************************

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.MDIntegrator import *
from _espressopp import integrator_VelocityVerletPI

class VelocityVerletPILocal(MDIntegratorLocal, integrator_VelocityVerletPI):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system, _verletlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerletPI, system, _verletlist)
    
    def setTimeStep(self, real):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setTimeStep(self, real)        
            
    def setmStep(self, int):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setmStep(self, int)

    def getmStep(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getmStep(self)
        
    def setsStep(self, int):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setsStep(self, int)

    def getsStep(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getsStep(self)
        
    def setNtrotter(self, int):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setNtrotter(self, int)

    def getNtrotter(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNtrotter(self)
        
    def setTemperature(self, real):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setTemperature(self, real)

    def getTemperature(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getTemperature(self)
        
    def setGamma(self, real):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setGamma(self, real)

    def getGamma(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getGamma(self)

    def setCMDparameter(self, real):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setCMDparameter(self, real)

    def getCMDparameter(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCMDparameter(self)

    def setPILElambda(self, real):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPILElambda(self, real)

    def getPILElambda(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPILElambda(self)
        
    def setClmassmultiplier(self, real):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setClmassmultiplier(self, real)

    def getClmassmultiplier(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getClmassmultiplier(self)
        
    def setSpeedup(self, bool):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setSpeedup(self, bool)

    def getSpeedup(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getSpeedup(self)
        
    def setKTI(self, bool):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setKTI(self, bool)
            
    def getKTI(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getKTI(self)

    def setCentroidThermostat(self, bool):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setCentroidThermostat(self, bool)
            
    def getCentroidThermostat(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCentroidThermostat(self)

    def setPILE(self, bool):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPILE(self, bool)
            
    def getPILE(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPILE(self)

    def setRealKinMass(self, bool):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setRealKinMass(self, bool)
            
    def getRealKinMass(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getRealKinMass(self)

    def setConstKinMass(self, bool):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setConstKinMass(self, bool)

    def getConstKinMass(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getConstKinMass(self)
        
    def setVerletList(self, _verletlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setVerletList(self, _verletlist)
            
    def getVerletList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)
    
    def addEigenvalues(self, evalueslist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            for eval in evalueslist: 
                self.cxxclass.addValues(self, eval)
                
    def computeRingEnergy(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeRingEnergy(self)
        
    def computeRingEnergyRaw(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeRingEnergyRaw(self)
        
    def computeKineticEnergy(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeKineticEnergy(self)
        
    def computePositionDrift(self, int):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computePositionDrift(self, int)
        
    def computeMomentumDrift(self, int):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeMomentumDrift(self, int)    
                
    def addEigenvectors(self, evlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            for ev in evlist: 
                for val in ev:
                    self.cxxclass.add(self, val)
                self.cxxclass.addEV(self)
            self.cxxclass.transp(self)
            
if pmi.isController :
    class VelocityVerletPI(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.integrator.VelocityVerletPILocal',
          pmiproperty = [ 'sStep', 'mStep', 'ntrotter', 'gamma', 'CMDparameter', 'PILElambda', 'temperature', 'clmassmultiplier', 'speedup', 'KTI', 'constkinmass', 'verletList', 'centroidthermostat', 'PILE', 'realkinmass' ],
          pmicall = ['resetTimers', 'setTimeStep', 'setmStep', 'setsStep', 'setNtrotter', 'setTemperature', 'setGamma', 'setCMDparameter', 'setPILElambda', 'setClmassmultiplier', 'setSpeedup', 'setKTI', 'setPILE', 'setRealKinMass', 'setCentroidThermostat', 'setConstKinMass', 'addEigenvectors', 'addEigenvalues', 'setVerletList', 'computeKineticEnergy', 'computeRingEnergy', 'computeRingEnergyRaw', 'computeMomentumDrift', 'computePositionDrift'],
          pmiinvoke = ['getTimeStep', 'getTimers', 'getmStep', 'getsStep', 'getNtrotter', 'getTemperature', 'getGamma', 'getCMDparameter', 'getPILElambda', 'getClmassmultiplier', 'getSpeedup', 'getKTI', 'getPILE', 'getRealKinMass', 'getCentroidThermostat', 'getConstKinMass', 'getVerletList' ]
        )
