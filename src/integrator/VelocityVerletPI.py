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
        
    def setTemperature(self, int):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setTemperature(self, int)

    def getTemperature(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getTemperature(self)
        
    def setGamma(self, int):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setGamma(self, int)

    def getGamma(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getGamma(self)
        
    def setSpeedup(self, bool):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setSpeedup(self, bool)

    def getSpeedup(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getSpeedup(self)
        
    def setVerletList(self, _verletlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setVerletList(self, _verletlist)
            
    def getVerletList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)
        
    def addEigenvectors(self, evlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            for ev in evlist: 
                for val in ev:
                    self.cxxclass.add(self, val)
                self.cxxclass.addEV(self);
            self.cxxclass.transp(self);
            
if pmi.isController :
    class VelocityVerletPI(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.integrator.VelocityVerletPILocal',
          pmiproperty = [ 'sStep', 'mStep', 'ntrotter', 'gamma', 'temperature', 'speedup', 'verletList' ],
          pmicall = ['resetTimers', 'setTimeStep', 'setmStep', 'setsStep', 'setNtrotter', 'setTemperature', 'setGamma', 'setSpeedup', 'addEigenvectors', 'setVerletList'],
          pmiinvoke = ['getTimeStep', 'getTimers', 'getmStep', 'getsStep', 'getNtrotter', 'getTemperature', 'getGamma', 'getSpeedup', 'getVerletList' ]
        )
