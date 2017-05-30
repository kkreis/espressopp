/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "python.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "ParticleDistribution.hpp"
#include "bc/BC.hpp"

#include <boost/serialization/map.hpp>

using namespace espressopp;
using namespace iterator;
using namespace std;

namespace espressopp {
  namespace analysis {

    python::list ParticleDistribution::computeArray(int bins) const {

      System& system = getSystemRef();
      Real3D Li = system.bc->getBoxL();

      real dr;
      if (geometry == 1) {dr = Li[0] / (real)bins;}
      else if (geometry == 2) {dr = Li[1] / (real)bins;}
      else {dr = Li[1] / (real)bins;}

      real * histogram = 0;
      histogram = new real[bins];
      for(int i=0;i<bins;i++) histogram[i]=0.0;

      CellList realCells = system.storage->getRealCells();

      real pos = 0.0;
      int bin = 0;

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        if(particlelist.count(cit->id())){
           if (geometry == 1) {pos = cit->position()[0];}
           else if (geometry == 2) {pos = cit->position()[1];}
           else {pos = cit->position()[2];}

           if (geometry == 1) {
             while(pos < 0.0) pos = pos + Li[0];
             while(pos > Li[0]) pos = pos - Li[0];
           }

           if (geometry == 2) {
             while(pos < 0.0) pos = pos + Li[1];
             while(pos > Li[1]) pos = pos - Li[1];
           }

           if (geometry == 3) {
             while(pos < 0.0) pos = pos + Li[2];
             while(pos > Li[2]) pos = pos - Li[2];
           }

           bin = floor (pos/dr);
           histogram[bin] += 1.0;
        }
      }

      real * totHistogram = 0;
      totHistogram = new real[bins];
      for(int i=0;i<bins;i++) totHistogram[i]=0.0;

      boost::mpi::all_reduce(*mpiWorld, histogram, bins, totHistogram, std::plus<real>());

      python::list pyli;
      for(int i=0; i < bins; i++){
        pyli.append( totHistogram[i] );
      }

      delete[] totHistogram;
      totHistogram = 0;
      delete[] histogram;
      histogram = 0;

      return pyli;
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real ParticleDistribution::compute() const {
      return -1.0;
    }

    void ParticleDistribution::registerPython() {
      using namespace espressopp::python;
      class_<ParticleDistribution, bases< Observable > >
        ("analysis_ParticleDistribution", init< shared_ptr< System >, int >())
        .def("addPID", &ParticleDistribution::addPID)
        .def("compute", &ParticleDistribution::computeArray)
      ;
    }
  }
}
