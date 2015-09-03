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
#include <cmath>
#include "ClassicalPressurePI.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "VerletList.hpp"

using namespace espressopp;
using namespace iterator;
using namespace interaction;

namespace espressopp {
  namespace analysis {
    real ClassicalPressurePI::compute() const {

      System& system = getSystemRef();
      mpi::communicator communic = *system.comm;

      // determine volume of the box
      Real3D Li = system.bc->getBoxL();
      real tripleV = 3.0 * Li[0] * Li[1] * Li[2];

      // compute the kinetic contriubtion (2/3 \sum 1/2mv^2)
      real p_kinetic;
      real v2sum = 0.0;
      real v2 = 0.0;

      CellList realCells = system.storage->getRealCells();
           
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        const Particle& p = *cit;
        v2 += p.mass() * (p.velocity() * p.velocity());
      }
          
            
      boost::mpi::all_reduce(*getSystem()->comm, v2, v2sum, std::plus<real>());
      //boost::mpi::all_reduce( communic, v2, v2sum, std::plus<real>());
      p_kinetic = v2sum;
      
      // compute the short-range nonbonded contribution
      real rij_dot_Fij = 0.0;
      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        rij_dot_Fij += srIL[j]->computeVirial();
        //std::cout << "srIL[" << j << "]: " << srIL[j]->computeVirial() << "\n";
      }

      real p_nonbonded = rij_dot_Fij;
      
      //DEBUG (This pressure calculation seems incorrect at least for liquid water.
      // Be careful when doing system with electrostatic interactions)
      //return (p_kinetic + p_nonbonded) / tripleV;   TEST FOR CONSTRAINTS
      /*std::cout << "tripleV: " << tripleV << "\n";
      std::cout << "p_kinetic: " << p_kinetic << "\n";      
      std::cout << "p_nonbonded: " << p_nonbonded << "\n";
      std::cout << "p_kinetic/tripleV: " << p_kinetic/tripleV << "\n";
      std::cout << "p_nonbonded/tripleV: " << p_nonbonded/tripleV << "\n";     */
      //DEBUG
      
      return (p_kinetic + p_nonbonded) / tripleV;     
      //return ((p_kinetic/tripleV)*(3.0/2.0) + p_nonbonded/tripleV);  //SETTLE CONSTAINTS
    }

    void ClassicalPressurePI::registerPython() {
      using namespace espressopp::python;
      class_<ClassicalPressurePI, bases< Observable > >
        ("analysis_ClassicalPressurePI", init< shared_ptr< System > >())
      ;
    }
  }
}