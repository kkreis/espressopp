/*
  Copyright (C) 2012,2013,2014,2015
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
#include "mpi.hpp"
#include "AdressPI.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "FixedTupleListAdress.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"
#include <iomanip>

namespace espressopp {

  namespace integrator {

    using namespace espressopp::iterator;

    AdressPI::AdressPI(shared_ptr<System> _system, shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList, int _ntrotter, real _clmassmultiplier, bool _KTI /*= false*/)
        : Extension(_system), verletList(_verletList), fixedtupleList(_fixedtupleList), ntrotter(_ntrotter), clmassmultiplier(_clmassmultiplier), KTI(_KTI){
        LOG4ESPP_INFO(theLogger, "construct AdressPI");
        type = Extension::Adress;

        // AdResS PI stuff
        dhy = verletList->getHy();
        pidhy2 = M_PI/(dhy * 2.0);
        dex = verletList->getEx();
        dex2 = dex * dex;
        dexdhy = dex + verletList->getHy();
        dexdhy2 = dexdhy * dexdhy;

        //communicateAdrPositions();
    }


    AdressPI::~AdressPI() {
      LOG4ESPP_INFO(theLogger, "~AdressPI");
      disconnect();
    }

    void AdressPI::disconnect(){
        _SetPosVel.disconnect();
        // _initForces.disconnect();
        _setweights.disconnect();
        _aftCalcF.disconnect();
    }

    void AdressPI::connect() {

        // connection to after runInit()
        _SetPosVel = integrator->runInit.connect(
                boost::bind(&AdressPI::SetPosVel, this), boost::signals2::at_front);

        // // connection to after initForces()
        // _initForces = integrator->aftInitF.connect(
        //         boost::bind(&AdressPI::initForces, this), boost::signals2::at_front);

        // connection to inside of integrate1()
        _setweights = integrator->recalc2.connect(
                boost::bind(&AdressPI::setweights, this), boost::signals2::at_front);

        // connection to after _aftCalcFPI()
        _aftCalcF = integrator->aftCalcFPI.connect(
                boost::bind(&AdressPI::aftCalcF, this), boost::signals2::at_back);
    }


    void AdressPI::SetPosVel(){

        System& system = getSystemRef();

        // Set the positions and velocity of CG particles & update weights.
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {

              Particle &vp = *cit;

              FixedTupleListAdress::iterator it3;
              it3 = fixedtupleList->find(&vp);

              if (it3 != fixedtupleList->end()) {

                  std::vector<Particle*> atList;
                  atList = it3->second;

                  // Compute center of mass
                  Real3D cmp(0.0, 0.0, 0.0); // center of mass position
                  Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
                  for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                      Particle &at = **it2;
                      cmp += at.position();
                      cmv += at.velocity();
                      // if(at.pib() == 1)
                      // {
                      //     cmv = (1.0/sqrt(ntrotter)) * at.velocity(); //?????
                      // }
                  }
                  cmv /= ntrotter;
                  cmp /= ntrotter;

                  // update (overwrite) the position and velocity of the VP
                  vp.position() = cmp;
                  vp.velocity() = cmv;

                  if (KTI == false) {
                      // calculate distance to nearest adress particle or center
                      std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                      Real3D pa = **it2; // position of adress particle
                      Real3D d1(0.0, 0.0, 0.0);
                      //Real3D d1 = vp.position() - pa;                                                      // X SPLIT VS SPHERE CHANGE
                      verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                      real min1sq;
                      //real d1 = vp.position()[0] - pa[0];                                                // X SPLIT VS SPHERE CHANGE
                      //real min1sq = d1.sqr();  // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                      //real min1sq = d1*d1;   // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                      if (verletList->getAdrRegionType()) { // spherical adress region
                        min1sq = d1.sqr(); // set min1sq before loop
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                             pa = **it2;
                             verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                             real distsq1 = d1.sqr();
                             if (distsq1 < min1sq) min1sq = distsq1;
                        }
                      }
                      else { //slab-type adress region
                        min1sq = d1[0]*d1[0];   // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                             pa = **it2;
                             //d1 = vp.position() - pa;                                                          // X SPLIT VS SPHERE CHANGE
                             //d1 = vp.position()[0] - pa[0];
                             verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);        // X SPLIT VS SPHERE CHANGE
                             //real distsq1 = d1.sqr();                                                          // X SPLIT VS SPHERE CHANGE
                             real distsq1 = d1[0]*d1[0];                                                           // X SPLIT VS SPHERE CHANGE
                             //std::cout << pa << " " << sqrt(distsq1) << "\n";
                             if (distsq1 < min1sq) min1sq = distsq1;

                        }
                      }

                      real w = weight(min1sq);
                      vp.lambda() = w;
                      real wDeriv = weightderivative(sqrt(min1sq));
                      vp.lambdaDeriv() = wDeriv;
                      vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
                  }

              }
              else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return;
              }


        }

    }


    // void AdressPI::initForces(){

    //     System& system = getSystemRef();

    //     // AT reals
    //     ParticleList& adrATparticles = system.storage->getAdrATParticles();
    //     for (std::vector<Particle>::iterator it = adrATparticles.begin();
    //             it != adrATparticles.end(); ++it) {
    //         it->force() = 0.0;
    //         it->forcem() = 0.0;
    //         // it->drift() = 0.0;
    //     }

    //     // AT ghosts
    //     typedef std::list<ParticleList> ParticleListAdr;
    //     ParticleListAdr& adrATparticlesG = system.storage->getAdrATParticlesG();
    //     for (ParticleListAdr::iterator it = adrATparticlesG.begin();
    //             it != adrATparticlesG.end(); ++it) {

    //         for (ParticleList::iterator it2 = it->begin();
    //                 it2 != it->end(); ++it2) {
    //             it2->force() = 0.0;
    //             it2->forcem() = 0.0;
    //             // it2->drift() = 0.0;
    //         }
    //     }
    // }

    void AdressPI::setweights(){

        System& system = getSystemRef();

        // Set the positions and velocity of CG particles & update weights.
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {


              Particle &vp = *cit;

              FixedTupleListAdress::iterator it3;
              it3 = fixedtupleList->find(&vp);

              if (it3 != fixedtupleList->end()) {

                  if (KTI == false) {

                      // calculate distance to nearest adress particle or center
                      std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                      Real3D pa = **it2; // position of adress particle
                      Real3D d1(0.0, 0.0, 0.0);
                      real min1sq;
                      //Real3D d1 = vp.position() - pa;                                                      // X SPLIT VS SPHERE CHANGE
                      verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                      //real d1 = vp.position()[0] - pa[0];                                                // X SPLIT VS SPHERE CHANGE
                      //real min1sq = d1.sqr();  // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                      //real min1sq = d1*d1;   // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                      if (verletList->getAdrRegionType()) { // spherical adress region
                        min1sq = d1.sqr(); // set min1sq before loop
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                             pa = **it2;
                             verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                             real distsq1 = d1.sqr();
                             if (distsq1 < min1sq) min1sq = distsq1;
                        }
                      }
                      else { //slab-type adress region
                        min1sq = d1[0]*d1[0];   // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                             pa = **it2;
                             //d1 = vp.position() - pa;                                                          // X SPLIT VS SPHERE CHANGE
                             //d1 = vp.position()[0] - pa[0];
                             verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);        // X SPLIT VS SPHERE CHANGE
                             //real distsq1 = d1.sqr();                                                          // X SPLIT VS SPHERE CHANGE
                             real distsq1 = d1[0]*d1[0];                                                           // X SPLIT VS SPHERE CHANGE
                             //std::cout << pa << " " << sqrt(distsq1) << "\n";
                             if (distsq1 < min1sq) min1sq = distsq1;
                        }
                      }


                      real w = weight(min1sq);
                      vp.lambda() = w;
                      real wDeriv = weightderivative(min1sq);
                      vp.lambdaDeriv() = wDeriv;
                      vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
                  }

              }
              else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return;
              }


        }
    }

    // AdResS Weighting function
    real AdressPI::weight(real distanceSqr){
        if (dex2 >= distanceSqr) return 1.0;
        else if (dexdhy2 <= distanceSqr) return 0.0;
        else {
            real argument = sqrt(distanceSqr) - dex;
            return pow(cos(pidhy2 * argument),2.0);
        }
    }
    // AdResS Weight Derivative function
    real AdressPI::weightderivative(real distanceSqr){
        if (dex2 >= distanceSqr) return 0.0;
        else if (dexdhy2 <= distanceSqr) return 0.0;
        else{
            real argument = sqrt(distanceSqr) - dex;
            return -1.0 * pidhy2 * 2.0 * cos(pidhy2*argument) * sin(pidhy2*argument);
        }
    }

    void AdressPI::aftCalcF(){
        System& system = getSystemRef();
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {

        Particle &vp = *cit;

        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);

        if (it3 != fixedtupleList->end()) {

            std::vector<Particle*> atList;
            atList = it3->second;

            // update force of AT particles belonging to a VP
            Real3D vpfm = vp.force();
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                                 it2 != atList.end(); ++it2) {
                Particle &at = **it2;
                at.force() += (1.0/ntrotter)*vpfm;
            }
        }
        else { // this should not happen
            std::cout << " particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
            std::cout << " (" << vp.position() << ")\n";
            exit(1);
            return;
        }
      }
    }



    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void AdressPI::registerPython() {
      using namespace espressopp::python;

      class_<AdressPI, shared_ptr<AdressPI>, bases<Extension> >
        ("integrator_AdressPI", init<shared_ptr<System>, shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress>, int, real, bool >())
        .def("connect", &AdressPI::connect)
        .def("disconnect", &AdressPI::disconnect)
        ;
    }

  }

}
