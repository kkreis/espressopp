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

//#include <iomanip>
#include "python.hpp"
#include "VelocityVerletPI.hpp"
#include <iomanip>
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"
#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "VerletListAdress.hpp"

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER( name)
#endif

namespace espressopp {
  using namespace std;
  namespace integrator {
    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    LOG4ESPP_LOGGER(VelocityVerletPI::theLogger, "VelocityVerletPI");

    VelocityVerletPI::VelocityVerletPI(shared_ptr< System > system, shared_ptr<VerletListAdress> _verletList)
    : MDIntegrator(system), verletList(_verletList)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerletPI");
      resortFlag = true;
      maxDist    = 0.0;
      dt2 = dt*sStep;
      dt3 = dt*mStep*sStep;
              
      // Initialize variables
      temperature = 0.0;
      gamma = 0.0;
      ntrotter = 0;
      CMDparameter = 1.0;
      PILElambda = 0.5;
      speedup = false;
      constkinmass = false;
      realkinmass = false;
      KTI = false;
      centroidthermostat = true;
      PILE = false;
      sStep = 0;
      mStep = 0;
      omega2 = 0.0;
      clmassmultiplier = 100.0;
      
      UpdateCounter = 0;
      
      // AdResS PI stuff
      dhy = verletList->getHy();
      pidhy2 = M_PI/(dhy * 2.0);
      dex = verletList->getEx();
      dex2 = dex * dex;
      dexdhy = dex + verletList->getHy();
      dexdhy2 = dexdhy * dexdhy;
      
      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;
    }

    VelocityVerletPI::~VelocityVerletPI()
    {
      LOG4ESPP_INFO(theLogger, "free VelocityVerletPI");
    }

    void VelocityVerletPI::run(int nsteps)
    {
    	/*std::cout << "PILE: " << PILE << "\n";
    	std::cout << "realkinmass: " << realkinmass << "\n";
    	std::cout << "centroidthermostat: " << centroidthermostat << "\n";
    	std::cout << "CMDparameter: " << CMDparameter << "\n\n";*/

      // Check if Eigenvalues start with zero:
      if(Eigenvalues.at(0) > 0.00000000001){
        std::cout << "Eigenvalues don't start with zero!\n";
        exit(1);
        return;
      }
      
      if(PILE == true && realkinmass == false){
        std::cout << "Using the path integral Langevin scheme (PILE) only makes sense when using real masses for kinetic masses!\n";
        exit(1);
        return;
      }

      System& system = getSystemRef();
      storage::Storage& storage = *system.storage;
      real skinHalf = 0.5 * system.getSkin();

      // signal
      runInit();
      
      // Before start make sure that particles are on the right processor
      if (resortFlag) {
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
      }
      
      transPos2(); // Update mode positions.
  
      // outer loop
      for (int i = 0; i < nsteps; i++) {

        // FULL PSEUDOCODE SCHEME
        // Initialize run from input (real AT positions): update VP centroid position and velocity (0.0), calculate weights. 
        // Update/initialize mode positions
        // Calculate all forces and transform them (don't forget the VP forces from CL region)
        // Calculate centroid forces from FEC - Extension!!!
        
            // main step loop
    
            // Propagate all (mode) velocities according to slow forces 1 step with dt3/2

                 // mstep loop

                 // Propagate all (mode) velocities according to medium forces 1 step with dt2/2

                       // sstep loop

                       // Propagate higher (mode) velocities according to fast forces 1 step with dt/4
                       // Propagate centroid (mode) velocities according to fast forces 1 step with dt/2
                       // Propagate higher (mode) velocities according to fast forces 1 step with dt/4

                       // Propagate all (mode) positions with velocities with dt/2

                       // Propagate all (mode) velocities according to Ornstein Uhlenbeck with dt

                       // Propagate all (mode) positions with velocities with dt/2

                       // Update real positions - don't forget the VP centroid particle (both velocity and position).
                       // Check the VerletList, maybe rebuild it.
                       // Recalculate weights and mass
                       
                       // Calculate fast forces.
                       // Propagate higher (mode) velocities according to fast forces 1 step with dt/4
                       // Propagate centroid (mode) velocities according to fast forces 1 step with dt/2
                       // Propagate higher (mode) velocities according to fast forces 1 step with dt/4

                       // sstep loop

                // Calculate medium forces and transform them (don't forget the VP forces from CL region)
                // Propagate all (mode) velocities according to medium forces 1 step with dt2/2

                // mstep loop

            // Calculate slow forces and transform them (don't forget the VP forces from CL region)
            // Calculate centroid forces from FEC - Extension!!!
            // Propagate all (mode) velocities according to slow forces 1 step with dt3/2
        
            // main step loop
        
        // (Update VP velocities - optional)
        // FULL PSEUDOCODE SCHEME
        
        if(i==0){ updateForces(3); }
        //updateForces(3);
        integrateV1(3);
        
        // mstep loop
        for (int j = 0; j < mStep; j++){
            
            if(j==0){ updateForces(2); }
            //updateForces(2);
            integrateV1(2);
            
            // lstep loop
            for (int k= 0; k < sStep; k++){
                
                if(k==0){ updateForces(1); }                
                //updateForces(1);
                integrateV2();
                
                integrateModePos();
                OUintegrate();
                integrateModePos();
                
                // Update real positions - don't forget the VP centroid particle (both velocity and position).
                // Recalculate weights and mass
                transPos1();
                // Check the VerletList, maybe rebuild it.
                if (maxDist > skinHalf) resortFlag = true;
                if (resortFlag) {
                    storage.decompose();
                    transPos2();
                    maxDist  = 0.0;
                    resortFlag = false;
                }
                
                updateForces(1);
                integrateV2();
            
            }
            
            updateForces(2);
            integrateV1(2);
        
        }
        
        updateForces(3);
        integrateV1(3);
        
        step++;
      }

    }


    using namespace boost::python;


    void VelocityVerletPI::integrateV1(int t)
    {   
        real half_dt = 0.0;
        if (t==2){
            half_dt = 0.5 * dt2;
        }
        else if (t==3){
            half_dt = 0.5 * dt3;
        }
        else{
            std::cout << "integrateV1 routine in VelocityVerletPI integrator received wrong integer.\n";
            exit(1);
            return;
        }
                
        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) { 
            
            Particle &vp = *cit;
            
            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);
            if (it3 != fixedtupleList->end()) {
                std::vector<Particle*> atList;
                atList = it3->second;

                for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                    Particle &at = **it2;

					if((at.pib() != 1) && (speedup == true) && (vp.lambda() < 0.000000001)){
						continue;
					}else{
						at.velocity() += half_dt * at.forcem();
						if(at.pib() == 1)
						{
							vp.velocity() = (1.0/sqrt(ntrotter)) * at.velocity();
						}
					}
                    
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
    
    
    void VelocityVerletPI::integrateV2()
    {     
        real half_dt = 0.5 * dt;
        real half_dt4 = 0.25 * dt;
                
        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();  //!!!!
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) { 
            
            Particle &vp = *cit;
            
            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);
            if (it3 != fixedtupleList->end()) {
                std::vector<Particle*> atList;
                atList = it3->second;

                if (constkinmass == false){ // kinetic mass also changes - hence, second decomposition level required

					// First half_dt4 loop without centroid mode
                	if((speedup == false) || (vp.lambda() > 0.000000001)){

						for (std::vector<Particle*>::iterator it2 = atList.begin();
											   it2 != atList.end(); ++it2) {
							Particle &at = **it2;

							if(at.pib() == 1){
								continue;
							}
							else if(at.pib() > 1 && at.pib() <= ntrotter){
									at.velocity() += half_dt4 * at.forcem();
							}
							else{
								 std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
								 exit(1);
								 return;
							}
						}

                	}

					// half_dt2 centroid mode
					for (std::vector<Particle*>::iterator it2 = atList.begin();
										   it2 != atList.end(); ++it2) {
						Particle &at = **it2;

						if(at.pib() == 1){

							if(vp.lambda()<1.0 && vp.lambda()>0.0){

								// Get the inner factor
								real xi = 0.0;
								for (std::vector<Particle*>::iterator it5 = atList.begin();
											 it5 != atList.end(); ++it5) {
									Particle &at2 = **it5;
									if(at2.pib() != 1){
										if(realkinmass == false){
											xi += at2.velocity().sqr()/(Eigenvalues.at(at2.pib()-1));
										}
										else{
											xi += at2.velocity().sqr();
										}
									}
								}

								// calculate distance to nearest adress particle or center
								std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
								Real3D pa = **it2; // position of adress particle
								Real3D mindriftforce(0.0, 0.0, 0.0);
								//Real3D mindriftforce = vp.position() - pa;                                                           // X SPLIT VS SPHERE CHANGE
								//real mindriftforce = vp.position()[0] - pa[0];                                                         // X SPLIT VS SPHERE CHANGE
								verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, vp.position(), pa);
								real min1sq = mindriftforce[0]*mindriftforce[0]; // mindriftforce.sqr(); // set min1sq before loop
								++it2;
								for (; it2 != verletList->getAdrPositions().end(); ++it2) {
									   pa = **it2;
									   Real3D driftforce(0.0, 0.0, 0.0);
									   //Real3D driftforce = vp.position() - pa;                                                          // X SPLIT VS SPHERE CHANGE
									   //real driftforce = vp.position()[0] - pa[0];                                                    // X SPLIT VS SPHERE CHANGE
									   verletList->getSystem()->bc->getMinimumImageVector(driftforce, vp.position(), pa);
									   //real distsq1 = driftforce.sqr();//driftforce*driftforce; //driftforce.sqr();                     // X SPLIT VS SPHERE CHANGE
									   real distsq1 = driftforce[0]*driftforce[0];                                                          // X SPLIT VS SPHERE CHANGE
									   //std::cout << pa << " " << sqrt(distsq1) << "\n";
									   if (distsq1 < min1sq) {
											min1sq = distsq1;
											mindriftforce = driftforce;
									   }
								}
								min1sq = sqrt(min1sq);   // distance to nearest adress particle or center
								real mindriftforceX = (1.0/min1sq)*mindriftforce[0];  // normalized driftforce vector

								mindriftforceX *= 0.5*(clmassmultiplier-1.0)*vp.mass()*xi*vp.lambdaDeriv()/(vp.varmass()*vp.varmass()*sqrt(ntrotter)*CMDparameter);                                                              // X SPLIT VS SPHERE CHANGE
								Real3D driftforceadd(mindriftforceX,0.0,0.0);

								at.velocity() += half_dt * (at.forcem() - driftforceadd);
							}
							else{
								at.velocity() += half_dt * at.forcem();
							}
							vp.velocity() = (1.0/sqrt(ntrotter)) * at.velocity();
							break;

						}
						else if(at.pib() > 1 && at.pib() <= ntrotter){
							continue;
						}
						else{
							 std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
							 exit(1);
							 return;
						}
					}

					// Second half_dt4 loop without centroid mode
					if((speedup == false) || (vp.lambda() > 0.000000001)){

						for (std::vector<Particle*>::iterator it2 = atList.begin();
											   it2 != atList.end(); ++it2) {
							Particle &at = **it2;

							if(at.pib() == 1){
								continue;
							}
							else if(at.pib() > 1 && at.pib() <= ntrotter){
								at.velocity() += half_dt4 * at.forcem();
							}
							else{
								 std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
								 exit(1);
								 return;
							}
						}

					}

                }
                else{ // constant kinetic mass - hence, no second decomposition level required

					for (std::vector<Particle*>::iterator it2 = atList.begin();
										   it2 != atList.end(); ++it2) {
						Particle &at = **it2;

						// MOMENTUM APPROACH:
						if((at.pib() != 1) && (speedup == true) && (vp.lambda() < 0.000000001)){
							continue;
						}else{
							at.velocity() += half_dt * at.forcem();
							if(at.pib() == 1)
							{
								vp.velocity() = (1.0/sqrt(ntrotter)) * at.velocity();
							}
						}
					}

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
    
    
    void VelocityVerletPI::integrateModePos()
    {

        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells(); //!!!!
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            Particle &vp = *cit;

            if (constkinmass == false){ // kinetic mass also changes - hence, second decomposition level required

				real half_dt4 = 0.25 * dt / (vp.mass());

				FixedTupleListAdress::iterator it3;
				it3 = fixedtupleList->find(&vp);
				if (it3 != fixedtupleList->end()) {
					std::vector<Particle*> atList;
					atList = it3->second;

					// First half_dt4 loop with only centroid mode
					for (std::vector<Particle*>::iterator it2 = atList.begin();
										   it2 != atList.end(); ++it2) {
						Particle &at = **it2;

						if(at.pib() == 1){
							at.modepos() += half_dt4 * at.velocity();
							if (KTI == false){
								// Update resolution and variable masses
								// calculate distance to nearest adress particle or center
								std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
								Real3D pa = **it2; // position of adress particle
								Real3D d1(0.0, 0.0, 0.0);
								real min1sq;
								verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);

								if (verletList->getAdrRegionType()) { // spherical adress region
								  min1sq = d1.sqr(); // set min1sq before loop
								  ++it2;
								  for (; it2 != verletList->getAdrPositions().end(); ++it2) {
									   pa = **it2;
									   verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
									   real distsq1 = d1.sqr();
									   if (distsq1 < min1sq) min1sq = distsq1;
								  }
								}
								else { //slab-type adress region
								  min1sq = d1[0]*d1[0];   // set min1sq before loop
								  ++it2;
								  for (; it2 != verletList->getAdrPositions().end(); ++it2) {
									   pa = **it2;
									   verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
									   real distsq1 = d1[0]*d1[0];
									   if (distsq1 < min1sq) min1sq = distsq1;
								  }
								}
								real w = weight(min1sq);
								vp.lambda() = w;
								real wDeriv = weightderivative(min1sq);
								vp.lambdaDeriv() = wDeriv;
								vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
							}

							break;
						}else{
							continue;
						}
					}

					// half_dt2 loop without centroid mode
					real half_dt = 0.5 * dt / (CMDparameter * vp.varmass());
					for (std::vector<Particle*>::iterator it2 = atList.begin();
										   it2 != atList.end(); ++it2) {
						Particle &at = **it2;

						if((at.pib() == 1) || ((speedup == true) && (vp.lambda() < 0.000000001) && (at.pib() != 1))){
							continue;
						}
						else{
							if(realkinmass == false){
								at.modepos() += half_dt * at.velocity() / (Eigenvalues.at(at.pib()-1));
							}
							else{
								at.modepos() += half_dt * at.velocity();
							}
						}
					}

					// Second half_dt4 loop with only centroid mode
					for (std::vector<Particle*>::iterator it2 = atList.begin();
										   it2 != atList.end(); ++it2) {
						Particle &at = **it2;

						if(at.pib() == 1){
							at.modepos() += half_dt4 * at.velocity();
							if (KTI == false){
								// Update resolution and variable masses
								// calculate distance to nearest adress particle or center
								std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
								Real3D pa = **it2;
								Real3D d1(0.0, 0.0, 0.0);
								real min1sq;
								verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);

								if (verletList->getAdrRegionType()) { // spherical adress region
								  min1sq = d1.sqr(); // set min1sq before loop
								  ++it2;
								  for (; it2 != verletList->getAdrPositions().end(); ++it2) {
									   pa = **it2;
									   verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
									   real distsq1 = d1.sqr();
									   if (distsq1 < min1sq) min1sq = distsq1;
								  }
								}
								else { //slab-type adress region
								  min1sq = d1[0]*d1[0];   // set min1sq before loop
								  ++it2;
								  for (; it2 != verletList->getAdrPositions().end(); ++it2) {
									   pa = **it2;
									   verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
									   real distsq1 = d1[0]*d1[0];
									   if (distsq1 < min1sq) min1sq = distsq1;
								  }
								}
								real w = weight(min1sq);
								vp.lambda() = w;
								real wDeriv = weightderivative(min1sq);
								vp.lambdaDeriv() = wDeriv;
								vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
							}

							break;
						}else{
							continue;
						}
					}

				}
				else { // this should not happen
					  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
					  std::cout << " (" << vp.position() << ")\n";
					  exit(1);
					  return;
				}

            }
            else { // constant kinetic mass - hence, no second decomposition level required

            	real half_dt = 0.5 * dt / (vp.mass());

				FixedTupleListAdress::iterator it3;
				it3 = fixedtupleList->find(&vp);
				if (it3 != fixedtupleList->end()) {
					std::vector<Particle*> atList;
					atList = it3->second;

					for (std::vector<Particle*>::iterator it2 = atList.begin();
										   it2 != atList.end(); ++it2) {
						Particle &at = **it2;

						if(at.pib() == 1){
							at.modepos() += half_dt * at.velocity();

							if (KTI == false){
								// Update resolution and variable masses
								// calculate distance to nearest adress particle or center
								std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
								Real3D pa = **it2; // position of adress particle
								Real3D d1(0.0, 0.0, 0.0);
								real min1sq;
								verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);

								if (verletList->getAdrRegionType()) { // spherical adress region
								  min1sq = d1.sqr(); // set min1sq before loop
								  ++it2;
								  for (; it2 != verletList->getAdrPositions().end(); ++it2) {
									   pa = **it2;
									   verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);
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
									   verletList->getSystem()->bc->getMinimumImageVector(d1, (1.0/sqrt(ntrotter))*at.modepos(), pa);   // X SPLIT VS SPHERE CHANGE
									   //real distsq1 = d1.sqr();                                                          // X SPLIT VS SPHERE CHANGE
									   real distsq1 = d1[0]*d1[0];                                                           // X SPLIT VS SPHERE CHANGE
									   if (distsq1 < min1sq) min1sq = distsq1;
								  }
								}
								real w = weight(min1sq);
								vp.lambda() = w;
								real wDeriv = weightderivative(min1sq);
								vp.lambdaDeriv() = wDeriv;
								vp.varmass() = vp.mass()*( w*(1.0-clmassmultiplier) + clmassmultiplier );
							}

						}else if(at.pib() > 1 && at.pib() <= ntrotter){
							if((speedup == true) && (vp.lambda() < 0.000000001)){
								continue;
							}else{
								if(realkinmass == false){
									at.modepos() += half_dt * at.velocity() / (CMDparameter * Eigenvalues.at(at.pib()-1));
								}
								else{
									at.modepos() += half_dt * at.velocity() / CMDparameter;
								}
							}
						}else{
							std::cout << "at.pib() outside of trotter range in integrateModePos routine (VelocityVerletPI) \n";
	                        exit(1);
	                        return;
						}
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
        
    }
    
    
    void VelocityVerletPI::OUintegrate()
    {   
        real prefac1 = exp(-gamma*dt);
        real prefac2 = sqrt(12.0*temperature*( 1.0-exp(-2.0*gamma*dt) ));
        // careful with 12... it's because the uniform distribution needs to have the same width as the Gaussian one.
        // also note that kb is inside temperature, we use gromacs units
        
        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        
        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) { 
            
            Particle &vp = *cit;
            
            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);
            if (it3 != fixedtupleList->end()) {
                std::vector<Particle*> atList;
                atList = it3->second;

                for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                    Particle &at = **it2;
                    
                    if(at.pib() == 1){
                    	if(centroidthermostat == true || vp.lambda() < 1.0){
							Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);
							at.velocity() = prefac1 * at.velocity() + prefac2*sqrt(vp.mass()) * ranval;
							vp.velocity() = (1.0/sqrt(ntrotter)) * at.velocity();
                    	}
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);
                        if(PILE == false){
							if(constkinmass == false){
								if(realkinmass == false){
									at.velocity() = prefac1 * at.velocity() + (prefac2*sqrt(CMDparameter*vp.varmass()*Eigenvalues.at(at.pib()-1))) * ranval;
								}
								else{
									at.velocity() = prefac1 * at.velocity() + (prefac2*sqrt(CMDparameter*vp.varmass())) * ranval;
								}
							}
							else{
								if(realkinmass == false){
									at.velocity() = prefac1 * at.velocity() + (prefac2*sqrt(CMDparameter*vp.mass()*Eigenvalues.at(at.pib()-1))) * ranval;
								}
								else{
									at.velocity() = prefac1 * at.velocity() + (prefac2*sqrt(CMDparameter*vp.mass())) * ranval;
								}
							}
                        }
                        else{
                        	real modegamma = 2.0 * PILElambda * sqrt(omega2) * Eigenvalues.at(at.pib()-1);
                        	real prefac1_PILE = exp(-modegamma*dt);
                        	real prefac2_PILE = sqrt(12.0*temperature*( 1.0-exp(-2.0*modegamma*dt) ));
							if(constkinmass == false){
								at.velocity() = prefac1_PILE * at.velocity() + (prefac2_PILE*sqrt(CMDparameter*vp.varmass())) * ranval;
							}
							else{
								at.velocity() = prefac1_PILE * at.velocity() + (prefac2_PILE*sqrt(CMDparameter*vp.mass())) * ranval;
							}
                        }
                    }
                    else{
                         std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
                         exit(1);
                         return;
                    }
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
    
    void VelocityVerletPI::transPos1(){
        real maxSqDist = 0.0;
        
        System& system = getSystemRef();
        CellList localCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

        for(CellListIterator cit(localCells); !cit.isDone(); ++cit){
            Particle &vp = *cit;
            FixedTupleListAdress::iterator it2;
            it2 = fixedtupleList->find(&vp);

            if(it2 != fixedtupleList->end()){
                  std::vector<Particle*> atList;
                  atList = it2->second;

                  Real3D oldpos = vp.position();

                  for (std::vector<Particle*>::iterator it3 = atList.begin();
                                       it3 != atList.end(); ++it3) {
                      Particle &at = **it3;

                      Real3D zero(0.0,0.0,0.0);
                      at.position()=zero;

                      for (std::vector<Particle*>::iterator it5 = atList.begin();
                                   it5 != atList.end(); ++it5) {
                          Particle &at2 = **it5;

                          if(at.pib() <= ntrotter){
                              at.position()+= at2.modepos()*Tvectors.at(at.pib()-1).at(at2.pib()-1);
                              if((at2.pib() == 1) && (at.pib() == 1)){
                                  vp.position() = (1.0/sqrt(ntrotter)) * at2.modepos();
                                  vp.velocity() = (1.0/sqrt(ntrotter)) * at2.velocity();
                              }
                          }
                          else{
                               std::cout << "at2.pib() outside of trotter range in transPos1 routine (VelocityVerletPI) \n";
                               exit(1);
                               return;
                          }

                      }

                  }
                  
                  real sqDist = (oldpos-vp.position()).sqr();
                  maxSqDist = std::max(maxSqDist, sqDist);
            }
            else{
                 std::cout << "VP particle not found in tuples (VelocityVerletPI): " << vp.id() << "-" <<
                         vp.ghost();
                 std::cout << " (" << vp.position() << ")\n";
                 exit(1);
                 return;      
            }
                
        }

        real maxAllSqDist;
        mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());
        maxDist+=sqrt(maxAllSqDist);
        
    }
    
    
    void VelocityVerletPI::transPos2(){
        System& system = getSystemRef();
        CellList localCells = system.storage->getRealCells(); //!!!!!!
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

        for(CellListIterator cit(localCells); !cit.isDone(); ++cit){
            Particle &vp = *cit;
            FixedTupleListAdress::iterator it2;
            it2 = fixedtupleList->find(&vp);

            if(it2 != fixedtupleList->end()){
                  std::vector<Particle*> atList;
                  atList = it2->second;
                  for (std::vector<Particle*>::iterator it3 = atList.begin();
                                       it3 != atList.end(); ++it3) {
                      Particle &at = **it3;
                      Real3D zero(0.0,0.0,0.0);
                      at.modepos()=zero;
                      
                      for (std::vector<Particle*>::iterator it5 = atList.begin();
                                   it5 != atList.end(); ++it5) {
                          Particle &at2 = **it5;

                          if(at.pib() <= ntrotter){
                              at.modepos()+= at2.position()*Eigenvectors.at(at.pib()-1).at(at2.pib()-1);
                          }
                          else{
                               std::cout << "at2.pib() outside of trotter range in transPos2 routine (VelocityVerletPI) \n";
                               exit(1);
                               return;
                          }

                    }
                  } 
            }
            else{
                 std::cout << "VP particle not found in tuples (VelocityVerletPI): " << vp.id() << "-" <<
                         vp.ghost();
                 std::cout << " (" << vp.position() << ")\n";
                 exit(1);
                 return;      
            }
                
        }
        
    }
    
    
    void VelocityVerletPI::transForces(){
        System& system = getSystemRef();
        CellList localCells = system.storage->getLocalCells(); //!!!!
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

        for(CellListIterator cit(localCells); !cit.isDone(); ++cit){
            Particle &vp = *cit;
            FixedTupleListAdress::iterator it2;
            it2 = fixedtupleList->find(&vp);

            if(it2 != fixedtupleList->end()){
                  std::vector<Particle*> atList;
                  atList = it2->second;
                  for (std::vector<Particle*>::iterator it3 = atList.begin();
                                       it3 != atList.end(); ++it3) {
                      Particle &at = **it3;

                      for (std::vector<Particle*>::iterator it5 = atList.begin();
                                   it5 != atList.end(); ++it5) {
                          Particle &at2 = **it5;

                          if(at.pib() == 1){
                        	  at.forcem()+= (1.0/sqrt(ntrotter))* at2.force();
                          }
                          else if(at.pib() > 1 && at.pib() <= ntrotter){
                        	  at.forcem()+= at2.force()*Eigenvectors.at(at.pib()-1).at(at2.pib()-1);
                          }
                          else{
                               std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
                               exit(1);
                               return;
                          }

                    } 

                  } 
            }
            else{
                 std::cout << "VP particle not found in tuples (VelocityVerletPI): " << vp.id() << "-" <<
                         vp.ghost();
                 std::cout << " (" << vp.position() << ")\n";
                 exit(1);
                 return;      
            }
                
        }
        
    }
    
    
    void VelocityVerletPI::calcForcesS()
    {
      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;

      for (size_t i = 0; i < srIL.size(); i++) {
        if(srIL[i]->bondType() == Nonbonded){
            srIL[i]->addForces();
        }
      }
    }
    
    
    void VelocityVerletPI::calcForcesM()
    {
      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;

      for (size_t i = 0; i < srIL.size(); i++) {
	    LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
        if(srIL[i]->bondType() == Pair || srIL[i]->bondType() == Angular){
            srIL[i]->addForces();
        }
      }

      // Signal // For FEC
      aftCalcF();
    }
        
    
    void VelocityVerletPI::calcForcesF()
    {
      System& system = getSystemRef();
      CellList localCells = system.storage->getRealCells(); //!!!!
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit){
                Particle &vp = *cit;
                FixedTupleListAdress::iterator it2;
                it2 = fixedtupleList->find(&vp);
                
                if(speedup==true){
                    
                    if(vp.lambda() < 0.000000001){
                        continue;
                    }
                    else{
                        if(it2 != fixedtupleList->end()){
                              std::vector<Particle*> atList;
                              atList = it2->second;
                              for (std::vector<Particle*>::iterator it3 = atList.begin();
                                                   it3 != atList.end(); ++it3) {
                                  Particle &at = **it3;

                                  if(at.pib() > 1 && at.pib() <= ntrotter){
                                    at.forcem() -= at.modepos()*omega2*vp.varmass()*Eigenvalues.at(at.pib()-1);
                                  }
                                  else if(at.pib() == 1){
                                	  if(vp.lambda()<1.0 && vp.lambda()>0.0){
                                        real xi = 0.0;
                                        for (std::vector<Particle*>::iterator it5 = atList.begin();
                                                     it5 != atList.end(); ++it5) {
                                            Particle &at2 = **it5;
                                            if(at2.pib() != 1){
                                            	xi += -0.5*(clmassmultiplier-1.0)*(1.0/sqrt(ntrotter))*vp.lambdaDeriv()*vp.mass()*omega2*at2.modepos().sqr()*Eigenvalues.at(at2.pib()-1);
                                            }
                                        }

                                        // calculate distance to nearest adress particle or center
                                        std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                                        Real3D pa = **it2; // position of adress particle
                                        Real3D mindriftforce(0.0, 0.0, 0.0);
                                        //Real3D mindriftforce = vp.position() - pa;                                                           // X SPLIT VS SPHERE CHANGE
                                        //real mindriftforce = vp.position()[0] - pa[0];                                                         // X SPLIT VS SPHERE CHANGE
                                        verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, vp.position(), pa);
                                        real min1sq = mindriftforce[0]*mindriftforce[0]; // mindriftforce.sqr(); // set min1sq before loop
                                        ++it2;
                                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                                               pa = **it2;
                                               Real3D driftforce(0.0, 0.0, 0.0);
                                               //Real3D driftforce = vp.position() - pa;                                                          // X SPLIT VS SPHERE CHANGE
                                               //real driftforce = vp.position()[0] - pa[0];                                                    // X SPLIT VS SPHERE CHANGE
                                               verletList->getSystem()->bc->getMinimumImageVector(driftforce, vp.position(), pa);
                                               //real distsq1 = driftforce.sqr();//driftforce*driftforce; //driftforce.sqr();                     // X SPLIT VS SPHERE CHANGE
                                               real distsq1 = driftforce[0]*driftforce[0];                                                          // X SPLIT VS SPHERE CHANGE
                                               if (distsq1 < min1sq) {
                                                    min1sq = distsq1;
                                                    mindriftforce = driftforce;
                                               }
                                        }
                                        min1sq = sqrt(min1sq);   // distance to nearest adress particle or center
                                        real mindriftforceX = (1.0/min1sq)*mindriftforce[0];  // normalized driftforce vector
                                        mindriftforceX *= xi;   // get the energy differences which were calculated previously and put in drift force

                                        //mindriftforceX *= vp.lambdaDeriv();  // multiplication with derivative of the weighting function
                                        //vp.force() += mindriftforce;   // add drift force to virtual particles                                                                    // X SPLIT VS SPHERE CHANGE
                                        Real3D driftforceadd(mindriftforceX,0.0,0.0);                                                                                            // X SPLIT VS SPHERE CHANGE
                                        //Real3D driftforceadd(0.0,0.0,0.0);   
                                        //vp.force() += driftforceadd;             // Goes in, if one wants to apply the "normal" drift force - also improve using [0] ...           // X SPLIT VS SPHERE CHANGE


                                        at.forcem() -= driftforceadd;
                                    }
                                  }
                                  else{
                                       std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
                                       exit(1);
                                       return;
                                  }
                              } 
                        }
                        else{
                             std::cout << "VP particle not found in tuples (VelocityVerletPI): " << vp.id() << "-" <<
                                     vp.ghost();
                             std::cout << " (" << vp.position() << ")\n";
                             exit(1);
                             return;      
                        }        
                    }                   
                    
                }
                else{
                    
                    if(it2 != fixedtupleList->end()){
                          std::vector<Particle*> atList;
                          atList = it2->second;
                          for (std::vector<Particle*>::iterator it3 = atList.begin();
                                               it3 != atList.end(); ++it3) {
                              Particle &at = **it3;

                              if(at.pib() > 1 && at.pib() <= ntrotter){
                                at.forcem() -= at.modepos()*omega2*vp.varmass()*Eigenvalues.at(at.pib()-1);
                              }
                              else if(at.pib() == 1){
                                if(vp.lambda()<1.0 && vp.lambda()>0.0)
                                {
                                    real xi = 0.0;
                                    for (std::vector<Particle*>::iterator it5 = atList.begin();
                                                 it5 != atList.end(); ++it5) {
                                        Particle &at2 = **it5;
                                        if(at2.pib() != 1){
                                            xi += at2.modepos().sqr()*Eigenvalues.at(at2.pib()-1);
                                        }
                                    }

                                    // calculate distance to nearest adress particle or center
                                    std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                                    Real3D pa = **it2; // position of adress particle
                                    Real3D mindriftforce(0.0, 0.0, 0.0);
                                    //Real3D mindriftforce = vp.position() - pa;                                                           // X SPLIT VS SPHERE CHANGE
                                    //real mindriftforce = vp.position()[0] - pa[0];                                                         // X SPLIT VS SPHERE CHANGE
                                    verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, vp.position(), pa);
                                    real min1sq = mindriftforce[0]*mindriftforce[0]; // mindriftforce.sqr(); // set min1sq before loop
                                    ++it2;
                                    for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                                           pa = **it2;
                                           Real3D driftforce(0.0, 0.0, 0.0);
                                           //Real3D driftforce = vp.position() - pa;                                                          // X SPLIT VS SPHERE CHANGE
                                           //real driftforce = vp.position()[0] - pa[0];                                                    // X SPLIT VS SPHERE CHANGE
                                           verletList->getSystem()->bc->getMinimumImageVector(driftforce, vp.position(), pa);
                                           //real distsq1 = driftforce.sqr();//driftforce*driftforce; //driftforce.sqr();                     // X SPLIT VS SPHERE CHANGE
                                           real distsq1 = driftforce[0]*driftforce[0];                                                          // X SPLIT VS SPHERE CHANGE
                                           if (distsq1 < min1sq) {
                                                min1sq = distsq1;
                                                mindriftforce = driftforce;
                                           }
                                    }
                                    min1sq = sqrt(min1sq);   // distance to nearest adress particle or center
                                    real mindriftforceX = (1.0/min1sq)*mindriftforce[0];  // normalized driftforce vector
                                    mindriftforceX *= -1.0*(clmassmultiplier-1.0)*0.5*vp.mass()*(1.0/sqrt(ntrotter))*omega2*xi*vp.lambdaDeriv();

                                    //mindriftforceX *= vp.lambdaDeriv();  // multiplication with derivative of the weighting function
                                    //vp.force() += mindriftforce;   // add drift force to virtual particles                                                                    // X SPLIT VS SPHERE CHANGE
                                    Real3D driftforceadd(mindriftforceX,0.0,0.0);                                                                                            // X SPLIT VS SPHERE CHANGE
                                    //Real3D driftforceadd(0.0,0.0,0.0);   
                                    //vp.force() += driftforceadd;             // Goes in, if one wants to apply the "normal" drift force - also improve using [0] ...           // X SPLIT VS SPHERE CHANGE


                                    at.forcem() -= driftforceadd;
                                }
                              }
                              else{
                                   std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
                                   exit(1);
                                   return;
                              }
                              
                          } 
                    }
                    else{
                         std::cout << "VP particle not found in tuples (VelocityVerletPI): " << vp.id() << "-" <<
                                 vp.ghost();
                         std::cout << " (" << vp.position() << ")\n";
                         exit(1);
                         return;      
                    }
                    
                }
                
      }
      
    }

    
    void VelocityVerletPI::updateForces(int f)
    {
      initForces();

      // signal
      aftInitF();  
        
      storage::Storage& storage = *getSystemRef().storage;

      if (f==2 || f==3){
        storage.updateGhosts();
      }
      recalc2();
      
      if (f==1){
        calcForcesF();
      }
      else if (f==2){
        calcForcesM();
      }
      else if (f==3){
        calcForcesS();
      }
      else {
        std::cout << "updateForces routine in VelocityVerletPI integrator received wrong integer.\n";
        exit(1);
        return;
      }    
      
      if (f==2 || f==3){ // Not necessary for internal ring forces
        storage.collectGhostForces();
      }
      
      if (f==2 || f==3){
        // signal 
        aftCalcFPI();
        transForces();
      }
    }

    void VelocityVerletPI::initForces()
    {
      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
        cit->forcem() = 0.0;
        cit->drift() = 0.0;   // Can in principle be commented, when drift is not used.
      }
    }
    
    
    real VelocityVerletPI::computeKineticEnergy()
    {
        real esum = 0.0;

        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        
        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) { 
            
            Particle &vp = *cit;
            
            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);
            if (it3 != fixedtupleList->end()) {
                std::vector<Particle*> atList;
                atList = it3->second;

                for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                    Particle &at = **it2;
                                        
                    if(at.pib() == 1){
                        esum += at.velocity().sqr() / vp.mass();
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        if((speedup == true) && (vp.lambda() < 0.000000001)){
                            continue;
                        }else{
                        	if(constkinmass == false){
                        		if(realkinmass == false){
                        			esum += at.velocity().sqr() / (vp.varmass()*CMDparameter*Eigenvalues.at(at.pib()-1));
                        		}
                        		else{
                        			esum += at.velocity().sqr() / (vp.varmass()*CMDparameter);
                        		}
                        	}else{
                        		if(realkinmass == false){
                        			esum += at.velocity().sqr() / (vp.mass()*CMDparameter*Eigenvalues.at(at.pib()-1));
                        		}
                        		else{
                        			esum += at.velocity().sqr() / (vp.mass()*CMDparameter);
                        		}
                        	}
                        }
                    }
                    else{
                         std::cout << "at.pib() outside of trotter range in computeKineticEnergy routine (VelocityVerletPI) \n";
                         exit(1);
                         return 0.0;
                    }
                    
                }
            
            }
            else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return 0.0;
            }
        }
        esum *= 0.5;

        real esumtotal;
        boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());

        return esumtotal;
    }


    real VelocityVerletPI::computeRingEnergy()
    {
        real esum = 0.0;

        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            Particle &vp = *cit;

            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);
            if (it3 != fixedtupleList->end()) {
                std::vector<Particle*> atList;
                atList = it3->second;

                for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                    Particle &at = **it2;

                    if(at.pib() == 1){
                        continue;
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        esum += at.modepos().sqr() * omega2 *  (vp.varmass()*Eigenvalues.at(at.pib()-1));
                    }
                    else{
                         std::cout << "at.pib() outside of trotter range in computeRingEnergy routine (VelocityVerletPI) \n";
                         exit(1);
                         return 0.0;
                    }

                }

            }
            else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return 0.0;
            }
        }
        esum *= 0.5;

        real esumtotal;
        boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());

        return esumtotal;
    }

    real VelocityVerletPI::computeRingEnergyRaw()
    {
        real esum = 0.0;

        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            Particle &vp = *cit;

            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);
            if (it3 != fixedtupleList->end()) {
                std::vector<Particle*> atList;
                atList = it3->second;

            	Real3D pos1(0.0,0.0,0.0);
            	Real3D posN(0.0,0.0,0.0);

                for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                    Particle &at = **it2;

                    if(at.pib() == 1){
                    	pos1 = at.position();
                    	posN = at.position();
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter-1){
                        esum += (at.position()-posN).sqr() * omega2 * vp.varmass();
                        posN = at.position();
                    }
                    else if(at.pib() == ntrotter){
                        esum += (at.position()-posN).sqr() * omega2 * vp.varmass();
                        esum += (at.position()-pos1).sqr() * omega2 * vp.varmass();
                    }
                    else{
                         std::cout << "at.pib() outside of trotter range in computeRingEnergy routine (VelocityVerletPI) \n";
                         exit(1);
                         return 0.0;
                    }

                }

            }
            else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return 0.0;
            }
        }
        esum *= 0.5;

        real esumtotal;
        boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());

        return esumtotal;
    }


    real VelocityVerletPI::computeMomentumDrift(int parttype)
    {
        real esum = 0.0;

        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            Particle &vp = *cit;

            if(vp.type() == parttype){

				FixedTupleListAdress::iterator it3;
				it3 = fixedtupleList->find(&vp);
				if (it3 != fixedtupleList->end()) {
					std::vector<Particle*> atList;
					atList = it3->second;

					for (std::vector<Particle*>::iterator it2 = atList.begin();
										   it2 != atList.end(); ++it2) {
						Particle &at = **it2;

						if(at.pib() == 1){
							continue;
						}
						else if(at.pib() > 1 && at.pib() <= ntrotter){
							if(realkinmass == false){
								esum += (clmassmultiplier-1.0)*0.5*vp.mass()*at.velocity().sqr()/(vp.varmass()*vp.varmass()*CMDparameter*Eigenvalues.at(at.pib()-1));
							}
							else{
								esum += (clmassmultiplier-1.0)*0.5*vp.mass()*at.velocity().sqr()/(vp.varmass()*vp.varmass()*CMDparameter);
							}
						}
						else{
							 std::cout << "at.pib() outside of trotter range in computeRingEnergy routine (VelocityVerletPI) \n";
							 exit(1);
							 return 0.0;
						}

					}

				}
				else { // this should not happen
					  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
					  std::cout << " (" << vp.position() << ")\n";
					  exit(1);
					  return 0.0;
				}

            }
        }

        real esumtotal;
        boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());

        return esumtotal;
    }


    real VelocityVerletPI::computePositionDrift(int parttype)
    {
        real esum = 0.0;

        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            Particle &vp = *cit;

            if(vp.type() == parttype){

				FixedTupleListAdress::iterator it3;
				it3 = fixedtupleList->find(&vp);
				if (it3 != fixedtupleList->end()) {
					std::vector<Particle*> atList;
					atList = it3->second;

					for (std::vector<Particle*>::iterator it2 = atList.begin();
										   it2 != atList.end(); ++it2) {
						Particle &at = **it2;

						if(at.pib() == 1){
							continue;
						}
						else if(at.pib() > 1 && at.pib() <= ntrotter){
							esum -= (clmassmultiplier-1.0)*0.5*vp.mass()*omega2*at.modepos().sqr()*Eigenvalues.at(at.pib()-1);
						}
						else{
							 std::cout << "at.pib() outside of trotter range in computeRingEnergy routine (VelocityVerletPI) \n";
							 exit(1);
							 return 0.0;
						}

					}

				}
				else { // this should not happen
					  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
					  std::cout << " (" << vp.position() << ")\n";
					  exit(1);
					  return 0.0;
				}

            }
        }

        real esumtotal;
        boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, esum, esumtotal, std::plus<real>());

        return esumtotal;
    }



    real VelocityVerletPI::weight(real distanceSqr){
        if (dex2 >= distanceSqr) return 1.0;
        else if (dexdhy2 <= distanceSqr) return 0.0;
        else {
            real argument = sqrt(distanceSqr) - dex;
            return pow(cos(pidhy2 * argument),2.0);
        }
    }

    real VelocityVerletPI::weightderivative(real distanceSqr){
        if (dex2 >= distanceSqr) return 0.0;
        else if (dexdhy2 <= distanceSqr) return 0.0;
        else{
            real argument = sqrt(distanceSqr) - dex;
            return -1.0 * pidhy2 * 2.0 * cos(pidhy2*argument) * sin(pidhy2*argument);
        }
    }


    void VelocityVerletPI::setTimeStep(real _dt)
    {
      dt = _dt;
      dt2 = dt*sStep;
      dt3 = dt*mStep*sStep;
    }

    void VelocityVerletPI::setmStep(int _mStep)
    {
      if (_mStep == 0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "mStep must be non-zero!";
        err.setException(msg.str());
      }
      mStep = _mStep;
    }

    void VelocityVerletPI::setsStep(int _sStep)
    {
      if (_sStep == 0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "lStep must be non-zero!";
        err.setException(msg.str());
      }
      sStep = _sStep;
    }

    void VelocityVerletPI::setNtrotter(int _ntrotter)
    {
      if (_ntrotter == 0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "ntrotter must be non-zero!";
        err.setException(msg.str());
      }
      ntrotter = _ntrotter;
    }

    void VelocityVerletPI::setSpeedup(bool _speedup)
    {
      speedup = _speedup;
    }

    void VelocityVerletPI::setKTI(bool _KTI)
    {
      KTI = _KTI;
    }

    void VelocityVerletPI::setPILE(bool _PILE)
    {
      PILE = _PILE;
    }

    void VelocityVerletPI::setRealKinMass(bool _realkinmass)
    {
    	realkinmass = _realkinmass;
    }

    void VelocityVerletPI::setCentroidThermostat(bool _centroidthermostat)
    {
    	centroidthermostat = _centroidthermostat;
    }

    void VelocityVerletPI::setConstKinMass(bool _constkinmass)
    {
    	constkinmass = _constkinmass;
    }

    void VelocityVerletPI::setTemperature(real _temperature)
    {
      if (_temperature < 0.0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "temperature must be larger or equal zero!";
        err.setException(msg.str());
      }
      real hbar = 1.054571726;
      real Nav = 6.02214129;
      real kb = 1.3806488;
      temperature = _temperature;
      omega2 = ntrotter*temperature*temperature*kb*kb*Nav*1.66053892/(hbar*hbar*0.00831451*0.00831451*1000.0);
    }

    void VelocityVerletPI::setGamma(real _gamma)
    {
      if (_gamma < 0.0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "gamma must be larger or equal zero!";
        err.setException(msg.str());
      }
      gamma = _gamma;
    }

    void VelocityVerletPI::setCMDparameter(real _CMDparameter)
    {
      if (_CMDparameter <= 0.0 || _CMDparameter > 1.0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "CMDparameter must be larger than zero and smaller or equal one!";
        err.setException(msg.str());
      }
      CMDparameter = _CMDparameter;
    }

    void VelocityVerletPI::setPILElambda(real _PILElambda)
    {
      if (_PILElambda < 0.0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "PILElambda must be larger or equal than zero!";
        err.setException(msg.str());
      }
      PILElambda = _PILElambda;
    }

    void VelocityVerletPI::setClmassmultiplier(real _clmassmultiplier)
    {
      if (_clmassmultiplier < 1.0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "clmassmultiplier must be larger or equal one!";
        err.setException(msg.str());
      }
      clmassmultiplier = _clmassmultiplier;
    }

    void VelocityVerletPI::setVerletList(shared_ptr<VerletListAdress> _verletList)
    {
      if (!_verletList) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "No Verletlist given in setVerletList";
        err.setException(msg.str());
      }
      verletList = _verletList;
    }

    void VelocityVerletPI::transp()
    {
        if (Tvectors.size() != ntrotter)
        {
            std::cout << "Number of Eigenvectors not equal to number of Trotter beads! \n";
            exit(1);
            return;
        }

        std::vector<real> tmpvec;

        for (size_t i = 0; i < ntrotter; ++i)
        {
            for (size_t j = 0; j < ntrotter; ++j)
            {
                tmpvec.push_back( Tvectors.at(j).at(i) );
            }

            Eigenvectors.push_back(tmpvec);
            tmpvec.clear();
        }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerletPI::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<VelocityVerletPI, bases<MDIntegrator>, boost::noncopyable >
        ("integrator_VelocityVerletPI", init< shared_ptr<System>, shared_ptr<VerletListAdress> >())
        .def("setTimeStep", &VelocityVerletPI::setTimeStep)
        .def("setmStep", &VelocityVerletPI::setmStep)
        .def("setsStep", &VelocityVerletPI::setsStep)
        .def("setNtrotter", &VelocityVerletPI::setNtrotter)
        .def("setTemperature", &VelocityVerletPI::setTemperature)
        .def("setGamma", &VelocityVerletPI::setGamma)
		.def("setCMDparameter", &VelocityVerletPI::setCMDparameter)
		.def("setPILElambda", &VelocityVerletPI::setPILElambda)
		.def("setClmassmultiplier", &VelocityVerletPI::setClmassmultiplier)
        .def("setSpeedup", &VelocityVerletPI::setSpeedup)
        .def("setKTI", &VelocityVerletPI::setKTI)
		.def("setCentroidThermostat", &VelocityVerletPI::setCentroidThermostat)
		.def("setPILE", &VelocityVerletPI::setPILE)
		.def("setRealKinMass", &VelocityVerletPI::setRealKinMass)
		.def("setConstKinMass", &VelocityVerletPI::setConstKinMass)
        .def("setVerletList", &VelocityVerletPI::setVerletList)
        //.def("setTimeStep", &VelocityVerletPI::setTimeStep)
        .def("add", &VelocityVerletPI::add)
        .def("addEV", &VelocityVerletPI::addEV)
        .def("addValues", &VelocityVerletPI::addValues)
        .def("transp", &VelocityVerletPI::transp)
        .def("computeKineticEnergy", &VelocityVerletPI::computeKineticEnergy)
        .def("computeRingEnergy", &VelocityVerletPI::computeRingEnergy)
		.def("computeRingEnergyRaw", &VelocityVerletPI::computeRingEnergyRaw)
        .def("computeMomentumDrift", &VelocityVerletPI::computeMomentumDrift)
        .def("computePositionDrift", &VelocityVerletPI::computePositionDrift)
        .add_property("mStep", &VelocityVerletPI::getmStep, &VelocityVerletPI::setmStep)
        .add_property("sStep", &VelocityVerletPI::getsStep, &VelocityVerletPI::setsStep)
        .add_property("ntrotter", &VelocityVerletPI::getNtrotter, &VelocityVerletPI::setNtrotter)
        .add_property("temperature", &VelocityVerletPI::getTemperature, &VelocityVerletPI::setTemperature)
        .add_property("gamma", &VelocityVerletPI::getGamma, &VelocityVerletPI::setGamma)
		.add_property("CMDparameter", &VelocityVerletPI::getCMDparameter, &VelocityVerletPI::setCMDparameter)
		.add_property("PILElambda", &VelocityVerletPI::getPILElambda, &VelocityVerletPI::setPILElambda)
        .add_property("speedup", &VelocityVerletPI::getSpeedup, &VelocityVerletPI::setSpeedup)
        .add_property("KTI", &VelocityVerletPI::getKTI, &VelocityVerletPI::setKTI)
        .add_property("centroidthermostat", &VelocityVerletPI::getCentroidThermostat, &VelocityVerletPI::setCentroidThermostat)
		.add_property("PILE", &VelocityVerletPI::getPILE, &VelocityVerletPI::setPILE)
		.add_property("realkinmass", &VelocityVerletPI::getRealKinMass, &VelocityVerletPI::setRealKinMass)
		.add_property("constkinmass", &VelocityVerletPI::getConstKinMass, &VelocityVerletPI::setConstKinMass)
        .add_property("verletList", &VelocityVerletPI::getVerletList, &VelocityVerletPI::setVerletList)
        ;
    }
  }
}
