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
      dt2 = sStep*dt;
      dt3 = mStep*sStep*dt;
              
      // Initialize variables
      temperature = 0.0;
      gamma = 0.0;
      ntrotter = 0;
      speedup = false;
      sStep = 0;
      mStep = 0;
      //verletList = NULL;
      
      //kb = 1.3806488;//*pow(10.0,-23); // FIX
      //real hbar = 1.054571726;//*pow(10.0,-34); // FIX
      //real Nav = 6.02214129;//*pow(10.0,23); // FIX
      //omega2 = ntrotter*temperature*temperature*kb*kb*Nav*1.66053892/(hbar*hbar*0.00831451*0.00831451);
      
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
      // Check if Eigenvalues start with zero:
      if(Eigenvalues.at(0) > 0.00000000001){
        std::cout << "Eigenvalues don't start with zero!\n";
        exit(1);
        return;
      }
      
          //DEBUG
          /*cout << "First Eigenvector: " << "( " << Eigenvectors.at(0).at(0) << ", "
                  << Eigenvectors.at(0).at(1) << ", " << Eigenvectors.at(0).at(2) << ", " << Eigenvectors.at(0).at(3) << " )\n";
          cout << "Second Eigenvector: " << "( " << Eigenvectors.at(1).at(0) << ", "
                  << Eigenvectors.at(1).at(1) << ", " << Eigenvectors.at(1).at(2) << ", " << Eigenvectors.at(1).at(3) << " )\n";
          cout << "Third Eigenvector: " << "( " << Eigenvectors.at(2).at(0) << ", "
                  << Eigenvectors.at(2).at(1) << ", " << Eigenvectors.at(2).at(2) << ", " << Eigenvectors.at(2).at(3) << " )\n";
          cout << "Fourth Eigenvector: " << "( " << Eigenvectors.at(3).at(0) << ", "
                  << Eigenvectors.at(3).at(1) << ", " << Eigenvectors.at(3).at(2) << ", " << Eigenvectors.at(3).at(3) << " )\n";
          
          cout << "First Tvector: " << "( " << Tvectors.at(0).at(0) << ", "
                  << Tvectors.at(0).at(1) << ", " << Tvectors.at(0).at(2) << ", " << Tvectors.at(0).at(3) << " )\n";
          cout << "Second Tvector: " << "( " << Tvectors.at(1).at(0) << ", "
                  << Tvectors.at(1).at(1) << ", " << Tvectors.at(1).at(2) << ", " << Tvectors.at(1).at(3) << " )\n";
          cout << "Third Tvector: " << "( " << Tvectors.at(2).at(0) << ", "
                  << Tvectors.at(2).at(1) << ", " << Tvectors.at(2).at(2) << ", " << Tvectors.at(2).at(3) << " )\n";
          cout << "Fourth Tvector: " << "( " << Tvectors.at(3).at(0) << ", "
                  << Tvectors.at(3).at(1) << ", " << Tvectors.at(3).at(2) << ", " << Tvectors.at(3).at(3) << " )\n";
          
          cout << "First Eigenvalue: " << Eigenvalues.at(0) << "\n";
          cout << "Second Eigenvalue: " << Eigenvalues.at(1) << "\n";
          cout << "Third Eigenvalue: " << Eigenvalues.at(2) << "\n";
          cout << "Fourth Eigenvalue: " << Eigenvalues.at(3) << "\n";*/
          //DEBUB
        
      //VT_TRACER("run");
      //int nResorts = 0;
      //real time;
      //timeIntegrate.reset();
      //resetTimers();
      System& system = getSystemRef();
      storage::Storage& storage = *system.storage;
      real skinHalf = 0.5 * system.getSkin();

      // signal
      runInit();

            //DEBUG
        /*//System& system = getSystemRef();
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
                
                    cout << "Position of AT particle with pib " << at.pib() << " of VP particle " << vp.id() << ":" << at.position() << "\n";
                    cout << "ModePosition of AT particle with pib " << at.pib() << " of VP particle " << vp.id() << ":" << at.modepos() << "\n";

                }
            
            }
        }
        cout << "1\n";*/
      //DEBUG
      
      // Before start make sure that particles are on the right processor
      //cout << "Resort Flag at start: " << resortFlag << "\n";
      if (resortFlag) {
        //cout << "DECOMPOSING at start!\n";
        //VT_TRACER("resort");
        // time = timeIntegrate.getElapsedTime();
        //LOG4ESPP_INFO(theLogger, "resort particles");
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
        // timeResort += timeIntegrate.getElapsedTime();
      }
      
            //DEBUG
        /*//System& system = getSystemRef();
        //CellList realCells = system.storage->getRealCells();
        //shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        
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
                
                    cout << "Position of AT particle with pib " << at.pib() << " of VP particle " << vp.id() << ":" << at.position() << "\n";
                    cout << "ModePosition of AT particle with pib " << at.pib() << " of VP particle " << vp.id() << ":" << at.modepos() << "\n";

                }
            
            }
        }
        cout << "2\n";*/
      //DEBUG
      
      transPos2(); // Update mode positions.

      //DEBUG
        /*//System& system = getSystemRef();
        //CellList realCells = system.storage->getRealCells();
        //shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        
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
                
                    cout << "Position of AT particle with pib " << at.pib() << " of VP particle " << vp.id() << ":" << at.position() << "\n";
                    cout << "ModePosition of AT particle with pib " << at.pib() << " of VP particle " << vp.id() << ":" << at.modepos() << "\n";

                }
            
            }
        }
        cout << "3\n";*/
      //DEBUG
      
      
      //bool recalcForces = true;  // TODO: more intelligent

      //if (recalcForces) {
        //LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

        // signal
        //recalc1();

        //updateForces(4);
        //if (LOG4ESPP_DEBUG_ON(theLogger)) {
            // printForces(false);   // forces are reduced to real particles
        //}

        // signal
        //recalc2();
      //}

      //LOG4ESPP_INFO(theLogger, "starting main integration loop (nsteps=" << nsteps << ")");
  
      // outer loop
      for (int i = 0; i < nsteps; i++) {
        //LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

        //saveOldPos(); // save particle positions needed for constraints
        //transPos2();
        
        
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
          //cout << "omega2: " << omega2 << "\n";
        
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
                //cout << "Resort Flag in loop: " << resortFlag << "\n";
                if (resortFlag) {
                    //cout << "DECOMPOSING in loop!\n";
                    storage.decompose();
                    //transPos2();
                    maxDist  = 0.0;
                    resortFlag = false;
                }               
                //transPos2();
                // CHECK SIGNALS BELOW!!!
               
                
                updateForces(1);
                integrateV2();
            
            }
            
            updateForces(2);
            integrateV1(2);
        
        }
        
        updateForces(3);
        integrateV1(3); 
        
        step++;
        
        
        // signal
        //befIntP();

        //time = timeIntegrate.getElapsedTime();
        //LOG4ESPP_INFO(theLogger, "updating positions and velocities")
        //maxDist += integrate1();
        //timeInt1 += timeIntegrate.getElapsedTime() - time;

        /*
        real cellsize = 1.4411685442;
        if (maxDist > 1.4411685442){
          cout<<"WARNING!!!!!! huge jump: "<<maxDist<<endl;
          exit(1);
        }*/
        
        // signal
        //aftIntP();

        //LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

        //transPos1(); // Update real positions in order to calculate forces afterwards and maybe trigger list update
        
        /*if (maxDist > skinHalf) resortFlag = true;
        
        if (resortFlag) {
            //VT_TRACER("resort1");
            //time = timeIntegrate.getElapsedTime();
            //LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
            storage.decompose();
            maxDist  = 0.0;
            resortFlag = false;
            //nResorts ++;
            //timeResort += timeIntegrate.getElapsedTime() - time;
        }*/

        //LOG4ESPP_INFO(theLogger, "updating forces")

        //updateForces(); FIX!!

        // signal
        //befIntV();

        //time = timeIntegrate.getElapsedTime();
        //integrate2();
        //timeInt2 += timeIntegrate.getElapsedTime() - time;

        // signal
        //aftIntV();
      }
        
        
        

              //DEBUG
        /*//System& system = getSystemRef();
        //CellList realCells = system.storage->getRealCells();
        //shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        
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
                
                    cout << "Position of AT particle with pib " << at.pib() << " of VP particle " << vp.id() << ":" << at.position() << "\n";
                    cout << "ModePosition of AT particle with pib " << at.pib() << " of VP particle " << vp.id() << ":" << at.modepos() << "\n";

                }
            
            }
        }
        cout << "4\n";*/
      //DEBUG
        
        
        
        
      /*timeRun = timeIntegrate.getElapsedTime();
      timeLost = timeRun - (timeForceComp[0] + timeForceComp[1] + timeForceComp[2] +
                 timeComm1 + timeComm2 + timeInt1 + timeInt2 + timeResort);

      LOG4ESPP_INFO(theLogger, "finished run");*/
    }

    /*void VelocityVerletPI::resetTimers() {
      timeForce  = 0.0;
      for(int i = 0; i < 100; i++)
        timeForceComp[i] = 0.0;
      timeComm1  = 0.0;
      timeComm2  = 0.0;
      timeInt1   = 0.0;
      timeInt2   = 0.0;
      timeResort = 0.0;
    }*/

    using namespace boost::python;

    /*static object wrapGetTimers(class VelocityVerletPI* obj) {
      real tms[10];
      obj->loadTimers(tms);
      return make_tuple(tms[0],
                        tms[1],
                        tms[2],
                        tms[3],
                        tms[4],
                        tms[5],
                        tms[6],
                        tms[7],
                        tms[8],
                        tms[9]);
    }*/

    /*void VelocityVerletPI::loadTimers(real t[10]) {
      t[0] = timeRun;
      t[1] = timeForceComp[0];
      t[2] = timeForceComp[1];
      t[3] = timeForceComp[2];
      t[4] = timeComm1;
      t[5] = timeComm2;
      t[6] = timeInt1;
      t[7] = timeInt2;
      t[8] = timeResort;
      t[9] = timeLost;
    }*/

    /*void VelocityVerletPI::printTimers() {

      using namespace std;
      real pct;

      cout << endl;
      cout << "run = " << setiosflags(ios::fixed) << setprecision(1) << timeRun << endl;
      pct = 100.0 * (timeForceComp[0] / timeRun);
      cout << "pair (%) = " << timeForceComp[0] << " (" << pct << ")" << endl;
      pct = 100.0 * (timeForceComp[1] / timeRun);
      cout << "FENE (%) = " << timeForceComp[1] << " (" << pct << ")" << endl;
      pct = 100.0 * (timeForceComp[2] / timeRun);
      cout << "angle (%) = " << timeForceComp[2] << " (" << pct << ")" << endl;
      pct = 100.0 * (timeComm1 / timeRun);
      cout << "comm1 (%) = " << timeComm1 << " (" << pct << ")" << endl;
      pct = 100.0 * (timeComm2 / timeRun);
      cout << "comm2 (%) = " << timeComm2 << " (" << pct << ")" << endl;
      pct = 100.0 * (timeInt1 / timeRun);
      cout << "int1 (%) = " << timeInt1 << " (" << pct << ")" << endl;
      pct = 100.0 * (timeInt2 / timeRun);
      cout << "int2 (%) = " << timeInt2 << " (" << pct << ")" << endl;
      pct = 100.0 * (timeResort / timeRun);
      cout << "resort (%) = " << timeResort << " (" << pct << ")" << endl;
      pct = 100.0 * (timeLost / timeRun);
      cout << "other (%) = " << timeLost << " (" << pct << ")" << endl;
      cout << endl;
    }*/

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
        //cout << "dt: " << dt << ", dt2: " << dt2 << ", dt3: " << dt3 <<"\n";
        //cout << "mStep: " << mStep << ", sStep: " << sStep <<"\n";
        //cout << "integrateV1(int t), t: " << t << " and half_df: " << half_dt << "\n";
                
        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        //cout << "IntegrateV1, half_dt: " << half_dt << " for t = " << t << "\n";
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
                        real dtfm = half_dt / vp.mass();
                        at.velocity() += dtfm * at.forcem();
                        //cout << "IntegrateV1, at.forcem() for pib " << at.pib() << ": " << at.forcem() << "\n";
                        //cout << "IntegrateV1, dtfm * at.forcem() for pib " << at.pib() << ": " << dtfm * at.forcem() << "\n";
                        vp.velocity() = at.velocity();
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        //real dtfm = half_dt / (vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter)));
                        real dtfm = half_dt / (vp.varmass()*Eigenvalues.at(at.pib()-1));        
                        at.velocity() += dtfm * at.forcem();
                        //cout << "IntegrateV1, at.forcem() for pib " << at.pib() << ": " << at.forcem() << "\n";
                        //cout << "IntegrateV1, dtfm * at.forcem() for pib " << at.pib() << ": " << dtfm * at.forcem() << "\n";
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
    
    
    void VelocityVerletPI::integrateV2()
    {     
        real half_dt = 0.5 * dt;
        real half_dt4 = 0.25 * dt;
                
        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        //cout << "IntegrateV2, half_dt: " << half_dt << " and half_dt4: " << half_dt4 << "\n";
        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) { 
            
            Particle &vp = *cit;
            
            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);
            if (it3 != fixedtupleList->end()) {
                std::vector<Particle*> atList;
                atList = it3->second;

                // First half_dt4 loop without centroid mode
                for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                    Particle &at = **it2;
                
                    if(at.pib() == 1){
                        continue;
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        //real dtfm = half_dt4 / (vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter)));
                        real dtfm = half_dt4 / (vp.varmass()*Eigenvalues.at(at.pib()-1));        
                        at.velocity() += dtfm * at.forcem();
                        //cout << "IntegrateV2 - 1, at.forcem() for pib " << at.pib() << ": " << at.forcem() << "\n";
                        //cout << "IntegrateV2 - 1, dtfm * at.forcem() for pib " << at.pib() << ": " << dtfm * at.forcem() << "\n";
                    }
                    else{
                         std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
                         exit(1);
                         return;
                    }
                }
                
                // half_dt2 centroid mode
                for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                    Particle &at = **it2;
                
                    if(at.pib() == 1){
                        
                        // Get the inner factor
                        real xi = 0.0;
                        for (std::vector<Particle*>::iterator it5 = atList.begin();
                                     it5 != atList.end(); ++it5) {
                            Particle &at2 = **it5;
                            if(at2.pib() != 1){
                                //real chi = 2.0*(1.0-cos(2.0*M_PI*(at2.pib()-1.0)/ntrotter))*at2.velocity().sqr();
                                real chi = Eigenvalues.at(at2.pib()-1)*at2.velocity().sqr();        
                                xi += chi;
                            }
                        }
                        
                        real dtfm = half_dt / (vp.mass());
                        
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
                        //mindriftforce *= weightderivative(min1sq);  // multiplication with derivative of the weighting function
                        //mindriftforceX *= 0.5;
                        mindriftforceX *= 49.5*xi*vp.lambdaDeriv();   // get the energy differences which were calculated previously and put in drift force

                        //vp.drift() += mindriftforceX ;//* vp.lambdaDeriv();    // USE ONLY LIKE THAT, IF DOING ITERATIVE FEC INCLUDING ITERATIVE PRESSURE FEC               

                        //mindriftforceX *= vp.lambdaDeriv();  // multiplication with derivative of the weighting function
                        //vp.force() += mindriftforce;   // add drift force to virtual particles                                                                    // X SPLIT VS SPHERE CHANGE
                        Real3D driftforceadd(mindriftforceX,0.0,0.0);                                                                                            // X SPLIT VS SPHERE CHANGE
                        //Real3D driftforceadd(0.0,0.0,0.0);   
                        //vp.force() += driftforceadd;             // Goes in, if one wants to apply the "normal" drift force - also improve using [0] ...           // X SPLIT VS SPHERE CHANGE
                        //std::cout << "Added Drift Force: " << driftforceadd << " for particle at pos(x).: " << vp.position()[0] << "\n";                       
                        at.velocity() += dtfm * (at.forcem() - driftforceadd);
                        
                        //cout << "IntegrateV2 - 2, at.forcem() - driftforceadd for pib " << at.pib() << ": " << at.forcem() - driftforceadd << "\n";
                        
                        vp.velocity() = at.velocity();
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
                for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                    Particle &at = **it2;
                
                    if(at.pib() == 1){
                        continue;
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        //real dtfm = half_dt4 / (vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter)));
                        real dtfm = half_dt4 / (vp.varmass()*Eigenvalues.at(at.pib()-1));
                        at.velocity() += dtfm * at.forcem();
                        //cout << "IntegrateV2 - 3, at.forcem() for pib " << at.pib() << ": " << at.forcem() << "\n";
                        //cout << "IntegrateV2 - 3, dtfm * at.forcem() for pib " << at.pib() << ": " << dtfm * at.forcem() << "\n";
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
    
    
    void VelocityVerletPI::integrateModePos()
    {   
        real half_dt = 0.5 * dt; 
                
        System& system = getSystemRef();
        CellList realCells = system.storage->getRealCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
        //cout << "IntegrateModePos, half_dt: " << half_dt << "\n";
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
                    
                    if((speedup == true) && (vp.lambda() < 0.000000001) && (at.pib() != 1)){
                        continue;
                    }else{
                        at.modepos() += half_dt * at.velocity();
                        //cout << "IntegrateModePos, half_dt * at.velocity() for pib " << at.pib() << ": " << half_dt * at.velocity() << "\n";
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
                        Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);
                        at.velocity() = prefac1 * at.velocity() + prefac2/sqrt(vp.mass()) * ranval;
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);
                        //at.velocity() = prefac1 * at.velocity() + (prefac2/sqrt((vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter))))) * ranval;
                        at.velocity() = prefac1 * at.velocity() + (prefac2/sqrt((vp.varmass()*Eigenvalues.at(at.pib()-1)))) * ranval;
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
    
    
    /*real VelocityVerletPI::integrate1()
    {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells
      int count = 0;
      real maxSqDist = 0.0; // maximal square distance a particle moves
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real sqDist = 0.0;
        LOG4ESPP_INFO(theLogger, "updating first half step of velocities and full step of positions")
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << 
                ", pos = " << cit->position() <<
                ", v = " << cit->velocity() << 
                ", f = " << cit->force());

        /* more precise for DEBUG:
        printf("Particle %d, pos = %16.12f %16.12f %16.12f, v = %16.12f, %16.12f %16.12f, f = %16.12f %16.12f %16.12f\n",
            cit->p.id, cit->r.p[0], cit->r.p[1], cit->r.p[2],
                cit->m.v[0], cit->m.v[1], cit->m.v[2],
            cit->f.f[0], cit->f.f[1], cit->f.f[2]);
        */

        /*real dtfm = 0.5 * dt / cit->mass();

        // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) 
        cit->velocity() += dtfm * cit->force();

        // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt) 
        Real3D deltaP = cit->velocity();
        
        deltaP *= dt;
        cit->position() += deltaP;
        sqDist += deltaP * deltaP;

        count++;

        maxSqDist = std::max(maxSqDist, sqDist);
      }
      
      // signal
      //inIntP(maxSqDist);

      real maxAllSqDist;
      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());

      LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
		    ", max move local = " << sqrt(maxSqDist) <<
		    ", global = " << sqrt(maxAllSqDist));
      
      return sqrt(maxAllSqDist);
    }*/

    /*void VelocityVerletPI::integrate2()
    {
      LOG4ESPP_INFO(theLogger, "updating second half step of velocities")
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells
      real half_dt = 0.5 * dt; 
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real dtfm = half_dt / cit->mass();
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        /*cit->velocity() += dtfm * cit->force();
      }
      
      step++;
    }*/

    /*void VelocityVerletPI::calcForces()
    {
      VT_TRACER("forces");

      LOG4ESPP_INFO(theLogger, "calculate forces");

      initForces();

      // signal
      aftInitF();

      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;

      for (size_t i = 0; i < srIL.size(); i++) {
	    LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
        real time;
        time = timeIntegrate.getElapsedTime();
        srIL[i]->addForces();
        timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
      }
    }*/
 
    
    void VelocityVerletPI::transPos1(){
        // Update real positions based on mode positions
        
        real maxSqDist = 0.0;
        
        System& system = getSystemRef();
        CellList localCells = system.storage->getLocalCells();
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
                      
                      Real3D oldpos = vp.position();
                      //cout << "oldpos: " << oldpos << "\n";
                      Real3D zero(0.0,0.0,0.0);
                      at.position()=zero;

                      for (std::vector<Particle*>::iterator it5 = atList.begin();
                                   it5 != atList.end(); ++it5) {
                          Particle &at2 = **it5;

                          if(at.pib() <= ntrotter){
                              at.position()+= sqrt(ntrotter)*at2.modepos()*Tvectors.at(at.pib()-1).at(at2.pib()-1);
                              
                              // Update VP particle (this thing is executed too often... but this shouldn't matter, might be fixed later)
                              if(at2.pib() == 1){
                                  //vp.position() = at.position();
                                  //vp.velocity() = at.velocity();
                                  vp.position() = at2.modepos();
                                  vp.velocity() = at2.velocity();
                              }
                          }
                          else{
                               std::cout << "at2.pib() outside of trotter range in transPos1 routine (VelocityVerletPI) \n";
                               exit(1);
                               return;
                          }

                      }
                      
                      //cout << "transPos1(), at.position() for pib: " << at.pib() << ", ghost: " << vp.ghost() << ", position: " << at.position() <<  ", modepos: " << at.modepos() << "\n";
                      
                      //cout << "newpos: " << vp.position() << "\n";
                      
                      real sqDist = (oldpos-vp.position())*(oldpos-vp.position());
                      maxSqDist = std::max(maxSqDist, sqDist);
                      //cout << "maxSqDist: " << maxSqDist << "\n"; 

                  }
                  
                  /*//cout << "transPos1(), at.position() for pib: " << at.pib() << ", ghost: " << vp.ghost() << ", position: " << at.position() <<  ", modepos: " << at.modepos() << "\n";

                  cout << "newpos: " << vp.position() << "\n";

                  real sqDist = (oldpos-vp.position())*(oldpos-vp.position());
                  maxSqDist = std::max(maxSqDist, sqDist);
                  cout << "maxSqDist: " << maxSqDist << "\n";*/ 
                      
            }
            else{
                 std::cout << "VP particle not found in tuples (VelocityVerletPI): " << vp.id() << "-" <<
                         vp.ghost();
                 std::cout << " (" << vp.position() << ")\n";
                 exit(1);
                 return;      
            }
                
        }
      
        // signal
        inIntP(maxSqDist);

        real maxAllSqDist;
        mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());
        //cout << "maxAllSqDist: " << maxAllSqDist << "\n";
        maxDist+=sqrt(maxAllSqDist);
        
    }
    
    
    void VelocityVerletPI::transPos2(){
        // Update mode positions based on real positions
        
        System& system = getSystemRef();
        CellList localCells = system.storage->getLocalCells();
        shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();

        for(CellListIterator cit(localCells); !cit.isDone(); ++cit){
            Particle &vp = *cit;
            FixedTupleListAdress::iterator it2;
            it2 = fixedtupleList->find(&vp);

            if(it2 != fixedtupleList->end()){
                  std::vector<Particle*> atList;
                  atList = it2->second;
                  int pib_check = 0;
                  for (std::vector<Particle*>::iterator it3 = atList.begin();
                                       it3 != atList.end(); ++it3) {
                      Particle &at = **it3;
                      Real3D zero(0.0,0.0,0.0);
                      at.modepos()=zero;
                      
                      for (std::vector<Particle*>::iterator it5 = atList.begin();
                                   it5 != atList.end(); ++it5) {
                          Particle &at2 = **it5;

                          if(at.pib() <= ntrotter){
                              at.modepos()+= sqrt(1.0/ntrotter)*at2.position()*Eigenvectors.at(at.pib()-1).at(at2.pib()-1);                             
                          }
                          else{
                               std::cout << "at2.pib() outside of trotter range in transPos2 routine (VelocityVerletPI) \n";
                               exit(1);
                               return;
                          }

                    }        
                    //cout << "transPos2(), at.position() for pib: " << at.pib() << ", ghost: " << vp.ghost() << ", position: " << at.position() <<  ", modepos: " << at.modepos() << ", velocity: " << at.velocity() << "\n";
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
        // Update mode forces based on real forces
        
        System& system = getSystemRef();
        CellList localCells = system.storage->getLocalCells();
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
                              at.forcem()+= at2.force(); // + (1.0/ntrotter)*vp.force();   // Is the (1.0/ntrotter)*vp.force() correct?
                          }
                          else if(at.pib() > 1 && at.pib() <= ntrotter){
                              at.forcem()+= sqrt(ntrotter)*at2.force()*Eigenvectors.at(at.pib()-1).at(at2.pib()-1);
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
      //VT_TRACER("forces");

      //LOG4ESPP_INFO(theLogger, "calculate forces");

      //initForces();

      // signal
      //aftInitF();

      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;

      for (size_t i = 0; i < srIL.size(); i++) {
	    //LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
        if(srIL[i]->bondType() == Nonbonded){
            //real time;
            //time = timeIntegrate.getElapsedTime();
            srIL[i]->addForces();
            //timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
        }
      }
      
      // Signal // For FEC
      aftCalcFslow();       
    }
    
    
    void VelocityVerletPI::calcForcesM()
    {
      //VT_TRACER("forces");

      //LOG4ESPP_INFO(theLogger, "calculate forces");

      //initForces();

      // signal
      //aftInitF();

      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;

      for (size_t i = 0; i < srIL.size(); i++) {
	    LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
        if(srIL[i]->bondType() == Pair || srIL[i]->bondType() == Angular){
            //real time;
            //time = timeIntegrate.getElapsedTime();
            srIL[i]->addForces();
            //timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
        }
      }      
    }
        
    
    void VelocityVerletPI::calcForcesF()
    {
      //VT_TRACER("forces");

      //LOG4ESPP_INFO(theLogger, "calculate forces");

      //initForces();

      // signal
      //aftInitF();

      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();
      shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
      //for (size_t i = 0; i < srIL.size(); i++) {
	//    LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
      //real time;
      //time = timeIntegrate.getElapsedTime();

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
                                    //at.forcem() -= at.modepos()*omega2*vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter));
                                    at.forcem() -= at.modepos()*omega2*vp.varmass()*Eigenvalues.at(at.pib()-1);
                                    cout << "calcForcesF, at.forcem(): " << at.forcem() << " for at.modepos(): " << at.modepos() << "\n";
                                  }
                                  else if(at.pib() == 1){
                                    if(vp.lambda()<0.999999999 && vp.lambda()>0.000000001){                          
                                        real xi = 0.0;
                                        for (std::vector<Particle*>::iterator it5 = atList.begin();
                                                     it5 != atList.end(); ++it5) {
                                            Particle &at2 = **it5;
                                            if(at2.pib() != 1){
                                                //real chi = 2.0*(1.0-cos(2.0*M_PI*(at2.pib()-1.0)/ntrotter));   
                                                //xi += -49.5*omega2*at2.modepos().sqr()*chi;
                                                xi += -49.5*omega2*at2.modepos().sqr()*Eigenvalues.at(at2.pib()-1);
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
                                        //mindriftforce *= weightderivative(min1sq);  // multiplication with derivative of the weighting function
                                        //mindriftforceX *= 0.5;
                                        mindriftforceX *= xi*vp.lambdaDeriv();   // get the energy differences which were calculated previously and put in drift force

                                        //vp.drift() += mindriftforceX ;//* vp.lambdaDeriv();    // USE ONLY LIKE THAT, IF DOING ITERATIVE FEC INCLUDING ITERATIVE PRESSURE FEC               

                                        //mindriftforceX *= vp.lambdaDeriv();  // multiplication with derivative of the weighting function
                                        //vp.force() += mindriftforce;   // add drift force to virtual particles                                                                    // X SPLIT VS SPHERE CHANGE
                                        Real3D driftforceadd(mindriftforceX,0.0,0.0);                                                                                            // X SPLIT VS SPHERE CHANGE
                                        //Real3D driftforceadd(0.0,0.0,0.0);   
                                        //vp.force() += driftforceadd;             // Goes in, if one wants to apply the "normal" drift force - also improve using [0] ...           // X SPLIT VS SPHERE CHANGE
                                        //std::cout << "Added Drift Force: " << driftforceadd << " for particle at pos(x).: " << vp.position()[0] << "\n";


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
                                //at.forcem() -= at.modepos()*omega2*vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter));
                                at.forcem() -= at.modepos()*omega2*vp.varmass()*Eigenvalues.at(at.pib()-1);
                                //std::cout << "Eigenvalue for " << at.pib() << " is: Eigenvalues.at(at.pib()-1): " << Eigenvalues.at(at.pib()-1) << "\n";
                              }
                              else if(at.pib() == 1){
                                if(vp.lambda()<0.999999999 && vp.lambda()>0.000000001){                          
                                    real xi = 0.0;
                                    for (std::vector<Particle*>::iterator it5 = atList.begin();
                                                 it5 != atList.end(); ++it5) {
                                        Particle &at2 = **it5;
                                        if(at2.pib() != 1){
                                            //real chi = 2.0*(1.0-cos(2.0*M_PI*(at2.pib()-1.0)/ntrotter));
                                            //xi += -49.5*omega2*at2.modepos().sqr()*chi;
                                            xi += -49.5*omega2*at2.modepos().sqr()*Eigenvalues.at(at2.pib()-1);
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
                                    //mindriftforce *= weightderivative(min1sq);  // multiplication with derivative of the weighting function
                                    //mindriftforceX *= 0.5;
                                    mindriftforceX *= xi*vp.lambdaDeriv();   // get the energy differences which were calculated previously and put in drift force

                                    //vp.drift() += mindriftforceX ;//* vp.lambdaDeriv();    // USE ONLY LIKE THAT, IF DOING ITERATIVE FEC INCLUDING ITERATIVE PRESSURE FEC               

                                    //mindriftforceX *= vp.lambdaDeriv();  // multiplication with derivative of the weighting function
                                    //vp.force() += mindriftforce;   // add drift force to virtual particles                                                                    // X SPLIT VS SPHERE CHANGE
                                    Real3D driftforceadd(mindriftforceX,0.0,0.0);                                                                                            // X SPLIT VS SPHERE CHANGE
                                    //Real3D driftforceadd(0.0,0.0,0.0);   
                                    //vp.force() += driftforceadd;             // Goes in, if one wants to apply the "normal" drift force - also improve using [0] ...           // X SPLIT VS SPHERE CHANGE
                                    //std::cout << "Added Drift Force: " << driftforceadd << " for particle at pos(x).: " << vp.position()[0] << "\n";


                                    at.forcem() -= driftforceadd;
                                }
                              }
                              else{
                                   std::cout << "at.pib() outside of trotter range in TransForce routine (VelocityVerletPI) \n";
                                   exit(1);
                                   return;
                              }
                              
                              //std::cout << "FastForce for pib " << at.pib() << " is: at.forcem(): " << at.forcem() << "\n"; 
                              
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
      
      //timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
      //}
    }

    
    void VelocityVerletPI::updateForces(int f)
    {
      initForces();

      // signal
      aftInitF();  
        
      //LOG4ESPP_INFO(theLogger, "update ghosts, calculate forces and collect ghost forces")
      //real time;
      storage::Storage& storage = *getSystemRef().storage;
      //time = timeIntegrate.getElapsedTime();
      { 
        //VT_TRACER("commF");
        storage.updateGhosts();
        transPos2();
      }
      //timeComm1 += timeIntegrate.getElapsedTime() - time;
      //time = timeIntegrate.getElapsedTime();
      
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
      
      //calcForces();
      //timeForce += timeIntegrate.getElapsedTime() - time;
      //time = timeIntegrate.getElapsedTime();
      {
        //VT_TRACER("commR");
        storage.collectGhostForces();
      }
      //timeComm2 += timeIntegrate.getElapsedTime() - time;

      // signal
      //aftCalcF();
      
      if (f==2 || f==3){
        // signal 
        aftCalcF();
        transForces();
      }
    }

    void VelocityVerletPI::initForces()
    {
      // forces are initialized for real + ghost particles

      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        //cout << "\n";
        //cout << "cit->force(): " << cit->force() << "\n";
        //cout << "cit->forcem(): " << cit->forcem() << "\n";
        cit->force() = 0.0;
        cit->forcem() = 0.0;   // Also initialize mode force to zero
        cit->drift() = 0.0;   // Can in principle be commented, when drift is not used.
        //cout << "cit->force(): " << cit->force() << "\n";
        //cout << "cit->forcem(): " << cit->forcem() << "\n";
        //cout << "\n";
      }
    }

    
    /*void VelocityVerletPI::printForces(bool withGhosts)
    {
      // print forces of real + ghost particles

      System& system = getSystemRef();
      CellList cells;

      if (withGhosts) {
	    cells = system.storage->getLocalCells();
	    LOG4ESPP_DEBUG(theLogger, "local forces");
      } else {
	    cells = system.storage->getRealCells();
	    LOG4ESPP_DEBUG(theLogger, "real forces");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
	    LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", force = " << cit->force());
      }
    }*/

    /*void VelocityVerletPI::printPositions(bool withGhosts)
    {
      // print positions of real + ghost particles

      System& system = getSystemRef();
      CellList cells;

      if (withGhosts) {
	    cells = system.storage->getLocalCells();
	    LOG4ESPP_DEBUG(theLogger, "local positions");
      } else {
	    cells = system.storage->getRealCells();
	    LOG4ESPP_DEBUG(theLogger, "real positions");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
	    LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", position = " << cit->position());
      }
    }*/
    
    
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
                        esum += at.velocity().sqr() * vp.mass();
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        if((speedup == true) && (vp.lambda() < 0.000000001)){
                            continue;
                        }else{
                            //esum += at.velocity().sqr() * (vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter)));
                            esum += at.velocity().sqr() * (vp.varmass()*Eigenvalues.at(at.pib()-1));        
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
                    
                    //esum += at.velocity().sqr() * (vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter)));
                    //esum += at.velocity().sqr() * (vp.varmass()*Eigenvalues.at(at.pib()-1));        
                    
                    if(at.pib() == 1){
                        continue;
                    }
                    else if(at.pib() > 1 && at.pib() <= ntrotter){
                        //esum += at.modepos().sqr() * omega2 * (vp.varmass()*2.0*(1.0-cos(2.0*M_PI*(at.pib()-1.0)/ntrotter)));
                        esum += at.modepos().sqr() * omega2 * (vp.varmass()*Eigenvalues.at(at.pib()-1));
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
    
    void VelocityVerletPI::setTimeStep(real _dt)
    {
      /*if (_mStep == 0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "mStep must be non-zero!";
        err.setException(msg.str());
      }*/

      dt = _dt;
      dt2 = sStep*dt;
      dt3 = mStep*sStep*dt;
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
      /*if (_ntrotter == 0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "ntrotter must be non-zero!";
        err.setException(msg.str());
      }*/

      speedup = _speedup;
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
      real hbar = 1.054571726;//*pow(10.0,-34); // FIX
      real Nav = 6.02214129;//*pow(10.0,23); // FIX
      real kb = 1.3806488;
      temperature = _temperature;
      //omega2 = ntrotter*temperature*temperature*kb*kb*Nav*1.66053892*pow(10.,-48)/(hbar*hbar*0.00831451*0.00831451);
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
        //cout << "YAY!!!\n";
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
    
    
    /*void VelocityVerletPI::addEV(tuple vals) {
        //bool returnVal = true;
        
        Eigenvectors.push_back(vals);
        
        System& system = storage->getSystemRef();
        esutil::Error err(system.comm);

        // ADD THE LOCAL PARTICLES (pointers)
        Particle* vp, *at;
        std::vector<Particle*> tmp; // temporary vector
        std::vector<longint> pidstmp; // temporary vector
        longint pidK; // the pid used as key

        tuple::iterator it = pids.begin();
        vp = storage->lookupRealParticle(*it);
        if (!vp) { // Particle does not exist here, return false
            //std::cout << "particle " << *it << " not found in localParticles \n";
            returnVal = false;
        }
        else{
          pidK = *it; // first pid is key
          //std::cout << "Add key: " << *it << "\n";

          for (++it; it != pids.end(); ++it) {

              at = storage->lookupAdrATParticle(*it);
              if (!at) { // Particle does not exist here, return false
                  std::stringstream msg;
                  msg << "ERROR: AT particle " << *it << " not found in localAdrATParticles \n";
                  err.setException( msg.str() );
                  returnVal = false;
                  break;
              }
              tmp.push_back(at);
              //std::cout << " add: " << *it << "\n";
              pidstmp.push_back(*it); // pidK is not in this vector
          }
        }
        err.checkException();
        
        if(returnVal){
            this->add(vp, tmp); // add to TupleList

            // ADD THE GLOBAL PARTICLES (ids)
            globalTuples.insert(make_pair(pidK, pidstmp));
        }
        LOG4ESPP_INFO(theLogger, "added fixed tuple to global tuples");

        tmp.clear();
        pids.clear();
        pidstmp.clear();

        //std::cout << "\n";

        return returnVal;
    }*/
    

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerletPI::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<VelocityVerletPI, bases<MDIntegrator>, boost::noncopyable >
        ("integrator_VelocityVerletPI", init< shared_ptr<System>, shared_ptr<VerletListAdress> >())
        //.def("getTimers", &wrapGetTimers)
        //.def("resetTimers", &VelocityVerletPI::resetTimers)
        .def("setTimeStep", &VelocityVerletPI::setTimeStep)
        .def("setmStep", &VelocityVerletPI::setmStep)
        .def("setsStep", &VelocityVerletPI::setsStep)
        .def("setNtrotter", &VelocityVerletPI::setNtrotter)
        .def("setTemperature", &VelocityVerletPI::setTemperature)
        .def("setGamma", &VelocityVerletPI::setGamma)
        .def("setSpeedup", &VelocityVerletPI::setSpeedup)
        //.def("addEigenvectors", &VelocityVerletPI::setNtrotter)
        .def("setVerletList", &VelocityVerletPI::setVerletList)
        .def("setTimeStep", &VelocityVerletPI::setTimeStep)
        .def("add", &VelocityVerletPI::add)
        .def("addEV", &VelocityVerletPI::addEV)
        .def("addValues", &VelocityVerletPI::addValues)
        .def("transp", &VelocityVerletPI::transp)
        .def("computeKineticEnergy", &VelocityVerletPI::computeKineticEnergy)
        .def("computeRingEnergy", &VelocityVerletPI::computeRingEnergy)
        .add_property("mStep", &VelocityVerletPI::getmStep, &VelocityVerletPI::setmStep)
        .add_property("sStep", &VelocityVerletPI::getsStep, &VelocityVerletPI::setsStep)
        .add_property("ntrotter", &VelocityVerletPI::getNtrotter, &VelocityVerletPI::setNtrotter)
        .add_property("temperature", &VelocityVerletPI::getTemperature, &VelocityVerletPI::setTemperature)
        .add_property("gamma", &VelocityVerletPI::getGamma, &VelocityVerletPI::setGamma)
        .add_property("speedup", &VelocityVerletPI::getSpeedup, &VelocityVerletPI::setSpeedup)
        .add_property("verletList", &VelocityVerletPI::getVerletList, &VelocityVerletPI::setVerletList)
        ;
    }
  }
}
