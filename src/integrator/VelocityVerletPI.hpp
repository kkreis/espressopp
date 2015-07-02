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

// ESPP_CLASS
#ifndef _INTEGRATOR_VELOCITYVERLETPI_HPP
#define _INTEGRATOR_VELOCITYVERLETPI_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"
//#include "esutil/Timer.hpp"
#include <boost/signals2.hpp>
#include "VerletListAdress.hpp"

namespace espressopp {
  namespace integrator {

    /** Velocity Verlet Integrator */
    class VelocityVerletPI : public MDIntegrator {

      public:
        //shared_ptr<VerletListAdress> verletList;
        
        VelocityVerletPI(shared_ptr<class espressopp::System> system, shared_ptr<VerletListAdress> _verletList);

        virtual ~VelocityVerletPI();

        shared_ptr<VerletListAdress> verletList;        
       
        /* TODO should be removed after signals will be tested
        void setLangevin(shared_ptr<class Langevin> langevin);
        shared_ptr<class Langevin> getLangevin() { return langevin; }

        void setStochasticVelocityRescaling(shared_ptr<class StochasticVelocityRescaling> stochasticVelocityRescaling);
        shared_ptr<class StochasticVelocityRescaling> getStochasticVelocityRescaling() { return stochasticVelocityRescaling; }

        void setIsokinetic(shared_ptr<class Isokinetic> isokinetic);
        shared_ptr<class Isokinetic> getIsokinetic() { return isokinetic; }
        
        // set & get barostat (Berendsen)
        void setBerendsenBarostat(shared_ptr<class BerendsenBarostat> berendsenBarostat);
        shared_ptr<class BerendsenBarostat> getBerendsenBarostat() { return berendsenBarostat; }
        
        // set & get thermostat (Berendsen)
        void setBerendsenThermostat(shared_ptr<class BerendsenThermostat> berendsenThermostat);
        shared_ptr<class BerendsenThermostat> getBerendsenThermostat() { return berendsenThermostat; }
        
        // set & get barostat (Langevin-Hoover)
        void setLangevinBarostat(shared_ptr<class LangevinBarostat> langevinBarostat);
        shared_ptr<class LangevinBarostat> getLangevinBarostat() { return langevinBarostat; }
        
        // set & get FixPositions
        void setFixPositions(shared_ptr <class FixPositions> _fixPositions);
        shared_ptr<class FixPositions> getFixPositions () {return fixPositions; }
        */

        void run(int nsteps);
        
        /** Load timings in array to export to Python as a tuple. */
        //void loadTimers(real t[10]);

        //void resetTimers();

        /** Setter routine for the fast TimeStep. */
        void setTimeStep(real _dt);
        /** Getter routine for the fast TimeStep. */
        real getTimeStep() { return dt; }
        
        /** Setter routine for the mStep. */
        void setmStep(int _mStep);
        /** Getter routine for the mStep. */
        real getmStep() { return mStep; }
        
        /** Setter routine for the sStep. */
        void setsStep(int _sStep);
        /** Getter routine for the sStep. */
        real getsStep() { return sStep; }
        
        /** Setter routine for ntrotter. */
        void setNtrotter(int _ntrotter);
        /** Getter routine for ntrotter. */
        real getNtrotter() { return ntrotter; }

        void add(real val) { tmpvals.push_back(val); } // add particle id (called from python)
        //void addEVs() { addEV(tmpvals); tmpvals.clear(); } // add tuple (called from python)
        void addEV() { Eigenvectors.push_back(tmpvals); tmpvals.clear();}
        
        void transp();
        
        /** Setter routine for temperature. */
        void setTemperature(real _temperature);
        /** Getter routine for temperature. */
        real getTemperature() { return temperature; }
        
        /** Setter routine for gamma. */
        void setGamma(real _gamma);
        /** Getter routine for gamma. */
        real getGamma() { return gamma; }
        
        /** Setter routine for speedup. */
        void setSpeedup(bool _speedup);
        /** Getter routine for speedup. */
        bool getSpeedup() { return speedup; }        
        
        /** Setter routine for verletList. */
        void setVerletList(shared_ptr<VerletListAdress> _verletList);
        /** Getter routine for verletList. */ 
        shared_ptr<VerletListAdress> getVerletList() { return verletList; }
        
        // Compute internal ring energies
        real computeRingEnergy();
        
        // Compute special kinetic energy
        real computeKineticEnergy();
        
        
        // signal used for constraints
        //boost::signals2::signal0 <void> saveOldPos;
		//boost::signals2::signal0 <void> applyPosConst;
		//boost::signals2::signal0 <void> applyVelConst;


        /* -- moved to superclass
        // signals to extend the integrator
        boost::signals2::signal0 <void> runInit; // initialization of run()
        boost::signals2::signal0 <void> recalc1; // inside recalc, before updateForces()
        boost::signals2::signal0 <void> recalc2; // inside recalc, after  updateForces()
        boost::signals2::signal0 <void> befIntP; // before integrate1()
        boost::signals2::signal1 <void, real&> inIntP; // inside end of integrate1()
        boost::signals2::signal0 <void> aftIntP; // after  integrate1()
        boost::signals2::signal0 <void> aftInitF; // after initForces()
        boost::signals2::signal0 <void> aftCalcF; // after calcForces()
        boost::signals2::signal0 <void> befIntV; // before integrate2()
        boost::signals2::signal0 <void> aftIntV; // after  integrate2()
        */

        //System& getSystem();


        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:

        int sStep; // medium time step = sStep * dt
        int mStep; // large time step = mStep * sStep * dt
        int ntrotter; // number of Trotter beads
        bool resortFlag;  //!< true implies need for resort of particles
        bool speedup; // Freeze rings in classical region?
        real maxDist;
        real dt2;
        real dt3;
        real kb;
                
        real omega2;
        
        real gamma;
        real temperature;

        //shared_ptr<VerletListAdress> verletList;
        
        std::vector< std::vector<real> > Eigenvectors;
        std::vector< std::vector<real> > Tvectors;
        
        //real maxCut;
        
        std::vector<real> tmpvals;
        
        shared_ptr< esutil::RNG > rng;
        
        //typedef std::vector<real> tuple;
        //tuple tmpvals;
        //void addEV(tuple vals); // add tuple
        
        /* TODO should be removed after signals will be tested
        shared_ptr< class Langevin > langevin;  //!< Langevin thermostat if available
        shared_ptr< class Isokinetic > isokinetic;  //!< Isokinetic thermostat if available
        shared_ptr< class StochasticVelocityRescaling > stochasticVelocityRescaling;  //!< Stochastic velocity rescaling thermostat if available
        shared_ptr< class BerendsenBarostat > berendsenBarostat;  //!< Berendsen barostat if available
        shared_ptr< class BerendsenThermostat> berendsenThermostat;  //!< Berendsen thermostat if available
        
        shared_ptr< class LangevinBarostat > langevinBarostat;  //!< Langevin-Hoover barostat if available
        
        shared_ptr< class FixPositions > fixPositions; // fix positions of a group of particles
        */

        /** Method updates particle positions and velocities.
            \return maximal square distance a particle has moved.
        */


        //real integrate1();

        //void integrate2();
        
        void integrateV1(int t);
        void integrateV2();
        void integrateModePos();
        void OUintegrate();

        void initForces();

        void updateForces(int f);

        //void calcForces();
        
        void calcForcesF();
        void calcForcesM();
        void calcForcesS();

        void transForces();
        void transPos1();
        void transPos2();
        
        //void printPositions(bool withGhost);

        //void printForces(bool withGhost);

        //void setUp();   //!< set up for a new run

        //void printTimers();

        //esutil::WallTimer timeIntegrate;  //!< used for timing

        // variables that keep time information about different phases
        /*real timeRun;
        real timeLost;
        real timeForce;
        real timeForceComp[100];
        real timeComm1;
        real timeComm2;
        real timeInt1;
        real timeInt2;
        real timeResort;*/

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
