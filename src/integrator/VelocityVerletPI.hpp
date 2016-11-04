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

        void run(int nsteps);

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

        void addEV() { Tvectors.push_back(tmpvals); tmpvals.clear();}

        void addValues(real val) { Eigenvalues.push_back(val); } // add particle id (called from python)

        void transp();

        /** Setter routine for temperature. */
        void setTemperature(real _temperature);
        /** Getter routine for temperature. */
        real getTemperature() { return temperature; }

        /** Setter routine for gamma. */
        void setGamma(real _gamma);
        /** Getter routine for gamma. */
        real getGamma() { return gamma; }

        /** Setter routine for CMDparameter. */
        void setCMDparameter(real _CMDparameter);
        /** Getter routine for CMDparameter. */
        real getCMDparameter() { return CMDparameter; }

        /** Setter routine for CMDparameter. */
        void setPILElambda(real _PILElambda);
        /** Getter routine for CMDparameter. */
        real getPILElambda() { return PILElambda; }

        /** Setter routine for clmassmultiplier. */
        void setClmassmultiplier(real _clmassmultiplier);
        /** Getter routine for clmassmultiplier. */
        real getClmassmultiplier() { return clmassmultiplier; }

        /** Setter routine for speedup. */
        void setSpeedup(bool _speedup);
        /** Getter routine for speedup. */
        bool getSpeedup() { return speedup; }

        /** Setter routine for KTI. */
        void setKTI(bool _KTI);
        /** Getter routine for KTI. */
        bool getKTI() { return KTI; }

        /** Setter routine for centroidthermostat. */
        void setCentroidThermostat(bool _centroidthermostat);
        /** Getter routine for centroidthermostat. */
        bool getCentroidThermostat() { return centroidthermostat; }

        /** Setter routine for _PILE. */
        void setPILE(bool _PILE);
        /** Getter routine for _PILE. */
        bool getPILE() { return PILE; }

        /** Setter routine for _realkinmass. */
        void setRealKinMass(bool _realkinmass);
        /** Getter routine for _realkinmass. */
        bool getRealKinMass() { return realkinmass; }

        /** Setter routine for speedup. */
        void setConstKinMass(bool _constkinmass);
        /** Getter routine for speedup. */
        bool getConstKinMass() { return constkinmass; }

        /** Setter routine for verletList. */
        void setVerletList(shared_ptr<VerletListAdress> _verletList);
        /** Getter routine for verletList. */
        shared_ptr<VerletListAdress> getVerletList() { return verletList; }

        /** Getter routine for verletlistBuilds. */
        int getVerletlistBuilds() { return verletlistBuilds; }

        // Compute internal ring energies
        real computeRingEnergy();
        real computeRingEnergyRaw();

        // Compute special kinetic energy
        real computeKineticEnergy();

        // Compute internal ring energies
        real computePositionDrift(int parttype);

        // Compute special kinetic energy
        real computeMomentumDrift(int parttype);

        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:

        int sStep; // medium time step = sStep * dt
        int mStep; // large time step = mStep * sStep * dt
        int ntrotter; // number of Trotter beads
        int verletlistBuilds; // number of Verlet list builds from beginning
        bool resortFlag;  // true implies need for resort of particles
        bool speedup; // freeze rings in classical region?
        bool KTI; // KTI-like simulation with constant and uniform resolution?
        bool constkinmass; // Use a constant (real) kinetic mass for all modes?
        bool realkinmass; // Use real masses for kinetic masses (useful especially in TRPMD)?
        bool centroidthermostat; // Thermostat the centroid mode?
        bool PILE; // PILE thermostating with different frictions for different modes (useful especially in TRPMD)? (Ceriotti et al, JCP 133, 124104 (2010))
        real maxDist;
        real dt2;
        real dt3;
        //real kb;

        real dhy;
        real pidhy2;
        real dex;
        real dex2;
        real dexdhy;
        real dexdhy2;

        real omega2;
        real clmassmultiplier;

        real CMDparameter;
        real PILElambda;

        real gamma;
        real temperature;

        std::vector< std::vector<real> > Eigenvectors;
        std::vector< std::vector<real> > Tvectors;
        std::vector< real > Eigenvalues;

        std::vector<real> tmpvals;

        shared_ptr< esutil::RNG > rng;

        void integrateV1(int t);
        void integrateV2();
        void integrateModePos();
        void OUintegrate();

        void initForces();

        void updateForces(int f);

        void calcForcesF();
        void calcForcesM();
        void calcForcesS();

        void transForces();
        void transPos1();
        void transPos2();
        void transMom1();
        void transMom2();

        real weight(real distanceSqr);
        real weightderivative(real distanceSqr);

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
