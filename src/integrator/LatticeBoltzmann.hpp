// ESPP_CLASS
#ifndef _INTEGRATOR_LATTICEBOLTZMANN_HPP
#define _INTEGRATOR_LATTICEBOLTZMANN_HPP

//#include "boost/serialization/vector.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "LatticeSite.hpp"

namespace espresso {
  namespace integrator {

    class LatticeBoltzmann : public Extension {
      /*
      *  Three Constructors may be used: expecting 9, 7 or 3 parameters.
      *  Obligatory parameters are the lattice size in 3D. Without further specifications,
      *  the D3Q19 model will be used with both lattice spacing and time set to 1.
      *  Furthermore, for the initialization the density of 1.0 and velocity (0.,0.,0.) at
      *  each lattice site will be set.
      *  In a more advanced constructor, user specifies lattice spacing and time AND
      *  initial density and velocity explicitly.
      *  The most advanced constructor takes also the lattice model as parameters. For this,
      *  user has to specify number of dimensions and velocities. Note that all the realisations
      *  in the current version of the code are aimed at D3Q19 model and if you want to use
      *  something else, please, add the corresponded code parts.
      *
      *  Originally, we planned this module to operate in 3D only, so if you need 2D version,
      *  a bit more tuning is needed. On the other hand, adding different 3D lattice models
      *  (such as D3Q15 or D3Q27) should be straightforward.
      *
      */
      public:
        LatticeBoltzmann (shared_ptr< System > _system, int _x, int _y, int _z,
            real _a, real _tau, real _rho0, Real3D _u0, int _numDims, int _numVels);
        LatticeBoltzmann (shared_ptr< System > _system, int _x, int _y, int _z,
            real _a, real _tau, real _rho0, Real3D _u0);
        LatticeBoltzmann (shared_ptr< System > _system, int _x, int _y, int _z);
        ~LatticeBoltzmann ();

        /* SET AND GET DECLARATION */
        void setNx(int _Nx);         	  // set lattice size in x-direction
        void setNy(int _Ny);         	  // set lattice size in y-direction
        void setNz(int _Nz);         	  // set lattice size in z-direction

        int getNx();              	      // get lattice size in x-direction
        int getNy();              	      // get lattice size in y-direction
        int getNz();                	    // get lattice size in z-direction

        void setA (real _a);             // set lattice spacing
        real getA ();                    // get lattice spacing

        void setTau (real _tau);         // set lattice timestep
        real getTau ();                  // get lattice timestep

        void setInitDen (real _rho0);	  // set initial density
        real getInitDen ();		          // get initial density

        void setInitVel (Real3D _u0);	  // set initial velocity
        Real3D getInitVel ();		        // get initial velocity

        void setNumDims (int _numDims);	// set number of dimensions
        int getNumDims ();		            // get number of dimensions

        void setNumVels (int _numVels);	// set number of velocities
        int getNumVels ();              // get number of velocities

        void setEqWeight (int _l, real _value); // set eq.weights
        real getEqWeight (int _l);      // get eq.weights

        void setCi (int _l, Real3D _vec); // set c_i's
        Real3D getCi (int _l);           // get c_i's

        void setInvBi (int _l, real _value);    // set inverse b_i's
        real getInvBi (int _l);         // get inverse b_i's
        /* END OF SET AND GET DECLARATION */

        /* FUNCTIONS DECLARATION */
        void initLatticeModel ();   // initialize lattice model (weights, cis)
        void initPopulations (real _rho0, Real3D _u0); // initialize populations
        void makeLBStep ();          // perform one step of LB
        void addPolyLBForces();     // add to polymers forces due to LBsites

        void collideStream ();      // use collide-stream scheme

        void streaming (int _i, int _j, int _k);  // streaming along the velocity vectors

        /* control functions */
        void computeDensity (int _i, int _j, int _k, int _numVels);
        void computeMomentum (int _i, int _j, int _k, int _numVels);
        /* END OF FUNCTIONS DECLARATION */

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        int numDims;          		  // number of dimensions
        int numVels;          		  // number of velocities
        real cs2;             		  // squared speed of sound
        real invCs2;          		  // inverse square of speed of sound
        real a;                     // lattice spacing
        real tau;                   // lattice timestep
        std::vector<real> eqWeight; // lattice weights
        std::vector<Real3D> c_i;    // velocity vectors
        std::vector<real> inv_b_i;  // back-transformation weights
        real rho0;	      		      // initial density
        Real3D u0;	      		      // initial velocity
        int Nx, Ny, Nz;     		    // lattice lengths in 3D
        int idX, idY, idZ, index;	  // indexes in 3D and aligned 1D index

        /* two lattices. lbfluid has f,m and meq. ghostlat has f only.
         * the latter one used for sake of simplicity during streaming
         * */
        std::vector< std::vector< std::vector<LBSite> > > lbfluid;
        std::vector< std::vector< std::vector<GhostLattice> > > ghostlat;

        boost::signals2::connection _befIntP, _befIntV;
        void connect();
        void disconnect();

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif