/** @file benchmark.cpp

    A simple benchmark of the all-pairs loop operating on an unsorted
    ParticleStorage.
 */

#define LOG4ESPP_LEVEL_WARN

#include <iostream>
#include <iomanip>
#include <vector>

#include "espresso_common.hpp"
#include "types.hpp"
#include "bc/PBC.hpp"
#include "particlestorage/ParticleStorage.hpp"
#include "particleset/All.hpp"
#include "pairs/All.hpp"
#include "pairs/PairForceComputer.hpp"
#include "interaction/LennardJones.hpp"
#include "esutil/Timer.hpp"

using namespace std;
using namespace espresso;

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif


/// number of particles in each dimension
const int N = 20;
/// dimension of the cubic simulation box
const real size = 5.0;

class TestEspresso {
public:
    typedef particlestorage::ParticleStorage Storage;
    Storage storage;
    size_t position, force;
    size_t npart;

    TestEspresso(size_t nparticles): npart(nparticles) {
        position = storage.addProperty<Real3D>();
        force = storage.addProperty<Real3D>();
    }

    void addParticle(const Real3D &pos);

    void calculateForces(real epsilon, real sigma, real cutoff) NOINLINE;

    void calculateForces2(real epsilon, real sigma, real cutoff) NOINLINE;

    void runEmptyPairLoop() NOINLINE;

    real calculateMinDist() NOINLINE;

    void runEmptyLoop() NOINLINE;

    real calculateAverage() NOINLINE;

    Real3D getForce(size_t i) {
        // HACK!
        Storage::reference ref = storage.getParticleByID(i + 1);
        return storage.getProperty<Real3D>(force)[ref];
    }
};

void TestEspresso::addParticle(const Real3D &pos)
{
    Storage::reference ref = storage.addParticle();
    Storage::PropertyTraits<Real3D>::Reference positionRef = storage.getProperty<Real3D>(position);
    Storage::PropertyTraits<Real3D>::Reference forceRef    = storage.getProperty<Real3D>(force);

    positionRef[ref] = pos;
    forceRef[ref] = 0.0;
}

void TestEspresso::calculateForces(real epsilon, real sigma, real cutoff) {
    bc::PBC pbc(size);
    particleset::All allset(&storage);
    pairs::All allpairs(pbc, allset, position);
    interaction::LennardJones ljint;
    ljint.setCutoff(cutoff);
    ljint.setEpsilon(epsilon);
    ljint.setSigma(sigma);
    pairs::PairForceComputer forcecompute(storage.getProperty<Real3D>(force), ljint);
    allpairs.foreach(forcecompute);
}

class MyLJComputer: public pairs::ParticlePairComputer {
public:
    typedef particlestorage::ParticleStorage::PropertyTraits<Real3D>::Reference ForceReference;

    ForceReference force;

    real sigma;
    real epsilon;
    real cutoffSqr;

    MyLJComputer(real _sigma, real _epsilon, real cutoff,
                 ForceReference _force)
        : force(_force), sigma(_sigma), epsilon(_epsilon), cutoffSqr(cutoff*cutoff) {}

    virtual void operator()(const Real3D &dist,
			    particlestorage::ParticleStorage::reference ref1,
			    particlestorage::ParticleStorage::reference ref2) {
        real distSqr = dist.sqr();

        if (distSqr < cutoffSqr) {
            Real3D f(0.0, 0.0, 0.0);
            real frac2 = sigma / distSqr;
            real frac6 = frac2 * frac2 * frac2;
            real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5 * frac6) * frac2;
            f = dist * ffactor;
            force[ref1] += f;
            force[ref2] -= f;
        }
    }
};

void TestEspresso::calculateForces2(real epsilon, real sigma, real cutoff) {
    Storage::PropertyTraits<Real3D>::Reference forceRef = storage.getProperty<Real3D>(force);

    bc::PBC pbc(size);
    particleset::All allset(&storage);
    pairs::All allpairs(pbc, allset, position);
    MyLJComputer ljc(sigma, epsilon, cutoff, forceRef);
    allpairs.foreach(ljc);
}

class EmptyPairComputer: public pairs::ConstParticlePairComputer {
public:
    EmptyPairComputer() {}

    virtual void operator()(const Real3D &dist,
			    particlestorage::ParticleStorage::const_reference ref1,
			    particlestorage::ParticleStorage::const_reference ref2) {
    }
};

void TestEspresso::runEmptyPairLoop() {
    bc::PBC pbc(size);
    particleset::All allset(&storage);
    pairs::All allpairs(pbc, allset, position);
    EmptyPairComputer ljc;
    allpairs.foreach(ljc);
}

class MinDistComputer: public pairs::ConstParticlePairComputer {
public:
    real min;

    MinDistComputer(): min(1e10) {}

    virtual void operator()(const Real3D &dist,
			    particlestorage::ParticleStorage::const_reference ref1,
			    particlestorage::ParticleStorage::const_reference ref2) {
	real d = dist.sqr();
	if (min > d) {
	    min = d;
	}
    }
};

real TestEspresso::calculateMinDist() {
    bc::PBC pbc(size);
    particleset::All allset(&storage);
    pairs::All allpairs(pbc, allset, position);
    MinDistComputer mincomp;
    allpairs.foreach(mincomp);
    return sqrt(mincomp.min);
}

class AverageComputer: public particlestorage::ConstParticleComputer {
public:
    typedef particlestorage::ParticleStorage::PropertyTraits<Real3D>::ConstReference Reference;

    Reference property;
    real average;

    AverageComputer(const Reference &_property): property(_property), average(0) {}

    virtual void operator()(particlestorage::ParticleStorage::const_reference ref) {
        const Real3D p = property[ref];
        average += sqrt(p.sqr());
    }
};

real TestEspresso::calculateAverage() {
    particleset::All allset(&storage);
    AverageComputer avgcompute(storage.getProperty<Real3D>(force));
    allset.foreach(avgcompute);
    return avgcompute.average;
}

class EmptyComputer: public particlestorage::ConstParticleComputer {
public:
    EmptyComputer() {}

    virtual void operator()(particlestorage::ParticleStorage::const_reference ref) {
    }
};

void TestEspresso::runEmptyLoop() {
    particleset::All allset(&storage);
    EmptyComputer avgcompute;
    allset.foreach(avgcompute);
}

class TestBasic {
public:
    vector<Real3D> position, force;
    size_t npart;

    TestBasic(size_t nparticles): npart(nparticles) {
        position.reserve(nparticles);
        force.reserve(nparticles);
    }

    void addParticle(const Real3D &pos);

    void calculateForces(real epsilon, real sigma, real cutoff) NOINLINE;

    real calculateMinDist() NOINLINE;

    real calculateAverage() NOINLINE;

    Real3D getForce(size_t i) {
        return force[i];
    }
};

void TestBasic::addParticle(const Real3D &pos) {
    position.push_back(pos);
    force.push_back(0.0);
}

void TestBasic::calculateForces(real epsilon, real sigma, real cutoff) {
    real cutoffSqr = cutoff*cutoff;
    bc::PBC pbc(size);

    for (size_t i = 0; i < npart; ++i) {
        for (size_t j = i+1; j < npart; ++j) {
            Real3D pos1 = position[i];
            Real3D pos2 = position[j];
            Real3D dist = pbc.getDist(pos1, pos2);

            Real3D f(0.0, 0.0, 0.0);
            {
                real   frac2;
                real   frac6;
                real distSqr = dist.sqr();
                
                if (distSqr < cutoffSqr) {
                    frac2 = sigma / distSqr;
                    frac6 = frac2 * frac2 * frac2;
                    real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5 * frac6) * frac2;
                    f = dist * ffactor;
                } 
            }
            force[i] += f;
            force[j] -= f;
        }
    }
}

real TestBasic::calculateMinDist() {
    real min = 1e10;

    bc::PBC pbc(size);
    for (size_t i = 0; i < npart; ++i) {
        for (size_t j = i+1; j < npart; ++j) {
            Real3D pos1 = position[i];
            Real3D pos2 = position[j];
            real d = pbc.getDist(pos1, pos2).sqr();
	    if (min > d) {
		min = d;
            }
        }
    }
    return sqrt(min);
}

real TestBasic::calculateAverage() {
    real average = 0;
    for (size_t i = 0; i < npart; ++i) {
        Real3D p = force[i];
        average += sqrt(p.sqr());
    }
    return average;
}

template<class Test>
void generateParticles(Test &test) {
    srand48(123);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) { 
            for (int k = 0; k < N; k++) {
      
                real r;
                r = 0.4 + 0.2 * drand48();
                Real3D pos = Real3D(
                    (i + r) / N * size,
                    (j + r) / N * size, 
                    (k + r) / N * size);

                test.addParticle(pos);
            }
        }
    }
}

/// this routine runs the tests as defined above
int main()
{
    IF_MPI(initMPI());
    esutil::WallTimer timer;
    TestEspresso espresso(N*N*N);
    cout << setw(20) << "setup Espresso: " << timer << endl;

    timer.reset();
    TestBasic basic(N*N*N);
    cout << setw(20) << "setup Basic: " << timer << endl;

    // generate particles in the particle storage

    timer.reset();
    generateParticles(espresso);
    cout << setw(20) << "generate Espresso: " << timer << endl;

    timer.reset();
    generateParticles(basic);
    cout << setw(20) << "generate Basic: " << timer << endl;

    cout << endl << "PARTICLE PAIR LOOPING TESTS" << endl;

    // check empty pair loop

    timer.reset();
    espresso.runEmptyPairLoop();
    cout << setw(30) << "empty pair loop: " << timer << endl << endl;

    // calculate forces

    timer.reset();
    basic.calculateForces(1.1, 1.2, 2.5);
    real basicCalcTime = timer.getElapsedTime();
    cout << setw(30) << "calc basic: " << timer << endl;

    timer.reset();
    espresso.calculateForces(1.1, 1.2, 2.5);
    cout << setw(30) << "calc Espresso: " << timer << endl;
    cout << "RATIO: " << (timer.getElapsedTime() / basicCalcTime) << endl;

    timer.reset();
    espresso.calculateForces2(1.1, 1.2, 2.5);
    cout << setw(30) << "calc Espresso/LJ-Computer: " << timer << endl;
    cout << "RATIO: " << (timer.getElapsedTime() / basicCalcTime) << endl;

    // calculate minimum distance

    timer.reset();
    real mine = espresso.calculateMinDist();
    real espressoMinTime = timer.getElapsedTime();
    cout << setw(30) << "min Espresso: " << timer << endl;

    timer.reset();
    real minb = basic.calculateMinDist();
    cout << setw(30) << "min Basic: " << timer << endl; 
    cout << "RATIO: " << (espressoMinTime / timer.getElapsedTime()) << endl;

    cout << endl << "PARTICLE LOOPING TESTS" << endl;

    // check empty loop

    timer.reset();
    for (size_t cnt = 0; cnt < 10000; ++cnt)
        espresso.runEmptyLoop();
    cout << setw(30) << "empty loop: " << timer << endl << endl;

    // calculate average

    timer.reset();
    real ave;
    for (size_t cnt = 0; cnt < 10000; ++cnt)
        ave = espresso.calculateAverage();
    real espressoAvgTime = timer.getElapsedTime();
    cout << setw(30) << "average Espresso: " << timer << endl;

    timer.reset();
    real avb;
    for (size_t cnt = 0; cnt < 10000; ++cnt)
        avb = basic.calculateAverage();
    cout << setw(30) << "average Basic: " << timer << endl;
    cout << "RATIO: " << (espressoAvgTime / timer.getElapsedTime()) << endl;

    // check consistency
    cout << "min dists: " << mine << " " << minb << endl;
    if (abs(mine-minb)/abs(mine) > 1e-5) {
        cerr << "ERROR: minima are different: " << mine << " != " << minb << endl;        
    }

    // take into account that Espresso calculates forces twice to test two algorithms
    cout << "average force: " << 0.5*ave << " " << avb << endl;
    if (abs(0.5*ave-avb)/abs(ave) > 1e-5) {
        cerr << "ERROR: averages are different: " << ave << " != " << avb << endl;        
    }

    for (size_t i = 0; i < N*N*N; ++i) {
        Real3D f1 = espresso.getForce(i);
        Real3D f2 = basic.getForce(i);
        real diff = sqrt((f1-f2).sqr());
        if (diff/abs(f1.sqr()) > 1e-5) {
            cerr << "ERROR: difference " << diff << " too big for particle " << i << endl;
            cerr << "ERROR: " << f1.getX() << " vs. " << f2.getX() << endl;
            cerr << "ERROR: " << f1.getY() << " vs. " << f2.getY() << endl;
            cerr << "ERROR: " << f1.getZ() << " vs. " << f2.getZ() << endl;
        }
    }
    IF_MPI(finalizeMPI());
}