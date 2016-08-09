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
#ifndef _INTERACTION_FIXEDPAIRLISTPIADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTPIADRESSINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "FixedPairListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "Interaction.hpp"
#include "types.hpp"
#include "FixedTupleListAdress.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _Potential >
    class FixedPairListPIadressInteractionTemplate: public Interaction, SystemAccess {
        
    protected:
      typedef _Potential Potential;
      
    public:
      FixedPairListPIadressInteractionTemplate
      (shared_ptr < System > system,
       shared_ptr < FixedPairList > _fixedpairList,
       shared_ptr <FixedTupleListAdress> _fixedtupleList,
       shared_ptr < Potential > _potential,
       int _ntrotter,
       bool _speedup)
        : SystemAccess(system), fixedpairList(_fixedpairList), fixedtupleList(_fixedtupleList),
          potential(_potential), ntrotter(_ntrotter), speedup(_speedup)
      {
        if (! potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      virtual ~FixedPairListPIadressInteractionTemplate() {};

      void
      setFixedPairList(shared_ptr < FixedPairList > _fixedpairList) {
        fixedpairList = _fixedpairList;
      }

      shared_ptr < FixedPairList > getFixedPairList() {
        return fixedpairList;
      }
      
      void
      setFixedTupleList(shared_ptr<FixedTupleListAdress> _fixedtupleList) {
          fixedtupleList = _fixedtupleList;
      }
      
      void
      setNTrotter(int _ntrotter) {
          ntrotter = _ntrotter;
      }
      
      void
      setSpeedup(bool _speedup) {
          speedup = _speedup;
      }

      void
      setPotential(shared_ptr < Potential> _potential) {
        if (_potential) {
          potential = _potential;
        } else {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      shared_ptr < Potential > getPotential() {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeEnergyAA();
      virtual real computeEnergyCG();
      virtual real computeEnergyAA(int atomtype);
      virtual real computeEnergyCG(int atomtype);
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins); 
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Pair; }

    protected:
      int ntypes;
      int ntrotter; // Trotter number
      bool speedup; // Choose whether to approximate rings in Classical region by single particles
      shared_ptr < FixedPairList > fixedpairList;
      shared_ptr<FixedTupleListAdress> fixedtupleList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedPairListPIadressInteractionTemplate < _Potential >::addForces() {
        //std::cout << "Adding Forces in FixedPairListPIadressInteractionTemplate\n";
      LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedPairList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      real ltMaxBondSqr = fixedpairList->getLongtimeMaxBondSqr();
      for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        
        //weights 
        real w1 = p1.lambda();               
        real w2 = p2.lambda(); 
        
        // Completely in classical region?
        if ( (w1 == 0.0) && (w2 == 0.0) ) {
            
            if(speedup == true){
                Real3D dist;
                bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
                Real3D force;
                real d = dist.sqr();
                if (d > ltMaxBondSqr) {
                        fixedpairList->setLongtimeMaxBondSqr(d);
                        ltMaxBondSqr = d;
                }
                if(potential->_computeForce(force, dist)) {
                  p1.force() += force;
                  p2.force() -= force;
                  LOG4ESPP_DEBUG(_Potential::theLogger, "p" << p1.id() << "(" << p1.position()[0] << "," << p1.position()[1] << "," << p1.position()[2] << ") "
                                                             << "p" << p2.id() << "(" << p2.position()[0] << "," << p2.position()[1] << "," << p2.position()[2] << ") "
                                                             << "dist=" << sqrt(dist*dist) << " "
                                                             << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
                }
            }
            else{
                // Get the corresponding tuples
                FixedTupleListAdress::iterator it3;
                FixedTupleListAdress::iterator it4;
                it3 = fixedtupleList->find(&p1);
                it4 = fixedtupleList->find(&p2);
                
                if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {
            
                    // Get the PI bead lists (i.e. the AdResS particles)
                    std::vector<Particle*> atList1;
                    std::vector<Particle*> atList2;
                    atList1 = it3->second;
                    atList2 = it4->second;

                    // Iterate the two iterators in a parallel fashion
                    std::vector<Particle*>::iterator itv2 = atList2.begin();
                    for (std::vector<Particle*>::iterator itv = atList1.begin();
                         itv != atList1.end(); ++itv) {

                         // they should be the same length... Total Trotter beads the same everywhere in the system
                         if (itv2 == atList2.end()){
                             std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                             p1.id() << "\n";
                             exit(1);
                             return;
                         }

                         // Get the individual PI beads
                         Particle &p3 = **itv;
                         Particle &p4 = **itv2;

                         // the beads we get should have the same Trotter bead number to interact with each other
                         if (p3.pib() != p4.pib()){
                             std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                             p3.id() << " and " << p4.id() << "\n";
                             exit(1);
                             return;
                         }

                         // Calculate forces
                         Real3D dist;
                         bc.getMinimumImageVectorBox(dist, p3.position(), p4.position());
                         Real3D force;
                         real d = dist.sqr();
                         if (d > ltMaxBondSqr) {
                                fixedpairList->setLongtimeMaxBondSqr(d);
                                ltMaxBondSqr = d;
                         }
                         if(potential->_computeForce(force, dist)) {
                           force *= 1.0/ntrotter;
                           p3.force() += force;
                           p4.force() -= force;
                           LOG4ESPP_DEBUG(_Potential::theLogger, "p" << p1.id() << "(" << p1.position()[0] << "," << p1.position()[1] << "," << p1.position()[2] << ") "
                                                                     << "p" << p2.id() << "(" << p2.position()[0] << "," << p2.position()[1] << "," << p2.position()[2] << ") "
                                                                     << "dist=" << sqrt(dist*dist) << " "
                                                                     << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
                         }

                         //Iterate the second iterator
                         ++itv2;

                    }               
                }
                else { // this should not happen
                   std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                         p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
                   std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
                   exit(1);
                   return;
                }  
            }
        }
        // Otherwise...
        else{ 
            // Get the corresponding tuples
            FixedTupleListAdress::iterator it3;
            FixedTupleListAdress::iterator it4;
            it3 = fixedtupleList->find(&p1);
            it4 = fixedtupleList->find(&p2);
            
            if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {
            
                // Get the PI bead lists (i.e. the AdResS particles)
                std::vector<Particle*> atList1;
                std::vector<Particle*> atList2;
                atList1 = it3->second;
                atList2 = it4->second;
            
                // Iterate the two iterators in a parallel fashion
                std::vector<Particle*>::iterator itv2 = atList2.begin();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {
                
                     // they should be the same length... Total Trotter beads the same everywhere in the system
                     if (itv2 == atList2.end()){
                         std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                         p1.id() << "\n";
                         exit(1);
                         return;
                     }

                     // Get the individual PI beads
                     Particle &p3 = **itv;
                     Particle &p4 = **itv2;

                     // the beads we get should have the same Trotter bead number to interact with each other
                     if (p3.pib() != p4.pib()){
                         std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                         p3.id() << " and " << p4.id() << "\n";
                         exit(1);
                         return;
                     }
                     
                     // Calculate forces
                     Real3D dist;
                     bc.getMinimumImageVectorBox(dist, p3.position(), p4.position());
                     Real3D force;
                     real d = dist.sqr();
                     if (d > ltMaxBondSqr) {
                            fixedpairList->setLongtimeMaxBondSqr(d);
                            ltMaxBondSqr = d;
                     }
                     if(potential->_computeForce(force, dist)) {
                       force *= 1.0/ntrotter;
                       p3.force() += force;
                       p4.force() -= force;
                       
                       LOG4ESPP_DEBUG(_Potential::theLogger, "p" << p1.id() << "(" << p1.position()[0] << "," << p1.position()[1] << "," << p1.position()[2] << ") "
                                                                 << "p" << p2.id() << "(" << p2.position()[0] << "," << p2.position()[1] << "," << p2.position()[2] << ") "
                                                                 << "dist=" << sqrt(dist*dist) << " "
                                                                 << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
                     }
                     
                     //Iterate the second iterator
                     ++itv2;
                
                }               
            }
            else { // this should not happen
               std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
               std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
               exit(1);
               return;
            }
        
        }
        
      }
    }
    
    template < typename _Potential > inline real
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeEnergy() {

      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPairList pairs");

      real e = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
	   it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        
        //weights 
        real w1 = p1.lambda();               
        real w2 = p2.lambda(); 
        
        // Completely in classical region?
        if ( (w1 == 0.0) && (w2 == 0.0) ) {
            
            if(speedup == true){
                Real3D r21;
                bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
                e += potential->_computeEnergy(r21);
            }
            else{
                // Get the corresponding tuples
                FixedTupleListAdress::iterator it3;
                FixedTupleListAdress::iterator it4;
                it3 = fixedtupleList->find(&p1);
                it4 = fixedtupleList->find(&p2);
                
                if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {
            
                    // Get the PI bead lists (i.e. the AdResS particles)
                    std::vector<Particle*> atList1;
                    std::vector<Particle*> atList2;
                    atList1 = it3->second;
                    atList2 = it4->second;

                    // Iterate the two iterators in a parallel fashion
                    std::vector<Particle*>::iterator itv2 = atList2.begin();
                    for (std::vector<Particle*>::iterator itv = atList1.begin();
                         itv != atList1.end(); ++itv) {

                         // they should be the same length... Total Trotter beads the same everywhere in the system
                         if (itv2 == atList2.end()){
                             std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                             p1.id() << "\n";
                             exit(1);
                             return 0.0;
                         }

                         // Get the individual PI beads
                         Particle &p3 = **itv;
                         Particle &p4 = **itv2;

                         // the beads we get should have the same Trotter bead number to interact with each other
                         if (p3.pib() != p4.pib()){
                             std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                             p3.id() << " and " << p4.id() << "\n";
                             exit(1);
                             return 0.0;
                         }

                         // Calculate energies
                         Real3D r21;
                         bc.getMinimumImageVectorBox(r21, p3.position(), p4.position());
                         e += (1.0/ntrotter)*potential->_computeEnergy(r21);

                         //Iterate the second iterator
                         ++itv2;
                    }               
                }
                else { // this should not happen
                   std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                         p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
                   std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
                   exit(1);
                   return 0.0;
                }
                
                
            }
        }        
        // Otherwise... 
        else{ 
            // Get the corresponding tuples
            FixedTupleListAdress::iterator it3;
            FixedTupleListAdress::iterator it4;
            it3 = fixedtupleList->find(&p1);
            it4 = fixedtupleList->find(&p2);
            
            if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {
            
                // Get the PI bead lists (i.e. the AdResS particles)
                std::vector<Particle*> atList1;
                std::vector<Particle*> atList2;
                atList1 = it3->second;
                atList2 = it4->second;
            
                // Iterate the two iterators in a parallel fashion
                std::vector<Particle*>::iterator itv2 = atList2.begin();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {
                
                     // they should be the same length... Total Trotter beads the same everywhere in the system
                     if (itv2 == atList2.end()){
                         std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                         p1.id() << "\n";
                         exit(1);
                         return 0.0;
                     }

                     // Get the individual PI beads
                     Particle &p3 = **itv;
                     Particle &p4 = **itv2;

                     // the beads we get should have the same Trotter bead number to interact with each other
                     if (p3.pib() != p4.pib()){
                         std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                         p3.id() << " and " << p4.id() << "\n";
                         exit(1);
                         return 0.0;
                     }
                
                     // Calculate energies
                     Real3D r21;
                     bc.getMinimumImageVectorBox(r21, p3.position(), p4.position());
                     e += (1.0/ntrotter)*potential->_computeEnergy(r21);
                     
                     //Iterate the second iterator
                     ++itv2;
                }               
            }
            else { // this should not happen
               std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
               std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
               exit(1);
               return 0.0;
            }
        
        }
        
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }
    
    template < typename _Potential > inline real
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in FixedPairListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }
    
    template < typename _Potential > inline real
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in FixedPairListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in FixedPairListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in FixedPairListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }
    
    template < typename _Potential >
    inline void
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
              LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");
              std::cout << "Warning! At the moment computeVirialX in FixedPairListPIadressInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
              exit(1);
              return;
    }

       
    template < typename _Potential > inline real
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeVirial() {
    	real w = 0.0;

    	const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
    	      real ltMaxBondSqr = fixedpairList->getLongtimeMaxBondSqr();
    	      for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it) {
    	        Particle &p1 = *it->first;
    	        Particle &p2 = *it->second;

    	        //weights
    	        real w1 = p1.lambda();
    	        real w2 = p2.lambda();

    	        // Completely in classical region?
    	        if ( (w1 < 0.000000001) && (w2 < 0.000000001) ) {

    	            if(speedup == true){
    	                Real3D dist;
    	                bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
    	                Real3D force;
    	                if(potential->_computeForce(force, dist)) {
    	              	  w += dist * force;
    	                  LOG4ESPP_DEBUG(_Potential::theLogger, "p" << p1.id() << "(" << p1.position()[0] << "," << p1.position()[1] << "," << p1.position()[2] << ") "
    	                                                             << "p" << p2.id() << "(" << p2.position()[0] << "," << p2.position()[1] << "," << p2.position()[2] << ") "
    	                                                             << "dist=" << sqrt(dist*dist) << " "
    	                                                             << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
    	                }
    	            }
    	            else{
    	                // Get the corresponding tuples
    	                FixedTupleListAdress::iterator it3;
    	                FixedTupleListAdress::iterator it4;
    	                it3 = fixedtupleList->find(&p1);
    	                it4 = fixedtupleList->find(&p2);

    	                if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

    	                    // Get the PI bead lists (i.e. the AdResS particles)
    	                    std::vector<Particle*> atList1;
    	                    std::vector<Particle*> atList2;
    	                    atList1 = it3->second;
    	                    atList2 = it4->second;

    	                    // Iterate the two iterators in a parallel fashion
    	                    std::vector<Particle*>::iterator itv2 = atList2.begin();
    	                    for (std::vector<Particle*>::iterator itv = atList1.begin();
    	                         itv != atList1.end(); ++itv) {

    	                         // they should be the same length... Total Trotter beads the same everywhere in the system
    	                         if (itv2 == atList2.end()){
    	                             std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " <<
    	                             p1.id() << "\n";
    	                             exit(1);
    	                             return 0.0;
    	                         }

    	                         // Get the individual PI beads
    	                         Particle &p3 = **itv;
    	                         Particle &p4 = **itv2;

    	                         // the beads we get should have the same Trotter bead number to interact with each other
    	                         if (p3.pib() != p4.pib()){
    	                             std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " <<
    	                             p3.id() << " and " << p4.id() << "\n";
    	                             exit(1);
    	                             return 0.0;
    	                         }

    	                         // Calculate forces
    	                         Real3D dist;
    	                         bc.getMinimumImageVectorBox(dist, p3.position(), p4.position());
    	                         Real3D force;
    	                         if(potential->_computeForce(force, dist)) {
    	                           force *= 1.0/ntrotter;
    	                           w += dist * force;
    	                           LOG4ESPP_DEBUG(_Potential::theLogger, "p" << p1.id() << "(" << p1.position()[0] << "," << p1.position()[1] << "," << p1.position()[2] << ") "
    	                                                                     << "p" << p2.id() << "(" << p2.position()[0] << "," << p2.position()[1] << "," << p2.position()[2] << ") "
    	                                                                     << "dist=" << sqrt(dist*dist) << " "
    	                                                                     << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
    	                         }

    	                         //Iterate the second iterator
    	                         ++itv2;

    	                    }
    	                }
    	                else { // this should not happen
    	                   std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
    	                         p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
    	                   std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
    	                   exit(1);
    	                   return 0.0;
    	                }
    	            }
    	        }
    	        else{
    	            // Get the corresponding tuples
    	            FixedTupleListAdress::iterator it3;
    	            FixedTupleListAdress::iterator it4;
    	            it3 = fixedtupleList->find(&p1);
    	            it4 = fixedtupleList->find(&p2);

    	            if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

    	                // Get the PI bead lists (i.e. the AdResS particles)
    	                std::vector<Particle*> atList1;
    	                std::vector<Particle*> atList2;
    	                atList1 = it3->second;
    	                atList2 = it4->second;

    	                // Iterate the two iterators in a parallel fashion
    	                std::vector<Particle*>::iterator itv2 = atList2.begin();
    	                for (std::vector<Particle*>::iterator itv = atList1.begin();
    	                     itv != atList1.end(); ++itv) {

    	                     // they should be the same length... Total Trotter beads the same everywhere in the system
    	                     if (itv2 == atList2.end()){
    	                         std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " <<
    	                         p1.id() << "\n";
    	                         exit(1);
    	                         return 0.0;
    	                     }

    	                     // Get the individual PI beads
    	                     Particle &p3 = **itv;
    	                     Particle &p4 = **itv2;

    	                     // the beads we get should have the same Trotter bead number to interact with each other
    	                     if (p3.pib() != p4.pib()){
    	                         std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " <<
    	                         p3.id() << " and " << p4.id() << "\n";
    	                         exit(1);
    	                         return 0.0;
    	                     }

    	                     // Calculate forces
    	                     Real3D dist;
    	                     bc.getMinimumImageVectorBox(dist, p3.position(), p4.position());
    	                     Real3D force;
    	                     if(potential->_computeForce(force, dist)) {
    	                       force *= 1.0/ntrotter;
    	                       w += dist * force;

    	                       LOG4ESPP_DEBUG(_Potential::theLogger, "p" << p1.id() << "(" << p1.position()[0] << "," << p1.position()[1] << "," << p1.position()[2] << ") "
    	                                                                 << "p" << p2.id() << "(" << p2.position()[0] << "," << p2.position()[1] << "," << p2.position()[2] << ") "
    	                                                                 << "dist=" << sqrt(dist*dist) << " "
    	                                                                 << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
    	                     }

    	                     //Iterate the second iterator
    	                     ++itv2;

    	                }
    	            }
    	            else { // this should not happen
    	               std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
    	                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
    	               std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
    	               exit(1);
    	               return 0.0;
    	            }

    	        }

    	      }
      
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

   template < typename _Potential > inline void
    FixedPairListPIadressInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");
      std::cout << "Warning! At the moment computeVirialTensor() in FixedPairListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
    }

    template < typename _Potential > inline void
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");
      std::cout << "Warning! At the moment computeVirialTensor() in FixedPairListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
    }
    
    template < typename _Potential > inline void
    FixedPairListPIadressInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");
      std::cout << "Warning! At the moment computeVirialTensor() in FixedPairListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
    }
    
    template < typename _Potential >
    inline real
    FixedPairListPIadressInteractionTemplate< _Potential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
