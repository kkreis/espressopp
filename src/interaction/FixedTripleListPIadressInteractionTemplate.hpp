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
#ifndef _INTERACTION_FIXEDTRIPLELISTPIADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDTRIPLELISTPIADRESSINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedTripleList.hpp"
#include "FixedTripleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"
#include "FixedTupleListAdress.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _AngularPotential >
    class FixedTripleListPIadressInteractionTemplate : public Interaction, SystemAccess {
        
    protected:
      typedef _AngularPotential Potential;
      
    public:
      FixedTripleListPIadressInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedTripleList > _fixedtripleList,
       shared_ptr <FixedTupleListAdress> _fixedtupleList,
       shared_ptr < Potential > _potential,
       int _ntrotter,
       bool _speedup)
        : SystemAccess(_system), fixedtripleList(_fixedtripleList), fixedtupleList(_fixedtupleList),
          potential(_potential),  ntrotter(_ntrotter), speedup(_speedup)
      {
          if (! potential) {
                LOG4ESPP_ERROR(theLogger, "NULL potential");
          }
      }

      virtual ~FixedTripleListPIadressInteractionTemplate() {};

      void
      setFixedTripleList(shared_ptr < FixedTripleList > _fixedtripleList) {
        fixedtripleList = _fixedtripleList;
      }

      shared_ptr < FixedTripleList > getFixedTripleList() {
        return fixedtripleList;
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
      virtual int bondType() { return Angular; }

    protected:
      int ntypes;
      int ntrotter; // Trotter number
      bool speedup; // Choose whether to approximate rings in Classical region by single particles
      shared_ptr<FixedTripleList> fixedtripleList;
      shared_ptr<FixedTupleListAdress> fixedtupleList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _AngularPotential > inline void
    FixedTripleListPIadressInteractionTemplate <_AngularPotential>::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        
        real w1 = p1.lambda();               
        real w2 = p2.lambda();        
        real w3 = p3.lambda();        
        
        // Completely in classical region?
        if ((w1 == 0.0) && (w2 == 0.0) && (w3 == 0.0)) {
            
            if(speedup == true){
                Real3D dist12, dist32;
                bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
                bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
                Real3D force12, force32;
                potential->_computeForce(force12, force32, dist12, dist32);
                p1.force() += force12;
                p2.force() -= force12 + force32;
                p3.force() += force32;
            }
            else{
                // Get the corresponding tuples
                FixedTupleListAdress::iterator it4;
                FixedTupleListAdress::iterator it5;
                FixedTupleListAdress::iterator it6;
                it4 = fixedtupleList->find(&p1);
                it5 = fixedtupleList->find(&p2);
                it6 = fixedtupleList->find(&p3);
                
                if (it4 != fixedtupleList->end() && it5 != fixedtupleList->end() && it6 != fixedtupleList->end()) {
            
                    // Get the PI bead lists (i.e. the AdResS particles)
                    std::vector<Particle*> atList1;
                    std::vector<Particle*> atList2;
                    std::vector<Particle*> atList3;
                    atList1 = it4->second;
                    atList2 = it5->second;
                    atList3 = it6->second;

                    // Iterate the two iterators in a parallel fashion
                    std::vector<Particle*>::iterator itv2 = atList2.begin();
                    std::vector<Particle*>::iterator itv3 = atList3.begin();
                    for (std::vector<Particle*>::iterator itv = atList1.begin();
                         itv != atList1.end(); ++itv) {

                         // they should be the same length... Total Trotter beads the same everywhere in the system
                         if (itv2 == atList2.end() || itv3 == atList3.end()){
                             std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                             p1.id() << "\n";
                             exit(1);
                             return;
                         }

                         // Get the individual PI beads
                         Particle &p4 = **itv;
                         Particle &p5 = **itv2;
                         Particle &p6 = **itv3;

                         // the beads we get should have the same Trotter bead number to interact with each other
                         if (p4.pib() != p5.pib() || p5.pib() != p6.pib()){
                             std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                             p4.id() << " and " << p5.id() << " and " << p6.id() << "\n";
                             exit(1);
                             return;
                         }

                         // Calculate forces
                         Real3D dist12, dist32;
                         bc.getMinimumImageVectorBox(dist12, p4.position(), p5.position());
                         bc.getMinimumImageVectorBox(dist32, p6.position(), p5.position());
                         Real3D force12, force32;
                         potential->_computeForce(force12, force32, dist12, dist32);
                         force12 *= 1.0/ntrotter;
                         force32 *= 1.0/ntrotter;
                         p4.force() += force12;
                         p5.force() -= force12 + force32;
                         p6.force() += force32;                

                         //Iterate the second and third iterator
                         ++itv2;
                         ++itv3;

                    }               
                }
                else { // this should not happen
                   std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                         p1.ghost() << ", " << p2.id() << "-" << p2.ghost() << ", " << p3.id() << "-" << p3.ghost();
                   std::cout << " (" << p1.position() << ") (" << p2.position() << ") (" << p3.position() << ") \n";
                   exit(1);
                   return;
                }              
                
            }
            
        }
        // Otherwise...
        else{ 
            // Get the corresponding tuples
            FixedTupleListAdress::iterator it4;
            FixedTupleListAdress::iterator it5;
            FixedTupleListAdress::iterator it6;
            it4 = fixedtupleList->find(&p1);
            it5 = fixedtupleList->find(&p2);
            it6 = fixedtupleList->find(&p3);
            
            if (it4 != fixedtupleList->end() && it5 != fixedtupleList->end() && it6 != fixedtupleList->end()) {
            
                // Get the PI bead lists (i.e. the AdResS particles)
                std::vector<Particle*> atList1;
                std::vector<Particle*> atList2;
                std::vector<Particle*> atList3;
                atList1 = it4->second;
                atList2 = it5->second;
                atList3 = it6->second;
            
                // Iterate the two iterators in a parallel fashion
                std::vector<Particle*>::iterator itv2 = atList2.begin();
                std::vector<Particle*>::iterator itv3 = atList3.begin();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {
                
                     // they should be the same length... Total Trotter beads the same everywhere in the system
                     if (itv2 == atList2.end() || itv3 == atList3.end()){
                         std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                         p1.id() << "\n";
                         exit(1);
                         return;
                     }

                     // Get the individual PI beads
                     Particle &p4 = **itv;
                     Particle &p5 = **itv2;
                     Particle &p6 = **itv3;

                     // the beads we get should have the same Trotter bead number to interact with each other
                     if (p4.pib() != p5.pib() || p5.pib() != p6.pib()){
                         std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                         p4.id() << " and " << p5.id() << " and " << p6.id() << "\n";
                         exit(1);
                         return;
                     }
                     
                     // Calculate forces
                     Real3D dist12, dist32;
                     bc.getMinimumImageVectorBox(dist12, p4.position(), p5.position());
                     bc.getMinimumImageVectorBox(dist32, p6.position(), p5.position());
                     Real3D force12, force32;
                     potential->_computeForce(force12, force32, dist12, dist32);
                     force12 *= 1.0/ntrotter;
                     force32 *= 1.0/ntrotter;
                     p4.force() += force12;
                     p5.force() -= force12 + force32;
                     p6.force() += force32;                
                     
                     //Iterate the second and third iterator
                     ++itv2;
                     ++itv3;
                     
                }               
            }
            else { // this should not happen
               std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost() << ", " << p3.id() << "-" << p3.ghost();
               std::cout << " (" << p1.position() << ") (" << p2.position() << ") (" << p3.position() << ") \n";
               exit(1);
               return;
            }
        
        }

        
      }
    }

    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real e = 0.0;
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        
        //weights 
        real w1 = p1.lambda();               
        real w2 = p2.lambda();        
        real w3 = p3.lambda();        
        
        // Completely in classical region?
        if ( (w1 == 0.0) && (w2 == 0.0) && (w3 == 0.0) ) {
            
            if(speedup == true){
                Real3D dist12 = bc.getMinimumImageVector(p1.position(), p2.position());
                Real3D dist32 = bc.getMinimumImageVector(p3.position(), p2.position());
                e += potential->_computeEnergy(dist12, dist32);
            }
            else{
                // Get the corresponding tuples
                FixedTupleListAdress::iterator it4;
                FixedTupleListAdress::iterator it5;
                FixedTupleListAdress::iterator it6;
                it4 = fixedtupleList->find(&p1);
                it5 = fixedtupleList->find(&p2);
                it6 = fixedtupleList->find(&p3);                
                
                
                if (it4 != fixedtupleList->end() && it5 != fixedtupleList->end() && it6 != fixedtupleList->end()) {

                    // Get the PI bead lists (i.e. the AdResS particles)
                    std::vector<Particle*> atList1;
                    std::vector<Particle*> atList2;
                    std::vector<Particle*> atList3;
                    atList1 = it4->second;
                    atList2 = it5->second;
                    atList3 = it6->second;

                    // Iterate the two iterators in a parallel fashion
                    std::vector<Particle*>::iterator itv2 = atList2.begin();
                    std::vector<Particle*>::iterator itv3 = atList3.begin();
                    for (std::vector<Particle*>::iterator itv = atList1.begin();
                         itv != atList1.end(); ++itv) {

                         // they should be the same length... Total Trotter beads the same everywhere in the system
                         if (itv2 == atList2.end() || itv3 == atList3.end()){
                             std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                             p1.id() << "\n";
                             exit(1);
                             return 0.0;
                         }

                         // Get the individual PI beads
                         Particle &p4 = **itv;
                         Particle &p5 = **itv2;
                         Particle &p6 = **itv3;

                         // the beads we get should have the same Trotter bead number to interact with each other
                         if (p4.pib() != p5.pib() || p5.pib() != p6.pib()){
                             std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                             p4.id() << " and " << p5.id() << " and " << p6.id() << "\n";
                             exit(1);
                             return 0.0;
                         }

                         // Calculate energies
                         Real3D dist12 = bc.getMinimumImageVector(p4.position(), p5.position());
                         Real3D dist32 = bc.getMinimumImageVector(p6.position(), p5.position());
                         e += (1.0/ntrotter)*potential->_computeEnergy(dist12, dist32);

                         // Iterate the second and third iterator
                         ++itv2;
                         ++itv3;
                    }               
                }
                else { // this should not happen
                   std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                         p1.ghost() << ", " << p2.id() << "-" << p2.ghost() << ", " << p3.id() << "-" << p3.ghost();
                   std::cout << " (" << p1.position() << ") (" << p2.position() << ") (" << p3.position() << ") \n";
                   exit(1);
                   return 0.0;
                } 
            }
        }
        // Otherwise...
        else{  
            // Get the corresponding tuples
            FixedTupleListAdress::iterator it4;
            FixedTupleListAdress::iterator it5;
            FixedTupleListAdress::iterator it6;
            it4 = fixedtupleList->find(&p1);
            it5 = fixedtupleList->find(&p2);
            it6 = fixedtupleList->find(&p3);
            
            if (it4 != fixedtupleList->end() && it5 != fixedtupleList->end() && it6 != fixedtupleList->end()) {
            
                // Get the PI bead lists (i.e. the AdResS particles)
                std::vector<Particle*> atList1;
                std::vector<Particle*> atList2;
                std::vector<Particle*> atList3;
                atList1 = it4->second;
                atList2 = it5->second;
                atList3 = it6->second;
            
                // Iterate the two iterators in a parallel fashion
                std::vector<Particle*>::iterator itv2 = atList2.begin();
                std::vector<Particle*>::iterator itv3 = atList3.begin();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {
                
                     // they should be the same length... Total Trotter beads the same everywhere in the system
                     if (itv2 == atList2.end() || itv3 == atList3.end()){
                         std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                         p1.id() << "\n";
                         exit(1);
                         return 0.0;
                     }

                     // Get the individual PI beads
                     Particle &p4 = **itv;
                     Particle &p5 = **itv2;
                     Particle &p6 = **itv3;

                     // the beads we get should have the same Trotter bead number to interact with each other
                     if (p4.pib() != p5.pib() || p5.pib() != p6.pib()){
                         std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                         p4.id() << " and " << p5.id() << " and " << p6.id() << "\n";
                         exit(1);
                         return 0.0;
                     }
                     
                     // Calculate energies
                     Real3D dist12 = bc.getMinimumImageVector(p4.position(), p5.position());
                     Real3D dist32 = bc.getMinimumImageVector(p6.position(), p5.position());
                     e += (1.0/ntrotter)*potential->_computeEnergy(dist12, dist32);
                     
                     // Iterate the second and third iterator
                     ++itv2;
                     ++itv3;
                }               
            }
            else { // this should not happen
               std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost() << ", " << p3.id() << "-" << p3.ghost();
               std::cout << " (" << p1.position() << ") (" << p2.position() << ") (" << p3.position() << ") \n";
               exit(1);
               return 0.0;
            }
        
        }
       
        //const Potential &potential = getPotential(p1.type(), p2.type());
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }
    
    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }
    
    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }
           
    template < typename _AngularPotential >
    inline void
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
        exit(1);
        return;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");

      real w = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;

        //weights
        real w1 = p1.lambda();
        real w2 = p2.lambda();
        real w3 = p3.lambda();

        // Completely in classical region?
        if ( (w1 < 0.000000001) && (w2 < 0.000000001) && (w3 < 0.000000001) ) {

            if(speedup == true){
                Real3D dist12, dist32;
                bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
                bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
                Real3D force12, force32;
                potential->_computeForce(force12, force32, dist12, dist32);
                w += dist12 * force12 + dist32 * force32;
            }
            else{
                // Get the corresponding tuples
                FixedTupleListAdress::iterator it4;
                FixedTupleListAdress::iterator it5;
                FixedTupleListAdress::iterator it6;
                it4 = fixedtupleList->find(&p1);
                it5 = fixedtupleList->find(&p2);
                it6 = fixedtupleList->find(&p3);

                if (it4 != fixedtupleList->end() && it5 != fixedtupleList->end() && it6 != fixedtupleList->end()) {

                    // Get the PI bead lists (i.e. the AdResS particles)
                    std::vector<Particle*> atList1;
                    std::vector<Particle*> atList2;
                    std::vector<Particle*> atList3;
                    atList1 = it4->second;
                    atList2 = it5->second;
                    atList3 = it6->second;

                    // Iterate the two iterators in a parallel fashion
                    std::vector<Particle*>::iterator itv2 = atList2.begin();
                    std::vector<Particle*>::iterator itv3 = atList3.begin();
                    for (std::vector<Particle*>::iterator itv = atList1.begin();
                         itv != atList1.end(); ++itv) {

                         // they should be the same length... Total Trotter beads the same everywhere in the system
                         if (itv2 == atList2.end() || itv3 == atList3.end()){
                             std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " <<
                             p1.id() << "\n";
                             exit(1);
                             return 0.0;
                         }

                         // Get the individual PI beads
                         Particle &p4 = **itv;
                         Particle &p5 = **itv2;
                         Particle &p6 = **itv3;

                         // the beads we get should have the same Trotter bead number to interact with each other
                         if (p4.pib() != p5.pib() || p5.pib() != p6.pib()){
                             std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " <<
                             p4.id() << " and " << p5.id() << " and " << p6.id() << "\n";
                             exit(1);
                             return 0.0;
                         }

                         // Calculate forces
                         Real3D dist12, dist32;
                         bc.getMinimumImageVectorBox(dist12, p4.position(), p5.position());
                         bc.getMinimumImageVectorBox(dist32, p6.position(), p5.position());
                         Real3D force12, force32;
                         potential->_computeForce(force12, force32, dist12, dist32);
                         force12 *= 1.0/ntrotter;
                         force32 *= 1.0/ntrotter;
                         w += dist12 * force12 + dist32 * force32;

                         //Iterate the second and third iterator
                         ++itv2;
                         ++itv3;

                    }
                }
                else { // this should not happen
                   std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                         p1.ghost() << ", " << p2.id() << "-" << p2.ghost() << ", " << p3.id() << "-" << p3.ghost();
                   std::cout << " (" << p1.position() << ") (" << p2.position() << ") (" << p3.position() << ") \n";
                   exit(1);
                   return 0.0;
                }

            }

        }
        // Otherwise...
        else{
            // Get the corresponding tuples
            FixedTupleListAdress::iterator it4;
            FixedTupleListAdress::iterator it5;
            FixedTupleListAdress::iterator it6;
            it4 = fixedtupleList->find(&p1);
            it5 = fixedtupleList->find(&p2);
            it6 = fixedtupleList->find(&p3);

            if (it4 != fixedtupleList->end() && it5 != fixedtupleList->end() && it6 != fixedtupleList->end()) {

                // Get the PI bead lists (i.e. the AdResS particles)
                std::vector<Particle*> atList1;
                std::vector<Particle*> atList2;
                std::vector<Particle*> atList3;
                atList1 = it4->second;
                atList2 = it5->second;
                atList3 = it6->second;

                // Iterate the two iterators in a parallel fashion
                std::vector<Particle*>::iterator itv2 = atList2.begin();
                std::vector<Particle*>::iterator itv3 = atList3.begin();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {

                     // they should be the same length... Total Trotter beads the same everywhere in the system
                     if (itv2 == atList2.end() || itv3 == atList3.end()){
                         std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " <<
                         p1.id() << "\n";
                         exit(1);
                         return 0.0;
                     }

                     // Get the individual PI beads
                     Particle &p4 = **itv;
                     Particle &p5 = **itv2;
                     Particle &p6 = **itv3;

                     // the beads we get should have the same Trotter bead number to interact with each other
                     if (p4.pib() != p5.pib() || p5.pib() != p6.pib()){
                         std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " <<
                         p4.id() << " and " << p5.id() << " and " << p6.id() << "\n";
                         exit(1);
                         return 0.0;
                     }

                     // Calculate forces
                     Real3D dist12, dist32;
                     bc.getMinimumImageVectorBox(dist12, p4.position(), p5.position());
                     bc.getMinimumImageVectorBox(dist32, p6.position(), p5.position());
                     Real3D force12, force32;
                     potential->_computeForce(force12, force32, dist12, dist32);
                     force12 *= 1.0/ntrotter;
                     force32 *= 1.0/ntrotter;
                     w += dist12 * force12 + dist32 * force32;

                     //Iterate the second and third iterator
                     ++itv2;
                     ++itv3;

                }
            }
            else { // this should not happen
               std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost() << ", " << p3.id() << "-" << p3.ghost();
               std::cout << " (" << p1.position() << ") (" << p2.position() << ") (" << p3.position() << ") \n";
               exit(1);
               return 0.0;
            }

        }


      }
        
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _AngularPotential > inline void
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");
      std::cout << "Warning! At the moment computeVirialTensor in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
    }

    template < typename _AngularPotential > inline void
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");
      std::cout << "Warning! At the moment IK computeVirialTensor in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
    }
    
    template < typename _AngularPotential > inline void
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor *w, int n) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");
      std::cout << "Warning! At the moment IK computeVirialTensor in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
    }
    
    template < typename _AngularPotential >
    inline real
    FixedTripleListPIadressInteractionTemplate< _AngularPotential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif