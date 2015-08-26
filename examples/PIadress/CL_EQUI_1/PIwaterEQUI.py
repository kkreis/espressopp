#!/usr/bin/env python
# -*- coding: utf-8 -*-

# relevant imports
import math
import sys
import time
import espressopp
import mpi4py.MPI as MPI
#import MPI4PY.MPI as MPI
import numpy as np
import random as rand
from numpy import linalg as la
from scipy.sparse import dia_matrix

from espressopp import Real3D, Int3D
from espressopp.tools import decomp
from espressopp.tools import timers

# integration steps, cutoff, skin, AdResS specifications
steps =50000#500
intervals = 500#500
timestep_short = 0.0005/4.0 #0.002/50.0
timestepdevide_medium = 2#5
timestepdevide_short = 2
rc = 0.98 #0.9 # cutoff coarse-grained potential
rca = 0.78 # cutoff atomistic potential (cutoff (2^(1/6)), WCA)
skin = 0.26 # skin

# Trotter number
nTrotter = 8#48

# Parameters for the thermostat
gamma = 2.0#2.0#0.1#1.0#0.0#0.01#1.0
temp = 2.494353 # 300 K* 0.00831451

# Parameters for size of AdResS dimensions
ex_size = 500#1.5#500.0
hy_size = 2.0#5.0

# Speedup in classical region for interatomic interactions?
speedupInterAtom = False
# Speedup in classical region by freezing rings?
speedupFreezeRings = False

# read equilibrated configuration file
#pid, types, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz("equilibrated_conf.xyz")
pid, types, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz("input.xyz")
#pid, types, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz("single_mol.xyz")

# Make masses
masses = []
for item in types:
    if item == 1:
        masses.append(15.9994)
    else:
        masses.append(1.008)

# Table for potentials
tabAngle = "POTS/tableESP_angle.dat"
tabBondHH = "POTS/tableESP_bondHH.dat"
tabBondOH = "POTS/tableESP_bondOH.dat"
tabHW_HW = "POTS/tableESP_HW_HW.dat"
tabHW_OW = "POTS/tableESP_HW_OW.dat"
tabOW_OW = "POTS/tableESP_OW_OW.dat"

# number of CG particles
num_particlesCG = len(x)

# number of AT particles
num_particlesAT = len(x)*nTrotter

# density, size
density = num_particlesCG / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

# System, boundary conditions, skin, communicator, node & cell grids
print 'Setting up simulation ...\n'
system = espressopp.System()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)

# Random Number Generator
xs = time.time()
seed = int(xs % int(xs) * 10000000000)
rng = espressopp.esutil.RNG()
rng.seed(seed)
system.rng = rng

# (H-)AdResS domain decomposition
system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

# add particles to the system and then decompose
props = ['id', 'pos', 'v', 'f', 'pib', 'type', 'mass', 'adrat']
allParticlesAT = []
allParticles = []
tuples = []

# prepare AT particles
for pidCG in range(num_particlesCG):
    for trot in range(nTrotter):
        allParticlesAT.append([pidCG*nTrotter + trot + 1, # add here these particles just temporarily
                            Real3D(x[pidCG]+ 0.01*(rand.random()-0.5), y[pidCG]+ 0.01*(rand.random()-0.5), z[pidCG]+ 0.01*(rand.random()-0.5)), # position
                            #Real3D(x[pidCG], y[pidCG], z[pidCG]), # position
                            Real3D(0, 0, 0), # velocity
                            Real3D(0, 0, 0), # force
                            trot+1, types[pidCG], masses[pidCG], 1]) # pib, type, mass, is AT particle

# create CG particles
for pidCG in range(num_particlesCG):

    # Preparation of tuples (tuples define, which atoms/trotter beads belong to which CG molecules/CG atoms)
    tmptuple = [pidCG+num_particlesAT+1]
    for pidAT2 in range(nTrotter):
        pid = pidCG*nTrotter+pidAT2
        tmptuple.append((allParticlesAT[pid])[0])
    firsParticleId=tmptuple[1]
    cmp=allParticlesAT[firsParticleId-1][1]

    #typeCG=max(types)+1
    # append CG particles
    allParticles.append([pidCG+num_particlesAT+1, # CG particle has to be added first!
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(0, 0, 0), # vel
                         Real3D(0, 0, 0), # force
                         0, types[pidCG], masses[pidCG], 0]) # pib, type, mass, is not AT particle
    # append AT particles
    for pidAT in range(nTrotter):
        pid = pidCG*nTrotter+pidAT
        allParticles.append([(allParticlesAT[pid])[0], # now the AT particles can be added
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # force
                            (allParticlesAT[pid])[4], # pib
                            (allParticlesAT[pid])[5], # type
                            (allParticlesAT[pid])[6], # mass
                            (allParticlesAT[pid])[7]]) # is AT particle
    # append tuple to tuplelist
    tuples.append(tmptuple)

# add particles to system
system.storage.addParticles(allParticles, *props)

# create FixedTupleList object
ftpl = espressopp.FixedTupleListAdress(system.storage)

# and add the tuples
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

# decompose
print "Added particles and tuples, decompose ..."
system.storage.decompose()
print "Decomposing done\n"

# create bond lists between CG particles
bondsOH = []
bondsHH = []
for part in range(num_particlesCG/3):
    bondsOH.append((num_particlesAT + 1 + 3*part, num_particlesAT + 1 + 3*part+1))
    bondsOH.append((num_particlesAT + 1 + 3*part, num_particlesAT + 1 + 3*part+2))
    bondsHH.append((num_particlesAT + 1 + 3*part+1, num_particlesAT + 1 + 3*part+2))

# create Verlet List exclusions
excl = list(bondsOH)
for item in bondsHH:
    excl.append(item)

# add bonds between CG particles
fplOH = espressopp.FixedPairList(system.storage)
fplHH = espressopp.FixedPairList(system.storage)
fplOH.addBonds(bondsOH)
fplHH.addBonds(bondsHH)

# create a verlet list
vl = espressopp.VerletListAdress(system, cutoff=rc, adrcut=rc,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=[Lx/2, Ly/2, Lz/2], exclusionlist=excl)

# create angles between CG particles
angles = []
for part in range(num_particlesCG/3):
    angles.append((num_particlesAT + 1 + 3*part+1, num_particlesAT + 1 + 3*part, num_particlesAT + 1 + 3*part+2))

# add angles between CG particles
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)

# decompose
print "Added bonds, angles, and VerletList, decompose ..."
system.storage.decompose()
print "Decomposing done\n"

# non-bonded potentials
interNB = espressopp.interaction.VerletListPIadressTabulated(vl, ftpl, nTrotter, speedupInterAtom) # Here we need specific PI AdResS interaction type
potOO = espressopp.interaction.Tabulated(itype=3, filename=tabOW_OW, cutoff=rca)
potHO = espressopp.interaction.Tabulated(itype=3, filename=tabHW_OW, cutoff=rca)
potHH = espressopp.interaction.Tabulated(itype=3, filename=tabHW_HW, cutoff=rca)
interNB.setPotentialQM(type1=1, type2=1, potential=potOO)
interNB.setPotentialCL(type1=1, type2=1, potential=potOO)
interNB.setPotentialQM(type1=1, type2=0, potential=potHO)
interNB.setPotentialCL(type1=1, type2=0, potential=potHO)
interNB.setPotentialQM(type1=0, type2=0, potential=potHH)
interNB.setPotentialCL(type1=0, type2=0, potential=potHH)
system.addInteraction(interNB)

# bonded potentials
# Quartic potential between AT particles
potBondHH = espressopp.interaction.Tabulated(itype=3, filename=tabBondHH)
potBondOH = espressopp.interaction.Tabulated(itype=3, filename=tabBondOH)
interBondedHH = espressopp.interaction.FixedPairListPIadressTabulated(system, fplHH, ftpl, potBondHH, nTrotter, speedupInterAtom)
interBondedOH = espressopp.interaction.FixedPairListPIadressTabulated(system, fplOH, ftpl, potBondOH, nTrotter, speedupInterAtom)
system.addInteraction(interBondedHH)
system.addInteraction(interBondedOH)

# angle potentials
potAngle = espressopp.interaction.TabulatedAngular(itype=3, filename=tabAngle)
interAngle = espressopp.interaction.FixedTripleListPIadressTabulatedAngular(system, ftl, ftpl, potAngle, nTrotter, speedupInterAtom)
system.addInteraction(interAngle)

# set up matrix for diagonalization
diag = np.zeros(nTrotter) + 2
#print diag
#print ''
off_diag = np.zeros(nTrotter) - 1
#print off_diag
#print ''
matrix = dia_matrix(([diag, off_diag, off_diag, off_diag, off_diag], [0, 1, -1, nTrotter-1, -1*(nTrotter-1)]), shape=(nTrotter, nTrotter)).toarray()
print "Matrix to be diagonalized for decomposition into normal modes:"
print matrix
print ''

# calculate Eigenvectors and -values
eigenvals, eigenvecs = la.eigh(matrix)
evals = eigenvals.tolist()
evecs = eigenvecs.tolist()

# Make sure, evals and evecs start with first centroid mode
zero_index = evals.index(min(evals))
if evals.index(min(evals)) != 0:
    evals[0], evals[zero_index] = evals[zero_index], evals[0]
    for vec in evecs:
        vec[0], vec[zero_index] = vec[zero_index], vec[0]
    #evecs[0], evecs[zero_index] = evecs[zero_index], evecs[0]
    print "Eigenvectors and Eigenvals switched for indices {0} and {1}".format(0,zero_index)

# print Eigenvectors and -values
print "Eigenvalues:\n", evals, "\n"
#print evals
#print ''
print "Transposed Eigenvectors:"
for item in evecs:
    print item
print ''

# VelocityVerlet integrator
integrator = espressopp.integrator.VelocityVerletPI(system, vl)
integrator.setmStep(timestepdevide_medium)
integrator.setsStep(timestepdevide_short)
integrator.setTimeStep(timestep_short)
integrator.setNtrotter(nTrotter)
integrator.setTemperature(temp)
integrator.setGamma(gamma)
integrator.setSpeedup(speedupFreezeRings)
integrator.addEigenvectors(evecs)
integrator.addEigenvalues(evals)

# add AdResS extension
adress = espressopp.integrator.AdressPI(system, vl, ftpl, nTrotter, KTI=True)
integrator.addExtension(adress)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass
espressopp.tools.AdressDecomp(system, integrator)

# Set lambdas and derivates to zero
print "Set lambdas, lambda derivatives to zero, and varmasses on large masses..."
for i in range(1, num_particlesAT + num_particlesCG + 1):
    system.storage.modifyParticle(i, 'lambda_adrd', 0.0)     ### TO ARRANGE!! ###
    system.storage.modifyParticle(i, 'lambda_adr', 0.0)     ### TO ARRANGE!! ###
    system.storage.modifyParticle(i, 'varmass', 100.0*system.storage.getParticle(i).mass)
system.storage.decompose()


print "System setup done, information:"

# system information
print ''
print 'AdResS Center =', [Lx/2, Ly/2, Lz/2]
print 'Size of high resolution region', ex_size
print 'Size of hybrid region', hy_size
print 'Trotter number =', nTrotter
print 'Number of total Trotter beads =', num_particlesAT
print 'Number of atomistic particles =', num_particlesCG
print 'Atomistic density = %.4f' % (density)
print ''
print 'rc =', rc
print 'rca =', rca
print 'skin =', system.skin
print ''
print 'Short timestep =', integrator.dt
print 'Medium timestep = ', integrator.dt * timestepdevide_medium
print 'Long timestep = ', integrator.dt * timestepdevide_medium * timestepdevide_short
print 'Outer steps =', steps
print 'Intervals = ', intervals
print ''
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print 'Temperature =', temp
print 'Gamma = ', gamma
print 'Using centers of mass in classical region for force calculations?', speedupInterAtom
print 'Freezing internal ring vibrations in classical region?', speedupFreezeRings
print ''
print ''
print "Starting the integration loop"
print ''

# analysis
temperature = espressopp.analysis.Temperature(system)

# output
outfile = open("esp.dat", "w")

# Output format for screen and file
print "step, time (ps), T, Eb, EAng, EPI, ELj, Ek, Etotal"
fmt='%8d %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g\n'

# initial configuration
T = temperature.compute() * 120.27267#* 1.5
Eb = interBondedHH.computeEnergy() + interBondedOH.computeEnergy()
EAng = interAngle.computeEnergy()
ELj= interNB.computeEnergy()
Ek = integrator.computeKineticEnergy()#0.5 * T * (3 * num_particles)
EPI = integrator.computeRingEnergy()
Etotal = Ek+Eb+EAng+ELj+EPI
outfile.write(fmt%(0, 0, T, Eb, EAng, EPI, ELj, Ek, Etotal))
print ((fmt%(0, 0, T, Eb, EAng, EPI, ELj, Ek, Etotal)))

# Timer, Steps
nsteps = steps / intervals

# write the start configuration to gro trajectory
dump_conf_gro = espressopp.io.DumpGRO(system, integrator, filename='trajCG.gro')
dump_conf_gro_adr = espressopp.io.DumpGROAdress(system, ftpl, integrator, filename='trajAT.gro')
dump_conf_gro.dump()
dump_conf_gro_adr.dump()

# integration and on the fly analysis
for s in range(1, intervals + 1):
    integrator.run(nsteps)
    step = nsteps * s
    time = step * timestep_short * timestepdevide_medium * timestepdevide_short
    #print "Step", step, "done."
    T = temperature.compute() * 120.27267 #* 1.5
    Eb = interBondedHH.computeEnergy() + interBondedOH.computeEnergy()
    EAng = interAngle.computeEnergy()
    ELj= interNB.computeEnergy()
    Ek = integrator.computeKineticEnergy() #0.5 * T * (3 * num_particles)
    EPI = integrator.computeRingEnergy()
    Etotal = Ek+Eb+EAng+ELj+EPI
    outfile.write(fmt%(step, time, T, Eb, EAng, EPI, ELj, Ek, Etotal))
    print (fmt%(step, time, T, Eb, EAng, EPI, ELj, Ek, Etotal))
    dump_conf_gro.dump()
    dump_conf_gro_adr.dump()

# Done, close file
outfile.close()

# Write the thermalized AdResS configuration
filenamexyz = "output.xyz"
espressopp.tools.writexyz(filenamexyz, system)
filenamepdb = "output.pdb"
espressopp.tools.pdbwrite(filenamepdb, system)


print "Successfully finished"
# Output format for screen and file
#print "i*timestep, T, P, Eb, EAng, ELj, EQQ, Ek, Ecorr, Etotal"
#fmt='%5.5f %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8f %15.8f\n'

# Timer and energy file
#start_time = time.clock()
#outfile = open("esp.dat", "w")

# Main MD loop
#for i in range(check):

    # Analysis
#    T = temperature.compute() * 120.27267 * 1.5
#    P = pressure.compute() * 16.6054
#    Eb = 0
#    EAng = 0
   # for bd in bondedinteractions.values(): Eb+=bd.computeEnergy()
    # for ang in angleinteractions.values(): EAng+=ang.computeEnergy()
    # ELj= ljinteraction.computeEnergy()
    # EQQ= qq_interactions.computeEnergy()
    # Ek = 0.5 * T * (2 * num_particles) / 120.27267
    # Ecorr = fec.computeCompEnergy()
    # Etotal = Ek+Eb+EAng+EQQ+ELj+Ecorr
    # dump_conf_gro.dump()
    # dump_conf_gro_adr.dump()
    # outfile.write(fmt%(i*steps/check*timestep, T, P, Eb, EAng, ELj, EQQ, Ek, Ecorr, Etotal))
    # print (fmt%(i*steps/check*timestep, T, P, Eb, EAng, ELj, EQQ, Ek, Ecorr, Etotal))


# simulation information
#end_time = time.clock()
#timers.show(integrator.getTimers(), precision=3)
#sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
#sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
#sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
#sys.stdout.write('Integration steps = %d\n' % integrator.step)
#sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

