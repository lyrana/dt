# Particles_Mesh

__version__ = 0.1
__author__ = 'T. P. Hughes'
__all__ = ['Particles_Mesh.create_initial_particles', 
           'Particles_Mesh.particle_species', ]

import ParticlesM as PS_M
import dolfin as df

class Particles_Mesh(object):
    """Particles interacting with a mesh.
    """

    def __init__(self, allSpecies):
        """ Process the list of particle species given by the user."""
        self.particle_species = PS_M.Particles(allSpecies)
        self.setup()

    def setup(self):
        x = 1


    def add_particle_species(self, species_name, initial_distribution):
# push the particle species onto a list
# use xylose storage
        y = 1

    def create_initial_particles(self, pmesh):
        """Puts particles on the given mesh according to a specified density and velocity distribution"""
        z = 1        
        
# loop on the mesh

        for c in df.cells(pmesh):
# Print the vertices of the mesh
#            for v in df.vertices(c):
#                x = v.x(0)
#                y = v.x(1)
#                print x, y
#                pass
            vol = c.volume()
            p = c.midpoint()
            print " volume = ", vol

# Loop on species to create initial particles
            for sp_name in self.particle_species.names:
                dist_funct = self.particle_species.initial_distribution[sp_name]
                num_per_cell = self.particle_species.number_per_cell[sp_name]
#                print "dist_funct = ", dist_funct
                # compute the density at the center of the cell
                if dist_funct != None:
                    n, v, T = dist_funct(p)
                    print "n = ", n, "v = ", v, "T = ", T
                    # Compute number of real particle in the cell
                    number_of_real_particles = n*vol
                    # Compute weight
                    weight = number_of_real_particles/num_per_cell
                    print "weight for ", sp_name, " is ", weight

# lay down species in this cell, according to the specified density

                    xp, yp, 
