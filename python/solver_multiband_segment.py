################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


from pytriqs.applications.impurity_solvers.cthyb_segment import *
from pytriqs.applications.dft.U_matrix import *
import pytriqs.utility.mpi as mpi
from types import *
import numpy
from solver_multiband import set_U_matrix

#########################################
#
#  Solver for the Multi-Band problem
#
#########################################

class SolverMultiBandSegment (Solver):
    """
          This is a solver for a multiband local Hamiltonian and segment picture.
              Calling arguments:
              beta = inverse temperature
              n_orb = Number of local orbitals
              U_interact = Average Coulomb interaction
              J_hund     = Hund coupling
              use_matrix: Use the interaction matrix calculated from the Slater integrals
              is use_matrix, you need also:
                  l: angular momentum of the orbital, l=2 is d
                  T: Transformation matrix for U vertex. If not present, use standard complex harmonics
              block_names should contain only blocs with dim. 1
              w_max = maximum frequency for Green functions

              dynamical_U = if dynamical interaction, you need to provide K(tau) later

    """

    def __init__(self, beta, n_orb, U_interact=None, J_hund=None, block_names=False, map=False,
                 use_matrix = True, l=2, T=None, dimreps=None, irep=None, deg_orbs = [], sl_int = None, w_max=200):

        # initialize Segment solver, parameters are to be redefined later
        Solver.__init__(self,
                    beta = beta,
                    block_names = block_names,
                    n_tau_delta= 10000,
                    n_tau_g = 10000,
                    n_w = int(w_max*beta/(2.*numpy.pi)),
                    #N_Frequencies_Accumulated = int(Omega_Max*Beta/(10.*numpy.pi)),
                    #Fitting_Frequency_Start = int(Omega_Max*Beta/(20.*numpy.pi)),
                    )
        print "use_matrix =",use_matrix
        #set up vertex
        U, Up, U4ind, offset = set_U_matrix(U_interact,J_hund,n_orb,l,use_matrix,T,sl_int,use_spinflip=False,dim_reps=dimreps,irep=irep)

        Umat = numpy.zeros([n_orb*2,n_orb*2])

        #if (numpy.shape(U)[0] == 2 * l +1 ):
        if (numpy.shape(U)[0] == n_orb):
            #size of U is 2*l+1, set from U, Up.
            Umat[0:n_orb,0:n_orb]=U
            Umat[n_orb:2*n_orb,n_orb:2*n_orb]=U
            Umat[0:n_orb,n_orb:2*n_orb]=Up
            Umat[n_orb:2*n_orb,0:n_orb]=Up

        else:

            Umat = U

        self.U = Umat 

    def solve(self, dynamical_U=False):
       Solver.solve( n_cycles = 10000,
            length_cycle = 100,
            n_warmup_cycles = 1000,
            dynamical_U = dynamical_U,
            U=self.U
         )
