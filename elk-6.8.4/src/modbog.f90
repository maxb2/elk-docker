
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modbog

! Bogoliubov equation eigenvalues
real(8), allocatable :: evaluv(:,:)
! U-norm for each state and k-point
real(8), allocatable :: unormk(:,:)
! average U-norm over all states and k-points
real(8) unorm
! Bogoliubov equation occupation numbers
real(8), allocatable :: occuv(:,:)
! Fermi energy adjustment step size
real(8) tauefm
! Fermi energy convergence tolerance
real(8) epsefm
! Hartree-Fock-Bogoliubov coupling constant
real(8) ehfb

end module

