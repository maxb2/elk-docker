
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevaluv
use modmain
use modbog
implicit none
! local variables
integer ik,ist
! write out the valence eigenvalues
open(50,file='EIGVALUV.OUT',form='FORMATTED')
write(50,'(I6," : nkpt")') nkpt
write(50,'(I6," : nstsv")') nstsv
do ik=1,nkpt
  write(50,*)
  write(50,'(I6,3G18.10," : k-point, vkl")') ik,vkl(:,ik)
  write(50,'(" (state, eigenvalue, occupancy, U-norm below)")')
  do ist=1,nstsv
    write(50,'(I6,3G18.10)') ist,evaluv(ist,ik),occuv(ist,ik),unormk(ist,ik)
  end do
  write(50,*)
end do
close(50)
end subroutine

