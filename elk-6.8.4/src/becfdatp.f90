
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine becfdatp
use modmain
use modphonon
use modtddft
implicit none
integer its,is,ia
! generate the time step grid
call gentimes
! write zero force to file for all time steps
open(50,file='FORCETOT_TD.OUT',form='FORMATTED')
write(50,'(I8," : number of time steps")') ntimes
do its=1,ntimes
  write(50,'(I8,G18.10)') its,times(its)
  do is=1,nspecies
    do ia=1,natoms(is)
      write(50,'(2I4,3G18.10)') is,ia,0.d0,0.d0,0.d0
    end do
  end do
end do
close(50)
! write the atomic displacements and velocities
datposc(:,:,:,:)=0.d0
datposc(ipph,1,iaph,isph)=deltaph/tstime
call writedatposc
end subroutine

