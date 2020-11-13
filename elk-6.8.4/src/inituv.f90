
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine inituv
use modmain
use modbog
use modrandom
implicit none
! local variables
integer ik,ist,jst
real(8), parameter :: rndevv=0.01d0
! allocatable arrays
complex(8), allocatable :: evecu(:,:),evecv(:,:)
! initialise the Bogoliubov eigenvectors and write to file if required
if (task.eq.800) then
  allocate(evecu(nstsv,nstsv),evecv(nstsv,nstsv))
  evecu(:,:)=0.d0
  do ist=1,nstsv
    evecu(ist,ist)=1.d0
  end do
  do jst=1,nstsv
    do ist=1,nstsv
      evecv(ist,jst)=rndevv*cmplx(randomu()-0.5d0,randomu()-0.5d0,8)
    end do
  end do
  do ik=1,nkpt
    call putevecuv(ik,evecu,evecv)
  end do
  deallocate(evecu,evecv)
end if
! allocate the Bogoliubov eigenvalue array
if (allocated(evaluv)) deallocate(evaluv)
allocate(evaluv(nstsv,nkpt))
! allocate the occupation number array
if (allocated(occuv)) deallocate(occuv)
allocate(occuv(nstsv,nkpt))
if (task.eq.800) then
! initialise the occupation numbers (equal to those of the normal state)
  occuv(:,:)=occsv(:,:)
else
! read the occupation numbers from file
  do ik=1,nkpt
    call getoccuv(ik,occuv(:,ik))
  end do
end if
! allocate the U-norm array
if (allocated(unormk)) deallocate(unormk)
allocate(unormk(nstsv,nkpt))
end subroutine

