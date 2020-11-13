
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevecfv(fext,ik,evecfv)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
! local variables
integer recl
character(256) fname
! find the record length
inquire(iolength=recl) vkl(:,ik),nmatmax,nstfv,nspnfv,evecfv
fname=trim(scrpath)//'EVECFV'//trim(fext)
!$OMP CRITICAL(u122)
open(122,file=trim(fname),form='UNFORMATTED',access='DIRECT',action='WRITE', &
 recl=recl)
write(122,rec=ik) vkl(:,ik),nmatmax,nstfv,nspnfv,evecfv
close(122)
!$OMP END CRITICAL(u122)
end subroutine

