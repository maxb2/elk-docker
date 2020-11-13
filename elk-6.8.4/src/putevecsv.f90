
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevecsv(fext,ik,evecsv)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer recl
character(256) fname
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,evecsv
fname=trim(scrpath)//'EVECSV'//trim(fext)
!$OMP CRITICAL(u126)
open(126,file=trim(fname),form='UNFORMATTED',access='DIRECT',action='WRITE', &
 recl=recl)
write(126,rec=ik) vkl(:,ik),nstsv,evecsv
close(126)
!$OMP END CRITICAL(u126)
end subroutine

