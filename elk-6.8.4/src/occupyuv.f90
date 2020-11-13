
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine occupyuv
use modmain
use modbog
use modmpi
implicit none
integer ik,ist
real(8) e,chg,x
real(8) t1,t2,t3
! external functions
real(8), external :: stheta
t1=1.d0/swidth
unorm=0.d0
chg=0.d0
do ik=1,nkpt
  do ist=1,nstsv
    e=evaluv(ist,ik)
    x=-e*t1
    t2=stheta(stype,x)
! Bogoliubov occupation number
    occuv(ist,ik)=occmax*t2
! U-norm for current state and k-point
    t3=unormk(ist,ik)
! add to the average U-norm
    unorm=unorm+wkpt(ik)*t3
! add electron and hole contributions to the total charge
    t2=occmax*(t2*t3+(1.d0-t2)*(1.d0-t3))
    chg=chg+wkpt(ik)*t2
  end do
end do
unorm=unorm/dble(nstsv)
! adjust the Fermi energy
efermi=efermi+tauefm*(chgval-chg)
if (mp_mpi) then
  if (abs(chg-chgval).gt.epschg) then
    write(*,*)
    write(*,'("Warning(occupyuv): incorrect charge : ",2G18.10)') chg,chgval
  end if
end if
end subroutine

