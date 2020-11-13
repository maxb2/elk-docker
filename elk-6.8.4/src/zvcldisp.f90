
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zvcldisp(zvclmt)
use modmain
implicit none
! arguments
complex(8), intent(inout) :: zvclmt(npmtmax,natmtot)
! local variables
integer is,ia,ias
integer nr,nri,np,i,j
real(8) t1,t2
! automatic arrays
complex(8) gzfmt1(npmtmax,3),gzfmt2(npmtmax,3)
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the gradient of the nuclear and core Coulomb potential
    call grzvclncr(ias,gzfmt1)
    do i=1,3
      t1=-datposc(i,0,ia,is)
! add the gradient to the Coulomb potential for first-order term
      zvclmt(1:np,ias)=zvclmt(1:np,ias)+t1*gzfmt1(1:np,i)
! compute the gradient of the gradient for second-order term
      call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),gzfmt1(:,i),npmtmax, &
       gzfmt2)
      do j=1,3
        t2=-0.5d0*t1*datposc(j,0,ia,is)
        zvclmt(1:np,ias)=zvclmt(1:np,ias)+t2*gzfmt2(1:np,j)
      end do
    end do
  end do
end do
end subroutine

