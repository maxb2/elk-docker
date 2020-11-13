
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradzvcln(is,gzfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: is
complex(8), intent(out) :: gzfmt(npmtmax,3)
! local variables
integer nr,nri,ir,i
! automatic arrays
complex(8) zvclmt(npmtmax)
nr=nrmt(is)
nri=nrmti(is)
! convert nuclear Coulomb potential to complex spherical harmonics expansion
zvclmt(1:npmt(is))=0.d0
i=1
do ir=1,nri
  zvclmt(i)=vcln(ir,is)
  i=i+lmmaxi
end do
do ir=nri+1,nr
  zvclmt(i)=vcln(ir,is)
  i=i+lmmaxo
end do
! compute the gradient of the potential
call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),zvclmt,npmtmax,gzfmt)
end subroutine

