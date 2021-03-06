
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine energykncr
use modmain
implicit none
integer ist,is,ias
integer nr,nri,ir,i
! automatic arrays
real(8) rfmt(npmtmax)
! external functions
real(8), external :: rfmtinp
! calculate the kinetic energy for core states
engykncr=0.d0
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! sum of core eigenvalues
  do ist=1,nstsp(is)
    if (spcore(ist,is)) engykncr=engykncr+occcr(ist,ias)*evalcr(ist,ias)
  end do
! core density
  rfmt(1:npmt(is))=0.d0
  i=1
  if (spincore) then
! spin-polarised core
    do ir=1,nri
      rfmt(i)=rhocr(ir,ias,1)+rhocr(ir,ias,2)
      i=i+lmmaxi
    end do
    do ir=nri+1,nr
      rfmt(i)=rhocr(ir,ias,1)+rhocr(ir,ias,2)
      i=i+lmmaxo
    end do
  else
! spin-unpolarised core
    do ir=1,nri
      rfmt(i)=rhocr(ir,ias,1)
      i=i+lmmaxi
    end do
    do ir=nri+1,nr
      rfmt(i)=rhocr(ir,ias,1)
      i=i+lmmaxo
    end do
  end if
  engykncr=engykncr-rfmtinp(nr,nri,wrmt(:,is),rfmt,vsmt(:,ias))
end do
end subroutine

