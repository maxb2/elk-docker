
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zftzf(ngp,jlgpr,ylmgp,ld,sfacgp,zfmt,zfir,zfgp)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: jlgpr(njcmax,nspecies,ngp)
complex(8), intent(in) :: ylmgp(lmmaxo,ngp)
integer, intent(in) :: ld
complex(8), intent(in) :: sfacgp(ld,natmtot)
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtc)
complex(8), intent(out) :: zfgp(ngp)
! local variables
integer is,ia,ias,ig
integer nrc,nrci,irco,irc
integer l,m,lm,i,j
real(8) t0,y0
complex(8) zsm,z1,z2
! automatic arrays
complex(8) ylm(lmmaxo),zfft(ngtc)
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
! multiply by coarse characteristic function
zfft(:)=zfir(:)*cfrc(:)
! Fourier transform to coarse G-grid
call zfftifc(3,ngdgc,-1,zfft)
zfgp(1:ngp)=zfft(igfc(1:ngp))
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
t0=fourpi/omega
y0=t0*y00
do ig=1,ngp
  lm=1
  do l=1,lmaxo
    z1=t0*zilc(l)
    do m=-l,l
      lm=lm+1
      ylm(lm)=z1*ylmgp(lm,ig)
    end do
  end do
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    irco=nrci+1
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      zsm=0.d0
      i=1
      j=1
! inner part of muffin-tin
      if (lmaxi.eq.1) then
        do irc=1,nrci
          zsm=zsm+wrcmt(irc,is)*(jlgpr(j,is,ig)*zfmt(i,ias)*y0+jlgpr(j+1,is,ig)&
           *(zfmt(i+1,ias)*ylm(2)+zfmt(i+2,ias)*ylm(3)+zfmt(i+3,ias)*ylm(4)))
          i=i+4
          j=j+2
        end do
      else
        do irc=1,nrci
          z1=jlgpr(j,is,ig)*zfmt(i,ias)*y0
          i=i+1
          j=j+1
          lm=1
          do l=1,lmaxi
            lm=lm+1
            z2=zfmt(i,ias)*ylm(lm)
            i=i+1
            do m=1-l,l
              lm=lm+1
              z2=z2+zfmt(i,ias)*ylm(lm)
              i=i+1
            end do
            z1=z1+jlgpr(j,is,ig)*z2
            j=j+1
          end do
          zsm=zsm+wrcmt(irc,is)*z1
        end do
      end if
! outer part of muffin-tin
      do irc=irco,nrc
        z1=jlgpr(j,is,ig)*zfmt(i,ias)*y0
        i=i+1
        j=j+1
        lm=1
        do l=1,lmaxo
          lm=lm+1
          z2=zfmt(i,ias)*ylm(lm)
          i=i+1
          do m=1-l,l
            lm=lm+1
            z2=z2+zfmt(i,ias)*ylm(lm)
            i=i+1
          end do
          z1=z1+jlgpr(j,is,ig)*z2
          j=j+1
        end do
        zsm=zsm+wrcmt(irc,is)*z1
      end do
      zfgp(ig)=zfgp(ig)+conjg(sfacgp(ig,ias))*zsm
    end do
  end do
end do
end subroutine

