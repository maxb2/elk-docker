
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvbmatk(vmt,vir,bmt,bir,ngp,igpig,wfmt,ld,wfgp,vbmat)
use modmain
use modomp
implicit none
! arguments
! the potential and field are multiplied by the radial integration weights in
! the muffin-tin and by the characteristic function in the interstitial region
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(8), intent(in) :: wfgp(ld,nspinor,nstsv)
complex(8), intent(out) :: vbmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,npc,igp,nthd
! automatic arrays
complex(8) wfmt2(npcmtmax,nspinor),z(ngkmax)
! allocatable arrays
complex(8), allocatable :: wfir1(:,:),wfir2(:,:)
! external functions
complex(8), external :: zdotc
! zero the matrix elements
vbmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt2,ias,is,npc,ist) &
!$OMP NUM_THREADS(nthd)
do jst=1,nstsv
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
! apply local potential and magnetic field to spinor wavefunction
    if (ncmag) then
! non-collinear case
      call vbmk1(npc,vmt(:,ias),bmt(:,ias,1),bmt(:,ias,2),bmt(:,ias,3), &
       wfmt(:,ias,1,jst),wfmt(:,ias,2,jst),wfmt2,wfmt2(:,2))
    else
! collinear case
      call vbmk2(npc,vmt(:,ias),bmt(:,ias,1),wfmt(:,ias,1,jst), &
       wfmt(:,ias,2,jst),wfmt2,wfmt2(:,2))
    end if
! compute the inner products
    do ist=1,jst
      vbmat(ist,jst)=vbmat(ist,jst) &
       +zdotc(npc,wfmt(:,ias,1,ist),1,wfmt2,1) &
       +zdotc(npc,wfmt(:,ias,2,ist),1,wfmt2(:,2),1)
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
!---------------------------!
!     interstitial part     !
!---------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,wfir2,z) &
!$OMP PRIVATE(ispn,jspn,igp,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfir1(ngtot,nspinor),wfir2(ngtot,nspinor))
!$OMP DO
do jst=1,nstsv
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    wfir1(:,ispn)=0.d0
    do igp=1,ngp(jspn)
      wfir1(igfft(igpig(igp,jspn)),ispn)=wfgp(igp,ispn,jst)
    end do
    call zfftifc(3,ngridg,1,wfir1(:,ispn))
  end do
! apply local potential and magnetic field to spinor wavefunction
  if (ncmag) then
! non-collinear case
    call vbmk1(ngtot,vir,bir,bir(:,2),bir(:,3),wfir1,wfir1(:,2),wfir2, &
     wfir2(:,2))
  else
! collinear case
    call vbmk2(ngtot,vir,bir,wfir1,wfir1(:,2),wfir2,wfir2(:,2))
  end if
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir2(:,ispn))
    do igp=1,ngp(jspn)
      z(igp)=wfir2(igfft(igpig(igp,jspn)),ispn)
    end do
    do ist=1,jst
      vbmat(ist,jst)=vbmat(ist,jst)+zdotc(ngp(jspn),wfgp(:,ispn,ist),1,z,1)
    end do
  end do
end do
!$OMP END DO
deallocate(wfir1,wfir2)
!$OMP END PARALLEL
call freethd(nthd)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vbmat(ist,jst)=conjg(vbmat(jst,ist))
  end do
end do
return

contains

pure subroutine vbmk1(n,v,b1,b2,b3,wf11,wf12,wf21,wf22)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: v(n),b1(n),b2(n),b3(n)
complex(8), intent(in) :: wf11(n),wf12(n)
complex(8), intent(out) :: wf21(n),wf22(n)
! local variables
integer i
complex(8) z1
do i=1,n
  z1=cmplx(b1(i),b2(i),8)
  wf21(i)=(v(i)+b3(i))*wf11(i)+conjg(z1)*wf12(i)
  wf22(i)=(v(i)-b3(i))*wf12(i)+z1*wf11(i)
end do
end subroutine

pure subroutine vbmk2(n,v,b,wf11,wf12,wf21,wf22)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: v(n),b(n)
complex(8), intent(in) :: wf11(n),wf12(n)
complex(8), intent(out) :: wf21(n),wf22(n)
! local variables
integer i
do i=1,n
  wf21(i)=(v(i)+b(i))*wf11(i)
  wf22(i)=(v(i)-b(i))*wf12(i)
end do
end subroutine

end subroutine

