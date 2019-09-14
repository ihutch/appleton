!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cold plasma dispersion relation. Appleton-Hartree display.
! Frequency on logarithmic axis, plotting N rather than N^2
! Variables:
! Frequencies normalized to omegape. omega. omegace (positive) (=>B).
! Species of mass Aj (times electron mass), density nj (times electron den)
! Charge Zj. Electrons have charge -1. 
! S = 1 - \sum_j n_j/(Aj*omega^2-Zj^2*omegace^2/Aj)
! D = \sum_j Zj*omegace/(Aj*omega) n_j/(Aj*omega^2-omegace^2/Aj)
! P = 1 - \sum_j n_j/(Aj*omega)^2.
! R = 1 - \sum_j n_j /(Aj*omega*(omega+omegace*Zj/Aj) = S+D
! L = 1 - \sum_j n_j /(Aj*omega*(omega-omegace*Zj/Aj) = S-D
! Propagation angle theta
! A = S*\sin^2 \theta + P*\cos^2 \theta
! B = R*L*\sin^2 \theta + P*S*(1 + \cos^2 \theta )
! C = P*R*L
! F^2 = (RL - PS)^2 \sin^4 \theta + 4 P^2 D^2 \cos^2 \theta
! is the discriminant, and is manifestly positive. 
! Solution is
! N^2 = B\pm F/2A. But it is better evaluated as N^2=2C/(B\mp F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module plasma
implicit none
integer, parameter :: nspecies=2,nomega=2000
integer :: ispecies,i,itheta,ntheta,iw,ipf
integer :: ntp,ntm,ntp2
real, dimension(nspecies) :: nj,Aj,Zj
real :: omegace,omegamax,domce   ! over omegape.
real :: N2min,N2max
real, dimension(nomega) :: S,D,P,R,L,A,B,C,F2,resfac
real, dimension(nomega) ::  N2P,N2M,omega,X,Y,Nplus,Nminus
real, dimension(nomega) :: EyExP,EzExP,EyExM,EzExM 
real :: theta,sintheta,costheta
character*3 sangle

! Control parameters:
character*8 string
character*20 argument
integer :: ikey=0
real :: izero,xnlower
logical :: lEz=.false.,logx=.true.,logy=.true.

data Aj/1.,1836./
data Zj/-1.,1./
data nj/1.,1./
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine evalN2
  S=1. ; D=0.; P=1.
  do ispecies=1,nspecies
     S=S-nj(ispecies)/&
          (Aj(ispecies)*omega**2-Zj(ispecies)**2*omegace**2/Aj(ispecies))
     D=D- (Zj(ispecies)*omegace/(Aj(ispecies)*omega))* nj(ispecies)/&
          (Aj(ispecies)*omega**2-Zj(ispecies)**2*omegace**2/Aj(ispecies))
     P=P- nj(ispecies)/(Aj(ispecies)*omega)**2     
  enddo
  R=S+D
  L=S-D
  A=S*sintheta**2+P*costheta**2
  B=R*L*sintheta**2+P*S*(1+costheta**2)
  C=P*R*L
  F2=(R*L-P*S)**2*sintheta**4+4*P**2*D**2*costheta**2
  ! Explicit Sign flips at omega= omegac
  if(.true.)then ! This version works well at low omega.
     N2P=2.*C/(B-sign(1.,(1.-omega/omegace)*(omega*Aj(2)/omegace-1))*sqrt(F2))
     N2M=2.*C/(B+sign(1.,(1.-omega/omegace)*(omega*Aj(2)/omegace-1))*sqrt(F2))
  else  ! This does not work well at low omega.
  N2P=(B+sign(1.,(1.-omega/omegace)*(omega*Aj(2)/omegace-1))*sqrt(F2))/(2.*A)
  N2M=(B-sign(1.,(1.-omega/omegace)*(omega*Aj(2)/omegace-1))*sqrt(F2))/(2.*A)
  endif
  Nplus=sign(1.,N2P)*sqrt(abs(N2P))
  Nminus=sign(1.,N2M)*sqrt(abs(N2M))
! Polarizations Ey/iEx and Ez/Ex.
  EyExP=D/(N2P-S)
  EyExM=D/(N2M-S)
  EzExP=N2P*sintheta*costheta/(N2P*sintheta**2-S)
  EzExM=N2M*sintheta*costheta/(N2M*sintheta**2-S)
end subroutine evalN2

subroutine initialize
  if(logx)then
     omega(1)=0.2*omegamax/nomega
  else
     omega(1)=omegamax/nomega
  endif
  izero=nomega*alog(omega(1))/alog(omegamax/omega(1))
  do i=2,nomega
     if(logx)then
        omega(i)=omega(1)*exp(alog(omegamax/omega(1))*(i-1.)/(nomega-1.))
     else
        omega(i)=i*omega(1)   !linear
     endif
  enddo
end subroutine initialize

subroutine markatomega(om,frac,trace,yinc,mychar)
  real om,frac,yinc
  real, dimension(nomega) :: trace
  character*(*) mychar
  real, external :: wx2nx,wy2ny
  if(logx)then  !This switch does not work without consistency
     i=1+int(alog(om/omega(1))/alog(omegamax/omega(1))*(nomega-1.))
  else
     i=int(om/omega(1))
  endif
  if(i.le.nomega.and.i.gt.0)call jdrwstr(wx2nx(omega(i)),wy2ny(trace(i))+yinc,mychar,frac)
end subroutine markatomega

subroutine savesettings
  open(13,file='coldplas.input',status='unknown')
  write(13,'(4f10.5,i4)')omegace,omegamax,N2min,N2max,ntheta
  write(13,*)ipf
  write(13,*)' #omegace,omegamax,N2min,N2max,ntheta, // ipf'
  write(13,*)' key    k,l y-range; s,d x-range; r,e B-value'
  close(13)
end subroutine savesettings

subroutine readsettings
  open(12,file='coldplas.input',status='old',err=2)
  read(12,*,end=2,err=2)omegace,omegamax,N2min,N2max,ntheta
  write(*,'(a,4f10.3,i4)')'Read parameters',omegace,omegamax,N2min,N2max,ntheta
  read(12,*,end=2,err=2)ipf
2 continue
  close(12)
  do i=1,iargc()
     call getarg(i,argument)
     ikey=ichar(argument(1:1))
     call uif(ikey)                    ! Perform possible key action
  enddo
  if(ipf.ne.0)write(*,*)'Nonstop plotting',ipf
  if(ikey.ne.0)call savesettings    ! If key was a valid command
end subroutine readsettings

subroutine uif(iw)
! User interface routine
  integer :: iw
if(iw.eq.ichar('r'))then
   if(domce.lt..034*omegace)domce=domce*10.
   omegace=omegace+domce
elseif(iw.eq.ichar('e'))then
   if(domce.gt..32*omegace)domce=domce/10.
   omegace=omegace-domce
elseif(iw.eq.ichar('p'))then
   call pfset(3)
elseif(iw.eq.ichar('h'))then
   write(*,*)'Usage coldplas [<key> ... ]'
   write(*,*)'Keyboard graph-control keys as follows'
   write(*,*)'y-range up/k down/l; x-range left/s right/d' 
   write(*,*)'B-value e=dec r=inc; p: print, v: save settings,'
   write(*,*)'Toggle polarization from Ey/Ex to Ez/Ex: z'
   write(*,*)'Toggle logarithmic axes: x, y'
   write(*,*)'w: wait, n: no-wait plotting, else quit'
elseif(iw.eq.65362.or.iw.eq.ichar('k'))then
   N2max=2.*N2max
   N2min=2.*N2min
elseif(iw.eq.65364.or.iw.eq.ichar('l'))then
   N2max=N2max/2.
   N2min=N2min/2.
elseif(iw.eq.65361.or.iw.eq.ichar('s'))then
   if(logx)then
      if(omegamax.gt.0.02*omegace)omegamax=omegamax/2.
   else
      if(omegamax.gt.0.0001*omegace)omegamax=omegamax/2.
   endif
elseif(iw.eq.65363.or.iw.eq.ichar('d'))then
   if(logx)then
      if(omegamax.lt.10.*omegace)omegamax=omegamax*2.
   else
      if(omegamax.lt.3.*omegace)omegamax=omegamax*2.
   endif
elseif(iw.eq.ichar('w'))then
   ipf=abs(ipf)
elseif(iw.eq.ichar('n'))then
   ipf=-3
elseif(iw.eq.ichar('z'))then
   lEz=.not.lEz
elseif(iw.eq.ichar('x'))then
   logx=.not.logx
elseif(iw.eq.ichar('y'))then
   if(.not.logy)then; N2min=.1;   N2max=2000.
   else; N2min=-1.; N2max=10.
   endif
   logy=.not.logy
elseif(iw.eq.ichar('c'))then
omegace=2.5
omegamax=5.0001
N2min=1.e-1
N2max=8000.
elseif(iw.eq.ichar('v'))then
   call savesettings
else
   iw=0
endif
end subroutine uif

subroutine annotations(thetadeg)
  implicit none
  real :: thetadeg,ynorm,oalfven
  real, external :: wx2nx,wy2ny
  ! Annotations
  call color(15)
  if(Nminus(ntm-40).lt.N2max*(1.-(100-thetadeg)*.002) &
        .and.ntm.lt.nomega)then
      if(nint(thetadeg).eq.90)then
         call jdrwstr(wx2nx(omega(ntm-1)), &
              wy2ny(N2max)-.002*(100-thetadeg),'!Aq!@='//sangle,0.)
         call charangl(90.)
         if(omega(ntm-1).lt.omegamax.and.wx2nx(omega(ntm-1)).gt.xnlower) &
              call jdrwstr(wx2nx(omega(ntm-1)), &
              wy2ny(N2max)-.002*(100-thetadeg),'Lower Hybrid',-1.5)
         call charangl(0.)
      else
         call jdrwstr(wx2nx(omega(ntm-1)), &
              wy2ny(N2max)-.002*(100-thetadeg),sangle,0.)
      endif
   endif
   if(nint(thetadeg).eq.90)then
      if(omegace.gt.1.02)call markatomega(1.03*omegace,0.,Nplus,-0.014,'E')
      call markatomega(omegace*(1+sqrt(1+4./omegace**2))/2.+.1,&
           0.,Nplus,.01,'E')
      call markatomega(.98*sqrt(omegace**2+1.),1.,N2P,0.0,'!Aw!@!duh!d')
      call markatomega(.98*sqrt(omegace**2+1.),1.,N2P,.02,'!A}!@')
      if(omegace.gt.1.02)call markatomega((2.*omegace+1.)/3.5,-0.1,N2M,-.005,'O')
      call markatomega(.85,0.,N2P,0.01,'E')
!      call markatomega(.85,0.,N2M,0.01,'O')
      call vecw(1.,N2min,0)
      call vecw(1.,N2max,1)
      call vecw(omega(1),0.,0)
      call vecw(omegamax,0.,1)
      ynorm=wy2ny(N2min)+.035
      call jdrwstr(wx2nx(omegace),ynorm,'!AW!@',0.)
      call jdrwstr(wx2nx(omegace),ynorm-.022,'!A{!@',0.)
      call jdrwstr(wx2nx(omegace),ynorm+.022,'!A}!@',0.)
      if(omegace/Aj(2).ge.omega(2))then
         call jdrwstr(wx2nx(omegace/Aj(2)),ynorm,'!AW!@!di!d',0.)
         call jdrwstr(wx2nx(omegace/Aj(2)),ynorm-.022,'!A{!@',0.)
         call jdrwstr(wx2nx(omegace/Aj(2)),ynorm+.022,'!A}!@',0.)
      endif
   elseif(abs(thetadeg).lt.1)then
      if(omegace.gt.1.02)call markatomega((omegace+1.3)/2.,-1.,Nplus,0.01,'R')
      call markatomega(1.+.1*omegace,-0.02,Nminus,-.002,'L')
      if(wx2nx(omegace/Aj(2)).gt.xnlower) &
           call markatomega(0.7*omegace/Aj(2),1.2,Nplus,0.,'L')
!      call markatomega(sqrt(omegace**2+1.),-2.,N2P,0.,'R')
      call markatomega(.85,0.,Nplus,-0.01,'L')
      call markatomega(.25,-1.,Nminus,-0.01,'R')
      call markatomega(omegace*(1+sqrt(1+4./omegace**2))/2.,&
           .75,N2P,.033,'!A{!@!Aw!@!dR!d')      
      call markatomega((-(1-Aj(1)/Aj(2))*omegace &
           +sqrt(omegace**2*(1+Aj(1)/Aj(2))**2+4.*(1+Aj(1)/Aj(2))))/2.,&
           -.75,N2P,.033,'!Aw!@!dL!d!A{!@') ! Left cut needs ion corrections
      call charangl(90.)
      if(omegace.ge.1.)then
         if(wy2ny(N2max).gt.wy2ny(2.)+.33) &
              call markatomega(.6,2.,Nminus,0.,'!Aw!@!dp!dcos!Aq!B{!@')
         call markatomega(.2,1.,Nminus,0.,' Helicon or')
         call markatomega(.3,1.,Nminus,0.,' Whistler Branch')
      else
         call markatomega(.5*omegace,1.,Nminus,0.,' Electron cyclotron waves')
      endif
      if(wx2nx(omegace/Aj(2)).gt.xnlower)then
         call markatomega(0.8*omegace/Aj(2),1.2,Nplus,0.,'Ion cyclotron waves')
      endif
      oalfven=omega(1)+.05*omegace/Aj(2)
      if(wx2nx(oalfven).gt.xnlower)then
         call markatomega(oalfven,-1.,Nminus,0.,'Alfven !A_!@')
         call markatomega(oalfven,1.5,Nminus,0.,'Shear Alfven')
      endif
      if(wx2nx(omegace/250).gt.xnlower)then
         call markatomega(omegace/250.,-1.,Nminus,0.,'Whistlers !A_!@')
      endif
      call winset(.false.)
      call winset(.true.)
      call charangl(0.)
   endif
 end subroutine annotations

 subroutine findresonances
   implicit none
! Ensure plotting of resonances
   ntm=nomega
   ntp=nomega
   ntp2=nomega
   do i=3,nomega ! Find the resonance gaps for plotting.
      if(ntp.eq.nomega.and.N2P(i-1).gt.0.and.N2P(i).lt.0.)ntp=i
      if(ntp.ne.nomega.and.N2P(i-1).gt.0.and.N2P(i).lt.0.)ntp2=i
      if(ntm.eq.nomega.and.N2M(i-1).gt.1.and.N2M(i).lt.1.)then
!         write(*,*)'A',A(i),' B',B(i),' C',C(i),' F2',F2(i)
!         write(*,*)'S',S(i),' P',P(i),' D',D(i),' o',omega(i)
!         write(*,*)i,ntm,N2M(i-1),N2M(i)
         ntm=i
      endif
   enddo
   if(ntp.ne.nomega)then
      Nplus(ntp-1)=N2max
      if(.not.logy)Nplus(ntp)=N2min
   endif
   if(ntp2.ne.nomega)then
      Nplus(ntp2-1)=N2max
      if(.not.logy)Nplus(ntp2)=N2min
   endif
   if(ntm.ne.nomega)then
      Nminus(ntm-1)=N2max
      if(.not.logy)Nminus(ntm)=N2min
   endif
 end subroutine findresonances

 subroutine plotthesolutions(EyExmax,EzExmax)
   real :: EyExmax,EzExmax
   if(lEz)then
      call polycolorgapline(omega,Nplus,nomega, &
           int(1+min(.99,abs(EzExP)/EzExmax)*240.),Nplus.ge.0.)   
      if(.not.logy)call polycolorgapline(omega,Nplus,nomega, &
           int(1+min(.99,abs(EzExP)/EzExmax)*240.),Nplus.le.0.5)   
      call polycolorgapline(omega,Nminus,nomega, &
           int(1+min(.99,abs(EzExM)/EzExmax)*240.),Nminus.ge.0)
      if(.not.logy)call polycolorgapline(omega,Nminus,nomega, &
           int(1+min(.99,abs(EzExM)/EzExmax)*240.),Nminus.le.0.5)
   else
      call polycolorgapline(omega,Nplus,nomega, &
           int(1+min(.99,abs(EyExP)/EyExmax)*240.),Nplus.ge.0.)   
      if(.not.logy)call polycolorgapline(omega,Nplus,nomega, &
           int(1+min(.99,abs(EyExP)/EyExmax)*240.),Nplus.le.0.5)   
      call polycolorgapline(omega,Nminus,nomega, &
           int(1+min(.99,abs(EyExM)/EyExmax)*240.),Nminus.ge.0)
      if(.not.logy)call polycolorgapline(omega,Nminus,nomega, &
           int(1+min(.99,abs(EyExM)/EyExmax)*240.),Nminus.le.0.5)
   endif
 end subroutine plotthesolutions
 
end module plasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Has to be after the plasma module in this file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program coldplas
use plasma

!Defaults
omegace=2.5
domce=0.1
omegamax=5.0001
ipf=0
EyExmax=2.8
EzExmax=2.8
ntheta=10

! Handle command line arguments and settings.
call readsettings
call pfset(ipf)
call brgwscaled(0.,.9)
if(.not.logx)then; omegamax=3.5; else; omegamax=5.0001; endif
if(logy)then
   N2min=1.e-1
   N2max=2000.
else
   N2min=-1.
   N2max=10.
endif

xnlower=.17
! Start of user interface loop.
1 continue
call initialize
call axregion(.15,.88,.1,.7)
call pltinit(0.,omega(nomega),N2min,N2max)
call charsize(.018,.018)
call gradlegend(0.,EyExmax,1.04,0.,1.04,1.,.01,.true.)
if(lEz)then
   call legendline(1.045,0.45,258,'|E!dz!d/E!dx!d|')
else
   call legendline(1.045,0.45,258,'|E!dy!d/E!dx!d|')
endif
call scalewn(omega(1),omega(nomega),N2min,N2max,logx,logy)
call axis()
call axis2
call axlabels('!Aw!@/!Aw!@!dp!d','N')
if(N2min.lt.0.)call jdrwstr(0.08,wy2ny(0.97*N2min+0.03*N2max), &
     '-|N!u2!u|!u1/2!u',0.)
call boxtitle('!Aq!@=0: !p!n-!n!qRight-hand, !p!n-!n!qLeft-hand; !Aq!@=90: !p!n-!n!qOrdinary, !p!n-!n!qExtraordinary') 
call winset(.true.)
do itheta=1,ntheta
   theta=3.1415926*(itheta-1.)/(ntheta-1.)/2.
   sintheta=sin(theta) ; costheta=cos(theta)
   thetadeg=180.*theta/3.1415926
   write(sangle,'(i2)')nint(thetadeg)
   call evalN2
   call findresonances
   call plotthesolutions(EyExmax,EzExmax)
   call color(15)
   call annotations(thetadeg)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call prtend('')
if(ipf.lt.0)then
   call savesettings
   call exit()
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! User interface.
call pfset(0)
call eye3d(iw)
call uif(iw)
if(iw.ne.0)goto 1

write(*,*)'  omegace, omegamax, N2min,   N2max,  ntheta      current values are'
write(*,'(4f10.4,i4)')omegace,omegamax,N2min,N2max,ntheta
!write(*,*)'Append a further integer parameter -3 to print plot without stopping.'
end program coldplas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
