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
! N^2 = B\pm F/2A.
module plasma
implicit none
integer, parameter :: nspecies=2,nomega=1000
integer :: ispecies,i,itheta,ntheta,iw,ipf
integer :: ntp,ntm
real, dimension(nspecies) :: nj,Aj,Zj
real :: omegace,omegamax,domce   ! over omegape.
real :: N2min,N2max
real, dimension(nomega) :: S,D,P,R,L,A,B,C,F2,resfac
real, dimension(nomega) ::  N2P,N2M,omega,X,Y,Nplus,Nminus
real :: theta,sintheta,costheta
character*8 string
character*3 sangle
character*1 :: key=' '
integer :: ikey=0
real :: izero

data Aj/1.,1836./
data Zj/-1.,1./
data nj/1.,1./
contains

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
  if(.false.)then
     resfac=(omegace-omega)*(omega*Aj(2)-omegace)
     N2P=(B*resfac+sqrt(F2*resfac**2))/(2.*A*resfac)
     N2M=(B*resfac-sqrt(F2*resfac**2))/(2.*A*resfac)
  endif
  if(.false.)then
! Direct Appleton-Hartree formula. Electrons only.
     X=1/omega**2   ;  Y=omegace/omega
     F2=(Y*sintheta)**4/4.+((1-X)*Y*costheta)**2
     F2=sqrt(F2)
     N2M=1-X*(1-X)/(1-X-(Y*sintheta)**2/2.+F2)
     N2P=1-X*(1-X)/(1-X-(Y*sintheta)**2/2.-F2)
  endif
! Explicit Sign flips at omega= omegac
  N2P=(B+sign(1.,(1.-omega/omegace)*(omega*Aj(2)/omegace-1))*sqrt(F2))/(2.*A)
  N2M=(B-sign(1.,(1.-omega/omegace)*(omega*Aj(2)/omegace-1))*sqrt(F2))/(2.*A)
  Nplus=sign(1.,N2P)*sqrt(abs(N2P))
  Nminus=sign(1.,N2M)*sqrt(abs(N2M))
end subroutine evalN2

subroutine initialize
  omega(1)=omegamax/nomega
  izero=nomega*alog(omega(1))/alog(omegamax/omega(1))
  do i=2,nomega
!     omega(i)=i*omega(1)   !linear
     omega(i)=omega(1)*exp(alog(omegamax/omega(1))*(i-1.)/(nomega-1.))
  enddo
end subroutine initialize

subroutine markatomega(om,frac,trace,yinc,mychar)
  real om,frac,yinc
  real, dimension(nomega) :: trace
  character*(*) mychar
  real, external :: wx2nx,wy2ny
!  i=int(nomega*om/omegamax)   !linear
  i=1+int(alog(om/omega(1))/alog(omegamax/omega(1))*(nomega-1.))
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
!  write(*,*)'iargc()=',iargc(),'key=',key
  open(12,file='coldplas.input',status='old',err=2)
  read(12,*,end=2,err=2)omegace,omegamax,N2min,N2max,ntheta
  write(*,'(a,4f10.3,i4)')'Read parameters',omegace,omegamax,N2min,N2max,ntheta
  read(12,*,end=2,err=2)ipf
2 continue
  close(12)
  do i=1,iargc()
     call getarg(i,key)
     ikey=ichar(key)
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
   write(*,*)'Usage coldplas [<key>]'
   write(*,*)'y-range up/k down/l; x-range left/s right/d' 
   write(*,*)'B-value e=dec r=inc; p: print, v: save settings,'
   write(*,*)'w: wait, n: no-wait plotting, else quit'
elseif(iw.eq.65362.or.iw.eq.ichar('k'))then
   N2max=2.*N2max
   N2min=2.*N2min
elseif(iw.eq.65364.or.iw.eq.ichar('l'))then
   N2max=N2max/2.
   N2min=N2min/2.
elseif(iw.eq.65361.or.iw.eq.ichar('s'))then
   omegamax=omegamax/2.
elseif(iw.eq.65363.or.iw.eq.ichar('d'))then
   omegamax=omegamax*2.
elseif(iw.eq.ichar('w'))then
   ipf=abs(ipf)
elseif(iw.eq.ichar('n'))then
   ipf=-3
elseif(iw.eq.ichar('c'))then
omegace=2.5
omegamax=4.001
N2min=-2.
N2max=8.
elseif(iw.eq.ichar('v'))then
   call savesettings
else
   iw=0
endif
end subroutine uif

end module plasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program coldplas
use plasma
!write(*,'(a,2f10.4)')'nj ',nj,'Aj ',Aj,'Zj ',Zj
!Defaults
omegace=2.5
domce=0.1
omegamax=4.001
N2min=-2.
N2max=8.
ntheta=10
ipf=0

call readsettings
call pfset(ipf)
! Start of user interface loop.
1 continue
call initialize
call axregion(.15,.9,.1,.7)
call pltinit(0.,omega(nomega),N2min,N2max)
call scalewn(omega(1),omega(nomega),N2min,N2max,.true.,.false.)
call charsize(.018,.018)
call axis()
call axlabels('!Aw!@/!Aw!@!dp!d','N')
call boxtitle('!Aq!@=0: Right-hand, Left-hand; !Aq!@=90: Ordinary, Extraordinary') 
call jdrwstr(0.08,wy2ny(0.9*N2min+0.1*N2max),'-|N!u2!u|!u1/2!u',0.)
call winset(.true.)
call vecw(1.,N2min,0)
call vecw(1.,N2max,1)
call vecw(omega(1),0.,0)
call vecw(omegamax,0.,1)
call jdrwstr(wx2nx(omegace),wy2ny(0.9*N2min+0.1*N2max),'!AW!@',0.)
call jdrwstr(wx2nx(omegace),wy2ny(0.9*N2min+0.1*N2max)-.022,'!A{!@',0.)
call jdrwstr(wx2nx(omegace),wy2ny(0.9*N2min+0.1*N2max)+.022,'!A}!@',0.)
if(omegace/Aj(2).ge.omega(2))then
   call jdrwstr(wx2nx(omegace/Aj(2)),wy2ny(0.9*N2min+0.1*N2max),'!AW!@!di!d',0.)
   call jdrwstr(wx2nx(omegace/Aj(2)),wy2ny(0.9*N2min+0.1*N2max)-.022,'!A{!@',0.)
   call jdrwstr(wx2nx(omegace/Aj(2)),wy2ny(0.9*N2min+0.1*N2max)+.022,'!A}!@',0.)
endif
do itheta=1,ntheta
   theta=3.1415926*(itheta-1.)/(ntheta-1.)/2.
   sintheta=sin(theta) ; costheta=cos(theta)
   thetadeg=180.*theta/3.1415926
   write(sangle,'(i2)')nint(thetadeg)
   call evalN2

   call color(15)
 ! Annotations
   if(abs(thetadeg-90.).lt.1)then
      if(omegace.gt.1.02)call markatomega(omegace,0.,Nplus,-0.014,'E')
      call markatomega(omegace*(1+sqrt(1+4./omegace**2))/2.+.1,&
           0.,Nplus,.01,'E')
      call markatomega(.98*sqrt(omegace**2+1.),1.,N2P,0.,'!Aw!@!duh!d')
      call markatomega(.98*sqrt(omegace**2+1.),1.,N2P,.02,'!A}!@')
      if(omegace.gt.1.02)call markatomega((2.*omegace+1.)/3.,0.,N2M,-.01,'O')
      call markatomega(.85,0.,N2P,0.01,'E')
!      call markatomega(.85,0.,N2M,0.01,'O')
   elseif(abs(thetadeg).lt.1)then
      if(omegace.gt.1.02)call markatomega((omegace+1.)/2.,-2.,Nplus,0.01,'R')
      call markatomega((2.*omegace+1.)/3.,0.,Nminus,.01,'L')
!      call markatomega(sqrt(omegace**2+1.),-2.,N2P,0.,'R')
      call markatomega(.85,0.,Nplus,-0.01,'L')
      call markatomega(.25,-1.,Nminus,-0.01,'R')
!      call markatomega(omegace*(1+sqrt(1+4./omegace**2))/2.,&
!           -.75,N2P,.013,'!Aw!@!dR!d!A{!@')
      call markatomega(omegace*(1+sqrt(1+4./omegace**2))/2.,&
           .75,N2P,.013,'!A{!@!Aw!@!dR!d')      
      call markatomega((-(1-Aj(1)/Aj(2))*omegace &
           +sqrt(omegace**2*(1+Aj(1)/Aj(2))**2+4.*(1+Aj(1)/Aj(2))))/2.,&
           -.75,N2P,.013,'!Aw!@!dL!d!A{!@') ! Left cut needs ion corrections
      call charangl(90.)
      if(omegace.ge.1.)then
         call markatomega(.5,3.,Nminus,0.,'!Aw!@!dp!dcos!Aq _!@')
         call markatomega(.3,1.,Nminus,0.,' Whistler Branch')
      else
         call markatomega(.5*omegace,1.,Nminus,0.,' Electron cyclotron waves')
      endif
      if(omega(1).lt.omegace/Aj(2))then
         call markatomega(0.6*omegace/Aj(2),1.2,Nplus,0.,'Ion cyclotron waves')
      endif
      if(3.*omega(1).lt.omegace/Aj(2))then
         call markatomega(3.*omega(1),-1.,Nminus,0.,'Alfven !A_!@')
         if(3.*omega(1).lt..3*omegace/Aj(2))call markatomega(3.*omega(1),1.5,Nminus,0.,'Shear Alfven')
      elseif(omega(1).lt.omegace/40.)then
         call markatomega(3.*omega(1),-1.,Nminus,0.,'Whistlers !A_!@')
      endif
      call winset(.false.)
!      call jdrwstr(wx2nx(-.15*omegamax),wy2ny(N2max),'Propagating',-1.3)
!      call jdrwstr(wx2nx(-.15*omegamax),wy2ny(0.),'Evanescent',-1.)
      call winset(.true.)
      call charangl(0.)
   endif
   ntp=nomega
   ntm=nomega
   do i=2,nomega ! Find the resonance gaps for plotting.
      if(ntp.eq.nomega.and.N2P(i-1).gt.0.and.N2P(i).lt.0.)ntp=i
      if(ntm.eq.nomega.and.N2M(i-1).gt.1.and.N2M(i).lt.1.)ntm=i
   enddo
   if(ntp.ne.nomega)Nplus(ntp-1)=N2max
   if(ntm.ne.nomega)Nminus(ntm-1)=N2max

   call color(1)
   if(.false.)then
      call polyline(omega,N2P,ntp-1)  ! plot skipping resonance
      call polyline(omega(ntp),N2P(ntp),nomega-ntp+1)
      call color(2)
      call polyline(omega,N2M,ntm-1)
      !   call fwrite(thetadeg,iw,0,string)
      !   call labeline(omega,N2M,ntm-1,string,iw-1)
      call polyline(omega(ntm),N2M(ntm),nomega-ntm+1)
   else
      call polyline(omega,Nplus,ntp-1)  ! plot skipping resonance
   if(ntp.ne.nomega)Nplus(ntp-1)=N2min
      call polyline(omega(ntp-1),Nplus(ntp-1),nomega-ntp+1)
      call color(2)
      call polyline(omega,Nminus,ntm-1)
   if(ntm.ne.nomega)Nminus(ntm)=N2min
      call polyline(omega(ntm),Nminus(ntm),nomega-ntm+1)
      nback=ntm/2
      if(Nminus(ntm-nback).lt.N2max*(1.-(100-thetadeg)*.005) &
           .and.ntm.lt.nomega)then
         call jdrwstr(wx2nx(omega(ntm-1)), &
              wy2ny(N2max)-.002*(100-thetadeg),sangle,0.)
         if(nint(thetadeg).eq.90)then
            call charangl(90.)
            call color(15)
            call jdrwstr(wx2nx(omega(ntm-1)), &
                 wy2ny(N2max)-.002*(100-thetadeg),'Lower Hybrid',-1.5)
            call charangl(0.)
         endif
      else
!         write(*,*)sangle,ntm,N2max,Nminus(ntm-nback)
      endif
      call color(15)
   endif
   
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
!write(*,*)'iw=',iw
if(iw.ne.0)goto 1

write(*,*)'  omegace, omegamax, N2min,   N2max,  ntheta      current values are'
write(*,'(4f10.4,i4)')omegace,omegamax,N2min,N2max,ntheta
!write(*,*)'Append a further integer parameter -3 to print plot without stopping.'
end program coldplas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

