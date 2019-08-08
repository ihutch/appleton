! Cold plasma dispersion relation. Appleton-Hartree display.
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
integer, parameter :: nspecies=2,nomega=400
integer :: ispecies,i,itheta,ntheta,iw,ipf
integer :: ntp,ntm
real, dimension(nspecies) :: nj,Aj,Zj
real :: omegace,omegamax   ! over omegape.
real :: N2min,N2max
real, dimension(nomega) :: S,D,P,R,L,A,B,C,F2,resfac
real, dimension(nomega) ::  N2P,N2M,omega,X,Y
real :: theta,sintheta,costheta
character*8 string

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

end subroutine evalN2

subroutine initialize
  do i=1,nomega
     omega(i)=i*omegamax/nomega
  enddo
end subroutine initialize

subroutine markatomega(om,frac,trace,yinc,mychar)
  real om,frac,yinc
  real, dimension(nomega) :: trace
  character*(*) mychar
  real, external :: wx2nx,wy2ny
  i=int(nomega*om/omegamax)
  if(i.le.nomega)call jdrwstr(wx2nx(omega(i)),wy2ny(trace(i))+yinc,mychar,frac)
end subroutine markatomega

end module plasma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program coldplas
use plasma
!write(*,'(a,2f10.4)')'nj ',nj,'Aj ',Aj,'Zj ',Zj
!Defaults
omegace=2.5
domce=0.1
omegamax=3.001
N2min=-2.
N2max=8.
ntheta=10
ipf=0

open(12,file='coldplas.input',status='old',err=2)
read(12,*,end=2,err=2)omegace,omegamax,N2min,N2max,ntheta
write(*,'(a,4f10.3,i4)')'Read parameters',omegace,omegamax,N2min,N2max,ntheta
read(12,*,end=2,err=2)ipf
write(*,*)'Nonstop plotting',ipf
2 continue

if(ipf.eq.0)then
   write(*,*)'up/k: increase, down/l: decrease y-range; right/s incr, left/d decr x-range' 
   write(*,*)'e: decrease B, r: increase B, p: print, else quit'
else
   call pfset(ipf)
endif
! Start of user interface loop.
1 continue
call initialize
call axregion(.15,.9,.1,.7)
call pltinit(0.,omega(nomega),N2min,N2max)
call charsize(.018,.018)
call axis()
call axlabels('!Aw!@/!Aw!@!dp!d','N!u2!u')
call boxtitle('!Aq!@=0: Right-hand, Left-hand; !Aq!@=90: Ordinary, Extraordinary') 
call winset(.true.)
call vecw(1.,-10.,0)
call vecw(1.,100.,1)
call vecw(0.,0.,0)
call vecw(omegamax,0.,1)
call jdrwstr(wx2nx(omegace),wy2ny(0.9*N2min+0.1*N2max),'!AW!@',0.)
call jdrwstr(wx2nx(omegace),wy2ny(0.9*N2min+0.1*N2max)-.022,'!A{!@',0.)
call jdrwstr(wx2nx(omegace),wy2ny(0.9*N2min+0.1*N2max)+.022,'!A}!@',0.)
if(omegace.ge..02*Aj(2)*omegamax)then
call jdrwstr(wx2nx(omegace/Aj(2)),wy2ny(0.9*N2min+0.1*N2max),'!AW!@!di!d',0.)
call jdrwstr(wx2nx(omegace/Aj(2)),wy2ny(0.9*N2min+0.1*N2max)-.022,'!A{!@',0.)
call jdrwstr(wx2nx(omegace/Aj(2)),wy2ny(0.9*N2min+0.1*N2max)+.022,'!A}!@',0.)
endif
do itheta=1,ntheta
   theta=3.14159*(itheta-1.)/(ntheta-1.)/2.
   sintheta=sin(theta) ; costheta=cos(theta)
   thetadeg=180.*theta/3.1415926
   call evalN2

   call color(15)
! Annotations
   if(abs(thetadeg-90.).lt.1)then
      if(omegace.gt.1.02)call markatomega(omegace,0.,N2P,-0.014,'E')
      call markatomega(omegace*(1+sqrt(1+4./omegace**2))/2.+.1,&
           0.,N2P,.01,'E')
      call markatomega(.98*sqrt(omegace**2+1.),1.,N2P,0.,'!Aw!@!duh!d')
      call markatomega(.98*sqrt(omegace**2+1.),1.,N2P,.02,'!A}!@')
      if(omegace.gt.1.02)call markatomega((2.*omegace+1.)/3.,0.,N2M,-.01,'O')
      call markatomega(.85,0.,N2P,0.01,'E')
      call markatomega(.85,0.,N2M,0.01,'O')
   elseif(abs(thetadeg).lt.1)then
      if(omegace.gt.1.02)call markatomega((omegace+1.)/2.,-2.,N2P,0.01,'R')
      call markatomega((2.*omegace+1.)/3.,0.,N2M,.01,'L')
      call markatomega(sqrt(omegace**2+1.),-2.,N2P,0.,'R')
      call markatomega(.85,0.,N2P,-0.01,'L')
      call markatomega(.25,-1.,N2M,-0.01,'R')
      call markatomega(omegace*(1+sqrt(1+4./omegace**2))/2.,&
           -.75,N2P,.013,'!Aw!@!dR!d!A{!@')
      call markatomega(omegace*(-1+sqrt(1+4./omegace**2))/2.,&
           -.75,N2P,.013,'!Aw!@!dL!d!A{!@')
      call charangl(90.)
      if(omegace.ge.1)call markatomega(.8,3.,N2M,0.,'!Aw!@!dp!dcos!Aq _!@')
      if(omegace.ge.1)call markatomega(.6,1.5,N2M,0.,'Whistler Branch')
      if(omegamax.le.4.)call markatomega(.05*omegamax,-1.,N2M,0.,'Whistlers !A_!@')
      call winset(.false.)
      call jdrwstr(wx2nx(-.15*omegamax),wy2ny(N2max),'Propagating',-1.3)
      call jdrwstr(wx2nx(-.15*omegamax),wy2ny(0.),'Evanescent',-1.)
      call winset(.true.)
      call charangl(0.)
   endif
   do i=2,nomega ! Find the resonance gaps for plotting.
      if(N2P(i-1).gt.0.and.N2P(i).lt.0.)ntp=i
      if(N2M(i-1).gt.0.and.N2M(i).lt.0.)ntm=i
   enddo
   !write(*,'(4f8.4)')(omega(i),S(i),N2P(i),N2M(i),i=1,nomega)
   call color(1)
   call polyline(omega,N2P,ntp-1)  ! plot skipping resonance
   call polyline(omega(ntp),N2P(ntp),nomega-ntp+1)
   call color(2)
   call polyline(omega,N2M,ntm-1)
!   call fwrite(thetadeg,iw,0,string)
!   call labeline(omega,N2M,ntm-1,string,iw-1)
   call polyline(omega(ntm),N2M(ntm),nomega-ntm+1)
   call color(15)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call prtend('')
if(ipf.ne.0)call exit()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! User interface.
call pfset(0)
call eye3d(iw)
if(iw.eq.ichar('r'))then
   if(domce.lt..034*omegace)domce=domce*10.
   omegace=omegace+domce
   goto 1
elseif(iw.eq.ichar('e'))then
   if(domce.gt..32*omegace)domce=domce/10.
   omegace=omegace-domce
   goto 1
elseif(iw.eq.ichar('p'))then
   call pfset(3)
   goto 1
elseif(iw.eq.65362.or.iw.eq.ichar('k'))then
   N2max=2.*N2max
   N2min=2.*N2min
   goto 1
elseif(iw.eq.65364.or.iw.eq.ichar('l'))then
   N2max=N2max/2.
   N2min=N2min/2.
   goto 1
elseif(iw.eq.65361.or.iw.eq.ichar('s'))then
   omegamax=omegamax/2.
   goto 1
elseif(iw.eq.65363.or.iw.eq.ichar('d'))then
   omegamax=omegamax*2.
   goto 1
endif

write(*,*)'To change defaults create file coldplas.input containing values'
write(*,*)'omegace,omegamax,N2min,N2max,ntheta  whose current values are'
!write(*,*)omegace,omegamax,N2min,N2max,ntheta
write(*,'(4f10.4,i4)')omegace,omegamax,N2min,N2max,ntheta
write(*,*)'Append a further integer parameter -3 to print plot without stopping.'
end program coldplas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

