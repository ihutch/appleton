!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cold plasma dispersion relation. Appleton-Hartree display.
! Frequency on logarithmic axis, plotting N rather than N^2
! for various angles of propagation theta.
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Alternative display of dispersion relation
! N_\perp^2 for various given N_\parallel^2
! Variables:
! Parallel refractive index squared N2a, Perp squared N2e
! A = S
! B = (S-N2a)*(S+P)-D^2
! C = P*[(S-N2a)^2 -D^2]
! F^2 = B^2 -4AC
! is the discriminant, and may be negative, giving complex N2e
! Solution is
! N2e^2 = (B\pm F)/2A.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Non-uniform plasma plot with given Npar uses the alternative display
! but the omegaf is fixed and the omegape varies. We arrive at the same
! calculation by substituting omega=omegaf/omegape, and
! omegace=omegace/omegape. Ensure that no element of omegape is zero!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module coldplasmawaves
implicit none
integer, parameter :: nspecies=2,nomega=2000
integer :: ispecies,i,itheta,ntheta,iw,ipf,iNpar,nNpar=10
integer :: ntp,ntm,ntp2
real, dimension(nspecies) :: nj,Aj,Zj
real :: omegace,omegamax,domce   ! over omegape.
real :: omegfac=0.2              ! Factor for log omega minumum.
real :: N2min,N2max,EyExmax,EzExmax,thetadeg
real, dimension(nomega) :: S,D,P,R,L,A,B,C,F2,resfac
real, dimension(nomega) ::  N2P,N2M,omega,X,Y,Nplus,Nminus
real, dimension(nomega) :: EyExP,EzExP,EyExM,EzExM 
real :: theta,sintheta,costheta
character*3 sangle
real, external :: wy2ny,wx2nx
! Coldperp and Nomegape extra variables
logical, dimension(nomega) :: logic
real, dimension(nomega) :: omegape=1.,omegapi
logical :: lperp=.false.
integer :: iplottype=1
real :: N2a,Npmax=5,Na,olh,ouh
integer :: colorfast=2,colorslow=1
complex, dimension(nomega) :: CN2P,CN2M
real :: omegaf   ! The fixed frequency for omegape plots, unnormalized
! Control parameters:
character*8 string
character*20 argument
integer :: ikey=0
real :: izero,xnlower=.17
logical :: lEz=.false.,logx=.true.,logy=.true.,lplotomega=.false.

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
! Polarizations Ey/iEx and Ez/Ex.
  EyExP=D/(N2P-S)
  EyExM=D/(N2M-S)
  EzExP=N2P*sintheta*costheta/(N2P*sintheta**2-S)
  EzExM=N2M*sintheta*costheta/(N2M*sintheta**2-S)
  if(lplotomega)then
     N2P=N2P*omega**2
     N2M=N2M*omega**2
  endif
  Nplus=sign(1.,N2P)*sqrt(abs(N2P))
  Nminus=sign(1.,N2M)*sqrt(abs(N2M))
end subroutine evalN2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine evalNperp2
  ! omega, and omegace here are normalized to omegap
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
  A=S
  B=(S-N2a)*(S+P)-D**2
  C=P*((S-N2a)**2-D**2)
  F2=B**2-4.*A*C
  resfac=1.
  CN2P=(B*resfac+sqrt(F2*resfac**2))/(2.*A*resfac)
  CN2M=(B*resfac-sqrt(F2*resfac**2))/(2.*A*resfac)
  where(F2.lt.0)
     N2P=0.
     N2M=0.
  elsewhere(real(CN2P).gt.real(CN2m))
     N2P=real(CN2P)
     N2M=real(CN2M)
  elsewhere
     N2P=real(CN2M)
     N2M=real(CN2P)
  endwhere
end subroutine evalNperp2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine evalNperpop
  ! omegape, and omegace here are unnormalized; so should be omegaf.
  S=1. ; D=0.; P=1.
  do ispecies=1,nspecies
     S=S-nj(ispecies)/ &
          (Aj(ispecies)*(omegaf/omegape)**2 &
             -Zj(ispecies)**2*(omegace/omegape)**2/Aj(ispecies))
     D=D- (Zj(ispecies)*omegace/(Aj(ispecies)*omegaf))* nj(ispecies)/&
          (Aj(ispecies)*(omegaf/omegape)**2 &
             -Zj(ispecies)**2*(omegace/omegape)**2/Aj(ispecies))
     P=P- nj(ispecies)/(Aj(ispecies)*(omegaf/omegape))**2     
  enddo
  R=S+D
  L=S-D
  A=S
  B=(S-N2a)*(S+P)-D**2
  C=P*((S-N2a)**2-D**2)
  F2=B**2-4.*A*C
  resfac=1.
  CN2P=(B*resfac+sqrt(F2*resfac**2))/(2.*A*resfac)
  CN2M=(B*resfac-sqrt(F2*resfac**2))/(2.*A*resfac)
  where(F2.lt.0)
     N2P=0.
     N2M=0.
  elsewhere(real(CN2P).gt.real(CN2m))
     N2P=real(CN2P)
     N2M=real(CN2M)
  elsewhere
     N2P=real(CN2M)
     N2M=real(CN2P)
  endwhere
end subroutine evalNperpop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Various utilities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize ! Only the scope and omega arrays.
  if(logx)then
     omega(1)=omegfac*omegamax/nomega
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initomegape
  ! When doing omegape plots, interpret omegamax as the
  ! maximum of omegapi/omegaf but omegapi is not normalized.
  omegaf=sqrt(Aj(1)/Aj(2))
  if(logx)then
     omegapi(1)=omegfac*omegamax*omegaf/nomega
  else
     omegapi(1)=omegamax*omegaf/nomega
  endif
  do i=2,nomega
     if(logx)then
        omegapi(i)=omegapi(1)  &
             *exp(alog(omegamax*omegaf/omegapi(1))*(i-1.)/(nomega-1.))
     else
        omegapi(i)=i*omegapi(1)   !linear
     endif
  enddo
  ! The calculation is done in terms of omegape.
  omegape=omegapi*sqrt(Aj(2)/Aj(1))                  ! Unnormalized.
end subroutine initomegape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine markatomega(om,frac,trace,yinc,mychar)
  real om,frac,yinc
  real, dimension(nomega) :: trace
  character*(*) mychar
  real, external :: wx2nx,wy2ny
  i=iatomega(om)
  if(i.le.nomega.and.i.gt.0)call jdrwstr(wx2nx(omega(i)),wy2ny(trace(i))+yinc,mychar,frac)
end subroutine markatomega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function iatomega(om)
  real om
  if(logx)then  !This switch does not work without consistency
     iatomega=1+int(alog(om/omega(1))/alog(omegamax/omega(1))*(nomega-1.))
  else
     iatomega=int(om/omega(1))
  endif  
end function iatomega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine savesettings
  open(13,file='coldplas.input',status='unknown')
  write(13,'(4f10.5,i4)')omegace,omegamax,N2min,N2max,ntheta
  write(13,*)ipf
  write(13,*)' #omegace,omegamax,N2min,N2max,ntheta, // ipf'
  write(13,*)' key    k,l y-range; s,d x-range; r,e B-value'
  close(13)
end subroutine savesettings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readsettings
  open(12,file='coldplas.input',status='old',err=2)
  read(12,*,end=2,err=2)omegace,omegamax,N2min,N2max,ntheta
  write(*,'(a,4f10.3,i4)')'Read parameters',omegace,omegamax,N2min,N2max,ntheta
  read(12,*,end=2,err=2)ipf
2 continue
  close(12)
! Parse commandline arguments, single characters passed to uif.
  do i=1,iargc()
     call getarg(i,argument)
     ikey=ichar(argument(1:1))
     call uif(ikey)                    ! Perform possible key action
  enddo
  if(ipf.ne.0)write(*,*)'Nonstop plotting',ipf
  if(ikey.ne.0)call savesettings    ! If key was a valid command
! Decide default plot ranges.
  if(.not.logx)then; omegamax=3.5;
  else; omegamax=5.0001 ; if(iplottype.eq.3)omegamax=20.00001
  endif
  if(logy)then; N2min=1.e-1; N2max=2000. ; if(iplottype.eq.3)N2max=200000
  else; N2min=-1. ; N2max=10. ;
  endif
end subroutine readsettings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine findresonances
   implicit none
! Ensure correct plotting of resonances by adjusting at gaps.
   ntm=nomega
   do i=3,nomega ! Find the resonance gaps for plotting.
      if(ntm.eq.nomega.and.N2M(i-1).gt.1.and.N2M(i).lt.1.)ntm=i
      if(Nminus(i-1).gt.0.and.Nminus(i).lt.0.)Nminus(i-1)=N2max
      if(Nminus(i-1).lt.N2min.and.Nminus(i).gt.N2min)Nminus(i)=N2min
      if(Nplus(i-1).gt.0.and.Nplus(i).lt.0)Nplus(i-1)=N2max
      if(Nplus(i-1).lt.N2min.and.Nplus(i).gt.N2min)Nplus(i)=N2min
   enddo
 end subroutine findresonances
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! User interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine uif(iw)
! User interface routine
  integer :: iw
call pfset(0) ! By default each new plot is not printed?
if(iw.eq.ichar('r'))then
   if(domce.lt..034*omegace)domce=domce*5.
   omegace=omegace+domce
elseif(iw.eq.ichar('e'))then
   if(domce.gt..32*omegace)domce=domce/5.
   omegace=omegace-domce
elseif(iw.eq.ichar('p'))then
   call pfset(3)
elseif(iw.eq.ichar('h'))then
   write(*,*)'Usage coldplas [<key> ... ]'
   write(*,*)'Keyboard and command line graph-control keys as follows'
   write(*,*)'Y-RANGE SHIFT: up/k down/l; X-RANGE: left/s right/d'
   write(*,*)'Y-RANGE EXPAND:u contract:i; X-RANGE (log) expnd:n cntrct:m d' 
   write(*,*)'B-VALUE e=decrs r=incrs; PRINT: p, SAVE settings: v,'
   write(*,*)'POLARIZATION toggle from Ey/Ex to Ez/Ex: z'
   write(*,*)'LOGARITHMIC axes toggle: x, y.',' NPERP, N plotting toggle: -'
   write(*,*)'ORDINATE plotting toggle from N to k: o'
   write(*,*)'WAIT at end of plotting: w, NO-WAIT: a.',' TELL parameters: t'
   write(*,*)'QUIT: q, return, or single left click in window.'
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
elseif(iw.eq.ichar('n'))then
   omegfac=omegfac/2.
elseif(iw.eq.ichar('m'))then
   omegfac=omegfac*2.
elseif(iw.eq.ichar('u'))then
   N2min=N2min/2.
elseif(iw.eq.ichar('i'))then
   N2min=N2min*2.
elseif(iw.eq.ichar('w'))then
   ipf=abs(ipf)
elseif(iw.eq.ichar('a'))then
   ipf=-3
elseif(iw.eq.ichar('z'))then
   lEz=.not.lEz
elseif(iw.eq.ichar('x'))then
   logx=.not.logx
elseif(iw.eq.ichar('o'))then
   lplotomega=.not.lplotomega
elseif(iw.eq.ichar('y'))then
   logy=.not.logy
   if(logy)then
      N2min=.1;   N2max=2000.; if(iplottype.eq.3)N2max=200000
   else
      if(lperp)then; N2max=600.;N2min=-60.
      else; N2min=-1.; N2max=10.;endif
   endif
elseif(iw.eq.ichar('c'))then
   omegace=2.5
   omegamax=5.0001; if(iplottype.eq.3)omegamax=20.00001
   N2min=1.e-1
   N2max=8000.
elseif(iw.eq.ichar('v'))then
   call savesettings
elseif(iw.eq.ichar('-').or.iw.eq.189)then   
   lperp=.not.lperp
   iplottype=mod(iplottype,3)+1
   if(iplottype.eq.3)omegamax=20.00001
elseif(iw.eq.ichar('q').or.iw.eq.65293)then
   iw=0
elseif(iw.eq.ichar('t'))then
   write(*,*)'  omegace, omegamax, N2min,   N2max,  ntheta  current values'
   write(*,'(4f10.4,i4)')omegace,omegamax,N2min,N2max,ntheta

else
   write(*,*)'key number=',iw,' not recognized'
endif
!   write(*,*)'iw number=',iw
end subroutine uif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plotting and Annotating 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine annotationN2
  implicit none
  real :: ynorm,oalfven
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
         if(logx.or.omegamax.lt.4.)then
            call charangl(70.)
            call markatomega(.1,1.,Nminus,0.,' Helicon or')
            call markatomega(.16,1.,Nminus,0.,' Whistler Branch')
            call charangl(90.)
         endif
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
!         call markatomega(omegace/250.,-1.,Nminus,0.,'Whistlers !A_!@')
         call charangl(0.)
         if(lplotomega)then
            call markatomega(omegace/200.,1.,Nminus,-0.01,'Whistlers')
         else
            call markatomega(omegace/200.,-1.,Nminus,-0.01,'Whistlers')
         endif
      endif
      call charangl(0.)
   endif
 end subroutine annotationN2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine annotationperp
   external wy2ny,wx2nx
   real wy2ny,wx2nx,ynorm
   ynorm=wy2ny(N2min)+.035
   ouh=sqrt(omegace**2+1.)
   call jdrwstr(wx2nx(ouh),ynorm,'!Aw!@!duh!d',1.2)
   call jdrwstr(wx2nx((-(1-Aj(1)/Aj(2))*omegace &
        +sqrt(omegace**2*(1+Aj(1)/Aj(2))**2+4.*(1+Aj(1)/Aj(2))))/2.),&
        ynorm,'!Aw!@!dL!d!A{!@',-.75)
   call jdrwstr(wx2nx(ouh),ynorm,'!A}!@',0.)
   call jdrwstr(wx2nx(omegace),ynorm,'!AW!@',-1.5)
   call jdrwstr(wx2nx(omegace),ynorm-0.012,'!A{!@',0.)
   call jdrwstr(wx2nx(omegace),ynorm+0.012,'!A}!@',0.)
   if(omegace/Aj(2).ge.omega(2))then
      call jdrwstr(wx2nx(omegace/Aj(2)),ynorm,'!AW!@!di!d',0.)
      call jdrwstr(wx2nx(omegace/Aj(2)),ynorm-.022,'!A{!@',0.)
      call jdrwstr(wx2nx(omegace/Aj(2)),ynorm+.022,'!A}!@',0.)
   endif
   ! Lower hybrid frequency for single ion species.
   olh=((1.+Aj(1)/Aj(2))+omegace**2*(1+(Aj(1)/Aj(2))**2))/2.
   olh=olh-sqrt((((1.-Aj(1)/Aj(2))+omegace**2*(1-(Aj(1)/Aj(2))**2))/2.)**2&
        +1./Aj(2))
   olh=sqrt(olh)
   call jdrwstr(wx2nx(olh),ynorm,'!A}!@',0)
   call drcstr('!Aw!@!dLH!d')
   call color(colorslow)
   call polyline((/ouh,ouh/),(/N2min,N2max/),2)
   call polyline((/olh,olh/),(/N2min,N2max/),2)
   call color(15)
 end subroutine annotationperp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine annotateNomegape
   external wy2ny,wx2nx
   real wy2ny,wx2nx,ynorm
   real ogc,oin,opeeqom
   real Rde,Rdi
   ynorm=wy2ny(N2min)+.035
   ogc=omegace/omegaf*sqrt(Aj(1)/Aj(2))
   oin=omegace/omegaf*Aj(1)/Aj(2)
   opeeqom=sqrt(Aj(1)/Aj(2))
   ! Lower hybrid frequency for single ion species. Solving for omegapi
   ! at given omega. omegapi^2=Rde*Rdi/(Rde+Aj(2)/Aj(1)*Rdi) where
   ! Rde,i = (omega^2-Omegae,i^2). Then one must normalize by omegaf.
   Rde=omegaf**2-omegace**2
   Rdi=omegaf**2-(omegace*Aj(1)/Aj(2))**2
   olh=Rde*Rdi/(Rde+Rdi*Aj(2)/Aj(1))
   olh=sqrt(olh)/omegaf ! This becomes nan and suppresses plotting when needed
   call color(colorslow)
   call polyline((/olh,olh/),(/N2min,N2max/),2)
   call color(15)
!   write(*,*)ogc,omegace,omegaf,olh
   call jdrwstr(wx2nx(ogc),ynorm,'!A}!@',0)
   call drcstr('(!AW!@!di!d!AW!@!de!d)!u1/2!u/!Aw!@')
   call jdrwstr(wx2nx(oin),ynorm,'!A}!@',0)
   call drcstr('!AW!@!di!d/!Aw!@')
   call jdrwstr(wx2nx(opeeqom),ynorm+.035,'!A{!@',0)
   call jdrwstr(wx2nx(opeeqom),ynorm+.035,'!Aw!@!dpe!d/!Aw!@=1',-1.1)
   call jdrwstr(wx2nx(olh),ynorm+.07,'!A{!@',0)
   call jdrwstr(wx2nx(olh),ynorm+.07,'!Aw!@!dLH!d/!Aw!@=1',-1.1)

  
 end subroutine annotateNomegape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 subroutine setupplots
   call axregion(.15,.88,.1,.7)
   call pltinit(0.,omega(nomega),N2min,N2max)
   call color(15)
   call charsize(.018,.018)
   if(.not.lperp)then
      call gradlegend(0.,EyExmax,1.04,0.,1.04,1.,.01,.true.)
      if(lEz)then
         call legendline(1.045,0.45,258,'|E!dz!d/E!dx!d|')
      else
         call legendline(1.045,0.45,258,'|E!dy!d/E!dx!d|')
      endif
   endif
   call scalewn(omega(1),omega(nomega),N2min,N2max,logx,logy)
   call axis()
   call axis2
   if(lplotomega)then
      if(lperp)then
         call axlabels('!Aw!@/!Aw!@!dp!d','k!d!A`!@!dc/!Aw!@!dp!d')
      else
         call axlabels('!Aw!@/!Aw!@!dp!d','kc/!Aw!@!dp!d')
      endif
      if(N2min.lt.0.)then
         call jdrwstr(0.08,wy2ny(0.97*N2min+0.03*N2max),'-|k!u2!u|!u1/2!u',0.)
         call jdrwstr(0.08,wy2ny(0.97*N2min+0.03*N2max)-.01,'_______',0.)
         call jdrwstr(0.08,wy2ny(0.97*N2min+0.03*N2max)-.033,'!Aw!@!dp!d/c',0.)
      endif
   else
      if(lperp)then
         call axlabels('!Aw!@/!Aw!@!dp!d','N!d!a`!@!d')
         if(N2min.lt.0.)call jdrwstr(0.15,wy2ny(0.)-.035, &
           '-|N!d!A`!@!d!u2!u|!u1/2!u',-1.1)               
      else
         call axlabels('!Aw!@/!Aw!@!dp!d','N')
         if(N2min.lt.0.)call jdrwstr(0.15,wy2ny(0.97*N2min+0.03*N2max), &
           '-|N!u2!u|!u1/2!u',0.-1.1)
      endif
   endif
 end subroutine setupplots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 subroutine setupomegape
   call axregion(.15,.88,.1,.7)
   call pltinit(0.,omegapi(nomega),N2min,N2max)
   call color(15)
   call charsize(.018,.018)
   call scalewn(omegapi(1),omegapi(nomega),N2min,N2max,logx,logy)
   call axis()
   call axis2
   call axlabels('!Aw!@!dpi!d/!Aw!@','N!d!a`!@!d')
   if(N2min.lt.0.)call jdrwstr(0.15,wy2ny(0.)-.035, &
        '-|N!d!A`!@!d!u2!u|!u1/2!u',-1.1)
 end subroutine setupomegape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine plotNs
  call setupplots
  call boxtitle('!Aq!@=0: !p!n-!n!qRight-hand, !p!n-!n!qLeft-hand;'// &
      '!Aq!@=90: !p!n-!n!qOrdinary, !p!n-!n!qExtraordinary')
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
     call annotationN2 ! Must be after findresonances
  enddo
  call prtend('')
end subroutine plotNs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotNperps
  real dxa,dya,xleg
  call setupplots
  call boxtitle('N!d!A`!@!d for N!d!A|!@!d given by line labels')
  call winset(.true.)
  call color(15)
  if(.not.logy)call polyline((/omega(1),omega(nomega)/),(/1.,1./),2)
  call annotationperp
  do iNpar=1,nNpar
     Na=(Npmax*(iNpar)/(nNpar))
     N2a=Na**2
     thetadeg=45.
     call evalNperp2
     call color(15)
     call fwrite(Na,iw,1,string)
     logic=.true.
     where(F2.le.0. .and. cshift(F2,1).gt.0.)
        N2P=0.5*(cshift(N2P,1)+cshift(N2M,1))
        N2M=N2P
     elsewhere(F2.le.0. .and. cshift(F2,-1).gt.0.)
        N2P=0.5*(cshift(N2P,-1)+cshift(N2M,-1))
        N2M=N2P
     elsewhere(F2.le.0. &
          .or. cshift(N2M,-1)*N2M.lt.-10. &
          .or. cshift(N2M,1)*N2M.lt.-10. &
          )
        logic=.false.
     end where
     call color(colorslow)
     call polygapline(omega,N2P,nomega,logic)
     i=(iatomega(olh)+nomega/5-(nomega/16)*mod(1+iNpar,2))
     if(Na.le.1)i=(iatomega(olh)-(nomega/20))
     i=max(min(i,nomega-1),1)
     dxa=wx2nx(omega(i+1))-wx2nx(omega(i))
     dya=wy2ny(N2P(i+1))-wy2ny(N2P(i))
     call charangl(atan2(dya,dxa)*180./3.1415926)
     call jdrwstr(wx2nx(omega(i)),wy2ny(N2P(i)),string,0.5)
     call charangl(0.)
     xleg=0.02
     if(.not.logx.or.olh/omega(1).lt.15.)xleg=0.7
     call legendline(xleg,.95,0,'Slow Wave')

     call color(colorfast)
     call dashset(colorfast)
     call legendline(xleg,.9,0,'Fast Wave')
     call polygapline(omega,N2M,nomega,logic)
     i=iatomega(olh)+nomega/20
     i=max(min(i,nomega-1),1)
     dxa=wx2nx(omega(i+1))-wx2nx(omega(i))
     dya=wy2ny(N2M(i+1))-wy2ny(N2M(i))
     call charangl(atan2(dya,dxa)*180./3.1415926)
     call jdrwstr(wx2nx(omega(i)),wy2ny(N2M(i)),string,0.5)
     call charangl(0.)
     call dashset(0)
  enddo
  call prtend('')
end subroutine plotNperps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotNomegape
  real dxa,dya,xleg
  integer iwidth
!  real vminus,vat,vplus
  omegapi=omegapi/omegaf  ! Scale omegapi for plotting.
  call setupomegape
  call fwrite(omegace/omegaf*sqrt(Aj(1)/Aj(2)),iwidth,3,string)
  call boxtitle('N!d!A`!@!d for N!d!A|!@!d given by line labels.' &
   //'  (!AW!@!di!d!AW!@!de!d)!u1/2!u/!Aw!@='//string(1:iwidth))
  call winset(.true.)
  call color(15)
  call annotateNomegape
  do iNpar=1,nNpar
     Na=(Npmax*(iNpar)/(nNpar))
     N2a=Na**2
     call evalNperpop
     call color(15)
     call fwrite(Na,iw,1,string)
     logic=.true.
     where(F2.le.0. .and. cshift(F2,1).gt.0.)
        N2P=0.5*(cshift(N2P,1)+cshift(N2M,1))
        N2M=N2P
     elsewhere(F2.le.0. .and. cshift(F2,-1).gt.0.)
        N2P=0.5*(cshift(N2P,-1)+cshift(N2M,-1))
        N2M=N2P
     elsewhere(F2.le.0. &
          .or. cshift(N2M,-1)*N2M.lt.-10. &
          .or. cshift(N2M,1)*N2M.lt.-10. &
          )
        logic=.false.
     end where

     xleg=0.02
     !     if(.not.logx.or.olh/omegapi(1).lt.15.)xleg=0.7
     call color(colorfast)
     call dashset(colorfast)
     call legendline(xleg,.9,0,'Fast Wave')
     call polygapline(omegapi,N2M,nomega,logic)
!     i=iatomega(olh)+nomega/20
!     i=max(min(i,nomega-1),1)
     i=nomega/3
     dxa=wx2nx(omegapi(i+1))-wx2nx(omegapi(i))
     dya=wy2ny(N2M(i+1))-wy2ny(N2M(i))
     call charangl(atan2(dya,dxa)*180./3.1415926)
     call jdrwstr(wx2nx(omegapi(i)),wy2ny(N2M(i)),string,0.5)
     call charangl(0.)
     call dashset(0)

! Attempts to improve slow wave plot don't work.     
!     where(N2P.gt.1.e5.and.cshift(N2P,1).lt.1.); N2P=100.; endwhere
!     where(N2P.lt.1..and.cshift(N2P,-1).gt.1.e5); N2P=cshift(N2P,1); endwhere
!     where(N2P.lt.1.e-5.and.cshift(N2P,1).gt.1) ; logic=.false.; endwhere
     call color(colorslow)
     call polygapline(omegapi,N2P,nomega,logic)
     i=nomega/2
     if(mod(iNpar,2).eq.0)i=nomega/3
     dxa=wx2nx(omegapi(i+1))-wx2nx(omegapi(i))
     dya=wy2ny(N2P(i+1))-wy2ny(N2P(i))
     call charangl(atan2(dya,dxa)*180./3.1415926)
     call jdrwstr(wx2nx(omegapi(i)),wy2ny(N2P(i)),string,0.5)
     call charangl(0.)
     call legendline(xleg,.95,0,'Slow Wave')
  enddo
  call prtend('')
  
end subroutine plotNomegape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module coldplasmawaves
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Has to be after the plasma module in this file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program coldplas
use coldplasmawaves

!Defaults
omegace=2.5
domce=0.1
ipf=0
EyExmax=2.8
EzExmax=2.8
ntheta=10
Npmax=5

! Handle command line arguments and settings.
call readsettings
! Plotting settings initialization
call pfset(ipf)
call brgwscaled(0.,.9)
call setiwarn(0) ! Silence range warnings

! Start user interface loop.
1 continue

call initialize
if(iplottype.eq.1)then
   call plotNs
elseif(iplottype.eq.2)then
   call plotNperps
elseif(iplottype.eq.3)then
   call initomegape
   call plotNomegape
endif

if(ipf.lt.0)then  ! If in no-wait plotting mode always save settings.
   call savesettings
   call exit()
endif

call eye3d(iw)
call uif(iw)
if(iw.ne.0)goto 1
! End of user interface loop

end program coldplas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
