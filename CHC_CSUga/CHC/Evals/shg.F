      subroutine shgval(myproc,icalc,nparams,funcval,psi)
c     *********************************************
c     
c     delta-j=+/- 1 only (no Q-transitions)
c     
c     This program calculates the time evolution of a 5 level system
c     using the dressed state picture
c     
      implicit none
      integer i,igen
      character*4 cpuid    
      character*75 prefix
      character*79 filename
      include 'hamilt.inc'

      integer myproc,icalc,nparams
      real*8  funcval,shg,psi(1:nparam)
c     
      real*8 funcold,time1,cclock
      complex*16 efield(-20000:20000)
c      write(cpuid,'(i4.4)') myproc
c      filename='debug2.dat.'//cpuid
c      open(unit=13,file=filename,status='unknown')
c      filename='debug4.dat.'//cpuid
c      open(unit=21,file=filename,status='unknown')
c      filename='phases.dat.'//cpuid
c      open(unit=31,file=filename,status='unknown')
c      write(31,1347)(psi(i),i=1,80)       
c      open(unit=66,file='amp.dat',status='unknown')
      
       if (nparam.ne.nparams) then
         write(6,'("nparam and nparams differ",I5,I5)'),nparam,nparams
       end if

c 1347 format(1f14.6)
      igen=1
      time0 = cclock()
c     Laser detuning between ground and excited electronic state (s**-1)
      delta1=1.d15
c     Bandwith limited duration of the laser pulse (without pulse shaping)
      dpulse1=10.0d-14
c     Rabi frequency between ground and excited electronic state (s**-1)
      rxa=160.d12
c     Integration timestep
      timeint=5.d-15
      istart=5.d0*dpulse1/timeint
c     Ground state rotational constant
      bg=5.d0
      dg=0.d0
c     Excited state rotational constant
      be=5.d0
      de=0.d0
c     Number of ground state rotational levels included in the calculation
      ng=nnrot
c     Number of excited state rotational levels included in the calculation
      ne=nnrot
c     
      factr(0)=1.d0
      factr(1)=1.d0
      do 299 i=2,100
         factr(i)=factr(i-1)*i
 299  continue
c     
c      write(6,*)(psi(i),i=1,80)
      call shaper(efield,psi)
      shg=0.0d0
      do 49 i=-nfour,nfour
         shg=shg+abs(efield(i))**4
 49   continue
      funcval=shg
      return
      end     
