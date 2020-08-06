      subroutine shaper(ef,psi)
c     
      implicit none
c     
      include 'hamilt.inc'
c     
      complex*16 ef(-20000:20000)
      real*8 psi(nparam)
c      integer nparam
c      nn = 500
c     
      real*8 fiomega(-20000:20000),eomega(-20000:20000),
     *     dat(1:4*nfour)
      real*8 taug,fwhmomega,deltat,deltaomega,fimin,fimax,
     *     time
      real*8 omeg(1:2*nparam),y2(1:2*nparam),omega,y,psibar(1:2*nparam)
      integer i,j,k,ii,is,parinter
c     
c      open(unit=27,file='pulse.dat',status='unknown')
c      open(unit=43,file='spectrum.dat',status='unknown')
c     
      taug=dpulse1*1.d15
      fwhmomega=2.355d0/taug
      deltat=100*taug/nfour
      deltaomega=2*pi/dble(nfour)/dble(deltat)
      parinter = nparam/2
c     Initialise fiomega and eomega (FD)
c     
      do 115 i=-20000,20000
         fiomega(i) = 0.0d0
         eomega(i)  = 0.0d0
 115  continue
c     
c      write(6,91)
c 91   format('Starting...')  
      do 110 ii=1,2*nparam
         omeg(ii)=(ii-nparam-0.5)*fwhmomega/24.
         if (ii.gt.parinter.and.ii.lt.3*parinter) then 
            psibar(ii)=psi(ii-parinter)
         else
            if (ii.lt.parinter+1) then 
               psibar(ii)=psi(1)
            else
               psibar(ii)=psi(nparam)
            end if
         end if
c      write(27,403)omeg(ii),psibar(ii)
c 403  format(4e14.5,i8) 
 110  continue
c 111  format(2e14.5)
c      close(27)
c      write(6,92)
c 92   format('Call Spline...')  
      call spline(omeg,psibar,2*nparam,1e32,1e32,y2)
c     
c      write(6,93)
c 93   format('Call Splint...')  
      do 120 i=-4000,4000
         omega=i*fwhmomega/12.d2
         call splint(omeg,psibar,y2,2*nparam,omega,y)
         fiomega(i)=y
         eomega(i)=exp(-((omega*taug)**2)/4.)
c         write(43,98)i,omega,fiomega(i),eomega(i)
c 98      format(i5,3e14.5)
 120  continue
      do 300 i=-nfour+1,nfour
         is=i+nfour
         dat(2*is-1)=eomega(i)*dcos(fiomega(i))
         dat(2*is)=eomega(i)*dsin(fiomega(i))
 300  continue
c      write(6,94)
c 94   format('Call Fourier...')  
      call four1(dat,2*nfour,1)
      do 400 is=0,nfour
         ef(is)=dcmplx(dat(2*is+1),dat(2*is+2))/dsqrt(3.2628d6)
c     normalisation so that ft limited for nfour=16384 is 1)
 400  continue
      do 401 is=1,nfour
         ef(-is)=dcmplx(dat(4*nfour-2*is+1),dat(4*nfour+2-2*is))/
     *        dsqrt(3.2628d6)
 401  continue
      do 405 j=-nfour/2+1,nfour/2
         time=j*deltat
         ef(j)=ef(j)*dcmplx(dcos((nfour/2.d0-1.d0)*deltaomega*time),
     *        dsin((nfour/2.d0-1.d0)*deltaomega*time))
 405  continue
c      write(6,95)
c 95   format('Done')   
c      do 402 i=-nfour+1,nfour
c         write(27,403)pi*2400.d0*i/2.d0/nfour/fwhmomega,abs(ef(i))**2,
c     *        ef(i),i
c 402  continue
c 403  format(4e14.5,i8)
     
c      close(27)
c      close(43)
      return
      end        
