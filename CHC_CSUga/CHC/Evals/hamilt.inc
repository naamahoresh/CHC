      integer nnrot,nfour,nparam
      parameter(nnrot=20,nfour=16384,nparam=64)
c   
      integer jmax
      parameter(jmax=0)
c
      logical lapack
      parameter (lapack=.true.)
c
      real*8 ftol,pi      
      parameter(ftol=1.d-2,pi=3.14159265358979323846264d0)
c
      real*8 timeint,dpulse1
      integer istart
      common / experiment / timeint,dpulse1,istart
c
      complex*16 crham(1:2*nnrot,1:2*nnrot)
      integer ng,ne
      common / hamilt / crham,ng,ne
c
      real*8 fi
      common / phase / fi
c
      real*8 delta1,rxa,bg,dg,be,de
      common / inputs / delta1,rxa,bg,dg,be,de
c
      real*8 factr(0:100)
      common / factor / factr
c
      complex*8 amp(1:2*nnrot)
      common / amplitudes / amp
c
      real*8 rnorm
      common / norm / rnorm
c
      real*8 pcos2(0:100,1:91,1:3500)
      common / align / pcos2
c
      real*8 arpcos2(0:100,1:40,1:40)
      common / arraypcos / arpcos2
c
      real*8 time0
      common /timecom/ time0
