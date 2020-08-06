C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.
c
c
c
        subroutine cludcmp(a,n,np,indx,d)
        implicit none
        integer n,np,indx(n),nmax
        real*8 d,TINY
        complex*16 a(np,np),cdum,sum
        parameter (nmax=500,tiny=1.0d-20)
        integer i,imax,j,k
        real*8 aamax,dum,vv(nmax)
        d=1.d0
        do 12 i=1,n
          aamax=0.d0
          do 11 j=1,n
            if (cdabs(a(i,j)).gt.aamax) aamax=cdabs(a(i,j))
11        continue
          if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
          vv(i)=1.d0/aamax
12      continue
        do 19 j=1,n
          do 14 i=1,j-1
            sum=a(i,j)
            do 13 k=1,i-1
              sum=sum-a(i,k)*a(k,j)
13          continue
            a(i,j)=sum
14        continue
          aamax=0.d0
          do 16 i=j,n
            sum=a(i,j)
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
            dum=vv(i)*cdabs(sum)
            if (dum.ge.aamax) then
              imax=i
              aamax=dum
            endif
16        continue
          if (j.ne.imax)then
            do 17 k=1,n
              cdum=a(imax,k)
              a(imax,k)=a(j,k)
              a(j,k)=cdum
17          continue
            d=-d
            vv(imax)=vv(j)
          endif
          indx(j)=imax
          if(a(j,j).eq.0.d0)a(j,j)=TINY
          if(j.ne.n)then
            cdum=1.d0/a(j,j)
            do 18 i=j+1,n
              a(i,j)=a(i,j)*cdum
18          continue
          endif
19      continue
        return
        END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.
C

c
c
c       ************************************************************
c
c
        SUBROUTINE cclubksb(aa,n,np,indx,b)
c
        implicit none
c
        include 'hamilt.inc'
c
        integer n,np,indx(n)
        complex*16 b(n),aa(2*nnrot,2*nnrot)
        integer i,ii,j,ll
        complex*16 sum
        ii=0
        do 12 i=1,n
          ll=indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if (ii.ne.0)then
            do 11 j=ii,i-1
              sum=sum-aa(i,j)*b(j)
11          continue
          else if (sum.ne.0.d0) then
            ii=i
          endif
          b(i)=sum
12      continue
        do 14 i=n,1,-1
          sum=b(i)
          do 13 j=i+1,n
            sum=sum-aa(i,j)*b(j)
13        continue
          b(i)=sum/aa(i,i)
14      continue
        return
        end
c
      subroutine four1(data,nn,isign)
c
      implicit none
c
      integer isign,nn
      real*8 data(2*nn)
      integer i,istep,j,m,mmax,n
      real*8 tempi,tempr
      real*8 theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
c
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.
c
c
	SUBROUTINE spline(x,y,n,yp1,ypn,y2)
        INTEGER n,NMAX
        REAL*8 yp1,ypn,x(n),y(n),y2(n)
        PARAMETER (NMAX=500)
        INTEGER i,k
        REAL*8 p,qn,sig,un,u(NMAX)
        if (yp1.gt..99e30) then
          y2(1)=0.
          u(1)=0.
        else
          y2(1)=-0.5
          u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        do 11 i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.
          y2(i)=(sig-1.)/p
          u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *  1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *  u(i-1))/p
11      continue
        if (ypn.gt..99e30) then
          qn=0.
          un=0.
        else
          qn=0.5
          un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
        do 12 k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
12      continue
        return
        END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.
c
c
        SUBROUTINE splint(xa,ya,y2a,n,x,y)
        INTEGER n
        REAL*8 x,y,xa(n),y2a(n),ya(n)
        INTEGER k,khi,klo
        REAL*8 a,b,h
        klo=1
        khi=n
1       if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
        goto 1
        endif
        h=xa(khi)-xa(klo)
        if (h.eq.0.) pause 'bad xa input in splint'
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *  2)/6.
        return
        END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]3"iA#.

