      FUNCTION zbrent(func,x1,x2,tol)
      IMPLICIT NONE
      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
c     Using Brent's method, find the root of a function func known to lie 
c     between x1 and x2. The root, returned as zbrent, will be refined until 
c     its accuracy is tol.
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) 
     &     write(*,*) 'root must be bracketed for zbrent'
      c=b
      fc=fb
      do iter=1,ITMAX
         if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
            c=a
            fc=fa
            d=b-a
            e=d
         end if
         if(abs(fc).lt.abs(fb))then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         end if
         tol1=2.*EPS*abs(b)+.5*tol
         xm=.5*(c-b)
         if(abs(xm).le.tol1 .or. fb.eq.0.)then
            zbrent=b
            return
         end if
         if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb))then
            s=fb/fa
            if(a.eq.c)then
               p=2.*xm*s
               q=1.-s
            else
               q=fa/fc
               r=fb/fc
               p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
               q=(q-1.)*(r-1.)*(s-1.)
            end if
            if(p.gt.0.) q=-q
            p=abs(p)
            if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q)))then
               e=d
               d=p/q
            else
               d=xm
               e=d
            end if
         else
            d=xm
            e=d
         end if
         a=b
         fa=fb
         if(abs(d).gt.tol1)then
            b=b+d
         else
            b=b+sign(tol1,xm)
         end if
         fb=func(b)
      end do
      write(*,*) 'zbrent exceeding maximum iterations'
      zbrent = b
      end

c$$$      FUNCTION zbrent_test_func(x)
c$$$      zbrent_test_func=x**2-1.
c$$$      end
c$$$      
c$$$      program zbrent_test
c$$$      EXTERNAL zbrent_test_func
c$$$c     FUNCTION zbrent(func,x1,x2,tol)
c$$$      write(*,*) zbrent(zbrent_test_func,0.,1.5,1e-4)
c$$$      end
