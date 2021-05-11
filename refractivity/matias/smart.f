c<---><---><---><---><---><---><--->*<---><---><---><---><---><---><--->
c
      subroutine smartexponentialneff(r1, r2, n1, dr1r2, neff)
c
c     Evaluating average exponential index of refraction along a
c     straight line that goes from point r1 to point r2.
c
c     The matter path is defined via
c
c                           point r2
c          neff = 1+ 1E-6 INTEGRAL'  ns exp(kr.zv) dl / dr1r2
c                           point r1
c
c     where ns exp(kr.z) is the exponential refractivity model, 
c     which depends only on zv, the vertical altitude. The prime in the
c     integral indicates that it is done along the specified axis.
c
c     The index of refraction model must be properly initialized before
c     invoking this routine.
c
c     Written by: S. J. Sciutto, La Plata 2020., M.J. Tueros, Paris 2020
c
c
c     Arguments:
c     =========
c
c     r1, r2.......... (input, double precision, array(3)) Coordinates
c                      of the initial point.
c     n1   ........... (output, double precision) n at the starting point (in g/ m.cm2).
c     dr1r2........... (output, double precision) Metric distance.  The returned value is -1 if the ground is touched.                      
c                      (in m) between the two points.
c     neff........... (output, double precision)  The corresponding average n (in g/ m.cm2).
c
c
      implicit none
c
c     Compilation parameters.
c
c     include 'constants.f'
c     include 'initpar.f'
c     include 'initcomm.f'      
c     include 'fieldcomm.f'      
c
c     Declaration of arguments.
c
      double precision  r1(3)
      double precision  r2(3)
      double precision  dr1r2, neff,n1
c
c     Declaration of parameters
c
      double precision ns, kr_m, rearth, groundz
      parameter (ns=315.d0, kr_m=-1.218d-4, rearth=6370949.d0)
      parameter (groundz=-1000.d0)

c     Declaration of internal variables and arrays.
c
      double precision  sum,ds,zv1,zv2,Ref1,Ref2,s
      double precision  uxp(3),r1e(3)
      integer           l,nint,i

c     double precision  adstymfromz,depthfromz
c     external          adstymfromz
c     external          depthfromz
c
c     FIRST EXECUTABLE STATEMENT
c
      neff = 0
      dr1r2 = 0
c
c     Straight line parameterization.
c
      dr1r2 = 0
      do l = 1, 3
        r1e(l) = r1(l)
        uxp(l) = r2(l) - r1(l)
        dr1r2  = dr1r2 + uxp(l) ** 2
      enddo
c
      if (dr1r2 .le. 0) then
       write(*,*) "starting distance is 0,this is worrisome"
       dr1r2=0d0
       return
      endif
c
      dr1r2   = sqrt(dr1r2)
      r1e(3)  = r1e(3) + rearth
c
c     Determine number of steps
c      
      nint=(dr1r2/20000)+1      
c      write(*,*) "nint = ",nint
      ds=1.0/nint
c
c     Get altitude of starting point
c      
      zv1=0.d0
      do l = 1, 3
        zv1 = zv1 + (r1e(l)) ** 2
      enddo
      zv1 = sqrt(zv1) - rearth
c
      if(zv1.lt.groundz) then
        write(*,*) "starting point is underground, this is wrong"
        dr1r2=-1
        return        
      endif
c
c     get the Refractivity at the starting point        
c
      Ref1=exp(kr_m*zv1)      
c
c     get the index of refraction at the start, for the output
c
      n1=1.d0+1.d-6*ns*Ref1
c           
c     now the integral
c         
      sum=0.d0      
      do i=1,nint      
         s=i*ds
c
c        get altitude of next point
c         
         zv2 = 0
         do l = 1, 3
           zv2 = zv2 + (r1e(l) + uxp(l) * s) ** 2
         enddo
c
         zv2 = sqrt(zv2) - rearth 
c
         if(zv2.lt.groundz) then
           write(*,*) "point is underground, shadow?"
           dr1r2=-1
           return        
         endif
c          
         if(dabs(zv1-zv2) .gt. 0.1) then
           Ref2=exp(kr_m*zv2)
           sum=sum+(Ref2-Ref1)/(kr_m*(zv2-zv1))
c           write(*,*) "vertical",sum,i,nint,zv1,zv2,(x1-x2)/(zv2-zv1)
c           write(*,*) "        ",(x1-x2),(zv2-zv1),zv1,zv2,x1,x2            
         else
           Ref2=Ref1
           sum=sum+Ref1
c           write(*,*) "horizont",sum,i,nint,zv1,zv2,adstymfromz(zv2,l)             
         endif          
c
c        now switch x1 with x2
         Ref1=Ref2
         zv1=zv2 
c               
      enddo      
c
c
      neff=1.d0+1.d-6*ns*sum/nint
c      write(*,*) "sum0",sum,nint,dr1r2
c
      return
      end
c     --- End of routine smartexponential
c
c<---><---><---><---><---><---><--->*<---><---><---><---><---><---><--->
