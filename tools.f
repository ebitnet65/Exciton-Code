
ctools.f,v 1.2 2004/09/29 00:48:13 bittner Exp
c2004/09/29 00:48:13
ctools.f,v
cRevision 1.2  2004/09/29 00:48:13  bittner
cMajor changes from original code written by Stoyan.
c
c	tools.for
c	Collection of subroutines and functions
c	Zerovar
c	ZeroSite
c	SCopy
c	SAdd
c	SSubt
c	SMult
c	NormVector
c	NormPop
c	SMin
c	SMax
c	IMax
c	SSum
c	lidx
c	ScalProd
c	DScalProd
c	TraceOfMatProd
c	Swapvec
c	SwapReal
c	CalcPMat
c	swapint
c	Sort
c	SpinState
c	EnergyToOmega
c	DLorentz
c	DBoseEinstein
c	IndexofE
c	EofIndex
c
      subroutine Zerovar(a,n)
c     Subroutine to initialize variables to zero.
c     a - array of variable
c     n - max number of a variables
      integer i, n
      real(kind(1D0)) a(n)
      do i = 1,n
         a(i)=0.0d0
      end do
      return
      end  
c     
c     
      subroutine ZeroSite(n)
c     Subroutine to initialize site parameters
c     param_site - specific site to initialize
      include 'sizes.fi'
      include 'structure.fi'
      integer n


      x(1,n) = 0.0d0
      x(2,n) = 0.0d0
      x(3,n) = 0.0d0
      tdipole(1,n)= 0.d0
      tdipole(2,n)= 0.d0
      tdipole(3,n)= 0.d0
      eelec(n)=0.0d0
      ehole(n)=0.0d0
      chain = 1
      return
      end  
c     
c     
      subroutine SCopy(source,dest,n)
c     Subroutine that copies values from source into dest
c     source - source
c     dest - destination
      integer i, n 
      real(kind(1D0)) source(n), dest(n)
c     101	format(' n=',i5)      
      do i=1,n
         dest(i)=source(i)
      end do
      return
      end
c     
c     
      subroutine SAdd(a,b,c,n,fact)
c     Subroutine that calculates the scalar sum of a and b and assign to c
      integer i, n
      real(kind(1D0)) a(n),b(n),c(n),fact
      do i=1,n
         c(i)=a(i)+b(i)*fact
      end do
      return
      end       
c     
c     
      subroutine SSubt(a,b,c,n)
      integer i, n
      real(kind(1D0)) a(n),b(n),c(n)
      do i=1,n
         c(i)=a(i)-b(i)
      end do
      return
      end       
c     
c     
      subroutine SMult(a,b,n,fact)
      integer i, n
      real(kind(1D0)) a(n),b(n),fact
      do i=1,n
         b(i)=a(i)*fact
      end do
      return
      end   
c     
c     
      subroutine NormVector(a,n)
      integer n
      real(kind(1D0)) a(n), s
      real(kind(1D0)) ScalProd 
      s=ScalProd(a,a,n)
      if (s.gt.1.0d-20) then
         s=1.0d0/dble(sqrt(s))
         call SMult(a,a,n,s)
      end if
      return
      end         
c     
c     
      subroutine NormPop(a,n,anorm)
      integer n, i
      real(kind(1D0)) a(n), s, anorm
      s=0.0d0
      do i=1,n
         s=s+dabs(a(i))
      end do
      if (s.gt.1.0d-5) then
         s=anorm/s
         do i=1,n
            a(i)=dabs(a(i))*s
         end do
      end if
      return
      end
c     
c     
      function SMin(a,n)
      integer i, n
      real*8 smin
      real(kind(1D0)) a(n), s   
      s=a(1)
      do i=1,n
         if (a(i).lt.s) s=a(i)
      end do
      Smin=s
      return
      end      
c     
c     
      real function SMax(a,n)
      integer n, i
      real(kind(1D0)) s, a(n)   
      s=a(1)
      do i=1,n
         if (a(i).gt.s) s=a(i)
      end do
      SMax=s
      return
      end
c     
c     
      integer function IMax(a,n)
      integer n, i, im
      real(kind(1D0)) a(n)   
      im=1
      do i=2,n
         if (a(i).gt.a(im)) im=i
      end do
      IMax=im
      return
      end
c     
c
      real(kind(1D0)) function SSum(a,n)
c     Returns the sum of all elements of array a
      integer i, n
      real(kind(1D0)) s, a(n)   
      s=0.0d0
      do i=1,n
         s=s+a(i)
      end do
      SSum=s
      return
      end      
c     
c     
      integer function lidx(i,j) 
      integer i,j
      if (i.gt.j) then 
         lidx=((i-1)*i)/2+j
      else
         lidx=((j-1)*j)/2+i
      endif    
      return
      end
c     
c     
      real(kind(1D0)) function ScalProd(a,b,n)
      
c     Function that calculates the scalar product of a and b
      integer i, n
      real(kind(1D0)) s, a(n), b(n)
      s=0.0d0
      do i=1, n
         s=s+a(i)*b(i)
      end do
      ScalProd=s
      return
      end
c     
c     
      real function DScalProd(a,b,c,n)
      integer i, n
      real(kind(1D0)) s, a(n), b(n), c(n)
      S=0.0d0
      do i=1, n
         S=S+a(i)*b(i)*c(i)
      end do         
      DScalProd=s
      return
      end   
c     
c     
      real function TraceOfMatProd(a,b,n, nevec)
      integer n, nevec, i
      real(kind(1D0)) ScalProd
      real(kind(1D0)) s, a(nevec,n), b(nevec,n)
      S=0.0d0
      do i=1, n
         S=S+ScalProd(a(1,i),b(1,i),n)
      end do         
      TraceOfMatProd=s
      return
      end
c     
c     
      subroutine Swapvec(a,b,n)
      integer i, n
      real(kind(1D0)) t, a(n), b(n)
      do i=1,n
         t=a(i)
         a(i)=b(i)
         b(i)=t
      end do
      return
      end
c     
c     
      subroutine SwapReal(a,b)
      real(kind(1D0)) a, b, t
      t=a
      a=b
      b=t
      return
      end
c     
c     
      subroutine CalcPMat(c,p,oc,n,nevec)
      integer n, nevec, i, j, k
      real(kind(1D0)) p(nevec,n),c(nevec,n),oc(n)
      call zerovar(p,nevec*n)
c     
      do k=1,n
         if (oc(k).gt.0.0d0) then
            do i=1,n
               do j=1,i
                  p(i,j)=p(i,j)+c(i,k)*c(j,k)*oc(k)
               end do
            end do
         endif 
      end do
c     
      do i=1,n
         do j=1,i
            p(j,i)=p(i,j)
         end do
      end do
c     
      return
      end
c     
c     
      subroutine swapint(i1,i2)
      integer i, i1, i2
      i=i1
      i1=i2
      i2=i
      return
      end
c     
c     
      subroutine Sort(a,num,n,iord)
c     Sorting subroutine
c     on input:	a(1:n) real array                              
c     iord > 0 increasing; iord < 0 decreasing 
c     on output:	num(1:n) such that a(num(i)) is a sorted array      
      integer n, iord, i, j, k, numx
      real(kind(1D0)) a(n)
      integer num(n)
      do i=1,n
         num(i)=i
      end do
      do i=1,n
         k=i
         do j=i+1,n
            if ((a(num(j))-a(num(k)))*iord.lt.0.0d0) k=j	
         end do
         numx=num(k)
         num(k)=num(i)
         num(i)=numx
      end do
c     write(6,101) (num(i),i=1,n)
c     101	format(' Sort',/,1000i3)       
      return
      end
c     
c     
      character function SpinState(i)
      include 'sizes.fi'
      integer i
      if (i.ge.nqvars) then
         SpinState='S'
      else
         if (i.le.nqdof.and.ispinstate.ne.1) then
            SpinState='S'
         else
            SpinState='T'
         endif 
      endif
      return
      end
c     
c     
      real function EnergyToOmega(e)
      real(kind(1D0)) e
      include 'constant.fi'
      EnergyToOmega= e/HBarMine
      return
      end
c     
c     
      real(kind(1D0)) function DLorentz(En, Gama)
      real(kind(1D0)) En, Gama
c     write(6,101) En, Gama
c     101	format(2x,2e10.3)
c     pause 
      DLorentz=Gama/(Gama**2+En**2)
      return
      end
c     
c     
      real(kind(1D0)) function DBoseEinstein(En,KbT)
      include 'sizes.fi'
      include 'constant.fi'
      real(kind(1D0)) En, KbT, E
      if (abs(En).gt.1.0d-4*abs(KbT)) then
         E=exp(En/kbT)
         DBoseEinstein=1.0d0/(E-1.0d0)
      else
         if (En*KbT.gt.0.0d0) then
            DBoseEinstein=1.0d4
         else
            DBoseEinstein=-1.0d4
         endif
      endif
      return
      end
c     
c     
      integer function IndexofE(e,e1,e2,numpoints)
      real(kind(1D0)) e,e1,e2
      integer numpoints
      IndexofE= nint(1.0d0+(e-e1)*dble(numpoints-1)/(e2-e1))
      return
      end
c     
c     
      real function EofIndex(i,e1,e2,numpoints)
      real(kind(1D0)) e1,e2
      integer i, numpoints
      EofIndex= e1 + (e2-e1)*dble(i-1)/dble(numpoints-1)
      return
      end
c     
