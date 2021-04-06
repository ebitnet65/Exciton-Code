cintegrate.f,v 1.2 2004/09/29 00:48:13 bittner Exp
c2004/09/29 00:48:13
cintegrate.f,v
cRevision 1.2  2004/09/29 00:48:13  bittner
cMajor changes from original code written by Stoyan.
c
c	integrate.for
c	rk4
c	inittime
c	thetime
c	roundtime
c	reporttime
c	derivs
c	Report
c
	subroutine rk4(y,dydx,n,x,h,yout)
c	Given values for the variables y(1:n) and their derivatives dydx(1:n)
c	known at x, use the fourth-order Runge-Kuta method to advance the 
c	solution over an interval h and return the incremented variables as 
c	yout(1:n), which need not be a distinct array from y. The user 
c	supplies the subroutine derivs(x,y,dydx) which returns the 
c	derivatives dydx at x.
	include 'sizes.fi'
	integer n, i
	real(kind(1D0)) h,x,dydx(n),y(n),yout(n)
	real(kind(1D0)) h6,hh,xh,dym(NQMAX),dyt(NQMAX),yt(NQMAX)
	hh=h*0.5
	h6=h/6.
	xh=x+hh
	do 11 i=1,n
	   yt(i)=y(i)+hh*dydx(i)
 11	continue
	call derivs(xh,yt,dyt,n)
	do 12 i=1,n
	   yt(i)=y(i)+hh*dyt(i)
 12	continue
	call derivs(xh,yt,dym,n)
	do 13 i=1,n
	   yt(i)=y(i)+h*dym(i)
	   dym(i)=dyt(i)+dym(i)
 13	continue
	call derivs(x+h,yt,dyt,n)
	do 14 i=1,n
	   yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
 14	continue
	return
	end
c       
c       
	subroutine inittime(xl,xr,steps,first)
	include 'sizes.fi'
	real(kind(1D0)) xl,xr
	integer steps,first
	real(kind(1D0)) x1,x2,nfac,nmode
	integer nstep,nfirst,kmode
	common /CYCLE/ x1, x2, nfac, nmode, nstep, nfirst, kmode
	x1=xl
	x2=xr
	nstep=steps
	nfirst=first
	nmode=nfirst
	nfac=2.0d0
	if (int(nstep/nfirst).gt.2) then
	   nfac=dble(nstep/nfirst)**(1.0d0/dble(NTIMEMAX-4))
	endif
	kmode=max(1,int(nstep/(x2-x1)))
	do while (nstep/kmode > 1000)
	   kmode=kmode*10
	end do
c       kmode=nfirst
	return
	end
c       
c       
	real(kind(1d0)) function thetime(k)
	real(kind(1D0)) x1,x2,nfac,nmode,thetime
	integer nstep,nfirst,kmode
	common /CYCLE/ x1, x2, nfac, nmode, nstep, nfirst, kmode
	integer k
	thetime=x1+(x2-x1)*dfloat(k)/dfloat(nstep)
	return
	end
c       
c       
	logical function roundtime(k)
	real(kind(1D0)) x1,x2,nfac,nmode
	integer nstep,nfirst,kmode
	common /CYCLE/ x1, x2, nfac, nmode, nstep, nfirst, kmode
	integer k
	roundtime = (mod(k,kmode).eq.0)
	return
	end
c       
c       
	logical function reporttime(k)
	real(kind(1D0)) x1,x2,nfac,nmode
	integer nstep,nfirst,kmode
	common /CYCLE/ x1, x2, nfac, nmode, nstep, nfirst, kmode
	integer k,kk
	nmode=nfirst
	do while (k.gt.int(nmode/kmode)*kmode) 
	   nmode=nmode*nfac
	end do
	reporttime = (k.eq.int(nmode/kmode)*kmode)
	return
	end
c       
c       
	subroutine rkdumb(iw,vstart,nvar,x1,x2,nstep)
	include 'sizes.fi'
	include 'elstates.fi'
	integer nstep,nvar,iw
	real(kind(1D0)) x1,x2,vstart(nvar)
	real(kind(1D0)) ScalProd, ExpectValue, thetime
	logical roundtime, reporttime
c       COMMON /path/ xx,y
c       USES rk4
	integer i,k, count
	real(kind(1D0)) step,x,dv(NQMAX),v(NQMAX)
	external ScalProd
	external DipoleElem, ExpectValue
	call NormPop(vstart,nvar,1.0d0)
	call Scopy(vstart,v,nvar)
	x=x1
	count=1
	call Report(count,x,v,nvar)
	write(iw) x,(v(i),i=1,nvar)

 101	format (1000f10.4)
	step=(x2-x1)/(nstep)
	write(6,102) 0
	do k=1,nstep
	   call derivs(x,v,dv,nvar)
	   call rk4(v,dv,nvar,x,step,v)
	   call NormPop(v,nvar,1.0d0)
	   call Sadd(v,vstart,vstart,nvar,1.0d0)
	   if(x+step.eq.x) pause 'stepsize not significant in rkdumb'
	   x=x1+(x2-x1)*dfloat(k)/dfloat(nstep)
c	   if (roundtime(k)) write(iw) x, (v(i), i= 1,nvar)
	   write(iw) x, (v(i), i= 1,nvar)
	   if (reporttime(k)) then
	      write(6,102) count, int(x), x
	      count=count+1
c       call Report(count,x,v,nvar)
	   endif
	end do
	call NormPop(vstart,nvar,1.0d0)
 102	format('Step ',2i5,f10.4)
	return
	end
c       
c       
	subroutine derivs(x,y,dy,nvar)
	include 'sizes.fi'
	include 'rates.fi'
	integer nvar, k
	real(kind(1D0)) x
	real(kind(1D0)) ScalProd
	real(kind(1D0)) y(NQMAX),dy(NQMAX)
	do k=1, nvar
	   dy(k)=ScalProd(rk(1,k),y,nvar)
	end do
	return
	end
c       
c       
	subroutine Report(count,x,v,nvar)
	include 'sizes.fi'
	include 'parameters.fi'
	include 'elstates.fi'
	include 'spectra.fi'
	integer nvar, count
	real(kind(1D0)) x, v(*)
	real(kind(1D0)) poffdiag
	real(kind(1D0)) DipoleElem
	real(kind(1D0)) ExpectValue
	integer i,j
	if (count.ge.1.and.count.le.NTIMEMAX) then
	   ntime=count
c       call ElVibEmit(v,spectrum(1,ntime),eleft,eright,nspec,ElGama)
c       call SpectralDensity(spectrum(1,ntime), eigs, v
c       +		,eleft, eright, ElGama, nspec, nvar, 0)
c       call DensityofPopStates
c       +		(v,spectrum(1,ntime),eleft,eright,nspec,ElGama)
	   Xtime(ntime)=x
	endif
c       write(6,102) 'time:',x
c       write(6,101) (v(i),i=1,nqdof)
c       do i=2,nqdof
c       write(6,101) (poffdiag(i,j,x),j=1,i-1)
 101	format(2000f10.4)
c       end do
c       write(6,102) 'dipole:',ExpectValue(v,x,DipoleElem)
 102	format(a8,f10.4)
c       pause
	return
	end
c       
	
