c	relaxation.for
c	CalcRateConst
c	CalcDepopulation
c	CalcGkl
c	ConvertToCoord
c	poffdiag
c	ExpectValue
c	ChargeDistribution
c	StrainEnergy
c	CalcRelaxedStates
c
c-----------------------------------------------------------------------------
c	Subroutine that calculates the rates of the one-photon interstate
c	transition
c-----------------------------------------------------------------------------

	subroutine CalcRateConst

	include 'sizes.fi'
	include 'parameters.fi'
        include 'elstates.fi'
	include 'couplings.fi' 
	include 'constant.fi'
	include 'vibmodes.fi'
	include 'coords.fi'
	include 'rates.fi'

	real(kind(1D0)) gkl(NQMAX2), gab , fe
	real(kind(1D0)) JPlus, JMinus, Okl, Okx, StatWeight, w , a
        real(kind(1D0)) gabelem, DLorentz, DBoseEinstein
        real(kind(1D0)) helem 
	integer lidx,istatemax
	integer ia, ib, kx, k, l, kl, i


	write(6,*)"In CalcRateConst 1"

	open(unit=10,file = "exci.gabxi")
	write(10,*)"## Diabatic coupling constants g(a,b,xi) "
	write(10,*)"#  note that g is dimensionalless"
	write(10,*)"#  gab = hder(kx,i)/sqrt(2*qomega(kx)**3*Hbar)"
        write(10,'(a5,a5,a5,5a12)')"a","b","xi","w_xi","g"
c	Open(unit=13,file ="couplings.mctdh")
c	open(unit=12,file ="parameters.mctdh")
c	open(unit=14,file="parameters.hopping")
c	open(unit=14,file="parameters.tclme")

	write(6,*)" in relaxation nckeep = ",nckeep
	write(6,*)" calculating rate constants for transitions"

	if(iex.eq.0)stop
	

c	basic cycle over normal modes

c	do kx = 1,ncdof
c	   write(14,'(" w[",i3,"] =  ",f12.6,";")') kx,qomega(kx)*HBarMine
c	enddo
c	do kx = 1,nqkeep
c	   write(14,'(" e[",i3,"] = ",f12.6,";")')kx,eigs(kx)
c	enddo


	do kx=1,ncdof
	   Okx=qomega(kx)*HBarMine 
	   wt = sqrt(2.0d0*qomega(kx)**3 *hbarmine)
	   call ZeroVar(gkl, NQMAX2)
	   do i=1,nder
	      ia=ider(i)
	      ib=jder(i)

	      gab=-hder(kx,i)/wt
	      

c	This below was for testing purposes
c       gab=gabelem(mic(1,ia),mic(2,ia),mic(1,ib),mic(2,ib),kx)
c       call ZeroVar(cc,ncdof)
c       fe=helem(mic(1,ia),mic(2,ia),mic(1,ib),mic(2,ib))
c       call Scopy(qmode(1,kx),cc,ncdof)
c       fe=helem(mic(1,ia),mic(2,ia),mic(1,ib),mic(2,ib))-fe
c       call ZeroVar(cc,ncdof)
c       fe=fe/sqrt(2.0d0*qomega(kx)**3*HbarMine) 
c       if (fe.ne.0.0d0) write(6,101) ia,ib,gab,fe
c       gab=fe
c       101			format(2i4,2f8.4)

	      w=2.0d0
c	for diagonal element

	      if (ia.eq.ib) w=1.0d0
	      
c	off-diagonal only
	      
	      if (gab.ne.0.0d0) then
		 do k=1,nqdof
		    do l=1,k
		       kl=lidx(k,l)
		       gkl(kl)=gkl(kl)+gab*h(ia,k)*h(ib,l)*w
		    end do
		 end do
	      endif 
	   end do
c       pause



c       G(k,l) matrix for each mode
c  for the mctdh calculation, we'll only consider dynamics within the lowest 9 electronic states
c  add more here if nesc.

	   

c	   do k = 1,nqkeep


c       do l = 1,nqkeep
c		 kl = lidx(k,l)
		 
c		 write(10,'(3i5,5f12.6)')k,l,kx,qomega(kx),gkl(kl)
c		 if(abs(gkl(kl)).gt.1e-10) then
c		 write(14,'("G[",i3,",",i3,",",i3,"] =  ",e15.7,"*",e15.7,";")')k,l,kx,gkl(kl),qomega(kx)*hbarmine
c		    if(kx.lt.10) then
c		       write(12,'(3(a1,i1)," = ",f12.6,", ev")')"g",k,"_",l,"_",kx,gkl(kl)
c
c		       write(13,'(3(a1,i1),5x,"|1  S",i1,"&",i1,3x,"|",i1,"  q")')"g",k,"_",l,"_",kx,k,l,kx
c		    endif
		    
c		    if((kx.ge.10).and.(kx.lt.100)) then
c		       write(12,'(2(a1,i1),(a1,i2)," = ",f12.6,", ev")')"g",k,"_",l,"_",kx,gkl(kl)
c		       write(13,'(2(a1,i1),(a1,i2),5x,"|1  S",i1,"&",i1,3x,"|",i2,"  q")')"g",k,"_",l,"_",kx,k,l,kx
c		    endif
		    
c		    if((kx.gt.100).and.(kx.le.999)) then
c		       
c		       write(12,'(2(a1,i1),(a1,i3)," = ",f12.6,", ev")')"g",k,"_",l,"_",kx,gkl(kl)
c		       write(13,'(2(a1,i1),(a1,i3),5x,"|1  S",i1,"&",i1,3x,"|",i3,"  q")')"g",k,"_",l,"_",kx,k,l,kx
c		       
c		    endif
c		 endif
		 
		 
c	      enddo
c	   enddo


c       call WritePackedMatrix(6,gkl,nqdof,1)
c       pause
c	now gkl is calculated
c	cycle over the transitions now
	   do k=1,nqdof
	      do l=1,nqdof
		 if (k.ne.l) then 	 
		    kl=lidx(k,l)
		    Okl=(eigs(k)-eigs(l))/HBarMine
		    JPlus=DLorentz(Okl-Okx,PhGama/HBarMine)
		    JMinus=DLorentz(-Okl-Okx,PhGama/HBarMine)  
		    JMinus=DLorentz(Okl+Okx,PhGama/HBarMine)  
		    StatWeight=1.0d0+DBoseEinstein(Okl*HBarMine,EKbT)

		    rk(k,l)=rk(k,l)+gkl(kl)**2*Okl**2*StatWeight*(JPlus-JMinus)*2.0*Pi

c		    rk(k,l)=rk(k,l)+(gkl(kl)*Okl)**2*StatWeight*(JPlus-JMinus)*2.0d0*Pi

		 endif
	      end do
	   end do   
	end do

	open(12,file="rateconstants.out")
	do i = 1,nqdof
	   do j = i+1,nqdof
	      write(12,'(2i6,3e15.6)')i,j,rk(i,j),rk(j,i),rk(i,j)/rk(j,i)
	   enddo
	enddo

c       104	format(a,100e12.4)
	close(12)
	close(13)
	close(14)

	return
	end
c       
	subroutine CalcDepopulation(rk,n,ndim)
	integer n,ndim,k,l
	real(kind(1D0)) rk(ndim,*)
	real(kind(1D0)) a
	do k=1,n
	   a=0.0d0
	   do l=1,n
	      a=a+rk(k,l)
	   end do  
	   rk(k,k)=-a
	end do
	return
	end
c       
c       
	subroutine CalcGkl(ck,cl,q)
c	Calculate the normal displacement q(kx) = b(kx) = g(kx,k,k) of the
c	electronic state k, given by vec 
c       ck,cl = electronic states
c       q = normal displacement vector
c       
	include 'sizes.fi'
	include 'elstates.fi'
	include 'couplings.fi' 
	include 'constant.fi'
	include 'vibmodes.fi'
c       
	real(kind(1D0)) q(*), ck(*), cl(*)
	real(kind(1D0)) gab, w,sum,ww,w3,hbar2
	integer kx, ia, ib, i
c       
        hbar2=sqrt(hbarmine)
	
	call ZeroVar(q, ncdof)
	do i=1,nder
	   ia=ider(i)
	   ib=jder(i)
	   if (ia.eq.ib) then
	      w=1.0d0
	   else
	      w=2.0d0
	   endif
c	   write(12,'(3a6,10a12)')"ia","ib","kx","gab","hder","qomega","ck","cl","w"

	   do kx=1,ncdof
c       gab=gabelem(mic(1,ia),mic(2,ia),mic(1,ib),mic(2,ib),kx)
	      ww = qomega(kx)
	      w3 = ww*ww*ww
	      
	      gab=hder(kx,i)/sqrt(2.0d0*w3) 
c              write(12,'(3i6,10f12.6)')ia, ib, kx,gab,hder(kx,i),qomega(kx),ck(ia),cl(ib),w

	      q(kx)=q(kx)+gab*ck(ia)*cl(ib)*w/hbar2
	   end do


c	   write(6,'("q = ",1000f12.6)')(q(kx),kx=1,ncdof)
c	   pause
	end do
	return
	end
c
        subroutine CalcGkl_old(ck,cl,q)
c       Calculate the normal displacement q(kx) = b(kx) = g(kx,k,k) of the
c       electronic state k, given by vec
c
        include 'sizes.fi'
        include 'elstates.fi'
        include 'couplings.fi'
        include 'constant.fi'
        include 'vibmodes.fi'
c
        real(kind(1D0)) q(*), ck(*), cl(*)
        real(kind(1D0)) gab, w
        integer kx, ia, ib, i
c

        call ZeroVar(q, ncdof)
        do i=1,nder
           ia=ider(i)
           ib=jder(i)
           if (ia.eq.ib) then
              w=1.0d0
           else
              w=2.0d0
           endif
c           if(debug)write(11,*)"--kx    gab     hder    qomega     ck     cl    w"
           do kx=1,ncdof



              gab=hder(kx,i)/sqrt(2.0d0*qomega(kx)**3*HbarMine)

c              if(debug)write(11,'(i6,10f12.6)')kx,gab,hder(kx,i),qomega(kx),ck(ia),cl(ib),w
              q(kx)=q(kx)+gab*ck(ia)*cl(ib)*w
           end do
        end do
        return
        end

c       
	subroutine ConvertToCoord(q, coord, fac)
c	Transform a vector in normal coordinates (q) 
c	into site coordinates (coord)
	include 'sizes.fi'
	include 'constant.fi'
	include 'vibmodes.fi'
	
	real(kind(1D0)) q(*), coord(*)
	real(kind(1D0)) w, fac
	integer kx, i
	call ZeroVar(coord, ncdof)
	do kx =1, ncdof
	   w=sqrt(2.0d0*HBarMine/qomega(kx))
	   do i= 1,ncdof
	      coord(i)=coord(i)+w*q(kx)*qmode(i,kx)*fac
	   end do
	end do
c	write(6,101) 'coord', (coord(i),i=1,ncdof)
c 101	format(a8,1000e10.3)
c	pause
	return
	end      
c       
c       
	real(kind(1D0)) function poffdiag(k,l,time)
	include 'sizes.fi'
	include 'elstates.fi' 
	include 'coords.fi'
	include 'rates.fi'
	include 'constant.fi'
	integer k,l
	real(kind(1D0)) time, omega, gama, p0
	omega=(eigs(k)-eigs(l))/HBarMine
	gama=0.5d0*(rk(k,k)+rk(l,l))
	p0=qcaux(k)*qcaux(l)
	poffdiag=p0*cos(omega*time)*exp(gama*time)
	return
	end
c       
c       
	real(kind(1D0)) function ExpectValue(pdiag, time, oper)
c	Calculate the expectation value = tr[Ro*Oper]
c       where	Ro is the density matrix
c       oper is an operator
c	The function oper(k,l) returns the matrix element of the operator
c	in eigenstate representation.
c	pdiag(1:nqdof) - state populations (the diagonal of the density 
c	matrix)
	include 'sizes.fi'
	external oper
	real(kind(1D0)) poffdiag
	real(kind(1D0)) time, EVal
	real(kind(1D0)) pdiag(*)
	real(kind(1D0)) oper
	integer i,j
	EVal=0.0d0
	do i=1, nqdof
	   do j=1, i-1
	      EVal=Eval+poffdiag(i,j,time)*oper(i,j)  
	   end do
	end do
	Eval=Eval*2.0d0
	do i=1, nqdof
	   Eval=Eval+pdiag(i)*oper(i,i)
	end do
	ExpectValue=Eval
	return
	end
c       
	subroutine ChargeDistribution(pdiag, time, hdist, edist)
	include 'sizes.fi'
	include 'elstates.fi'
	real(kind(1D0)) hdist(*), edist(*)
	real(kind(1D0)) pdiag(*)
	real(kind(1D0)) time, poffdiag, weight
	integer k, l, i, hi, ei
	call ZeroVar(hdist, nsites)
	call ZeroVar(edist, nsites)
	do k=1,nqdof
	   do l=1,k-1   
	      weight=poffdiag(k,l,time)*2.0d0
	      if (dabs(weight).gt.1.0d-4) then
		 do i=1,nqdof 
		    hi=mic(1,i)
		    ei=mic(2,i)
		    hdist(hi)=hdist(hi)+h(i,k)*h(i,l)*weight
		    edist(ei)=edist(ei)-h(i,k)*h(i,l)*weight
		 end do
	      endif
	   end do
	   weight=pdiag(k)
	   do i=1,nqdof 
	      hi=mic(1,i)
	      ei=mic(2,i)
	      hdist(hi)=hdist(hi)+h(i,k)*h(i,k)*weight
	      edist(ei)=edist(ei)-h(i,k)*h(i,k)*weight
	   end do
	end do
	return
	end
c       
	real(kind(1D0)) function StrainEnergy(g)
c	Function that calculates the strain energy
	include 'sizes.fi'
	include 'constant.fi'
	include 'vibmodes.fi'
	real(kind(1D0)) g(*), e
	integer kx
	e=0.0d0
	do kx=1, ncdof
	   e=e+HBarMine*qomega(kx)*g(kx)**2
	end do
	StrainEnergy=e 
	return
	end


	real(kind(1d0)) function force(g)
c	function that calculates the derivative of the strain energy
	include 'sizes.fi'
	include 'constant.fi'
	include 'vibmodes.fi'
        real(kind(1D0)) g(*), f
        integer kx

        f=0.0d0
        do kx=1, ncdof
           f=f+HBarMine*2.0*qomega(kx)*g(kx)
        end do
        force = f
        return
        end
c       
	subroutine CalcRelaxedStates(iw, Gx, Ex, W, EDiab, Vdiab, spin)
c	Subroutine that calculates the relaxed (adiabatic) eigenstates
c       
	include 'sizes.fi'
	include 'constant.fi'
	include 'couplings.fi'
	include 'coords.fi'
	include 'elstates.fi'
	include 'vibmodes.fi'
	include 'parameters.fi'
c       
	integer spin
	real(kind(1D0)) Gx(NCMAX,*), Ex(NQMAX)
	real(kind(1D0)) W(NQMAX,NQMAX)
	real(kind(1D0)) E(NQMAX)
	real(kind(1D0)) EDiab(*)
	real(kind(1D0)) Vdiab(NQMAX,*)
	real(kind(1D0)) wk
	real(kind(1D0)) StrainEnergy, ScalProd
	real(kind(1D0)) eold
	real(kind(1d0)) sum,detol
	real(kind(1d0)) wtol
	
	integer it
	integer nc
	integer nk
	integer kxi
	integer iw
	integer max
	integer maxsteps
	logical haltcode	

	
	maxsteps = max(maxits,100)


	haltcode = .false.

	detol = 1.d-4

	wtol = 0.4d0
	
	write(6,*)" iw = ",iw

	open(unit=iw,err=990,file="relaxed_density.out")

	write(6,*)"In Relaxation.f---- 1"

	write(6,'("  Number of adiabatic states requested: ",i5)') min(nqdof,maxrelax)
	write(6,'("  Max number of states: nqdof = ",i5)')nqdof
	write(6,'("Max number of SCF iterations requested: ",i5)')maxsteps
	if(min(nqdof,maxrelax).eq.0) then
	   write(6,*)"---- ERROR ------"
	   call stoppg(1)
	endif

c       
	do nc=1,min(nqdof,maxrelax)

	   call CalcGkl(Vdiab(1,nc),Vdiab(1,nc),Gx(1,nc))
	   write(6,'(i3,a10,1000f12.6)')nc,"gx: ",(gx(k,nc),k=1,ncdof)
	   call ConvertToCoord(Gx(1,nc),cc,-1.0d0)
	   write(6,'(i3,a10,1000f12.6)')nc,"cc: ",(cc(k),k=1,ncdof)



           Ex(nc)=EDiab(nc)
	   write(6,*)"Ex(nc) = ",ex(nc)
           eold=Ex(nc)
           do kxi=1, ncdof
              Ex(nc)=Ex(nc)-HbarMine*qomega(kxi)*Gx(kxi,nc)**2
	      write(6,*)"kxi = ",kxi
	      write(6,*)"Ex(nc) = ",ex(nc)
	      write(6,*)"Gx(kxi,nc) = ",Gx(kxi,nc)
           end do


	   call ZeroVar(cc,ncdof)
	   it=0 
	   wk=1.0d0
	   tol = -1.0d-4
	   do while ((Ex(nc)-eold.lt.tol).and.(it.le.maxsteps))

	      eold=Ex(nc)
	      it=it+1
	      call ConvertToCoord(Gx(1,nc),cc,-1.0d0)
	      call HSparseCalc(hsp,spin)
	      call Hbisect(hsp,nqdof,min(nc+1,nqdof),E,W,NQMAX)
	      nk=nc
	      wk=dabs(ScalProd(Vdiab(1,nc),W(1,nk),nqdof))
	      if (wk.lt.0.5d0) then
		 do j= min(nc+1,nqdof),1,-1
		    if (dabs(ScalProd(Vdiab(1,nc),W(1,j),nqdof)).gt.wk) then
		       nk=j
		       wk=dabs(ScalProd(Vdiab(1,nc),W(1,nk),nqdof))
		    endif
		 end do
	      endif
	      if (wk.lt.0.4d0) go to 200
	      Ex(nc)=E(nk)+StrainEnergy(Gx(1,nc))


	      call CalcGkl(W(1,nk),W(1,nk),Gx(1,nc))

	   end do 
c       if (nc.eq.1) call SCopy(W(1,1),xtrelax,nqdof)
c       if (nc.eq.nct) call SCopy(W(1,nk),ctrelax,nqdof)

	   write(6,103) nc, nk, Ex(nc)-EDiab(nc), it, wk, Ex(nc)
	   call ConvertToCoord(Gx(1,nc),cc,-1.0d0)
	   write(6,*)"strain energy = ",strainenergy(gx(1,nc))
	   write(6,*)"strain force = ",Force(Gx(1,nc))
	   write(6,'(a12,100f12.6)')" normal:",(gx(k,nc),k=1,ncdof)
	   write(6,'(a12,100f12.6)')"lattice:",(cc(k),k=1,ncdof)


	   go to 500
 200	   write(6,104) nc, wk
 500	   continue 
	   write(iw,*)"Relaxed Eigenstate ",nc
	   write(iw,*)"Adiabtic Energy  = ",ex(nc)
	   
	   call EHDensity(iw,W(1,nk))
	   
	   

	end do   

 110	format(a10,i6)

c$$$c  	nc = diabatic
c$$$c       nk = relaxed state index.

c$$$
c$$$c add statement inquiring of "Stop command" has been issued
c$$$	   inquire(file="EXCITON_STOP",exist=haltcode)
c$$$	   if(haltcode) call stoppg(1)




 103	format(' State:',i4,i4,' Shift:',f10.4,' Iter:',i4,' at:',3f10.4)
 104	format(' State:',i4,' missed at',f10.4)
	call ZeroVar(cc,ncdof)
	close(iw)
	close(iw+1)
	return

 990	Continue
	write(6,*)"problem opening unit = ",iw
	
	end
	




	subroutine CalcRelaxedStates2(iw, Gx, Ex, W, EDiab, Vdiab, spin)
c	Subroutine that calculates the relaxed (adiabatic) eigenstates
c       
	include 'sizes.fi'
	include 'constant.fi'
	include 'couplings.fi'
	include 'coords.fi'
	include 'elstates.fi'
	include 'vibmodes.fi'
	include 'parameters.fi'
c       
c       
	integer spin
	real(kind(1D0)) Gx(NCMAX,*), Ex(NQMAX)
	real(kind(1D0)) W(NQMAX,NQMAX)
	real(kind(1D0)) E(NQMAX)
	real(kind(1D0)) EDiab(*)
	real(kind(1D0)) Vdiab(NQMAX,*)
	real(kind(1D0)) wk
	real(kind(1D0)) StrainEnergy, ScalProd
	real(kind(1D0)) eold
	integer it, nc, nk, j, kxi, iw
c       
	do nc=1,nqdof
	   call CalcGkl(Vdiab(1,nc),Vdiab(1,nc),Gx(1,nc))
	   Ex(nc)=EDiab(nc)
	   eold=Ex(nc)
	   do kxi=1, ncdof
	      Ex(nc)=Ex(nc)-HbarMine*qomega(kxi)*Gx(kxi,nc)**2
	   end do
	   call ZeroVar(cc,ncdof)
	   it=0 
	   wk=1.0d0
	   do while (Ex(nc)-eold.lt.-1.0d-4)
	      eold=Ex(nc)
	      it=it+1
	      call ConvertToCoord(Gx(1,nc),cc,-1.0d0)
	      call HSparseCalc(hsp,spin)
	      call Hbisect(hsp,nqdof,min(nc+1,nqdof),E,W,NQMAX)
	      nk=nc
	      wk=dabs(ScalProd(Vdiab(1,nc),W(1,nk),nqdof))
	      if (wk.lt.0.5d0) then
		 do j= min(nc+1,nqdof),1,-1
		    if (dabs(ScalProd(Vdiab(1,nc),W(1,j),nqdof)).gt.wk) then
		       nk=j
		       wk=dabs(ScalProd(Vdiab(1,nc),W(1,nk),nqdof))
		    endif
		 end do
	      endif
	      if (wk.lt.0.4d0) go to 200
	      Ex(nc)=E(nk)+StrainEnergy(Gx(1,nc))
	      call CalcGkl(W(1,nk),W(1,nk),Gx(1,nc))
	   end do 
c       if (nc.eq.1) call SCopy(W(1,1),xtrelax,nqdof)
c       if (nc.eq.nct) call SCopy(W(1,nk),ctrelax,nqdof)
 110	   format(a10,i6)
	   write(iw,103) nc, nk, Ex(nc)-EDiab(nc), it, wk, Ex(nc)
	   go to 500
 200	   write(iw,104) nc, wk
 500	   continue 
	end do   
 103	format(' State:',i4,i4,' Shift:',f10.4,' Iter:',i4,' at:',3f10.4)
 104	format(' State:',i4,' missed at',f10.4)
	call ZeroVar(cc,ncdof)
	return
	end
	
