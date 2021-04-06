cprepare,v 1.2 2004/09/29 00:48:13 bittner Exp
c2004/09/29 00:48:13
cprepare.f,v
cRevision 1.2  2004/09/29 00:48:13  bittner
cMajor changes from original code written by Stoyan.
c
c	prepare.f
c	prepare
c	lowct
c	PrepareElStates
c	PrepareModes
c	PrepareCouplings
c	mdist
c	CoulombInt
c	ExchangeInt
c	DipoleInt
c	felectron1
c	fhole1
c	edirect
c	exchange
c	eindirect
c	helem
c	HSparseCalc
c	TransitionDipole
c
      subroutine Prepare
c	Subroutine for setting up the energy levels:  calculating the 
c     vibrational modes, diabatic and adiabatic eigenstates
c     
      include 'sizes.fi' 
      include 'parameters.fi'
      include 'elstates.fi'      
      include 'couplings.fi'
      include 'coords.fi'	 
      include 'vibmodes.fi'
      include 'rates.fi'
      include 'constant.fi'
c     
      real(kind(1D0)) ehsym, RadRate
      real(kind(1d0))DBoseEinstein
      integer i, spin
      
      
      
c     Initialize coords
c     reminder:
c     qc and qcaux store the diagonal of the reduced (electronic) density
c     matrix, which have the meaning of eigenstate populations
c     cc and ccaux store the "nuclear" coordinates, i.e. the deviations of
c     the local harmonic oscillators from equilibrium
c     
      
      call ZeroVar(cc,NCMAX)
      call ZeroVar(qc,NQMAX*2)    
      call ZeroVar(ccaux,NCMAX)
      call ZeroVar(qcaux,NQMAX*2)
c     Initialize coupling
      
      
      
      call ZeroVar(Gkk, NCMAX*NQMAX*2)  
      
      
c     Initialize elstates
      call ZeroVar(transdip,NQMAX)
      
      
      call ZeroVar(sradrate,NQMAX)
      
      
      call ZeroVar(eigs,NQMAX*2)
      
      
      call ZeroVar(eadiab, NQMAX*2)
      
      
c     Initialize rates
      call ZeroVar(rk, NQMAX*NQMAX*4)
      
      
c     Generate the monoexcited configurations
      call PrepareElStates(mic)
      
      
c     Calculate the vibrational modes
      call PrepareModes(qomega,qmode)
      
      
      write(6,*)" NORMAL MODE FREQUENCIES (units = eV) "
      write(6,'(3a12)')"mode #"," (rad/ps)","ev"
      do i = 1,ncdof
         write(6,'(i12,2f12.6)')i,qomega(i),qomega(i)*6.58E-4
      enddo
      
c     Calculate el-ph couplings in electronic-site vibration-eigenmode
c     representation
      
      
      call HSparseCalc(hsp,-1)
      call PrepareCouplings
      spin=max(0,ispinstate)
      nqexc=nqdof
c     Generate the electronic hamiltonian and diagonalize it
      write(6,*)" Spin = ",spin
      call HSparseCalc(hsp,spin)
      
c     should be able to use LAPACK routine here. 
      
      call Hbisect(hsp,nqdof,nqdof,eigs,h,NQMAX)		
      call SCopy(eigs, eadiab, nqdof)
      
c     Calculate the squared one-photon transition dipoles and
c     fluorescence rates
      
c     if (spin.eq.0) then


      call TransitionDipole(h,transdip)
      do i=1,nqdof
         
         Okl=eigs(i)/HBarMine
         Coupling=transdip(i)/(6.0*EpsZeroMine*HBarMine**2)
         
c     Photon density of the vacuum
         Density=2.0d0*HBarMine*Okl**3/(Pi*CLightMine**3)
c     Equilibrium population factor
         
         StatWeight=1.0d0+DBoseEinstein(Okl*HBarMine,EKbT)
        

         sradrate(i)=Coupling*StatWeight*Density
      end do
      if(spin.eq.1) then
         write(6,*)"TRIPLET SPIN STATES have been requested"
         write(6,*)"all radiative rates are 0 by spin symmetry"
      endif
      
      write(6,*)"Lowest E(in eV) and RADIATIVE RATES for CI STATES"
      write(6,'(4a12)')"state#","Energy","Rad (1/ns)","<0|x.mu|n>"
      
      DO i = 1,min(nqdof,100)
         write(6,'(i12,10f12.6)')i,eigs(i),sradrate(i),transdip(i)
      enddo
      
      
      
c     Calculated the relaxed (adiabatic) eigenstates
      
      if ( relax) then
         
         write(6,*)" COMPUTING RELAXED STATES"

         call CalcRelaxedStates(60, Gkk, eadiab, work, eigs, h, spin)
         
      endif
      
c     Do the same for triplets, if requested
      
      if (ispinstate.lt.0) then
         spin=1
         nqexc=nqdof*2
         
c     USE LAPACK ROUTINES HERE 
         
         call HSparseCalc(hsp,spin)
         call Hbisect(hsp,nqdof,nqdof,eigs(nqdof+1),h(1,nqdof+1),NQMAX)
         call SCopy(eigs(nqdof+1), eadiab(nqdof+1),nqdof)
C     
         
         
         if ( relax)  call CalcRelaxedStates(6, Gkk(1,nqdof+1), 
     *        eadiab(nqdof+1), work, eigs(nqdof+1), h(1,nqdof+1), spin)
      endif
      
c     Determine e-h symmetry of states
      do i=1,nqexc
         if (ehsym(h(1,i)).ge.0) then
            iehpar(i)=1
         else
            iehpar(i)=-1
         endif
      end do
      return
c     110	format(a10,i6)
c     200	format(a10,3i10)
      end     
c     
c     
      integer function lowct(iehpar,n)
      integer n, iehpar(*)
      integer nct
      nct=1
      do while (nct.le.n.and.iehpar(nct).ge.0)
         nct=nct+1
      end do
      if (nct.gt.n) nct=1
      lowct=nct
      return
      end	      
c     
c     
      subroutine PrepareElStates(p_mic)
c     Subroutine that initialized the monoexcited electronic states
c     nq - i-th monoexcited configuration (counter)
c     idisp - displacement from i-th site
c     p_mic(1,i) - i-th site localized configuration with hole
c     p_mic(2,i) - i-th site localized configuration with electron
      include 'sizes.fi'
      integer p_mic(2,*)
      integer nq, i, idisp
      
      nq=0
      
      
c     monoexcited configurations with both electron and hole on i-th site
      do i=1,nsites
         nq=nq+1
         p_mic(1,nq)=i
         p_mic(2,nq)=i
      end do
      
c     monoexcited configurations with electon and hole idisp sites apart
      
      do idisp=1,nsites-1
         do i=1,nsites-idisp
            nq=nq+1
            p_mic(1,nq)=i
            p_mic(2,nq)=i+idisp
            nq=nq+1
            p_mic(1,nq)=i+idisp
            p_mic(2,nq)=i
         end do
      end do
      
      nqdof=nq
      return
      end
c     
c     
      subroutine PrepareModes(pqomega,pqmode)
c     Prepares the vibrational modes
      include 'sizes.fi'
      include 'constant.fi'
      include 'parameters.fi'
      include 'coords.fi'
c     include 'vibmodes.fi'
      
      real(kind(1D0)) K, G, o0, o1
      real(kind(1d0)) pqmode(NCMAX,NCMAX),pqomega(*)
      integer i,j,info
      integer lwork
      real(kind(1.d0)) scratch(3*NCDOF)
      lwork = 3*NCDOF
      
c     Determine number of degrees of freedom based on presence of a
c     second phonon coupling
      
      
      if (fediag2.ne.0.0d0) then
         ncdof=nsites+nsites
      else
         ncdof=nsites
      endif
      
      call ZeroVar(pqmode,NCMAX*NCMAX)
      
      
      o0=forcek1/HBarMine
      o1=forcek1/HBarMine-forceo1/HBarMine
      K=(o0**2+o1**2)/2.0d0
      G=(o0**2-o1**2)/4.0d0
      
c     this is for a homogeneous chain.  
c     rewrite for a chain with arb. couplings. 
c     the tbond matrix determines whether or not we have 
c     phonon coupling between sites. 
      
      do i=1,nsites
         do j=i,nsites
            pqmode(i,j) = G*tbond1(i,j)
            pqmode(j,i) = G*tbond1(j,i)
         enddo
         pqmode(i,i)=K
      enddo
      
      
      
c     do i=1, nsites-1
c     j=i+1
c     pqmode(i,j)=G
c     pqmode(j,i)=G
c     end do
c     Second phonon branch if present
      if (ncdof.ne.nsites) then
         o0=forcek2/HBarMine
         o1=forcek2/HBarMine-forceo2/HBarMine
         K=(o0**2+o1**2)/2
         G=(o0**2-o1**2)/4
         
c     do i=nsites+1,ncdof
c     pqmode(i,i)=K 
c     end do  
         
c     do i=nsites+1, ncdof-1
c     j=i+1
c     pqmode(i,j)=G
c     pqmode(j,i)=G
c     end do
         
         
         do i=nsites+1,ncdof
            do j=nsites+1,ncdof
               pqmode(i,j) =G*tbond2(i-nsites,j-nsites)
               pqmode(j,i) =G*tbond2(j-nsites,i-nsites)
            enddo
            pqmode(i,i)=K
         enddo
      endif
      
      
      if(debug) then
         write(6,*)"Phonon Hessian Matrix"
         do i = 1,ncdof
            write(6,'(100f12.5)')(pqmode(i,j),j=1,ncdof)
         enddo
      endif
      
      
c     Obtain the eigenvectors and values (normal modes and frequencies)
c     (eigval.for)
      
C     Use Lapack routine here. 
      
C     Default Routine
      call diagw(pqmode,ncdof,pqomega,NCMAX)
C     LAPACK ROUTINE
c     call dsyev('V','U',ncdof,pqmode,NCMAX,pqomega,scratch,3*ncmax,info)
      
      
      do i=1,ncdof
         pqomega(i)=sqrt(dabs(pqomega(i)))
      end do
      
      if(debug) then
         write(6,*)"---first phonon branch ---- "
c     -first phonon branch	
         do i = 1,nsites
            write(6,'(1000f12.6)')(pqmode(i,j+nsites),j=1,nsites) 
         enddo
         
c     -second phonon branch
         write(6,*)"----- second phonon branch -----"
         do i = 1,nsites
            write(6,'(1000f12.6)')(pqmode(i+nsites,j),j=1,nsites)
         enddo
      endif
      
      
      return
      end
c     
c     
      subroutine PrepareCouplings
c     Calculate el-ph couplings in electronic-site vibration-eigenmode
c     representation between the normal mode, hsp(kx), and the element
c     of the electronic Hamiltonian in site representation, H(ia,ib).  The
c     couplings are stored in hder.
      include 'sizes.fi' 
      include 'parameters.fi'
      include  'constant.fi'
      include 'couplings.fi'
      include 'coords.fi'
      include 'vibmodes.fi'
      include 'elstates.fi'
c     real(kind(1D0)) helem
      integer ia, ib, i
      integer kx
      integer h1,h2,e1,e2
      integer lidx
      real(kind(1d0)) ham1,ham2,helem
      call ZeroVar(hder,NQMAX*3)
      do i=1,ncdof
         cc(i)=1.0d0
      end do
      nder=0
      do ia=1,nqdof
         do ib=1,ia
            h1 = mic(1,ia)
            e1 = mic(2,ia)
            h2 = mic(1,ib)
            e2 = mic(2,ib)
            ham1 = helem(h1,e1,h2,e2)
            ham2 = hsp(lidx(ia,ib))
            f1 = felectron1(e1,e2)
            f2 = fhole1(h1,h2)
            
            if (ham1-ham2.ne.0.0d0) then
c     write(6,'(7i6,4f12.6)')h1,e1,h2,e2,ia,ib,lidx(ia,ib),
c     1	      ham1,ham2,f1,f2
               nder=nder+1
               ider(nder)=ia
               jder(nder)=ib
            endif
         end do
      end do
      
      
      if(debug) then
         write(6,*)" generating debug info on prepare.f in fort.15"
         write(15,*)" generated in prepare.f"
         
         write(15,'(7a4,10a12)')"i","ia","ib","h1",
     >        "e1","h2","e2","hder","ho","hsp"
         write(15,*)"-----------------------------------------------"
      endif
      do kx=1,ncdof
         call SCopy(qmode(1,kx),cc,ncdof)
         do i=1,nder
            ia=ider(i)
            ib=jder(i)
            h1 = mic(1,ia)
            e1 = mic(2,ia)
            h2 = mic(1,ib)
            e2 = mic(2,ib)
            ho  = helem(h1,e1,h2,e2)
            hder(kx,i) = ho - hsp(lidx(ia,ib))
            
c     hder(kx,i)=helem(mic(1,ia),mic(2,ia),mic(1,ib),mic(2,ib))
c     1	   -hsp(lidx(ia,ib))
            
            if(debug)write(15,101)i,ia,ib,h1,e1,h2,e2,
     $           hder(kx,i),ho,hsp(lidx(ia,ib))
            
 101        format(7i4,10f12.6)
         end do
      end do
      call ZeroVar(cc,ncdof)
      
      
      return
      end 
c     
c     
      integer function mdist(m1,m2)
      integer m1,m2
      include 'sizes.fi'
      if (mcyclic.ne.0) then
         mdist=iabs(m1-m2)
      else
         mdist=iabs(m1-m2)
      end if
      return
      end
c     
c       
      real(kind(1D0)) function CoulombInt(r)
c     Calculates coulomb integral. See edirect for addn'l explanation.
      include 'parameters.fi'
      real(kind(1D0)) r,denom
      denom = 1.d0 + r/coulombr
      
      if(coulombr.eq.0d0)Write(6,*)"----PROBLEM!!!---"
      CoulombInt=coulombe/(1.0d0+r/coulombr)
      
      return
      end         
c     
c     
      function ExchangeInt(r)
c     Calculates exchange integral. See eindirect for addn'l explanation.
      real*8 exchangeint,r
      include 'parameters.fi'
      
      ExchangeInt=potexe*exp(-r/potexr)
      return
      end
c     
c     
      real(kind(1D0)) function DipoleInt(r)
c     Calculates dipole-dipole coupling. See eindirect for addn'l explantn.
      include 'parameters.fi'
      real(kind(1D0)) r
c     if (r.gt.0.0.and.r.lt.4.0d0) then
c     DipoleInt=dipole/(r**3/dipolrc**3)
c     else
c     DipoleInt=0.0
c     endif
      if(r.gt.0.0) then
         DipoleInt = dipole/((r/dipolr)**3)
      else
         DipoleInt = 0.d0
      endif		
      return
      end
c     
c     
c$$$  real(kind(1D0)) function felectron1(m1,m2)
c$$$  c	Gives the <m1|f|m2> where m1 and m2 are the electon (CB) states
c$$$  c	at sites m1 and m2, and f is the one-electron hamiltonian
c$$$  c       
c$$$  
c$$$  include 'parameters.fi' 
c$$$  include 'sizes.fi'
c$$$  include 'coords.fi'
c$$$  include 'structure.fi'
c$$$  
c$$$  c       
c$$$  integer m1, m2
c$$$  real(kind(1D0)) fval, fcor,HuangRhysS
c$$$  c
c$$$	fval = 0.0d0
c$$$	fcor = 0.0d0
c$$$
c$$$	deltaelec=(forcem-1.0d0)/(forcem+1.0d0)
c$$$
c$$$c	if (thop(m1,m2).ne.0.0d0) then
c$$$	if(iabs(m1-m2).eq.1) then
c$$$	   fval= thop(m1,m2)
c$$$	1	+fediag1*0.40d0*(cc(m1)+cc(m2)) 
c$$$	2	+fediag2*0.40d0*(cc(m1+nsites)+cc(m2+nsites))
c$$$
c$$$	   fcor=fval
c$$$
c$$$	elseif (m1.eq.m2) then
c$$$
c$$$	   fval=sites(m1).eelec-fediag1*cc(m1)-fediag2*cc(m1+nsites)
c$$$
c$$$	   fcor=(sites(m1).eelec-sites(m1).ehole)/2.0
c$$$	1	-fediag1*cc(m1)-fediag2*cc(m1+nsites)
c$$$
c$$$	endif
c$$$
c$$$	felectron1 = fval + fcor*deltaelec
c$$$	return
c$$$	end
c$$$
c$$$
c$$$c
c$$$
c$$$c       
c$$$c       
c$$$	real(kind(1D0)) function fhole1(m1,m2)
c$$$c	Returns <m1|f|m2> where m1 and m2 are the hole (VB) states
c$$$c	at sites m1 and m2, and f is the one-particle hamiltonian
c$$$c       
c$$$	include 'parameters.fi' 
c$$$	include 'sizes.fi'
c$$$	include 'coords.fi'
c$$$	include 'structure.fi'
c$$$c       
c$$$	integer m1, m2
c$$$	real(kind(1D0)) fval, fcor
c$$$c       
c$$$	fcor = 0.0d0
c$$$	fval = 0.0d0
c$$$	deltahole=(1.0d0-forcem)/(forcem+1.0d0)
c$$$
c$$$	if (thop(m1,m2).ne.0.0d0) then
c$$$c	if(iabs(m1-m2).eq.1) then
c$$$	   fval=thop(m1,m2)
c$$$	1	+fediag1*0.4d0*(cc(m1)+cc(m2))
c$$$	2	+fediag2*0.4d0*(cc(m1+nsites)+cc(m2+nsites))
c$$$
c$$$	   fcor=fval
c$$$	   
c$$$	elseif (m1.eq.m2) then
c$$$
c$$$	   fval=-sites(m1).ehole-fediag1*cc(m1)-fediag2*cc(m1+nsites)
c$$$
c$$$	   fcor=(sites(m1).eelec-sites(m1).ehole)/2.0
c$$$	1	-fediag1*cc(m1)-fediag2*cc(m1+nsites)
c$$$
c$$$	endif
c$$$	fhole1=fval+fcor*deltahole
c$$$	return
c$$$	end

c

c
      real(kind(1D0)) function felectron1(m1,m2)
c     Gives the <m1|f|m2> where m1 and m2 are the electon (CB) states
c     at sites m1 and m2, and f is the one-electron hamiltonian
c     
      include 'parameters.fi'
      include 'sizes.fi'
      include 'coords.fi'
      include 'structure.fi'
      
c     
      integer m1, m2
      real(kind(1D0)) fval, fcor
c     
      ewidth = 0.536
      deltaelec=(forcem-1.0d0)/(forcem+1.0d0)
      
c     write(6,*)'deltaelec = ',deltaelec
      
      
c...  orig
c     if (iabs(m1-m2).eq.1) then
c     fval=ewidth+fediag1*0.4d0*(cc(m1)+cc(m2))
c---  mod
c     why the 0.4d0????
c     first term modulates the electronic hopping integral between sites m1 and m2

      if(thop(m1,m2).ne.0.d0) then
         fval = thop(m1,m2)+fediag1*0.4d0*(cc(m1)+cc(m2))
     +        +fediag2*0.4d0*(cc(m1+nsites)+cc(m2+nsites))
         fcor=fval
         felectron1=fval+fcor*deltaelec
         
         
c     this term modulates the site energy by the phonon 
      else
         if (m1.eq.m2) then
            fval=eelec(m1)-fediag1*cc(m1)-fediag2*cc(m1+nsites)
            fcor=(eelec(m1)-ehole(m1))/2.0
     >           -fediag1*cc(m1)-fediag2*cc(m1+nsites)
            felectron1=fval+fcor*deltaelec
         else
            felectron1=0.0d0
         endif
      endif
      return
      end
      

      real(kind(1D0)) function fhole1(m1,m2)
c     Returns <m1|f|m2> where m1 and m2 are the hole (VB) states
c     at sites m1 and m2, and f is the one-particle hamiltonian
c     
      include 'parameters.fi'
      include 'sizes.fi'
      include 'coords.fi'
      include 'structure.fi'
c     
      integer m1, m2
      real(kind(1D0)) fval, fcor
      ewidth = 0.536
      deltahole=(1.0d0-forcem)/(forcem+1.0d0)	
c     write(6,*)'deltahole = ',deltahole
c     pause
c     
c     if (iabs(m1-m2).eq.1) then
c     fval=ewidth+fediag1*0.4d0*(cc(m1)+cc(m2))
      if(thop(m1,m2).ne.0.0d0) then
         fval=thop(m1,m2)+fediag1*0.4d0*(cc(m1)+cc(m2))
     +        +fediag2*0.4d0*(cc(m1+nsites)+cc(m2+nsites))
         fcor=fval
         fhole1=fval+fcor*deltahole
      else
         if (m1.eq.m2) then
            fval=-ehole(m1)-fediag1*cc(m1)-fediag2*cc(m1+nsites)
            fcor=(eelec(m1)-ehole(m1))/2.0
     $           -fediag1*cc(m1)-fediag2*cc(m1+nsites)
            fhole1=fval+fcor*deltahole
         else
            fhole1=0.0d0
         endif
      endif
      return
      end
c
c
c       
c----------------------------------------c
c     
      double precision function exchange(e1,h2,h1,e2)
      integer e1, h1, e2, h2
      include 'sizes.fi'
      include 'coords.fi'
      double precision ExchangeInt, r 
      
      
      r = distmat(e1,e2)
      if (e1.eq.h1.and.e2.eq.h2) then
c     exchange=ExchangeInt(dabs(dfloat(h1-h2)))
         exchange=ExchangeInt(r)
      else
         exchange=0.0d0
      endif
      return
      end
      
      
c----------------------------------------c       
c     real(kind(1D0)) function edirect(m1,n1,m2,n2)
c     Function that calculates direct coulombic interaction
c     J|m-n| = <m1n1|v|m2n2> = Jo (1 + r/ro)^-1
c     where	m1 & m2, and n1 & n2 are orthonormal
c     J|m-n| = CoulombInt
c     |m-n| = r (# of sites separating electron(m) & hole(n))
c     
      double precision function edirect(m1,n1,m2,n2)
      integer m1, n1, m2, n2
      include 'sizes.fi'
      include 'coords.fi'
      real*8  CoulombInt,r
      if (m1.eq.m2.and.n1.eq.n2) then
         r = distmat(m1,n1)
         
c     edirect=CoulombInt(dabs(dfloat(m1-n1)))
         
         edirect=CoulombInt(r)
         
      else
         edirect=0.0d0
      endif
      return
      end
      

       
c----------------------------------------c       
c     eindirect(m1,n1,n2,m2)
c     Calculates indirect coulombic interaction, short-range spin-exchange,
c     K, or same unit transition dipole-dipole coupling, D.
c     K = <m1n1|v|n2m2> = Ko exp(-r/ro)
c     where	m1 & m2, and n1 & n2 are orthonormal
c     r = |m1-n1|
c       Ko = potexe
c     ro = potexr
c     D = <m1n1|v|m2n2> = Do (r/ro)^-3
c     where	m1 & n2, and m2 & n1 are orthonormal and m1 <> m2
c       r = |m1-m2| site separation of electron 1 and 2
c     Do = dipole
c     r0 = dipolr
c     
      double precision function eindirect(m1,n1,n2,m2)
      integer m1, n1, m2, n2
      double precision ExchangeInt, DipoleInt
      double precision ei,r
      include 'sizes.fi'
      include 'coords.fi'
      
c     
      if (m1.eq.m2.and.n1.eq.n2) then
         r = distmat(m1,n1)
c     ei=ExchangeInt(dabs(dfloat(m1-n1)))
         ei=ExchangeInt(r)
      else
         if (m1.eq.n2.and.m2.eq.n1.and.m1.ne.m2) then
            r = distmat(m1,m2)
c     ei=DipoleInt(dabs(dfloat(m1-m2)))
            ei=DipoleInt(r)
         else
            ei=0.0d0
         endif   
      endif 	  
      eindirect=ei
      return
      end     
c--------------------------------------------


c       
c       
      real(kind(1D0)) function helem(h1,e1,h2,e2)
c     Calculate the one-particle CI matrix element of the configurations
c	with a hole at monomer h1(h2) and electron at monomer e1(e2)
c     (Fij) = <h1|h2><e1|f|e2> + <e1|e2><h1|f|h2>
c     where	h1 & h2, and e1 & e2 are orthonormal
c     f is the one-electron hamiltonian
c     include 'sizes.fi' 
      include 'parameters.fi'
      integer h1,e1,h2,e2      
      real(kind(1D0)) f
      real(kind(1D0)) felectron1, fhole1
c     *	, CoulombInt, DipoleInt, edirect, eindirect
      f=0.0d0
      
      if (h1.eq.h2) f=felectron1(e1,e2)
      if (e1.eq.e2) f=f+fhole1(h1,h2)
        
      
      helem=f
c     101	format(f10.4)
      return
      end                                     
c     
c     
      
	
	
      subroutine HSparseCalc(hs,spinstate)
c     Calculate the Hel matrix elements, monoexcited CI in localized basis
c     Hel = sum(Fij + Vij(Singlet/Triplet)) |i><j|
c     where	Fij = CI matrix elements of the one-particle
c     band-structure operator f (calculated by helem)
c     Vij = CI of singlet and triplet monoexcitations
c     (claculated by edirect and eindirect)
c     h1 = site of hole for i-th config
c     e1 = site of electron for i-th config
c     h2 = site of hole for j-th config
c     e2 = site of electron for j-th config
      include 'sizes.fi'
      include 'elstates.fi'
      include 'parameters.fi'
      integer h1, e1, h2, e2, i, j, k
      integer spinstate
      real(kind(1D0)) hel, hs(*)
      real*8 helem, edirect, eindirect,eJ,ho,eK
      call ZeroVar(hs,NQMAX2)
      k=0
      
      if(spinstate.ne.-1) then
         write(6,*)" in HSParseCalc, spinstate = ",spinstate
         if(debug) then
            write(6,*)" Generating DEBUG info on Hamiltonian Matrix elements -> fort.11"
            write(11,*)"  Hamiltonian Matrix Element "
            write(11,*)"  Spinstate = ",spinstate
            write(11,'(6a5,10a12)')"i","j","e1","h1","e2","h2","ho","J","K","H"
            write(11,*)"-----------------------------------------------------------"
         endif
      endif
      
      do i=1,nqdof
         h1=mic(1,i)
         e1=mic(2,i)
         do j=1,i
            h2=mic(1,j)
            e2=mic(2,j) 
            k=k+1
            
c     -- this part need to be checked ---
            ho = helem(h1,e1,h2,e2)
c     eJ = edirect(e1,h2,e2,h1)
c     eK = eindirect(e1,h2,h1,e2)
            eJ = edirect(h1,e1,h2,e2)
            eK = eindirect(h1,e1,h2,e2)
            hel = ho
            
            
c     if (spinstate.eq.1) hel=ho-edirect(e1,h2,e2,h1)
c     if (spinstate.eq.0) hel=ho-edirect(e1,h2,e2,h1)+2.0*eindirect(e1,h2,h1,e2)
            if(spinstate.eq.-1)hel = ho
c     if(spinstate.eq.0)hel = ho - eJ + 2.0*eK
c     if(spinstate.eq.1)hel = ho - eJ 
            if(spinstate.eq.0)hel = ho - eJ + 2.0*eK
            if(spinstate.eq.1)hel = ho - eJ 
c     -- ----------------------------------	      
            hs(k)=hel
            
            if(debug.and.spinstate.ne.-1)write(11,'(6i5,4f12.6)')i,j,e1,h1,e2,h2,ho,eJ,eK,hel
            
         end do
      end do  
      
      return
      end
c     
c     
      subroutine TransitionDipole(v,dip)
c     Subroutine that calculates the squared one-photon transition dipoles
      include 'sizes.fi'
      include 'elstates.fi'
      include 'coords.fi'
      real(kind(1D0)) v(NQMAX,*), dip(*),rnm
      integer nq, noff0, noff1, noff2
      real(kind(1D0)) dx, dy
      real(kind(1D0)) SSum
      real(kind(1D0)) TDip(2,3),y1,y2

c     TDip - components of the transition dipoles of the elementary
c     e-h pairs in the XY conjugation plane as obtained for PPV. 
c     The three 2D vectors corresponds to e-h pairs with separations 
c     zero, one and two.
c     See table 1 in JCP 118 4291 (2003).
      
c     v contains the CI eigenvectors
      
      
      integer m,n
      
      data TDip /-1.1612d0, 0.2285d0,
     $     0.5838d0, 0.0324d0,
     $     -0.1493d0, 0.0166d0 /
      pi = 3.14159d0
      noff0=nsites
      noff1=noff0+max((nsites-1)*2,0)
      noff2=noff1+max((nsites-2)*2,0)
      do nq=1,nqdof
         dx = 0.0d0
         dy = 0.0d0
c         y1 = coords(nq,2)
         do nq1 = 1,nqdof
            m = mic(2,nq1)
            n = mic(1,nq1)
            inm = iabs(n-m)
c            y2 = coords(nq1,2)
c        this is a functional approx to the numbers above.
            rnm = distmat(n,m)
c	      if(y1.eq.y2) then   ! consider only intramolecular contributions
            dx = dx + tdip(1,1)*exp(-0.8d0*rnm)*cos(Pi*rnm)*v(nq1,nq)
            dy = dy + tdip(2,1)*exp(-1.9d0)*rnm*v(nq1,nq)
c     endif
         enddo
         dip(nq) = dx*dx + dy*dy
c$$$  dx =  SSum(v(1,nq),noff0)*TDip(1,1)
c$$$  1	+SSum(v(noff0+1,nq),noff1-noff0)*TDip(1,2)
c$$$  2	+SSum(v(noff1+1,nq),noff2-noff1)*TDip(1,3)
c$$$  dy =  SSum(v(1,nq),noff0)*TDip(2,1)
c$$$  1	+SSum(v(noff0+1,nq),noff1-noff0)*TDip(2,2)
c$$$  2	+SSum(v(noff1+1,nq),noff2-noff1)*TDip(2,3)
c$$$  dip(nq)=(dx*dx+dy*dy)
      end do
c     110	format(i4,3f10.4)
      return
      end
      
