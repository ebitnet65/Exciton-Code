c Copyright (C) 2004  Eric R. Bittner
c
c   This program is free software; you can redistribute it and/or
c   modify it under the terms of the GNU General Public License
c   as published by the Free Software Foundation; either version 2
c   of the License, or (at your option) any later version.

c   This program is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c   GNU General Public License for more details.
c
c   You should have received a copy of the GNU General Public License
c   along with this program; if not, write to the Free Software
c   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
c       
c       input.f,v 1.2 2004/09/29 00:48:13 bittner Exp
c       2004/09/29 00:48:13
c       input.f,v
c       Revision 1.2  2004/09/29 00:48:13  bittner
c       Major changes from original code written by Stoyan.
c       
c       Revision 1.1  2004/09/24 23:05:58  bittner
c       Initial revision
c       
c	Input.for
c
	subroutine Input
c
c	something about input parameters
c       ispinstate=0 - singlet; > - triplet
c       
c	Electron model
c	egap:	e-level (the energy difference of the band centers / 2)
c       h-level = -e-level
c	ewidth: the intersite transfer parameter (negative) for e-levels
c       the band dispersion is 2*ewidth*cos(k), so than the bandwidth
c       is actually 4*ewidth
c       
c	Phonon model 
c	forcek: on-site spring constant: U(x)=0.5*forcek*x^2
c	forceo: off-site spring constant: U(x1-x2)=0.5*forceo*(x1-x2)^2
c	forcem: reduced mass of the on-site oscillators
c       
c	Electron-phonon coupling
c	fediag: on-site electron-phonon coupling constant
c	On-site one-electron levels are modulated by
c       e(x) = e-level - fediag*x (+/-:occupied/unoccupied)
c       
c	Electron-electron interaction 
c	Coulomb interaction: diagonal in site representation
c       U(h,e)= -coulombe(1+|h-e|/coulombr)^-1
c       where h and e are the position of the hole and electron
c	Dipole-dipole interaction: effective for excitonic states (h=e=m)
c       D(m1,m2)= dipole/(1+|m1-m2|^3/dipolr^3)^-1  
c       
c	Envionment
c	EKbT - thermodynamic temparature in energy units
c	PhGama - phonon line-width 
c	ElGama - electron line width
c       
	include 'sizes.fi' 
	include 'parameters.fi'
	include 'constant.fi'
	include 'structure.fi'
	include 'coords.fi'
        include 'elstates.fi'
c       
c	File Descriptors
	character*100 keywrd
	integer count,i,j
	integer jbonds(4)
	real(kind(1D0)) ave,xx,yy,x1,x2,y1,y2
	CHARACTER LINE*80,LINEOLD*80,basefilename
	logical nmlist,disorder,help

	namelist/electrons/nsites,ispinstate,maxrelax,maxits,iex
	namelist/keywds/disorder,mctdh,relax,evolve,exciton,
	1    chargesep,absorb,pdensity,radrates,debug,iolevel
	namelist/params/sigma_beta,coulombe,coulombr,dipole,dipolr,potexe,potexr,EKbT,PhGama,ElGama,knonrad
	namelist/phonons/forcek1,forcek2,forceo1,forceo2,gediag1,gediag2,forcem
	namelist/mctdh_control/nqkeep,nckeep
	namelist/evolve_control/timall,timstep,timshot


c       
c	Parameters initialization
c       234567
C       EKbT=298.0d0*BoltzmanKMine
C       PhGama=0.3d0
C       ElGama=0.05d0
C       potexe=0.0d0
C       potexr=1.0d0
C       forcek2=0.0d0
C       forceo2=0.0d0
C       fediag2=0.0d0
C       gediag2=0.0d0
C       gediag1=0.0d0
C       mcyclic=0
C      key.relax = .false.
C       key.absorb = .false.
C       key.evolve = .false.
C       key.pdensity = .true.
C       key.pemission = .false.
C       key.radrates = .false.
C       key.knonrad=0.0d0
C       
	
c       default parameters for electron and phonon models





c23456	
c       - parameters from table 1 in JCP 118 4291 (2003)
	
	egap = 2.75
	ewidth = 0.5360
	forcem = 1.1d0
	
	forcek1 = 0.2
	forcek2 = 0.0125
	forceo1 = 0.02
	forceo2 = 0.0013
	gediag1 = 0.60
	gediag2 = 4.0
	coulombe = 3.0921
	coulombr = 0.6840
	dipole = -0.0321
	dipolr = 1.0
	potexe = 1.0574
	potexr = 0.4743
	EKbT  =     0.0250   
	PhGama =   0.0100
	ElGama=    0.0300

c default parameters
c  this is for a simple 4-site ppv model	

c$$$	timall = 1000.0
c$$$	timstep = 1.0
c$$$	tmshot = 20.0
c$$$	sigma_beta = 0.0
c$$$	DISORDER =  .false.
c$$$	MCTDH =  .False.
c$$$	RELAX =  .true.
c$$$	EVOLVE =  .false.
c$$$	EXCITON =  .false.
c$$$	CHARGESEP =  .false.
c$$$	ABSORB =  .false.
c$$$	PDENSITY =  .false.
c$$$	RADRATES =  .false.
c$$$	NSITES =   4
c$$$	ISPINSTATE =   0
c$$$	MAXRELAX =   16
c$$$	MAXITS =   100
c$$$	IEX =   1
c$$$	IOLEVEL = 3
c$$$	debug = .false.
c$$$
	
c23456

	do count =1,NSMAX
	   call ZeroSite(count)
	end do
	
	
C       ==--------------------------------------------------------------==
C       ==  OPEN INPUT FILE (THE FILE NAME IS GIVEN ON THE COMMAND LINE)==
C       ==--------------------------------------------------------------==
	ICARG=IARGC()
	IF(ICARG.LT.1) THEN
	   WRITE(6,*) ' NO INPUT FILE NAME SPECIFIED '
	   CALL STOPGM('CONTROL',' ')
	ENDIF
	CALL GETARG(1,line)
	if(index(line,"-generate_input").ne.0) then
	   write(6,*)"#-----BEGIN SAMPLE INPUT FOR 4-SITE MODEL------"
	   write(6,nml=electrons)
	   write(6,nml=keywds)
	   write(6,nml=params)
	   write(6,nml=evolve_control)
	   write(6,nml=mctdh_control)

           write(6,*)"    1  -2.7500 2.7500 0.0   0.0"
           write(6,*)"    2  -2.7500 2.7500 1.0   0.0"
	   write(6,*)"    3  -2.7500 2.7500 2.0   0.0"
	   write(6,*)"    4  -2.7500 2.7500 3.0   0.0"
	   write(6,*)"    1 2  0.536  1.0  1.0"
	   write(6,*)"    2 3  0.536  1.0  1.0"
	   write(6,*)"    3 4  0.536  1.0  1.0"
	   write(6,*)"    999 999 1.0 1.0 1.0"
	   write(6,*)"#-----END SAMPLE INPUT FOR 4-SITE MODEL------"
	   call stoppg(0)
	endif
	if(index(line,"-help").ne.0) then
	   help = .true.
	   write(6,*)"#______DEFAULT NAMELIST INPUT BLOCK______"
	   write(6,*)"#Copy and paste this into your own input file"
	   write(6,nml=electrons)
	   write(6,nml=keywds)
	   write(6,nml=params)
	   write(6,nml=evolve_control)
	   write(6,nml=mctdh_control)

	   write(6,*)"#---This is then followed by the site energies and hopping integrals---"
	   write(6,*)"#---For example for a 4 site model"
           write(6,*)"    1  -2.7500 2.7500 0.0   0.0"
           write(6,*)"    2   -2.7500 2.7500 1.0   0.0"
	   write(6,*)"    3  -2.7500 2.7500 2.0   0.0"
	   write(6,*)"    4  -2.7500 2.7500 3.0   0.0"
	   write(6,*)"    1 2  0.536  1.0 1.0"
	   write(6,*)"    2 3  0.536  1.0 1.0"
	   write(6,*)"    3 4  0.536  1.0 1.0"
	   write(6,*)"    999 999 1.0 1.0 1.0"
	   write(6,*)"#---------------"
	   write(6,*)" note: the first 4 lines are the site numbers, site energies for h and e and the x/y coord."
	   write(6,*)" the second set are the hopping integrals beween sites i & j.  This is terminated with 999 999 1"
	   write(6,*)" to indicate the end of the input info."
	   write(6,*)" " 
	   call information(6)
	   call stoppg(0)
	endif
	write(6,'(" Reading from File name:   ",A80)')line
	OPEN(UNIT=5,FILE=LINE,STATUS='OLD')

	iline = 0

C 	read(UNIT=5,fmt='(a80)')line
C 	if(index(line,"COMMENT_BEGIN").ne.0)then
C 	   write(6,'(a80)')line
C  11	   continue
C 	   read(5,fmt='(a80)')line
C 	   write(6,'(a80)')line
C 	   if(index(line,"COMMENT_END").eq.0)goto 11
C 	endif

C First, a series of namelists for the input options
c This replaces the key-word parsing 
	write(6,*)"#  reading NML input"
 90	read(5,nml=electrons)
 91	read(5,nml=keywds)
 92	read(5,nml=params)

 94	read(5,nml=evolve_control)
 95	read(5,nml=mctdh_control)

 96	continue

	write(6,nml=electrons)
	write(6,nml=keywds)
	write(6,nml=params)
	write(6,nml=phonons)
	write(6,nml=evolve_control)
	write(6,nml=mctdh_control)


	write(6,*)"#----NML info read"




C once the namelists have been parsed, input the site model 
	write(6,*)"#reading in site energies"
	do count = 1,nsites
	   read(5,*)i,eh,ee,xx,yy
	   eelec(count)=ee
	   ehole(count)=eh
	   coords(i,1) =xx
	   coords(i,2) =yy
c	   coords(i,3) =zz
	   write(6,'(i12,5f12.6)')count,eelec(count),ehole(count),xx,yy
	enddo
c       after the site energies, we read in the connectivity 
	write(6,*)"#finished reading in site energies"

	call zerovar(thop,nsmax*nsmax)
	call zerovar(tbond1,nsmax*nsmax)
	call zerovar(tbond2,nsmax*nsmax)

c--additional option to have phonons couple 
c  differently

	do while(.true.)
	   read(5,*,end=46)i,j,ee,tb1,tb2
	   write(6,*)i,j,ee,tb1,tb2
	   if(i.ne.999)then 
	      thop(i,j) = ee
	      thop(j,i) = ee
	      tbond1(i,j) = tb1
	      tbond1(j,i) = tb1
	      tbond2(j,i) = tb2
	      tbond2(i,j) = tb2
	   endif
	enddo
 46     Continue
	iw = 4
	if(debug.or.(nsites.lt.12)) then
	   iw = 6
	   write(iw,*)"huckle matrix"
	   do i = 1,nsites
	      write(iw,'(12f12.6)')(thop(i,j),j=1,nsites)
	   enddo
	endif


	
c       
c       generate distance matrix
	write(6,*)"Generating Distance Matrix"
	do i = 1,nsites
	   x1 = coords(i,1)
	   y1 = coords(i,2)
c	   z1 = coords(i,3)
	   do j = i,nsites
	      x2 = coords(j,1)
	      y2 = coords(j,2)
c	      y2 = coords(j,3)
	      distmat(i,j) = dsqrt((x1-x2)**2 + (y1-y2)**2 )
	      distmat(j,i) = distmat(i,j)
	   enddo
	enddo

	
c       nqkeep = min(nqkeep,nqdof)
c       nckeep = min(nckeep,ncdof)

	   
c bail outs
	if((evolve).and.(( exciton.or. chargesep).eqv.(.false.)))then
	   write(6,*)"-----"
	   write(6,*)"EVOLVE keyword selected"
	   write(6,*)" evolve = ", evolve
	   write(6,*)" exciton = ", exciton
	   write(6,*)" chargesep = ", chargesep
	   write(6,*)"In order to perform this calculation, "
	   write(6,*)"you need to select your initial configuration"
	   write(6,*)"EXCITON"
	   write(6,*)"NN  <----  initial state"
	   write(6,*)"or use the ``CHARGE_SEPARATED_STATE'' "
	   write(6,*)"keyword to inject an electron and hole on opposite"
	   write(6,*)"ends of the molecule."
	   write(6,*)" nontheless..."
	   write(6,*)"...you've been terminated!!!"
	   stop
	endif
	if( relax.and.maxrelax.eq.0) maxrelax = nqdof

	
	write(6,*)"FINISHED READING INPUT "
	write(6,'("READ ",i3," INPUT FLAGS")')iline
	close(5)
	if(forcek1.eq.0.0d0) then
	   write(6,*)"problem reading PHONON_MODEL section"
	   print *,forcek1
	   stop
	endif
c       - parameters from table 1 in JCP 118 4291 (2003)
	
	egap = 2.75
	ewidth = 0.5360
	forcem = 1.1d0
	
	forcek1 = 0.2
	forcek2 = 0.0125
	forceo1 = 0.02
	forceo2 = 0.0013
c- huang rhys parameters
	gediag1 = 0.60
	gediag2 = 4.0
c
	coulombe = 3.0921
	coulombr = 0.6840
	dipole = -0.0321
	dipolr = 1.0
	potexe = 1.0574
	potexr = 0.4743
	EKbT  =     0.0250   
	PhGama =   0.0100
	ElGama=    0.0300


	
	
	write(6,111)	'nsites:',nsites,'ispinst:',ispinstate,'mcyclic:', mcyclic
	write(6,122)	'egap:',egap,'ewidth:',ewidth,'forcem:',forcem
	write(6,122)'forcek1:',forcek1,'forceo1:',forceo1,'gediag1:',gediag1
	write(6,122)'forcek2:',forcek2,'forceo2:',forceo2,'gediag2:',gediag2
	write(6,122)	'coulombe:',coulombe,'coulombr:',coulombr
	write(6,122)	'dipole:',dipole,'dipolr:',dipolr
	write(6,122)	'potexe:',potexe,'potexr:',potexr
	write(6,122)	'EKbT:',EKbT,'PhGama:',PhGama,'ElGama:',ElGama
	write(6,122)	'timall:',timall,'timstep:',timstep,'timeshot:',timshot
     
	fediag1=sqrt(0.5d0*forcek1**3*gediag1/(HBarMine**2))
	if (gediag2.gt.0.0d0) then
	   fediag2=sqrt(0.5*forcek2**3*gediag2/(HBarMine**2))
	else
	   fediag2=0.0d0
	endif
	ncdof=nsites
        
c	Write output for site energies 
	
	write(6,*)"---- SITE ENERGIES FOR ELECTRONS AND HOLES ---- "
	if(disorder) then
	   write(6,*)"Off diagonal disorder will be imposed"
	   write(6,*)"   sigma  = ",sigma_beta

c--idum simply must be some large negative number to initiate the sequence. 	   
	   idum = -9999
	   do i = 1,nsites
	      do j = i,nsites
		 thop(i,j) = thop(i,j)+ gasdev(idum)*sigma_beta
		 thop(j,i) = thop(i,j)
	      enddo
	   enddo
	endif

	if(nsites.le.10) then 
	   write(6,*)"connectivity matrix 1"
	   do i = 1,nsites
	      write(6,'(10f12.0)')(tbond1(i,j),j=1,nsites)
	   enddo
	endif

	write(6,'(a12,a12,a12)')"site #","e-level","h-level"
	do count = 1,nsites
	   write(6,'(i12,2f12.6)')count,eelec(count),ehole(count)
	enddo
	


	write(6,*)"------------------------------------------------"
	write(6,nml=keywds)
        write(6,*)"  KEY WORDS SELECTED" 
	write(6,*)" relax = ", relax
	write(6,*)" evolve = ", evolve
	write(6,*)" pdensity = ", pdensity
	write(6,*)" radrates = ", radrates
	write(6,*)"------------------------------------------------"
	
	if( evolve)then
	   write(6,*)" proceeding with density matrix evaluation" 
	   write(6,'(3f12.6)')timall, timstep, timshot 
	endif
c       
 100	continue
c       
 101	format(10i5)
 102	format(10f10.4)   
 103	format(a100)
 104	format(a10,1000f8.4)
 105	format(a10,1000i8)
 111	format(10(a10,i10))
 122	format(5(a10,f10.4))
 133	format(a10,a100)
	
	
	
	return
 99	continue
	WRITE(6,*) ' CONTROL: ERROR IN READING INPUT FILE'
	WRITE(6,*) ' THE LAST TWO LINES READ WERE '
	WRITE(6,*) LINEOLD
	WRITE(6,*) LINE
	call stoppg(0)

 200	continue
	write(6,nml=electrons)
	write(6,nml=keywds)
	write(6,nml=params)
	write(6,nml=phonons)
	write(6,nml=evolve_control)
	write(6,nml=mctdh_control)
	call stoppg(0)
	
	end
	
	
C       @(#)input.f,v 1.2 2004/09/29 00:48:13 bittner Exp 
	
C       ==================================================================
	FUNCTION INSCAN(IUNIT,LABEL)
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*(*) LABEL
	CHARACTER LINE*80
	LOGICAL EXISTS
C       ==--------------------------------------------------------------==
C       ==  Scans unit 'iunit' for the section header 'label'           ==
C       ==--------------------------------------------------------------==
	INQUIRE(UNIT=IUNIT,EXIST=EXISTS)
	IF(.NOT.EXISTS) THEN
	   WRITE(6,'(A,I3,A)') ' INSCAN: Unit',iunit,' does not exist'
	   CALL STOPGM('INSCAN',' ')
	ENDIF
	INSCAN=1
	REWIND(IUNIT)
C       ==--------------------------------------------------------------==
	CALL XSTRING(LABEL,IA,IE)
 10	CONTINUE
	READ(IUNIT,END=20,ERR=20,FMT='(A80)') LINE
	IFIND=INDEX(LINE,LABEL(IA:IE))
	IF(IFIND.NE.0) THEN
	   INSCAN=0
	   GOTO 20
	ENDIF
	GOTO 10
 20	CONTINUE
C       ==--------------------------------------------------------------==
	RETURN
	END
	
	
	
C       @(#)input.f,v 1.2 2004/09/29 00:48:13 bittner Exp 
	
C       ==================================================================
	SUBROUTINE STOPGM(A,B)
	CHARACTER*(*) A,B
	WRITE(6,'(//,A,A,A,A)') ' PROGRAM STOPS IN SUBROUTINE ',A,'| ',B
	
	
	STOP 999
	
	END
C       ==================================================================
	
	
C       ==================================================================
	SUBROUTINE XSTRING(STRING,IA,IE)
	CHARACTER STRING*(*)
	SLEN=LEN(STRING)
	IA=1
	DO I=1,SLEN
	   IF(STRING(I:I).NE.' ') THEN
	      IA=I
	      GOTO 10
	   ENDIF
	ENDDO
 10	CONTINUE
	DO I=IA,SLEN
	   IF(STRING(I:I).EQ.' ') THEN
	      IE=I-1
	      GOTO 20
	   ENDIF
	ENDDO
	IE=SLEN
 20	CONTINUE
	RETURN
	END
C       ==================================================================

	FUNCTION RAN1(IDUM0)
	implicit real*8(a-h,o-z)
	real*8 R(97)
        parameter (m1=259200,ia1=7147,ic1=54773,rm1=1./m1)
	parameter (m2=134456,ia2=8182,ic2=28411,rm2=1./m2)
	parameter (m3=243000,ia3=4561,ic3=51349)
	data iff/0/
	if(idum.lt.0.or.iff.eq.0) then
	   iff=1
	   ix1=mod(ic1-idum,m1)
	   ix1=mod(ia1*ix1+ic1,m1)
	   ix2=mod(ix1,m2)
	   ix1=mod(ia1*ix1+ic1,m1)
	   ix3=mod(ix1,m3)
	   do 11 j = 1,97
	      ix1=mod(ia1*ix1+ic1,m1)
	      ix2=mod(ia2*ix2+ic2,m2)
	      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
 11	   continue
	   idum=1
	endif
	ix1=mod(ia1*ix1+ic1,m1)
	ix2=mod(ia2*ix2+ic2,m2)
	ix3=mod(ia3*ix3+ic3,m3)
	j = 1+(97*ix3)/m3
	if(j.gt.97.or.j.lt.1) pause
	ran1 = r(j)
	r(j) = (float(ix1)+float(ix2)*rm2)*rm1
	return
	end
	function gasdev(idum)
	implicit real*8(a-h,o-z)
	data iset/0/
	if(iset.eq.0) then
 1	   v1 = 2.*ran1(idum) - 1
	   v2 = 2.*ran1(idum) - 1
	   r = v1*v1 + v2*v2
	   if(r.ge.1.0.or.r.eq.0.0)go to 1
	   fac = sqrt(-2.*log(r)/r)
	   gset = v*fac
	   gasdev = v2*fac
	   iset = 1
	else
	   gasdev = gset
	   iset = 0
	endif
	return
	end
	
	
	
	subroutine information(iw)
	integer iw
	write(6,*)
	write(6,*)
	write(6,*)
	write(6,*)"  Keyword List and default Values for Exciton Code"
	write(6,*)"---------------------------------------------------------"
	write(6,*)"NML = ELECTRONS "
	write(6,*)"   NSITES      = TOTAL NUMBER OF WANNIER FUNCTIONS"
	write(6,*)"   ISPINSTATE  = Spin state (0= singlet, 1  = triplet)"
	write(6,*)"   MAXRELAX    = Total number of states to relax during"
	write(6,*)"                 adiabatic relaxation step"
	write(6,*)"                 MAXRELAX < NSITES*NSITES"
	write(6,*)"   MAXITS      = Max number of scf iterations during adiabatic relaxation"
	write(6,*)"   IEX         = Initial CI eigenstate for relaxation calculation"
	write(6,*)"NML = KEYWDS"
	write(6,*)"   DISORDER    = T/F  Introduce off-diagonal disorder"
	write(6,*)"   MCTDH       = T/F  Generate input Hamiltonian for MCTDH calculation later"
	write(6,*)"   EVOLVE      = T/F  Evolve Pauli Master Eq. starting from either a CT state"
	write(6,*)"                 where the e & h are located on the end sites of the chain, or " 
	write(6,*)"                 from the IEX-th CI state allowing one to start from an excitonic "
	write(6,*)"                 eigenstate of the CI hamiltonian"
	write(6,*)"   RELAX       = T/F  Compute adiabatic states by relaxing the lattice."
	write(6,*)"   EXCITON     = T/F  If EVOLVE=T, then start off in the IEX state."
	write(6,*)"   CHARGESEP   = T/F  If EVOLVE=T, start off in a CT state. "
	write(6,*)"  note:  If EVOLVE=T then one and only one of EXCITON or CHARGESEP may be T as well."
	write(6,*)"   ABSORB      = T/F  Compute Absorption spectrum (output to ________)"
	write(6,*)"   PDENSITY    = T/F  Output density"
	write(6,*)"   RADRATES    = T/F  Compute Radiative Rates."
	write(6,*)"NML=PARAMS"
	write(6,*)"   SIGMA_BETA  = fluctuations in hopping integral (if DISORDER = T)"
	write(6,*)"   COULOMBE & COULOMBR :  parameters for direct Coulomb interaction"
	write(6,*)"   DIPOLE  & DIPOLR   :  parameters for dipole-dipole interaction"
	write(6,*)"   POTEXE  & POTEXR    :  parameters for exchange interaction"
	write(6,*)"   EKBT = Boltzmann constant"
	write(6,*)"   PhGama = phonon linewidth"
	write(6,*)"   ElGama = Spectral linewidth"
	write(6,*)"   KNONRAD = Non-radiative rate"
	write(6,*)"NML=PHONONS"
	write(6,*)"   FORCEK1 & FORCEK2:  force constants for 2 phonon branches"
	write(6,*)"   FORCEO1 & FORCEO2:  "
	write(6,*)"   GEDIAG1 & GEDIAG2: electron/phonon couplings"
	write(6,*)"NML=MCTDH_CONTROL"
	write(6,*)"   NQKEEP & NCKEEP : number of quantum and classical modes to "
	write(6,*)"          consider in generating MCTDH input Hamiltonian"
	write(6,*)"NML=EVOLVE_CONTROL"
	write(6,*)"   TIMALL  = Maximum time in ps"
	write(6,*)"   TIMSTEP = time-step (ps)"
	write(6,*)"   TMSHOT  = output interval"
	write(6,*)"-------------------------------------------------"
	return
	end
	
	

