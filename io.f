c       io.f,v 1.2 2004/09/29 00:48:13 bittner Exp
c       2004/09/29 00:48:13
c       io.f,v
c       Revision 1.2  2004/09/29 00:48:13  bittner
c       Major changes from original code written by Stoyan.
c
c	io.for
c	input/outpute subroutines
c	WritePackedMatrix
c	WriteSparseMatrix
c	WriteMatrix
c	WriteVec
c	WriteMat
c	Output
c	WriteSpectra
c	getlowstate
c	writedtinfo
c       
	subroutine WritePackedMatrix(iw,a,n,mode)
	integer iw,n,mode, k, l
	integer lidx
	real(kind(1D0)) a(*)
	if (mode.le.0) then
	   do k=1,n
	      write(iw,102) k,(a(lidx(l,k)),l=1,k)
	   end do
	else
	   write(iw,103) (l,l=1,n)
	   do k=1,n
	      write(iw,102) k,(a(lidx(l,k)),l=1,n)
	   end do
	endif

 102	format(i6,2000f8.3)
 103	format(a6,2000f8.3)
 101	format(6x,2000i8)
 105	format(a,2f10.5)
	return
	end
c       
	subroutine WriteSparseMatrix(iw,a,n,nevec,threshold)
	integer, parameter :: icols=6
	integer iw, n, nevec
	real(kind(1D0)) a(nevec,n), threshold
	real(kind(1D0)) stack(icols)
	integer ik(icols), jk(icols), ncols, k, l, j, kk, ll
c       integer idx(*)
c       write(6,110) threshold
c       pause
 110	format('Treshold=',e12.5) 
	ncols=0
	do k=1,n
c       kk=idx(k)
	   kk=k
	   do l=1,n
c       ll=idx(l)
	      ll=l
	      if (dabs(a(kk,ll)).gt.threshold) then
		 ncols=ncols+1
		 stack(ncols)=a(kk,ll)
		 ik(ncols)=kk
		 jk(ncols)=ll
	      end if
	      if (ncols.ge.icols) then
		 write(iw, 101) (ik(j),jk(j),stack(j), j=1,ncols)
		 ncols=0
	      end if 
	   end do
	end do
 101	format(20(i5,i3,f10.4))
	if (ncols.ge.0) then
	   write(iw, 101) (ik(j),jk(j),stack(j), j=1,ncols)
	   ncols=0
	end if 
	return
	end
c       
c       
	subroutine WriteMatrix(iw,a,n,nevec)
	integer iw, n, nevec
	real(kind(1D0)) a(nevec,n)
	integer  k, l
	write(iw,103) (l,l=1,n)
	do k=1,n
	   write(iw,102) k,(a(k,l),l=1,n)
	end do
c  102	format(i4,2000f10.5)
 102	format(i4,2000e15.6)
 103	format(4x,2000i10)
	return
	end
c       
c       
	subroutine WriteVec(iw,a,n)
	integer iw, n, nevec
	real(kind(1D0)) a(n)
	integer l
	write(iw,102) (a(l),l=1,n)
 102	format(2000f10.5)
	return
	end
c       
c       
	subroutine WriteMat(iw,a,n,nevec)
	integer iw, n, nevec
	real(kind(1D0)) a(nevec,n)
	integer  k, l
	do k=1,n
	   write(iw,102) (a(k,l),l=1,n)
	end do
 102	format(2000e15.6)
	return
	end
c       
c       
	subroutine OutPut(iw,level,numstates)
c	Subroutine that writes data and results to exci.out
	include 'sizes.fi'
	include 'elstates.fi'
	include 'vibmodes.fi'
	include 'rates.fi'
	include 'constant.fi'
	include 'couplings.fi'
	include 'coords.fi'
	include 'parameters.fi'
	integer numstates
	integer iw, level, i, k, l, j
	integer getlowstate
	real(kind(1D0)) ehsym
	real(kind(1D0)) RadRate
	character*80 line

	write(iw,*)"Exciton Output File"
	write(iw,*)"Verbose level = ",level

	if (level.ge.4) then
	   write(iw,*) '   Normal Modes and Huang-Rhys factors'
	   write(iw,101) (i, i=1,ncdof)
	   write(iw, 103) '  E/eV',(qomega(i)*HBarMine, i=1,ncdof)
	   write(iw, 103) '      ', (Gkk(i,1)**2, i=1,ncdof) 
	   if (level.ge.5) then
	      do j=1,ncdof
		 write(iw,102) j, (qmode(j,i),i=1,ncdof)
	      end do
	   endif
	   write(iw,*)
	endif
	if (level.ge.2) then
	   k=getlowstate(0,0)
	   if (k.gt.0) then
	      write(iw,105)
	      write(iw,105) ' Singlets'
	      write(iw,103) ' EVert', (eigs(i),i=k,k+nqdof-1)
	      write(iw,103) ' EAdia', (eadiab(i),i=k,k+nqdof-1)
	      write(iw,103) ' Dipol', (transdip(i),i=k,k+nqdof-1)
	      write(iw,103) ' Rns-1', (RadRate(eigs(i),transdip(i))*1.0d03
	1 	   ,i=k,k+nqdof-1)
	   endif
	   k=getlowstate(1,0)
	   if (k.gt.0) then
	      write(iw,105)
	      write(iw,105) ' Triplets'
	      write(iw,103) ' EVert', (eigs(i),i=k,k+nqdof-1)
	      write(iw,103) ' EAdia', (eadiab(i),i=k,k+nqdof-1)
	   endif
	endif
	write(iw,105)
	if (level.ge.6.and.ispinstate.eq.0) then
	   write(iw,*) ' Huang-Rhys factors for transitions from electronic state n'
	   write(iw,*) 'n\        phonon mode      -> '
	   do j=1,ncdof
	      write(iw,102) j, (Gkk(j,i)**2 ,i=1,nqdof)
	   end do
	endif


	write(iw,103)
	if (level.ge.4) then
	   write(iw,*) ' Eigen Vectors'
	   do j=1,nqdof
	      write(iw,102) j, (h(j,i),i=1,nqdof)
	   end do
	   write(iw,*)
	endif

C 	if (level.ge.3) then
C 	   k=getlowstate(0,1)
C 	   if (k.gt.0) then
C 	      write(iw,105) ' Sxt', ehsym(h(1,k)), eigs(k)
C 	      call EHDensity(iw,h(1,k))
C 	   endif
C 	   k=getlowstate(0,-1)
C 	   if (k.gt.0) then
C 	      write(iw,105) ' Sct', ehsym(h(1,k)), eigs(k)
C 	      call EHDensity(iw,h(1,k))
C 	   endif
C 	   k=getlowstate(1,1)
C 	   if (k.gt.0) then
C 	      write(iw,105) ' Txt', ehsym(h(1,k)), eigs(k)
C 	      call EHDensity(iw,h(1,k))
C 	   endif
C 	   k=getlowstate(1,-1)
C 	   if (k.gt.0) then
C 	      write(iw,105) ' Tct', ehsym(h(1,k)), eigs(k)
C 	      call EHDensity(iw,h(1,k))
C 	   endif
C 	endif
	


	if (level.ge.2) then 
	   write(iw,*) ' Rate constants'
	   call WriteMatrix(iw,rk,nqvars,NQMAX*2)
	endif

	if(mctdh) then 

	   write(6,*)" MCTDH operator input data --> exci-mctdh.op "
	   open(unit = 21,file = "exci-mctdh.op")
	   open(unit =22, file = "Hamiltonian_op.inp")




	   write(21,*) 'Hamiltonian Operator input for MCTDH Calculation'
	   write(21,*) '#------------cut here ----------------------------'
	   write(21,*) 'OP_DEFINE-SECTION'
	   write(21,*) 'title'
	   write(21,*) 'input generated by exciton.x for MCTDH'
	   write(21,*) '    '
	   write(21,*) 'end-title'
	   write(21,*) 'end-op_define-section'
	   write(21,*) '    '
	   write(21,*) 'PARAMETER-SECTION'

	   write(22,'("c------- Hamiltonian Terms for Surface Hopping Calc----")')
	   write(22,'("       nc = ",i3)')nckeep
	   write(22,'("       nq = ",i3)')nqkeep
	   write(22,'("c-- Parameters  ")')
	   write(22,'("c    phonon frequencies (in rad/ps) ")')
	   do kx = 1,nckeep
	      write(22,'("              w(",i2,")  = ",f12.6)')kx,qomega(kx)
	   enddo
	   write(22,'("c    vertical energies (in eV)")')
	   do i=1,nqkeep
	      write(22,'("              e(",i2,")  = ",f12.6)')i,eigs(i)
	   enddo
	   write(22,'("c  gab(xi) couplings-- these are hder elements in the diab. rep")')
	   
c need to generate unique mode names
c note, that the qomega's are in rad/ps
c need to convert to ev using hbar =  6.58*10^-4

c  also, I've explicitly hardwired here and in relaxation.f
c  a cutoff for only printing the lowest 9 electronic states 
c the mctdh input generation part. 

c==eigs(i) should be the adiabatic energies.


c    for the mctdh calculation the gab elements are in dimensionless units as
c    per given in May and Kuhn's book.

	   do kx = 1,min(nckeep,9)
	      write(21,'("w",i1,"  = ",f12.6 , ", ev")')kx,qomega(kx)*6.58E-4
	   enddo
	   do kx = 10,min(nckeep,99)
	      write(21,'("w",i2,"  = ",f12.6 , ", ev")')kx,qomega(kx)*6.58E-4
	   enddo
	   do i = 1,min(nqkeep,9)
	      write(21,'("e",i1,"  = ",f12.6 , ", ev")')i,eigs(i)
	   enddo
c	   do i = 10,min(nqkeep,99)
c	      write(21,'("e",i2,"  = ",f12.6 , ", ev")')i,eigs(i)
c	   enddo
	   open(unit=12,file="parameters.mctdh")
	   open(unit=13,file="parameters.hopping")
 30	   continue
	   read(12,'(A80)',end=40)line
	   write(21,*)line
 40	   continue	   
	   read(13,'(A80)',end=45)line
	   write(22,*)line
	   goto 30
 45	   continue
	   close(12)
	   close(22)
	   

	   write(21,*) 'end-parameter-section'
	   write(21,*)
   	   write(21,*) 'hamiltonian-section'
	   write(21,'("modes  |  el")')
	   do i = 1,min(9,nckeep)
	      write(21,'("modes   |  v",i1)')i
	   enddo
	   do i = 10,min(99,nckeep)
	      write(21,'("modes   |  v",i2)')i
	   enddo
	   do i = 1,min(9,nqkeep)
	      write(21,'("1.0*e",i1,10x,"|1  S",i1,"&",i1)')i,i,i       
	   enddo
	   do i=1,min(9,nckeep)
	      write(21,'("0.5*w",i1,10x,"|",i1,"   q^2")')i,i+1
	      write(21,'("0.5*w",i1,10x,"|",i1,"   p^2")')i,i+1
	   enddo
	   do i=10,min(99,nckeep)
	      write(21,'("0.5*w",i2,10x,"|",i2,"   q^2")')i,i+1
	      write(21,'("0.5*w",i2,10x,"|",i2,"   p^2")')i,i+1
	   enddo

	   open(unit=13,file ="couplings.mctdh")
 10	   continue
	   read(13,'(A80)',end=20)line
	   write(21,*)line
	   goto 10
 20	   continue
	   close(13)

c also need to include the energy shifts by completing the square in the 
c phonon displacement
	   
   	   write(21,*) 'end-hamiltonian-section'
	endif
	close(21)

 102	format(i6,2000f8.3)
 103	format(a6,2000f8.3)
 101	format(6x,2000i8)
 105	format(a,2f10.5)
	return
	end
c       
c       
	subroutine WriteSpectra(iw)
	include 'sizes.fi'
	include 'spectra.fi'
	integer iw, i, it
	real(kind(1D0)) e, EofIndex
 102	format(1000es15.6)
	do i=1, nspec
	   e = EofIndex(i,eleft,eright,nspec)
c       e = eleft+(eright-eleft)*dble(i-1)/dble(nspec-1)
	   write(iw,102) e,(spectrum(i,it),it=1,ntime)
	end do 
	return
	end
c       
c       
	integer function getlowstate(spin,ehpar)
	include 'sizes.fi'
	include 'elstates.fi'
	integer spin,ehpar
	integer nstart, i
	nstart=0
	getlowstate=0
	if (ispinstate.lt.0) then
	   if (spin.eq.1) nstart=nqdof
	else
	   if (spin.ne.ispinstate) return
	endif
	do i=nstart+1,nstart+nqdof
	   if (iehpar(i)*ehpar.ge.0) then
	      getlowstate=i
	      return
	   endif
	end do
	return
	end
c       
c       
	subroutine writedtinfo(iw)
	include 'genparam.fi'
	integer i
	write(iw,100)	''
	write(iw,100)	'Start:  ',
	1   sdt(5),':',
	2   sdt(6),':',
	3   sdt(7),' ',
	4   sdt(3),'/',
	5   sdt(2),'/',
	6   sdt(1),'(hh:mm:ss DD:MM:YYYY)'
	write(iw,100)	'End:  ',
	1   edt(5),':',
	2   edt(6),':',
	3   edt(7),' ',
	4   edt(3),'/',
	5   edt(2),'/',
	6   edt(1),'(hh:mm:ss DD:MM:YYYY)'
	do i=1,8 
	   jdt(i) = edt(i)-sdt(i)
	end do
c	Determine seconds elapsed (index 7)
	if (jdt(7).lt.0) then
	   jdt(7) = 60+jdt(7)
	   jdt(6) = jdt(6)-1
	endif
c	Determine minutes elapsed (index 6)
	if (jdt(6).lt.0) then
	   jdt(6) = 60+jdt(6)
	   jdt(5) = jdt(5)-1
	endif
c	Determine hours elapsed (index 5)
	if (jdt(5).lt.0) then
	   jdt(5) = 24+jdt(5)
	   jdt(3) = jdt(3)-1
	endif
c	Determine days elapsed (index 3)
	if (jdt(3).lt.0) then
	   select case (sdt(3))
	   case (JAN,MAR,MAY,JUL,AUG,OCT,DEC);
	      jdt(3)=31+jdt(3)
	   case (APR,JUN,SEP,NOV);
	      jdt(3)=31+jdt(3)
	   case (FEB);
	      if (mod(sdt(1),4).eq.0) then
		 jdt(3)=29+jdt(3)
	      else
		 jdt(3)=28+jdt(3)
	      endif
	   end select
	   jdt(2) = jdt(2)-1
	endif
c	Determine months elapsed (index 2)
	if (jdt(2).lt.0) then
	   jdt(2) = 12+jdt(2)
	   jdt(1) = jdt(1)-1
	endif
	write(iw,100)	'Job Time:  '
	1   , jdt(5),':'
	2   , jdt(6),':'
	3   , jdt(7),' '
	4   , jdt(3),'/'
	5   , jdt(2),'/'
	6   , jdt(1),'(hh:mm:ss DD:MM:YYYY)'
 100	format(a15,5(i3,a1),i4,a25)
	end
c       
	
