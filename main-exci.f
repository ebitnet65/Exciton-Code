c Copyright (C) 2004  Eric R. Bittner
c
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
c
c     exci.f,v 1.2 2004/09/25 02:42:12 bittner Exp
c     2004/09/25 02:42:12
c     exci.f,v
c     Revision 1.2  2004/09/25 02:42:12  bittner
c     Input part works.
c     Crashes when we generate the redfield tensor
c         solution is to set the stacksize  to unlimited
c     
c     Revision 1.1  2004/09/24 23:07:54  bittner
c     Initial revision
c     
c     exci.for
c     Exci
c     CalcRedfield
c     TimeAverage
c     EHInjection
c     WriteEvolution
c     
      program Exci
      
      
c     Input dimensions
      include 'sizes.fi' 
      
c     Other input parameters
      include 'parameters.fi'
      
c     Electonic eigenstates
      include 'elstates.fi'
      
c     Vibrational coordinates and eigenstate populations
      include 'coords.fi'
      
c     Storage of different and time-dependent spectral distributions,
c     subject to output
      include 'spectra.fi'
      
c     El-ph coupling constants
      include 'couplings.fi'
      
c     Vibrational frequencies and modes
      include 'vibmodes.fi'
      
c     Basic units and constants
      include 'constant.fi'
      
c     Interstate transition rates and Redfield matrix
      include 'rates.fi'
      
      
c     General Parameters
      include 'genparam.fi'

c     
      integer nsteps, nshots, i,level
      
      character*10 user,hostname,cwd
      integer procid
      logical locked
      
      locked = .false.
      inquire(file="RUNNING",exist=locked)
      open(999,file="RUNNING")
c      if(locked)call stoppg(2)


      call Header

c     Determine machine limit on precision (eigval.for)



      call GetMachEps

c     Get Start Time
      
      call date_and_time(sdat,stim,szon,sdt)
      
      write(6,*)"calling input"
      
      call Input
      
      
      
c     Setup the energy levels:
c     Calculate vibrational modes, diabatic and adiabatic eigenstates
      
      
      write(6,*)"calling Prepare"
      

      call Prepare

      open(unit=10,file = "diabatic_densities.out")
      if(maxrelax.eq.0)maxrelax = 20
      do i=1,min(maxrelax,nqdof)
         write(10,*)"Eigenstate ",i
         write(10,*)"Diabatic Energy = ",eigs(i)
         call EHDensity(10,h(1,i))
      enddo
      close(10)
      
      write(6,*)"We are prepared...."
      
c     Determine the spectral range and the times at which
c     time-dependent quantities and distributions will be stored for 
c     subsequent output.
      
      nspec=NSPECMAX
      eleft=min(eigs(1)-1.0d0,1.0d0)
      eright=min(eigs(nqdof)+0.5d0,10.0d0)
      nqvars=nqexc+1
c     Calculate interstate rates and construct the Redfield matrix

c      if(( mctdh).and.( evolve.eq..false.)) then
c         write(6,*)"Computing  couplings for later use in MCTDH code"
c         write(6,*)"the rest is up to you"
c         call CalcRedfield    
c      endif



c-----  master eq. relaxation module --------

      if (evolve) then
         surfacehop = .false.
         write(6,*)"Computing rates and constructing Redfield matrix"
c         write(6,*)" OK to here"
         
         call CalcRedfield    
c         write(6,*)" OK to here"
c         if(debug.and.(nqdof.le.16)) then 
c         write(6,*)" Writing  REDFIELD TENSOR to External File"
c         write(6,*)" redfield.txt "
c         open(101,file="redfield.txt")
c         write(6,'(5x,16(i5,7x))')(i,i=1,nqdof)
c         do j = 1,nqdof
c            write(101,*)(rk(j,i),i=1,nqdof)
c         enddo
c         close(101)
c      endif
               

                  

         
         
c     unit 9 contains the time-dependent populations (from integrate.for)
         write(6,*)"Writing the time-dep pop to exci.tmp (unformatted)"
         open(unit=9, file='exci.tmp',
     >        status='replace',access='sequential',form='unformatted')
         
         
         ntime=0
         
         
         write(6,*)" Initializing initial e/h pair"
         nsteps=int(timall/timstep)
         nshots=max(int(timshot/timstep),10)
         write(6,*)"timshot = ",timshot
         write(6,*)"timstep = ",timstep
         write(6,*)"nsteps = ",nsteps
         write(6,*)"nshots = ",nshots

         call inittime(0.0,timall,nsteps,nshots)
         
c     Initialize the populations qc with an injected e-h pair 
         
         if( chargesep) then 
            write(6,*)" Generating initial charge separated state "
            call EHInjection(qc)
            exciton = .false.

         endif
         

         if( exciton)then
            write(6,'(" Initial state is the ",i2,"th CI eigenstate")')iex
            do i = 1,nqdof
               qc(i) = 0.0d0
            enddo
            qc(iex)  = 1.0d0
            write(6,*)" Writing Initial EH density-> exci.init"
            
            open(unit=8, file='exci.init',status='replace',
     *           access='sequential',form='formatted')                 
            
            call EHDensity(8,h(1,iex))
            close(8)
         endif
         

        
         write(6,*)"Integrating Redfield Matrix"
c     Integrate the Redfield equation

         call rkdumb(9,qc,nqvars,0.0d0,timall,nsteps)





c     qcav contains the time-average populations
c     qcav(iord(i)), i=1..nqvar, maps them in decreasing order 
         
         write(6,*)"Performing Time-averages"
         call TimeAverage(9,qcav)
         call Sort(qcav,iord,nqvars,-1)

         write(6,*)"Generating output files"
c     Calculate and output (exci.csv) time-dependent population density
         if ( pdensity) then      
            write(6,*)"exci.csv-- time-dependent population density"
            call ScanDensity(9,nsteps,1,nqvars)  

            open(unit=8, file='exci.csv',
     >           status='replace',
     *           access='sequential',
     >           form='formatted')                 
            write(8,'("#     time dependent population density")')
            write(8,120)'#  Energy | time =',(xtime(i), i=1,ntime)
 120        format(a20,1000f8.3)   
            call WriteSpectra(8)
            close(8)
         endif
         
c     Calculate and output (exci.lum) time-dependent fluorescence
         if ( pemission) then
            write(6,*)"exci.lum-- time-dependent fluorescence data"
            call ScanLuminescence(9,nsteps)      
            open(unit=8, 
     >           file='exci.lum',
     >           status='replace',
     *           access='sequential',
     >           form='formatted')             
            write(8,'("#    exci.lum-- time-dep fluorescence data")')
            write(8,120)'#  Energy | time =',(xtime(i), i=1,ntime)
            call WriteSpectra(8)
            close(8)
         endif
c     Output (exci.dat) the time-dependent populations of the
c     essential excited states

         write(6,*)"exci.dat-- time-dependent populations data"
         open(unit=8, file='exci.dat',status='replace',             
     *        access='sequential',form='formatted')
         write(8,'("## exci.dat-- time-dependent populations data")')
c         write(6,*)(iord(i),i=1,nqvars)

         do i = 1,nqdof
              iord(i) = i
         enddo
         call WriteEvolution(9,8,min(nqdof,maxrelax),iord)
         close(8)
      endif


c     Calculate and output time-independent or time-integrated
c     spectral distributions


      open(unit=8, file='exci.txt',status='replace',
     *     access='sequential',form='formatted')
      write(8,*)"#    exci.txt-- time-dependent populations data (text)"

 110  format(1000a15)
 111  format("#",I5,A15,5x,A40)
      ntime = 0
      write(6,*)" Generating Spectral Distributions"
      write(6,*)" Output file = exci.txt "
      write(6,*)"             Column Headings" 
      write(6,'(A5,A15,A20)')"#"," Title ","Description"
      write(6,*)"---------------------------------------------------"


      write(8,'("#",A5,A15,A20)')"#"," Title ","Description"
      write(8,'("#--------------------------------------------------")')
      
      ntime=ntime+1
      
c     Density of states
      
      ctitle(ntime) ='      dDOS'
      write(6,111)ntime,ctitle(ntime),"Density of States"
      write(8,111)ntime,ctitle(ntime),"Density of States"

      call DensityOfStates(spectrum(1,ntime),eigs,iehpar,0,eleft
     *     ,eright,ElGama,nspec,nqexc)
      ntime=ntime+1
c     Initial population
      
      
      ctitle(ntime)='      iPOP'
      write(6,111)ntime,ctitle(ntime),"Initial Population" 
      write(8,111)ntime,ctitle(ntime),"Initial Population" 
      call EHInjection(qcaux)
      call   SpectralDensity(spectrum(1,ntime),eigs,qcaux
     *     ,eleft,eright,ElGama,nspec,nqexc,0)

      
      
c     Diabatic (vertical) absorption spectrum
      if (ispinstate.le.0) then
         ntime=ntime+1
         ctitle(ntime)='      dABS'
         write(6,111)ntime,ctitle(ntime),"Vertical Absorb. Spectrum" 
         write(8,111)ntime,ctitle(ntime),"Vertical Absorb. Spectrum" 
         
         call   SpectralDensity(spectrum(1,ntime),eigs,transdip
     *        ,eleft,eright,ElGama,nspec,nqdof,0)
      endif
      if ( evolve.and.ispinstate.le.0) then
         ntime=ntime+1
c     Time-integrated luminescence
         ctitle(ntime) ='    intLUM'
         write(6,111)ntime,ctitle(ntime),"Time Integrated Luminescence"
         write(8,111)ntime,ctitle(ntime),"Time Integrated Luminescence"
         call ElVibEmit(qcav,spectrum(1,ntime),eleft,eright,nspec,ElGama)
         ntime=ntime+1
         call ZeroVar(qcaux,NQMAX)
         qcaux(1)=1.0d0
c     Luminescence from S1
         ctitle(ntime) ='   relxLUM'
         write(6,111)ntime,ctitle(ntime),"Luminescence from S1"
         write(8,111)ntime,ctitle(ntime),"Luminescence from S1"
         call ElVibEmit(qcaux,spectrum(1,ntime),eleft,eright,nspec,ElGama)

      endif
c     Adiabatic (electron-vibrational) absorption
      if ( relax) then
         ntime=ntime+1
         ctitle(ntime) ='      aDOS'
         write(6,111)ntime,ctitle(ntime),"Adiabatic Density of States"
         write(8,111)ntime,ctitle(ntime),"Adiabatic Density of States"
         call DensityOfStates(spectrum(1,ntime),eadiab,iehpar,0
     *        ,eleft,eright,ElGama,nspec,nqdof)
         if (ispinstate.eq.0) then
            ntime=ntime+1
            ctitle(ntime) ='      aABS'
         write(6,111)ntime,ctitle(ntime),"Adiabatic Absorption Spectum"
         write(8,111)ntime,ctitle(ntime),"Adiabatic Absorption Spectum"
            call ElVibAbsorb(spectrum(1,ntime),
     *           eleft,eright,nspec,ElGama)
         endif


      endif
c------------ end of master equation evolution section ------------





      write(8,110) '#Energy', (ctitle(i),i=1,ntime)
      call WriteSpectra(8)

      close(8)

c     level is apparently a flag for 
c     how verbose of output is desired.  
c     set this to be high for lots of info. 
      level = 10

      write(6,*)"exci.out-- accumulated output (text)"
      open(unit=8, file='exci.out',status='replace')
      call Output(8,iolevel,min(nqvars,6))
      close(8)

c     Get End Time
c     pause
      close(5)
      close(7)       
      

      call stoppg(0)

      end 


      subroutine stoppg(i)
      include "genparam.fi"
      integer i
      call date_and_time(edat,etim,ezon,edt)
      if(i.eq.0)then 
         write(6,'("Program terminated normally ")')
         call writedtinfo(6)
         close(999,status="delete")
      endif
      if(i.eq.1) then
         write(6,*) "::::::::::  INTERRUPT SIGNAL CAUGHT :::::::"
         write(6,*)"::::::::::  PROGRAM TERMINATED ::::::::::::"
         open(unit = 10,file ="EXCITON_STOP")
         close(unit=10,status="delete")
         call writedtinfo(6)
         close(999,status="delete")
      endif
      if(i.eq.2) then
         write(999,*)":::::: DIRECTORY LOCKED :::::::::::"
         write(999,*)"  exciton.x is running in this directory"
         write(999,*)"  terminating program"
         call writedtinfo(999)
      endif

      stop
      end
c     
c     
      subroutine CalcRedfield
c     Subroutine that constructs the Redfield matrix and caluculate
c     interstate rates
      include 'sizes.fi' 
      include 'parameters.fi'
      include 'elstates.fi'      
      include 'couplings.fi'
      include 'coords.fi'    
      include 'vibmodes.fi'
      include 'rates.fi'
      integer k,l,nt,istatemax
c     real(kind(1D0)) ScalProd, DLorentz
      real(kind(1D0)) overf, overe, spino
      

      call ZeroVar(rk,NQMAX*NQMAX*4)
      

      
c     Calculate the rates of the one-phonon interstate transition
c     (relaxation.for)

      
      write(6,*)"CalcRF: OK to here"       
      call CalcRateConst
      write(6,*)"CalcRF: OK to here"      

      
c      if (ispinstate.lt.0) then
c         nt=nqdof+1
c         call CalcRateConst(rk(nt,nt),h(1,nt),eigs(nt))
c      endif
      

      l=nqexc+1
c     Include the radiative transition rates to the ground state,
c     if requested
      if ( radrates) then
         do k=1,nqdof
            rk(k,l)=sradrate(k)
         end do
      endif

c     Include the nonradiative transition rates to the ground state
      do k=1,nqexc
         rk(k,l)=rk(k,l)+ knonrad
      end do

c     Include the intersystem crossing rates: testing phase yet (tools.for)
c      if (ispinstate.lt.0) then
c         spino=1.0d-2
c         do k=1,nqdof
c            do l=nt,nqdof*2
c               overf=ScalProd(h(1,k),h(1,l),nqdof)
c               overe=DLorentz(eigs(k)-eigs(l), ElGama)*ElGama
c               rk(k,l)=spino*overf*overf*overe
c               rk(l,k)=rk(k,l)
c           end do
c        end do
c      endif
C     Here rk stores the transition rates. To transform it into the 
c     Redfield matrix, the resultant depopulation rates must enter 
c     the diagonal. (relaxation.for)
      call CalcDepopulation(rk,nqexc+1,NQMAX*2)
      return
      end     
c     
c     
      subroutine TimeAverage(ir,abuff)
c     Subroutine that calculates the time-average populations
      include 'sizes.fi'
      integer ir,i
      real(kind(1D0)) ave
      real(kind(1D0)) buff(NQMAX*2),abuff(NQMAX*2),tim
      call ZeroVar(abuff,nqvars)
      rewind(ir)
      ave = 0.0d0
      do while (.true.)
         read(ir,end=100)tim,(buff(i),i=1,nqvars)
         ave=ave+1
         call SAdd(buff, abuff, abuff, nqvars, 1.0d0)
      end do
 100  continue
      do i=1,nqvars
         abuff(i)=abuff(i)/ave
      end do
      return
      end           
c     
c     
      subroutine EHInjection(qc)
c     Subroutine that initializes the populations qc with an injected
c     e-h pair
      include 'sizes.fi'
      include 'elstates.fi'
      real(kind(1D0)) qc(*)
      integer i

c  The matrix h at this point contains the CI 
c  eigenstates. (i believe)
c  are relaxing from the  highest energy eigenstate????
      do i=1,nqdof
         qc(i)=h(nqdof,i)**2
      end do
      if (ispinstate.lt.0) then
         do i=nqdof+1,nqdof*2
            qc(i)=3.0*h(nqdof,i)**2
         end do
         do i=1,nqdof*2
            qc(i)=qc(i)/3.0d0
         end do
      endif
      return
      end
c     
c     
      subroutine WriteEvolution(ir,iw,nmax,ibuff)
c     Subroutine to write the time-dependent populations of the essential
c     excited states
      include 'sizes.fi'
      include 'parameters.fi'
      include 'elstates.fi'
      integer ir,iw,nmax
      real(kind(1D0)) buff(NQMAX*2)
      integer ibuff(*)  
      real(kind(1D0)) SMin, SSum,  ScalProd, ExpectValue
      external DipoleElem, ExpectValue
      real(kind(1D0)) tim, yield, energy, dipolemom
      integer i
      character SpinState
      rewind(ir)
      write(iw,120) (SpinState(ibuff(i)),eigs(ibuff(i)),i=1,nmax)
      do while (.true.)
         read(ir,end=100)  tim, (buff(i),i=1,nqvars)
         energy=ScalProd(buff,eigs,nqvars)
         yield=ScalProd(buff,sradrate,nqvars)
         dipolemom=ExpectValue(buff,tim,DipoleElem)
         write(iw,110) tim, energy, dipolemom, yield*1.0D3,
     *        (buff(ibuff(i)),i=1,nmax) 
      end do
 100  continue
 110  format(1000f10.4)       
 120  format('      Time',
     *     '    Energy',
     *     '    Dipole',   
     *     '     Yield',
     *     1000('   ',a2,f5.3))
      return
      end           
      
c--------
      subroutine header

      write(6,*)"-----------------------------------------------------------------------"
      write(6,*)"-                                                                     -"
      write(6,*)"-                                                                     -"
      write(6,*)"- EEEEEEEE   X      X      CCCC     II   TTTTTTTT    OOO     N      NN-"
      write(6,*)"- EE          X    X      CC        II      TT     OO   OO   NN     NN-"
      write(6,*)"- EE           X  X       CC        II      TT     OO   OO   NNN    NN-"
      write(6,*)"- EEEEE         XX        CC        II      TT     OO   OO   NN N   NN-"
      write(6,*)"- EE           X  X       CC        II      TT     OO   OO   NN  N  NN-"
      write(6,*)"- EE          X    X      CC        II      TT     OO   OO   NN   N NN-"
      write(6,*)"- EEEEEEE    X      X      CCCC     II      TT       OOO     NN     NN-"
      write(6,*)"-                                                                     -"
      write(6,*)"-                                                                     -"
      write(6,*)"-                             EXCITON'O4                              -"
      write(6,*)"-                     (-)------(+) -> -----(*)-----                   -"
      write(6,*)"-                                                                     -"
      write(6,*)"-                         Bittner Research Group                      -"
      write(6,*)"-                            Univ of HOUSTON                          -"
      write(6,*)"-                                                                     -"
      write(6,*)"-                             exciton.v1.0                            -"
      write(6,*)"-   References:                                                       -"
      write(6,*)"-  1.) Polaron-excitons and electron-vibrational band shapes          -"
      write(6,*)"-      in conjugated polymers, S. Karabunarliev and E. R. Bittner,    -"
      write(6,*)"-      J. Chem. Phys.  118 , 4291 (2003).                             -"
      write(6,*)"-  2.)Spin dependent electron-hole recombination kinetics in          -"
      write(6,*)"-     luminescent organic polymers,                                   -"
      write(6,*)"-     S. Karabunarliev and E. R. Bittner,                             -"
      write(6,*)"-     Phys. Rev. Lett.  90 , 057402 (2003).                           -"
      write(6,*)"-  3.) Energy relaxation dynamics and universal scaling laws in       -"
      write(6,*)"-      organic light emitting diodes,                                 -"
      write(6,*)"-      Eric R. Bittner and Stoyan Karabunarliev,                      -"
      write(6,*)"-      Int. J. Quant. Chem. 95 , 521 (2003).                          -"
      write(6,*)"-                                                                     -"
      write(6,*)"- Original Source code is Copyright (c) 2004 by Eric R. Bittner       -"
      write(6,*)"-                                                                     -"
      write(6,*)"-  exciton.x, a program for computing electron/hole dynamics in       -"
      write(6,*)"-             conjugated polymers.                                    -"
      write(6,*)"-                                                                     -"
      write(6,*)"-    This program is free software; you can redistribute it and/or    -"
      write(6,*)"-    modify it under the terms of the GNU General Public License      -"
      write(6,*)"-    as published by the Free Software Foundation; either version 2   -"
      write(6,*)"-    of the License, or (at your option) any later version.           -"
      write(6,*)"-    This program is distributed in the hope that it will be useful,  -"
      write(6,*)"-    but WITHOUT ANY WARRANTY; without even the implied warranty of   -"
      write(6,*)"-    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    -"
      write(6,*)"-    GNU General Public License for more details.                     -"
      write(6,*)"-                                                                     -"
      write(6,*)"-    You should have received a copy of the GNU General Public License-"
      write(6,*)"-    along with this program; if not, write to the Free Software      -"
      write(6,*)"-    Foundation, Inc., 59 Temple Place                                -"
      write(6,*)"     Suite 330, Boston, MA  02111-1307, USA.                          -"
      write(6,*)"-                                                                     -" 
      write(6,*)"-  Contact information                                                -"
      write(6,*)"-   Prof. Eric R. Bittner                                             -"
      write(6,*)"-  Department of Chemistry                                            -"
      write(6,*)"-  University of Houston                                              -"
      write(6,*)"-  Houston, TX  77204                                                 -"
      write(6,*)"- email:  bittner@uh.edu                                              -"
      write(6,*)"-   url:  http://k2.chem.uh.edu                                       -"
      write(6,*)"-----------------------------------------------------------------------"
      RETURN
      end
