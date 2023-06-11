	!********************************************************************************************************************
	! WF1DSAT
	!********************************************************************************************************************
	! TITLE         : WF1DSAT Main program
	! PROJECT       : WF1DSAT
	! MODULE        : -
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	!> <B>Main program to calculate the watertable on a multilayer system due to an infiltration defined in time</B>
	!>
	!> The program with one parameter equal to fileinput path: WF1DSAT.exe  xxxx
	!>
	!> In the same path than the input:
	!> xxxx.wfsinp: Input textfile containing information of parameters (after heading BLOCK A), materials (after heading BLOCK B),
	!> layers (after heading BLOCK C), mesh (after heading BLOCK D)
	!>
	!> xxxx.wfbound: Input textfile containing the information about boundary values in time (table: t(T),h(L), or t(T),q(L3·T-1)
	!>
	!> As a result, five output *.csv files are generated:
	!>
	!> xxxx.outnods.csv: Values at the nodes at each defined print time
	!>
	!> xxxx.outelms.csv: Values at the elements at each defined print time
	!
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	PROGRAM WF1DSAT
	use sat_mod_ty_model,only:ty_sat_model
	use sat_mod_ty_elements,only:ty_sat_elements
	use com_mod_ty_elements,only:ty_com_elements

	IMPLICIT NONE
	INCLUDE 'inc_precision.fi'

	!---- 00 ---- Create instances (sat, elemcom, elem), of classes (ty_sat_model, ty_com_elements, ty_sat_elements)
	type(ty_sat_model),target::sat
	!CHARACTER*400,parameter::FILEINPUT='INPUTS\sat_inputs.txt'
	!character*400,parameter::fileboundary				='INPUTS\flow1d_bound_q.txt'
	!CHARACTER*400,parameter::FILEOUTPUT_NODES		='OUTPUTS\sat_outputs_nodes.csv'
	!CHARACTER*400,parameter::FILEOUTPUT_ELEMENTS='OUTPUTS\sat_outputs_elements.csv'
	character*400::fileinput
	character*400::namefile !This will be the name without extension
	character*400::fileboundary
	character*400::FILEOUTPUT_NODES
	character*400::FILEOUTPUT_ELEMENTS
	character*400::FILEOUTPUT_STATS

	logical,parameter::IS_NOT_TIME_DEPENDANT=.false.,IS_TIME_DEPENDANT=.true.
	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDERHOLD=2
	class(ty_com_elements),pointer::elemcom
	class(ty_sat_elements),pointer::elem
	real(kind=dpd),pointer::t

	integer::iterconverg, npos, nargs
	logical::isconverged,isconvergedall=.true.

	!For statistics:
	real(kind=8)	  ::cputimeinit=0.0_8,cputimeend=0.0_8,systimeinit=0.0_8,systimeend=0.0_8,totaldt=0.0_8,meandt=0.0_8,maxdt=0.0_8,mindt=1.0E10_8
	integer(8)			::itertotal=0,itersteps=0,timemilisec=0,countrate=1000

	!---- 01 ---- Read file inputs
	call CPU_TIME(cputimeinit)
	call SYSTEM_CLOCK(timemilisec,countrate)
	systimeinit=real(timemilisec,8)/real(countrate,8)

	nargs = COMMAND_ARGUMENT_COUNT()
	if (nargs==0) then
		WRITE(*,*) 'Pease include input file as the first argument argument'
		goto 999 !If there is no argument then end programm
	end if

	call GET_COMMAND_ARGUMENT(1,fileinput) !Input file is the first argument of the program
	write(*,*) 'number of arguments',nargs
	npos = scan(trim(fileinput),".", BACK= .true.) !To remove extension
	if (npos>0) namefile = trim(fileinput(1:npos-1)) !Now without extension

	fileboundary				=trim(namefile)//'.wfbound'
	FILEOUTPUT_NODES		=trim(namefile)//'.outnods'
	FILEOUTPUT_ELEMENTS	=trim(namefile)//'.outelms'
	FILEOUTPUT_STATS	  =trim(namefile)//'.outstat'

	call sat%readfileinput(FILEINPUT,fileboundary)

	!---- 02 ---- Allocate and construct all instances
	call sat%construct()

	!---- 03 ---- Open output files
	OPEN (31,FILE=FILEOUTPUT_NODES,     STATUS='unknown')
	OPEN (32,FILE=FILEOUTPUT_ELEMENTS,     STATUS='unknown')

	!---- 04 ---- Set initial conditions
	call sat%calc%set_to_initial() !Set time to initial and hnew, htemp and hold to hinit
	call sat%print_timestep(31,32,1)

	!---- 05 ---- build matrixes of linear system that don’t depend on time
	call sat%calc%build_linearsystem(IS_NOT_TIME_DEPENDANT,CONSIDER_HTEMP) !Construct all matrices that do not depend on time

	elemcom=> sat%calc%elements
	select type (elemcom)
	type is (ty_sat_elements)
		elem => elemcom
	end select

	t=> sat%time%t

	!---- 06 ---- TIMESTEPPING: do while time < time to model
	TIMESTEPPING: do while(sat%time%t < sat%calc%parameters%Tmax)

		!---- 06.01 ---- increase t with dt
		call sat%time%increase_time() !Also if tprint is between t and t+dt then t=tprint and checkprint=true

		!---- 06.02 ----	Estimate initial hnew for current timestep
		call sat%calc%estimate_hnew_for_new_timestep() !Estimate hnew to begin iterations (hnew=hold+slopeold·dt with slopeold = (hold_n-1-hold_n)/dtold, in the first step slopeold is 0.0)

		! ---- 06.03 ---- Update Newman boundary conditions at top and apply
		if (sat%calc%layers%topboundbyq) then
			sat%mesh%vmod_qent =  sat%calc%boundary%get_qbound_file(sat%time%t)
		end if

		! ---- 06.04 ---- Update Dirichlet boundary conditions from constraints to nodes
		call sat%set_qent_from_constraint_to_nodes()

		!Begin iterations stepping (until converged or dt below dtmin with max iterations)

		! ---- 06.05 ---- Iteration counter to 0
		iterconverg = 0
		! ---- 06.06 ---- Initialize isconverged to false
		isconverged = .false.

		! ---- 06.07 ---- CONVERGENCE: do while not(isconverged or dt<dtmin and max iteration reached))
		CONVERGENCE:do while(.not.(isconverged .or. (sat%time%dt<=sat%calc%parameters%dtmin .and. iterconverg==(sat%calc%parameters%it_max-1))))
			! ---- 06.07.01 ---- Increase iteration counter
			iterconverg = iterconverg+1
			! ---- 06.07.02 ---- Update wc_temp from h_temp

			! ---- 06.07.02 ---- Call subroutine iterate
			! ---- 06.07.02.01 ----	Update h_temp with previous h_new
			! ---- 06.07.02.02 ----	Build linear system
			! ---- 06.07.02.03 ----	Apply Dirichlet conditions to linear system
			! ---- 06.07.02.04 ----	Solve linear system
			! ---- 06.07.02.05 ---- Update h_new from solution
			! ---- 06.07.02.06 ----	Calculate error epsh
			! ---- 06.07.02.07 ----	Isconverged = epsh<epsh_tolerance
			! ---- 06.07.02.08 ----	End subroutine iterate
			!call sat%calc%iterate(isconverged,isconvergedall)
			itertotal=itertotal+1
			call sat%calc%iterate(isconverged)

			if (iterconverg<sat%calc%parameters%it_min) isconverged=.false.
			! ---- 06.07.03 ----	If (not isconverged and max iteration reached)
			if((.not.isconverged .and. iterconverg==sat%calc%parameters%it_max).or.(.not.isconvergedall)) then
				sat%time%checkprint = .false.
				! ---- 06.07.03.01 ---- Reduce timestep
				call sat%calc%time%factor_timestep(1/3.0_dpd) !Update timestep with boundaries on tmin and tmax (dt=factor·dtold)
				! ---- 06.07.03.02 ---- Revert to old
				call sat%calc%revert_to_old()
				! ---- 06.07.03.03 ---- Exit convergence iterations
				iterconverg = 0
				isconvergedall=.true.
			end if

		end do CONVERGENCE

		!Iteration converged or min dtmin reached

		! ---- 06.08 ---- Get results of timestep
		call sat%calc%get_results_elements()
		call sat%calc%get_results_nodes()
		! ---- 06.09 ---- Write outputs	 if print times reached
		WRITE(*,'("Iter: ",i6," t: ", f10.3," dt: ",E10.3," h(0): ",E10.3," qent(0): ",E10.3," incvoldt(l): ",E10.3," dqhordx_all: ",E10.3," dqhordx_from_inc_all: ",E10.3)') iterconverg, sat%time%t, sat%time%dt, sat%calc%nodes%hnew(1), sum(elem%results_qent), sum(elem%results_incvoldt), sum(elem%results_dqhordx_all), sum(elem%results_dqhordx_from_incvoldt_all)
		IF (sat%time%checkprint) call sat%print_timestep(31,32,1) !Check if we are in printstep and print

		! ---- 06.10 ---- Update timestep dt
		call sat%time%update_dt(iterconverg) !update dt, decreasing if>it_dec_dt or increasing if < it_inc_dt

		! ---- 06.11 ---- set old values from new values for next timestep
		call sat%calc%set_old() !All new values passes now to old, and Slopeold=(hnew-hold)/(t-told), told=t, dtold=dt, hold=hnew

		!Log statistics
		itersteps=itersteps+1
		totaldt=totaldt+sat%time%dt
		mindt=min(mindt,sat%time%dt)
		maxdt=max(maxdt,sat%time%dt)
	end do TIMESTEPPING

	! ---- 07 ---- Close output files	and deallocate all
	close(31)
	close(32)

	!PRINT STATS:
	!Log statistics
	call CPU_TIME(cputimeend)
	call SYSTEM_CLOCK(timemilisec,countrate)
	systimeend=real(timemilisec,8)/real(countrate,8)
	meandt=totaldt/real(itersteps,8)

	open (33,file=FILEOUTPUT_STATS,status='unknown')
	write(33,'("itertotal,nsteps,cputime_s,systime_s,mintimestep,meantimestep,maxtimestep")')
	write(33,'(i6,",",i6,",",F10.3,",",F10.3,",",E10.3,",",E10.3,",",E10.3)') itertotal,itersteps,cputimeend-cputimeinit,systimeend-systimeinit,mindt,meandt,maxdt
	close(33)


999 END PROGRAM WF1DSAT