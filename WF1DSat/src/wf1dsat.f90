	PROGRAM sat_flow1d
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
			sat%mesh%vmod_qent =  sat%calc%boundary%get_qbound_file(sat%time%t) !CHECK if qtop has to be multiplied by 2.
		end if

		! ---- 06.04 ---- Update Dirichlet boundary conditions from constraints to nodes
		call sat%set_qent_from_constraint_to_nodes()
		!call sat%calc%update_load(sat%mesh)

		!!Update loads in each node (flow entering) for actual timestep
		!if (sat%time%t<sat%calc%parameters%tmax/2) then !Update here values of qent in nvmod
		!	!sat%mesh%nvmod_qent = sat%mesh%nvmod_qent !In this line the flow enter data is updated for this timestep.
		!else
		!	sat%mesh%vmod_qent = 0.0
		!end if
		!!Waterflow entering is updated with f_q_hor_infiltr(nodes%x,mesh) it does an interpolation from the nvmod to each node given x values of nvmod and x values of nodes) (Can be optimized)
		!call sat%calc%update_load(sat%mesh)

		!Begin iterations stepping (until converged or dt below dtmin with max iterations)

		! ---- 06.05 ---- Iteration counter to 0
		iterconverg = 0
		! ---- 06.06 ---- Initialize isconverged to false
		isconverged = .false.

		! ---- 06.07 ---- CONVERGENCE: do while not(isconverged or dt<dtmin and max iteration reached))
		CONVERGENCE:do while(.not.(isconverged .or. (sat%time%dt<=sat%calc%parameters%dtmin .and. iterconverg==(sat%calc%parameters%it_max-1))))
			! ---- 06.07.01 ---- Increase iteration counter
			iterconverg = iterconverg+1
			! ---- 06.07.02 ---- Update wc_temp from h_temp (CHECK)

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
		!For printing purposes (CHECK) (This is not Ok, reverse the Solution to get a result.
		!WRITE(*,'("Iter: ",i6," t: ", f10.3," dt: ",E10.3," h(0): ",f10.3," h(l): ",f10.3," Qent: ",E10.3," QSal: ",E10.3)') iterconverg, sat%calc%t, sat%calc%dt, sat%calc%nodes%hnew(1), sat%calc%nodes%hnew(sat%calc%nodes%count), sat%mesh%vmod_qent(1), 0.0
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
	

999	END PROGRAM