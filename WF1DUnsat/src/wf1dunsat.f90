	!------------------------------------------------------------------------------
	!        FLOW1DUNSAT
	!------------------------------------------------------------------------------
	! TITLE         : FLOW1DUNSAT: One dimensional flow in variably saturated media
	! PROJECT       : An 1.5D simulation for a variably saturated multilayer system
	! MODULE        : FLOW1DUNSAT
	! URL           : ...
	! AFFILIATION   : The University of Nottingham
	! DATE          : 23th of October of 2019
	! REVISION      : 0.0
	!> @author
	!> Iván Campos-Guereta Díez
	!
	! DESCRIPTION:
	!> This is the main program to simulate one dimensional variably saturated
	!> materials.
	!------------------------------------------------------------------------------

	program unsat_flow1d
	use unsat_mod_ty_model,only:ty_unsat_model
	use unsat_mod_ty_elements,only:ty_unsat_elements
	use com_mod_ty_elements,only:ty_com_elements

	implicit none
	include 'inc_precision.fi'

	!---- 00 ---- Create instances (unsat, elemcom, elem), of classes (ty_unsat_model, ty_com_elements, ty_sat_elements)
	logical,parameter						::IS_NOT_TIME_DEPENDANT=.false.,IS_TIME_DEPENDANT=.true.
	integer,parameter						::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDER_HOLD=2, CONSIDER_ALL=3
	!character*400,parameter			::FILEINPUT						='INPUTS\unsat_inputs.txt'
	!character*400,parameter			::FILEBOUNDARY				='INPUTS\flow1d_bound_h.txt'
	!character*400,parameter			::FILEOUTPUT_NODES		='OUTPUTS\unsat_outputs_nodes.csv'
	!character*400,parameter			::FILEOUTPUT_ELEMENTS	='OUTPUTS\unsat_outputs_elements.csv'
	character*400::fileinput
	character*400::namefile !This will be the name without extension
	character*400::fileboundary						
	character*400::FILEOUTPUT_NODES			
	character*400::FILEOUTPUT_ELEMENTS	

	type(ty_unsat_model),target				::unsat
	class(ty_com_elements),pointer		::elemcom
	class(ty_unsat_elements),pointer	::elem

	real(kind=dpd)	::htop,qtop
	logical					::isconverged=.false.,isconvergedall=.true.
	integer					::iterconverg, npos, nargs

	!---- 01 ---- Read file inputs
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
	FILEOUTPUT_NODES		=trim(namefile)//'.outnodu'
	FILEOUTPUT_ELEMENTS	=trim(namefile)//'.outelmu'	
	
	call unsat%readfileinput(FILEINPUT,fileboundary)

	!---- 02 ---- Allocate and construct all instances
	call unsat%construct()

	!---- 03 ---- Open output files
	open (31,file=FILEOUTPUT_NODES,		status='unknown')
	open (32,file=FILEOUTPUT_ELEMENTS,status='unknown')

	!Construct a pointer to activate extension of element derived type>
	elemcom => unsat%calc%elements
	select type (elemcom)
	type is (ty_unsat_elements)
		elem=>elemcom
	end select

	!---- 04 ---- Set initial conditions
	call unsat%calc%set_to_initial() !Read initial conditions

	!---- 05 ---- Update all interface values for initial conditions
	call unsat%constraints%update_all(unsat%calc,unsat%calc%time%dt)
	!call unsat%constraints%update_all(unsat%calc,unsat%calc%time%dt,1.0_dpd)

	!---- 06 ---- Update wc_old, wc_temp, wc_new (water content)
	call unsat%calc%update_th_from_h(CONSIDER_ALL)
	call unsat%print_timestep(31,32,1)

	!---- 07 ---- assign Dirichlet conditions for phreatic on the bottom
	call unsat%calc%set_dirichlet_to_node(1,unsat%calc%nodes%hnew(1))

	!---- 08 ---- build matrixes of linear system that don’t depend on time
	call unsat%calc%build_linearsystem(IS_NOT_TIME_DEPENDANT,CONSIDER_HTEMP)

	!---- 09 ---- TIMESTEPPING: do while time < time to model
	TIMESTEPPING: do while(unsat%time%t < unsat%calc%parameters%Tmax)

		!---- 09.01 ---- increase t with dt
		call unsat%time%increase_time() !Increase timestep t=t+dt (and if tprint is between t and t+dt then t=tprint and checkprint=true)

		!---- 09.02 ----	Estimate initial hnew for current timestep
		call unsat%calc%estimate_hnew_for_new_timestep() !Estimate hnew to begin iterations (hnew=hold+slopeold·dt with slopeold = (hold_n-1-hold_n)/dtold, in the first step slopeold is 0.0)

		! ---- 09.03 ---- Update Dirichlet boundary conditions from constraints to nodes
		if (unsat%calc%layers%topboundbyh) then
			htop =  unsat%calc%boundary%get_hbound_file(unsat%time%t)
			call unsat%calc%set_dirichlet_to_node(unsat%calc%nodes%count,htop)
		end if

		!Include here boundary at the bottom if it changes

		! ---- 09.04 ---- Update Newman boundary conditions at top and apply
		if (unsat%calc%layers%topboundbyq) then
			qtop =  unsat%calc%boundary%get_qbound_file(unsat%time%t) !Check if qtop has to be multiplied by 2.
			unsat%calc%colbound(unsat%calc%nodes%count)= qtop !Put qtop on colbound vector.
		end if

		! ---- 09.05 ---- Iteration counter to 0
		iterconverg = 0

		! ---- 09.06 ---- Initialize isconverged to false
		isconverged = .false.

		! ---- 09.07 ---- CONVERGENCE: do while not(isconverged or dt<dtmin and max iteration reached))
		CONVERGENCE:do while(.not.(isconverged .or. (unsat%time%dt<=unsat%calc%parameters%dtmin .and. iterconverg==(unsat%calc%parameters%it_max-1))))
			! ---- 09.07.01 ---- Increase iteration counter
			iterconverg = iterconverg+1

			! ---- 09.07.02 ---- Update wc_temp from h_temp (CHECK)
			call unsat%calc%update_th_from_h(CONSIDER_HNEW)

			! ---- 09.07.03 ---- Call subroutine iterate
			! ---- 09.07.03.01 ----	Update h_temp with previous h_new
			! ---- 09.07.03.02 ----	Build linear system
			! ---- 09.07.03.03 ----	Apply Dirichlet conditions to linear system
			! ---- 09.07.03.04 ----	Solve linear system
			! ---- 09.07.03.05 ---- Update h_new from solution
			! ---- 09.07.03.06 ----	Calculate error epsh
			! ---- 09.07.03.07 ----	Isconverged = epsh<epsh_tolerance
			! ---- 09.07.03.08 ----	End subroutine iterate
			!call unsat%calc%iterate(isconverged,isconvergedall)
			call unsat%calc%iterate(isconverged)

			! ---- 09.07.04 ----	Update th_new values from h_new values
			call unsat%calc%update_th_from_h(CONSIDER_HNEW)


			if (iterconverg<unsat%calc%parameters%it_min) isconverged=.false.
			! ---- 09.07.05 ----	If (not isconverged and max iteration reached)
			if((.not.isconverged.and. iterconverg==unsat%calc%parameters%it_max).or.(.not.isconvergedall)) then
				unsat%time%checkprint = .false.
				! ---- 09.07.05.01 ---- Reduce timestep
				call unsat%calc%time%factor_timestep(1/3.0_dpd) !Update timestep with boundaries on tmin and tmax (dt=factor·dtold)
				! ---- 09.07.05.02 ---- Revert to old
				call unsat%calc%revert_to_old()
				! ---- 09.07.05.03 ---- Exit convergence iterations
				iterconverg = 0
				isconvergedall=.true.
			end if

		end do CONVERGENCE

		!Iterantion converged or min dtmin reached
		!For printing purposes (Check) (This is not Ok, reverse the Solution to get a result.
		!WRITE(*,'("Iter: ",i6," t: ", f10.3," dt: ",E10.3," h(0): ",f10.3," h(l): ",f10.3," Qent: ",E10.3," QSal: ",E10.3)') iterconverg, sat%calc%t, sat%calc%dt, sat%calc%nodes%hnew(1), sat%calc%nodes%hnew(sat%calc%nodes%count), sat%mesh%vmod_qent(1), 0.0

		! ---- 09.08 ---- Get results of timestep
		call unsat%calc%get_results_elements()
		!call unsat%calc%get_results_nodes()

		! ---- 09.09 ---- Write outputs	 if print times reached
		!WRITE(*,'("Iter: ",i6," t: ", f10.3," dt: ",E10.3," h(0): ",E10.3," qent(0): ",E10.3," qincvol(l): ",E10.3," qhor: ",E10.3)') iterconverg, unsat%calc%t, unsat%calc%dt, unsat%calc%nodes%hnew(1), sum(elem%results_qent), sum(elem%results_incvol), sum(elem%results_incqhor)
		WRITE(*,'("Iter: ",i6," t: ", f10.3," dt: ",E10.3," h(0): ",E10.3," h(end): ",E10.3," h(end/2): ",E10.3," qhor: ",E10.3)') iterconverg, unsat%time%t, unsat%time%dt, unsat%calc%nodes%hnew(1), unsat%calc%nodes%hnew(unsat%calc%nodes%count), unsat%calc%nodes%hnew(int(unsat%calc%nodes%count/2))
		IF (unsat%time%checkprint) call unsat%print_timestep(31,32,1) !Check if we are in printstep and print

		! ---- 09.10 ---- Update timestep dt
		call unsat%time%update_dt(iterconverg) !update dt, decreasing if>it_dec_dt or increasing if < it_inc_dt

		! ---- 09.11 ---- Set old values from new values
		call unsat%calc%set_old() !All new values passes now to old, and slopeold=(hnew-hold)/(t-told), told=t, dtold=dt, hold=hnew
	end do TIMESTEPPING

	! ---- 10 ---- Close output files	and deallocate all
	close(31)
	close(32)


999	end program unsat_flow1d