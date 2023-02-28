	!********************************************************************************************************************
	! WF15DSATUNSAT
	!********************************************************************************************************************
	! TITLE         : WF15DSATUNSAT Main program
	! PROJECT       : WF15DSATUNSAT
	! MODULE        : -
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	!> <B>Main program to calculate the water flow on a multilayer section by using an asembly of 1D submodels</B>
	!>
	!> The program with one parameter equal to fileinput path: WF15DSATUNSAT.exe  xxxx
	!>
	!> In the same path than the input:
	!> xxxx.wfuinp: Input textfile containing information of parameters (after heading BLOCK A), materials (after heading BLOCK B),
	!> layers (after heading BLOCK C), mesh (after heading BLOCK D)
	!>
	!> xxxx.wfbound: Input textfile containing the information about boundary values in time (table: t(T),h(L), or t(T),q(L3·T-1)
	!>
	!> As a result, five output *.csv files are generated:
	!>
	!> xxxx.outcons.csv: Output at the constraints at all timesteps
	!>
	!> xxxx.outstat.csv: Statistics about the run
	!>
	!> xxxx.outnods.csv: Values at the nodes on WF1DSAT submodels at each defined print time
	!>
	!> xxxx.outelms.csv: Values at the elements on WF1DSAT submodels at each defined print time
	!>
	!> xxxx.outnodu.csv: Values at the nodes on WF1DUNSAT submodels at each defined print time
	!>
	!> xxxx.outelmu.csv: Values at the elements on WF1DUNSAT submodels at each defined print time
	!
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************


	program WF15DSATUNSAT
	use model_mod_ty_model,only:ty_model

	implicit none
	include 'inc_precision.fi'

	!---- 00 ---- Create instances:model of class:ty_unsat_model, and define some parameters.
	logical,parameter::IS_NOT_TIME_DEPENDANT=.false.,IS_TIME_DEPENDANT=.true.
	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDER_HOLD=2, CONSIDER_ALL=3
	integer,parameter::MIN_ITERATIONS_IN_MODEL=1,MAX_ITERATIONS_IN_MODEL=20
	real(kind=dpd),parameter::MAX_ERROR_HSATHUNSAT=1E-3,MAX_ERROR_QSATQUNSAT=1E-6

	type(ty_model),target::model
	character*400::fileinput
	character*400::namefile !This will be the name without extension

	character*400::fileboundary
	character*400::fileoutput_nodes_v
	character*400::fileoutput_elements_v
	character*400::fileoutput_nodes_h
	character*400::fileoutput_elements_h
	character*400::fileoutput_constraints
	character*400::FILEOUTPUT_STATS
	character*100::debugmsg_txt=''
	character*1::printchar = ' '

	real(kind=dpd)::htop,qtop,epshsat,inchsat,maxhsat=0.0_dpd,maxhsatold=0.0_dpd,maxhsattemp=0.0_dpd,epshsathunsat,epsqsatunsat
	logical::isconverged=.false.
	integer::iterconverg,itermodel,iterunsat,itersat,iternrel,iterconvergmax
	logical::isfirsthsatloop=.false.
	logical,allocatable::ishsatover0(:),isfirsthsat(:)
	integer::nunsat,nsat,iu,is, npos, nargs
	real(kind=dpd),pointer::t
	logical::isconvergedall,ismaxitreached

	!For statistics:
	real(kind=8)	  ::cputimeinit=0.0_8,cputimeend=0.0_8,systimeinit=0.0_8,systimeend=0.0_8,totaldt=0.0_8,meandt=0.0_8,maxdt=0.0_8,mindt=1.0E10_8
	integer(8)			::itertotal=0,itersteps=0,timemilisec=0,countrate=1000


	!---- 01 ---- read input files
	call CPU_TIME(cputimeinit)
	call SYSTEM_CLOCK(timemilisec,countrate)
	systimeinit=real(timemilisec,8)/real(countrate,8)

	nargs = COMMAND_ARGUMENT_COUNT()
	if (nargs==0) then
		write(*,*) 'C Pease include input filename as the first argument'
		goto 999 !If there is no argument then end programm
	end if

	call GET_COMMAND_ARGUMENT(1,fileinput) !Input file is the first argument of the program
	!write(*,*) 'number of arguments',nargs
	npos = scan(trim(fileinput),".", BACK= .true.) !To remove extension
	if (npos>0) namefile = trim(fileinput(1:npos-1)) !Now without extension

	fileboundary					=trim(namefile)//'.wfbound'
	fileoutput_nodes_v		=trim(namefile)//'.outnodu.csv'
	fileoutput_elements_v	=trim(namefile)//'.outelmu.csv'
	fileoutput_nodes_h		=trim(namefile)//'.outnods.csv'
	fileoutput_elements_h	=trim(namefile)//'.outelms.csv'
	fileoutput_constraints	=trim(namefile)//'.outcons.csv'
	fileoutput_stats			=trim(namefile)//'.outstat.csv'

	call model%readfileinput(fileinput,fileboundary)

	!---- 02 ---- construct all instances
	call model%construct()
	!---- 03 ---- open output files
	open (31,file=fileoutput_nodes_v,     status='unknown')
	open (32,file=fileoutput_elements_v,  status='unknown')
	open (33,file=fileoutput_nodes_h,     status='unknown')
	open (34,file=fileoutput_elements_h,  status='unknown')
	open (51,file=fileoutput_constraints,  status='unknown')


	!Assign some variables and pointers
	nunsat=model%mesh%vmod_count
	if(model%layers%bottombyphl) then
		nsat = model%layers%count-1
	else
		nsat = model%layers%count
	end if
	t=>model%time%t
	iterunsat = 0
	itersat = 0
	itermodel = 0

	!---- 04 ---- Set initial conditions to WF1DUNSAT
	do iu=1,nunsat
		!Set time to initial and hnew, htemp and hold to hinit
		call model%unsat(iu)%calc%set_to_initial()
		call model%unsat(iu)%calc%update_th_from_h(CONSIDER_ALL)

		!Set Dirichlet into first node
		if (model%layers%bottombyphl) call model%unsat(iu)%calc%set_dirichlet_to_node(1,model%unsat(iu)%calc%nodes%hnew(1))
	end do

	!---- 05 ---- Set initial conditions to WF1DSAT
	do is=1,nsat
		!Set time to initial and hnew, htemp and hold to hinit
		call model%sat(is)%calc%set_to_initial()
		call model%sat(is)%calc%get_results_elements()
		call model%sat(is)%calc%get_results_nodes()
	end do

	!Update results...
	do iu=1,nunsat
		call model%unsat(iu)%calc%get_results_elements()
		call model%unsat(iu)%calc%get_results_nodes()
		call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt)
	end do

	!---- 06 ---- build WF1DUNSAT linear systems arrays that don’t change on time
	do iu=1,nunsat
		call model%unsat(iu)%calc%build_linearsystem(IS_NOT_TIME_DEPENDANT,CONSIDER_HTEMP)
	end do

	!---- 07 ---- build WF1DSAT linear systems arrays that don’t change in time on time
	do is=1,nsat
		call model%sat(is)%calc%build_linearsystem(IS_NOT_TIME_DEPENDANT,CONSIDER_HTEMP)
	end do

	!Initialize some variables...
	if(.not.allocated(ishsatover0)) allocate(ishsatover0(nsat))
	if(.not.allocated(isfirsthsat)) allocate(isfirsthsat(nsat))
	isfirsthsat = .false.
	ishsatover0 = .false.
	model%time%iter_total = 0

	!---- 08 ---- print first timestep to files
	do iu=1,nunsat
		call model%unsat(iu)%print_timestep(31,32,iu)
	end do
	do is=1,nsat
		call model%sat(is)%print_timestep(33,34,is)
	end do
	call model%print_alltimes(51)

	!---- 09 ---- TIMESTEPPING: do while t < tmax
	TIMESTEPPING: do while(model%time%t < model%parameters%Tmax)
		!---- 09.01 ---- increase time
		call model%time%increase_time()

		!---- 09.02 ---- read boundary properties on top of 1DUNSAT elements (but don’t apply yet)
		if (model%layers%topboundbyh) htop =  model%boundary%get_hbound_file(model%time%t)
		if (model%layers%topboundbyq) qtop =  model%boundary%get_qbound_file(model%time%t)

		!---- 09.03 ---- run WF1DUNSAT: with imposed source/sink terms qs (or Neumann boundary conditions) from WF1DSAT (old)

		!---- 09.03.01 ---- estimate hnew for new timestep
		do iu=1,nunsat
			call model%unsat(iu)%calc%estimate_hnew_for_new_timestep()
		end do

		!---- 09.03.02 ---- impose source/sink terms qs in WF1DUNSAT from WF1DSAT (as Neumann boundary condition)
		call model%constraints%set_newmann_from_sat_h(1.0_dpd)
		do iu=1,nunsat
			!boundary properties on top
			if (model%layers%topboundbyh) call model%unsat(iu)%calc%set_dirichlet_to_node(model%unsat(iu)%calc%nodes%count,htop)
			if (model%layers%topboundbyq) call model%unsat(iu)%calc%set_newmann_to_node(model%unsat(iu)%calc%nodes%count,qtop)
		end do

		!---- 09.03.03 ---- impose boundary on top and bottom of WF1DUNSAT
		do iu=1,nunsat
			call model%unsat(iu)%calc%assign_newmann()
		end do

		!---- 09.03.04 ---- CONVERGENCE_WF1DUNSAT: iterate WF1DUNSAT while not converged or min dt or max iterations reached --------------
		!PARALLEL COSIMULATION (CAN BE DONE IN PARALLEL) ----------------------------------------------------------------------------------------
		isconvergedall = .true.
		ismaxitreached = .false.
		iterconvergmax = 0
		iterunsat = 0
		WFUNSAT1: do iu=1,nunsat
			isconverged = .false.
			iterconverg = 0
			do while(.not.(isconverged .or. (model%time%dt<=model%parameters%dtmin .and. iterconverg==(model%parameters%it_max-1))))
				model%time%iter_total = model%time%iter_total + 1
				iterconverg = iterconverg+1
				call model%unsat(iu)%calc%update_th_from_h(CONSIDER_HNEW)
				itertotal=itertotal+1
				call model%unsat(iu)%calc%iterate(isconverged)
				if (iterconverg<model%parameters%it_min) isconverged=.false.
				if(.not.isconverged .and. iterconverg==model%parameters%it_max) exit
			end do
			iterunsat = max(iterconverg,iterunsat)
			if (.not.isconverged) isconvergedall=.false.
			if (.not.isconverged .and. iterconverg==model%parameters%it_max) ismaxitreached = .true.
			iterconvergmax=max(iterconverg,iterconvergmax)
		end do WFUNSAT1

		!---- 09.04 ---- estimate hnew for new timestep
		do is=1,nsat
			call model%sat(is)%calc%estimate_hnew_for_new_timestep()
		end do

		!---- 09.05 ---- check if WF1DSAT is activated (a water-table appears)
		!---- 09.06 ---- if first time activated-> set initial conditions to 0 in WF1DSAT
		do is=1,nsat
			if(model%constraints%get_sum_hsatv_on_is(is)>0.0_dpd) then
				if  (.not.ishsatover0(is)) then
					!write(*,*) "C Calculating saturated layers"
					model%sat(is)%calc%nodes%hnew=0.0_dpd
					model%sat(is)%calc%nodes%hold=0.0_dpd
					model%sat(is)%calc%nodes%qent=0.0_dpd
					call model%sat(is)%calc%get_results_nodes()
					call model%sat(is)%calc%get_results_elements()
					isfirsthsat(is)=.true.
				else
					isfirsthsat(is)=.false.
				end if
				ishsatover0(is)=.true.	!1DSAT(is): activated
			else
				ishsatover0(is)=.false. !1DSAT(is): deactivated
				isfirsthsat(is)=.false.
			end if
		end do

		!---- 09.06 ---- run WF1DSAT: With qvtb_sat and nrel from WF1DUNSAT
		!---- 09.06.01 ----  put nrel and qvtb from constraints to nodes
		do is=1,nsat
			call model%sat(is)%set_qent_from_constraint_to_nodes() 
			call model%sat(is)%set_nrel_from_constraint_to_nodes()
		end do

		!---- 09.06.02 ----  SAT_CONVERGENCE: do while not converge or dtminor itmax not reached
		!(CAN BE DONE IN PARALLEL)
		itersat=0
		do is=1,nsat
			if(ishsatover0(is)) then
				isconverged = .false.
				iterconverg = 0
				do while(.not.(isconverged .or. (model%time%dt<=model%parameters%dtmin .and. iterconverg==(model%parameters%it_max-1))))
					model%time%iter_total = model%time%iter_total + 1
					iterconverg = iterconverg+1
					itertotal=itertotal+1
					call model%sat(is)%calc%iterate(isconverged)
					if(.not.isconverged .and. iterconverg==model%parameters%it_max) exit
				end do
				itersat = max(iterconverg,itersat)
				if (.not.isconverged) isconvergedall=.false.
				if (.not.isconverged .and. iterconverg==model%parameters%it_max) ismaxitreached = .true.
				iterconvergmax=max(iterconverg,iterconvergmax)
			end if
		end do

		!---- 09.08 ----  If not converged reduce timestep, revert to old and restart in 09
		if (.not.isconvergedall) then
			write(*,*) 'C REVERT TO OLD: Not converge...'
			model%time%checkprint = .false.
			call model%time%factor_timestep(1/3.0_dpd) !Update timestep with boundaries on tmin and tmax (dt=factor·dtold)
			do iu=1,nunsat
				call model%unsat(iu)%calc%revert_to_old()
			end do
			do is=1,nsat
				call model%sat(is)%calc%revert_to_old()
				if (isfirsthsat(is)==.true.) then
					isfirsthsat(is)= .false.
					ishsatover0(is)=.false.
				end if
			end do
			cycle
		end if

		!---- 09.09 ---- run WF1DUNSAT: With imposed watertable height from WF1DSAT
		!---- 09.09.01 ---- impose the water-table height on constraints from WF1DSAT (as Dirichlet boundary)
		call model%constraints%set_dirichlet_from_hsat_h(1.0_dpd)

		!---- 09.09.02 ---- UNSAT_CONVERGENCE: do while not converge or dtminor itmax not reached
		!(CAN BE DONE IN PARALLEL)
		isconvergedall = .true.
		ismaxitreached = .false.
		do iu=1,nunsat
			isconverged = .false.
			iterconverg = 0
			do while(.not.(isconverged .or. (model%time%dt<=model%parameters%dtmin .and. iterconverg==(model%parameters%it_max-1))))
				model%time%iter_total = model%time%iter_total + 1
				iterconverg = iterconverg+1
				iterunsat = iterunsat+1
				call model%unsat(iu)%calc%update_th_from_h(CONSIDER_HNEW)
				itertotal=itertotal+1
				call model%unsat(iu)%calc%iterate(isconverged)
				if (iterconverg<model%parameters%it_min) isconverged=.false.
				if(.not.isconverged .and. iterconverg==model%parameters%it_max) exit
			end do
			if (.not.isconverged) isconvergedall=.false.
			if(.not.isconverged .and. iterconverg==model%parameters%it_max) ismaxitreached = .true.
			iterconvergmax=max(iterconverg,iterconvergmax)
		end do

		!---- 09.10 ---- Check if increment on hsat is within the tolerance. If not set converged to false
		maxhsat = maxval(model%constraints%get_hsatv_mat())
		inchsat = abs(maxhsat-maxhsatold)
		if(inchsat>model%parameters%maxhsatinc) isconvergedall=.false.

		!---- 09.11 ---- If not converged reduce timestep, revert to old and restart in 09
		if (.not.isconvergedall) then
			if(inchsat>model%parameters%maxhsatinc) write(*,*) 'C REVERT TO OLD: inchsat>incmax...'
			model%time%checkprint = .false.
			call model%time%factor_timestep(1/3.0_dpd) !Update timestep with boundaries on tmin and tmax (dt=factor·dtold)
			do iu=1,nunsat
				call model%unsat(iu)%calc%revert_to_old()
			end do
			do is=1,nsat
				call model%sat(is)%calc%revert_to_old()
			end do
			cycle
		end if

		!---- 09.12 ---- Model CONVERGES!: Update results in nodes and elements in WF1DUNSAT
		do iu=1,nunsat
			call model%unsat(iu)%calc%get_results_nodes()
			call model%unsat(iu)%calc%get_results_elements()
			call model%unsat(iu)%calc%get_results_nodes()
			!call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt)
		end do
		!---- 09.13 ---- Model CONVERGES!: Update results in nodes and elements in WF1DSAT
		do is=1,nsat
			call model%sat(is)%calc%get_results_elements()
			call model%sat(is)%calc%get_results_nodes()
			call model%sat(is)%calc%get_results_elements()
			call model%sat(is)%calc%get_results_nodes()
			!call model%sat(is)%put_results_in_constraints()
		end do

		!---- 09.14 ---- Update the value of hsat in constraints
		do iu=1,nunsat
			call model%unsat(iu)%constraints%update_hsat(model%unsat(iu)%calc%elements)
		end do

		!---- 09.15 ---- Update the value of nrel from WF1DUNSAT
		do iu=1,nunsat
			call model%constraints%update_nrel_to_fit_inc_hunsat(model%time%dt)
		end do
		!---- 09.16 ---- Update the value of nrel from WF1DUNSAT
		do iu=1,nunsat
			call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt)
		end do
		do is=1,nsat
			call model%sat(is)%put_results_in_constraints()
		end do

		!---- 09.17 ---- If t=tprint: Print output
		if (model%time%checkprint) then
			printchar ='P'
			do is=1,nsat
				call model%sat(is)%print_timestep(33,34,is)
			end do
			do iu=1,nunsat
				call model%unsat(iu)%print_timestep(31,32,iu)
			end do
			model%time%checkprint = .false.
		else
			printchar = '·'
		end if
		call model%print_alltimes(51)

		!Write console...
		!					printchar		,itermodel		,iterunsat,		,itersat			,t,							,dt							,hnew_u																															,hsat_u															,hsat_s														,qv																	,dqhdx_s
		!WRITE(*,'(A2					," IM: ",i3.3	," IU: ",i3.3	," IS: ",i3.3	," t: ", f10.4	," dt: "	,E10.3	," hu(1,2): ",E10.3																								," hsatv(1,2): ",E10.3							," hsath(1,2): ",E10.3						," qverv(1,2): ",E10.3							," dqdxh: ",E10.3,A40)') &
		!	&				printchar		,itermodel		,iterunsat		,itersat			, model%time%t	,model%time%dt	,model%unsat(2)%calc%nodes%hnew(model%constraints%get_idnodev(1,2))	, model%constraints%get_hsatv(1,2)	, model%constraints%get_hsath(1,2), model%constraints%get_qverv(1,2)	, model%constraints%get_dqhordxh(1,2), trim(debugmsg_txt)
		WRITE(*,'(A2					," IM: ",i4.4	," IU: ",i4.4	," IS: ",i4.4	," t: ", E10.3	," dt: "	,E10.3	)') &
		& printchar		,itermodel		,iterunsat		,itersat			, model%time%t	,model%time%dt	
		
		
		!---- 09.18 ---- Update timestep dt, and set old values from calculated new values
		call model%time%update_dt(iterconvergmax)
		!set all values to old
		do iu=1,nunsat
			call model%unsat(iu)%calc%set_old()
		end do
		do is=1,nsat
			call model%sat(is)%calc%set_old()
		end do
		maxhsatold = maxval(model%constraints%get_hsatv_mat())

		!Log statistics
		itersteps=itersteps+1
		totaldt=totaldt+model%time%dt
		mindt=min(mindt,model%time%dt)
		maxdt=max(maxdt,model%time%dt)

	end do TIMESTEPPING

	close(31)
	close(32)
	close(33)
	close(34)
	close(51)

	!Print and log statistics...
	call CPU_TIME(cputimeend)
	call SYSTEM_CLOCK(timemilisec,countrate)
	systimeend=real(timemilisec,8)/real(countrate,8)
	meandt=totaldt/real(itersteps,8)

	open (33,file=FILEOUTPUT_STATS,status='unknown')
	write(33,'("itertotal,nsteps,cputime_s,systime_s,mintimestep,meantimestep,maxtimestep")')
	write(33,'(i6,",",i6,",",F10.3,",",F10.3,",",E10.3,",",E10.3,",",E10.3)') itertotal,itersteps,cputimeend-cputimeinit,systimeend-systimeinit,mindt,meandt,maxdt
	close(33)

999 end program WF15DSATUNSAT