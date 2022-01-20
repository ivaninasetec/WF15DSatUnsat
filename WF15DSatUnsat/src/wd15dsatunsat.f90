	program model_flow1d
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

	!---- 01 ---- read input files
	nargs = COMMAND_ARGUMENT_COUNT()
	if (nargs==0) then
		WRITE(*,*) 'Pease include input file as the first argument argument'
		goto 999 !If there is no argument then end programm
	end if

	call GET_COMMAND_ARGUMENT(1,fileinput) !Input file is the first argument of the program
	!write(*,*) 'number of arguments',nargs
	npos = scan(trim(fileinput),".", BACK= .true.) !To remove extension
	if (npos>0) namefile = trim(fileinput(1:npos-1)) !Now without extension

	fileboundary					=trim(namefile)//'.wfbound'
	fileoutput_nodes_v		=trim(namefile)//'.outnodu'
	fileoutput_elements_v	=trim(namefile)//'.outelmu'
	fileoutput_nodes_h		=trim(namefile)//'.outnods'
	fileoutput_elements_h	=trim(namefile)//'.outelms'
	fileoutput_constraints	=trim(namefile)//'.outcons'


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

	!*** set intial conditions...
	do iu=1,nunsat
		!Set time to initial and hnew, htemp and hold to hinit
		call model%unsat(iu)%calc%set_to_initial()
		call model%unsat(iu)%calc%update_th_from_h(CONSIDER_ALL)

		!Set Dirichlet into first node
		if (model%layers%bottombyphl) call model%unsat(iu)%calc%set_dirichlet_to_node(1,model%unsat(iu)%calc%nodes%hnew(1))
		!Construct matrix not depending on time
		call model%unsat(iu)%calc%build_linearsystem(IS_NOT_TIME_DEPENDANT,CONSIDER_HTEMP)
	end do

	do is=1,nsat
		!Set time to initial and hnew, htemp and hold to hinit
		call model%sat(is)%calc%set_to_initial()
		call model%sat(is)%calc%get_results_elements()
		call model%sat(is)%calc%get_results_nodes()
		!Construct matrix not depending on time
		call model%sat(is)%calc%build_linearsystem(IS_NOT_TIME_DEPENDANT,CONSIDER_HTEMP)
	end do

	do iu=1,nunsat
		call model%unsat(iu)%calc%get_results_elements()
		call model%unsat(iu)%calc%get_results_nodes()
		call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt) !CHECK
		!call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt,10.0_dpd) !CHECK
	end do

	if(.not.allocated(ishsatover0)) allocate(ishsatover0(nsat))
	if(.not.allocated(isfirsthsat)) allocate(isfirsthsat(nsat))
	isfirsthsat = .false.
	ishsatover0 = .false.
	model%time%iter_total = 0

	!*** print initial time
	do iu=1,nunsat
		call model%unsat(iu)%print_timestep(31,32,iu)
	end do
	do is=1,nsat
		call model%sat(is)%print_timestep(33,34,is)
	end do
	call model%print_alltimes(51)

	!BEGIN TIMESTEPPING -------------------------------------------------------------------------------------------------

	TIMESTEPPING: do while(model%time%t < model%parameters%Tmax)
		!*** increase time from told (t=t+dt or tprint or tmax) (checkprint=true if tprint)
		call model%time%increase_time()

		!*** update boundaries properties on top (dont apply yet)
		if (model%layers%topboundbyh) htop =  model%boundary%get_hbound_file(model%time%t)
		if (model%layers%topboundbyq) qtop =  model%boundary%get_qbound_file(model%time%t)

		!Estimate hnew to start iterations...
		do iu=1,nunsat
			call model%unsat(iu)%calc%estimate_hnew_for_new_timestep()
		end do

		do is=1,nsat
			call model%sat(is)%calc%estimate_hnew_for_new_timestep()
		end do

		!*** 1.5DSATUNSAT: maxhsattemp
		!maxhsattemp = maxval(model%constraints%get_hsatv_mat())

		!*** 1DUNSAT: apply Newmann from interfaces and deactivate Dirichlet...
		!boundary properties on interfaces
		call model%constraints%set_newmann_from_sat_h(1.0_dpd)
		do iu=1,nunsat
			!boundary properties on top
			if (model%layers%topboundbyh) call model%unsat(iu)%calc%set_dirichlet_to_node(model%unsat(iu)%calc%nodes%count,htop)
			if (model%layers%topboundbyq) call model%unsat(iu)%calc%set_newmann_to_node(model%unsat(iu)%calc%nodes%count,qtop)
		end do
		
		!apply Newmann on nodes with isnewmann
		do iu=1,nunsat
			call model%unsat(iu)%calc%assign_newmann()
		end do
		
		!*** CHECK IF HSAT IS ACTIVATED (hsat begin to appear in unsat and put initial conditions)...
		do is=1,nsat
			if(model%constraints%get_sum_hsatv_on_is(is)>0.0_dpd) then
				if  (.not.ishsatover0(is)) then
					!write(*,*) "Calculating saturated layers"
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

		!*** assign infiltration load from interfaces on 1DSAT
		do is=1,nsat
			call model%sat(is)%set_qent_from_constraint_to_nodes() !Check the value of qent in sat(is)%constraints%qent(iu) and interpolate linearly along x in 1DSAT 
			call model%sat(is)%set_nrel_from_constraint_to_nodes() !Check the value of nrel in sat(is)%constraints%qent(iu) and interpolate linearly along x in 1DSAT 
		end do

		!PARALLEL COSIMULATION ----------------------------------------------------------------------------------------
		!*** 1DFLOW_CONVERGENCE: 1DUNSAT...(CAN BE DONE IN PARALLEL)
		isconvergedall = .true.
		ismaxitreached = .false.
		iterconvergmax = 0
		iterunsat = 0
		do iu=1,nunsat
			isconverged = .false.
			iterconverg = 0
			do while(.not.(isconverged .or. (model%time%dt<=model%parameters%dtmin .and. iterconverg==(model%parameters%it_max-1))))
				model%time%iter_total = model%time%iter_total + 1
				iterconverg = iterconverg+1
				call model%unsat(iu)%calc%update_th_from_h(CONSIDER_HNEW)
				!call model%unsat(iu)%calc%iterate(isconverged,isconvergedall)
				call model%unsat(iu)%calc%iterate(isconverged)
				if (iterconverg<model%parameters%it_min) isconverged=.false.
				if(.not.isconverged .and. iterconverg==model%parameters%it_max) exit
			end do
			iterunsat = max(iterconverg,iterunsat)
			if (.not.isconverged) isconvergedall=.false.
			if (.not.isconverged .and. iterconverg==model%parameters%it_max) ismaxitreached = .true.
			iterconvergmax=max(iterconverg,iterconvergmax)
		end do

		!!*** 1DFLOW_CONVERGENCE: 1DSAT...	(CAN BE DONE IN PARALLEL)
		itersat=0	
		do is=1,nsat
			if(ishsatover0(is)) then
				isconverged = .false.
				iterconverg = 0
				do while(.not.(isconverged .or. (model%time%dt<=model%parameters%dtmin .and. iterconverg==(model%parameters%it_max-1))))
					model%time%iter_total = model%time%iter_total + 1
					iterconverg = iterconverg+1
					!call model%sat(is)%calc%iterate(isconverged,isconvergedall)
					call model%sat(is)%calc%iterate(isconverged)
					if(.not.isconverged .and. iterconverg==model%parameters%it_max) exit
					!if (iterconverg<model%parameters%it_min) isconverged=.false.
				end do
				itersat = max(iterconverg,itersat)
				if (.not.isconverged) isconvergedall=.false.
				if (.not.isconverged .and. iterconverg==model%parameters%it_max) ismaxitreached = .true.
				iterconvergmax=max(iterconverg,iterconvergmax)
			end if
		end do

		!*** NOT CONVERGED: reduce timestep, revert to old and cycle
		if (.not.isconvergedall) then
			write(*,*) 'REVERT TO OLD: Not converge...'
			model%time%checkprint = .false.
			call model%time%factor_timestep(1/3.0_dpd) !Update timestep with boundaries on tmin and tmax (dt=factor·dtold)
			do iu=1,nunsat
				call model%unsat(iu)%calc%revert_to_old()
				!call model%unsat(iu)%calc%get_results_nodes()
				!call model%unsat(iu)%calc%get_results_elements()
				!call model%unsat(iu)%calc%get_results_nodes()
				!call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt)
			end do
			do is=1,nsat
				call model%sat(is)%calc%revert_to_old()
				if (isfirsthsat(is)==.true.) then
					isfirsthsat(is)= .false.
					ishsatover0(is)=.false.
				end if
				!call model%sat(is)%calc%get_results_nodes()
				!call model%sat(is)%calc%get_results_elements() !Need nodes to be calculated previously
				!call model%sat(is)%put_results_in_constraints() !Update: hsat_mean, qent_mean, incvoldt_mean, dqhordx_mean,  dqhordx_all_mean
			end do
			cycle
		end if

		!ADJUST hsat,u to hsat,s ----------------------------------------------------------------------------------------
		!*** 1DUNSAT: set dirichlet from hsat,s (and deactivate newmann)
		call model%constraints%set_dirichlet_from_hsat_h(1.0_dpd)
  
		!*** 1DFLOW_CONVERGENCE: 1DUNSAT...(CAN BE DONE IN PARALLEL)
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
				call model%unsat(iu)%calc%iterate(isconverged)
				!call model%unsat(iu)%calc%iterate(isconverged,isconvergedall)
				if (iterconverg<model%parameters%it_min) isconverged=.false.
				if(.not.isconverged .and. iterconverg==model%parameters%it_max) exit
			end do
			if (.not.isconverged) isconvergedall=.false.
			if(.not.isconverged .and. iterconverg==model%parameters%it_max) ismaxitreached = .true.
			iterconvergmax=max(iterconverg,iterconvergmax)
		end do

		!Check if increment on hsat is not too much...
		maxhsat = maxval(model%constraints%get_hsatv_mat())
		!epshsat = abs(maxhsat-maxhsattemp)
		inchsat = abs(maxhsat-maxhsatold)
		if(inchsat>model%parameters%maxhsatinc) isconvergedall=.false.

		!*** NOT CONVERGED: reduce timestep, revert to old and cycle
		if (.not.isconvergedall) then
			if(inchsat>model%parameters%maxhsatinc) write(*,*) 'REVERT TO OLD: inchsat>incmax...'
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

		! CONVERGED ----------------------------------------------------------------------------------
		!update results in nodes and elements:
		do iu=1,nunsat
			call model%unsat(iu)%calc%get_results_nodes()
			call model%unsat(iu)%calc%get_results_elements()
			call model%unsat(iu)%calc%get_results_nodes()
			!call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt)
		end do
		!on 1DSAT
		do is=1,nsat
			call model%sat(is)%calc%get_results_elements()
			call model%sat(is)%calc%get_results_nodes()
			call model%sat(is)%calc%get_results_elements()
			call model%sat(is)%calc%get_results_nodes()
			!call model%sat(is)%put_results_in_constraints()
		end do
		
		!Update hsat in constraints...
		do iu=1,nunsat
		call model%unsat(iu)%constraints%update_hsat(model%unsat(iu)%calc%elements)
		end do
		
		!Get new value for nrel... !This is nrel=(q'v,u_old-Q'newman,u_new)/inchunsat_new
		do iu=1,nunsat
		call model%constraints%update_nrel_to_fit_inc_hunsat(model%time%dt) !This is nrel=(q'u-Q'u)/inchunsat
		end do
		!Update values on constraints
		do iu=1,nunsat
			call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt)
			!call model%unsat(iu)%constraints%update_all(model%unsat(iu)%calc,model%time%dt,1.0_dpd) !CHECK
		end do
		do is=1,nsat
			call model%sat(is)%put_results_in_constraints()
		end do
				
		
		

		!Print...
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
			printchar = ' '
		end if
		call model%print_alltimes(51)

		!!Write console...
		!!					printchar		,itermodel		,iterunsat,		,itersat			,t,							,dt							,hnew_u																															,hsat_u															,hsat_s														,qv																	,dqhdx_s
		!WRITE(*,'(A1					,"IM:",i3.3		," IU:",i3.3	," IS:",i3.3	," t: ", f10.4	," dt: "	,E10.3	," hu(1,2): ",E10.3																								," hsatv(1,2): ",E10.3							," hsath(1,2): ",E10.3						," qverv(1,2): ",E10.3							," dqdxh: ",E10.3)') &
		!	&				printchar		,itermodel		,iterunsat		,itersat			, model%time%t	,model%time%dt	,model%unsat(2)%calc%nodes%hnew(model%constraints%get_idnodev(1,2))	, model%constraints%get_hsatv(1,2)	, model%constraints%get_hsath(1,2), model%constraints%get_qverv(1,2)	, model%constraints%get_dqhordxh(1,2)

		!Write console...
		!					printchar		,itermodel		,iterunsat,		,itersat			,t,							,dt							,hnew_u																															,hsat_u															,hsat_s														,qv																	,dqhdx_s
		WRITE(*,'(A2					," IM: ",i3.3	," IU: ",i3.3	," IS: ",i3.3	," t: ", f10.4	," dt: "	,E10.3	," hu(1,2): ",E10.3																								," hsatv(1,2): ",E10.3							," hsath(1,2): ",E10.3						," qverv(1,2): ",E10.3							," dqdxh: ",E10.3,A40)') &
			&				printchar		,itermodel		,iterunsat		,itersat			, model%time%t	,model%time%dt	,model%unsat(2)%calc%nodes%hnew(model%constraints%get_idnodev(1,2))	, model%constraints%get_hsatv(1,2)	, model%constraints%get_hsath(1,2), model%constraints%get_qverv(1,2)	, model%constraints%get_dqhordxh(1,2), trim(debugmsg_txt)

		
		!*** update dt depending on iterations:
		call model%time%update_dt(iterconvergmax)

		!*** set all values to old
		do iu=1,nunsat
			call model%unsat(iu)%calc%set_old()
		end do
		do is=1,nsat
			call model%sat(is)%calc%set_old()
		end do
		maxhsatold = maxval(model%constraints%get_hsatv_mat())

	end do TIMESTEPPING

	close(31)
	close(32)
	close(33)
	close(34)
	close(51)

999	end program model_flow1d