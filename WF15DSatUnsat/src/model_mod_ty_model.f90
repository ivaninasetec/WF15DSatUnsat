	!********************************************************************************************************************
	!        CLASS THAT INCLUDE THE WHOLE MODEL
  !********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : mod_sat_ty_nodes
	! URL           : ...
	! AFFILIATION   : ...
	! DATE          : ...
	! REVISION      : ... V 0.0
	! LICENSE				: This software is copyrighted 2019(C)
	!> @author
	!> Iván Campos-Guereta Díez
  !  MSc Civil Engineering by Polytechnic University of Madrid                                                     *
  !  PhD Student by University of Nottingham                                                                       *
  !  eMBA by International Institute San Telmo in Seville                                                          *
  !  ivan.camposguereta@nottingham.ac.uk  
	! DESCRIPTION:
	!> Class for the collection of nodes-classes in the saturated model            
  !********************************************************************************************************************
	
	module model_mod_ty_model
		use com_mod_ty_time,					only: ty_com_time
		use com_mod_ty_material,			only: ty_com_material
		use com_mod_ty_parameters,		only: ty_com_parameters
		use com_mod_ty_boundary,			only: ty_com_boundary
		use com_mod_ty_layers,				only: ty_com_layers
		
		use model_mod_ty_mesh,				only: ty_model_mesh
		use model_mod_ty_constraints,	only: ty_model_constraints
		use unsat_mod_ty_model,				only: ty_unsat_model
		use sat_mod_ty_model,					only: ty_sat_model
	
  implicit none
  include 'inc_precision.fi'

  private
	
	  type,public::ty_model
			type(ty_com_parameters)::						parameters
			type(ty_com_material),allocatable::	material(:)
			type(ty_com_layers)::layers
			type(ty_model_mesh)::mesh
			type(ty_com_boundary)::boundary
			type(ty_model_constraints)::constraints
			
			type(ty_com_time)::time
			type(ty_unsat_model)	,allocatable::unsat(:)
			type(ty_sat_model)		,allocatable::sat(:)
		contains
			procedure,public:: readfileinput=> s_model_model_readfileinput
			procedure,public:: construct=> s_model_model_construct
			!procedure,public:: pass_unsat_to_constraint_to_sat => s_model_model_pass_unsat_to_constraints_to_sat
			!procedure,public:: pass_sat_to_constraint_to_unsat => s_model_model_pass_sat_to_constraints_to_unsat
			procedure,public:: print_alltimes => s_model_model_print_alltimes
			!procedure,public:: print_timestep=> s_unsat_model_print_timestep
			
		end type ty_model
		
	contains

	!---------------------------------------------------------------------------------------------------------------------
	! READFILEINPUT
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to construct the model from the inputs.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------
	
	subroutine s_model_model_readfileinput(this,fileinput,fileboundary)	
	use model_mod_inputs, only: s_model_inputs_parameters,s_model_inputs_materials,s_model_inputs_layers,s_model_inputs_mesh,s_model_inputs_boundary

	class(ty_model),target::this
	character*400,intent(in)::fileinput
	character*400,intent(in)::fileboundary
	
	call s_model_inputs_parameters(this%parameters,fileinput)
	call s_model_inputs_materials(this%material,fileinput)
	call s_model_inputs_layers(this%layers,fileinput,this%material)
	call s_model_inputs_mesh(this%mesh,fileinput,this%layers)
	call s_model_inputs_boundary(fileboundary,this%boundary,this%layers)
		
	end subroutine s_model_model_readfileinput	
	
	!---------------------------------------------------------------------------------------------------------------------
	! CONSTRUCTOR
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to construct the model from the inputs.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------
	
	subroutine s_model_model_construct(this)
	use sat_mod_inputs, only:s_sat_inputs_nodes,s_sat_inputs_elements
	use unsat_mod_inputs, only:s_unsat_inputs_nodes,s_unsat_inputs_elements
	
	class(ty_model),target::this
	integer::i,nh,nv,l,nsat,isat1,j
		
	!Construct all instances in model from input file...
	if(this%layers%bottombyphl) then
	nsat = this%layers%count-1
	isat1 = 2
	else
	nsat = this%layers%count
	isat1 = 1
	end if
	
	if(.not.allocated(this%sat))		allocate(this%sat		(nsat))
	if(.not.allocated(this%unsat))	allocate(this%unsat	(this%mesh%vmod_count))
	

	
	!FILL INPUTS IN SATURATED ELEMENTS:
	do i=1, nsat
		this%sat(i)%parameters = this%parameters
		if(.not.allocated(this%sat(i)%material)) allocate(this%sat(i)%material(size(this%material)))
		this%sat(i)%material = this%material
		!this%sat(i)%boundary= this%boundary !This is not a must
		
		!Fill layers:
		call this%sat(i)%layers%allocateall(nsat-i)
		this%sat(i)%layers%count = nsat-i+1
		this%sat(i)%layers%height = this%layers%height(i+isat1-1:this%layers%count)
		this%sat(i)%layers%width = this%layers%width
		this%sat(i)%layers%material = this%layers%material(i+isat1-1:this%layers%count)
		this%sat(i)%layers%htop = this%layers%htop(i+isat1-1:this%layers%count)-this%layers%hbottom(i+isat1-1)
		this%sat(i)%layers%hbottom = this%layers%hbottom(i+isat1-1:this%layers%count)-this%layers%hbottom(i+isat1-1)
		this%sat(i)%layers%slope = this%layers%slope(i+isat1-1:this%layers%count)
		this%sat(i)%layers%topboundbyh = this%layers%topboundbyh
		this%sat(i)%layers%topboundbyq = this%layers%topboundbyq

		this%sat(i)%calc%parameters => this%parameters
		this%sat(i)%calc%time => this%time
		this%sat(i)%calc%layers => this%sat(i)%layers
		!this%sat(i)%calc%boundary =>	this%boundary
		
		!FILL MESH ON SAT:
		call this%sat(i)%mesh%allocateall(this%mesh%vmod_count)
		this%sat(i)%mesh%nn = this%mesh%nn
		this%sat(i)%mesh%nc = this%mesh%nc
		this%sat(i)%mesh%nnc = this%mesh%nnc
		this%sat(i)%mesh%nelemh_count = this%mesh%nelemh_count
		this%sat(i)%mesh%nelemh = this%mesh%nelemh
		this%sat(i)%mesh%nnodh_count = this%mesh%nnodh_count
		this%sat(i)%mesh%nnodclassh_count = this%mesh%nnodclassh_count
		this%sat(i)%mesh%width = this%mesh%width
		this%sat(i)%mesh%vmod_count = this%mesh%vmod_count
		this%sat(i)%mesh%vmod_x = this%mesh%vmod_x
		this%sat(i)%mesh%vmod_idnod = this%mesh%vmod_idnod
		!this%sat(i)%mesh%vmod_hsat =	0.0_dpd
		!this%sat(i)%mesh%vmod_hsattemp =	0.0_dpd
		!this%sat(i)%mesh%vmod_hsatold =	0.0_dpd
		!this%sat(i)%mesh%vmod_qent =	0.0_dpd
		!this%sat(i)%mesh%vmod_qenttemp =	0.0_dpd
		!this%sat(i)%mesh%vmod_qentold =	0.0_dpd
		
		call this%sat(i)%calc%allocateall()
	
		call s_sat_inputs_nodes(this%sat(i)%calc%nodes,this%sat(i)%mesh,this%sat(i)%layers)	
		call s_sat_inputs_elements(this%sat(i)%calc%elements,this%sat(i)%mesh,this%sat(i)%calc%nodes,this%sat(i)%layers%material(1))
		call this%sat(i)%calc%construct()
		!call this%sat(i)%mesh%construct_interpmatrix(this%sat(i)%calc%nodes)
		call this%sat(i)%constraints%allocateall(this%mesh%vmod_count,this%mesh%nnodclassh_count)
		call this%sat(i)%constraints%construct(this%mesh%vmod_idnod,this%sat(i)%calc%nodes)
		
	end do
	
	!FILL INPUTS IN UNSATURATED ELEMENTS:
	do i=1, this%mesh%vmod_count
		this%unsat(i)%parameters = this%parameters
		if(.not.allocated(this%unsat(i)%material)) allocate(this%unsat(i)%material(size(this%material)))
		this%unsat(i)%material = this%material
		!this%sat(i)%boundary= this%boundary !This is not a must
		
		!Fill layers:
		!this%unsat(i)%calc%layers => this%layers
		call this%unsat(i)%layers%allocateall(this%layers%count)
		this%unsat(i)%layers%height = this%layers%height
		this%unsat(i)%layers%width = this%layers%width
		this%unsat(i)%layers%material = this%layers%material
		this%unsat(i)%layers%htop = this%layers%htop
		this%unsat(i)%layers%hbottom = this%layers%hbottom
		this%unsat(i)%layers%zphr = this%layers%zphr
		this%unsat(i)%layers%topboundbyh = .false.
		this%unsat(i)%layers%topboundbyq = .true.

		this%unsat(i)%calc%parameters => this%parameters
		this%unsat(i)%calc%time => this%time
		this%unsat(i)%calc%layers => this%unsat(i)%layers
		this%unsat(i)%calc%boundary =>	this%boundary
		
		!Fill unsaturated mesh...
		call this%unsat(i)%mesh%allocateall(this%layers%count)
		this%unsat(i)%mesh%nn = this%mesh%nn
		this%unsat(i)%mesh%nc = this%mesh%nc
		this%unsat(i)%mesh%nnc = this%mesh%nnc
		this%unsat(i)%mesh%nelemv_count = this%mesh%nelemv_count
		this%unsat(i)%mesh%nelemv = this%mesh%nelemv
		this%unsat(i)%mesh%nnodv_count = this%mesh%nnodv_count
		this%unsat(i)%mesh%nnodclassv_count = this%mesh%nnodclassv_count
		this%unsat(i)%mesh%height = this%mesh%height

		
		!this%unsat(i)%mesh%hsat = 0.0_dpd
		!this%unsat(i)%mesh%hsatold = 0.0_dpd
		!this%unsat(i)%mesh%idnodsat = 0
		!this%unsat(i)%mesh%idnodlayer(1)=1
		!do l=2,this%layers%count
		!this%unsat(i)%mesh%idnodlayer(l)=this%unsat(i)%mesh%idnodlayer(l-1)+this%unsat(i)%mesh%nelemv(l)	
		!end do
		
	
	!Construct mesh:
	do j=1, size(this%mesh%nelemv)
		this%unsat(i)%mesh%idnode_bottom(j) = sum(this%mesh%nelemv(1:j-1))+1
		this%unsat(i)%mesh%idnode_top(j) = sum(this%mesh%nelemv(1:j))
	end do

	do j=1, size(this%mesh%nelemv) 
		this%unsat(i)%mesh%idnode_bottom(j) = sum(this%mesh%nelemv(1:j-1))+1
		this%unsat(i)%mesh%idnode_top(j) = sum(this%mesh%nelemv(1:j))+1
	end do

	do j=1, size(this%mesh%nelemv)
		this%unsat(i)%mesh%idelem_bottom(j) = sum(this%mesh%nelemv(1:j-1))+1
		this%unsat(i)%mesh%idelem_top(j) = sum(this%mesh%nelemv(1:j))
	end do
	
		!this%unsat(i)%mesh%incvsat=0.0_dpd
		!this%unsat(i)%mesh%qvsat=0.0_dpd
		!do l=1,this%layers%count
		!this%unsat(i)%mesh%nmean=this%layers%material(l)%thsat-this%layers%material(l)%thres	
		!end do	
		!this%unsat(i)%mesh%nrel=1.0_dpd
		!this%unsat(i)%mesh%dqhordx=0.0_dpd
		
		!Calculat nodes of each layer
		call this%unsat(i)%calc%allocateall()
	
		call s_unsat_inputs_nodes(this%unsat(i)%calc%nodes,this%unsat(i)%mesh,this%unsat(i)%layers)	
		call s_unsat_inputs_elements(this%unsat(i)%calc%elements,this%unsat(i)%mesh,this%unsat(i)%calc%nodes,this%unsat(i)%layers)
		call this%unsat(i)%calc%construct()
		
		call this%unsat(i)%constraints%allocateall(nsat,this%unsat(i)%calc%nodes%count)
		!call this%unsat(i)%constraints%construct(this%unsat(i)%mesh,this%unsat(i)%calc%nodes)
    call this%unsat(i)%constraints%construct(this%unsat(i)%mesh%idnode_bottom(isat1:),this%unsat(i)%calc%nodes)
    
		
	end do	
	
	!Update layers idelements
	do l=1, size(this%mesh%nelemv)
		this%layers%id_elem_bottom(l) = sum(this%mesh%nelemv(1:l-1))+1
		this%layers%id_elem_top(l) = sum(this%mesh%nelemv(1:l))
	end do
	
	!constraints module
	call this%constraints%construct(this%unsat,this%sat)
	
	!!Allocate constraints and initialize
	!call this%constraints%allocateall(this%layers%count-1,this%mesh%vmod_count)
	!
	!do nh=1,this%layers%count-1
	!	do nv=1,this%mesh%vmod_count
	!		this%constraints%idnode_h(nh,nv) = this%sat(1)%mesh%vmod_idnod(nv) !Id node for the 
	!		this%constraints%idnode_v(nh,nv) = sum(this%unsat(1)%mesh%nelemv(1:nh))+1
	!	end do
	!end do
	!
	!this%constraints%hsat_h = 0.0_dpd
	!this%constraints%hsat_v = 0.0_dpd
	!this%constraints%dqhordx_h = 0.0_dpd
	!this%constraints%qver_v = 0.0_dpd
	!
	!this%constraints%unsat => this%unsat
	!this%constraints%sat => this%sat
	
	this%time%parameters => this%parameters
	
	end subroutine s_model_model_construct		
	
	
	!!---------------------------------------------------------------------------------------------------------------------
	!! PASS_UNSAT_TO_CONSTRAINTS_TO_SAT
	!!---------------------------------------------------------------------------------------------------------------------
	!!> @author Iván Campos-Guereta Díez
	!!> @brief
	!!> Subroutine to construct the model from the inputs.
	!!> @param[in] ntype
	!!> @param[inout] mx
	!!> @param[in] option
	!!---------------------------------------------------------------------------------------------------------------------
	!
	!subroutine s_model_model_pass_unsat_to_constraints_to_sat(this)
	!
	!class(ty_model),target::this
	!
	!integer::iu,is,nunsat,nsat,i
	!
	!!Number of unsat and sat elements:
	!nunsat=this%mesh%vmod_count
	!nsat = this%layers%count-1
	!
	!!Fill constraint class with unsat data:
	!do iu=1,nunsat
	!call this%unsat(iu)%set_hsat_to_mesh()
	!call this%unsat(iu)%set_satvalues_to_mesh()
	!
	!this%constraints%hsat_v(1:nsat,iu)			= this%unsat(iu)%mesh%hsat(2:nsat+1)
	!this%constraints%qver_v(1:nsat,iu)			= this%unsat(iu)%mesh%qvsat(2:nsat+1)
	!this%constraints%nrel(1:nsat,iu)				= this%unsat(iu)%mesh%nrel(2:nsat+1)	
	!end do
 !
	!!Set vertical flow to each vertical module in the sat model:
	!do is=1,nsat
	!	do iu=1,nunsat
	!	this%sat(is)%mesh%vmod_qent(iu)=this%constraints%qver_v(is,iu)
	!	this%sat(is)%mesh%vmod_nrel(iu)=this%constraints%nrel(is,iu)
	!	end do
	!end do
	!
	!end subroutine s_model_model_pass_unsat_to_constraints_to_sat
 !
	!!---------------------------------------------------------------------------------------------------------------------
	!! PASS_UNSAT_TO_CONSTRAINTS_TO_SAT
	!!---------------------------------------------------------------------------------------------------------------------
	!!> @author Iván Campos-Guereta Díez
	!!> @brief
	!!> Subroutine to construct the model from the inputs.
	!!> @param[in] ntype
	!!> @param[inout] mx
	!!> @param[in] option
	!!---------------------------------------------------------------------------------------------------------------------
	!
	!subroutine s_model_model_pass_sat_to_constraints_to_unsat(this)
	!
	!class(ty_model),target::this
	!
	!integer::iu,is,nunsat,nsat,i
	!
	!!Number of unsat and sat elements:
	!nunsat=this%mesh%vmod_count
	!nsat = this%layers%count-1
	!
	!!Fill constraint class with sat data:
	!do is=1,nsat
	!	call this%sat(is)%calc%get_results_nodes()
	!	call this%sat(is)%set_sat_to_mesh()
	!
	!this%constraints%hsat_h(is,:)		 = this%sat(is)%mesh%vmod_hsat(:)
	!this%constraints%dqhordx_h(is,:) = this%sat(is)%mesh%vmod_dqhordx(:)
	!	
	!end do
 !
	!!Set horizontal dqhordx_h to unsaturated mesh:
	!do iu=1,nunsat		
	!	do is=1,nsat
	!	this%unsat(iu)%mesh%dqhordx(is+1)=-this%constraints%dqhordx_h(is,iu)
	!	end do
	!end do
	!
	!end subroutine s_model_model_pass_sat_to_constraints_to_unsat	

	!!---------------------------------------------------------------------------------------------------------------------
	!! PASS_UNSAT_TO_CONSTRAINTS_TO_SAT
	!!---------------------------------------------------------------------------------------------------------------------
	!!> @author Iván Campos-Guereta Díez
	!!> @brief
	!!> Subroutine to construct the model from the inputs.
	!!> @param[in] ntype
	!!> @param[inout] mx
	!!> @param[in] option
	!!---------------------------------------------------------------------------------------------------------------------
	!
	!subroutine s_model_model_print_alltimes(this,findexconstraint)
	!use unsat_mod_ty_constraints,only:ty_unsat_constraints
	!use sat_mod_ty_constraints,only:ty_sat_constraints
	!
	!class(ty_model),target::this
	!integer,intent(in)::findexconstraint
	!class(ty_unsat_constraints),pointer::const_u
	!class(ty_sat_constraints),pointer::const_s
	!
	!integer::iu,is
	!
	!if(this%time%tprint == this%parameters%tinit) then
	!	write(findexconstraint,'(A)') 't,iteration,iu,is,v_idelem_hsat,v_hsat,v_qver,v_qver,v_qinf,v_incvoldt,v_qnewmann,v_qnewmann_all,nrel,h_width,h_hsat_mean, h_qent_mean,h_incvoldt_mean,h_dqhordx_mean,h_dqhordx_all_mean	'
	!end if	
	!
	!do iu=1,size(this%constraints%unsat)
	!	do is=1,size(this%constraints%sat)
	!		const_u => this%constraints%unsat(iu)%constraints
	!		const_s => this%constraints%sat(is)%constraints
	!		!														t						,iteration						,iu			,is			,v_idelem_hsat						,v_hsat							,v_qver						,v_qver						,v_qinf						,v_incvoldt						,v_qnewmann							,v_qnewmann_all						,nrel								,h_width							,h_hsat_mean						, h_qent_mean						,h_incvoldt_mean						,h_dqhordx_mean						,h_dqhordx_all_mean																	'
	!		write	(	findexconstraint, '(E15.6E3,","		,I10,","							,I2,","	,I2,","	,E15.6E3,","								,E15.6E3,","					,E15.6E3,","				,E15.6E3,","				,E15.6E3,","				,E15.6E3,","						,E15.6E3,","							,E15.6E3,","								,E15.6E3,","					,E15.6E3,","						,E15.6E3,","							,E15.6E3,","							,E15.6E3,","									,E15.6E3,","								,E15.6E3)')&
	!			!&													this%time%t	,this%time%iter_total	,iu			,is			,const_u%idelem_hsat(is)	,const_u%hsat(is)		,const_u%qver(is)	,const_u%qver(is)	,const_u%qinf(is)	,const_u%incvoldt(is)	,const_u%qnewmann(is)%p ,const_u%qnewmann_all(is)	,const_s%nrel(iu)%p	,const_s%width(iu)		,const_s%hsat_mean(iu)	,const_s%qent_mean(iu)	,const_s%incvoldt_mean(iu)	,const_s%dqhordx_mean(iu)	,const_s%dqhordx_all_mean(iu)	
	!			&													this%time%t	,this%time%iter_total	,iu			,is			,const_u%idelem_hsat(is)	,const_u%hsat(is)		,const_u%qver(is)	,const_u%qver(is)	,const_u%qinf(is)	,const_u%incvoldt(is)	,const_u%qnewmann(is)%p ,const_u%qnewmann_all(is)	,const_s%nrel(iu)%p	,const_s%width(iu)		,const_s%hsat_mean(iu)	,const_s%qent_mean(iu)	,const_s%incvoldt_mean(iu)	,const_s%dqhordx_mean(iu)	,const_s%dqhordx_all_mean(iu)	
 !       end do
	!end do
	!
	!
	!
	!end subroutine s_model_model_print_alltimes		
	
	!---------------------------------------------------------------------------------------------------------------------
	! PASS_UNSAT_TO_CONSTRAINTS_TO_SAT
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to construct the model from the inputs.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------
	
	subroutine s_model_model_print_alltimes(this,findexconstraint)
	use unsat_mod_ty_constraints,only:ty_unsat_constraints
	use sat_mod_ty_constraints,only:ty_sat_constraints
	
	class(ty_model),target::this
	integer,intent(in)::findexconstraint
	class(ty_unsat_constraints),pointer::const_u
	class(ty_sat_constraints),pointer::const_s
	
	integer::iu,is
	
	if(this%time%tprint == this%parameters%tinit) then
				write(findexconstraint,'(A)') 'idconstraint,t,iteration,iu,is,v_idelem_hsat,v_hnew,v_thnew,v_hsat,v_qvtb,v_qvt,v_qvb,v_incvoldt,v_Qs_layer,v_Qs_all,nrel,h_width,h_hsat_mean, h_qent_mean,h_incvoldt_mean,h_dqhordx_mean,h_dqhordx_all_mean,inchnew_mean	'
				!write(findexconstraint,'(A)') 'idconstraint,t,iteration,iu,is,v_idelem_hsat,v_hnew,v_thnew,v_hsat,v_qvtb,v_qvt,v_qvb,v_incvoldt,v_Qs_layer,v_Qs_all,nrel,h_width,h_hsat_mean, h_qent_mean,h_incvoldt_mean,h_dqhordx_mean,h_dqhordx_all_mean'
	end if	
	
	do iu=1,size(this%constraints%unsat)
		do is=1,size(this%constraints%sat)
			const_u => this%constraints%unsat(iu)%constraints
			const_s => this%constraints%sat(is)%constraints
			this%constraints%idprint = this%constraints%idprint+1
			!														idconstraint										t						,iteration						,iu			,is			,v_idelem_hsat					,v_hnew,						,v_thnew							,v_hsat							,v_qvtb						,v_qvt						,v_qvb						,v_incvoldt						,v_qnewmann							,v_qnewmann_all						,nrel								,h_width							,h_hsat_mean						, h_qent_mean						,h_incvoldt_mean						,h_dqhordx_mean						,h_dqhordx_all_mean				,h_dqhordx_all_mean				
			write	(	findexconstraint, '(I10,",",											E15.6E3,","		,I10,","							,I2,","	,I2,","	,E15.6E3,","				,E15.6E3,","						,E15.6E3,","					,E15.6E3,","				,E15.6E3,","			,E15.6E3,","			,E15.6E3,","			,E15.6E3,","					,E15.6E3,","						,E15.6E3,","							,E15.6E3,","				,E15.6E3,","					,E15.6E3,","						,E15.6E3,","						,E15.6E3,","								,E15.6E3,","							,E15.6E3,","									,E15.6E3)')&
				&													this%constraints%idprint,	this%time%t				,this%time%iter_total	,iu			,is			,const_u%idelem_hsat(is),const_u%hnew(is)%p	,const_u%thnew(is)%p	,const_u%hsat(is)		,const_u%qver(is)	,const_u%qsat(is)	,const_u%qinf(is)	,const_u%incvoldt(is)	,const_u%qnewmann(is)%P ,const_u%qnewmann_all(is)	,const_s%nrel(iu)%p	,const_s%width(iu)		,const_s%hsat_mean(iu)	,const_s%qent_mean(iu)	,const_s%incvoldt_mean(iu)	,const_s%dqhordx_mean(iu)	,const_s%dqhordx_all_mean(iu)	, const_s%inchnew_mean(iu)
			!write	(	findexconstraint, '(I10,",",											E15.6E3,","		,I10,","							,I2,","	,I2,","	,E15.6E3,","				,E15.6E3,","						,E15.6E3,","					,E15.6E3,","				,E15.6E3,","			,E15.6E3,","			,E15.6E3,","			,E15.6E3,","					,E15.6E3,","						,E15.6E3,","							,E15.6E3,","				,E15.6E3,","					,E15.6E3,","						,E15.6E3,","						,E15.6E3,","								,E15.6E3,","							,E15.6E3)')&
			!	&													this%constraints%idprint,	this%time%t				,this%time%iter_total	,iu			,is			,const_u%idelem_hsat(is),const_u%hnew(is)%p	,const_u%thnew(is)%p	,const_u%hsat(is)		,const_u%qver(is)	,const_u%qsat(is)	,const_u%qinf(is)	,const_u%incvoldt(is)	,const_u%qnewmann(is)%P ,const_u%qnewmann_all(is)	,const_s%nrel(iu)%p	,const_s%width(iu)		,const_s%hsat_mean(iu)	,const_s%qent_mean(iu)	,const_s%incvoldt_mean(iu)	,const_s%dqhordx_mean(iu)	,const_s%dqhordx_all_mean(iu)
		end do
	end do
	
	
	
	end subroutine s_model_model_print_alltimes		
	
	end module model_mod_ty_model