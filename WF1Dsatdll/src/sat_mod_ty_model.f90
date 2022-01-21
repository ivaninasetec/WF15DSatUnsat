	!********************************************************************************************************************
	!        CLASS FOR THE COLLECTION OF NODES-CLASSES IN THE SATURATED MODEL
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
	!> Iv�n Campos-Guereta D�ez
	!  MSc Civil Engineering by Polytechnic University of Madrid                                                     *
	!  PhD Student by University of Nottingham                                                                       *
	!  eMBA by International Institute San Telmo in Seville                                                          *
	!  ivan.camposguereta@nottingham.ac.uk
	! DESCRIPTION:
	!> Class for the collection of nodes-classes in the saturated model
	!********************************************************************************************************************

	module sat_mod_ty_model
	use com_mod_ty_parameters,only: ty_com_parameters
	use com_mod_ty_material,	only: ty_com_material
	use sat_mod_ty_mesh,			only: ty_sat_mesh
	use com_mod_ty_layers,		only: ty_com_layers
	use com_mod_ty_boundary,	only: ty_com_boundary
	use com_mod_ty_time,			only: ty_com_time
	use sat_mod_ty_calc,			only: ty_sat_calc
	use sat_mod_ty_constraints,only:ty_sat_constraints


	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_sat_model
		type(ty_com_parameters)::						parameters
		type(ty_com_material),allocatable::	material(:)
		type(ty_sat_mesh)::									mesh
		type(ty_com_layers)::								layers
		type(ty_com_time)::									time
		type(ty_com_boundary)::							boundary
		type(ty_sat_calc)::									calc
		type(ty_sat_constraints)::					constraints
	contains
	procedure,public:: readfileinput		=> s_sat_model_readfileinput
	procedure,public:: construct				=> s_sat_model_construct
	procedure,public:: print_timestep		=> s_sat_model_print_timestep
	!procedure,public:: set_sat_to_mesh	=> s_sat_model_set_hsat_to_mesh
	!procedure,public:: set_qent_from_mesh_to_nodes	=> s_sat_model_put_qent_from_mesh
	!procedure,public:: set_nrel_from_mesh_to_nodes	=> s_sat_model_put_nrel_from_mesh
	procedure,public:: set_qent_from_constraint_to_nodes	=> s_sat_model_put_qent_from_constraint
	procedure,public:: set_nrel_from_constraint_to_nodes	=> s_sat_model_put_nrel_from_constraint
	procedure,public:: put_results_in_constraints	=> s_sat_model_put_results_in_constraints
	end type ty_sat_model

	contains

	subroutine s_sat_model_readfileinput(this,fileinput,fileboundary)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_model_readfileinput
	!DEC$ endif
	use sat_mod_inputs, only: s_sat_inputs_parameters,s_sat_inputs_materials,s_sat_inputs_layers,s_sat_inputs_mesh,s_sat_inputs_boundary

	class(ty_sat_model)::this
	character*400,intent(in)::fileinput
	character*400,intent(in)::fileboundary

	call s_sat_inputs_parameters(this%parameters,fileinput)
	call s_sat_inputs_materials(this%material,fileinput)
	call s_sat_inputs_layers(this%layers,fileinput,this%material)
	call s_sat_inputs_mesh(this%mesh,fileinput,this%layers)
	call s_sat_inputs_boundary(fileboundary,this%boundary,this%layers)

	end subroutine s_sat_model_readfileinput


	subroutine s_sat_model_construct(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_model_construct
	!DEC$ endif
	use sat_mod_inputs, only: s_sat_inputs_elements,s_sat_inputs_nodes

	class(ty_sat_model),target::this

	!Construct all instances in model from input file...
	this%calc%parameters => this%parameters
	this%calc%time => this%time
	this%calc%layers => this%layers
	this%calc%boundary => this%boundary
	this%time%parameters => this%parameters
	call this%calc%allocateall()

	call s_sat_inputs_nodes(this%calc%nodes,this%mesh,this%layers)
	call s_sat_inputs_elements(this%calc%elements,this%mesh,this%calc%nodes,this%layers%material(1))

	call this%calc%construct()

	call this%constraints%allocateall(size(this%mesh%vmod_idnod),this%mesh%nnodclassh_count)
	call this%constraints%construct(this%mesh%vmod_idnod,this%calc%nodes,this%mesh%vmod_qent,this%mesh%vmod_nrel)

	end subroutine s_sat_model_construct

	subroutine s_sat_model_print_timestep(this,findexnode,findexelem,is)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_model_print_timestep
	!DEC$ endif
	use sat_mod_ty_nodes,only:ty_sat_nodes
	use sat_mod_ty_elements,only:ty_sat_elements
	use com_mod_ty_nodes,only:ty_com_nodes
	use com_mod_ty_elements,only:ty_com_elements

	class(ty_sat_model),intent(in),target::this
	integer,intent(in)::findexnode
	integer,intent(in)::findexelem
	integer,intent(in)::is
	class(ty_com_nodes),pointer::nodescom
	class(ty_sat_nodes),pointer::nodes
	class(ty_com_elements),pointer::elementscom
	class(ty_sat_elements),pointer::elements
	type(ty_sat_nodes)::satnodes

	integer::i

	nodescom=>this%calc%nodes
	select type(nodescom)
	type is(ty_sat_nodes)
		nodes =>nodescom
	end select

	elementscom=>this%calc%elements
	select type(elementscom)
	type is(ty_sat_elements)
		elements => elementscom
	end select

	!Check if we are in printstep
	if((this%calc%time%tprint == this%calc%parameters%tinit).and.(is==1)) then
		write(findexnode,'(A)') 't,is,x,z,head,qent,qincvol,qhor,dqhordx,dqhordx_from_incvol,qhor_all,dqhordx_all,dqhordx_from_incvol_all'
		write(findexelem,'(A)') 't,is,ie,x0,x1,z0,z1,h0,h1,h,dhdx,dhzdx,qent,incvoldt,dqhordx,dqhordx_all,dqhordx_from_incvoldt,dqhordx_from_incvoldt_all,q,q_all'
	end if


	!this%time%checkprint = .false.
	!this%time%tprint = min(this%time%tprint+this%calc%parameters%tprintinc,this%calc%parameters%tmax)
	! Here all the code that we want for when the time is iqual to print time

	satnodes = nodes
	!write(*,'("Printing... t= ",E10.3)') this%calc%time%t
	do i=1,this%calc%nodes%count
		!										t									is				x						z							head						qent										qincvol											qhor										dqhordx											dqhordx_from_incvol												qhor_all										dqhordx_all											dqhordx_from_incvol_all
		write(findexnode,'( E15.6E3,","				,I2,","		,E15.6E3,",", E15.6E3,","	,E15.6E3,","		, E15.6E3,","						,E15.6E3,","								,	E15.6E3,","						, E15.6E3,","								, E15.6E3,","															, E15.6E3,","								, E15.6E3,","									  , E15.6E3										)', ADVANCE = 'YES')&
			&									this%calc%time%t, is				, nodes%x(i), nodes%z(i)	,nodes%hnew(i)	, nodes%results_qent(i)	,nodes%results_incvoldt(i)	, nodes%results_qhor(i)	, nodes%results_dqhordx(i)	, nodes%results_dqhordx_from_incvoldt(i)	, nodes%results_qhor_all(i)	, nodes%results_dqhordx_all(i)	, nodes%results_dqhordx_from_incvoldt_all(i)
	end do


	do i=1,this%calc%elements%count
		!										t 								,is			,ie			,x0							,x1							,z0							,z1							,h0							,h1							,h								,dhdx							,dhzdx							,qent											,incvoldt											,dqhordx											,dqhordx_all											,dqhordx_from_incvoldt											,dqhordx_from_incvoldt_all										,q							,q_all'
		write(findexelem,'(	E15.6,","					,I2,","	,I3,","	,E15.6,","			,E15.6,","			,E15.6,","			,E15.6,","			,E15.6,","			,E15.6,","			,E15.6,","				,E15.6,","				,E15.6,","					,E15.6,","								,E15.6,","										,E15.6,","										,E15.6,","												,E15.6,","																	,E15.6,","																		,E15.6,","			,E15.6											)', ADVANCE = 'YES') &
			&									this%calc%time%t	,is			,i			,elements%x0(i)	,elements%x1(i)	,elements%z0(i)	,elements%z1(i)	,elements%h0(i)	,elements%h1(i)	,elements%hnew(i)	,elements%dhdx(i)	,elements%dhzdx(i)	,elements%results_qent(i)	,elements%results_incvoldt(i)	,elements%results_dqhordx(i)	,elements%results_dqhordx_all(i)	,elements%results_dqhordx_from_incvoldt(i)	,elements%results_dqhordx_from_incvoldt_all(i),	elements%q(i)	,elements%q_all(i)
	end do


	end subroutine s_sat_model_print_timestep


	!!******************************************************************************************************************
	!! Sub: s_sat_model_set_hsat_to_mesh(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!!******************************************************************************************************************
	!
	!subroutine s_sat_model_set_hsat_to_mesh(this)
	!use com_mod_ty_nodes, only: ty_com_nodes
	!use sat_mod_ty_nodes, only: ty_sat_nodes
	!class(ty_sat_model),intent(inout),target::this
	!class(ty_com_nodes),pointer::nodes
	!
	!integer::opt
	!integer::iv
	!real(kind=dpd)::htemp,hsattemp,hsatmax
	!
	!nodes => this%calc%nodes
	!
	!!opt = 0
	!!
	!!if(present(option)) opt = option
	!!
	!do iv=1,this%mesh%vmod_count
	!!!get hsat from the bottom node of the layer
	!!select case(opt)
	!!case (1)
	!!	htemp = this%calc%nodes%hold	(this%mesh%vmod_idnode(iv))
	!!case (2)
	!!	htemp = this%calc%nodes%htemp	(this%mesh%vmod_idnode(iv))
	!!case default
	!!	htemp = this%calc%nodes%hnew	(this%mesh%vmod_idnode(iv))
	!!end select
	!
	!this%mesh%vmod_hsat(iv)			= this%calc%nodes%hnew(this%mesh%vmod_idnod(iv))
	!this%mesh%vmod_hsattemp(iv) = this%calc%nodes%htemp(this%mesh%vmod_idnod(iv))
	!this%mesh%vmod_hsatold(iv)	= this%calc%nodes%hold(this%mesh%vmod_idnod(iv))
	!
	!select type(nodes)
	!	type is(ty_sat_nodes)
	!	this%mesh%vmod_dqhordx(iv)	= nodes%results_dqhordx(this%mesh%vmod_idnod(iv))
	!end select
	!
	!end do
	!
	!end subroutine s_sat_model_set_hsat_to_mesh

	!!******************************************************************************************************************
	!! Sub: s_sat_model_set_hsat_to_mesh(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!!******************************************************************************************************************
	!
	!subroutine s_sat_model_put_qent_from_mesh(this)
	!
	!class(ty_sat_model),intent(inout),target::this
	!
	!this%calc%nodes%qent = matmul(this%mesh%vmod_intep_matrix,this%mesh%vmod_qent)
	!
	!end subroutine s_sat_model_put_qent_from_mesh

	!******************************************************************************************************************
	! Sub: s_sat_model_set_hsat_to_mesh(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_sat_model_put_qent_from_constraint(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_model_put_qent_from_constraint
	!DEC$ endif

	class(ty_sat_model),intent(inout),target::this
	real(kind=dpd)::qent(size(this%constraints%idnode))

	integer::i

	do i=1,size(this%constraints%idnode)
		qent(i) = this%constraints%qent(i)%p
	end do

	!!CHECK: I am doing the mean with nodes from 2 to end because the first unsaturated node gives unstabilities and to smooth the whole model:
	!qent = sum(qent(2:size(this%constraints%idnode)))/real(size(this%constraints%idnode),dpd)

	this%calc%nodes%qent = matmul(this%constraints%intep_matrix,qent)

	end subroutine s_sat_model_put_qent_from_constraint

	!******************************************************************************************************************
	! Sub: s_sat_model_set_hsat_to_mesh(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_sat_model_put_nrel_from_constraint(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_model_put_nrel_from_constraint
	!DEC$ endif

	class(ty_sat_model),intent(inout),target::this
	real(kind=dpd)::nrel(size(this%constraints%idnode))

	integer::i

	do i=1,size(this%constraints%idnode)
		nrel(i) = min(1.0_dpd,this%constraints%nrel(i)%p) !CHECK: I put this because I was getting nrel>1
	end do
	!!CHECK: I am doing the mean with nodes from 2 to end because the first unsaturated node gives unstabilities and to smooth the whole model:
	!nrel = sum(nrel(2:size(this%constraints%idnode)))/real(size(this%constraints%idnode),dpd)

	this%calc%nrel = matmul(this%constraints%intep_matrix,nrel)

	end subroutine s_sat_model_put_nrel_from_constraint

	!!******************************************************************************************************************
	!! Sub: s_sat_model_set_hsat_to_mesh(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!!******************************************************************************************************************
	!
	!subroutine s_sat_model_put_nrel_from_mesh(this)
	!
	!class(ty_sat_model),intent(inout),target::this
	!
	!this%calc%nrel = matmul(this%mesh%vmod_intep_matrix,this%mesh%vmod_nrel)
	!
	!end subroutine s_sat_model_put_nrel_from_mesh

	!--------------------------------------------------------------------------------------------------------------------
	! F: s_sat_model_put_results_in_constraints(is,iu)
	!--------------------------------------------------------------------------------------------------------------------
	!> @brief
	!> Put results from elements in constraints from results in elements.
	!> In this case results are taken from the mean values of the segments, values updated in this%constraints% are:
	!> - hsat_mean: mean elems%hnew in interval.
	!> - qent_mean: mean elems%results_qent in interval.
	!> - incvoldt_mean: mean qent in interval.
	!> @param[inout] this This class.
	!> @param[in] is Layer where to measure hsat.
	!> @param[in] iu UNSAT model where to measure hsat.
	!> @return	array nxxnu with the values of dqhordx in the SAT models at each (is,iu).
	!> @author Iv�n Campos-Guereta D�ez

	subroutine s_sat_model_put_results_in_constraints(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_model_put_results_in_constraints
	!DEC$ endif
	use com_mod_ty_elements, only: ty_com_elements
	use sat_mod_ty_elements, only: ty_sat_elements

	class(ty_sat_model),intent(inout),target::this
	class(ty_com_elements),pointer::elem
	class(ty_sat_elements),pointer::elems

	integer::i,idelem0,idelem1

	elem=>this%calc%elements

	select type(elem)
	type is(ty_sat_elements)
		elems=>elem
	end select

	do i=1,size(this%constraints%idnode)
		idelem0 = this%constraints%idelem0(i)
		idelem1 = this%constraints%idelem1(i)
		this%constraints%hsat_mean(i)= dot_product(elems%hnew(idelem0:idelem1),elems%lenght(idelem0:idelem1))/this%constraints%width(i)
		this%constraints%qent_mean(i)= dot_product(elems%results_qent(idelem0:idelem1),elems%lenght(idelem0:idelem1))/this%constraints%width(i)
		this%constraints%incvoldt_mean(i)= dot_product(elems%results_incvoldt(idelem0:idelem1),elems%lenght(idelem0:idelem1))/this%constraints%width(i)
		this%constraints%dqhordx_mean(i)= dot_product(elems%results_dqhordx_from_incvoldt(idelem0:idelem1),elems%lenght(idelem0:idelem1))/this%constraints%width(i)
		this%constraints%dqhordx_all_mean(i)= dot_product(elems%results_dqhordx_from_incvoldt_all(idelem0:idelem1),elems%lenght(idelem0:idelem1))/this%constraints%width(i)
	end do
	idelem0 = this%constraints%idelem0(1)
	idelem1 = this%constraints%idelem1(1)
	this%constraints%qhor0(1) = elems%q0(idelem0)
	this%constraints%qhor1(1) = (elems%q1(idelem1-1)+elems%q0(idelem1))/2.0_dpd
	do i=2,size(this%constraints%idnode)-1
		idelem0 = this%constraints%idelem0(i)
		idelem1 = this%constraints%idelem1(i)
		this%constraints%qhor0(i) = (elems%q1(idelem0-1)+elems%q0(idelem0))/2.0_dpd
		this%constraints%qhor1(i) = (elems%q1(idelem1-1)+elems%q0(idelem1))/2.0_dpd
	end do
	idelem0 = this%constraints%idelem0(size(this%constraints%idnode))
	idelem1 = this%constraints%idelem1(size(this%constraints%idnode))
	this%constraints%qhor0(size(this%constraints%idnode)) = (elems%q1(idelem0-1)+elems%q0(idelem0))/2.0_dpd
	this%constraints%qhor1(size(this%constraints%idnode)) = elems%q1(idelem1)

	end subroutine s_sat_model_put_results_in_constraints

	end module sat_mod_ty_model