	!********************************************************************************************************************
	!        CLASS FOR THE COLLECTION OF CLASSES IN THE UNSATURATED MODEL
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

	module unsat_mod_ty_model
	use com_mod_ty_parameters,	only: ty_com_parameters
	use com_mod_ty_material,		only: ty_com_material
	use com_mod_ty_time,				only: ty_com_time
	use com_mod_ty_boundary,		only: ty_com_boundary
	use com_mod_ty_layers,			only: ty_com_layers

	use unsat_mod_ty_mesh,			only: ty_unsat_mesh
	use unsat_mod_ty_calc,			only: ty_unsat_calc
	use unsat_mod_ty_constraints, only: ty_unsat_constraints

	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_unsat_model
		type(ty_com_parameters)::parameters
		type(ty_com_material),allocatable::material(:)
		type(ty_unsat_mesh)::mesh
		type(ty_com_layers)::layers
		type(ty_com_boundary)::boundary
		type(ty_com_time)::time
		type(ty_unsat_constraints)::constraints

		type(ty_unsat_calc)::calc
	contains
	procedure,public:: readfileinput => s_unsat_model_readfileinput
	procedure,public:: construct=> s_unsat_model_construct
	procedure,public:: print_timestep=> s_unsat_model_print_timestep
	!procedure,public:: set_hsat_to_mesh => s_unsat_model_set_hsat_to_mesh
	!procedure,public:: set_satvalues_to_mesh => s_unsat_model_set_satvalues_to_mesh
	!procedure,public:: set_newmann_from_mesh_hsat => s_unsat_model_set_Newmann_from_mesh_hsat
	end type ty_unsat_model

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

	subroutine s_unsat_model_readfileinput(this,fileinput,fileboundary)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_model_readfileinput" :: s_unsat_model_readfileinput
	!DEC$ endif
	use unsat_mod_inputs, only: s_unsat_inputs_parameters,s_unsat_inputs_materials,s_unsat_inputs_layers,s_unsat_inputs_mesh,s_unsat_inputs_boundary

	class(ty_unsat_model),target::this
	character*400,intent(in)::fileinput
	character*400,intent(in)::fileboundary

	call s_unsat_inputs_parameters(this%parameters,fileinput)
	call s_unsat_inputs_materials(this%material,fileinput)
	call s_unsat_inputs_layers(this%layers,fileinput,this%material)
	call s_unsat_inputs_mesh(this%mesh,fileinput,this%layers)
	call s_unsat_inputs_boundary(fileboundary,this%boundary,this%layers)

	end subroutine s_unsat_model_readfileinput

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

	subroutine s_unsat_model_construct(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_model_construct" :: s_unsat_model_construct
	!DEC$ endif
	use unsat_mod_inputs, only: s_unsat_inputs_elements,s_unsat_inputs_nodes

	class(ty_unsat_model),target::this
	integer::i

	!Construct all instances in model from input file...
	this%calc%parameters=>this%parameters
	this%calc%time=>this%time
	this%calc%boundary=>this%boundary
	this%calc%layers => this%layers
	this%calc%material => this%material
	this%time%parameters => this%parameters

	call this%calc%allocateall()

	call s_unsat_inputs_nodes(this%calc%nodes,this%mesh,this%layers)
	call s_unsat_inputs_elements(this%calc%elements,this%mesh,this%calc%nodes,this%layers)

	call this%calc%construct()

	!Construct mesh:
	do i=1, size(this%mesh%nelemv) !CHECK: This doesnt consider multiple nodes in an element
		this%mesh%idnode_bottom(i) = sum(this%mesh%nelemv(1:i-1))+1
		this%mesh%idnode_top(i) = sum(this%mesh%nelemv(1:i))
	end do

	do i=1, size(this%mesh%nelemv) !CHECK: This doesnt consider multiple nodes in an element
		this%mesh%idnode_bottom(i) = sum(this%mesh%nelemv(1:i-1))+1
		this%mesh%idnode_top(i) = sum(this%mesh%nelemv(1:i))+1
	end do

	do i=1, size(this%mesh%nelemv)
		this%mesh%idelem_bottom(i) = sum(this%mesh%nelemv(1:i-1))+1
		this%mesh%idelem_top(i) = sum(this%mesh%nelemv(1:i))
	end do

	call this%constraints%allocateall(size(this%mesh%nelemv)-1,this%calc%nodes%count)

	!call this%constraints%construct(this%mesh,this%calc%nodes,this%mesh%dqhordx)
	call this%constraints%construct(this%mesh%idnode_bottom(2:),this%calc%nodes,this%mesh%dqhordx)

	end subroutine s_unsat_model_construct

	!---------------------------------------------------------------------------------------------------------------------
	! PRINT_TIMESTEP
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to print a time step in the output file.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------


	subroutine s_unsat_model_print_timestep(this,findexnode,findexelem,iu)
	!DEC$ if defined(_DLL)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_model_print_timestep" :: s_unsat_model_print_timestep
	!DEC$ endif
	!DEC$ endif
	use unsat_mod_ty_nodes,only:ty_unsat_nodes
	use unsat_mod_ty_elements,only:ty_unsat_elements
	use com_mod_ty_nodes,only:ty_com_nodes
	use com_mod_ty_elements,only:ty_com_elements

	class(ty_unsat_model),target::this
	class(ty_com_nodes),pointer::nodescom
	class(ty_unsat_nodes),pointer::nodes
	class(ty_com_elements),pointer::elementscom
	class(ty_unsat_elements),pointer::elements
	type(ty_unsat_nodes)::unsatnodes
	integer,intent(in)::findexnode,findexelem
	integer,intent(in)::iu

	integer::i

	nodescom=>this%calc%nodes
	elementscom=>this%calc%elements
	select type(nodescom)
	type is(ty_unsat_nodes)
		nodes=>nodescom
	end select

	select type(elementscom)
	type is(ty_unsat_elements)
		elements=>elementscom
	end select

	!Printing...
	if((this%calc%time%tprint == this%calc%parameters%tinit).and.(iu==1)) then
		write(findexnode,'(A)') 't,iu,x,z,h,th'
		write(findexelem,'(A)') 't,iu,nelem,x0,x1,hnew,hold,thnew,thold,qent,qsal,qmed,incvoldt,kmed,h0,h1,th0,th1,k0,k1,dhdx0,dhdx1,dhxdx0,dhxdx1,dhdxmed,dhxdxmed'
	end if


	!this%time%checkprint = .false.
	!this%time%tprint = min(this%time%tprint+this%calc%parameters%tprintinc,this%calc%parameters%tmax)
	! Here all the code that we want for when the time is iqual to print time

	unsatnodes = nodes
	!write(*,'("Printing... t= ",E10.3)') this%calc%time%t
	do i=1,nodes%count
		write(findexnode,'(E15.5E3,",",I2,",",E15.5E3,",",E15.5E3,",",E15.5E3,",",E15.5E3)', ADVANCE = 'YES') this%calc%time%t,iu, nodes%x(i), nodes%z(i),nodes%hnew(i),nodes%thnew(i)
	end do



	do i=1,elements%count
		!									'	t								,iu			,nelem	,x0																				,x1																			,hnew							,hold							,thnew							,thold							,qent											,qsal											,qmed											,incvoldt											,kmed											,h0							,h1							,th0						,th1							,k0							,k1							,dhdx0						,dhdx1						,dhxdx0							,dhxdx1							,dhdxmed										,dhxdxmed'
		write(findexelem,'(	E15.5,","				,I2,","	,I3,","	,E15.5,","																,E15.5,","															,E15.5,","				,E15.5,","				,E15.5,","					,E15.5,","					,E15.5,","								,E15.5,","								,E15.5,","								,E15.5,","										,E15.5,","								,E15.5,","			,E15.5,","			,E15.5,","			,E15.5,","				,E15.5,","			,E15.5,","			,E15.5,","				,E15.5,","				,E15.5,","					,E15.5,","					,E15.5,","									,E15.5)', ADVANCE = 'YES')&
			&									this%calc%time%t,iu			,i			,this%calc%nodes%x(elements%idnode(i,1))	,this%calc%nodes%x(elements%idnode(i,2)),elements%hnew(i)	,elements%hold(i), elements%thnew(i)	,elements%thold(i)	,elements%results_qent(i)	,elements%results_qsal(i)	,elements%results_qmed(i)	,elements%results_incvoldtperm(i)	,elements%results_kmed(i)	,elements%h0(i)	,elements%h1(i)	,elements%th0(i),elements%th1(i)	,elements%k0(i)	,elements%k1(i)	,elements%dhdx0(i),elements%dhdx1(i),elements%dhxdx0(i)	,elements%dhxdx1(i)	,elements%results_dhdxmed(i),elements%results_dhxdxmed(i)
	end do

	end subroutine s_unsat_model_print_timestep


	!!******************************************************************************************************************
	!! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!!******************************************************************************************************************
	!
	!subroutine s_unsat_model_set_hsat_to_mesh(this,option)
	!use unsat_mod_ty_layers,only:ty_unsat_layers
	!
	!class(ty_unsat_model),intent(inout)::this
	!integer,intent(in),optional::option
	!
	!integer::opt
	!integer::i,l,idnode
	!real(kind=dpd)::htemp,hsattemp,hsatmax
	!
	!opt = 0
	!
	!if(present(option)) opt = option
	!
	!do l=1,this%layers%count
	!idnode = sum(this%mesh%nelemv(1:l-1))+1
	!!get hsat from the bottom node of the layer
	!select case(opt)
	!case (1)
	!	htemp = this%calc%nodes%hold(idnode)
	!case (2)
	!	htemp = this%calc%nodes%htemp(idnode)
	!case default
	!	htemp = this%calc%nodes%hnew(idnode)
	!end select
	!
	!if (htemp > 0.0_dpd) then
	!	this%mesh%hsat(l) = htemp
	!	this%mesh%idnodsat(l) = idnode
	!	this%mesh%idnodlayer(l) = idnode
	!else
	!	this%mesh%hsat(l) = 0.0_dpd
	!	this%mesh%idnodsat(l) = 0
	!	this%mesh%idnodlayer(l) = idnode
	! end if
	!
	!this%mesh%hsatold(l) = max(0.0_dpd,this%calc%nodes%hold(idnode))
	!	end do
	!
	!!this is another method to get hsat (from the top node)
	!!this%hsat = -0.0_dpd
	!!hsatmax = -0.0
	!!this%nnodsat = 0
	!!
	!!do i=this%nnodes,1,-1
	!!	select case(opt)
	!!	case (1)
	!!		htemp = this%nodeptr(i)%node%hold
	!!	case (2)
	!!		htemp = this%nodeptr(i)%node%htemp
	!!		case default
	!!		htemp = this%nodeptr(i)%node%hnew
	!!	end select
	!!		!hsatmax is going to be the max of hsats on nodes below
	!!		htemp = max(hsatmax-(this%nodeptr(i)%node%x-this%hbottom),htemp)
	!!
	!!	if (htemp > 1.0e-8) then
	!!		this%nnodsat = this%nnodsat + 1
	!!		hsatmax = htemp+(this%nodeptr(i)%node%x-this%hbottom)
	!!	else
	!!		exit
	!!	end if
	!!end do
	!!
	!!this%hsat = hsatmax
	!
	!end subroutine s_unsat_model_set_hsat_to_mesh





	!!******************************************************************************************************************
	!! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!!******************************************************************************************************************
	!
	!subroutine s_unsat_model_set_satvalues_to_mesh(this,option)
	!use unsat_mod_ty_layers,only:ty_unsat_layers
	!use com_mod_linear_interpolation,only:linearinterpolation_ordered
	!use com_mod_fem_intelement, only:intelement_abs_1don_sca
	!
	!integer,parameter::INTEGRATION_ORDER=10
	!class(ty_unsat_model),intent(inout)::this
	!integer,intent(in),optional::option
	!real(kind=dpd)::xsat(2)
	!
	!integer::opt
	!integer::i,l,idnode
	!real(kind=dpd)::htemp,hsattemp,hsatmax
	!
	!do l=1,this%layers%count
	!	if (this%mesh%hsat(l).ne.0.0_dpd) then
	!		xsat(2) = this%calc%nodes%x(this%mesh%idnodlayer(l))+this%mesh%hsat(l) !This is x for hsatnew
	!		xsat(1) = this%calc%nodes%x(this%mesh%idnodlayer(l))+this%mesh%hsatold(l)	!!This is x for hsatold
	!     if (xsat(2)==xsat(1)) then
	!		this%mesh%incvsat(l)=0.0_dpd
	!		this%mesh%qvsat(l)	= this%mesh%incvsat(l)/this%calc%time%dt+this%mesh%dqhordx(l)	!CHECK for when more than 1 layer. Vertical waterflow that will enter to saturated model due to the increment of water from hsatold to hsat
	!     !this%mesh%nmean(l) = this%mesh%nmean(l)
	!     !this%mesh%nrel(l) = this%mesh%nrel(l)
	!     else
	!     this%mesh%incvsat(l)=intelement_abs_1don_sca(f_th_in_x,xsat,INTEGRATION_ORDER)			!Increment of water content from xsatold to xsatnew
	!     this%mesh%qvsat(l)	= this%mesh%incvsat(l)/this%calc%time%dt+this%mesh%dqhordx(l)	!CHECK for when more than 1 layer. Vertical waterflow that will enter to saturated model due to the increment of water from hsatold to hsat
	!     this%mesh%nmean(l)	= abs(this%mesh%incvsat(l)/(xsat(2)-xsat(1)))										!Mean non-wetting porosity available between xsatold and xsatnew (to pass it to sat model in this layer)
	!		this%mesh%nrel(l)		= this%mesh%nmean(l)/(this%calc%nodes%material(this%calc%get_idelement_from_x(xsat(2)))%thsat-this%calc%nodes%material(this%calc%get_idelement_from_x(xsat(2)))%thres) !CHECK because n maybe has to be integrated from hsat to hsatold
	!     end if
	!
	!	else
	!		this%mesh%incvsat(l)=0.0_dpd	!CHECK for the case when h goes from hsatold to 0.0
	!		this%mesh%qvsat(l)=0.0_dpd		!In this case no vertical flow is passed to horizontal model as there is no saturated layer.
	!		this%mesh%nmean(l)=this%layers%material(l)%thsat-this%layers%material(l)%thres
	!		this%mesh%nrel(l)=1.0_dpd
	!	end if
	!end do
	!
	!contains
	!function f_th_in_x(x)
	!use com_mod_linear_interpolation,only:linearinterpolation_ordered
	!use com_mod_hyd,only:f_hyd_th_h_sca
	!
	!real(kind=dpd),intent(in)::x
	!real(kind=dpd)::f_th_in_x
	!!CHECK: Check this because we are integrating with the material in node, but the equations are being done with material in element.
	!f_th_in_x =		f_hyd_th_h_sca(linearinterpolation_ordered(x,this%calc%nodes%x,this%calc%nodes%hnew),this%calc%nodes%material(this%calc%get_idelement_from_x(x)))-&
	!						&	f_hyd_th_h_sca(linearinterpolation_ordered(x,this%calc%nodes%x,this%calc%nodes%hold),this%calc%nodes%material(this%calc%get_idelement_from_x(x)))
	!
	!end function f_th_in_x
	!
	!end subroutine s_unsat_model_set_satvalues_to_mesh

	!!******************************************************************************************************************
	!! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!!******************************************************************************************************************
	!
	!subroutine s_unsat_model_set_Newmann_from_mesh_hsat(this)
	!class(ty_unsat_model),intent(inout)::this
	!
	!integer::l
	!
	!do l=2,this%layers%count
	!if (this%mesh%idnodsat(l)>0) then
	!	this%calc%nodes%qnewmann(this%mesh%idnodsat(l)) = this%mesh%dqhordx(l)
	!	this%calc%nodes%isnewmann(this%mesh%idnodsat(l)) = .true.
	!else
	!	this%calc%nodes%qnewmann(this%mesh%idnodlayer(l)) = 0.0_dpd
	!	this%calc%nodes%isnewmann(this%mesh%idnodlayer(l)) = .false.
	!end if
	!end do
	!
	!end subroutine s_unsat_model_set_Newmann_from_mesh_hsat

	end module unsat_mod_ty_model