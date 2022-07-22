	!********************************************************************************************************************
	! TITLE         : UNSAT_MOD_TY_MODEL: DERIVED TYPE TO DEFINE PROPERTIES AND METHODS SPECIFIC TO THE MODEL WF1DUNSAT
	! PROJECT       : WF1DUNSATDLL
	! MODULE        : UNSAT_MOD_TY_MODEL
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Derived type to define properties and methods specific to the model wf1dunsat
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
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
	end type ty_unsat_model

	contains

	!---------------------------------------------------------------------------------------------------------------------
	! READFILEINPUT
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to construct the model from the inputs.
	!> @param[in] fileinput
	!> @param[in] fileboundary
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
	!> Assign asignable properties to the class instance.
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
	do i=1, size(this%mesh%nelemv)
		this%mesh%idnode_bottom(i) = sum(this%mesh%nelemv(1:i-1))+1
		this%mesh%idnode_top(i) = sum(this%mesh%nelemv(1:i))
	end do

	do i=1, size(this%mesh%nelemv)
		this%mesh%idnode_bottom(i) = sum(this%mesh%nelemv(1:i-1))+1
		this%mesh%idnode_top(i) = sum(this%mesh%nelemv(1:i))+1
	end do

	do i=1, size(this%mesh%nelemv)
		this%mesh%idelem_bottom(i) = sum(this%mesh%nelemv(1:i-1))+1
		this%mesh%idelem_top(i) = sum(this%mesh%nelemv(1:i))
	end do

	call this%constraints%allocateall(size(this%mesh%nelemv)-1,this%calc%nodes%count)

	call this%constraints%construct(this%mesh%idnode_bottom(2:),this%calc%nodes,this%mesh%dqhordx)

	end subroutine s_unsat_model_construct

	!---------------------------------------------------------------------------------------------------------------------
	! PRINT_TIMESTEP
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to print a time step in the output file.
	!> @param[in] findexnode
	!> @param[in] findexelem
	!> @param[in] iu
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

	unsatnodes = nodes

	do i=1,nodes%count
		write(findexnode,'(E15.5E3,",",I2,",",E15.5E3,",",E15.5E3,",",E15.5E3,",",E15.5E3)', ADVANCE = 'YES') this%calc%time%t,iu, nodes%x(i), nodes%z(i),nodes%hnew(i),nodes%thnew(i)
	end do



	do i=1,elements%count
		!									'	t								,iu			,nelem	,x0																				,x1																			,hnew							,hold							,thnew							,thold							,qent											,qsal											,qmed											,incvoldt											,kmed											,h0							,h1							,th0						,th1							,k0							,k1							,dhdx0						,dhdx1						,dhxdx0							,dhxdx1							,dhdxmed										,dhxdxmed'
		write(findexelem,'(	E15.5,","				,I2,","	,I3,","	,E15.5,","																,E15.5,","															,E15.5,","				,E15.5,","				,E15.5,","					,E15.5,","					,E15.5,","								,E15.5,","								,E15.5,","								,E15.5,","										,E15.5,","								,E15.5,","			,E15.5,","			,E15.5,","			,E15.5,","				,E15.5,","			,E15.5,","			,E15.5,","				,E15.5,","				,E15.5,","					,E15.5,","					,E15.5,","									,E15.5)', ADVANCE = 'YES')&
			&									this%calc%time%t,iu			,i			,this%calc%nodes%x(elements%idnode(i,1))	,this%calc%nodes%x(elements%idnode(i,2)),elements%hnew(i)	,elements%hold(i), elements%thnew(i)	,elements%thold(i)	,elements%results_qent(i)	,elements%results_qsal(i)	,elements%results_qmed(i)	,elements%results_incvoldtperm(i)	,elements%results_kmed(i)	,elements%h0(i)	,elements%h1(i)	,elements%th0(i),elements%th1(i)	,elements%k0(i)	,elements%k1(i)	,elements%dhdx0(i),elements%dhdx1(i),elements%dhxdx0(i)	,elements%dhxdx1(i)	,elements%results_dhdxmed(i),elements%results_dhxdxmed(i)
	end do

	end subroutine s_unsat_model_print_timestep

	end module unsat_mod_ty_model