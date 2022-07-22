	!********************************************************************************************************************
	! TITLE         : SAT_MOD_TY_MESH: DERIVED TYPE TO DEFINE POROPERTIES AND METHODS OF THE WF1DSAT MESH
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_CALC
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Derived type to define poroperties and methods of the wf1dsat mesh
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module sat_mod_ty_mesh

	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_sat_mesh
		integer::nn		!<number of nodes per element
		integer::nc		!<number of clases of nodes (0=c0, 1=c1,...)
		integer::nnc	!<number of clases and nodes in elements (0=c0, 1=c1,...)
		integer::nelemh_count	!<number of elements
		integer,allocatable					::nelemh(:) !<Number of elements between vertical modules (first always at x=0)
		integer::nnodh_count	!<number of nodes
		integer::nnodclassh_count	!<number of nodes and classes
		real(kind=dps)::width	!<total width of the model


		integer::vmod_count !<number of vertical modules
		real(kind=dps),allocatable	::vmod_x(:) !<x coord of every vertical module
		integer,allocatable					::vmod_idnod(:) !<node number of every vertical module (before numnod)

		real(kind=dps),allocatable	::vmod_qent(:) !especific flow in each vertical module node

		real(kind=dps),allocatable	::vmod_nrel(:) !Relarive porosity of non-wetting voids [int(th,hold,hnew)/(|Hnew-hold|·(thsat-thres))]

	contains
	procedure,public:: allocateall		=> s_sat_mesh_allocateall
	procedure,public:: deallocateall	=> s_sat_mesh_deallocateall
	procedure,public:: set_z					=> s_sat_mesh_set_z

	end type ty_sat_mesh

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_sat_mesh_allocateall(this,nvmod)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_mesh_allocateall
	!DEC$ endif
	class(ty_sat_mesh),intent(inout)::this
	integer,intent(in)::nvmod
	integer::nnodes

	nnodes = this%nnodclassh_count
	this%vmod_count = nvmod

	if(.not.allocated(this% nelemh))		allocate(this% nelemh (nvmod))

	if(.not.allocated(this% vmod_x))				allocate(this% vmod_x (nvmod))
	if(.not.allocated(this% vmod_idnod))		allocate(this% vmod_idnod (nvmod))

	if(.not.allocated(this% vmod_qent))			allocate(this% vmod_qent (nvmod))

	if(.not.allocated(this% vmod_nrel))			allocate(this% vmod_nrel (nvmod))

	this% vmod_nrel = 1.0_dpd
	end subroutine s_sat_mesh_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate all arrays
	!---------------------------------------------------------------------------------------------------------------------


	subroutine s_sat_mesh_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_mesh_deallocateall
	!DEC$ endif
	class(ty_sat_mesh),intent(inout)::this

	if(allocated(this% nelemh))		deallocate(this% nelemh )

	if(allocated(this% vmod_x))				deallocate(this% vmod_x )
	if(allocated(this% vmod_idnod))		deallocate(this% vmod_idnod )

	if(allocated(this% vmod_qent))		deallocate(this% vmod_qent)

	if(allocated(this% vmod_nrel))			deallocate(this% vmod_nrel)

	end subroutine s_sat_mesh_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Set coord z at each node
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_mesh_set_z(this,nodes,layers)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_mesh_set_z
	!DEC$ endif
	use sat_mod_ty_nodes, only:ty_sat_nodes
	use com_mod_ty_layers, only:ty_com_layers

	class(ty_sat_mesh),intent(in)::this
	class(ty_sat_nodes),intent(inout)::nodes
	class(ty_com_layers),intent(in)::layers

	nodes%z = (this%width-nodes%x)*layers%slope(1)

	end subroutine	s_sat_mesh_set_z

	end module sat_mod_ty_mesh