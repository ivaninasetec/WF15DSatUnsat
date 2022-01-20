	!********************************************************************************************************************
	!        MODULE: CLASS WITH MESH CHARACTERISTICS OF SATURATED MODEL
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : sat_mod_ty_mesh
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
	!> Class of collection of layers in horizontal saturated model
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
		!
		!  real(kind=dps),allocatable	::vmod_hsat(:)
		!real(kind=dps),allocatable	::vmod_hsattemp(:)
		!real(kind=dps),allocatable	::vmod_hsatold(:)
		!
		real(kind=dps),allocatable	::vmod_qent(:) !especific flow in each vertical module node
		!real(kind=dps),allocatable	::vmod_qenttemp(:) !especific flow in each vertical module node
		!real(kind=dps),allocatable	::vmod_qentold(:) !especific flow in each vertical module node
		!
		!real(kind=dps),allocatable	::vmod_dqhordx(:) !especific flow in each vertical module node
		real(kind=dps),allocatable	::vmod_nrel(:) !Relarive porosity of non-wetting voids [int(th,hold,hnew)/(|Hnew-hold|·(thsat-thres))]
		!real(kind=dps),allocatable	::vmod_intep_matrix(:,:) !Interpolation matrix to get values in nodes from values in vertical modules.

	contains
	procedure,public:: allocateall		=> s_sat_mesh_allocateall
	procedure,public:: deallocateall	=> s_sat_mesh_deallocateall
	procedure,public:: set_z					=> s_sat_mesh_set_z
	!procedure,public:: construct_interpmatrix	=> s_sat_mesh_construct_interpmatrix

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
	!
	!if(.not.allocated(this% vmod_hsat))			allocate(this% vmod_hsat (nvmod))
	!if(.not.allocated(this% vmod_hsattemp)) allocate(this% vmod_hsattemp (nvmod))
	!if(.not.allocated(this% vmod_hsatold))	allocate(this% vmod_hsatold (nvmod))
	if(.not.allocated(this% vmod_qent))			allocate(this% vmod_qent (nvmod))
	!if(.not.allocated(this% vmod_qenttemp)) allocate(this% vmod_qenttemp (nvmod))
	!if(.not.allocated(this% vmod_qentold))	allocate(this% vmod_qentold (nvmod))
	!if(.not.allocated(this% vmod_dqhordx))	allocate(this% vmod_dqhordx (nvmod))
	if(.not.allocated(this% vmod_nrel))			allocate(this% vmod_nrel (nvmod))
	!if(.not.allocated(this% vmod_intep_matrix))			allocate(this% vmod_intep_matrix (nvmod,nnodes))
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
	!
	!if(allocated(this% vmod_hsat))		deallocate(this% vmod_hsat )
	!if(allocated(this% vmod_hsattemp))deallocate(this% vmod_hsattemp)
	!if(allocated(this% vmod_hsatold))	deallocate(this% vmod_hsatold )
	if(allocated(this% vmod_qent))		deallocate(this% vmod_qent)
	!if(allocated(this% vmod_qenttemp))deallocate(this% vmod_qenttemp)
	!if(allocated(this% vmod_qentold))	deallocate(this% vmod_qentold)
	if(allocated(this% vmod_nrel))			deallocate(this% vmod_nrel)
	!if(allocated(this% vmod_intep_matrix))			deallocate(this% vmod_intep_matrix)
	end subroutine s_sat_mesh_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_nodes. Stablish soil shape (coord z) from slope.
	!> Modify this subroutine to define other soil shape.
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

	!!---------------------------------------------------------------------------------------------------------------------
	!!> @author Iván Campos-Guereta Díez
	!!> @brief
	!!> Procedure inside the class ty_sat_nodes. Stablish soil shape (coord z) from slope.
	!!> Modify this subroutine to define other soil shape.
	!!---------------------------------------------------------------------------------------------------------------------
	!
	!	subroutine s_sat_mesh_construct_interpmatrix(this,nodes)
	!	use com_mod_ty_nodes, only:ty_com_nodes
	!	use com_mod_fem_shapefunctions, only:calc_linear_interpolation_matrix
	!
	!	class(ty_sat_mesh),intent(inout)::this
	!	class(ty_com_nodes),intent(in)::nodes
	!	real(kind=dpd)::xvmod(this%vmod_count+1)
	!	real(kind=dpd)::tempmatrix(nodes%count,this%vmod_count+1)
	!
	!	xvmod(1:this%vmod_count) = this%vmod_x
	!	xvmod(this%vmod_count+1) = this%width
	!	tempmatrix = calc_linear_interpolation_matrix(nodes%x,xvmod)
	!
	!	this%vmod_intep_matrix = tempmatrix(:,1:this%vmod_count)
	!	this%vmod_intep_matrix(:,this%vmod_count) = tempmatrix(:,this%vmod_count)+tempmatrix(:,this%vmod_count+1) !This is to interpolate over the final nvmod point.
	!
	!	end subroutine	s_sat_mesh_construct_interpmatrix


	end module sat_mod_ty_mesh