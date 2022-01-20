	!********************************************************************************************************************
	!        MODULE: CLASS WITH MESH CHARACTERISTICS OF THE WHOLE MODEL
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : unsat_mod_ty_mesh
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

	module model_mod_ty_mesh

	implicit none
	include 'inc_precision.fi'

	private
	
  type,public::ty_model_mesh
    integer::nn		!<number of nodes per element
    integer::nc		!<number of clases of nodes (0=c0, 1=c1,...)
    integer::nnc	!<number of clases and nodes in elements (0=c0, 1=c1,...)
		!integer,pointer::layers_count	!<number of layers
		integer::nelemv_count !Total vertical elements
		integer::nelemh_count !Total horizontal elements
		integer,allocatable::nelemv(:) !Number of vertical elements in each layer.
		integer,allocatable::nelemh(:) !Number of vertical elements per vertical modules
		integer::nnodh_count			!<number of nodes (horizontal)
		integer::nnodclassh_count	!<number of nodes and classes (horizontal)
		integer::nnodv_count			!<number of nodes (vertical)
		integer::nnodclassv_count	!<number of nodes and classes (vertical)
		real(kind=dps),pointer::width		!<total width of the model (point to layers%width
		real(kind=dps)::height	!<total height of the model (sum of layers%height)
		
    integer::vmod_count !<number of vertical modules		
    real(kind=dps),allocatable	::vmod_x(:) !<x coord of every vertical module		
    integer,allocatable					::vmod_idnod(:) !<node number of every vertical module (before numnod)
		!integer,allocatable					::vmod_nelem(:) !<Number of elements between vertical modules (first always at x=0)
  !  real(kind=dps),allocatable	::vmod_hsat(:)
		!real(kind=dps),allocatable	::vmod_hsattemp(:)		
		!real(kind=dps),allocatable	::vmod_hsatold(:)
  !
  !  real(kind=dps),allocatable	::vmod_qent(:) !especific flow in each vertical module node
		!real(kind=dps),allocatable	::vmod_qenttemp(:) !especific flow in each vertical module node
		!real(kind=dps),allocatable	::vmod_qentold(:) !especific flow in each vertical module node 		
		
	contains
		procedure,public:: allocateall		=> s_model_mesh_allocateall
		procedure,public:: deallocateall	=> s_model_mesh_deallocateall
		
	end type ty_model_mesh	

	contains
		
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_model_mesh_allocateall(this,nlayers,nvmod)
	class(ty_model_mesh),intent(inout)::this
	integer,intent(in)::nlayers
	integer,intent(in)::nvmod
	
	if(.not.allocated(this% nelemv))		allocate(this% nelemv(nlayers))
	if(.not.allocated(this% nelemh))		allocate(this% nelemh(nvmod+1))
	this%vmod_count = nvmod
	if(.not.allocated(this% vmod_x))				allocate(this% vmod_x (nvmod))
	if(.not.allocated(this% vmod_idnod))		allocate(this% vmod_idnod (nvmod))
	!if(.not.allocated(this% vmod_nelem))		allocate(this% vmod_nelem (nvmod))
	
	end subroutine s_model_mesh_allocateall
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	
	subroutine s_model_mesh_deallocateall(this)
	class(ty_model_mesh),intent(inout)::this
	
	if(allocated(this% nelemv))				deallocate(this% nelemv )
	if(allocated(this% nelemh))				deallocate(this% nelemh )
	if(allocated(this% vmod_x))				deallocate(this% vmod_x )
	if(allocated(this% vmod_idnod))		deallocate(this% vmod_idnod )
	!if(allocated(this% vmod_nelem))		deallocate(this% vmod_nelem )
	
	end subroutine s_model_mesh_deallocateall
	
	end module model_mod_ty_mesh