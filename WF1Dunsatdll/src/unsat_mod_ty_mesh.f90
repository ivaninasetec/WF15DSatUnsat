	!********************************************************************************************************************
	!        MODULE: CLASS WITH MESH CHARACTERISTICS OF UNSATURATED MODEL
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

	module unsat_mod_ty_mesh

	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_unsat_mesh
		integer::nn		!<number of nodes per element
		integer::nc		!<number of clases of nodes (0=c0, 1=c1,...)
		integer::nnc	!<number of clases and nodes in elements (0=c0, 1=c1,...)
		integer::nelemv_count	!<number of elements
		integer,allocatable::nelemv(:)
		integer::nnodv_count			!<number of nodes
		integer::nnodclassv_count	!<number of nodes and classes
		real(kind=dps)::height		!<total leght of the model

		!real(kind=dps),allocatable::	hsat(:)			!<Value of hsat for the node (possitive pressure at node)
		!real(kind=dps),allocatable::	hsatold(:)	!<Value of hsat for the node (old value) (possitive pressure at node)
		!integer,allocatable::					idnodsat(:)	!<Node at which nsat reaches
		integer,allocatable::					idnode_bottom(:)	!<Node at wich layer start
		integer,allocatable::					idnode_top(:)	!<Node at wich layer start
		integer,allocatable::					idelem_bottom(:)	!<Node at wich layer start
		integer,allocatable::					idelem_top(:)	!<Node at wich layer start
		!real(kind=dps),allocatable::	incvsat(:)	!<Increment of water content	from hsat_old to hsat_new = Int(th,hsatold,hsat)
		!real(kind=dps),allocatable::	qvsat(:)		!<Vertical water flow to pass from hsatold to hsat = incvsat/Dt
		!real(kind=dps),allocatable::	nmean(:)		!<Mean value for the non wetting available porosity in between hsatold to hsat =incvsat/abs(hsat-hsatold)
		!real(kind=dps),allocatable::	nrel(:)			!<relative_n=nmean/n
		real(kind=dps),allocatable::	dqhordx(:)	!<Water leaving the watertable (dqh/dt)

	contains
	procedure,public:: allocateall		=> s_unsat_mesh_allocateall
	procedure,public:: deallocateall	=> s_unsat_mesh_deallocateall

	end type ty_unsat_mesh

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_unsat_mesh_allocateall(this,nlayers)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_mesh_allocateall" :: s_unsat_mesh_allocateall
	!DEC$ endif
	class(ty_unsat_mesh),intent(inout)::this
	integer,intent(in)::nlayers

	if(.not.allocated(this% nelemv))			allocate(this% nelemv(nlayers))
	!if(.not.allocated(this% hsat))				allocate(this% hsat(nlayers))
	!if(.not.allocated(this% hsatold))			allocate(this% hsatold(nlayers))
	!if(.not.allocated(this% idnodsat))		allocate(this% idnodsat(nlayers))
	if(.not.allocated(this% idnode_bottom))	allocate(this% idnode_bottom(nlayers))
	if(.not.allocated(this% idnode_top))		allocate(this% idnode_top(nlayers))
	if(.not.allocated(this% idelem_bottom))	allocate(this% idelem_bottom(nlayers))
	if(.not.allocated(this% idelem_top))		allocate(this% idelem_top(nlayers))
	!if(.not.allocated(this% incvsat))			allocate(this% incvsat(nlayers))
	!if(.not.allocated(this% qvsat))				allocate(this% qvsat(nlayers))
	!if(.not.allocated(this% nmean))				allocate(this% nmean(nlayers))
	!if(.not.allocated(this% nrel))				allocate(this% nrel(nlayers))
	if(.not.allocated(this% dqhordx))			allocate(this% dqhordx(nlayers))
	this% dqhordx = 0.0_dpd

	end subroutine s_unsat_mesh_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate all arrays
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_mesh_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_mesh_deallocateall" :: s_unsat_mesh_deallocateall
	!DEC$ endif
	class(ty_unsat_mesh),intent(inout)::this

	if(allocated(this% nelemv))			deallocate(this% nelemv )
	!if(allocated(this% hsat))				deallocate(this% hsat )
	!if(allocated(this% hsatold))		deallocate(this% hsatold )
	!if(allocated(this% idnodsat))		deallocate(this% idnodsat )
	if(allocated(this% idnode_bottom))	deallocate(this% idnode_bottom )
	!if(allocated(this% incvsat))		deallocate(this% incvsat )
	!if(allocated(this% qvsat))			deallocate(this% qvsat )
	!if(allocated(this% nmean))			deallocate(this% nmean )
	!if(allocated(this% nrel))				deallocate(this% nrel )
	if(allocated(this% dqhordx))		deallocate(this% dqhordx )

	end subroutine s_unsat_mesh_deallocateall

	end module unsat_mod_ty_mesh