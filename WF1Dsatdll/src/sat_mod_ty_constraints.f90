	!********************************************************************************************************************
	!        MODULE: CLASS OF COLLECTION OF CONSTRAINTS CALCULATIONS FOR SATURATED MODEL
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : mod_sat_ty_layers
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

	module sat_mod_ty_constraints
	use com_mod_ty_layers, only: ty_com_layers
	use com_mod_ty_pointers, only: ty_com_pointer_real,ty_com_pointer_integer

	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_sat_constraints !< CLASS: layers collection in sat model
		integer::count !<number of vertical modules
		type(ty_com_pointer_real),allocatable			::x(:)			!<x coord of every vertical module | point to x(idnode)
		type(ty_com_pointer_real),allocatable			::z(:)			!<z coord of every vertical module | point to z(idnode)
		type(ty_com_pointer_integer),allocatable	::idnode(:) !<node number of every vertical module | from input point to mesh%vmod_idnod(1:nu)

		type(ty_com_pointer_real),allocatable			::hsat(:)			!<hsat			| points to hnew	(idnode)
		type(ty_com_pointer_real),allocatable			::hsattemp(:)	!<hsattemp	| points to htemp	(idnode)
		type(ty_com_pointer_real),allocatable			::hsatold(:)	!<hsatold		| points to hold	(idnode)

		type(ty_com_pointer_real),allocatable			::qent(:) !<Inflow into water-table | point to constructor in 1DSAT | point to unsat(iu)%constratints%qver(is) in 15DSATUNSAT.

		type(ty_com_pointer_real),allocatable			::dqhordx(:)				!< Increment of qh respect to x | points to this%dqhordx_mean(1:nu) (so mean value has to be updated)
		type(ty_com_pointer_real),allocatable			::nrel(:)						!< nrel= non-wetting porosity/porosity | points to mesh%vmod_nrel(1:nu) in 1DSAT | points to  unsat(iu)%constratints%nrel(is) in 15DSATUNSAT
		real(kind=dps),allocatable								::intep_matrix(:,:)	!< Interpolation matrix to get values in nodes from values in vertical modules.

		integer,allocatable												::idelem0(:)	!<id of start sat element associated to iu
		integer,allocatable												::idelem1(:)	!<id of final sat element associated to iu
		integer,allocatable												::idnode0(:)	!<id of start sat node associated to iu
		integer,allocatable												::idnode1(:)	!<id of final sat node associated to iu
		real(kind=dpd),allocatable								::width(:)		!<width associated to iu (length from idnode0 to idnode1)

		real(kind=dpd),allocatable			::hsat_mean(:)				!<	hsat mean between idnode0 and idnode 1 associated to iu				| updated with: sat(is)%model%put_results_in_constraints()
		real(kind=dpd),allocatable			::qent_mean(:)				!<	qent mean between idnode0 and idnode 1 associated to iu				|	updated with: sat(is)%model%put_results_in_constraints()
		real(kind=dpd),allocatable			::incvoldt_mean(:)		!<	incvoldt mean between idnode0 and idnode 1 associated to iu		| updated with: sat(is)%model%put_results_in_constraints()
		real(kind=dpd),allocatable			::dqhordx_mean(:)			!<	dqhordx mean between idnode0 and idnode 1 associated to iu		| updated with: sat(is)%model%put_results_in_constraints()
		real(kind=dpd),allocatable			::dqhordx_all_mean(:)	!<	dqhordx_al mean between idnode0 and idnode 1 associated to iu | updated with: sat(is)%model%put_results_in_constraints()
		real(kind=dpd),allocatable			::qhor0(:)						!<	qh on the first element associated to iu | updated with: sat(is)%model%put_results_in_constraints()
		real(kind=dpd),allocatable			::qhor1(:)						!<	qh on the final element associated to iu | updated with: sat(is)%model%put_results_in_constraints()

	contains
	procedure,public:: allocateall		=> s_sat_constraints_allocateall		!< allocateall(nvmod,nnodes): Allocate all arrays in this instance.
	procedure,public:: deallocateall	=> s_sat_constraints_deallocateall	!< deallocateall(): deallocate all arrays in this instance.
	procedure,public:: construct			=> s_sat_constraints_construct			!< construct(idnode(:),nodesarg(:),qent(:),nrel(:)): Allocate and point idnode(:) to arg idnode, and point x(:),z(:),hsat(:),hsattemp(:),hsatold(:) to nodes%x..., and point qent(:) to arg qent and nrel(:) to arg nrel

	end type ty_sat_constraints

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Initalize count and allocate layer(:) (Override procedure)
	!> @param[in] nlayers
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_constraints_allocateall(this,nvmod,nnodes)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_constraints_allocateall
	!DEC$ endif
	class(ty_sat_constraints),intent(inout)::this
	integer,intent(in)::nvmod
	integer,intent(in)::nnodes

	this%count = nvmod

	if(.not.allocated(this% x))				allocate(this% x (nvmod))
	if(.not.allocated(this% z))				allocate(this% z (nvmod))
	if(.not.allocated(this% idnode))		allocate(this% idnode (nvmod))

	if(.not.allocated(this% hsat))		allocate(this% hsat (nvmod))
	if(.not.allocated(this% hsattemp))allocate(this% hsattemp (nvmod))
	if(.not.allocated(this% hsatold))	allocate(this% hsatold (nvmod))
	if(.not.allocated(this% qent))		allocate(this% qent (nvmod))
	!if(.not.allocated(this% qenttemp)) allocate(this% qenttemp (nvmod))
	!if(.not.allocated(this% qentold))	allocate(this% qentold (nvmod))
	if(.not.allocated(this% dqhordx))	allocate(this% dqhordx (nvmod))
	if(.not.allocated(this% nrel))		allocate(this% nrel (nvmod))
	if(.not.allocated(this% intep_matrix))			allocate(this% intep_matrix (nvmod,nnodes))

	if(.not.allocated(this% idelem0))		allocate(this% idelem0 (nvmod))
	if(.not.allocated(this% idelem1))		allocate(this% idelem1 (nvmod))
	if(.not.allocated(this% idnode0))		allocate(this% idnode0 (nvmod))
	if(.not.allocated(this% idnode1))		allocate(this% idnode1 (nvmod))
	if(.not.allocated(this% width))			allocate(this% width (nvmod))

	if(.not.allocated(this% hsat_mean))			allocate(this% hsat_mean (nvmod))
	if(.not.allocated(this% qent_mean))			allocate(this% qent_mean (nvmod))
	if(.not.allocated(this% incvoldt_mean))	allocate(this% incvoldt_mean (nvmod))
	if(.not.allocated(this% dqhordx_mean))	allocate(this% dqhordx_mean (nvmod))
	if(.not.allocated(this% dqhordx_all_mean))	allocate(this% dqhordx_all_mean (nvmod))
	if(.not.allocated(this% qhor0))	allocate(this% qhor0 (nvmod))
	if(.not.allocated(this% qhor1))	allocate(this% qhor1 (nvmod))

	end subroutine s_sat_constraints_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate all arrays
	!---------------------------------------------------------------------------------------------------------------------


	subroutine s_sat_constraints_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_constraints_deallocateall
	!DEC$ endif
	class(ty_sat_constraints),intent(inout)::this

	if(allocated(this% x))				deallocate(this% x )
	if(allocated(this% z))				deallocate(this% z )
	if(allocated(this% idnode))		deallocate(this% idnode )

	if(allocated(this% hsat))			deallocate(this% hsat )
	if(allocated(this% hsattemp))	deallocate(this% hsattemp)
	if(allocated(this% hsatold))	deallocate(this% hsatold )
	if(allocated(this% qent))			deallocate(this% qent)
	!if(allocated(this% qenttemp))	deallocate(this% qenttemp)
	!if(allocated(this% qentold))	deallocate(this% qentold)
	if(allocated(this% dqhordx))	deallocate(this% dqhordx)
	if(allocated(this% nrel))			deallocate(this% nrel)
	if(allocated(this% intep_matrix))	deallocate(this% intep_matrix)

	if(allocated(this% idelem0))		deallocate(this% idelem0 )
	if(allocated(this% idelem1))		deallocate(this% idelem1 )
	if(allocated(this% idnode0))		deallocate(this% idnode0 )
	if(allocated(this% idnode1))		deallocate(this% idnode1 )
	if(allocated(this% width))			deallocate(this% width )

	if(allocated(this% hsat_mean))			deallocate(this% hsat_mean )
	if(allocated(this% qent_mean))			deallocate(this% qent_mean )
	if(allocated(this% incvoldt_mean))	deallocate(this% incvoldt_mean )
	if(allocated(this% dqhordx_mean))	deallocate(this% dqhordx_mean )
	if(allocated(this% dqhordx_all_mean))	deallocate(this% dqhordx_all_mean )
	if(allocated(this% qhor0))	deallocate(this% qhor0 )
	if(allocated(this% qhor1))	deallocate(this% qhor1 )

	end subroutine s_sat_constraints_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate and point idnode(:) to argument idnode. Point x(:) to nodesarg%x(:), z to nodesarg%z(:),
	!> hsat(:) to nodesarg%hnew(:), hsattemp(:) and hsatold(:) the same. qent(:) points to argument qent,
	!> and nrel(:) points to argument nrel.
	!> @param[inout] idnode
	!> @param[inout] nodesarg
	!> @param[inout] qent
	!> @param[inout] nrel
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_constraints_construct(this,idnode,nodesarg,qent,nrel)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_constraints_construct
	!DEC$ endif
	use com_mod_ty_nodes, only:ty_com_nodes
	use sat_mod_ty_nodes, only:ty_sat_nodes

	class(ty_sat_constraints),intent(inout),target::this
	integer,intent(inout),target::idnode(:)
	class(ty_com_nodes),intent(inout),target::nodesarg
	real(kind=dpd),intent(inout),target,optional::qent(:)
	real(kind=dpd),intent(inout),target,optional::nrel(:)

	class(ty_com_nodes),pointer::nodescom
	class(ty_sat_nodes),pointer::nodes
	integer::i

	nodescom => nodesarg
	select type (nodescom)
	type is(ty_sat_nodes)
		nodes => nodescom
	end select

	this%count = size(idnode)
	if(.not.allocated(this%idnode)) allocate(this%idnode(this%count))


	do i=1,this%count
		this%idnode(i)%p	=> idnode(i)
		this%x(i)%p		=> nodes%x(this%idnode(i)%p)
		this%z(i)%p		=> nodes%z(this%idnode(i)%p)

		this%hsat(i)%p			=> nodes%hnew(this%idnode(i)%p)
		this%hsattemp(i)%p	=> nodes%htemp(this%idnode(i)%p)
		this%hsatold(i)%p		=> nodes%hold(this%idnode(i)%p)
		if(present(qent)) then
			this%qent(i)%p		=> qent(i)
		end if

		if(present(nrel)) then
			this%nrel(i)%p		=> nrel(i)
		end if


		this%dqhordx(i)%p		=> nodes%results_dqhordx_from_incvoldt(this%idnode(i)%p) !Take results_dqhordx_from_incvol because expect that calc is more precise than results_dqhordx
		!this%dqhordx(i)%p		=> this%dqhordx_mean(i)	 !Take results_dqhordx_from_incvol because expect that calc is more precise than results_dqhordx


	end do

	this%idnode0(1) = 1
	this%idnode1(1) = ((this%idnode0(1)+this%idnode(2)%p)/2)
	do i=2,this%count-1
		this%idnode0(i) = this%idnode1(i-1)
		this%idnode1(i) = ((this%idnode0(i)+this%idnode(i+1)%p)/2)
	end do
	this%idnode0(this%count) = this%idnode1(this%count-1)
	this%idnode1(this%count) = nodes%count

	this%idelem0 = this%idnode0
	this%idelem1 = this%idnode1-1

	do i=1,this%count
		this%width(i)=nodes%x(this%idnode1(i))-nodes%x(this%idnode0(i))
	end do

	call s_sat_constraints_construct_interpmatrix(this,nodes)

	this%hsat_mean = 0.0_dpd
	this%qent_mean = 0.0_dpd
	this%incvoldt_mean = 0.0_dpd
	this%dqhordx_mean = 0.0_dpd
	this%dqhordx_all_mean = 0.0_dpd
	this%qhor0 = 0.0_dpd
	this%qhor1 = 0.0_dpd

	end subroutine s_sat_constraints_construct



	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_nodes. Stablish soil shape (coord z) from slope.
	!> Modify this subroutine to define other soil shape.
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_constraints_construct_interpmatrix(this,nodes)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_constraints_construct_interpmatrix
	!DEC$ endif
	use com_mod_ty_nodes, only:ty_com_nodes
	use com_mod_fem_shapefunctions, only:calc_linear_interpolation_matrix

	class(ty_sat_constraints),intent(inout)::this
	class(ty_com_nodes),intent(in)::nodes
	real(kind=dpd)::xvmod(this%count+1)
	real(kind=dpd)::tempmatrix(nodes%count,this%count+1)

	integer::i

	do i=1,this%count
		xvmod(i) = nodes%x(this%idnode(i)%p)
	end do
	xvmod(this%count+1) = nodes%x(nodes%count)
	tempmatrix = calc_linear_interpolation_matrix(nodes%x,xvmod)

	this%intep_matrix = tempmatrix(:,1:this%count)
	this%intep_matrix(:,this%count) = tempmatrix(:,this%count)+tempmatrix(:,this%count+1) !This is to interpolate over the final nvmod point.

	end subroutine	s_sat_constraints_construct_interpmatrix

	end module sat_mod_ty_constraints