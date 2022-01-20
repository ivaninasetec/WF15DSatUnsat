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

	module unsat_mod_ty_constraints
	use com_mod_ty_pointers, only: ty_com_pointer_real,ty_com_pointer_integer

	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_unsat_constraints !< CLASS: layers collection in sat model
		integer::count !< number of vertical modules
		type(ty_com_pointer_real),allocatable			::x(:) !<x(is) x coord where is 1DSAT cross | points to nodes%x(idnode_layer_bottom)
		type(ty_com_pointer_integer),allocatable	::idnode_layer_bottom(:) !<id of node of bottom of layer assigned to 1DSAT | pointer is assigned in construct, inside model%construct in 15DSAT | points to mesh%idnode_bottom(is) in 1DUNSAT | points to mesh%vmod_idnod in 15DSATUNSAT
		type(ty_com_pointer_integer),allocatable	::idelem_layer_bottom(:) !<Same as idnode_layer_bottom but id of element
		type(ty_com_pointer_integer),allocatable	::idelem_layer_top(:) !<Same as idnode_layer_bottom but id of element on top of layer is
				type(ty_com_pointer_real),allocatable	::hnew(:)			!<hsat on 1DUNSAT | need to be updated with 'update_hsat(elements)', or with update_all(calc,dt)
		type(ty_com_pointer_real),allocatable	::thnew(:)			!<hsat on 1DUNSAT | need to be updated with 'update_hsat(elements)', or with update_all(calc,dt)

		integer,allocatable::idelem_hsat(:)			!<id of element where watertable is located for 1DSAT simulation is
		integer,allocatable::idelem_top(:)			!<id of element where watertable is located for 1DSAT or top of layer if water-table exceed.
		integer,allocatable::idelem_hsatold(:)	!<id of element where old watertable is located.

		real(kind=dps),allocatable		::hsat(:)			!<hsat on 1DUNSAT | need to be updated with 'update_hsat(elements)', or with update_all(calc,dt)
		!real(kind=dps),allocatable		::hsattemp(:)	!<hsattemp on 1DUNSAT
		real(kind=dps),allocatable		::hsatold(:)	!<hsatold on 1DUNSAT | need to be updated with 'update_hsat(elements)', or with update_all(calc,dt)
		real(kind=dps),allocatable		::xhsat(:)		!<location in z of hsat on 1DUNSAT | need to be updated with 'update_hsat(elements)'
		!real(kind=dps),allocatable		::xhsattemp(:)	!<location in z of hsat on 1DUNSAT | need to be updated with 'update_hsat(elements)'
		real(kind=dps),allocatable		::xhsatold(:) !<location in z of hsatold on 1DUNSAT | need to be updated with 'update_hsat(elements)'
		real(kind=dps),allocatable		::chi_hsat(:)	!<relative coord chi inside the element where hsat is located in 1DUNSAT | need to be updated with 'update_hsat(elements)'

		real(kind=dps),allocatable		::qver(:)			!<Water entering on the top of water-table minus leaving on layer below | updated with update_qsupinf() or with update_all(calc,dt)
		real(kind=dps),allocatable		::incvoldt(:) !<Increment of volume filled by watertable in the timestep divided by timestep | updated with update_incvoldt(calc,dt) or update_all(calc,dt)
		type(ty_com_pointer_real),allocatable		::dqhordx(:)	!<Points to sat dqhordx, would be the flow leaving 1DSAT because of horizontal flow | always pointing ot 1DSAT(is)%contraints%dqhordx(iu)
		real(kind=dps),allocatable		::qsat(:) !<Vertical flow just on the top of hsat | updated with update_qsat(calc) or from update_all(calc,dt)
		real(kind=dps),allocatable		::qtop(:) !<Vertical flow just on the top of hsat or top of layer (if wt is over layer) | updated with update_qtop(calc) (CHECK: not used in this moment)
		real(kind=dps),allocatable		::qinf(:) !<Vertical flow just below the layer (on the other layer element) | Updated with update_qinf(calc) or with update_all(calc,dt)
		real(kind=dps),allocatable		::qsupinf(:) !<Vertical flow just below the layer (on the other layer element) | Updated with update_qsupinf(calc) or with update_all(calc,dt)
		type(ty_com_pointer_real),allocatable		::qnewmann(:) !This is the result water leaving when Dirichlet conditions are imposed | pointing to nodes%results_qnewmann(idnode_layer_bottom) (nodes have to be updated
		real(kind=dps),allocatable		::qnewmann_all(:) !Sum all qnewmann until hsat from the sum of qnewmann | updated with update_qnewmann_all()

		!real(kind=dps),allocatable	::qenttemp(:) !especific flow in each vertical module node
		!real(kind=dps),allocatable	::qentold(:) !especific flow in each vertical module node


		real(kind=dps),allocatable	::nrel(:)								!Relative porosity of non-wetting voids [int(th,hold,hnew)/(|Hnew-hold|·(thsat-thres))]

	contains
	procedure,public:: allocateall		=> s_unsat_constraints_allocateall
	procedure,public:: deallocateall	=> s_unsat_constraints_deallocateall
	procedure,public:: construct			=> s_unsat_constraints_construct

	procedure,public:: update_hsat			=> s_unsat_model_update_hsat_from_hbottom
	procedure,public:: update_incvoldt	=> s_unsat_constraints_update_incvoldt
	procedure,public:: update_qsat			=> s_unsat_constraints_update_qsat_from_qmed
	procedure,public:: update_qtop			=> s_unsat_constraints_update_qtop_from_qmed
	procedure,public:: update_qinf			=> s_unsat_constraints_update_qinf_from_qtop
	procedure,public:: update_qsupinf		=> s_unsat_constraints_update_qsupinf
	procedure,public:: update_qnewmann_all	=> s_unsat_constraints_update_qnewman_all
	procedure,public:: update_all				=> s_unsat_model_updateall	!<update_all(calc,dt): Updates hsat,hsatold,xhsat,xhsatold,chi_hsat,qver,incvoldt,qsat,qtop,qinf,qsuping,qnewmann,qnewmannall.

	procedure,public:: set_newmann_from_dqhordx	=> s_unsat_constraint_set_newmann_from_dqhordx
	!
	end type ty_unsat_constraints

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Initalize count and allocate layer(:) (Override procedure)
	!> @param[in] nlayers
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_constraints_allocateall(this,nvmod,nnodes)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_allocateall" :: s_unsat_constraints_allocateall
	!DEC$ endif
	class(ty_unsat_constraints),intent(inout)::this
	integer,intent(in)::nvmod
	integer,intent(in)::nnodes

	this%count = nvmod
	if(.not.allocated(this% x))					allocate(this% x (nvmod))
	if(.not.allocated(this% idnode_layer_bottom))		allocate(this% idnode_layer_bottom (nvmod))
	if(.not.allocated(this% idelem_layer_bottom))		allocate(this% idelem_layer_bottom (nvmod))
	if(.not.allocated(this% idelem_layer_top))			allocate(this% idelem_layer_top (nvmod))

	!if(.not.allocated(this% idnode_hsat))		allocate(this% idnode_hsat (nvmod))
	!if(.not.allocated(this% idnode_hsattemp))		allocate(this% idnode_hsattemp (nvmod))
	!if(.not.allocated(this% idnode_hsatold))		allocate(this% idnode_hsatold (nvmod))

	if(.not.allocated(this% idelem_hsat))		allocate(this% idelem_hsat (nvmod))
	if(.not.allocated(this% idelem_top))		allocate(this% idelem_top (nvmod))
	if(.not.allocated(this% idelem_hsatold))		allocate(this% idelem_hsatold (nvmod))

	if(.not.allocated(this% hsat))			allocate(this% hsat (nvmod))
	if(.not.allocated(this% hnew))			allocate(this% hnew (nvmod))
	if(.not.allocated(this% thnew))			allocate(this% thnew (nvmod))
	
	if(.not.allocated(this% chi_hsat))			allocate(this% chi_hsat (nvmod))
	!if(.not.allocated(this% hsattemp))	allocate(this% hsattemp (nvmod))
	if(.not.allocated(this% hsatold))		allocate(this% hsatold (nvmod))
	if(.not.allocated(this% xhsat))			allocate(this% xhsat (nvmod))
	!if(.not.allocated(this% xhsattemp))	allocate(this% xhsattemp (nvmod))
	if(.not.allocated(this% xhsatold))	allocate(this% xhsatold (nvmod))
	if(.not.allocated(this% incvoldt))	allocate(this% incvoldt (nvmod))
	if(.not.allocated(this% qver))			allocate(this% qver (nvmod))
	if(.not.allocated(this% dqhordx))		allocate(this% dqhordx (nvmod))
	if(.not.allocated(this% nrel))			allocate(this% nrel (nvmod))
	if(.not.allocated(this% qsat))			allocate(this% qsat (nvmod))
	if(.not.allocated(this% qtop))			allocate(this% qtop (nvmod))
	if(.not.allocated(this% qinf))			allocate(this% qinf (nvmod))
	if(.not.allocated(this% qsupinf))			allocate(this% qsupinf (nvmod))
	if(.not.allocated(this% qnewmann))			allocate(this% qnewmann (nvmod))
	if(.not.allocated(this% qnewmann_all))			allocate(this% qnewmann_all (nvmod))

	end subroutine s_unsat_constraints_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate all arrays
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_constraints_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_deallocateall" :: s_unsat_constraints_deallocateall
	!DEC$ endif
	class(ty_unsat_constraints),intent(inout)::this

	if(allocated(this% x))						deallocate(this% x )
	if(allocated(this% idnode_layer_bottom))				deallocate(this% idnode_layer_bottom )
	!if(allocated(this% idnode_hsat))	deallocate(this% idnode_hsat )
	!if(allocated(this% idnode_hsattemp))	deallocate(this% idnode_hsattemp )
	!if(allocated(this% idnode_hsatold))	deallocate(this% idnode_hsatold )
	if(allocated(this% idelem_hsat))	deallocate(this% idelem_hsat )
	if(allocated(this% idelem_hsatold))	deallocate(this% idelem_hsatold )

	if(allocated(this% hsat))			deallocate(this% hsat )
	if(allocated(this% hsat))			deallocate(this% hsat )
	if(allocated(this% hnew))			deallocate(this% hnew )
	if(allocated(this% chi_hsat))			deallocate(this% chi_hsat )
	!if(allocated(this% hsattemp))	deallocate(this% hsattemp)
	if(allocated(this% hsatold))	deallocate(this% hsatold )
	if(allocated(this% xhsat))		deallocate(this% xhsat )
	!if(allocated(this% xhsattemp))deallocate(this% xhsattemp )
	if(allocated(this% xhsatold))	deallocate(this% xhsatold )
	if(allocated(this% incvoldt))	deallocate(this% incvoldt)
	if(allocated(this% qver))			deallocate(this% qver)
	if(allocated(this% dqhordx))	deallocate(this% dqhordx)
	if(allocated(this% nrel))			deallocate(this% nrel)
	if(allocated(this% qsat))			deallocate(this% qsat)
	if(allocated(this% qtop))			deallocate(this% qtop)
	if(allocated(this% qinf))			deallocate(this% qinf)
	if(allocated(this% qsupinf))	deallocate(this% qsupinf)
	if(allocated(this% qnewmann))	deallocate(this% qnewmann)
	if(allocated(this% qnewmann_all))	deallocate(this% qnewmann_all)

	end subroutine s_unsat_constraints_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate all arrays
	!---------------------------------------------------------------------------------------------------------------------


	subroutine s_unsat_constraints_construct(this,idnode,nodes,dqhordx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_construct" :: s_unsat_constraints_construct
	!DEC$ endif
	use com_mod_ty_nodes, only:ty_com_nodes
	use unsat_mod_ty_nodes, only:ty_unsat_nodes

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_com_nodes),intent(inout),target::nodes
	integer,intent(inout),target::idnode(:)
	real(kind=dpd),intent(inout),target,optional::dqhordx(:)

	class(ty_com_nodes),pointer::nodescom
	class(ty_unsat_nodes),pointer::nodesunsat

	integer::i

	this%count = size(idnode)
	if(.not.allocated(this%idnode_layer_bottom)) allocate(this%idnode_layer_bottom(this%count))

	nodescom=>nodes
	select type(nodescom)
	type is(ty_unsat_nodes)
		nodesunsat=> nodescom
	end select

	do i=1,this%count
		this%idnode_layer_bottom(i)%p	=> idnode(i)
		this%idelem_layer_bottom(i)%p	=> idnode(i)
		if (i==this%count) then
			this%idelem_layer_top(i)%p		=> nodes%id(this%count-1)
		else
			this%idelem_layer_top(i)%p		=> nodes%id(idnode(i+1)-1)
		end if

		!this%idnode_hsat(i) = this%idnode_layer_bottom(i)%p
		this%idelem_hsat(i) = this%idelem_layer_bottom(i)%p
		!this%idelem_top(i)	= this%idelem_hsat(i)

		this%x(i)%p				=> nodes%x(this%idnode_layer_bottom(i)%p)
		this%xhsat(i) = this%x(i)%p
		this%hnew(i)%p				=> nodes%hnew(this%idnode_layer_bottom(i)%p)
		this%thnew(i)%p				=> nodes%thnew(this%idnode_layer_bottom(i)%p)
		if(present(dqhordx)) then
			this%dqhordx(i)%p	=> dqhordx(i)
		end if
		this%qnewmann(i)%p => nodesunsat%results_qnewmann(this%idnode_layer_bottom(i)%p)
	end do

	!this%idnode_hsattemp = this%idnode_hsat
	!this%idnode_hsatold	= this%idnode_hsat
	this%hsat			= 0.0_dpd
	!this%hsattemp = 0.0_dpd
	this%hsatold	= 0.0_dpd
	!this%xhsattemp = this%xhsat
	this%xhsatold = this%xhsat
	this%incvoldt = 0.0_dpd
	this%qver			= 0.0_dpd
	this%nrel			= 0.1_dpd
	this%qnewmann_all			= 0.0_dpd

	end subroutine s_unsat_constraints_construct

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_model_update_hsat_from_hbottom(this,elementsarg)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_model_update_hsat_from_hbottom" :: s_unsat_model_update_hsat_from_hbottom
	!DEC$ endif
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_com_elements),intent(in),target::elementsarg
	class(ty_com_elements),pointer::elementscom
	class(ty_unsat_elements),pointer::elements

	integer::ic,ielem,ielemlybot
	real(kind=dpd)::xsat0,xsat1,xlay,hsat0,hsat,hsat1,p0,p1,psat

	elementscom => elementsarg
	select type(elementscom)
	type is (ty_unsat_elements)
		elements => elementscom
	end select

	!This is more precise (need to be precise), it searches when there is an unsat node over the saturated and then calculat hsat
	!For hnew:
	do ic=1,this%count
		do ielem=this%idelem_layer_bottom(ic)%p,elements%count
			if(elements%h0(ielem)<0.0_dpd) exit
		end do
		if(ielem==this%idelem_layer_bottom(ic)%p) then !Not saturated
			this%idelem_hsat(ic)=ielem
			!this%idelem_top(ic)=this%idelem_hsat(ic)
			this%hsat(ic)		= 0.0_dpd
			this%xhsat(ic)	= elements%x0(ielem)
			this%chi_hsat(ic)	 = -1.0_dpd
		else !Saturated
			ielem = ielem-1
			ielemlybot = this%idelem_layer_bottom(ic)%p
			this%idelem_hsat(ic)=ielem
			!this%idelem_top(ic)=min(this%idelem_hsat(ic),this%idelem_layer_top(ic)%p)
			this%hsat(ic)		= elements%h0(ielem)+(elements%x0(ielem)-elements%x0(ielemlybot))
			this%xhsat(ic)	= elements%h0(ielem)+elements%x0(ielem)
			this%chi_hsat(ic)	 = min(1.0_dpd,2.0_dpd*elements%h0(ielem)/(elements%x1(ielem)-elements%x0(ielem))-1.0_dpd)
		end if
	end do

	!For hold:
	do ic=1,this%count
		do ielem=this%idelem_layer_bottom(ic)%p,elements%count
			if(elements%h0old(ielem)<0.0_dpd) exit
		end do
		if(ielem==this%idelem_layer_bottom(ic)%p) then !Not saturated
			this%idelem_hsatold(ic)=ielem
			this%hsatold(ic)		= 0.0_dpd
			this%xhsatold(ic)	= elements%x0(ielem)
			!this%chi_hsatold(ic)	 = -1.0_dpd
		else !Saturated
			ielem = ielem-1
			ielemlybot = this%idelem_layer_bottom(ic)%p
			this%idelem_hsatold(ic)=ielem
			this%hsatold(ic)		= elements%h0old(ielem)+(elements%x0(ielem)-elements%x0(ielemlybot))
			this%xhsatold(ic)	= elements%h0old(ielem)+elements%x0(ielem)
			!this%chi_hsatold(ic)	 = min(1.0_dpd,2.0_dpd*elements%hold(ielem)/(elements%x1(ielem)-elements%x0(ielem))-1.0_dpd)
		end if
	end do

	end subroutine s_unsat_model_update_hsat_from_hbottom



	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************
	!CHECK: UPDATE Hsat: From elements_hmean is not working.
	subroutine s_unsat_model_update_hsat_from_hmean(this,elements)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_model_update_hsat_from_hmean" :: s_unsat_model_update_hsat_from_hmean
	!DEC$ endif
	use com_mod_ty_elements,only:ty_com_elements

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_com_elements),intent(in)::elements

	integer::ic,ielem,ielemlybot
	real(kind=dpd)::xsat0,xsat1,xlay,hsat0,hsat,hsat1,p0,p1,psat


	!This is more precise (need to be precise), it searches when there is an unsat node over the saturated and then calculat hsat
	!For hnew:
	do ic=1,this%count
		do ielem=this%idelem_layer_bottom(ic)%p,elements%count
			if(elements%hnew(ielem)<0.0_dpd) exit
		end do
		if(ielem==this%idelem_layer_bottom(ic)%p) then !Not saturated
			this%idelem_hsat(ic)=ielem
			!this%idelem_top(ic)=this%idelem_hsat(ic)
			this%hsat(ic)		= 0.0_dpd
			this%xhsat(ic)	= elements%x0(ielem)
			this%chi_hsat(ic)	 = -1.0_dpd
		else !Saturated
			ielem = ielem-1
			ielemlybot = this%idelem_layer_bottom(ic)%p
			this%idelem_hsat(ic)=ielem
			!this%idelem_top(ic)=min(this%idelem_hsat(ic),this%idelem_layer_top(ic)%p)
			this%hsat(ic)		= elements%hnew(ielem)+(elements%x0(ielem)-elements%x0(ielemlybot))
			this%xhsat(ic)	= elements%hnew(ielem)+elements%x0(ielem)
			this%chi_hsat(ic)	 = min(1.0_dpd,2.0_dpd*elements%hnew(ielem)/(elements%x1(ielem)-elements%x0(ielem))-1.0_dpd)
		end if
	end do

	!For hold:
	do ic=1,this%count
		do ielem=this%idelem_layer_bottom(ic)%p,elements%count
			if(elements%hold(ielem)<0.0_dpd) exit
		end do
		if(ielem==this%idelem_layer_bottom(ic)%p) then !Not saturated
			this%idelem_hsatold(ic)=ielem
			this%hsatold(ic)		= 0.0_dpd
			this%xhsatold(ic)	= elements%x0(ielem)
			!this%chi_hsatold(ic)	 = -1.0_dpd
		else !Saturated
			ielem = ielem-1
			ielemlybot = this%idelem_layer_bottom(ic)%p
			this%idelem_hsatold(ic)=ielem
			this%hsatold(ic)		= elements%hold(ielem)+(elements%x0(ielem)-elements%x0(ielemlybot))
			this%xhsatold(ic)	= elements%hold(ielem)+elements%x0(ielem)
			!this%chi_hsatold(ic)	 = min(1.0_dpd,2.0_dpd*elements%hold(ielem)/(elements%x1(ielem)-elements%x0(ielem))-1.0_dpd)
		end if
	end do

	end subroutine s_unsat_model_update_hsat_from_hmean


	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_model_update_hsat_from_nodes(this,nodes)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_model_update_hsat_from_nodes" :: s_unsat_model_update_hsat_from_nodes
	!DEC$ endif
	use com_mod_ty_nodes,only:ty_com_nodes

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_com_nodes),intent(in)::nodes

	integer::ic,inod
	real(kind=dpd)::xsat0,xsat1,xlay,hsat0,hsat,hsat1,p0,p1,psat

	!do ic=1,this%count
	!	this%hsat(ic)			= max(0.0_dpd,nodes%hnew(this%idnode(ic)%p))
	!	this%hsattemp(ic)	= max(0.0_dpd,nodes%htemp(this%idnode(ic)%p))
	!	this%hsatold(ic)		= max(0.0_dpd,nodes%hold(this%idnode(ic)%p))
	!	this%xhsat(ic)			= nodes%x(this%idnode(ic)%p)
	!	this%idnode_hsat(ic)= this%idnode(ic)%p
	!end do

	!This is more precise (need to be precise), it searches when there is an unsat node over the saturated and then calculat hsat
	!For hnew:
	do ic=1,this%count
		this%idelem_hsat(ic) = this%idnode_layer_bottom(ic)%p
		do inod=this%idnode_layer_bottom(ic)%p,nodes%count
			if(nodes%hnew(inod)<0.0_dpd) exit
		end do
		if(inod==this%idnode_layer_bottom(ic)%p) then !Not saturated
			this%idelem_hsat(ic)=inod
			this%hsat(ic)		= 0.0_dpd
			this%xhsat(ic)	= nodes%x(inod)
			this%chi_hsat(ic)	 = -1.0_dpd
		else if (inod>=nodes%count) then
			inod=nodes%count
			if(nodes%hnew(inod)>0.0_dpd) then !Saturation surpasses the top surface
				xsat0 = nodes%x(inod)
				xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
				hsat = nodes%hnew(inod)+xsat0-xlay

				this%idelem_hsat(ic)=inod
				this%hsat(ic)		= hsat
				this%xhsat(ic)	= xlay+hsat
				this%chi_hsat(ic)	 = 1.0_dpd

			else !Saturation just below top surface
				xsat0 = nodes%x(inod-1)
				xsat1 = nodes%x(inod)
				xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
				p0 = nodes%hnew(inod-1)
				p1 = nodes%hnew(inod)
				hsat0 = p0+xsat0-xlay
				hsat1 = p1+xsat1-xlay
				hsat = 0.5_dpd*(hsat1+hsat0)
				psat = 0.0_dpd

				this%idelem_hsat(ic)=inod-1
				this%hsat(ic)		= hsat
				this%xhsat(ic)	= xsat0+(xsat1-xsat0)*(psat-p0)/(p1-p0)
				this%chi_hsat(ic)	= (this%xhsat(ic)-xsat0)/(xsat1-xsat0)
			end if
		else
			xsat0 = nodes%x(inod-1)
			xsat1 = nodes%x(inod)
			xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
			p0 = nodes%hnew(inod-1)
			p1 = nodes%hnew(inod)
			hsat0 = p0+xsat0-xlay
			hsat1 = p1+xsat1-xlay
			hsat = 0.5_dpd*(hsat1+hsat0)
			psat = 0.0_dpd

			this%idelem_hsat(ic)=inod-1
			this%hsat(ic)		= hsat
			this%xhsat(ic)	= xsat0+(xsat1-xsat0)*(psat-p0)/(p1-p0)
			this%chi_hsat(ic)	= (this%xhsat(ic)-xsat0)/(xsat1-xsat0)
		end if
	end do



	!!For htemp:
	!do ic=1,this%count
	!	this%idnode_hsattemp(ic) = this%idnode_layer_bottom(ic)%p
	!	do inod=this%idnode_layer_bottom(ic)%p,nodes%count
	!		if(nodes%htemp(inod)<0.0_dpd) exit
	!	end do
	!	if(inod==this%idnode_layer_bottom(ic)%p) then !Not saturated
	!		this%idnode_hsattemp(ic)=inod
	!		this%hsattemp(ic)		= 0.0_dpd
	!		this%xhsattemp(ic)	= nodes%x(inod)
	!	else if (inod>=nodes%count) then
	!		inod=nodes%count
	!		if(nodes%htemp(inod)>0.0_dpd) then !Saturation surpasses the top surface
	!		xsat0 = nodes%x(inod)
	!		xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
	!		hsat = nodes%htemp(inod)+xsat0-xlay
	!
	!		this%idnode_hsattemp(ic)=inod
	!		this%hsattemp(ic)		= hsat
	!		this%xhsattemp(ic)	= xlay+hsat
	!		else !Saturation just below top surface
	!		xsat0 = nodes%x(inod-1)
	!		xsat1 = nodes%x(inod)
	!		xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
	!		p0 = nodes%htemp(inod-1)
	!		p1 = nodes%htemp(inod)
	!		hsat0 = p0+xsat0-xlay
	!		hsat1 = p1+xsat1-xlay
	!		hsat = 0.5_dpd*(hsat1+hsat0)
	!		psat = 0.0_dpd
	!
	!		this%idnode_hsattemp(ic)=inod-1
	!		this%hsattemp(ic)		= hsat
	!		this%xhsattemp(ic)	= xsat0+(xsat1-xsat0)*(psat-p0)/(p1-p0)
	!		end if
	!	else
	!		xsat0 = nodes%x(inod-1)
	!		xsat1 = nodes%x(inod)
	!		xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
	!		p0 = nodes%htemp(inod-1)
	!		p1 = nodes%htemp(inod)
	!		hsat0 = p0+xsat0-xlay
	!		hsat1 = p1+xsat1-xlay
	!		hsat = 0.5_dpd*(hsat1+hsat0)
	!		psat = 0.0_dpd
	!
	!		this%idnode_hsattemp(ic)=inod-1
	!		this%hsattemp(ic)		= hsat
	!		this%xhsattemp(ic)	= xsat0+(xsat1-xsat0)*(psat-p0)/(p1-p0)
	!	end if
	!end do

	!For hold:
	do ic=1,this%count
		this%idelem_hsatold(ic) = this%idnode_layer_bottom(ic)%p
		do inod=this%idnode_layer_bottom(ic)%p,nodes%count
			if(nodes%hold(inod)<0.0_dpd) exit
		end do
		if(inod==this%idnode_layer_bottom(ic)%p) then !Not saturated
			this%idelem_hsatold(ic)=inod
			this%hsatold(ic)		= 0.0_dpd
			this%xhsatold(ic)	= nodes%x(inod)
		else if (inod>=nodes%count) then
			inod=nodes%count
			if(nodes%hold(inod)>0.0_dpd) then !Saturation surpasses the top surface
				xsat0 = nodes%x(inod)
				xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
				hsat = nodes%hold(inod)+xsat0-xlay

				this%idelem_hsatold(ic)=inod
				this%hsatold(ic)		= hsat
				this%xhsatold(ic)	= xlay+hsat
			else !Saturation just below top surface
				xsat0 = nodes%x(inod-1)
				xsat1 = nodes%x(inod)
				xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
				p0 = nodes%hold(inod-1)
				p1 = nodes%hold(inod)
				hsat0 = p0+xsat0-xlay
				hsat1 = p1+xsat1-xlay
				hsat = 0.5_dpd*(hsat1+hsat0)
				psat = 0.0_dpd

				this%idelem_hsatold(ic)=inod-1
				this%hsatold(ic)		= hsat
				this%xhsatold(ic)	= xsat0+(xsat1-xsat0)*(psat-p0)/(p1-p0)
			end if
		else
			xsat0 = nodes%x(inod-1)
			xsat1 = nodes%x(inod)
			xlay = nodes%x(this%idnode_layer_bottom(ic)%p)
			p0 = nodes%hold(inod-1)
			p1 = nodes%hold(inod)
			hsat0 = p0+xsat0-xlay
			hsat1 = p1+xsat1-xlay
			hsat = 0.5_dpd*(hsat1+hsat0)
			psat = 0.0_dpd

			this%idelem_hsatold(ic)=inod-1
			this%hsatold(ic)		= hsat
			this%xhsatold(ic)	= xsat0+(xsat1-xsat0)*(psat-p0)/(p1-p0)
		end if
	end do

	!!this%idelem_hsat		= this%idnode_hsat
	!do ic=1,size(this%idelem_hsat)
	!this%idelem_top			= min(this%idelem_hsat(ic),this%idelem_layer_top(ic)%p)
	!end do
	!this%idelem_hsatold = this%idnode_hsatold

	end subroutine s_unsat_model_update_hsat_from_nodes




	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_incvoldt(this,calc,dt)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_incvoldt" :: s_unsat_constraints_update_incvoldt
	!DEC$ endif
	use com_mod_fem_intelement, only:intelement_abs_1don_sca
	use unsat_mod_ty_calc, only:ty_unsat_calc

	integer,parameter::INTEGRATION_ORDER=40
	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in)::calc
	real(kind=dpd),intent(in)::dt
	real(kind=dpd)::incvol
	integer::i

	do i=1,this%count
		this%incvoldt(i) = calc%get_thnewold_from_x0_to_x1(this%xhsatold(i),this%xhsat(i),INTEGRATION_ORDER)/dt
	end do

	end subroutine s_unsat_constraints_update_incvoldt

	!!******************************************************************************************************************
	!! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!!******************************************************************************************************************
	!
	!subroutine s_unsat_constraints_update_qsup_from_qmed(this,calc)
	!use unsat_mod_ty_calc, only:ty_unsat_calc
	!use com_mod_ty_elements,only:ty_com_elements
	!use unsat_mod_ty_elements,only:ty_unsat_elements
	!
	!class(ty_unsat_constraints),intent(inout)::this
	!class(ty_unsat_calc),intent(in),target::calc
	!class(ty_com_elements),pointer::elements
	!
	!integer::i
	!
	!elements => calc%elements
	!!This gets the mean flow at element where it is hsat and elements below layer
	!	do i=1,this%count
	!	select type (elements)
	!		type is (ty_unsat_elements)
	!			this%qsup(i) = elements%results_qmed(max(this%idnode_hsat(i),elements%count))
	!		end select
	!	end do
	!
	!end subroutine s_unsat_constraints_update_qsup_from_qmed
	!
	!!******************************************************************************************************************
	!! Sub: s_unsat_constraints_update_qsup_from_qsat(calc):	Returns the water flow just in the watertable.
	!!******************************************************************************************************************
	!
	!subroutine s_unsat_constraints_update_qsup_from_qsat(this,calc)
	!use unsat_mod_ty_calc, only:ty_unsat_calc
	!use com_mod_ty_elements,only:ty_com_elements
	!use unsat_mod_ty_elements,only:ty_unsat_elements
	!
	!class(ty_unsat_constraints),intent(inout)::this
	!class(ty_unsat_calc),intent(in),target::calc
	!class(ty_com_elements),pointer::elements
	!real(kind=dpd)::ksat
	!integer::i,idelement
	!
	!elements => calc%elements
	!!This gets the mean flow at element where it is hsat and elements below layer
	!	do i=1,this%count
	!	select type (elements)
	!	type is (ty_unsat_elements)
	!		if(this%idnode_hsat(i)>elements%count) then
	!			ksat = 1.0_dpd
	!			this%qsup(i) = calc%nodes%qnewmann(calc%nodes%count)
	!		else
	!			idelement = this%idnode_hsat(i)
	!			this%qsup(i) = elements%material(idelement)%ksat*elements%results_dhxdxmed(idelement)
	!		end if
	!		end select
	!	end do
	!
	!end subroutine s_unsat_constraints_update_qsup_from_qsat
	!
	!!******************************************************************************************************************
	!! Sub: s_unsat_constraints_update_qsup_from_qsat(calc):	Returns the water flow just in the watertable.
	!!******************************************************************************************************************
	!
	!subroutine s_unsat_constraints_update_qsup_from_qtopinf(this,calc)
	!use unsat_mod_ty_calc, only:ty_unsat_calc
	!use com_mod_ty_elements,only:ty_com_elements
	!use unsat_mod_ty_elements,only:ty_unsat_elements
	!
	!class(ty_unsat_constraints),intent(inout)::this
	!class(ty_unsat_calc),intent(in),target::calc
	!class(ty_com_elements),pointer::elements
	!real(kind=dpd)::ksat,xpercent
	!integer::i,idelement
	!
	!elements => calc%elements
	!!This gets the mean flow at element where it is hsat and elements below layer
	!	do i=1,this%count
	!	select type (elements)
	!	type is (ty_unsat_elements)
	!		if(this%idnode_hsat(i)>elements%count) then
	!			ksat = 1.0_dpd
	!			this%qsup(i) = calc%nodes%qnewmann(calc%nodes%count)
	!		else
	!			idelement = this%idnode_hsat(i)
	!			xpercent = (this%idnode_hsat(i)-elements%h0(idelement))/(elements%h1(idelement)-elements%h0(idelement))
	!			this%qsup(i) =  xpercent*elements%results_qsal(idelement)+(1-xpercent)*elements%results_qent(idelement)
	!		end if
	!		end select
	!	end do
	!
	!end subroutine s_unsat_constraints_update_qsup_from_qtopinf

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qsat_from_qmed(this,calc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qsat_from_qmed" :: s_unsat_constraints_update_qsat_from_qmed
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in),target::calc
	class(ty_com_elements),pointer::elementcom
	class(ty_unsat_elements),pointer::elements

	integer::i

	elementcom => calc%elements
	select type (elementcom)
	type is (ty_unsat_elements)
		elements => elementcom
	end select


	!This gets the mean flow at element where it is hsat and elements below layer
	!do i=1,this%count
	this%qsat = elements%results_qmed(max(this%idelem_hsat,elements%count))
	!this%qsat = elements%results_qmed(this%idelem_hsat)
	!end do

	end subroutine s_unsat_constraints_update_qsat_from_qmed

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qsat_from_qtop(this,calc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qsat_from_qtop" :: s_unsat_constraints_update_qsat_from_qtop
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in),target::calc
	class(ty_com_elements),pointer::elementcom
	class(ty_unsat_elements),pointer::elements

	integer::i

	elementcom => calc%elements
	select type (elementcom)
	type is (ty_unsat_elements)
		elements => elementcom
	end select


	!This gets the mean flow at element where it is hsat and elements below layer
	!do i=1,this%count
	this%qsat = elements%results_qent(this%idelem_hsat)
	!end do

	end subroutine s_unsat_constraints_update_qsat_from_qtop

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qtop_from_qmed(this,calc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qtop_from_qmed" :: s_unsat_constraints_update_qtop_from_qmed
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in),target::calc
	class(ty_com_elements),pointer::elementcom
	class(ty_unsat_elements),pointer::elements
	!integer::i

	elementcom => calc%elements
	select type (elementcom)
	type is (ty_unsat_elements)
		elements => elementcom
	end select

	!This gets the mean flow at element where it is hsat and elements below layer
	!do i=1,this%count
	this%qtop = elements%results_qmed(this%idelem_top)
	!end do

	end subroutine s_unsat_constraints_update_qtop_from_qmed

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qtop_from_qtop(this,calc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qtop_from_qtop" :: s_unsat_constraints_update_qtop_from_qtop
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in),target::calc
	class(ty_com_elements),pointer::elementcom
	class(ty_unsat_elements),pointer::elements
	!integer::i

	elementcom => calc%elements
	select type (elementcom)
	type is (ty_unsat_elements)
		elements => elementcom
	end select

	!This gets the mean flow at element where it is hsat and elements below layer
	!do i=1,this%count
	this%qtop = elements%results_qent(this%idelem_top)
	!end do

	end subroutine s_unsat_constraints_update_qtop_from_qtop

	!******************************************************************************************************************
	! Sub: s_unsat_constraints_update_qsup_from_qsat(calc):	Returns the water flow just in the watertable.
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qsat_from_satpoint(this,calc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qsat_from_satpoint" :: s_unsat_constraints_update_qsat_from_satpoint
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements
	use com_mod_fem_shapefunctions,only:interp_on_element

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in),target::calc
	class(ty_com_elements),pointer::elementcom
	class(ty_unsat_elements),pointer::elements
	real(kind=dpd)::ksat
	integer::i,idelement

	elementcom => calc%elements
	select type (elementcom)
	type is (ty_unsat_elements)
		elements => elementcom
	end select

	!This gets the mean flow at element where it is hsat and elements below layer

	do i=1,this%count
		idelement = this%idelem_hsat(i)
		this%qsat(i) = interp_on_element(this%chi_hsat(i),(/elements%results_qsal(idelement),elements%results_qent(idelement)/))
	end do


	end subroutine s_unsat_constraints_update_qsat_from_satpoint



	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qinf_from_qmed(this,calc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qinf_from_qmed" :: s_unsat_constraints_update_qinf_from_qmed
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in),target::calc
	class(ty_com_elements),pointer::elements

	integer::i

	elements => calc%elements
	!This gets the mean flow at element where it is hsat and elements below layer
	do i=1,this%count
		select type (elements)
		type is (ty_unsat_elements)
			if ((this%idnode_layer_bottom(i)%p-1)==0) then
				this%qinf(i)=0.0_dpd
			else
				this%qinf(i) = elements%results_qmed(this%idnode_layer_bottom(i)%p-1) !On the top of the element below
			end if
		end select
	end do

	end subroutine s_unsat_constraints_update_qinf_from_qmed

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qinf_from_qtop(this,calc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qinf_from_qtop" :: s_unsat_constraints_update_qinf_from_qtop
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in),target::calc
	class(ty_com_elements),pointer::elements

	integer::i

	elements => calc%elements
	!This gets the mean flow at element where it is hsat and elements below layer
	do i=1,this%count
		select type (elements)
		type is (ty_unsat_elements)
			if ((this%idnode_layer_bottom(i)%p-1)==0) then
				this%qinf(i)=0.0_dpd
			else
				this%qinf(i) = elements%results_qent(this%idnode_layer_bottom(i)%p-1) !On the top of the element below
			end if
		end select
	end do

	end subroutine s_unsat_constraints_update_qinf_from_qtop

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qsupinf(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qsupinf" :: s_unsat_constraints_update_qsupinf
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc

	class(ty_unsat_constraints),intent(inout)::this

	integer::i

	do i=1,this%count
		this%qsupinf(i) = this%qsat(i)-this%qinf(i)
		this%qver(i) = this%qsupinf(i)
	end do

	end subroutine s_unsat_constraints_update_qsupinf

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraints_update_qnewman_all(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraints_update_qnewman_all" :: s_unsat_constraints_update_qnewman_all
	!DEC$ endif
	use unsat_mod_ty_calc, only:ty_unsat_calc

	class(ty_unsat_constraints),intent(inout)::this

	integer::i,j

	do i=1,this%count
		this%qnewmann_all(i)=0.0_dpd
		do j=1,this%count
			if(this%idelem_layer_bottom(j)%p<=this%idelem_hsat(i)) this%qnewmann_all(i) = this%qnewmann_all(i)+this%qnewmann(j)%p
		end do
	end do

	end subroutine s_unsat_constraints_update_qnewman_all

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_model_updateall(this,calc,dt)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_model_updateall" :: s_unsat_model_updateall
	!DEC$ endif
	use unsat_mod_ty_calc,only:ty_unsat_calc
	use com_mod_ty_time,only:ty_com_time

	class(ty_unsat_constraints),intent(inout)::this
	class(ty_unsat_calc),intent(in),target::calc
	real(kind=dpd),intent(in)::dt

	integer::i

	!call this%update_hsat(calc%nodes)
	call this%update_hsat(calc%elements)
	call this%update_incvoldt(calc,dt)
	call this%update_qsat(calc)
	call this%update_qinf(calc)
	call this%update_qsupinf()
	call this%update_qnewmann_all()

	end subroutine s_unsat_model_updateall

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_unsat_constraint_set_newmann_from_dqhordx(this,nodes)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_constraint_set_newmann_from_dqhordx" :: s_unsat_constraint_set_newmann_from_dqhordx
	!DEC$ endif
	use com_mod_ty_nodes,only:ty_com_nodes
	!CHECK: set_newmann_from_dqhordx: This procedure put the boundary only in the bottom node of the layer, it would be better to distribute over
	!idnode to idnode_hsat

	class(ty_unsat_constraints),intent(in)::this
	class(ty_com_nodes),intent(inout)::nodes

	integer::l

	do l=1,this%count
		if ((this%hsat(l)>0.00_dpd).and.(this%dqhordx(l)%p>0.0_dpd)) then !CHECK: Check this if I have to put hsat>0 or dqhordx>0 as when hsat=0 there has to be no horizontal flow.
			nodes%qnewmann	(this%idnode_layer_bottom(l)%p)	= -this%dqhordx(l)%p
			nodes%isnewmann	(this%idnode_layer_bottom(l)%p)	= .true.
		else
			nodes%qnewmann	(this%idnode_layer_bottom(l)%p)	= 0.0_dpd
			nodes%isnewmann	(this%idnode_layer_bottom(l)%p)	= .false.
		end if
	end do

	end subroutine s_unsat_constraint_set_newmann_from_dqhordx


	end module unsat_mod_ty_constraints