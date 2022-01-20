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
	!> Iván Campos-Guereta Díez
	!  MSc Civil Engineering by Polytechnic University of Madrid                                                     *
	!  PhD Student by University of Nottingham                                                                       *
	!  eMBA by International Institute San Telmo in Seville                                                          *
	!  ivan.camposguereta@nottingham.ac.uk
	! DESCRIPTION:
	!> Class for the collection of nodes-classes in the saturated model
	!********************************************************************************************************************

	module com_mod_ty_nodes
	use com_mod_ty_material, only: ty_com_material
	!use com_mod_ty_calc, only: ty_com_calc

	implicit none
	include 'inc_precision.fi'

	private

	type,abstract,public::ty_com_nodes
		integer::count !<Total number of nodes-classes
		integer::nn		 !<Number of nodes per element
		integer::nc		 !<Number of classes
		!type(ty_com_calc),allocatable::calc

		integer,allocatable::id(:) !<ID of node-class
		integer,allocatable::idnode(:) !ID of node
		integer,allocatable::ncnode(:) !class number of node
		integer,allocatable::ne(:) !number of elements on node
		integer,allocatable:: idelement(:,:)	!Id of elements on node (count,ne)
		real(kind=dps),allocatable::x(:)			!<Coordinate of node
		real(kind=dps),allocatable::z(:)			!<Z coordinate of node [L]

		real(kind=dps),allocatable::qent(:)		!<Imposed water flow at entrance [L/T]
		!real(kind=dps),allocatable::bound(:)	!<Imposed dirichlet [L]
		real(kind=dps),allocatable::minmeshlenght(:)

		real(kind=dps),allocatable::hinit(:)
		real(kind=dps),allocatable::dhinit(:)

		real(kind=dps),allocatable::hnew(:)
		real(kind=dps),allocatable::htemp(:)
		real(kind=dps),allocatable::hold(:)

		real(kind=dps),allocatable::dhnew(:)
		real(kind=dps),allocatable::dhtemp(:)
		real(kind=dps),allocatable::dhold(:)

		real(kind=dps),allocatable::thnew(:)
		real(kind=dps),allocatable::thtemp(:)
		real(kind=dps),allocatable::thold(:)

		!boundary conditions on node
		logical,allocatable::isdirichlet(:)
		logical,allocatable::isnewmann(:)
		real(kind=dps),allocatable::hdirichlet(:)
		real(kind=dps),allocatable::qnewmann(:)

		type(ty_com_material),allocatable::material(:)

	contains
	procedure,public:: allocateall					=> f_com_nodes_allocateall		!<Initializes nn,nc, and count and allocate all matrix
	procedure,public:: deallocateall				=> f_com_nodes_deallocateall	!<Deallocate all arrays
	procedure,public:: get_id_from_x_bottom	=> f_nodes_get_id_from_x

	end type ty_com_nodes

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_com_nodes. Allocate all vectors and create instance of calc.
	!> @param[in] ne
	!> @param[in] nn
	!> @param[in] nc
	!---------------------------------------------------------------------------------------------------------------------

	subroutine f_com_nodes_allocateall(this,nnc,nn,nc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_com_nodes_allocateall" :: f_com_nodes_allocateall
	!DEC$ endif

	integer,parameter::MAXNELEMENTS_PER_NODE=2

	class(ty_com_nodes),intent(inout)::this
	integer, intent(in)::nn !< Number of nodes per element
	integer, intent(in)::nc !< Class of nodes
	integer, intent(in)::nnc !Total number of nodes and classes

	this%nn = nn
	this%nc = nc
	this%count = nnc

	if(.not.allocated(this%id))					allocate(this%id(nnc))
	if(.not.allocated(this%idnode))			allocate(this%idnode(nnc))
	if(.not.allocated(this%ncnode))			allocate(this%ncnode(nnc))
	if(.not.allocated(this%ne))					allocate(this%ne(nnc))

	if(.not.allocated(this%idelement))	 allocate(this%idelement(nnc,MAXNELEMENTS_PER_NODE))
	if(.not.allocated(this%x))				   allocate(this%x(nnc))
	if(.not.allocated(this%z))				   allocate(this%z(nnc))
	if(.not.allocated(this%qent))				 allocate(this%qent(nnc))
	if(.not.allocated(this%minmeshlenght))allocate(this%minmeshlenght(nnc))
	if(.not.allocated(this%hinit))			 allocate(this%hinit(nnc))
	if(.not.allocated(this%dhinit))			 allocate(this%dhinit(nnc))
	if(.not.allocated(this%hnew))				 allocate(this%hnew(nnc))
	if(.not.allocated(this%htemp))			 allocate(this%htemp(nnc))
	if(.not.allocated(this%hold))				 allocate(this%hold(nnc))
	if(.not.allocated(this%dhnew))			 allocate(this%dhnew(nnc))
	if(.not.allocated(this%dhtemp))			 allocate(this%dhtemp(nnc))
	if(.not.allocated(this%dhold))			 allocate(this%dhold(nnc))
	if(.not.allocated(this%thnew))				 allocate(this%thnew(nnc))
	if(.not.allocated(this%thtemp))			 allocate(this%thtemp(nnc))
	if(.not.allocated(this%thold))				 allocate(this%thold(nnc))

	if(.not.allocated(this%isdirichlet)) allocate(this%isdirichlet(nnc))
	if(.not.allocated(this%isnewmann))	 allocate(this%isnewmann(nnc))
	if(.not.allocated(this%hdirichlet))	 allocate(this%hdirichlet(nnc))
	if(.not.allocated(this%qnewmann))		 allocate(this%qnewmann(nnc))

	if(.not.allocated(this%material))		 allocate(this%material(nnc))

	!call this%calc%init(nnc,sparsity)

	this%isdirichlet = .false.
	this%isnewmann = .false.
	this%hdirichlet = 0.0_dpd
	this%qnewmann = 0.0_dpd

	end subroutine f_com_nodes_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_nodes. Deallocate all vectors and instance of calc.
	!---------------------------------------------------------------------------------------------------------------------


	subroutine f_com_nodes_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_com_nodes_deallocateall" :: f_com_nodes_deallocateall
	!DEC$ endif

	class(ty_com_nodes),intent(inout)::this

	if(allocated(this%id))						deallocate(this%id)
	if(allocated(this%idnode))				deallocate(this%idnode)
	if(allocated(this%ncnode))				deallocate(this%ncnode)
	if(allocated(this%ne))						deallocate(this%ne)

	if(allocated(this%idelement))			deallocate(this%idelement)
	if(allocated(this%x))							deallocate(this%x)
	if(allocated(this%z))							deallocate(this%z)
	if(allocated(this%hinit))					deallocate(this%hinit)
	if(allocated(this%qent))					deallocate(this%qent)
	if(allocated(this%minmeshlenght))	deallocate(this%minmeshlenght)
	if(allocated(this%hnew))					deallocate(this%hnew)
	if(allocated(this%htemp))					deallocate(this%htemp)
	if(allocated(this%hold))					deallocate(this%hold)
	if(allocated(this%thnew))					deallocate(this%thnew)
	if(allocated(this%thtemp))				deallocate(this%thtemp)
	if(allocated(this%thold))					deallocate(this%thold)

	if(allocated(this%isdirichlet))		deallocate(this%isdirichlet)
	if(allocated(this%isnewmann))			deallocate(this%isnewmann)
	if(allocated(this%hdirichlet))		deallocate(this%hdirichlet)
	if(allocated(this%qnewmann))			deallocate(this%qnewmann)

	if(allocated(this%material))			deallocate(this%material)

	!call this%calc%deallocateall()

	end subroutine f_com_nodes_deallocateall


	!---------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns the id of the layer given the height x from the bottom.
	!> @param[in] x
	!---------------------------------------------------------------------------

	function f_nodes_get_id_from_x(this,x)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_nodes_get_id_from_x" :: f_nodes_get_id_from_x
	!DEC$ endif


	class(ty_com_nodes),intent(in)::this
	real(kind=dpd),intent(in)::x(:)
	integer::f_nodes_get_id_from_x(size(x))
	integer::i,findtemp(1)
	logical::mask(this%count)

	DO i=1, size(x)
		mask = this%x<=x(i)
		findtemp = FINDLOC(mask,.true.,BACK=.true.)
		f_nodes_get_id_from_x(i) = findtemp(1)
	end do

	end function f_nodes_get_id_from_x

	end module com_mod_ty_nodes