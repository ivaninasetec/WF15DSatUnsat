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

	module com_mod_ty_calc

	use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded
	use com_mod_ty_nodes, only:ty_com_nodes
	use com_mod_ty_elements, only:ty_com_elements
	use com_mod_ty_parameters, only:ty_com_parameters
	use com_mod_ty_boundary, only:ty_com_boundary
	use com_mod_ty_time, only:ty_com_time


	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_com_calc
		class(ty_com_parameters),pointer::parameters
		class(ty_com_time),pointer::time
		!real(kind=dpd):: t			!<Actual time
		!real(kind=dpd):: told		!<Time at previous step
		!real(kind=dpd):: dt			!<Actual time
		!real(kind=dpd):: dtold	!<Time at previous step
		!real(kind=dpd):: tprint	!<Actual printtime
		!logical				:: CheckPrint	!<If true then print

		real(kind=dpd)::epsh		!<Actual value of iteration error in h	(hnew-htemp)
		real(kind=dpd)::epsth		!<Actual value of iteration error in th (hnew-htemp)

		
		class(ty_com_nodes)			,allocatable::nodes
		class(ty_com_elements)	,allocatable::elements
		type(ty_com_boundary)		,pointer::boundary

		!Right hand side of the system
		real(kind=dps),allocatable::rhs(:)			!<Calculated right hand side vector

		real(kind=dps),allocatable::solution(:)	!<Vector with the solution of the linear system
		real(kind=dps),allocatable::slopeold(:) !<Equal to (hnew-hold)/dt for previous step in order to estimate new h.

		!Left han side fo the system
		class(ty_mx),allocatable::mx			!<Matrix of the FEM linear system

		!Others
		integer				,allocatable::perm(:) !<Permutation vector for a better conditioned system.


	contains
	procedure,public:: allocateall					=> s_com_calc_allocateall !<Initializes nn,nc, and count and allocate all matrix
	procedure,public:: deallocateall				=> s_com_calc_deallocateall !<Initializes nn,nc, and count and allocate all matrix
	procedure,public:: shapematrix					=> s_com_calc_shapematrix
	procedure,public:: allocate_matrix			=> s_com_calc_allocate_other_matrix !<Initializes nn,nc, and count and allocate all matrix
	procedure,public:: buildcoefmatrix			=> s_com_calc_buildcoefmatrix
	procedure,public:: build_linearsystem		=>	s_com_buildlinearsystem
	procedure,public:: iterate							=>	s_com_calc_iteration
	procedure,public:: set_to_initial				=>	s_com_calc_set_to_initial
	procedure,public:: revert_to_old				=>	s_com_calc_revert_to_old
	procedure,public:: set_old							=>	s_com_calc_set_old
	procedure,public:: estimate_hnew_for_new_timestep	=> s_com_calc_estimate_hnew_for_new_timestep
	procedure,public:: estimate_htemp_for_new_iteration	=> s_com_calc_estimate_htemp_for_new_iteration
	procedure,public:: update_hnew_from_solution	=> s_com_calc_update_hnew_from_solution
	procedure,public:: assign_dirichlet			=> s_vector_assign_dirichlet_from_nodes
	procedure,public:: set_dirichlet_to_node 	=> s_com_calc_set_dirichlet_to_node_sca,s_com_calc_set_dirichlet_to_node_vec
	procedure,public:: set_newmann_to_node 	=> s_com_calc_set_newmann_to_node_sca,s_com_calc_set_newmann_to_node_vec
	procedure,public:: get_idelement_from_x 	=> f_com_calc_get_idelement_from_x
	
	end type ty_com_calc

	interface
	subroutine s_com_buildlinearsystem(this,IsTimeDependant,option)
	import::ty_com_calc,dpd
	class(ty_com_calc),intent(inout),target::this
	!real(kind=dpd),intent(in)::dt
	logical,intent(in)::IsTimeDependant
	integer,intent(in),optional::option
	end subroutine s_com_buildlinearsystem
	end interface

	contains
	
	!---------------------------------------------------------------------------------------------------------------------
	! Set Dirichlet to node
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_model. Build coefficient matrix for actual head pressures.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------
	
	subroutine s_com_calc_set_dirichlet_to_node_sca(this,idnode,hdirichlet)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_set_dirichlet_to_node_sca" :: s_com_calc_set_dirichlet_to_node_sca
	!DEC$ endif

	class(ty_com_calc),intent(inout)::this
	integer,intent(in)::idnode
	real(kind=dpd)::hdirichlet
	
	this%nodes%isdirichlet(idnode) = .true.
	this%nodes%hdirichlet(idnode) = hdirichlet
		this%nodes%isnewmann(idnode) = .false.

	end subroutine s_com_calc_set_dirichlet_to_node_sca
	
	
	subroutine s_com_calc_set_dirichlet_to_node_vec(this,idnode,hdirichlet)
		!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_set_dirichlet_to_node_vec" :: s_com_calc_set_dirichlet_to_node_vec
	!DEC$ endif

	class(ty_com_calc),intent(inout)::this
	integer,intent(in)::idnode(:)
	real(kind=dpd)::hdirichlet(:)
	
	this%nodes%isdirichlet(idnode) = .true.
	this%nodes%hdirichlet(idnode) = hdirichlet
		this%nodes%isnewmann(idnode) = .false.

	end subroutine s_com_calc_set_dirichlet_to_node_vec

	!---------------------------------------------------------------------------------------------------------------------
	! Set Dirichlet to node
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_model. Build coefficient matrix for actual head pressures.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------
	
	subroutine s_com_calc_set_newmann_to_node_sca(this,idnode,qnewmann)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_set_newmann_to_node_sca" :: s_com_calc_set_newmann_to_node_sca
	!DEC$ endif

	class(ty_com_calc),intent(inout)::this
	integer,intent(in)::idnode
	real(kind=dpd)::qnewmann
	
	this%nodes%isnewmann(idnode)	= .true.
	this%nodes%qnewmann(idnode)		= qnewmann
	this%nodes%isdirichlet(idnode) = .false.

	end subroutine s_com_calc_set_newmann_to_node_sca
	
	subroutine s_com_calc_set_newmann_to_node_vec(this,idnode,qnewmann)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_set_newmann_to_node_vec" :: s_com_calc_set_newmann_to_node_vec
	!DEC$ endif

	class(ty_com_calc),intent(inout)::this
	integer,intent(in)::idnode(:)
	real(kind=dpd)::		qnewmann(:)
	
	this%nodes%isnewmann(idnode)		= .true.
	this%nodes%qnewmann(idnode)			= qnewmann
	this%nodes%isdirichlet(idnode)	= .false.

	end subroutine s_com_calc_set_newmann_to_node_vec	
	

	!---------------------------------------------------------------------------------------------------------------------
	! ALLOCATE ALL
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_model. Build coefficient matrix for actual head pressures.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_com_calc_allocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_allocateall" :: s_com_calc_allocateall
	!DEC$ endif

	!integer,parameter::KU=1,KL=1

	class(ty_com_calc),intent(inout)::this
	integer::nnc !<Number of nodes-classes

	nnc=this%nodes%count

	!Allocate all vectors
	if(.not.allocated(this%slopeold)) allocate(this%slopeold(nnc))
	if(.not.allocated(this%rhs))			allocate(this%rhs(nnc))
	if(.not.allocated(this%solution)) allocate(this%solution(nnc))
	!if(.not.allocated(this%perm))			allocate(this%perm(nnc)) !If perm is not allocated then it is not used.

	!allocate matrix
	select case(this%parameters%typematrixstorage)
	case(1) !dense matrix-------------------------------------------------------------------------------------------------
		allocate(ty_mx_dense::this%mx)
	case(2) !csr matrix-------------------------------------------------------------------------------------------------
		allocate(ty_mx_csr::this%mx)
	case(3) !banded matrix-------------------------------------------------------------------------------------------------
		allocate(ty_mx_banded::this%mx)
	end select
	call this%shapematrix()

	end subroutine s_com_calc_allocateall


	!---------------------------------------------------------------------------------------------------------------------
	! DEALLOCATE ALL
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_model. Build coefficient matrix for actual head pressures.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_deallocateall" :: s_com_calc_deallocateall
	!DEC$ endif

	class(ty_com_calc),intent(inout)::this

	!if(allocated(this%colbound))	deallocate(this%colbound)
	!if(allocated(this%colqent))	deallocate(this%colqent)
	if(allocated(this%slopeold))	deallocate(this%slopeold)
	if(allocated(this%rhs))				deallocate(this%rhs)
	if(allocated(this%solution))	deallocate(this%solution)
	if(allocated(this%perm))			deallocate(this%perm)

	if(allocated(this%mx))				deallocate(this%mx)

	end subroutine s_com_calc_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_nodes. Shape the matrix from the distribution of nodes and elements.
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_shapematrix(this)
		!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_shapematrix" :: s_com_calc_shapematrix
	!DEC$ endif
	use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded

	class(ty_com_calc),target,intent(inout)::this

	class(ty_mx),pointer::mx,mxmass,mxstiff,mxload,mxbound
	integer::i,j,ic,jc,n,m,e,ncount,kl,ku
	integer::nelem,numnod,numclass,numnodclass,numnp,nstartelem,nendelem

	nelem = this%elements%count
	numnod=	this%elements%nn !Nodes on element
	numclass = this%elements%nc+1 !Classes on element
	numnodclass = this%elements%nnc !Nodes and classes on element
	kl = numclass
	ku = numclass

	numnp=this%nodes%count !Total number of nodes

	!Pointers to each array
	mx			=>this%mx
	!mxd			=>this%mxd
	!mxmass	=>this%mxmass
	!mxstiff	=>this%mxstiff
	!mxload	=>this%mxload
	!mxbound =>this%mxbound

	!initialize matrix.....
	select type(mx)
	type is(ty_mx_dense) !dense matrix---------------------------------------------------------------------------------------

		!initiallize...
		call mx       %init(numnp,numnp)
		!call mxmass   %init(numnp,numnp)
		!call mxstiff  %init(numnp,numnp)
		!call mxload   %init(numnp,numnp)
		!call mxbound  %init(numnp,numnp)

	type is(ty_mx_csr) !csr matrix---------------------------------------------------------------------------------------------------

		!calculate number of values...
		if (.not.allocated(mx%pos)) allocate(mx%pos(numnp,numnp))

		mx%pos = 0.0_dpd
		ncount = 0

		!first construct pos(r,c) matrix and ncount
		do e=1,this%elements%count
			nstartelem = this%elements%idnode(e,1) !id of nodeclass that start element
			nendelem = this%elements%idnode(e,numnodclass) !id of nodeclass that end element

			do n=nstartelem,nendelem !rows
				do m=nstartelem,nendelem !columns
					if(mx%pos(n,m)==0) then
						ncount=ncount+1
						mx%pos(n,m) = ncount
					end if
				end do
			end do
		end do

		!initiallize csr matrix...
		call mx       %init(ncount,numnp)
		!call mxmass   %init(ncount,numnp)
		!call mxstiff  %init(ncount,numnp)
		!call mxload   %init(ncount,numnp)
		!call mxbound  %init(ncount,numnp)

		!fills csrows, csrcols and this%pos...
		ncount=0

		mx%csrrows(1) = 1 !first element of csrrows
		do n=1,numnp
			do m=1,numnp
				if (mx%pos(n,m).ne.0) then
					ncount=ncount+1
					!select type(mx)
					!	type is (ty_mx_csr)
					mx%csrcols(ncount) = m
					!end select
					mx%pos(n,m)=ncount
				end if
			end do
			mx%csrrows(n+1) = ncount+1
		end do



	type is(ty_mx_banded) !banded matrix------------------------------------------------------------------------------------------------

		!initialize matrix...
		call mx      %init(numnp,kl,ku)
		!call mxmass  %init(numnp,kl,ku)
		!call mxstiff %init(numnp,kl,ku)
		!call mxload  %init(numnp,kl,ku)
		!call mxbound %init(numnp,kl,ku)

	end select

	!!update csrrows and csrcols for all csr matrixs in the model...
	!
	!select type(mx)
	!type is(ty_mx_csr)
	!	select type(mxmass)
	!	type is(ty_mx_csr)
	!		mxmass%csrrows  = mx%csrrows
	!		mxmass%csrcols  = mx%csrcols
	!		mxmass%pos  = mx%pos
	!	end select
	!
	!	select type(mxstiff)
	!	type is(ty_mx_csr)
	!		mxstiff%csrrows  = mx%csrrows
	!		mxstiff%csrcols  = mx%csrcols
	!		mxstiff%pos  = mx%pos
	!	end select
	!
	!	select type(mxload)
	!	type is(ty_mx_csr)
	!		mxload%csrrows  = mx%csrrows
	!		mxload%csrcols  = mx%csrcols
	!		mxload%pos  = mx%pos
	!	end select
	!
	!	select type(mxbound)
	!	type is(ty_mx_csr)
	!		mxbound%csrrows  = mx%csrrows
	!		mxbound%csrcols  = mx%csrcols
	!		mxbound%pos  = mx%pos
	!	end select
	!
	!end select


	end subroutine s_com_calc_shapematrix

	!---------------------------------------------------------------------------------------------------------------------
	!	ALLOCATE OTHER MATRIX
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_nodes. Shape the matrix from the distribution of nodes and elements.
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_allocate_other_matrix(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_allocate_other_matrix" :: s_com_calc_allocate_other_matrix
	!DEC$ endif
	use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded

	class(ty_com_calc),target,intent(inout)::this
	class(ty_mx),allocatable,intent(inout)::mx

	class(ty_mx),pointer::mxbase
	integer::ncount,numnp,kl,ku

	mxbase => this%mx

	!!allocate matrix.....
	select type(mxbase)
	type is(ty_mx_dense) !dense matrix---------------------------------------------------------------------------------------
		if(.not.allocated(mx)) allocate(ty_mx_dense::mx)
	type is(ty_mx_csr) !dense matrix---------------------------------------------------------------------------------------
		if(.not.allocated(mx)) allocate(ty_mx_csr::mx)
	type is(ty_mx_banded) !dense matrix---------------------------------------------------------------------------------------
		if(.not.allocated(mx)) allocate(ty_mx_banded::mx)
	end select

	!Get values of number of nodes, rows and columns
	select type(mxbase)
	type is(ty_mx_dense)
		numnp = size(mxbase%mx,1)
	type is(ty_mx_csr)
		numnp = size(mxbase%csrrows-1)
		ncount = size(mxbase%csrcols)
	type is(ty_mx_banded)
		kl = mxbase%kl
		ku = mxbase%ku
		numnp = size(mxbase%mx,2)
	end select

	!initialize matrix.....
	select type(mx)
	type is(ty_mx_dense)
		call mx       %init(numnp,numnp)
	type is(ty_mx_csr)
		call mx       %init(ncount,numnp)
	type is(ty_mx_banded)
		call mx       %init(numnp,kl,ku)
	end select

	!update csrrows and csrcols for all csr matrixs in the model...
	select type(mxbase)
	type is(ty_mx_csr)
		select type(mx)
		type is(ty_mx_csr)
			mx%csrrows		= mxbase%csrrows
			mx%csrcols		= mxbase%csrcols
			mx%pos				= mxbase%pos
		end select
	end select

	end subroutine s_com_calc_allocate_other_matrix

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> This procedure build the coefficient matrix mx given the function f(h).
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_buildcoefmatrix(this,fun_quad,mx,masslump,optionbasis)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_buildcoefmatrix" :: s_com_calc_buildcoefmatrix
	!DEC$ endif
	use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded
	use com_mod_fem_intelement,only:intelement_rel_1don

	interface
	function fun_quad(chi,e)
	import::dps
	real(kind=dps),intent(in)::chi(:)
	integer,intent(in)::e
	real(kind=dps)::fun_quad(size(chi))
	end function
	end interface

	integer,parameter::COEF_MASS=1,COEF_STIFF=2,COEF_LOAD=3,COEF_BOUND=4
	integer,parameter::OPTIONBASIS_H=1,OPTIONBASIS_DH=2,OPTIONBASIS_H_H=3,OPTIONBASIS_H_DH=4,OPTIONBASIS_DH_H=5,OPTIONBASIS_DH_DH=6

	class(ty_com_calc),intent(in)::this


	!procedure(auxfun2)::fun_quad !Function to integrate, with arguments chi and the element number e.

	class(ty_mx),intent(inout)::mx
	integer,intent(in)::optionbasis
	logical,intent(in)::masslump
	real(kind=dpd)::henodes(this%elements%nn*(this%elements%nc+1))
	real(kind=dpd)::henodesold(this%elements%nn*(this%elements%nc+1))
	real(kind=dpd)::hlayer

	integer::n,m,i,e,ishape,jshape
	integer::nelem,numnod,numclass,numnodclass,nstartelem,nendelem

	nelem=this%elements%count
	numnod=this%elements%nn
	numclass = this%elements%nc+1
	numnodclass = numnod*numclass
	
	mx=0.0_dpd

	do e=1,nelem
		nstartelem	= this%elements%idnode(e,1)
		nendelem		= this%elements%idnode(e,numnodclass)

		do n=nstartelem,nendelem
			ishape = n-nstartelem+1
			do m=nstartelem,nendelem
				jshape = m-nstartelem+1
				if (masslump) then
					!This is the integral on the element of the function auxfun(chi)
					!mx(n,n)=mx(n,n)+int(f,xnodonelem,no)
					call mx%set(n,n,mx%get(n,n)+intelement_rel_1don(auxfun,this%nodes%x(this%elements%idnode(e,:)),this%parameters%quadratureorder))
				else
					call mx%set(n,m,mx%get(n,m)+intelement_rel_1don(auxfun,this%nodes%x(this%elements%idnode(e,:)),this%parameters%quadratureorder))
				end if
			end do
		end do
	end do

	contains
	!This function returns the value to integrate in element e. Given relative coord chi and element number e


	function auxfun(chi)
	!use com_mod_fem_shapefunctions,only:shape1d,dshape1d
	use com_mod_fem_shapefunctions,only:shape1d,dphi1d
	use com_mod_fem_shapefunctions,only:jacobian1d
	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::auxfun(size(chi)),xne(numnodclass)
	!real(kind=dps)::shapefunresult(size(chi),numnodclass),dshapefunresult(size(chi),numnodclass)
	real(kind=dps)::shapefunresult(size(chi),numnodclass),dphifunresult(size(chi),numnodclass)
	integer::nstartelem,nendelem

	select case(optionbasis)
	case(OPTIONBASIS_H) !f(chi,e)·basis_i(chi)
		shapefunresult = shape1d(chi,this%elements%chi(e,:))
		auxfun = fun_quad(chi,e)*shapefunresult(:,ishape)
	case(OPTIONBASIS_DH) !f(chi,e)·dbasis_i(chi)
		nstartelem	= this%elements%idnode(e,1)
		nendelem		= this%elements%idnode(e,numnodclass)
		xne = this%nodes%x(nstartelem : nendelem)
		!dshapefunresult = dshape1d(chi,this%elements%chi(e,:))
		!auxfun = fun_quad(chi,e)*dshapefunresult(:,ishape)/Jacobian1D(chi,xne)
		dphifunresult = dphi1d(chi,xne)
		auxfun = fun_quad(chi,e)*dphifunresult(:,ishape)
	case(OPTIONBASIS_H_H) !f(chi,e)·basis_i(chi)·basis_j(chi)
		shapefunresult = shape1d(chi,this%elements%chi(e,:))
		auxfun = fun_quad(chi,e)*shapefunresult(:,ishape)*shapefunresult(:,jshape)
	case(OPTIONBASIS_H_DH) !f(chi,e)·basis_i(chi)·dbasis_j(chi)
		nstartelem	= this%elements%idnode(e,1)
		nendelem		= this%elements%idnode(e,numnodclass)
		xne = this%nodes%x(nstartelem : nendelem)
		shapefunresult = shape1d(chi,this%elements%chi(e,:))
		!dshapefunresult = dshape1d(chi,this%elements%chi(e,:))
		!auxfun = fun_quad(chi,e) * shapefunresult(:,ishape) * dshapefunresult(:,jshape)/Jacobian1D(chi,xne)
		dphifunresult = dphi1d(chi,xne)
		auxfun = fun_quad(chi,e)*shapefunresult(:,ishape)*dphifunresult(:,jshape)
	case(OPTIONBASIS_DH_H) !f(chi,e)·dbasis_i(chi)·basis_j(chi)
		nstartelem	= this%elements%idnode(e,1)
		nendelem		= this%elements%idnode(e,numnodclass)
		xne = this%nodes%x(nstartelem : nendelem)
		shapefunresult = shape1d(chi,this%elements%chi(e,:))
		!dshapefunresult = dshape1d(chi,this%elements%chi(e,:))
		!auxfun = fun_quad(chi,e) * dshapefunresult(:,ishape)/Jacobian1D(chi,xne)	 * shapefunresult(:,jshape)
		dphifunresult = dphi1d(chi,xne)
		auxfun = fun_quad(chi,e)*dphifunresult(:,ishape)*shapefunresult(:,jshape)
	case(OPTIONBASIS_DH_DH)	!f(chi,e)·dbasis_i(chi)·dbasis_j(chi)
		nstartelem	= this%elements%idnode(e,1)
		nendelem		= this%elements%idnode(e,numnodclass)
		xne = this%nodes%x(nstartelem : nendelem)
		!dshapefunresult = dshape1d(chi,this%elements%chi(e,:))
		!auxfun = fun_quad(chi,e) * dshapefunresult(:,ishape)/Jacobian1D(chi,xne) * dshapefunresult(:,jshape)/Jacobian1D(chi,xne)
		dphifunresult = dphi1d(chi,xne)
		auxfun = fun_quad(chi,e)*dphifunresult(:,ishape)*dphifunresult(:,jshape)
	end select

	end function auxfun

	end subroutine s_com_calc_buildcoefmatrix

	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_estimate_htemp_for_new_iteration(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_estimate_htemp_for_new_iteration" :: s_com_calc_estimate_htemp_for_new_iteration
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this
	
	this%nodes%htemp = this%nodes%hnew
	
	!this%nodes%hnew	= this%nodes%hold+this%slopeold*this%dt
	end subroutine s_com_calc_estimate_htemp_for_new_iteration

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_iteration(this,IsConverged)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_iteration" :: s_com_calc_iteration
	!DEC$ endif
	use com_mod_mx_solvers,only: solve_linalg
	use com_mod_ty_layers,only:ty_com_layers

	logical,parameter::IS_TIME_DEPENDANT=.true.,IS_NOT_TIME_DEPENDANT=.false.

	INTEGER,PARAMETER::COEF_MASS=1,COEF_STIFF=2,COEF_LOAD=3,COEF_BOUND=4
	INTEGER,PARAMETER::USE_H_TEMP=1,USE_H_NEW=0

	class(ty_com_calc),intent(inout)::this
	!real(kind=dpd),intent(in)::dt
	logical,intent(inout)::IsConverged

	integer::numnp
	real(kind=dps)::ksal

	numnp = this%nodes%count

	call this%estimate_htemp_for_new_iteration()
	!this%nodes%htemp = this%nodes%hnew

	!build linear system
	call this%build_linearsystem(IS_TIME_DEPENDANT,USE_H_TEMP)

	!Assign Dirichlet conditions
	call this%assign_dirichlet()
	
	! solve linear algebra system
	call solve_linalg(this%mx,this%rhs,this%solution,this%parameters%typesolver,this%perm)

	!update hnew wiith overrelaxation: (hnew = htemp+crelax·(hnew-htemp)
	call this%update_hnew_from_solution()

	!calculate error in h (in this case error in nodes):
	!this%epsh = maxval(abs(this%solution-this%nodes%htemp))
	
	this%epsh = maxval(abs(this%nodes%hnew-this%nodes%htemp)) !Chechk (If we are updating hnew to only positive, this is the one that has to be put.

	isconverged = this%epsh<this%parameters%epsh_tol
	


	end subroutine s_com_calc_iteration

	!---------------------------------------------------------------------------------------------------------------------
	! UPDATE HNEW FROM SOLUTION AND GET EPSH
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_update_hnew_from_solution(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_update_hnew_from_solution" :: s_com_calc_update_hnew_from_solution
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this

	this%nodes%hnew = this%nodes%htemp + this%parameters%crelax*(this%solution-this%nodes%htemp)

	end subroutine s_com_calc_update_hnew_from_solution

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_set_to_initial(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_set_to_initial" :: s_com_calc_set_to_initial
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this
	this%time%t			= this%parameters%tinit
	this%time%dt			= this%parameters%dtinit
	this%time%told		= this%parameters%tinit
	this%time%dtold	= this%parameters%dtinit
	this%time%checkprint = .true.

	!All values to initial conditions
	this%nodes%hnew		= this%nodes%hinit
	this%nodes%htemp	= this%nodes%hinit
	this%nodes%hold		= this%nodes%hinit
	!All values to initial conditions
	this%nodes%dhnew	= this%nodes%dhinit
	this%nodes%dhtemp	= this%nodes%dhinit
	this%nodes%dhold	= this%nodes%dhinit
	

	end subroutine s_com_calc_set_to_initial

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_revert_to_old(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_revert_to_old" :: s_com_calc_revert_to_old
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this
	this%time%t			= this%time%told
	!call this%time%increase_time() !It is done in the timestepping
	!this%dt			= this%dtold

	!All values to initial conditions
	this%nodes%hnew		= this%nodes%hold
	this%nodes%htemp	= this%nodes%hold

	!All values to initial conditions
	this%nodes%dhnew	= this%nodes%dhold
	this%nodes%dhtemp	= this%nodes%dhold

	end subroutine s_com_calc_revert_to_old





	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_set_old(this)
		!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_set_old" :: s_com_calc_set_old
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this

	!this%slopeold = 0.0_dpd
	
	if(this%time%told.ne.this%time%t) then
	this%time%dtold	= (this%time%t-this%time%told)
	this%time%told	= this%time%t
	!this%time%dtold	= this%time%dt
	end if
	this%slopeold = min(1.0_dpd,(this%nodes%hnew-this%nodes%hold)/this%time%dtold)
	this%nodes%hold = this%nodes%hnew

	end subroutine s_com_calc_set_old

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_estimate_hnew_for_new_timestep(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_estimate_hnew_for_new_timestep" :: s_com_calc_estimate_hnew_for_new_timestep
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this
	
	this%nodes%hnew	= this%nodes%hold+0.5_dpd*this%slopeold*this%time%dt !CHECK
	!this%nodes%hnew	= this%nodes%hold
	
	end subroutine s_com_calc_estimate_hnew_for_new_timestep
	
	!******************************************************************************************************************
	! Sub: s_vector_assign_dirichlet_from_nodes(this,NODS):   Modify matrix and RHS with dirichlet conditions.
	!******************************************************************************************************************		
	
	SUBROUTINE s_vector_assign_dirichlet_from_nodes(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_vector_assign_dirichlet_from_nodes" :: s_vector_assign_dirichlet_from_nodes
	!DEC$ endif
	CLASS(ty_com_calc),INTENT(INOUT)::this
	
	!CLASS(TY_MX),POINTER::MX
	!CLASS(TY_MX),POINTER::MXD
	
	integer::n,m
	real(kind=dps)::hnod

	!MX => THIS%MX
	!MXD => THIS%MXD
	!
	!MXD = MX
	!this%RHSD = this%RHS
 !
 ! SELECT TYPE (MX)
 ! TYPE IS (ty_mx_dense)
 !   NUMNODCLASES = SIZE(MX%MX,1)
 ! TYPE IS (ty_mx_csr)
 !   NUMNODCLASES = SIZE(MX%CSRROWS)-1
 ! TYPE IS (ty_mx_banded)
 !   NUMNODCLASES = SIZE(MX%MX,2)
 ! END SELECT

  DO n = 1, this%nodes%count
    !DO NC = 0, NODS%NODE(nd)%nclass

      !--------
      IF (this%nodes%Isdirichlet(n)) THEN !Assign dirichlet pressure condition to node
        !IDNOD = (nd-1)*(NODS%NODE(nd)%nclass+1)+(nc+1)
        !HNOD = NODS%NODE(nd)%HDirichlet

        DO m=1,this%nodes%count
          IF(this%MX%Get(n,m).NE.0) this%RHS(m) = this%RHS(m)-this%nodes%HDirichlet(n)*this%MX%Get(n,m)
          IF(this%MX%Get(n,m).NE.0) CALL this%MX%Set(n,m,0.0_dpd)
          IF(this%MX%Get(m,n).NE.0) CALL this%MX%Set(m,n,0.0_dpd)
        END DO
        CALL this%MX%Set(n,n,1.0_dpd)
        this%rhs(n) = this%nodes%HDirichlet(n)
      END IF

      !!------------
      !
      !!INCLUDE DIRICHLET CONDITIONS OF WATER FLOWS WHEN DIFFERENTIABILITY OF NODES IS 1
      !IF (NODS%NODE(nd)%IsNewman.AND.nc==1) THEN
      !  IDNOD = (nd-1)*(NODS%NODE(nd)%nclass+1)+(nc+1)
      !  HNOD = NODS%NODE(nd)%QNewman
      !
      !  DO i=1,NUMNODCLASES
      !    IF(MX%Get(IDNOD,i).NE.0) this%RHSD(i) = this%RHSD(i)-HNOD*MX%Get(IDNOD,i)
      !    IF(MX%Get(IDNOD,i).NE.0) CALL MXD%Set(IDNOD,i,0.0_dpd)
      !    IF(MX%Get(i,IDNOD).NE.0) CALL MXD%Set(i,IDNOD,0.0_dpd)
      !  END DO
      !  CALL MXD%Set(IDNOD,IDNOD,1.0_dpd)
      !  this%RHSD(IDNOD) = HNOD
      !END IF

    END DO	
	
	END SUBROUTINE s_vector_assign_dirichlet_from_nodes

	!******************************************************************************************************************
	! Sub: s_vector_assign_dirichlet_from_nodes(this,NODS):   Modify matrix and RHS with dirichlet conditions.
	!******************************************************************************************************************		
	
	function f_com_calc_get_idelement_from_x(this,x)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_com_calc_get_idelement_from_x" :: f_com_calc_get_idelement_from_x
	!DEC$ endif
	class(ty_com_calc),intent(in)::this
	real(kind=dpd),intent(in)::x
	integer::f_com_calc_get_idelement_from_x
	logical::isascending
	integer::idelem(1)
	
	if(this%nodes%x(1)<this%nodes%x(this%nodes%count)) then
		isascending = .true.
	else
		isascending = .false.
	end if
	
	select case(isascending)
	case(.true.)
		if(x<=this%nodes%x(1)) then
			f_com_calc_get_idelement_from_x=1
		else if (x>=this%nodes%x(this%nodes%count)) then
			f_com_calc_get_idelement_from_x=this%nodes%count
		else
			idelem = minloc(x-this%nodes%x,(x-this%nodes%x)>=0.0_dpd)
			f_com_calc_get_idelement_from_x=idelem(1)
		end if
		
	case(.false.)
		if(x>=this%nodes%x(1)) then
			f_com_calc_get_idelement_from_x=1
		else if (x<=this%nodes%x(this%nodes%count)) then
			f_com_calc_get_idelement_from_x=this%nodes%count
		else
			idelem = minloc(this%nodes%x-x,(this%nodes%x-x)>=0.0_dpd)
			f_com_calc_get_idelement_from_x=idelem(1)
		end if
	
	end select	
	
	end function f_com_calc_get_idelement_from_x
	

	end module com_mod_ty_calc
	
	