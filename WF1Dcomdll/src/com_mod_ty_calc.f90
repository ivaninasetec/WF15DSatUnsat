	!********************************************************************************************************************
	! TITLE         : COM_MOD_TY_CALC: DERIVED TYPE THAT INCLUDE NODES AND ELEMENTS AND COMMON METHODS TO PERFORM ON THEM
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_CALC
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Derived type that include nodes and elements and common methods to perform on them
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
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

		real(kind=dpd)::epsh		!<Actual value of iteration error in h	(hnew-htemp)
		real(kind=dpd)::epsth		!<Actual value of iteration error in th (hnew-htemp)

		class(ty_com_nodes)			,allocatable::nodes
		class(ty_com_elements)	,allocatable::elements
		type(ty_com_boundary)		,pointer::boundary

		real(kind=dps),allocatable::rhs(:)			!<Calculated right hand side vector

		real(kind=dps),allocatable::solution(:)	!<Vector with the solution of the linear system
		real(kind=dps),allocatable::slopeold(:) !<Equal to (hnew-hold)/dt for previous step in order to estimate new h.

		class(ty_mx),allocatable::mx			!<Matrix of the FEM linear system (left side of the linear system)

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
	!> Set hdirichlet to the node idnode
	!> @param[in] idnode			ID of the node
	!> @param[in] hdirichlet	Piezometric pressure to impose at the node
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
	! Set Newmann condition to node
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Method to set the waterflow entering/leving the model on the node
	!> @param[in] idnode	  ID of the node
	!> @param[in] qnewman		Waterflow to impose on the node as a boundary condition
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
	!> Method to allocate all allocatable properties of the derived type(class)
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_com_calc_allocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_allocateall" :: s_com_calc_allocateall
	!DEC$ endif

	class(ty_com_calc),intent(inout)::this
	integer::nnc !<Number of nodes-classes

	nnc=this%nodes%count

	!Allocate all vectors
	if(.not.allocated(this%slopeold)) allocate(this%slopeold(nnc))
	if(.not.allocated(this%rhs))			allocate(this%rhs(nnc))
	if(.not.allocated(this%solution)) allocate(this%solution(nnc))

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
	!> Method to deallocate all allocatable properties of the derived type
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_deallocateall" :: s_com_calc_deallocateall
	!DEC$ endif

	class(ty_com_calc),intent(inout)::this

	if(allocated(this%slopeold))	deallocate(this%slopeold)
	if(allocated(this%rhs))				deallocate(this%rhs)
	if(allocated(this%solution))	deallocate(this%solution)
	if(allocated(this%perm))			deallocate(this%perm)

	if(allocated(this%mx))				deallocate(this%mx)

	end subroutine s_com_calc_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	! SHAPEMATRIX
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

	mx			=>this%mx

	!initialize matrix.....
	select type(mx)
	type is(ty_mx_dense) !dense matrix---------------------------------------------------------------------------------------

		!initiallize...
		call mx       %init(numnp,numnp)

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

		!fills csrows, csrcols and this%pos...
		ncount=0

		mx%csrrows(1) = 1 !first element of csrrows
		do n=1,numnp
			do m=1,numnp
				if (mx%pos(n,m).ne.0) then
					ncount=ncount+1
					mx%csrcols(ncount) = m
					mx%pos(n,m)=ncount
				end if
			end do
			mx%csrrows(n+1) = ncount+1
		end do

	type is(ty_mx_banded) !banded matrix------------------------------------------------------------------------------------------------

		!initialize matrix...
		call mx      %init(numnp,kl,ku)

	end select

	end subroutine s_com_calc_shapematrix

	!---------------------------------------------------------------------------------------------------------------------
	!	ALLOCATE OTHER MATRIX
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate matrix mx
	!> @param[inout] mx			Matrix to be allocated
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
	!	BUILD COEFMATRIX
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> This procedure build the coefficient matrix mx given the function fun_quad(chi,e) which multiplies basis functions.
	!> @param[in] fun_quad	Function to integrate (with arguments chi and the element e) muliplying basis functions and used to fill matrix rows and columns
	!> @param[inout] mx			Matrix to be build
	!> @param[in] masslump	If true, the matrix is lumped diagonally (is diagonal)
	!> @param[in] optionbasis Define if the basis is (1)h, (2)dh, (3)h·h, (4)h·dh, (5)dh·h, (6)dh·dh
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

	use com_mod_fem_shapefunctions,only:shape1d,dphi1d
	use com_mod_fem_shapefunctions,only:jacobian1d
	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::auxfun(size(chi)),xne(numnodclass)
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
		dphifunresult = dphi1d(chi,xne)
		auxfun = fun_quad(chi,e)*shapefunresult(:,ishape)*dphifunresult(:,jshape)
	case(OPTIONBASIS_DH_H) !f(chi,e)·dbasis_i(chi)·basis_j(chi)
		nstartelem	= this%elements%idnode(e,1)
		nendelem		= this%elements%idnode(e,numnodclass)
		xne = this%nodes%x(nstartelem : nendelem)
		shapefunresult = shape1d(chi,this%elements%chi(e,:))
		dphifunresult = dphi1d(chi,xne)
		auxfun = fun_quad(chi,e)*dphifunresult(:,ishape)*shapefunresult(:,jshape)
	case(OPTIONBASIS_DH_DH)	!f(chi,e)·dbasis_i(chi)·dbasis_j(chi)
		nstartelem	= this%elements%idnode(e,1)
		nendelem		= this%elements%idnode(e,numnodclass)
		xne = this%nodes%x(nstartelem : nendelem)
		dphifunresult = dphi1d(chi,xne)
		auxfun = fun_quad(chi,e)*dphifunresult(:,ishape)*dphifunresult(:,jshape)
	end select

	end function auxfun

	end subroutine s_com_calc_buildcoefmatrix

	!---------------------------------------------------------------------------------------------------------------------
	!	ESTIMATE HTEMP FOR NEW ITERATION
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Set htemp from previouly calculated hnew
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_estimate_htemp_for_new_iteration(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_estimate_htemp_for_new_iteration" :: s_com_calc_estimate_htemp_for_new_iteration
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this

	this%nodes%htemp = this%nodes%hnew

	end subroutine s_com_calc_estimate_htemp_for_new_iteration

	!---------------------------------------------------------------------------------------------------------------------
	!	CALCULATE ITERATION
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Perform one iteration
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
	logical,intent(inout)::IsConverged

	integer::numnp
	real(kind=dps)::ksal

	numnp = this%nodes%count

	call this%estimate_htemp_for_new_iteration()

	!build linear system
	call this%build_linearsystem(IS_TIME_DEPENDANT,USE_H_TEMP)

	!Assign Dirichlet conditions
	call this%assign_dirichlet()

	! solve linear algebra system
	call solve_linalg(this%mx,this%rhs,this%solution,this%parameters%typesolver,this%perm)

	!update hnew wiith overrelaxation: (hnew = htemp+crelax·(hnew-htemp)
	call this%update_hnew_from_solution()

	!calculate error in h (in this case error in nodes):
	this%epsh = maxval(abs(this%nodes%hnew-this%nodes%htemp)) !Chechk (If we are updating hnew to only positive, this is the one that has to be put.

	isconverged = this%epsh<this%parameters%epsh_tol

	end subroutine s_com_calc_iteration

	!---------------------------------------------------------------------------------------------------------------------
	! UPDATE HNEW FROM SOLUTION
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Set the value of hnew from the value of the solution of the linear system using the coefficient of relaxation
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_update_hnew_from_solution(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_update_hnew_from_solution" :: s_com_calc_update_hnew_from_solution
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this

	this%nodes%hnew = this%nodes%htemp + this%parameters%crelax*(this%solution-this%nodes%htemp)

	end subroutine s_com_calc_update_hnew_from_solution

	!---------------------------------------------------------------------------------------------------------------------
	! SET HNEW FROM INITIAL VALUES
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Set the value of hnew from the definition of the initial conditions
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
	! REVERT TO OLD
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Revert values to the ones in the previous timestep
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_revert_to_old(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_revert_to_old" :: s_com_calc_revert_to_old
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this
	this%time%t			= this%time%told

	!All values to initial conditions
	this%nodes%hnew		= this%nodes%hold
	this%nodes%htemp	= this%nodes%hold

	!All values to initial conditions
	this%nodes%dhnew	= this%nodes%dhold
	this%nodes%dhtemp	= this%nodes%dhold

	end subroutine s_com_calc_revert_to_old


	!---------------------------------------------------------------------------------------------------------------------
	! SET OLD VALUES
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Set old values from current calculation of the timestep
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_set_old(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_set_old" :: s_com_calc_set_old
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this

	if(this%time%told.ne.this%time%t) then
		this%time%dtold	= (this%time%t-this%time%told)
		this%time%told	= this%time%t
	end if
	this%slopeold = min(1.0_dpd,(this%nodes%hnew-this%nodes%hold)/this%time%dtold)
	this%nodes%hold = this%nodes%hnew

	end subroutine s_com_calc_set_old

	!---------------------------------------------------------------------------------------------------------------------
	! ESTIMATE HNEW FOR NEW TIMESTEP
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Estimate the value of hew to begin iterations in the current timestep
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_estimate_hnew_for_new_timestep(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_estimate_hnew_for_new_timestep" :: s_com_calc_estimate_hnew_for_new_timestep
	!DEC$ endif
	class(ty_com_calc),intent(inout)::this

	this%nodes%hnew	= this%nodes%hold+0.5_dpd*this%slopeold*this%time%dt

	end subroutine s_com_calc_estimate_hnew_for_new_timestep

	!---------------------------------------------------------------------------------------------------------------------
	! ASSIGN DIRECHLET FROM NODES
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Modify matrix and RHS with Dirichlet conditions
	!---------------------------------------------------------------------------------------------------------------------

	SUBROUTINE s_vector_assign_dirichlet_from_nodes(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_vector_assign_dirichlet_from_nodes" :: s_vector_assign_dirichlet_from_nodes
	!DEC$ endif
	CLASS(ty_com_calc),INTENT(INOUT)::this

	integer::n,m
	real(kind=dps)::hnod

	DO n = 1, this%nodes%count

		!--------
		IF (this%nodes%Isdirichlet(n)) THEN !Assign dirichlet pressure condition to node
			DO m=1,this%nodes%count
				IF(this%MX%Get(n,m).NE.0) this%RHS(m) = this%RHS(m)-this%nodes%HDirichlet(n)*this%MX%Get(n,m)
				IF(this%MX%Get(n,m).NE.0) CALL this%MX%Set(n,m,0.0_dpd)
				IF(this%MX%Get(m,n).NE.0) CALL this%MX%Set(m,n,0.0_dpd)
			END DO
			CALL this%MX%Set(n,n,1.0_dpd)
			this%rhs(n) = this%nodes%HDirichlet(n)
		END IF

	END DO

	END SUBROUTINE s_vector_assign_dirichlet_from_nodes

	!---------------------------------------------------------------------------------------------------------------------
	! GET IDELEMENT FROM X
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Get the element ID given the absolute coordinate x over the model
	!---------------------------------------------------------------------------------------------------------------------

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

