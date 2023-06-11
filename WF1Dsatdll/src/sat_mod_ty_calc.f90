	!********************************************************************************************************************
	! TITLE         : SAT_MOD_TY_CALC: EXTENDED DERIVED TYPE OF COM_MOD_TY_CALC TO INCLUDE PROPERTIES AND METHODS OF WF1DSAT
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_CALC
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Extended derived type of com_mod_ty_calc to include properties and methods of wf1dsat
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module sat_mod_ty_calc
	use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded
	use com_mod_ty_calc,			only: ty_com_calc
	use com_mod_ty_layers,		only: ty_com_layers
	use sat_mod_ty_nodes,			only: ty_sat_nodes
	use sat_mod_ty_elements,	only: ty_sat_elements

	implicit none
	include 'inc_precision.fi'

	private

	type,extends(ty_com_calc),public::ty_sat_calc	!< CLASS: Definition of the layer in saturated model
		type(ty_com_layers)	,pointer::layers

		real(kind=dps),allocatable::nrel(:)	!<Factor to correct the relative porosity when coupling with the unsaturated model on an increasing waterlevel.

		class(ty_mx),allocatable::mxmass	!<Mass matrix: of the FEM linear system
		class(ty_mx),allocatable::mxstiff !<Stiffness matrix: (depend on results but not on time)
		class(ty_mx),allocatable::mxload	!<Loads matrix: Loads on the FEM
		class(ty_mx),allocatable::mxbound !<Boundary matrix: Boundary conditions of FEM



	contains
	procedure,public:: allocateall 							=> s_sat_calc_allocateall
	procedure,public:: construct 								=> s_sat_calc_construct
	procedure,public:: deallocateall 						=> s_sat_calc_deallocateall
	procedure,public:: build_linearsystem 			=> s_sat_calc_build_linearsystem
	procedure,public:: update_hnew_from_solution	=> s_sat_calc_update_hnew_from_solution
	procedure,public:: estimate_htemp_for_new_iteration	=> s_sat_calc_estimate_htemp_for_new_iteration
	procedure,public:: get_results_elements			=> s_sat_calc_results_in_elements
	procedure,public:: get_results_nodes				=> s_sat_calc_results_in_nodes
	end type ty_sat_calc

	contains

	!---------------------------------------------------------------------------------------------------------------------
	! ALLOCATE ALL (OVERRIDE)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all allocatable properties of the class
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_sat_calc_allocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_calc_allocateall
	!DEC$ endif
	use sat_mod_ty_nodes, only:ty_sat_nodes
	use sat_mod_ty_elements, only:ty_sat_elements
	use com_mod_ty_parameters, only:ty_com_parameters

	class(ty_sat_calc),intent(inout)::this

	if(.not.allocated(this%nodes))allocate(ty_sat_nodes::this%nodes)
	if(.not.allocated(this%elements))allocate(ty_sat_elements::this%elements)

	end subroutine s_sat_calc_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	! CONSTRUCTOR (OVERRIDE)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Construct the class.
	!> @param[in] ntype
	!> @param[inout] mx
	!> @param[in] option
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_sat_calc_construct(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_calc_construct
	!DEC$ endif
	use com_mod_ty_elements,	only: ty_com_elements
	use sat_mod_ty_elements,	only: ty_sat_elements
	use com_mod_ty_nodes,	only: ty_com_nodes
	use sat_mod_ty_nodes,	only: ty_sat_nodes
	integer,parameter::KU=1,KL=1

	class(ty_sat_calc),intent(inout),target::this
	class(ty_com_elements),pointer::elemcom
	class(ty_sat_elements),pointer::elem
	class(ty_com_nodes),pointer::nodescom
	class(ty_sat_nodes),pointer::nodes
	integer::nnc !<Number of nodes-classes

	nnc=this%nodes%count
	elemcom=>this%elements
	select type(elemcom)
	type is(ty_sat_elements)
		elem => elemcom
	end select

	nodescom=>this%nodes
	select type(nodescom)
	type is(ty_sat_nodes)
		nodes => nodescom
	end select

	!Allocate all vectors
	if(.not.allocated(this%nrel))	allocate(this%nrel(nnc))
	this%nrel = 1.0_dpd
	if(.not.allocated(this%slopeold)) allocate(this%slopeold(nnc))
	if(.not.allocated(this%rhs))			allocate(this%rhs(nnc))
	if(.not.allocated(this%solution)) allocate(this%solution(nnc))
	if(.not.allocated(this%perm))			allocate(this%perm(nnc)) !If perm is not allocated then it is not used.

	if(.not.allocated(elem%results_qent))		allocate(elem%results_qent(this%elements%count))
	if(.not.allocated(elem%results_qsat))		allocate(elem%results_qsat(this%elements%count))
	if(.not.allocated(elem%results_qinf))		allocate(elem%results_qinf(this%elements%count))
	if(.not.allocated(elem%results_inchnew))		allocate(elem%results_inchnew(this%elements%count))
	if(.not.allocated(elem%results_incvoldt)) allocate(elem%results_incvoldt(this%elements%count))
	if(.not.allocated(elem%results_dqhordx))allocate(elem%results_dqhordx(this%elements%count))
	if(.not.allocated(elem%results_dqhordx_from_incvoldt))allocate(elem%results_dqhordx_from_incvoldt(this%elements%count))
	if(.not.allocated(elem%results_dqhordx_all))allocate(elem%results_dqhordx_all(this%elements%count))
	if(.not.allocated(elem%results_dqhordx_from_incvoldt_all))allocate(elem%results_dqhordx_from_incvoldt_all(this%elements%count))
	if(.not.allocated(elem%x0))		allocate(elem%x0(this%elements%count))
	if(.not.allocated(elem%x1))		allocate(elem%x1(this%elements%count))
	if(.not.allocated(elem%z))		allocate(elem%z(this%elements%count))
	if(.not.allocated(elem%hz))		allocate(elem%hz(this%elements%count))
	if(.not.allocated(elem%hzold))		allocate(elem%hzold(this%elements%count))
	if(.not.allocated(elem%dhdx))		allocate(elem%dhdx(this%elements%count))
	if(.not.allocated(elem%dhzdx))		allocate(elem%dhzdx(this%elements%count))
	if(.not.allocated(elem%h0))		allocate(elem%h0(this%elements%count))
	if(.not.allocated(elem%h1))		allocate(elem%h1(this%elements%count))
	if(.not.allocated(elem%z0))		allocate(elem%z0(this%elements%count))
	if(.not.allocated(elem%z1))		allocate(elem%z1(this%elements%count))
	if(.not.allocated(elem%dhdx0))		allocate(elem%dhdx0(this%elements%count))
	if(.not.allocated(elem%dhdx1))		allocate(elem%dhdx1(this%elements%count))
	if(.not.allocated(elem%dhzdx0))		allocate(elem%dhzdx0(this%elements%count))
	if(.not.allocated(elem%dhzdx1))		allocate(elem%dhzdx1(this%elements%count))
	if(.not.allocated(elem%q0))		allocate(elem%q0(this%elements%count))
	if(.not.allocated(elem%q1))		allocate(elem%q1(this%elements%count))
	if(.not.allocated(elem%q0_all))		allocate(elem%q0_all(this%elements%count))
	if(.not.allocated(elem%q1_all))		allocate(elem%q1_all(this%elements%count))
	if(.not.allocated(elem%q))		allocate(elem%q(this%elements%count))
	if(.not.allocated(elem%q_all))		allocate(elem%q_all(this%elements%count))
	if(.not.allocated(elem%k0))		allocate(elem%k0(this%elements%count))
	if(.not.allocated(elem%k1))		allocate(elem%k1(this%elements%count))
	if(.not.allocated(elem%k0_all))		allocate(elem%k0_all(this%elements%count))
	if(.not.allocated(elem%k1_all))		allocate(elem%k1_all(this%elements%count))
	if(.not.allocated(elem%k))		allocate(elem%k(this%elements%count))

	elem%results_qent=0.0_dpd
	elem%results_qsat=0.0_dpd
	elem%results_qinf=0.0_dpd
	elem%results_inchnew=0.0_dpd	
	elem%results_incvoldt=0.0_dpd
	elem%results_dqhordx=0.0_dpd
	elem%results_dqhordx_from_incvoldt=0.0_dpd
	elem%results_dqhordx_all=0.0_dpd
	elem%results_dqhordx_from_incvoldt_all=0.0_dpd

	if(.not.allocated(nodes%results_qent))		allocate(nodes%results_qent(nodes%count))
	if(.not.allocated(nodes%results_qsat))		allocate(nodes%results_qsat(nodes%count))
	if(.not.allocated(nodes%results_qinf))		allocate(nodes%results_qinf(nodes%count))
	if(.not.allocated(nodes%results_inchnew))		allocate(nodes%results_inchnew(nodes%count))
	if(.not.allocated(nodes%results_incvoldt))	allocate(nodes%results_incvoldt(nodes%count))
	if(.not.allocated(nodes%results_qhor))		allocate(nodes%results_qhor(nodes%count))
	if(.not.allocated(nodes%results_dqhordx))	allocate(nodes%results_dqhordx(nodes%count))
	if(.not.allocated(nodes%results_qhor_all))		allocate(nodes%results_qhor_all(nodes%count))
	if(.not.allocated(nodes%results_dqhordx_all))	allocate(nodes%results_dqhordx_all(nodes%count))
	if(.not.allocated(nodes%results_dqhordx_from_incvoldt))	allocate(nodes%results_dqhordx_from_incvoldt(nodes%count))
	if(.not.allocated(nodes%results_dqhordx_from_incvoldt_all))	allocate(nodes%results_dqhordx_from_incvoldt_all(nodes%count))
	nodes%results_qent								=0.0_dpd
	nodes%results_qsat								=0.0_dpd
	nodes%results_qinf								=0.0_dpd
	nodes%results_inchnew							=0.0_dpd
	nodes%results_incvoldt						=0.0_dpd
	nodes%results_incvoldt							=0.0_dpd
	nodes%results_qhor								=0.0_dpd
	nodes%results_dqhordx							=0.0_dpd
	nodes%results_qhor_all						=0.0_dpd
	nodes%results_dqhordx_all					=0.0_dpd
	nodes%results_dqhordx_from_incvoldt_all	=0.0_dpd

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
	!Allocate other matrix with same spartity and shape than mx...
	call this%allocate_matrix(this%mxmass)
	call this%allocate_matrix(this%mxstiff)
	call this%allocate_matrix(this%mxload)
	call this%allocate_matrix(this%mxbound)

	this%slopeold = 0.0_dpd

	end subroutine s_sat_calc_construct


	!---------------------------------------------------------------------------------------------------------------------
	! DEALLOCATE ALL (OVERRIDE)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate all allocatable properties of the class
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_calc_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_calc_deallocateall
	!DEC$ endif

	class(ty_sat_calc),intent(inout)::this

	if(allocated(this%slopeold))		deallocate(this%slopeold)
	if(allocated(this%rhs))				deallocate(this%rhs)
	if(allocated(this%solution))		deallocate(this%solution)
	if(allocated(this%perm))				deallocate(this%perm)
	if(allocated(this%nrel))				deallocate(this%nrel)

	if(allocated(this%mx))					deallocate(this%mx)
	if(allocated(this%mxmass))			deallocate(this%mxmass)
	if(allocated(this%mxstiff))			deallocate(this%mxstiff)
	if(allocated(this%mxload))			deallocate(this%mxload)
	if(allocated(this%mxbound))			deallocate(this%mxbound)

	if(allocated(this%elements)) then
		call this%elements%deallocateall()
		deallocate(this%elements)
	end if

	if(allocated(this%nodes)) then
		call this%nodes%deallocateall()
		deallocate(this%nodes)
	end if

	end subroutine s_sat_calc_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	! BUILD LINEAR SYSTEM
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Build the linear system of the FEM methods
	!> @param[in] IsTimeDependant If true, the matrix depend on time
	!> @param[in] option OPTIONBASIS_H=1,OPTIONBASIS_DH=2,OPTIONBASIS_H_H=3,OPTIONBASIS_H_DH=4,OPTIONBASIS_DH_H=5,OPTIONBASIS_DH_DH=6
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_calc_build_linearsystem(this,IsTimeDependant,option)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_calc_build_linearsystem
	!DEC$ endif
	use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded

	integer,parameter::COEF_MASS=1,COEF_STIFF=2,COEF_LOAD=3,COEF_BOUND=4
	integer,parameter::OPTIONBASIS_H=1,OPTIONBASIS_DH=2,OPTIONBASIS_H_H=3,OPTIONBASIS_H_DH=4,OPTIONBASIS_DH_H=5,OPTIONBASIS_DH_DH=6
	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDERHOLD=2
	logical,parameter::IS_EQ_FIRST_VERSION=.TRUE.
	class(ty_sat_calc),intent(inout),target::this
	logical,intent(in)::IsTimeDependant
	integer,intent(in),optional::option
	integer::opt,nelem,numnod,numclass,numnodclass

	real(kind=dpd)::henodes(this%elements%nn*(this%elements%nc+1))
	real(kind=dpd)::dhzenodes(this%elements%nn*(this%elements%nc+1))
	real(kind=dpd)::henodesold(this%elements%nn*(this%elements%nc+1))
	real(kind=dpd)::nrelnodes(this%elements%nn*(this%elements%nc+1))

	opt=0
	if(present(option)) opt = option

	nelem				=this%elements%count
	numnod			=this%elements%nn
	numclass		=this%elements%nc+1
	numnodclass =numnod*numclass

	!Mass matrix:
	!tex:\[\left[ {{M_{MASSi,j}}} \right] = \left[ {\mathop \smallint \limits_{\rm{\Omega }} {{\overline {{n_{nw}}} }^{k,n + 1}}{\phi _j}{\phi _i}d{\rm{\Omega }}} \right] = \left[ {{M_{MASSi,j}}} \right] = \left[ {\mathop \smallint \limits_{\rm{\Omega }} {r_{rel}}\cdot\frac{{\sum\limits_l {({\theta _{sat,l}} - {\theta _{res,l}})\Delta {h_l}} }}{{\Delta h}}{\phi _j}{\phi _i}d{\rm{\Omega }}} \right]{\overline {{n_{nw}}} ^{k,n + 1}}\]
	if(IsTimeDependant)				call this%buildcoefmatrix(f_mass,this%mxmass,this%parameters%masslump,OPTIONBASIS_H_H)	!Build mass matrix
	!Stiffness matrix:


	!1st version:
	!tex:
	!$\left[ {{M_{Ki,j}}} \right] = \left[ {\int\limits_\Omega  {k_{eff}^{k,n + 1} \cdot {h^{k,n + 1}}\cdot\frac{{\partial {\phi _j}}}{{\partial x}}\frac{{\partial {\phi _i}}}{{\partial x}}} d\Omega } \right]$

	!2nd version (not working):
	!tex:
	!\[\left[ {{M_{Ki,j}}} \right] = \left[ {\int\limits_\Omega  {k_{eff}^{k,n + 1} \cdot \left( {\frac{{\partial {h^{k,n + 1}}}}{{\partial x}} + \frac{{\partial z}}{{\partial x}}} \right) \cdot {\phi _j}\frac{{\partial {\phi _i}}}{{\partial x}}} d\Omega } \right]\]

	if(IS_EQ_FIRST_VERSION.and.IsTimeDependant)				call this%buildcoefmatrix(f_stiffA,this%mxstiff,.false.,OPTIONBASIS_DH_DH)
	if(.not.IS_EQ_FIRST_VERSION.and.IsTimeDependant)	call this%buildcoefmatrix(f_stiffB,this%mxstiff,.false.,OPTIONBASIS_H_DH)								!Build stiff matrix


	!Identity matrix for the waterflow:
	!tex:$\left[ {{I_{i,j}}} \right] = \left[ {\mathop \smallint \limits_{\rm{\Omega }} {\phi _j}{\phi _i}{\rm{\;}}d{\rm{\Omega }}} \right]$
	if(.not.IsTimeDependant)	call this%buildcoefmatrix(f_identity,this%mxload,.false.,OPTIONBASIS_H_H)								!Build load matrix
	!Boundary matrix:
	!tex:
	!	$\left[ {{M_{boundi,j}}} \right] = \left[ {\begin{array}{*{20}{c}}$
	!0&0&0\\
	!0& \ddots & \vdots \\
	!0& \cdots &{C·k_{eff}^{k,n + 1}}
	!\end{array}} \right]$
	if(IsTimeDependant)				call s_aux_buildcoefbound(this,this%mxbound)																						!Build boundary matrix (in this case pressure pass from h to 0 in a length of kmul·h)

	if (isTimeDependant) then
		!1st version:
		!tex: \[\left( {\frac{1}{{\Delta t}}\left[ {{M_{MASSi,j}}} \right] + \left[ {{M_{Ki,j}}} \right] + \left[ {{M_{boundi,j}}} \right]} \right)\left\{ {\psi _j^{k + 1,n + 1}} \right\} = \frac{1}{{\Delta t}}\left[ {{M_{MASSi,j}}} \right]\left\{ {\psi _j^n} \right\} - \left[ {{M_{Ki,j}}} \right]\left\{ {{z_j}} \right\} + \left[ {{I_{i,j}}} \right]\left\{ {{{\rm{q}}_{{\rm{vtb}},{\rm{j}}}}} \right\}\]

		!2nd version:
		!tex:
		!$\left( {\frac{1}{{\Delta t}}\left[ {{M_{MASSi,j}}} \right] + \left[ {{M_{Ki,j}}} \right] + \left[ {{M_{boundi,j}}} \right]} \right)\left\{ {\psi _j^{k + 1,n + 1}} \right\} = \frac{1}{{\Delta t}}\left[ {{M_{MASSi,j}}} \right]\left\{ {\psi _j^n} \right\} + \left[ {{I_{i,j}}} \right]\left\{ {{{\rm{q}}_{{\rm{vtb}},{\rm{i}}}}} \right\}$

		!create mx matrix: [mx]=[mass]/dt+[stiff]+[bound]
		this%mx = this%mxmass/this%time%dt+this%mxstiff+this%mxbound

		!1dt version:
		if(IS_EQ_FIRST_VERSION)				this%rhs = this%mxmass*this%nodes%hold/this%time%dt-this%mxstiff*this%nodes%z+this%mxload*this%nodes%qent
		!2nd version:
		if(.not.IS_EQ_FIRST_VERSION)	this%rhs = this%mxmass*this%nodes%hold/this%time%dt												   +this%mxload*this%nodes%qent

	end if

	contains
	!--------------------------------------------------------------------------
	!Functions to integrate in each kind of matrix (for [LOAD] matrix)
	function f_identity(chi,e)
	real(kind=dps),intent(in)::chi(:)
	integer,intent(in)::e
	real(kind=dps)::f_identity(size(chi))
	f_identity = 1.0_dps
	end function f_identity
	!-----------------------------------------------------------------------------------------
	!Returns the increment of water content from hold to htemp in coord chi of element e (for [MASS] matrix)
	function f_mass(chi,e)
	use com_mod_fem_shapefunctions,only:interp_on_element

	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDER_HOLD=2
	real(kind=dps),intent(in)::chi(:)
	integer,intent(in)::e
	real(kind=dps)::f_mass(size(chi))
	real(kind=dps)::nreltemp(size(chi)),headtemp(size(chi)),headold(size(chi))

	integer::nstartelem,nendelem

	nstartelem	= this%elements%idnode(e,1)
	nendelem		= this%elements%idnode(e,numnodclass)
	select case (opt)
	case (CONSIDER_HNEW)
		henodes = this%nodes%hnew(nstartelem : nendelem)
	case (CONSIDER_HTEMP)
		henodes = this%nodes%htemp(nstartelem : nendelem)
	case (CONSIDER_HOLD)
		henodes = this%nodes%hold(nstartelem : nendelem)
		case default
		stop('Matrix cannot be builded from old values in nodes')
	end select

	henodesold	= this%nodes%hold(nstartelem : nendelem)
	nrelnodes		= this%nrel(nstartelem : nendelem)

	headtemp	= interp_on_element(chi,henodes)		!Returns htemp in each chi
	headold		= interp_on_element(chi,henodesold)	!Returns hold in each chi
	nreltemp	= interp_on_element(chi,nrelnodes)	!Return relative porosity in chi. (non-wettin voids divided by total voids as calculated in the unsaturated vertical model)

	!tex:
	!$\left[ {{M_{MASSi,j}}} \right] = \left[ {\mathop \smallint \limits_{\rm{\Omega }} {{\overline {{n_{nw}}} }^{k,n + 1}}{\phi _j}{\phi _i}d{\rm{\Omega }}} \right]$
	!tex:
	!\[{\overline {{n_{nw}}} ^{k,n + 1}} = {r_{rel}}\cdot\overline {({\theta _{sat,l}} - {\theta _{res,l}})}  = {r_{rel}}\cdot\frac{{\sum\limits_l {({\theta _{sat,l}} - {\theta _{res,l}})\Delta {h_l}} }}{{\Delta h}}\]
	f_mass = nreltemp*this%layers%get_water_inc_med(headold,headtemp)

	end function f_mass

	!-----------------------------------------------------------------------------------------
	!Returns ksat·head... (for [stiff] matrix)
	function f_stiffA(chi,e)
	use com_mod_fem_shapefunctions,only:interp_on_element

	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDER_HOLD=2
	real(kind=dps),intent(in)::chi(:)
	integer,intent(in)::e
	real(kind=dps)::f_stiffA(size(chi))
	integer::nstartelem,nendelem
	real(kind=dps)::head(size(chi),this%layers%count+1),ksat(this%layers%count+1),headonchi(size(chi)),kmed(size(chi)),headoldonchi(size(chi))

	nstartelem	= this%elements%idnode(e,1)
	nendelem		= this%elements%idnode(e,numnodclass)
	select case (opt)
	case (CONSIDER_HNEW)
		henodes = this%nodes%hnew(nstartelem : nendelem)
	case (CONSIDER_HTEMP)
		henodes = this%nodes%htemp(nstartelem : nendelem)
	case (CONSIDER_HOLD)
		henodes = this%nodes%hold(nstartelem : nendelem)
		case default
		stop('Matrix cannot be builded from old values in nodes')
	end select
	!Get a vector of ksat for every layer and another ksat=1.0 over the top layer
	ksat(1:this%layers%count) = this%layers%material(:)%ksat
	ksat(this%layers%count+1) = 1.0_dpd !Assumed that permeatility over the top layer is 1m3/2

	!Get a matrix for the height of water in each layer (for each value of chi) and over the top layer
	headonchi = interp_on_element(chi,henodes)
	!headoldonchi = interp_on_element(chi,this%nodes%hold(nstartelem : nendelem))
	head = this%layers%get_inc_h_from0_vec(headonchi)
	kmed = matmul(head,ksat)
	where (headonchi<=0)
		kmed = ksat(1)
	else where
		kmed = abs(kmed/headonchi)
	end where

	f_stiffA = kmed * headonchi

	end function f_stiffA

	!-----------------------------------------------------------------------------------------
	!Returns ksat·head... (for [stiff] matrix)
	function f_stiffB(chi,e)
	use com_mod_fem_shapefunctions,only:interp_on_element
	use com_mod_fem_shapefunctions,only:dphi1d

	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDER_HOLD=2
	real(kind=dps),intent(in)::chi(:)
	integer,intent(in)::e
	real(kind=dps)::f_stiffB(size(chi))
	integer::nstartelem,nendelem
	real(kind=dps)::head(size(chi),this%layers%count+1),ksat(this%layers%count+1),headonchi(size(chi)),kmed(size(chi)),headoldonchi(size(chi)),dhzonchi(size(chi))

	nstartelem	= this%elements%idnode(e,1)
	nendelem		= this%elements%idnode(e,numnodclass)
	select case (opt)
	case (CONSIDER_HNEW)
		henodes = this%nodes%hnew(nstartelem : nendelem)
		dhzenodes= matmul(dphi1d(chi,this%nodes%x(nstartelem : nendelem)),this%nodes%hnew(nstartelem : nendelem)+this%nodes%z(nstartelem : nendelem))
	case (CONSIDER_HTEMP)
		henodes = this%nodes%hnew(nstartelem : nendelem)
		dhzenodes= matmul(dphi1d(chi,this%nodes%x(nstartelem : nendelem)),this%nodes%htemp(nstartelem : nendelem)+this%nodes%z(nstartelem : nendelem))
	case (CONSIDER_HOLD)
		henodes = this%nodes%hnew(nstartelem : nendelem)
		dhzenodes= matmul(dphi1d(chi,this%nodes%x(nstartelem : nendelem)),this%nodes%hold(nstartelem : nendelem)+this%nodes%z(nstartelem : nendelem))
		case default
		stop('Matrix cannot be builded from old values in nodes')
	end select

	!Get a vector of ksat for every layer and another ksat=1.0 over the top layer
	ksat(1:this%layers%count) = this%layers%material(:)%ksat
	ksat(this%layers%count+1) = 1.0_dpd !Assumed that permeatility over the top layer is 1m3/2

	!Get a matrix for the height of water in each layer (for each value of chi) and over the top layer
	headonchi = interp_on_element(chi,henodes)
	dhzonchi =  interp_on_element(chi,dhzenodes)
	head = this%layers%get_inc_h_from0_vec(headonchi)

	kmed = matmul(head,ksat)
	where (headonchi<=0)
		kmed = ksat(1)
	else where
		kmed = abs(kmed/headonchi)
	end where

	!tex:${k_{eff}^{k,n + 1}\cdot\left( {\frac{{\partial {h^{k,n + 1}}}}{{\partial x}} + \frac{{\partial z}}{{\partial x}}} \right)\cdot}$
	f_stiffB = kmed * dhzonchi

	end function f_stiffB

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Aux subroutine to build boundary matrix
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_aux_buildcoefbound(this,mx)
	use com_mod_mx_datatypes, only:ty_mx

	class(ty_sat_calc),intent(in)::this
	class(ty_mx),intent(inout)::mx

	integer::nnodes
	real(kind=dpd)::ksal,multksal
	real(kind=dpd)::ktemp(this%layers%count+1),htemp(this%layers%count+1)

	multksal = min(this%parameters%multksal,(this%nodes%hold(this%nodes%count-1)/2.0_dpd+this%nodes%z(this%nodes%count-1)-this%nodes%z(this%nodes%count))/abs(this%nodes%x(this%nodes%count)-this%nodes%x(this%nodes%count-1))) !As boundary dh/dx is the min between multk and (0.5·hi-1+(zi-1-zi))/(xi-1-xi)
	nnodes = this%nodes%count
	if (this%nodes%htemp(this%nodes%count)==0.0_dpd) then
		ksal= 0.0_dpd
	else
		ktemp(1:this%layers%count) = this%layers%material(:)%ksat
		ktemp(this%layers%count+1) = 1.0_dpd
		htemp = this%layers%get_inc_h_from0_sca(this%nodes%htemp(this%nodes%count))
		ksal = dot_product(ktemp,htemp)/this%nodes%htemp(this%nodes%count)*multksal
	end if
	!tex:
	!	$\left[ {{M_{boundi,j}}} \right] = \left[ {\begin{array}{*{20}{c}}
	!0&0&0\\
	!0& \ddots & \vdots \\
	!0& \cdots &{k_{eff}^{k,n + 1}}
	!\end{array}} \right]$
	mx = 0.0_dpd
	call mx%set(nnodes,nnodes,ksal)

	end subroutine s_aux_buildcoefbound

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Aux subroutine to build boundary matrix (In this case the permeability increase to the boundary).
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------


	subroutine s_aux_buildcoefbound_multk(this,mx)
	use com_mod_mx_datatypes, only:ty_mx
	use com_mod_fem_shapefunctions,only:shape1d,dphi1d
	use com_mod_fem_shapefunctions,only:jacobian1d

	class(ty_sat_calc),intent(in)::this
	class(ty_mx),intent(inout)::mx
	integer::nnodes
	real(kind=dpd)::ksal
	real(kind=dps)::dphifunresult(this%elements%nn*(this%elements%nc+1)),xne(this%elements%nn*(this%elements%nc+1))
	integer::nstartelem,nendelem,i


	nelem=this%elements%count
	numnod=this%elements%nn
	numclass = this%elements%nc+1
	numnodclass = numnod*numclass

	nstartelem	= this%elements%idnode(nelem,1)
	nendelem		= this%elements%idnode(nelem,numnodclass)

	xne = this%nodes%x(nstartelem : nendelem)
	dphifunresult = dphi1d(1.0_dpd,xne) !Derivatives at the end of the element

	mx = 0.0_dpd
	ksal = -this%layers%material(1)%ksat*min(this%parameters%multksal,(this%nodes%htemp(numnodclass-1)+this%nodes%z(numnodclass-1)-0.75_dpd*this%nodes%htemp(numnodclass)-this%nodes%z(numnodclass))/abs(this%nodes%x(numnodclass)-this%nodes%x(numnodclass-1)))
	do i=1,numnodclass
		call mx%set(nelem+1,nstartelem+i-1,ksal*this%nodes%htemp(numnodclass)*dphifunresult(i))
	end do

	end subroutine s_aux_buildcoefbound_multk

	end subroutine s_sat_calc_build_linearsystem

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Update value of hnew from the solution
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_calc_update_hnew_from_solution(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_calc_update_hnew_from_solution
	!DEC$ endif
	class(ty_sat_calc),intent(inout)::this

	this%nodes%hnew = max(0.0_dpd,this%nodes%htemp + 0.25*this%parameters%crelax*(this%solution-this%nodes%htemp))

	end subroutine s_sat_calc_update_hnew_from_solution

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Estimate htemp for the new iteration
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_calc_estimate_htemp_for_new_iteration(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_calc_estimate_htemp_for_new_iteration
	!DEC$ endif
	class(ty_sat_calc),intent(inout)::this

	integer::i

	this%nodes%htemp = this%nodes%hnew

	end subroutine s_sat_calc_estimate_htemp_for_new_iteration

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to calculate results in nodes:
	!> - n%results_qent = n%qent
	!> - n%results_incvoldt = nrel·waterincmed·hsatinc/dtinc
	!> - n%results_qhor_all = kmed·hnew·dhnew/dx
	!> - n%results_dqhordx_all = get_derivatives(n%results_qhor_all)
	!> - n%results_dqhordx_from_incvoldt_all = n%results_qent-n%results_incvoldt
	!> - n%results_qhor = k1·h1·dhnew/dx
	!> - n%results_dqhordx = get_derivatives(n%results_qhor_all)
	!> - n%results_dqhordx_from_incvoldt = (h1·ksat1/(kmed·hnew))·n%results_dqhordx_from_incvoldt_all (it is the ponderated by permeability)
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_calc_results_in_nodes(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_calc_results_in_nodes
	!DEC$ endif
	use com_mod_ty_nodes,only:ty_com_nodes
	use sat_mod_ty_nodes,only:ty_sat_nodes
	use com_other_custom_functions,only:f_get_derivatives
	class(ty_sat_calc),intent(inout),target::this
	real(kind=dpd)::ksatlayers(this%layers%count+1),ksatlayers1(this%layers%count+1)

	class(ty_com_nodes),pointer::nodescom
	type(ty_sat_nodes),pointer::nodes
	integer:: i

	nodescom => this%nodes
	select type(nodescom)
	type is (ty_sat_nodes)
		nodes => nodescom
	end select


	ksatlayers(1:this%layers%count) = this%layers%material%ksat
	ksatlayers(this%layers%count+1) = 1.0_dpd
	ksatlayers1 = 0.0_dpd
	ksatlayers1(1) = this%layers%material(1)%ksat

	nodes%results_qent		= nodes%qent
	nodes%results_incvoldt	= this%nrel*this%layers%get_water_inc_med(nodes%hold,nodes%hnew)*(nodes%hnew-nodes%hold)/this%time%dt
	nodes%results_qhor_all		= -matmul(this%layers%get_inc_h_from0_vec(nodes%hnew),ksatlayers)*f_get_derivatives(nodes%x,nodes%hnew+nodes%z) !kmed·hnew·dhnew/dx
	nodes%results_dqhordx_all = f_get_derivatives(nodes%x,nodes%results_qhor_all)
	
	!Calculate the pressure drom in each node due to waterflow qsal and waterflow leaving each layer.
	!inch(i,l) = qsal·hl/kl+1/2·ql/(hl·kl)·(hl^2)  | ql = -kl·hl·dh/dx
	nodes%results_inchnew = nodes%results_qinf*matmul(this%layers%get_inc_h_from0_vec(nodes%hnew),1.0_dpd/ksatlayers)-1.0_dpd/2.0_dpd*matmul(this%layers%get_inc_h_from0_vec(nodes%hnew)**2.0_dpd,ksatlayers/ksatlayers)*f_get_derivatives(nodes%x,nodes%hnew+nodes%z) 

	nodes%results_qhor				= -matmul(this%layers%get_inc_h_from0_vec(nodes%hnew),ksatlayers1)*f_get_derivatives(nodes%x,nodes%hnew+nodes%z) !Only for layer 1
	nodes%results_dqhordx			= f_get_derivatives(nodes%x,nodes%results_qhor)

	nodes%results_dqhordx_from_incvoldt_all = nodes%results_qent-nodes%results_incvoldt
	where(nodes%hnew==0.0_dpd)
		nodes%results_dqhordx_from_incvoldt = nodes%results_dqhordx_from_incvoldt_all
	else where
		nodes%results_dqhordx_from_incvoldt = matmul(this%layers%get_inc_h_from0_vec(nodes%hnew),ksatlayers1)/matmul(this%layers%get_inc_h_from0_vec(nodes%hnew),ksatlayers)*nodes%results_dqhordx_from_incvoldt_all
	end where

	end subroutine s_sat_calc_results_in_nodes


	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to calculate results in elements
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_calc_results_in_elements(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_calc_results_in_elements
	!DEC$ endif
	use com_mod_fem_intelement,only:intelement_rel_1don
	use com_mod_ty_elements,only:ty_com_elements
	use sat_mod_ty_elements,only:ty_sat_elements
	use com_mod_ty_nodes,only:ty_com_nodes
	use sat_mod_ty_nodes,only:ty_sat_nodes
	use com_mod_fem_shapefunctions,only:dphi1d
	class(ty_sat_calc),intent(inout),target::this
	class(ty_com_elements),pointer::elemcom
	class(ty_sat_elements),pointer::elem
	class(ty_com_nodes),pointer::nodescom
	class(ty_sat_nodes),pointer::nodes
	type(ty_sat_elements)::tempelem
	integer::e
	integer,parameter::QUADRATURE_ORDER=40
	real(kind=dpd)::ksat(this%layers%count+1),ksat1(this%layers%count+1)
	real(kind=dpd)::dhzdx1,dhzdx0

	ksat(1:this%layers%count) = this%layers%material(:)%ksat
	ksat(this%layers%count+1) = 1.0_dpd !Assumed that permeatility over the top layer is 1m3/2
	ksat1=0.0_dpd
	ksat1(1) = this%layers%material(1)%ksat

	elemcom => this%elements
	nodescom => this%nodes

	select type(elemcom)
	type is (ty_sat_elements)
		elem => elemcom
	end select

	select type(nodescom)
	type is (ty_sat_nodes)
		nodes => nodescom
	end select
	do e=1,this%elements%count
		!Vertical flow entering element: Int(qv·dx)
		elem%results_qent(e)			=intelement_rel_1don(f_qent_chi		,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_qsat(e)			=intelement_rel_1don(f_qsat_chi		,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_qinf(e)			=intelement_rel_1don(f_qinf_chi		,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_inchnew(e)			=intelement_rel_1don(f_inchnew_chi		,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_incvoldt(e)	=intelement_rel_1don(f_incvoldt_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_dqhordx(e)		=intelement_rel_1don(f_dqhordx_chi	,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_dqhordx_from_incvoldt(e)				=intelement_rel_1don(f_dqhordx_from_incvoldt_chi			,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_dqhordx_all(e)									=intelement_rel_1don(f_dqhordx_all_chi								,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_dqhordx_from_incvoldt_all(e)		=intelement_rel_1don(f_dqhordx_from_incvoldt_all_chi	,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)

		elem%hnew(e)	=intelement_rel_1don(f_h_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%hold(e)	=intelement_rel_1don(f_hold_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%z(e)			=intelement_rel_1don(f_z_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%hz(e)		=intelement_rel_1don(f_hz_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%hzold(e)	=intelement_rel_1don(f_hzold_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%dhdx(e)	=intelement_rel_1don(f_dhdx_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%dhzdx(e)	=intelement_rel_1don(f_dhzdx_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)

		elem%x0(e)	= nodes%x(elem%idnode(e,1))
		elem%x1(e)	= nodes%x(elem%idnode(e,2))
		elem%h0(e)	= nodes%hnew(elem%idnode(e,1))
		elem%h1(e)	= nodes%hnew(elem%idnode(e,2))
		elem%z0(e)	= nodes%z(elem%idnode(e,1))
		elem%z1(e)	= nodes%z(elem%idnode(e,2))
		elem%dhdx0(e)	= dot_product(dphi1d(-1.0_dpd,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:)))
		elem%dhdx1(e)	= dot_product(dphi1d(1.0_dpd,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:)))
		elem%dhzdx0(e)	= dot_product(dphi1d(-1.0_dpd,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:))+nodes%z(elem%idnode(e,:)))
		elem%dhzdx1(e)	= dot_product(dphi1d(1.0_dpd,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:))+nodes%z(elem%idnode(e,:)))

		elem%q_all(e)	=intelement_rel_1don(f_q_all_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%q0_all(e)	= sum(f_q_all_chi((/-1.0_dpd/)))
		elem%q1_all(e)	= sum(f_q_all_chi((/1.0_dpd/)))

		elem%q(e)	=intelement_rel_1don(f_q_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%q0(e)	= sum(f_q_chi((/-1.0_dpd/)))
		elem%q1(e)	= sum(f_q_chi((/1.0_dpd/)))

		elem%k0(e)	= sum(f_k_chi((/-1.0_dpd/)))
		elem%k1(e)	= sum(f_k_chi((/1.0_dpd/)))
		elem%k(e)	=intelement_rel_1don(f_k_chi,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
	end do


	contains
	!----- qent
	function f_qent_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_qent_chi(size(chi))

	!Int(qent·phi)
	f_qent_chi = interp_on_element(chi,nodes%qent(elem%idnode(e,:)))
	end function f_qent_chi
		!----- qsat
	function f_qsat_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_qsat_chi(size(chi))

	!Int(qent·phi)
	f_qsat_chi = interp_on_element(chi,nodes%results_qsat(elem%idnode(e,:)))
	end function f_qsat_chi
	!----- qinf
	function f_qinf_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_qinf_chi(size(chi))

	!Int(qent·phi)
	f_qinf_chi = interp_on_element(chi,nodes%results_qinf(elem%idnode(e,:)))
	end function f_qinf_chi
	!----- inchnew
	function f_inchnew_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_inchnew_chi(size(chi))

	!Int(qent·phi)
	f_inchnew_chi = interp_on_element(chi,nodes%results_inchnew(elem%idnode(e,:)))
	end function f_inchnew_chi
	!----- h
	function f_h_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_h_chi(size(chi))
	!Int(qent·phi)
	f_h_chi = interp_on_element(chi,nodes%hnew(elem%idnode(e,:)))
	end function f_h_chi
	!----- hold
	function f_hold_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_hold_chi(size(chi))
	!Int(qent·phi)
	f_hold_chi = interp_on_element(chi,nodes%hold(elem%idnode(e,:)))
	end function f_hold_chi
	!----- z
	function f_z_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_z_chi(size(chi))
	!Int(qent·phi)
	f_z_chi = interp_on_element(chi,nodes%z(elem%idnode(e,:)))
	end function f_z_chi
	!----- hz
	function f_hz_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_hz_chi(size(chi))
	!Int(qent·phi)
	f_hz_chi = f_h_chi(chi)+f_z_chi(chi)
	end function f_hz_chi
	!------ hzold
	function f_hzold_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_hzold_chi(size(chi))
	!Int(qent·phi)
	f_hzold_chi = f_hold_chi(chi)+f_z_chi(chi)
	end function f_hzold_chi
	!----- dhdx
	function f_dhdx_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_dhdx_chi(size(chi))
	!Int(qent·phi)
	f_dhdx_chi = matmul(dphi1d(chi,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:)))
	end function f_dhdx_chi
	!----- dhzdx
	function f_dhzdx_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_dhzdx_chi(size(chi))
	!Int(qent·phi)
	f_dhzdx_chi = matmul(dphi1d(chi,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:))+nodes%z(elem%idnode(e,:)))
	end function f_dhzdx_chi
	!-----	incvoldt
	function f_incvoldt_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_incvoldt_chi(size(chi))
	real(kind=dps)::hnew_temp(size(chi)),hold_temp(size(chi)),nrel_temp(size(chi)),thsatres_temp(size(chi))
	!Int(qent·phi)

	hnew_temp = f_h_chi(chi)
	hold_temp = f_hold_chi(chi)
	nrel_temp = interp_on_element(chi,this%nrel(elem%idnode(e,:)))
	thsatres_temp = this%layers%get_water_inc_med(hold_temp,hnew_temp)
	f_incvoldt_chi = (hnew_temp-hold_temp)*nrel_temp*thsatres_temp/this%time%dt

	end function f_incvoldt_chi

	!------ dqhord (need nodes be calculated previously)
	function f_dqhordx_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_dqhordx_chi(size(chi))

	f_dqhordx_chi = interp_on_element(chi,nodes%results_dqhordx(elem%idnode(e,:)))

	end function f_dqhordx_chi

	!------ dqhordx_all (need nodes be calculated previously)
	function f_dqhordx_all_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_dqhordx_all_chi(size(chi))
	!Int(qent·phi)
	f_dqhordx_all_chi = interp_on_element(chi,nodes%results_dqhordx_all(elem%idnode(e,:)))
	end function f_dqhordx_all_chi

	!------ dqhordx_from_incvoldt_all
	function f_dqhordx_from_incvoldt_all_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_dqhordx_from_incvoldt_all_chi(size(chi))
	f_dqhordx_from_incvoldt_all_chi = f_qent_chi(chi)-f_incvoldt_chi(chi)
	end function f_dqhordx_from_incvoldt_all_chi

	!------ dqhordx_from_incvoldt
	function f_dqhordx_from_incvoldt_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_dqhordx_from_incvoldt_chi(size(chi))
	real(kind=dps)::flux_ratio_for_layer1(size(chi)),head(size(chi),this%layers%count+1)
	real(kind=dps)::incvoldt_chi(size(chi)),qent(size(chi))
	integer::i

	head = this%layers%get_inc_h_from0_vec(interp_on_element(chi,this%nodes%hnew(elem%idnode(e,:)))) !This is a matrix [nodesxlayers+1]

	do i=1,size(chi)
		if (sum(head(i,:))<1E-10_dpd) then
			flux_ratio_for_layer1(i) = 1.0_dpd
		else
			flux_ratio_for_layer1(i) = dot_product(head(i,:),ksat1)/dot_product(head(i,:),ksat)
		end if
	end do
	incvoldt_chi = f_incvoldt_chi(chi)
	qent = f_qent_chi(chi)
	f_dqhordx_from_incvoldt_chi = flux_ratio_for_layer1*(qent-incvoldt_chi)

	end function f_dqhordx_from_incvoldt_chi

	!------ k
	function f_k_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_k_chi(size(chi))
	real(kind=dpd)::head(size(chi),this%layers%count+1)
	integer::i

	head = this%layers%get_inc_h_from0_vec(interp_on_element(chi,this%nodes%hnew(elem%idnode(e,:))))

	do i=1,size(chi)
		if (sum(head(i,:))<1E-10_dpd) then
			f_k_chi(i) = ksat(1)
		else
			f_k_chi(i) = dot_product(head(i,:),ksat)/sum(head(i,:))
		end if
	end do
	end function f_k_chi

	!------ q_all
	function f_q_all_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_q_all_chi(size(chi))
	real(kind=dpd)::head(size(chi),this%layers%count+1)

	head = this%layers%get_inc_h_from0_vec(interp_on_element(chi,this%nodes%hnew(elem%idnode(e,:))))
	f_q_all_chi = -matmul(head,ksat)*f_dhzdx_chi(chi)

	end function f_q_all_chi

	!------ q
	function f_q_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_q_chi(size(chi))
	real(kind=dpd)::head(size(chi),this%layers%count+1)

	head = this%layers%get_inc_h_from0_vec(interp_on_element(chi,this%nodes%hnew(elem%idnode(e,:))))
	f_q_chi = -matmul(head,ksat1)*f_dhzdx_chi(chi)

	end function f_q_chi

	end subroutine s_sat_calc_results_in_elements

	end module sat_mod_ty_calc