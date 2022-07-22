	!********************************************************************************************************************
	! TITLE         : UNSAT_MOD_TY_CALC: EXTENDED DERIVED TYPE OF COM_MOD_TY_CALC TO INCLUDE PROPERTIES AND METHODS OF WF1DUNSAT
	! PROJECT       : WF1DUNSATDLL
	! MODULE        : COM_MOD_TY_CALC
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Extended derived type of com_mod_ty_elements to include properties and methods of wf1dunsat
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module unsat_mod_ty_calc
	use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded
	use com_mod_ty_calc,			only: ty_com_calc
	use com_mod_ty_material,	only: ty_com_material
	use com_mod_ty_layers,		only: ty_com_layers
	use unsat_mod_ty_elements, only:ty_unsat_elements
	use unsat_mod_ty_nodes	, only:ty_unsat_nodes

	implicit none
	include 'inc_precision.fi'

	private

	type,extends(ty_com_calc),public::ty_unsat_calc	!< CLASS: Definition of the layer in saturated model
		type(ty_com_layers)	,pointer::layers
		type(ty_com_material) ,pointer::material(:)
		type(ty_unsat_elements),pointer::elementsptr
		type(ty_unsat_nodes),pointer::nodesptr

		class(ty_mx),allocatable::mxmass1	!<Mass matrix: of the FEM linear system
		class(ty_mx),allocatable::mxmass2	!<Mass matrix: of the FEM linear system
		class(ty_mx),allocatable::mxstiff !<Stiffness matrix: (depend on results but not on time)
		class(ty_mx),allocatable::mxsink	!<Loads matrix: Loads on the FEM

		class(ty_mx),allocatable::mxinit	!<initial mx before changing by dirichtlet conditions
		real(kind=dpd),allocatable::rhsinit(:) !<initial rhs before changing dirichlet conditions

		real(kind=dpd),allocatable::colsink(:)
		real(kind=dpd),allocatable::colbound(:)


	contains
	procedure,public:: allocateall 							=> s_unsat_calc_allocateall
	procedure,public:: construct 								=> s_unsat_calc_construct
	procedure,public:: deallocateall 						=> s_unsat_calc_deallocateall
	procedure,public:: assign_newmann 					=> s_unsat_calc_assign_newmann
	procedure,public:: get_qnewmann 						=> s_unsat_get_qnewmann
	procedure,public:: update_th_from_h					=> s_unsat_calc_update_th_from_h
	procedure,public:: build_linearsystem 			=> s_unsat_calc_build_linearsystem
	procedure,public:: estimate_hnew_for_new_timestep	=> s_unsat_calc_estimate_hnew_for_new_timestep !Overriden
	procedure,public:: estimate_htemp_for_new_iteration	=> s_unsat_calc_estimate_htemp_for_new_iteration !Overriden
	procedure,public:: revert_to_old						=>	s_unsat_calc_revert_to_old !Overriden
	procedure,public:: set_old									=>	s_unsat_calc_set_old !Overriden
	procedure,public:: update_hnew_from_solution	=> s_unsat_calc_update_hnew_from_solution !Overriden
	procedure,public:: get_results_elements	=> s_unsat_calc_results_in_elements
	procedure,public:: get_results_nodes	=> s_unsat_calc_results_in_nodes
	procedure,public:: q_in_x => f_unsat_calc_q_in_x
	procedure,public:: q_in_element_chi => f_unsat_elem_q_in_chi_elem
	procedure,public:: get_thsatres_mean_from_x0_to_x1 => f_unsat_calc_get_thsatres_mean_from_x0_to_x1_sca
	procedure,public:: get_th_from_x0_to_x1 => f_unsat_calc_get_th_from_x0_to_x1_sca
	procedure,public:: get_thnewold_from_x0_to_x1 => f_unsat_calc_get_thnewold_from_x0_to_x1_sca
	end type ty_unsat_calc

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Get the value of thsat-thres on x0 to x1 using an integration order
	!> @param[inout] x0	First coordinate
	!> @param[inout] x1	Second coordinate
	!> @param[inout] integration_order	Quadrature order for the integral
	!---------------------------------------------------------------------------------------------------------------------

	function f_unsat_calc_get_thsatres_mean_from_x0_to_x1_sca(this,x0,x1,integration_order)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_unsat_calc_get_thsatres_mean_from_x0_to_x1_sca" :: f_unsat_calc_get_thsatres_mean_from_x0_to_x1_sca
	!DEC$ endif
	use com_mod_fem_intelement, only:intelement_abs_1don_sca

	class(ty_unsat_calc),intent(in)::this
	real(kind=dpd)::x0
	real(kind=dpd)::x1
	integer,intent(in)::integration_order
	real(kind=dpd)::f_unsat_calc_get_thsatres_mean_from_x0_to_x1_sca

	if (x0==x1) then
		f_unsat_calc_get_thsatres_mean_from_x0_to_x1_sca = f_thsatres_in_x(x0)
	else
		f_unsat_calc_get_thsatres_mean_from_x0_to_x1_sca = abs(intelement_abs_1don_sca(f_thsatres_in_x,(/x0,x1/)	,integration_order)/(x1-x0))
	end if

	!----------------------------------------------------------
	contains

	!This function returns thsat-thres in the height x
	function f_thsatres_in_x(x)

	real(kind=dpd),intent(in)::x
	real(kind=dpd)::f_thsatres_in_x
	integer::elementid_in_x

	elementid_in_x	= this%elements%id_from_x_sca(x)
	f_thsatres_in_x =	this%elements%material(elementid_in_x)%thsat-this%elements%material(elementid_in_x)%thres

	end function f_thsatres_in_x

	end function f_unsat_calc_get_thsatres_mean_from_x0_to_x1_sca


	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Get the value of th from x0 to x1 using an integration order
	!> @param[inout] x0	First coordinate
	!> @param[inout] x1	Second coordinate
	!> @param[inout] integration_order	Quadrature order for the integral
	!> @param[inout] option	0 use hnew, 1 use htemp, 2 use hold
	!---------------------------------------------------------------------------------------------------------------------

	function f_unsat_calc_get_th_from_x0_to_x1_sca(this,x0,x1,integration_order,option)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_unsat_calc_get_th_from_x0_to_x1_sca" :: f_unsat_calc_get_th_from_x0_to_x1_sca
	!DEC$ endif
	use com_mod_fem_intelement, only:intelement_abs_1don_sca

	class(ty_unsat_calc),intent(in)::this
	real(kind=dpd)::x0
	real(kind=dpd)::x1
	integer,intent(in)::integration_order
	integer,intent(in)::option
	real(kind=dpd)::f_unsat_calc_get_th_from_x0_to_x1_sca

	if (x0==x1) then
		f_unsat_calc_get_th_from_x0_to_x1_sca = 0.0_dpd
	else
		f_unsat_calc_get_th_from_x0_to_x1_sca = intelement_abs_1don_sca(f_th_in_x,(/x0,x1/)	,integration_order)
	end if

	!----------------------------------------------------------
	contains

	!This function returns th in the height x
	function f_th_in_x(x)
	use com_mod_linear_interpolation,only:linearinterpolation_ordered

	real(kind=dpd),intent(in)::x
	real(kind=dpd)::f_th_in_x
	integer::elementid_in_x
	real(kind=dpd)::hpres

	elementid_in_x = this%get_idelement_from_x(x)

	select case (option)
	case (0)
		hpres = linearinterpolation_ordered(x,this%nodes%x,this%nodes%hnew)
	case(1)
		hpres = linearinterpolation_ordered(x,this%nodes%x,this%nodes%htemp)
	case(2)
		hpres = linearinterpolation_ordered(x,this%nodes%x,this%nodes%hold)
	end  select

	f_th_in_x = this%elements%material(elementid_in_x)%get_th_sca(hpres)

	end function f_th_in_x

	end function f_unsat_calc_get_th_from_x0_to_x1_sca

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Get the value of thnew-thold (current timestep - previous) from x0 to x1 using an integration order
	!> @param[inout] x0	First coordinate
	!> @param[inout] x1	Second coordinate
	!> @param[inout] integration_order	Quadrature order for the integral
	!---------------------------------------------------------------------------------------------------------------------

	function f_unsat_calc_get_thnewold_from_x0_to_x1_sca(this,x0,x1,integration_order)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_unsat_calc_get_thnewold_from_x0_to_x1_sca" :: f_unsat_calc_get_thnewold_from_x0_to_x1_sca
	!DEC$ endif
	use com_mod_fem_intelement, only:intelement_abs_1don_sca

	class(ty_unsat_calc),intent(in)::this
	real(kind=dpd)::x0
	real(kind=dpd)::x1
	integer,intent(in)::integration_order
	real(kind=dpd)::f_unsat_calc_get_thnewold_from_x0_to_x1_sca

	if (x0==x1) then
		f_unsat_calc_get_thnewold_from_x0_to_x1_sca = 0.0_dpd
	else
		f_unsat_calc_get_thnewold_from_x0_to_x1_sca = intelement_abs_1don_sca(f_thnewold_in_x,(/x0,x1/)	,integration_order)
	end if

	!----------------------------------------------------------
	contains

	!This function returns thsat-thres in the height x
	function f_thnewold_in_x(x)
	use com_mod_linear_interpolation,only:linearinterpolation_ordered

	real(kind=dpd),intent(in)::x
	real(kind=dpd)::f_thnewold_in_x
	integer::elementid_in_x
	real(kind=dpd)::hnew,hold,thsatres_in_x

	elementid_in_x = this%get_idelement_from_x(x)
	hnew = linearinterpolation_ordered(x,this%nodes%x,this%nodes%hnew)
	hold = linearinterpolation_ordered(x,this%nodes%x,this%nodes%hold)

	if (elementid_in_x>this%elements%count) then
		f_thnewold_in_x = 1.0_dpd
	else
		thsatres_in_x =	this%elements%material(elementid_in_x)%thsat-this%elements%material(elementid_in_x)%thres

		f_thnewold_in_x = thsatres_in_x*this%elements%material(elementid_in_x)%get_incs_h1_to_h2_sca(hold,hnew)
	end if

	end function f_thnewold_in_x

	end function f_unsat_calc_get_thnewold_from_x0_to_x1_sca

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Estimate hnew for the new timestep
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_estimate_hnew_for_new_timestep(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_estimate_hnew_for_new_timestep" :: s_unsat_calc_estimate_hnew_for_new_timestep
	!DEC$ endif

	class(ty_unsat_calc),intent(inout)::this
	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDERHOLD=2,CONSIDER_ALL=3

	this%nodes%hnew	= this%nodes%hold
	call this%update_th_from_h(CONSIDER_HNEW)

	end subroutine s_unsat_calc_estimate_hnew_for_new_timestep

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Estimate htemp for the new iteration
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_estimate_htemp_for_new_iteration(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_estimate_htemp_for_new_iteration" :: s_unsat_calc_estimate_htemp_for_new_iteration
	!DEC$ endif
	class(ty_unsat_calc),intent(inout)::this
	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDERHOLD=2,CONSIDER_ALL=3

	this%nodes%htemp = this%nodes%hnew
	call this%update_th_from_h(CONSIDER_HTEMP)
	end subroutine s_unsat_calc_estimate_htemp_for_new_iteration

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Revert values to previous timestep (usually when the current timestep didn't converge)
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_revert_to_old(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_revert_to_old" :: s_unsat_calc_revert_to_old
	!DEC$ endif
	class(ty_unsat_calc),intent(inout)::this
	this%time%t			= this%time%told

	!All values to initial conditions
	this%nodes%hnew		= this%nodes%hold
	this%nodes%htemp	= this%nodes%hold

	this%nodes%thnew	= this%nodes%thold
	this%nodes%thtemp	= this%nodes%thold

	!All values to initial conditions
	this%nodes%dhnew	= this%nodes%dhold
	this%nodes%dhtemp	= this%nodes%dhold

	end subroutine s_unsat_calc_revert_to_old

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Set old values from current calculation (usually when iteration process has finished and converged)
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_set_old(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_set_old" :: s_unsat_calc_set_old
	!DEC$ endif
	class(ty_unsat_calc),intent(inout)::this

	if(this%time%told.ne.this%time%t) then
		this%time%dtold = (this%time%t-this%time%told)
		this%time%told		= this%time%t
	end if
	this%slopeold = min(1.0_dpd,(this%nodes%hnew-this%nodes%hold)/this%time%dtold)

	this%nodes%hold = this%nodes%hnew
	this%nodes%thold = this%nodes%thnew

	end subroutine s_unsat_calc_set_old

	!---------------------------------------------------------------------------------------------------------------------
	! UPDATE HNEW FROM SOLUTION AND GET EPSH
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Set value of hnew from the value of solution and coeficient of relaxation
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_update_hnew_from_solution(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_update_hnew_from_solution" :: s_unsat_calc_update_hnew_from_solution
	!DEC$ endif
	class(ty_unsat_calc),intent(inout)::this
	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDERHOLD=2,CONSIDER_ALL=3

	this%nodes%hnew = this%nodes%htemp + this%parameters%crelax*(this%solution-this%nodes%htemp)
	call this%update_th_from_h(CONSIDER_HTEMP)

	end subroutine s_unsat_calc_update_hnew_from_solution

	!---------------------------------------------------------------------------------------------------------------------
	! ASSIGN_NEWMAN
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Assing Newmman boundaries (waterflow on nodes)
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_assign_newmann(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_assign_newmann" :: s_unsat_calc_assign_newmann
	!DEC$ endif

	class(ty_unsat_calc),intent(inout)::this

	where(this%nodes%isNewmann)
		this%colbound = this%nodes%qnewmann
	else where
		this%colbound = 0.0_dpd
	end where

	end subroutine s_unsat_calc_assign_newmann

	!---------------------------------------------------------------------------------------------------------------------
	! GET QNEWMANN DUE TO DIRICHLET
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Get the values of the flow due to the imposition of the dirichlet boundary conditions
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_get_qnewmann(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_get_qnewmann" :: s_unsat_get_qnewmann
	!DEC$ endif

	class(ty_unsat_calc),intent(inout)::this

	this%nodes%qnewmann = this%mx*this%nodes%hnew-this%rhs+this%colbound

	end subroutine s_unsat_get_qnewmann

	!---------------------------------------------------------------------------------------------------------------------
	! GET TH FROM H
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Get the water content on the nodes knowing the value of hnew, htemp or hold
	!> @param[in] option 0 use hnew, 1 use htemp, 2 use hold
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_update_th_from_h(this,option)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_update_th_from_h" :: s_unsat_calc_update_th_from_h
	!DEC$ endif
	use unsat_mod_ty_nodes,only:ty_unsat_nodes
	use com_mod_ty_nodes,only:ty_com_nodes

	class(ty_unsat_calc),intent(inout),target::this
	integer,intent(in),optional::option
	class(ty_com_nodes),pointer::nodescom
	class(ty_unsat_nodes),pointer::nodes
	integer,parameter						::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDER_HOLD=2,CONSIDER_ALL=3

	integer::opt,i

	nodescom=>this%nodes
	opt = 0

	if(present(option)) opt = option

	select type(nodescom)
	type is(ty_unsat_nodes)
		nodes=>nodescom
	end select

	select case (opt)
	case(CONSIDER_HOLD)
		!$OMP DO
		do i=1,size(nodes%thold)
			nodes%thold(i)	= this%nodes%material(i)%Get_th_sca(this%nodes%hold(i))
		end do
		!$OMP END DO
	case(CONSIDER_HTEMP)
		!$OMP DO
		do i=1,size(nodes%thtemp)
			nodes%thtemp(i)	= this%nodes%material(i)%Get_th_sca(this%nodes%htemp(i))
		end do
		!$OMP END DO
	case(CONSIDER_ALL)
		!$OMP DO
		do i=1,size(nodes%thnew)
			nodes%thold(i)	= this%nodes%material(i)%Get_th_sca(this%nodes%hold(i))
			nodes%thtemp(i)	= this%nodes%material(i)%Get_th_sca(this%nodes%htemp(i))
			nodes%thnew(i)	= this%nodes%material(i)%Get_th_sca(this%nodes%hnew(i))
		end do
		!$OMP END DO
		case default
		!$OMP DO
		do i=1,size(nodes%thnew)
			nodes%thnew(i)	= this%nodes%material(i)%Get_th_sca(this%nodes%hnew(i))
		end do
		!$OMP END DO
	end select

	end subroutine s_unsat_calc_update_th_from_h

	!---------------------------------------------------------------------------------------------------------------------
	! ALLOCATE ALL (OVERRIDE)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all allocatable properties of the parent class
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_unsat_calc_allocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_allocateall" :: s_unsat_calc_allocateall
	!DEC$ endif
	use unsat_mod_ty_nodes, only:ty_unsat_nodes
	use unsat_mod_ty_elements, only:ty_unsat_elements
	use com_mod_ty_parameters, only:ty_com_parameters

	class(ty_unsat_calc),intent(inout)::this

	if(.not.allocated(this%nodes))			allocate(ty_unsat_nodes::this%nodes)
	if(.not.allocated(this%elements))		allocate(ty_unsat_elements::this%elements)

	end subroutine s_unsat_calc_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	! DEALLOCATE ALL (OVERRIDE)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate all allocatable properties of the parent class
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_deallocateall" :: s_unsat_calc_deallocateall
	!DEC$ endif

	class(ty_unsat_calc),intent(inout)::this

	if(allocated(this%slopeold))		deallocate(this%slopeold)
	if(allocated(this%rhs))					deallocate(this%rhs)
	if(allocated(this%solution))		deallocate(this%solution)
	if(allocated(this%perm))				deallocate(this%perm)

	if(allocated(this%mx))					deallocate(this%mx)
	if(allocated(this%mxmass1))			deallocate(this%mxmass1)
	if(allocated(this%mxmass2))			deallocate(this%mxmass2)
	if(allocated(this%mxstiff))			deallocate(this%mxstiff)
	if(allocated(this%mxsink))			deallocate(this%mxsink)

	if(allocated(this%elements)) then
		call this%elements%deallocateall()
		deallocate(this%elements)
	end if

	if(allocated(this%nodes)) then
		call this%nodes%deallocateall()
		deallocate(this%nodes)
	end if

	end subroutine s_unsat_calc_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	! CONSTRUCTOR (OVERRIDE)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Construct the class (asigning all asignable properties)
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_unsat_calc_construct(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_construct" :: s_unsat_calc_construct
	!DEC$ endif
	use com_mod_ty_elements,	only: ty_com_elements
	use unsat_mod_ty_elements,	only: ty_unsat_elements
	use com_mod_ty_nodes,	only: ty_com_nodes
	use unsat_mod_ty_nodes,	only: ty_unsat_nodes
	integer,parameter::KU=1,KL=1

	class(ty_unsat_calc),intent(inout),target::this
	class(ty_com_elements),pointer::elemcom
	class(ty_unsat_elements),pointer::elem
	class(ty_com_nodes),pointer::nodescom
	class(ty_unsat_nodes),pointer::nodes
	integer::nnc !<Number of nodes-classes

	nnc=this%nodes%count
	elemcom=>this%elements
	select type(elemcom)
	type is(ty_unsat_elements)
		elem => elemcom
	end select

	nodescom=>this%nodes
	select type(nodescom)
	type is(ty_unsat_nodes)
		nodes => nodescom
	end select

	this%nodesptr=> nodes
	this%elementsptr=> elem

	!Allocate all vectors
	if(.not.allocated(this%colbound)) allocate(this%colbound(nnc))
	if(.not.allocated(this%colsink)) allocate(this%colsink(nnc))
	if(.not.allocated(this%slopeold)) allocate(this%slopeold(nnc))
	if(.not.allocated(this%rhs))			allocate(this%rhs(nnc))
	if(.not.allocated(this%solution)) allocate(this%solution(nnc))
	if(.not.allocated(this%perm))			allocate(this%perm(nnc)) !If perm is not allocated then it is not used.
	if(.not.allocated(this%rhsinit))			allocate(this%rhsinit(nnc))

	if(.not.allocated(elem%thnew))		allocate(elem%thnew(this%elements%count))
	if(.not.allocated(elem%thtemp))		allocate(elem%thtemp(this%elements%count))
	if(.not.allocated(elem%thold))		allocate(elem%thold(this%elements%count))
	if(.not.allocated(elem%results_qent))			allocate(elem%results_qent(this%elements%count))
	if(.not.allocated(elem%results_qsal))			allocate(elem%results_qsal(this%elements%count))
	if(.not.allocated(elem%results_incvoldtperm)) allocate(elem%results_incvoldtperm(this%elements%count))
	if(.not.allocated(elem%results_qmed))			allocate(elem%results_qmed(this%elements%count))
	if(.not.allocated(elem%results_kmed))			allocate(elem%results_kmed(this%elements%count))
	if(.not.allocated(elem%results_dhdxmed))			allocate(elem%results_dhdxmed(this%elements%count))
	if(.not.allocated(elem%results_dhxdxmed))			allocate(elem%results_dhxdxmed(this%elements%count))
	if(.not.allocated(elem%h0))			allocate(elem%h0(this%elements%count))
	if(.not.allocated(elem%h1))			allocate(elem%h1(this%elements%count))
	if(.not.allocated(elem%h0old))			allocate(elem%h0old(this%elements%count))
	if(.not.allocated(elem%h1old))			allocate(elem%h1old(this%elements%count))
	if(.not.allocated(elem%th0))			allocate(elem%th0(this%elements%count))
	if(.not.allocated(elem%th1))			allocate(elem%th1(this%elements%count))
	if(.not.allocated(elem%k0))			allocate(elem%k0(this%elements%count))
	if(.not.allocated(elem%k1))			allocate(elem%k1(this%elements%count))
	if(.not.allocated(elem%dhdx0))			allocate(elem%dhdx0(this%elements%count))
	if(.not.allocated(elem%dhdx1))			allocate(elem%dhdx1(this%elements%count))
	if(.not.allocated(elem%dhxdx0))			allocate(elem%dhxdx0(this%elements%count))
	if(.not.allocated(elem%dhxdx1))			allocate(elem%dhxdx1(this%elements%count))
	elem%results_qent=0.0_dpd
	elem%results_qsal=0.0_dpd
	elem%results_qmed=0.0_dpd
	elem%results_incvoldtperm=0.0_dpd

	if(.not.allocated(nodes%thnew))		allocate(nodes%thnew(nodes%count))
	if(.not.allocated(nodes%thtemp))	allocate(nodes%thtemp(nodes%count))
	if(.not.allocated(nodes%thold))		allocate(nodes%thold(nodes%count))
	if(.not.allocated(nodes%results_qnewmann)) allocate(nodes%results_qnewmann(nodes%count))
	nodes%results_qnewmann = 0.0_dpd

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
	call this%allocate_matrix(this%mxmass1)
	call this%allocate_matrix(this%mxmass2)
	call this%allocate_matrix(this%mxstiff)
	call this%allocate_matrix(this%mxsink)
	call this%allocate_matrix(this%mxinit)

	this%slopeold = 0.0_dpd

	end subroutine s_unsat_calc_construct

	!---------------------------------------------------------------------------------------------------------------------
	! BUILD LINEAR SYSTEM
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Build coefficient matrix for actual head pressures.
	!> @param[in] IsTimeDependant If true, the matrix is time dependant.
	!> @param[in] option CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDER_HOLD=2,CONSIDER_ALL=3
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_build_linearsystem(this,IsTimeDependant,option)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_build_linearsystem" :: s_unsat_calc_build_linearsystem
	!DEC$ endif
	use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded
	use com_mod_ty_nodes, only:ty_com_nodes
	use unsat_mod_ty_nodes, only: ty_unsat_nodes

	integer,parameter::COEF_MASS=1,COEF_STIFF=2,COEF_LOAD=3,COEF_BOUND=4
	integer,parameter::OPTIONBASIS_H=1,OPTIONBASIS_DH=2,OPTIONBASIS_H_H=3,OPTIONBASIS_H_DH=4,OPTIONBASIS_DH_H=5,OPTIONBASIS_DH_DH=6
	integer,parameter::CONSIDER_HNEW=0,CONSIDER_HTEMP=1,CONSIDER_HOLD=2,CONSIDER_ALL=3
	class(ty_unsat_calc),intent(inout),target::this

	logical,intent(in)::IsTimeDependant
	integer,intent(in),optional::option
	integer::opt,nelem,numnod,numclass,numnodclass
	class(ty_com_nodes),pointer::nodescom
	class(ty_unsat_nodes),pointer::nodes

	real(kind=dpd)::henodes(this%elements%nn*(this%elements%nc+1))
	real(kind=dpd)::henodesold(this%elements%nn*(this%elements%nc+1))

	nodescom => this%nodes
	select type(nodescom)
	type is(ty_unsat_nodes)
		nodes => nodescom
	end select


	opt=CONSIDER_HNEW
	if(present(option)) opt = option

	nelem				=this%elements%count
	numnod			=this%elements%nn
	numclass		=this%elements%nc+1
	numnodclass =numnod*numclass

	if (this%parameters%isModifiedPicard) call this%update_th_from_h(CONSIDER_ALL)
	!tex:Mass1 matrix:
	!$\mathop \smallint \limits_{\rm{\Omega }} {\left. {\frac{{d\theta }}{{d\psi }}} \right|^{k,n + 1}}{\phi _j}{\phi _i}{\rm{\;}}d{\rm{\Omega }}$
	if(IsTimeDependant)				call this%buildcoefmatrix(f_mass1,this%mxmass1,this%parameters%masslump,OPTIONBASIS_H_H)				!Build [MASS1]
	!tex:Mass2 matrix:
	!$\mathop \smallint \limits_{\rm{\Omega }} {\phi _j}{\phi _i}{\rm{\;}}d{\rm{\Omega }}$]
	if(.not.IsTimeDependant)	call this%buildcoefmatrix(f_identity,this%mxmass2,this%parameters%masslump,OPTIONBASIS_H_H)			!Build [MASS2]
	!tex:Stiffness matrix:
	!$\mathop \smallint \limits_{\rm{\Omega }} {k^{k,n + 1}}\frac{{\partial {\phi _j}}}{{\partial x}}\frac{{\partial {\phi _i}}}{{\partial x}}{\rm{\;}}d{\rm{\Omega }}$
	if(IsTimeDependant)				call this%buildcoefmatrix(f_stiff,this%mxstiff,.false.,OPTIONBASIS_DH_DH)												!Build [STIFF]
	!tex:Sink matrix:
	!$\mathop \smallint \limits_{\rm{\Omega }} {\phi _j}{\phi _i}{\rm{\;}}d{\rm{\Omega }}$]
	if(.not.IsTimeDependant)	call this%buildcoefmatrix(f_identity,this%mxsink,.false.,OPTIONBASIS_H_H)			!Build [SINK]

	if (isTimeDependant) then
		!create mx matrix: [mx]=[mass]/dt+[stiff]+[bound]
		!call construct_mx_dt(dt,mx,param%typesolution)

		if (this%parameters%isModifiedPicard) then
			!Linerization: Celia, 1990 Modified Picard linearization...
			![MX] = [MASS1]/dt+[STIFF]
			!{rhs}=[MASS1]{htemp}/dt-[MASS2]{thtemp-thold}/dt-[STIFF]{z}+[SINK]{sink}+{COLBOUND}

			!tex: Modified Picard linearization:
			!$\left( {\frac{1}{{{\rm{\Delta }}t}}\left[ {{M_{MASSi,j}}} \right] + \left[ {{M_{Ki,j}}} \right]} \right)\left\{ {\psi _j^{k + 1,n + 1}} \right\} = \frac{1}{{{\rm{\Delta }}t}}\left[ {{M_{MASSi,j}}} \right]\left\{ {\psi _j^{k,n + 1}} \right\} - \frac{1}{{{\rm{\Delta }}t}}\left[ {{I_{i,j}}} \right]\left\{ {\theta _j^{k,n + 1} - \theta _j^n} \right\} - \left[ {{M_{Ki,j}}} \right]\left\{ {{z_j}} \right\} + \left\{ {{Q_{bound,i}}} \right\}$

			this%mx = this%mxmass1/this%time%dt+this%mxstiff
			this%rhs=	this%mxmass1*nodes%htemp/this%time%dt-this%mxmass2*(nodes%thtemp-nodes%thold)/this%time%dt-this%mxstiff*nodes%x+this%mxsink*this%colsink+this%colbound


		ELSE
			!Linearization: simple Picard linearization:
			![MX] = [MASS1]/dt+[STIFF]
			!{RHS}=[MASS1]·{hold}/dt-[STIFF]{z}+[SINK]{sink}+{colbound}

			!tex: Picard linearization
			!$\left( {\frac{1}{{{\rm{\Delta }}t}}\left[ {{M_{MASSi,j}}} \right] + \left[ {{M_{,j}}} \right]} \right)\left\{ {\psi _j^{k + 1,n + 1}} \right\} = \frac{1}{{{\rm{\Delta }}t}}\left[ {{M_{MASSi,j}}} \right]\left\{ {\psi _j^n} \right\} - \left[ {{M_{Ki,j}}} \right]\left\{ {{z_j}} \right\} + \left\{ {{Q_{bound,i}}} \right\}$
			this%mx = this%mxmass1/this%time%dt+this%mxstiff
			this%rhs=this%mxmass1*this%nodes%hold/this%time%dt-this%mxstiff*this%nodes%x+this%mxsink*this%colsink+this%colbound

		end if

	end if

	this%mxinit = this%mx
	this%rhsinit = this%rhs

	contains
	!--------------------------------------------------------------------------
	!Functions to integrate in each kind of matrix (for [LOAD] matrix)
	function f_identity(chi,e) !For example for mass2
	real(kind=dps),intent(in)::chi(:)
	integer,intent(in)::e
	real(kind=dps)::f_identity(size(chi))
	f_identity = 1.0_dps
	end function f_identity
	!-----------------------------------------------------------------------------------------
	!Returns the increment of water content from hold to htemp in coord chi of element e (for [MASS] and [SINK] matrixes)
	function f_mass1(chi,e)
	use com_mod_fem_shapefunctions,only:interp_on_element
	real(kind=dps),intent(in)::chi(:)
	integer,intent(in)::e
	real(kind=dps)::f_mass1(size(chi))
	real(kind=dps)::headtemp(size(chi)),material(size(chi))

	integer::nstartelem,nendelem

	nstartelem	= this%elements%idnode(e,1)
	nendelem		= this%elements%idnode(e,numnodclass)
	select case (opt)
	case (CONSIDER_HNEW)
		henodes = this%nodes%hnew(nstartelem : nendelem)
	case (CONSIDER_HTEMP)
		henodes = this%nodes%htemp(nstartelem : nendelem)
		case default
		stop('Matrix cannot be builded from old values in nodes')
	end select

	!tex:
	!$\mathop \smallint \limits_{\rm{\Omega }} {\left. {\frac{{d\theta }}{{d\psi }}} \right|^{k,n + 1}}{\phi _j}{\phi _i}{\rm{\;}}d{\rm{\Omega }}$
	headtemp = interp_on_element(chi,henodes)			!Returns htemp in each chi
	f_mass1 = this%elements%material(e)%get_cap_vec(headtemp)

	end function f_mass1

	!-----------------------------------------------------------------------------------------
	!Returns the increment of water content from hold to htemp in coord chi of element e (for [MASS] matrix)
	function f_stiff(chi,e)
	use com_mod_fem_shapefunctions,only:interp_on_element
	real(kind=dps),intent(in)::chi(:)
	integer,intent(in)::e
	real(kind=dps)::f_stiff(size(chi))
	real(kind=dps)::headtemp(size(chi)),material(size(chi))

	integer::nstartelem,nendelem

	nstartelem	= this%elements%idnode(e,1)
	nendelem		= this%elements%idnode(e,numnodclass)
	select case (opt)
	case (CONSIDER_HNEW)
		henodes = this%nodes%hnew(nstartelem : nendelem)
	case (CONSIDER_HTEMP)
		henodes = this%nodes%htemp(nstartelem : nendelem)
		case default
		stop('Matrix cannot be builded from old values in nodes')
	end select

	!tex:Stiffneww matrix:
	!$\mathop \smallint \limits_{\rm{\Omega }} {k^{k,n + 1}}\frac{{\partial {\phi _j}}}{{\partial x}}\frac{{\partial {\phi _i}}}{{\partial x}}{\rm{\;}}d{\rm{\Omega }}$
	headtemp = interp_on_element(chi,henodes)			!Returns htemp in each chi

	f_stiff = this%elements%material(e)%get_k_vec(headtemp)

	end function f_stiff

	end subroutine s_unsat_calc_build_linearsystem

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure calculate the results in the nodes
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_results_in_nodes(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_results_in_nodes" :: s_unsat_calc_results_in_nodes
	!DEC$ endif
	use com_mod_ty_nodes,only:ty_com_nodes
	use unsat_mod_ty_nodes,only:ty_unsat_nodes

	class(ty_unsat_calc),intent(inout),target::this
	class(ty_com_nodes),pointer::nodescom
	class(ty_unsat_nodes),pointer::nodes

	nodescom => this%nodes

	select type(nodescom)
	type is(ty_unsat_nodes)
		nodes=>nodescom
	end select

	nodes%results_qnewmann = -(this%mxinit*nodes%hnew-(this%rhsinit-this%colbound)) !Use init, before changing mx and rhs by dirichlet conditions.


	end subroutine s_unsat_calc_results_in_nodes

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure to calculate the results in the elements
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_calc_results_in_elements(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_calc_results_in_elements" :: s_unsat_calc_results_in_elements
	!DEC$ endif
	use com_mod_fem_intelement,only:intelement_rel_1don
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements
	use com_mod_ty_nodes,only:ty_com_nodes
	use unsat_mod_ty_nodes,only:ty_unsat_nodes
	use com_mod_fem_shapefunctions,only:dphi1d

	class(ty_unsat_calc),intent(inout),target::this
	class(ty_com_elements),pointer::elemcom
	class(ty_unsat_elements),pointer::elem
	class(ty_com_nodes),pointer::nodescom
	class(ty_unsat_nodes),pointer::nodes
	integer::e
	real(kind=dpd)::dhzdx1,dhzdx0
	integer,parameter::QUADRATURE_ORDER=40


	elemcom => this%elements
	select type(elemcom)
	type is (ty_unsat_elements)
		elem => elemcom
	end select

	nodescom => this%nodes
	select type(nodescom)
	type is (ty_unsat_nodes)
		nodes => nodescom
	end select

	do e=1,this%elements%count
		elem%results_qent(e) = elem%material(e)%get_k_sca(nodes%hnew(elem%idnode(e,2)))*dot_product(nodes%hnew(elem%idnode(e,:))+nodes%x(elem%idnode(e,:)),dphi1d(1.0_dpd,nodes%x(elem%idnode(e,:))))
		elem%results_qsal(e) = elem%material(e)%get_k_sca(nodes%hnew(elem%idnode(e,1)))*dot_product(nodes%hnew(elem%idnode(e,:))+nodes%x(elem%idnode(e,:)),dphi1d(-1.0_dpd,nodes%x(elem%idnode(e,:))))

		elem%results_qmed(e)			=intelement_rel_1don(f_qmed_chi		,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)

		elem%hnew(e)							=intelement_rel_1don(f_hnew_chi	,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%htemp(e)							=intelement_rel_1don(f_htemp_chi	,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%hold(e)							=intelement_rel_1don(f_hold_chi	,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)

		elem%thnew(e)							=intelement_rel_1don(f_thnew_chi	,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%thtemp(e)						=intelement_rel_1don(f_thtemp_chi	,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%thold(e)							=intelement_rel_1don(f_thold_chi	,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_incvoldtperm(e)	= (elem%thnew(e)-elem%thold(e))/this%time%dt
		elem%h0(e)							= nodes%hnew(elem%idnode(e,1))
		elem%h1(e)							= nodes%hnew(elem%idnode(e,2))
		elem%h0old(e)						= nodes%hold(elem%idnode(e,1))
		elem%h1old(e)						= nodes%hold(elem%idnode(e,2))
		elem%th0(e)							= elem%material(e)%get_th_sca(nodes%hnew(elem%idnode(e,1)))
		elem%th1(e)							= elem%material(e)%get_th_sca(nodes%hnew(elem%idnode(e,2)))
		elem%k0(e)							= elem%material(e)%get_k_sca(nodes%hnew(elem%idnode(e,1)))
		elem%k1(e)							= elem%material(e)%get_k_sca(nodes%hnew(elem%idnode(e,2)))
		elem%dhdx0(e)						= dot_product(dphi1d(1.0_dpd,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:)))
		elem%dhdx1(e)						= dot_product(dphi1d(-1.0_dpd,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:)))
		elem%dhxdx0(e)					= dot_product(dphi1d(1.0_dpd,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:))+nodes%x(elem%idnode(e,:)))
		elem%dhxdx1(e)					= dot_product(dphi1d(-1.0_dpd,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:))+nodes%x(elem%idnode(e,:)))
		elem%results_dhdxmed(e)		=intelement_rel_1don(f_dhdxmed_chi		,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_dhxdxmed(e)	=intelement_rel_1don(f_dhxdxmed_chi		,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
		elem%results_kmed(e)			=intelement_rel_1don(f_kmed_chi		,this%nodes%x(elem%idnode(e,:))	, QUADRATURE_ORDER)/elem%lenght(e)
	end do

	contains
	!-----
	function f_qmed_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element,dphi1d

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_qmed_chi(size(chi))
	real(kind=dps)::ktemp(size(chi)),hnewtemp(size(chi)),dphi1dtemp(size(chi),size(elem%idnode(e,:))),dhdxtemp(size(chi)),hxtemp(size(elem%idnode(e,:)))

	hnewtemp = interp_on_element(chi,nodes%hnew(elem%idnode(e,:)))
	ktemp = 	elem%material(e)%get_k_vec(hnewtemp)
	dphi1dtemp = dphi1d(chi,nodes%x(elem%idnode(e,:)))
	hxtemp = nodes%hnew(elem%idnode(e,:))+nodes%x(elem%idnode(e,:))
	dhdxtemp = matmul(dphi1dtemp,hxtemp)
	f_qmed_chi = ktemp*dhdxtemp

	end function f_qmed_chi
	!-----
	function f_kmed_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_kmed_chi(size(chi))

	f_kmed_chi = elem%material(e)%get_k_vec(interp_on_element(chi,nodes%hnew(elem%idnode(e,:))))

	end function f_kmed_chi

	!-----
	function f_dhdxmed_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element,dphi1d

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_dhdxmed_chi(size(chi))

	f_dhdxmed_chi = matmul(dphi1d(chi,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:)))

	end function f_dhdxmed_chi
	!-----
	function f_dhxdxmed_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element,dphi1d

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_dhxdxmed_chi(size(chi))

	f_dhxdxmed_chi = matmul(dphi1d(chi,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:))+nodes%x(elem%idnode(e,:)))

	end function f_dhxdxmed_chi
	!------
	function f_thnew_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_thnew_chi(size(chi))

	f_thnew_chi = elem%material(e)%get_th_vec(interp_on_element(chi,nodes%hnew(elem%idnode(e,:))))

	end function f_thnew_chi
	!------
	function f_thtemp_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_thtemp_chi(size(chi))

	f_thtemp_chi = elem%material(e)%get_th_vec(interp_on_element(chi,nodes%htemp(elem%idnode(e,:))))

	end function f_thtemp_chi

	!------
	function f_thold_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_thold_chi(size(chi))

	f_thold_chi = elem%material(e)%get_th_vec(interp_on_element(chi,nodes%hold(elem%idnode(e,:))))

	end function f_thold_chi

	!---------------

	function f_hnew_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_hnew_chi(size(chi))

	f_hnew_chi = interp_on_element(chi,nodes%hnew(elem%idnode(e,:)))

	end function f_hnew_chi
	!------
	function f_htemp_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_htemp_chi(size(chi))

	f_htemp_chi = interp_on_element(chi,nodes%htemp(elem%idnode(e,:)))

	end function f_htemp_chi

	!------
	function f_hold_chi(chi)
	use com_mod_fem_shapefunctions,only:interp_on_element

	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::f_hold_chi(size(chi))

	f_hold_chi = interp_on_element(chi,nodes%hold(elem%idnode(e,:)))

	end function f_hold_chi



	end subroutine s_unsat_calc_results_in_elements

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Calculate the waterflow in absolute coordinate x
	!> @param[in] x coordinate where to calculate wht waterflow
	!> @param[in] el optional parameter to indicate the number of the element where to calculate the waterflow (where x is located)
	!---------------------------------------------------------------------------------------------------------------------

	function f_unsat_calc_q_in_x(this,x,el)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_unsat_calc_q_in_x" :: f_unsat_calc_q_in_x
	!DEC$ endif
	use com_mod_ty_elements,only:ty_com_elements
	use unsat_mod_ty_elements,only:ty_unsat_elements
	use com_mod_ty_nodes,only:ty_com_nodes
	use unsat_mod_ty_nodes,only:ty_unsat_nodes
	use com_mod_fem_shapefunctions,only:dphi1d
	use com_mod_fem_shapefunctions,only:interp_on_element

	class(ty_unsat_calc),intent(in),target::this
	real(kind=dpd),intent(in)::x
	integer,intent(in),optional::el
	real(kind=dpd)::f_unsat_calc_q_in_x
	class(ty_com_elements),pointer::elemcom
	class(ty_unsat_elements),pointer::elem
	class(ty_com_nodes),pointer::nodescom
	class(ty_unsat_nodes),pointer::nodes
	real(kind=dpd)::chi,dhzdx,k
	integer::e

	elemcom => this%elements
	nodescom => this%nodes

	select type(elemcom)
	type is (ty_unsat_elements)
		elem=> elemcom
	end select

	select type(nodescom)
	type is (ty_unsat_nodes)
		nodes=> nodescom
	end select

	if(present(el)) then
		e= el
	else
		e = this%get_idelement_from_x(x)
	end if

	chi = 2.0_dpd*(x-elem%x0(e))/(elem%x1(e)-elem%x0(e))-1.0_dpd

	dhzdx = dot_product(dphi1d(chi,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:))+nodes%x(elem%idnode(e,:)))

	k = elem%material(e)%get_k_sca(interp_on_element(chi,nodes%hnew(elem%idnode(e,:))))
	f_unsat_calc_q_in_x = k*dhzdx


	end function f_unsat_calc_q_in_x

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Return the waterflow at the relative coodinate chi of the element e
	!> @param[inout] chi	Relative coordinate where to get the waterflow
	!> @param[inout] e		Element where to calculate
	!---------------------------------------------------------------------------------------------------------------------

	function f_unsat_elem_q_in_chi_elem(this,chi,e)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_unsat_elem_q_in_chi_elem" :: f_unsat_elem_q_in_chi_elem
	!DEC$ endif
	use com_mod_ty_elements,only:ty_com_elements
	use com_mod_ty_nodes,only:ty_com_nodes
	use com_mod_fem_shapefunctions,only:dphi1d
	use com_mod_fem_shapefunctions,only:interp_on_element

	class(ty_unsat_calc),intent(in),target::this
	real(kind=dpd),intent(in)::chi
	integer,intent(in)::e !id of element
	class(ty_com_nodes),pointer::nodes
	class(ty_com_elements),pointer::elem
	real(kind=dpd)::f_unsat_elem_q_in_chi_elem

	real(kind=dpd)::dhzdx,k

	elem => this%elements
	nodes => this%nodes

	dhzdx =  dot_product(dphi1d(chi,nodes%x(elem%idnode(e,:))),nodes%hnew(elem%idnode(e,:))+nodes%x(elem%idnode(e,:)))
	k = elem%material(e)%get_k_sca(interp_on_element(chi,nodes%hnew(elem%idnode(e,:))))
	f_unsat_elem_q_in_chi_elem = k*dhzdx

	end function f_unsat_elem_q_in_chi_elem


	end module unsat_mod_ty_calc