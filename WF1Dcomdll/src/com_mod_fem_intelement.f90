	!********************************************************************************************************************
	! TITLE         : COM_MOD_FEM_INTELEMENT: FUNCTIONS TO INTEGRATE IN AN ELEMENT USING GAUSS-LEGENDRE QUARATURE
	! PROJECT       : WF1DCOMDLL
	! MODULE        : MOD_COM_QUADRATURE
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> this module include the function that integrate other function (func1d) in an element with nodes with absolute
	!> coordinates  xne(nnodes). Function func1d can be in absolute coordinates or in relative coordinates to the element
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module com_mod_fem_intelement

	implicit none
	private
	include 'inc_precision.fi'

	public:: intelement_abs_1don,intelement_rel_1don,intelement_abs_1don_sca

	contains

	!********************************************************************************************************************
	! F: INTELEMENT_ABS_1DON(FUNC1D,XNE,NO)
	!--------------------------------------------------------------------------------------------------------------------
	!		Function that returns the integral in the element. In this case the function is in absolute coords to the element.
	!   func1d(x):   Function to integrate (function of absolute coordinates)
	!   xne(nnodes): Is the absolute coodinates of the nodes of element.
	!********************************************************************************************************************

	function intelement_abs_1don_sca(func1d_x,xne,no)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"intelement_abs_1don_sca" :: intelement_abs_1don_sca
	!DEC$ endif

	use com_mod_fem_quadrature,only:        quadrature1d_sca
	use com_mod_fem_shapefunctions,only:          jacobian1d

	interface
	function func1d_x(x)
	import::dps
	real(kind=dps),intent(in)::x
	real(kind=dps)::func1d_x
	end function func1d_x
	end interface

	real(kind=dps),intent(in)::xne(:) !absolute x on nodes
	integer,intent(in)::no !quadrature order

	integer, parameter:: ndim=1				!number of dimensions=1
	integer, parameter:: nne=2				!number of nodes in element=2
	real(kind=dps)::intelement_abs_1don_sca

	intelement_abs_1don_sca = quadrature1d_sca(func_in_rel,no)

	contains
	function func_in_rel(chi) !function func1d·jacobian in relative coords (chi)
	use com_mod_fem_shapefunctions, only:   interp_on_element
	real(kind=dps),intent(in)::chi
	real(kind=dps)::func_in_rel

	func_in_rel = func1d_x(interp_on_element(chi,xne))*jacobian1d(chi,xne) !in this case jacobian returns 2.
	end function func_in_rel

	end function intelement_abs_1don_sca

	!********************************************************************************************************************
	! F: INTELEMENT_ABS_1DON(FUNC1D,XNE,NO)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the integral in the element. In this case the function is in relative coords to the element.
	!   fun1d(chi):  Function to integrate (function of absolute coordinates)
	!   xne(nnodes): Is the absolute coodinates of the nodes of element.
	!********************************************************************************************************************

	function intelement_abs_1don(func1d_x,xne,no)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"intelement_abs_1don_sca" :: intelement_abs_1don_sca
	!DEC$ endif

	use com_mod_fem_quadrature,only:        quadrature1d
	use com_mod_fem_shapefunctions,only:          jacobian1d

	interface
	function func1d_x(x)
	import::dps
	real(kind=dps),intent(in)::x(:)
	real(kind=dps)::func1d_x(size(x))
	end function func1d_x
	end interface

	real(kind=dps),intent(in)::xne(:) !absolute x on nodes
	integer,intent(in)::no !quadrature order

	integer, parameter:: ndim=1				!number of dimensions=1
	integer, parameter:: nne=2				!number of nodes in element=2
	real(kind=dps)::intelement_abs_1don

	intelement_abs_1don = quadrature1d(func_in_rel,no)

	contains
	function func_in_rel(chi) !function func1d·jacobian in relative coords (chi)
	use com_mod_fem_shapefunctions, only:   interp_on_element
	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::func_in_rel(size(chi))

	func_in_rel = func1d_x(interp_on_element(chi,xne))*jacobian1d(chi,xne) !in this case jacobian returns 2.
	end function func_in_rel

	end function intelement_abs_1don

	!********************************************************************************************************************
	! F: INTELEMENT_REL_1DON(FUNC1D,XNE,NO)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the integral in the element. In this case the function is in relative coords to the element.
	!   fun1d(chi):  Function to integrate (function of relative coordinates to the element)
	!   xne(nnodes): Is the absolute coodinates of the nodes of element.
	!********************************************************************************************************************

	function intelement_rel_1don(func1d,xne,no)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"intelement_rel_1don" :: intelement_rel_1don
	!DEC$ endif


	use com_mod_fem_quadrature,only:        quadrature1d
	use com_mod_fem_shapefunctions,only:          jacobian1d
	use com_mod_fem_shapefunctions, only:   interp_on_element

	interface
	function func1d(chi)
	import::dps
	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::func1d(size(chi))
	end function func1d
	end interface

	real(kind=dps),intent(in)::xne(:) !absolute x on nodes
	integer,intent(in)::no !quadrature order

	integer, parameter:: ndim=1				!number of dimensions=1
	integer, parameter:: nne=2				!number of nodes in element=2
	real(kind=dps)::intelement_rel_1don

	intelement_rel_1don = quadrature1d(func_in_rel,no)

	contains
	function func_in_rel(chi) !function func1d·jacobian in relative coords (chi)
	real(kind=dps),intent(in)::chi(:)
	real(kind=dps)::func_in_rel(size(chi))

	func_in_rel = func1d(chi)*jacobian1d(chi,xne) !in this case jacobian returns (x1-x0)/2.
	end function func_in_rel

	end function intelement_rel_1don

	end module com_mod_fem_intelement