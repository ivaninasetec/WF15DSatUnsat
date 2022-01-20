	!********************************************************************************************************************
	!*                                                                                                                  *
	!*                                 MODULE: MOD_COM_INTELEMENT                                                       *
	!*                                                                                                                  *
	!********************************************************************************************************************
	!* THIS MODULE INCLUDE THE FUNCTION THAT INTEGRATE OTHER FUNCTION (func1D) IN AN ELEMENT WITH NODES xne(nnodes)     *
	!* USING GAUSS-LEGENDRE QUADRATURE OF ORDER NO.                                                                     *
	!*                                                                                                                  *
	!*    Jacobian1D(chi,xne):   Calculate jacobian on the relative coord chi (from -1 to 1), of an element with        *
	!*                           absolute coords in nodes xne(nnodes).                                                  *
	!*                           The Jacobian is the determinant of the matrix of the derivatives of absolute coords,   *
	!*                           With respect the relative coords (chi) =|dx_i(chi)/dchi_j|.                            *
	!*                                                                                                                  *
	!* It uses the modules:                                                                                             *
	!*    MOD_COM_QUADRATURE                                                                                            *
	!*    MOD_COM_JACOBIAN                                                                                              *
	!*    MOD_COM_SHAPEFUNCTIONS                                                                                        *
	!*                                                                                                                  *
	!*    Iván Campos-Guereta Díez                                                                                      *
	!*    MSc Civil Engineering by Polytechnic University of Madrid                                                     *
	!*    PhD Student by University of Nottingham                                                                       *
	!*    eMBA by International Institute San Telmo in Seville                                                          *
	!*    ivan.camposguereta@nottingham.ac.uk                                                                           *
	!*                                                                                                                  *
	!*    This software is copyrighted 2019(C)                                                                          *
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
	!   func1d(x):         Is the relative coord at wich calculate the jacobian
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
	!   fun1d(chi):         Is the relative coord at wich calculate the jacobian
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
	!   fun1d(chi):         Is the relative coord at wich calculate the jacobian
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

	!
	!procedure(func_in_rel)::func1d

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