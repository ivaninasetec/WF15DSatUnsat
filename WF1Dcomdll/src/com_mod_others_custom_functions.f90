	!********************************************************************************************************************
	!        CLASS TO INCLUDE CUSTOM FUNCTIONS TO USE IN CALCULATIONS
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

	module com_other_custom_functions

	implicit none
	include 'inc_precision.fi'

	private


	public::f_get_derivatives
	contains

	function f_get_derivatives (x,y)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_get_derivatives" :: f_get_derivatives
	!DEC$ endif

	real(kind=dpd),intent(in)::x(:)
	real(kind=dpd),intent(in)::y(:)
	real(kind=dpd)::f_get_derivatives(size(x))
	integer::n

	n=size(x)

	f_get_derivatives(1)= (y(2)-y(1))/(x(2)-x(1)) + (y(2)-y(3))/(x(3)-x(2))+(y(3)-y(1))/(x(3)-x(1))
	f_get_derivatives(n)= (y(n-2)-y(n-1))/(x(n-1)-x(n-2)) + (y(n)-y(n-1))/(x(n)-x(n-1))+(y(n)-y(n-2))/(x(n)-x(n-2))

	f_get_derivatives(2:n-1) = ( (y(3:n)-y(2:n-1)) *(x(2:n-1)-x(1:n-2))**2+(y(2:n-1)-y(1:n-2))*(x(3:n)-x(2:n-1))**2 )/((x(2:n-1)-x(1:n-2))*(x(3:n)-x(2:n-1))*(x(3:n)-x(1:n-2)))


	end function f_get_derivatives




	end module com_other_custom_functions