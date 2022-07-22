	!********************************************************************************************************************
	! TITLE         : LIBRARY WITH CUSTOM FUNCTIONS TO USE IN CALCULATIONS
	! PROJECT       : WF1DCOMDLL
	! MODULE        : mod_sat_ty_nodes
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Class for the collection of nodes-classes in the saturated model
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
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