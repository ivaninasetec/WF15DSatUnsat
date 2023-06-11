	!********************************************************************************************************************
	! TITLE         : LIBRARY OF FUNCTIONS TO PERFORM LINEAR INTERPOLATIONS
	! PROJECT       : WF1DCOMDLL
	! MODULE        : COM_MOD_LINEAR_INTERPOLATION
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Module with a function to perform a linear interpolation
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************


	module com_mod_linear_interpolation

	implicit none
	private
	include 'inc_precision.fi'

	public:: linearinterpolation_ordered

	contains

	!********************************************************************************************************************
	! F: linearinterpolation_ordered(X,XVAL(:),YVAL(:))
	!--------------------------------------------------------------------------------------------------------------------
	!		Linear interpolation at x given values of points xval(:) and yval(:)
	!   x:					real value where to get the interpolated value
	!   XVAL(:)  Coordinate x of known points
	!   YVAL(:)  Known function values at x
	!********************************************************************************************************************

	function linearinterpolation_ordered(x,xval,yval)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"linearinterpolation_ordered" :: linearinterpolation_ordered
	!DEC$ endif

	!x is sorted from lower to mac

	real(kind=dpd),intent(in)::x
	real(kind=dpd),intent(in)::xval(:)
	real(kind=dpd),intent(in)::yval(:)
	real(kind=dpd)::linearinterpolation_ordered
	integer::idx0(1),idx1(1)
	real(kind=dpd)::x0,x1
	real(kind=dpd)::y0,y1

	logical::isascending

	if(xval(1)<xval(size(xval))) then
		isascending = .true.
	else
		isascending = .false.
	end if

	select case(isascending)
	case(.true.)
		if(x<=xval(1)) then
			linearinterpolation_ordered=yval(1)
		else if (x>=xval(size(xval))) then
			linearinterpolation_ordered=yval(size(xval))
		else
			idx0 = minloc(x-xval,(x-xval)>=0.0_dpd)
			idx1 = idx0+1
			x0 = xval(idx0(1))
			x1 = xval(idx1(1))
			y0 = yval(idx0(1))
			y1 = yval(idx1(1))
			linearinterpolation_ordered = y0+(y1-y0)*(x-x0)/(x1-x0)
		end if

	case(.false.)
		if(x>=xval(1)) then
			linearinterpolation_ordered=yval(1)
		else if (x<=xval(size(xval))) then
			linearinterpolation_ordered=yval(size(xval))
		else
			idx1 = minloc(xval-x,(xval-x)>=0.0_dpd)
			idx0 = idx0-1
			x0 = xval(idx0(1))
			x1 = xval(idx1(1))
			y0 = yval(idx0(1))
			y1 = yval(idx1(1))
			linearinterpolation_ordered = y0+(y1-y0)*(x-x0)/(x1-x0)
		end if

	end select

	end function linearinterpolation_ordered


	end module com_mod_linear_interpolation