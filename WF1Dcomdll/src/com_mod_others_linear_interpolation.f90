	!********************************************************************************************************************
	!        INCLUDE FUNCTIONS FOR LINEAR INTERPOLATION
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
	!  MSc Civil Engineering by Polytechnic University of Madrid
	!  PhD Student by University of Nottingham
	!  eMBA by International Institute San Telmo in Seville
	!  ivan.camposguereta@nottingham.ac.uk
	! DESCRIPTION:
	!> Class for the collection of nodes-classes in the saturated model
	!********************************************************************************************************************


	module com_mod_linear_interpolation

	implicit none
	private
	include 'inc_precision.fi'

	public:: linearinterpolation_ordered

	contains

	!********************************************************************************************************************
	! F: LINEARINTERPOLATION(FUNC1D,XVAL(:),YVAL(:))
	!--------------------------------------------------------------------------------------------------------------------
	!		Function that returns the integral in the element. In this case the function is in absolute coords to the element.
	!   func1d(x):         Is the relative coord at wich calculate the jacobian
	!   xne(nnodes): Is the absolute coodinates of the nodes of element.
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