	!********************************************************************************************************************
	! COM_MOD_FEM_BASIS
	!********************************************************************************************************************
	! TITLE         : LIBRARY OF COMMON FUNCTIONS AND CLASSES TO USE IN FLOW1DSAT, FLOW1DUNSAT AND FLOW15DSATUNSAT.
	! PROJECT       : FLOW15DSATUNSAT
	! MODULE        : MOD_COM_SHAPEFUNCTIONS
	! URL           : ...
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2020
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2019(C)
	!
	!> <B>RETURN THE VALUE AND DERIVATIVE OF THE BASIS FUNCTION GIVEN ABS COORD. ALSO INTERPOLATE WITH BASIS</B>
	!>
	!> This module include 3 functions that returns the value of the piecewise  basis function in an absolute coord
	!> its derivative and also do the linear interpolation given the values on nodes using those basis functions.
	!>
	!> Include the following 3 functions:
	!>
	!>- <UL><B>\link basis_i \endlink (x,i,xnod):</B>	  Returns the value of the  \b i basis (piecewise linear function),
	!>	in the coord \b x given the absolute coords on the nodes \b xnod.</UL>
	!>
	!>- <UL><B>\link dbasis_i \endlink (x,i,xnod):</B>  Returns the derivative of the \b i basis (piecewise linear
	!>	function), in the coord \b x given the absolute coords on the nodes \b xnod.</UL>
	!>
	!>- <UL><B>\link interpolate_with_basis \endlink (x,nodevalues,xpoints):</B> Interpolate at \b x using the basis
	!>  functions with values at nodes \b nodevalues, and with the position of the nodes \b xpoints.</UL>
	!>
	!>- <UL><B>\link dinterpolate_with_basis \endlink (x,nodevalues,xpoints):</B> Interpolate the derivative at \b x
	!>  using basis functions with values at nodes \b nodevalues, and with the position of the nodes \b xpoints.</UL>
	!
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module com_mod_fem_basis

	implicit none
	public
	include 'inc_precision.fi'

	contains

	!--------------------------------------------------------------------------------------------------------------------
	! F: BASIS_I(X,I,XNOD)
	!--------------------------------------------------------------------------------------------------------------------
	!> @brief
	!> Returns the value of the '\b i' basis (piecewise linear function), at the coord '\b x' given the
	!> absolute coords on the nodes '\b xnod'.
	!>
	!> For the case of first element (i=1):
	!>
	!> \f$xrel2 = \frac{{{xnod_{i + 1}} - x}}{{{xnod_{i + 1}} - {xnod_i}}}\f$ \n
	!> \f$basis\_i = \begin{cases} xrel2 & 0 \leq xrel2 \leq 1 \\ 0 & true \end{cases}\f$
	!>
	!> For the last element (i=nnod):
	!>
	!> \f$xrel1 = \frac{{x - xnod_{i - 1}}}{{xno{d_i} - xnod_{i - 1}}}\f$ \n
	!> \f$basis\_i = \begin{cases} xrel1 & 0 \leq xrel1 \leq 1 \\ 0 & true \end{cases}\f$
	!>
	!> For internal elements (1<i<nnod):
	!>
	!> \f$\begin{array}{*{20}{c}} {xrel1 = \frac{{x - xnod_{i - 1}}}{{xno{d_i} - xnod_{i - 1}}}}		&{xrel2 = \frac{{{xnod_{i + 1}} - x}}{{{xnod_{i + 1}} - {xnod_i}}}}  \end{array}\f$ \n
	!> \f$basis\_i = \max \left({xrel1,xrel2}\right)\f$
	!>
	!> @param[in] x(:)				Coord where to get value of basis
	!> @param[in] i						Index of basis function
	!> @param[in] xnod(nnod)	Coord of nodes
	!> @return The values of basis function 'i' at 'x(:)'.
	!> @author Iván Campos-Guereta Díez

	pure function basis_i(x,i,xnod)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"basis_i" :: basis_i
	!DEC$ endif

	real(kind=dps),intent(in)::x(:)
	integer,intent(in)::i
	real(kind=dps),intent(in)::xnod(:)

	real(kind=dpd)::basis_i(size(x))
	real(kind=dpd)::xrel1(size(x))
	real(kind=dpd)::xrel2(size(x))

	integer, parameter:: ndim=1 !At this moment number of dimensions=1
	integer, parameter:: nne=2	!At this moment number of nodes in element=2
	integer::nnod

	nnod = size(xnod) !number of nodes

	!for piecewise linear (2 nodes per element, class.0):
	if (i==1) then !for the case of first element:
		xrel2 = (xnod(i+1)-x)/(xnod(i+1)-xnod(i))
		basis_i = merge(xrel2,0.0_dpd,(xrel2>=0).and.(xrel2<=1))
		!basis_i = (xrel2>=0)*(xrel2<=1)*xrel2
		!basis_i = min(1.0_dpd,max(0.0_dpd, xrel2))
	else if (i==nnod) then !for the case of last element:
		xrel1 = (x-xnod(i-1))/(xnod(i)-xnod(i-1))
		basis_i = merge(xrel1,0.0_dpd,(xrel1>=0).and.(xrel1<=1))
		!basis_i = (xrel1>=0)*(xrel1<=1)*xrel1
		!basis_i = min(1.0_dpd,max(0.0_dpd, xrel1))
	else !interior element
		xrel2 = (xnod(i+1)-x)/(xnod(i+1)-xnod(i))
		!xrel2 = (xrel2>=0)*(xrel2<=1)*xrel2
		xrel1 = (x-xnod(i-1))/(xnod(i)-xnod(i-1))
		!xrel1 = (xrel1>=0)*(xrel1<=1)*xrel1
		!basis_i = max(0.0_dpd,max(xrel1,xrel2))
		basis_i = max(merge(xrel1,0.0_dpd,(xrel1>=0).and.(xrel1<=1)),merge(xrel2,0.0_dpd,(xrel2>=0).and.(xrel2<=1)))
	end if

	end function basis_i

	!--------------------------------------------------------------------------------------------------------------------
	! F: DBASIS_I(X,I,XNOD)
	!--------------------------------------------------------------------------------------------------------------------
	!> @brief
	!> Returns the first derivative of the '\b i' basis (piecewise linear function), at the coord '\b x' given the
	!> absolute coords on the nodes '\b xnod'.
	!>
	!> For the case of first element (i=1):
	!>
	!> \f$xrel2 = \frac{{{xnod_{i + 1}} - x}}{{{xnod_{i + 1}} - {xnod_i}}}\f$ \n
	!> \f$dbasis\_i = \begin{cases} \frac{{-1}}{{{xnod_{i + 1}} - {xnod_i}}} & 0 \leq xrel2 < 1 \\ 0 & true \end{cases}\f$
	!>
	!> For the last element (i=nnod):
	!>
	!> \f$xrel1 = \frac{{x - xnod_{i - 1}}}{{xnod_i - xnod_{i - 1}}}\f$ \n
	!> \f$dbasis\_i = \begin{cases} \frac{{1}}{{xnod_i - xnod_{i - 1}}} & 0 < xrel1 \leq 1 \\ 0 & true \end{cases}\f$
	!>
	!> For internal elements (1<i<nnod):
	!>
	!> 	\f$\begin{array}{*{20}{c}} {dbasis\_{i_{left}}= \begin{cases} \frac{{1}}{{xnod_i - xnod_{i - 1}}} & 0 \leq xrel1 < 1 \\ 0 & true \end{cases}}
	!>                            &{dbasis\_i_{right} = \begin{cases} \frac{{-1}}{{{xnod_{i + 1}} - {xnod_i}}} & 0 < xrel2 \leq 1 \\ 0 & true \end{cases}}
	!>  \end{array}\f$ \n
	!> \f$dbasis\_i = \max \left({dbasis\_i_{left}  , dbasis\_i_{right} }\right)\f$
	!>
	!> @param[in] x(:)				Coord where to get value of basis
	!> @param[in] i						Index of basis function
	!> @param[in] xnod(nnod)	Coord of nodes
	!> @return The values of basis function 'i' at 'x(:)'.
	!> @author Iván Campos-Guereta Díez
	!> @param[in] x(:)		Coords where to get derivative of basis
	!> @param[in] i				Index of basis function
	!> @param[in] xnod(nnod)	Coord of nodes
	!> @return The derivatives of basis function 'i' at 'x(:)'.
	!> @author Iván Campos-Guereta Díez

	pure function dbasis_i(x,i,xnod)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"dbasis_i" :: dbasis_i
	!DEC$ endif

	real(kind=dps),intent(in)::x(:)
	integer,intent(in)::i
	real(kind=dps),intent(in)::xnod(:)

	real(kind=dpd)::dbasis_i(size(x)),xrel1(size(x)),xrel2(size(x)),dbasis_i_left(size(x)),dbasis_i_right(size(x))

	integer, parameter:: ndim=1 !number of dimensions=1
	integer, parameter:: nne=2	!number of nodes in element=2
	integer::nnod

	nnod = size(xnod)

	!case of piecewise linear (2 nodes per element, class 0).
	if (i==1) then !for the case of first element:
		xrel2 = (xnod(i+1)-x)/(xnod(i+1)-xnod(i))
		dbasis_i = merge(-1.0_dpd/(xnod(i+1)-xnod(i)),0.0_dpd,(xrel2>=0).and.(xrel2<1))
	else if (i==nnod) then !for the case of last element:
		xrel1 = (x-xnod(i-1))/(xnod(i)-xnod(i-1))
		dbasis_i = merge(1.0_dpd/(xnod(i)-xnod(i-1)),0.0_dpd,(xrel1>0).and.(xrel1<=1))
	else  !interior element, in this case is undefined and depend on the element where is included
		xrel2 = (xnod(i+1)-x)/(xnod(i+1)-xnod(i))
		!xrel2 = (xrel2>0)*(xrel2<1)*xrel2
		xrel1 = (x-xnod(i-1))/(xnod(i)-xnod(i-1))
		!xrel1 = (xrel1>0)*(xrel1<1)*xrel1
		dbasis_i_left = merge(1.0_dpd/(xnod(i+1)-xnod(i)), 0.0_dpd,(xrel1>=0).and.(xrel1<1))
		dbasis_i_right = merge(-1.0_dpd/(xnod(i+1)-xnod(i)), 0.0_dpd,(xrel2>0).and.(xrel2<=1))

		dbasis_i = merge(dbasis_i_left,dbasis_i_right,(xrel1>=0).and.(xrel1<1));

		!dbasis_i = max(dbasis_i_left,dbasis_i_right)
	end if

	end function dbasis_i

	!--------------------------------------------------------------------------------------------------------------------
	! F: interpolate_with_basis(x,nodevalues,xpoints)
	!--------------------------------------------------------------------------------------------------------------------
	!> @brief
	!> Interpolate at coord \b x given the values on nodes (\b nodevalues) and coords of nodes (\b xpoints), unsing basis
	!> functions.
	!>
	!>  \f$interpolate\_with\_basis = \sum\limits_{i = 1}^{nnod} {nodevalues_i \cdot basis\_i \left( {x,i,xpoints} \right)}\f$
	!>
	!> @param[in] x(:)							Coords where to get the interpolated values.
	!> @param[in] nodevalues(nnod)	Values at the nodes
	!> @param[in] xpoints(nnod)			Coordinates of the nodes.
	!> @return Interpolated values at \b x.
	!> @author Iván Campos-Guereta Díez

	pure function interpolate_with_basis(x,nodevalues,xpoints)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"interpolate_with_basis" :: interpolate_with_basis
	!DEC$ endif




	real(kind=dps),intent(in)::nodevalues(:),xpoints(:)!values and coordinate values of nodes
	real(kind=dps),intent(in)::x(:)			!x value where to get value interpolated with basis
	integer::i							!number of the basis function
	real(kind=dpd)::interpolate_with_basis(size(x))

	interpolate_with_basis = 0.0_dps
	do i=1, size(nodevalues)
		interpolate_with_basis = interpolate_with_basis + nodevalues(i)*basis_i(x,i,xpoints)
	end do

	end function interpolate_with_basis

	!--------------------------------------------------------------------------------------------------------------------
	! F: interpolate_with_basis(x,nodevalues,xpoints)
	!--------------------------------------------------------------------------------------------------------------------
	!> @brief
	!> Interpolate at coord \b x given the values on nodes (\b nodevalues) and coords of nodes (\b xpoints), unsing basis
	!> functions.
	!>
	!>  \f$dinterpolate\_with\_basis = \sum\limits_{i = 1}^{nnod} {nodevalues_i \cdot dbasis\_i \left( {x,i,xpoints} \right)}\f$
	!>
	!> @param[in] x(:)							Coords where to get the interpolated values.
	!> @param[in] nodevalues(nnod)	Values at the nodes
	!> @param[in] xpoints(nnod)			Coordinates of the nodes.
	!> @return Interpolated values at \b x.
	!> @author Iván Campos-Guereta Díez

	pure function dinterpolate_with_basis(x,nodevalues,xpoints)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"dinterpolate_with_basis" :: dinterpolate_with_basis
	!DEC$ endif




	real(kind=dps),intent(in)::x(:)
	real(kind=dps),intent(in)::nodevalues(:)
	real(kind=dps),intent(in)::xpoints(:)
	real(kind=dpd)::dinterpolate_with_basis(size(x))

	integer::i

	dinterpolate_with_basis = 0.0_dps
	do i=1, size(nodevalues)
		dinterpolate_with_basis = dinterpolate_with_basis + nodevalues(i)*dbasis_i(x,i,xpoints)
	end do

	end function dinterpolate_with_basis

	end module com_mod_fem_basis