	!********************************************************************************************************************
	! TITLE         : MOD_COM_QUADRATURE: FUNCTIONS TO CALCULATE THE GAUSS LEGENDRE QUADRATURE OF ORDER 1 TO 60
	! PROJECT       : WF1DCOMDLL
	! MODULE        : MOD_COM_QUADRATURE
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> This module include a function to calculate the gauss legendre quadrature of order no (1 to 60) for a function
	!> of the relative coord of an element
	!>
	!>    Quadrature1D_sca(func,NO): calculate the Gauss-Legendre quadrature of order NO (1 to 60) for the function
	!>                               func(chi). Being chi the relative coordinate on the element from -1 to 1
	!>
	!>    Quadrature1D(func,NO):   calculate the Gauss-Legendre quadrature of order NO (1 to 60) for the function
	!>                               func(chi). In this case chi is defined as a vector
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module com_mod_fem_quadrature

	implicit none

	include 'inc_precision.fi'

	!Weights and roots of legendre polinomials:
	include 'inc_legendre_parameters.fi'
	include 'inc_legendre_weights.fi'

	public::quadrature1d,quadrature1d_sca

	contains
	!DEFINE ALL THE OPERATORS FOR A PIECEWISE LINEAR BASIS ***************************************************************

	!---------------------------------------------------------------------------------------------------------------------
	! quadrature1d(func,no)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Calculate the Gauss-Legendre quadrature in 1D element for a function func fo an order 'no' between 1 and 61
	!> @param[in] func  real function (function to integrate) of one real vector argument
	!> @param[in] no    quadrature order of the integral
	!> @return A vector with the values of the integral on the element
	!---------------------------------------------------------------------------------------------------------------------

	real(kind=dpd) function quadrature1d(func,no)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"quadrature1d" :: quadrature1d
	!DEC$ endif

	interface
	function func(ch)
	include 'inc_precision.fi'
	real(kind=dpd),intent(in)::ch(:)
	real(kind=dpd)::func(size(ch))
	end function func
	end interface

	integer,intent(in)::no !<quadrature order (number of points to be integrated)

	select case (no) !select case depending on quadrature order (number of integration points)
	case(1)
		quadrature1d = dot_product(wi01,func(chii01))
	case(2)
		quadrature1d = dot_product(wi02,func(chii02))
	case(3)
		quadrature1d = dot_product(wi03,func(chii03))
	case(4)
		quadrature1d = dot_product(wi04,func(chii04))
	case(5)
		quadrature1d = dot_product(wi05,func(chii05))
	case(6)
		quadrature1d = dot_product(wi06,func(chii06))
	case(7)
		quadrature1d = dot_product(wi07,func(chii07))
	case(8)
		quadrature1d = dot_product(wi08,func(chii08))
	case(9)
		quadrature1d = dot_product(wi09,func(chii09))
	case(10)
		quadrature1d = dot_product(wi10,func(chii10))
	case(11)
		quadrature1d = dot_product(wi11,func(chii11))
	case(12)
		quadrature1d = dot_product(wi12,func(chii12))
	case(13)
		quadrature1d = dot_product(wi13,func(chii13))
	case(14)
		quadrature1d = dot_product(wi14,func(chii14))
	case(15)
		quadrature1d = dot_product(wi15,func(chii15))
	case(16)
		quadrature1d = dot_product(wi16,func(chii16))
	case(17)
		quadrature1d = dot_product(wi17,func(chii17))
	case(18)
		quadrature1d = dot_product(wi18,func(chii18))
	case(19)
		quadrature1d = dot_product(wi19,func(chii19))
	case(20)
		quadrature1d = dot_product(wi20,func(chii20))
	case(21)
		quadrature1d = dot_product(wi21,func(chii21))
	case(22)
		quadrature1d = dot_product(wi22,func(chii22))
	case(23)
		quadrature1d = dot_product(wi23,func(chii23))
	case(24)
		quadrature1d = dot_product(wi24,func(chii24))
	case(25)
		quadrature1d = dot_product(wi25,func(chii25))
	case(26)
		quadrature1d = dot_product(wi26,func(chii26))
	case(27)
		quadrature1d = dot_product(wi27,func(chii27))
	case(28)
		quadrature1d = dot_product(wi28,func(chii28))
	case(29)
		quadrature1d = dot_product(wi29,func(chii29))
	case(30)
		quadrature1d = dot_product(wi30,func(chii30))
	case(31)
		quadrature1d = dot_product(wi31,func(chii31))
	case(32)
		quadrature1d = dot_product(wi32,func(chii32))
	case(33)
		quadrature1d = dot_product(wi33,func(chii33))
	case(34)
		quadrature1d = dot_product(wi34,func(chii34))
	case(35)
		quadrature1d = dot_product(wi35,func(chii35))
	case(36)
		quadrature1d = dot_product(wi36,func(chii36))
	case(37)
		quadrature1d = dot_product(wi37,func(chii37))
	case(38)
		quadrature1d = dot_product(wi38,func(chii38))
	case(39)
		quadrature1d = dot_product(wi39,func(chii39))
	case(40)
		quadrature1d = dot_product(wi40,func(chii40))
	case(41)
		quadrature1d = dot_product(wi41,func(chii41))
	case(42)
		quadrature1d = dot_product(wi42,func(chii42))
	case(43)
		quadrature1d = dot_product(wi43,func(chii43))
	case(44)
		quadrature1d = dot_product(wi44,func(chii44))
	case(45)
		quadrature1d = dot_product(wi45,func(chii45))
	case(46)
		quadrature1d = dot_product(wi46,func(chii46))
	case(47)
		quadrature1d = dot_product(wi47,func(chii47))
	case(48)
		quadrature1d = dot_product(wi48,func(chii48))
	case(49)
		quadrature1d = dot_product(wi49,func(chii49))
	case(50)
		quadrature1d = dot_product(wi50,func(chii50))
	case(51)
		quadrature1d = dot_product(wi51,func(chii51))
	case(52)
		quadrature1d = dot_product(wi52,func(chii52))
	case(53)
		quadrature1d = dot_product(wi53,func(chii53))
	case(54)
		quadrature1d = dot_product(wi54,func(chii54))
	case(55)
		quadrature1d = dot_product(wi55,func(chii55))
	case(56)
		quadrature1d = dot_product(wi56,func(chii56))
	case(57)
		quadrature1d = dot_product(wi57,func(chii57))
	case(58)
		quadrature1d = dot_product(wi58,func(chii58))
	case(59)
		quadrature1d = dot_product(wi59,func(chii59))
	case(60)
		quadrature1d = dot_product(wi60,func(chii60))

		case default
		quadrature1d = dot_product(wi20,func(chii20))
	end select

	end function quadrature1d

	!---------------------------------------------------------------------------------------------------------------------
	! quadrature1d(func,no)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Calculate the Gauss-Legendre quadrature in 1D element for a function func fo an order 'no' between 1 and 61
	!> @param[in] func  real function (function to integrate) of one real scalar argument
	!> @param[in] no    quadrature order of the integral
	!> @return					A scalar with the values of the integral on the element
	!---------------------------------------------------------------------------------------------------------------------

	real(kind=dpd) function quadrature1d_sca(func,no)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"quadrature1d_sca" :: quadrature1d_sca
	!DEC$ endif

	interface
	function func(ch)
	include 'inc_precision.fi'
	real(kind=dpd),intent(in)::ch
	real(kind=dpd)::func
	end function func
	end interface

	integer,intent(in)::no !<quadrature order (number of points to be integrated)

	select case (no) !select case depending on quadrature order (number of integration points)
	case(1)
		quadrature1d_sca = dot_product(wi01,func_mat(chii01))
	case(2)
		quadrature1d_sca = dot_product(wi02,func_mat(chii02))
	case(3)
		quadrature1d_sca = dot_product(wi03,func_mat(chii03))
	case(4)
		quadrature1d_sca = dot_product(wi04,func_mat(chii04))
	case(5)
		quadrature1d_sca = dot_product(wi05,func_mat(chii05))
	case(6)
		quadrature1d_sca = dot_product(wi06,func_mat(chii06))
	case(7)
		quadrature1d_sca = dot_product(wi07,func_mat(chii07))
	case(8)
		quadrature1d_sca = dot_product(wi08,func_mat(chii08))
	case(9)
		quadrature1d_sca = dot_product(wi09,func_mat(chii09))
	case(10)
		quadrature1d_sca = dot_product(wi10,func_mat(chii10))
	case(11)
		quadrature1d_sca = dot_product(wi11,func_mat(chii11))
	case(12)
		quadrature1d_sca = dot_product(wi12,func_mat(chii12))
	case(13)
		quadrature1d_sca = dot_product(wi13,func_mat(chii13))
	case(14)
		quadrature1d_sca = dot_product(wi14,func_mat(chii14))
	case(15)
		quadrature1d_sca = dot_product(wi15,func_mat(chii15))
	case(16)
		quadrature1d_sca = dot_product(wi16,func_mat(chii16))
	case(17)
		quadrature1d_sca = dot_product(wi17,func_mat(chii17))
	case(18)
		quadrature1d_sca = dot_product(wi18,func_mat(chii18))
	case(19)
		quadrature1d_sca = dot_product(wi19,func_mat(chii19))
	case(20)
		quadrature1d_sca = dot_product(wi20,func_mat(chii20))
	case(21)
		quadrature1d_sca = dot_product(wi21,func_mat(chii21))
	case(22)
		quadrature1d_sca = dot_product(wi22,func_mat(chii22))
	case(23)
		quadrature1d_sca = dot_product(wi23,func_mat(chii23))
	case(24)
		quadrature1d_sca = dot_product(wi24,func_mat(chii24))
	case(25)
		quadrature1d_sca = dot_product(wi25,func_mat(chii25))
	case(26)
		quadrature1d_sca = dot_product(wi26,func_mat(chii26))
	case(27)
		quadrature1d_sca = dot_product(wi27,func_mat(chii27))
	case(28)
		quadrature1d_sca = dot_product(wi28,func_mat(chii28))
	case(29)
		quadrature1d_sca = dot_product(wi29,func_mat(chii29))
	case(30)
		quadrature1d_sca = dot_product(wi30,func_mat(chii30))
	case(31)
		quadrature1d_sca = dot_product(wi31,func_mat(chii31))
	case(32)
		quadrature1d_sca = dot_product(wi32,func_mat(chii32))
	case(33)
		quadrature1d_sca = dot_product(wi33,func_mat(chii33))
	case(34)
		quadrature1d_sca = dot_product(wi34,func_mat(chii34))
	case(35)
		quadrature1d_sca = dot_product(wi35,func_mat(chii35))
	case(36)
		quadrature1d_sca = dot_product(wi36,func_mat(chii36))
	case(37)
		quadrature1d_sca = dot_product(wi37,func_mat(chii37))
	case(38)
		quadrature1d_sca = dot_product(wi38,func_mat(chii38))
	case(39)
		quadrature1d_sca = dot_product(wi39,func_mat(chii39))
	case(40)
		quadrature1d_sca = dot_product(wi40,func_mat(chii40))
	case(41)
		quadrature1d_sca = dot_product(wi41,func_mat(chii41))
	case(42)
		quadrature1d_sca = dot_product(wi42,func_mat(chii42))
	case(43)
		quadrature1d_sca = dot_product(wi43,func_mat(chii43))
	case(44)
		quadrature1d_sca = dot_product(wi44,func_mat(chii44))
	case(45)
		quadrature1d_sca = dot_product(wi45,func_mat(chii45))
	case(46)
		quadrature1d_sca = dot_product(wi46,func_mat(chii46))
	case(47)
		quadrature1d_sca = dot_product(wi47,func_mat(chii47))
	case(48)
		quadrature1d_sca = dot_product(wi48,func_mat(chii48))
	case(49)
		quadrature1d_sca = dot_product(wi49,func_mat(chii49))
	case(50)
		quadrature1d_sca = dot_product(wi50,func_mat(chii50))
	case(51)
		quadrature1d_sca = dot_product(wi51,func_mat(chii51))
	case(52)
		quadrature1d_sca = dot_product(wi52,func_mat(chii52))
	case(53)
		quadrature1d_sca = dot_product(wi53,func_mat(chii53))
	case(54)
		quadrature1d_sca = dot_product(wi54,func_mat(chii54))
	case(55)
		quadrature1d_sca = dot_product(wi55,func_mat(chii55))
	case(56)
		quadrature1d_sca = dot_product(wi56,func_mat(chii56))
	case(57)
		quadrature1d_sca = dot_product(wi57,func_mat(chii57))
	case(58)
		quadrature1d_sca = dot_product(wi58,func_mat(chii58))
	case(59)
		quadrature1d_sca = dot_product(wi59,func_mat(chii59))
	case(60)
		quadrature1d_sca = dot_product(wi60,func_mat(chii60))

		case default
		quadrature1d_sca = dot_product(wi20,func_mat(chii20))
	end select

	contains

	function func_mat(chi)
	real(kind=dpd),intent(in)::chi(:)
	real(kind=dpd)::func_mat(size(chi))
	integer::i

	do i=1,size(chi)
		func_mat(i) = func(chi(i))
	end do


	end function func_mat

	end function quadrature1d_sca

	end module com_mod_fem_quadrature