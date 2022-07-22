	!********************************************************************************************************************
	! TITLE         : COM_MOD_TY_MATERIAL: DERIVED TYPE THAT DEFINES COMMON PROPERTIES AND METHODS OF THE MATERIALS
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_ELEMENTS
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Derived type that defines properties and methods of the materials
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module com_mod_ty_material
	implicit none
	include 'inc_precision.fi'

	private

	!******************************************************************************************************************
	! TY_COM_MATERIAL
	! Derived type that defines a material
	!------------------------------------------------------------------------------------------------------------------
	!	 |ID              [integer]:        ID number
	!	 |KindMat         [integer]:        Kind of material 1: Van Genuchten-Mualem, 2:Exp hyd functions
	!	 |ThSat           [real]:           Saturated water content (L3/L3)
	!	 |ThRes           [real]:           Residual water content  (L3/L3)
	!	 |ksat            [real]:           Saturated permeability  (L/T)
	!	 |a               [real]:           Parameter a of Van-Genuchten (L-1) or a for Exp hyd func (L-1)
	!	 |n               [real]:           Parameter n of Van-Genuchten (-) or n for Exp hyd func (-)
	!	 |m               [real]:           Parameter m of Van-Genuchten (-)
	!	 |l               [real]:           Parameter l of Mualem-Van-Genuchten (-)
	!******************************************************************************************************************

	type,public::ty_com_material !< CLASS: Material (common description)
		integer				::id			!< ID number
		integer				::kindmat !< Kind of material (1: Van Genuchten-Mualem, 2:Exp hyd functions)
		real(kind=dps)::thsat		!< Saturated water content [L3/L3]
		real(kind=dps)::thres		!< Residual water content  [L3/L3]
		real(kind=dps)::ksat		!< Saturated permeability  [L/T]
		real(kind=dps)::a				!< Parameter alfa of Van-Genuchten [L-1] or a for Exp hyd func [L-1]
		real(kind=dps)::n				!< Parameter n of Van-Genuchten (-) or n for Exp hyd func (-)
		real(kind=dps)::m				!< Parameter m of Van-Genuchten (-)
		real(kind=dps)::l				!< Parameter l of Mualem-Van-Genuchten (-)

	contains
	procedure,public:: get_incs_h1_to_h2_sca		=> f_hyd_incs_h1_to_h2_sca
	procedure,public:: get_s_sca		=> f_hyd_s_h_sca
	procedure,public:: get_s_vec		=> f_hyd_s_h_vec2
	procedure,public:: get_th_sca		=> f_hyd_th_h_sca
	procedure,public:: get_th_vec		=> f_hyd_th_h_vec2
	procedure,public:: get_kr_sca		=> f_hyd_kr_h_sca
	procedure,public:: get_kr_vec		=> f_hyd_kr_h_vec2
	procedure,public:: get_k_sca		=> f_hyd_k_h_sca
	procedure,public:: get_k_vec		=> f_hyd_k_h_vec2
	procedure,public:: get_dk_sca		=> f_hyd_dk_h_sca
	procedure,public:: get_dk_vec		=> f_hyd_dk_h_vec2
	procedure,public:: get_cap_sca	=> f_hyd_cap_h_sca
	procedure,public:: get_cap_vec	=> f_hyd_cap_h_vec2
	procedure,public:: get_cmax	=> f_hyd_cmax

	end type ty_com_material

	contains

	!Include here distinct hydraulic functions
	include 'com_mod_hyd_vg.fi' !Van-Genuchten
	include 'com_mod_hyd_pw.fi' !Power function in Hayek,2016
	include 'com_mod_hyd_bc.fi'	!Brooks and Corey
	include 'com_mod_hyd_ba.fi'	!Clean ballast: wrc as Brooks and Corey, kr=Sr^(2+l)

	!********************************************************************************************************************
	! F: F_HYD_S_H_PW_SCA(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_s_h_sca(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_s_h_sca" :: f_hyd_s_h_sca
	!DEC$ endif
	!returns specific saturation from pore pressure using kindmat function
	real(kind=dps),intent(in)::h
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_s_h_pw_sca(material,h)
	case(3)
		rout = f_hyd_s_h_bc_sca(material,h)
	case(4)
		rout = f_hyd_s_h_ba_sca(material,h)
		case default
		rout = f_hyd_s_h_vg_sca(material,h)
	end select

	end function f_hyd_s_h_sca

	!********************************************************************************************************************
	! F: f_hyd_S_h_vec(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_s_h_vec2(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_s_h_vec2" :: f_hyd_s_h_vec2
	!DEC$ endif
	!returns specific saturation from pore pressure using kindmat function

	real(kind=dps),intent(in)::h(:)
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_s_h_pw_vec2(material,h)
	case(3)
		rout = f_hyd_s_h_bc_vec2(material,h)
	case(4)
		rout = f_hyd_s_h_ba_vec2(material,h)
		case default
		rout = f_hyd_s_h_vg_vec2(material,h)
	end select

	end function f_hyd_s_h_vec2

	!********************************************************************************************************************
	! F: f_hyd_th_h_sca(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_th_h_sca(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_th_h_sca" :: f_hyd_th_h_sca
	!DEC$ endif
	!returns water content from pore pressure using kindmat function

	real(kind=dps),intent(in)::h
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	select case (material%kindmat)
	case(2)
		rout = f_hyd_th_h_pw_sca(material,h)
	case(3)
		rout = f_hyd_th_h_bc_sca(material,h)
	case(4)
		rout = f_hyd_th_h_ba_sca(material,h)
		case default
		rout = f_hyd_th_h_vg_sca(material,h)
	end select

	END FUNCTION f_hyd_th_h_sca


	!********************************************************************************************************************
	! F: f_hyd_th_h_vec(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f_hyd_th_h_vec2(material,h) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_th_h_vec2" :: f_hyd_th_h_vec2
	!DEC$ endif
	!Returns water content from pore pressure using kindmat function

	REAL(KIND=dps),INTENT(IN)::h(:)
	class(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	select case (material%kindmat)
	case(2)
		Rout = f_hyd_th_h_pw_vec2(material,h)
	case(3)
		Rout = f_hyd_th_h_bc_vec2(material,h)
	case(4)
		Rout = f_hyd_th_h_ba_vec2(material,h)
		case default
		Rout = f_hyd_th_h_vg_vec2(material,h)
	end select

	END FUNCTION f_hyd_th_h_vec2


	!********************************************************************************************************************
	! F: f_hyd_kr_h_sca(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_kr_h_sca(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_kr_h_sca" :: f_hyd_kr_h_sca
	!DEC$ endif
	!returns relative permeability from pore pressure using kindmat function

	real(kind=dps),intent(in)::h
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_kr_h_pw_sca(material,h)
	case(3)
		rout = f_hyd_kr_h_bc_sca(material,h)
	case(4)
		rout = f_hyd_kr_h_ba_sca(material,h)
		case default
		rout = f_hyd_kr_h_vg_sca(material,h)
	end select

	end function f_hyd_kr_h_sca

	!********************************************************************************************************************
	! F: f_hyd_kr_h_vec(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_kr_h_vec2(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_kr_h_vec2" :: f_hyd_kr_h_vec2
	!DEC$ endif
	!returns relative permeability from pore pressure using kindmat function

	real(kind=dps),intent(in)::h(:)
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_kr_h_pw_vec2(material,h)
	case(3)
		rout = f_hyd_kr_h_bc_vec2(material,h)
	case(4)
		rout = f_hyd_kr_h_ba_vec2(material,h)
		case default
		rout = f_hyd_kr_h_vg_vec2(material,h)
	end select

	end function f_hyd_kr_h_vec2

	!********************************************************************************************************************
	! F: f_hyd_k_h_sca(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the permeability (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_k_h_sca(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_k_h_sca" :: f_hyd_k_h_sca
	!DEC$ endif
	!returns permeability from pore pressure using kindmat function

	real(kind=dps),intent(in)::h
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_k_h_pw_sca(material,h)
	case(3)
		rout = f_hyd_k_h_bc_sca(material,h)
	case(4)
		rout = f_hyd_k_h_ba_sca(material,h)
		case default
		rout = f_hyd_k_h_vg_sca(material,h)
	end select

	end function f_hyd_k_h_sca

	!********************************************************************************************************************
	! F: f_hyd_k_h_vec(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the permeability (vector) from pressure head (vector pressure) using the hyraulic funcion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f_hyd_k_h_vec2(material,h) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_k_h_vec2" :: f_hyd_k_h_vec2
	!DEC$ endif
	!Returns permeability from pore pressure using kindmat function

	REAL(KIND=dps),INTENT(IN)::h(:)
	class(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	select case (material%kindmat)
	case(2)
		Rout = f_hyd_k_h_pw_vec2(material,h)
	case(3)
		Rout = f_hyd_k_h_bc_vec2(material,h)
	case(4)
		Rout = f_hyd_k_h_ba_vec2(material,h)

		case default
		Rout = f_hyd_k_h_VG_vec2(material,h)
	END select

	END FUNCTION f_hyd_k_h_vec2

	!********************************************************************************************************************
	! F: f_hyd_Cap_h_sca(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the water capacity (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_cap_h_sca(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_cap_h_sca" :: f_hyd_cap_h_sca
	!DEC$ endif
	!returns water capacity from pore pressure using kindmat function

	real(kind=dps),intent(in)::h
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_cap_h_pw_sca(material,h)
	case(3)
		rout = f_hyd_cap_h_bc_sca(material,h)
	case(4)
		rout = f_hyd_cap_h_ba_sca(material,h)
		case default
		rout = f_hyd_cap_h_vg_sca(material,h)
	end select

	end function f_hyd_cap_h_sca

	!********************************************************************************************************************
	! F: f_hyd_Cap_h_vec(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the water capacity (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f_hyd_Cap_h_vec2(material,h) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_Cap_h_vec2" :: f_hyd_Cap_h_vec2
	!DEC$ endif
	!Returns water capacity from pore pressure using kindmat function

	REAL(KIND=dps),INTENT(IN)::h(:)
	class(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	select case (material%kindmat)
	case(2)
		Rout = f_hyd_Cap_h_PW_vec2(material,h)
	case(3)
		Rout = f_hyd_Cap_h_bc_vec2(material,h)
	case(4)
		Rout = f_hyd_Cap_h_ba_vec2(material,h)
		case default
		Rout = f_hyd_cap_h_VG_vec2(material,h)
	END select

	END FUNCTION f_hyd_Cap_h_vec2

	!********************************************************************************************************************
	! F: f_hyd_dkr_h_sca(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the derivative of relative permeability (scalar) respect pressure head (scalar pressure)
	! using the hyraulicfuncion defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f_hyd_dkr_h_sca(material,h) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_dkr_h_sca" :: f_hyd_dkr_h_sca
	!DEC$ endif
	!Returns relative permeability from pore pressure using kindmat function

	REAL(KIND=dps),INTENT(IN)::h
	class(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f_hyd_dkr_h_PW_sca(material,h)
	CASE(3)
		Rout = f_hyd_dkr_h_bc_sca(material,h)
	CASE(4)
		Rout = f_hyd_dkr_h_ba_sca(material,h)
		CASE DEFAULT
		Rout = f_hyd_dkr_h_VG_sca(material,h)
	END SELECT

	END FUNCTION f_hyd_dkr_h_sca

	!********************************************************************************************************************
	! F: f_hyd_dkr_h_vec(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the derivative of relative permeability (vector) respect pressure head (vector pressure)
	! using the hyraulicfuncion defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_dkr_h_vec2(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_dkr_h_vec2" :: f_hyd_dkr_h_vec2
	!DEC$ endif
	!returns relative permeability from pore pressure using kindmat function

	real(kind=dps),intent(in)::h(:)
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	!select proper hydraulic function depending on kindmat
	SELECT CASE (material%KindMat)
	CASE(2)
		rout = f_hyd_dkr_h_pw_vec2(material,h)
	CASE(3)
		rout = f_hyd_dkr_h_bc_vec2(material,h)
	CASE(4)
		rout = f_hyd_dkr_h_ba_vec2(material,h)
		case default
		rout = f_hyd_dkr_h_vg_vec2(material,h)
	end select

	end function f_hyd_dkr_h_vec2

	!********************************************************************************************************************
	! F: f_hyd_dk_h_sca(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the derivative of permeability (scalar) respect pressure head (scalar pressure)
	! using the hyraulicfuncion defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_dk_h_sca(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_dk_h_sca" :: f_hyd_dk_h_sca
	!DEC$ endif
	!returns permeability from pore pressure using kindmat function

	real(kind=dps),intent(in)::h
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_dk_h_pw_sca(material,h)
	case(3)
		rout = f_hyd_dk_h_bc_sca(material,h)
	case(4)
		rout = f_hyd_dk_h_ba_sca(material,h)
		case default
		rout = f_hyd_dk_h_vg_sca(material,h)
	end select

	end function f_hyd_dk_h_sca

	!********************************************************************************************************************
	! F: f_hyd_dk_h_vec(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the derivative of permeability (vector) respect pressure head (vector pressure)
	! using the hyraulicfuncion defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************
	pure function f_hyd_dk_h_vec2(material,h) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_dk_h_vec2" :: f_hyd_dk_h_vec2
	!DEC$ endif
	!returns permeability from pore pressure using kindmat function

	real(kind=dps),intent(in)::h(:)
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_dk_h_pw_vec2(material,h)
	case(3)
		rout = f_hyd_dk_h_bc_vec2(material,h)
	case(4)
		rout = f_hyd_dk_h_ba_vec2(material,h)
		case default
		rout = f_hyd_dk_h_vg_vec2(material,h)
	end select

	end function f_hyd_dk_h_vec2


	!********************************************************************************************************************
	! F: f_hyd_incs_h1_to_h2_sca(material,h1,h2)
	!--------------------------------------------------------------------------------------------------------------------
	! Return the icrement of water content from h1 to h2
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f_hyd_incs_h1_to_h2_sca(material,h1,h2) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_hyd_incs_h1_to_h2_sca" :: f_hyd_incs_h1_to_h2_sca
	!DEC$ endif

	!Returns specific saturation from pore pressure using kindmat function
	REAL(KIND=dps),INTENT(IN)::h1,h2
	class(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f_hyd_incs_h1_to_h2_pw_sca(material,h1,h2)
	CASE(3)
		Rout = f_hyd_incs_h1_to_h2_bc_sca(material,h1,h2)
	CASE(4)
		Rout = f_hyd_incs_h1_to_h2_ba_sca(material,h1,h2)
		CASE DEFAULT
		Rout = f_hyd_incs_h1_to_h2_vg_sca(material,h1,h2)
	END SELECT

	END FUNCTION f_hyd_incs_h1_to_h2_sca

	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_SCA(material,h)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	pure function f_hyd_Cmax(material) result(rout)
	!dec$ if defined(_dll)
	!dec$ attributes dllexport, alias:"f_hyd_Cmax" :: f_hyd_Cmax
	!dec$ endif

	!returns specific saturation from pore pressure using kindmat function
	class(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	!select proper hydraulic function depending on kindmat
	select case (material%kindmat)
	case(2)
		rout = f_hyd_Cmax_pw(material)
	case(3)
		rout = f_hyd_Cmax_bc(material)
	case(4)
		rout = f_hyd_Cmax_ba(material)
		case default
		rout = f_hyd_Cmax_vg(material)
	end select

	end function f_hyd_Cmax


	end module com_mod_ty_material