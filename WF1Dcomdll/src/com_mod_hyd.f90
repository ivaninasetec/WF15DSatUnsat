	!******************************************************************************************************************
	!*                                                                                                                *
	!*                                 MODULE: MOD_UNSAT_HYDRAULICFUN                                                 *
	!*                                                                                                                *
	!******************************************************************************************************************
	!* THIS MODULE AGGREGATE ALL VARIABLY SATURATED HYDRAULIC FUNCTIONS																								*
	!* The hydraulic function is selected with the parameter material%kindmad (1 default: Van Genuchten, 2: Exp func  *
	!*    f2_hyd_S_h:   Returns the specific saturation from pressure head                                             *
	!*    f2_hyd_th_h:  Returns water content from pressure head                                                       *
	!*    f2_hyd_kr_h:  Returns relative permeability from pressure head                                               *
	!*    f2_hyd_k_h:   Returns permeability from pressure head                                                        *
	!*    f2_hyd_cap_h: Returns water capacity from pressure head                                                      *
	!*                                                                                                                *
	!*    Iván Campos-Guereta Díez                                                                                    *
	!*    MSc Civil Engineering by Polytechnic University of Madrid                                                   *
	!*    PhD Student by University of Nottingham                                                                     *
	!*    eMBA by International Institute San Telmo in Seville                                                        *
	!*    ivan.camposguereta@nottingham.ac.uk                                                                         *
	!*                                                                                                                *
	!*    This software is copyrighted 2019(C)                                                                        *
	!******************************************************************************************************************

	MODULE com_mod_hyd

	IMPLICIT NONE
	PRIVATE
	include 'inc_precision.fi'

	! !********************************************************************************************************************
	! ! POLIMORPHIC INTERFACES FOR THE FOLLOWING HYDRAULIC FUNCIONS
	! !--------------------------------------------------------------------------------------------------------------------
	! ! f2_hyd_S_h		:Returns Equivalent saturation (sca or vec) from pressure head (sca or vec)
	! ! f2_hyd_th_h	:Returns water content (sca or vec)					from pressure head (sca or vec)
	! ! f2_hyd_kr_h	:Returns relative permeability (sca or vec) from pressure head (sca or vec)
	! ! f2_hyd_k_h		:Returns permeability (sca or vec)					from pressure head (sca or vec)
	! ! f2_hyd_cap_h	:Returns water capacity (sca or vec)				from pressure head (sca or vec)
	! !********************************************************************************************************************
	!
	!  interface f2_hyd_S_h
	!  module procedure f2_hyd_S_h_sca
	!  module procedure f2_hyd_S_h_vec
	!module procedure f2_hyd_S_h_vec2
	!  end interface f2_hyd_S_h
	!
	!  interface f2_hyd_th_h
	!  module procedure f2_hyd_th_h_sca
	!  module procedure f2_hyd_th_h_vec
	!module procedure f2_hyd_th_h_vec2
	!  end interface f2_hyd_th_h
	!
	!  interface f2_hyd_kr_h
	!  module procedure f2_hyd_kr_h_sca
	!  module procedure f2_hyd_kr_h_vec
	!module procedure f2_hyd_kr_h_vec2
	!  end interface f2_hyd_kr_h
	!
	!  interface f2_hyd_k_h
	!  module procedure f2_hyd_k_h_sca
	!  module procedure f2_hyd_k_h_vec
	!module procedure f2_hyd_k_h_vec2
	!  end interface f2_hyd_k_h
	!
	!  interface f2_hyd_cap_h
	!  module procedure f2_hyd_cap_h_sca
	!  module procedure f2_hyd_cap_h_vec
	!module procedure f2_hyd_cap_h_vec2
	!end interface f2_hyd_cap_h
	!
	!interface f2_hyd_dkr_h
	!  module procedure f2_hyd_dkr_h_sca
	!  module procedure f2_hyd_dkr_h_vec
	!module procedure f2_hyd_dkr_h_vec2
	!  end interface f2_hyd_dkr_h
	!
	!  interface f2_hyd_dk_h
	!  module procedure f2_hyd_dk_h_sca
	!  module procedure f2_hyd_dk_h_vec
	!module procedure f2_hyd_dk_h_vec2
	!  end interface f2_hyd_dk_h


	PUBLIC::  f2_hyd_S_h_sca,  f2_hyd_S_h_vec,  f2_hyd_S_h_vec2, &
		f2_hyd_th_h_sca, f2_hyd_th_h_vec, f2_hyd_th_h_vec2,  &
		f2_hyd_kr_h_sca, f2_hyd_kr_h_vec, f2_hyd_kr_h_vec2, &
		f2_hyd_k_h_sca,  f2_hyd_k_h_vec,  f2_hyd_k_h_vec2,  &
		f2_hyd_cap_h_sca,f2_hyd_cap_h_vec,f2_hyd_cap_h_vec2, &
		f2_hyd_dkr_h_sca,f2_hyd_dkr_h_vec,f2_hyd_dkr_h_vec2, &
		f2_hyd_dk_h_sca, f2_hyd_dk_h_vec, f2_hyd_dk_h_vec2, &
		f2_hyd_incs_h1_to_h2_sca

	CONTAINS

	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_SCA(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_S_h_sca(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_S_h_sca" :: f2_hyd_S_h_sca
	!DEC$ endif

	!Returns specific saturation from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_s_h_pw_sca
	USE com_mod_hyd_vg, ONLY: f2_hyd_S_h_VG_sca
	USE com_mod_hyd_bc, ONLY: f2_hyd_S_h_bc_sca
	REAL(KIND=dps),INTENT(IN)::h
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f2_hyd_S_h_PW_sca(h,material)
	CASE(3)
		Rout = f2_hyd_S_h_bc_sca(h,material)
		CASE DEFAULT
		Rout = f2_hyd_S_h_VG_sca(h,material)
	END SELECT

	END FUNCTION f2_hyd_S_h_sca

	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_SCA(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_incs_h1_to_h2_sca(h1,h2,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_incs_h1_to_h2_sca" :: f2_hyd_incs_h1_to_h2_sca
	!DEC$ endif

	!Returns specific saturation from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_incs_h1_to_h2_pw_sca
	USE com_mod_hyd_vg, ONLY: f2_hyd_incs_h1_to_h2_vg_sca
	USE com_mod_hyd_bc, ONLY: f2_hyd_incs_h1_to_h2_bc_sca
	REAL(KIND=dps),INTENT(IN)::h1,h2
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f2_hyd_incs_h1_to_h2_pw_sca(h1,h2,material)
	CASE(3)
		Rout = f2_hyd_incs_h1_to_h2_bc_sca(h1,h2,material)
		CASE DEFAULT
		Rout = f2_hyd_incs_h1_to_h2_vg_sca(h1,h2,material)
	END SELECT

	END FUNCTION f2_hyd_incs_h1_to_h2_sca


	!********************************************************************************************************************
	! F: f2_hyd_S_h_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_S_h_vec(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_S_h_vec" :: f2_hyd_S_h_vec
	!DEC$ endif

	!Returns specific saturation from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_S_h_PW_vec
	USE com_mod_hyd_vg, ONLY: f2_hyd_S_h_VG_vec
	USE com_mod_hyd_bc, ONLY: f2_hyd_S_h_bc_vec
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material(:)
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	WHERE (material%KindMat==2)
		Rout = f2_hyd_S_h_PW_vec(h,material)
	ELSE WHERE (material%KindMat==3)
		Rout = f2_hyd_S_h_bc_vec(h,material)
		ELSE WHERE
		Rout = f2_hyd_S_h_VG_vec(h,material)
	END WHERE

	END FUNCTION f2_hyd_S_h_vec

	!------

	PURE FUNCTION f2_hyd_S_h_vec2(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_S_h_vec2" :: f2_hyd_S_h_vec2
	!DEC$ endif

	!Returns specific saturation from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_S_h_PW_vec2
	USE com_mod_hyd_vg, ONLY: f2_hyd_S_h_VG_vec2
	USE com_mod_hyd_bc, ONLY: f2_hyd_S_h_bc_vec2
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	if (material%KindMat==2) then
		Rout = f2_hyd_S_h_PW_vec2(h,material)
	ELSE if (material%KindMat==3) then
		Rout = f2_hyd_S_h_BC_vec2(h,material)
		ELSE
		Rout = f2_hyd_S_h_VG_vec2(h,material)
	END if

	END FUNCTION f2_hyd_S_h_vec2


	!********************************************************************************************************************
	! F: f2_hyd_th_h_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_th_h_sca(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_sca" :: f2_hyd_th_h_sca
	!DEC$ endif

	!Returns water content from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_th_h_PW_sca
	USE com_mod_hyd_vg, ONLY: f2_hyd_th_h_VG_sca
	USE com_mod_hyd_bc, ONLY: f2_hyd_th_h_bc_sca
	REAL(KIND=dps),INTENT(IN)::h
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f2_hyd_th_h_PW_sca(h,material)
	CASE(3)
		Rout = f2_hyd_th_h_BC_sca(h,material)
		CASE DEFAULT
		Rout = f2_hyd_th_h_VG_sca(h,material)
	END SELECT

	END FUNCTION f2_hyd_th_h_sca


	!********************************************************************************************************************
	! F: f2_hyd_th_h_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_th_h_vec(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_vec" :: f2_hyd_th_h_vec
	!DEC$ endif

	!Returns water content from pore pressure using kindmat function

	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_th_h_PW_vec
	USE com_mod_hyd_vg, ONLY: f2_hyd_th_h_VG_vec
	USE com_mod_hyd_bc, ONLY: f2_hyd_th_h_bc_vec

	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material(:)
	REAL(KIND=dps)::Rout(SIZE(h))
	integer::i

	WHERE (material%KindMat==2)
		Rout = f2_hyd_th_h_PW_vec(h,material)
	ELSE WHERE (material%KindMat==3)
		Rout = f2_hyd_th_h_BC_vec(h,material)
	ELSE WHERE
		Rout = f2_hyd_th_h_VG_vec(h,material)
	END WHERE
	!do i=1,size(h)
	!	Rout(i)=material(i)%get_th_sca(h(i))
	!end do
	

	END FUNCTION f2_hyd_th_h_vec

	!--------------

	PURE FUNCTION f2_hyd_th_h_vec2(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_vec2" :: f2_hyd_th_h_vec2
	!DEC$ endif

	!Returns water content from pore pressure using kindmat function

	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_th_h_PW_vec2
	USE com_mod_hyd_vg, ONLY: f2_hyd_th_h_VG_vec2
	USE com_mod_hyd_bc, ONLY: f2_hyd_th_h_bc_vec2

	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	if (material%KindMat==2) then
		Rout = f2_hyd_th_h_PW_vec2(h,material)
	ELSE IF (material%KindMat==3) then
		Rout = f2_hyd_th_h_bc_vec2(h,material)
	ELSE
		Rout = f2_hyd_th_h_VG_vec2(h,material)
	END if

	END FUNCTION f2_hyd_th_h_vec2


	!********************************************************************************************************************
	! F: f2_hyd_kr_h_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_kr_h_sca(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_sca" :: f2_hyd_kr_h_sca
	!DEC$ endif

	!Returns relative permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_kr_h_PW_sca
	USE com_mod_hyd_vg, ONLY: f2_hyd_kr_h_VG_sca
	USE com_mod_hyd_bc, ONLY: f2_hyd_kr_h_bc_sca
	
	REAL(KIND=dps),INTENT(IN)::h
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f2_hyd_kr_h_PW_sca(h,material)
	CASE(3)
		Rout = f2_hyd_kr_h_BC_sca(h,material)
		CASE DEFAULT
		Rout = f2_hyd_kr_h_VG_sca(h,material)
	END SELECT

	END FUNCTION f2_hyd_kr_h_sca


	!********************************************************************************************************************
	! F: f2_hyd_kr_h_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_kr_h_vec(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_vec" :: f2_hyd_kr_h_vec
	!DEC$ endif

	!Returns relative permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_kr_h_PW_vec
	USE com_mod_hyd_vg, ONLY: f2_hyd_kr_h_VG_vec
	USE com_mod_hyd_bc, ONLY: f2_hyd_kr_h_bc_vec
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material(:)
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	WHERE (material%KindMat==2)
		Rout = f2_hyd_kr_h_PW_vec(h,material)
	ELSE WHERE (material%KindMat==3)
		Rout = f2_hyd_kr_h_BC_vec(h,material)
		ELSE WHERE
		Rout = f2_hyd_kr_h_VG_vec(h,material)
	END WHERE

	END FUNCTION f2_hyd_kr_h_vec

	!----------------

	PURE FUNCTION f2_hyd_kr_h_vec2(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_vec2" :: f2_hyd_kr_h_vec2
	!DEC$ endif

	!Returns relative permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_kr_h_PW_vec2
	USE com_mod_hyd_vg, ONLY: f2_hyd_kr_h_VG_vec2
	USE com_mod_hyd_bc, ONLY: f2_hyd_kr_h_bc_vec2
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	if (material%KindMat==2) then
		Rout = f2_hyd_kr_h_PW_vec2(h,material)
	ELSE IF(material%KindMat==3) then
		Rout = f2_hyd_kr_h_BC_vec2(h,material)
		ELSE
		Rout = f2_hyd_kr_h_VG_vec2(h,material)
	END if

	END FUNCTION f2_hyd_kr_h_vec2


	!********************************************************************************************************************
	! F: f2_hyd_k_h_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the permeability (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_k_h_sca(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_sca" :: f2_hyd_k_h_sca
	!DEC$ endif

	!Returns permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_k_h_PW_sca
	USE com_mod_hyd_vg, ONLY: f2_hyd_k_h_VG_sca
	USE com_mod_hyd_bc, ONLY: f2_hyd_k_h_bc_sca
	
	REAL(KIND=dps),INTENT(IN)::h
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f2_hyd_k_h_PW_sca(h,material)
	CASE(3)
		Rout = f2_hyd_k_h_BC_sca(h,material)		
		CASE DEFAULT
		Rout = f2_hyd_k_h_VG_sca(h,material)
	END SELECT

	END FUNCTION f2_hyd_k_h_sca


	!********************************************************************************************************************
	! F: f2_hyd_k_h_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns sthe permeability (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_k_h_vec(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_vec" :: f2_hyd_k_h_vec
	!DEC$ endif

	!Returns permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_k_h_PW_vec
	USE com_mod_hyd_vg, ONLY: f2_hyd_k_h_VG_vec
	USE com_mod_hyd_bc, ONLY: f2_hyd_k_h_bc_vec
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material(:)
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	WHERE (material%KindMat==2)
		Rout = f2_hyd_k_h_PW_vec(h,material)
	ELSE WHERE (material%KindMat==3)
		Rout = f2_hyd_k_h_BC_vec(h,material)		
	ELSE WHERE
		Rout = f2_hyd_k_h_VG_vec(h,material)
	END WHERE

	END FUNCTION f2_hyd_k_h_vec

	!-------------

	PURE FUNCTION f2_hyd_k_h_vec2(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_vec2" :: f2_hyd_k_h_vec2
	!DEC$ endif

	!Returns permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_k_h_PW_vec2
	USE com_mod_hyd_vg, ONLY: f2_hyd_k_h_VG_vec2
	USE com_mod_hyd_bc, ONLY: f2_hyd_k_h_bc_vec2
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	if (material%KindMat==2) then
		Rout = f2_hyd_k_h_PW_vec2(h,material)
	ELSE IF (material%KindMat==3) then
		Rout = f2_hyd_k_h_BC_vec2(h,material)
	ELSE
		Rout = f2_hyd_k_h_VG_vec2(h,material)
	END if


	END FUNCTION f2_hyd_k_h_vec2


	!********************************************************************************************************************
	! F: f2_hyd_Cap_h_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the water capacity (scalar) from pressure head (scalar pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_Cap_h_sca(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_Cap_h_sca" :: f2_hyd_Cap_h_sca
	!DEC$ endif

	!Returns water capacity from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_cap_h_PW_sca
	USE com_mod_hyd_vg, ONLY: f2_hyd_cap_h_VG_sca
	USE com_mod_hyd_bc, ONLY: f2_hyd_cap_h_bc_sca
	
	REAL(KIND=dps),INTENT(IN)::h
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f2_hyd_Cap_h_PW_sca(h,material)
	CASE(3)
		Rout = f2_hyd_Cap_h_BC_sca(h,material)
		CASE DEFAULT
		Rout = f2_hyd_cap_h_VG_sca(h,material)
	END SELECT

	END FUNCTION f2_hyd_Cap_h_sca


	!********************************************************************************************************************
	! F: f2_hyd_Cap_h_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the water capacity (vector) from pressure head (vector pressure) using the hyraulicfuncion
	! defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_Cap_h_vec(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_Cap_h_vec" :: f2_hyd_Cap_h_vec
	!DEC$ endif

	!Returns water capacity from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_cap_h_PW_vec
	USE com_mod_hyd_vg, ONLY: f2_hyd_cap_h_VG_vec
	USE com_mod_hyd_bc, ONLY: f2_hyd_cap_h_bc_vec
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material(:)
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	WHERE (material%KindMat==2)
		Rout = f2_hyd_Cap_h_PW_vec(h,material)
	ELSE WHERE (material%KindMat==3)
		Rout = f2_hyd_Cap_h_BC_vec(h,material)		
	ELSE WHERE
		Rout = f2_hyd_cap_h_VG_vec(h,material)
	END WHERE

	END FUNCTION f2_hyd_Cap_h_vec

	!-------

	PURE FUNCTION f2_hyd_Cap_h_vec2(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_Cap_h_vec2" :: f2_hyd_Cap_h_vec2
	!DEC$ endif

	!Returns water capacity from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_cap_h_PW_vec2
	USE com_mod_hyd_vg, ONLY: f2_hyd_cap_h_VG_vec2
	USE com_mod_hyd_bc, ONLY: f2_hyd_cap_h_bc_vec2
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	if (material%KindMat==2) then
		Rout = f2_hyd_Cap_h_PW_vec2(h,material)
	ELSE IF (material%KindMat==3) then
		Rout = f2_hyd_cap_h_BC_vec2(h,material)
	ELSE
		Rout = f2_hyd_cap_h_VG_vec2(h,material)
	END if

	END FUNCTION f2_hyd_Cap_h_vec2


	!********************************************************************************************************************
	! F: f2_hyd_dkr_h_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the derivative of relative permeability (scalar) respect pressure head (scalar pressure)
	! using the hyraulicfuncion defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_dkr_h_sca(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_sca" :: f2_hyd_dkr_h_sca
	!DEC$ endif

	!Returns relative permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_dkr_h_PW_sca
	USE com_mod_hyd_vg, ONLY: f2_hyd_dkr_h_VG_sca
	USE com_mod_hyd_bc, ONLY: f2_hyd_dkr_h_bc_sca
	
	REAL(KIND=dps),INTENT(IN)::h
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f2_hyd_dkr_h_PW_sca(h,material)
	CASE(3)
		Rout = f2_hyd_dkr_h_BC_sca(h,material)
		CASE DEFAULT
		Rout = f2_hyd_dkr_h_VG_sca(h,material)
	END SELECT

	END FUNCTION f2_hyd_dkr_h_sca


	!********************************************************************************************************************
	! F: f2_hyd_dkr_h_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the derivative of relative permeability (vector) respect pressure head (vector pressure)
	! using the hyraulicfuncion defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_dkr_h_vec(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_vec" :: f2_hyd_dkr_h_vec
	!DEC$ endif

	!Returns relative permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_dkr_h_PW_vec
	USE com_mod_hyd_vg, ONLY: f2_hyd_dkr_h_VG_vec
	USE com_mod_hyd_bc, ONLY: f2_hyd_dkr_h_bc_vec
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material(:)
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	WHERE (material%KindMat==2)
		Rout = f2_hyd_dkr_h_PW_vec(h,material)
	ELSE WHERE (material%KindMat==3)
		Rout = f2_hyd_dkr_h_BC_vec(h,material)
	ELSE WHERE
		Rout = f2_hyd_dkr_h_VG_vec(h,material)
	END WHERE

	END FUNCTION f2_hyd_dkr_h_vec

	!-------

	PURE FUNCTION f2_hyd_dkr_h_vec2(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_vec2" :: f2_hyd_dkr_h_vec2
	!DEC$ endif

	!Returns relative permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_dkr_h_PW_vec2
	USE com_mod_hyd_vg, ONLY: f2_hyd_dkr_h_VG_vec2
	USE com_mod_hyd_bc, ONLY: f2_hyd_dkr_h_bc_vec2
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	if (material%KindMat==2) then
		Rout = f2_hyd_dkr_h_PW_vec2(h,material)
	ELSE if (material%KindMat==3) then
		Rout = f2_hyd_dkr_h_BC_vec2(h,material)
	ELSE
		Rout = f2_hyd_dkr_h_VG_vec2(h,material)
	END if

	END FUNCTION f2_hyd_dkr_h_vec2


	!********************************************************************************************************************
	! F: f2_hyd_dk_h_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the derivative of permeability (scalar) respect pressure head (scalar pressure)
	! using the hyraulicfuncion defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_dk_h_sca(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_sca" :: f2_hyd_dk_h_sca
	!DEC$ endif

	!Returns permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_dk_h_PW_sca
	USE com_mod_hyd_vg, ONLY: f2_hyd_dk_h_VG_sca
	USE com_mod_hyd_bc, ONLY: f2_hyd_dk_h_bc_sca
	
	REAL(KIND=dps),INTENT(IN)::h
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout

	!Select proper hydraulic function depending on KindMat
	SELECT CASE (material%KindMat)
	CASE(2)
		Rout = f2_hyd_dk_h_PW_sca(h,material)
	CASE(3)
		Rout = f2_hyd_dk_h_BC_sca(h,material)
		CASE DEFAULT
		Rout = f2_hyd_dk_h_VG_sca(h,material)
	END SELECT

	END FUNCTION f2_hyd_dk_h_sca


	!********************************************************************************************************************
	! F: f2_hyd_dk_h_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns the derivative of permeability (vector) respect pressure head (vector pressure)
	! using the hyraulicfuncion defined with parameter MATERIAL%KINDMAT:
	!	1: Mualem-Van Genuchten (Default)
	!	2: Exp function
	!********************************************************************************************************************

	PURE FUNCTION f2_hyd_dk_h_vec(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_vec" :: f2_hyd_dk_h_vec
	!DEC$ endif

	!Returns permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_dk_h_PW_vec
	USE com_mod_hyd_vg, ONLY: f2_hyd_dk_h_VG_vec
	USE com_mod_hyd_bc, ONLY: f2_hyd_dk_h_bc_vec
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material(:)
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	WHERE (material%KindMat==2)
		Rout = f2_hyd_dk_h_PW_vec(h,material)
	ELSE WHERE (material%KindMat==3)
		Rout = f2_hyd_dk_h_BC_vec(h,material)
	ELSE WHERE
		Rout = f2_hyd_dk_h_VG_vec(h,material)
	END WHERE


	END FUNCTION f2_hyd_dk_h_vec

	!------

	PURE FUNCTION f2_hyd_dk_h_vec2(h,material) RESULT(Rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_vec2" :: f2_hyd_dk_h_vec2
	!DEC$ endif

	!Returns permeability from pore pressure using kindmat function
	USE com_mod_ty_material,ONLY:TY_COM_MATERIAL
	USE com_mod_hyd_pw, ONLY: f2_hyd_dk_h_PW_vec2
	USE com_mod_hyd_vg, ONLY: f2_hyd_dk_h_VG_vec2
	USE com_mod_hyd_bc, ONLY: f2_hyd_dk_h_bc_vec2
	
	REAL(KIND=dps),INTENT(IN)::h(:)
	TYPE(TY_COM_MATERIAL),INTENT(IN)::material
	REAL(KIND=dps)::Rout(SIZE(h))

	!Select proper hydraulic function depending on KindMat
	if (material%KindMat==2) then
		Rout = f2_hyd_dk_h_PW_vec2(h,material)
	ELSE if (material%KindMat==3) then
		Rout = f2_hyd_dk_h_BC_vec2(h,material)
	ELSE
		Rout = f2_hyd_dk_h_VG_vec2(h,material)
	END if


	END FUNCTION f2_hyd_dk_h_vec2

	END MODULE com_mod_hyd