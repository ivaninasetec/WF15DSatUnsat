	!********************************************************************************************************************
	!*                                                                                                                  *
	!*                                 MODULE: MOD_FUN_HYDRAULIC_PW                                                     *
	!*                                                                                                                  *
	!********************************************************************************************************************
	!* THIS MODULE INCLUDE VARIABLY SATURATED HYDRAULIC FUNCTIONS: SIMPLIFIED POWER FUNCITONS	(Paper Hayek,2016)        *
	!*    f2_hyd_S_h_PW:   Returns the specific saturation from pressure head                                            *
	!*    f2_hyd_th_h_PW:  Returns water content from pressure head                                                      *
	!*    f2_hyd_kr_h_PW:  Returns relative permeability from pressure head                                              *
	!*    f2_hyd_k_h_PW:   Returns permeability from pressure head                                                       *
	!*    f2_hyd_cap_h_PW: Returns water capacity from pressure head                                                     *
	!*                                                                                                                  *
	!*    Iván Campos-Guereta Díez                                                                                      *
	!*    MSc Civil Engineering by Polytechnic University of Madrid                                                     *
	!*    PhD Student by University of Nottingham                                                                       *
	!*    eMBA by International Institute San Telmo in Seville                                                          *
	!*    ivan.camposguereta@nottingham.ac.uk                                                                           *
	!*                                                                                                                  *
	!*    This software is copyrighted 2019(C)                                                                          *
	!********************************************************************************************************************

	module com_mod_hyd_pw

	implicit none
	private
	include 'inc_precision.fi'

	!********************************************************************************************************************
	! POLIMORPHIC INTERFACES FOR THE FOLLOWING HYDRAULIC FUNCIONS
	!--------------------------------------------------------------------------------------------------------------------
	! f2_hyd_S_h_PW	:Returns Equivalent saturation (sca or vec) from pressure head (sca or vec)
	! f2_hyd_th_h_PW	:Returns water content (sca or vec)					from pressure head (sca or vec)
	! f2_hyd_kr_h_PW	:Returns relative permeability (sca or vec) from pressure head (sca or vec)
	! f2_hyd_k_h_PW	:Returns permeability (sca or vec)					from pressure head (sca or vec)
	! f2_hyd_cap_h_PW:Returns water capacity (sca or vec)				from pressure head (sca or vec)
	!********************************************************************************************************************

	!interface f2_hyd_s_h_pw
	!module procedure f2_hyd_s_h_pw_sca
	!module procedure f2_hyd_s_h_pw_vec
	!module procedure f2_hyd_s_h_pw_vec2
	!end interface f2_hyd_s_h_pw
	!
	!interface f2_hyd_th_h_pw
	!module procedure f2_hyd_th_h_pw_sca
	!module procedure f2_hyd_th_h_pw_vec
	!module procedure f2_hyd_th_h_pw_vec2
	!end interface f2_hyd_th_h_pw
	!
	!interface f2_hyd_kr_h_pw
	!module procedure f2_hyd_kr_h_pw_sca
	!module procedure f2_hyd_kr_h_pw_vec
	!module procedure f2_hyd_kr_h_pw_vec2
	!end interface f2_hyd_kr_h_pw
	!
	!interface f2_hyd_k_h_pw
	!module procedure f2_hyd_k_h_pw_sca
	!module procedure f2_hyd_k_h_pw_vec
	!module procedure f2_hyd_k_h_pw_vec2
	!end interface f2_hyd_k_h_pw
	!
	!interface f2_hyd_cap_h_pw
	!module procedure f2_hyd_cap_h_pw_sca
	!module procedure f2_hyd_cap_h_pw_vec
	!module procedure f2_hyd_cap_h_pw_vec2
	!end interface f2_hyd_cap_h_pw
	!
	!interface f2_hyd_dkr_h_pw
	!module procedure f2_hyd_dkr_h_pw_sca
	!module procedure f2_hyd_dkr_h_pw_vec
	!module procedure f2_hyd_dkr_h_pw_vec2
	!end interface f2_hyd_dkr_h_pw
	!
	!interface f2_hyd_dk_h_pw
	!module procedure f2_hyd_dk_h_pw_sca
	!module procedure f2_hyd_dk_h_pw_vec
	!module procedure f2_hyd_dk_h_pw_vec2
	!end interface f2_hyd_dk_h_pw
	!
	!public::f2_hyd_s_h_pw,   f2_hyd_th_h_pw,  f2_hyd_kr_h_pw,  f2_hyd_k_h_pw,  f2_hyd_cap_h_pw,  f2_hyd_dkr_h_pw,  f2_hyd_dk_h_pw

	PUBLIC::  f2_hyd_S_h_pw_sca,  f2_hyd_S_h_pw_vec,  f2_hyd_S_h_pw_vec2, &
		f2_hyd_th_h_pw_sca, f2_hyd_th_h_pw_vec, f2_hyd_th_h_pw_vec2,  &
		f2_hyd_kr_h_pw_sca, f2_hyd_kr_h_pw_vec, f2_hyd_kr_h_pw_vec2, &
		f2_hyd_k_h_pw_sca,  f2_hyd_k_h_pw_vec,  f2_hyd_k_h_pw_vec2,  &
		f2_hyd_cap_h_pw_sca,f2_hyd_cap_h_pw_vec,f2_hyd_cap_h_pw_vec2, &
		f2_hyd_dkr_h_pw_sca,f2_hyd_dkr_h_pw_vec,f2_hyd_dkr_h_pw_vec2, &
		f2_hyd_dk_h_pw_sca, f2_hyd_dk_h_pw_vec, f2_hyd_dk_h_pw_vec2, &
		f2_hyd_incs_h1_to_h2_pw_sca

	contains

	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_SCA(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		Se=EXP(a*h/n)
	!********************************************************************************************************************

	pure function f2_hyd_s_h_pw_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_pw_sca" :: f2_hyd_s_h_pw_sca
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	if (h<0.0) then
		rout = merge(exp(material%a*h/material%n),1.0_dps,h<0.0_dps)
	else
		rout = 1.0_dps
	end if

	end function f2_hyd_s_h_pw_sca

	!********************************************************************************************************************
	! F: f2_hyd_incs_h1_to_h2_vg_sca(H1,H2,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		Se1=EXP(a*h1/n) Se2=EXP(a*h2/n) -> Out Se2-Se1
	!********************************************************************************************************************

	pure function f2_hyd_incs_h1_to_h2_pw_sca(h1,h2,material)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_incs_h1_to_h2_pw_sca" :: f2_hyd_incs_h1_to_h2_pw_sca
	!DEC$ endif

	!return the specific saturation from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h1,h2
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::f2_hyd_incs_h1_to_h2_pw_sca
	real(kind=16)::a,hb1,hb2,n,m,rout,rout1,rout2

	a = real(material%a,16)
	hb1 = real(h1,16)
	hb2 = real(h2,16)
	n= real(material%n,16)
	m= real(material%m,16)

	if (hb1<0.0_16) then
		rout1 = merge(exp(material%a*h1/material%n),1.0_dps,h1<0.0_dps)
	else
		rout1 =  1.0_16
	end if

	if (hb2<0.0_16) then
		rout2 = merge(exp(material%a*h2/material%n),1.0_dps,h2<0.0_dps)
	else
		rout2 =  1.0_16
	end if

	rout = rout2-rout1

	f2_hyd_incs_h1_to_h2_pw_sca = real(rout,dpd)

	end function f2_hyd_incs_h1_to_h2_pw_sca



	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_VEC(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (vector) from pressure head (vector pressure), with simplified power
	! function:
	!		Se=EXP(a*h/n)
	!********************************************************************************************************************

	pure function f2_hyd_s_h_pw_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_pw_vec" :: f2_hyd_s_h_pw_vec
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = merge(exp(material%a*h/material%n),1.0_dps,h<0.0_dps)
	else where
		rout = 1.0_dps
	end where

	end function f2_hyd_s_h_pw_vec

	!-----------------

	pure function f2_hyd_s_h_pw_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_pw_vec2" :: f2_hyd_s_h_pw_vec2
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = merge(exp(material%a*h/material%n),1.0_dps,h<0.0_dps)
	else where
		rout = 1.0_dps
	end where

	end function f2_hyd_s_h_pw_vec2


	!********************************************************************************************************************
	!	F: f2_hyd_th_h_PW_sca(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		th = thres+(thsat-thres)*Se
	!********************************************************************************************************************

	pure function f2_hyd_th_h_pw_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_pw_sca" :: f2_hyd_th_h_pw_sca
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	if (h<0.0) then
		rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_pw_sca(h,material)
	else
		rout = material%thsat
	end if

	end function f2_hyd_th_h_pw_sca


	!********************************************************************************************************************
	!	F: f2_hyd_th_h_PW_vec(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (vector) from pressure head (vector pressure), with simplified power
	! function:
	!		th = thres+(thsat-thres)*Se
	!********************************************************************************************************************

	pure function f2_hyd_th_h_pw_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_pw_vec" :: f2_hyd_th_h_pw_vec
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_pw_vec(h,material)
	else where
		rout = material%thsat
	end where

	end function f2_hyd_th_h_pw_vec

	!---------

	pure function f2_hyd_th_h_pw_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_pw_vec2" :: f2_hyd_th_h_pw_vec2
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_pw_vec2(h,material)

	!where (h<0.0)
	!	rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_pw_sca(h,material)
	!else where
	!	rout = material%thsat
	!end where

	end function f2_hyd_th_h_pw_vec2

	!********************************************************************************************************************
	!	F: f2_hyd_kr_h_PW_sca(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		Kr = EXP(a·h)
	!********************************************************************************************************************

	pure function f2_hyd_kr_h_pw_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_pw_sca" :: f2_hyd_kr_h_pw_sca
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	if (h<0.0) then
		rout = exp(material%a*h)
	else
		rout = 1.0_dps
	end if

	end function f2_hyd_kr_h_pw_sca


	!********************************************************************************************************************
	!	F: f2_hyd_kr_h_PW_vec(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (vector) from pressure head (vector pressure), with simplified power
	! function:
	!		Kr = EXP(a·h)
	!********************************************************************************************************************

	pure function f2_hyd_kr_h_pw_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_pw_vec" :: f2_hyd_kr_h_pw_vec
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	!kr = exp(a·h)
	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = exp(material%a*h)
	else where
		rout = 1.0_dps
	end where

	end function f2_hyd_kr_h_pw_vec

	!--------

	pure function f2_hyd_kr_h_pw_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_pw_vec2" :: f2_hyd_kr_h_pw_vec2
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	!kr = exp(a·h)
	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = exp(material%a*h)
	else where
		rout = 1.0_dps
	end where

	end function f2_hyd_kr_h_pw_vec2

	!********************************************************************************************************************
	!	F: f2_hyd_k_h_PW_sca(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns permeability (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		k = ksat·kr
	!********************************************************************************************************************

	pure function f2_hyd_k_h_pw_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_pw_sca" :: f2_hyd_k_h_pw_sca
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	rout = material%ksat*f2_hyd_kr_h_pw_sca(h,material)

	end function f2_hyd_k_h_pw_sca


	!********************************************************************************************************************
	!	F: f2_hyd_k_h_PW_vec(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns permeability (vector) from pressure head (vector pressure), with simplified power
	! function:
	!		k = ksat·kr
	!********************************************************************************************************************

	pure function f2_hyd_k_h_pw_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_pw_vec" :: f2_hyd_k_h_pw_vec
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dps)::rout(size(h))

	rout = material%ksat*f2_hyd_kr_h_pw_vec(h,material)

	end function f2_hyd_k_h_pw_vec

	!---------

	pure function f2_hyd_k_h_pw_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_pw_vec2" :: f2_hyd_k_h_pw_vec2
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	rout = material%ksat*f2_hyd_kr_h_pw_vec2(h,material)

	end function f2_hyd_k_h_pw_vec2

	!********************************************************************************************************************
	!	F: f2_hyd_cap_h_PW_sca(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water capacity (dth/dh) (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		Cap=(thsat-thres)·(a/n)*EXP(a·h/n)
	!********************************************************************************************************************

	pure function f2_hyd_cap_h_pw_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_pw_sca" :: f2_hyd_cap_h_pw_sca
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	if (h<0.0) then
		rout = (material%thsat-material%thres)*(material%a/material%n)*exp(material%a*h/material%n)
	else
		rout = (material%thsat-material%thres)*(material%a/material%n)
	end if

	end function f2_hyd_cap_h_pw_sca


	!********************************************************************************************************************
	!	F: f2_hyd_cap_h_PW_sca(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water capacity (dth/dh) (vector) from pressure head (vector pressure), with simplified power
	! function:
	!		Cap=(thsat-thres)·(a/n)*EXP(a·h/n)
	!********************************************************************************************************************

	pure function f2_hyd_cap_h_pw_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_pw_vec" :: f2_hyd_cap_h_pw_vec
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = (material%thsat-material%thres)*(material%a/material%n)*exp(material%a*h/material%n)
	else where
		rout = (material%thsat-material%thres)*(material%a/material%n)
	end where

	end function f2_hyd_cap_h_pw_vec

	!------

	pure function f2_hyd_cap_h_pw_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_pw_vec2" :: f2_hyd_cap_h_pw_vec2
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = (material%thsat-material%thres)*(material%a/material%n)*exp(material%a*h/material%n)
	else where
		rout = (material%thsat-material%thres)*(material%a/material%n)
	end where

	end function f2_hyd_cap_h_pw_vec2

	!********************************************************************************************************************
	!	F: f2_hyd_dkr_h_PW_sca(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of relative permeability (scalar) over pressure head (scalar pressure)
	! with simplified power function:
	!		dKr/dh = a·EXP(a·h)
	!********************************************************************************************************************

	pure function f2_hyd_dkr_h_pw_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_pw_sca" :: f2_hyd_dkr_h_pw_sca
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	if (h<0.0) then
		rout = material%a*exp(material%a*h)
	else
		rout = 1.0_dps
	end if

	end function f2_hyd_dkr_h_pw_sca


	!********************************************************************************************************************
	!	F: f2_hyd_dkr_h_PW_vec(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of relative permeability (vector) over pressure head (vector pressure)
	! with simplified power function:
	!		dKr/dh = a·EXP(a·h)
	!********************************************************************************************************************

	pure function f2_hyd_dkr_h_pw_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_pw_vec" :: f2_hyd_dkr_h_pw_vec
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	!kr = exp(a·h)
	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = material%a*exp(material%a*h)
	else where
		rout = 1.0_dps
	end where

	end function f2_hyd_dkr_h_pw_vec

	!-------------------

	pure function f2_hyd_dkr_h_pw_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_pw_vec2" :: f2_hyd_dkr_h_pw_vec2
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	!kr = exp(a·h)
	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = material%a*exp(material%a*h)
	else where
		rout = 1.0_dps
	end where

	end function f2_hyd_dkr_h_pw_vec2

	!********************************************************************************************************************
	!	F: f2_hyd_dk_h_PW_sca(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of permeability (scalar) over pressure head (scalar pressure)
	! with simplified power function:
	!		dK/dh = ksat·dKr/dh
	!********************************************************************************************************************

	pure function f2_hyd_dk_h_pw_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_pw_sca" :: f2_hyd_dk_h_pw_sca
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	rout = material%ksat*f2_hyd_dkr_h_pw_sca(h,material)

	end function f2_hyd_dk_h_pw_sca


	!********************************************************************************************************************
	!	F: f2_hyd_dk_h_PW_vec(h,material)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of permeability (vector) over pressure head (vector pressure)
	! with simplified power function:
	!		dK/dh = ksat·dKr/dh
	!********************************************************************************************************************

	pure function f2_hyd_dk_h_pw_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_pw_vec" :: f2_hyd_dk_h_pw_vec
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dps)::rout(size(h))

	rout = material%ksat*f2_hyd_dkr_h_pw_vec(h,material)

	end function f2_hyd_dk_h_pw_vec

	!------------------

	pure function f2_hyd_dk_h_pw_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_pw_vec2" :: f2_hyd_dk_h_pw_vec2
	!DEC$ endif


	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	rout = material%ksat*f2_hyd_dkr_h_pw_vec2(h,material)

	end function f2_hyd_dk_h_pw_vec2


	end module com_mod_hyd_pw