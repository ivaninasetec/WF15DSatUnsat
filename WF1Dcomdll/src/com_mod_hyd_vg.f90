	!********************************************************************************************************************
	!*                                                                                                                  *
	!*                                 MODULE: MOD_UNSAT_HYDRAULICFUN_VG                                                *
	!*                                                                                                                  *
	!********************************************************************************************************************
	!* THIS MODULE INCLUDE VARIABLY SATURATED HYDRAULIC FUNCTIONS BY VAN GENUCHTEN AND MUALEM                           *
	!*    f2_hyd_S_h_VG:   Returns the specific saturation from pressure head                                            *
	!*    f2_hyd_th_h_VG:  Returns water content from pressure head                                                      *
	!*    f2_hyd_kr_h_VG:  Returns relative permeability from pressure head                                              *
	!*    f2_hyd_k_h_VG:   Returns permeability from pressure head                                                       *
	!*    f2_hyd_cap_h_VG: Returns water capacity from pressure head                                                     *
	!*                                                                                                                  *
	!*    Iv�n Campos-Guereta D�ez                                                                                      *
	!*    MSc Civil Engineering by Polytechnic University of Madrid                                                     *
	!*    PhD Student by University of Nottingham                                                                       *
	!*    eMBA by International Institute San Telmo in Seville                                                          *
	!*    ivan.camposguereta@nottingham.ac.uk                                                                           *
	!*                                                                                                                  *
	!*    This software is copyrighted 2019(C)                                                                          *
	!********************************************************************************************************************

	module com_mod_hyd_vg

	implicit none
	private
	include 'inc_precision.fi'

	! !********************************************************************************************************************
	! ! POLIMORPHIC INTERFACES FOR THE FOLLOWING HYDRAULIC FUNCIONS
	! !--------------------------------------------------------------------------------------------------------------------
	! ! f2_hyd_S_h_VG	:Returns Equivalent saturation (sca or vec) from pressure head (sca or vec)
	! ! f2_hyd_th_h_VG	:Returns water content (sca or vec)					from pressure head (sca or vec)
	! ! f2_hyd_kr_h_VG	:Returns relative permeability (sca or vec) from pressure head (sca or vec)
	! ! f2_hyd_k_h_VG	:Returns permeability (sca or vec)					from pressure head (sca or vec)
	! ! f2_hyd_cap_h_VG:Returns water capacity (sca or vec)				from pressure head (sca or vec)
	! !********************************************************************************************************************
	!
	! interface f2_hyd_s_h_vg
	! module procedure f2_hyd_s_h_vg_sca
	! module procedure f2_hyd_s_h_vg_vec
	! module procedure f2_hyd_s_h_vg_vec2
	! end interface f2_hyd_s_h_vg
	!
	! interface f2_hyd_th_h_vg
	! module procedure f2_hyd_th_h_vg_sca
	! module procedure f2_hyd_th_h_vg_vec
	!module procedure f2_hyd_th_h_vg_vec2
	! end interface f2_hyd_th_h_vg
	!
	! interface f2_hyd_kr_h_vg
	! module procedure f2_hyd_kr_h_vg_sca
	! module procedure f2_hyd_kr_h_vg_vec
	!module procedure f2_hyd_kr_h_vg_vec2
	! end interface f2_hyd_kr_h_vg
	!
	! interface f2_hyd_k_h_vg
	! module procedure f2_hyd_k_h_vg_sca
	! module procedure f2_hyd_k_h_vg_vec
	!module procedure f2_hyd_k_h_vg_vec2
	! end interface f2_hyd_k_h_vg
	!
	! interface f2_hyd_cap_h_vg
	! module procedure f2_hyd_cap_h_vg_sca
	! module procedure f2_hyd_cap_h_vg_vec
	!module procedure f2_hyd_cap_h_vg_vec2
	!end interface f2_hyd_cap_h_vg
	!
	! interface f2_hyd_dkr_h_vg
	! module procedure f2_hyd_dkr_h_vg_sca
	! module procedure f2_hyd_dkr_h_vg_vec
	!module procedure f2_hyd_dkr_h_vg_vec2
	!end interface f2_hyd_dkr_h_vg
	!
	!interface f2_hyd_dk_h_vg
	! module procedure f2_hyd_dk_h_vg_sca
	! module procedure f2_hyd_dk_h_vg_vec
	!module procedure f2_hyd_dk_h_vg_vec2
	! end interface f2_hyd_dk_h_vg
	!
	! public::f2_hyd_s_h_vg,   f2_hyd_th_h_vg,  f2_hyd_kr_h_vg,  f2_hyd_k_h_vg,  f2_hyd_cap_h_vg,  f2_hyd_dkr_h_vg,  f2_hyd_dk_h_vg, f2_hyd_s0_h_vg_sca, f2_hyd_incs_h1_to_h2_vg_sca

	PUBLIC::  f2_hyd_S_h_vg_sca,  f2_hyd_S_h_vg_vec,  f2_hyd_S_h_vg_vec2, &
		f2_hyd_th_h_vg_sca, f2_hyd_th_h_vg_vec, f2_hyd_th_h_vg_vec2,  &
		f2_hyd_kr_h_vg_sca, f2_hyd_kr_h_vg_vec, f2_hyd_kr_h_vg_vec2, &
		f2_hyd_k_h_vg_sca,  f2_hyd_k_h_vg_vec,  f2_hyd_k_h_vg_vec2,  &
		f2_hyd_cap_h_vg_sca,f2_hyd_cap_h_vg_vec,f2_hyd_cap_h_vg_vec2, &
		f2_hyd_dkr_h_vg_sca,f2_hyd_dkr_h_vg_vec,f2_hyd_dkr_h_vg_vec2, &
		f2_hyd_dk_h_vg_sca, f2_hyd_dk_h_vg_vec, f2_hyd_dk_h_vg_vec2, &
		f2_hyd_incs_h1_to_h2_vg_sca



	contains


	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_SCA(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		Se=(1+(-a�h))^n)^(-m) for h<0 | Se=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_s_h_vg_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_vg_sca" :: f2_hyd_s_h_vg_sca
	!DEC$ endif

	!return the specific saturation from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	if (h<0.0_dps) then
		rout = (1.0_dps+(-material%a*h)**material%n)**(-material%m)
	else
		rout =  1.0_dps
	end if

	end function f2_hyd_s_h_vg_sca

	!!********************************************************************************************************************
	!! F: f2_hyd_S_H_PW_SCA(H,MATERIAL)
	!!--------------------------------------------------------------------------------------------------------------------
	!! Function that returns specific saturation (scalar) from pressure head (scalar pressure), with simplified power
	!! function:
	!!		Se=(1+(-a�h))^n)^(-m) for h<0 | Se=1 for h>=0
	!!********************************************************************************************************************
 !
	!pure function f2_hyd_s0_h_vg_sca(h,material)
	!!DEC$ if defined(_DLL)
	!!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s0_h_vg_sca" :: f2_hyd_s0_h_vg_sca
	!!DEC$ endif
 !
	!!return the specific saturation from pressure head(scalar)
 !
	!use com_mod_ty_material,only:ty_com_material
 !
	!real(kind=dps),intent(in)::h
	!type(ty_com_material),intent(in)::material
	!real(kind=dpd)::f2_hyd_s0_h_vg_sca
	!real(kind=16)::a,hb,n,m,rout
 !
	!a = real(material%a,16)
	!hb = real(h,16)
	!n= real(material%n,16)
	!m= real(material%m,16)
 !
	!if (hb<0.0) then
	!	rout = (1.0_16+(-a*hb)**n)**(-m)-1.0_16
	!	f2_hyd_s0_h_vg_sca = real(rout,dpd)
	!else
	!	f2_hyd_s0_h_vg_sca =  0.0_dpd
	!end if
 !
	!end function f2_hyd_s0_h_vg_sca

	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_SCA(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		Se=(1+(-a�h))^n)^(-m) for h<0 | Se=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_incs_h1_to_h2_vg_sca(h1,h2,material)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_incs_h1_to_h2_vg_sca" :: f2_hyd_incs_h1_to_h2_vg_sca
	!DEC$ endif

	!return the specific saturation from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h1,h2
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::f2_hyd_incs_h1_to_h2_vg_sca
	real(kind=16)::a,hb1,hb2,n,m,rout,rout1,rout2

	a = real(material%a,16)
	hb1 = real(h1,16)
	hb2 = real(h2,16)
	n= real(material%n,16)
	m= real(material%m,16)

	if (hb1<0.0_16) then
		rout1 = (1.0_16+(-a*hb1)**n)**(-m)
	else
		rout1 =  1.0_16
	end if

	if (hb2<0.0_16) then
		rout2 = (1.0_16+(-a*hb2)**n)**(-m)
	else
		rout2 =  1.0_16
	end if

	rout = rout2-rout1

	f2_hyd_incs_h1_to_h2_vg_sca = real(rout,dpd)

	end function f2_hyd_incs_h1_to_h2_vg_sca


	!********************************************************************************************************************
	! F: f2_hyd_S_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (vector) from pressure head (vector pressure), with Van Genuchten
	!	expressions:
	!		Se=(1+(-a�h))^n)^(-m) for h<0 | Se=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_s_h_vg_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_vg_vec" :: f2_hyd_s_h_vg_vec
	!DEC$ endif

	!return the specific saturation from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	withsuction:where (h<0.0_dps)
		rout = (1.0_dps+(-material%a*h)**material%n)**(-material%m)
	elsewhere
		rout =  1.0_dps
	end where withsuction

	end function f2_hyd_s_h_vg_vec

	!-------

	pure function f2_hyd_s_h_vg_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_vg_vec2" :: f2_hyd_s_h_vg_vec2
	!DEC$ endif

	!return the specific saturation from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	withsuction:where (h<0.0_dps)
		rout = (1.0_dps+(-material%a*h)**material%n)**(-material%m)
	elsewhere
		rout =  1.0_dps
	end where withsuction

	end function f2_hyd_s_h_vg_vec2

	!********************************************************************************************************************
	! F: f2_hyd_th_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (scalar) from pressure head (vector pressure), with Van Genuchten
	! expressions:
	!		th = thres+(thsat-thres)*Se
	!********************************************************************************************************************

	pure function f2_hyd_th_h_vg_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_vg_sca" :: f2_hyd_th_h_vg_sca
	!DEC$ endif

	!returns the moisture content from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_vg_sca(h,material)

	end function f2_hyd_th_h_vg_sca


	!********************************************************************************************************************
	! F: f2_hyd_th_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (scalar) from pressure head (vector pressure), with Van Genuchten
	! expressions:
	!		th = thres+(thsat-thres)*Se
	!********************************************************************************************************************

	pure function f2_hyd_th_h_vg_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_vg_vec" :: f2_hyd_th_h_vg_vec
	!DEC$ endif

	!returns the moisture content from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_vg_vec(h,material)

	end function f2_hyd_th_h_vg_vec

	!------------

	pure function f2_hyd_th_h_vg_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_vg_vec2" :: f2_hyd_th_h_vg_vec2
	!DEC$ endif

	!returns the moisture content from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_vg_vec2(h,material)

	end function f2_hyd_th_h_vg_vec2

	!********************************************************************************************************************
	! F: f2_hyd_kr_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (scalar) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		kr= Se^l�(1-(1-Se^(1/m))^m)^2 | Kr=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_kr_h_vg_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_vg_sca" :: f2_hyd_kr_h_vg_sca
	!DEC$ endif

	!returns relative permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout,se

	if (h<0.0_dps) then
		se = f2_hyd_s_h_vg_sca(h,material)
		rout = (se**material%l)*(1.0_dpd-(1.0_dpd-se**(1.0_dpd/material%m))**material%m)**2
	else
		rout =  1.0_dpd
	end if


	end function f2_hyd_kr_h_vg_sca


	!********************************************************************************************************************
	! F: f2_hyd_kr_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		kr= Se^l�(1-(1-Se^(1/m))^m)^2 | Kr=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_kr_h_vg_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_vg_vec" :: f2_hyd_kr_h_vg_vec
	!DEC$ endif

	!returns relative permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h)),se(size(h))

	withsuction:where (h<0.0_dps)
		se = f2_hyd_s_h_vg_vec(h,material)
		rout = (se**material%l)*(1.0_dpd-(1.0_dpd-se**(1.0_dpd/material%m))**material%m)**2
	elsewhere
		rout =  1.0_dps
	end where withsuction

	end function f2_hyd_kr_h_vg_vec

	!--------------

	pure function f2_hyd_kr_h_vg_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_vg_vec2" :: f2_hyd_kr_h_vg_vec2
	!DEC$ endif

	!returns relative permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h)),se(size(h))

	withsuction:where (h<0.0_dps)
		se = f2_hyd_s_h_vg_vec2(h,material)
		rout = (se**material%l)*(1.0_dpd-(1.0_dpd-se**(1.0_dpd/material%m))**material%m)**2
	elsewhere
		rout =  1.0_dps
	end where withsuction

	end function f2_hyd_kr_h_vg_vec2

	!********************************************************************************************************************
	! F: f2_hyd_k_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns permeability (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		k= ksat�kr
	!********************************************************************************************************************

	pure function f2_hyd_k_h_vg_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_vg_sca" :: f2_hyd_k_h_vg_sca
	!DEC$ endif

	!returns permeability from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	rout = material%ksat*f2_hyd_kr_h_vg_sca(h,material)

	end function f2_hyd_k_h_vg_sca


	!********************************************************************************************************************
	! F: f2_hyd_k_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns permeability (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		k= ksat�kr
	!********************************************************************************************************************

	pure function f2_hyd_k_h_vg_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_vg_vec" :: f2_hyd_k_h_vg_vec
	!DEC$ endif

	!returns permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	rout = material%ksat*f2_hyd_kr_h_vg_vec(h,material)

	end function f2_hyd_k_h_vg_vec

	!-----------

	pure function f2_hyd_k_h_vg_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_vg_vec2" :: f2_hyd_k_h_vg_vec2
	!DEC$ endif

	!returns permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	rout = material%ksat*f2_hyd_kr_h_vg_vec2(h,material)

	end function f2_hyd_k_h_vg_vec2

	!********************************************************************************************************************
	! F: f2_hyd_cap_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water capacity (scalar) from pressure head (scalar pressure), with
	!	Mualem-Van Genuchten expressions:
	!		Cap = -(thsat-thres)�m�n�(-a�h)^n�(1+(-h�a)^n)^(-1-m)/h (Problems when very near to saturation)
	!		Cap = -(thsat-thres)�m�n�a^n�(-h)^(n-1)�(1+(-h�a)^n)^(-1-m) (Avoid problems with n>1)
	!********************************************************************************************************************

	pure function f2_hyd_cap_h_vg_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_vg_sca" :: f2_hyd_cap_h_vg_sca
	!DEC$ endif

	!returns water capacity from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	if (h<0.0_dps) then
		!need to check the sign of cap!!!!!!!!!!!!!!!!!!!!
		rout = (material%thsat-material%thres)*material%m*material%n*(material%a**material%n)*(-h)**(material%n-1)&
			&*(1.0_dpd+(-material%a*h)**material%n)**(-1.0_dpd-material%m)
		!rout = ((1.0_dps+(-material%a*h)**material%n)**(-material%m-1.0_dps))*(material%thsat-material%thres)&
		!*material%m*material%n*(material%a**material%n)*(-h)**(material%n-1)
	else
		rout =  0.0_dps
	end if


	end function f2_hyd_cap_h_vg_sca


	!********************************************************************************************************************
	! F: f2_hyd_cap_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water capacity (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		Cap = -(thsat-thres)�m�n�(-a�h)^n�(1+(-h�a)^n)^(-1-m)/h (Problems when very near to saturation)
	!		Cap = -(thsat-thres)�m�n�a^n�(-h)^(n-1)�(1+(-h�a)^n)^(-1-m) (Avoid problems with n>1)
	!********************************************************************************************************************

	pure function f2_hyd_cap_h_vg_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_vg_vec" :: f2_hyd_cap_h_vg_vec
	!DEC$ endif

	!returns water capacity from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	withsuction:where (h<0.0_dps)
		!need to check the sign of cap!!!!!!!!!!!!!!!!!!!!
		rout = (material%thsat-material%thres)*material%m*material%n*(material%a**material%n)*(-h)**(material%n-1)&
			&*(1.0_dpd+(-material%a*h)**material%n)**(-1.0_dpd-material%m)
	elsewhere
		rout =  0.0_dps
	end where withsuction

	end function f2_hyd_cap_h_vg_vec

	!-------

	pure function f2_hyd_cap_h_vg_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_vg_vec2" :: f2_hyd_cap_h_vg_vec2
	!DEC$ endif

	!returns water capacity from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	withsuction:where (h<0.0_dps)
		!need to check the sign of cap!!!!!!!!!!!!!!!!!!!!
		rout = (material%thsat-material%thres)*material%m*material%n*(material%a**material%n)*(-h)**(material%n-1)&
			&*(1.0_dpd+(-material%a*h)**material%n)**(-1.0_dpd-material%m)
	elsewhere
		rout =  0.0_dps
	end where withsuction

	end function f2_hyd_cap_h_vg_vec2

	!********************************************************************************************************************
	! F: f2_hyd_dkrcap_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of relative permeability over h
	!	Mualem-Van Genuchten expressions:
	!		dkr/dh = m�n�x^[-1-(1+l)�m]�{l�x^m+(2+l-l�x)�y^(m-1)}�[1-(y/x)^m]�a^n�(-h)^(n-1)
	!   x= 1+(-a�h)^n
	!   y= x-1
	!********************************************************************************************************************

	pure function f2_hyd_dkr_h_vg_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_vg_sca" :: f2_hyd_dkr_h_vg_sca
	!DEC$ endif

	!returns water capacity from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout,x,y

	if (h<0.0_dps) then
		x = 1.0_dpd+(-material%a*h)**material%n
		y = x-1.0_dpd
		rout = material%m*material%n*x**(-1.0_dpd-(1+material%l)*material%m)*(material%l*x**material%m+(2.0_dpd+material%l-material%l*x)*&
			&    y**(material%m-1.0_dpd))*(1.0_dpd-(y/x)**material%m)*material%a**material%n*(-h)**(material%n-1)
	else
		rout =  0.0_dps
	end if


	end function f2_hyd_dkr_h_vg_sca

	!********************************************************************************************************************
	! F: f2_hyd_dkrcap_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of relative permeability over h (vector)
	!	Mualem-Van Genuchten expressions:
	!		dkr/dh = m�n�x^[-1-(1+l)�m]�{l�x^m+(2+l-l�x)�y^(m-1)}�[1-(y/x)^m]�a^n�(-h)^(n-1)
	!   x= 1+(-a�h)^n
	!   y= x-1
	!********************************************************************************************************************

	pure function f2_hyd_dkr_h_vg_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_vg_vec" :: f2_hyd_dkr_h_vg_vec
	!DEC$ endif

	!returns water capacity from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h)),x(size(h)),y(size(h))

	withsuction:where (h<0.0_dps)
		x = 1+(-material%a*h)**material%n
		y = x-1
		rout = material%m*material%n*x**(-1-(1+material%l)*material%m)*(material%l*x**material%m+(2+material%l-material%l*x)*&
			&    y**(material%m-1))*(1-(y/x)**material%m)*material%a**material%n*(-h)**(material%n-1)
	elsewhere
		rout =  0.0_dps
	end where withsuction

	end function f2_hyd_dkr_h_vg_vec

	!------------------

	pure function f2_hyd_dkr_h_vg_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_vg_vec2" :: f2_hyd_dkr_h_vg_vec2
	!DEC$ endif

	!returns water capacity from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h)),x(size(h)),y(size(h))

	withsuction:where (h<0.0_dps)
		x = 1+(-material%a*h)**material%n
		y = x-1
		rout = material%m*material%n*x**(-1-(1+material%l)*material%m)*(material%l*x**material%m+(2+material%l-material%l*x)*&
			&    y**(material%m-1))*(1-(y/x)**material%m)*material%a**material%n*(-h)**(material%n-1)
	elsewhere
		rout =  0.0_dps
	end where withsuction

	end function f2_hyd_dkr_h_vg_vec2

	!********************************************************************************************************************
	! F: f2_hyd_dk_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns permeability (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		k= ksat�kr
	!********************************************************************************************************************

	pure function f2_hyd_dk_h_vg_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_vg_sca" :: f2_hyd_dk_h_vg_sca
	!DEC$ endif

	!returns permeability from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	rout = material%ksat*f2_hyd_dkr_h_vg_sca(h,material)

	end function f2_hyd_dk_h_vg_sca


	!********************************************************************************************************************
	! F: f2_hyd_dk_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of permeability (vector) respect h from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		k= ksat�kr
	!********************************************************************************************************************

	pure function f2_hyd_dk_h_vg_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_vg_vec" :: f2_hyd_dk_h_vg_vec
	!DEC$ endif

	!returns permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	rout = material%ksat*f2_hyd_dkr_h_vg_vec(h,material)

	end function f2_hyd_dk_h_vg_vec

	!-------------

	pure function f2_hyd_dk_h_vg_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_vg_vec2" :: f2_hyd_dk_h_vg_vec2
	!DEC$ endif

	!returns permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	rout = material%ksat*f2_hyd_dkr_h_vg_vec2(h,material)

	end function f2_hyd_dk_h_vg_vec2


	end module com_mod_hyd_vg