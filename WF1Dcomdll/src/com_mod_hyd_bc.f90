	!********************************************************************************************************************
	!*                                                                                                                  *
	!*                                 MODULE: MOD_UNSAT_HYDRAULICFUN_BC                                                *
	!*                                                                                                                  *
	!********************************************************************************************************************
	!* THIS MODULE INCLUDE VARIABLY SATURATED HYDRAULIC FUNCTIONS BY VAN GENUCHTEN AND MUALEM                           *
	!*    f2_hyd_S_h_BC:   Returns the specific saturation from pressure head                                            *
	!*    f2_hyd_th_h_BC:  Returns water content from pressure head                                                      *
	!*    f2_hyd_kr_h_BC:  Returns relative permeability from pressure head                                              *
	!*    f2_hyd_k_h_BC:   Returns permeability from pressure head                                                       *
	!*    f2_hyd_cap_h_BC: Returns water capacity from pressure head                                                     *
	!*                                                                                                                  *
	!*    Iván Campos-Guereta Díez                                                                                      *
	!*    MSc Civil Engineering by Polytechnic University of Madrid                                                     *
	!*    PhD Student by University of Nottingham                                                                       *
	!*    eMBA by International Institute San Telmo in Seville                                                          *
	!*    ivan.camposguereta@nottingham.ac.uk                                                                           *
	!*                                                                                                                  *
	!*    This software is copyrighted 2019(C)                                                                          *
	!********************************************************************************************************************

	module com_mod_hyd_bc

	implicit none
	private
	include 'inc_precision.fi'

	! !********************************************************************************************************************
	! ! POLIMORPHIC INTERFACES FOR THE FOLLOWING HYDRAULIC FUNCIONS
	! !--------------------------------------------------------------------------------------------------------------------
	! ! f2_hyd_S_h_BC	:Returns Equivalent saturation (sca or vec) from pressure head (sca or vec)
	! ! f2_hyd_th_h_BC	:Returns water content (sca or vec)					from pressure head (sca or vec)
	! ! f2_hyd_kr_h_BC	:Returns relative permeability (sca or vec) from pressure head (sca or vec)
	! ! f2_hyd_k_h_BC	:Returns permeability (sca or vec)					from pressure head (sca or vec)
	! ! f2_hyd_cap_h_BC:Returns water capacity (sca or vec)				from pressure head (sca or vec)
	! !********************************************************************************************************************
	!
	! interface f2_hyd_s_h_bc
	! module procedure f2_hyd_s_h_bc_sca
	! module procedure f2_hyd_s_h_bc_vec
	! module procedure f2_hyd_s_h_bc_vec2
	! end interface f2_hyd_s_h_bc
	!
	! interface f2_hyd_th_h_bc
	! module procedure f2_hyd_th_h_bc_sca
	! module procedure f2_hyd_th_h_bc_vec
	!module procedure f2_hyd_th_h_bc_vec2
	! end interface f2_hyd_th_h_bc
	!
	! interface f2_hyd_kr_h_bc
	! module procedure f2_hyd_kr_h_bc_sca
	! module procedure f2_hyd_kr_h_bc_vec
	!module procedure f2_hyd_kr_h_bc_vec2
	! end interface f2_hyd_kr_h_bc
	!
	! interface f2_hyd_k_h_bc
	! module procedure f2_hyd_k_h_bc_sca
	! module procedure f2_hyd_k_h_bc_vec
	!module procedure f2_hyd_k_h_bc_vec2
	! end interface f2_hyd_k_h_bc
	!
	! interface f2_hyd_cap_h_bc
	! module procedure f2_hyd_cap_h_bc_sca
	! module procedure f2_hyd_cap_h_bc_vec
	!module procedure f2_hyd_cap_h_bc_vec2
	!end interface f2_hyd_cap_h_bc
	!
	! interface f2_hyd_dkr_h_bc
	! module procedure f2_hyd_dkr_h_bc_sca
	! module procedure f2_hyd_dkr_h_bc_vec
	!module procedure f2_hyd_dkr_h_bc_vec2
	!end interface f2_hyd_dkr_h_bc
	!
	!interface f2_hyd_dk_h_bc
	! module procedure f2_hyd_dk_h_bc_sca
	! module procedure f2_hyd_dk_h_bc_vec
	!module procedure f2_hyd_dk_h_bc_vec2
	! end interface f2_hyd_dk_h_bc
	!
	! public::f2_hyd_s_h_bc,   f2_hyd_th_h_bc,  f2_hyd_kr_h_bc,  f2_hyd_k_h_bc,  f2_hyd_cap_h_bc,  f2_hyd_dkr_h_bc,  f2_hyd_dk_h_bc, f2_hyd_s0_h_bc_sca, f2_hyd_incs_h1_to_h2_bc_sca

	PUBLIC::  f2_hyd_S_h_bc_sca,  f2_hyd_S_h_bc_vec,  f2_hyd_S_h_bc_vec2, &
		f2_hyd_th_h_bc_sca, f2_hyd_th_h_bc_vec, f2_hyd_th_h_bc_vec2,  &
		f2_hyd_kr_h_bc_sca, f2_hyd_kr_h_bc_vec, f2_hyd_kr_h_bc_vec2, &
		f2_hyd_k_h_bc_sca,  f2_hyd_k_h_bc_vec,  f2_hyd_k_h_bc_vec2,  &
		f2_hyd_cap_h_bc_sca,f2_hyd_cap_h_bc_vec,f2_hyd_cap_h_bc_vec2, &
		f2_hyd_dkr_h_bc_sca,f2_hyd_dkr_h_bc_vec,f2_hyd_dkr_h_bc_vec2, &
		f2_hyd_dk_h_bc_sca, f2_hyd_dk_h_bc_vec, f2_hyd_dk_h_bc_vec2, &
		f2_hyd_incs_h1_to_h2_bc_sca



	contains


	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_SCA(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		Se=(1+(-a·h))^n)^(-m) for h<0 | Se=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_s_h_bc_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_bc_sca" :: f2_hyd_s_h_bc_sca
	!DEC$ endif

	!return the specific saturation from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	if (h<0.0_dps) then
		rout = (h/(-material%a))**(-material%n)
	else
		rout =  1.0_dps
	end if

	end function f2_hyd_s_h_bc_sca

	!********************************************************************************************************************
	! F: f2_hyd_S_H_PW_SCA(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (scalar) from pressure head (scalar pressure), with simplified power
	! function:
	!		Se=(1+(-a·h))^n)^(-m) for h<0 | Se=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_incs_h1_to_h2_bc_sca(h1,h2,material)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_incs_h1_to_h2_bc_sca" :: f2_hyd_incs_h1_to_h2_bc_sca
	!DEC$ endif

	!CHECK:IMPORTANT: In this case the increment S is not between h1 and h2 but between h1+psi_b and h2+psi_b, as the watertable in reality begin in psi_b (material%a)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h1,h2
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::f2_hyd_incs_h1_to_h2_bc_sca
	real(kind=16)::a,hb1,hb2,n,m,rout,rout1,rout2

	a = real(material%a,16)
	hb1 = real(h1,16)
	hb2 = real(h2,16)
	n= real(material%n,16)
	m= real(material%m,16)

	if (hb1<0.0_16) then
		rout1 = (hb1/(-material%a))**(-material%n) !CHECK: If this need to be hb1+material%a instead of hb1
	else
		rout1 =  1.0_16
	end if

	if (hb2<0.0_16) then
		rout2 = (hb2/(-material%a))**(-material%n)
	else
		rout2 =  1.0_16
	end if

	rout = rout2-rout1

	f2_hyd_incs_h1_to_h2_bc_sca = real(rout,dpd)

	end function f2_hyd_incs_h1_to_h2_bc_sca


	!********************************************************************************************************************
	! F: f2_hyd_S_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns specific saturation (vector) from pressure head (vector pressure), with Van Genuchten
	!	expressions:
	!		Se=(1+(-a·h))^n)^(-m) for h<0 | Se=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_s_h_bc_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_bc_vec" :: f2_hyd_s_h_bc_vec
	!DEC$ endif

	!return the specific saturation from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	withsuction:where (h<0.0_dps)
		rout = (h/(-material%a))**(-material%n)
	elsewhere
		rout =  1.0_dps
	end where withsuction

	end function f2_hyd_s_h_bc_vec

	!-------

	pure function f2_hyd_s_h_bc_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_s_h_bc_vec2" :: f2_hyd_s_h_bc_vec2
	!DEC$ endif

	!return the specific saturation from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	withsuction:where (h<0.0_dps)
		rout = (h/(-material%a))**(-material%n)
	elsewhere
		rout =  1.0_dps
	end where withsuction

	end function f2_hyd_s_h_bc_vec2

	!********************************************************************************************************************
	! F: f2_hyd_th_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (scalar) from pressure head (vector pressure), with Van Genuchten
	! expressions:
	!		th = thres+(thsat-thres)*Se
	!********************************************************************************************************************

	pure function f2_hyd_th_h_bc_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_bc_sca" :: f2_hyd_th_h_bc_sca
	!DEC$ endif

	!returns the moisture content from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_bc_sca(h,material)

	end function f2_hyd_th_h_bc_sca


	!********************************************************************************************************************
	! F: f2_hyd_th_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water content (scalar) from pressure head (vector pressure), with Van Genuchten
	! expressions:
	!		th = thres+(thsat-thres)*Se
	!********************************************************************************************************************

	pure function f2_hyd_th_h_bc_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_bc_vec" :: f2_hyd_th_h_bc_vec
	!DEC$ endif

	!returns the moisture content from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_bc_vec(h,material)

	end function f2_hyd_th_h_bc_vec

	!------------

	pure function f2_hyd_th_h_bc_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_th_h_bc_vec2" :: f2_hyd_th_h_bc_vec2
	!DEC$ endif

	!returns the moisture content from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	rout = material%thres+(material%thsat-material%thres)*f2_hyd_s_h_bc_vec2(h,material)

	end function f2_hyd_th_h_bc_vec2

	!********************************************************************************************************************
	! F: f2_hyd_kr_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (scalar) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		kr= Se^l·(1-(1-Se^(1/m))^m)^2 | Kr=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_kr_h_bc_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_bc_sca" :: f2_hyd_kr_h_bc_sca
	!DEC$ endif

	!returns relative permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout,se

	if (h<0.0_dps) then
		se = f2_hyd_s_h_bc_sca(h,material)
		rout = se**(2.0_dpd+material%l+2.0_dpd/material%n)
	else
		rout =  1.0_dpd
	end if


	end function f2_hyd_kr_h_bc_sca


	!********************************************************************************************************************
	! F: f2_hyd_kr_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns relative permeability (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		kr= Se^l·(1-(1-Se^(1/m))^m)^2 | Kr=1 for h>=0
	!********************************************************************************************************************

	pure function f2_hyd_kr_h_bc_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_bc_vec" :: f2_hyd_kr_h_bc_vec
	!DEC$ endif

	!returns relative permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h)),se(size(h))

	withsuction:where (h<0.0_dps)
		se = f2_hyd_s_h_bc_vec(h,material)
		rout = se**(2.0_dpd+material%l+2.0_dpd/material%n)
	elsewhere
		rout =  1.0_dps
	end where withsuction

	end function f2_hyd_kr_h_bc_vec

	!--------------

	pure function f2_hyd_kr_h_bc_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_kr_h_bc_vec2" :: f2_hyd_kr_h_bc_vec2
	!DEC$ endif

	!returns relative permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h)),se(size(h))

	withsuction:where (h<0.0_dps)
		se = f2_hyd_s_h_bc_vec2(h,material)
		rout = se**(2.0_dpd+material%l+2.0_dpd/material%n)
	elsewhere
		rout =  1.0_dps
	end where withsuction

	end function f2_hyd_kr_h_bc_vec2

	!********************************************************************************************************************
	! F: f2_hyd_k_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns permeability (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		k= ksat·kr
	!********************************************************************************************************************

	pure function f2_hyd_k_h_bc_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_bc_sca" :: f2_hyd_k_h_bc_sca
	!DEC$ endif

	!returns permeability from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	rout = material%ksat*f2_hyd_kr_h_bc_sca(h,material)

	end function f2_hyd_k_h_bc_sca


	!********************************************************************************************************************
	! F: f2_hyd_k_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns permeability (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		k= ksat·kr
	!********************************************************************************************************************

	pure function f2_hyd_k_h_bc_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_bc_vec" :: f2_hyd_k_h_bc_vec
	!DEC$ endif

	!returns permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	rout = material%ksat*f2_hyd_kr_h_bc_vec(h,material)

	end function f2_hyd_k_h_bc_vec

	!-----------

	pure function f2_hyd_k_h_bc_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_k_h_bc_vec2" :: f2_hyd_k_h_bc_vec2
	!DEC$ endif

	!returns permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	rout = material%ksat*f2_hyd_kr_h_bc_vec2(h,material)

	end function f2_hyd_k_h_bc_vec2

	!********************************************************************************************************************
	! F: f2_hyd_cap_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water capacity (scalar) from pressure head (scalar pressure), with
	!	Mualem-Van Genuchten expressions:
	!		Cap = -(thsat-thres)·m·n·(-a·h)^n·(1+(-h·a)^n)^(-1-m)/h (Problems when very near to saturation)
	!		Cap = -(thsat-thres)·m·n·a^n·(-h)^(n-1)·(1+(-h·a)^n)^(-1-m) (Avoid problems with n>1)
	!********************************************************************************************************************

	pure function f2_hyd_cap_h_bc_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_bc_sca" :: f2_hyd_cap_h_bc_sca
	!DEC$ endif

	!returns water capacity from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	if (h<0.0_dps) then
		!need to check the sign of cap!!!!!!!!!!!!!!!!!!!!
		rout =((material%thres-material%thsat)*material%n*(h/(-material%a))**(-material%n))/h
	else
		rout =  0.0_dps
	end if


	end function f2_hyd_cap_h_bc_sca


	!********************************************************************************************************************
	! F: f2_hyd_cap_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns water capacity (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		Cap = -(thsat-thres)·m·n·(-a·h)^n·(1+(-h·a)^n)^(-1-m)/h (Problems when very near to saturation)
	!		Cap = -(thsat-thres)·m·n·a^n·(-h)^(n-1)·(1+(-h·a)^n)^(-1-m) (Avoid problems with n>1)
	!********************************************************************************************************************

	pure function f2_hyd_cap_h_bc_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_bc_vec" :: f2_hyd_cap_h_bc_vec
	!DEC$ endif

	!returns water capacity from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	withsuction:where (h<0.0_dps)
		!need to check the sign of cap!!!!!!!!!!!!!!!!!!!!
		rout = ((material%thres-material%thsat)*material%n*(h/(-material%a))**(-material%n))/h
	elsewhere
		rout =  0.0_dps
	end where withsuction

	end function f2_hyd_cap_h_bc_vec

	!-------

	pure function f2_hyd_cap_h_bc_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_cap_h_bc_vec2" :: f2_hyd_cap_h_bc_vec2
	!DEC$ endif

	!returns water capacity from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	withsuction:where (h<0.0_dps)
		!need to check the sign of cap!!!!!!!!!!!!!!!!!!!!
		rout = ((material%thres-material%thsat)*material%n*(h/(-material%a))**(-material%n))/h
	elsewhere
		rout =  0.0_dps
	end where withsuction

	end function f2_hyd_cap_h_bc_vec2

	!********************************************************************************************************************
	! F: f2_hyd_dkrcap_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of relative permeability over h
	!	Mualem-Van Genuchten expressions:
	!		dkr/dh = m·n·x^[-1-(1+l)·m]·{l·x^m+(2+l-l·x)·y^(m-1)}·[1-(y/x)^m]·a^n·(-h)^(n-1)
	!   x= 1+(-a·h)^n
	!   y= x-1
	!********************************************************************************************************************

	pure function f2_hyd_dkr_h_bc_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_bc_sca" :: f2_hyd_dkr_h_bc_sca
	!DEC$ endif

	!returns water capacity from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout,x,y

	if (h<0.0_dps) then
		x = 1.0_dpd+(-material%a*h)**material%n
		y = x-1.0_dpd
		rout = -((2.0_dpd+(2.0_dpd+material%l)*material%n)*((h/(-material%a))**(-material%n))**(2.0_dpd+material%l+2.0_dpd/material%n))/h		
	else
		rout =  0.0_dps
	end if


	end function f2_hyd_dkr_h_bc_sca

	!********************************************************************************************************************
	! F: f2_hyd_dkrcap_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of relative permeability over h (vector)
	!	Mualem-Van Genuchten expressions:
	!		dkr/dh = m·n·x^[-1-(1+l)·m]·{l·x^m+(2+l-l·x)·y^(m-1)}·[1-(y/x)^m]·a^n·(-h)^(n-1)
	!   x= 1+(-a·h)^n
	!   y= x-1
	!********************************************************************************************************************

	pure function f2_hyd_dkr_h_bc_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_bc_vec" :: f2_hyd_dkr_h_bc_vec
	!DEC$ endif

	!returns water capacity from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h)),x(size(h)),y(size(h))

	withsuction:where (h<0.0_dps)
		x = 1+(-material%a*h)**material%n
		y = x-1
		rout = -((2.0_dpd+(2.0_dpd+material%l)*material%n)*((h/(-material%a))**(-material%n))**(2.0_dpd+material%l+2.0_dpd/material%n))/h
	elsewhere
		rout =  0.0_dps
	end where withsuction

	end function f2_hyd_dkr_h_bc_vec

	!------------------

	pure function f2_hyd_dkr_h_bc_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dkr_h_bc_vec2" :: f2_hyd_dkr_h_bc_vec2
	!DEC$ endif

	!returns water capacity from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h)),x(size(h)),y(size(h))

	withsuction:where (h<0.0_dps)
		x = 1+(-material%a*h)**material%n
		y = x-1
		rout = -((2.0_dpd+(2.0_dpd+material%l)*material%n)*((h/(-material%a))**(-material%n))**(2.0_dpd+material%l+2.0_dpd/material%n))/h
	elsewhere
		rout =  0.0_dps
	end where withsuction

	end function f2_hyd_dkr_h_bc_vec2

	!********************************************************************************************************************
	! F: f2_hyd_dk_h_VG_sca(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns permeability (vector) from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		k= ksat·kr
	!********************************************************************************************************************

	pure function f2_hyd_dk_h_bc_sca(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_bc_sca" :: f2_hyd_dk_h_bc_sca
	!DEC$ endif

	!returns permeability from pressure head(scalar)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout

	rout = material%ksat*f2_hyd_dkr_h_bc_sca(h,material)

	end function f2_hyd_dk_h_bc_sca


	!********************************************************************************************************************
	! F: f2_hyd_dk_h_VG_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns derivative of permeability (vector) respect h from pressure head (vector pressure), with
	!	Mualem-Van Genuchten expressions:
	!		k= ksat·kr
	!********************************************************************************************************************

	pure function f2_hyd_dk_h_bc_vec(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_bc_vec" :: f2_hyd_dk_h_bc_vec
	!DEC$ endif

	!returns permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material(:)
	real(kind=dpd)::rout(size(h))

	rout = material%ksat*f2_hyd_dkr_h_bc_vec(h,material)

	end function f2_hyd_dk_h_bc_vec

	!-------------

	pure function f2_hyd_dk_h_bc_vec2(h,material) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f2_hyd_dk_h_bc_vec2" :: f2_hyd_dk_h_bc_vec2
	!DEC$ endif

	!returns permeability from pressure head(vector)

	use com_mod_ty_material,only:ty_com_material

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dpd)::rout(size(h))

	rout = material%ksat*f2_hyd_dkr_h_bc_vec2(h,material)

	end function f2_hyd_dk_h_bc_vec2


	end module com_mod_hyd_bc