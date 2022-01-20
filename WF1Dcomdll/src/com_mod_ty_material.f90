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
	
	end type ty_com_material
	
	contains
	
	!********************************************************************************************************************
	! F: F_HYD_S_H_PW_SCA(H,MATERIAL)
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
       rout = f_hyd_s_h_pw_sca() 
       case default
       rout = f_hyd_s_h_vg_sca()
			 end select
			 
		contains
		!-------

		pure function f_hyd_s_h_pw_sca() result(rout)
		real(kind=dps)::rout

		if (h<0.0) then
			rout = merge(exp(material%a*h/material%n),1.0_dps,h<0.0_dps)
		else
			rout = 1.0_dps
		end if

		end function f_hyd_s_h_pw_sca

		!-----------

		pure function f_hyd_s_h_vg_sca() result(rout)
		!return the specific saturation from pressure head(scalar)
		real(kind=dpd)::rout

		if (h<0.0_dps) then
			rout = (1.0_dps+(-material%a*h)**material%n)**(-material%m)
		else
			rout =  1.0_dps
		end if

		end function f_hyd_s_h_vg_sca


		end function f_hyd_s_h_sca
		
	!********************************************************************************************************************
	! F: f_hyd_S_h_vec(H,MATERIAL)
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
		if (material%kindmat==2) then
			rout = f_hyd_s_h_pw_vec2()
		else
			rout = f_hyd_s_h_vg_vec2()
		end if

		contains

		pure function f_hyd_s_h_pw_vec2() result(rout)
		real(kind=dps)::rout(size(h))

		where (h<0.0)
			rout = merge(exp(material%a*h/material%n),1.0_dps,h<0.0_dps)
		else where
			rout = 1.0_dps
		end where

		end function f_hyd_s_h_pw_vec2

		!----------------

		pure function f_hyd_s_h_vg_vec2() result(rout)
		!return the specific saturation from pressure head(vector)

		real(kind=dpd)::rout(size(h))

		withsuction:where (h<0.0_dps)
			rout = (1.0_dps+(-material%a*h)**material%n)**(-material%m)
		elsewhere
			rout =  1.0_dps
		end where withsuction

		end function f_hyd_s_h_vg_vec2

		end function f_hyd_s_h_vec2		
		
	!********************************************************************************************************************
	! F: f_hyd_th_h_sca(H,MATERIAL)
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
			rout = f_hyd_th_h_pw_sca()
			case default
			rout = f_hyd_th_h_vg_sca()
		end select

		contains
		
		!---------------------
		
	pure function f_hyd_th_h_pw_sca()
	real(kind=dps)::f_hyd_th_h_pw_sca

	if (h<0.0) then
		f_hyd_th_h_pw_sca = material%thres+(material%thsat-material%thres)*material%get_s_sca(h)
	else
		f_hyd_th_h_pw_sca = material%thsat
	end if

	end function f_hyd_th_h_pw_sca
	
	!-------------
	
  pure function f_hyd_th_h_vg_sca()
  !returns the moisture content from pressure head(scalar)
  real(kind=dpd)::f_hyd_th_h_vg_sca

  f_hyd_th_h_vg_sca = material%thres+(material%thsat-material%thres)*material%get_s_sca(h)
	
	!--------

	end function f_hyd_th_h_vg_sca	
			
		
	
		END FUNCTION f_hyd_th_h_sca	
	
	
	!********************************************************************************************************************
	! F: f_hyd_th_h_vec(H,MATERIAL)
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
    
    if (material%KindMat==2) then
       Rout = f_hyd_th_h_pw_vec2(h,material)
       ELSE
       Rout = f_hyd_th_h_vg_vec2(h,material)
			 END if
			 
		contains
		
		pure function f_hyd_th_h_pw_vec2(h,material) result(rout)
 
	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))
 
	where (h<0.0)
		rout = material%thres+(material%thsat-material%thres)*material%get_s_vec(h)
	else where
		rout = material%thsat
	end where
 
		end function f_hyd_th_h_pw_vec2
		
		!----
		
	pure function f_hyd_th_h_vg_vec2(h,material) result(rout)
  !returns the moisture content from pressure head(vector)
 
  real(kind=dps),intent(in)::h(:)
  type(ty_com_material),intent(in)::material
  real(kind=dpd)::rout(size(h))
 
  rout = material%thres+(material%thsat-material%thres)*material%get_s_vec(h)
 
	end function f_hyd_th_h_vg_vec2
	
 
		END FUNCTION f_hyd_th_h_vec2
	
		
	!********************************************************************************************************************
	! F: f_hyd_kr_h_sca(H,MATERIAL)
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
			rout = f_hyd_kr_h_pw_sca(h,material)
			case default
			rout = f_hyd_kr_h_vg_sca(h,material)
		end select

		contains

		!-------

		pure function f_hyd_kr_h_pw_sca(h,material) result(rout)

		real(kind=dps),intent(in)::h
		type(ty_com_material),intent(in)::material
		real(kind=dps)::rout

		if (h<0.0) then
			rout = exp(material%a*h)
		else
			rout = 1.0_dps
		end if

		end function f_hyd_kr_h_pw_sca

		!-------

		pure function f_hyd_kr_h_vg_sca(h,material) result(rout)
		!returns relative permeability from pressure head(vector)

		real(kind=dps),intent(in)::h
		type(ty_com_material),intent(in)::material
		real(kind=dpd)::rout,se

		if (h<0.0_dps) then
			se = material%get_s_sca(h)
			rout = (se**material%l)*(1.0_dpd-(1.0_dpd-se**(1.0_dpd/material%m))**material%m)**2
		else
			rout =  1.0_dpd
		end if

		end function f_hyd_kr_h_vg_sca

		end function f_hyd_kr_h_sca

	!********************************************************************************************************************
	! F: f_hyd_kr_h_vec(H,MATERIAL)
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
		if (material%kindmat==2) then
			rout = f_hyd_kr_h_pw_vec2(h,material)
		else
			rout = f_hyd_kr_h_vg_vec2(h,material)
		end if

		contains

		!-----

		pure function f_hyd_kr_h_pw_vec2(h,material) result(rout)

		!kr = exp(a·h)
		real(kind=dps),intent(in)::h(:)
		type(ty_com_material),intent(in)::material
		real(kind=dps)::rout(size(h))

		where (h<0.0)
			rout = exp(material%a*h)
		else where
			rout = 1.0_dps
		end where

		end function f_hyd_kr_h_pw_vec2

		!---------

		pure function f_hyd_kr_h_vg_vec2(h,material) result(rout)
		!returns relative permeability from pressure head(vector)

		real(kind=dps),intent(in)::h(:)
		type(ty_com_material),intent(in)::material
		real(kind=dpd)::rout(size(h)),se(size(h))

		withsuction:where (h<0.0_dps)
			!se = material%get_s_vec(h) !PREVIOUS
			se = f_hyd_s_h_vec2(material,h)
			rout = (se**material%l)*(1.0_dpd-(1.0_dpd-se**(1.0_dpd/material%m))**material%m)**2
		elsewhere
			rout =  1.0_dps
		end where withsuction

		end function f_hyd_kr_h_vg_vec2

		end function f_hyd_kr_h_vec2
		
	!********************************************************************************************************************
	! F: f_hyd_k_h_sca(H,MATERIAL)
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
			rout = f_hyd_k_h_pw_sca(h,material)
			case default
			rout = f_hyd_k_h_vg_sca(h,material)
		end select

		contains

		!----

		pure function f_hyd_k_h_pw_sca(h,material) result(rout)

		real(kind=dps),intent(in)::h
		type(ty_com_material),intent(in)::material
		real(kind=dps)::rout

		rout = material%ksat*material%get_kr_sca(h)

		end function f_hyd_k_h_pw_sca


		!-----

		pure function f_hyd_k_h_vg_sca(h,material) result(rout)
		!returns permeability from pressure head(scalar)

		real(kind=dps),intent(in)::h
		type(ty_com_material),intent(in)::material
		real(kind=dpd)::rout

		rout = material%ksat*material%get_kr_sca(h)

		end function f_hyd_k_h_vg_sca

		end function f_hyd_k_h_sca
		
	!********************************************************************************************************************
	! F: f_hyd_k_h_vec(H,MATERIAL)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that returns sthe permeability (vector) from pressure head (vector pressure) using the hyraulicfuncion
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
    if (material%KindMat==2) then
       Rout = f_hyd_k_h_PW_vec2(h,material)
       ELSE
       Rout = f_hyd_k_h_VG_vec2(h,material)
    END if    
    
		contains
		
		!-------
		
			pure function f_hyd_k_h_pw_vec2(h,material) result(rout)

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	rout = material%ksat*material%get_kr_vec(h)

			end function f_hyd_k_h_pw_vec2
			
			!------
			
	pure function f_hyd_k_h_vg_vec2(h,material) result(rout)
  !returns permeability from pressure head(vector)

  real(kind=dps),intent(in)::h(:)
  type(ty_com_material),intent(in)::material
  real(kind=dpd)::rout(size(h))

  rout = material%ksat*material%get_kr_vec(h)

	end function f_hyd_k_h_vg_vec2
	
		END FUNCTION f_hyd_k_h_vec2		
		
	!********************************************************************************************************************
	! F: f_hyd_Cap_h_sca(H,MATERIAL)
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
			rout = f_hyd_cap_h_pw_sca(h,material)
			case default
			rout = f_hyd_cap_h_vg_sca(h,material)
		end select

		contains
		!-------

		pure function f_hyd_cap_h_pw_sca(h,material) result(rout)

		real(kind=dps),intent(in)::h
		type(ty_com_material),intent(in)::material
		real(kind=dps)::rout

		if (h<0.0) then
			rout = (material%thsat-material%thres)*(material%a/material%n)*exp(material%a*h/material%n)
		else
			rout = (material%thsat-material%thres)*(material%a/material%n)
		end if

		end function f_hyd_cap_h_pw_sca

		!-----

		pure function f_hyd_cap_h_vg_sca(h,material) result(rout)
		!returns water capacity from pressure head(scalar)

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

		end function f_hyd_cap_h_vg_sca


		end function f_hyd_cap_h_sca
		
	!********************************************************************************************************************
	! F: f_hyd_Cap_h_vec(H,MATERIAL)
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
    if (material%KindMat==2) then
       Rout = f_hyd_Cap_h_PW_vec2(h,material)
       ELSE 
       Rout = f_hyd_cap_h_VG_vec2(h,material)
			 END if       
    
		contains
			 !-------
		
	pure function f_hyd_cap_h_pw_vec2(h,material) result(rout)

	real(kind=dps),intent(in)::h(:)
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout(size(h))

	where (h<0.0)
		rout = (material%thsat-material%thres)*(material%a/material%n)*exp(material%a*h/material%n)
	else where
		rout = (material%thsat-material%thres)*(material%a/material%n)
	end where

	end function f_hyd_cap_h_pw_vec2			 
	
	!-------
	
	pure function f_hyd_cap_h_vg_vec2(h,material) result(rout)
  !returns water capacity from pressure head(vector)

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

	end function f_hyd_cap_h_vg_vec2	
			 
		END FUNCTION f_hyd_Cap_h_vec2		
		
	!********************************************************************************************************************
	! F: f_hyd_dkr_h_sca(H,MATERIAL)
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
       Rout = f_hyd_dkr_h_PW_sca(h,material)
       CASE DEFAULT
       Rout = f_hyd_dkr_h_VG_sca(h,material)
			 END SELECT   
		
		contains
		!----------------
		
	pure function f_hyd_dkr_h_pw_sca(h,material) result(rout)

	real(kind=dps),intent(in)::h
	type(ty_com_material),intent(in)::material
	real(kind=dps)::rout

	if (h<0.0) then
		rout = material%a*exp(material%a*h)
	else
		rout = 1.0_dps
	end if

	end function f_hyd_dkr_h_pw_sca
	
	!---------------
	
  pure function f_hyd_dkr_h_vg_sca(h,material) result(rout)
  !returns water capacity from pressure head(scalar)

  real(kind=dps),intent(in)::h
  type(ty_com_material),intent(in)::material
  real(kind=dpd)::rout,x,y

  if (h<0.0_dps) then
    x = 1+(-material%a*h)**material%n
		y = x-1
    rout = material%m*material%n*x**(-1-(1+material%l)*material%m)*(material%l*x**material%m+(2+material%l-material%l*x)*&
			&    y**(material%m-1))*(1-(y/x)**material%m)*material%a**material%n*(-h)**(material%n-1)
  else
    rout =  0.0_dps
  end if

	end function f_hyd_dkr_h_vg_sca	

		END FUNCTION f_hyd_dkr_h_sca	
		
	!********************************************************************************************************************
	! F: f_hyd_dkr_h_vec(H,MATERIAL)
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
		if (material%kindmat==2) then
			rout = f_hyd_dkr_h_pw_vec2(h,material)
		else
			rout = f_hyd_dkr_h_vg_vec2(h,material)
		end if

		contains

		pure function f_hyd_dkr_h_pw_vec2(h,material) result(rout)

		!kr = exp(a·h)
		real(kind=dps),intent(in)::h(:)
		type(ty_com_material),intent(in)::material
		real(kind=dps)::rout(size(h))

		where (h<0.0)
			rout = material%a*exp(material%a*h)
		else where
			rout = 1.0_dps
		end where

		end function f_hyd_dkr_h_pw_vec2

		!-----------

		pure function f_hyd_dkr_h_vg_vec2(h,material) result(rout)
		!returns water capacity from pressure head(scalar)

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

		end function f_hyd_dkr_h_vg_vec2

		end function f_hyd_dkr_h_vec2
		
	!********************************************************************************************************************
	! F: f_hyd_dk_h_sca(H,MATERIAL)
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
			rout = f_hyd_dk_h_pw_sca(h,material)
			case default
			rout = f_hyd_dk_h_vg_sca(h,material)
		end select

		contains

		!----

		pure function f_hyd_dk_h_pw_sca(h,material) result(rout)

		real(kind=dps),intent(in)::h
		type(ty_com_material),intent(in)::material
		real(kind=dps)::rout

		rout = material%ksat*f_hyd_dkr_h_sca(material,h)

		end function f_hyd_dk_h_pw_sca

		!-----

		pure function f_hyd_dk_h_vg_sca(h,material) result(rout)
		!returns permeability from pressure head(scalar)

		real(kind=dps),intent(in)::h
		type(ty_com_material),intent(in)::material
		real(kind=dpd)::rout

		rout = material%ksat*f_hyd_dkr_h_sca(material,h)

		end function f_hyd_dk_h_vg_sca

		end function f_hyd_dk_h_sca
		
	!********************************************************************************************************************
	! F: f_hyd_dk_h_vec(H,MATERIAL)
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
		if (material%kindmat==2) then
			rout = f_hyd_dk_h_pw_vec2(h,material)
		else
			rout = f_hyd_dk_h_vg_vec2(h,material)
		end if

		contains

		!------

		pure function f_hyd_dk_h_pw_vec2(h,material) result(rout)

		real(kind=dps),intent(in)::h(:)
		type(ty_com_material),intent(in)::material
		real(kind=dps)::rout(size(h))

		rout = material%ksat*f_hyd_dkr_h_vec2(material,h)

		end function f_hyd_dk_h_pw_vec2

		!----------

		pure function f_hyd_dk_h_vg_vec2(h,material) result(rout)
		!returns permeability from pressure head(vector)

		real(kind=dps),intent(in)::h(:)
		type(ty_com_material),intent(in)::material
		real(kind=dpd)::rout(size(h))

		rout = material%ksat*f_hyd_dkr_h_vec2(material,h)

		end function f_hyd_dk_h_vg_vec2

		end function f_hyd_dk_h_vec2
		

	end module com_mod_ty_material