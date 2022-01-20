  module com_mod_ty_materials
	use com_mod_ty_material, only: ty_com_material
	
  implicit none
  include 'inc_precision.fi'
	
  private
	
	!******************************************************************************************************************
  ! TY_COM_MATERIALS
  ! Includes the collection of materials for the model
  !------------------------------------------------------------------------------------------------------------------
  !	 |Count           [integer]:        Number of different materials
  !  |Material(:)     [ty_com_material]:Array of objects [ty_com_material]
  !******************************************************************************************************************

  type,public::ty_com_materials !CLASS: Materials collection (common description)
    integer::count !< total number of materials
    type(ty_com_material),allocatable::material(:) !< material(nmat) nmat is number of materials.
		
	contains
		procedure,public:: allocateall	 => s_com_materials_allocateall
		procedure,public:: deallocateall	 => s_com_materials_deallocateall
		procedure,public:: get_s	 => f_com_materials_get_s_vec,f_com_materials_get_s_sca
		procedure,public:: get_th	 => f_com_materials_get_th_vec,f_com_materials_get_th_sca
		procedure,public:: get_k	 => f_com_materials_get_k_vec,f_com_materials_get_k_sca
		procedure,public:: get_cap => f_com_materials_get_cap_vec,f_com_materials_get_cap_sca
	end type ty_com_materials
	
	contains
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_com_materials. Initiallize count and allocate material.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_materials_allocateall(this,nmat)
	class(ty_com_materials),intent(inout)::this
	integer,intent(in)::nmat
	
	this%count = nmat
	if(.not.allocated(this%material)) allocate(this%material(nmat))
	
	end subroutine s_com_materials_allocateall
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_com_materials. deallocate material.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------
	
	subroutine s_com_materials_deallocateall(this)
	class(ty_com_materials),intent(inout)::this

	if(allocated(this%material)) deallocate(this%material)
	
	end subroutine s_com_materials_deallocateall
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns saturation at pore pressures h
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------
	function f_com_materials_get_s_vec(this,h,nmat)
	use com_mod_hyd,only:f_hyd_S_h_vec
	class(ty_com_materials),intent(inout)::this
	real(kind=dpd),intent(in)::h(:)
	integer,intent(in)::nmat(:)
	real(kind=dpd)::f_com_materials_get_s_vec(size(h))
	
	f_com_materials_get_s_vec = f_hyd_S_h_vec(h,this%material(nmat))
	
	end function f_com_materials_get_s_vec
	
	function f_com_materials_get_s_sca(this,h,nmat)
	use com_mod_hyd,only:f_hyd_S_h_sca
	class(ty_com_materials),intent(inout)::this
	real(kind=dpd),intent(in)::h
	integer,intent(in)::nmat
	real(kind=dpd)::f_com_materials_get_s_sca
	
	f_com_materials_get_s_sca = f_hyd_S_h_sca(h,this%material(nmat))
	
	end function f_com_materials_get_s_sca
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns water content given pore pressures h
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------
	function f_com_materials_get_th_vec(this,h,nmat)
	use com_mod_hyd,only:f_hyd_th_h_vec
	class(ty_com_materials),intent(inout)::this
	real(kind=dpd),intent(in)::h(:)
	integer,intent(in)::nmat(:)
	real(kind=dpd)::f_com_materials_get_th_vec(size(h))
	
	f_com_materials_get_th_vec = f_hyd_th_h_vec(h,this%material(nmat))
	
	end function f_com_materials_get_th_vec
	
	function f_com_materials_get_th_sca(this,h,nmat)
	use com_mod_hyd,only:f_hyd_th_h_sca
	class(ty_com_materials),intent(inout)::this
	real(kind=dpd),intent(in)::h
	integer,intent(in)::nmat
	real(kind=dpd)::f_com_materials_get_th_sca
	
	f_com_materials_get_th_sca = f_hyd_th_h_sca(h,this%material(nmat))
	
	end function f_com_materials_get_th_sca	

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns unsat permeability given pore pressures h
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------
	function f_com_materials_get_k_vec(this,h,nmat)
	use com_mod_hyd,only:f_hyd_k_h_vec
	class(ty_com_materials),intent(inout)::this
	real(kind=dpd),intent(in)::h(:)
	integer,intent(in)::nmat(:)
	real(kind=dpd)::f_com_materials_get_k_vec(size(h))
	
	f_com_materials_get_k_vec = f_hyd_k_h_vec(h,this%material(nmat))
	
	end function f_com_materials_get_k_vec
	
	function f_com_materials_get_k_sca(this,h,nmat)
	use com_mod_hyd,only:f_hyd_k_h_sca
	class(ty_com_materials),intent(inout)::this
	real(kind=dpd),intent(in)::h
	integer,intent(in)::nmat
	real(kind=dpd)::f_com_materials_get_k_sca
	
	f_com_materials_get_k_sca = f_hyd_k_h_sca(h,this%material(nmat))
	
	end function f_com_materials_get_k_sca
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns unsat permeability given pore pressures h
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------
	function f_com_materials_get_cap_vec(this,h,nmat)
	use com_mod_hyd,only:f_hyd_cap_h_vec
	class(ty_com_materials),intent(inout)::this
	real(kind=dpd),intent(in)::h(:)
	integer,intent(in)::nmat(:)
	real(kind=dpd)::f_com_materials_get_cap_vec(size(h))
	
	f_com_materials_get_cap_vec = f_hyd_cap_h_vec(h,this%material(nmat))
	
	end function f_com_materials_get_cap_vec
	
	function f_com_materials_get_cap_sca(this,h,nmat)
	use com_mod_hyd,only:f_hyd_cap_h_sca
	class(ty_com_materials),intent(inout)::this
	real(kind=dpd),intent(in)::h
	integer,intent(in)::nmat
	real(kind=dpd)::f_com_materials_get_cap_sca
	
	f_com_materials_get_cap_sca = f_hyd_cap_h_sca(h,this%material(nmat))
	
	end function f_com_materials_get_cap_sca	
	
	end module com_mod_ty_materials