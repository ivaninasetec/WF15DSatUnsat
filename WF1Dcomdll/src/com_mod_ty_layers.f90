	!********************************************************************************************************************
	! TITLE         : COM_MOD_TY_LAYERS: DERIVED TYPE THAT DEFINES COMMON PROPERTIES OF ALL LAYERS
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_ELEMENTS
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Derived type that defines common properties and methods of the layers
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module com_mod_ty_layers
	use com_mod_ty_material, only: ty_com_material

	implicit none
	include 'inc_precision.fi'


	private

	!******************************************************************************************************************
	! TY_COM_LAYERS
	! Derived type that defines properties and methods of the layers
	!******************************************************************************************************************

	type,public::ty_com_layers !< Class: layer (common definition)

		integer::count !< Number of layers

		real(kind=dps),allocatable::height(:) !< Height of the layer [L]
		real(kind=dps)::width !< Widht of the layer  [L]
		type(ty_com_material),allocatable:: material(:) !< Pointer to the material of the layer
		real(kind=dps),allocatable::htop(:)		!<total height in the top of the layer
		real(kind=dps),allocatable::hbottom(:) !<total height in the bottom of the layer
		integer,allocatable::id_elem_top(:)
		integer,allocatable::id_elem_bottom(:)

		real(kind=dpd)::zphr	!<Phreatic level (z coordinate in [L]) from the top left of the model.
		logical::topboundbyh !<If true, water piezometric head boundary condition defined in the top with flow1d_bound_h.txt
		logical::topboundbyq !<If true, water flow boundary contidion defined in the top with flow1d_bound_q.txt
		logical::bottombyphl !<If true, then the bottom  will be a dirichlet condition defined by phreatic level, otherwise bottom is newmann with q=0
		real(kind=dps),allocatable::slope(:)		!<Bank of the bottom of the layer

	contains

	procedure,public:: allocateall		=> s_layers_allocateall
	procedure,public:: deallocateall	=> s_layers_deallocateall

	procedure,public:: get_id_from_h	=> f_layers_get_id_from_h

	procedure,public:: get_h_from0_sca		=> f_layers_get_h_from0_sca
	procedure,public:: get_h_from0_vec		=> f_layers_get_h_from0_vec
	procedure,public:: get_h_from_h_to_h_sca	=> f_layers_get_h_from_h_to_h_sca
	procedure,public:: get_h_from_h_to_h_vec	=> f_layers_get_h_from_h_to_h_vec

	procedure,public:: get_inc_h_sca					=> f_layers_get_inc_h_sca
	procedure,public:: get_inc_h_vec					=> f_layers_get_inc_h_vec
	procedure,public:: get_inc_h_from0_sca		=> f_layers_get_inc_h_from0_sca
	procedure,public:: get_inc_h_from0_vec		=> f_layers_get_inc_h_from0_vec

	procedure,public:: get_water_inc_med	=> f_layers_get_water_inc_med_constant_simple

	end type ty_com_layers

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Get a vector with heights until a given h, begining from h=0 to given h
	!> @param[in] nlayers
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_h_from0_sca(this,h)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_h_from0_sca" :: f_layers_get_h_from0_sca
	!DEC$ endif

	class(ty_com_layers),intent(inout)::this
	real(kind=dpd),intent(in)::h
	real(kind=dpd)::f_layers_get_h_from0_sca(this%count+1)
	integer::n

	n=this%count

	f_layers_get_h_from0_sca(:)=max(0.0_dpd,min(this%htop(:),h)-this%hbottom(1))

	end function f_layers_get_h_from0_sca


	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Get a vector with heights until a given h, begining from h=0 to given h
	!> @param[in] nlayers
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_h_from0_vec(this,h)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_h_from0_vec" :: f_layers_get_h_from0_vec
	!DEC$ endif

	class(ty_com_layers),intent(inout)::this
	real(kind=dpd),intent(in)::h(:)
	real(kind=dpd)::f_layers_get_h_from0_vec(size(h),this%count+1)
	integer::n,i

	n=this%count

	do i=1,n+1
		f_layers_get_h_from0_vec(:,i)		=max(0.0_dpd,min(this%htop(i),h)-this%hbottom(1)) !Difference from 0
	end do

	end function f_layers_get_h_from0_vec

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Initalize count and allocate layer(:)
	!> @param[in] nlayers
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_layers_allocateall(this,nlayers)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_layers_allocateall" :: s_layers_allocateall
	!DEC$ endif

	class(ty_com_layers),intent(inout)::this
	integer,intent(in)::nlayers

	this%count = nlayers
	if(.not.allocated(this%height)) allocate(this%height(nlayers))
	if(.not.allocated(this%material)) allocate(this%material(nlayers))
	if(.not.allocated(this%htop)) allocate(this%htop(nlayers))
	if(.not.allocated(this%hbottom)) allocate(this%hbottom(nlayers))
	if(.not.allocated(this%id_elem_top)) allocate(this%id_elem_top(nlayers))
	if(.not.allocated(this%id_elem_bottom)) allocate(this%id_elem_bottom(nlayers))

	if(.not.allocated(this%slope)) allocate(this%slope(nlayers+1))
	end subroutine s_layers_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Deallocate layer(:)
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_layers_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_layers_deallocateall" :: s_layers_deallocateall
	!DEC$ endif

	class(ty_com_layers),intent(inout)::this

	if(allocated(this%height)) deallocate(this%height)
	if(allocated(this%material)) deallocate(this%material)
	if(allocated(this%htop)) deallocate(this%htop)
	if(allocated(this%hbottom)) deallocate(this%hbottom)
	if(allocated(this%id_elem_top)) deallocate(this%id_elem_top)
	if(allocated(this%id_elem_bottom)) deallocate(this%id_elem_bottom)
	if(allocated(this%slope)) deallocate(this%slope)

	end subroutine s_layers_deallocateall

	!---------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns the id of the layer given the height h.
	!> @param[in] h
	!---------------------------------------------------------------------------

	function f_layers_get_id_from_h(this,h)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_id_from_h" :: f_layers_get_id_from_h
	!DEC$ endif

	class(ty_com_layers),intent(in)::this
	real(kind=dpd),intent(in)::h(:)
	integer::f_layers_get_id_from_h(size(h))
	integer::i,findtemp(1)
	logical::mask(size(this%hbottom))


	DO i=1, size(h)
		mask = this%hbottom<=h(i)
		findtemp = FINDLOC(mask,.true.,BACK=.true.)
		f_layers_get_id_from_h(i) = max(findtemp(1),1)
	end do

	end function f_layers_get_id_from_h

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> This function returns the increment of height from the layer h0 to the top of the layer or to h1
	!> (the minimum of both) in each of the layers of the collection. The output is a double matrix
	!>	with dimension nlayers+1 to account for water on top of upper layer.
	!> @param[in] h0
	!> @param[in] h1
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_inc_h_sca(this,h0,h1) 
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_inc_h_sca" :: f_layers_get_inc_h_sca
	!DEC$ endif

	class(ty_com_layers),intent(in)::this
	real(kind=dpd),intent(in)::h0 !< Height of first  point from the bottom of first layer
	real(kind=dpd),intent(in)::h1 !< Height of second point from the bottom of first layer
	real(kind=dpd)::f_layers_get_inc_h_sca(this%count+1)

	if (h1>h0) then
		f_layers_get_inc_h_sca(1:this%count) = max(0.0_dpd,min(this%htop,h1)-max(this%hbottom,h0))
		f_layers_get_inc_h_sca(this%count+1) = max(0.0_dpd,h1-max(this%htop(this%count),h0))
	else
		f_layers_get_inc_h_sca(1:this%count) = max(0.0_dpd,min(this%htop,h0)-max(this%hbottom,h1))
		f_layers_get_inc_h_sca(this%count+1) = max(0.0_dpd,h0-max(this%htop(this%count),h1))
	end if

	end function f_layers_get_inc_h_sca

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> This function returns the increment of height from the layer h0 to the top of the layer or to h1
	!> (the minimum of both) in each of the layers of the collection. The output is a double matrix
	!>	with dimension nlayers+1 to account for water on top of upper layer.
	!> @param[in] h0
	!> @param[in] h1
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_inc_h_vec(this,h0,h1)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_inc_h_vec" :: f_layers_get_inc_h_vec
	!DEC$ endif

	class(ty_com_layers),intent(in)::this
	real(kind=dpd),intent(in)::h0(:) !< Height of first  point from the bottom of first layer
	real(kind=dpd),intent(in)::h1(:) !< Height of second point from the bottom of first layer
	real(kind=dpd)::f_layers_get_inc_h_vec(size(h0),this%count+1)

	integer::l

	do l=1,this%count
		where (h1>h0)
			f_layers_get_inc_h_vec(:,l) = max(0.0_dpd,min(this%htop(l),h1(:))-max(this%hbottom(l),h0(:)))
		else where
			f_layers_get_inc_h_vec(:,l) = max(0.0_dpd,min(this%htop(l),h0(:))-max(this%hbottom(l),h1(:)))
		end where
	end do

	!For final layer
	where (h1>h0)
		f_layers_get_inc_h_vec(:,this%count+1) = max(0.0_dpd,h1(:)-max(this%htop(this%count),h0(:)))
	else where
		f_layers_get_inc_h_vec(:,this%count+1) = max(0.0_dpd,h0(:)-max(this%htop(this%count),h1(:)))
	end where

	end function f_layers_get_inc_h_vec

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> This function returns the increment of height from the layer h0 to the top of the layer or to h1
	!> (the minimum of both) in each of the layers of the collection. The output is a double matrix
	!>	with dimension nlayers+1 to account for water on top of upper layer.
	!> @param[in] h0
	!> @param[in] h1
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_inc_h_from0_sca(this,h1)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_inc_h_from0_sca" :: f_layers_get_inc_h_from0_sca
	!DEC$ endif

	class(ty_com_layers),intent(in)::this
	real(kind=dpd),intent(in)::h1 !< Height of second point from the bottom of first layer
	real(kind=dpd)::f_layers_get_inc_h_from0_sca(this%count+1)

	if (h1>0.0_dpd) then
		f_layers_get_inc_h_from0_sca(1:this%count) = max(0.0_dpd,min(this%htop,h1)-max(this%hbottom,0.0_dpd))
		f_layers_get_inc_h_from0_sca(this%count+1) = max(0.0_dpd,h1-max(this%htop(this%count),0.0_dpd))
	else
		f_layers_get_inc_h_from0_sca = 0.0_dpd
		f_layers_get_inc_h_from0_sca(1) = h1
	end if

	end function f_layers_get_inc_h_from0_sca

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> This function returns the increment of height from the layer h0 to the top of the layer or to h1
	!> (the minimum of both) in each of the layers of the collection. The output is a double matrix
	!>	with dimension nlayers+1 to account for water on top of upper layer.
	!> @param[in] h0
	!> @param[in] h1
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_inc_h_from0_vec(this,h1)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_inc_h_from0_vec" :: f_layers_get_inc_h_from0_vec
	!DEC$ endif

	class(ty_com_layers),intent(in)::this
	real(kind=dpd),intent(in)::h1(:) !< Height of first  point from the bottom of first layer
	real(kind=dpd)::f_layers_get_inc_h_from0_vec(size(h1),this%count+1)

	integer::l

	do l=1,this%count
		if (l==1) then
			where (h1>0.0_dpd)
				f_layers_get_inc_h_from0_vec(:,l) = max(0.0_dpd,min(this%htop(l),h1(:))-max(this%hbottom(l),0.0_dpd))
			else where
				f_layers_get_inc_h_from0_vec(:,l) = h1(:)
			end where
		else
			where (h1>0.0_dpd)
				f_layers_get_inc_h_from0_vec(:,l) = max(0.0_dpd,min(this%htop(l),h1(:))-max(this%hbottom(l),0.0_dpd))
			else where
				f_layers_get_inc_h_from0_vec(:,l) = 0.0_dpd
			end where
		end if
	end do

	!For final layer
	where (h1>0.0_dpd)
		f_layers_get_inc_h_from0_vec(:,this%count+1) = max(0.0_dpd,h1(:)-max(this%htop(this%count),0.0_dpd))
	else where
		f_layers_get_inc_h_from0_vec(:,l) = 0.0_dpd
	end where

	end function f_layers_get_inc_h_from0_vec


	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns height in each layer from h0, to h1. So in the second layer is added the height from 0, from the first layer.
	!> @param[in] h0
	!> @param[in] h1
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_h_from_h_to_h_sca(this,h0,h1)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_h_from_h_to_h_sca" :: f_layers_get_h_from_h_to_h_sca
	!DEC$ endif

	class(ty_com_layers),intent(in)::this
	real(kind=dpd),intent(in)::h0 !< Height of first  point from the bottom of first layer
	real(kind=dpd),intent(in)::h1 !< Height of second point from the bottom of first layer
	real(kind=dpd)::f_layers_get_h_from_h_to_h_sca(this%count+1)

	integer::nlayers,i

	if (h1>h0) then
		f_layers_get_h_from_h_to_h_sca(1:this%count) = max(0.0_dpd,min(this%htop,h1)-max(this%hbottom,h0))
		f_layers_get_h_from_h_to_h_sca(this%count+1) = max(0.0_dpd,h1-max(this%htop(this%count),h0))
	else
		f_layers_get_h_from_h_to_h_sca(1:this%count) = - max(0.0_dpd,min(this%htop,h0)-max(this%hbottom,h1))
		f_layers_get_h_from_h_to_h_sca(this%count+1) = - max(0.0_dpd,h0-max(this%htop(this%count),h1))
	end if

	end function f_layers_get_h_from_h_to_h_sca

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns height in each layer from h0, to h1. So in the second layer is added the height from 0, from the first layer.
	!> @param[in] h0
	!> @param[in] h1
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_h_from_h_to_h_vec(this,h0,h1)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_h_from_h_to_h_vec" :: f_layers_get_h_from_h_to_h_vec
	!DEC$ endif

	class(ty_com_layers),intent(in)::this
	real(kind=dpd),intent(in)::h0(:) !< Height of first  point from the bottom of first layer
	real(kind=dpd),intent(in)::h1(:) !< Height of second point from the bottom of first layer
	real(kind=dpd)::f_layers_get_h_from_h_to_h_vec(size(h0),this%count+1)

	integer::l

	do l=1,this%count
		where (h1>h0)
			f_layers_get_h_from_h_to_h_vec(:,l) = max(0.0_dpd,min(this%htop(l),h1(:))-max(this%hbottom(l),h0(:)))

		else where
			f_layers_get_h_from_h_to_h_vec(:,l) = - max(0.0_dpd,min(this%htop(l),h0(:))-max(this%hbottom(l),h1(:)))
		end where
	end do

	!Over the top layer:
	where (h1>h0)
		f_layers_get_h_from_h_to_h_vec(:,this%count+1) = max(0.0_dpd,h1(:)-max(this%htop(this%count),h0(:)))
	else where
		f_layers_get_h_from_h_to_h_vec(:,this%count+1) = - max(0.0_dpd,h0(:)-max(this%htop(this%count),h1(:)))
	end where

	end function f_layers_get_h_from_h_to_h_vec

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> This function returns the ammount of water needed to fill from h0 (or h1) to h1 (or h0) in the
	!> unsaturated multilayered system considering van-genuchten unsaturated materials
	!> In this case we consider the unsaturated material over water table completely dry (not real) and return the mean
	! thsat-thres between h0 and h1.
	!> @param[in] h0
	!> @param[in] h1
	!---------------------------------------------------------------------------------------------------------------------

	function f_layers_get_water_inc_med_constant_simple(this,h0,h1)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_layers_get_water_inc_med_constant_simple" :: f_layers_get_water_inc_med_constant_simple
	!DEC$ endif

	integer,parameter::NDIVISIONS=10
	real(kind=dpd), parameter::MIN_H0_TO_H1 = 1.0E-5_dpd;
	class(ty_com_layers),intent(in)::this
	real(kind=dpd),intent(in)::h0(:),h1(size(h0))
	real(kind=dpd)::f_layers_get_water_inc_med_constant_simple(size(h0))

	real(kind=dpd)::thsatrestemp(this%count+1)

	thsatrestemp(1:this%count) = this%material%thsat-this%material%thres
	thsatrestemp(this%count+1) = 1.0_dpd

	where (abs(h0-h1)<MIN_H0_TO_H1)
		f_layers_get_water_inc_med_constant_simple = this%material(this%get_id_from_h(h1))%thsat-this%material(this%get_id_from_h(h1))%thres
	else where
		!tex:
		!\[\overline {({\theta _{sat,l}} - {\theta _{res,l}})}  = \frac{{\sum\limits_l {({\theta _{sat,l}} - {\theta _{res,l}})\Delta {h_l}} }}{{\Delta h}}\]
		f_layers_get_water_inc_med_constant_simple = abs(matmul(this%get_inc_h_vec(h0,h1),thsatrestemp)/(h1-h0))
	end where

	end function f_layers_get_water_inc_med_constant_simple


	end module com_mod_ty_layers