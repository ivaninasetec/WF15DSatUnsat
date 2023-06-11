	!********************************************************************************************************************
	!        CLASS THAT INCLUDE THE WHOLE MODEL
  !********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : mod_sat_ty_nodes
	! URL           : ...
	! AFFILIATION   : ...
	! DATE          : ...
	! REVISION      : ... V 0.0
	! LICENSE				: This software is copyrighted 2019(C)
	!> @author
	!> Iván Campos-Guereta Díez
  !  MSc Civil Engineering by Polytechnic University of Madrid                                                     *
  !  PhD Student by University of Nottingham                                                                       *
  !  eMBA by International Institute San Telmo in Seville                                                          *
  !  ivan.camposguereta@nottingham.ac.uk  
	! DESCRIPTION:
	!> Class for the collection of nodes-classes in the saturated model            
  !********************************************************************************************************************
	
	module model_mod_ty_constraints
		use unsat_mod_ty_model,				only: ty_unsat_model
		use sat_mod_ty_model,					only: ty_sat_model
	
  implicit none
  include 'inc_precision.fi'

  private
		!real(dpd),parameter::MIN_NREL=1.0E-10,MAX_NREL=1.0
	
	  type,public::ty_model_constraints
			type(ty_unsat_model),pointer::unsat(:)
			type(ty_sat_model),pointer::sat(:)
			integer(KIND=4)::idprint=0							!<printed line in OUTPUTS\model_constraints.cs
			!integer,allocatable::idnode_h(:,:)
			!integer,allocatable::idnode_v(:,:)
			!integer,allocatable::idnode_hsat_v(:,:)
			!
			!real(kind=dpd),allocatable::hsat_h(:,:)
			!real(kind=dpd),allocatable::hsat_v(:,:)
			!real(kind=dpd),allocatable::hsattemp_v(:,:)
			!real(kind=dpd),allocatable::hsatold_v(:,:)
			!real(kind=dpd),allocatable::dqhordx_h(:,:)
			!real(kind=dpd),allocatable::qver_v(:,:)
			!real(kind=dpd),allocatable::nrel(:,:)
			!
			!type(ty_unsat_model)	,pointer::unsat(:)
			!type(ty_sat_model)		,pointer::sat(:)
		contains
		procedure,public:: construct		=> s_model_constraints_construct
		procedure,public:: get_hsatv_mat=> f_model_constraints_get_hsat_v
		procedure,public:: get_hsatv		=> f_model_constraints_get_hsat_v_i_j
		procedure,public:: get_sum_hsatv_on_is	=> f_model_constraints_get_sum_hsatv_on_is
		procedure,public:: get_hsath_mat=> f_model_constraints_get_hsat_h
		procedure,public:: get_hsath		=> f_model_constraints_get_hsat_h_i_j
		procedure,public:: get_hsath_mean_mat=> f_model_constraints_get_hsat_h_mean
		procedure,public:: get_hsath_mean		=> f_model_constraints_get_hsat_h_mean_i_j
		procedure,public:: get_hnewh_mean_mat=> f_model_constraints_get_hnew_h_mean
		procedure,public:: get_hnewh_mean		=> f_model_constraints_get_hnew_h_mean_i_j
		procedure,public:: get_idnodeh	=> f_model_constraints_get_idnode_h_i_j
		procedure,public:: get_idnodev	=> f_model_constraints_get_idnode_v_i_j
		
		procedure,public:: get_qverv_mat=> f_model_constraints_get_qver_v
		procedure,public:: get_qverv		=> f_model_constraints_get_qver_v_i_j
		procedure,public:: get_qnewmannv_mat=> f_model_constraints_get_qnewmann_v
		procedure,public:: get_qnewmannv	=> f_model_constraints_get_qnewmann_v_i_j
		procedure,public:: get_dqhordxh_mat	=> f_model_constraints_get_dqhordx_h
		procedure,public:: get_dqhordxh	=> f_model_constraints_get_dqhordx_h_i_j
		procedure,public:: get_dqhordxh_mean_mat	=> f_model_constraints_get_dqhordx_mean_h
		procedure,public:: get_dqhordxh_mean	=> f_model_constraints_get_dqhordx_mean_h_i_j
		procedure,public:: update_nrel_to_fit_inc_hunsat	=> s_model_constraints_calc_nrel_to_fit_inc_hunsat
		procedure,public:: update_nrel_to_fit_inc_hsat		=> s_model_constraints_calc_nrel_to_fit_inc_hsat
		
		procedure,public:: set_dirichlet_from_hsat_h	=> s_model_constraint_set_dirichlet_from_hsath_mean
		procedure,public:: set_newmann_from_sat_h	=> s_model_constraint_set_newmann_from_sath_with_relaxation
		
		procedure,public:: get_error_in_hsat => f_get_error_in_hsat
		procedure,public:: get_error_in_q => f_get_error_in_q
		end type ty_model_constraints
		
	contains
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	subroutine s_model_constraints_construct(this,unsat,sat)
		use unsat_mod_ty_model,				only: ty_unsat_model
		use sat_mod_ty_model,					only: ty_sat_model
	
	class(ty_model_constraints),intent(inout)::this
	type(ty_unsat_model),intent(in),target::unsat(:)
	type(ty_sat_model),intent(in),target::sat(:)
	
	integer::ih,iv
	
	this%unsat	=> unsat
	this%sat		=> sat
	
	do ih=1,size(this%sat)
		do iv=1,size(this%unsat)
			this%sat(ih)%constraints%qent(iv)%p				=> this%unsat(iv)%constraints%qver(ih)
			this%sat(ih)%constraints%qsat(iv)%p				=> this%unsat(iv)%constraints%qsat(ih)
			this%sat(ih)%constraints%qinf(iv)%p				=> this%unsat(iv)%constraints%qinf(ih)
			this%sat(ih)%constraints%nrel(iv)%p				=> this%unsat(iv)%constraints%nrel(ih)
			
			this%unsat(iv)%constraints%dqhordx(ih)%p	=> this%sat(ih)%constraints%dqhordx(iv)%p
		end do
	end do
	
	end subroutine s_model_constraints_construct
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_hsat_v(this)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_hsat_v(size(this%sat),size(this%unsat))
	integer::iu,is
	do is=1,size(this%sat)
		do iu=1,size(this%unsat)
			f_model_constraints_get_hsat_v(is,iu)=this%unsat(iu)%constraints%hsat(is)
		end do
	end do
	end function f_model_constraints_get_hsat_v
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_sum_hsatv_on_is(this,is)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_sum_hsatv_on_is
	integer,intent(in)::is
	integer::iu
	
	f_model_constraints_get_sum_hsatv_on_is = 0.0_dpd
		do iu=1,size(this%unsat)
			f_model_constraints_get_sum_hsatv_on_is = f_model_constraints_get_sum_hsatv_on_is+this%unsat(iu)%constraints%hsat(is)
		end do
	end function f_model_constraints_get_sum_hsatv_on_is
	
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_hsat_v_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_hsat_v_i_j

			f_model_constraints_get_hsat_v_i_j = this%unsat(iu)%constraints%hsat(is)

	end function f_model_constraints_get_hsat_v_i_j

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_hsat_h_mean(this)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_hsat_h_mean(size(this%sat),size(this%unsat))
	integer::iu,is
	do is=1,size(this%sat)
		do iu=1,size(this%unsat)
			f_model_constraints_get_hsat_h_mean(is,iu)=this%sat(is)%constraints%hsat_mean(iu)
		end do
	end do
	end function f_model_constraints_get_hsat_h_mean
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_hsat_h_mean_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_hsat_h_mean_i_j

			f_model_constraints_get_hsat_h_mean_i_j = this%sat(is)%constraints%hsat_mean(iu)

	end function f_model_constraints_get_hsat_h_mean_i_j	
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_hnew_h_mean(this)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_hnew_h_mean(size(this%sat),size(this%unsat))
	integer::iu,is
	do is=1,size(this%sat)
		do iu=1,size(this%unsat)
			f_model_constraints_get_hnew_h_mean(is,iu)=this%sat(is)%constraints%hsat_mean(iu)-this%sat(is)%constraints%inchnew_mean(iu)
		end do
	end do
	end function f_model_constraints_get_hnew_h_mean
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_hnew_h_mean_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_hnew_h_mean_i_j

			f_model_constraints_get_hnew_h_mean_i_j = this%sat(is)%constraints%hsat_mean(iu)-this%sat(is)%constraints%inchnew_mean(iu)

	end function f_model_constraints_get_hnew_h_mean_i_j	
	
	
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_hsat_h(this)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_hsat_h(size(this%sat),size(this%unsat))
	integer::iu,is
	do is=1,size(this%sat)
		do iu=1,size(this%unsat)
			f_model_constraints_get_hsat_h(is,iu)=this%sat(is)%constraints%hsat(iu)%p
		end do
	end do
	end function f_model_constraints_get_hsat_h
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_hsat_h_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_hsat_h_i_j

			f_model_constraints_get_hsat_h_i_j = this%sat(is)%constraints%hsat(iu)%p

	end function f_model_constraints_get_hsat_h_i_j

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_idnode_h_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_idnode_h_i_j

			f_model_constraints_get_idnode_h_i_j = this%sat(is)%constraints%idnode(iu)%p

	end function f_model_constraints_get_idnode_h_i_j
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_idnode_v_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_idnode_v_i_j
 
			!f_model_constraints_get_idnode_v_i_j = this%unsat(iu)%constraints%idnode_layer_bottom(is)%p
      f_model_constraints_get_idnode_v_i_j = this%unsat(iu)%constraints%idnode_layer_bottom(is)%p
 
	end function f_model_constraints_get_idnode_v_i_j	
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_qver_v(this)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_qver_v(size(this%sat),size(this%unsat))
	integer::iu,is
	do is=1,size(this%sat)
		do iu=1,size(this%unsat)
			f_model_constraints_get_qver_v(is,iu)=this%unsat(iu)%constraints%qver(is)
		end do
	end do
	end function f_model_constraints_get_qver_v
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_qver_v_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_qver_v_i_j

			f_model_constraints_get_qver_v_i_j = this%unsat(iu)%constraints%qver(is)

	end function f_model_constraints_get_qver_v_i_j
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_qnewmann_v(this)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_qnewmann_v(size(this%sat),size(this%unsat))
	integer::iu,is
	do is=1,size(this%sat)
		do iu=1,size(this%unsat)
			f_model_constraints_get_qnewmann_v(is,iu)=this%unsat(iu)%constraints%qnewmann(is)%p
		end do
	end do
	end function f_model_constraints_get_qnewmann_v
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_qnewmann_v_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_qnewmann_v_i_j

			f_model_constraints_get_qnewmann_v_i_j = this%unsat(iu)%constraints%qnewmann(is)%p

	end function f_model_constraints_get_qnewmann_v_i_j
	
		!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_dqhordx_h(this)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_dqhordx_h(size(this%sat),size(this%unsat))
	integer::iu,is
	do is=1,size(this%sat)
		do iu=1,size(this%unsat)
			f_model_constraints_get_dqhordx_h(is,iu)=this%sat(is)%constraints%dqhordx(iu)%p
		end do
	end do
	end function f_model_constraints_get_dqhordx_h
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_dqhordx_h_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_dqhordx_h_i_j

			f_model_constraints_get_dqhordx_h_i_j = this%sat(is)%constraints%dqhordx(iu)%p

	end function f_model_constraints_get_dqhordx_h_i_j
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_dqhordx_mean_h(this)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_model_constraints_get_dqhordx_mean_h(size(this%sat),size(this%unsat))
	integer::iu,is
	do is=1,size(this%sat)
		do iu=1,size(this%unsat)
			f_model_constraints_get_dqhordx_mean_h(is,iu)=this%sat(is)%constraints%dqhordx_mean(iu)
		end do
	end do
	end function f_model_constraints_get_dqhordx_mean_h
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------
	function f_model_constraints_get_dqhordx_mean_h_i_j(this,is,iu)
	class(ty_model_constraints),intent(inout)::this
	integer,intent(in)::is,iu
	real(kind=dpd)::f_model_constraints_get_dqhordx_mean_h_i_j

			f_model_constraints_get_dqhordx_mean_h_i_j = this%sat(is)%constraints%dqhordx_mean(iu)

	end function f_model_constraints_get_dqhordx_mean_h_i_j

	
	!!---------------------------------------------------------------------------------------------------------------------
	!!> @author Iván Campos-Guereta Díez
	!!> @brief
	!!> Allocate all arrays
	!!---------------------------------------------------------------------------------------------------------------------
	!function f_model_constraints_get_max_hsat_v(this)
	!class(ty_model_constraints),intent(inout)::this
	!real(kind=dpd)::f_model_constraints_get_max_hsat_v
	!integer::iu,is
	!
	!f_model_constraints_get_max_hsat_v = 0.0_dpd
	!do is=1,size(this%sat)
	!	do iu=1,size(this%unsat)
	!		f_model_constraints_get_max_hsat_v = max(f_model_constraints_get_max_hsat_v,this%unsat(iu)%constraints%hsat(is))
	!	end do
	!end do
	!end function f_model_constraints_get_max_hsat_v	
	!
	
	subroutine s_model_constraints_calc_nrel_to_fit_inc_hunsat(this,dt)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd),intent(in)::dt
	real(kind=dpd)::hsat_u,hsatold_u,qver,qnewmann,sat_waterinc_med(1),nrel_temp(1)
	
	
	integer::is,iu,nsat,nunsat
	
	nsat		= size(this%sat)
	nunsat	= size(this%unsat)
	
	!This it to get the nrel needed to increase in sat model from hsat_u to hsatold_u with given qver entering-leaving.
	!dt·(qent-qsal-qnewmann)=nrel·(thsat-threl)·(hsat_u-hsatold_u)
	do is=1,nsat
		do iu=1,nunsat
				hsat_u= this%unsat(iu)%constraints%hsat(is)
				hsatold_u= this%unsat(iu)%constraints%hsatold(is)
				qver = this%unsat(iu)%constraints%qver(is)
				!qnewmann = this%unsat(iu)%constraints%qnewmann(is)%p
				qnewmann = this%unsat(iu)%constraints%qnewmann_all(is)
				sat_waterinc_med = this%sat(is)%layers%get_water_inc_med((/hsatold_u/),(/hsat_u/)) !This returns the mean of (thsat-thres) between hsatold and hsat
				if ((abs(sat_waterinc_med(1)*(hsat_u-hsatold_u))<1.0E-20_dpd).or.(abs(qver-qnewmann)<1.0E-20_dpd)) then
					!Maintain previous nrel.
					!this%unsat(iu)%constraints%nrel(is)=1.0_dpd !CHECK
				else
          !this%unsat(iu)%constraints%nrel(is) = 0.01_dpd
					this%unsat(iu)%constraints%nrel(is) = max(this%unsat(iu)%parameters%nrelmin,min(this%unsat(iu)%parameters%nrelmax,abs((qver-qnewmann)*dt/(sat_waterinc_med(1)*(hsat_u-hsatold_u))))) 
				end if
		end do 
	end do
	
!!This is getting nrel automatically by using the hypergeometric function:	
!		do is=1,nsat
!		do iu=1,nunsat
!			hsat_u= this%unsat(iu)%constraints%hsat(is)
!			hsatold_u= this%unsat(iu)%constraints%hsatold(is)
!			if(abs(hsat_u-hsatold_u)>0.0_dpd) then
!			!nrel_temp = this%sat(is)%layers%get_nrel_from_hsat_to_h((/hsatold_u/),(/hsat_u/))
!			nrel_temp = 0.1_dpd !Forcing to 0.1
!			this%unsat(iu)%constraints%nrel(is) = nrel_temp(1)
!			end if
!			end do 
!	end do
	
	end subroutine s_model_constraints_calc_nrel_to_fit_inc_hunsat
	

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Allocate all arrays
	!---------------------------------------------------------------------------------------------------------------------

	
	subroutine s_model_constraints_calc_nrel_to_fit_inc_hsat(this,epshsathunsat,krelax)
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd),intent(in),optional::krelax
	real(kind=dpd),intent(inout)::epshsathunsat
	real(kind=dpd)::hsat_s,hsatold_s,hsat_u,hsatold_u,qver,qnewmann,sat_waterinc_med_u(1),sat_waterinc_med_s(1),nrel_temp(1)
	
	integer::is,iu,nsat,nunsat
	real(kind=dpd)::krelx
		
	if(present(krelax)) then
		krelx = krelax
	else
		krelx = 1.0_dpd
	end if
	
	nsat		= size(this%sat)
	nunsat	= size(this%unsat)
	
	!This it to get the nrel needed to increase in sat model from hsat_u to hsatold_u with given qver entering-leaving.
	!dt·(qent-qsal-qnewmann)=nrel·(thsat-threl)·(hsat_u-hsatold_u)
	epshsathunsat = 0.0_dpd
	do is=1,nsat
		do iu=1,nunsat
				hsat_u= this%unsat(iu)%constraints%hsat(is)
				hsatold_u= this%unsat(iu)%constraints%hsatold(is)
				hsat_s= this%sat(is)%constraints%hsat(iu)%p
				hsatold_s= this%sat(is)%constraints%hsatold(iu)%p
				epshsathunsat = max(epshsathunsat,abs(hsat_u-hsat_s))
				sat_waterinc_med_u = this%sat(is)%layers%get_water_inc_med((/hsatold_s/),(/hsat_u/)) !This returns the mean of (thsat-thres) between hsatold and hsat
				sat_waterinc_med_s = this%sat(is)%layers%get_water_inc_med((/hsatold_s/),(/hsat_s/)) !This returns the mean of (thsat-thres) between hsatold and hsat

				if(sat_waterinc_med_u(1)*(hsat_u-hsatold_s)>1.0E-20_dpd) then
					this%unsat(iu)%constraints%nrel(is) = max(this%unsat(iu)%parameters%nrelmin,min(this%unsat(iu)%parameters%nrelmax,(1-krelax)*this%unsat(iu)%constraints%nrel(is)+krelax*this%unsat(iu)%constraints%nrel(is)*sat_waterinc_med_s(1)*(hsat_s-hsatold_s)/(sat_waterinc_med_u(1)*(hsat_u-hsatold_s))))
				end if
		end do 
	end do
	
	end subroutine s_model_constraints_calc_nrel_to_fit_inc_hsat
	
	
	
	
	
	
	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_model_constraint_set_dirichlet_from_hsat_h(this)
	use com_mod_ty_nodes,only:ty_com_nodes
	!idnode to idnode_hsat
	
	class(ty_model_constraints),intent(inout)::this
	
	integer::is,iu,nsat,nunsat,idnode_u
	
	nsat = size(this%sat)
	nunsat = size(this%unsat)
	
	do is=1,nsat
	do iu=1,nunsat
		!idnode_u = this%unsat(iu)%constraints%idnode_layer_bottom(is)%p
    idnode_u = this%unsat(iu)%constraints%idnode_layer_bottom(is)%p
	if ((this%get_hsath(is,iu)>1E-10_dpd)) then 
		this%unsat(iu)%calc%nodes%hdirichlet	(idnode_u)	= this%get_hsath(is,iu) !******We put as dirichlet hsath, not hsath_mean
		this%unsat(iu)%calc%nodes%isdirichlet	(idnode_u)	= .true.
	else
		this%unsat(iu)%calc%nodes%hdirichlet	(idnode_u)	= 0.0_dpd
		this%unsat(iu)%calc%nodes%isdirichlet	(idnode_u)	= .false.
	end if
	end do
	end do
	
	end subroutine s_model_constraint_set_dirichlet_from_hsat_h

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_model_constraint_set_dirichlet_from_hsath_mean(this,krelax)
	use com_mod_ty_nodes,only:ty_com_nodes
	
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd),optional::krelax	
	integer::is,iu,nsat,nunsat,idnode_u
	real(kind=dpd)::krelx
		
	nsat = size(this%sat)
	nunsat = size(this%unsat)
	
	if(present(krelax)) then
		krelx = krelax
	else
		krelx = 1.0_dpd
	end if
	
	do is=1,nsat
	do iu=1,nunsat
		!idnode_u = this%unsat(iu)%constraints%idnode_layer_bottom(is)%p
    idnode_u = this%unsat(iu)%constraints%idnode_layer_bottom(is)%p
		this%unsat(iu)%calc%nodes%qnewmann	(idnode_u)		= 0.0_dpd
		this%unsat(iu)%calc%nodes%isnewmann	(idnode_u)		= .false.
	!if (this%get_hsath_mean(is,iu)>1.0E-3_dpd) then !Only acttivate dirichlet on watertables over 1E-3_dpd
	if ((this%get_hsath_mean(is,iu)>1E-10_dpd)) then 
		this%unsat(iu)%calc%nodes%hdirichlet	(idnode_u)	= (1.0_dpd-krelx)*this%get_hsatv(is,iu)+krelx*this%get_hsath_mean(is,iu) !******We put as dirichlet the mean, not hsat
		!this%unsat(iu)%calc%nodes%hdirichlet	(idnode_u)	= (1.0_dpd-krelx)*this%get_hsath(is,iu)+krelx*this%get_hnewh_mean(is,iu) !CHECK:******We put as dirichlet the mean, not hsat
		this%unsat(iu)%calc%nodes%isdirichlet	(idnode_u)	= .true.
	else
		this%unsat(iu)%calc%nodes%hdirichlet	(idnode_u)	= 0.0_dpd
		this%unsat(iu)%calc%nodes%isdirichlet	(idnode_u)	= .false.
	end if
	end do
	end do
	
	end subroutine s_model_constraint_set_dirichlet_from_hsath_mean	
	
	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	subroutine s_model_constraint_set_newmann_from_sath_with_relaxation(this,krelax)
	use com_mod_ty_nodes,only:ty_com_nodes
	
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd),optional::krelax
	integer::is,iu,nsat,nunsat,idnode_u
	real(kind=dpd)::krelx
	
	nsat = size(this%sat)
	nunsat = size(this%unsat)
	
	if(present(krelax)) then
		krelx = krelax
	else
		krelx = 1.0_dpd
	end if
	
	do is=1,nsat
	do iu=1,nunsat
    idnode_u = this%unsat(iu)%constraints%idnode_layer_bottom(is)%p
		this%unsat(iu)%calc%nodes%hdirichlet	(idnode_u)	= 0.0_dpd
		this%unsat(iu)%calc%nodes%isdirichlet	(idnode_u)	= .false.		
	if ((this%get_hsath_mean(is,iu)>1E-10_dpd)) then 
		this%unsat(iu)%calc%nodes%qnewmann	(idnode_u)		= (1.0_dpd-krelx)*this%get_qnewmannv(is,iu)+krelx*(this%get_dqhordxh(is,iu))
		!this%unsat(iu)%calc%nodes%qnewmann	(idnode_u)		= krelx*(-this%get_dqhordxh_mean(is,iu))+(1.0_dpd-krelx)*(-this%get_dqhordxh(is,iu)) !CHECK.Changed (sign?)
		this%unsat(iu)%calc%nodes%isnewmann	(idnode_u)		= .true.		
	else
		this%unsat(iu)%calc%nodes%qnewmann	(idnode_u)		= 0.0_dpd
		this%unsat(iu)%calc%nodes%isnewmann	(idnode_u)		= .false.
	end if
	end do
	end do
	
	end subroutine s_model_constraint_set_newmann_from_sath_with_relaxation		
	

	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	function f_get_error_in_hsat(this)
	use com_mod_ty_nodes,only:ty_com_nodes
	
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_get_error_in_hsat
	
	integer::is,iu,nsat,nunsat
	
	nsat = size(this%sat)
	nunsat = size(this%unsat)
	
	f_get_error_in_hsat = 0.0_dpd
	
	do is=1,nsat
	do iu=1,nunsat
		f_get_error_in_hsat = max(f_get_error_in_hsat,abs(this%get_hsatv(is,iu)-this%get_hsath(is,iu)))
	end do
	end do
	
	end function f_get_error_in_hsat		
	
	!******************************************************************************************************************
	! Sub: s_layer_get_hsat(layer,option):	Update hsat in the layer (saturation height from nodes with h>0)
	!******************************************************************************************************************

	function f_get_error_in_q(this)
	use com_mod_ty_nodes,only:ty_com_nodes
	
	class(ty_model_constraints),intent(inout)::this
	real(kind=dpd)::f_get_error_in_q
	
	integer::is,iu,nsat,nunsat
	
	nsat = size(this%sat)
	nunsat = size(this%unsat)
	
	f_get_error_in_q = 0.0_dpd
	
	do is=1,nsat
	do iu=1,nunsat
		f_get_error_in_q = max(f_get_error_in_q,abs(this%get_qnewmannv(is,iu)-(this%get_dqhordxh(is,iu))))
	end do
	end do
	
	end function f_get_error_in_q		
	
	end module model_mod_ty_constraints