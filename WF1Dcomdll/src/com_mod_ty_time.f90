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

	module com_mod_ty_time
	use com_mod_ty_parameters,		only: ty_com_parameters

	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_com_time
		real(kind=dpd):: t			!<Actual time
		real(kind=dpd):: told		!<Time at previous step
		real(kind=dpd):: dt			!<Actual time
		real(kind=dpd):: dtold	!<Time at previous step
		real(kind=dpd):: tprint	!<Actual printtime
		integer(8)::iter_total
		logical				:: CheckPrint	!<If true then print
		type(ty_com_parameters),pointer::parameters
	contains
	procedure,public:: increase_time				=>	s_com_time_increase_time
	procedure,public:: factor_timestep			=>	s_com_time_factor_timestep
	procedure,public:: update_dt						=>	s_com_calc_update_dt
	!procedure,public:: construct=> s_unsat_model_construct
	!procedure,public:: print_timestep=> s_unsat_model_print_timestep

	end type ty_com_time

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_time_increase_time(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_time_increase_time" :: s_com_time_increase_time
	!DEC$ endif

	class(ty_com_time),intent(inout)::this

	real(kind=dpd)::ttemp
	ttemp = this%t
	this%t			= this%t+this%dt
	this%tprint = real(int(ttemp/this%parameters%tprintinc),dpd)*this%parameters%tprintinc+this%parameters%tprintinc

	!Adjust t to print time or increase with dt
	if(ttemp<this%tprint.and.this%t>=this%tprint) then
		this%t = this%tprint
		this%checkprint = .true.
		!this%tprint  = min(this%parameters%tmax,this%tprint+this%parameters%tprintinc) !Update: tprint= tprintinc·int(ttemp/tprintinc)+tprintinc
		!else if (this%t==this%tprint) then
		!  this%checkprint = .true.
	else
		this%checkprint = .false.
	end if
	this%t = min(this%t,this%parameters%tmax)

	end subroutine s_com_time_increase_time

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_time_factor_timestep(this,factor)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_time_factor_timestep" :: s_com_time_factor_timestep
	!DEC$ endif

	class(ty_com_time),intent(inout)::this
	real(kind=dpd),intent(in)::factor


	this%dt			= min(this%parameters%dtmax,max(factor*this%dt,this%parameters%dtmin))

	end subroutine s_com_time_factor_timestep


	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside ty_sat_model. Perform one iteration
	!> @param[inout] mx
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_calc_update_dt(this,iterconverg)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_calc_update_dt" :: s_com_calc_update_dt
	!DEC$ endif

	class(ty_com_time),intent(inout)::this
	integer,intent(in)::iterconverg

	if (iterconverg<this%parameters%it_inc_dt)  this%dt = min(this%parameters%dtinc*this%dt,this%parameters%dtmax)
	if (iterconverg>this%parameters%it_dec_dt)  this%dt = max(this%parameters%dtdec*this%dt,this%parameters%dtmin)

	end subroutine s_com_calc_update_dt

	end module com_mod_ty_time