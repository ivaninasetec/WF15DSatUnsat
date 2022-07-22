	!********************************************************************************************************************
	! TITLE         : COM_MOD_TY_TIME: DERIVED TYPE THAT DEFINE COMMON PARAMETERS RELATED TO TIME AND TIMESSTEPING
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_CALC
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Derived type that define common parameters related to time and timessteping
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
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

	end type ty_com_time

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Increase time and update time related parameters
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
	else
		this%checkprint = .false.
	end if
	this%t = min(this%t,this%parameters%tmax)

	end subroutine s_com_time_increase_time

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Increase Delta t by multiplying by factor, having into consideration dtmin and dtmax
	!> @param[inout] factor	Factor to multiply Dt
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
	!> Update increment of t depending on the number of iterations performed (dynmic timestepping)
	!> @param[inout] iterconvergence	Number of iterations in the previous convergence process
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