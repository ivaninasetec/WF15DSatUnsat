	!********************************************************************************************************************
	! TITLE         : SAT_MOD_TY_NODES: EXTENDED DERIVED TYPE OF COM_MOD_TY_NODES TO INCLUDE PROPERTIES AND METHODS OF WF1DSAT
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_CALC
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Extended derived type of com_mod_ty_elements to include properties and methods of wf1dsat
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module sat_mod_ty_nodes
	use com_mod_ty_nodes, only: ty_com_nodes

	implicit none
	include 'inc_precision.fi'

	private

	type,extends(ty_com_nodes),public::ty_sat_nodes	!< CLASS: Definition of the layer in saturated model
		real(kind=dpd),allocatable::results_qent(:)
		real(kind=dpd),allocatable::results_incvoldt(:)
		real(kind=dpd),allocatable::results_qhor(:)
		real(kind=dpd),allocatable::results_dqhordx(:)
		real(kind=dpd),allocatable::results_dqhordx_from_incvoldt(:)
		real(kind=dpd),allocatable::results_qhor_all(:)
		real(kind=dpd),allocatable::results_dqhordx_all(:)
		real(kind=dpd),allocatable::results_dqhordx_from_incvoldt_all(:)

	contains
	procedure,public:: set_hinit			=> s_sat_set_hinit				!<Procedure to set initial h
	procedure,public:: set_z					=> s_sat_set_z						!<Procedure to set initial h

	end type ty_sat_nodes

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to define the initial h.
	!> @param[in] ne
	!> @param[in] nn
	!> @param[in] nc
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_set_hinit(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_set_hinit
	!DEC$ endif
	class(ty_sat_nodes),intent(inout)::this

	this%hinit = 0.0_dpd
	this%dhinit = 0.0_dpd

	end subroutine s_sat_set_hinit

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure to update the coordinate z from slope and lenght of the layer
	!> @param[in] slope
	!> @param[in] lenght
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_sat_set_z(this,slope,lenght)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT :: s_sat_set_z
	!DEC$ endif
	class(ty_sat_nodes),intent(inout)::this
	real(kind=dps),intent(in)::slope
	real(kind=dps),intent(in)::lenght

	integer::i

	do i=1, this%count
		this%z(i) = (lenght-this%x(i))*slope
	end do

	end subroutine s_sat_set_z

	end module sat_mod_ty_nodes