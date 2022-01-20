	!********************************************************************************************************************
	!        EXTENSION OF CLASS TY_COM_NODES FOR VERTICUAL UNSATURATED MODEL
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : MOD_SAT_TY_LAYER
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
	!> Class for horizontal saturated layer. Extend common class of layers.
	!********************************************************************************************************************

	module unsat_mod_ty_nodes
	use com_mod_ty_nodes, only: ty_com_nodes

	implicit none
	include 'inc_precision.fi'

	private

	type,extends(ty_com_nodes),public::ty_unsat_nodes	!< CLASS: Definition of the layer in saturated model
		real(kind=dpd),allocatable::results_qnewmann(:)
		!real(kind=dpd),allocatable::results_qent(:)
		!real(kind=dpd),allocatable::results_incvol(:)
		!real(kind=dpd),allocatable::results_qhor(:)
		!real(kind=dpd),allocatable::results_dqhordx(:)

	contains
	procedure,public:: set_hinit			=> s_unsat_set_hinit				!<Procedure to set initial h
	!procedure,public:: set_z					=> s_sat_set_z						!<Procedure to set initial h

	end type ty_unsat_nodes

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_nodes. Subroutine to define the initial h.
	!> @param[in] ne
	!> @param[in] nn
	!> @param[in] nc
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_set_hinit(this,zphr)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_set_hinit" :: s_unsat_set_hinit
	!DEC$ endif
	class(ty_unsat_nodes),intent(inout)::this
	real(kind=dpd)::zphr

	this%hinit = -this%x+zphr
	this%dhinit = 0.0_dpd

	end subroutine s_unsat_set_hinit

	end module unsat_mod_ty_nodes