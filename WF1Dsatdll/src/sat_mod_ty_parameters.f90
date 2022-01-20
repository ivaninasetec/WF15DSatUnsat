	!********************************************************************************************************************
	!        EXTENSION OF CLASS TY_COM_PARAMETERS FOR SATURATED MODEL
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
	
	module sat_mod_ty_parameters
	use com_mod_ty_parameters, only: ty_com_parameters
	
  implicit none
  include 'inc_precision.fi'

  private
	
	type,extends(ty_com_parameters),public::ty_sat_parameters	!< CLASS: Definition of the layer in saturated model

		
	end type ty_sat_parameters
	
	end module sat_mod_ty_parameters