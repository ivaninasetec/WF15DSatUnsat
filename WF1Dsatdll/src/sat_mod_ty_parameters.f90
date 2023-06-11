	!********************************************************************************************************************
	! TITLE         : SAT_MOD_TY_PARAMETERS: EXTENDED DERIVED TYPE OF COM_MOD_TY_PARAMETERS TO INCLUDE PROPERTIES AND METHODS OF WF1DSAT
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_CALC
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Extended derived type of com_mod_ty_parameters to include properties and methods of wf1dsat
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************
	
	module sat_mod_ty_parameters
	use com_mod_ty_parameters, only: ty_com_parameters
	
  implicit none
  include 'inc_precision.fi'

  private
	
	type,extends(ty_com_parameters),public::ty_sat_parameters	!< CLASS: Definition of the layer in saturated model

		
	end type ty_sat_parameters
	
	end module sat_mod_ty_parameters