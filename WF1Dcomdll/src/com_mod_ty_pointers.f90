	!********************************************************************************************************************
	! TITLE         : COM_MOD_TY_POINTERS: DERIVED TYPES TO USE POINTERS
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : COM_MOD_TY_PARAMETERS
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Derived types to be able to use pointers (real or integers) in other derived types.
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************
	

	module com_mod_ty_pointers

	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_com_pointer_real
		real(kind=dpd),pointer:: p
	end type ty_com_pointer_real

	type,public::ty_com_pointer_integer
		integer,pointer:: p
	end type ty_com_pointer_integer

	end module com_mod_ty_pointers