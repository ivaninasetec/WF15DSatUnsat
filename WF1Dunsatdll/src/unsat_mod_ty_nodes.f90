	!********************************************************************************************************************
	! TITLE         : UNSAT_MOD_TY_NODES: EXTENDED DERIVED TYPE OF UNSAT_MOD_TY_NODES TO INCLUDE PROPERTIES AND METHODS OF WF1DUNSAT
	! PROJECT       : WF1DUNSATDLL
	! MODULE        : UNSAT_MOD_TY_NODES
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Extended derived type of unsat_mod_ty_nodes to include properties and methods of wf1dunsat
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module unsat_mod_ty_nodes
	use com_mod_ty_nodes, only: ty_com_nodes

	implicit none
	include 'inc_precision.fi'

	private

	type,extends(ty_com_nodes),public::ty_unsat_nodes	!< CLASS: Definition of the layer in saturated model
		real(kind=dpd),allocatable::results_qnewmann(:)

	contains
	procedure,public:: set_hinit			=> s_unsat_set_hinit				!<Procedure to set initial h

	end type ty_unsat_nodes

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Subroutine to define the initial h from the definition of the phreatic level.
	!> @param[in] zphr  Phreatic level
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