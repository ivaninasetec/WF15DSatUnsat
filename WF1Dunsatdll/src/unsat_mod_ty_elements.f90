	!********************************************************************************************************************
	! TITLE         : UNSAT_MOD_TY_ELEMENTS: EXTENDED DERIVED TYPE OF COM_MOD_TY_ELEMENTS TO INCLUDE PROPERTIES AND METHODS OF WF1DUNSAT
	! PROJECT       : WF1DUNSATDLL
	! MODULE        : COM_MOD_TY_CALC
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Extended derived type of com_mod_ty_elements to include properties and methods of wf1dunsat
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module unsat_mod_ty_elements
	use com_mod_ty_elements, only: ty_com_elements

	implicit none
	include 'inc_precision.fi'

	private

	type,extends(ty_com_elements),public::ty_unsat_elements	!< CLASS: Definition of the layer in saturated model
		real(kind=dps),allocatable::thnew(:)	!<Current water content integrated on element (n+1,k+1)
		real(kind=dps),allocatable::thtemp(:)	!<Current water content integrated on element for previous iteration (n+1,k)
		real(kind=dps),allocatable::thold(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::results_qent(:)
		real(kind=dps),allocatable::results_qsal(:)
		real(kind=dps),allocatable::results_incvoldtperm(:)
		real(kind=dps),allocatable::results_qmed(:)
		real(kind=dps),allocatable::results_kmed(:)
		real(kind=dps),allocatable::results_dhdxmed(:)
		real(kind=dps),allocatable::results_dhxdxmed(:)
		real(kind=dps),allocatable::h0(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::h1(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::h0old(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::h1old(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::th0(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::th1(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::k0(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::k1(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::dhdx0(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::dhdx1(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::dhxdx0(:)	!<Water content integrated on element on the previous timestep (n)
		real(kind=dps),allocatable::dhxdx1(:)	!<Water content integrated on element on the previous timestep (n)

	end type ty_unsat_elements

	end module unsat_mod_ty_elements