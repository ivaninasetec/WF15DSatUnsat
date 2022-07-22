	!********************************************************************************************************************
	! TITLE         : SAT_MOD_TY_ELEMENTS: EXTENDED DERIVED TYPE OF COM_MOD_TY_ELEMENTS TO INCLUDE PROPERTIES AND METHODS OF WF1DSAT
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
	
	module sat_mod_ty_elements
	use com_mod_ty_elements, only: ty_com_elements
	
  implicit none
  include 'inc_precision.fi'

  private
	
	type,extends(ty_com_elements),public::ty_sat_elements	!< CLASS: Definition of the layer in saturated model
	real(kind=dps),allocatable::results_qent(:)
	real(kind=dps),allocatable::results_incvoldt(:)		!dV/dt in element
	real(kind=dps),allocatable::results_dqhordx(:) !dqhor/dx in 1st layer in element
	real(kind=dps),allocatable::results_dqhordx_from_incvoldt(:) !dqhor/dx=qent-dV/dt in 1st layer 
	real(kind=dps),allocatable::results_dqhordx_all(:) !dqhor/dx in all layers in node
	real(kind=dps),allocatable::results_dqhordx_from_incvoldt_all(:) !dqhor/dx=qent-dV/dt in all layers 
		real(kind=dpd),allocatable::z(:)
		real(kind=dpd),allocatable::hz(:)
		real(kind=dpd),allocatable::hzold(:)
		real(kind=dpd),allocatable::dhdx(:)
		real(kind=dpd),allocatable::dhzdx(:)
		real(kind=dpd),allocatable::h0(:)
		real(kind=dpd),allocatable::h1(:)
		real(kind=dpd),allocatable::z0(:)
		real(kind=dpd),allocatable::z1(:)
		real(kind=dpd),allocatable::dhdx0(:)
		real(kind=dpd),allocatable::dhdx1(:)
		real(kind=dpd),allocatable::dhzdx0(:)
		real(kind=dpd),allocatable::dhzdx1(:)
		real(kind=dpd),allocatable::q0(:)
		real(kind=dpd),allocatable::q1(:)
		real(kind=dpd),allocatable::q0_all(:)
		real(kind=dpd),allocatable::q1_all(:)
		real(kind=dpd),allocatable::q(:)
		real(kind=dpd),allocatable::q_all(:)
		real(kind=dpd),allocatable::k0(:)
		real(kind=dpd),allocatable::k1(:)
		real(kind=dpd),allocatable::k0_all(:)
		real(kind=dpd),allocatable::k1_all(:)
		real(kind=dpd),allocatable::k(:)
	end type ty_sat_elements
	
	end module sat_mod_ty_elements