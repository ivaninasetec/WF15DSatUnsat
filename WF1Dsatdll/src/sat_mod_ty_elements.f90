	!********************************************************************************************************************
	!        EXTENSION OF CLASS TY_COM_ELEMENTS FOR SATURATED MODEL
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
		!real(kind=dpd),allocatable::x0(:)
		!real(kind=dpd),allocatable::x1(:)
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