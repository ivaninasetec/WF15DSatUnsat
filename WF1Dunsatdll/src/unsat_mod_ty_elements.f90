	!********************************************************************************************************************
	!        EXTENSION OF CLASS TY_COM_ELEMENTS FOR UNSATURATED MODEL
  !********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : UNSAT_MOD_TY_ELEMENTS
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
		!real(kind=dps),allocatable::results_incvol(:)
		!real(kind=dps),allocatable::results_incqhor(:)
		
	end type ty_unsat_elements
	


	
	
	end module unsat_mod_ty_elements