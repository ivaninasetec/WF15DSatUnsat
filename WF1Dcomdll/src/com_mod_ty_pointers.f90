	!********************************************************************************************************************
	!        CLASS TO INCLUDE POINTER (THIS ALLOW TO HAVE AN ARRAY OF POINTERS WHEN INCLUDING IT IN OTHER CLASSES)
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : mod_com_ty_pointers
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
	!> Class for the collection of nodes-classes in the saturated model
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