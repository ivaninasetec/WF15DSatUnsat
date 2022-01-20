	module com_mod_ty_boundary

	implicit none
	include 'inc_precision.fi'

	private

	!******************************************************************************************************************
	! TY_COM_BOUND
	! Derived type that defines properties of boundaries
	!------------------------------------------------------------------------------------------------------------------
	!	 |TopBoundbyH     [logical]:        If true, the top boundary is Dirichlet, defined by pore pressure in time.
	!	 |TopBoundbyQ     [logical]:        If true, the top boundary is Newmann, defined by waterflow in time.
	!	 |timestepscount  [integer]:        Total number of steps in the boundary table
	!	 |index           [integer]:        Index of the actual time step (to lookup in boundary table)
	!	 |timebound (:)   [real]:           Times of the boundary table t(index) is time at beginning of step: index
	!	 |hbound (:)      [real]:           Value of Dirichlet boundary pressure of the boundary table.
	!  |                                  hbound(index) is pressure at top at step index.
	!	 |qbound (:)      [real]:           Value of Newmann boundary water flow of the boundary table.
	!  |                                  qbound(index) is water flow at top at step index.
	!******************************************************************************************************************

	type, public::ty_com_boundary !< CLASS: boundary (common definition)
		logical::										topboundbyh=.false. !< If true, the top boundary is Dirichlet, defined by pore pressure in time.
		logical::										topboundbyq=.false. !< If true, the top boundary is Newmann, defined by waterflow in time.
		integer::										timestepscount !< Total number of steps in the boundary table

		integer::										index=1 !< Index of the actual time step (to lookup in boundary table)
		real(kind=dps),allocatable::timebound(:)!< Times of the boundary table t(index) is time at beginning of step: index
		real(kind=dps),allocatable::hbound(:) !< Value of Dirichlet boundary pressure of the boundary table. hbound(index,idnode) is pressure at top at step index.
		real(kind=dps),allocatable::qbound(:) !< Value of Newmann boundary water flow of the boundary table. qbound(index,idnode) is water flow at top at step index.

	contains
	procedure,public:: get_hbound_file => f_h_top_file
	procedure,public:: get_qbound_file => f_q_top_file

	end type ty_com_boundary

	contains

	!********************************************************************************************************************
	! F: f_h_top_file(T,BD)
	!--------------------------------------------------------------------------------------------------------------------
	! This function returns the head pressure at the boundary at time t. Discrete values of water pressures
	! has been previously readed from file and include in the class BD (of type: TY_COM_BOUNDARY), in BD%hbound(:)
	!
	!	This subROUTine check actual index (BD%index) and if t falls between t(index) and t(index+1) and increase index in
	! case necessary, modifying actual BD%index until t(index)<=t<t(index+1).
	! Then return q(index) (already stored in BD%hbound
	!********************************************************************************************************************

	function f_h_top_file(this,t)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_h_top_file" :: f_h_top_file
	!DEC$ endif

	!returns the head at boundary at time t from hbound(:)

	real(kind=dps),intent(in)::t
	class(ty_com_boundary),intent(inout)::this
	real(kind=dpd)::f_h_top_file

	!real(kind=dpd)::t0,t1,q0,q1

	!adjust index for actual time
	if(this%index>1) then
		do while(t<this%timebound(this%index).or.this%index==1)
			this%index = max(1,this%index-1)
		end do
	end if

	if (this%index<this%timestepscount) then
		do while (t>this%timebound(min(this%index+1,this%timestepscount)).and.this%index.ne.this%timestepscount)
			this%index = min(this%index + 1,this%timestepscount)
		end do
	end if

	!return value of table for actual time step
	f_h_top_file = this%hbound(this%index)

	end function f_h_top_file

	!********************************************************************************************************************
	! F: f_q_top_file(T,BD)
	!--------------------------------------------------------------------------------------------------------------------
	! This function returns the waterflow at the boundary at time t. Discrete values of water pressures
	! has been previously readed from file and include in the class BD (of type: TY_COM_BOUNDARY), in BD%qbound(:)
	!
	!	This subROUTine check actual index (BD%index) and if t falls between t(index) and t(index+1) and increase index in
	! case necessary, modifying actual BD%index until t(index)<=t<t(index+1).
	! Then return q(index) (already stored in BD%qbound
	!********************************************************************************************************************

	function f_q_top_file(this,t)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_q_top_file" :: f_q_top_file
	!DEC$ endif


	real(kind=dps),intent(in)::t
	class(ty_com_boundary),intent(inout)::this
	real(kind=dpd)::f_q_top_file

	if (this%index<this%timestepscount) then
		do while (t>this%timebound(min(this%index,this%timestepscount)).and.this%index.ne.this%timestepscount)
			this%index = min(this%index + 1,this%timestepscount)
		end do
	end if

	!return value of table for actual time step
	f_q_top_file = this%qbound(this%index)

	end function f_q_top_file

	end module com_mod_ty_boundary