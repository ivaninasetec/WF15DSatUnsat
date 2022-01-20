	module com_mod_ty_parameters

	implicit none
	include 'inc_precision.fi'

	private

	!******************************************************************************************************************
	! TY_COM_PARAMETERS
	! Derived type that defines calculation options
	!------------------------------------------------------------------------------------------------------------------
	!	 |epsh_tol        [real]:           Abs Tolerance on pore pressure for iteration convergence (L)
	!	 |epsth_tol       [real]:           Abs Tolerance on water content for iteration convergence (L3/L3)
	!	 |it_inc_dt       [integer]:        If iterations are above this number, time step is increased
	!	 |it_dec_dt       [integer]:        If iterations are below this number, time step is decreased
	!	 |it_max          [integer]:        If iterations ar above this number, then the actual time step is
	!  |                                  considered that dont converge, iterations start again with a
	!  |                                  decreased time step.
	!	 |crelax          [real]:           Coeficient of relaxation for Newton Method (>0) (1 means no relaxation)
	!	 |masslump        [logical]:        If true, then the mass coefficients are lumped in the diagonal, this
	!  |                                  is convenient to reduce spureus oscillations but can be less accurate.
	!	 |ErrorOnNode     [logical]:        If true, the Modified Picard is used, if false, the diferential equation
	!  |                                  is soved in h and the water content error control is done
	!  |                                  on elements with the L2 space-norm.
	!  |QuadratureOrder [integer]:        Number of points of elements integration for Legendre quadrature. (1-60)
	!	 |typesolution    [integer]:        Kind of matrix storage, solver and preconditioners used for the solution
	!  |                                  of the Linear algebra system:
	!  |                                    0: Dense matrix, direct gauss solver, no preconditioner
	!  |                                    1: Sparse CSR matrix, Direct Sparse Solver (Intel DSS),no preconditioner
	!  |                                    2: Sparse CSR matrix, FGMRES solver, ILUT preconditioner
	!  |                                    3: Banded matrix, banded direct solver, no preconditioner
	!	 |ccourant        [real]:           Coefficient of Courant to limit the timestepping with changes on solution
	!	 |Multksal        [real]:           Factor to multiply k on the seepage boundary (for saturated model)
	!	 |NODESPERELEMENT [integer]:        Number of nodes per element (at this moment 2).
	!	 |CLASSESOFNODES  [integer]:        Max class number of the node (at this moment 0, continuos function on C0)
	!	 |KL              [integer]:        Number of lower diagonals for banded matrix (calculated)
	!	 |KU              [integer]:        Number of upper diagonals for banded matrix (calculated)
	!******************************************************************************************************************

	type,public::ty_com_parameters !< Class with parameters for model calculation
		real(kind=dps)	::epsh_tol	!< Tolerance in pressure head [L]
		real(kind=dps)	::epsth_tol !< Tolerance in watercontent [-]
		real(kind=dps)	::epshsat_tol !< Tolerance in watercontent [-]
		integer					::it_inc_dt !< Iterations below which dt is increased
		integer					::it_dec_dt !< Iterations over which dt is decreased
		integer					::it_min		!< Minimum number of iterations on each time loop even if it converged for less iterations
		integer					::it_max		!< Max number of iterations, over which time step restart decreased
		real(kind=dps)	::crelax		!< Relaxation coefficient in the update pressure head in each iteration
		logical					::masslump	!< Mass lumping used? (.true. or .false.)
		logical					::erroronnode=.true. !< Is used the error on node or error on element? (.true.= error on node)
		integer					::quadratureorder=5	 !< Quadrature order for integration inside element
		integer					::typesolution !< Type of matrix solver (0 for dense, 1 csr-dss, 2 csr-fgmres, 3 banded-direct)
		integer					::typematrixstorage !< Type of matrix sparsity: 1 for dense, 2 sparse csr, 3 banded
		integer					::typepreconditioner !< Preconditioner for the solver (0 for none, 1 jacobi, 2 ilu0, 2 ilut)
		integer					::typesolver !< Type of solver (0 for gauss, 1 direct dss, 2 direct paradiso, 3 fgmres)
		real(kind=dps)	::ccourant = 1.0 !< Courant number
		real(kind=dps)	::multksal !< Permeability at the seepage surface compared to permeability on element.
		real(kind=dps)	::maxhsatinc !< Max saturated height.
		real(kind=dps)	::nrelmin !< Max saturated height.
		real(kind=dps)	::nrelmax !< Max saturated height.
		real(kind=dps)	::crelax_nrel !< Max saturated height.
		real(kind=dps)	::crelax_q !< Max saturated height.
		integer					::n_repetitions_satunsat
		!integer					::nodesperelement !< Nodes per element
		!integer					::classesofnodes  !< Classes per node
		!integer					::kl							!< Number of lower bands when band storage is used
		!integer					::ku							!< Number of upper bands when band storage is used

		!******************************************************************************************************************
		! TY_COM_TIME
		! Derived type that defines time options
		!------------------------------------------------------------------------------------------------------------------
		!	 |tinit           [real]:           Time at the begining of the calculations.
		!	 |dtinit          [real]:           Initial time step 'dt'.
		!	 |tmax            [integer]:        Final time of simulation
		!	 |dtinc           [integer]:        Factor to multiply 'dt' when 'it_inc_dt' is reached
		!	 |dtdec           [integer]:        Factor to multiply 'dt' when 'it_dec_dt' is reached
		!	 |dtmax           [real]:           Max value for time step 'dt'.
		!	 |dtmin           [logical]:        Min value for time step 'dt', below this the solution is explicit.
		!	 |tprintinc       [logical]:        Increment of time in which output record will be written.
		!******************************************************************************************************************

		real(kind=dps)	::tinit			!<Time at the begining of the calculations.
		real(kind=dps)	::dtinit		!<Initial time step 'dt'.
		real(kind=dps)	::tmax			!< Final time of simulation
		real(kind=dps)	::dtinc			!< Factor to multiply 'dt' when 'it_inc_dt' is reached
		real(kind=dps)	::dtdec			!< Factor to multiply 'dt' when 'it_dec_dt' is reached
		real(kind=dps)	::dtmax			!< Max value for time step 'dt'.
		real(kind=dps)	::dtmin			!< Min value for time step 'dt', below this the solution is explicit.
		real(kind=dps)	::tprintinc !< Increment of time in which output record will be written.
		logical					::isrestarttime=.false. !<If true then the model began to restart at time trestart from previous output
		real(kind=dps)	::trestart !< Increment of time in which output record will be written.

	end type ty_com_parameters

	end module com_mod_ty_parameters