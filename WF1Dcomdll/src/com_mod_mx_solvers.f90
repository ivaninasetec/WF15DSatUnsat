	!********************************************************************************************************************
	!*                                                                                                                  *
	!*                                 MODULE: MOD_LINALG_SOLVERS                                                       *
	!*                                                                                                                  *
	!********************************************************************************************************************
	!* THIS MODULE INCLUDE VARIABLY SATURATED HYDRAULIC FUNCTIONS BY VAN GENUCHTEN                                      *
	!*    SOLVE_LINALG(MX,SOLVER):   Calculate linear system								                                            *
	!*                                                                                                                  *
	!*    SOLVER = 0 Dense Matrix                                                                                       *
	!*    SOLVER = 1 CSR Matrix solver by Intel Direct Sparse Solver (Intel DSS)                                        *
	!*    SOLVER = 2 CSR MATRIX solved by FMGRES Iterative solver with ILUT Preconditioning                             *
	!*    SOLVER = 3 BANDED MATRIX solved by direct gaussian method                                                     *
	!*                                                                                                                  *
	!*    Iván Campos-Guereta Díez                                                                                      *
	!*    MSc Civil Engineering by Polytechnic University of Madrid                                                     *
	!*    PhD Student by University of Nottingham                                                                       *
	!*    eMBA by International Institute San Telmo in Seville                                                          *
	!*    ivan.camposguereta@nottingham.ac.uk                                                                           *
	!*                                                                                                                  *
	!*    This software is copyrighted 2019(C)                                                                          *
	!********************************************************************************************************************

	include 'mkl_dss.f90' ! include the standard dss "header file."
	include 'mkl_spblas.f90' !Sparse blas
	!!include 'mkl_pardiso.f90'
	!
	!include	'lapack.f90'

	!include	'blas.f90'

	module com_mod_mx_solvers

	implicit none

	private

	include 'inc_precision.fi'

	!********************************************************************************************************************
	! POLIMORPHIC INTERFACES FOR THE FOLLOWING LINEAR ALGEBRA SOLVER (select solver depending on arguments)
	!--------------------------------------------------------------------------------------------------------------------
	! SOLVE_DENSE(MX,RHS,SOLUTION,SOLVER)		:Dense storage
	! SOLVE_BANDED(MXBANDED,RHS,SOLUTION,SOLVER,MXBANDEDLS,(klarg),(kuarg)): Banded storage
	! SOLVE_CSR (MX,RHS,SOLUTION,SOLVER,CSRROWS,CSRCOLS,PERM)	: Sparse CSR storage
	!********************************************************************************************************************

	!interface SOLVE_LINALG
	!module procedure SOLVE_DENSE !MATRIX VECTOR MULTIPLICATION OF BANDED MATRIX
	!module procedure SOLVE_BANDED !MATRIX VECTOR MULTIPLICATION OF BANDED MATRIX
	!module procedure SOLVE_CSR !MATRIX VECTOR MULTIPLICATION OF BANDED MATRIX
	!end interface SOLVE_LINALG

	public:: solve_linalg

	contains

	!SUBROUTINE SOLVE_LINALG(MX,SOLVER)
	!!Solve a línear system given matrix in a given format
	!USE DATATYPES,ONLY:TYLINALG
	!
	!TYPE(TYLINALG),INTENT(INOUT)::MX
	!INTEGER,INTENT(IN)::SOLVER
	!!REAL(KIND=dpd),INTENT(IN)::MX(:),RHS(:)
	!!INTEGER,INTENT(IN),OPTIONAL::CSRROWS(:),CSRCOLS(:)
	!!INTEGER,INTENT(INOUT),OPTIONAL,ALLOCATABLE::PERM(:)
	!
	!!REAL(KIND=dpd),INTENT(OUT)::SOLUTION(:)
	!
	!SELECT CASE (SOLVER)
	!CASE (0) !General dense matrix (general matrix, dense storage)
	!     CALL SOLVE_DENSE(MX%MXDENSE,MX%RHS,MX%SOLUTION)
	!CASE (1) !CSR Sparse Direct DSS solver (general matrix, CSR storage)
	!		 CALL SOLVE_DSS(MX%MXCSR,MX%CSRROWS,MX%CSRCOLS,MX%RHS,MX%SOLUTION,MX%PERM)
	!CASE (2) !CSR iterative FMGRES with ILUT preconditioning (general matrix, CSR storage)
	!	   CALL SOLVE_FGMRES(MX%MXCSR,MX%CSRROWS,MX%CSRCOLS,MX%RHS,MX%SOLUTION)
	!CASE (3) !BANDED iterative FMGRES with preconditioning (general matrix, CSR storage)
	!	   CALL SOLVE_BANDED(MX%MXBANDED,MX%RHS,MX%SOLUTION,MX%MXBANDEDLS,MX%kl,MX%ku)
	!CASE DEFAULT
	!		STOP('No proper Linear Algebra Solver selected')
	!END SELECT
	!
	!END SUBROUTINE SOLVE_LINALG


	!********************************************************************************************************************
	! S: SOLVE_DENSE(MX,RHS,SOLUTION,SOLVER)
	!--------------------------------------------------------------------------------------------------------------------
	! Solves linear algebra system when dense matrix storage. For the following cases (solver):
	!	solver=0: lapack dgesv (solution for general unsymetric matrices)
	!	solver=other: STOP
	!
	! Parameters:
	! MX(nrows,ncols): Left hand side  matrix in dense format (input)
	! RHS(ncols):  Right hand side vector (output)
	! SOLUTION(ncols): Vector with solution
	! solver: 0 for dense, other=error
	!********************************************************************************************************************

	subroutine solve_linalg(mx,rhs,solution,solver,perm)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"solve_linalg" :: solve_linalg
	!DEC$ endif

	use com_mod_mx_datatypes

	class(ty_mx),intent(in)::mx
	real(kind=dpd),intent(in)::rhs(:)
	real(kind=dpd),intent(inout)::solution(:)
	integer,intent(in)::solver
	integer,intent(inout),allocatable,optional::perm(:)
	real(kind=dpd),allocatable,save::mxbandedls(:,:)

	select type (mx)
	type is (ty_mx_dense)
		call solve_dense(mx%mx,rhs,solution,solver)
	type is (ty_mx_csr)
		if(present(perm)) call solve_csr(mx%mxval,rhs,solution,solver,mx%csrrows,mx%csrcols,perm)
		if(.not.present(perm)) call solve_csr(mx%mxval,rhs,solution,solver,mx%csrrows,mx%csrcols)
	type is (ty_mx_banded)
		if(.not.allocated(mxbandedls)) allocate(mxbandedls(1+mx%ku+2*mx%kl,size(mx%mx,2)))
		call solve_banded(mx%mx,rhs,solution,solver,mxbandedls,mx%kl,mx%ku)
	end select
	end subroutine

	!********************************************************************************************************************

	subroutine solve_dense(mx,rhs,solution,solver)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"solve_dense" :: solve_dense
	!DEC$ endif

	real(kind=dpd),intent(in)::mx(:,:),rhs(:)
	real(kind=dpd),intent(inout)::solution(:)
	integer,intent(in)::solver

	select case (solver)
	case (0) !general dense matrix (general matrix, dense storage)
		call solve_dense_case_0(mx,rhs,solution)
		case default
		stop("arguments dont match for dense matrix linear algebra solver")
	end select

	contains

	subroutine solve_dense_case_0(mx,rhs,solution)

	use f95_precision, only: wp => sp
	use lapack95, only: gesv


	real(kind=dpd),intent(in)::mx(:,:),rhs(:)
	real(kind=dpd),intent(inout)::solution(:)

	integer::numnod,info
	integer:: ipiv(size(rhs))


	numnod = size(rhs)

	!call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
	solution = rhs

	!call gesv(  mx, solution, ipiv, info )

	call dgesv( numnod, 1, mx, numnod, ipiv, solution, numnod, info ) !v_rhs is overwritten with solution

	end subroutine solve_dense_case_0

	end subroutine	solve_dense


	!********************************************************************************************************************
	! S: SOLVE_BANDED(MXBANDED,RHS,SOLUTION,SOLVER,MXBANDEDLS,KLARG,KUARG)
	!--------------------------------------------------------------------------------------------------------------------
	! Solves linear algebra system when banded matrix storage. For the following cases (solver):
	!	solver=3: lapack dgbsv (solution for general unsymetric banded matrices)
	!	solver=other: STOP
	!
	! Parameters:
	! MXBANDED(ndiagonals,rcols): Left hand side  matrix in banded format (input) (the diagonals are in rows)
	! RHS(ncols):  Right hand side vector (output)
	! SOLUTION(ncols): Vector with solution
	! solver: 3 for banded, other=error
	! MXBANDED(nupperrows+1+ldiagonals,rcols): Include MX in banded storage but with the whole upper diagonals, needed
	!                                          for calculation purposes.
	! KLARG (optional): Number of lower diagonals
	! KuARG (optional): Number of upper diagonals
	!********************************************************************************************************************

	subroutine solve_banded(mxbanded,rhs,solution,solver,mxbandedls,klarg,kuarg)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"solve_banded" :: solve_banded
	!DEC$ endif


	real(kind=dpd),intent(in)::		mxbanded(:,:),rhs(:)
	real(kind=dpd),intent(inout)::	solution(:)
	integer,intent(in)::solver
	real(kind=dpd),intent(inout)::mxbandedls(:,:)
	integer, intent(in), optional:: klarg,kuarg

	integer::kl,ku

	if(present(klarg)) then
		kl = klarg
	else
		kl = 1
	end if

	if(present(kuarg)) then
		ku = kuarg
	else
		ku = 1
	end if

	select case (solver)
	case (0) !Direct solver for banded matrix (general matrix, dense storage)
		call solve_banded_direct_none(mxbanded,rhs,solution,mxbandedls,kl,ku)
		case default
		stop("arguments dont match for banded matrix linear algebra solver")
	end select

	contains

	subroutine solve_banded_direct_none(mxbanded,rhs,solution,mxbandedls,kl,ku)

	real(kind=dpd),intent(in)::		mxbanded(:,:),rhs(:)
	real(kind=dpd),intent(inout)::	solution(:)
	real(kind=dpd),intent(inout)::mxbandedls(:,:)
	integer, intent(in):: kl,ku


	integer::numnod,info
	integer:: ipiv(size(rhs))

	numnod = size(rhs)

	mxbandedls = 0.0_dpd
	mxbandedls(kl+1:2*kl+ku+1,:) = mxbanded
	solution = rhs
	!call dgbsv( n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info )
	call dgbsv(numnod,kl,ku,1,mxbandedls,1+2*kl+ku,ipiv,solution,numnod,info)

	end subroutine solve_banded_direct_none


	end subroutine solve_banded


	!********************************************************************************************************************
	! S: SOLVE_CSR(MX,RHS,SOLUTION,SOLVER,CSRROWS,CSRCOLS,PERM)
	!--------------------------------------------------------------------------------------------------------------------
	! Solves linear algebra system when sparse CSR matrix storage. For the following cases (solver):
	!	solver=1: Direct sparse solver without preconditioning (Intel DSS)
	!	solver=2: Iterative FGMRES solver with ILUT preconditioning (lapack dfgmre and dcsrilut)
	!
	! Parameters:
	! MX(nvalues): Left hand side  matrix in CSR format (input) (only contain values of mx with whole diagonal)
	! RHS(ncols):  Right hand side vector (output)
	! SOLUTION(ncols): Vector with solution
	! solver: 1 for DSS solver, 2 for FGMRES with ILUT, other: error
	! CSRROWS(nrows+1): Specify at wich column start the first non cero MX value in that row (end with total values+1)
	! CSRCOLS(nvalues): Specify the column of each MX sparse value
	! PERM(optional): Permutation matrix, can be included for not to recalcuate each time
	!********************************************************************************************************************

	subroutine solve_csr(mx,rhs,solution,solver,csrrows,csrcols,perm)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"solve_csr" :: solve_csr
	!DEC$ endif


	use mkl_dss

	implicit none

	real(kind=dpd),intent(in)::mx(:),rhs(:)
	real(kind=dpd),intent(inout)::solution(:)
	integer,intent(in)::solver
	integer,intent(in)::csrrows(:),csrcols(:)
	integer,intent(inout),allocatable,optional::perm(:)


	select case (solver)
	case (1) !general dense matrix (general matrix, dense storage)
		call solve_csr_dss(mx,rhs,solution,csrrows,csrcols,perm)
	case (3) !general dense matrix (general matrix, dense storage)
		call solve_csr_fgmres_ilut(mx,rhs,solution,csrrows,csrcols)
		case default
		stop("arguments dont match for banded matrix linear algebra solver")
	end select

	contains

	subroutine solve_csr_dss(mx,rhs,solution,csrrows,csrcols,perm)

	use mkl_dss

	implicit none


	real(kind=dpd),intent(in)::mx(:),rhs(:)
	integer,intent(in)::csrrows(:),csrcols(:)
	integer,intent(inout),allocatable,optional::perm(:)
	real(kind=dpd2),intent(inout)::solution(:)
	!type(tymatrix),intent(inout)::	mx
	!type(mkl_pardiso_handle)::pt(64)

	integer :: error,i,optcreate,nrhs,optstructure,optfactor,optsolve, opt, numnp,ncount
	! define the data arrays and the solution and rhs vectors.
	type(mkl_dss_handle) :: handle ! allocate storage for the solver handle.
	character*15::statin


	numnp = size(csrrows)-1
	opt = mkl_dss_defaults

	!mkl_dss_defaults: mkl_dss_refinement_on
	optcreate = mkl_dss_defaults

	!mkl_dss_non_symmetric, mkl_dss_symmetric, mkl_dss_symmetric_structure
	optstructure = mkl_dss_non_symmetric

	!mkl_dss_positive_definite, mkl_dss_indefinite, mkl_dss_hermitian_positive_definite,mkl_dss_hermitian_indefinite
	optfactor = mkl_dss_defaults

	!mkl_dss_forward_solve, mkl_dss_diagonal_solve,mkl_dss_backward_solve | mkl_dss_refinement_off, mkl_dss_refinement_on
	optsolve = mkl_dss_defaults
	nrhs = 1
	! initialize the solver.
	error = dss_create(handle, optcreate)
	if (error /= mkl_dss_success) goto 999
	! define the non-zero structure of the matrix.
	ncount = size(csrcols)
	error = dss_define_structure(handle, optstructure, csrrows, numnp, numnp, csrcols, ncount)
	if (error /= mkl_dss_success) goto 999
	! reorder the matrix.

	if(present(perm)) then
		if (.not.allocated(perm)) then
			allocate(perm(numnp))
			error = dss_reorder(handle, mkl_dss_get_order, perm)
		else if (perm(1)==0) then
			error = dss_reorder(handle, mkl_dss_get_order, perm)
		else
			error = dss_reorder(handle, mkl_dss_my_order, perm)
		end if
	end if

	if (error /= mkl_dss_success) goto 999
	! factor the matrix.
	error = dss_factor_real(handle, optfactor, mx)
	if (error /= mkl_dss_success) goto 999
	! allocate the solution vector and solve the problem.
	error = dss_solve_real(handle, optsolve, rhs, nrhs, solution)
	if (error /= mkl_dss_success) goto 999
	! deallocate solver storage and various local arrays.
	error = dss_delete(handle, opt)
	if (error /= mkl_dss_success) goto 999

	goto 1000
	! print an error message and exit
999 write(*,*) "solver returned error code ", error
	stop
1000 continue

	end subroutine solve_csr_dss

	!----------------------------------------------------------------------------------------------------

	subroutine solve_csr_fgmres_ilut(mx,rhs,solution,csrrows,csrcols)

	use mkl_spblas

	implicit none

	real(kind=dpd),parameter::relative_tolerance=1.0d-6 !default value 10e-6
	real(kind=dpd),parameter::absolute_tolerance=1.0d-6 !default value 10e-0
	real(kind=dpd),parameter::residual_tolerance=1.0d-8 !initially 10e-3 but it is too high
	real(kind=dpd),parameter::tol_norm_gen_vector= 1.d-12!max norm of generated vector (coef hk+1,k of hessenberg matrix). initially 1.0d-12 (let that)
	real(kind=dpd),parameter::tol_zero_norm= 1.0d-12 !tolerance for zero norm, has to be lower than tol_norm_gen_vector
	integer,parameter::restartiterations=10
	integer,parameter::maxiterations=100
	!maximum fill-in, which is half of the preconditioner bandwidth. the number of non-zero elements
	!in the rows of the preconditioner cannot exceed (2*maxfil+1).
	integer,parameter::maxfil=5 !the max number of non zero elements in a row will be 2·maxfil+1

	!tolerance fhor threshold criterion for the resulting entries of the preconditioner
	real(kind=dpd),parameter::ilut_tol=1.0e-10 !tolerance of ilut


	real(kind=dpd),intent(in)::mx(:),rhs(:)
	integer,intent(in)::csrrows(:),csrcols(:)
	real(kind=dpd2),intent(inout)::solution(:)

	!type(tymatrix),intent(inout)::	mx
	!type(tynod),intent(inout)::	no

	type(sparse_matrix_t)::csra,csrl  ! structure with sparse matrix
	type(matrix_descr)::descra,descrl,descru    ! sparse matrix descriptor

	integer, parameter::sizep=128
	integer::i,info,rci_request,numnp,ierr,numilu
	real(kind=dpd2)::alpha, beta, eps,brhs(size(csrrows)-1),residual(size(csrrows)-1),b(size(csrrows)-1),trvec(size(csrrows)-1),tol
	real(kind=dpd2),save,allocatable::tmp(:)
	real(kind=dpd)::expected_solution(size(csrrows)-1),computed_solution(size(csrrows)-1),dpar(sizep)
	integer:: ipar(sizep)
	integer:: itercount
	integer:: sizetmp
	real(kind=dpd):: dnrm2,dvar
	real(kind=dpd2)::bilut((2*maxfil+1)*(size(csrrows)-1)-maxfil*(maxfil+1)+1)
	integer::ibilut(size(csrrows)),jbilut((2*maxfil+1)*(size(csrrows)-1)-maxfil*(maxfil+1)+1)
	!real(kind=dpd2)::bilut((2*maxfil+1)*(size(csrrows)-1)-maxfil*(maxfil+1)+1)
	!integer::ibilut(size(csrrows)),jbilut((2*maxfil+1)*(size(csrrows)-1)-maxfil*(maxfil+1)+1)

	!numilu = (2*maxfil+1)*(size(csrrows)-1)-maxfil*(maxfil+1)+1
	numnp = size(csrrows)-1
	!real(kind=dpd): alpha, beta, eps
	!double precision alpha, beta, eps
	alpha = 1.0_dpd
	beta  = 0.0_dpd
	!   ! parameters ipar defined in: https://software.intel.com/en-us/mkl-developer-reference-c-fgmres-interface-description#rci_id_commonparameters
	ipar(15) = restartiterations !restart after x iterations
	!ipar(4) = 1000 !maximum iterations
	!ipar(6)=0 !0 Means no warning messages
	!ipar(7)=0 !0 Means no warning messages
	!ipar(8) = 1 !if 0 do not do the stopping test for the maximal number of iterations. stopping is when ipar[3]≤ipar[4] (ipar(3) is current iteration number)
	!ipar(11) = 0 !if the value is 0.0, the dfgmres routine runs the non-preconditioned version of the fgmres method. otherwise preconditioned. and need step setting rci_request=3.
	!dpar(1) = relative_tolerance !specifies the relative tolerance. the default value is 1.0e-6.
	!!dpar(2) = 1.0d-10 !specifies the absolute tolerance. the default value is 0.0e+0.

	!if(.not.allocated(tmp)) allocate(tmp(no%numnp*(2*no%numnp+1)+(no%numnp*(no%numnp+9))/2+1))
	!sizetmp=(2*ipar(15) + 1)*numnp + ipar(15)*(ipar(15) + 9)/2 + 1
	sizetmp=max((2*ipar(15) + 1)*numnp + ipar(15)*(ipar(15) + 9)/2 + 1,(numnp*(2*numnp+1)+(numnp*(numnp+9))/2+1))
	if(.not.allocated(tmp)) allocate(tmp(sizetmp))

	!   create matrix descriptor
	descra % type = sparse_matrix_type_general
	descra % mode = sparse_fill_mode_upper
	descra % diag = sparse_diag_non_unit


	!   create csr matrix
	i = mkl_sparse_d_create_csr(csra,sparse_index_base_one,numnp,numnp,csrrows(1:numnp),csrrows(2:numnp+1),csrcols,mx)
	!rhs = mx%bx
	residual = 0.0_dpd
	tmp = 0.0_dpd2
	expected_solution=solution
	computed_solution=solution

	!---------------------------------------------------------------------------
	! initialize variables and the right hand side through matrix-vector product
	!---------------------------------------------------------------------------
	solution=rhs
	info = mkl_sparse_d_mv(sparse_operation_non_transpose,alpha,csra,descra,expected_solution,beta,solution)
	!solution=rhs
	!rhs = mx%bx
	!---------------------------------------------------------------------------
	! initialize the solver
	!---------------------------------------------------------------------------


	call dfgmres_init(numnp, computed_solution, rhs, rci_request, ipar, dpar, tmp)


	if (rci_request .ne. 0) go to 999


	!---------------------------------------------------------------------------
	! calculate ilut preconditioner.
	!                      !attention!
	! dcsrilut routine uses some ipar, dpar set by dfgmres_init routine.
	! important for dcsrilut default entries set by dfgmres_init are
	! ipar(2) = 6 - output of error messages to the screen,
	! ipar(6) = 1 - allow output of error messages,
	! ipar(31)= 0 - abort dcsrilut calculations if routine meets zero diagonal element.
	! ipar(7) = 1 - output warn messages if any and continue
	!
	! if ilut is going to be used out of intel(r) mkl fgmres context, than the values
	! of ipar(2), ipar(6), ipar(31), and dpar(31), should be user
	! provided before the dcsrilut routine call.
	!
	! in this example, specific for dcsrilut entries are set in turn:
	! ipar(31)= 1 - change small diagonal value to that given by dpar(31),
	! dpar(31)= 1.d-5  instead of the default value set by dfgmres_init.
	!                  it is the target value of the diagonal value if it is
	!                  small as compared to given tolerance multiplied
	!                  by the matrix row norm and the routine should
	!                  change it rather than abort dcsrilut calculations.
	!---------------------------------------------------------------------------

	ipar(31) = 1 !if the value of diagonal is 0 then set the value to dpar(31)
	dpar(31) = 10.0*ilut_tol !need to be equal or greater than ilut_tol
	dpar(8) = tol_zero_norm !tolerance for zero norm
	!tol = 1.d-20

	!ilut preconditioner based on the incomplete lu factorization with a threshold of a sparse matrix.
	call dcsrilut(numnp, mx, csrrows, csrcols, bilut, ibilut, jbilut, ilut_tol, maxfil, ipar, dpar, ierr)
	if (ierr.ne.0) then
		write(*,*) 'error: ',ierr
		stop
	end if

	info = mkl_sparse_d_create_csr(csrl,sparse_index_base_one,numnp,numnp,ibilut,ibilut(2),jbilut,bilut)
	!nrm2 = dnrm2(matsize, bilut, incx)
	!---------------------------------------------------------------------------
	! set the desired parameters:
	! do the restart after 2 iterations
	! logical parameters:
	! do not do the stopping test for the maximal number of iterations
	! do the preconditioned iterations of fgmres method
	! double precision parameters
	! set the relative tolerance to 1.0d-3 instead of default value 1.0d-6
	!---------------------------------------------------------------------------
	!! parameters ipar defined in: https://software.intel.com/en-us/mkl-developer-reference-c-fgmres-interface-description#rci_id_commonparameters
	!        ipar(15) = 5000 !restart after x iterations
	!        ipar(8) = 0 !if 0 do not do the stopping test for the maximal number of iterations. stopping is when par[3]≤ipar[4]
	!        ipar(11) = 0 !if the value is 0.0, the dfgmres routine runs the non-preconditioned version of the fgmres method. otherwise preconditioned. and need step setting rci_request=3.
	!        dpar(1) = 1.0d-20 !specifies the relative tolerance. the default value is 1.0e-6.
	!        !dpar(2) = 1.0d-10 !specifies the absolute tolerance. the default value is 0.0e+0.

	! parameters ipar defined in: https://software.intel.com/en-us/mkl-developer-reference-c-fgmres-interface-description#rci_id_commonparameters
	!ipar(15) = restartiterations !restart after x iterations (if deactivated then the non restarted version of fgmres is used) (unchecking this can cause inconsistencies, warnings and the stop of the program)

	ipar(5) = maxiterations !this is the maximum number of iterations
	ipar(8) = 0 !if 0 do not do the stopping test for the maximal number of iterations. stopping is when ipar[3]≤ipar[4] (ipar(3) is current iteration number) (default 0)
	ipar(11) = 1 !if the value is 0.0, the dfgmres routine runs the non-preconditioned version of the fgmres method. otherwise preconditioned. and need step setting rci_request=3.
	dpar(1) = relative_tolerance !specifies the relative tolerance. the default value is 1.0e-6.
	dpar(2) = absolute_tolerance !specifies the absolute tolerance. the default value is 0.0e+0.

	!---------------------------------------------------------------------------
	! check the correctness and consistency of the newly set parameters
	!---------------------------------------------------------------------------

	call dfgmres_check(numnp, computed_solution, rhs, rci_request, ipar, dpar, tmp)
	if (rci_request .ne. 0) go to 999
	!---------------------------------------------------------------------------
	! compute the solution by rci (p)fgmres solver with preconditioning
	! reverse communication starts here
	!---------------------------------------------------------------------------
1	call dfgmres(numnp, computed_solution, rhs, rci_request, ipar, dpar, tmp)
	!---------------------------------------------------------------------------
	! if rci_request=0.0, then the solution was found with the required precision
	!---------------------------------------------------------------------------
	if (rci_request .eq. 0) go to 3
	!---------------------------------------------------------------------------
	! if rci_request=1, then compute the vector a*tmp(ipar[22]-1:ipar[22]+n-2)
	! and put the result in vector tmp(ipar(23)-1:ipar[23]+n-2)
	!---------------------------------------------------------------------------
	if (rci_request .eq. 1) then
		info = mkl_sparse_d_mv(sparse_operation_non_transpose,alpha,csra,descra,tmp(ipar(22)),beta,tmp(ipar(23)))
		!call mkl_dcsrgemv('n', numnp, mx, csrrows, csrcols, tmp(ipar(22):ipar(22)+numnp-1), tmp(ipar(23):ipar(23)+numnp-1)) !doesnt work
		go to 1
	end if
	!---------------------------------------------------------------------------
	! if rci_request=2, then do the user-defined stopping test
	! the residual stopping test for the computed solution is performed here
	!---------------------------------------------------------------------------
	! note: from this point vector b(n) is no longer containing the right-hand
	! side of the problem! it contains the current fgmres approximation to the
	! solution. if you need to keep the right-hand side, save it in some other
	! vector before the call to dfgmres routine. here we saved it in vector
	! rhs(n). the vector b is used instead of rhs to preserve the original
	! right-hand side of the problem and guarantee the proper restart of fgmres
	! method. vector b will be altered when computing the residual stopping
	! criterion!
	!---------------------------------------------------------------------------
	if (rci_request .eq. 2) then
		! request to the dfgmres_get routine to put the solution into b(n) via ipar(13)
		ipar(13) = 1 !if 0 updates de solution to vector x, if positive write the solution to b
		! get the current fgmres solution in the vector b(n)
		call dfgmres_get(numnp, computed_solution, solution, rci_request, ipar,dpar, tmp, itercount)
		! compute the current true residual via intel(r) mkl (sparse) blas routines
		info = mkl_sparse_d_mv(sparse_operation_non_transpose,alpha,csra,descra,solution,beta,residual)
		!call mkl_dcsrgemv('n', numnp, mx, csrrows, csrcols, solution, residual)
		!call daxpy(numnp, -1.0d0.0, rhs, 1, residual, 1)
		residual = residual-rhs
		dvar = dnrm2(numnp, residual, 1)
		if (dvar .lt. residual_tolerance) then !initially 1.0e-3
			go to 3 !go to 3 if success
		else
			go to 1 !start new iteration
		end if
	end if
	!---------------------------------------------------------------------------
	! if rci_request=3, then apply the preconditioner on the vector
	! tmp(ipar(22)-1:ipar(22)+n-2) and put the result in vector tmp(ipar(23)-1:ipar(23)+n-2)
	!---------------------------------------------------------------------------
	if (rci_request .eq. 3) then
		!ilut preconditioner
		!mkl_?csrtrsv = triangular solver with simplified interface for sparse matrix in csr format
		!with one based indexing
		!call mkl_dcsrtrsv(uplo, transa, diag, m, a, ia, ja, x, y)
		call mkl_dcsrtrsv('l','n','u', numnp, bilut, ibilut, jbilut, tmp(ipar(22)),trvec)
		call mkl_dcsrtrsv('u','n','n', numnp, bilut, ibilut, jbilut, trvec, tmp(ipar(23)))
		!descrl % type = sparse_matrix_type_triangular
		!descrl % mode = sparse_fill_mode_lower
		!descrl % diag = sparse_diag_unit
		!info = mkl_sparse_d_trsv(sparse_operation_non_transpose,alpha,csrl,descrl,tmp(ipar(22)),trvec)
		!
		!descrl % mode = sparse_fill_mode_upper
		!descrl % diag = sparse_diag_non_unit
		!info = mkl_sparse_d_trsv(sparse_operation_non_transpose,alpha,csrl,descrl,trvec,tmp(ipar(23)))

		go to 1
	end if
	!---------------------------------------------------------------------------
	! if rci_request=4, then check if the norm of the next generated vector is
	! not zero up to rounding and computational errors. the norm is contained
	! in dpar(7) parameter
	!---------------------------------------------------------------------------
	if (rci_request .eq. 4) then
		if (dpar(7) .le. tol_norm_gen_vector) then
			go to 3
		else
			go to 1
		end if
		!---------------------------------------------------------------------------
		! if rci_request=anything else, then dfgmres subroutine failed
		! to compute the solution vector: computed_solution(n)
		!---------------------------------------------------------------------------
	else
		go to 999
	end if
	!---------------------------------------------------------------------------
	! reverse communication ends here
	! get the current iteration number and the fgmres solution. (do not forget to
	! call dfgmres_get routine as computed_solution is still containing
	! the initial guess!). request to dfgmres_get to put the solution into
	! vector computed_solution(n) via ipar(13)
	!---------------------------------------------------------------------------
3	ipar(13) = 0 !if 0 updates the solution to vector x, positive to the rhs vector b.
	call dfgmres_get(numnp, computed_solution, rhs, rci_request, ipar, dpar, tmp, itercount)
	!---------------------------------------------------------------------------
	! release internal intel(r) mkl memory that might be used for computations
	! note: it is important to call the routine below to avoid memory leaks
	! unless you disable intel(r) mkl memory manager
	!---------------------------------------------------------------------------
	call mkl_free_buffers

	!dvar = dnrm2(n,expected_solution,1) !euclidean norm

	solution=computed_solution

	go to 1000

	!---------------------------------------------------------------------------
	! release internal intel(r) mkl memory that might be used for computations
	! note: it is important to call the routine below to avoid memory leaks
	! unless you disable intel(r) mkl memory manager
	!---------------------------------------------------------------------------
999 write( *,'(a,a,i5)') 'this example failed as the solver has', &
		&  ' returned the error code', rci_request
	info = mkl_sparse_destroy(csra)
	info = mkl_sparse_destroy(csrl)
	call mkl_free_buffers
	stop 1

1000 end subroutine solve_csr_fgmres_ilut



	end subroutine solve_csr



	!---------------------------------------------------------------------------------------------------------------






	! *******************************************************************************************************************************



	!	!*************************************************************************************************************************************
	!
	!	SUBROUTINE SOLVE_GMRES_SCALED(NO,MX)
	!
	!	USE MKL_SPBLAS
	!
	!	IMPLICIT NONE
	!
	!
	!	TYPE(TYMATRIX),INTENT(INOUT)::	MX
	!	TYPE(TYNOD),INTENT(INOUT)::	NO
	!
	!	TYPE(SPARSE_MATRIX_T)::csrA,csrD,csrAScaled ! Structure with sparse matrix
	!	TYPE(MATRIX_DESCR)::descrA,descrD,descrAScaled    ! Sparse matrix descriptor
	!
	!	INTEGER, PARAMETER::SIZE=128
	!	INTEGER::i,j,INFO,RCI_REQUEST
	!	REAL(KIND=dpd2)::alpha, beta, rhs(NO%NumNP),residual(NO%NumNP),B(NO%NumNP)
	!	REAL(KIND=dpd2),ALLOCATABLE::tmp(:)
	!	REAL(KIND=dpd)::EXPECTED_SOLUTION(NO%NumNP),COMPUTED_SOLUTION(NO%NumNP),DPAR(SIZE)
	!	REAL(KIND=dpd)::D(NO%NumNP)
	!	INTEGER::DROW(NO%NumNP+1)
	!	INTEGER::DCOL(NO%NumNP)
	!
	!	INTEGER:: IPAR(SIZE)
	!	INTEGER:: ITERCOUNT
	!	INTEGER:: SIZETMP
	!	REAL(KIND=dpd):: DNRM2,DVAR
	!
	!	alpha = 1.0_dpd
	!	beta  = 0.0_dpd
	!	DROW = (/(i,i=1,NO%NumNP+1)/)
	!	DCOL = (/(i,i=1,NO%NumNP)/)
	!
	!
	!	! Parameters IPAR defined in: https://software.intel.com/en-us/mkl-developer-reference-c-fgmres-interface-description#RCI_ID_COMMONPARAMETERS
	!	IPAR(15) = 2 !Restart after X iterations
	!	ipar(4) =  100 !Maximum iterations
	!	IPAR(8) = 0 !if 0 do not do the stopping test for the maximal number of iterations. Stopping is when ipar[3]≤ipar[4] (ipar(3) is current iteration number)
	!	IPAR(11) = 0 !if the value is 0.0, the dfgmres ROUTine runs the non-preconditioned version of the FGMRES method. Otherwise preconditioned. And need step setting RCI_request=3.
	!	DPAR(1) = 1.0D-20 !Specifies the relative tolerance. The default value is 1.0e-6.
	!	!DPAR(2) = 1.0D-8 !Specifies the absolute tolerance. The default value is 0.0e+0.
	!
	!	!IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(NO%NumNp*(2*NO%NumNp+1)+(NO%NumNp*(NO%NumNp+9))/2+1))
	!	SIZETMP=NO%NumNp*(2*NO%NumNp+1)+(NO%NumNp*(NO%NumNp+9))/2+1
	!	IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(SIZETMP))
	!
	!	!   Create matrix descriptor
	!	descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
	!	descrA % MODE = SPARSE_FILL_MODE_UPPER
	!	descrA % DIAG = SPARSE_DIAG_NON_UNIT
	!	descrD % TYPE = SPARSE_MATRIX_TYPE_GENERAL
	!	descrD % MODE = SPARSE_FILL_MODE_UPPER
	!	descrD % DIAG = SPARSE_DIAG_NON_UNIT
	!	descrAScaled % TYPE = SPARSE_MATRIX_TYPE_GENERAL
	!	descrAScaled % MODE = SPARSE_FILL_MODE_UPPER
	!	descrAScaled % DIAG = SPARSE_DIAG_NON_UNIT
	!
	!
	!
	!
	!
	!
	!	!!Also modifie AR matrix by diagonally multiply for D: (Preconditioning scale)
	!	!DO i=1,NO%NumNP
	!	!	DO j=1,NO%NumNP
	!	!			IF(MX%AR_POS(i,j) .NE. 0) MX%AR_values(mx%AR_POS(i,j))=MX%AR_values(mx%AR_POS(i,j))*D(i)
	!	!	END DO
	!	!END DO
	!
	!	!!   Create CSR matrix
	!	i = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,NO%numnp,NO%numnp,MX%AR_CSR3_rowindex(1:NO%NumNP),MX%AR_CSR3_rowindex(2:NO%NumNP+1),MX%AR_CSR3_col,MX%AR_values)
	!
	!	!Scaling matrix is a diagonal matrix with the inverse of NO%HNew, doing that for all expected solution are aproximated to 1. Except for 0 values in which scale is not performed.
	!
	!	!Scale preconditioner:
	!	!D = MERGE(ABS(NO%HNew),1.0_dpd,ABS(NO%HNew)>1.0E-10_dpd)
	!
	!
	!	!JACOBI PRECONDITIONER:-------------------------------------------------------
	!	!DO i=1,NO%NumNP
	!	!D = MERGE(1/MX%AR_values(mx%AR_POS(i,i)),1.0_dpd,abs(MX%AR_values(mx%AR_POS(i,i)))>1.0E-10_dpd)
	!	!END DO
	!
	!	D = 1.0_dpd
	!
	!	descrA % TYPE = SPARSE_MATRIX_TYPE_DIAGONAL
	!	i = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,1.0_dps,csrA,descrA,D,0.0_dps,D)
	!	descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
	!
	!
	!	D = MERGE(D,1.0_dpd,ABS(D)>1.0E-10)
	!	D = 1/D
	!
	!	i = MKL_SPARSE_D_CREATE_CSR(csrD,SPARSE_INDEX_BASE_ONE,NO%numnp,NO%numnp,DROW(1:NO%NumNP),DROW(2:NO%NumNP+1),DCOL,D)
	!
	!
	!	!MX%BX = MX%BX/D !Scaled RHS as preconditioner
	!	residual = 0.0_dpd
	!	tmp = 0.0_dpd
	!
	!
	!	!Multiply sparse matrix.
	!	!stat = mkl_sparse_spmm (operation, A, B, C)
	!	i = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE,csrA,csrD,csrAScaled)
	!
	!	!descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
	!	!descrA % MODE = SPARSE_FILL_MODE_UPPER
	!	!descrA % DIAG = SPARSE_DIAG_NON_UNIT
	!	!
	!	!descrD % TYPE = SPARSE_MATRIX_TYPE_DIAGONAL
	!	!descrD % MODE = SPARSE_FILL_MODE_UPPER
	!	!descrD % DIAG = SPARSE_DIAG_NON_UNIT
	!	!i = mkl_sparse_sp2m (SPARSE_OPERATION_NON_TRANSPOSE, descrA, csrA, SPARSE_OPERATION_NON_TRANSPOSE, descrD, csrD, SPARSE_STAGE_FULL_MULT, csrAScaled)
	!	!IF (i .NE. SPARSE_STATUS_SUCCESS) STOP
	!
	!
	!	EXPECTED_SOLUTION=NO%HNew/D !Expected to be 1.0 except for values near to 0.0
	!	!COMPUTED_SOLUTION=NO%HNew
	!	COMPUTED_SOLUTION=EXPECTED_SOLUTION
	!	!---------------------------------------------------------------------------
	!	! Initialize variables and the right hand side through matrix-vector product
	!	!---------------------------------------------------------------------------
	!	!info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,EXPECTED_SOLUTION,beta,RHS)
	!	B=MX%BX !Temporal right hand side
	!	!rhs = MX%BX
	!	!---------------------------------------------------------------------------
	!	! Initialize the solver
	!	!---------------------------------------------------------------------------
	!	CALL DFGMRES_INIT(NO%NumNP, COMPUTED_SOLUTION, MX%BX, RCI_REQUEST, IPAR, DPAR, TMP)
	!	IF (RCI_REQUEST .NE. 0) GO TO 999
	!	!---------------------------------------------------------------------------
	!	! Set the desired parameters:
	!	! do the restart after 2 iterations
	!	! LOGICAL parameters:
	!	! do not do the stopping test for the maximal number of iterations
	!	! do the Preconditioned iterations of FGMRES method
	!	! DOUBLE PRECISION parameters
	!	! set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
	!	!---------------------------------------------------------------------------
	!	!! Parameters IPAR defined in: https://software.intel.com/en-us/mkl-developer-reference-c-fgmres-interface-description#RCI_ID_COMMONPARAMETERS
	!	!        IPAR(15) = 5000 !Restart after X iterations
	!	!        IPAR(8) = 0 !if 0 do not do the stopping test for the maximal number of iterations. Stopping is when par[3]≤ipar[4]
	!	!        IPAR(11) = 0 !if the value is 0.0, the dfgmres ROUTine runs the non-preconditioned version of the FGMRES method. Otherwise preconditioned. And need step setting RCI_request=3.
	!	!        DPAR(1) = 1.0D-20 !Specifies the relative tolerance. The default value is 1.0e-6.
	!	!        !DPAR(2) = 1.0D-10 !Specifies the absolute tolerance. The default value is 0.0e+0.
	!
	!
	!	!---------------------------------------------------------------------------
	!	! Check the correctness and consistency of the newly set parameters
	!	!---------------------------------------------------------------------------
	!	CALL DFGMRES_CHECK(NO%NumNP, COMPUTED_SOLUTION, MX%BX, RCI_REQUEST, IPAR, DPAR, TMP)
	!	IF (RCI_REQUEST .NE. 0) GO TO 999
	!	!---------------------------------------------------------------------------
	!	! Compute the solution by RCI (P)FGMRES solver with preconditioning
	!	! Reverse Communication starts here
	!	!---------------------------------------------------------------------------
	!1	CALL DFGMRES(NO%NumNP, COMPUTED_SOLUTION, MX%BX, RCI_REQUEST, IPAR, DPAR, TMP)
	!	!---------------------------------------------------------------------------
	!	! If RCI_REQUEST=0.0, then the SOLUTION WAS FOUND with the required precision
	!	!---------------------------------------------------------------------------
	!	IF (RCI_REQUEST .EQ. 0) GO TO 3
	!	!---------------------------------------------------------------------------
	!	! If RCI_REQUEST=1, then compute the vector A*TMP(ipar[22]-1:ipar[22]+n-2)
	!	! and put the result in vector TMP(IPAR(23)-1:ipar[23]+n-2)
	!	!---------------------------------------------------------------------------
	!	IF (RCI_REQUEST .EQ. 1) THEN
	!		info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrAScaled,descrA,TMP(IPAR(22)),beta,TMP(IPAR(23)))
	!		GO TO 1
	!	END IF
	!	!---------------------------------------------------------------------------
	!	! If RCI_request=2, then do the user-defined stopping test
	!	! The RESIDUAL STOPPING TEST for the computed solution is performed here
	!	!---------------------------------------------------------------------------
	!	! NOTE: from this point vector B(N) is no longer containing the right-hand
	!	! side of the problem! It contains the current FGMRES approximation to the
	!	! solution. If you need to keep the right-hand side, save it in some other
	!	! vector before the call to DFGMRES ROUTine. Here we saved it in vector
	!	! RHS(N). The vector B is used instead of RHS to preserve the original
	!	! right-hand side of the problem and guarantee the proper restart of FGMRES
	!	! method. Vector B will be altered when computing the residual stopping
	!	! criterion!
	!	!---------------------------------------------------------------------------
	!	IF (RCI_REQUEST .EQ. 2) THEN
	!		! Request to the DFGMRES_GET ROUTine to put the solution into B(N) via IPAR(13)
	!		IPAR(13) = 1 !if 0 updates de solution to vector x, if positive write the solution to B
	!		! Get the current FGMRES solution in the vector B(N)
	!		CALL DFGMRES_GET(NO%NumNP, COMPUTED_SOLUTION, B, RCI_REQUEST, IPAR,DPAR, TMP, ITERCOUNT)
	!		! Compute the current true residual via Intel(R) MKL (Sparse) BLAS ROUTines
	!		info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrAScaled,descrA,B,beta,RESIDUAL)
	!		CALL DAXPY(NO%NumNP, -1.0D0.0, MX%BX, 1, RESIDUAL, 1)
	!		DVAR = DNRM2(NO%NumNP, RESIDUAL, 1)
	!		IF (DVAR .LT. 1.0E-20) THEN !Initially 1.0E-3
	!			GO TO 3 !Go to 3 if Success
	!		ELSE
	!			GO TO 1 !Start new iteration
	!		END IF
	!	END IF
	!	!---------------------------------------------------------------------------
	!	! If RCI_REQUEST=3, then apply the preconditioner on the vector
	!	! TMP(IPAR(22)-1:ipar(22)+n-2) and put the result in vector TMP(IPAR(23)-1:ipar(23)+n-2)
	!	!---------------------------------------------------------------------------
	!	IF (RCI_REQUEST .EQ. 3) THEN
	!		IF (IPAR(4) .EQ. 3) THEN !ipar[4] is current iteration number
	!			TMP(IPAR(23)+0)=-2.0D0 !ipar(23) specifies memory location which second vector (output start)
	!			TMP(IPAR(23)+1)= 0.08519601586107672D0
	!			TMP(IPAR(23)+2)=-1.1590871369607090D0
	!			TMP(IPAR(23)+3)=-0.65791939687456980D0
	!			TMP(IPAR(23)+4)= 0.75660051476696133D0
	!		ELSE
	!			IF(IPAR(4).EQ.4) THEN
	!				TMP(IPAR(23)+0)= 0.0D0
	!				TMP(IPAR(23)+1)= 0.0D0
	!				TMP(IPAR(23)+2)= 0.0D0
	!				TMP(IPAR(23)+3)= 1.0D0
	!				TMP(IPAR(23)+4)=-1.0D0
	!			ELSE
	!				DO I = 0.0, NO%NumNP-1
	!					TMP(IPAR(23)+I) = I * TMP(IPAR(22)+I)
	!				END DO
	!			END IF
	!		END IF
	!		GO TO 1
	!	END IF
	!	!---------------------------------------------------------------------------
	!	! If RCI_REQUEST=4, then check if the norm of the next generated vector is
	!	! not zero up to rounding and computational errors. The norm is contained
	!	! in DPAR(7) parameter
	!	!---------------------------------------------------------------------------
	!	IF (RCI_REQUEST .EQ. 4) THEN
	!		IF (DPAR(7) .LT. 1.0D-12) THEN
	!			GO TO 3
	!		ELSE
	!			GO TO 1
	!		END IF
	!		!---------------------------------------------------------------------------
	!		! If RCI_REQUEST=anything else, then DFGMRES subROUTine failed
	!		! to compute the solution vector: COMPUTED_SOLUTION(N)
	!		!---------------------------------------------------------------------------
	!	ELSE
	!		GO TO 999
	!	END IF
	!	!---------------------------------------------------------------------------
	!	! Reverse Communication ends here
	!	! Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
	!	! call DFGMRES_GET ROUTine as computed_solution is still containing
	!	! the initial guess!). Request to DFGMRES_GET to put the solution into
	!	! vector COMPUTED_SOLUTION(N) via IPAR(13)
	!	!---------------------------------------------------------------------------
	!3	IPAR(13) = 0 !IF 0 updates the solution to vector x, positive to the RHS vector B.
	!	CALL DFGMRES_GET(NO%NumNP, COMPUTED_SOLUTION, MX%BX, RCI_REQUEST, IPAR, DPAR, TMP, ITERCOUNT)
	!	!---------------------------------------------------------------------------
	!	! Release internal Intel(R) MKL memory that might be used for computations
	!	! NOTE: It is important to call the ROUTine below to avoid memory leaks
	!	! unless you disable Intel(R) MKL Memory Manager
	!	!---------------------------------------------------------------------------
	!
	!    info = MKL_SPARSE_DESTROY(csrA)
	!	info = MKL_SPARSE_DESTROY(csrD)
	!	info = MKL_SPARSE_DESTROY(csrAScaled)
	!
	!	CALL MKL_FREE_BUFFERS
	!
	!	!DVAR = DNRM2(N,EXPECTED_SOLUTION,1) !Euclidean Norm
	!
	!	MX%BSol=COMPUTED_SOLUTION*D
	!
	!	GO TO 1000
	!
	!	!---------------------------------------------------------------------------
	!	! Release internal Intel(R) MKL memory that might be used for computations
	!	! NOTE: It is important to call the ROUTine below to avoid memory leaks
	!	! unless you disable Intel(R) MKL Memory Manager
	!	!---------------------------------------------------------------------------
	!999 WRITE( *,'(A,A,I5)') 'This example FAILED as the solver has', &
	!		&  ' returned the ERROR code', RCI_REQUEST
	!	info = MKL_SPARSE_DESTROY(csrD)
	!	info = MKL_SPARSE_DESTROY(csrAScaled)
	!	CALL MKL_FREE_BUFFERS
	!	STOP 1
	!
	!1000 END SUBROUTINE SOLVE_GMRES_SCALED

	end module com_mod_mx_solvers