	!------------------------------------------------------------------------------
	!        MOD_COM_DATATYPES_MATRIX
	!------------------------------------------------------------------------------
	! TITLE         : project name
	! PROJECT       : 1.5D MULTILAYER FLOW MODEL
	! MODULE        : MOD_COM_DATATYPES_MATRIX
	! URL           : ...
	! AFFILIATION   : ...
	! DATE          : ...
	! REVISION      : ... V 0.15
	! @author
	! Iván Campos-Guereta Díez
	!
	! DESCRIPTION:
	! Module to include matrix datatype and operations.
	!------------------------------------------------------------------------------

	!********************************************************************************************************************
	! COM_MOD_MX_DATATYPES
	!********************************************************************************************************************
	! TITLE         : LIBRARY OF COMMON FUNCTIONS AND CLASSES TO USE IN FLOW1DSAT, FLOW1DUNSAT AND FLOW15DSATUNSAT.
	! PROJECT       : FLOW15DSATUNSAT
	! MODULE        : COM_MOD_MX_DATATYPES
	! URL           : ...
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2020
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2019(C)
	!
	!> <B>MODULE TO INCLUDE MATRIX DATATYPE AND OPERATIONS</B>
	!> This module define an abstract class \b ty_mx, that need to be subclassed in one of three possible classes:
	!> <UL><B>\link ty_mx_dense \endlink:</B> Class for dense matrix.</UL>
	!> <UL><B>\link ty_mx_CSR \endlink:</B> Class for sparse CSR matrix.</UL>
	!> <UL><B>\link ty_mx_dense \endlink:</B> Class for sparse banded matrix.</UL>
	!>
	!> An example code to create a class:
	!> \code{.f90}
	!> program test_mx_datatypes
	!> use test_com_mod_mx_datatypes,only:test_mod_com_basis
	!>
	!> call test_mod_com_basis()
	!> end program
	!>
	!> module test_com_mod_mx_datatypes
	!>
	!> contains
	!> subroutine test_mod_com_basis
	!> use com_mod_mx_datatypes, only:ty_mx,ty_mx_dense,ty_mx_csr,ty_mx_banded
	!>
	!> class(ty_mx),allocatable::mxsample1,mxsample2,mxresult
	!> class(ty_mx),pointer::mx1,mx2,mx3
	!> integer::numnp,ncount
	!>
	!> allocate(ty_mx_csr::mxsample1,ty_mx_csr::mxsample2,ty_mx_csr::mxsample3) !In this case ty_mx_csr is dynamically allocated as CSR
	!> mx1 => mxsample1
	!> numnp = 3	!Number of nodes
	!> ncount = 5 !Number of non zero values
	!> kl = 1 !Number of lower diagonals
	!> ku = 1 !Number of upper diagonals
	!>
	!> select type(mx1)
	!>  type is(ty_mx_dense)
	!>   mx1%init(numnp,numnp)
	!>  type is(ty_mx_csr)
	!>   mx1%init(ncount,numnp)
	!>  type is(ty_mx_csr)
	!>   call mx1%init(numnp,kl,ku)
	!> end select
	!>
	!> end subroutine test_mod_com_basis
	!> end module test_com_mod_mx_datatypes
	!> \endcode
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************


	module com_mod_mx_datatypes

	implicit none
	include 'inc_precision.fi'

	private

	!abstract class---------
	type,abstract,public :: ty_mx
		real(kind=dpd),public,allocatable::mx(:,:)
		real(kind=dpd),public,allocatable::mxval(:)
		integer,public,allocatable::csrrows(:),csrcols(:)

	contains
	private

	procedure(genmulmv),pass(this),deferred::mul_mv
	procedure(genmul_msca4),pass(this),deferred::mul_msca4
	procedure(genmul_msca8),pass(this),deferred::mul_msca8
	procedure(genmul_msca16),pass(this),deferred::mul_msca16
	procedure(genmul_sca4m),pass(this),deferred::mul_sca4m
	procedure(genmul_sca8m),pass(this),deferred::mul_sca8m
	procedure(genmul_sca16m),pass(this),deferred::mul_sca16m
	generic, public :: operator(*) => mul_mv,mul_msca4,mul_msca8,mul_msca16,mul_sca4m,mul_sca8m,mul_sca16m

	procedure(gendiv_msca4),pass(this),deferred::div_msca4
	procedure(gendiv_msca8),pass(this),deferred::div_msca8
	procedure(gendiv_msca16),pass(this),deferred::div_msca16
	generic, public :: operator(/) => div_msca4,div_msca8,div_msca16

	procedure(gensum_mm), pass(this),deferred::sum_mm
	generic, public :: operator(+) => sum_mm

	procedure(gensub_mm), pass(this),deferred::sub_mm
	generic, public :: operator(-) => sub_mm

	procedure(geninit),public,pass(this),deferred:: init

	procedure(genget),public,pass(this),deferred:: get
	procedure(genset),public,pass(this),deferred:: set

	procedure(genassign_matrix_to_array), pass(this),deferred::assign_matrix_to_array
	procedure(genassign_array_to_matrix), pass(this),deferred::assign_array_to_matrix
	procedure(genassign_array_to_array), pass(this),deferred::assign_array_to_array
	procedure(genassign_array_to_sca4), pass(this),deferred::assign_array_to_sca4
	procedure(genassign_array_to_sca8), pass(this),deferred::assign_array_to_sca8
	procedure(genassign_array_to_sca16), pass(this),deferred::assign_array_to_sca16
	generic, public :: assignment(=) => assign_matrix_to_array, assign_array_to_matrix,assign_array_to_array,assign_array_to_sca4,assign_array_to_sca8,assign_array_to_sca16

	end type

	!abstract interfaces -----------
	abstract interface
	subroutine geninit(this,int1,int2,int3)
	import::ty_mx,dpd
	class(ty_mx),intent(inout)::this
	integer,intent(in)::int1
	integer,intent(in)::int2
	integer,intent(in),optional::int3
	end subroutine geninit

	function genmulmv(this,v)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real(kind=dpd),intent(in)::v(:)
	real(kind=dpd)::genmulmv(size(v))
	end function genmulmv

	function genmul_msca4(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*4,intent(in)::sca
	class(ty_mx),allocatable::genmul_msca4
	end function genmul_msca4

	function genmul_msca8(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*8,intent(in)::sca
	class(ty_mx),allocatable::genmul_msca8
	end function genmul_msca8

	function genmul_msca16(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*16,intent(in)::sca
	class(ty_mx),allocatable::genmul_msca16
	end function genmul_msca16

	function genmul_sca4m(sca,this)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*4,intent(in)::sca
	class(ty_mx),allocatable::genmul_sca4m
	end function genmul_sca4m

	function genmul_sca8m(sca,this)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*8,intent(in)::sca
	class(ty_mx),allocatable::genmul_sca8m
	end function genmul_sca8m

	function genmul_sca16m(sca,this)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*16,intent(in)::sca
	class(ty_mx),allocatable::genmul_sca16m
	end function genmul_sca16m

	function gendiv_msca4(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*4,intent(in)::sca
	class(ty_mx),allocatable::gendiv_msca4
	end function gendiv_msca4

	function gendiv_msca8(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*8,intent(in)::sca
	class(ty_mx),allocatable::gendiv_msca8
	end function gendiv_msca8

	function gendiv_msca16(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real*16,intent(in)::sca
	class(ty_mx),allocatable::gendiv_msca16
	end function gendiv_msca16

	function gensum_mm(this,mx)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::gensum_mm
	end function gensum_mm

	function gensub_mm(this,mx)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::gensub_mm
	end function gensub_mm

	subroutine genassign_array_to_matrix(this,mx)
	import::ty_mx,dpd
	class(ty_mx),intent(inout)::this
	real(kind=dpd),intent(in)::mx(:,:)
	end subroutine genassign_array_to_matrix

	subroutine genassign_matrix_to_array(mx,this)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	real(kind=dpd),intent(inout)::mx(:,:)
	end subroutine genassign_matrix_to_array

	subroutine genassign_array_to_array(this,mx)
	import::ty_mx,dpd
	class(ty_mx),intent(inout)::this
	class(ty_mx),intent(in)::mx
	end subroutine genassign_array_to_array

	subroutine genassign_array_to_sca4(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(inout)::this
	real(kind=4),intent(in)::sca
	end subroutine genassign_array_to_sca4

	subroutine genassign_array_to_sca8(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(inout)::this
	real(kind=8),intent(in)::sca
	end subroutine genassign_array_to_sca8

	subroutine genassign_array_to_sca16(this,sca)
	import::ty_mx,dpd
	class(ty_mx),intent(inout)::this
	real(kind=16),intent(in)::sca
	end subroutine genassign_array_to_sca16

	function genget(this,nrow,ncol)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	integer,intent(in)::nrow
	integer,intent(in)::ncol
	real(kind=dpd)::genget
	end function genget

	subroutine genset(this,nrow,ncol,val)
	import::ty_mx,dpd
	class(ty_mx),intent(inout)::this
	integer,intent(in)::nrow
	integer,intent(in)::ncol
	real(kind=dpd),intent(in)::val
	end subroutine genset

	function mul_mm(this,mx)
	import::ty_mx,dpd
	class(ty_mx),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::mul_mm
	end function mul_mm
	end interface

	!---------------------------------------------------------------------------------------------------------------------
	! TY_MX_DENSE
	!---------------------------------------------------------------------------------------------------------------------

	type,extends(ty_mx),public::ty_mx_dense

	contains
	procedure,pass(this):: mul_mv => mulmv_dns
	procedure,public,pass(this):: init => init_dense
	procedure,pass(this):: mul_msca4 => mul_msca4_dns
	procedure,pass(this):: mul_msca8 => mul_msca8_dns
	procedure,pass(this):: mul_msca16 => mul_msca16_dns
	procedure,pass(this):: mul_sca4m => mul_sca4m_dns
	procedure,pass(this):: mul_sca8m => mul_sca8m_dns
	procedure,pass(this):: mul_sca16m => mul_sca16m_dns
	procedure,pass(this):: div_msca4 => div_msca4_dns
	procedure,pass(this):: div_msca8 => div_msca8_dns
	procedure,pass(this):: div_msca16 => div_msca16_dns
	procedure,pass(this):: sum_mm => sum_mm_dns
	procedure,pass(this):: sub_mm => sub_mm_dns
	procedure, pass(this)::assign_matrix_to_array => assign_matrix_to_array_dns
	procedure, pass(this)::assign_array_to_matrix => assign_array_to_matrix_dns
	procedure, pass(this)::assign_array_to_array  => assign_array_to_array_dns
	procedure, pass(this)::assign_array_to_sca4   => assign_array_to_sca4_dns
	procedure, pass(this)::assign_array_to_sca8   => assign_array_to_sca8_dns
	procedure, pass(this)::assign_array_to_sca16  => assign_array_to_sca16_dns
	procedure, public,pass(this):: get => get_dns
	procedure, public,pass(this):: set => set_dns


	!PROCEDURES IMPLEMENTED ONLY IN DENSE MATRIX
	procedure,public,pass(this)::mul_mm => mul_mm_dns

	end type

	type,extends(ty_mx),public::ty_mx_banded
		integer,public::kl
		integer,public::ku
	contains
	procedure,pass(this):: mul_mv => mulmv_bnd
	procedure,public,pass(this):: init => init_banded
	procedure,pass(this):: mul_msca4 => mul_msca4_bnd
	procedure,pass(this):: mul_msca8 => mul_msca8_bnd
	procedure,pass(this):: mul_msca16 => mul_msca16_bnd
	procedure,pass(this):: mul_sca4m => mul_sca4m_bnd
	procedure,pass(this):: mul_sca8m => mul_sca8m_bnd
	procedure,pass(this):: mul_sca16m => mul_sca16m_bnd
	procedure,pass(this):: div_msca4 => div_msca4_bnd
	procedure,pass(this):: div_msca8 => div_msca8_bnd
	procedure,pass(this):: div_msca16 => div_msca16_bnd
	procedure,pass(this):: sum_mm => sum_mm_bnd
	procedure,pass(this):: sub_mm => sub_mm_bnd
	procedure, pass(this)::assign_matrix_to_array => assign_matrix_to_array_bnd
	procedure, pass(this)::assign_array_to_matrix => assign_array_to_matrix_bnd
	procedure, pass(this)::assign_array_to_array  => assign_array_to_array_bnd
	procedure, pass(this)::assign_array_to_sca4   => assign_array_to_sca4_bnd
	procedure, pass(this)::assign_array_to_sca8   => assign_array_to_sca8_bnd
	procedure, pass(this)::assign_array_to_sca16  => assign_array_to_sca16_bnd
	procedure, public,pass(this):: get => get_bnd
	procedure, public,pass(this):: set => set_bnd

	end type

	type,extends(ty_mx),public::ty_mx_csr
		integer,public,allocatable::POS(:,:)
	contains
	procedure,pass(this):: mul_mv => mulmv_csr
	procedure,public,pass(this):: init => init_csr
	procedure,pass(this):: mul_msca4 => mul_msca4_csr
	procedure,pass(this):: mul_msca8 => mul_msca8_csr
	procedure,pass(this):: mul_msca16 => mul_msca16_csr
	procedure,pass(this):: mul_sca4m => mul_sca4m_csr
	procedure,pass(this):: mul_sca8m => mul_sca8m_csr
	procedure,pass(this):: mul_sca16m => mul_sca16m_csr
	procedure,pass(this):: div_msca4 => div_msca4_csr
	procedure,pass(this):: div_msca8 => div_msca8_csr
	procedure,pass(this):: div_msca16 => div_msca16_csr
	procedure,pass(this):: sum_mm => sum_mm_csr
	procedure,pass(this):: sub_mm => sub_mm_csr
	procedure, pass(this)::assign_matrix_to_array => assign_matrix_to_array_csr
	procedure, pass(this)::assign_array_to_matrix => assign_array_to_matrix_csr
	procedure, pass(this)::assign_array_to_array  => assign_array_to_array_csr
	procedure, pass(this)::assign_array_to_sca4   => assign_array_to_sca4_csr
	procedure, pass(this)::assign_array_to_sca8   => assign_array_to_sca8_csr
	procedure, pass(this)::assign_array_to_sca16  => assign_array_to_sca16_csr
	procedure, public,pass(this):: get => get_csr
	procedure, public,pass(this):: set => set_csr

	end type

	contains

	!********************************************************************************************************************
	! F: PERMUTE ROWS (DENSE)
	!********************************************************************************************************************

	function permute_rows_dns(this,p)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"permute_rows_dns" :: permute_rows_dns
	!DEC$ endif


	class(ty_mx_dense),intent(in)::this
	integer,intent(in)::p(:)
	class(ty_mx),allocatable::permute_rows_dns

	allocate(ty_mx_dense::permute_rows_dns)

	select type (permute_rows_dns)
	type is (ty_mx_dense)
		call permute_rows_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
		permute_rows_dns%mx = this%mx(p,:)
	end select

	end function permute_rows_dns

	!********************************************************************************************************************
	! F: PERMUTE COLUMNS (DENSE)
	!********************************************************************************************************************

	function permute_cols_dns(this,p)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"permute_cols_dns" :: permute_cols_dns
	!DEC$ endif


	class(ty_mx_dense),intent(in)::this
	integer,intent(in)::p(:)
	class(ty_mx),allocatable::permute_cols_dns

	allocate(ty_mx_dense::permute_cols_dns)

	select type (permute_cols_dns)
	type is (ty_mx_dense)
		call permute_cols_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
		permute_cols_dns%mx = this%mx(:,p)
	end select

	end function permute_cols_dns

	!********************************************************************************************************************
	! F: PERMUTE ROWS AND COLUMNS (DENSE)
	!********************************************************************************************************************

	function permute_rows_cols_dns(this,p)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"permute_rows_cols_dns" :: permute_rows_cols_dns
	!DEC$ endif


	class(ty_mx_dense),intent(in)::this
	integer,intent(in)::p(:)
	class(ty_mx),allocatable::permute_rows_cols_dns

	allocate(ty_mx_dense::permute_rows_cols_dns)

	select type (permute_rows_cols_dns)
	type is (ty_mx_dense)
		call permute_rows_cols_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
		permute_rows_cols_dns%mx = this%mx(p,:)
		permute_rows_cols_dns%mx = permute_rows_cols_dns%mx(:,p)
	end select

	end function permute_rows_cols_dns

	!********************************************************************************************************************
	! F: MATRIX VECTOR MULTIPLICATION FUNCTIONS
	!********************************************************************************************************************

	function mulmv_dns(this,v) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mulmv_dns" :: mulmv_dns
	!DEC$ endif


	class(ty_mx_dense),intent(in)::this
	real(kind=dpd),intent(in)::v(:)
	real(kind=dpd)::rout(size(v))

	!allocate(rout(size(this%mx,1)))
	rout = matmul(this%mx,v)

	end function mulmv_dns

	!------------------------------------------------------------

	function mulmv_bnd(this,v) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mulmv_bnd" :: mulmv_bnd
	!DEC$ endif


	real(kind=dpd),parameter::alpha = 1.0_dpd, beta = 0.0_dpd
	class(ty_mx_banded),intent(in)::this
	real(kind=dpd),intent(in)::v(:)
	real(kind=dpd)::rout(size(v))
	integer::numnp

	numnp = size(this%mx,2)
	call dgbmv('N',numnp,numnp, this%kl, this%ku, alpha,this%mx, 1+this%kl+this%ku, v,1,beta,rout,1)

	end function mulmv_bnd

	!------------------------------------------------------------

	function mulmv_csr(this,v) result(rout)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mulmv_csr" :: mulmv_csr
	!DEC$ endif


	real(kind=dpd),parameter::alpha = 1.0_dpd, beta = 0.0_dpd
	class(ty_mx_csr),intent(in)::this
	real(kind=dpd),intent(in)::v(:)
	real(kind=dpd)::rout(size(v))
	integer::numnp,i,n1,n2

	numnp = size(this%csrrows,1)-1

	!$OMP PARALLEL SHARED(rout) PRIVATE(n1,n2,i)
	!$OMP DO
	do i=1, numnp
		n1 = this%csrrows(i)
		n2 = this%csrrows(i+1)-1
		rout(i) = dot_product(this%mxval(n1:n2),v(this%csrcols(n1:n2)))
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	end function mulmv_csr


	function mul_msca4_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca4_dns" :: mul_msca4_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::mul_msca4_dns

	allocate(ty_mx_dense::mul_msca4_dns)
	call mul_msca4_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (mul_msca4_dns)
	type is (ty_mx_dense)
		mul_msca4_dns%mx = sca*this%mx
	end select


	end function mul_msca4_dns


	function mul_msca4_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca4_bnd" :: mul_msca4_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::mul_msca4_bnd

	allocate(ty_mx_banded::mul_msca4_bnd)
	call mul_msca4_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (mul_msca4_bnd)
	type is (ty_mx_banded)
		mul_msca4_bnd%mx = sca*this%mx
	end select

	end function mul_msca4_bnd

	function mul_msca4_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca4_csr" :: mul_msca4_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::mul_msca4_csr

	allocate(ty_mx_csr::mul_msca4_csr)
	call mul_msca4_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (mul_msca4_csr)
	type is (ty_mx_csr)
		mul_msca4_csr%mxval = sca*this%mxval
	end select

	end function mul_msca4_csr

	!-------------------
	function mul_msca8_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca8_dns" :: mul_msca8_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::mul_msca8_dns

	allocate(ty_mx_dense::mul_msca8_dns)
	call mul_msca8_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (mul_msca8_dns)
	type is (ty_mx_dense)
		mul_msca8_dns%mx = sca*this%mx
	end select


	end function mul_msca8_dns


	function mul_msca8_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca8_bnd" :: mul_msca8_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::mul_msca8_bnd

	allocate(ty_mx_banded::mul_msca8_bnd)
	call mul_msca8_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (mul_msca8_bnd)
	type is (ty_mx_banded)
		mul_msca8_bnd%mx = sca*this%mx
	end select

	end function mul_msca8_bnd

	function mul_msca8_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca8_csr" :: mul_msca8_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::mul_msca8_csr

	allocate(ty_mx_csr::mul_msca8_csr)
	call mul_msca8_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (mul_msca8_csr)
	type is (ty_mx_csr)
		mul_msca8_csr%mxval = sca*this%mxval
	end select

	end function mul_msca8_csr
	!
	!	!-------------------
	!
	function mul_msca16_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca16_dns" :: mul_msca16_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::mul_msca16_dns

	allocate(ty_mx_dense::mul_msca16_dns)
	call mul_msca16_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (mul_msca16_dns)
	type is (ty_mx_dense)
		mul_msca16_dns%mx = sca*this%mx
	end select


	end function mul_msca16_dns


	function mul_msca16_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca16_bnd" :: mul_msca16_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::mul_msca16_bnd

	allocate(ty_mx_banded::mul_msca16_bnd)
	call mul_msca16_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (mul_msca16_bnd)
	type is (ty_mx_banded)
		mul_msca16_bnd%mx = sca*this%mx
	end select

	end function mul_msca16_bnd

	function mul_msca16_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_msca16_csr" :: mul_msca16_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::mul_msca16_csr

	allocate(ty_mx_csr::mul_msca16_csr)
	call mul_msca16_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (mul_msca16_csr)
	type is (ty_mx_csr)
		mul_msca16_csr%mxval = sca*this%mxval
	end select

	end function mul_msca16_csr

	!--------------------------------------------------------------------
	function mul_sca4m_dns(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca4m_dns" :: mul_sca4m_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::mul_sca4m_dns

	allocate(ty_mx_dense::mul_sca4m_dns)
	call mul_sca4m_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (mul_sca4m_dns)
	type is (ty_mx_dense)
		mul_sca4m_dns%mx = sca*this%mx
	end select

	end function mul_sca4m_dns


	function mul_sca4m_bnd(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca4m_bnd" :: mul_sca4m_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::mul_sca4m_bnd

	allocate(ty_mx_banded::mul_sca4m_bnd)
	call mul_sca4m_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (mul_sca4m_bnd)
	type is (ty_mx_banded)
		mul_sca4m_bnd%mx = sca*this%mx
	end select

	end function mul_sca4m_bnd

	function mul_sca4m_csr(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca4m_csr" :: mul_sca4m_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::mul_sca4m_csr

	allocate(ty_mx_csr::mul_sca4m_csr)
	call mul_sca4m_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (mul_sca4m_csr)
	type is (ty_mx_csr)
		mul_sca4m_csr%mxval = sca*this%mxval
	end select

	end function mul_sca4m_csr

	!	!-------------------
	function mul_sca8m_dns(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca8m_dns" :: mul_sca8m_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::mul_sca8m_dns

	allocate(ty_mx_dense::mul_sca8m_dns)
	call mul_sca8m_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (mul_sca8m_dns)
	type is (ty_mx_dense)
		mul_sca8m_dns%mx = sca*this%mx
	end select


	end function mul_sca8m_dns


	function mul_sca8m_bnd(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca8m_bnd" :: mul_sca8m_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::mul_sca8m_bnd

	allocate(ty_mx_banded::mul_sca8m_bnd)
	call mul_sca8m_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (mul_sca8m_bnd)
	type is (ty_mx_banded)
		mul_sca8m_bnd%mx = sca*this%mx
	end select

	end function mul_sca8m_bnd

	function mul_sca8m_csr(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca8m_csr" :: mul_sca8m_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::mul_sca8m_csr

	allocate(ty_mx_csr::mul_sca8m_csr)
	call mul_sca8m_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (mul_sca8m_csr)
	type is (ty_mx_csr)
		mul_sca8m_csr%mxval = sca*this%mxval
	end select

	end function mul_sca8m_csr
	!
	!	!-------------------
	!
	function mul_sca16m_dns(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca16m_dns" :: mul_sca16m_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::mul_sca16m_dns

	allocate(ty_mx_dense::mul_sca16m_dns)
	call mul_sca16m_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (mul_sca16m_dns)
	type is (ty_mx_dense)
		mul_sca16m_dns%mx = sca*this%mx
	end select


	end function mul_sca16m_dns


	function mul_sca16m_bnd(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca16m_bnd" :: mul_sca16m_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::mul_sca16m_bnd

	allocate(ty_mx_banded::mul_sca16m_bnd)
	call mul_sca16m_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (mul_sca16m_bnd)
	type is (ty_mx_banded)
		mul_sca16m_bnd%mx = sca*this%mx
	end select

	end function mul_sca16m_bnd

	function mul_sca16m_csr(sca,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_sca16m_csr" :: mul_sca16m_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::mul_sca16m_csr

	allocate(ty_mx_csr::mul_sca16m_csr)
	call mul_sca16m_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (mul_sca16m_csr)
	type is (ty_mx_csr)
		mul_sca16m_csr%mxval = sca*this%mxval
	end select

	end function mul_sca16m_csr

	! !********************************************************************************************************************
	!	! f: division matrix by scalar
	!	!********************************************************************************************************************


	function div_msca4_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca4_dns" :: div_msca4_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::div_msca4_dns

	allocate(ty_mx_dense::div_msca4_dns)
	call div_msca4_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (div_msca4_dns)
	type is (ty_mx_dense)
		div_msca4_dns%mx = this%mx/sca
	end select


	end function div_msca4_dns


	function div_msca4_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca4_bnd" :: div_msca4_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::div_msca4_bnd

	allocate(ty_mx_banded::div_msca4_bnd)
	call div_msca4_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (div_msca4_bnd)
	type is (ty_mx_banded)
		div_msca4_bnd%mx = this%mx/sca
	end select

	end function div_msca4_bnd

	function div_msca4_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca4_csr" :: div_msca4_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=4),intent(in)::sca
	class(ty_mx),allocatable::div_msca4_csr

	allocate(ty_mx_csr::div_msca4_csr)
	call div_msca4_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (div_msca4_csr)
	type is (ty_mx_csr)
		div_msca4_csr%mxval = this%mxval/sca
	end select

	end function div_msca4_csr

	!	!-------------------
	function div_msca8_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca8_dns" :: div_msca8_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::div_msca8_dns

	allocate(ty_mx_dense::div_msca8_dns)
	call div_msca8_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (div_msca8_dns)
	type is (ty_mx_dense)
		div_msca8_dns%mx = this%mx/sca
	end select


	end function div_msca8_dns


	function div_msca8_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca8_bnd" :: div_msca8_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::div_msca8_bnd

	allocate(ty_mx_banded::div_msca8_bnd)
	call div_msca8_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (div_msca8_bnd)
	type is (ty_mx_banded)
		div_msca8_bnd%mx = this%mx/sca
	end select

	end function div_msca8_bnd

	function div_msca8_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca8_csr" :: div_msca8_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=8),intent(in)::sca
	class(ty_mx),allocatable::div_msca8_csr

	allocate(ty_mx_csr::div_msca8_csr)
	call div_msca8_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (div_msca8_csr)
	type is (ty_mx_csr)
		div_msca8_csr%mxval = this%mxval/sca
	end select

	end function div_msca8_csr
	!
	!	!-------------------
	!
	function div_msca16_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca16_dns" :: div_msca16_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::div_msca16_dns

	allocate(ty_mx_dense::div_msca16_dns)
	call div_msca16_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
	select type (div_msca16_dns)
	type is (ty_mx_dense)
		div_msca16_dns%mx = this%mx/sca
	end select


	end function div_msca16_dns


	function div_msca16_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca16_bnd" :: div_msca16_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::div_msca16_bnd

	allocate(ty_mx_banded::div_msca16_bnd)
	call div_msca16_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor

	select type (div_msca16_bnd)
	type is (ty_mx_banded)
		div_msca16_bnd%mx = this%mx/sca
	end select

	end function div_msca16_bnd

	function div_msca16_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"div_msca16_csr" :: div_msca16_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=16),intent(in)::sca
	class(ty_mx),allocatable::div_msca16_csr

	allocate(ty_mx_csr::div_msca16_csr)
	call div_msca16_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor

	select type (div_msca16_csr)
	type is (ty_mx_csr)
		div_msca16_csr%mxval = this%mxval/sca
	end select

	end function div_msca16_csr

	! !********************************************************************************************************************
	!	! F: MULTIPLICATION OF TWO MATRIX (DENSE)
	!	!********************************************************************************************************************

	function mul_mm_dns(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"mul_mm_dns" :: mul_mm_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::mul_mm_dns

	select type (mx)
	type is (ty_mx_dense)
		allocate(ty_mx_dense::mul_mm_dns)
		!allocate(sum_mm,source=ty_mx_dense())
		call mul_mm_dns%init(size(this%mx,1),size(mx%mx,2)) !constructor
		select type(mul_mm_dns)
		type is (ty_mx_dense)
			mul_mm_dns%mx = matmul(this%mx,mx%mx)
		end select
	class default
		stop('error: multiply between different storages kinds')
	end select

	end function mul_mm_dns


	! !********************************************************************************************************************
	!	! F: TRANSPOSE OF A MATRIX (DENSE)
	!	!********************************************************************************************************************

	function transpose_dns(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"transpose_dns" :: transpose_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	class(ty_mx),allocatable::transpose_dns

	allocate(ty_mx_dense::transpose_dns)

	select type(transpose_dns)
	type is (ty_mx_dense)
		call transpose_dns%init(size(this%mx,2),size(this%mx,1)) !constructor
		transpose_dns%mx = transpose(this%mx)
	class default
		stop('error: multiply between different storages kinds')
	end select

	end function transpose_dns








	!  !********************************************************************************************************************
	!	! F: TWO MATRIX SUM OF THE SAME KIND AND SAME SPARSITY
	!	!********************************************************************************************************************
	! !


	function sum_mm_dns(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"sum_mm_dns" :: sum_mm_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::sum_mm_dns

	select type (mx)
	type is (ty_mx_dense)
		allocate(ty_mx_dense::sum_mm_dns)
		!allocate(sum_mm,source=ty_mx_dense())
		call sum_mm_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
		select type(sum_mm_dns)
		type is (ty_mx_dense)
			sum_mm_dns%mx = this%mx+mx%mx
		end select
	class default
		stop('error: sum between different storages kinds')
	end select

	end function sum_mm_dns



	function sum_mm_bnd(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"sum_mm_bnd" :: sum_mm_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::sum_mm_bnd

	select type (mx)
	type is (ty_mx_banded)
		allocate(ty_mx_banded::sum_mm_bnd)
		!allocate(sum_mm,source=ty_mx_banded(kl=this%kl,ku=this%ku))
		call sum_mm_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor
		select type(sum_mm_bnd)
		type is (ty_mx_banded)
			sum_mm_bnd%mx = this%mx+mx%mx
		end select
	class default
		stop('error: sum between different storages kinds')
	end select

	end function sum_mm_bnd



	function sum_mm_csr(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"sum_mm_csr" :: sum_mm_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::sum_mm_csr

	select type (mx)
	type is (ty_mx_csr)
		allocate(ty_mx_csr::sum_mm_csr)
		!allocate(sum_mm,source=ty_mx_csr())
		call sum_mm_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor
		select type(sum_mm_csr)
		type is (ty_mx_csr)
			sum_mm_csr%mxval = this%mxval+mx%mxval
		end select
	class default
		stop('error: sum between different storages kinds')
	end select

	end function sum_mm_csr


	!
	!  !********************************************************************************************************************
	!	! F: TWO MATRIX SUBTRACT OF THE SAME KIND AND SAME SPARSITY
	!	!********************************************************************************************************************
	! !
	function sub_mm_dns(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"sub_mm_dns" :: sub_mm_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::sub_mm_dns

	select type (mx)
	type is (ty_mx_dense)
		allocate(ty_mx_dense::sub_mm_dns)
		!allocate(sub_mm,source=ty_mx_dense())
		call sub_mm_dns%init(size(this%mx,1),size(this%mx,2)) !constructor
		select type(sub_mm_dns)
		type is (ty_mx_dense)
			sub_mm_dns%mx = this%mx-mx%mx
		end select
	class default
		stop('error: sub between different storages kinds')
	end select

	end function sub_mm_dns

	function sub_mm_bnd(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"sub_mm_bnd" :: sub_mm_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::sub_mm_bnd

	select type (mx)
	type is (ty_mx_banded)
		allocate(ty_mx_banded::sub_mm_bnd)
		!allocate(sub_mm,source=ty_mx_banded(kl=this%kl,ku=this%ku))
		call sub_mm_bnd%init(size(this%mx,2),this%kl,this%ku) !constructor
		select type(sub_mm_bnd)
		type is (ty_mx_banded)
			sub_mm_bnd%mx = this%mx-mx%mx
		end select
	class default
		stop('error: sub between different storages kinds')
	end select

	end function sub_mm_bnd


	function sub_mm_csr(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"sub_mm_csr" :: sub_mm_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	class(ty_mx),intent(in)::mx
	class(ty_mx),allocatable::sub_mm_csr

	select type (mx)
	type is (ty_mx_csr)
		allocate(ty_mx_csr::sub_mm_csr)
		!allocate(sub_mm,source=ty_mx_csr())
		call sub_mm_csr%init(size(this%mxval),size(this%csrrows)-1) !constructor
		select type(sub_mm_csr)
		type is (ty_mx_csr)
			sub_mm_csr%mxval = this%mxval-mx%mxval
		end select
	class default
		stop('error: sub between different storages kinds')
	end select

	end function sub_mm_csr


	!
	!********************************************************************************************************************
	! s: MATRIXCONSTRUCTORS
	!********************************************************************************************************************

	subroutine init_dense(this,int1,int2,int3)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"init_dense" :: init_dense
	!DEC$ endif

	class(ty_mx_dense),intent(inout)::this
	integer,intent(in)::int1 !nrows
	integer,intent(in)::int2 !ncols
	integer,intent(in),optional::int3 !not used

	if(.not.allocated(this%mx)) allocate(this%mx(int1,int2))

	end subroutine init_dense

	!-----------------------------------------------------------------

	subroutine init_banded(this,int1,int2,int3)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"init_banded" :: init_banded
	!DEC$ endif

	class(ty_mx_banded),intent(inout)::this
	integer,intent(in)::int1 !ncols
	integer,intent(in)::int2 !kl
	integer,intent(in),optional::int3 !ku
	integer::ku,kl,ncols

	ncols = int1
	kl = int2
	if(present(int3)) then
		ku=int3
	else
		ku=int2
	end if


	if(.not.allocated(this%mx)) allocate(this%mx(1+kl+ku,int1))
	this%kl=kl
	this%ku=ku

	end subroutine init_banded

	!-----------------------------------------------------------------

	subroutine init_csr(this,int1,int2,int3)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"init_csr" :: init_csr
	!DEC$ endif

	class(ty_mx_csr),intent(inout)::this
	integer,intent(in)::int1 !nval
	integer,intent(in)::int2 !nrows
	integer,intent(in),optional::int3 !not used

	if(.not.allocated(this%mxval)) allocate(this%mxval(int1))

	if(.not.allocated(this%csrcols)) allocate(this%csrcols(int1))

	if(.not.allocated(this%csrrows)) allocate(this%csrrows(int2+1))

	if(.not.allocated(this%pos)) allocate(this%pos(int2,int2))

	end subroutine init_csr

	!
	!  !********************************************************************************************************************
	!	! s: ASSIGNMENT OVERRIDE
	!	!********************************************************************************************************************

	! MATRIX CLASS TO ARRAY

	subroutine assign_matrix_to_array_dns(mx,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_matrix_to_array_dns" :: assign_matrix_to_array_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	real(kind=dpd),intent(inout)::mx(:,:)



	mx = this%mx

	end subroutine assign_matrix_to_array_dns


	subroutine assign_matrix_to_array_bnd(mx,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_matrix_to_array_bnd" :: assign_matrix_to_array_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	real(kind=dpd),intent(inout)::mx(:,:)

	integer::ncol,i

	mx = 0.0_dpd
	!it only fills the sparse values
	!$OMP PARALLEL SHARED(mx) PRIVATE(ncol,i)
	!$OMP DO
	do ncol=1,size(this%mx,2)
		do i=-this%kl,this%ku
			if (ncol+i>=1.and.ncol+i<=1) then
				mx(ncol+i,ncol)=this%mx(i,ncol)
			else
				mx(i,ncol)=0.0_dpd
			end if
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL


	end subroutine assign_matrix_to_array_bnd

	subroutine assign_matrix_to_array_csr(mx,this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_matrix_to_array_csr" :: assign_matrix_to_array_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	real(kind=dpd),intent(inout)::mx(:,:)

	integer::nrow,ncol,i

	mx =  0.0_dpd
	!$OMP PARALLEL SHARED(mx) PRIVATE(nrow,ncol,i)
	!$OMP DO
	do nrow=1,size(this%csrrows)-1
		do i=this%csrrows(nrow),this%csrrows(nrow+1)-1
			ncol = this%csrcols(i)
			mx(nrow,ncol)=this%mxval(i)
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	end subroutine assign_matrix_to_array_csr






	! ARRAY TO MATRIX CLASS

	subroutine assign_array_to_matrix_dns(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_matrix_dns" :: assign_array_to_matrix_dns
	!DEC$ endif

	class(ty_mx_dense),intent(inout)::this
	real(kind=dpd),intent(in)::mx(:,:)



	this%mx=mx

	end subroutine assign_array_to_matrix_dns


	subroutine assign_array_to_matrix_bnd(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_matrix_bnd" :: assign_array_to_matrix_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(inout)::this
	real(kind=dpd),intent(in)::mx(:,:)

	integer::ncol,i

	!it only fills the sparse values
	!$OMP PARALLEL SHARED(this) PRIVATE(ncol,i)
	!$OMP DO
	do ncol=1,size(this%mx,2)
		do i=-this%kl,this%ku
			if (ncol+i>=1.and.ncol+i<=1) then
				this%mx(i,ncol)=mx(ncol+i,ncol)
			else
				this%mx(i,ncol)=0.0_dpd
			end if
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	end subroutine assign_array_to_matrix_bnd

	subroutine assign_array_to_matrix_csr(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_matrix_csr" :: assign_array_to_matrix_csr
	!DEC$ endif

	class(ty_mx_csr),intent(inout)::this
	real(kind=dpd),intent(in)::mx(:,:)

	integer::nrow,ncol,i

	!$OMP PARALLEL SHARED(this) PRIVATE(nrow,ncol,i)
	!$OMP DO
	do nrow=1,size(this%csrrows)-1
		do i=this%csrrows(nrow),this%csrrows(nrow+1)-1
			ncol = this%csrcols(i)
			this%mxval(i)=mx(nrow,ncol)
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	end subroutine assign_array_to_matrix_csr


	!MATRIX CLASS TO MATRIX CLASS....

	subroutine assign_array_to_array_dns(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_array_dns" :: assign_array_to_array_dns
	!DEC$ endif

	class(ty_mx_dense),intent(inout)::this
	class(ty_mx),intent(in)::mx

	select type (mx)
	type is (ty_mx_dense)
		this%mx=mx%mx
	class default
		stop('error: assignment between different storages kinds')
	end select

	end subroutine assign_array_to_array_dns

	subroutine assign_array_to_array_bnd(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_array_bnd" :: assign_array_to_array_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(inout)::this
	class(ty_mx),intent(in)::mx

	select type (mx)
	type is (ty_mx_banded)
		this%mx=mx%mx
	class default
		stop('error: assignment between different storages kinds')
	end select

	end subroutine assign_array_to_array_bnd


	subroutine assign_array_to_array_csr(this,mx)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_array_csr" :: assign_array_to_array_csr
	!DEC$ endif

	class(ty_mx_csr),intent(inout)::this
	class(ty_mx),intent(in)::mx

	select type (mx)
	type is (ty_mx_csr)
		this%mxval=mx%mxval
	class default
		stop('error: assignment between different storages kinds')
	end select

	end subroutine assign_array_to_array_csr


	!SCALAR TO MATRIX CLASS....
	!
	!  SUBROUTINE assign_array_to_sca4(this,sca)
	! CLASS(ty_mx),INTENT(INOUT)::this
	! REAL(KIND=4),INTENT(IN)::sca
	!
	! SELECT TYPE (this)
	! TYPE IS (ty_mx_dense)
	! this%mx=sca
	! TYPE IS (ty_mx_banded)
	! this%mx=sca
	! TYPE IS (ty_mx_csr)
	! this%mxval=sca
	! END SELECT
	!  END SUBROUTINE assign_array_to_sca4

	subroutine assign_array_to_sca4_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca4_dns" :: assign_array_to_sca4_dns
	!DEC$ endif

	class(ty_mx_dense),intent(inout)::this
	real(kind=4),intent(in)::sca

	this%mx=sca
	end subroutine assign_array_to_sca4_dns

	subroutine assign_array_to_sca4_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca4_bnd" :: assign_array_to_sca4_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(inout)::this
	real(kind=4),intent(in)::sca

	this%mx=sca
	end subroutine assign_array_to_sca4_bnd

	subroutine assign_array_to_sca4_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca4_csr" :: assign_array_to_sca4_csr
	!DEC$ endif

	class(ty_mx_csr),intent(inout)::this
	real(kind=4),intent(in)::sca

	this%mxval=sca
	end subroutine assign_array_to_sca4_csr
	! !----------------------------------------------------------------------
	!
	subroutine assign_array_to_sca8_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca8_dns" :: assign_array_to_sca8_dns
	!DEC$ endif

	class(ty_mx_dense),intent(inout)::this
	real(kind=8),intent(in)::sca

	this%mx=sca
	end subroutine assign_array_to_sca8_dns

	subroutine assign_array_to_sca8_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca8_bnd" :: assign_array_to_sca8_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(inout)::this
	real(kind=8),intent(in)::sca

	this%mx=sca
	end subroutine assign_array_to_sca8_bnd

	subroutine assign_array_to_sca8_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca8_csr" :: assign_array_to_sca8_csr
	!DEC$ endif

	class(ty_mx_csr),intent(inout)::this
	real(kind=8),intent(in)::sca

	this%mxval=sca
	end subroutine assign_array_to_sca8_csr
	!
	!   !----------------------------------------------------------------------
	!
	subroutine assign_array_to_sca16_dns(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca16_dns" :: assign_array_to_sca16_dns
	!DEC$ endif

	class(ty_mx_dense),intent(inout)::this
	real(kind=16),intent(in)::sca

	this%mx=sca
	end subroutine assign_array_to_sca16_dns

	subroutine assign_array_to_sca16_bnd(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca16_bnd" :: assign_array_to_sca16_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(inout)::this
	real(kind=16),intent(in)::sca

	this%mx=sca
	end subroutine assign_array_to_sca16_bnd

	subroutine assign_array_to_sca16_csr(this,sca)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"assign_array_to_sca16_csr" :: assign_array_to_sca16_csr
	!DEC$ endif

	class(ty_mx_csr),intent(inout)::this
	real(kind=16),intent(in)::sca

	this%mxval=sca
	end subroutine assign_array_to_sca16_csr


	!********************************************************************************************************************
	! s: GET VALUES AT GIVEN ROW AND COLUMN IN THE MATRIX
	!********************************************************************************************************************

	function get_dns(this,nrow,ncol)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"get_dns" :: get_dns
	!DEC$ endif

	class(ty_mx_dense),intent(in)::this
	integer,intent(in)::nrow
	integer,intent(in)::ncol
	real(kind=dpd)::get_dns



	get_dns = this%mx(nrow,ncol)

	end function get_dns

	!--------------------------------------
	function get_bnd(this,nrow,ncol)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"get_bnd" :: get_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(in)::this
	integer,intent(in)::nrow
	integer,intent(in)::ncol
	real(kind=dpd)::get_bnd



	if ((1+this%ku+nrow-ncol).ge.1.and.(1+this%ku+nrow-ncol).le.(1+this%ku+this%kl)) then
		get_bnd = this%mx(1+this%ku+nrow-ncol,ncol)
	else
		get_bnd = 0.0_dpd
	end if

	end function get_bnd

	!-------------------------
	function get_csr(this,nrow,ncol)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"get_csr" :: get_csr
	!DEC$ endif

	class(ty_mx_csr),intent(in)::this
	integer,intent(in)::nrow
	integer,intent(in)::ncol
	real(kind=dpd)::get_csr



	if(this%pos(nrow,ncol).ne.0) then
		get_csr = this%mxval(this%pos(nrow,ncol))
	else
		get_csr = 0.0_dpd
	end if


	!npos = FINDLOC(this%csrcols(this%csrrows(nrow):this%csrrows(nrow+1)-1),ncol)
	!if (npos(1).NE.0) then
	!  npos(1) = this%csrrows(nrow) + npos(1)-1
	!  get_csr = this%mxval(npos(1))
	!else
	!  get_csr = 0.0_dpd
	!end if

	end function get_csr
	!
	!********************************************************************************************************************
	! s: SET VALUES AT GIVEN ROW AND COLUMN IN THE MATRIX
	!********************************************************************************************************************
	!
	!  SUBROUTINE set_generic(this,nrow,ncol,val)
	!  CLASS(ty_mx),INTENT(INOUT)::this
	!  INTEGER,INTENT(IN)::nrow
	!  INTEGER,INTENT(IN)::ncol
	!  REAL(dpd),intent(in)::val
	!
	!  INTEGER::npos(1),ku,kl
	!
	!  SELECT TYPE (this)
	!  type IS (ty_mx_dense)
	!      this%mx(nrow,ncol)=val
	!	type IS (ty_mx_banded)
	!		ku = this%ku
	!		kl = this%kl
	!    if ((1+this%ku+nrow-ncol).GE.1.AND.(1+this%ku+nrow-ncol).LE.(1+this%ku+this%kl)) THEN
	!      this%mx(1+this%ku+nrow-ncol,ncol)=val
	!    ELSE
	!      STOP('TRIED TO FILL OUTSIDE SPARSITY OF THE MATRIX')
	!    END IF
	!  type IS (ty_mx_csr)
	!    npos = FINDLOC(this%csrcols(this%csrrows(nrow):this%csrrows(nrow+1)-1),ncol)
	!    if (npos(1).NE.0) then
	!      npos(1) = this%csrrows(nrow) + npos(1)-1
	!      this%mxval(npos(1)) = val
	!    else
	!      STOP('TRIED TO FILL OUTSIDE SPARSITY OF THE MATRIX')
	!    end if
	!  END SELECT
	!
	!  END SUBROUTINE set_generic

	subroutine set_dns(this,nrow,ncol,val)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"set_dns" :: set_dns
	!DEC$ endif

	class(ty_mx_dense),intent(inout)::this
	integer,intent(in)::nrow
	integer,intent(in)::ncol
	real(dpd),intent(in)::val



	this%mx(nrow,ncol)=val

	end subroutine set_dns


	subroutine set_bnd(this,nrow,ncol,val)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"set_bnd" :: set_bnd
	!DEC$ endif

	class(ty_mx_banded),intent(inout)::this
	integer,intent(in)::nrow
	integer,intent(in)::ncol
	real(dpd),intent(in)::val

	integer::ku,kl

	ku = this%ku
	kl = this%kl
	if ((1+this%ku+nrow-ncol).ge.1.and.(1+this%ku+nrow-ncol).le.(1+this%ku+this%kl)) then
		this%mx(1+this%ku+nrow-ncol,ncol)=val
	else
		stop('tried to fill outside sparsity of the matrix')
	end if

	end subroutine set_bnd


	subroutine set_csr(this,nrow,ncol,val)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"set_csr" :: set_csr
	!DEC$ endif

	class(ty_mx_csr),intent(inout)::this
	integer,intent(in)::nrow
	integer,intent(in)::ncol
	real(dpd),intent(in)::val



	this%mxval(this%pos(nrow,ncol)) = val
	!npos = findloc(this%csrcols(this%csrrows(nrow):this%csrrows(nrow+1)-1),ncol)
	!if (npos(1).ne.0) then
	!  npos(1) = this%csrrows(nrow) + npos(1)-1
	!  this%mxval(npos(1)) = val
	!else
	!  stop('tried to fill outside sparsity of the matrix')
	!end if

	end subroutine set_csr

	end module com_mod_mx_datatypes