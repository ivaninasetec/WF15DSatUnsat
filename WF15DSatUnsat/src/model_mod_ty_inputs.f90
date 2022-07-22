	!********************************************************************************************************************
	!        MODULE TO INPUT ALL DATA IN THE WHOLE MODEL AND POLULATE ALL CLASS INSTANCES
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : model_mod_inputs
	! URL           : ...
	! AFFILIATION   : ...
	! DATE          : ...
	! REVISION      : ... V 0.0
	! LICENSE				: This software is copyrighted 2019(C)
	!> @author
	!> Iván Campos-Guereta Díez
	!  MSc Civil Engineering by Polytechnic University of Madrid
	!  PhD Student by University of Nottingham
	!  eMBA by International Institute San Telmo in Seville
	!  ivan.camposguereta@nottingham.ac.uk
	! DESCRIPTION:
	!> Class for horizontal saturated layer. Extend common class of layers.
	!********************************************************************************************************************

	module model_mod_inputs

	implicit none
	include 'inc_precision.fi'

	private

	public::  s_model_inputs_parameters,s_model_inputs_materials,s_model_inputs_layers,s_model_inputs_mesh, s_model_inputs_boundary,s_model_inputs_hinit

	contains
	!---------------------------------------------------------------------------------------------------------------------
	! s_unsat_inputs_parameters(materials)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_materials. Input all data in class from file.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_model_inputs_hinit(nodes,parameters,FileINPUT,tinit)
	use com_mod_ty_nodes, only: ty_com_nodes
	use com_mod_ty_parameters, only: ty_com_parameters
	
	class(ty_com_nodes),intent(inout)::nodes
	class(ty_com_parameters),intent(inout)::parameters
	character*400,intent(in)::fileinput
	real(kind=dpd),intent(in)::tinit
	
	integer::readerror,i

	character*400:: chblock

	real(kind=dpd)::tini,t,dumb
	
	!reading file, block a: parameters...
	open(990, file=fileinput, status='old')
	read (990,'(A)', iostat=readerror) chblock
	
	tini = 0.0_dpd
	
	readerror = 0

	do while (t<tini)
		read (990, iostat=readerror) t
		if (readerror.ne.0) then
			write(*,*) 'Specified time does not exist on previous output file'
			stop
		end if
		parameters%tinit=tini
		tini=t
	end do
	BACKSPACE(990)
	do i=1,nodes%count
		BACKSPACE(990)
	end do

	do i=1,nodes%count
		read(990,*) t,dumb,dumb,nodes%hinit(i)
	end do
		
	close(990)
	
	end subroutine s_model_inputs_hinit
	
	!---------------------------------------------------------------------------------------------------------------------
	! s_unsat_inputs_parameters(materials)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_materials. Input all data in class from file.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_model_inputs_parameters(parameters,FileINPUT)
	use com_mod_ty_parameters, only: ty_com_parameters
	use com_mod_inputs, only:s_com_inputs_locateblock,s_com_inputs_nextrecord
	
	integer,parameter::STORAGE_DENSE=1,STORAGE_CSR=2,STORAGE_BANDED=3
	integer,parameter::PRECONDITIONER_NONE=0,PRECONDITIONER_JACOBI=1,PRECONDITIONER_ILU0=2,PRECONDITIONER_ILUT=3
	integer,parameter::SOLVER_GAUSS=0,SOLVER_DSS=1,SOLVER_PARADISO=2,SOLVER_FMGRES=3
	class(ty_com_parameters),intent(inout)::parameters
	character*400,intent(in)::fileinput
	
	integer::readerror	
	
		!Temporal variables to store inputs...
	real(kind=dps)::param_epsh_tol !tolerance in h
	real(kind=dps)::param_epsth_tol !tolerance in th
	real(kind=dps)::param_epsthhsat_tol
	integer::       param_itinc_dt
	integer::       param_itdec_dt
	integer::       param_itmin
	integer::       param_itmax
	real(kind=dps)::param_crelax
	logical::       param_masslump
	logical::       param_isModifiedPicard
	integer::       param_quadratureorder
	integer::       param_typesolution !0 for dense, 1 csr-dss, 2 csr-fgmres, 3 banded-direct
	integer::       param_typematrixstorage !type of matrix sparsity: 1 for dense, 2 sparse csr, 3 banded
	integer::       param_typepreconditioner !preconditioner for the solver (0 for none, 1 jacobi, 2 ilu0, 3 ilut)
	integer::       param_typesolver !type of solver (0 for gauss, 1 direct dss, 2 direct paradiso, 3 fgmres)

	real(kind=dps)::param_ccourant !coefficient of courant for implicit euler (very high to avoid courant)
	real(kind=dps)::param_mulksal  ! this is a fixed slope (d(h+z)/dx defined at exit
	real(kind=dps)::param_max_hsatinc !Max increment in hsat.

	real(kind=dps)::param_tinit  
	real(kind=dps)::param_dtinit 
	real(kind=dps)::param_tmax   
	real(kind=dps)::param_dtinc  
	real(kind=dps)::param_dtdec  
	real(kind=dps)::param_dtmax  
	real(kind=dps)::param_dtmin  
	real(kind=dps)::param_tprintinc
	logical				::param_isrestarttime 
	real(kind=dps)::param_trestart 
	
	!reading file, block a: parameters...
	open(40, file=fileinput, status='old')
	call s_com_inputs_locateblock(40, 'A', readerror)
	
	!read(40,*)
	!read(40,*)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_epsh_tol, param_epsth_tol, param_epsthhsat_tol
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_itmin,param_itinc_dt, param_itdec_dt,param_itmax
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_crelax
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_masslump, param_isModifiedPicard
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_quadratureorder	
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_typesolution
	!read(40,*)
	!read(40,*) param_typematrixstorage
	!read(40,*)
	!read(40,*) param_typepreconditioner
	!read(40,*)
	!read(40,*) param_typesolver
	!read(40,*)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_ccourant
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_mulksal
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_max_hsatinc
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_tinit, param_dtinit, param_tmax, param_dtmin, param_dtmax
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_dtinc,	param_dtdec
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_tprintinc
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_isrestarttime
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_trestart
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) parameters%nrelmin,parameters%nrelmax
	
	select case (param_typesolution)
		
	case (1) !dense      gauss    none
		param_typematrixstorage = STORAGE_DENSE
		param_typesolver = SOLVER_GAUSS
		param_typepreconditioner = PRECONDITIONER_NONE
	case (2) !dense      gauss    jacobi
		param_typematrixstorage = STORAGE_DENSE
		param_typesolver = SOLVER_GAUSS
		param_typepreconditioner = PRECONDITIONER_JACOBI
			case (3) !csr        DSS      none
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_DSS
		param_typepreconditioner = PRECONDITIONER_NONE
			case (4) !csr        DSS      jacobi
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_DSS
		param_typepreconditioner = PRECONDITIONER_JACOBI
			case (5) !csr        DSS      ilu0
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_DSS
		param_typepreconditioner = PRECONDITIONER_ILU0
			case (6) !csr        DSS      ilut
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_DSS
		param_typepreconditioner = PRECONDITIONER_ILUT
			case (7) !csr        PARADISO none
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_PARADISO
		param_typepreconditioner = PRECONDITIONER_NONE
			case (8) !csr        PARADISO jacobi
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_PARADISO
		param_typepreconditioner = PRECONDITIONER_JACOBI
			case (9) !csr        PARADISO ilu0
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_PARADISO
		param_typepreconditioner = PRECONDITIONER_ILU0
			case (10) !csr        PARADISO ilut
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_PARADISO
		param_typepreconditioner = PRECONDITIONER_ILUT
			case (11) !csr        FGMRES   none
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_FMGRES
		param_typepreconditioner = PRECONDITIONER_NONE
			case (12) !csr        FGMRES   jacobi
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_FMGRES
		param_typepreconditioner = PRECONDITIONER_JACOBI
			case (13) !csr        FGMRES   ilu0
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_FMGRES
		param_typepreconditioner = PRECONDITIONER_ILU0
			case (14) !csr        FGMRES   ilut
		param_typematrixstorage = STORAGE_CSR
		param_typesolver = SOLVER_FMGRES
		param_typepreconditioner = PRECONDITIONER_ILUT
			case (15) !banded     gauss	  none
		param_typematrixstorage = STORAGE_BANDED
		param_typesolver = SOLVER_GAUSS
		param_typepreconditioner = PRECONDITIONER_NONE
			case (16) !banded     gauss	  jacobi
		param_typematrixstorage = STORAGE_BANDED
		param_typesolver = SOLVER_GAUSS
		param_typepreconditioner = PRECONDITIONER_JACOBI
		end select
	
	
	CLOSE(40)
	
	!Assing readings to parameters instance...
	parameters%epsh_tol		= PARAM_EPSH_TOL ! Tolerance in pressure head [L]
	parameters%epsth_tol	= PARAM_EPSTH_TOL ! Tolerance in watercontent [-]
	parameters%epshsat_tol	= param_epsthhsat_tol
	parameters%it_inc_dt	= PARAM_ITINC_DT ! Iterations below which dt is increased
	parameters%it_dec_dt	= PARAM_ITDEC_DT ! Iterations over which dt is decreased
	parameters%it_min			= PARAM_ITMIN
	parameters%it_max			= PARAM_ITMAX ! Max number of iterations, over which time step restart decreased
	parameters%crelax			= PARAM_CRELAX ! Relaxation coefficient in the update pressure head in each iteration
	parameters%masslump		= PARAM_MASSLUMP ! Mass lumping used? (.true. or .false.)
	parameters%isModifiedPicard				= PARAM_ISMODIFIEDPICARD ! Is used the error on node or error on element? (.true.= error on node)
	parameters%quadratureorder		=	PARAM_QUADRATUREORDER ! Quadrature order for integration inside element
	parameters%typesolution				= PARAM_TYPESOLUTION ! Type of matrix solver (0 for dense, 1 csr-dss, 2 csr-fgmres, 3 banded-direct)
	parameters%typematrixstorage	= PARAM_TYPEMATRIXSTORAGE ! Type of matrix sparsity: 1 for dense, 2 sparse csr, 3 banded
	parameters%typepreconditioner = PARAM_TYPEPRECONDITIONER! Preconditioner for the solver (0 for none, 1 jacobi, 2 ilu0, 2 ilut)
	parameters%typesolver					= PARAM_TYPESOLVER! Type of solver (0 for gauss, 1 direct dss, 2 direct paradiso, 3 fgmres)
	parameters%ccourant						= PARAM_CCOURANT ! Courant number
	parameters%multksal						= PARAM_MULKSAL ! Permeability at the seepage surface compared to permeability on element
	parameters%maxhsatinc					= param_max_hsatinc
	
	
	
	
	parameters%tinit = param_tinit
	parameters%dtinit = param_dtinit
	parameters%tmax = param_tmax
	parameters%dtinc = param_dtinc
	parameters%dtdec = param_dtdec
	parameters%dtmax = param_dtmax
	parameters%dtmin = param_dtmin
	parameters%tprintinc = param_tprintinc
	parameters%isrestarttime			= param_isrestarttime
	parameters%trestart			= param_trestart
	
	end subroutine s_model_inputs_parameters
		
	!---------------------------------------------------------------------------------------------------------------------
	! s_sat_inputs_materials(materials)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_materials. Input all data in class from file.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_model_inputs_materials(material,fileinput)
	use com_mod_ty_material, only: ty_com_material
	use com_mod_inputs, only:s_com_inputs_locateblock,s_com_inputs_nextrecord

	type(ty_com_material),allocatable,intent(inout)::material(:)
	character*400,intent(in)::fileinput
	
	integer::readerror,i

	!temporal parameters to store inputs...
	integer:: param_materials_count
	integer,allocatable::				param_kindmat(:)
	real(kind=dps),allocatable::param_thsat(:)
	real(kind=dps),allocatable::param_thres(:)
	real(kind=dps),allocatable::param_ksat(:)
	real(kind=dps),allocatable::param_a(:)
	real(kind=dps),allocatable::param_n(:)
	real(kind=dps),allocatable::param_m(:)
	real(kind=dps),allocatable::param_l(:)

	!reading file, block b: materials...
	open(40, file=fileinput, status='old')
	call s_com_inputs_locateblock(40, 'B', readerror)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_materials_count
	
	allocate(param_kindmat(param_materials_count))
	allocate(param_thsat(param_materials_count))
	allocate(param_thres(param_materials_count))
	allocate(param_ksat(param_materials_count))	
	allocate(param_a(param_materials_count))
	allocate(param_n(param_materials_count))
	allocate(param_m(param_materials_count))
	allocate(param_l(param_materials_count))

	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) (param_kindmat(i),i=1,param_materials_count)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) (param_thsat(i),i=1,param_materials_count)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) (param_thres(i),i=1,param_materials_count)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) (param_ksat(i),i=1,param_materials_count)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) (param_a(i),i=1,param_materials_count)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) (param_n(i),i=1,param_materials_count)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) (param_m(i),i=1,param_materials_count)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) (param_l(i),i=1,param_materials_count)
	
	close(40)

	
	!Assing values to materials instance...
	if(.not.allocated(material)) allocate(material(param_materials_count))
	!call materials%allocateall(PARAM_MATERIALS_COUNT)

	do i=1, param_materials_count
		material(i)%id = i
	end do

	material(:)%kindmat	= PARAM_KINDMAT
	material(:)%thsat		= PARAM_THSAT
	material(:)%thres		= PARAM_THRES
	material(:)%ksat		= PARAM_KSAT
	material(:)%a				= PARAM_A
	material(:)%n				= PARAM_N
	material(:)%m				= PARAM_M
	material(:)%l				= PARAM_L
	
	deallocate(PARAM_KINDMAT)
	deallocate(PARAM_THSAT)
	deallocate(PARAM_Thres)
	deallocate(PARAM_Ksat)	
	deallocate(PARAM_a)
	deallocate(PARAM_n)
	deallocate(PARAM_m)
	deallocate(PARAM_l)

	end subroutine s_model_inputs_materials	
	
	!---------------------------------------------------------------------------------------------------------------------
	! s_model_inputs_layers(layers)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Build the instance layers from the data in file.
	!> @param[in] fileinput
	!> @param[in] materials
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_model_inputs_layers(layers,fileinput,material)
	use com_mod_ty_layers, only: ty_com_layers
	use com_mod_ty_material,only:ty_com_material
	use com_mod_inputs, only:s_com_inputs_locateblock,s_com_inputs_nextrecord
		
	class(ty_com_layers),intent(inout)::layers
	character*400,intent(in)::fileinput
	type(ty_com_material),intent(inout)::material(:)
	integer::i,readerror

	!temporal variables to store inputs...
	integer::				 param_mod_layers_count	
	real(kind=dps):: param_mod_layers_width

	real(kind=dps),allocatable::param_mod_layer_h(:)
	real(kind=dps),allocatable::param_mod_layer_slope(:)
	integer,allocatable::				param_mod_layer_nmaterial(:)

	!Start reading file...
	open(40, file=fileinput, status='old')
	call s_com_inputs_locateblock(40, 'C', readerror)

	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) param_mod_layers_count !READ: count

	allocate(param_mod_layer_h(param_mod_layers_count))
	allocate(param_mod_layer_slope(param_mod_layers_count+1))
	allocate(param_mod_layer_nmaterial(param_mod_layers_count))
	
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) param_mod_layers_width	!READ: width

	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) (param_mod_layer_h(i),i=1,param_mod_layers_count) !READ: height(:)
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) (param_mod_layer_slope(i),i=1,param_mod_layers_count+1) !READ: slope(:)
	call s_com_inputs_nextrecord(40,readerror) 
	READ(40,*) (param_mod_layer_nmaterial(i),i=1,param_mod_layers_count) !READ: material(:)
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) layers%zphr	!READ, ASSIGN: zphr
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) layers%topboundbyh !READ, ASSIGN: layers%topboundbyh
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) layers%topboundbyq	!READ, ASSIGN: layers%topboundbyq
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) layers%bottombyphl	!READ, ASSIGN: layers%topboundbyq
	close(40)
	
	!Build layers instance with readings...
	
	layers%count = param_mod_layers_count !ASSIGN: count
	call layers%allocateall(param_mod_layers_count)
	layers%width		= param_mod_layers_width	!ASSIGN: width
	layers%height		= param_mod_layer_h !ASSIGN: height(:)
	layers%slope		= param_mod_layer_slope !ASSIGN: slope(:)

	do i=1, PARAM_MOD_LAYERS_COUNT
		layers%material(i)		= material(param_mod_layer_nmaterial(i)) !ASSIGN: material(:)
	end do

	layers%hbottom(1) = 0.0
	layers%htop(1) = layers%height(1)
	do i=2,PARAM_MOD_LAYERS_COUNT
		layers%hbottom(i) = layers%htop(i-1)	!CALC: htop(:)
		layers%htop(i)		= layers%hbottom(i)+layers%height(i) !CALC: hbottom(:)
	end do

	deallocate(param_mod_layer_h)
	deallocate(param_mod_layer_slope)
	deallocate(param_mod_layer_nmaterial)
	
	end subroutine s_model_inputs_layers
	
	!---------------------------------------------------------------------------------------------------------------------
	! s_model_inputs_mesh(mesh)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Build the instance layers from the data in file.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_model_inputs_mesh(mesh,fileinput,layers)
	use model_mod_ty_mesh, only: ty_model_mesh
	use com_mod_ty_layers, only: ty_com_layers
	use com_mod_inputs, only:s_com_inputs_locateblock,s_com_inputs_nextrecord
	
	class(ty_model_mesh),intent(inout)::mesh
	character*400,intent(in)::fileinput
	class(ty_com_layers),intent(in),target::layers
	
	integer::i,e,idnode,readerror
	real(kind=dps)::incx
	
	!temporal variables to store inputs...
	integer::PARAM_MESH_NODES_PER_ELEMENT
	integer::PARAM_MESH_CLASS_OF_NODE
	integer:: PARAM_MOD_NMODVER
	integer:: PARAM_VMOD_COUNT
	integer,ALLOCATABLE:: PARAM_NELEMH(:),PARAM_NELEMV(:)
	real(kind=dps),allocatable::PARAM_X_VMOD(:)
	
	!Start reading file...
	open(40, file=fileinput, status='old')
	call s_com_inputs_locateblock(40, 'D', readerror)
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) PARAM_MESH_NODES_PER_ELEMENT, PARAM_MESH_CLASS_OF_NODE	!READ: nn, nc
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) PARAM_VMOD_COUNT !READ: vmod_count
	
	!mesh%layers_count => layers%count !ASSIGN: layers_count
	
	allocate(PARAM_NELEMH(PARAM_VMOD_COUNT),PARAM_X_VMOD(PARAM_VMOD_COUNT),PARAM_NELEMV(layers%count))
	
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) (PARAM_NELEMH(i),i=1,PARAM_VMOD_COUNT) !READ: nelemh(:)
	
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) (PARAM_X_VMOD(i),i=1,PARAM_VMOD_COUNT) !READ: vmod_x(:)	
	
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) (PARAM_NELEMV(i),i=1,layers%count) !READ: nelemv(:)

	
	!Build mesh instance from inputs...
	
	!MESH------
	mesh%nelemv_count		= sum(PARAM_NELEMV)					!CALC: nelemv_count	
	mesh%nelemh_count		= sum(PARAM_NELEMH)					!CALC: nelemh_count	

	mesh%nn					= PARAM_MESH_NODES_PER_ELEMENT	!ASSIGN:	nn
	mesh%nc					= PARAM_MESH_CLASS_OF_NODE			!ASSIGN:	nc
	mesh%nnc				= mesh%nn*(mesh%nnc+1)					!CALC: nnc
	mesh%nnodh_count			= mesh%nelemh_count*(mesh%nn-1)+1			!CALC: nnodh_count (horizontal)
	mesh%nnodclassh_count	= mesh%nnodh_count*(mesh%nc+1)				!CALC: nnodclassh_count (horizontal)
	mesh%nnodv_count			= mesh%nelemv_count*(mesh%nn-1)+1			!CALC: nnodv_count (vertical)
	mesh%nnodclassv_count	= mesh%nnodv_count*(mesh%nc+1)				!CALC: nnodclassv_count (vertical)

	mesh%height =sum(layers%height)									!CALC:	height
	mesh%width	=> layers%width											!ASSIGN:width

	call mesh%allocateall(layers%count,PARAM_VMOD_COUNT)

	mesh%nelemv = PARAM_NELEMV											!ASSIGN:nelemv(:)
	mesh%nelemh = PARAM_NELEMH											!ASSIGN:nelemh(:)
	mesh%vmod_x = PARAM_X_VMOD											!ASSIGN:vmod_x(:)
	
	
	mesh%vmod_idnod(1)=1 !First vertical module on node 1
	do i=2, PARAM_VMOD_COUNT
		mesh%vmod_idnod(i)=mesh%vmod_idnod(i-1)+PARAM_NELEMH(i)	!CALC: vmod_idnod
	end do	
	
	deallocate(PARAM_NELEMH)
	deallocate(PARAM_NELEMV)
	deallocate(PARAM_X_VMOD)
	
	end subroutine s_model_inputs_mesh
	
	!---------------------------------------------------------------------------------------------------------------------
	!> @author iván campos-guereta díez
	!> @brief
	!> This subroutine reads boundary condition file and fill the class "ty_com_boundary".
	!> @param[in] ne
	!> @param[in] nn
	!> @param[in] nc
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_model_inputs_boundary(fileboundary,boundary,layers)
	use com_mod_ty_boundary, only:ty_com_boundary
	use com_mod_ty_layers, only:ty_com_layers
	
	character*400,intent(in)::fileboundary
	class(ty_com_boundary),intent(inout)::boundary
	class(ty_com_layers),intent(in)::layers
	
	integer:: readerror,i
	real(kind=dpd)::dumb
	
	!Read data from file: flow1d_bound_h  (top point dirichlet boundary condition) *****
	boundary%topboundbyh = layers%topboundbyh

	if (boundary%topboundbyh ) then
		open (55,  file=fileboundary, status='old')

		boundary%timestepscount=0
		readerror=0
		read(55,'(a)')
		do while (readerror==0)
			read (40,*, iostat=readerror) dumb,dumb
			boundary%timestepscount=boundary%timestepscount+1
			!write(*,*) mo%bd%timestepscount
		end do
		boundary%timestepscount=boundary%timestepscount-1

		!write(*,*) 'finish read bound_h'
		write(*,*) 'allocating boundary array...'
		allocate(boundary%timebound(boundary%timestepscount))
		allocate(boundary%hbound(boundary%timestepscount))

		rewind(55)
		read(55,'(a)')
		!write(*,*) 'rewinded'
		!write(*,*) 'timestep count: ', mo%bd%timestepscount
		do i=1, boundary%timestepscount
			read (40,*) boundary%timebound(i), boundary%hbound(i)
			!write(*,*) i
		end do

		close(55)

	end if
	
	!Read data from file: flow1d_bound_q  (top point dirichlet boundary condition) *****
	boundary%TopBoundbyQ = layers%TopBoundbyq
	
	if (boundary%topboundbyq) then
		open (55, file=fileboundary, status='old')

		boundary%timestepscount=0
		readerror=0
		read(55,'(a)')
		do while (readerror==0)
			read (55,*, iostat=readerror) dumb,dumb
			boundary%timestepscount=boundary%timestepscount+1
		end do
		boundary%timestepscount=boundary%timestepscount-1
		!boundary%timestepscount=boundary%timestepscount-2

		write(*,*) 'allocating boundary array q...'
		if (allocated(boundary%timebound)) deallocate(boundary%timebound)
		allocate(boundary%timebound(boundary%timestepscount))
		allocate(boundary%qbound(boundary%timestepscount))

		rewind(55)
		read(55,'(a)')
		do i=1, boundary%timestepscount
			read (55,*) boundary%timebound(i), boundary%qbound(i)
		end do

		close(55)

	end if
	
	end subroutine s_model_inputs_boundary	
	
	
	
	end module model_mod_inputs