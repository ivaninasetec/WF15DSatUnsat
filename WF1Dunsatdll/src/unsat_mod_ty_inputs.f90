	!********************************************************************************************************************
	!        MODULE TO INPUT ALL DATA IN THE SATURATED MODEL AND POLULATE ALL CLASS INSTANCES
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : MOD_SAT_TY_LAYER
	! URL           : ...
	! AFFILIATION   : ...
	! DATE          : ...
	! REVISION      : ... V 0.
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

	module unsat_mod_inputs

	implicit none
	include 'inc_precision.fi'

	private

	public::  s_unsat_inputs_parameters,s_unsat_inputs_materials,s_unsat_inputs_layers,s_unsat_inputs_mesh,s_unsat_inputs_elements,s_unsat_inputs_nodes,s_unsat_inputs_boundary

	contains

	!---------------------------------------------------------------------------------------------------------------------
	! s_unsat_inputs_parameters(materials)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_materials. Input all data in class from file.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_inputs_parameters(parameters,FileINPUT)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_inputs_parameters" :: s_unsat_inputs_parameters
	!DEC$ endif
	use com_mod_ty_parameters, only: ty_com_parameters
	!use com_mod_inputs, only:s_com_inputs_locateblock,s_com_inputs_nextrecord

	integer,parameter::STORAGE_DENSE=1,STORAGE_CSR=2,STORAGE_BANDED=3
	integer,parameter::PRECONDITIONER_NONE=0,PRECONDITIONER_JACOBI=1,PRECONDITIONER_ILU0=2,PRECONDITIONER_ILUT=3
	integer,parameter::SOLVER_GAUSS=0,SOLVER_DSS=1,SOLVER_PARADISO=2,SOLVER_FMGRES=3
	class(ty_com_parameters),intent(inout)::parameters
	character*400,intent(in)::fileinput

	integer::readerror

	!Temporal variables to store inputs...
	real(kind=dps)::param_epsh_tol !tolerance in h
	real(kind=dps)::param_epsth_tol !tolerance in th
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

	real(kind=dps)::param_tinit
	real(kind=dps)::param_dtinit
	real(kind=dps)::param_tmax
	real(kind=dps)::param_dtinc
	real(kind=dps)::param_dtdec
	real(kind=dps)::param_dtmax
	real(kind=dps)::param_dtmin
	real(kind=dps)::param_tprintinc

	!reading file, block a: parameters...
	open(40, file=fileinput, status='old')
	call s_com_inputs_locateblock(40, 'A', readerror)

	!read(40,*)
	!read(40,*)
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_epsh_tol, param_epsth_tol
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
	read(40,*) param_tinit, param_dtinit, param_tmax, param_dtmin, param_dtmax
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_dtinc,	param_dtdec
	call s_com_inputs_nextrecord(40,readerror)
	read(40,*) param_tprintinc

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
	parameters%it_inc_dt	= PARAM_ITINC_DT ! Iterations below which dt is increased
	parameters%it_dec_dt	= PARAM_ITDEC_DT ! Iterations over which dt is decreased
	parameters%it_min			= PARAM_ITMIN
	parameters%it_max			= PARAM_ITMAX ! Max number of iterations, over which time step restart decreased
	parameters%crelax			= PARAM_CRELAX ! Relaxation coefficient in the update pressure head in each iteration
	parameters%masslump		= PARAM_MASSLUMP ! Mass lumping used? (.true. or .false.)
	parameters%isModifiedPicard				= PARAM_isModifiedPicard ! Is used the error on node or error on element? (.true.= error on node)
	parameters%quadratureorder		=	PARAM_QUADRATUREORDER ! Quadrature order for integration inside element
	parameters%typesolution				= PARAM_TYPESOLUTION ! Type of matrix solver (0 for dense, 1 csr-dss, 2 csr-fgmres, 3 banded-direct)
	parameters%typematrixstorage	= PARAM_TYPEMATRIXSTORAGE ! Type of matrix sparsity: 1 for dense, 2 sparse csr, 3 banded
	parameters%typepreconditioner = PARAM_TYPEPRECONDITIONER! Preconditioner for the solver (0 for none, 1 jacobi, 2 ilu0, 2 ilut)
	parameters%typesolver					= PARAM_TYPESOLVER! Type of solver (0 for gauss, 1 direct dss, 2 direct paradiso, 3 fgmres)
	parameters%ccourant						= PARAM_CCOURANT ! Courant number
	parameters%multksal						= PARAM_MULKSAL ! Permeability at the seepage surface compared to permeability on element


	parameters%tinit = param_tinit
	parameters%dtinit = param_dtinit
	parameters%tmax = param_tmax
	parameters%dtinc = param_dtinc
	parameters%dtdec = param_dtdec
	parameters%dtmax = param_dtmax
	parameters%dtmin = param_dtmin
	parameters%tprintinc = param_tprintinc


	end subroutine s_unsat_inputs_parameters


	!---------------------------------------------------------------------------------------------------------------------
	! s_sat_inputs_materials(materials)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_materials. Input all data in class from file.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_inputs_materials(material,fileinput)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_inputs_materials" :: s_unsat_inputs_materials
	!DEC$ endif
	use com_mod_ty_material, only: ty_com_material
	!use com_mod_inputs, only:s_com_inputs_locateblock,s_com_inputs_nextrecord

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

	end subroutine s_unsat_inputs_materials

	!---------------------------------------------------------------------------------------------------------------------
	! s_sat_inputs_layers(layers)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Build the instance layers from the data in file.
	!> @param[in] fileinput
	!> @param[in] materials
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_inputs_layers(layers,fileinput,material)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_inputs_layers" :: s_unsat_inputs_layers
	!DEC$ endif
	use com_mod_ty_layers, only: ty_com_layers
	use com_mod_ty_material,only:ty_com_material
	!use com_mod_inputs, only:s_com_inputs_locateblock,s_com_inputs_nextrecord

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
	READ(40,*) param_mod_layers_count

	allocate(param_mod_layer_h(param_mod_layers_count))
	allocate(param_mod_layer_slope(param_mod_layers_count))
	allocate(param_mod_layer_nmaterial(param_mod_layers_count))

	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) param_mod_layers_width

	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) (param_mod_layer_h(i),i=1,param_mod_layers_count)
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) (param_mod_layer_nmaterial(i),i=1,param_mod_layers_count)
	!call s_com_inputs_nextrecord(40,readerror)
	!READ(40,*) (param_mod_layer_slope(i),i=1,param_mod_layers_count)
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) layers%zphr
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) layers%topboundbyh
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) layers%topboundbyq
	close(40)

	!Build layers instance with readings...

	layers%count = param_mod_layers_count
	call layers%allocateall(param_mod_layers_count)
	layers%width		= param_mod_layers_width
	layers%height		= param_mod_layer_h
	!layers%slope		= param_mod_layer_slope

	do i=1, PARAM_MOD_LAYERS_COUNT
		layers%material(i)		= material(param_mod_layer_nmaterial(i))
	end do

	!htop and hbottom need a review in order to use different slopes.
	layers%hbottom(1) = 0.0
	layers%htop(1) = layers%height(1)
	do i=2,PARAM_MOD_LAYERS_COUNT
		layers%hbottom(i) = layers%htop(i-1)
		layers%htop(i)		= layers%hbottom(i)+layers%height(i)
	end do

	deallocate(param_mod_layer_h)
	deallocate(param_mod_layer_slope)
	deallocate(param_mod_layer_nmaterial)

	end subroutine s_unsat_inputs_layers

	!---------------------------------------------------------------------------------------------------------------------
	! s_sat_inputs_mesh(mesh)
	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Build the instance layers from the data in file.
	!> @param[in] nmat
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_inputs_mesh(mesh,fileinput,layers)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_inputs_mesh" :: s_unsat_inputs_mesh
	!DEC$ endif
	use unsat_mod_ty_mesh, only: ty_unsat_mesh
	use com_mod_ty_layers, only: ty_com_layers
	!use com_mod_inputs, only:s_com_inputs_locateblock,s_com_inputs_nextrecord

	class(ty_unsat_mesh),intent(inout)::mesh
	character*400,intent(in)::fileinput
	class(ty_com_layers),intent(in)::layers

	integer::i,e,idnode,readerror
	real(kind=dps)::incx

	!temporal variables to store inputs...
	integer::PARAM_MESH_NODES_PER_ELEMENT
	integer::PARAM_MESH_CLASS_OF_NODE
	integer:: PARAM_NELEMENTS
	integer:: PARAM_MOD_NMODVER

	!real(kind=dps),parameter:: total_length = 10.0
	!real(kind=dps),parameter:: mesh_slope = 0.03

	integer:: PARAM_NUMELEM_LAYER(layers%count)
	!integer,allocatable:: NUMELEM_VMOD(:)
	!real(kind=dps),allocatable::PARAM_X_VMOD(:) !checkthis: the boundaries definition will be in an specific function
	!real(kind=dps),allocatable::PARAM_QENTVMOD(:) !checkthis: the boundaries definition will be in an specific function


	!Start reading file...
	open(40, file=fileinput, status='old')
	call s_com_inputs_locateblock(40, 'D', readerror)
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) (PARAM_NUMELEM_LAYER(i),i=1,layers%count)
	call s_com_inputs_nextrecord(40,readerror)
	READ(40,*) PARAM_MESH_NODES_PER_ELEMENT, PARAM_MESH_CLASS_OF_NODE
	!call s_com_inputs_nextrecord(40,readerror)
	!READ(40,*) PARAM_MOD_NMODVER
	!
	!allocate(NUMELEM_VMOD(PARAM_MOD_NMODVER))
	!allocate(PARAM_X_VMOD(PARAM_MOD_NMODVER))
	!allocate(PARAM_QENTVMOD(PARAM_MOD_NMODVER))
	!
	!call s_com_inputs_nextrecord(40,readerror)
	!read(40,*) (NUMELEM_VMOD(i),i=1,param_mod_nmodver)
	!call s_com_inputs_nextrecord(40,readerror)
	!read(40,*) (PARAM_X_VMOD(i),i=1,param_mod_nmodver)
	!call s_com_inputs_nextrecord(40,readerror)
	!read(40,*) (PARAM_QENTVMOD(i),i=1,param_mod_nmodver)
	!close(40)


	!Build mesh instance from inputs...

	!MESH------
	mesh%nelemv_count				= sum(PARAM_NUMELEM_LAYER)
	mesh%nelemv		= PARAM_NUMELEM_LAYER

	mesh%nn					= PARAM_MESH_NODES_PER_ELEMENT
	mesh%nc					= PARAM_MESH_CLASS_OF_NODE
	mesh%nnc				= mesh%nn*(mesh%nnc+1)
	mesh%nnodv_count				= mesh%nelemv_count	*(mesh%nn-1)+1
	mesh%nnodclassv_count	= mesh%nnodv_count*(mesh%nc+1)

	mesh%height =sum(layers%height)
	!mesh%length = layers%width
	!mesh%slope	= layers%slope(1)

	call mesh%allocateall(layers%count)

	!mesh%nelemlayer = PARAM_NUMELEM_LAYER
	!
	!mesh%vmod_nelem = NUMELEM_VMOD
	!mesh%vmod_idnod(1)=1
	!do i=2, PARAM_MOD_NMODVER
	!	mesh%vmod_idnod(i)=mesh%vmod_idnod(i-1)+NUMELEM_VMOD(i)
	!end do
	!
	!mesh%vmod_x=PARAM_X_VMOD
	!
	!mesh%vmod_hsat =0.0_dpd
	!mesh%vmod_hsattemp =0.0_dpd
	!mesh%vmod_hsatold =0.0_dpd
	!
	!mesh%vmod_qent			= PARAM_QENTVMOD
	!mesh%vmod_qenttemp	= PARAM_QENTVMOD
	!mesh%vmod_qentold		= PARAM_QENTVMOD

	!deallocate(NUMELEM_VMOD)
	!deallocate(PARAM_X_VMOD)
	!deallocate(PARAM_QENTVMOD)

	end subroutine s_unsat_inputs_mesh

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Fill all data needed for the ty_unsat_nodes from the mesh data.
	!> @param[in] mesh
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_inputs_nodes(nodes,mesh,layers)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_inputs_nodes" :: s_unsat_inputs_nodes
	!DEC$ endif
	use com_mod_ty_nodes,only:ty_com_nodes
	use unsat_mod_ty_nodes,only:ty_unsat_nodes
	use unsat_mod_ty_mesh,only:ty_unsat_mesh
	use com_mod_ty_layers,only:ty_com_layers

	class(ty_com_nodes),intent(inout)	::nodes
	class(ty_unsat_mesh),intent(in)		::mesh		!< Class with the properties of the unsat mesh.
	class(ty_com_layers),intent(in)	::layers	!< Class with the properties of the unsat layers.

	integer::n,nc,id,e,idnode,ncount,i
	real(kind=dps)::incx

	!Number of nodes, nodes on element, class of nodes and allocate all arrays
	call nodes%allocateall(mesh%nnodclassv_count,mesh%nn,mesh%nc)

	!assign id at each node, and count number of nodes and classes of nodes
	do n=1, nodes%count
		do nc=1,mesh%nc+1
			idnode = 1+(n-1)*nc
			nodes%id(idnode) = idnode
			nodes%idnode(idnode) = n
			nodes%ncnode(idnode) = nc-1
		end do
	end do

	!Assign x (z coord) of nodes (meshing):
	ncount=0
	do i=1,layers%count
		incx = layers%height(i)/real(mesh%nelemv(i),dps) !Assuming uniform distribution of x

		do n=1, mesh%nelemv(i)+1
			ncount = sum(mesh%nelemv(1:i-1))+n
			do nc=1,mesh%nc+1
				idnode = 1+(ncount-1)*nc
				nodes%x(idnode)=real(n-1,dps)*incx+layers%hbottom(i)
				nodes%minmeshlenght(idnode) = abs(incx) !Define min mesh lenght on node
			end do
		end do
	end do

	nodes%x=nodes%x-sum(layers%height)

	!!Assign x of nodes (meshing):
	!ncount=0
	!!On top layer
	!incx = -layers%height(layers%count)/real(mesh%nelemv(layers%count),dps)
	!
	!do n=1, mesh%nelemv(layers%count)+1
	!	ncount = n
	!	do nc=1,mesh%nc+1
	!		idnode = 1+(ncount-1)*nc
	!		nodes%x(idnode)=real(n-1,dps)*incx
	!		nodes%minmeshlenght(idnode) = abs(incx) !Define min mesh lenght on node
	!	end do
	!end do
	!
	!
	!do i=layers%count-1,1,-1
	!incx = -layers%height(i)/real(mesh%nelemv(i),dps)
	!
	!do n=1, mesh%nelemv(i)+1
	!	ncount = sum(mesh%nelemv(i+1:layers%count))+n
	!	do nc=1,mesh%nc+1
	!		idnode = 1+(ncount-1)*nc
	!		nodes%x(idnode)=real(n-1,dps)*incx+(layers%htop(i)-layers%htop(layers%count))
	!		nodes%minmeshlenght(idnode) = abs(incx) !Define min mesh lenght on node
	!	end do
	!end do
	!end do

	!Assign hinit to nodes...
	select type(nodes)
	type is(ty_unsat_nodes)
		call nodes%set_hinit(layers%zphr)
	end select

	!Identify elements on node
	nodes%ne=1
	do e=1,mesh%nelemv_count
		do n=1,mesh%nn
			do nc=1,mesh%nc+1
				!        Nodeclass start on element +  nodeclassnumber inside element
				idnode = 1+(mesh%nn-1)*(e-1)        +  (n-1)*(mesh%nc+1)+(nc-1)
				if(n==1) then
					if (e>1) then
						nodes%ne(idnode)=2
						nodes%idelement(idnode,1)=e-1
						nodes%idelement(idnode,2)=e
					else
						nodes%ne(idnode)=1
						nodes%idelement(idnode,1)=e
					end if

				elseif (n==mesh%nn) then
					if (e<mesh%nn) then
						nodes%ne(idnode)=2
						nodes%idelement(idnode,1)=e
						nodes%idelement(idnode,2)=e+1
					else
						nodes%ne(idnode)=1
						nodes%idelement(idnode,1)=e
					end if
				else
					nodes%ne(idnode)=1
					nodes%idelement(idnode,1)=e
				end if
			end do
		end do
	end do

	!populateboundary conditions
	nodes%isdirichlet		= .false.
	nodes%isnewmann			= .false.
	nodes%qnewmann			= 0.0_dpd
	nodes%hdirichlet		=	0.0_dpd

	!Initial waterflow equal cero
	nodes%qent = 0.0_dpd


	end subroutine s_unsat_inputs_nodes

	!---------------------------------------------------------------------------------------------------------------------
	!> @author iván campos-guereta díez
	!> @brief
	!> procedure inside the class ty_com_elements. fill all elements data. first nodes has to be filled.
	!> @param[in] ne
	!> @param[in] nn
	!> @param[in] nc
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_inputs_elements(elements,mesh,nodes,layers)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_inputs_elements" :: s_unsat_inputs_elements
	!DEC$ endif
	use com_mod_ty_elements, only:ty_com_elements
	use com_mod_ty_nodes, only:ty_com_nodes
	use unsat_mod_ty_mesh, only:ty_unsat_mesh
	use com_mod_ty_layers, only:ty_com_layers

	class(ty_com_elements),intent(inout)::elements
	class(ty_unsat_mesh),intent(inout)::mesh
	class(ty_com_nodes),intent(inout)::nodes
	class(ty_com_layers),intent(inout)::layers

	integer::ne !number of elements
	integer::nn !number of nodes per element
	integer::nc !class number of nodes


	integer::nnc,i,j,l,e,nelem

	ne = mesh%nnodv_count-1
	nn = mesh%nn
	nc = mesh%nc

	nnc=nn*nc

	!allocate all vectors and update count, nn, nc, and nnc
	call elements%allocateall(ne,nn,nc)

	do i=1,ne
		elements%id(i) = i
		do j=1,nn
			elements%idnode(i,j) = j+(i-1)*(nn-1)
		end do
		elements%x0(i) = nodes%x(elements%idnode(i,1))
		elements%x1(i) = nodes%x(elements%idnode(i,nn))
		elements%lenght(i) = abs(elements%x1(i)-elements%x0(i))
		do j=1,nn
			elements%chi(i,j) = min(1.0_dpd,max(-1.0_dpd,elements%chi_from_x(i,nodes%x(elements%idnode(i,j)))))
		end do
	end do

	!Assign material to elements
	nelem=0
	do l=1,layers%count
		do e=1,mesh%nelemv(l)
			nelem=nelem+1
			elements%material(nelem) = layers%material(l)
		end do
	end do


	!Update material in nodes
	call elements%update_materials_to_nodes(nodes)

	end subroutine s_unsat_inputs_elements


	!---------------------------------------------------------------------------------------------------------------------
	!> @author iván campos-guereta díez
	!> @brief
	!> This subroutine reads boundary condition file and fill the class "ty_com_boundary".
	!> @param[in] ne
	!> @param[in] nn
	!> @param[in] nc
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_unsat_inputs_boundary(fileboundary,boundary,layers)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_unsat_inputs_boundary" :: s_unsat_inputs_boundary
	!DEC$ endif
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
			read (55,*, iostat=readerror) dumb,dumb
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
			read (55,*) boundary%timebound(i), boundary%hbound(i)
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
		boundary%timestepscount=boundary%timestepscount-2

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


	end subroutine s_unsat_inputs_boundary

	!********************************************************************************************************************
	! S: LOCATEBLOCK(FileIndex, Bloq, Readerror)
	!--------------------------------------------------------------------------------------------------------------------
	! This subROUTine rewind the file with index "FileIndex" and search for the line with the word: 'Block '+Bloq
	! that will define the initial line for the input of that block.
	! In case of error returns Readerror<>0
	!********************************************************************************************************************

	subroutine s_com_inputs_locateblock(fileindex, bloq, readerror)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_inputs_locateblock" :: s_com_inputs_locateblock
	!DEC$ endif
	! search for 'block x' in a file and let the file pointer at that row

	integer, intent(in):: fileindex
	character*1, intent(in)::bloq
	integer, intent(out)::readerror

	character*400:: chblock
	logical::checkblock

	rewind fileindex
	readerror = 0
	checkblock=.false.
	do while (.not.checkblock)
		read (fileindex,'(A)', iostat=readerror) chblock
		!read (fileindex,'(A)') chblock
		if (readerror.ne.0) then
			write(*,*) 'BLOCK '//bloq//' not found on input file'
			stop
		end if
		if(index(chblock,'BLOCK '//bloq) .ne. 0) checkblock=.true.
	end do

	end subroutine s_com_inputs_locateblock

	!********************************************************************************************************************
	! S: NEXTRECORD(FileIndex, Readerror)
	!--------------------------------------------------------------------------------------------------------------------
	! This subROUTine rewind the file with index "FileIndex" and search for the line with the word: 'Block '+Bloq
	! that will define the initial line for the input of that block.
	! In case of error returns Readerror<>0
	!********************************************************************************************************************

	subroutine s_com_inputs_nextrecord(fileindex, readerror)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_inputs_nextrecord" :: s_com_inputs_nextrecord
	!DEC$ endif
	! search for 'block x' in a file and let the file pointer at that row

	integer, intent(in):: fileindex
	integer, intent(out)::readerror

	character:: chkline
	logical::checkblock

	readerror = 0
	checkblock=.false.
	chkline='!'
	do while(chkline=='!')
		read (fileindex,'(A)', iostat=readerror) chkline
	end do
	backspace(fileindex)

	end subroutine s_com_inputs_nextrecord


	end module unsat_mod_inputs