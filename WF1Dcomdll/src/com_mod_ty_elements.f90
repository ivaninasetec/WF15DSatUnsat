	!********************************************************************************************************************
	!        CLASS FOR THE COLLECTION OF ELEMENTS IN THE COMMON LIBRARY
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D COMMON MODEL LIBRARIES
	! MODULE        : MOD_COM_TY_LAYER
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

	module com_mod_ty_elements
	use com_mod_ty_material, only: ty_com_material

	implicit none
	include 'inc_precision.fi'

	private

	type,public::ty_com_elements
		integer::count    !number of elements
		integer,allocatable::nn    !number of nodes per element
		integer,allocatable::nc    !number of clases of nodes (0=c0, 1=c1,...)
		integer::nnc							 !nodes and clases on elements

		integer,allocatable::id(:)     !id of element

		real(kind=dps),allocatable::x0(:) !<x coord of first node
		real(kind=dps),allocatable::x1(:)	!<x coord of last node
		class(ty_com_material),allocatable::material(:) !<Material of the element
		integer,allocatable::idnode(:,:) !<id to the material of the element (ne, nn)
		real(kind=dps),allocatable::chi(:,:) !<relative coord of each node inside element (2·(x-x0)/(x1-x0)-1)
		real(kind=dps),allocatable::lenght(:) !<Lenght of element
		real(kind=dps),allocatable::hnew(:)	  !<Current pressure head integrated on elment (n+1,k+1)
		real(kind=dps),allocatable::htemp(:)	!<Current pressure head integrated on element for previous iteration (n+1,k)
		real(kind=dps),allocatable::hold(:)	  !<Pressure head on the integrated on element previous coverged timestep (n)

	contains
	procedure,public:: allocateall	 => s_com_elements_allocateall
	procedure,public:: deallocateall => s_com_elements_deallocateall
	!procedure,public:: chi_from_nod => f_chi_from_nod
	procedure,public:: chi_from_x => f_com_element_chi_from_x
	procedure,public:: h_from_x => f_com_element_h_from_x
	procedure,public:: dh_from_x => f_com_element_dh_from_x
	procedure,public:: id_from_x => f_com_element_dh_from_x
	procedure,public:: id_from_x_sca => f_com_elements_idelement_from_x_sca
	!procedure,public:: h_on_nodes => f_h_on_nodes
	procedure,public:: update_materials_to_nodes	 => s_com_elements_update_materias_to_nodes


	end type 	ty_com_elements

	contains

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_elements. Initiallize ne, nn, nc, and nnc and allocate all arrays
	!> @param[in] ne
	!> @param[in] nn
	!> @param[in] nc
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_elements_allocateall(this,ne,nn,nc)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_elements_allocateall" :: s_com_elements_allocateall
	!DEC$ endif

	class(ty_com_elements),intent(inout)::this
	integer, intent(in)::ne !< Number of elements
	integer, intent(in)::nn !< Number of nodes on element
	integer, intent(in)::nc !< Class of nodes

	this%count	= ne
	this%nn			= nn
	this%nc			= nc
	this%nnc		= nn*(nc+1)

	if(.not.allocated(this%id))				allocate(this%id(ne))
	if(.not.allocated(this%x0))			allocate(this%x0(ne))
	if(.not.allocated(this%x1))			allocate(this%x1(ne))
	if(.not.allocated(this%material))	allocate(this%material(ne))
	if(.not.allocated(this%idnode))	  allocate(this%idnode(ne,nn))
	if(.not.allocated(this%chi))			allocate(this%chi(ne,nn))
	if(.not.allocated(this%lenght))		allocate(this%lenght(ne))
	if(.not.allocated(this%hnew))			allocate(this%hnew(ne))
	if(.not.allocated(this%htemp))		allocate(this%htemp(ne))
	if(.not.allocated(this%hold))			allocate(this%hold(ne))
	end subroutine s_com_elements_allocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_elements. Deallocate all arrays
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_elements_deallocateall(this)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_elements_deallocateall" :: s_com_elements_deallocateall
	!DEC$ endif

	class(ty_com_elements),intent(inout)::this

	if(allocated(this%id))				deallocate(this%id)
	if(allocated(this%x0))			deallocate(this%x0)
	if(allocated(this%x1))			deallocate(this%x1)
	if(allocated(this%material))deallocate(this%material)
	if(allocated(this%idnode))	  deallocate(this%idnode)
	if(allocated(this%chi))				deallocate(this%chi)
	if(allocated(this%lenght))		deallocate(this%lenght)
	if(allocated(this%hnew))			deallocate(this%hnew)
	if(allocated(this%htemp))			deallocate(this%htemp)
	if(allocated(this%hold))			deallocate(this%hold)


	end subroutine s_com_elements_deallocateall

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_elements. Construct all known values from nodes data.
	!> Update: x0,x1,idnode,chi,lenght
	!> @param[in] nodes
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_elements_update_from_nodes(this,nodes)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_elements_update_from_nodes" :: s_com_elements_update_from_nodes
	!DEC$ endif


	use com_mod_ty_nodes, only:ty_com_nodes

	class(ty_com_elements),intent(inout)::this
	class(ty_com_nodes),intent(in)::nodes

	real(kind=dpd)::xi,xf,xc,x
	integer::idnode,i,j

	do i=1, this%count
		this%id(i)=i
		xi=nodes%x((i-1)*(this%nn-1)+1)
		xf=nodes%x(i*(this%nn-1)+1)
		xc=(xi+xf)/2.0_dps

		this%x0(i)=xi
		this%x1(i)=xf

		this%lenght(i)=abs(xf-xi)

		do j=1, this%nn
			idnode=j+(i-1)*(this%nn-1)
			this%idnode(i,j)= idnode
			x = nodes%x(idnode)
			this%chi(i,j) = 2.0_dpd*(x-xc)/(xf-xi)
		end do

	end do

	end subroutine s_com_elements_update_from_nodes

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_elements. Return chi value from x in the element
	!> Update: x0,x1,idnode,chi,lenght
	!> @param[in] idelem
	!> @param[in] x
	!---------------------------------------------------------------------------------------------------------------------

	function f_com_element_chi_from_x(this,idelem,x)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_com_element_chi_from_x" :: f_com_element_chi_from_x
	!DEC$ endif

	class(ty_com_elements),intent(inout)::this
	integer,intent(in)::idelem
	real(kind=dpd),intent(in)::x
	REAL(KIND=dps)::f_com_element_chi_from_x
	REAL(KIND=dps)::xi,xf,xc

	xi = this%x0(idelem)
	xf = this%x1(idelem)
	xc = (xi+xf)/2.0_dps

	f_com_element_chi_from_x = 2.0_dps*(x-xc)/(xf-xi)

	end function f_com_element_chi_from_x

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_com_elements. Return value of h in element in x (opt 0: hnew(default), 1: htemp, 2:hold)
	!> Update: xini,xend,idnode,chi,lenght
	!> @param[in] idelem
	!> @param[in] x
	!---------------------------------------------------------------------------------------------------------------------

	function f_com_element_h_from_x(this,idelem,nodes,x,option)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_com_element_h_from_x" :: f_com_element_h_from_x
	!DEC$ endif

	use com_mod_fem_shapefunctions,only:shape1d
	use com_mod_ty_nodes,only:ty_com_nodes

	class(ty_com_elements),intent(inout)::this
	class(ty_com_nodes),intent(inout)::nodes

	integer,intent(in)::idelem
	real(kind=dpd),intent(in)::x
	REAL(KIND=dps)::f_com_element_h_from_x

	REAL(KIND=dps)::htemp(this%nn),shapeonnodestemp(this%nn)
	INTEGER,INTENT(IN),OPTIONAL::option
	INTEGER::opt,i,idnodeini,idnodeend

	opt=0

	IF(PRESENT(option)) opt = option

	idnodeini = this%idnode(idelem,1)
	idnodeend = this%idnode(idelem,this%nn)

	SELECT CASE (opt)
	CASE (1)
		htemp = nodes%htemp(idnodeini:idnodeend)
	CASE (2)
		htemp = nodes%hold(idnodeini:idnodeend)
		CASE DEFAULT
		htemp = nodes%hnew(idnodeini:idnodeend)
	END SELECT

	shapeonnodestemp = Shape1D(this%chi_from_x(idelem,x),this%chi(idelem,:))

	f_com_element_h_from_x = DOT_PRODUCT(htemp,shapeonnodestemp)

	end function f_com_element_h_from_x

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_com_elements. Return value of dh in element in x (opt 0: hnew(default), 1: htemp, 2:hold)
	!> Update: xini,xend,idnode,chi,lenght
	!> @param[in] idelem
	!> @param[in] x
	!---------------------------------------------------------------------------------------------------------------------

	function f_com_element_dh_from_x(this,idelem,nodes,x,option)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_com_element_dh_from_x" :: f_com_element_dh_from_x
	!DEC$ endif

	use com_mod_fem_shapefunctions,only:dphi1d
	use com_mod_ty_nodes,only:ty_com_nodes

	class(ty_com_elements),intent(inout)::this
	class(ty_com_nodes),intent(inout)::nodes

	integer,intent(in)::idelem
	real(kind=dpd),intent(in)::x
	real(kind=dps)::f_com_element_dh_from_x

	real(kind=dps)::htemp(this%nn),dphionnodestemp(this%nn)
	integer,intent(in),optional::option
	integer::opt,i,idnodeini,idnodeend

	opt=0

	if(present(option)) opt = option

	idnodeini = this%idnode(idelem,1)
	idnodeend = this%idnode(idelem,this%nn)

	select case (opt)
	case (1)
		htemp = nodes%htemp(idnodeini:idnodeend)
	case (2)
		htemp = nodes%hold(idnodeini:idnodeend)
		case default
		htemp = nodes%hnew(idnodeini:idnodeend)
	end select


	!shapeonnodestemp = dshape1d(this%chi_from_x(idelem,x),this%chi(idelem,:))
	dphionnodestemp = dphi1d(this%chi_from_x(idelem,x),nodes%x(idnodeini:idnodeend))

	f_com_element_dh_from_x = dot_product(htemp,dphionnodestemp)

	end function f_com_element_dh_from_x

	!!---------------------------------------------------------------------------------------------------------------------
	!!> @author Iván Campos-Guereta Díez
	!!> @brief
	!!> Procedure inside the class ty_com_elements. Return value of id of element given x (opt 0: hnew(default), 1: htemp, 2:hold)
	!!> Update: xini,xend,idnode,chi,lenght
	!!> @param[in] idelem
	!!> @param[in] x
	!!---------------------------------------------------------------------------------------------------------------------
	!
	!function f_com_elements_idelement_from_x(this,x) result(idelement) !default hnew, 1 hold, 2 htemp
	!
	!class(ty_com_elements),intent(inout)::this
	!
	!real(kind=dps),intent(in)::x
	!real(kind=dps)::idelement
	!real(kind=dps)::chi
	!
	!integer::e
	!
	!idelement = 1
	!do e=1,this%count
	!	if(((x>=this%xini(e)).and.(x<=this%xend(e))).or.((x<=this%xini(e)).and.(x>=this%xend(e)))) then
	!		idelement = e
	!		exit
	!	end if
	!end do
	!
	!end function f_com_elements_idelement_from_x


	!---------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Returns the id of the layer given the height h.
	!> @param[in] h
	!---------------------------------------------------------------------------

	function f_com_elements_idelement_from_x(this,x)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_com_elements_idelement_from_x" :: f_com_elements_idelement_from_x
	!DEC$ endif

	class(ty_com_elements),intent(in)::this
	real(kind=dpd),intent(in)::x(:)
	integer::f_com_elements_idelement_from_x(size(x))
	integer::i,findtemp(1)
	logical::mask(this%count)


	DO i=1, size(x)
		mask = this%x0<=x(i)
		findtemp = FINDLOC(mask,.true.,BACK=.true.)
		f_com_elements_idelement_from_x(i) = findtemp(1)
	end do

	end function f_com_elements_idelement_from_x

	function f_com_elements_idelement_from_x_sca(this,x)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"f_com_elements_idelement_from_x_sca" :: f_com_elements_idelement_from_x_sca
	!DEC$ endif


	class(ty_com_elements),intent(in)::this
	real(kind=dpd),intent(in)::x
	integer::f_com_elements_idelement_from_x_sca
	integer::i,findtemp(1)
	logical::mask(this%count)

	mask = this%x0<=x
	findtemp = FINDLOC(mask,.true.,BACK=.true.)
	f_com_elements_idelement_from_x_sca = findtemp(1)

	end function f_com_elements_idelement_from_x_sca

	!---------------------------------------------------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez
	!> @brief
	!> Procedure inside the class ty_sat_nodes. Deallocate all vectors and instance of calc.
	!---------------------------------------------------------------------------------------------------------------------

	subroutine s_com_elements_update_materias_to_nodes(this,nodes)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_elements_update_materias_to_nodes" :: s_com_elements_update_materias_to_nodes
	!DEC$ endif

	use com_mod_ty_nodes,only:ty_com_nodes

	class(ty_com_nodes),intent(inout)::nodes
	class(ty_com_elements),intent(inout)::this

	integer:: idnode,ne
	!assign material on node (as the mean of parameter when materials are the same or as first material when elements are different)
	do idnode=1,nodes%count
		do ne=1, nodes%ne(idnode) !number of elements on node
			if (ne == 1) then
				nodes%material(idnode) = this%material(nodes%idelement(idnode,ne))
			else
				if (nodes%material(idnode)%kindmat.eq.this%material(nodes%idelement(idnode,ne))%kindmat) then

					nodes%material(idnode)% thsat = nodes%material(idnode)% thsat+this%material(nodes%idelement(idnode,ne))% thsat
					nodes%material(idnode)% thres = nodes%material(idnode)% thres+this%material(nodes%idelement(idnode,ne))% thres
					nodes%material(idnode)% ksat = nodes%material(idnode)% ksat+this%material(nodes%idelement(idnode,ne))% ksat
					nodes%material(idnode)% a = nodes%material(idnode)% a+this%material(nodes%idelement(idnode,ne))% a
					nodes%material(idnode)% n = nodes%material(idnode)% n+this%material(nodes%idelement(idnode,ne))% n
					nodes%material(idnode)% m = nodes%material(idnode)% m+this%material(nodes%idelement(idnode,ne))% m
					nodes%material(idnode)% l = nodes%material(idnode)% l+this%material(nodes%idelement(idnode,ne))% l
				end if
			end if

		end do
		nodes%material(idnode)% thsat = nodes%material(idnode)% thsat / real(nodes%ne(idnode),dps)
		nodes%material(idnode)% thres = nodes%material(idnode)% thres / real(nodes%ne(idnode),dps)
		nodes%material(idnode)% ksat = nodes%material(idnode)% ksat / real(nodes%ne(idnode),dps)
		nodes%material(idnode)% a = nodes%material(idnode)% a / real(nodes%ne(idnode),dps)
		nodes%material(idnode)% n = nodes%material(idnode)% n / real(nodes%ne(idnode),dps)
		nodes%material(idnode)% m = nodes%material(idnode)% m / real(nodes%ne(idnode),dps)
		nodes%material(idnode)% l = nodes%material(idnode)% l / real(nodes%ne(idnode),dps)
	end do

	end subroutine s_com_elements_update_materias_to_nodes




	end module com_mod_ty_elements