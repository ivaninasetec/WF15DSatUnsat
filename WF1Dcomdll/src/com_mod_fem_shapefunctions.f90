	!********************************************************************************************************************
	!        SHAPE FUNCTIONS SUBROUTINES FOR FINITE ELEMENT METHOD
  !********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D COMMON LIBRARIES
	! MODULE        : MOD_COM_SHAPEFUNCTIONS
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
	!> This module include shapefunctions for 1d elements of nne nodes of class 0, and their derivatives               
	!> 
	!> At this moment it is only implemented for the case of n nodes of class 0  (not upper classes)
	!>  
	!>    Shape1D(chi):   Calculate value of shape functions of 2 nodes at relative coord 2. Chi can be scalar or vector
	!>                    If chi is scalar, the result is a vector with the value of every each shape function          
	!>                    (1 per node in element) at relative coord chi.                                                
	!>                    If chi is a vector chi(npt), the result is a matrix (nbasis,npt) with the values of shape     
	!>                    functions in every relative coord chi. Chi goes from -1 to 1.                                 
	!>                                                                                                                  
	!>    Shape1D(chi,nne): The same but with nne nodes in element                                                      
	!>                                                                                                                  
	!>    dShape1D(chi):  The same that Shape1D(chi) but returns the derivative at chi for each shape function          
	!>                                                                                                                  
	!>    dShape1D(chi,nne): The same that Shape1D(chi,nne) but returns the derivative at chi for each shape function   
  !********************************************************************************************************************
  
  
  module com_mod_fem_shapefunctions


	implicit none
  
  !********************************************************************************************************************
	! POLIMORPHIC INTERFACES FOR THE FOLLOWING FUNCTIONS 
	!--------------------------------------------------------------------------------------------------------------------
	! Shape1D (chi, [chinod])		: Returns shape function value at relative coord chi (sca or vec), with relative coord on
  !                             nodes equal to chinod(nne). If chinod not specified, then is supposed 2 nodes element 
  !                             with relative coord -1 and 1.
	! dShape1D (chi, [chinod])	: Returns shape function derivative at relative coord chi (sca or vec), with relative 
  !                             coord on nodes equal to chinod(nne). If chinod not specified, then is supposed 2 nodes
  !                             element with relative coord -1 and 1.
	! Interp_onelement(chi,val_on_nodes): Interpolate on element at relative coords of element equal to chi, with value
  !                             on nodes equal to val_on_nodes(nne). (only implemented with 2 nodes at this moment)
	!********************************************************************************************************************

	!INTERFACE FOR SHAPE FUNCTIONS
	interface shape1d
	module procedure shape1d2n_sca !shape1d(chi) -> out: (nne)
	module procedure shape1d2n_vec !shape1d(chi(:)) -> out: (nne,npt)
	module procedure shape1dnn_sca !shape1d(chi,chinod) -> out: (nne)
	module procedure shape1dnn_vec !shape1d(chi(:),chinod) -> out: (nne,npt)
	end interface shape1d
 
	!INTERFACE FOR DERIVATIVE OF SHAPE FUNCTIONS (In relative coords) (it is not the derivative of basis function at chi)
	interface dshape1d
	module procedure dshape1d2n_sca !shape1d(chi) -> out: (nne)
	module procedure dshape1d2n_vec !shape1d(chi(:)) -> out: (nne,npt)
	module procedure dshape1dnn_sca !shape1d(chi,chinod) -> out: (nne)
	module procedure dshape1dnn_vec !shape1d(chi(:),chinod) -> out: (nne,npt)
	end interface dshape1d
	
	!INTERFACE FOR DERIVATIVE OF BASIS FUNCTION AT RELATIVE CHI (In relative coords)
	interface dphi1d
	module procedure dphi1d2n_sca !shape1d(chi) -> out: (nne)
	module procedure dphi1d2n_vec !shape1d(chi(:)) -> out: (nne,npt)
	end interface dphi1d
	
	!INTERFACE INTERPOLATE ON ELEMENT
	interface interp_on_element
	module procedure interp_onelement_sca !shape1d(chi) -> out: sca
	module procedure interp_onelement_vec !shape1d(chi(:)) -> out: (npt)
	end interface interp_on_element

	!INTERFACE FOR JACOBIAN
	interface jacobian1d
  module procedure jacobian1d2nod_sca		!jacobian1d(chi,xne) jacobian for element with 2 nodes (chi scalar)
  module procedure jacobian1d2nod_vec !jacobian for element with 2 nodes (chi vector
  end interface jacobian1d
	
	private
	public::shape1d,dshape1d,interp_on_element,calc_linear_interpolation_matrix,jacobian1d,dphi1d

	include 'inc_precision.fi'

  contains
  
  !********************************************************************************************************************
	! F: SHAPE1D2N_SCA(CHI)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function value at relative coord chi (sca). In this case element has only two 
  ! nodes at chi= -1 and at chi= 1.
  ! 
  ! Parameters:
  !   chi: Relative coord (sca) at wicht calculate values of each Shape Function  (element goes from -1 to 1)
  !
  ! Output:
  !   [nbasis]: Values at each basis (total nbasis)
  !   N[1] = (1-chi)/2
  !   N[2] = (1+chi)/2
	!********************************************************************************************************************

	pure function shape1d2n_sca(chi)
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"shape1d2n_sca" :: shape1d2n_sca
	!DEC$ endif


	integer, parameter:: nne=2		!nbasis

	real(kind=dpd),intent(in)::chi
	real(kind=dpd)::shape1d2n_sca(nne)

	shape1d2n_sca(1) = (1.0_dpd-chi)/2.0_dpd
	shape1d2n_sca(2) = (chi+1.0_dpd)/2.0_dpd

  end function shape1d2n_sca
  

  !********************************************************************************************************************
	! F: SHAPE1D2N_VEC(CHI)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function value at relative coord chi (vec). In this case element has only two 
  ! nodes at chi= -1 and at chi= 1.
  ! 
  ! Parameters:
  !   chi(npt): Relative coords (vec) at wicht calculate values of each Shape Function (element goes from -1 to 1)
  !
  ! Output:
  !   N[nbasis,npt]: Values at each basis (total nbasis) for every point defined at relative coord chi.
  !   N[1,:] = (1-chi(:))/2
  !   N[2,:] = (1+chi(:))/2
	!********************************************************************************************************************

	pure function shape1d2n_vec(chi)
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"shape1d2n_vec" :: shape1d2n_vec
	!DEC$ endif


	integer, parameter:: nne=2		!nbasis

	real(kind=dpd),intent(in)::chi(:)
	real(kind=dpd)::shape1d2n_vec(size(chi),nne)

	shape1d2n_vec(:,1) = (1.0_dpd-chi)/2.0_dpd
	shape1d2n_vec(:,2) = (chi+1.0_dpd)/2.0_dpd

  end function shape1d2n_vec


  !********************************************************************************************************************
	! F: SHAPE1DNN_SCA(CHI,CHINOD)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function value at relative coord chi (sca), with relative coord on 
  ! nodes equal to chinod(nne). 
  ! 
  ! Parameters:
  !   chi: Relative coords (sca) at wicht calculate values of each Shape Function (element goes from -1 to 1)
  !   chinod(nne): Relative coord of each of the nodes on element.
  !
  ! Output:
  !   [nbasis]: Values at each basis (total nbasis)
  !   N[i:1-nbasis] = 
  !         Mul_j:1-nbasis (and j<>i)
  !           (chi-chinod(j))/(chinod(i)-chinod(j))  
	!********************************************************************************************************************

	pure function shape1dnn_sca(chi,chinod) !lineal shape functions for 1 dimension and nn nodes
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"shape1dnn_sca" :: shape1dnn_sca
	!DEC$ endif


	real(kind=dpd),intent(in)::chi
    real(kind=dpd),intent(in)::chinod(:)
	real(kind=dpd)::shape1dnn_sca(size(chinod))

	integer::i,j,nne
  nne = size(chinod)

	do i=1, nne
		shape1dnn_sca(i) = 1.0
		do j=1, nne
			if (j.ne.i) shape1dnn_sca(i) = shape1dnn_sca(i)*(chi-chinod(j))/(chinod(i)-chinod(j))
		end do
	end do

  end function shape1dnn_sca


  !********************************************************************************************************************
	! F: SHAPE1DNN_VEC(CHI,CHINOD)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function value at relative coord chi (vec), with relative coord on 
  ! nodes equal to chinod(nne). 
  ! 
  ! Parameters:
  !   chi(npt): Relative coords (vec) at wicht calculate values of each Shape Function (element goes from -1 to 1)
  !   chinod(nne): Relative coord of each of the nodes on element.
  !
  ! Output:
  !   N[nbasis,npt]: Values at each basis (total nbasis)
  !   N[i:1-nbasis,1-npt] = 
  !         Mul_j:1-nbasis (and j<>i)
  !           (chi(1:npt)-chinod(j))/(chinod(i)-chinod(j))  
	!********************************************************************************************************************

	pure function shape1dnn_vec(chi,chinod) !lineal shape functions for 1 dimension and nn nodes
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"shape1dnn_vec" :: shape1dnn_vec
	!DEC$ endif

	
	real(kind=dpd),intent(in)::chi(:)
	real(kind=dpd),intent(in)::chinod(:)
    real(kind=dpd)::shape1dnn_vec(size(chi),size(chinod))

    integer::nne
	integer::i,j

    nne = size(chinod)

	!do i=1, nne
	!	chinod(i) = real(i-1)/real(nn-1)
	!end do

	do i=1, nne
		shape1dnn_vec(:,i) = 1.0
		do j=1, nne
			if (j.ne.i) shape1dnn_vec(:,i) = shape1dnn_vec(:,i)*(chi-chinod(j))/(chinod(i)-chinod(j))
		end do
	end do

  end function shape1dnn_vec


  !********************************************************************************************************************
	! F: DSHAPE1D2N_SCA(CHI)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function derivative at relative coord chi (sca). In this case element has only two 
  ! nodes at chi= -1 and at chi= 1.
  ! 
  ! Parameters:
  !   chi: Relative coord (sca) at wicht calculate values of each Shape Function  (element goes from -1 to 1)
  !
  ! Output:
  !   dN[nbasis=2]: Derivatives at each basis (total nbasis)
  !   dN[1] = -1.0_dps
  !   dN[2] =  1.0_dps
	!********************************************************************************************************************

	pure function dshape1d2n_sca(chi) !derivative of shape function for 1 dimension and 2 nodes per element
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"dshape1d2n_sca" :: dshape1d2n_sca
	!DEC$ endif


	integer, parameter:: ndim=1 !number of dimensions=1
	integer, parameter:: nne=2		!number of nodes in element=2

	real(kind=dpd),intent(in)::chi
	real(kind=dpd)::dshape1d2n_sca(nne)

	!this is dn/dchi...
	dshape1d2n_sca(1) = -1.0_dps/2.0_dps
	dshape1d2n_sca(2) = 1.0_dps/2.0_dps

  end function dshape1d2n_sca


  !********************************************************************************************************************
	! F: DSHAPE1D2N_VEC(CHI)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function derivative at relative coord chi (vec). In this case element has only two 
  ! nodes at chi= -1 and at chi= 1.
  ! 
  ! Parameters:
  !   chi(npt): Relative coords (vec) at wicht calculate values of each Shape Function (element goes from -1 to 1)
  !
  ! Output:
  !   dN[nbasis=2,npt]: Derivatives at each basis (total nbasis) for every point defined at relative coord chi(npt).
  !   dN[1,1-npt] = -1.0_dps
  !   dN[2,1-npt] =  1.0_dps
  !********************************************************************************************************************

	pure function dshape1d2n_vec(chi) !derivative of shape function for 1 dimension and 2 nodes per element
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"dshape1d2n_vec" :: dshape1d2n_vec
	!DEC$ endif

	
	integer, parameter:: ndim=1 !number of dimensions=1
	integer, parameter:: nne=2		!number of nodes in element=2

	real(kind=dpd),intent(in)::chi(:)
	real(kind=dpd)::dshape1d2n_vec(size(chi),nne)

	dshape1d2n_vec(:,1) = -1.0_dps/2.0_dps
	dshape1d2n_vec(:,2) = 1.0_dps/2.0_dps

	end function dshape1d2n_vec

	  !********************************************************************************************************************
	! F: DSHAPE1D2N_SCA(CHI)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function derivative at relative coord chi (sca). In this case element has only two 
  ! nodes at chi= -1 and at chi= 1.
  ! 
  ! Parameters:
  !   chi: Relative coord (sca) at wicht calculate values of each Shape Function  (element goes from -1 to 1)
  !
  ! Output:
  !   dN[nbasis=2]: Derivatives at each basis (total nbasis)
  !   dN[1] = -1.0_dps
  !   dN[2] =  1.0_dps
	!********************************************************************************************************************

	pure function dphi1d2n_sca(chi,xnod) !derivative of shape function for 1 dimension and 2 nodes per element
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"dphi1d2n_sca" :: dphi1d2n_sca
	!DEC$ endif

	
	integer, parameter:: ndim=1 !number of dimensions=1
	integer, parameter:: nne=2		!number of nodes in element=2

	real(kind=dpd),intent(in)::chi
	real(kind=dpd),intent(in)::xnod(:) !Absolute coordinate values on nodes
	real(kind=dpd)::dphi1d2n_sca(size(xnod))

	!this is dn/dchi...
	dphi1d2n_sca(1) = -(1.0_dps/2.0_dps)/Jacobian1D(chi,xnod)
	dphi1d2n_sca(2) = (1.0_dps/2.0_dps)/Jacobian1D(chi,xnod)

  end function dphi1d2n_sca


  !********************************************************************************************************************
	! F: DSHAPE1D2N_VEC(CHI)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function derivative at relative coord chi (vec). In this case element has only two 
  ! nodes at chi= -1 and at chi= 1.
  ! 
  ! Parameters:
  !   chi(npt): Relative coords (vec) at wicht calculate values of each Shape Function (element goes from -1 to 1)
  !
  ! Output:
  !   dN[nbasis=2,npt]: Derivatives at each basis (total nbasis) for every point defined at relative coord chi(npt).
  !   dN[1,1-npt] = -1.0_dps
  !   dN[2,1-npt] =  1.0_dps
  !********************************************************************************************************************

	pure function dphi1d2n_vec(chi,xnod) !derivative of shape function for 1 dimension and 2 nodes per element
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"dphi1d2n_vec" :: dphi1d2n_vec
	!DEC$ endif

	
	integer, parameter:: ndim=1 !number of dimensions=1
	integer, parameter:: nne=2		!number of nodes in element=2

	real(kind=dpd),intent(in)::chi(:)	!Relative chi at which we want to calculate derivative
	real(kind=dpd),intent(in)::xnod(:) !Absolute coordinate values on nodes
	real(kind=dpd)::dphi1d2n_vec(size(chi),size(xnod)) !Returns the derivative in chi

	dphi1d2n_vec(:,1) = -(1.0_dps/2.0_dps)/Jacobian1D(chi,xnod)
	dphi1d2n_vec(:,2) = (1.0_dps/2.0_dps)/Jacobian1D(chi,xnod)

  end function dphi1d2n_vec
	

  !********************************************************************************************************************
	! F: DSHAPE1DNN_SCA(CHI,CHINOD)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function derivatives at relative coord chi (sca), with relative coord on 
  ! nodes equal to chinod(nne). 
  ! 
  ! Parameters:
  !   chi: Relative coords (sca) at wicht calculate values of each Shape Function (element goes from -1 to 1)
  !   chinod(nne): Relative coord of each of the nodes on element.
  !
  ! Output:
  !   N[nbasis]: Derivatives at each basis (total nbasis)
  !   N[i:1-nbasis] = 
  !      Sum_j:1-nbasis (and i<>j)
  !         Mul_m:1-nbasis (and m<>i, and m<>j)
  !           (chi-chinod(m))/(chinod(j)-chinod(m))
	!********************************************************************************************************************

	pure function dshape1dnn_sca(chi,chinod) !derivative of shape function for 1 dimension and 2 nodes per element
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"dshape1dnn_sca" :: dshape1dnn_sca
	!DEC$ endif


	integer, parameter:: ndim=1 !number of dimensions=1

	real(kind=dpd),intent(in)::chi
  real(kind=dpd),intent(in)::chinod(:)
	real(kind=dpd)::dshape1dnn_sca(size(chinod))
	real(kind=dpd)::rtemp

	integer::nne
	integer::i,j,m

  nne = size(chinod)
    
	do i=1, nne !i is the number of shape function for the node
		dshape1dnn_sca(i) = 0.0
		do j=1, nne
			if (i.ne.j) then
				rtemp= 1.0_dpd
				do m=1,nne
					if (m.ne.i.and.m.ne.j) rtemp = rtemp*(chi-chinod(m))/(chinod(j)-chinod(m))
				end do
				dshape1dnn_sca(i) = dshape1dnn_sca(i)+rtemp/(chinod(j)-chinod(i))
			end if
		end do
	end do

  end function dshape1dnn_sca


  !********************************************************************************************************************
	! F: DSHAPE1DNN_VEC(CHI,CHINOD)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns shape function derivatives at relative coord chi (vec), with relative coord on 
  ! nodes equal to chinod(nne). 
  ! 
  ! Parameters:
  !   chi(npt): Relative coords (vec) at wicht calculate values of each Shape Function (element goes from -1 to 1)
  !   chinod(nne): Relative coord of each of the nodes on element.
  !
  ! Output:
  !   dN[nbasis,npt]: Derivatives at each basis (total nbasis) in every point defined with relative coords chi(npt)
  !   dN[i:1-nbasis,npt] = 
  !      Sum_j:1-nbasis (and i<>j)
  !         Mul_m:1-nbasis (and m<>i, and m<>j)
  !           (chi(1:npt)-chinod(m))/(chinod(j)-chinod(m))
	!********************************************************************************************************************

	pure function dshape1dnn_vec(chi,chinod) !derivative of shape function for 1 dimension and 2 nodes per element
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"dshape1dnn_vec" :: dshape1dnn_vec
	!DEC$ endif


	integer, parameter:: ndim=1 !number of dimensions=1

	real(kind=dpd),intent(in)::chi(:)
	real(kind=dpd),intent(in)::chinod(:)
    real(kind=dpd)::dshape1dnn_vec(size(chi),size(chinod))
	real(kind=dpd)::rtemp(size(chi))

    integer::nne
    integer::i,j,m
    
    nne = size(chinod)

	do i=1, nne !i: id of shape function
		dshape1dnn_vec(:,i) = 0.0
		do j=1, nne
			if (i.ne.j) then
				rtemp(:)= 1.0_dpd
				do m=1,nne
					if (m.ne.i.and.m.ne.j) rtemp(:) = rtemp(:)*(chi(:)-chinod(m))/(chinod(j)-chinod(m))
				end do
				dshape1dnn_vec(:,i) = dshape1dnn_vec(:,i)+rtemp(:)/(chinod(j)-chinod(i))
			end if
		end do
	end do

	end function dshape1dnn_vec





  !********************************************************************************************************************
	! F: INTERP_ONELEMENT_SCA(CHI,VAL_ON_NODES)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns the interpolated value at relative coord chi in element with values on nodes (only 2 nodes)
  ! equal to VAL_ON_NODES, using shape functions of 2 nodes Shape1D(chi)
  ! 
  ! Parameters:
  !   chi: Relative coord (sca) at wicht calculate value.
  !   val_on_nodes(nne=2): Values of the function on each of the 2 nodes.
  !
  ! Output:
  !   Interp(scalar): Interpolated value at chi.
  !   Interp={vi(nbasis)}'·{Ni(nbasis)}
	!********************************************************************************************************************

	pure function interp_onelement_sca(chi,val_on_nodes)
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"interp_onelement_sca" :: interp_onelement_sca
	!DEC$ endif

	!returns absolute coordinate given relative x=p(chi)

	!integer, parameter:: ndim=1				!number of dimensions=1
	!integer, parameter:: nne=2				!number of nodes in element=2
	real(kind=dps),intent(in)::chi		!relative coord at with we need to compute absolute coord
	real(kind=dps),intent(in)::val_on_nodes(:) !absolute coordinates of nodes of the element
	real(kind=dps)::interp_onelement_sca

	!quadrature1d = dot_product(xne,shape1d2n(chi)) !the absolute coordinate given relative x=sum(xi·ni)
  if (chi>=-1.0_dps.or.chi<=1.0_dps) then
	interp_onelement_sca = dot_product(val_on_nodes,shape1d(chi)) !the absolute coordinate given relative x=sum(xi·ni)
  else
    interp_onelement_sca=0.0_dps
  end if

	end function interp_onelement_sca

  
  !********************************************************************************************************************
	! F: INTERP_ONELEMENT_VEC(CHI,VAL_ON_NODES)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns the interpolated value at relative coord chi in element with values on nodes (only 2 nodes)
  ! equal to VAL_ON_NODES, using shape functions of 2 nodes Shape1D(chi)
  ! 
  ! Parameters:
  !   chi: Relative coord (sca) at wicht calculate value.
  !   val_on_nodes(nne=2): Values of the function on each of the 2 nodes.
  !
  ! Output:
  !   Interp[npt]: Interpolated value at each point defined in relative coords chi(npt).
  !   Interp={vi(nbasis)}'·[Ni(nbasis,npt)]
	!********************************************************************************************************************

	!---------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez - Univisersity of Nottingham
	!> @brief  
	!> Function that Returns the interpolated value at relative coord chi in element with values on nodes (only 2 nodes)
	!> equal to VAL_ON_NODES, using shape functions of 2 nodes Shape1D(chi)
	!>
  !> Output:
  !>   Interp[npt]: Interpolated value at each point defined in relative coords chi(npt).
  !>   Interp={vi(nbasis)}'·[Ni(nbasis,npt)]
	!> @param[in] chi(:)
	!> @param[in] val_on_nodes(nne)
	!---------------------------------------------------------------------------
	
	pure function interp_onelement_vec(chi,val_on_nodes)
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"interp_onelement_vec" :: interp_onelement_vec
	!DEC$ endif

	!returns absolute coordinate given relative x=p(chi)
	
	!integer, parameter:: ndim=1				!number of dimensions=1
	!integer, parameter:: nne=2				!number of nodes in element=2
	
	
	real(kind=dps),intent(in)::chi(:)								!<relative coord at with we need to compute absolute coord
	real(kind=dps),intent(in)::val_on_nodes(:)		!<absolute coordinates of nodes of the element
	real(kind=dps)::interp_onelement_vec(size(chi)) !<output: interpolated values over element.

  where (chi>=-1.0_dps.or.chi<=1.0_dps)
	interp_onelement_vec = matmul(shape1d(chi),val_on_nodes)
  else where
  interp_onelement_vec = 0.0_dps
  end where
  
	end function interp_onelement_vec
	
	
	!---------------------------------------------------------------------------
	!> @author Iván Campos-Guereta Díez - Univisersity of Nottingham
	!> @brief  
	!> Create interpolation matrix (given values at Xi, and shapefunctions return the interpolation matrix)
	!> such as {xi}=[R(i,j)]·{Xj}
	!>
  !> Output:
  !>   Interp[npt]: Interpolated value at each point defined in relative coords chi(npt).
  !>   Interp={vi(nbasis)}'·[Ni(nbasis,npt)]
	!> @param[in] chi(:)
	!> @param[in] val_on_nodes(nne)
	!---------------------------------------------------------------------------
	
	pure function calc_linear_interpolation_matrix(xi,Xj) !CHECK: To include class 0)
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"calc_linear_interpolation_matrix" :: calc_linear_interpolation_matrix
	!DEC$ endif

	
	real(kind=dpd),intent(in)::xi(:)
	real(kind=dpd),intent(in)::Xj(:)
	real(kind=dpd)::calc_linear_interpolation_matrix(size(xi),size(Xj))
	
	integer::nnodes,nvmod,nelem,ivmod,inod
	real(kind=dpd)::x,x0,xinf,xsup,chisup,chiinf
	
	nnodes	= size(xi)
	nvmod		= size(Xj)
	nelem = nvmod-1
	
	do ivmod=1,nvmod
		do inod=1,nnodes
			if (ivmod==1) then
			x = xi(inod)
			x0 = Xj(ivmod)
			xsup=Xj(ivmod+1)
			chisup = 1.0_dpd-(x-x0)/(xsup-x0) !From 0 to 1
			if (chisup>1.0_dpd) chisup = 0.0_dpd
			if (chisup<0.0_dpd) chisup = 0.0_dpd
			calc_linear_interpolation_matrix(inod,ivmod) = min(1.0_dpd,max(0.0_dpd,chisup))
			else if (ivmod==nvmod) then
			x = xi(inod)
			xinf=Xj(ivmod-1)
			x0 = Xj(ivmod)
			chiinf = (x-xinf)/(x0-xinf)
			if (chiinf>1.0_dpd) chiinf = 0.0_dpd
			if (chiinf<0.0_dpd) chiinf = 0.0_dpd
			calc_linear_interpolation_matrix(inod,ivmod) = min(1.0_dpd,max(0.0_dpd,chiinf))
			else
			x = xi(inod)
			xinf=Xj(ivmod-1)
			x0 = Xj(ivmod)
			xsup=Xj(ivmod+1)
			chisup = 1.0_dpd-(x-x0)/(xsup-x0) !From -1 to 1
			chiinf = (x-xinf)/(x0-xinf)
			if (chiinf>1.0_dpd) chiinf = 0.0_dpd
			if (chiinf<0.0_dpd) chiinf = 0.0_dpd
			if (chisup>1.0_dpd) chisup = 0.0_dpd
			if (chisup<0.0_dpd) chisup = 0.0_dpd			
			calc_linear_interpolation_matrix(inod,ivmod) = min(1.0_dpd,max(0.0_dpd,chisup+chiinf))
			end if
		end do
	end do
	
	end function calc_linear_interpolation_matrix
	
	!********************************************************************************************************************
	! F: JACOBIAN1D2NOD_SCA(CHI,XNE)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns the jacobian for 1D element of 2 nodes (scalar), in each chi relative coord
  !   chi:         Is the relative coord at wich calculate the jacobian
  !   xne(nnodes): Is the absolute coodinates of the nodes of element.  
	!********************************************************************************************************************
 
  pure function jacobian1d2nod_sca(chi,xne)
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"jacobian1d2nod_sca" :: jacobian1d2nod_sca
	!DEC$ endif

  !in this case chi is not needed but it is included for more general cases
 
  !integer, parameter:: ndim=1			!number of dimensions=1
  !integer, parameter:: nne=2				!number of nodes in element=2
  real(kind=dps),intent(in)::chi		!<relative coordinates of nodes of the element
  real(kind=dps),intent(in)::xne(:) !<absolute coordinates of nodes of the element
  real(kind=dps)::jacobian1d2nod_sca
	integer::nne
	nne=size(xne)
	
  if (nne==2) then
	jacobian1d2nod_sca = abs(xne(2)-xne(1))/2.0_dpd !this is dx/dchi
	else
		block
			real(kind=dps)::chinod(size(xne))
			jacobian1d2nod_sca = abs(dot_product(dshape1d(chi,chinod),xne)) !|dn_i(chi)/dchi·xnod_i|
		end block
	end if
	
	
  end function jacobian1d2nod_sca
  
  
	!********************************************************************************************************************
	! F: JACOBIAN1D2NOD_VEC(CHI,XNE)
	!--------------------------------------------------------------------------------------------------------------------
	! Function that Returns the jacobian for 1D element of 2 nodes (vector), in each chi relative coord
  !   chi(npt):    Is the relative coords at wich calculate the jacobian
  !   xne(nnodes): Is the absolute coodinates of the nodes of element.  
	!********************************************************************************************************************
  !--------------------------------------------------------------------------------------------------------
 
  pure function jacobian1d2nod_vec(chi,xne)
	!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"jacobian1d2nod_vec" :: jacobian1d2nod_vec
	!DEC$ endif

  !in this case chi is not needed but it is included for more general cases
  
  !integer, parameter:: ndim=1				!number of dimensions=1
  !integer, parameter:: nne=2				!number of nodes in element=2
  real(kind=dps),intent(in)::chi(:) !absolute coordinates of nodes of the element
  real(kind=dps),intent(in)::xne(:) !absolute coordinates of nodes of the element
  real(kind=dps)::jacobian1d2nod_vec(size(chi))
	integer::nne
	
	nne=size(xne)
	
	if (nne==2) then
		jacobian1d2nod_vec = abs(xne(2)-xne(1))/2.0_dpd !this is dx/dchi
	else
		block
			real(kind=dps)::chinod(size(xne))
			jacobian1d2nod_vec = abs(matmul(dshape1d(chi,chinod),xne)) !|dn_i(chi)/dchi·xnod_i|
		end block
	end if
	
	
	end function	jacobian1d2nod_vec

 ! !--------------------------------------------------------------------------------------------------------
 !
 !
 ! FUNCTION Jacobian1DNN(chi,xne)
 !!DEC$ if defined(_DLL)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"Jacobian1DNN" :: Jacobian1DNN
	!DEC$ endif

 ! 
 ! !INTEGER, PARAMETER:: ndim=1				!Number of dimensions=1
 ! !INTEGER, PARAMETER:: nne=2				!Number of nodes in element=2
 ! REAL(KIND=dps),INTENT(IN)::chi !Absolute coordinates of nodes of the element
 ! REAL(KIND=dps),INTENT(IN)::xne(:) !Absolute coordinates of nodes of the element
 ! REAL(KIND=dps)::Jacobian1DNN
 ! INTEGER nne
	!REAL(KIND=dps)::chinod(size(xne))
 ! 
 ! nne=SIZE(xne)
 ! 
	!chinod = 2*xne/(xne(nne)-xne(1))-1
	!
 ! Jacobian1D2Nod_vec = ABS(DOT_PRODUCT(xne,dshape1dnn(chi,chinod))) !|xnod_i·dN_i/dchi|
 ! 
 ! END FUNCTION Jacobian1DN

	end module com_mod_fem_shapefunctions