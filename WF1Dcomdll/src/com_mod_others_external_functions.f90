	module com_mod_external_functions

	implicit none
	public
	include 'inc_precision.fi'

	public::hyp

	contains
	!***********************************************************************
	!
	!   program :  hypergeometric function
	!
	!   notation :  F(a,b;c;z)
	!
	!   reference:  see article `Computing the hypergeometric function'
	!               by R.C. Forrey, J. Comp. Phys. 137, 79-100 (1997).
	!
	!   send comments to:
	!
	!        Robert C. Forrey
	!        Institute for Theoretical Atomic and Molecular Physics
	!        Harvard-Smithsonian Center for Astrophysics
	!        60 Garden Street Cambridge, MA 02138
	!        rforrey@cfa.harvard.edu
	!
	!***********************************************************************
	!
	!  subroutine name    - hyp
	!
	!  computation
	!  performed          - calculates the hypergeometric function
	!
	!  usage              - call hyp(z,a,b,c,re,im)
	!
	!  arguments
	!                  z  - the independent variable of the hypergeometric
	!                       function (must be real).
	!
	!               a,b,c - real parameters of the hypergeometric function.
	!
	!               re,im - the real and imaginary parts of the
	!                       hypergeometric function.
	!
	!  precision          - double
	!
	!  language           - fortran
	!
	!***********************************************************************

	FUNCTION HYP(Z,AINP,BINP,CINP) RESULT(Re)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"hyp" :: hyp
	!DEC$ endif





	IMPLICIT NONE
	!*--HYP45
	!*** Start of declarations inserted by SPAG

	REAL(KIND=DPD),INTENT(IN)::AINP , BINP , CINP , Z
	REAL(KIND=DPD)::Re,Im
	INTEGER k , m , n , nmax
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) ZERO , ONE , TWO , HALF
	PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,HALF=0.5D0)
	INTEGER flag , flag2 , neps
	REAL(KIND=DPD):: A , B , C , w , f1 , f2 , tol , test , pi ,     &
		& machep , BINom , eps , x1 , x2 , x3 , x4 ,       &
		& coeff1 , coeff2 , coeff3 , coeff4 , temp1 , temp2 , a1 ,   &
		& b1 , c1 , a2 , b2 , c2
	LOGICAL fix
	COMMON /BCOEFF/ BINom(5151)

	!  tabulate the binomial coefficients and set the defaults

	A = AINP
	B = BINP
	C = CINP

	fix = .FALSE.
	CALL BINOMC
	tol = .1D0
	Im = ZERO
	nmax = 100
	n = 5

	CALL GETEPS(machep,neps)
	pi = DACOS(-1.D0)
	!     write(20,'(/a,i3/)') ' machine epsilon is machep = (1/2)**',neps
	!     write(20,'(a,d23.15//)') ' to machine accuracy, pi = ',pi

	!  handle the special case when z=1

	IF ( Z==ONE ) THEN
		Re = GAMM(C)*GAMM(C-A-B)/GAMM(C-A)/GAMM(C-B)
		RETURN
	ENDIF

	!  transform to a new variable w which lies between 0 and 1/2

	IF ( Z<-ONE ) THEN
		a1 = A
		b1 = C - B
		c1 = A - B + 1
		a2 = B
		b2 = C - A
		c2 = B - A + 1
		w = ONE/(ONE-Z)
		flag = 1
	ELSEIF ( (-ONE<=Z) .AND. (Z<ZERO) ) THEN
		a1 = A
		b1 = C - B
		c1 = C
		a2 = a1
		b2 = b1
		c2 = c1
		w = Z/(Z-ONE)
		flag = 2
	ELSEIF ( (ZERO<=Z) .AND. (Z<=HALF) ) THEN
		a1 = A
		b1 = B
		c1 = C
		a2 = a1
		b2 = b1
		c2 = c1
		w = Z
		flag = 3
	ELSEIF ( (HALF<Z) .AND. (Z<=ONE) ) THEN
		a1 = A
		b1 = B
		c1 = A + B - C + 1
		a2 = C - A
		b2 = C - B
		c2 = C - A - B + 1
		w = ONE - Z
		flag = 4
	ELSEIF ( (ONE<Z) .AND. (Z<=TWO) ) THEN
		a1 = A
		b1 = A - C + 1
		c1 = A + B - C + 1
		a2 = C - A
		b2 = 1 - A
		c2 = C - A - B + 1
		w = ONE - (ONE/Z)
		flag = 5
	ELSEIF ( TWO<Z ) THEN
		a1 = A
		b1 = A - C + 1
		c1 = A - B + 1
		a2 = B - C + 1
		b2 = B
		c2 = B - A + 1
		w = ONE/Z
		flag = 6
	ENDIF

	!  compute the hypergeometric function of z via the transformation
	!  theory

	IF ( flag==1 ) THEN
		k = NINT(DBLE(A-B))
		test = A - B - DBLE(k)
		IF ( DABS(test)<tol ) THEN
			fix = .TRUE.
			flag2 = 0
			IF ( A<B ) THEN
				temp1 = A
				temp2 = B
				B = temp1
				A = temp2
				flag2 = 1
			ENDIF
			k = NINT(DBLE(A-B))
			eps = A - B - DBLE(k)
			CALL FIX1(A,B,C,n,k,f1,w,machep,eps)
			DO m = n + 5 , nmax , 5
				CALL FIX1(A,B,C,m,k,f2,w,machep,eps)
				test = DABS(f1-f2)
				IF ( test<=machep ) GOTO 20
				f1 = f2
			ENDDO
			WRITE (*,*) 'fix1 warning: not completely converged'
20		Re = f2
			IF ( flag2==1 ) THEN
				A = temp1
				B = temp2
			ENDIF
		ELSE
			CALL HYPER(w,a1,b1,c1,f1,machep)
			CALL HYPER(w,a2,b2,c2,f2,machep)

			x1 = B
			coeff1 = ONE
40		IF ( x1<ONE ) THEN
				coeff1 = coeff1*x1
				x1 = x1 + ONE
				GOTO 40
			ENDIF

			x2 = C - A
			coeff2 = ONE
60		IF ( x2<ONE ) THEN
				coeff2 = coeff2*x2
				x2 = x2 + ONE
				GOTO 60
			ENDIF

			x3 = A
			coeff3 = ONE
80		IF ( x3<ONE ) THEN
				coeff3 = coeff3*x3
				x3 = x3 + ONE
				GOTO 80
			ENDIF

			x4 = C - B
			coeff4 = ONE
100		IF ( x4<ONE ) THEN
				coeff4 = coeff4*x4
				x4 = x4 + ONE
				GOTO 100
			ENDIF

			Re = (w**A)*GAMM(C)*GAMM(B-A)*coeff1*coeff2/GAMM(x1)        &
				& /GAMM(x2)*f1 + (w**B)*GAMM(C)*GAMM(A-B)                &
				& *coeff3*coeff4/GAMM(x3)/GAMM(x4)*f2

		ENDIF
	ELSEIF ( flag==2 ) THEN
		CALL HYPER(w,a1,b1,c1,f1,machep)
		Re = ((ONE-w)**A)*f1
	ELSEIF ( flag==3 ) THEN
		CALL HYPER(w,a1,b1,c1,f1,machep)
		Re = f1
	ELSEIF ( flag==4 ) THEN
		k = NINT(DBLE(C-A-B))
		test = C - A - B - DBLE(k)
		IF ( DABS(test)<tol ) THEN
			fix = .TRUE.
			IF ( k>=ZERO ) THEN
				eps = C - A - B - DBLE(k)
				CALL FIX4A(A,B,C,n,k,f1,w,machep,eps)
				DO m = n + 5 , nmax , 5
					CALL FIX4A(A,B,C,m,k,f2,w,machep,eps)
					test = DABS(f1-f2)
					IF ( test<=machep ) GOTO 110
					f1 = f2
				ENDDO
				WRITE (*,*) 'fix4a warning: not completely converged'
110			Re = f2
			ELSE
				k = -k
				eps = C - A - B + DBLE(k)
				CALL FIX4B(A,B,C,n,k,f1,w,machep,eps)
				DO m = n + 5 , nmax , 5
					CALL FIX4B(A,B,C,m,k,f2,w,machep,eps)
					test = DABS(f1-f2)
					IF ( test<=machep ) GOTO 120
					f1 = f2
				ENDDO
				WRITE (*,*) 'fix4b warning: not completely converged'
120			Re = f2
			ENDIF
		ELSE
			CALL HYPER(w,a1,b1,c1,f1,machep)
			CALL HYPER(w,a2,b2,c2,f2,machep)

			x1 = C - A
			coeff1 = ONE
140		IF ( x1<ONE ) THEN
				coeff1 = coeff1*x1
				x1 = x1 + ONE
				GOTO 140
			ENDIF

			x2 = C - B
			coeff2 = ONE
160		IF ( x2<ONE ) THEN
				coeff2 = coeff2*x2
				x2 = x2 + ONE
				GOTO 160
			ENDIF

			x3 = A
			coeff3 = ONE
180		IF ( x3<ONE ) THEN
				coeff3 = coeff3*x3
				x3 = x3 + ONE
				GOTO 180
			ENDIF

			x4 = B
			coeff4 = ONE
200		IF ( x4<ONE ) THEN
				coeff4 = coeff4*x4
				x4 = x4 + ONE
				GOTO 200
			ENDIF

			Re = GAMM(C)*GAMM(C-A-B)*coeff1*coeff2/GAMM(x1)/GAMM(x2)    &
				& *f1 + w**(C-A-B)*GAMM(C)*GAMM(A+B-C)                   &
				& *coeff3*coeff4/GAMM(x3)/GAMM(x4)*f2

		ENDIF
	ELSEIF ( flag==5 ) THEN
		k = NINT(DBLE(C-A-B))
		test = C - A - B - DBLE(k)
		IF ( DABS(test)<tol ) THEN
			fix = .TRUE.
			IF ( k>=ZERO ) THEN
				eps = C - A - B - DBLE(k)
				CALL FIX5A(A,B,C,n,k,f1,Im,w,machep,eps,pi)
				DO m = n + 5 , nmax , 5
					CALL FIX5A(A,B,C,m,k,f2,Im,w,machep,eps,pi)
					test = DABS(f1-f2)
					IF ( test<=machep ) GOTO 210
					f1 = f2
				ENDDO
				WRITE (*,*) 'fix5a warning: not completely converged'
210			Re = f2
			ELSE
				k = -k
				eps = C - A - B + DBLE(k)
				CALL FIX5B(A,B,C,n,k,f1,Im,w,machep,eps,pi)
				DO m = n + 5 , nmax , 5
					CALL FIX5B(A,B,C,m,k,f2,Im,w,machep,eps,pi)
					test = DABS(f1-f2)
					IF ( test<=machep ) GOTO 220
					f1 = f2
				ENDDO
				WRITE (*,*) 'fix5b warning: not completely converged'
220			Re = f2
			ENDIF
		ELSE
			CALL HYPER(w,a1,b1,c1,f1,machep)
			CALL HYPER(w,a2,b2,c2,f2,machep)

			x1 = C - A
			coeff1 = ONE
240		IF ( x1<ONE ) THEN
				coeff1 = coeff1*x1
				x1 = x1 + ONE
				GOTO 240
			ENDIF

			x2 = C - B
			coeff2 = ONE
260		IF ( x2<ONE ) THEN
				coeff2 = coeff2*x2
				x2 = x2 + ONE
				GOTO 260
			ENDIF

			x3 = A
			coeff3 = ONE
280		IF ( x3<ONE ) THEN
				coeff3 = coeff3*x3
				x3 = x3 + ONE
				GOTO 280
			ENDIF

			x4 = B
			coeff4 = ONE
300		IF ( x4<ONE ) THEN
				coeff4 = coeff4*x4
				x4 = x4 + ONE
				GOTO 300
			ENDIF

			Re = (ONE-w)**A*GAMM(C)*GAMM(C-A-B)*coeff1*coeff2/GAMM(x1)  &
				& /GAMM(x2)*f1 + w**(C-A-B)*(ONE-w)**B*DCOS(pi*(C-A-B))  &
				& *GAMM(C)*GAMM(A+B-C)*coeff3*coeff3/GAMM(x3)/GAMM(x4)*f2

		ENDIF
	ELSEIF ( flag==6 ) THEN
		k = NINT(DBLE(A-B))
		test = A - B - DBLE(k)
		IF ( DABS(test)<tol ) THEN
			fix = .TRUE.
			flag2 = 0
			IF ( A<B ) THEN
				temp1 = A
				temp2 = B
				B = temp1
				A = temp2
				flag2 = 1
			ENDIF
			k = NINT(DBLE(A-B))
			eps = A - B - DBLE(k)
			CALL FIX6(A,B,C,n,k,f1,Im,w,machep,eps,pi)
			DO m = n + 5 , nmax , 5
				CALL FIX6(A,B,C,m,k,f2,Im,w,machep,eps,pi)
				test = DABS(f1-f2)
				IF ( test<=machep ) GOTO 320
				f1 = f2
			ENDDO
			WRITE (*,*) 'fix6 warning: not completely converged'
320		Re = f2
			IF ( flag2==1 ) THEN
				A = temp1
				B = temp2
			ENDIF
		ELSE
			CALL HYPER(w,a1,b1,c1,f1,machep)
			CALL HYPER(w,a2,b2,c2,f2,machep)

			x1 = B
			coeff1 = ONE
340		IF ( x1<ONE ) THEN
				coeff1 = coeff1*x1
				x1 = x1 + ONE
				GOTO 340
			ENDIF

			x2 = C - A
			coeff2 = ONE
360		IF ( x2<ONE ) THEN
				coeff2 = coeff2*x2
				x2 = x2 + ONE
				GOTO 360
			ENDIF

			x3 = A
			coeff3 = ONE
380		IF ( x3<ONE ) THEN
				coeff3 = coeff3*x3
				x3 = x3 + ONE
				GOTO 380
			ENDIF

			x4 = C - B
			coeff4 = ONE
400		IF ( x4<ONE ) THEN
				coeff4 = coeff4*x4
				x4 = x4 + ONE
				GOTO 400
			ENDIF

			Re = w**A*DCOS(pi*A)*GAMM(C)*GAMM(B-A)                      &
				& *coeff1*coeff2/GAMM(x1)/GAMM(x2)*f1 + w**B*DCOS(pi*B)  &
				& *GAMM(C)*GAMM(A-B)*coeff3*coeff4/GAMM(x3)/GAMM(x4)*f2

		ENDIF
	ENDIF

	!     if(fix)then
	!       write(20,'(2(a6,2x,i3,2x))')'case=',flag,'m=',m
	!     endif

	END
	!*==HYPER.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019


	!***********************************************************************
	!
	!  subroutine name     - hyper
	!
	!  computation
	!  performed           - calculates the hypergeometric function,
	!                        f(a,b;c;w), from its power series for 0<w<.5.
	!
	!  usage               - call hyper(w,a,b,c,f)
	!
	!  arguments
	!                   w  - the transformed independent variable.
	!
	!                a,b,c - the parameters of the hypergeometric function.
	!
	!                   f  - the computed value of the hypergeometric
	!                        function which is returned to the caller.
	!
	!  precision           - double
	!
	!  language            - fortran
	!
	!***********************************************************************

	SUBROUTINE HYPER(W,A,B,C,F,Machep)
	IMPLICIT NONE
	!*--HYPER462
	INTEGER i , m , n , NMAX , k , k0 , k1
	PARAMETER (NMAX=100)
	REAL*8 A , B , C , W , F , alpha0 , alpha1 , rn , term , Machep , &
		& BINom
	COMMON /BCOEFF/ BINom(5151)

	!  compute the number of sums needed to get good convergence

	alpha1 = A + B - C
	k1 = NINT(DBLE(alpha1))

	DO n = 1 , NMAX

		rn = 0.D0
		alpha0 = (A+n+1)*(B+n+1)/(C+n+1) - (n+1)
		k0 = NINT(DBLE(alpha0))
		k = MAX(k0,k1)
		IF ( k<=1 ) k = 1
		IF ( n+k>=100 ) THEN
			WRITE (*,*)                                                 &
				&'error in hyp:  binomial coefficient routine                      &
				&                no longer valid'
			RETURN
		ENDIF

		DO m = 0 , k
			rn = rn + BINom((n+k+1)*(n+k+2)/2+m+1)
		ENDDO

		term = 1.D0
		DO i = 1 , n + 1
			term = (A+i-1)*(B+i-1)/(C+i-1)/(k+i)*term
		ENDDO

		rn = rn*term*(W**(n+1))/(1-W)
		IF ( DABS(rn)<Machep ) GOTO 100

	ENDDO

	WRITE (*,*) 'error in hyp:  nmax not large enough'
	RETURN

	!     write(20,'(2(a6,i3,5x))')'n=',n

	!  evaluate the hypergeometric function of w from its power series

100 term = 1.D0
	F = 1.D0
	DO k = 1 , n
		term = term*(A+k-1)*(B+k-1)/(C+k-1)*W/k
		F = F + term
	ENDDO

	END
	!*==GAMM.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!***********************************************************************
	!
	!  function name      - gamm
	!
	!  computation
	!  performed          - calculates the gamma function
	!
	!  usage              - gamm(x)
	!
	!  argument         x - any real number (excluding negative integers).
	!
	!  precision          - double
	!
	!  language           - fortran
	!
	!***********************************************************************

	FUNCTION GAMM(X)
	IMPLICIT NONE
	!*--GAMM538

	REAL(KIND=DPD) ZERO , ONE
	PARAMETER (ZERO=0.D0,ONE=1.D0)
	REAL(KIND=DPD) X , xx , coeff , GAMM

	!  scale input variable and change it's name
	!  so that it does not get destroyed

	xx = X - ONE
	coeff = ONE

100 IF ( (ZERO<=xx) .AND. (xx<=ONE) ) THEN

		GAMM = coeff*G(xx)
	ELSE

		IF ( xx<ZERO ) THEN
			xx = xx + ONE
			coeff = coeff/xx
		ELSE
			coeff = xx*coeff
			xx = xx - ONE
		ENDIF
		GOTO 100
	ENDIF

	END
	!*==G.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!***********************************************************************
	!
	!  function name     - g
	!
	!  computation
	!  performed         - calculates gamma(xx+1) for xx in the interval
	!                      [0,1] using clenshaw's recurrence formula with
	!                      tchebychev polynomials and the tabulated
	!                      expansion coefficients.
	!
	!  usage             - g(xx)
	!
	!  argument       xx - scaled value for 'x' in 'gamm(x)'.
	!
	!  precision         - double
	!
	!  language          - fortran
	!
	!***********************************************************************

	FUNCTION G(Xx)
	IMPLICIT NONE
	!*--G590
	!*** Start of declarations inserted by SPAG
	INTEGER k
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) ZERO , ONE , TWO
	PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0)
	REAL(KIND=DPD) c(0:41) , Xx , y , y1 , y2 , G

	!  use clenshaw recurrence scheme with tchebychev polynomials
	!  and the expansion coefficients tabulated in 'cheb' for 0<xx<1 .

	CALL CHEB(c,41,1)

	y1 = ZERO
	y2 = ZERO

	DO k = 41 , 1 , -1
		y = TWO*(TWO*Xx-ONE)*y1 - y2 + c(k)
		y2 = y1
		y1 = y
	ENDDO

	G = -y2 + (TWO*Xx-ONE)*y1 + c(0)

	END
	!*==FIX1.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!**********************************************************************
	!
	!  subroutine name    - fix1
	!
	!  computation
	!  performed          - calculates the hypergeometric function for
	!                       z less than -1 when a-b is near an integer.
	!
	!  usage              - call fix1(a,b,c,n,k,f,w,machep,eps)
	!
	!  arguments    a,b,c - parameters of the hypergeometric function.
	!
	!                  n  - the upper limit of the finite series expansion
	!                       of the hypergeometric function.
	!
	!                  k  - equals the nearest integer of a-b.
	!
	!                  f  - computed value of the hypergeometric function.
	!
	!                  w  - transformed independent variable.
	!
	!              machep - equals machine epsilon.
	!
	!                eps  - equals a-b-k.
	!
	!  precision          - double
	!
	!  language           - fortran
	!
	!***********************************************************************

	SUBROUTINE FIX1(A,B,C,N,K,F,W,Machep,Eps)
	IMPLICIT NONE
	!*--FIX1651
	!*** Start of declarations inserted by SPAG
	INTEGER i , K , N , NMAX
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) ZERO , ONE , TWO , FOUR , EIGHTH
	PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FOUR=4.D0,EIGHTH=1.D0/8.D0,&
		& NMAX=100)
	REAL(KIND=DPD) A , B , C , W , F , Eps , Machep , test , arg , rn ,&
		& sum , et1 , et2 , term1 , term2 , term3 , term4 , term5 ,  &
		& term6 , term7 , temp , temp1 , temp2 , temp6 , coeff ,     &
		& coeff1 , coeff2 , coeff3 , coeff4 , x , x1 , x2 , x3 , x4 ,&
		& t1(0:80) , t2(0:80) , t3(0:80) , t4(0:80) , c1(0:80) ,     &
		& c2(0:80) , c3(0:80) , c4(0:80) , f1(0:80) , f2(0:80) ,     &
		& f3(0:80) , f4(0:80) , g1(0:NMAX) , g2(0:NMAX) , g3(0:NMAX) &
		& , g4(0:NMAX) , g5(0:NMAX) , fff1(0:NMAX) , ff1(0:NMAX) ,   &
		& fff2(0:NMAX) , ff2(0:NMAX) , ff3(0:NMAX) , ff4(0:NMAX) ,   &
		& poch1(0:NMAX) , poch2(0:NMAX) , e1(0:NMAX) , e2(0:NMAX) ,  &
		& e3(0:NMAX) , e4(0:NMAX)

	INTEGER flag

	!  calculate the extra terms

	x = B - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
100 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 100
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
150	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 150
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 100
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c1,41,1)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - ONE
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 41
			t1(i) = (FOUR*(x+Eps)-TWO)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-TWO)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c1,55,2)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps)
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 55
			t1(i) = FOUR*(x+Eps)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + FOUR*x*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c1,34,3)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - TWO
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 34
			t1(i) = (FOUR*(x+Eps)-FOUR)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-FOUR)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = B
		coeff1 = ONE
200	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 200
		ENDIF
		x2 = B + Eps
		coeff2 = ONE
250	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 250
		ENDIF
		temp = sum + coeff*temp
		et1 = -temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
300	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 300
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
350	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 350
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et1 = sum + coeff*temp
	ENDIF
	et1 = -et1
	!     write(10,*)et1,(one/gamm(a-dble(k)-eps)-one/gamm(a-dble(k)))
	!    #                                           /eps
	x = C - A + DBLE(K) - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
400 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 400
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
450	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 450
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 400
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c2,41,1)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps) - ONE
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 41
			t2(i) = (FOUR*(x+Eps)-TWO)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + (FOUR*x-TWO)*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c2,55,2)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps)
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 55
			t2(i) = FOUR*(x+Eps)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + FOUR*x*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c2,34,3)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps) - TWO
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 34
			t2(i) = (FOUR*(x+Eps)-FOUR)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + (FOUR*x-FOUR)*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = C - A + DBLE(K)
		coeff1 = ONE
500	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 500
		ENDIF
		x2 = C - A + DBLE(K) + Eps
		coeff2 = ONE
550	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 550
		ENDIF
		temp2 = sum + coeff*temp2
		et2 = -temp2*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
600	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 600
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
650	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 650
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et2 = sum + coeff*temp2
	ENDIF

	!     write(10,*)et2,(one/gamm(c-a+dble(k)+eps)-one/gamm(c-a+dble(k)))
	!    #                                                         /eps

	!  calculate the f-functions

	x1 = A
	coeff1 = ONE
700 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 700
	ENDIF

	x2 = A - DBLE(K)
	coeff2 = ONE
800 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 800
	ENDIF

	x3 = C - A
	coeff3 = ONE
900 IF ( x3<ONE ) THEN
		coeff3 = x3*coeff3
		x3 = x3 + ONE
		GOTO 900
	ENDIF

	x4 = C - A + DBLE(K)
	coeff4 = ONE
1000 IF ( x4<ONE ) THEN
		coeff4 = x4*coeff4
		x4 = x4 + ONE
		GOTO 1000
	ENDIF

	coeff = ONE
	arg = -Eps - DBLE(K)
1100 IF ( arg<-Eps ) THEN
		coeff = coeff/arg
		arg = arg + ONE
		GOTO 1100
	ENDIF

	fff1(0) = ONE
	fff2(0) = ONE
	ff1(0) = ONE
	ff2(0) = ONE
	DO i = 1 , K
		fff1(0) = (B+DBLE(i-1))*fff1(0)
		ff2(0) = (C-A+DBLE(i-1))*ff2(0)
	ENDDO

	fff1(0) = fff1(0)*coeff1/GAMM(x1)
	fff2(0) = coeff3/GAMM(x3)
	ff1(0) = coeff2/GAMM(x2)
	ff2(0) = ff2(0)*coeff4/GAMM(x4)
	ff3(0) = (-1)**(K+1)*GAMM(ONE-Eps)*coeff
	ff4(0) = GAMM(ONE+Eps)

	!     do 26 i=1,n
	!       fff1(i)=(b+dble(k+i-1))*fff1(i-1)
	!       fff2(i)=(c-b+dble(i-1))*fff2(i-1)
	!       ff1(i)=(a+dble(i-1))*ff1(i-1)
	!       ff2(i)=(c-a+dble(k+i-1))*ff2(i-1)
	!       ff3(i)=ff3(i-1)/(eps+dble(i)+dble(k))
	!       ff4(i)=ff4(i-1)/(dble(i)-eps)
	!26   continue

	!     do 27 i=0,n
	!       write(10,*)'fff1=',fff1(i),gamm(b+i+k)/gamm(a)/gamm(b)
	!       write(10,*)'ff1=',ff1(i),gamm(a+i)/gamm(a)/gamm(a-k)
	!       write(10,*)'fff2=',fff2(i),gamm(c-b+i)/gamm(c-a)/gamm(c-b)
	!       write(10,*)'ff2=',ff2(i),gamm(c-a+i+k)/gamm(c-a)/gamm(c-a+k)
	!       write(10,*)'ff3=',ff3(i),(-1)**(k+i)*eps*gamm(-k-i-eps)
	!       write(10,*)'ff4=',ff4(i),(-1)**i*eps*gamm(eps-i)
	!27   continue

	!   calculate  g1,g2

	x1 = A
	coeff1 = ONE
1200 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1200
	ENDIF

	x2 = C - A
	coeff2 = ONE
1300 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 1300
	ENDIF

	g1(0) = ZERO
	g2(0) = ZERO
	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , K
		g1(0) = g1(0)*(A-Eps+DBLE(i-K-1)) - poch1(0)
		poch1(0) = poch1(0)*(A+DBLE(i-K-1))
	ENDDO

	g1(0) = g1(0)*coeff1/GAMM(x1)
	g2(0) = g2(0)*coeff2/GAMM(x2)
	poch1(0) = poch1(0)*coeff1/GAMM(x1)
	poch2(0) = poch2(0)*coeff2/GAMM(x2)
	DO i = 1 , N
		poch1(i) = (A+DBLE(i-1))*poch1(i-1)
		poch2(i) = (C-A+DBLE(K+i-1))*poch2(i-1)
		g1(i) = g1(i-1)*(A-Eps+DBLE(i-1)) - poch1(i-1)
		g2(i) = g2(i-1)*(C-A+Eps+DBLE(K+i-1)) + poch2(i-1)
	ENDDO

	!     do 104 i=0,n
	!       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
	!       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps
	!104  continue

	!  calculate  g3,g4,g5

	x3 = ZERO
	CALL CHEB(c3,55,2)
	t3(0) = ONE
	t3(1) = TWO*(x3-Eps)
	f3(0) = ZERO
	f3(1) = -TWO
	g3(0) = c3(1)*f3(1)

	x4 = ZERO
	CALL CHEB(c4,55,2)
	t4(0) = ONE
	t4(1) = TWO*(x4+Eps)
	f4(0) = ZERO
	f4(1) = TWO
	g4(0) = c4(1)*f4(1)

	DO i = 2 , 55
		t3(i) = FOUR*(x3-Eps)*t3(i-1) - t3(i-2)
		t4(i) = FOUR*(x4+Eps)*t4(i-1) - t4(i-2)
		f3(i) = -FOUR*t3(i-1) + FOUR*x3*f3(i-1) - f3(i-2)
		f4(i) = FOUR*t4(i-1) + FOUR*x4*f4(i-1) - f4(i-2)
		g3(0) = g3(0) + c3(i)*f3(i)
		g4(0) = g4(0) + c4(i)*f4(i)
	ENDDO

	g3(0) = -g3(0)
	DO i = -K , -1
		g3(0) = (g3(0)+ONE/GAMM(DBLE(K+i+2)))/(DBLE(K+i+1)+Eps)
	ENDDO

	test = DABS(Eps*DLOG(W))
	temp = -DLOG(W)
	IF ( Eps<=ZERO ) THEN
		IF ( test>=EIGHTH ) THEN
			temp = (ONE-EXP(test))/Eps
		ELSE
			i = 1
1320	rn = (Eps**(i)*(DLOG(W))**(i+1))/GAMM(DBLE(i+2))
			IF ( DABS(rn)>=Machep ) THEN
				temp = temp - rn
				i = i + 1
				GOTO 1320
			ENDIF
		ENDIF
		g5(0) = W**(A-Eps)*temp
	ELSE
		IF ( test>=EIGHTH ) THEN
			temp = (EXP(test)-ONE)/Eps
		ELSE
			i = 1
1340	rn = (Eps**(i)*(-DLOG(W))**(i+1))/GAMM(DBLE(i+2))
			IF ( DABS(rn)>=Machep ) THEN
				temp = temp + rn
				i = i + 1
				GOTO 1340
			ENDIF
		ENDIF
		g5(0) = (W**A)*temp
	ENDIF

	!     write(10,*)g3(0),(-1)**k*gamm(-k-eps)+one/eps/gamm(dble(k+1))
	!     write(10,*)g4(0),gamm(eps)-one/eps
	!     write(10,*)g5(0),w**(a-eps)/eps-w**a/eps

	e1(0) = ONE/GAMM(DBLE(K+1))
	e2(0) = ONE
	DO i = 1 , N
		e1(i) = e1(i-1)/DBLE(K+i)
		e2(i) = e2(i-1)/DBLE(i)
		g3(i) = (g3(i-1)+e1(i))/(DBLE(K+i)+Eps)
		g4(i) = (g4(i-1)+e2(i))/(DBLE(i)-Eps)
		g5(i) = W*g5(i-1)
	ENDDO

	e1(0) = ONE
	e2(0) = ONE
	e3(0) = ONE
	e4(0) = ONE
	DO i = 1 , K
		e2(0) = (C-A+DBLE(i-1))/DBLE(i)*e2(0)
		e4(0) = e4(0)/DBLE(i)
	ENDDO

	!     do 140 i=1,n
	!       e1(i)=(a+dble(i-1))/dble(i)*e1(i-1)
	!       e2(i)=(c-a+dble(k+i-1))/dble(k+i)*e2(i-1)
	!       e3(i)=e3(i-1)/dble(i)
	!       e4(i)=e4(i-1)/dble(k+i)
	!140  continue

	!   put everything back together again

	term1 = -GAMM(C)*W**A*e3(0)*fff2(0)*ff3(0)*(-1)**K
	term2 = GAMM(C)*W**A*e3(0)*fff1(0)*ff3(0)*(-1)**K
	term3 = GAMM(C)*W**A*e3(0)*fff1(0)*ff2(0)*(-1)**K
	term4 = GAMM(C)*W**A*e4(0)*fff1(0)*ff2(0)*(-1)**K
	term5 = GAMM(C)*e4(0)*fff1(0)*ff2(0)*ff4(0)*(-1)**K
	term6 = GAMM(C)*W**A*e1(0)*fff2(0)*ff3(0)*(-1)**K
	term7 = GAMM(C)*W**(A-Eps)*e2(0)*fff1(0)*ff4(0)*(-1)**K

	temp = g1(0)*term1 + g2(0)*term2 + g3(0)*term3 + g4(0)            &
		& *term4 + g5(0)*term5
	temp1 = term6
	temp2 = term7
	DO i = 1 , N
		term1 = term1*W*(C-B+DBLE(i-1))/(Eps+DBLE(i+K))/DBLE(i)
		term2 = term2*W*(B+DBLE(K+i-1))/(Eps+DBLE(i+K))/DBLE(i)
		term3 = term3*W*(B+DBLE(K+i-1))*(C-A+DBLE(K+i-1))/DBLE(i)
		term4 = term4*W*(B+DBLE(K+i-1))*(C-A+DBLE(K+i-1))/DBLE(K+i)
		term5 = term5*(B+DBLE(K+i-1))*(C-A+DBLE(K+i-1))/DBLE(K+i)      &
			& /(DBLE(i)-Eps)
		term6 = term6*W*(A+DBLE(i-1))/DBLE(i)*(C-B+DBLE(i-1))          &
			& /(Eps+DBLE(i+K))
		term7 = term7*W*(C-A+DBLE(K+i-1))/DBLE(K+i)*(B+DBLE(K+i-1))    &
			& /(DBLE(i)-Eps)
		temp = temp + g1(i)*term1 + g2(i)*term2 + g3(i)*term3 + g4(i)  &
			& *term4 + g5(i)*term5
		temp1 = temp1 + term6
		temp2 = temp2 + term7
	ENDDO

	!  calculate the finite series term

	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , K - 1
		poch1(i) = (B+DBLE(i-1))*poch1(i-1)
		poch2(i) = (C-A+DBLE(i-1))*poch2(i-1)
	ENDDO

	temp6 = ZERO
	DO i = 0 , K - 1
		temp6 = temp6 + poch1(i)*poch2(i)*GAMM(DBLE(K-i)+Eps)          &
			& /GAMM(DBLE(i+1))*(-W)**i
	ENDDO

	x1 = A
	coeff1 = ONE
1400 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1400
	ENDIF

	x2 = C - B
	coeff2 = ONE
1500 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 1500
	ENDIF

	F = temp + et1*temp1 + et2*temp2 + coeff1*coeff2/GAMM(x1)/GAMM(x2)&
		& *GAMM(C)*W**(A-Eps-DBLE(K))*temp6

	!  alternative method  (must also turn on the individual functions)

	!     temp3=zero
	!     temp4=zero
	!     temp5=zero
	!     do 200 i=0,n
	!       term1=-gamm(c)*w**(a+dble(i))*e3(i)*g1(i)*fff2(i)*ff3(i)*(-1)**k
	!       term2=gamm(c)*w**(a+dble(i))*e3(i)*fff1(i)*g2(i)*ff3(i)*(-1)**k
	!       term3=gamm(c)*w**(a+dble(i))*e3(i)*fff1(i)*ff2(i)*g3(i)*(-1)**k
	!       term4=gamm(c)*w**(a+dble(i))*e4(i)*fff1(i)*ff2(i)*g4(i)*(-1)**k
	!       term5=gamm(c)*fff1(i)*e4(i)*ff2(i)*ff4(i)*g5(i)*(-1)**k
	!       temp3=temp3+term1+term2+term3+term4+term5
	!       temp4=temp4+gamm(c)*w**(a+dble(i))*e1(i)*fff2(i)*ff3(i)*(-1)**k
	!       temp5=temp5+gamm(c)*w**(a+dble(i)-eps)*e2(i)*fff1(i)*ff4(i)
	!    #                                                       *(-1)**k
	!200  continue
	!     write(10,*)'temp=',temp,temp3
	!     write(10,*)'temp1=',temp1,temp4
	!     write(10,*)'temp2=',temp2,temp5

	!     x=temp3+et1*temp4+et2*temp5+coeff1*coeff2/gamm(x1)/gamm(x2)
	!    #                             *gamm(c)*w**(a-eps-dble(k))*temp6
	!     write(10,*)'f=',f,x

	END
	!*==FIX4A.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!**********************************************************************
	!
	!  subroutine name    - fix4a
	!
	!  computation
	!  performed          - calculates the hypergeometric function for z
	!                       in the interval (.5,1) when c-a-b is near a
	!                       positive integer.
	!
	!  usage              - call fix4a(a,b,c,n,k,f,w,machep,eps)
	!
	!  arguments    a,b,c - parameters of the hypergeometric function.
	!
	!                  n  - the upper limit of the finite series expansion
	!                       of the hypergeometric function.
	!
	!                  k  - equals the nearest integer of c-a-b.
	!
	!                  f  - computed value of the hypergeometric function.
	!
	!                  w  - transformed independent variable.
	!
	!              machep - equals machine epsilon.
	!
	!                eps  - equals c-a-b-k.
	!
	!  precision          - double
	!
	!  language           - fortran
	!
	!***********************************************************************

	SUBROUTINE FIX4A(A,B,C,N,K,F,W,Machep,Eps)
	IMPLICIT NONE
	!*--FIX4A1213
	!*** Start of declarations inserted by SPAG
	INTEGER i , K , N , NMAX
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) ZERO , ONE , TWO , FOUR , EIGHTH
	PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FOUR=4.D0,EIGHTH=1.D0/8.D0,&
		& NMAX=100)
	REAL(KIND=DPD) A , B , C , W , F , Eps , Machep , test , arg , rn ,&
		& sum , et1 , et2 , term1 , term2 , term3 , term4 , term5 ,  &
		& term6 , term7 , term8 , temp , temp1 , temp2 , temp3 ,     &
		& temp4 , coeff , coeff1 , coeff2 , x , x1 , x2 , x3 , x4 ,  &
		& t1(0:80) , t2(0:80) , t3(0:80) , t4(0:80) , c1(0:80) ,     &
		& c2(0:80) , c3(0:80) , c4(0:80) , f1(0:80) , f2(0:80) ,     &
		& f3(0:80) , f4(0:80) , g1(0:NMAX) , g2(0:NMAX) , g3(0:NMAX) &
		& , g4(0:NMAX) , g5(0:NMAX) , fff1(0:NMAX) , ff1(0:NMAX) ,   &
		& fff2(0:NMAX) , ff2(0:NMAX) , ff3(0:NMAX) , ff4(0:NMAX) ,   &
		& poch1(0:NMAX) , poch2(0:NMAX) , e1(0:NMAX) , e2(0:NMAX) ,  &
		& e3(0:NMAX) , e4(0:NMAX)

	INTEGER flag

	!  calculate the extra terms

	x = A + DBLE(K) - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
100 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 100
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
150	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 150
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 100
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c1,41,1)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - ONE
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 41
			t1(i) = (FOUR*(x+Eps)-TWO)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-TWO)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c1,55,2)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps)
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 55
			t1(i) = FOUR*(x+Eps)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + FOUR*x*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c1,34,3)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - TWO
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 34
			t1(i) = (FOUR*(x+Eps)-FOUR)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-FOUR)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = A + DBLE(K)
		coeff1 = ONE
200	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 200
		ENDIF
		x2 = A + DBLE(K) + Eps
		coeff2 = ONE
250	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 250
		ENDIF
		temp = sum + coeff*temp
		et1 = -temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
300	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 300
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
350	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 350
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et1 = sum + coeff*temp
	ENDIF

	!     write(10,*)et1,(one/gamm(a+k+eps)-one/gamm(a+k))/eps

	x = B + DBLE(K) - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
400 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 400
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
450	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 450
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 400
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c2,41,1)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps) - ONE
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 41
			t2(i) = (FOUR*(x+Eps)-TWO)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + (FOUR*x-TWO)*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c2,55,2)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps)
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 55
			t2(i) = FOUR*(x+Eps)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + FOUR*x*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c2,34,3)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps) - TWO
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 34
			t2(i) = (FOUR*(x+Eps)-FOUR)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + (FOUR*x-FOUR)*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = B + DBLE(K)
		coeff1 = ONE
500	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 500
		ENDIF
		x2 = B + DBLE(K) + Eps
		coeff2 = ONE
550	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 550
		ENDIF
		temp2 = sum + coeff*temp2
		et2 = -temp2*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
600	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 600
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
650	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 650
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et2 = sum + coeff*temp2
	ENDIF

	!     write(10,*)et2,(one/gamm(b+k+eps)-one/gamm(b+k))/eps

	!  calculate the f-functions

	x1 = A
	coeff1 = ONE
700 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 700
	ENDIF

	x2 = B
	coeff2 = ONE
800 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 800
	ENDIF

	coeff = ONE
	arg = -Eps - DBLE(K)
900 IF ( arg<-Eps ) THEN
		coeff = coeff/arg
		arg = arg + ONE
		GOTO 900
	ENDIF

	ff1(0) = coeff1/GAMM(x1)
	ff2(0) = coeff2/GAMM(x2)
	fff1(0) = coeff1/GAMM(x1)
	fff2(0) = coeff2/GAMM(x2)
	ff3(0) = (-1)**(K+1)*GAMM(ONE-Eps)*coeff
	ff4(0) = GAMM(ONE+Eps)

	!     do 23 i=1,n
	!       ff1(i)=(a+dble(k+i-1))*ff1(i-1)
	!       ff2(i)=(b+dble(k+i-1))*ff2(i-1)
	!       fff1(i)=(c-b+dble(i-1))*fff1(i-1)
	!       fff2(i)=(c-a+dble(i-1))*fff2(i-1)
	!       ff3(i)=ff3(i-1)/(eps+dble(i)+dble(k))
	!       ff4(i)=ff4(i-1)/(dble(i)-eps)
	! 23  continue

	!     do 24 i=0,n
	!       write(10,*)'fff1=',fff1(i),gamm(c-b+i)/gamm(a)/gamm(c-b)
	!       write(10,*)'ff1=',ff1(i),gamm(a+k+i)/gamm(a)/gamm(a+k)
	!       write(10,*)'fff2=',fff2(i),gamm(c-a+i)/gamm(b)/gamm(c-a)
	!       write(10,*)'ff2=',ff2(i),gamm(b+k+i)/gamm(b)/gamm(b+k)
	!       write(10,*)'ff3=',ff3(i),(-1)**(k+i)*eps*gamm(-k-i-eps)
	!       write(10,*)'ff4=',ff4(i),(-1)**i*eps*gamm(eps-i)
	! 24  continue

	!   calculate  g1,g2

	g1(0) = ZERO
	g2(0) = ZERO
	poch1(0) = coeff1/GAMM(x1)
	poch2(0) = coeff2/GAMM(x2)
	DO i = 1 , N
		poch1(i) = (A+DBLE(K+i-1))*poch1(i-1)
		poch2(i) = (B+DBLE(K+i-1))*poch2(i-1)
		g1(i) = g1(i-1)*(A+Eps+DBLE(K+i-1)) + poch1(i-1)
		g2(i) = g2(i-1)*(B+Eps+DBLE(K+i-1)) + poch2(i-1)
	ENDDO

	!     do 101 i=0,n
	!       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
	!       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps
	!101  continue

	!  calculate  g3,g4,g5

	x3 = ZERO
	CALL CHEB(c3,55,2)
	t3(0) = ONE
	t3(1) = TWO*(x3-Eps)
	f3(0) = ZERO
	f3(1) = -TWO
	g3(0) = c3(1)*f3(1)

	x4 = ZERO
	CALL CHEB(c4,55,2)
	t4(0) = ONE
	t4(1) = TWO*(x4+Eps)
	f4(0) = ZERO
	f4(1) = TWO
	g4(0) = c4(1)*f4(1)

	DO i = 2 , 55
		t3(i) = FOUR*(x3-Eps)*t3(i-1) - t3(i-2)
		t4(i) = FOUR*(x4+Eps)*t4(i-1) - t4(i-2)
		f3(i) = -FOUR*t3(i-1) + FOUR*x3*f3(i-1) - f3(i-2)
		f4(i) = FOUR*t4(i-1) + FOUR*x4*f4(i-1) - f4(i-2)
		g3(0) = g3(0) + c3(i)*f3(i)
		g4(0) = g4(0) + c4(i)*f4(i)
	ENDDO

	g3(0) = -g3(0)
	DO i = -K , -1
		g3(0) = (g3(0)+ONE/GAMM(DBLE(K+i+2)))/(DBLE(K+i+1)+Eps)
	ENDDO

	test = Eps*DLOG(W)
	temp = DLOG(W)
	IF ( DABS(test)>=EIGHTH ) THEN
		temp = (EXP(test)-ONE)/Eps
	ELSE
		i = 1
950	rn = (Eps**(i)*(DLOG(W))**(i+1))/GAMM(DBLE(i+2))
		IF ( DABS(rn)>=Machep ) THEN
			temp = temp + rn
			i = i + 1
			GOTO 950
		ENDIF
	ENDIF
	g5(0) = (W**K)*temp

	!     write(10,*)g3(0),(-1)**k*gamm(-k-eps)+one/eps/gamm(dble(k+1))
	!     write(10,*)g4(0),gamm(eps)-one/eps
	!     write(10,*)g5(0),w**(k+eps)/eps-w**k/eps

	DO i = 1 , N
		g3(i) = (g3(i-1)+ONE/GAMM(DBLE(K+i+1)))/(DBLE(K+i)+Eps)
		g4(i) = (g4(i-1)+ONE/GAMM(DBLE(i+1)))/(DBLE(i)-Eps)
		g5(i) = W*g5(i-1)
	ENDDO

	e1(0) = ONE
	e2(0) = ONE
	e3(0) = -ONE
	e4(0) = ONE
	DO i = 1 , K
		e1(0) = (A+DBLE(i-1))*e1(0)
		e2(0) = (B+DBLE(i-1))*e2(0)
		e3(0) = e3(0)/DBLE(i)
	ENDDO

	!     do 140 i=1,n
	!       e1(i)=(a+dble(k+i-1))*e1(i-1)
	!       e2(i)=(b+dble(k+i-1))*e2(i-1)
	!       e3(i)=e3(i-1)/dble(k+i)
	!       e4(i)=e4(i-1)/dble(i)
	!140  continue

	!  put everything back together again

	term1 = GAMM(C)*(-1)**K*fff2(0)*ff3(0)*e4(0)*W**(C-A-B)
	term2 = GAMM(C)*(-1)**K*ff1(0)*ff3(0)*e4(0)*W**(C-A-B)
	term3 = GAMM(C)*(-1)**K*ff1(0)*ff2(0)*e4(0)*W**(C-A-B)
	term4 = -GAMM(C)*(-1)**K*ff1(0)*ff2(0)*e3(0)*W**(C-A-B)
	term5 = GAMM(C)*(-1)**K*ff1(0)*ff2(0)*e3(0)*ff4(0)
	term6 = -GAMM(C)*(-W)**K*et1*e1(0)*ff2(0)*e3(0)*ff4(0)
	term7 = -GAMM(C)*(-W)**K*et2*ff1(0)*e2(0)*e3(0)*ff4(0)
	term8 = -GAMM(C)*(-W)**K*Eps*et1*et2*e1(0)*e2(0)*e3(0)*ff4(0)

	temp = g1(0)*term1 + g2(0)*term2 + g3(0)*term3 + g4(0)            &
		& *term4 + g5(0)*term5
	temp1 = term6
	temp2 = term7
	temp3 = term8

	DO i = 1 , N
		term1 = term1*W*(B+Eps+DBLE(K+i-1))/(Eps+DBLE(i+K))/DBLE(i)
		term2 = term2*W*(A+DBLE(K+i-1))/(Eps+DBLE(i+K))/DBLE(i)
		term3 = term3*W*(A+DBLE(K+i-1))*(B+DBLE(K+i-1))/DBLE(i)
		term4 = term4*W*(A+DBLE(K+i-1))*(B+DBLE(K+i-1))/DBLE(K+i)
		term5 = term5*(A+DBLE(K+i-1))*(B+DBLE(K+i-1))/DBLE(K+i)        &
			& /(DBLE(i)-Eps)
		term6 = term6*W*(A+DBLE(K+i-1))*(B+DBLE(K+i-1))/DBLE(K+i)      &
			& /(DBLE(i)-Eps)
		term7 = term7*W*(A+DBLE(K+i-1))*(B+DBLE(K+i-1))/DBLE(K+i)      &
			& /(DBLE(i)-Eps)
		term8 = term8*W*(A+DBLE(K+i-1))*(B+DBLE(K+i-1))/DBLE(K+i)      &
			& /(DBLE(i)-Eps)
		temp = temp + g1(i)*term1 + g2(i)*term2 + g3(i)*term3 + g4(i)  &
			& *term4 + g5(i)*term5
		temp1 = temp1 + term6
		temp2 = temp2 + term7
		temp3 = temp3 + term8
	ENDDO

	!  calculate the finite series term

	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , K - 1
		poch1(i) = (A+DBLE(i-1))*poch1(i-1)
		poch2(i) = (B+DBLE(i-1))*poch2(i-1)
	ENDDO

	temp4 = ZERO
	DO i = 0 , K - 1
		temp4 = temp4 + poch1(i)*poch2(i)*GAMM(Eps+DBLE(K-i))*(-W)     &
			& **i/GAMM(DBLE(i+1))
	ENDDO

	x1 = C - A
	coeff1 = ONE
1000 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1000
	ENDIF

	x2 = C - B
	coeff2 = ONE
1100 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 1100
	ENDIF

	temp4 = temp4*GAMM(C)*coeff1*coeff2/GAMM(x1)/GAMM(x2)

	F = temp + temp1 + temp2 + temp3 + temp4

	!  alternative method  (must also turn on the individual functions)

	!     temp5=zero
	!     temp6=zero
	!     temp7=zero
	!     temp8=zero
	!     do 200 i=0,n
	!       term1=gamm(c)*(-1)**k*e4(i)*g1(i)*fff2(i)*ff3(i)
	!    #                                     *w**(dble(k+i)+eps)
	!       term2=gamm(c)*(-1)**k*e4(i)*g2(i)*ff1(i)*ff3(i)
	!    #                                     *w**(dble(k+i)+eps)
	!       term3=gamm(c)*(-1)**k*e4(i)*g3(i)*ff1(i)*ff2(i)
	!    #                                     *w**(dble(k+i)+eps)
	!       term4=-gamm(c)*(-1)**k*e3(i)*g4(i)*ff1(i)*ff2(i)
	!    #                                     *w**(dble(k+i)+eps)
	!       term5=gamm(c)*(-1)**k*e3(i)*g5(i)*ff1(i)*ff2(i)*ff4(i)

	!       temp5=temp5+term1+term2+term3+term4+term5
	!       temp6=temp6-gamm(c)*e3(i)*(-1)**k*et1*e1(i)*ff2(i)*ff4(i)
	!    #                                            *w**(k+i)
	!       temp7=temp7-gamm(c)*e3(i)*(-1)**k*et2*e2(i)*ff1(i)*ff4(i)
	!    #                                            *w**(k+i)
	!       temp8=temp8-gamm(c)*e3(i)*(-1)**k*eps*et1*et2*e1(i)*e2(i)*ff4(i)
	!    #                                            *w**(k+i)
	!200  continue
	!     write(10,*)'temp=',temp,temp5
	!     write(10,*)'temp1=',temp1,temp6
	!     write(10,*)'temp2=',temp2,temp7
	!     write(10,*)'temp3=',temp3,temp8

	!     x=temp5+temp6+temp7+temp8+temp4
	!     write(10,*)'f=',f,x

	END
	!*==FIX4B.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!**********************************************************************
	!
	!  subroutine name    - fix4b
	!
	!  computation
	!  performed          - calculates the hypergeometric function for z
	!                       in the interval (.5,1) when c-a-b is near a
	!                       negative integer.
	!
	!  usage              - call fix4b(a,b,c,n,k,f,w,machep,eps)
	!
	!  arguments    a,b,c - parameters of the hypergeometric function.
	!
	!                  n  - the upper limit of the finite series expansion
	!                       of the hypergeometric function.
	!
	!                  k  - equals the nearest integer of a+b-c.
	!
	!                  f  - computed value of the hypergeometric function.
	!
	!                  w  - transformed independent variable.
	!
	!              machep - equals machine epsilon.
	!
	!                eps  - equals c-a-b+k.
	!
	!  precision          - double
	!
	!  language           - fortran
	!
	!***********************************************************************

	SUBROUTINE FIX4B(A,B,C,N,K,F,W,Machep,Eps)
	IMPLICIT NONE
	!*--FIX4B1722
	!*** Start of declarations inserted by SPAG
	INTEGER i , j , K , N , NMAX
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) ZERO , ONE , TWO , FOUR , EIGHTH
	PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FOUR=4.D0,EIGHTH=1.D0/8.D0,&
		& NMAX=100)
	REAL(KIND=DPD) A , B , C , W , F , Eps , Machep , test , arg , rn ,&
		& sum , et1 , et2 , term1 , term2 , term3 , term4 , term5 ,  &
		& term6 , term7 , term8 , temp , temp1 , temp2 , temp3 ,     &
		& temp4 , coeff , coeff1 , coeff2 , coeff3 , coeff4 , x ,    &
		& x1 , x2 , x3 , x4 , t1(0:80) , t2(0:80) , t3(0:80) ,       &
		& t4(0:80) , c1(0:80) , c2(0:80) , c3(0:80) , c4(0:80) ,     &
		& f1(0:80) , f2(0:80) , f3(0:80) , f4(0:80) , g1(0:NMAX) ,   &
		& g2(0:NMAX) , g3(0:NMAX) , g4(0:NMAX) , g5(0:NMAX) ,        &
		& fff1(0:NMAX) , ff1(0:NMAX) , fff2(0:NMAX) , ff2(0:NMAX) ,  &
		& ff3(0:NMAX) , ff4(0:NMAX) , poch1(0:NMAX) , poch2(0:NMAX) ,&
		& e1(0:NMAX) , e2(0:NMAX) , e3(0:NMAX) , e4(0:NMAX)

	INTEGER flag

	!  calculate the extra terms

	x = A - DBLE(K) - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
100 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 100
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
150	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 150
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 100
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c1,41,1)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - ONE
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 41
			t1(i) = (FOUR*(x+Eps)-TWO)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-TWO)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c1,55,2)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps)
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 55
			t1(i) = FOUR*(x+Eps)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + FOUR*x*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c1,34,3)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - TWO
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 34
			t1(i) = (FOUR*(x+Eps)-FOUR)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-FOUR)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = A - DBLE(K)
		coeff1 = ONE
200	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 200
		ENDIF
		x2 = A - DBLE(K) + Eps
		coeff2 = ONE
250	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 250
		ENDIF
		temp = sum + coeff*temp
		et1 = -temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
300	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 300
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
350	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 350
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et1 = sum + coeff*temp
	ENDIF

	!     write(10,*)et1,(one/gamm(a-k+eps)-one/gamm(a-k))/eps

	x = B - DBLE(K) - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
400 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 400
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
450	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 450
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 400
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c2,41,1)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps) - ONE
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 41
			t2(i) = (FOUR*(x+Eps)-TWO)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + (FOUR*x-TWO)*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c2,55,2)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps)
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 55
			t2(i) = FOUR*(x+Eps)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + FOUR*x*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c2,34,3)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps) - TWO
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 34
			t2(i) = (FOUR*(x+Eps)-FOUR)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + (FOUR*x-FOUR)*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = B - DBLE(K)
		coeff1 = ONE
500	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 500
		ENDIF
		x2 = B - DBLE(K) + Eps
		coeff2 = ONE
550	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 550
		ENDIF
		temp2 = sum + coeff*temp2
		et2 = -temp2*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
600	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 600
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
650	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 650
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et2 = sum + coeff*temp2
	ENDIF

	!     write(10,*)et2,(one/gamm(b-k+eps)-one/gamm(b-k))/eps

	!  calculate the f-functions

	x1 = A
	coeff1 = ONE
700 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 700
	ENDIF

	x2 = B
	coeff2 = ONE
800 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 800
	ENDIF

	x3 = A - DBLE(K)
	coeff3 = ONE
900 IF ( x3<ONE ) THEN
		coeff3 = x3*coeff3
		x3 = x3 + ONE
		GOTO 900
	ENDIF

	x4 = B - DBLE(K)
	coeff4 = ONE
1000 IF ( x4<ONE ) THEN
		coeff4 = x4*coeff4
		x4 = x4 + ONE
		GOTO 1000
	ENDIF

	coeff = ONE
	arg = Eps - DBLE(K)
1100 IF ( arg<Eps ) THEN
		coeff = coeff/arg
		arg = arg + ONE
		GOTO 1100
	ENDIF

	fff1(0) = ONE
	fff2(0) = ONE
	ff1(0) = ONE
	ff2(0) = ONE
	DO i = 1 , K
		fff1(0) = (C-B+DBLE(i-1))*fff1(0)
		fff2(0) = (C-A+DBLE(i-1))*fff2(0)
	ENDDO

	fff1(0) = fff1(0)*coeff1/GAMM(x1)
	fff2(0) = fff2(0)*coeff2/GAMM(x2)
	ff1(0) = ff1(0)*coeff3/GAMM(x3)
	ff2(0) = ff2(0)*coeff4/GAMM(x4)
	ff3(0) = -GAMM(ONE-Eps)
	ff4(0) = (-1)**K*GAMM(ONE+Eps)*coeff

	!     do 26 i=1,n
	!       fff1(i)=(c-b+dble(k+i-1))*fff1(i-1)
	!       fff2(i)=(c-a+dble(k+i-1))*fff2(i-1)
	!       ff1(i)=(a+dble(i-1))*ff1(i-1)
	!       ff2(i)=(b+dble(i-1))*ff2(i-1)
	!       ff3(i)=ff3(i-1)/(eps+dble(i))
	!       ff4(i)=ff4(i-1)/(dble(k+i)-eps)
	! 26  continue

	!     do 27 i=0,n
	!       write(10,*)'fff1=',fff1(i),gamm(a+eps+i)/gamm(a)/gamm(c-b)
	!       write(10,*)'fff2=',fff2(i),gamm(b+eps+i)/gamm(b)/gamm(c-a)
	!       write(10,*)'ff1=',ff1(i),gamm(a+i)/gamm(a)/gamm(a-k)
	!       write(10,*)'ff2=',ff2(i),gamm(b+i)/gamm(b)/gamm(b-k)
	!       write(10,*)'ff3=',ff3(i),(-1)**i*eps*gamm(-i-eps)
	!       write(10,*)'ff4=',ff4(i),(-1)**(k+i)*eps*gamm(eps-k-i)
	! 27  continue

	!   calculate  g1,g2

	x1 = A
	coeff1 = ONE
1200 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1200
	ENDIF

	x2 = B
	coeff2 = ONE
1300 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 1300
	ENDIF

	g1(0) = ZERO
	g2(0) = ZERO
	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , K
		g1(0) = g1(0)*(A+Eps+DBLE(i-K-1)) + poch1(0)
		g2(0) = g2(0)*(B+Eps+DBLE(i-K-1)) + poch2(0)
		poch1(0) = poch1(0)*(A+DBLE(i-K-1))
		poch2(0) = poch2(0)*(B+DBLE(i-K-1))
	ENDDO

	g1(0) = g1(0)*coeff1/GAMM(x1)
	g2(0) = g2(0)*coeff2/GAMM(x2)
	poch1(0) = poch1(0)*coeff1/GAMM(x1)
	poch2(0) = poch2(0)*coeff2/GAMM(x2)
	DO i = 1 , N
		poch1(i) = (A+i-1)*poch1(i-1)/i
		poch2(i) = (B+i-1)*poch2(i-1)/i
		g1(i) = (g1(i-1)*(A+Eps+i-1)+poch1(i-1))/i
		g2(i) = (g2(i-1)*(B+Eps+i-1)+poch2(i-1))/i
	ENDDO

	!     do 104 i=0,n
	!       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps/gamma(i+1.0)
	!       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps/gamma(i+1.0)
	!104  continue

	!  calculate  g3,g4,g5

	x3 = ZERO
	CALL CHEB(c3,55,2)
	t3(0) = ONE
	t3(1) = TWO*(x3-Eps)
	f3(0) = ZERO
	f3(1) = -TWO
	g3(0) = c3(1)*f3(1)

	x4 = ZERO
	CALL CHEB(c4,55,2)
	t4(0) = ONE
	t4(1) = TWO*(x4+Eps)
	f4(0) = ZERO
	f4(1) = TWO
	g4(0) = c4(1)*f4(1)

	DO i = 2 , 55
		t3(i) = FOUR*(x3-Eps)*t3(i-1) - t3(i-2)
		t4(i) = FOUR*(x4+Eps)*t4(i-1) - t4(i-2)
		f3(i) = -FOUR*t3(i-1) + FOUR*x3*f3(i-1) - f3(i-2)
		f4(i) = FOUR*t4(i-1) + FOUR*x4*f4(i-1) - f4(i-2)
		g3(0) = g3(0) + c3(i)*f3(i)
		g4(0) = g4(0) + c4(i)*f4(i)
	ENDDO

	g3(0) = -g3(0)
	DO i = -K , -1
		g4(0) = (g4(0)+ONE/GAMM(DBLE(K+i+2)))/(DBLE(K+i+1)-Eps)
	ENDDO

	test = Eps*DLOG(W)
	temp = DLOG(W)
	IF ( DABS(test)>=EIGHTH ) THEN
		temp = (EXP(test)-ONE)/Eps
	ELSE
		i = 1
1350 rn = (Eps**(i)*(DLOG(W))**(i+1))/GAMM(DBLE(i+2))
		IF ( DABS(rn)>=Machep ) THEN
			temp = temp + rn
			i = i + 1
			GOTO 1350
		ENDIF
	ENDIF
	g5(0) = temp

	!     write(10,*)g3(0),gamm(-eps)+one/eps
	!     write(10,*)g4(0),(-1)**k*gamm(eps-dble(k))-one/eps/gamm(dble(k+1))
	!     write(10,*)g5(0),w**eps/eps-one/eps

	DO i = 1 , N
		temp = ONE/GAMM(DBLE(K+1))
		DO j = 1 , i
			temp = temp*DBLE(j)/DBLE(K+j)
		ENDDO
		g3(i) = (g3(i-1)*DBLE(i)+ONE)/(DBLE(i)+Eps)
		g4(i) = (g4(i-1)*DBLE(i)+temp)/(DBLE(K+i)-Eps)
		g5(i) = W*g5(i-1)
	ENDDO

	e1(0) = ONE
	e2(0) = ONE
	e3(0) = -ONE
	e4(0) = ONE
	DO i = 1 , K
		e4(0) = e4(0)/DBLE(i)
	ENDDO

	!     do 140 i=1,n
	!       e1(i)=(a+dble(i-1))*e1(i-1)
	!       e2(i)=(b+dble(i-1))*e2(i-1)
	!       e3(i)=e3(i-1)/dble(i)
	!       e4(i)=e4(i-1)/dble(k+i)
	!140  continue

	!  put everything back together again

	term1 = GAMM(C)*(-1)**K*fff2(0)*ff3(0)*e4(0)*W**Eps
	term2 = GAMM(C)*(-1)**K*ff1(0)*ff3(0)*e4(0)*W**Eps
	term3 = GAMM(C)*(-1)**K*ff1(0)*ff2(0)*e4(0)*W**Eps
	term4 = -GAMM(C)*(-1)**K*ff1(0)*ff2(0)*e3(0)*W**Eps
	term5 = GAMM(C)*(-1)**K*ff1(0)*ff2(0)*e3(0)*ff4(0)
	term6 = -GAMM(C)*(-1)**K*et1*e1(0)*ff2(0)*e3(0)*ff4(0)
	term7 = -GAMM(C)*(-1)**K*et2*ff1(0)*e2(0)*e3(0)*ff4(0)
	term8 = -GAMM(C)*(-1)**K*Eps*et1*et2*e1(0)*e2(0)*e3(0)*ff4(0)

	temp = g1(0)*term1 + g2(0)*term2 + g3(0)*term3 + g4(0)            &
		& *term4 + g5(0)*term5
	temp1 = term6
	temp2 = term7
	temp3 = term8

	DO i = 1 , N
		term1 = term1*W*(B+Eps+DBLE(i-1))/(Eps+DBLE(i))*DBLE(i)        &
			& /DBLE(K+i)
		term2 = term2*W*(A+DBLE(i-1))/(Eps+DBLE(i))*DBLE(i)/DBLE(K+i)
		term3 = term3*W*(A+DBLE(i-1))/DBLE(i)*(B+DBLE(i-1))/DBLE(K+i)
		term4 = term4*W*(A+DBLE(i-1))/DBLE(i)*(B+DBLE(i-1))/DBLE(i)
		term5 = term5*(A+DBLE(i-1))/DBLE(i)*(B+DBLE(i-1))              &
			& /(DBLE(K+i)-Eps)
		term6 = term6*W*(A+DBLE(i-1))/DBLE(i)*(B+DBLE(i-1))            &
			& /(DBLE(K+i)-Eps)
		term7 = term7*W*(A+DBLE(i-1))/DBLE(i)*(B+DBLE(i-1))            &
			& /(DBLE(K+i)-Eps)
		term8 = term8*W*(A+DBLE(i-1))/DBLE(i)*(B+DBLE(i-1))            &
			& /(DBLE(K+i)-Eps)
		temp = temp + g1(i)*term1 + g2(i)*term2 + g3(i)*term3 + g4(i)  &
			& *term4 + g5(i)*term5
		temp1 = temp1 + term6
		temp2 = temp2 + term7
		temp3 = temp3 + term8
	ENDDO

	!  calculate the finite series term

	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , K - 1
		poch1(i) = (C-A+DBLE(i-1))*poch1(i-1)
		poch2(i) = (C-B+DBLE(i-1))*poch2(i-1)
	ENDDO

	temp4 = ZERO
	DO i = 0 , K - 1
		temp4 = temp4 + poch1(i)*poch2(i)*GAMM(-Eps+DBLE(K-i))*(-1)    &
			& **i*W**(Eps+DBLE(i-K))/GAMM(DBLE(i+1))
	ENDDO

	x1 = A
	coeff1 = ONE
1400 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1400
	ENDIF

	x2 = B
	coeff2 = ONE
1500 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 1500
	ENDIF

	temp4 = temp4*GAMM(C)*coeff1*coeff2/GAMM(x1)/GAMM(x2)

	F = temp + temp1 + temp2 + temp3 + temp4

	!  alternative method (must also turn on the individual functions)

	!     temp5=zero
	!     temp6=zero
	!     temp7=zero
	!     temp8=zero
	!     do 200 i=0,n
	!       term1=w**(dble(i)+eps)/gamm(dble(k+i+1))
	!    #                      *g1(i)*fff2(i)*ff3(i)
	!       term2=w**(dble(i)+eps)/gamm(dble(k+i+1))
	!    #                      *g2(i)*ff1(i)*ff3(i)
	!       term3=w**(dble(i)+eps)/gamm(dble(k+i+1))
	!    #                      *g3(i)*ff1(i)*ff2(i)
	!       term4=w**(dble(i)+eps)/gamm(dble(i+1))
	!    #                      *g4(i)*ff1(i)*ff2(i)
	!       term5=-ff1(i)/gamm(dble(i+1))*ff2(i)*ff4(i)*g5(i)
	!
	!       temp5=temp5+term1+term2+term3+term4+term5
	!
	!       temp6=temp6+ff1(i)*et2/gamm(dble(i+1))*ff4(i)*e2(i)
	!    #                                 *w**(dble(i))
	!       temp7=temp7+ff2(i)*et1/gamm(dble(i+1))*ff4(i)*e1(i)
	!    #                                 *w**(dble(i))
	!       temp8=temp8+e1(i)*et1*et2*eps/gamm(dble(i+1))*ff4(i)
	!    #                                 *w**(dble(i))*e2(i)
	!200  continue
	!     write(10,*)'temp=',temp,temp5*gamm(c)*(-1)**k
	!     write(10,*)'temp1=',temp1,temp7*gamm(c)*(-1)**k
	!     write(10,*)'temp2=',temp2,temp6*gamm(c)*(-1)**k
	!     write(10,*)'temp3=',temp3,temp8*gamm(c)*(-1)**k
	!
	!     x=(-1)**k*gamm(c)*(temp5+temp6+temp7+temp8)+temp4
	!
	!     write(10,*)'f=',f,x

	END
	!*==FIX5A.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!**********************************************************************
	!
	!  subroutine name    - fix5a
	!
	!  computation
	!  performed          - calculates the hypergeometric function for z
	!                       in the interval (1,2) when c-a-b is near a
	!                       positive integer.
	!
	!  usage              - call fix5a(a,b,c,n,k,re,im,w,machep,eps,pi)
	!
	!  arguments    a,b,c - parameters of the hypergeometric function.
	!
	!                  n  - the upper limit of the finite series expansion
	!                       of the hypergeometric function.
	!
	!                  k  - equals the nearest integer of c-a-b.
	!
	!               re,im - computed values for the real and imaginary parts
	!                       of the hypergeometric function.
	!
	!                  w  - transformed independent variable.
	!
	!              machep - equals machine epsilon.
	!
	!                eps  - equals c-a-b-k.
	!
	!                 pi  - equals 3.1415... to machine accuracy.
	!
	!  precision          - double
	!
	!  language           - fortran
	!
	!***********************************************************************

	SUBROUTINE FIX5A(A,B,C,N,K,Re,Im,W,Machep,Eps,Pi)
	IMPLICIT NONE
	!*--FIX5A2291
	!*** Start of declarations inserted by SPAG
	REAL(KIND=DPD) flag
	INTEGER i , K , N , NMAX
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) ZERO , ONE , TWO , FOUR , EIGHTH
	PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FOUR=4.D0,EIGHTH=1.D0/8.D0,&
		& NMAX=100)
	REAL(KIND=DPD) A , B , C , W , Re , Im , temp , temp2 , g1(0:NMAX) &
		& , g2 , g3(0:NMAX) , g4(0:NMAX) , g5(0:NMAX) , x , x1 , x2 ,&
		& x3 , x4 , rn , t1(0:80) , t2(0:80) , t3(0:80) , t4(0:80) , &
		& test , Machep , Pi , f1(0:80) , f2(0:80) , f3(0:80) ,      &
		& f4(0:80) , ff3(0:NMAX) , Eps , ff4(0:NMAX) , coeff1 ,      &
		& coeff2 , c1(0:80) , c2(0:80) , c3(0:80) , c4(0:80) , sum , &
		& term1 , term2 , term3 , term4 , term5 , poch1(0:NMAX) ,    &
		& coeff , temp1 , et1 , et2 , e1(0:NMAX) , e2(0:NMAX) ,      &
		& e3(0:NMAX) , ff1(0:NMAX) , fff1(0:NMAX) , coeff3 , coeff4 ,&
		& f(0:NMAX) , error , poch2(0:NMAX)

	x3 = ZERO
	CALL CHEB(c3,55,2)
	t3(0) = ONE
	t3(1) = TWO*(x3+Eps)
	f3(0) = ZERO
	f3(1) = TWO
	g3(0) = c3(1)*f3(1)

	x4 = ZERO
	CALL CHEB(c4,55,2)
	t4(0) = ONE
	t4(1) = TWO*(x4-Eps)
	f4(0) = ZERO
	f4(1) = -TWO
	g4(0) = c4(1)*f4(1)

	DO i = 2 , 55
		t3(i) = FOUR*(x3+Eps)*t3(i-1) - t3(i-2)
		t4(i) = FOUR*(x4-Eps)*t4(i-1) - t4(i-2)
		f3(i) = FOUR*t3(i-1) + FOUR*x3*f3(i-1) - f3(i-2)
		f4(i) = -FOUR*t4(i-1) + FOUR*x4*f4(i-1) - f4(i-2)
		g3(0) = g3(0) + c3(i)*f3(i)
		g4(0) = g4(0) + c4(i)*f4(i)
	ENDDO

	g4(0) = -g4(0)
	DO i = -K , -1
		g4(0) = (g4(0)+ONE/GAMM(DBLE(K+i+2)))/(DBLE(K+i+1)+Eps)
	ENDDO

	test = Eps*DLOG(W)
	temp = DLOG(W)
	IF ( DABS(test)>=EIGHTH ) THEN
		temp = (EXP(test)-ONE)/Eps
	ELSE
		i = 1
50	rn = (Eps**(i)*(DLOG(W))**(i+1))/GAMM(DBLE(i+2))
		IF ( DABS(rn)>=Machep ) THEN
			temp = temp + rn
			i = i + 1
			GOTO 50
		ENDIF
	ENDIF
	g5(0) = temp*W**K

	!     write(10,*)g3(0),gamm(-eps)+one/eps
	!     write(10,*)g4(0),(-1)**k*gamm(eps-dble(k))-one/eps/gamm(dble(k+1))
	!     write(10,*)g5(0),w**eps/eps-one/eps

	DO i = 1 , N
		g3(i) = (g3(i-1)+ONE/GAMM(DBLE(i+1)))/(DBLE(i)-Eps)
		g4(i) = (g4(i-1)+ONE/GAMM(DBLE(K+i+1)))/(DBLE(K+i)+Eps)
		g5(i) = W*g5(i-1)
	ENDDO

	DO i = 0 , N
		ff3(i) = Eps*g3(i) + ONE/GAMM(DBLE(i+1))
		ff4(i) = Eps*g4(i) - ONE/GAMM(DBLE(K+i+1))
	ENDDO

	!  calculate the extra terms

	x = A - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
100 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 100
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
150	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 150
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 100
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c1,41,1)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - ONE
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 41
			t1(i) = (FOUR*(x+Eps)-TWO)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-TWO)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c1,55,2)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps)
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 55
			t1(i) = FOUR*(x+Eps)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + FOUR*x*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c1,34,3)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - TWO
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 34
			t1(i) = (FOUR*(x+Eps)-FOUR)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-FOUR)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = A
		coeff1 = ONE
200	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 200
		ENDIF
		x2 = A + Eps
		coeff2 = ONE
250	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 250
		ENDIF
		temp = sum + coeff*temp
		et1 = -temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
300	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 300
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
350	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 350
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et1 = sum + coeff*temp
	ENDIF
	et1 = -et1
	!     write(10,*)et1,(one/gamm(c-b-dble(k)-eps)-one/gamm(c-b-dble(k)))
	!    #                                                  /eps

	x = B - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
400 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 400
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
450	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 450
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 400
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c2,41,1)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps) - ONE
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 41
			t2(i) = (FOUR*(x+Eps)-TWO)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + (FOUR*x-TWO)*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c2,55,2)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps)
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 55
			t2(i) = FOUR*(x+Eps)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + FOUR*x*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c2,34,3)
		t2(0) = ONE
		t2(1) = TWO*(x+Eps) - TWO
		f2(0) = ZERO
		f2(1) = TWO
		temp2 = c2(1)*f2(1)
		DO i = 2 , 34
			t2(i) = (FOUR*(x+Eps)-FOUR)*t2(i-1) - t2(i-2)
			f2(i) = FOUR*t2(i-1) + (FOUR*x-FOUR)*f2(i-1) - f2(i-2)
			temp2 = temp2 + c2(i)*f2(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = B
		coeff1 = ONE
500	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 500
		ENDIF
		x2 = B + Eps
		coeff2 = ONE
550	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 550
		ENDIF
		temp2 = sum + coeff*temp2
		et2 = -temp2*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
600	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 600
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
650	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 650
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et2 = sum + coeff*temp2
	ENDIF

	!     write(10,*)et2,(one/gamm(b+eps)-one/gamm(b))/eps

	fff1(0) = ONE
	DO i = 1 , K
		fff1(0) = (A+DBLE(i-1))*fff1(0)
	ENDDO

	ff1(0) = ONE
	DO i = 1 , N
		fff1(i) = (A+DBLE(K+i-1))*fff1(i-1)
		ff1(i) = (C-B+DBLE(i-1))*ff1(i-1)
	ENDDO

	x1 = C - B
	coeff1 = ONE
700 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 700
	ENDIF

	x2 = C - B - DBLE(K)
	coeff2 = ONE
800 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 800
	ENDIF

	DO i = 0 , N
		x3 = B + Eps - DBLE(i)
		coeff3 = ONE
850	IF ( x3<ONE ) THEN
			coeff3 = x3*coeff3
			x3 = x3 + ONE
			GOTO 850
		ENDIF

		x4 = B - DBLE(i)
		coeff4 = ONE
900	IF ( x4<ONE ) THEN
			coeff4 = x4*coeff4
			x4 = x4 + ONE
			GOTO 900
		ENDIF
		f(i) = ff1(i)*coeff4/GAMM(x4)
		fff1(i) = fff1(i)*coeff1*coeff3/GAMM(x1)/GAMM(x3)
		ff1(i) = ff1(i)*coeff2*coeff4/GAMM(x2)/GAMM(x4)
		!       write(10,*)'fff1=',fff1(i),gamm(c-b-eps+dble(i))
		!    #            /gamm(a)/gamm(b+eps-dble(i))/gamm(c-b)
		!       write(10,*)'ff1=',ff1(i),gamm(c-b+dble(i))
		!    #            /gamm(a+eps)/gamm(b-dble(i))/gamm(c-b)
	ENDDO

	!   calculate  g1

	e1(0) = ZERO
	poch1(0) = ONE
	DO i = 1 , K
		e1(0) = e1(0)*(C-B-Eps+DBLE(i-K-1)) - poch1(0)
		poch1(0) = poch1(0)*(C-B+DBLE(i-K-1))
	ENDDO
	DO i = 1 , N
		poch1(i) = (C-B+DBLE(i-1))*poch1(i-1)
		e1(i) = e1(i-1)*(C-B-Eps+DBLE(i-1)) - poch1(i-1)
	ENDDO

	DO i = 0 , N
		x1 = B - DBLE(i)
		coeff1 = ONE
950	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 950
		ENDIF
		e1(i) = e1(i)*coeff1/GAMM(x1)
	ENDDO

	e2(0) = et2
	DO i = 1 , N
		x1 = B - DBLE(i-1)
		coeff1 = ONE
1000 IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 1000
		ENDIF
		e2(i) = e2(i-1)*(B+Eps-DBLE(i)) + coeff1/GAMM(x1)
	ENDDO

	e3(0) = ONE
	DO i = 1 , K
		e3(0) = (A+DBLE(i-1))*e3(0)
	ENDDO

	DO i = 1 , N
		e3(i) = (A+DBLE(K+i-1))*e3(i-1)
	ENDDO

	x1 = C - B
	coeff1 = ONE
1100 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1100
	ENDIF

	DO i = 0 , N
		g1(i) = (e2(i)*e3(i)+e1(i))*coeff1/GAMM(x1)
		!       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
	ENDDO

	!  calculate g2

	g2 = ZERO
	IF ( DABS(Eps)<.1D0 ) THEN
		i = 1
1150 rn = (-1)**i*Pi**(i+i)*Eps**(i+i-1)/GAMM(DBLE(i+i+1))
		IF ( DABS(rn)>=Machep ) THEN
			g2 = g2 + rn
			i = i + 1
			GOTO 1150
		ENDIF
	ELSE
		g2 = (COS(Pi*Eps)-ONE)/Eps
	ENDIF
	!     write(10,*)'g2=',g2,(cos(pi*eps)-one)/eps

	temp = ZERO
	temp1 = ZERO
	DO i = 0 , N
		term1 = -g1(i)*COS(Pi*Eps)/GAMM(DBLE(i+1))*ff4(i)              &
			& *W**(Eps+DBLE(i+K))
		term2 = fff1(i)*g2/GAMM(DBLE(i+1))*ff4(i)*W**(Eps+DBLE(i+K))
		term3 = -fff1(i)*g3(i)*ff4(i)*W**(Eps+DBLE(i+K))
		term4 = fff1(i)*ff3(i)*g4(i)*W**(Eps+DBLE(i+K))
		term5 = -fff1(i)*ff3(i)/GAMM(DBLE(K+i+1))*g5(i)
		temp = temp + (term1+term2+term3+term4+term5)*(-1)**i
		temp1 = temp1 + (et1*f(i)*COS(Pi*Eps)/GAMM(DBLE(i+1))*ff4(i)   &
			& *W**(DBLE(i+K)+Eps))*(-1)**i
	ENDDO

	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , K - 1
		poch1(i) = (A+DBLE(i-1))*poch1(i-1)
		poch2(i) = (ONE-C+A+DBLE(i-1))*poch2(i-1)
	ENDDO

	x1 = C - A
	coeff1 = ONE
1200 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1200
	ENDIF

	x2 = C - B
	coeff2 = ONE
1300 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 1300
	ENDIF

	temp2 = ZERO
	DO i = 0 , K - 1
		temp2 = temp2 + poch1(i)*poch2(i)*coeff1*coeff2/GAMM(x1)       &
			& /GAMM(x2)*GAMM(DBLE(K-i)+Eps)/GAMM(DBLE(i+1))*(-W)**i
	ENDDO

	!     term1=zero
	!     do 81 i=0,k-1
	!       term1=term1+gamm(a+dble(i))/gamm(a)*gamm(a-c+dble(1+i))
	!    #        /gamm(a-c+one)*gamm(eps+dble(k-i))*(-w)**i/gamm(dble(i+1))
	!    #        /gamm(c-a)/gamm(c-b)
	! 81  continue
	!     write(10,*)temp2,term1

	Re = (ONE-W)**A*GAMM(C)*(temp+temp1+temp2)

	!  calculate the imaginary part

	x1 = A
	coeff1 = ONE
1400 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1400
	ENDIF

	temp = ZERO
	DO i = 0 , N
		temp = temp + (-1)**i*f(i)/GAMM(DBLE(i+1))*ff4(i)              &
			& *W**(Eps+DBLE(i+K))*coeff1/GAMM(x1)
	ENDDO

	IF ( DABS(Eps)<.1D0 ) THEN
		temp1 = ONE
		i = 1
1450 temp2 = temp1 + (-1)**i*(Pi*Eps)**(i+i)/GAMM(DBLE(i+i+2))
		error = (temp2-temp1)/temp2
		IF ( DABS(error)>=Machep ) THEN
			i = i + 1
			temp1 = temp2
			GOTO 1450
		ENDIF
	ELSE
		temp2 = SIN(Pi*Eps)/Pi/Eps
	ENDIF
	!     write(10,*)temp2,sin(pi*eps)/pi/eps

	Im = -Pi*temp2*temp

	END
	!*==FIX5B.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!**********************************************************************
	!
	!  subroutine name    - fix5b
	!
	!  computation
	!  performed          - calculates the hypergeometric function for z
	!                       in the interval (1,2) when c-a-b is near a
	!                       negative integer.
	!
	!  usage              - call fix5b(a,b,c,n,k,re,im,w,machep,eps,pi)
	!
	!  arguments    a,b,c - parameters of the hypergeometric function.
	!
	!                  n  - the upper limit of the finite series expansion
	!                       of the hypergeometric function.
	!
	!                  k  - equals the nearest integer of a+b-c.
	!
	!               re,im - computed values for the real and imaginary parts
	!                       of the hypergeometric function.
	!
	!                  w  - transformed independent variable.
	!
	!              machep - equals machine epsilon.
	!
	!                eps  - equals c-a-b+k.
	!
	!                 pi  - equals 3.1415... to machine accuracy.
	!
	!  precision          - double
	!
	!  language           - fortran
	!
	!***********************************************************************

	SUBROUTINE FIX5B(A,B,C,N,K,Re,Im,W,Machep,Eps,Pi)
	IMPLICIT NONE
	!*--FIX5B2826
	!*** Start of declarations inserted by SPAG
	INTEGER i , K , N , NMAX
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) ZERO , ONE , TWO , FOUR , EIGHTH
	PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FOUR=4.D0,EIGHTH=1.D0/8.D0,&
		& NMAX=100)
	REAL(KIND=DPD) A , B , C , W , Re , Im , temp , temp2 , g1(0:NMAX) &
		& , g2(0:NMAX) , g3(0:NMAX) , g4(0:NMAX) , g5(0:NMAX) , x ,  &
		& x1 , x2 , x3 , x4 , rn , t1(0:80) , t3(0:80) , t4(0:80) ,  &
		& test , Machep , Pi , f1(0:80) , f3(0:80) , f4(0:80) ,      &
		& ff3(0:NMAX) , Eps , ff4(0:NMAX) , coeff1 , coeff2 ,        &
		& c1(0:80) , c3(0:80) , c4(0:80) , sum , term1 , term2 ,     &
		& term3 , term4 , term5 , term6 , coeff , temp1 , et1 , et2 ,&
		& e1 , e2(0:NMAX) , coeff3 , coeff4 , fff1(0:NMAX) ,         &
		& fff2(0:NMAX) , ff1(0:NMAX) , ff2(0:NMAX) , poch1(0:NMAX) , &
		& poch2(0:NMAX) , ttest , error

	INTEGER flag

	x3 = ZERO
	CALL CHEB(c3,55,2)
	t3(0) = ONE
	t3(1) = TWO*(x3-Eps)
	f3(0) = ZERO
	f3(1) = -TWO
	g3(0) = c3(1)*f3(1)

	x4 = ZERO
	CALL CHEB(c4,55,2)
	t4(0) = ONE
	t4(1) = TWO*(x4+Eps)
	f4(0) = ZERO
	f4(1) = TWO
	g4(0) = c4(1)*f4(1)

	DO i = 2 , 55
		t3(i) = FOUR*(x3-Eps)*t3(i-1) - t3(i-2)
		t4(i) = FOUR*(x4+Eps)*t4(i-1) - t4(i-2)
		f3(i) = -FOUR*t3(i-1) + FOUR*x3*f3(i-1) - f3(i-2)
		f4(i) = FOUR*t4(i-1) + FOUR*x4*f4(i-1) - f4(i-2)
		g3(0) = g3(0) + c3(i)*f3(i)
		g4(0) = g4(0) + c4(i)*f4(i)
	ENDDO

	g3(0) = -g3(0)
	DO i = -K , -1
		g4(0) = (g4(0)+ONE/GAMM(DBLE(K+i+2)))/(DBLE(K+i+1)-Eps)
	ENDDO

	test = Eps*DLOG(W)
	temp = DLOG(W)
	IF ( DABS(test)>=EIGHTH ) THEN
		temp = (EXP(test)-ONE)/Eps
	ELSE
		i = 1
50	rn = (Eps**(i)*(DLOG(W))**(i+1))/GAMM(DBLE(i+2))
		IF ( DABS(rn)>=Machep ) THEN
			temp = temp + rn
			i = i + 1
			GOTO 50
		ENDIF
	ENDIF
	g5(0) = temp

	!     write(10,*)g3(0),gamm(-eps)+one/eps
	!     write(10,*)g4(0),(-1)**k*gamm(eps-dble(k))-one/eps/gamm(dble(k+1))
	!     write(10,*)g5(0),w**eps/eps-one/eps

	DO i = 1 , N
		g3(i) = (g3(i-1)+ONE/GAMM(DBLE(i+1)))/(DBLE(i)+Eps)
		g4(i) = (g4(i-1)+ONE/GAMM(DBLE(K+i+1)))/(DBLE(K+i)-Eps)
		g5(i) = W*g5(i-1)
	ENDDO

	DO i = 0 , N
		ff3(i) = Eps*g3(i) - ONE/GAMM(DBLE(i+1))
		ff4(i) = Eps*g4(i) + ONE/GAMM(DBLE(K+i+1))
	ENDDO

	!  calculate the extra terms

	x = A - DBLE(K) - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
100 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 100
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
150	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 150
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 100
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c1,41,1)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - ONE
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 41
			t1(i) = (FOUR*(x+Eps)-TWO)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-TWO)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c1,55,2)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps)
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 55
			t1(i) = FOUR*(x+Eps)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + FOUR*x*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c1,34,3)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - TWO
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 34
			t1(i) = (FOUR*(x+Eps)-FOUR)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-FOUR)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = A - DBLE(K)
		coeff1 = ONE
200	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 200
		ENDIF
		x2 = A - DBLE(K) + Eps
		coeff2 = ONE
250	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 250
		ENDIF
		temp = sum + coeff*temp
		et1 = -temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
300	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 300
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
350	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 350
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et1 = sum + coeff*temp
	ENDIF

	!     write(10,*)et1,(one/gamm(a-dble(k)+eps)-one/gamm(a-dble(k)))
	!    #                                                  /eps

	x = B - DBLE(K) - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
400 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 400
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
450	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 450
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 400
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c1,41,1)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - ONE
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 41
			t1(i) = (FOUR*(x+Eps)-TWO)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-TWO)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c1,55,2)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps)
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 55
			t1(i) = FOUR*(x+Eps)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + FOUR*x*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c1,34,3)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - TWO
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 34
			t1(i) = (FOUR*(x+Eps)-FOUR)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-FOUR)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = B - DBLE(K)
		coeff1 = ONE
500	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 500
		ENDIF
		x2 = B - DBLE(K) + Eps
		coeff2 = ONE
550	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 550
		ENDIF
		temp = sum + coeff*temp
		et2 = -temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
600	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 600
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
650	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 650
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et2 = sum + coeff*temp
	ENDIF

	!     write(10,*)et2,(one/gamm(b-dble(k)+eps)-one/gamm(b-dble(k)))
	!    #                                                  /eps

	fff1(0) = ONE
	DO i = 1 , K
		fff1(0) = (C-B+DBLE(i-1))*fff1(0)
	ENDDO

	ff1(0) = ONE
	e2(0) = ONE
	DO i = 1 , N
		fff1(i) = (C-B+DBLE(K+i-1))*fff1(i-1)
		ff1(i) = (A+DBLE(i-1))*ff1(i-1)
		e2(i) = ff1(i)
	ENDDO

	x1 = A
	coeff1 = ONE
700 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 700
	ENDIF

	x2 = A - DBLE(K)
	coeff2 = ONE
800 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 800
	ENDIF

	DO i = 0 , N
		x3 = B + Eps - DBLE(i+K)
		coeff3 = ONE
850	IF ( x3<ONE ) THEN
			coeff3 = x3*coeff3
			x3 = x3 + ONE
			GOTO 850
		ENDIF

		x4 = B - DBLE(i+K)
		coeff4 = ONE
900	IF ( x4<ONE ) THEN
			coeff4 = x4*coeff4
			x4 = x4 + ONE
			GOTO 900
		ENDIF
		fff1(i) = fff1(i)*coeff1/GAMM(x1)
		ff1(i) = ff1(i)*coeff2/GAMM(x2)
		fff2(i) = coeff3/GAMM(x3)
		ff2(i) = coeff4/GAMM(x4)
		!       write(10,*)'fff1=',fff1(i),gamm(c-b+dble(i+k))/gamm(a)/gamm(c-b)
		!       write(10,*)'ff1=',ff1(i),gamm(a+dble(i))/gamm(a)/gamm(a-dble(k))
		!       write(10,*)'fff2=',fff2(i),one/gamm(b+eps-dble(k+i))
		!       write(10,*)'ff2=',ff2(i),one/gamm(b-dble(k+i))
	ENDDO

	!   calculate  g1

	g1(0) = ZERO
	poch1(0) = ONE
	DO i = 1 , K
		g1(0) = g1(0)*(A+Eps+DBLE(i-K-1)) + poch1(0)
		poch1(0) = poch1(0)*(A+DBLE(i-K-1))
	ENDDO
	DO i = 1 , N
		poch1(i) = (A+DBLE(i-1))*poch1(i-1)
		g1(i) = g1(i-1)*(A+Eps+DBLE(i-1)) + poch1(i-1)
	ENDDO

	x1 = A
	coeff1 = ONE
1000 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1000
	ENDIF
	DO i = 0 , N
		g1(i) = g1(i)*coeff1/GAMM(x1)
		!       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
	ENDDO

	!   calculate  g2

	g2(0) = et2
	DO i = 1 , N
		x1 = B - DBLE(K+i-1)
		coeff1 = ONE
1050 IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 1050
		ENDIF
		g2(i) = g2(i-1)*(B+Eps-DBLE(i+K)) + coeff1/GAMM(x1)
		!       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps
	ENDDO

	!  calculate  e1

	e1 = ZERO
	IF ( DABS(Eps)<.1D0 ) THEN
		i = 1
1100 rn = (-1)**i*Pi**(i+i)*Eps**(i+i-1)/GAMM(DBLE(i+i+1))
		IF ( DABS(rn)>=Machep ) THEN
			e1 = e1 + rn
			i = i + 1
			GOTO 1100
		ENDIF
	ELSE
		e1 = (COS(Pi*Eps)-ONE)/Eps
	ENDIF
	!     write(10,*)'e1=',e1,(cos(pi*eps)-one)/eps

	!  put everything back together again

	ttest = ZERO
	temp = ZERO
	temp1 = ZERO
	DO i = 0 , N
		term1 = g1(i)/GAMM(DBLE(K+i+1))*ff2(i)*ff3(i)*COS(Pi*Eps)      &
			& *W**(DBLE(i)+Eps)
		term2 = -ff1(i)/GAMM(DBLE(K+i+1))*g2(i)*ff3(i)*COS(Pi*Eps)     &
			& *W**(DBLE(i)+Eps)
		term3 = ff1(i)/GAMM(DBLE(K+i+1))*fff2(i)*g3(i)*COS(Pi*Eps)     &
			& *W**(DBLE(i)+Eps)
		term4 = ff1(i)/GAMM(DBLE(i+1))*fff2(i)*g4(i)*COS(Pi*Eps)       &
			& *W**(DBLE(i)+Eps)
		term5 = -ff1(i)/GAMM(DBLE(i+1))*ff4(i)*fff2(i)*COS(Pi*Eps)     &
			& *g5(i)
		term6 = -ff1(i)/GAMM(DBLE(i+1))*ff4(i)*fff2(i)*W**(DBLE(i))*e1
		temp = temp + (term1+term2+term3+term4+term5+term6)*(-1)**(K+i)
		temp1 = temp1 + e2(i)/GAMM(DBLE(i+1))*et1*fff2(i)*ff4(i)       &
			& *W**(DBLE(i))*(-1)**(K+i)
		!       ttest=ttest+(-1)**(k+i)*(cos(pi*eps)*fff1(i)*ff2(i)*ff3(i)
		!    #    /gamm(dble(k+i+1))*w**(dble(i)+eps)+ff1(i)*fff2(i)
		!    #    /gamm(dble(i+1))*ff4(i)*w**(dble(i)))/eps
		!       write(10,*)temp,ttest
		!       ttest=ttest+(-1)**(k+i)*gamm(a+dble(i))/gamm(a)*fff2(i)
		!    #    /gamm(dble(i+1))*ff4(i)*w**(dble(i))*et1
		!       write(10,*)temp1,ttest
		!       ttest=ttest+gamm(dble(i+k+1)-b)/gamm(one-b)*gamm(c-b+dble(i+k))
		!    #    /gamm(c-b)/gamm(a)/gamm(b)*gamm(-eps-dble(i))*(-1)**(i+k)
		!    #    *w**(eps+dble(i))/gamm(dble(i+k+1))*cos(pi*(eps-dble(k)))
		!    #    +gamm(a+dble(i))/gamm(a)*gamm(a-c+dble(i+1))/gamm(a-c+one)
		!    #    /gamm(c-a)/gamm(c-b)*gamm(eps-dble(k+i))*(-w)**i
		!    #    /gamm(dble(i+1))
		!       write(10,*)temp+temp1,ttest
	ENDDO

	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , K - 1
		poch1(i) = (C-B+DBLE(i-1))*poch1(i-1)
		poch2(i) = (DBLE(i)-B)*poch2(i-1)
	ENDDO

	x1 = A
	coeff1 = ONE
1200 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1200
	ENDIF

	x2 = B
	coeff2 = ONE
1300 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 1300
	ENDIF

	temp2 = ZERO
	DO i = 0 , K - 1
		temp2 = temp2 + coeff1*coeff2/GAMM(x1)/GAMM(x2)*poch1(i)       &
			& *poch2(i)*GAMM(DBLE(K-i)-Eps)/GAMM(DBLE(i+1))          &
			& *W**(Eps+DBLE(i-K))*(-1)**i
	ENDDO

	!     term1=zero
	!     do 81 i=0,k-1
	!       term1=term1+gamm(dble(i+1)-b)/gamm(a)*gamm(c-b+dble(i))/gamm(b)
	!    #        /gamm(one-b)/gamm(c-b)*(-1)**i/gamm(dble(i+1))
	!    #        *gamm(dble(k-i)-eps)*w**(eps+dble(i-k))
	! 81  continue
	!     write(10,*)temp2,term1

	Re = GAMM(C)*(ONE-W)**A*(temp+temp1+temp2*COS(Pi*(Eps-DBLE(K))))
	!     write(10,*)re,(ttest+term1*cos(pi*(eps-dble(k))))
	!    #             *gamm(c)*(one-w)**a

	!  calculate the imaginary part

	Im = temp2*SIN(Pi*(Eps-DBLE(K)))

	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , K
		poch1(0) = (C-B+DBLE(i-1))*poch1(0)
		poch2(0) = (DBLE(i)-B)*poch2(0)
	ENDDO
	DO i = 1 , N
		poch1(i) = (C-B+DBLE(K+i-1))*poch1(i-1)
		poch2(i) = (DBLE(K+i)-B)*poch2(i-1)
	ENDDO

	temp = ZERO
	DO i = 0 , N
		temp = temp + poch1(i)/GAMM(DBLE(K+i+1))*ff3(i)                &
			& *W**(Eps+DBLE(i))*poch2(i)
	ENDDO

	x1 = A
	coeff1 = ONE
1400 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1400
	ENDIF

	x2 = B
	coeff2 = ONE
1500 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 1500
	ENDIF

	IF ( DABS(Eps)<.1D0 ) THEN
		temp1 = ONE
		i = 1
1550 temp2 = temp1 + (-1)**i*(Pi*Eps)**(i+i)/GAMM(DBLE(i+i+2))
		error = (temp2-temp1)/temp2
		IF ( DABS(error)>=Machep ) THEN
			i = i + 1
			temp1 = temp2
			GOTO 1550
		ENDIF
	ELSE
		temp2 = SIN(Pi*Eps)/Pi/Eps
	ENDIF
	!     write(10,*)temp2,sin(pi*eps)/pi/eps

	Im = Im + temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)*Pi*temp2
	Im = -Im

	END
	!*==FIX6.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!**********************************************************************
	!
	!  subroutine name    - fix6
	!
	!  computation
	!  performed          - calculates the hypergeometric function for z
	!                       greater than 2 when a-b is near an integer.
	!
	!  usage              - call fix6(a,b,c,n,k,re,im,w,machep,eps,pi)
	!
	!  arguments    a,b,c - parameters of the hypergeometric function.
	!
	!                  n  - the upper limit of the finite series expansion
	!                       of the hypergeometric function.
	!
	!                  k  - equals the nearest integer of a-b.
	!
	!               re,im - computed values for the real and imaginary parts
	!                       of the hypergeometric function.
	!
	!                  w  - transformed independent variable.
	!
	!              machep - equals machine epsilon.
	!
	!                eps  - equals a-b-k.
	!
	!                 pi  - equals 3.1415... to machine accuracy.
	!
	!  precision          - double
	!
	!  language           - fortran
	!
	!***********************************************************************

	SUBROUTINE FIX6(A,B,C,N,K,Re,Im,W,Machep,Eps,Pi)
	IMPLICIT NONE
	!*--FIX63394
	!*** Start of declarations inserted by SPAG
	INTEGER i , K , N , NMAX
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) ZERO , ONE , TWO , FOUR , EIGHTH
	PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FOUR=4.D0,EIGHTH=1.D0/8.D0,&
		& NMAX=100)
	REAL(KIND=DPD) A , B , C , W , Re , Im , temp , temp2 , g1(0:NMAX) &
		& , g2(0:NMAX) , g3(0:NMAX) , g4(0:NMAX) , g5(0:NMAX) , x ,  &
		& x1 , x2 , x3 , x4 , rn , t1(0:80) , t3(0:80) , t4(0:80) ,  &
		& test , Machep , Pi , f1(0:80) , f3(0:80) , f4(0:80) ,      &
		& ff3(0:NMAX) , Eps , ff4(0:NMAX) , coeff1 , coeff2 ,        &
		& c1(0:80) , c3(0:80) , c4(0:80) , sum , term1 , term2 ,     &
		& term3 , term4 , term5 , et1 , et2 , error , term6 , temp1 ,&
		& coeff , coeff3 , coeff4 , fff1(0:NMAX) , ff1(0:NMAX) ,     &
		& fff2(0:NMAX) , ff2(0:NMAX) , poch1(0:NMAX) , poch2(0:NMAX) &
		& , e1

	INTEGER flag

	x3 = ZERO
	CALL CHEB(c3,55,2)
	t3(0) = ONE
	t3(1) = TWO*(x3+Eps)
	f3(0) = ZERO
	f3(1) = TWO
	g3(0) = c3(1)*f3(1)

	x4 = ZERO
	CALL CHEB(c4,55,2)
	t4(0) = ONE
	t4(1) = TWO*(x4-Eps)
	f4(0) = ZERO
	f4(1) = -TWO
	g4(0) = c4(1)*f4(1)

	DO i = 2 , 55
		t3(i) = FOUR*(x3+Eps)*t3(i-1) - t3(i-2)
		t4(i) = FOUR*(x4-Eps)*t4(i-1) - t4(i-2)
		f3(i) = FOUR*t3(i-1) + FOUR*x3*f3(i-1) - f3(i-2)
		f4(i) = -FOUR*t4(i-1) + FOUR*x4*f4(i-1) - f4(i-2)
		g3(0) = g3(0) + c3(i)*f3(i)
		g4(0) = g4(0) + c4(i)*f4(i)
	ENDDO

	g4(0) = -g4(0)
	DO i = -K , -1
		g4(0) = (g4(0)+ONE/GAMM(DBLE(K+i+2)))/(DBLE(K+i+1)+Eps)
	ENDDO

	test = -Eps*DLOG(W)
	temp = -DLOG(W)
	IF ( DABS(test)>=EIGHTH ) THEN
		temp = (EXP(test)-ONE)/Eps
	ELSE
		i = 1
50	rn = (Eps**(i)*(-DLOG(W))**(i+1))/GAMM(DBLE(i+2))
		IF ( DABS(rn)>=Machep ) THEN
			temp = temp + rn
			i = i + 1
			GOTO 50
		ENDIF
	ENDIF
	g5(0) = temp*W**A

	!     write(10,*)g3(0),gamm(eps)-one/eps
	!     write(10,*)g4(0),(-1)**k*gamm(-eps-k)+one/eps/gamm(dble(k+1))
	!     write(10,*)g5(0),w**(a-eps)/eps-w**a/eps

	DO i = 1 , N
		g3(i) = (g3(i-1)+ONE/GAMM(DBLE(i+1)))/(DBLE(i)-Eps)
		g4(i) = (g4(i-1)+ONE/GAMM(DBLE(K+i+1)))/(DBLE(K+i)+Eps)
		g5(i) = W*g5(i-1)
	ENDDO

	DO i = 0 , N
		ff3(i) = Eps*g3(i) + ONE/GAMM(DBLE(i+1))
		ff4(i) = Eps*g4(i) - ONE/GAMM(DBLE(K+i+1))
	ENDDO

	!  calculate the extra terms

	x = B - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
100 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 100
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
150	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 150
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 100
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c1,41,1)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - ONE
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 41
			t1(i) = (FOUR*(x+Eps)-TWO)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-TWO)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c1,55,2)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps)
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 55
			t1(i) = FOUR*(x+Eps)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + FOUR*x*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c1,34,3)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - TWO
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 34
			t1(i) = (FOUR*(x+Eps)-FOUR)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-FOUR)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = B
		coeff1 = ONE
200	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 200
		ENDIF
		x2 = B + Eps
		coeff2 = ONE
250	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 250
		ENDIF
		temp = sum + coeff*temp
		et1 = -temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
300	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 300
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
350	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 350
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et1 = sum + coeff*temp
	ENDIF
	et1 = -et1
	!     write(10,*)et1,(one/gamm(a-dble(k)-eps)-one/gamm(a-dble(k)))
	!    #                                                  /eps

	x = C - A + K - ONE
	sum = ZERO
	coeff = ONE
	flag = 0
400 IF ( x>ONE ) THEN
		sum = sum + coeff*GAMM(x+Eps)
		coeff = coeff*x
		x = x - ONE
		GOTO 400
	ELSEIF ( x<ZERO ) THEN
		x1 = x + Eps + TWO
		coeff1 = ONE
450	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 450
		ENDIF
		sum = sum + coeff*coeff1/GAMM(x1)
		coeff = coeff*(x+ONE)
		x = x + ONE
		flag = 1
		GOTO 400
	ENDIF

	IF ( (x>=.25D0) .AND. (x<=.75D0) ) THEN
		CALL CHEB(c1,41,1)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - ONE
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 41
			t1(i) = (FOUR*(x+Eps)-TWO)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-TWO)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>=0.D0) .AND. (x<.25D0) ) THEN
		CALL CHEB(c1,55,2)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps)
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 55
			t1(i) = FOUR*(x+Eps)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + FOUR*x*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ELSEIF ( (x>.75D0) .AND. (x<=1.D0) ) THEN
		CALL CHEB(c1,34,3)
		t1(0) = ONE
		t1(1) = TWO*(x+Eps) - TWO
		f1(0) = ZERO
		f1(1) = TWO
		temp = c1(1)*f1(1)
		DO i = 2 , 34
			t1(i) = (FOUR*(x+Eps)-FOUR)*t1(i-1) - t1(i-2)
			f1(i) = FOUR*t1(i-1) + (FOUR*x-FOUR)*f1(i-1) - f1(i-2)
			temp = temp + c1(i)*f1(i)
		ENDDO
	ENDIF

	IF ( flag==0 ) THEN
		x1 = C - A + DBLE(K)
		coeff1 = ONE
500	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 500
		ENDIF
		x2 = C - A + DBLE(K) + Eps
		coeff2 = ONE
550	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 550
		ENDIF
		temp = sum + coeff*temp
		et2 = -temp*coeff1*coeff2/GAMM(x1)/GAMM(x2)
	ELSEIF ( flag==ONE ) THEN
		x1 = x + ONE
		coeff1 = ONE
600	IF ( x1<ONE ) THEN
			coeff1 = x1*coeff1
			x1 = x1 + ONE
			GOTO 600
		ENDIF
		x2 = x + ONE + Eps
		coeff2 = ONE
650	IF ( x2<ONE ) THEN
			coeff2 = x2*coeff2
			x2 = x2 + ONE
			GOTO 650
		ENDIF
		coeff = -coeff*coeff1*coeff2/GAMM(x1)/GAMM(x2)
		et2 = sum + coeff*temp
	ENDIF
	et2 = -et2
	!     write(10,*)et2,(one/gamm(c-b-eps)-one/gamm(c-b))/eps
	!
	fff1(0) = ONE
	fff2(0) = ONE
	ff2(0) = ONE
	DO i = 1 , K
		fff1(0) = (B+DBLE(i-1))*fff1(0)
		fff2(0) = (B-C+Eps+DBLE(i))*fff2(0)
		ff2(0) = (B-C+DBLE(i))*ff2(0)
	ENDDO

	ff1(0) = ONE
	DO i = 1 , N
		fff1(i) = (B+DBLE(K+i-1))*fff1(i-1)
		fff2(i) = (B-C+Eps+DBLE(K+i))*fff2(i-1)
		ff1(i) = (A+DBLE(i-1))*ff1(i-1)
		ff2(i) = (B-C+DBLE(K+i))*ff2(i-1)
	ENDDO

	x1 = A
	coeff1 = ONE
700 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 700
	ENDIF

	x2 = A - DBLE(K)
	coeff2 = ONE
800 IF ( x2<ONE ) THEN
		coeff2 = x2*coeff2
		x2 = x2 + ONE
		GOTO 800
	ENDIF

	x3 = C - B - Eps
	coeff3 = ONE
900 IF ( x3<ONE ) THEN
		coeff3 = x3*coeff3
		x3 = x3 + ONE
		GOTO 900
	ENDIF

	x4 = C - B
	coeff4 = ONE
1000 IF ( x4<ONE ) THEN
		coeff4 = x4*coeff4
		x4 = x4 + ONE
		GOTO 1000
	ENDIF

	DO i = 0 , N
		fff1(i) = fff1(i)*coeff1/GAMM(x1)
		ff1(i) = ff1(i)*coeff2/GAMM(x2)
		fff2(i) = fff2(i)*coeff3/GAMM(x3)
		ff2(i) = ff2(i)*coeff4/GAMM(x4)
		!       write(10,*)'fff1=',fff1(i),gamm(b+dble(i+k))/gamm(a)/gamm(b)
		!       write(10,*)'ff1=',ff1(i),gamm(a+dble(i))/gamm(a)/gamm(a-dble(k))
		!       write(10,*)'fff2=',fff2(i),gamm(a-c+dble(i+1))/gamm(c-b-eps)
		!    #                               /gamm(one-c+b+eps)
		!       write(10,*)'ff2=',ff2(i),gamm(b-c+dble(i+k+1))/gamm(c-b)
		!    #                               /gamm(one-c+b)
	ENDDO

	!   calculate  g1

	g1(0) = ZERO
	poch1(0) = ONE
	DO i = 1 , K
		g1(0) = g1(0)*(A-Eps+DBLE(i-K-1)) - poch1(0)
		poch1(0) = poch1(0)*(A+DBLE(i-K-1))
	ENDDO
	DO i = 1 , N
		poch1(i) = (A+DBLE(i-1))*poch1(i-1)
		g1(i) = g1(i-1)*(A-Eps+DBLE(i-1)) - poch1(i-1)
	ENDDO

	x1 = A
	coeff1 = ONE
1100 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1100
	ENDIF
	DO i = 0 , N
		g1(i) = g1(i)*coeff1/GAMM(x1)
		!       write(10,*)'g1=',g1(i),(fff1(i)-ff1(i))/eps
	ENDDO

	!   calculate  g2

	g2(0) = ZERO
	poch2(0) = ONE
	DO i = 1 , K
		g2(0) = g2(0)*(B-C+Eps+DBLE(i)) + poch2(0)
		poch2(0) = poch2(0)*(B-C+DBLE(i))
	ENDDO
	DO i = 1 , N
		poch2(i) = (B-C+DBLE(i+K))*poch2(i-1)
		g2(i) = g2(i-1)*(B-C+Eps+DBLE(i+K)) + poch2(i-1)
	ENDDO

	x1 = C - B
	coeff1 = ONE
1200 IF ( x1<ONE ) THEN
		coeff1 = x1*coeff1
		x1 = x1 + ONE
		GOTO 1200
	ENDIF

	poch2(0) = ONE
	DO i = 1 , K
		poch2(0) = (B-C+Eps+DBLE(i))*poch2(0)
	ENDDO
	DO i = 1 , N
		poch2(i) = (B-C+Eps+DBLE(i+K))*poch2(i-1)
	ENDDO

	DO i = 0 , N
		g2(i) = et2*poch2(i) + g2(i)*coeff1/GAMM(x1)
		!       write(10,*)'g2=',g2(i),(fff2(i)-ff2(i))/eps
	ENDDO

	!  calculate  e1

	e1 = ZERO
	IF ( DABS(Eps)<.1D0 ) THEN
		i = 1
1250 rn = (-1)**i*Pi**(i+i)*Eps**(i+i-1)/GAMM(DBLE(i+i+1))
		IF ( DABS(rn)>=Machep ) THEN
			e1 = e1 + rn
			i = i + 1
			GOTO 1250
		ENDIF
	ELSE
		e1 = (COS(Pi*Eps)-ONE)/Eps
	ENDIF
	!     write(10,*)'e1=',e1,(cos(pi*eps)-one)/eps

	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , N
		poch1(i) = (A+DBLE(i-1))*poch1(i-1)
		poch2(i) = (A-C+DBLE(i))*poch2(i-1)
	ENDDO

	!  put everything back together again

	temp = ZERO
	temp1 = ZERO
	temp2 = ZERO
	DO i = 0 , N
		term1 = -g1(i)/GAMM(DBLE(i+1))*fff2(i)*ff4(i)*W**(A+DBLE(i))   &
			& *COS(Pi*Eps)
		term2 = fff1(i)/GAMM(DBLE(i+1))*g2(i)*ff4(i)*W**(A+DBLE(i))    &
			& *COS(Pi*Eps)
		term3 = -fff1(i)*ff2(i)*g3(i)*ff4(i)*W**(A+DBLE(i))*COS(Pi*Eps)
		term4 = fff1(i)*ff2(i)*ff3(i)*g4(i)*W**(A+DBLE(i))*COS(Pi*Eps)
		term5 = fff1(i)/GAMM(DBLE(K+i+1))*ff2(i)*ff3(i)*g5(i)          &
			& *COS(Pi*Eps)
		term6 = -fff1(i)/GAMM(DBLE(K+i+1))*ff2(i)*ff3(i)               &
			& *W**(A-Eps+DBLE(i))*e1
		temp = temp + term1 + term2 + term3 + term4 + term5 + term6
		temp1 = temp1 + poch1(i)/GAMM(DBLE(i+1))*poch2(i)*ff4(i)       &
			& *W**(A+DBLE(i))
		temp2 = temp2 + poch1(i)/GAMM(DBLE(i+1))*fff2(i)*ff4(i)        &
			& *W**(A+DBLE(i))
	ENDDO

	x1 = B
	coeff1 = ONE
1300 IF ( x1<ONE ) THEN
		coeff1 = coeff1*x1
		x1 = x1 + ONE
		GOTO 1300
	ENDIF

	x2 = C - A
	coeff2 = ONE
1400 IF ( x2<ONE ) THEN
		coeff2 = coeff2*x2
		x2 = x2 + ONE
		GOTO 1400
	ENDIF

	term1 = temp*GAMM(C)*COS(Pi*B)*(-1)**K
	term2 = -temp1*GAMM(C)*SIN(Pi*B)*coeff1/GAMM(x1)*coeff2/GAMM(x2)
	term3 = temp2*GAMM(C)*COS(Pi*B)*COS(Pi*Eps)*(-1)**K*et1
	term4 = temp*GAMM(C)*SIN(Pi*B)*(-1)**K
	term5 = temp1*GAMM(C)*COS(Pi*B)*coeff1/GAMM(x1)*coeff2/GAMM(x2)
	term6 = term6*GAMM(C)*SIN(Pi*B)*COS(Pi*Eps)*(-1)**K*et1

	IF ( DABS(Eps)<.1D0 ) THEN
		temp1 = ONE
		i = 1
1450 temp2 = temp1 + (-1)**i*(Pi*Eps)**(i+i)/GAMM(DBLE(i+i+2))
		error = (temp2-temp1)/temp2
		IF ( DABS(error)>=Machep ) THEN
			i = i + 1
			temp1 = temp2
			GOTO 1450
		ENDIF
	ELSE
		temp2 = SIN(Pi*Eps)/Pi/Eps
	ENDIF
	!     write(10,*)temp2,sin(pi*eps)/(pi*eps)

	term2 = term2*Pi*temp2
	term5 = term5*Pi*temp2

	Re = term1 + term2 + term3
	Im = term4 + term5 + term6

	!  calculate the finite series contribution

	poch1(0) = ONE
	poch2(0) = ONE
	DO i = 1 , N
		poch1(i) = (B+DBLE(i-1))*poch1(i-1)
		poch2(i) = (B-C+DBLE(i))*poch2(i-1)
	ENDDO

	temp = ZERO
	DO i = 0 , K - 1
		temp = temp + poch1(i)*poch2(i)/GAMM(DBLE(i+1))                &
			& *GAMM(Eps+DBLE(K-i))*(-1)**i*W**(B+DBLE(i))
	ENDDO

	x1 = A
	coeff1 = ONE
1500 IF ( x1<ONE ) THEN
		coeff1 = coeff1*x1
		x1 = x1 + ONE
		GOTO 1500
	ENDIF

	x2 = C - B
	coeff2 = ONE
1600 IF ( x2<ONE ) THEN
		coeff2 = coeff2*x2
		x2 = x2 + ONE
		GOTO 1600
	ENDIF

	temp = temp*GAMM(C)*coeff1/GAMM(x1)*coeff2/GAMM(x2)

	Re = Re + temp*COS(Pi*B)
	Im = Im + temp*SIN(Pi*B)

	END
	!*==GETEPS.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!**********************************************************************
	!
	!   subroutine name     - geteps
	!
	!   computation
	!   performed           - compute the smallest number machep such that
	!                           machep+1 is not equal to 1 in the finite
	!                           precision arithmetic used by the computer.
	!
	!   usage               - call geteps(machep,neps)
	!
	!   argument     machep - double precision (output).  the smallest
	!                           number such that machep+1 is not equal to 1
	!                           in the finite precision arithmetic used by
	!                           the computer.
	!                  neps - integer (output).  machine epsilon is machep =
	!                           (1/2)**neps
	!
	!   precision           - double
	!
	!   language            - fortran 77
	!
	!***********************************************************************
	!
	SUBROUTINE GETEPS(Machep,Neps)
	IMPLICIT NONE
	!*--GETEPS3954
	!
	REAL(KIND=DPD) Machep , ONE , TWO , temp
	INTEGER Neps
	PARAMETER (ONE=1.0D0,TWO=2.0D0)
	Machep = ONE
	Neps = 0
100 Machep = Machep/TWO
	Neps = Neps + 1
	temp = Machep + ONE
	IF ( temp/=ONE ) GOTO 100
	!
	Machep = TWO*Machep
	Neps = Neps - 1
	!
	END
	!*==BINOMC.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019
	!
	!***********************************************************************

	SUBROUTINE BINOMC
	IMPLICIT NONE
	!*--BINOMC3976
	!*** Start of declarations inserted by SPAG
	INTEGER i , ii , ij , imax , j , jmax , maxnll
	!*** End of declarations inserted by SPAG

	!       a
	!     (   ) = binom(a*(a+1)/2+b+1)
	!       b

	DOUBLE PRECISION BINom , one
	COMMON /BCOEFF/ BINom(5151)

	maxnll = 100

	IF ( maxnll>=0 ) THEN
		IF ( maxnll<=100 ) THEN
			one = 1.0D0
			BINom(1) = one
			IF ( maxnll/=0 ) THEN
				BINom(2) = one
				BINom(3) = one
				IF ( maxnll/=1 ) THEN
					ij = 4
					imax = maxnll + 1
					DO i = 3 , imax
						ii = ((i-1)*(i-2))/2
						BINom(ij) = one
						ij = ij + 1
						jmax = i - 1
						DO j = 2 , jmax
							BINom(ij) = BINom(ii+j-1) + BINom(ii+j)
							ij = ij + 1
						ENDDO
						BINom(ij) = one
						ij = ij + 1
					ENDDO
				ENDIF
			ENDIF
		ENDIF
	ENDIF
	!
	!
	END
	!*==CHEB.spg  processed by SPAG 6.72Dc at 13:31 on  9 Jul 2019

	!***********************************************************************
	!
	!  subroutine name    - cheb
	!
	!  computation
	!  performed          - tabulates the tchebychev coefficients which
	!                       were computed by the program 'tcheb2'.  the
	!                       three sets of coefficients correspond to
	!                       the three gamma function expansions shown in
	!                       equations (4.35),(4.36), and (4.37). see
	!                       'tcheb2' for additional documentation.
	!
	!  usage              - call cheb(c,n,flag)
	!
	!  arguments       c  - the array (output) which contains the
	!                       tchebychev coefficients.
	!
	!                  n  - the dimension (input) of the array 'c'.
	!
	!                flag - the parameter (input) which tells the sub-
	!                       routine which tchebychev coefficients to
	!                       return to the caller.
	!
	!  precision          - double (although the coefficients are
	!                               accurate to quadruple)
	!
	!  language           - fortran 77
	!
	!***********************************************************************

	SUBROUTINE CHEB(C,N,Flag)
	IMPLICIT NONE
	!*--CHEB4053
	!*** Start of declarations inserted by SPAG
	INTEGER N
	!*** End of declarations inserted by SPAG

	REAL(KIND=DPD) C(0:N)
	INTEGER Flag

	IF ( Flag/=1 ) THEN
		IF ( Flag==2 ) THEN

			!  tchebychev expansion coefficients for the range,  -.5<x<.5

			C(0) = 0.11528686913857579339872890819003657D+01
			C(1) = -0.39836641427188668813550502856567435D+00
			C(2) = 0.16381491849746834445969671065563396D+00
			C(3) = -0.41349972584595838242416447164595642D-01
			C(4) = 0.11739888104509743948748485834561229D-01
			C(5) = -0.31509159742825717845846783104528302D-02
			C(6) = 0.85084809366682540330028115184077086D-03
			C(7) = -0.22845443192182297253614554810213881D-03
			C(8) = 0.61296656896858907270916323759970391D-04
			C(9) = -0.16433766723011959082591541534833589D-04
			C(10) = 0.44046701847148520660258125028242579D-05
			C(11) = -0.11803851479587223345492859134791582D-05
			C(12) = 0.31630339312403588488305625683201151D-06
			C(13) = -0.84755796666686117564957022251013564D-07
			C(14) = 0.22710572677209079780536954678987573D-07
			C(15) = -0.60853209609268373214751556259951644D-08
			C(16) = 0.16305620921375867864482570008163625D-08
			C(17) = -0.43690846345047718022878883179027790D-09
			C(18) = 0.11706935476739890379554689241357534D-09
			C(19) = -0.31368649843198552351255033209421610D-10
			C(20) = 0.84052057618382692960217222664957228D-11
			C(21) = -0.22521682699590609081199019088965996D-11
			C(22) = 0.60346669123807723976181127096882828D-12
			C(23) = -0.16169841538137032176079290114309245D-12
			C(24) = 0.43326960175123609635570088625382667D-13
			C(25) = -0.11609424034675431553315176322024985D-13
			C(26) = 0.31107358004300087572452155428660087D-14
			C(27) = -0.83351914632193111475558815401948979D-15
			C(28) = 0.22334078222557889355389486422061460D-15
			C(29) = -0.59843982246058550382747881611851515D-16
			C(30) = 0.16035146716190080240936859943115090D-16
			C(31) = -0.42966046133076898235808019603294715D-17
			C(32) = 0.11512717363557431988678458870224873D-17
			C(33) = -0.30848233202835882015258583966299712D-18
			C(34) = 0.82657591746540727258216017499064442D-19
			C(35) = -0.22148034956862123422799663231945171D-19
			C(36) = 0.59345480806145642339133686333296721D-20
			C(37) = -0.15901573656881585725893714030807897D-20
			C(38) = 0.42608138203898096080539369435375448D-21
			C(39) = -0.11416816226321087557458906349840213D-21
			C(40) = 0.30591266842950015571055286508657438D-22
			C(41) = -0.81969053674548061989664444282339330D-23
			C(42) = 0.21963543471485197662543467891802004D-23
			C(43) = -0.58851140572211577956963471197095354D-24
			C(44) = 0.15769121438531798083082131134888596D-24
			C(45) = -0.42253211944581570323425035302537635D-25
			C(46) = 0.11321706791574145306428072576766804D-25
			C(47) = -0.30335842761477973373797446515125892D-26
			C(48) = 0.81281383350578045680446098123885346D-27
			C(49) = -0.21782407988772728568103833180457024D-27
			C(50) = 0.58395544064782062129754390403734767D-28
			C(51) = -0.15729062977489325257494410942884130D-28
			C(52) = 0.42390612257722955199550993363196147D-29
			C(53) = -0.11242203351086692027388616387423238D-29
			C(54) = 0.27892280419588143241883200553486195D-30
			C(55) = -0.75766427928255356179910217971637866D-31
			RETURN
		ELSEIF ( Flag==3 ) THEN

			!  tchebychev expansion coefficients for the range,  .5<x<1.5

			C(0) = 0.10532770878177862619534128247576828D+01
			C(1) = 0.21902166104535936497306369004840667D+00
			C(2) = 0.53885821783347712865216341722976574D-01
			C(3) = 0.25387290658986838596948519579519148D-02
			C(4) = 0.61466596479014144199820446583715941D-03
			C(5) = -0.32319247384294465724865638122474435D-05
			C(6) = 0.60054921157267140200751871810266970D-05
			C(7) = -0.41824428090189489334617924547407754D-06
			C(8) = 0.74607235650174366232051332482639985D-07
			C(9) = -0.84349526185192483560074198183789434D-08
			C(10) = 0.11322169721817117406057072389666464D-08
			C(11) = -0.14175349900034682206860980369914924D-09
			C(12) = 0.18156967683771854495445069753509525D-10
			C(13) = -0.23052163748763990586386231147733255D-11
			C(14) = 0.29327030584105892891631030300077869D-12
			C(15) = -0.37268590170679729030689484336505900D-13
			C(16) = 0.47360432581610222494078892575939043D-14
			C(17) = -0.60172423075766780010690060490450222D-15
			C(18) = 0.76443979970650480527157880770622904D-16
			C(19) = -0.97108892590783757664936380167684001D-17
			C(20) = 0.12335488659810502174628042595177563D-17
			C(21) = -0.15668997427423797214874298423999374D-18
			C(22) = 0.19902969432180950952170748993213290D-19
			C(23) = -0.25280701093316992983208535829903356D-20
			C(24) = 0.32111217127088658654008440525466587D-21
			C(25) = -0.40787027055654288157193053732139852D-22
			C(26) = 0.51806681115442807351458062924762066D-23
			C(27) = -0.65803415226414646040514708695329147D-24
			C(28) = 0.83581632724068042390791744946381128D-25
			C(29) = -0.10616267321620223331012310816058461D-25
			C(30) = 0.13484159784261929973156667845986312D-26
			C(31) = -0.17130640476670792317750095910458264D-27
			C(32) = 0.21720215147689502411187819143753676D-28
			C(33) = -0.27633054946463729557612727034555572D-29
			C(34) = 0.26664265210535867308016959008022142D-30
			GOTO 99999
		ENDIF
	ENDIF

	!  tchebychev expansion coefficients for the range, 0<x<1

	C(0) = 0.94178559779549466571096003120435196D+00
	C(1) = 0.44153813248410067571913157711414607D-02
	C(2) = 0.56850436815993633786326645888162378D-01
	C(3) = -0.42198353964185605010125001866024699D-02
	C(4) = 0.13268081812124602205840067963889683D-02
	C(5) = -0.18930245297988804325239470239464680D-03
	C(6) = 0.36069253274412452565780822094442805D-04
	C(7) = -0.60567619044608642184855483216922771D-05
	C(8) = 0.10558295463022833447318234541645507D-05
	C(9) = -0.18119673655423840482918555144273961D-06
	C(10) = 0.31177249647153222777902517006137963D-07
	C(11) = -0.53542196390196871408740949118221475D-08
	C(12) = 0.91932755198595889468877475468573503D-09
	C(13) = -0.15779412802883397617671187106425584D-09
	C(14) = 0.27079806229349545432695717700017206D-10
	C(15) = -0.46468186538257301439531283506784063D-11
	C(16) = 0.79733501920074196555512936759234830D-12
	C(17) = -0.13680782098309160264738694164685656D-12
	C(18) = 0.23473194865638006534799539031857605D-13
	C(19) = -0.40274326149490669507857892267787757D-14
	C(20) = 0.69100517473721009958174457696435176D-15
	C(21) = -0.11855845002219929396593062972684083D-15
	C(22) = 0.20341485424963760969383490105975402D-16
	C(23) = -0.34900543417173691101844936408331408D-17
	C(24) = 0.59879938564842634972645168624438135D-18
	C(25) = -0.10273780578716378747008169519685451D-18
	C(26) = 0.17627028160574041125936108594612916D-19
	C(27) = -0.30243206536626379817809691872233988D-20
	C(28) = 0.51889146600668142375785699199940389D-21
	C(29) = -0.89027708392150216484577040964212789D-22
	C(30) = 0.15274740724470977041487116294681806D-22
	C(31) = -0.26207312865170684216151526387496724D-23
	C(32) = 0.44964644619824783627762340991300087D-24
	C(33) = -0.77147147879836211531329396406348717D-25
	C(34) = 0.13236365808260955301316348853544449D-25
	C(35) = -0.22709797413377406198008958539204735D-26
	C(36) = 0.38966913277073699893252807432563276D-27
	C(37) = -0.66795989154793901466615113245736539D-28
	C(38) = 0.11456694360946249087722449327564468D-28
	C(39) = -0.20956088513945987438866120550893160D-29
	C(40) = 0.34345153487326051089311279207743562D-30
	C(41) = -0.74448389617685196161619686887550341D-31
	RETURN
99999 END

	END MODULE com_mod_external_functions