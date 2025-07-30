!--------------------------------------------------------
!  SOLVE THE TRIDIAGONAL MATRIX 
!    [b1  c1           ]
!    [a2  b2  c2       ]
!    [...              ]
!    [...              ]
!    [           aN  bN]
!--------------------------------------------------------

   SUBROUTINE tdma(a,b,c,r,u,n1,n2)

    USE global_parameters    
    IMPLICIT NONE

    INTEGER                           , INTENT (IN)  :: n1,n2
    REAL(KIND=CGREAL), DIMENSION(0:n2), INTENT (IN)  :: a,b,c,r
    REAL(KIND=CGREAL), DIMENSION(0:n2), INTENT (OUT) :: u

    REAL(KIND=CGREAL)                             :: binv
    REAL(KIND=CGREAL), DIMENSION(0:n2)            :: gam
    INTEGER                          :: j

    binv  = 1.0_CGREAL/b(n1)
    u(n1) = r(n1)*binv
    DO j=n1+1,n2
      gam(j)  = c(j-1)*binv
      binv    = 1.0_CGREAL/(b(j)-a(j)*gam(j))
      u(j)    = (r(j)-a(j)*u(j-1))*binv
    ENDDO
    DO j=n2-1,n1,-1
      u(j) = u(j) - gam(j+1)*u(j+1)
    ENDDO

   END SUBROUTINE tdma
!--------------------------------------------------------
   SUBROUTINE tdma2(a,b,c,r,u,ix1,ix2,n1,n2)
 
    USE global_parameters
    USE flow_parameters
 
    IMPLICIT NONE
 
    INTEGER          , INTENT (IN)  :: ix1,ix2,n1,n2
    REAL(KIND=CGREAL), INTENT (IN)  :: a(0:nx+1,0:n2)
    REAL(KIND=CGREAL), INTENT (IN)  :: b(0:nx+1,0:n2)
    REAL(KIND=CGREAL), INTENT (IN)  :: c(0:nx+1,0:n2)
    REAL(KIND=CGREAL), INTENT (IN)  :: r(0:nx+1,0:n2)
    REAL(KIND=CGREAL), INTENT (OUT) :: u(0:nx+1,0:n2)
 
    REAL(KIND=CGREAL) :: binv(ix1:ix2)
    REAL(KIND=CGREAL) :: gam(ix1:ix2,0:n2)
    INTEGER           :: ix,k
 
    do ix = ix1, ix2
      binv(ix)  = 1.0_CGREAL/b(ix,n1)
      u(ix,n1) = r(ix,n1)*binv(ix)
    end do 
    DO k=n1+1,n2 
      do ix = ix1, ix2 
        gam(ix,k)  = c(ix,k-1)*binv(ix)
        binv(ix)    = 1.0_CGREAL/(b(ix,k)-a(ix,k)*gam(ix,k))
        u(ix,k)    = (r(ix,k)-a(ix,k)*u(ix,k-1))*binv(ix)
      end do
    ENDDO 
    DO k=n2-1,n1,-1 
      do ix = ix1, ix2
        u(ix,k) = u(ix,k) - gam(ix,k+1)*u(ix,k+1)
      end do
    ENDDO 
       
   END SUBROUTINE tdma2
 

!******************************************************************
!...Subroutine to solve Ax=b where A is a cyclic tridiagonal matrix
!    [b c   Be]
!    [a b c   ]
!    [  a b c ]
!    [Al  a b ]
!******************************************************************
    SUBROUTINE cyclic_tdma(a,b,c,r,u,alp,bet,n1,n2)

      USE global_parameters    
      IMPLICIT NONE

      INTEGER                           , INTENT (IN)  :: n1,n2
      REAL(KIND=CGREAL), DIMENSION(0:n2), INTENT (IN)  :: a,b,c,r
      REAL(KIND=CGREAL)                 , INTENT (IN)  :: alp,bet
      REAL(KIND=CGREAL), DIMENSION(0:n2), INTENT (OUT) :: u

      REAL(KIND=CGREAL)                                :: gam,fact
      REAL(KIND=CGREAL), DIMENSION(0:n2)               :: x,y,z,bb
      INTEGER                                          :: j

      gam      = -b(n1)
      bb(n1)   =  b(n1)-gam
      bb(n2)   =  b(n2)-alp*bet/gam
 
      DO J=n1+1,n2-1
        bb(j) = b(j)
      ENDDO
 
      CALL tdma(a,bb,c,r,x,n1,n2) 

      y(n1) = gam
      y(n2) = alp
 
      DO j=n1+1,n2-1
        y(j) = 0.0_CGREAL
      ENDDO
 
      CALL tdma(a,bb,c,y,z,n1,n2) 

      fact =         ( x(n1) + bet*x(n2)/gam ) &
               /( 1.0_CGREAL + z(n1) + bet*z(n2)/gam )
 
      DO j=n1,n2
         x(j) = x(j) - fact*z(j)
         u(j) = x(j)
      ENDDO
 
   END SUBROUTINE cyclic_tdma
     
