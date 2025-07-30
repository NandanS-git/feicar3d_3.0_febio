
!-------------------------------------------------------------------
! Crank--Nicolson scheme for both convective and diffusion terms
!
!          n    dt {    n  }             {    n         n       n  }
!  RHS  = u   - ---{  NL   } + dt(1/2Re) {am u    + ac u  + ap u   }
!                2 {    i  }             {  i i-1     i i     i i+1} 
!-------------------------------------------------------------------

   SUBROUTINE rhs_advec_diff() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE MPI_module

    IMPLICIT NONE

    INTEGER :: i,j,k,n
!***********************************************

!------------------------------
!   Non-linear convection term
!------------------------------

    CALL rhs_advec(uOld,nlu)   
    CALL rhs_advec(vOld,nlv)    
    CALL rhs_advec(wOld,nlw)    

!------------------------------    
!   Diffusion terms 
!------------------------------

    CALL rhs_diff(uOld,nlu)    
    CALL rhs_diff(vOld,nlv)    
    CALL rhs_diff(wOld,nlw)    

    !PRINT*, 'Min-Max of NLU =',minval(nlu(1:nx-1,1:ny-1,1:nz-1)),maxval(nlu(1:nx-1,1:ny-1,1:nz-1))
    !PRINT*, 'Min-Max of NLV =',minval(nlv(1:nx-1,1:ny-1,1:nz-1)),maxval(nlv(1:nx-1,1:ny-1,1:nz-1))
    !PRINT*, 'Min-Max of NLW =',minval(nlw(1:nx-1,1:ny-1,1:nz-1)),maxval(nlw(1:nx-1,1:ny-1,1:nz-1))

    !store the r.h.s. terms at the hybrid nodes
    DO k = zc_start,zc_end   !1,nz-1
    DO j = yc_start,yc_end
    DO i = 1,nx-1
       if(abs(ghostCellMark(i,j,k)) == 1 ) then
          n = ghostCellIndex(i,j,k)
          ghostNScoeff(8 ,n) = nlu(i,j,k)
          ghostNScoeff(9 ,n) = nlv(i,j,k)
          ghostNScoeff(10,n) = nlw(i,j,k)
       endif
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k 

   END SUBROUTINE  rhs_advec_diff 

!------------------------------------------------------------------------------

    SUBROUTINE rhs_advec(vel,nlvel) 

      USE global_parameters
      USE flow_parameters
      USE flow_arrays
      USE grid_arrays
      USE boundary_arrays
      USE MPI_module

      IMPLICIT NONE

!... parameters
   
      REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(IN)  :: vel
      REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(OUT) :: nlvel

!... loop  variables

      INTEGER           :: i,j,k
      INTEGER           :: imm,jmm,kmm         ! Added by Rupesh
      INTEGER           :: ipp,jpp,kpp         ! used in 2nd Upwinding


      REAL(KIND=CGREAL) :: edxWeightm1, edxWeightm2, edxWeightp1, edxWeightp2
      REAL(KIND=CGREAL) :: wdxWeightm1, wdxWeightm2, wdxWeightp1, wdxWeightp2

      REAL(KIND=CGREAL) :: ndxWeightm1, ndxWeightm2, ndxWeightp1, ndxWeightp2
      REAL(KIND=CGREAL) :: sdxWeightm1, sdxWeightm2, sdxWeightp1, sdxWeightp2

      REAL(KIND=CGREAL) :: fdxWeightm1, fdxWeightm2, fdxWeightp1, fdxWeightp2
      REAL(KIND=CGREAL) :: bdxWeightm1, bdxWeightm2, bdxWeightp1, bdxWeightp2

      REAL(KIND=CGREAL) :: wsign, wsignP, Usign, UsignP, Vsign, VsignP

      REAL(KIND=CGREAL) :: vele,velw,veln,vels,velb,velf
      REAL(KIND=CGREAL) :: vele_Up,velw_Up,veln_Up,vels_Up,velb_Up,velf_Up ! Added by Rupesh 
      REAL(KIND=CGREAL) :: tempvel

      !Ye,hybrid
      REAL(KIND=CGREAL) :: vele_Up1,velw_Up1,veln_Up1,vels_Up1,velb_Up1,velf_Up1
      REAL(KIND=CGREAL) :: vele_Up2,velw_Up2,veln_Up2,vels_Up2,velb_Up2,velf_Up2

!******************************************************************************

!---------------------------------------------------------------------------------------------------
! Convective terms:
! nlvel = d(vel*U)/dx + d(vel*V)/dy 
!---------------------------------------------------------------------------------------------------
!   The following decription was added by Rupesh and pertains to 2nd Upwind scheme for 
!   convected face velocities.  
!    ___________________________________________________
!   |          |         |         |         |          |
!   |          |         |         |         |          |
!   |    o     |    o   w+    o    +e   o    |    o     |
!   |   WW     |    W    |    P    |    E    |    EE    |
!   |__________|_________|_________|_________|__________|
!
!   |<--S_uu-->|<--S_u-->|<--S_c-->|<--S_d-->|<--S_dd-->|  1.Flow: Left --> Right
!                                           
!   |<--S_dd-->|<--S_d-->|<--S_c-->|<--S_u-->|<--S_uu-->|  2.Flow: Left <-- Right 
!
!   Subscripts 'u' and 'd' in S refer to upstream and downstream resply.
!
!   LOGIC: 1. If the flow is from left to right across a face (say, face e), then
!
!                 |     { S_u + 2S_c}      {    S_c    }
!              u_e|   = {-----------}u_P - {-----------}u_W
!                 |up   { S_u + S_c }      { S_u + S_c } 
!
!          2. If the flow is from right to left across a face (say, face e), then
!
!                 |     { S_uu + 2S_u}      {     S_u    }
!              u_e|   = {------------}u_E - {------------}u_EE
!                 |up   { S_uu + S_u }      { S_uu + S_u } 
!        
!   - It should be noted that for u_w, the above formulae are still valid, provided the stencil  
!     is offset by one cell.
!   - These formulae are derived from:  u_face = u + (Grad u)(dot)(dS), where 
!
!     'u' and 'Grad u' are cellcentered value and its gradient in the upstream cell resply, 
!      and dS is the displacement vector from the upstream cell-centroid to the face centroid.
!      'Grad u' is approximated by upwind differencing based on the direction of the wind.
!
!----------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! x-direction 
!------------------------------------------------------------------------------
      DO k = zc_start,zc_end  !1,nz-1
      DO j = yc_start,yc_end
      DO i = 1,nx-1

       alfa = alfa_upwind

       imm = MAX(i-2,0)     
       ipp = MIN(i+2,nx)

       USign  = SIGN(1.0_CGREAL,face_u(i,  j,k)) 
       USignP = SIGN(1.0_CGREAL,face_u(i+1,j,k))

       vele =              fx(i+1)  *vel(i+1,j,k)                            &
            + ( 1.0_CGREAL-fx(i+1)) *vel(i,  j,k) 

       velw = (             fx(i)     *vel(i,j,k)                            &
            +  ( 1.0_CGREAL-fx(i)   ) *vel(i-1,j,k) )
!
!      2nd Upwind differencing added by Rupesh
!
       edxWeightm1 = ( dxinv(i-1)+2.0_CGREAL*dxinv(i) ) /   &
                                ( dxinv(i-1)+dxinv(i) )
       edxWeightm2 = dxinv(i) / ( dxinv(i-1)+dxinv(i) )
       edxWeightp1 = ( dxinv(ipp)+2.0_CGREAL*dxinv(i+1) ) / &
                                ( dxinv(ipp)+dxinv(i+1) )
       edxWeightp2 = dxinv(i+1) / ( dxinv(ipp)+dxinv(i+1) )

       wdxWeightm1 = ( dxinv(imm)+2.0_CGREAL*dxinv(i-1) ) / &
                                ( dxinv(imm)+dxinv(i-1) )
       wdxWeightm2 = dxinv(i-1) / ( dxinv(imm)+dxinv(i-1) )
       wdxWeightp1 = ( dxinv(i+1)+2.0_CGREAL*dxinv(i) ) /   &
                                ( dxinv(i+1)+dxinv(i) )
       wdxWeightp2 = dxinv(i) / ( dxinv(i+1)+dxinv(i) )


       IF (i==nx-1 .and. (UsignP < -0.9d0))then
           vele_Up2 =  vele
       else 
           vele_Up2 = ( ( 1.0_CGREAL+UsignP ) * ( edxWeightm1*vel(i,j,k)   - edxWeightm2*vel(i-1,j,k) )   &
                      + ( 1.0_CGREAL-UsignP ) * ( edxWeightp1*vel(i+1,j,k) - edxWeightp2*vel(ipp,j,k) )   &
                      )*0.5_CGREAL
       endif  

       IF (i==1    .and. (Usign  > 0.9d0)) then
           velw_Up2 =  velw
       ELSE
           velw_Up2 = ( ( 1.0_CGREAL+Usign )  * ( wdxWeightm1*vel(i-1,j,k) - wdxWeightm2*vel(imm,j,k) )   &
                      + ( 1.0_CGREAL-Usign )  * ( wdxWeightp1*vel(i,j,k)   - wdxWeightp2*vel(i+1,j,k) )   &
                      )*0.5_CGREAL
       ENDIF

!
!      1st order Upwind differencing added by Ye-------------------------------
!
       IF (i==nx-1 .and. (UsignP < -0.9d0))then
           vele_Up1 =  vele
       else
           vele_Up1 = ( ( 1.0_CGREAL+UsignP ) * ( vel(i,j,k) )   &
                      + ( 1.0_CGREAL-UsignP ) * ( vel(i+1,j,k) )   &
                      )*0.5_CGREAL
       endif

       IF (i==1    .and. (Usign  > 0.9d0)) then
           velw_Up1 =  velw
       ELSE
           velw_Up1 = ( ( 1.0_CGREAL+Usign )  * ( vel(i-1,j,k) )   &
                      + ( 1.0_CGREAL-Usign )  * ( vel(i,j,k) )   &
                      )*0.5_CGREAL
       ENDIF

       !Ye,hybrid
       IF ( (i.ge.25).and.(i.le.375) )then
          velw_Up = velw_Up2
          vele_Up = vele_Up2
       ELSE
          velw_Up = velw_Up1
          vele_Up = vele_Up1
       ENDIF

!-------------------------------------------------------------------------

         nlvel(i,j,k) = (  ( (1.0_CGREAL - alfa)*vele + alfa*vele_Up )*          &
                                                       face_u(i+1,j,k)           & 
                         - ( (1.0_CGREAL - alfa)*velw + alfa*velw_Up )*          &
                                                       face_u(i,j,k)             &
                        )*dxinv(i)           ! Hybrid definition added by Rupesh
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

!------------------------------------------------------------------------------
! y-direction
!------------------------------------------------------------------------------
     DO k = zc_start,zc_end  !1,nz-1
     DO j = yc_start,yc_end
     DO i = 1,nx-1

 
       alfa = alfa_upwind
       jmm = MAX(j-2,0)     
       jpp = MIN(j+2,ny)

       VSign  = SIGN(1.0_CGREAL,face_v(i,j,k)) 
       VSignP = SIGN(1.0_CGREAL,face_v(i,j+1,k))

       veln = (             fy(j+1)   *vel(i,j+1,k)                            &
            +  ( 1.0_CGREAL-fy(j+1) ) *vel(i,j,k)   )

       vels = (             fy(j)     *vel(i,j,k)                              &
            +  ( 1.0_CGREAL-fy(j)   ) *vel(i,j-1,k) )
!
!      2nd Upwind differencing added by Rupesh
!
       ndxWeightm1 = ( dyinv(j-1)+2.0_CGREAL*dyinv(j) ) /   &
                                ( dyinv(j-1)+dyinv(j) )
       ndxWeightm2 = dyinv(j) / ( dyinv(j-1)+dyinv(j) )
       ndxWeightp1 = ( dyinv(jpp)+2.0_CGREAL*dyinv(j+1) ) / &
                                ( dyinv(jpp)+dyinv(j+1) )
       ndxWeightp2 = dyinv(j+1) / ( dyinv(jpp)+dyinv(j+1) )

       sdxWeightm1 = ( dyinv(jmm)+2.0_CGREAL*dyinv(j-1) ) / &
                                ( dyinv(jmm)+dyinv(j-1) )
       sdxWeightm2 = dyinv(j-1) / ( dyinv(jmm)+dyinv(j-1) )
       sdxWeightp1 = ( dyinv(j+1)+2.0_CGREAL*dyinv(j) ) /   &
                                ( dyinv(j+1)+dyinv(j) )
       sdxWeightp2 = dyinv(j) / ( dyinv(j+1)+dyinv(j) )


       IF (j==ny-1 .and. (VsignP < -0.9d0))  then
          veln_Up2 =  veln
       ELSE
        veln_Up2 = ( ( 1.0_CGREAL+VsignP )  * ( ndxWeightm1*vel(i,j,k)   - ndxWeightm2*vel(i,j-1,k) )   &
                   + ( 1.0_CGREAL-VsignP )  * ( ndxWeightp1*vel(i,j+1,k) - ndxWeightp2*vel(i,jpp,k) )   &
                   )*0.5_CGREAL
       ENDIF


       IF (j==1    .and. (Vsign  > 0.9d0))   then
          vels_Up2 =  vels
       ELSE
        vels_Up2 = ( ( 1.0_CGREAL+Vsign )   * ( sdxWeightm1*vel(i,j-1,k) - sdxWeightm2*vel(i,jmm,k) )    &
                   + ( 1.0_CGREAL-Vsign )   * ( sdxWeightp1*vel(i,j,k)   - sdxWeightp2*vel(i,j+1,k) )    &
                   )*0.5_CGREAL
       ENDIF

!
!      1st Upwind differencing added by Ye-------------------------------
!
       IF (j==ny-1 .and. (VsignP < -0.9d0))  then
            veln_Up1 =  veln
       ELSE
        veln_Up1 = ( ( 1.0_CGREAL+VsignP )  * ( vel(i,j,k) )   &
                   + ( 1.0_CGREAL-VsignP )  * ( vel(i,j+1,k) )   &
                   )*0.5_CGREAL
       ENDIF


       IF (j==1    .and. (Vsign  > 0.9d0))   then
          vels_Up1 =  vels
       ELSE
        vels_Up1 = ( ( 1.0_CGREAL+Vsign )   * ( vel(i,j-1,k) )    &
                   + ( 1.0_CGREAL-Vsign )   * ( vel(i,j,k) )    &
                   )*0.5_CGREAL
       ENDIF

       !Ye,hybrid
       IF ( (i.ge.25).and.(i.le.375) )then
          veln_Up = veln_Up2
          vels_Up = vels_Up2
       ELSE
          veln_Up = veln_Up1
          vels_Up = vels_Up1
       ENDIF

!-------------------------------------------------------------------------

        nlvel(i,j,k) =nlvel(i,j,k) + (  ( (1.0_CGREAL - alfa)*veln + alfa*veln_Up )*          &
                                                        face_v(i,j+1,k)                       & 
                            - ( (1.0_CGREAL - alfa)*vels + alfa*vels_Up )*                    &
                                                        face_v(i,j,k)                         &
                                     )*dyinv(j)           ! Hybrid definition added by Rupesh
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

!=========== z direction
     DO k = zc_start,zc_end  !1,nz-1
     DO j = yc_start,yc_end
     DO i = 1,nx-1

       alfa = alfa_upwind

       kmm = MAX(k-2,0)     
       kpp = MIN(k+2,nz)

       WSign  = SIGN(1.0_CGREAL,face_w(i,j,k)) 
       WSignP = SIGN(1.0_CGREAL,face_w(i,j,k+1))

       velf = (             fz(k+1)   *vel(i,j,k+1)                            &
            +  ( 1.0_CGREAL-fz(k+1) ) *vel(i,j,k)   )

       velb = (             fz(k)     *vel(i,j,k)                              &
            +  ( 1.0_CGREAL-fz(k)   ) *vel(i,j,k-1) )
!
!      2nd Upwind differencing added by Rupesh
!
       fdxWeightm1 = ( dzinv(k-1)+2.0_CGREAL*dzinv(k) ) /   &
                                ( dzinv(k-1)+dzinv(k) )
       fdxWeightm2 = dzinv(k) / ( dzinv(k-1)+dzinv(k) )
       fdxWeightp1 = ( dzinv(kpp)+2.0_CGREAL*dzinv(k+1) ) / &
                                ( dzinv(kpp)+dzinv(k+1) )
       fdxWeightp2 = dzinv(k+1) / ( dzinv(kpp)+dzinv(k+1) )
!

       bdxWeightm1 = ( dzinv(kmm)+2.0_CGREAL*dzinv(k-1) ) / &
                               ( dzinv(kmm)+dzinv(k-1) )
       bdxWeightm2 = dzinv(k-1) / ( dzinv(kmm)+dzinv(k-1) )
       bdxWeightp1 = ( dzinv(k+1)+2.0_CGREAL*dzinv(k) ) /   &
                                ( dzinv(k+1)+dzinv(k) )
       bdxWeightp2 = dzinv(k) / ( dzinv(k+1)+dzinv(k) )
!


       IF      (k==nz-1 .and. (WsignP < -0.9d0))   then
          velf_Up2 =  velf
       ELSE  
        velf_Up2 = ( ( 1.0_CGREAL+wsignP )  * ( fdxWeightm1*vel(i,j,k)   - fdxWeightm2*vel(i,j,k-1) )   &
                   + ( 1.0_CGREAL-wsignP )  * ( fdxWeightp1*vel(i,j,k+1) - fdxWeightp2*vel(i,j,kpp) )   &
                   )*0.5_CGREAL

       ENDIF


       IF     (k==1    .and. (Wsign  > 0.9d0))   then
          velb_Up2 =  velb
       ELSE

        velb_Up2 = ( ( 1.0_CGREAL+wsign )   * ( bdxWeightm1*vel(i,j,k-1) - bdxWeightm2*vel(i,j,kmm) )    &
                   + ( 1.0_CGREAL-wsign )   * ( bdxWeightp1*vel(i,j,k)   - bdxWeightp2*vel(i,j,k+1) )    &
                   )*0.5_CGREAL
       ENDIF

!
!      1st Upwind differencing added by Ye
!
       IF      (k==nz-1 .and. (WsignP < -0.9d0))   then
          velf_Up1 =  velf
       ELSE
        velf_Up1 = ( ( 1.0_CGREAL+wsignP )  * ( vel(i,j,k) )   &
                   + ( 1.0_CGREAL-wsignP )  * ( vel(i,j,k+1) )   &
                   )*0.5_CGREAL

       ENDIF


       IF     (k==1    .and. (Wsign  > 0.9d0))   then
          velb_Up1 =  velb
       ELSE

        velb_Up1 = ( ( 1.0_CGREAL+wsign )   * ( vel(i,j,k-1) )    &
                   + ( 1.0_CGREAL-wsign )   * ( vel(i,j,k) )    &
                   )*0.5_CGREAL
       ENDIF      

       !Ye,hybrid
       IF ( (i.ge.25).and.(i.le.375) )then
          velf_Up = velf_Up2
          velb_Up = velb_Up2
       ELSE
          velf_Up = velf_Up1
          velb_Up = velb_Up1
       ENDIF
 
!----------------------------------------------------------------------------- 
        nlvel(i,j,k) =nlvel(i,j,k) + (  ( (1.0_CGREAL - alfa)*velf + alfa*velf_Up )*          &
                                                        face_w(i,j,k+1)                       & 
                            - ( (1.0_CGREAL - alfa)*velb + alfa*velb_Up )*                    &
                                                        face_w(i,j,k)                         &
                                     )*dzinv(k)           ! Hybrid definition added by Rupesh
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

    END SUBROUTINE rhs_advec

!------------------------------------------------------------------------------

    SUBROUTINE rhs_diff(vel,nlvel) 
    
      USE global_parameters
      USE flow_parameters
      USE flow_arrays
      USE grid_arrays
      USE boundary_arrays
      USE solver_arrays
      USE MPI_module

      IMPLICIT NONE

!... parameters
   
      REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(IN)    :: vel     
      REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(INOUT) :: nlvel

!... loop  variables

      INTEGER              :: i,j,k

!... local variables

      REAL(KIND=CGREAL)    :: aacx,aamx,aapx,aacy,aamy,aapy,aacz,aamz,aapz
      REAL(KIND=CGREAL)    :: diffvel
      REAL(KIND=CGREAL)    :: nuE,nuW,nuS,nuN,nuB,nuF

!******************************************************************************

!------------------------------------------------------------------------------
! Diffusive terms
! nlvel = d/dx_j [(1/Re+nut) d(vel)/dx_j ]
!
! Note: Since amx_ad coefficients have already a negative sign in them (-1/2 dt)
!       diffvel needs only to be subtracted in computing nlvel term
!------------------------------------------------------------------------------
     
      DO k = zc_start,zc_end  !1,nz-1
      DO j = yc_start,yc_end
      DO i = 1,nx-1
        nuE = (             fx(i+1)   *viscTot(i+1,j,k)                            &
            +  ( 1.0_CGREAL-fx(i+1) ) *viscTot(i,j,k)   )

        nuW = (             fx(i)     *viscTot(i,j,k)                              &
            +  ( 1.0_CGREAL-fx(i)   ) *viscTot(i-1,j,k) )

        nuN = (             fy(j+1)   *viscTot(i,j+1,k)                            &
            +  ( 1.0_CGREAL-fy(j+1) ) *viscTot(i,j,k)   )

        nuS = (             fy(j)     *viscTot(i,j,k)                              &
            +  ( 1.0_CGREAL-fy(j)   ) *viscTot(i,j-1,k) )

        nuf = (             fz(k+1)   *viscTot(i,j,k+1)                            &
            +  ( 1.0_CGREAL-fz(k+1) ) *viscTot(i,j,k)   )

        nub = (             fz(k)     *viscTot(i,j,k)                              &
            +  ( 1.0_CGREAL-fz(k)   ) *viscTot(i,j,k-1) )
      

        aamx =  dxcinv(i)*dxinv(i) 
        aapx =  dxcinv(i+1)*dxinv(i) 
        aamx =  aamx*nuW
        aapx =  aapx*nuE
        aacx = - ( aamx + aapx )

        aamy =  dycinv(j)*dyinv(j) 
        aapy =  dycinv(j+1)*dyinv(j) 
        aamy =  aamy * nuS
        aapy =  aapy * nuN
        aacy = - ( aamy + aapy )

        aamz =  dzcinv(k)*dzinv(k) 
        aapz =  dzcinv(k+1)*dzinv(k) 
        aamz =  aamz * nub
        aapz =  aapz * nuf
        aacz = - ( aamz + aapz )

        diffvel  = aamx*vel(i-1,j,k) + aacx*vel(i,j,k) + aapx*vel(i+1,j,k)  &
                 + aamy*vel(i,j-1,k) + aacy*vel(i,j,k) + aapy*vel(i,j+1,k)  &
                 + aamz*vel(i,j,k-1) + aacz*vel(i,j,k) + aapz*vel(i,j,k+1)
                    
        nlvel(i,j,k ) = vel(i,j,k)  + 0.5_CGREAL*dt*(-nlvel(i,j,k) + diffvel)
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k 

    END SUBROUTINE rhs_diff
!------------------------------------------------------------------------------

    SUBROUTINE rhs_vankan
    
      USE global_parameters
      USE flow_parameters
      USE flow_arrays
      USE grid_arrays
      USE boundary_arrays
      USE multiuse_arrays
      USE MPI_module

      IMPLICIT NONE

!... loop  variables

      INTEGER           :: i,j,k

!... local variables

      REAL(KIND=CGREAL) :: pf,pb,pe,pw,pn,ps,pgx,pgy,pgz
!******************************************************************************

!------------------------------------------------------------------------------      
! add pressure gradient
!------------------------------------------------------------------------------
    DO k = zc_start,zc_end   !1,nz-1 
    DO j = 1,ny-1
    DO i = 1,nx-1
      pe = ( fx(i+1)*p(i+1,j,k ) + (1.0_CGREAL-fx(i+1))*p(i,j,k )   )
      pw = ( fx(i)  *p(i,j,k )   + (1.0_CGREAL-fx(i))  *p(i-1,j,k ) ) 
      pn = ( fy(j+1)*p(i,j+1,k ) + (1.0_CGREAL-fy(j+1))*p(i,j,k )   ) 
      ps = ( fy(j)  *p(i,j,k )   + (1.0_CGREAL-fy(j))  *p(i,j-1,k ) ) 
      pf = ( fz(k+1)*p(i,j,k+1)  + (1.0_CGREAL-fz(k+1))*p(i,j ,k)   ) 
      pb = ( fz(k)  *p(i,j,k )   + (1.0_CGREAL-fz(k))  *p(i,j,k-1 ) ) 
       
      pgx= (pe-pw)*dxinv(i)
      pgy= (pn-ps)*dyinv(j)
      pgz= (pf-pb)*dzinv(k)
      
      nlu(i,j,k ) = nlu(i,j,k ) - dt*pgx*(1.0_CGREAL-REAL(iblank(i,j,k ),KIND=CGREAL))
      nlv(i,j,k ) = nlv(i,j,k ) - dt*pgy*(1.0_CGREAL-REAL(iblank(i,j,k ),KIND=CGREAL))
      nlw(i,j,k ) = nlw(i,j,k ) - dt*pgz*(1.0_CGREAL-REAL(iblank(i,j,k ),KIND=CGREAL))      
    ENDDO
    ENDDO
    ENDDO
   END SUBROUTINE  rhs_vankan
!------------------------------------------------------------------------------


!******************************************************************************
!
! Purpose: generalized kernel to compute the coefficients 
!          for the diffusion term
!
! Description: none.
!
! Input: grid size, total viscosity
!
! Output: amx_ad, amy_ad, amz_ad, apx_ad, apy_ad, apz_ad. 
!
!
!******************************************************************************
!------------------------------------------------------------------------------
   SUBROUTINE set_solve_ad()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE GCM_arrays
    USE flow_arrays
    USE MPI_module

    IMPLICIT NONE

!... Loop variables 
    
    INTEGER :: i,j,n,iRow,k
    
!... Local variables
   
    INTEGER :: iFr,jFr
    REAL(KIND=CGREAL) :: rFreshCell, rnDim
    REAL(KIND=CGREAL) :: nuE,nuW,nuS,nuN,nuB,nuF
!******************************************************************************

!------------------------------------------------------------------------------ 
! Initialize coefficients 
!------------------------------------------------------------------------------ 

    amx_ad = 0.0_CGREAL
    apx_ad = 0.0_CGREAL
    
    amy_ad = 0.0_CGREAL
    apy_ad = 0.0_CGREAL
    
    amz_ad = 0.0_CGREAL
    apz_ad = 0.0_CGREAL

    !nuE = (             fx(i+1)   *viscTot(i+1,j,k)                            &
    !    +  ( 1.0_CGREAL-fx(i+1) ) *viscTot(i,j,k)   )

    !nuW = (             fx(i)     *viscTot(i,j,k)                              &
    !    +  ( 1.0_CGREAL-fx(i)   ) *viscTot(i-1,j,k) )

    !nuN = (             fy(j+1)   *viscTot(i,j+1,k)                            &
    !    +  ( 1.0_CGREAL-fy(j+1) ) *viscTot(i,j,k)   )

    !nuS = (             fy(j)     *viscTot(i,j,k)                              &
    !    +  ( 1.0_CGREAL-fy(j)   ) *viscTot(i,j-1,k) )
      
    !nuF = (             fz(k+1)   *viscTot(i,j,k+1)                            &
    !    +  ( 1.0_CGREAL-fz(k+1) ) *viscTot(i,j,k)   )

    !nuB = (             fz(k)     *viscTot(i,j,k)                              &
    !    +  ( 1.0_CGREAL-fz(k)   ) *viscTot(i,j,k-1) )
             
    nuE = reinv
    nuW = reinv 
    nuN = reinv
    nuS = reinv
    nuF = reinv
    nuB = reinv
     
    do i=1,nx-1  
       amx_ad(i) = ( dxcinv(i)  )*dxinv(i)
       apx_ad(i) = ( dxcinv(i+1))*dxinv(i)

       amx_ad(i) =- (0.50_CGREAL*dt*nuW)*amx_ad(i)
       apx_ad(i) =- (0.50_CGREAL*dt*nuE)*apx_ad(i)
    enddo
      
    do j=yc_start,yc_end
       amy_ad(j) =  dycinv(j)*dyinv(j)
       apy_ad(j) =  dycinv(j+1)*dyinv(j)

       amy_ad(j) =- (0.50_CGREAL*dt*nuS)*amy_ad(j)
       apy_ad(j) =- (0.50_CGREAL*dt*nuN)*apy_ad(j)
    enddo
    
    do k=zc_start,zc_end  !1,nz-1         
       amz_ad(k) =  dzcinv(k)*dzinv(k)
       apz_ad(k) =  dzcinv(k+1)*dzinv(k)

       amz_ad(k) =- (0.50_CGREAL*dt*nuB)*amz_ad(k)
       apz_ad(k) =- (0.50_CGREAL*dt*nuF)*apz_ad(k)
    enddo
 
  END SUBROUTINE set_solve_ad   
!------------------------------------------------------------------------------
