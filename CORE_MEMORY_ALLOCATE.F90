!--------------------------------------------------
!   SUBROUTINE allocate_memory()
!   SUBROUTINE deallocate_memory()
!--------------------------------------------------

   SUBROUTINE allocate_memory()

    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE grid_arrays
    USE solver_arrays
    USE stat_arrays
    USE stat_vort_arrays
    USE MPI_module

    IMPLICIT NONE
    
    INTEGER :: nmax

    ! Note that each processor has the full grid
    ALLOCATE(x(0:nx+1)     ,y(0:ny+1)     ,z(0:nz+1)     )
    ALLOCATE(xc(0:nx+1)    ,yc(0:ny+1)    ,zc(0:nz+1)    )
    ALLOCATE(dx(0:nx+1)    ,dy(0:ny+1)    ,dz(0:nz+1)    )
    ALLOCATE(dxc(0:nx+1)   ,dyc(0:ny+1)   ,dzc(0:nz+1)   )
    ALLOCATE(dxinv(0:nx+1) ,dyinv(0:ny+1) ,dzinv(0:nz+1) )
    ALLOCATE(dxcinv(0:nx+1),dycinv(0:ny+1),dzcinv(0:nz+1))
    ALLOCATE(fx(0:nx+1)    ,fy(0:ny+1)    ,fz(0:nz+1)    )

    ! Note that each processor has only a subdomain, or slab, 
    ! including ghost slices
    ALLOCATE(u(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(v(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(w(0:nx+1,yb1:yb2,zb1:zb2))

    ALLOCATE(face_u(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(face_v(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(face_w(0:nx+1,yb1:yb2,zb1:zb2))

    ALLOCATE(p(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(pPrime(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(p_temp(0:nx+1,yb1:yb2,zb1:zb2))
 
    !ALLOCATE(dpdx(0:nx+1,0:ny+1,zb1:zb2))
    !ALLOCATE(dpdy(0:nx+1,0:ny+1,zb1:zb2))
    !ALLOCATE(dpdz(0:nx+1,0:ny+1,zb1:zb2))
    
    ALLOCATE(pgradx1(0:ny+1,0:nz+1))  !full length in z
    ALLOCATE(pgradx2(0:ny+1,0:nz+1))  ! 
    ALLOCATE(pgrady1(0:nx+1,0:nz+1))  !
    ALLOCATE(pgrady2(0:nx+1,0:nz+1))  !
    ALLOCATE(pgradz1(0:nx+1,0:ny+1))
    ALLOCATE(pgradz2(0:nx+1,0:ny+1))

    
    ALLOCATE(nlu(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(nlv(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(nlw(0:nx+1,yb1:yb2,zb1:zb2))


    !ALLOCATE(nluold(0:nx+1,0:ny+1,zb1:zb2))
    !ALLOCATE(nlvold(0:nx+1,0:ny+1,zb1:zb2))
    !ALLOCATE(nlwold(0:nx+1,0:ny+1,zb1:zb2))

    ALLOCATE(bcxu(0:1,   yb1:yb2,zb1:zb2))
    ALLOCATE(bcxv(0:1,   yb1:yb2,zb1:zb2))
    ALLOCATE(bcxw(0:1,   yb1:yb2,zb1:zb2))
    !ALLOCATE(bcxp(0:1,   0:ny+1,zb1:zb2))
    ALLOCATE(bcxp(0:1,   0:ny+1, 0:nz+1))  !full length in z

    ALLOCATE(bcyu(0:nx+1,0:1,   zb1:zb2))
    ALLOCATE(bcyv(0:nx+1,0:1,   zb1:zb2))
    ALLOCATE(bcyw(0:nx+1,0:1,   zb1:zb2))
    !ALLOCATE(bcyp(0:nx+1,0:1,   zb1:zb2))
    ALLOCATE(bcyp(0:nx+1,0:1,    0:nz+1))  !full length in z

    ALLOCATE(bczu(0:nx+1,yb1:yb2,0:1   ))
    ALLOCATE(bczv(0:nx+1,yb1:yb2,0:1   ))
    ALLOCATE(bczw(0:nx+1,yb1:yb2,0:1   ))
    ALLOCATE(bczp(0:nx+1,0:ny+1,0:1   )) !full length in y


    ALLOCATE(div(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(lmd(0:nx+1,yb1:yb2,zb1:zb2))

    ALLOCATE(amx(0:nx+1),acx(0:nx+1),apx(0:nx+1))
    ALLOCATE(amy(0:ny+1),acy(0:ny+1),apy(0:ny+1))
    !The following solver arrays still have full size
    ALLOCATE(amz(0:nz+1),acz(0:nz+1),apz(0:nz+1)) 

    nmax = max(max(nx,ny),nz)   !nx+ny+nz
    ALLOCATE(rhs(0:nmax+1),dummy(0:nmax+1))
    !ALLOCATE(face1(0:nmax+1),face2(0:nmax+1))

    ALLOCATE(amx_ad(0:nx+1))
    ALLOCATE(apx_ad(0:nx+1))
    
    ALLOCATE(amy_ad(0:ny+1))
    ALLOCATE(apy_ad(0:ny+1))

    ALLOCATE(amz_ad(zb1:zb2))
    ALLOCATE(apz_ad(zb1:zb2))

    ALLOCATE(uOld(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(vOld(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(wOld(0:nx+1,yb1:yb2,zb1:zb2))
    
!----------------------------------------------
!   Arrays pertinent to viscosity
!----------------------------------------------
    ALLOCATE(viscTot(0:nx+1,yb1:yb2,zb1:zb2),STAT=iErr)

    x       = 0.0_CGREAL
    y       = 0.0_CGREAL
    z       = 0.0_CGREAL
    
    dx      = 0.0_CGREAL
    dy      = 0.0_CGREAL
    dz      = 0.0_CGREAL 
    
    dxc     = 0.0_CGREAL
    dyc     = 0.0_CGREAL
    dzc     = 0.0_CGREAL 
    
    dxinv   = 0.0_CGREAL
    dyinv   = 0.0_CGREAL
    dzinv   = 0.0_CGREAL 
    
    dxcinv  = 0.0_CGREAL
    dycinv  = 0.0_CGREAL
    dzcinv  = 0.0_CGREAL 
    
    fx      = 0.0_CGREAL
    fy      = 0.0_CGREAL
    fz      = 0.0_CGREAL
       
    !  iblank      = 0 
    !  dead_cell   = 0
    
    u       = 0.0_CGREAL
    v       = 0.0_CGREAL
    w       = 0.0_CGREAL

    face_u  = 0.0_CGREAL
    face_v  = 0.0_CGREAL
    face_w  = 0.0_CGREAL

    p       = 0.0_CGREAL
    pPrime  = 0.0_CGREAL
    
    pgradx1 = 0.0_CGREAL
    pgradx2 = 0.0_CGREAL
    pgrady1 = 0.0_CGREAL
    pgrady2 = 0.0_CGREAL
    pgradz1 = 0.0_CGREAL
    pgradz2 = 0.0_CGREAL

    nlu     = 0.0_CGREAL
    nlv     = 0.0_CGREAL
    nlw     = 0.0_CGREAL

    !nluold  = 0.0_CGREAL
    !nlvold  = 0.0_CGREAL
    !nlwold  = 0.0_CGREAL

    uOld    = 0.0_CGREAL
    vOld    = 0.0_CGREAL
    wOld    = 0.0_CGREAL

    bcxu    = 0.0_CGREAL
    bcyu    = 0.0_CGREAL
    bczu    = 0.0_CGREAL
    
    bcxv    = 0.0_CGREAL
    bcyv    = 0.0_CGREAL
    bczv    = 0.0_CGREAL
    
    bcxw    = 0.0_CGREAL
    bcyw    = 0.0_CGREAL
    bczw    = 0.0_CGREAL
    
    amx     = 0.0_CGREAL
    apx     = 0.0_CGREAL
    acx     = 0.0_CGREAL
    
    amy     = 0.0_CGREAL
    apy     = 0.0_CGREAL
    acy     = 0.0_CGREAL
    
    amz     = 0.0_CGREAL
    apz     = 0.0_CGREAL
    acz     = 0.0_CGREAL
    
    rhs     = 0.0_CGREAL
    dummy   = 0.0_CGREAL
    
    !    face1   = 0.0_CGREAL
    !    face2   = 0.0_CGREAL
    
    viscTot = 0.0_CGREAL

  END SUBROUTINE allocate_memory
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE deallocate_memory()

    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE grid_arrays
    USE solver_arrays
    USE stat_arrays
    USE stat_vort_arrays
    USE mg_module
    USE usr_module 
    USE MPI_module
    
    IMPLICIT NONE

    DEALLOCATE(x  ,y  ,z)
    DEALLOCATE(xc ,yc ,zc)
    DEALLOCATE(dx ,dy ,dz)
    DEALLOCATE(dxc,dyc,dzc)
    DEALLOCATE(dxinv ,dyinv ,dzinv)
    DEALLOCATE(dxcinv,dycinv,dzcinv)
    DEALLOCATE(fx    ,fy    ,fz)

    DEALLOCATE(iblank)
    DEALLOCATE(dead_cell)

    DEALLOCATE(u)
    DEALLOCATE(v)
    DEALLOCATE(w)

    DEALLOCATE(face_u)
    DEALLOCATE(face_v)
    DEALLOCATE(face_w)

    DEALLOCATE(p)
    DEALLOCATE(pPrime)
    
    DEALLOCATE(pgradx1)
    DEALLOCATE(pgradx2)
    DEALLOCATE(pgrady1)
    DEALLOCATE(pgrady2)
    DEALLOCATE(pgradz1)
    DEALLOCATE(pgradz2)


    DEALLOCATE(nlu)
    DEALLOCATE(nlv)
    DEALLOCATE(nlw)

    !DEALLOCATE(nluold)
    !DEALLOCATE(nlvold)
    !DEALLOCATE(nlwold)
    
    DEALLOCATE(div )
    DEALLOCATE(lmd )

    DEALLOCATE(bcxu)
    DEALLOCATE(bcxv)
    DEALLOCATE(bcxw)
    DEALLOCATE(bcyu)
    DEALLOCATE(bcyv)
    DEALLOCATE(bcyw)
    DEALLOCATE(bczu)
    DEALLOCATE(bczv)
    DEALLOCATE(bczw)

    
    DEALLOCATE(rhs,dummy)
    !DEALLOCATE(face1,face2)
    DEALLOCATE(amx,acx,apx)
    DEALLOCATE(amy,acy,apy)
    DEALLOCATE(amz,acz,apz)
   
 
    DEALLOCATE(amx_ad)
    DEALLOCATE(apx_ad)    
    DEALLOCATE(amy_ad)
    DEALLOCATE(apy_ad)
    DEALLOCATE(amz_ad)
    DEALLOCATE(apz_ad)

    DEALLOCATE(uOld)
    DEALLOCATE(vOld)
    DEALLOCATE(wOld)
    
    IF (it_solver_type == IT_SOLVER_TYPE_MG) THEN
      DEALLOCATE (mgrid_I, mgrid_J, mgrid_K)
      DEALLOCATE (mg_array)
    END IF ! it_solver_type

  END SUBROUTINE deallocate_memory
!------------------------------------------------------------------------------
