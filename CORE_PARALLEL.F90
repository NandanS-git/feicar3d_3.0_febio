 !-------------------------------------------------- 
  SUBROUTINE MPI_Initialize

    USE MPI_module

    IMPLICIT NONE
    integer:: i 
    character*(MPI_MAX_PROCESSOR_NAME) hostname
    integer  :: result_len

!................... Initialize the MPI environment ............
 
    call MPI_INIT(ierr)

    ! Create MPI world include processors for both fluid and solid solvers.
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nProcs_g,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,lProc_g, ierr)    
    call MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierr)
  
    ! read input for MPI setup  
    open(unit=1,file='mpi_setup.dat')
    read(1,*)
    read(1,*)nProcY, nProcZ
    close(1)

    if (lProc_g.eq.1)then
       write(*,*)'nProcY/Z = ',nProcY,nProcZ
    endif
  
    CALL read_inputs_body() !read from 'Canonical_Body_In.dat'

    ! Assign last processors, from (nProcs-nElastic) to (nProcs-1), for the elastic bodies
    do i=1,nElastic
       !proc_s(i) = nProcs_g-nElastic + i-1 
       proc_s(i) = i-1
    enddo 

    ! create a group and communicator for flow solver
    ! call MPI_GROUP_EXCL(world_group, 1, proc_s, flow_group, ierr)
    
    call MPI_GROUP_EXCL(world_group, nElastic, proc_s, flow_group, ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD, flow_group, flow_comm, ierr)
 
    !if (lProc_g .lt. nProcs_g-nElastic) then
    if (lProc_g .gt. nElastic-1) then

       ! in the flow group, nProcs is number of processors, and lProc is the procssor ID.   
       call MPI_COMM_SIZE (flow_comm,  nProcs,  ierr)
       call MPI_GROUP_RANK(flow_group, lProc, ierr)

       call mpi_get_processor_name(hostname, result_len, ierr)
       write(*,'(a,I5,a,I5,a,I5,a,a30)')"Flow: process ", lProc, " of ", nProcs_g,  &
               " (lProc_g=",lProc_g, ") running on node: ", hostname
    else
       lProc = -1  !lProc_g
       call mpi_get_processor_name(hostname, result_len, ierr)
       write(*,'(a,I5,a,I5,a,a30)')'Solid: process ',lProc_g, " of ", nProcs_g, ' running on node ', hostname
    endif

    PROC_M=0  ! the master processor in the flow group.
    if(lProc==PROC_M) print*,'PROC_M: lProc_g=',lProc_g, nProcs
 
    !Ye, for uniform slices distribution only
    jProc = MOD(lProc,nProcY)
    kProc = lProc/nProcY

    !The leftmost and rightmost processors in CPU array. These can change for
    !for multigrids
    lProc_leftmost  =  0
    lProc_rightmost =  nProcY-1

    !The topest and bottomest processors in CPU array. These can change for for
    !multigrids
    lProc_btommost = 0
    lProc_topmost  = nProcZ-1

    !The Left/Right/Top/Bottom processors.
    lProc_left  =  kProc*nProcY+(jProc-1)
    lProc_right =  kProc*nProcY+(jProc+1)
    lProc_top   = (kProc+1)*nProcY+jProc
    lProc_bottom= (kProc-1)*nProcY+jProc

    !Debug
    !write(800+lProc,*)lProc,jProc,kProc
    !write(800+lProc,*)lProc_leftmost,lProc_left,lProc_right,lProc_rightmost
    !write(800+lProc,*)lProc_btommost,lProc_bottom,lProc_top,lProc_topmost

  END SUBROUTINE MPI_Initialize

!--------------------------------------------------
  SUBROUTINE slice_assignment

    USE flow_parameters
    USE MPI_module

    IMPLICIT NONE
    
    INTEGER :: i, id, izc
    !********************

    !Ye, store number of sclieces in 2D decomposition
    ALLOCATE(numSlicesY(1:nProcs))
    ALLOCATE(numSlicesZ(1:nProcs))
    
    if(y_slab_unif == UNIFORM_SLAB) then
       !Create uniform slabs in y
       jSlices  = (ny-1) / nProcY
       DO i=1,nProcs
          numSlicesY(i) = jSlices
       ENDDO

       yc_start = jProc * jSlices + 1
       yc_end   = yc_start + jSlices - 1

    else
       ! Read numSlices from y_slices.dat
       OPEN(UNIT=11,FILE='y_slices.dat')
       DO i=1,nProcs
          read(11,*) id, izc, numSlicesY(i)
          if(i == lProc+1) then      ! for the current processor
             yc_start= izc
             jSlices = numSlicesY(i)
             yc_end  = yc_start+jSlices-1
          endif
       ENDDO
    endif

    if(z_slab_unif == UNIFORM_SLAB) then
       !Create uniform slabs in z
       kSlices  = (nz-1) / nProcZ
       DO i=1,nProcs
          numSlicesZ(i) = kSlices
       ENDDO
      
       zc_start = kProc * kSlices + 1
       zc_end   = zc_start + kSlices - 1

    else       
       ! Read numSlices from z_slices.dat
       OPEN(UNIT=11,FILE='z_slices.dat')       
       DO i=1,nProcs
          read(11,*) id, izc, numSlicesZ(i)
          if(i == lProc+1) then      ! for the current processor
             zc_start= izc
             kSlices = numSlicesZ(i)
             zc_end  = zc_start+kSlices-1
          endif   
       ENDDO
    endif
    
    ! Cell face start/end
    y_start  = yc_start
    y_end    = yc_end + 1

    z_start  = zc_start
    z_end    = zc_end + 1

    write(*,20) lProc, yc_start, yc_end, numSlicesY(lProc+1)
20  FORMAT(4X,'lProc #',I4,' Y-starts from ',I5,', ends at ', I5, ' with ', I5, ' slices.')

    call MPI_BARRIER(flow_comm,ierr)

    !print*, 'lProc #',lProc,'starts from ', zc_start,', ends at ', zc_end, ' with ', kSlices, 'slices.'
    write(*,30) lProc, zc_start, zc_end, numSlicesZ(lProc+1)    
30  FORMAT(4X,'lProc #',I4,' Z-starts from ',I5,', ends at ', I5, ' with ', I5, ' slices.')

    call MPI_BARRIER(flow_comm,ierr)

    !yb1/yb2   
    if (jProc.eq.0)then
       yb1 = 0
       yb2 = yc_end+2
    else
       yb1 = yc_start-2
       yb2 = yc_end+2
    endif

    !zb1/zb2   
    if (kProc.eq.0)then
       zb1 = 0
       zb2 = zc_end+2
    else
       zb1 = zc_start-2
       zb2 = zc_end+2
    endif

    write(*,40) lProc, jProc, yb1, yb2
40  FORMAT(4X,'lProc # ',I5,', jProc # ', I5,', yb1 =', I5,', yb2 =', I5)

    call MPI_BARRIER(flow_comm,ierr)

    write(*,50) lProc, kProc, zb1, zb2
50  FORMAT(4X,'lProc # ',I5,', kProc # ', I5,', zb1 =', I5,', zb2 =', I5)

    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !STOP

    ! Array allocation for each processor, including buffer slices
    !if(lProc .eq. 0) then      ! first slab
    !   zb1 = 0
    !   zb2 = zc_end+2
    !elseif(lProc .eq. nProcs-1) then   ! last slab
    !   zb1 = zc_start-2
    !   zb2 = zc_end+2
    !elseif(lProc .gt. 0 .and. lProc .lt. nProcs-1) then    ! middle slabs
    !   zb1 = zc_start-2
    !   zb2 = zc_end+2       
    !endif

  END SUBROUTINE slice_assignment
!---------------------------------------------------------------------
!-----------------------------------------------------------------------
! Use an arbitrary matrix to test data exchange between adjacent subdomains.
  SUBROUTINE test_send_receive()
    USE flow_parameters
    USE grid_arrays
    USE MPI_module

    IMPLICIT NONE

    INTEGER :: i,j,k,id,ia
    REAL(KIND=CGREAL),  DIMENSION(:,:,:), ALLOCATABLE :: A0
    INTEGER(KIND=INT_K),DIMENSION(:,:,:), ALLOCATABLE :: M0
    
    !Use a small matrix size
    nx=7
    ALLOCATE(A0(0:nx+1,yb1:yb2,zb1:zb2),STAT=iErr)
    ALLOCATE(M0(0:nx+1,yb1:yb2,zb1:zb2),STAT=iErr)
    
    IF ( iErr /= ERR_NONE ) THEN
       WRITE(STDOUT,*) &
            'test_send_receive: Memory Allocation Error for A0'
       STOP
    ENDIF ! ierr

    !Fill with arbitray data
    id = 100 + jProc*10 + kProc  !lProc

    A0(0:nx+1,yc_start:yc_end,zc_start:zc_end) = id
    A0(0:nx+1,yb1:yb1,zb1:zb1) = id - 0.2
    A0(0:nx+1,yb1+1:yb1+1,zb1+1:zb1+1) = id - 0.1
    A0(0:nx+1,yb2-1:yb2-1,zb2-1:zb2-1) = id + 0.1
    A0(0:nx+1,yb2:yb2,zb2:zb2) = id + 0.2

    write(6,*)'lProc',lProc, ' writing initial data to the disk ...'
    write(id,*)'VARIABLES="X","Y","Z","A0"'
    write(id,*)'ZONE F=POINT, I=',nx+2,', J=',yb2-yb1+1,' K=',zb2-zb1+1
    DO j=yb1,yb2
      write(id,*)'j = ', j 
      DO k=zb1,zb2
      !DO i=0,nx+1   
         !write(100+id,123)xc(i),yc(j),zc(k),A0(i,j,k)
         write(id,123) (A0(i,j,k),i=0,nx+1)
      !ENDDO
      ENDDO
    ENDDO

    ! Exchange data slices
    call send_receive_slices_real_y(A0, 0,nx+1,yb1,yb2,zb1,zb2,1)
      
    write(6,*)'lProc',lProc, ' writing exchanged data to the disk ...'
    write(100+id,*)'VARIABLES="X","Y","Z","A0"'
    write(100+id,*)'ZONE F=POINT, I=',nx+2,', J=',yb2-yb1+1,' K=',zb2-zb1+1
    DO j=yb1,yb2
      write(100+id,*)'j = ', j
      DO k=zb1,zb2
      !DO i=0,nx+1   
        !write(200+id,123)xc(i),yc(j),zc(k),A0(i,j,k)
         write(100+id,123) (A0(i,j,k),i=0,nx+1)
      !ENDDO
      ENDDO
    ENDDO

123  format(20(2x,F7.2))
    
    close(id)
    close(100+id)

    DEALLOCATE(A0)

    !---------
    ! test the integer exchange subroutine
    !---------
    !Fill with arbitray data
    ia = 100 + jProc*10 + kProc
    M0(0:nx+1,yc_start:yc_end,zc_start:zc_end) = ia
    !M0(0:nx+1,yb1:yb1,zb1:zb1) = ia - 2
    !M0(0:nx+1,yb1+1:yb1+1,zb1+1:zb1+1) = ia - 1
    !M0(0:nx+1,yb2-1:yb2-1,zb2-1:zb2-1) = ia + 1
    !M0(0:nx+1,yb2:yb2,zb2:zb2) = ia + 2
    !id = ia

    write(6,*)'lProc',lProc, ' writing initial data to the disk ...'
    write(200+id,*)'VARIABLES="X","Y","Z","M0"'
    write(200+id,*)'ZONE F=POINT, I=',nx+2,', J=',yb2-yb1+1,' K=',zb2-zb1+1
    DO j=yb1,yb2
      write(200+id,*)'j = ', j
      DO k=zb1,zb2
      !DO i=0,nx+1   
         write(200+id,124) (M0(i,j,k),i=0,nx+1)
      !ENDDO
      ENDDO
    ENDDO

    ! Exchange data slices
    call send_receive_slices_integer_y(M0, 0,nx+1,yb1,yb2,zb1,zb2,2)
      
    write(6,*)'lProc',lProc, ' writing exchanged data to the disk ...'
    write(300+id,*)'VARIABLES="X","Y","Z","A0"'
    write(300+id,*)'ZONE F=POINT, I=',nx+2,', J=',yb2-yb1+1,' K=',zb2-zb1+1
    DO j=yb1,yb2
      write(300+id,*)'j = ', j
      DO k=zb1,zb2
      !DO i=0,nx+1   
        !write(200+id,123)xc(i),yc(j),zc(k),A0(i,j,k)
         write(300+id,124) (M0(i,j,k),i=0,nx+1)
      !ENDDO
      ENDDO
    ENDDO

124  format(20(2x,I7))
    
    close(200+id)
    close(300+id)
    
    call MPI_Finalize(ierr)
    stop

  END SUBROUTINE test_send_receive
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------    
  SUBROUTINE write_subdomain() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE
   
    INTEGER                   :: i,j,k,iBody,n,m
    CHARACTER*13              :: fname1

    PRINT*,'Writing out subdomain'
    IF ( lProc >= 0 .AND.  lProc .le. 9  )           &
        write(fname1,341) lProc 
    IF ( lProc >= 10 .AND. lProc .le. 99 )           &
        write(fname1,342) lProc
    IF ( lProc >= 100 .AND. lProc .le. 999 )         &
         write(fname1,343) lProc

341     format('q00',i1,'.dat')
342     format('q0', i2,'.dat')
343     format('q',  i3,'.dat')

    OPEN(UNIT=70,FILE=fname1,STATUS='UNKNOWN')
    
      write(70,*)'VARIABLES="X","Y","Z","U","V","W","P","xb","yb","zb","IBLANK","GHOST","DEAD","BodyNum"'
      write(70,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,' K=',kSlices+2
      do k=zc_start-1,zc_end+1  !1, nz-1
      do j=yc_start-1,yc_end+1
      do i=1,nx-1
         write(70,123)xc(i),yc(j),zc(k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k) &
                     ,xBItable(i,j,k),yBItable(i,j,k),zBItable(i,j,k) &
                     ,iblank(i,j,k),ghostcellmark(i,j,k),dead_cell(i,j,k),bodyNum(i,j,k)
      enddo
      enddo
      enddo

      close(70)
123   format(10(2x,e14.7),4(2x,i2))
      call sleep(5)
  END SUBROUTINE write_subdomain
!***********************************************************************
  SUBROUTINE send_receive_slices_integer_y(f,i1,i2,j1,j2,k1,k2,ks)

!  Send and received ks slices of data between two adjacent subdomains.
!  Each slice is of size [i1:i2, j1:j2] -- note that this has to be the entire
!  slice since
!  the data are contiguous.  The subdomain includes 2 ghost slices on each side.   
!  Does not wait for the receives to finish (call MPI_wait to
!  achieve this).

    USE flow_parameters
    USE MPI_module

    USE boundary_arrays
  
    IMPLICIT NONE
    
    INTEGER,INTENT (IN) :: i1,i2,j1,j2,k1,k2,ks
    INTEGER(KIND=INT_K), DIMENSION(i1:i2,j1:j2,k1:k2),INTENT (IN OUT) :: f

    INTEGER :: ndatay, l1, l2
    INTEGER :: i,j,k
   
    !Data needs to be continuous in memory for MPI passing.
    !Use a temporary array to hold the matrix transpose; only ghost slices are
    !needed.
    INTEGER(KIND=INT_K), DIMENSION(i1:i2,k1:k2,0:ks-1) :: g_tr1, g_tr2
    INTEGER(KIND=INT_K), DIMENSION(i1:i2,k1:k2,0:ks-1) :: g_tr3, g_tr4

    !left buffer transpose
    if (jProc > lProc_leftmost)then
    l1 = j1+2
    l2 = l1+ks-1
    do i=i1,i2
    do j=l1,l2
    do k=k1,k2
       g_tr1(i,k,j-l1) = f(i,j,k)
    enddo
    enddo
    enddo
    endif

    !right buffer transpose
    if (jProc < lProc_rightmost)then
    l2 = j2-2
    l1 = l2-ks+1
    do i=i1,i2
    do j=l1,l2
    do k=k1,k2
       g_tr2(i,k,j-l1) = f(i,j,k)
    enddo
    enddo
    enddo
    endif

    ! number of data to be exchanged, where ks is either 1 or 2.
    ndatay = (i2-i1+1)*ks*(k2-k1+1)
    !print*, 'ndata = ',ndata

    !Send to the left; provide the 1st element location and number of data
    !if (jProc > 0)  &
    if (jProc > lProc_leftmost) &
         call MPI_ISEND(g_tr1(i1,k1,0), ndatay, MPI_INTEGER1, lProc_left, &
                        toleft, FLOW_COMM, isend_rq_toleft, ierr)
    if (ierr /= 0) print*,'Warning: Isend toleft says', ierr

    !Send to the right
    if (jProc < lProc_rightmost)  &
         call MPI_ISEND(g_tr2(i1,k1,0), ndatay, MPI_INTEGER1, lProc_right, &
                        torght, FLOW_COMM, isend_rq_torght, ierr)
    if (ierr /= 0) print*,'Warning: Isend torght says', ierr

    !Receiving from right; provide the 1st element location and number of data
    if (jProc < lProc_rightmost)  &
         call MPI_IRECV(g_tr3(i1,k1,0), ndatay, MPI_INTEGER1, lProc_right, &
                        toleft, FLOW_COMM, irecv_rq_fromrght, ierr)
    if (ierr /= 0) print*,'Warning: Irecv toleft (fromrght) says', ierr

    !Receiving from left
    if (jProc > lProc_leftmost) &
         call MPI_IRECV(g_tr4(i1,k1,0), ndatay, MPI_INTEGER1, lProc_left, &
                        torght, FLOW_COMM, irecv_rq_fromleft, ierr)
    if (ierr /= 0) print*,'Warning: Irecv torght (fromleft) says', ierr
  
    !Wait for sending and receiving to complete
    !if (jProc > 0) &
    if (jProc > lProc_leftmost) &
         call MPI_WAIT(isend_rq_toleft, isend_stat_tl, ierr)
    if (ierr /= 0) print*,'Warning: Wait toleft says', ierr

    !if (jProc < nProcY-1)  &
    if (jProc < lProc_rightmost)  &
         call MPI_WAIT(isend_rq_torght, isend_stat_tr, ierr)
    if (ierr /= 0) print*,'Warning: Wait torght says', ierr

    !if (jProc < nProcY-1)  &
    if (jProc < lProc_rightmost)  &
         call MPI_WAIT(irecv_rq_fromrght, irecv_stat_fr, ierr)
    if (ierr /= 0) print*,'Warning: Wait toleft (fromrght) says', ierr
    
    !if (jProc > 0) &
    if (jProc > lProc_leftmost) &
         call MPI_WAIT(irecv_rq_fromleft, irecv_stat_fl, ierr)
    if (ierr /= 0) print*,'Warning: Wait torght (fromleft) says', ierr

    !transfer left buffer back
    if (jProc > lProc_leftmost)then
    do i=i1,i2
    do k=k1,k2
    do j=0,ks-1
       f(i,j1+2-ks+j,k) = g_tr4(i,k,j)
    enddo
    enddo
    enddo
    endif

    !transfer right buffer back
    if (jProc < lProc_rightmost)then
    do i=i1,i2
    do k=k1,k2
    do j=0,ks-1
      f(i,j2-1+j,k) = g_tr3(i,k,j)
    enddo
    enddo
    enddo
    endif

  END SUBROUTINE send_receive_slices_integer_y
!***********************************************************************
  SUBROUTINE send_receive_slices_integer_z(f,i1,i2,j1,j2,k1,k2,ks)

!  Send and received ks slices of data between two adjacent subdomains.
!  Each slice is of size [i1:i2, j1:j2] -- note that this has to be the entire
!  slice since
!  the data are contiguous.  The subdomain includes 2 ghost slices on each side.   
!  Does not wait for the receives to finish (call MPI_wait to
!  achieve this).

    USE flow_parameters
    USE MPI_module

    USE boundary_arrays
  
    IMPLICIT NONE
    
    INTEGER,INTENT (IN) :: i1,i2,j1,j2,k1,k2,ks
    INTEGER(KIND=INT_K), DIMENSION(i1:i2,j1:j2,k1:k2),INTENT (IN OUT) :: f

    INTEGER :: ndataz, l1, l2
    INTEGER :: i,j,k
   
    ! number of data to be exchanged, where ks is either 1 or 2.
    ndataz = (i2-i1+1)*(j2-j1+1)*ks
    !print*, 'ndata = ',ndata

    !Send to bottom
    l1 = k1+2
    l2 = l1+ks-1
    !if (kProc > 0)  &
    if (kProc > lProc_btommost) &
         call MPI_ISEND(f(i1,j1,l1), ndataz, MPI_INTEGER1, lProc_bottom, &
                        tobtom, FLOW_COMM, isend_rq_tobtom, ierr)
    if (ierr /= 0) print*,'Warning: Isend tobtom says', ierr

    !Send to the top
    l2 = k2-2
    l1 = l2-ks+1
    !if (kProc < nProcZ-1)  &
    if (kProc < lProc_topmost)  &
         call MPI_ISEND(f(i1,j1,l1), ndataz, MPI_INTEGER1, lProc_top, &
                        totop, FLOW_COMM, isend_rq_totop, ierr)
    if (ierr /= 0) print*,'Warning: Isend totop says', ierr

    !Receiving from top
    l1 = k2-1
    l2 = l1+ks-1
    !if (kProc < nProcZ-1)  &
    if (kProc < lProc_topmost)  &
         call MPI_IRECV(f(i1,j1,l1), ndataz, MPI_INTEGER1, lProc_top, &
                        tobtom, FLOW_COMM, irecv_rq_fromtop, ierr)
    if (ierr /= 0) print*,'Warning: Irecv tobtom (fromtop) says', ierr

    !Receiving from bottom
    l2 = k1+1
    l1 = l2-ks+1
    !if (kProc > 0) &
    if (kProc > lProc_btommost) &
         call MPI_IRECV(f(i1,j1,l1), ndataz, MPI_INTEGER1, lProc_bottom, &
                        totop, FLOW_COMM, irecv_rq_frombtom, ierr)
    if (ierr /= 0) print*,'Warning: Irecv totop (frombtom) says', ierr

    !Wait for sending and receiving to complete
    !if (kProc > 0) &
    if (kProc > lProc_btommost) &
         call MPI_WAIT(isend_rq_tobtom, isend_stat_tb, ierr)
    if (ierr /= 0) print*,'Warning: Wait totop says', ierr

    !if (kProc < nProcZ-1)  &
    if (kProc < lProc_topmost)  &
         call MPI_WAIT(isend_rq_totop, isend_stat_tt, ierr)
    if (ierr /= 0) print*,'Warning: Wait tobtom says', ierr

    !if (kProc < nProcZ-1)  &
    if (kProc < lProc_topmost)  &
         call MPI_WAIT(irecv_rq_fromtop, irecv_stat_ft, ierr)
    if (ierr /= 0) print*,'Warning: Wait tobtom (fromtop) says', ierr

    !if (kProc > 0) &
    if (kProc > lProc_btommost) &
         call MPI_WAIT(irecv_rq_frombtom, irecv_stat_fb, ierr)
    if (ierr /= 0) print*,'Warning: Wait totop (frombtom) says', ierr

  END SUBROUTINE send_receive_slices_integer_z
!***********************************************************************
  SUBROUTINE send_receive_slices_real_y(f,i1,i2,j1,j2,k1,k2,ks)

!  Send and received ks slices of data between two adjacent subdomains.
!  Each slice is of size [i1:i2, j1:j2] -- note that this has to be the entire
!  slice since
!  the data are contiguous.  The subdomain includes 2 ghost slices on each side.   
!  Does not wait for the receives to finish (call MPI_wait to
!  achieve this).

    USE flow_parameters
    USE MPI_module
   
    USE boundary_arrays

    IMPLICIT NONE
    
    INTEGER,INTENT (IN) :: i1,i2,j1,j2,k1,k2,ks
    REAL(KIND=CGREAL), DIMENSION(i1:i2,j1:j2,k1:k2),INTENT (IN OUT) :: f

    INTEGER :: ndatay, l1, l2
    INTEGER :: i,j,k

    !Data needs to be continuous in memory for MPI passing.
    !Use a temporary array to hold the matrix transpose; only ghost slices are
    !needed.
    REAL(KIND=CGREAL), DIMENSION(i1:i2,k1:k2,0:ks-1) :: h_tr1, h_tr2
    REAL(KIND=CGREAL), DIMENSION(i1:i2,k1:k2,0:ks-1) :: h_tr3, h_tr4

    !left buffer transpose
    if (jProc > lProc_leftmost)then
    l1 = j1+2
    l2 = l1+ks-1
    do i=i1,i2
    do j=l1,l2
    do k=k1,k2
       h_tr1(i,k,j-l1) = f(i,j,k)
    enddo
    enddo
    enddo
    endif

    !right buffer transpose
    if (jProc < lProc_rightmost)then
    l2 = j2-2
    l1 = l2-ks+1
    do i=i1,i2
    do j=l1,l2
    do k=k1,k2
       h_tr2(i,k,j-l1) = f(i,j,k)
    enddo
    enddo
    enddo
    endif

    ! number of data to be exchanged, where ks is either 1 or 2.
    ndatay = (i2-i1+1)*(k2-k1+1)*ks
    !write(*,'(a,3I5,I12)')'ndata = ',i2,j2,ks,ndata

    !Send to the left; provide the 1st element location and number of data
    if (jProc > lProc_leftmost) &
         call MPI_ISEND(h_tr1(i1,k1,0), ndatay, MPI_DOUBLE_PRECISION, lProc_left, &
                        toleft, FLOW_COMM, isend_rq_toleft, ierr)
    if (ierr /= 0) print*,'Warning: Isend toleft says', ierr

    !Send to the right
    if (jProc < lProc_rightmost)  &
         call MPI_ISEND(h_tr2(i1,k1,0), ndatay, MPI_DOUBLE_PRECISION, lProc_right, &
                        torght, FLOW_COMM, isend_rq_torght, ierr)
    if (ierr /= 0) print*,'Warning: Isend torght says', ierr

    !Receiving from right; provide the 1st element location and number of data
    if (jProc < lProc_rightmost)  &
         call MPI_IRECV(h_tr3(i1,k1,0), ndatay, MPI_DOUBLE_PRECISION, lProc_right, &
                        toleft, FLOW_COMM, irecv_rq_fromrght, ierr)
    if (ierr /= 0) print*,'Warning: Irecv toleft (fromrght) says', ierr

    !Receiving from left
    if (jProc > lProc_leftmost) &
         call MPI_IRECV(h_tr4(i1,k1,0), ndatay, MPI_DOUBLE_PRECISION, lProc_left, &
                        torght, FLOW_COMM, irecv_rq_fromleft, ierr)
    if (ierr /= 0) print*,'Warning: Irecv torght (fromleft) says', ierr
  

    !Wait for sending and receiving to complete
    !if (jProc > 0) &
    if (jProc > lProc_leftmost) &
         call MPI_WAIT(isend_rq_toleft, isend_stat_tl, ierr)
    if (ierr /= 0) print*,'Warning: Wait toleft says', ierr

    !if (jProc < nProcY-1)  &
    if (jProc < lProc_rightmost)  &
         call MPI_WAIT(isend_rq_torght, isend_stat_tr, ierr)
    if (ierr /= 0) print*,'Warning: Wait torght says', ierr

    !if (jProc < nProcY-1)  &
    if (jProc < lProc_rightmost)  &
         call MPI_WAIT(irecv_rq_fromrght, irecv_stat_fr, ierr)
    if (ierr /= 0) print*,'Warning: Wait toleft (fromrght) says', ierr
    
    !if (jProc > 0) &
    if (jProc > lProc_leftmost) &
         call MPI_WAIT(irecv_rq_fromleft, irecv_stat_fl, ierr)
    if (ierr /= 0) print*,'Warning: Wait torght (fromleft) says', ierr

    !transfer left buffer back
    if (jProc > lProc_leftmost)then
    do i=i1,i2
    do k=k1,k2
    do j=0,ks-1
       f(i,j1+2-ks+j,k) = h_tr4(i,k,j)
    enddo
    enddo
    enddo
    endif

    !transfer right buffer back
    if (jProc < lProc_rightmost)then
    do i=i1,i2
    do k=k1,k2
    do j=0,ks-1
      f(i,j2-1+j,k) = h_tr3(i,k,j)
    enddo
    enddo
    enddo
    endif

  END SUBROUTINE send_receive_slices_real_y
!---------------------------------------------------------------------
!***********************************************************************
  SUBROUTINE send_receive_slices_real_z(f,i1,i2,j1,j2,k1,k2,ks)

!  Send and received ks slices of data between two adjacent subdomains.
!  Each slice is of size [i1:i2, j1:j2] -- note that this has to be the entire
!  slice since
!  the data are contiguous.  The subdomain includes 2 ghost slices on each side.   
!  Does not wait for the receives to finish (call MPI_wait to
!  achieve this).

    USE flow_parameters
    USE MPI_module

    USE boundary_arrays

    IMPLICIT NONE
    
    INTEGER,INTENT (IN) :: i1,i2,j1,j2,k1,k2,ks
    REAL(KIND=CGREAL), DIMENSION(i1:i2,j1:j2,k1:k2),INTENT (IN OUT) :: f

    INTEGER :: ndataz, l1, l2
    INTEGER :: i,j,k


    ! number of data to be exchanged, where ks is either 1 or 2.
    ndataz = (i2-i1+1)*(j2-j1+1)*ks
    !write(*,'(a,3I5,I12)')'ndata = ',i2,j2,ks,ndata

    !Send to the bottom; provide the 1st element location and number of data
    l1 = k1+2
    l2 = l1+ks-1    
    !if (kProc > 0)  &
    if (kProc > lProc_btommost) &
         call MPI_ISEND(f(i1,j1,l1), ndataz, MPI_DOUBLE_PRECISION, lProc_bottom, &
                        tobtom, FLOW_COMM, isend_rq_tobtom, ierr)
    if (ierr /= 0) print*,'Warning: Isend tobtom says', ierr

    !Send to the top
    l2 = k2-2
    l1 = l2-ks+1    
    !if (kProc < nProcZ-1)  &
    if (kProc < lProc_topmost)  &
         call MPI_ISEND(f(i1,j1,l1), ndataz, MPI_DOUBLE_PRECISION, lProc_top, &
                        totop, FLOW_COMM, isend_rq_totop, ierr)
    if (ierr /= 0) print*,'Warning: Isend totop says', ierr

    !Receiving from top
    l1 = k2-1
    l2 = l1+ks-1
    !if (kProc < nProcZ-1)  &
    if (kProc < lProc_topmost)  &
         call MPI_IRECV(f(i1,j1,l1), ndataz, MPI_DOUBLE_PRECISION, lProc_top, &
                        tobtom, FLOW_COMM, irecv_rq_fromtop, ierr)
    if (ierr /= 0) print*,'Warning: Irecv tobtom (fromtop) says', ierr

    !Receiving from bottom
    l2 = k1+1
    l1 = l2-ks+1
    !if (kProc > 0) &
    if (kProc > lProc_btommost) &
         call MPI_IRECV(f(i1,j1,l1), ndataz, MPI_DOUBLE_PRECISION, lProc_bottom, &
                        totop, FLOW_COMM, irecv_rq_frombtom, ierr)
    if (ierr /= 0) print*,'Warning: Irecv totop (frombtom) says', ierr
  

    !Wait for sending and receiving to complete
    !if (kProc > 0) &
    if (kProc > lProc_btommost) &
         call MPI_WAIT(isend_rq_tobtom, isend_stat_tb, ierr)
    if (ierr /= 0) print*,'Warning: Wait totop says', ierr

    !if (kProc < nProcZ-1)  &
    if (kProc < lProc_topmost)  &
         call MPI_WAIT(isend_rq_totop, isend_stat_tt, ierr)
    if (ierr /= 0) print*,'Warning: Wait tobtom says', ierr

    !if (kProc < nProcZ-1)  &
    if (kProc < lProc_topmost)  &
         call MPI_WAIT(irecv_rq_fromtop, irecv_stat_ft, ierr)
    if (ierr /= 0) print*,'Warning: Wait tobtom (fromtop) says', ierr
    
    !if (kProc > 0) &
    if (kProc > lProc_btommost) &
         call MPI_WAIT(irecv_rq_frombtom, irecv_stat_fb, ierr)
    if (ierr /= 0) print*,'Warning: Wait totop (frombtom) says', ierr

  END SUBROUTINE send_receive_slices_real_z
!---------------------------------------------------------------------
