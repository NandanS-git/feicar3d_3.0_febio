c
c     Modified for BIG TIME LOOP (Ye, Apr-2015)
c4400
c     NON-liner 3D EXPlicit incremental analysis 
      subroutine non_3D_exp_Tian(wk,pload, disp, fmag,
     &                        dispt , 
     &                   mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt)

          use contact

          implicit real*8 (a-h,o-z)
          include 'commons.std'

         !--------------------------- F.-B. Tian
         include "mpif.h"

#ifdef _OPENMP
         include 'omp_lib.h'
#endif

         INTEGER ::  IERR,IPES,ierror,npes,nrank,proc_m
         integer ::  istatus(MPI_STATUS_SIZE)
         integer ::  icvg
         !---------------------------

         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=45000)
         integer nmap [allocatable] (:)
	 real*8  xyz_s [allocatable] (:)
	 real*8  xyz_f [allocatable] (:),contact_f [allocatable] (:)
	 real*8  xyz_f_temp [allocatable] (:)
         real*8  A0 [allocatable] (:),r0 [allocatable] (:)

         integer nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer inod(ireadmax,2),mnod(iconelem)

         real*8 tforce1(maxforce,2)
         real*8 velout(ireadmax)
         real*8 mass(maxmass),wk(neq)
     &          ,pload(neq ),pload_tmp(neq),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(maxxyz)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 grav(maxxyz)
         real*8 dispt(neq)
         real*8 dispfult(maxnode,3),dispfult_s(maxnode,3)
         real*8 disptm(neq),disptp(neq),dispt_s(neq)
         real*8 disptm_s(neq)
         real*8 relax_FSIp
         integer FSI_on
!=================for vel calculation
!	 real*8 dispt_1(neq),dispt_0(neq),dispt2(neq),dispt1(neq)
!=========F.-B. Tian
         real*8 dispti(neq),dispfulti(maxnode,3),disptmi(neq)
	 real*8 xyzti(maxxyz),xyzt_s(maxxyz)
!===================
c
         real*8 dload(neq),floadt(3*maxnode)

           !Ye, contact
           real*8 contact_fx,contact_fy,contact_fz
           real*8  cnt_fx [allocatable] (:)
           real*8  cnt_fy [allocatable] (:)
           real*8  cnt_fz [allocatable] (:)
           real*8  cnt_force [allocatable] (:)           

           real*8  cnt_fx12 [allocatable] (:)
           real*8  cnt_fx13 [allocatable] (:)
           real*8  cnt_fx21 [allocatable] (:)
           real*8  cnt_fx23 [allocatable] (:)
           real*8  cnt_fx31 [allocatable] (:)
           real*8  cnt_fx32 [allocatable] (:)

           real*8  cnt_fy12 [allocatable] (:)
           real*8  cnt_fy13 [allocatable] (:)
           real*8  cnt_fy21 [allocatable] (:)
           real*8  cnt_fy23 [allocatable] (:)
           real*8  cnt_fy31 [allocatable] (:)
           real*8  cnt_fy32 [allocatable] (:)

           real*8  cnt_fz12 [allocatable] (:)
           real*8  cnt_fz13 [allocatable] (:)
           real*8  cnt_fz21 [allocatable] (:)
           real*8  cnt_fz23 [allocatable] (:)
           real*8  cnt_fz31 [allocatable] (:)
           real*8  cnt_fz32 [allocatable] (:)

           integer  ElemNum, ii
           real*8 cnt_d, csi

           integer clock0, clock1, clock_rate
           integer clock2, clock3,clock4,clock5

           integer :: nthreads,tid
           real*8 rho_f

           real*8 cf_min, cf_max
           real   time0
           INTEGER,DIMENSION(1)     :: cfmin_loc(1)
           INTEGER,DIMENSION(1)     :: cfmax_loc(1)

         maxdim=maxnode
         maxdof=3*maxdim
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         write(ilog,*)'@@  non_3D_exp ver270  November 2007 '

!        algorithm parameters
!          write(*,*) maxstiff,maxnode,neq,' ms'
!          stop
         kref=2
         csi = 1.0E-10
         iter_fsi = 0

!        Get things ready for time integration loop
!         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
!     &                                        ' snap count '
!         call zzwrt('  -->  ')
!         read(ikbd,*) deltat,npt,iprcnt,isnap
         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
     &                       ' snap count | is_res | irestart | isub'
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,iprcnt,isnap, is_res,irestart,isub ! F.-B. Tian
!-------------------dt,    Nstep, nout, nscreen, res, ires===============

         !Ye,contact info input============================
       ! total # of surf tri elems
       totNumElem = 20952
       ! total number of surf pts
       nPtsBM = 10482

       allocate (cnt_d12(nPtsBM))
       allocate (cnt_d21(nPtsBM))
       allocate (cnt_d13(nPtsBM))
       allocate (cnt_d31(nPtsBM))
       allocate (cnt_d23(nPtsBM))
       allocate (cnt_d32(nPtsBM))

       allocate (ElemNum12(nPtsBM))
       allocate (ElemNum21(nPtsBM))
       allocate (ElemNum13(nPtsBM))
       allocate (ElemNum31(nPtsBM))
       allocate (ElemNum23(nPtsBM))
       allocate (ElemNum32(nPtsBM))

       !allocate (cnt_d(nPtsBM))
       allocate (BM_num(nPtsBM))
       allocate (mark(nPtsBM))
       !allocate (ElemNum(nPtsBM))
       allocate (ElemNeig(3,totNumElem))
       allocate (BM_x(nPtsBM),BM_y(nPtsBM),BM_z(nPtsBM))
       allocate (ElemCentx(totNumElem),ElemCenty(totNumElem),
     &           ElemCentz(totNumElem))
       allocate (ElemNorm_x(totNumElem),ElemNorm_y(totNumElem),
     &           ElemNorm_z(totNumElem))

       allocate (xNormBM(nPtsBM),yNormBM(nPtsBM),zNormBM(nPtsBM))

       cnt_d = 1.0E+20

       !Ye,read from formatted input file
       open(unit=1,file='leaflets_file.dat')
       read(1,*)
       !LEAFLET 1,mark=1 for contact detection
       do i=1,nPtsBM/3
          read(1,*)j,BM_num(i),mark(i)!solid node bianhao
       enddo
       do i=1,totNumElem/3
          read(1,*)ElemNeig(1,i),ElemNeig(2,i),ElemNeig(3,i)!xuhao
       enddo
       read(1,*)
       !LEAFLET 2
       do i=nPtsBM/3+1,2*nPtsBM/3
          read(1,*)j,BM_num(i),mark(i)
       enddo
       do i=totNumElem/3+1,2*totNumElem/3
          read(1,*)ElemNeig(1,i),ElemNeig(2,i),ElemNeig(3,i)
       enddo
       read(1,*)
       !LEAFLET 3
       do i=2*nPtsBM/3+1,nPtsBM
          read(1,*)j,BM_num(i),mark(i)
       enddo
       do i=2*totNumElem/3+1,totNumElem
          read(1,*)ElemNeig(1,i),ElemNeig(2,i),ElemNeig(3,i)
       enddo
       close(1)
       write(*,*)'finish reading leaflet_file.dat!!'

       pointOutsideX = 0.0d0
       pointOutsideY = 0.0d0
       pointOutsideZ = 0.0d0

       !Ye,READ 'COMPACT_CONTACT.DAT'
       !STORE CONTACT XUHAO
       !cnt_n = 2490
       allocate (cnt_xh(cnt_n))
       open(unit=1,file='compact_contact.dat')
       do i=1,cnt_n
          read(1,*)cnt_xh(i)
       enddo
       close(1)

       allocate (cnt_fx12(cnt_n/3),cnt_fy12(cnt_n/3)
     &          ,cnt_fz12(cnt_n/3))
       allocate (cnt_fx13(cnt_n/3),cnt_fy13(cnt_n/3)
     &          ,cnt_fz13(cnt_n/3))
       allocate (cnt_fx21(cnt_n/3),cnt_fy21(cnt_n/3)
     &          ,cnt_fz21(cnt_n/3))
       allocate (cnt_fx23(cnt_n/3),cnt_fy23(cnt_n/3)
     &          ,cnt_fz23(cnt_n/3))
       allocate (cnt_fx31(cnt_n/3),cnt_fy31(cnt_n/3)
     &          ,cnt_fz31(cnt_n/3))
       allocate (cnt_fx32(cnt_n/3),cnt_fy32(cnt_n/3)
     &          ,cnt_fz32(cnt_n/3))

       write(*,*)'finish reading compact_contact.dat!!'
       !contact===========================================

       !store triElem #s containing each surface point
       ALLOCATE (ElemC(nPtsBM,8))
       ElemC = 0
       do i=1,nPtsBM
          ii=1
          do m=1,totNumElem
             IF ( (ElemNeig(1,m).eq.i).or.(ElemNeig(2,m).eq.i).or.
     &            (ElemNeig(3,m).eq.i) ) THEN
                ElemC(i,ii)=m
                ii = ii + 1
             ENDIF
          enddo
          IF (ii.gt.9) THEN
             WRITE(*,*)'ERROR in storing TriElems!'
          ENDIF
        enddo

       open(unit=1,file='nmap.dat')
       read(1,*)nshare
       allocate(nmap(nshare))
       allocate(xyz_s(3*nshare))
       allocate(xyz_f(3*nshare))
       allocate(contact_f(nshare))
       allocate(xyz_f_temp(3*nshare))
       allocate(xyz_f_old(3*nshare))
       allocate(A0(nshare))
       allocate(r0(3*nshare))
       allocate(cnt_fx(nshare))
       allocate(cnt_fy(nshare))
       allocate(cnt_fz(nshare))
       allocate(cnt_force(nshare))
       xyz_f = 0.0d0
       xyz_f_temp = 0.0d0
       xyz_f_old = 0.0d0
       do i=1,nshare
          read(1,*)j,nmap(i),A0(i),r0(i)
       enddo
       close(1)

         if(is_res.eq.1)then
         call nond_3D_exp_restart_read(time,itime0,tf1,tf1s,k1old, 
     & maxnode, maxxyz, xyzt,dispt,disptm, dispfult, vel, acc)
!      call restart_read
!xyzt(),gsize,itime0-->itime, mass,dispt,vel,acc,disptm,time,k1old,tf1s,tf1
         else
          xn=0.0
          yn=0.0
          zn=0.0
          do i=1,nnp
              xyzt(i)=xyz0(i)
              xyzt(maxnode+i)=xyz0(maxnode+i)
              xyzt(maxnode*2+i)=xyz0(maxnode*2+i)
              xn=xn+(xyzt(          i)-xyzt(          1))**2
              yn=yn+(xyzt(maxnode  +i)-xyzt(maxnode  +1))**2
              zn=zn+(xyzt(maxnode*2+i)-xyzt(maxnode*2+1))**2
          enddo
          gsize=( (xn+yn+zn)/neq )**0.5
          write(ilog,*)'@@ gsize: ',gsize
          time = 0.0
          itime0=1
         endif

       do i=1,nshare
          xyz_s(i)=xyzt(nmap(i))
          xyz_s(nshare+i)=xyzt(maxnode+nmap(i))
          xyz_s(2*nshare+i)=xyzt(2*maxnode+nmap(i))
       enddo

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         dt=deltat
!==========================================================================
! Communicate with the flow solver

         irankf = 1     ! the processor from the flow side

         call MPI_SEND(deltat, 1, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(is_res, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(irestart, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(npt, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)
 
         !send the initial share coordinates
         call MPI_SEND(nshare, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)  
         call MPI_SEND(xyz_s(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare+1  ), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare*2+1), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)
!==========================================================================

         if (isnap .eq. 0) isnap=10000

         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'

         nsnapmax=npt
         write(*,*)' '
         write(*,*)'CHOOSE nodal output:'
         nout=0
27       continue
         nout=nout+1
         write(*,*)'TYPE:  node# |  DoF  | rate     <0 0 0 to end>'
         call zzwrt('  -->  ')
         read(ikbd,*) inode,idof,irate
         write(ilog,*) inode,idof,irate,' ::node  dof '
         if (inode .ne. 0) then
             inod(nout,1)=idbc(inode,idof)
             inod(nout,2)=irate
             write(ilog,*) '@@ eqn # ',inod(nout,1)
             goto 27
         endif
         nout=nout-1
         write(*,'(a)')' @@ '

!           form lumped mass
            ilump=1
!           write(*,*) maxstiff,maxnode,neq,' ms'
            write(ilog,*)'@@ FORMing lumped mass matrix '
!           call zztimer(ilog,'--> FORMmass')
            call formmass_Tian(mass(1), npijkm,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      idbc, maxnode,maxelem,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp)
            z1=0.0
            do i=1,neq
               z1=z1+mass(i)
            enddo
            write(iout,*)'Total Mass: ',z1/3
!           write(iout,86) (mass(j),j=1,neq)
!           call zztimer(ilog,'<-- FORMmass')
!           call zzflush(ilog)

            if (igravity_on .eq. 1) then
                call body_force_grav( iglobal,
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                   prop,nmprop,ippp,
     &                   wk(1)
     &                          )
                z1=0.0
                do i=1,neq
                   grav(i)=wk(i)
                   z1=z1+grav(i)
                enddo
                write(iout,*)'Total gravity force: ',z1
!               write(iout,86) (grav(j),j=1,neq)
            else
                do i=1,neq
                   grav(i)=0.0  
                enddo
            endif
!           stop
         rewind(ilod)
         do i=1,neq
            read(ilod) fmag(i)
         enddo
         write(*   ,*) '@@ Reloaded  {P}   OK'
         write(ilog,*) '@@ Reloaded  {P}   OK'

!        INPUT load history from a file and interpolate 
         call getload_hist(tforce1,maxforce,npt,deltat)
         write(*,*)'                    {P}  | {Grav} |     |     '
         write(*,*)'TYPE force scales:   P1  |  P2    |  P3 | P4  '
         call zzwrt('  -->  ')
         read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'

!        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
         write(ilog ,*)'@@  Beginning increment analysis   '

!       Initialize  times, disp
        if(is_res.eq.1)then


        else
           write(*  ,*) '@@ INITIAL disps set to zero'
           do i= 1, neq  
              dispt(i)=0.0
              vel(i)=0.0
              acc(i)=-dmm*vel(i)
              disptm(i) = dispt(i)-dt*vel(i)+0.5*dt*dt*acc(i)
           enddo
              call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))
        endif

           do 71 n=1,nout
              iprnode=inod(n,1)
              write(ilog,*)'@@ iprnode ',n,iprnode
              if (iprnode .eq. 0) then
                  velout(n) = 0.0
                  goto 71
              endif
              if (inod(n,2) .eq. 0) then
                  velout(n)=dispt(iprnode)
              endif
 71        continue
           rewind(idyn)

        if(is_res.eq.1)then
        istart0=itime0+1
        iend0=itime0+npt-1
        else
           k1old=1
           tf1s=tforce1(1,2)*scale_p1
           tf1=tforce1(1,2)
           write(idyn,182) time,tf1,(velout(n), n=1,nout ),tf1

           rewind(isnp)
           write(isnp) time, neq, nsnapmax 
           do i= 1, neq  
              write(isnp) dispt(i)
           enddo

           istart0=2
           iend0=npt
        endif

         kount = 0
         ksnap=1
         ksnap=0
         inc=1
         inc=iprcnt
         icarriage=0

         call zzwrt(' @@ ')

         call write_coord(xyzt,maxelem,                ! Fangbao Tian
     &                     npijkm,maxnode,maxxyz,0)


!        BIG TIME LOOP
!============================
         ramp=2.0

         do 200 itime= istart0,iend0
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
            dt0=deltat/isub
            
            !save the current itime results
            xyzt_s    =xyzt
            dispfult_s=dispfult
            disptm_s  =disptm
            dispt_s   =dispt

2013        continue

            time = time0 + real(itime-itime0-1)*deltat
            ksub=0 

!------------------------------------------------------
! Communicate with the flow solver to receive the load
! for each node.  

        if(ksub .eq. 0) then

          xyz_f_temp = xyz_f  ! dai
          xyz_f = 0.0

          call MPI_SEND(itime-1, 1, MPI_INTEGER, irankf, 1,  
     &                 MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(relax_FSIp, 1, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(FSI_on, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(iter_fsi, 1, MPI_INTEGER, irankf, 1,
     &                 MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(rho_f, 1, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(xyz_f(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(xyz_f(nshare+1  ), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(xyz_f(nshare*2+1), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          resMax_FSIp = 0.0d0
          do i=1,nshare*3
            resP = xyz_f_temp(i)-xyz_f(i)
            xyz_f(i) = xyz_f(i)+ resP*relax_FSIp
            resMax_FSIp  = max(resMax_FSIp,   abs(resP))
          enddo

          call MPI_SEND(resMax_fsip, 1, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          ! initialize the state for substepping
          xyzt    =xyzt_s
          dispfult=dispfult_s
          disptm  =disptm_s
          dispt   =dispt_s

         endif!(ksub .eq. 0)

130        continue ! return point for substep
            ksub=ksub+1
            time = time0 + real(itime-itime0-1)*deltat+dt0*ksub
            pt1=1.0/(2.0*dt0)
            pt2=1.0/dt0**2
!            call zzcount_2(itime,inc,icarriage)

!           set load increment
!           Fmag says where the load is applied
!           Done this way in case distributed load applied
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            do i= 1, neq  
               pload(i) =  tf1*(fmag(i)*scale_p1 + grav(i)*scale_p2) 
               pload_tmp(i) = pload(i)
            enddo 
!============================
      ! initialize the current substep
      xyzti=xyzt
      dispfulti=dispfult
      disptmi=disptm
      dispti=dispt

      !Ye,use the current position to calc triElem norm===
      do i=1,nPtsBM!xuhao
         BM_x(i)=xyzt(BM_num(i))
         BM_y(i)=xyzt(maxnode+BM_num(i))
         BM_z(i)=xyzt(2*maxnode+BM_num(i))
      enddo

      call calculate_arclength_norm_cnt(itime,istart0,ksub)
!=============================
        !reassign the distributed force to nodal load
       sc  = 1.0d0 - exp(-ramp*time) ! load ramp factor

       !Ye, contact force initialization
       do i=1,nPtsBM
          cnt_fx(i)=0.0d0
          cnt_fy(i)=0.0d0
          cnt_fz(i)=0.0d0
       enddo

       do i=1,nPtsBM/3
          cnt_fx12(i)=0.0d0
          cnt_fy12(i)=0.0d0
          cnt_fz12(i)=0.0d0
          cnt_fx13(i)=0.0d0
          cnt_fy13(i)=0.0d0
          cnt_fz13(i)=0.0d0
          cnt_fx21(i)=0.0d0
          cnt_fy21(i)=0.0d0
          cnt_fz21(i)=0.0d0
          cnt_fx23(i)=0.0d0
          cnt_fy23(i)=0.0d0
          cnt_fz23(i)=0.0d0
          cnt_fx31(i)=0.0d0
          cnt_fy31(i)=0.0d0
          cnt_fz31(i)=0.0d0
          cnt_fx32(i)=0.0d0
          cnt_fy32(i)=0.0d0
          cnt_fz32(i)=0.0d0
       enddo

       call system_clock(clock4, clock_rate)

       IF (ksub.eq.1) THEN
       call system_clock(clock2, clock_rate)
       ENDIF
       !write(*,*)'Nonstad: entering contact force.'

#ifdef _OPENMP
       call    OMP_SET_DYNAMIC(.FALSE.)
       call    OMP_SET_NUM_THREADS(8)
#endif

!!$OMP PARALLEL PRIVATE(nthreads, tid)
!!Obtain thread number
!      IF (ksub.eq.1) THEN
!         tid = OMP_GET_THREAD_NUM()
!         write(*,*) 'Hello from sub: tid=', tid

!!Only master thread does this
!         IF (tid .EQ. 0) THEN
!            nthreads = OMP_GET_NUM_THREADS()
!            write(*,*)'Number of threads = ', nthreads
!         ENDIF
!      ENDIF
!!$OMP END PARALLEL


!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum,tid)
!$OMP DO
            !Ye,Contact Btw LEAFLET_1 & LEAFLET_2
            !calculate cnt_d(i)
            do i=1,cnt_n/3
               ii = cnt_xh(i)
               !tid = OMP_GET_THREAD_NUM()
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+20

               !IF (ksub.eq.1) THEN
               !   write(*,'(a,I5,a,I5)'),'i=', i, ' by tid ', tid
               !ENDIF

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii
     &         ,nPtsBM/3+1,2*nPtsBM/3,cnt_d,ElemNum)

               cnt_d12(ii) = cnt_d
               ElemNum12(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               !calculate total external force at node i
               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -ElemNorm_x(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -ElemNorm_y(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -ElemNorm_z(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = ElemNorm_x(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = ElemNorm_y(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = ElemNorm_z(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,interp the contact force
               cnt_fx12(i)=contact_fx*tf1
               cnt_fy12(i)=contact_fy*tf1
               cnt_fz12(i)=contact_fz*tf1
               !=====================================
               !IF (ksub.eq.1) THEN
               !   write(*,'(a,I5,a,I5)'),'i=', i, ' by tid ', tid
               !ENDIF
            enddo
!$OMP END PARALLEL

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum)
!$OMP DO
           do i=cnt_n/3+1,2*cnt_n/3
               ii = cnt_xh(i)
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+20

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii
     &              ,1,nPtsBM/3,cnt_d,ElemNum)

               cnt_d21(ii) = cnt_d
               ElemNum21(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -ElemNorm_x(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -ElemNorm_y(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -ElemNorm_z(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = ElemNorm_x(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = ElemNorm_y(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = ElemNorm_z(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,debug
!               if ((itime.eq.101).and.(ii.eq.5443)!.and.(iter_fsi.eq.1)
!     &             .and.(ksub.eq.1))then
!                  write(502,'(4F12.7)')xyz_f(ii),xyz_f(nshare+ii),
!     &                                 xyz_f(nshare*2+ii),f_ext
!                  write(502,'(3F12.7)')ElemNorm_x(ElemNum),
!     &                                 ElemNorm_y(ElemNum),
!     &                                 ElemNorm_z(ElemNum)
!                  write(502,*)ElemNum,cnt_d,cnt_c,cnt_h
!                  write(502,'(3F12.7)')contact_fz,contact_fy,contact_fz
!               endif

               !Ye,interp the contact force
               cnt_fx21(i-nPtsBM/3)=contact_fx*tf1
               cnt_fy21(i-nPtsBM/3)=contact_fy*tf1
               cnt_fz21(i-nPtsBM/3)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum)
!$OMP DO
            !Ye,Contact Btw LEAFLET_1 & LEAFLET_3
            !calculate cnt_d(i)
            do i=1,cnt_n/3
               ii = cnt_xh(i)
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+20

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii
     &         ,2*nPtsBM/3+1,nPtsBM,cnt_d,ElemNum)

               cnt_d13(ii) = cnt_d
               ElemNum13(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -ElemNorm_x(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -ElemNorm_y(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -ElemNorm_z(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = ElemNorm_x(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = ElemNorm_y(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = ElemNorm_z(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,interp the contact force
               cnt_fx13(i)=contact_fx*tf1
               cnt_fy13(i)=contact_fy*tf1
               cnt_fz13(i)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum)
!$OMP DO
           do i=2*cnt_n/3+1,cnt_n
               ii = cnt_xh(i)
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+20

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii
     &              ,1,nPtsBM/3,cnt_d,ElemNum)

               cnt_d31(ii) = cnt_d
               ElemNum31(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -ElemNorm_x(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -ElemNorm_y(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -ElemNorm_z(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = ElemNorm_x(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = ElemNorm_y(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = ElemNorm_z(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,interp the contact force
               cnt_fx31(i-2*nPtsBM/3)=contact_fx*tf1
               cnt_fy31(i-2*nPtsBM/3)=contact_fy*tf1
               cnt_fz31(i-2*nPtsBM/3)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum)
!$OMP DO
            !Ye,Contact Btw LEAFLET_3 & LEAFLET_2
            !calculate cnt_d(i)
           do i=2*cnt_n/3+1,cnt_n
               ii = cnt_xh(i)
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+20

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii
     &         ,nPtsBM/3+1,2*nPtsBM/3,cnt_d,ElemNum)

               cnt_d32(ii) = cnt_d
               ElemNum32(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -ElemNorm_x(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -ElemNorm_y(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -ElemNorm_z(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = ElemNorm_x(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = ElemNorm_y(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = ElemNorm_z(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,interp the contact force
               cnt_fx32(i-2*nPtsBM/3)=contact_fx*tf1
               cnt_fy32(i-2*nPtsBM/3)=contact_fy*tf1
               cnt_fz32(i-2*nPtsBM/3)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum)
!$OMP DO
            do i=cnt_n/3+1,2*cnt_n/3
               ii = cnt_xh(i)
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+20

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii
     &         ,2*nPtsBM/3+1,nPtsBM,cnt_d,ElemNum)

               cnt_d23(ii) = cnt_d
               ElemNum23(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -ElemNorm_x(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -ElemNorm_y(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -ElemNorm_z(ElemNum)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = ElemNorm_x(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = ElemNorm_y(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = ElemNorm_z(ElemNum)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,interp the contact force
               cnt_fx23(i-nPtsBM/3)=contact_fx*tf1
               cnt_fy23(i-nPtsBM/3)=contact_fy*tf1
               cnt_fz23(i-nPtsBM/3)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL
            !=================================================
            !write(*,*)'Nonstad:contact force calculated !!!'
            IF(ksub.eq.isub)THEN
            call system_clock(clock3, clock_rate)
            WRITE(*,*) ' Total time for contact force is:',
     &            REAL(clock3-clock2)/REAL(clock_rate)
            ENDIF

!            IF(ksub.le.isub)THEN
!            call system_clock(clock5, clock_rate)
!            WRITE(*,*) ' Sub time for contact force is:',
!     &            REAL(clock5-clock4)/REAL(clock_rate)
!            ENDIF

            !Ye,tansform the cusp-wise contact force
            do i=1,cnt_n/3
               ii = cnt_xh(i)
               cnt_fx(ii)=cnt_fx12(i)+cnt_fx13(i)
               cnt_fy(ii)=cnt_fy12(i)+cnt_fy13(i)
               cnt_fz(ii)=cnt_fz12(i)+cnt_fz13(i)
            enddo

            do i=cnt_n/3+1,2*cnt_n/3
               ii = cnt_xh(i)
               cnt_fx(ii)=cnt_fx21(i-cnt_n/3)+cnt_fx23(i-cnt_n/3)
               cnt_fy(ii)=cnt_fy21(i-cnt_n/3)+cnt_fy23(i-cnt_n/3)
               cnt_fz(ii)=cnt_fz21(i-cnt_n/3)+cnt_fz23(i-cnt_n/3)
            enddo

            do i=2*cnt_n/3+1,cnt_n
               ii = cnt_xh(i)
               cnt_fx(ii)=cnt_fx31(i-2*cnt_n/3)+cnt_fx32(i-2*cnt_n/3)
               cnt_fy(ii)=cnt_fy31(i-2*cnt_n/3)+cnt_fy32(i-2*cnt_n/3)
               cnt_fz(ii)=cnt_fz31(i-2*cnt_n/3)+cnt_fz32(i-2*cnt_n/3)
            enddo
            !=======================================

       IF (ksub.eq.isub) THEN
       do i=1,nPtsBM
          cnt_force(i)=sqrt(cnt_fx(i)**2+cnt_fy(i)**2+cnt_fz(i)**2)
       enddo

       cf_min=minval(cnt_force(1:nPtsBM))
       cf_max=maxval(cnt_force(1:nPtsBM))
       cfmin_loc=minloc(cnt_force(1:nPtsBM))
       cfmax_loc=maxloc(cnt_force(1:nPtsBM))

       write(*,'(a,2E12.4)')'Min-Max of CntForce:',cf_min,cf_max
       write(*,*)'Min-Max of CntForce at:',cfmin_loc(1),cfmax_loc(1)
       ENDIF!ksub.eq.isub

102    CONTINUE

       do i=1,nshare
          ii=nmap(i)
          ij=idbc(ii,1)
          if (ij.ne.0) then
          pload(ij)=pload_tmp(ij)+xyz_f(i)*tf1+cnt_fx(i)
          endif
          ij=idbc(ii,2)
          if (ij.ne.0) then
          pload(ij)=pload_tmp(ij)+xyz_f(nshare+i)*tf1+cnt_fy(i)
          endif
          ij=idbc(ii,3)
          if (ij.ne.0) then
          pload(ij)=pload_tmp(ij)+xyz_f(nshare*2+i)*tf1+cnt_fz(i)
          endif
       enddo
!============================


!============================
!              forces from body stresses
               call body_stress_rub_tian(
     &                   dispfulti,
     &                   iglobal,
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                   prop,nmprop,ippp,
     &                   floadt,maxdim,maxdof
     &                )

              do i= 1, neq  
                 aterm =  (pt2)*(2*dispti(i)-disptmi(i))
                 vterm =  (pt1)*disptmi(i)
                 dload(i) = pload(i) - floadt(i) 
     &                      + mass(i)*aterm + dmm*mass(i)*vterm 
                 coeff = pt2 + 1.0*dmm*pt1
                 disptp(i) = dload(i)/( coeff*mass(i))
              enddo
              do i= 1, neq
                 acc(i) = (disptp(i)-2*dispti(i)+disptmi(i))*pt2
                 vel(i) = (disptp(i)-disptmi(i))*pt1
              enddo
              do i= 1, neq
                 disptmi(i) = dispti(i)
                 dispti(i)  = disptp(i)
              enddo
             call update_geom_3D(disptp, dispfulti,maxdim, idbc,maxnode,
     &                     xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                     xyzti(1),xyzti(maxnode+1),xyzti(maxnode*2+1))
 

!================TEST FSI LOOP=================
                do 1722 n=1,nout
                   iprnode=inod(n,1)
                   if (inod(n,2) .eq. 0) then
                       velout(n)=dispti(iprnode)
                   elseif (inod(n,2) .eq. 1) then
                       velout(n)=vel(iprnode)
                   elseif (inod(n,2) .eq. 2) then
                       velout(n)=acc(iprnode)
                   endif
 1722            continue
                write(1111,'(15F15.8)') time,(velout(n), n=1,nout)
!===============================================

        xyzt    =xyzti
        dispfult=dispfulti
        disptm  =disptmi
        dispt   =dispti

      if(ksub.lt.isub)then
         go to 130

      else  ! ksub = isub

         print*, 'Nonstad sending ... '

         ! pack and send the node position
          do i=1,nshare
            xyz_s(i)=xyzti(nmap(i))
            xyz_s(nshare+i)=xyzti(maxnode+nmap(i))
            xyz_s(2*nshare+i)=xyzti(2*maxnode+nmap(i))
          enddo
         
         call MPI_SEND(xyz_s(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare+1  ), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare*2+1), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          ! pack the velocity data and send to the flow solver
             do i=1,nshare
               ii=nmap(i)
               ij=idbc(ii,1)
               xyz_s(i)=vel(ij)
               ij=idbc(ii,2)
               xyz_s(nshare+i)=vel(ij)
               ij=idbc(ii,3)
               xyz_s(2*nshare+i)=vel(ij)
             enddo

         call MPI_SEND(xyz_s(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare+1  ), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare*2+1), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         !write(*,666) 'xvel in solv: ',xyz_v(1:10)
         !write(*,666) 'yvel in solv: ',xyz_v(maxnode+1: maxnode+10)
         !write(*,666) 'zvel in solv: ',xyz_v(maxnode*2+1:maxnode*2+10)

         ! Receive signal from the flow solver about the FSI convergence.
         call MPI_RECV(icvg, 1, MPI_INTEGER, irankf, 1,
     &                 MPI_COMM_WORLD,istatus,ierr)

         ! If not converged, go to back to load update
         if(icvg .eq. 0) then
            goto 2013
         else
            write(*,*) 'NONSTAD: FSI iteration done.'
         endif

      endif  ! end of (ksub .lt. isub)
!=====================================

            !====================================
            !Ye,send contact force to flow slover
            CALL MPI_SEND(nshare, 1, MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_fx, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_fy, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_fz, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_d12, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_d21, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_d13, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_d31, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_d23, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_d32, nshare,
     &           MPI_DOUBLE_PRECISION,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(ElemNum12, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(ElemNum21, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(ElemNum13, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(ElemNum31, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(ElemNum23, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(ElemNum32, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            write(*,*)'Nonstad:contact force sent !'

!          Print out results of interest at each time step
!           do 174 n=1,nout
!              iprnode=inod(n,1)
!              irate  =inod(n,2)
!              if (iprnode .eq. 0) then
!                 velout(n) = 0.0
!                 goto 174
!              endif
!              if (inod(n,2) .eq. 0) then
!                 velout(n)=dispt(iprnode)
!              elseif (inod(n,2) .eq. 1) then
!                 velout(n)=vel(iprnode)
!              elseif (inod(n,2) .eq. 2) then
!                 velout(n)=acc(iprnode)
!              endif
! 174          continue
!              tf1s=tf1*scale_p1
!              write(idyn,182) time,tf1s,(velout(n), n=1,nout ),tf1


!           Print out results of interest for this time loop
            if (iprcnt .eq. kount) then
                kount = 0
                nout1=nout-1
                do 172 n=1,nout
                   iprnode=inod(n,1)
                   irate  =inod(n,2)
                   if (iprnode .eq. 0) then
                       velout(n) = 0.0
                       goto 172
                   endif
                   if (inod(n,2) .eq. 0) then
                       velout(n)=dispt(iprnode)
                   elseif (inod(n,2) .eq. 1) then
                       velout(n)=vel(iprnode)
                   elseif (inod(n,2) .eq. 2) then
                       velout(n)=acc(iprnode)
                   endif
 172            continue
                tf1s=tf1*scale_p1
!               write(idyn,182) time,tf1s,(velout(n), n=1,nout ),tf1
!                write(*,*)'Itime=',itime

!               call write_coordHEX20(xyzt,                ! Fangbao Tian
!     &                     npijkm,maxnode,maxelem,maxxyz,itime-1)!
                call write_coord(xyzt,maxelem,                ! Fangbao Tian
     &                     npijkm,maxnode,maxxyz,itime-1)
             endif
             if (isnap .eq. ksnap) then
                 write(isnp) time, neq, nsnapmax 
                 ksnap=0
             endif

!  Fang-Bao Tian
! Writing the restart file
       if(mod(itime-1,irestart) .eq. 0) then
        call nond_3D_exp_restart_write(time,itime,tf1,tf1s,k1old,
     &  maxnode, maxxyz, xyzt,dispt,disptm, dispfult, vel,acc)
       endif
 200  continue
 210  continue
      write(*,'(a)')' @@'
c     END BIG TIME LOOP
c
 81        format(1x,a,1x,90(1x,g12.6))
 82        format(1x,90(1x,g12.6))
 83        format(1x,4(1x,g12.6),1x,a)
 84        format(1x,a,1x,2(1x,g12.6),1x,2(i5,1x))
 85        format(1x,90(1x,g12.6))
 86        format(1x, 6(1x,g12.6))
 182       format(1x,20(g12.6,1x),/,20(g12.6,1x),/,20(g12.6,1x))
 183       format(1x,i5,4x,3(g12.6,1x))

      return
      end


      subroutine body_stress_rub_tian(
     &                   dispful,
     &                   iglobal,
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xord0,yord0,zord0,
     &                   prop,nmprop,ippp,
     &                   gforce,maxdim,maxdof
     &                )

          implicit real*8 (a-h,o-z)
             include 'commons.std'

         integer npijkm(maxelem,21)
         integer idbc(maxnode,3)
         real*8  xord0(nnp), yord0(nnp), zord0(nnp)
!         real*8  dispful(900,3)

         real*8 prop(ippp,10),dd(6,6)
         integer nmprop(nel)

         real*8  stress(6,27), force20(3,20)
!         real*8  gforce(3000)
         real*8  dispful(maxdim,3)
         real*8  gforce(maxdof)
         real*8  uvw(3,20),xyz(3,20)
cccccccccccccccccc plastic  / rubber
         real*8  stresst(maxelem,6,27)
ccccccccccccccc


      call update_stress_rub_tian( 
     &                   stresst,
     &                   dispful, 
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xord0,yord0,zord0,
     &                   prop,nmprop,ippp,maxdim)

*        write(ilog,*)'@@ << in body_stress_rub >>'

         do i=1,neq
            gforce(i)=0.0
         enddo

*        rewind(igeo)

!          For each element, calculate the strain, stress at centroid
           do 50 n=1,nel
              mat=nmprop(n)
              neltype=npijkm(n,1)

!                 plate
                  e0=prop(mat,1)
                  g0=prop(mat,2)
                  t0=prop(mat,3)
                  r0=prop(mat,4)
                 pl0=prop(mat,5)
*               call dmat(e0,g0,dd)
c
c
*             write(iout,*)' [D] ',i
*             do ii=1,6
*                write(iout,82) (dd(ii,jj), jj=1,6)
*             enddo
c
              if (neltype .eq. 5) then
c                 tetra
c
              elseif (neltype .eq. 9) then
c                 hex8
c
              elseif (neltype .eq. 21) then
c                 hex20
                  kmax=neltype-1
                  do k=1,kmax
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
                      uvw(1,k)=dispful(node,1)
                      uvw(2,k)=dispful(node,2)
                      uvw(3,k)=dispful(node,3)
*                    duvw(1,k)=du  ful(node,1)
*                    duvw(2,k)=du  ful(node,2)
*                    duvw(3,k)=du  ful(node,3)
                  enddo
c
                  do i=1,6
                     do k=1,27
                        stress(i,k)=stresst(n,i,k)
                     enddo
                  enddo
c
               if (n .eq. 24) then
                  write(iout,*) '@@ stress elem ',n
                  do k=1,27  
                     write(iout,82) (stresst(n,i,k), i=1,6)
                  enddo
               endif
                  ired_int=3
                  ired_int=ri3
                  call fstress_rub_HEX20(neltype,ired_int, 
     &                            xyz,uvw,
     &                              stress,
     &                              force20)
*                 write(iout,*) '@@ force elem ',n
*                 do ii=1,20
*                    write(iout,82) (force20(jj,ii), jj=1,3)
*                 enddo
                  call assembFOR_HEX(gforce,force20,idbc,maxnode
     &                          ,npijkm,maxelem,n)
c
              endif
 50        continue
c          end of loop over elements
*             write(iout,*)'@@ gFORCE: '
*             write(iout,86) (gforce(k),k=1,neq)
c
 82        format(1x,40(g12.6,1x))
 83        format(1x,6(g12.6,1x))
 86        format(1x,6(g12.6,1x))
c
 999  continue
      return
      end

      subroutine update_stress_rub_tian( 
     &                   stresst_ip,
     &                   dispful, 
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xord0,yord0,zord0,
     &                   prop,nmprop,ippp,maxdim)
c
          implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3)
         real*8  xord0(nnp), yord0(nnp), zord0(nnp)
         real*8  dispful(maxdim,3)
c
         real*8 prop(ippp,10),dd(6,6)
         integer nmprop(nel)
c
         real*8 uvw(3,20),xyz(3,20)
         real*8 stresst_ip(maxelem,6,27)
         real*8 straint_ip(maxelem,6,27),stressc_ip(maxelem,6,27)
ccccccccccccccc
c
         real*8 uu(60),ee(6),ud(9),ss(6)
         real*8 skk(3,3),def(3,3),scc(3,3),cc(3,3)
         real*8 BE(6,60),Bd(9,60)
         real*8 xg(4,4),wgt(4,4)
         real*8 force(3,20),ff(60)
      nall=-1
      iextra=0
ccccccccccccccc
c
c        gauss-legendre sampling points
         data xg/ 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     &       -.5773502691896d0, 
     &        .5773502691896d0, 0.0d0, 0.0d0, -.7745966692415d0, 0.0d0,
     &        .7745966692415d0,0.0d0, -.8611363115941d0,
     &       -.3399810435849d0,  .3399810435849d0,  .8611363115941d0 /
c        gauss-legendre weights
         data wgt/ 2.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0,
     &        0.0d0, 0.0d0, .5555555555556d0, .8888888888889d0,
     &      .5555555555556d0, 0.0d0, .3478548451375d0, .6521451548625d0,
     &      .6521451548625d0, .3478548451375d0 /
c
c
*          write(ilog,*)'@@ << in update_stress_rub >>'
c
c          For each element, calculate the strain, stress at centroid
           if (nall .lt. 0) then
               nel1=1
               nel2=nel
           else
               nel1=nall
               nel2=nall
           endif
*          do 50 n=1,nel
           do 50 n=nel1,nel2
              mat=nmprop(n)
              neltype=npijkm(n,1)
c
c                 plate
                  e0=prop(mat,1)
                  g0=prop(mat,2)
                  t0=prop(mat,3)
                  r0=prop(mat,4)
                 pl0=prop(mat,5)
c
c
              if (neltype .eq. 5) then
c                 tetra
c
              elseif (neltype .eq. 9) then
c                 hex8
c
              elseif (neltype .eq. 21) then
c                 hex20
                  kmax=neltype-1
                  do k=1,kmax
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
                      uvw(1,k)=dispful(node,1)
                      uvw(2,k)=dispful(node,2)
                      uvw(3,k)=dispful(node,3)
*                    duvw(1,k)=du  ful(node,1)
*                    duvw(2,k)=du  ful(node,2)
*                    duvw(3,k)=du  ful(node,3)
*                 if (node .eq. 177 .OR. node .eq. 144) then
*                    write(iout,* ) ' u 177 ',node         
*                    write(iout,183) node, (dispful(node,m), m=1,3)
*                 endif
                  enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 update at each integration point
c
                  nn=0
                  do i=1,kmax
                     do j=1,3
                        nn=nn+1
                        uu(nn)=uvw(j,i)
                     enddo
                  enddo
                  if (n .eq. 24) then
                     write(iout,* ) ' uu '         
                     write(iout,182) (uu(k), k=1,60)
                     do i=1,3
                        write(iout,182) (uvw(i,k), k=1,20)
                     enddo
                  endif
c
                  do i=1,kmax*3
                     ff(i) = 0.0 
                  enddo
                  ired_int=3
                  ired_int=ri3
                  nint=ired_int
                  ip_n=0
                        do 80 lz=1,nint
                     do 80 ly=1,nint
                  do 80 lx=1,nint
                     ip_n=ip_n+1
                           ri=xg(lx,nint)
                           si=xg(ly,nint)
                           ti=xg(lz,nint)
c
c                    deriv oper and jacobian
                     call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
c                    displacement gradient
                     do i=1,9
                        sum=0.0
                        do k=1,kmax*3
                           sum=sum+Bd(i,k)*uu(k)
                        enddo
                        ud(i)=sum
                     enddo
                    iconstit=1
                     call rub_mat(ud,e0,g0,ss,dd,iconstit)
*                    if (ip_n  .eq. 27) then
*                        write(iout,*) ' [D] ',n,lx,ly,lz
*                        do i=1,6
*                           write(iout,82) (dd(i,j), j=1,6)
*                        enddo
*                     write(iout,* ) ' ud ss '         
*                     write(iout,82) (ud(k), k=1,9)
c*                    write(iout,82) (ee(k), k=1,6)
*                     write(iout,82) (ss(k), k=1,6)
*                     endif
c
*              iconstit=2
*                    if (iconstit .eq. 1) then
*                        call dmat(e0,g0,dd)
*                        call elast_lin(ud,dd,ss)
*                    elseif (iconstit .eq. 2) then
*                        a10=e0
*                        a01=g0
*                        call elast_rub(ud,e0,g0,ss)
*                    endif
*                    write(iout,* ) ' ES '         
*                    write(iout,82) (ud(k), k=1,9)
*                    write(iout,82) (ee(k), k=1,6)
*                    write(iout,82) (ss(k), k=1,6)
c
                     do i=1,6
                        stresst_ip(n,i,ip_n)=ss(i)
                     enddo
c
c             extra stuff
              if (iextra .eq. 1) then
c                 add contib to force
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
                  do 70 j=1,kmax*3
                     sum = 0.0
                     do k=1,6
                        sum = sum + ss(k)*BE(k,j)
                     enddo
                     ff(j) = ff(j) + sum*wt
 70               continue
c
c                 Lagrange strain
         ee(1)= ud(1) + ( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
         ee(2)= ud(5) + ( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
         ee(3)= ud(9) + ( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
         ee(4)= ud(2)+ud(4) + ud(1)*ud(2)+ud(4)*ud(5) 
     &                        + ud(7)*ud(8)
         ee(5)= ud(6)+ud(8) + ud(2)*ud(3)+ud(5)*ud(6) 
     &                        + ud(8)*ud(9)
         ee(6)= ud(3)+ud(7) + ud(1)*ud(3)+ud(4)*ud(6) 
     &                        + ud(7)*ud(9)
c
                     do i=1,6
                        straint_ip(n,i,ip_n)=ee(i)
                     enddo
c
c                     Cauchy stress
                      skk(1,1)=ss(1)
                      skk(2,2)=ss(2)
                      skk(3,3)=ss(3)
                      skk(1,2)=ss(4)
                      skk(1,3)=ss(6)      !
                      skk(2,3)=ss(5)      !
                      skk(2,1)=skk(1,2)
                      skk(3,1)=skk(1,3)
                      skk(3,2)=skk(2,3)
c
                      do j=1,3
                         def(1,j)=ud(0+j)
                         def(2,j)=ud(3+j)
                         def(3,j)=ud(6+j)
                      enddo
                         def(1,1)=def(1,1)+1.0
                         def(2,2)=def(2,2)+1.0
                         def(3,3)=def(3,3)+1.0
                zjj = def(1,1)*(def(2,2)*def(3,3)-def(3,2)*def(2,3))
     &               -def(1,2)*(def(2,1)*def(3,3)-def(3,1)*def(2,3))
     &               +def(1,3)*(def(2,1)*def(3,2)-def(3,1)*def(2,2))
c                     write(ilog,*)'@@ J: ', zjj
                      do i=1,3
                         do j=1,3
                            sum=0.0
                            do k=1,3
                               sum=sum + skk(i,k)*def(j,k)
                            enddo
                            cc(i,j)=sum
                         enddo
                      enddo
                      do i=1,3
                         do j=1,3
                            sum=0.0
                            do k=1,3
                               sum=sum + def(i,k)*cc(k,j)
                            enddo
*                           scc(i,j)=sum/det
                            scc(i,j)=sum/zjj
                         enddo
                       enddo
                         stressc_ip(n,1,ip_n)=scc(1,1)
                         stressc_ip(n,2,ip_n)=scc(2,2)
                         stressc_ip(n,3,ip_n)=scc(3,3)
                         stressc_ip(n,4,ip_n)=scc(1,2)
                         stressc_ip(n,6,ip_n)=scc(1,3)  !
                         stressc_ip(n,5,ip_n)=scc(2,3)  !
              endif
 80             continue
c               bottom loop over IP
                  nn=0
                  do j=1,kmax
                     do i=1,3
                        nn=nn+1
                        force(i,j)=ff(nn)
                     enddo
                  enddo
*                 write(iout,*)' force'
*                 do i=1,3
*                 write(iout,182)(force(i,j), j=1,10)
*                 enddo
c
              endif
c             bottom diff elems
 50        continue
c          end of loop over elements
*          write(iout,*) ' [D] ',lx,ly,lz
*          do i=1,6
*             write(iout,82) (dd(i,j), j=1,6)
*          enddo
c
 82        format(1x,40(g12.6,1x))
 83        format(1x,6(g12.6,1x))
 86        format(1x,6(g12.6,1x))
 866       format(1x,6(g12.6,1x),a)
 182       format(1x,90(g12.6,1x))
 183       format(1x,i5,1x,90(g12.6,1x))
c
 999  continue
      return
      end
