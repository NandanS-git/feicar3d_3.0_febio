c
c     Modified for BIG TIME LOOP (Ye, Mar-2015)
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

         integer nmap [allocatable] (:)
         real*8  xyz_s [allocatable] (:)
         real*8  xyz_f [allocatable] (:),contact_f [allocatable] (:)
         real*8  xyz_f_temp [allocatable] (:)
         !real*8  xyz_f_old [allocatable] (:)
         real*8  A0 [allocatable] (:),r0 [allocatable] (:)
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=45000)

         integer nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer inod(ireadmax,2),mnod(iconelem)
c
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
 
         !integer FSI_on
!=================for vel calculation
!	   real*8 dispt_1(neq),dispt_0(neq),dispt2(neq),dispt1(neq)
!==========F.-B. Tian
           real*8 dispti(neq),dispfulti(maxnode,3),disptmi(neq)
	   real*8 xyzti(maxxyz),xyzt_s(maxxyz)
!===================

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

           real*8  cnt_d12 [allocatable] (:)
           real*8  cnt_d21 [allocatable] (:)
           real*8  cnt_d13 [allocatable] (:)
           real*8  cnt_d31 [allocatable] (:)
           real*8  cnt_d23 [allocatable] (:)
           real*8  cnt_d32 [allocatable] (:)

           integer ElemNum12 [allocatable] (:)
           integer ElemNum21 [allocatable] (:)
           integer ElemNum13 [allocatable] (:)
           integer ElemNum31 [allocatable] (:)
           integer ElemNum23 [allocatable] (:)
           integer ElemNum32 [allocatable] (:)

           integer clock0, clock1, clock_rate
           integer clock2, clock3,clock4,clock5,clock6,clock7

           integer :: nthreads,tid

           real*8 cf_min, cf_max
           real   time0
           INTEGER,DIMENSION(1)     :: cfmin_loc(1)
           INTEGER,DIMENSION(1)     :: cfmax_loc(1)

           !Ye,added Mass
           !real*8 vel_old(neq), vel_s(neq), mass2(neq)
           !real*8 dv(neq), dp(neq), acc_s(neq)
           !real*8 b_AddMass(neq), f_AddMass(neq)
           !real*8 rho_f
           !INTEGER iter_fsi

           !real*8 beta
 
           CHARACTER*20    :: fname2

       !Ye,allocate space for 2nd PK stress output at 27 interpolation nodes
       allocate (pks2(nel,27,6))
       allocate (strn(nel,27,6))
       allocate (cauchy(nel,27,6))

         maxdim=maxnode
         maxdof=3*maxdim
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          write(ilog,*)'@@  non_3D_exp ver270  November 2007 '
 
         kref = 2
         csi = 1.0E-10
         iter_fsi = 0
         beta = 1.0d0

!        Get things ready for time integration loop
!        write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
!     &                                        ' snap count '
!        call zzwrt('  -->  ')
!        read(ikbd,*) deltat,npt,iprcnt,isnap
         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
     &                       ' snap count | is_res | irestart | isub'
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,iprcnt,isnap, is_res,irestart,isub ! F.-B. Tian
!-------------------dt,    Nstep, nout, nscreen, res, ires===============

       !Ye,contact info input============================
       ! total # of surf tri elems
       totNumElem = 21168
       ! total number of surf pts
       nPtsBM = 10590

       allocate (BI_normX(nPtsBM))
       allocate (BI_normY(nPtsBM))
       allocate (BI_normZ(nPtsBM))

       allocate (cnt_flg12(nPtsBM))
       allocate (cnt_flg21(nPtsBM))
       allocate (cnt_flg13(nPtsBM))
       allocate (cnt_flg31(nPtsBM))
       allocate (cnt_flg23(nPtsBM))
       allocate (cnt_flg32(nPtsBM))

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

       allocate (BM_num(nPtsBM))
       allocate (mark(nPtsBM))

       allocate (ElemNeig(3,totNumElem))
       allocate (BM_x(nPtsBM),BM_y(nPtsBM),BM_z(nPtsBM))
       allocate (BM_x0(nPtsBM),BM_y0(nPtsBM),BM_z0(nPtsBM))
       allocate (BM_x1(nPtsBM),BM_y1(nPtsBM),BM_z1(nPtsBM))
       allocate (BM_x2(nPtsBM),BM_y2(nPtsBM),BM_z2(nPtsBM))
       allocate (BM_x3(nPtsBM),BM_y3(nPtsBM),BM_z3(nPtsBM))
       allocate (ElemCentx(totNumElem),ElemCenty(totNumElem),
     &           ElemCentz(totNumElem))
       allocate (ElemNorm_x(totNumElem),ElemNorm_y(totNumElem),
     &           ElemNorm_z(totNumElem))

       allocate (xNormBM(nPtsBM),yNormBM(nPtsBM),zNormBM(nPtsBM))

       cnt_flg12 = 0
       cnt_flg21 = 0
       cnt_flg13 = 0
       cnt_flg31 = 0
       cnt_flg23 = 0
       cnt_flg32 = 0

       cnt_d = 1.0E+5

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
       !cnt_n = 2709
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
       !allocate(xyz_f_old(3*nshare))
       allocate(A0(nshare))
       allocate(r0(3*nshare))
       allocate(cnt_fx(nshare))
       allocate(cnt_fy(nshare))
       allocate(cnt_fz(nshare))
       allocate(cnt_force(nshare))
       xyz_f = 0.0d0
       xyz_f_temp = 0.0d0
       !xyz_f_old = 0.0d0
       do i=1,nshare
          read(1,*)j,nmap(i),A0(i),r0(i)
       enddo
       close(1)

         if(is_res.eq.1)then
         call nond_3D_exp_restart_read(time,itime0,tf1,tf1s,k1old, 
     & maxnode,maxxyz,nPtsBM,xyzt,dispt,disptm,dispfult,vel,acc,xyz_f)
         !Ye,debug=========
         !do i=1,10482
         !   write(500,'(1(1X,I6),3(1X,1PE13.5))')i,xyz_f(i),
!     &           xyz_f(nshare+i),xyz_f(2*nshare+i)
!         enddo
         !=================
         time0 = time
         !write(*,*)'!!',time,itime0
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
           time0 = 0.0
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
         print*,'Send to flow proc ID:', irankf

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
        
         call MPI_SEND(xyz_s(1), nshare, 
     &        MPI_DOUBLE_PRECISION,
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

            !form lumped mass
            ilump=1
            write(ilog,*)'@@ FORMing lumped mass matrix '
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
c
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

!        Initialize  times, disp
         if(is_res.eq.1)then


         else
           write(*  ,*) '@@ INITIAL disps set to zero'
           do i= 1, neq
              dispt(i)=0.0
              vel(i)=0.0
              acc(i)=-dmm*vel(i)!dmm=C1/solid_density
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
           write(*,*)'!!',istart0,iend0
        else

           k1old=1
           tf1s=tforce1(1,2)*scale_p1
           tf1=tforce1(1,2)
           write(idyn,182) time,tf1s,(velout(n), n=1,nout ),tf1
 
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

!        call write_coordHEX20(xyzt,                ! Fangbao Tian
!     &                     npijkm,maxnode,maxelem,maxxyz,0)!
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

            !Save the current itime results
            xyzt_s    =xyzt
            dispfult_s=dispfult
            disptm_s  =disptm
            dispt_s   =dispt 

            !Save the vel & force at beginning of timestep
            !vel_old = vel
            !xyz_f_old = xyz_f

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

!          call MPI_RECV(FSI_on, 1, MPI_INTEGER,
!     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

!          call MPI_RECV(iter_fsi, 1, MPI_INTEGER, irankf, 1,
!     &                 MPI_COMM_WORLD,istatus,ierr)

!          call MPI_RECV(rho_f, 1, MPI_DOUBLE_PRECISION,
!     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(xyz_f(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(xyz_f(nshare+1  ), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          call MPI_RECV(xyz_f(nshare*2+1), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)


          !Ye,added mass============
          !if (itime.gt.itime0) then
          ! do i=1,nshare
          !   ii=nmap(i)
          !   ij=idbc(ii,1)
          !   if (ij.ne.0) then
          !      dp(ij) = xyz_f(i) - xyz_f_old(i)
          !   endif
          !   ij=idbc(ii,2)
          !   if (ij.ne.0) then
          !      dp(ij) = xyz_f(nshare+i) - xyz_f_old(nshare+i)
          !   endif
          !   ij=idbc(ii,3)
          !   if (ij.ne.0) then
          !      dp(ij) = xyz_f(nshare*2+i) - xyz_f_old(nshare*2+i)
          !   endif
          ! enddo
          !endif
          !=========================

          resMax_FSIp = 0.0d0
          do i=1,nshare*3
            !Ye,debug
            !if (i.eq.4180) then
            !   write(44,'(1E16.8,1X,i6,1X,7(E16.8,1X))')time,i,
!     &              relax_FSIp,
!     &              xyz_f(i),xyz_f(nshare+i),xyz_f(2*nshare+i),
!     &              xyz_f_temp(i),xyz_f_temp(nshare+i),
!     &              xyz_f_temp(2*nshare+i)
            !endif
            resP = xyz_f(i)-xyz_f_temp(i)
            xyz_f(i) = xyz_f_temp(i)+ resP*relax_FSIp
            resMax_FSIp  = max(resMax_FSIp,   abs(resP))
          enddo

          !Ye,debug
!          write(729,'(1(E16.8,1X),i6,1X,i4,1X,3(E16.8,1X))')time,itime,
!     &         iter_fsi,xyz_f(934),xyz_f(nshare+934),xyz_f(nshare*2+934)

          call MPI_SEND(resMax_fsip, 1, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

          ! initialize the state for substepping
          xyzt    =xyzt_s
          dispfult=dispfult_s
          disptm  =disptm_s
          dispt   =dispt_s

         endif!(ksub .eq. 0)

         call system_clock(clock0) !start the clock for the time step

 130     continue  ! return point for substep
            ksub =ksub+1
            !time = real(itime-2)*deltat+dt0*ksub
            time = time0 + real(itime-itime0-1)*deltat+dt0*ksub
            pt1=1.0/(2.0*dt0)
            pt2=1.0/dt0**2
!           call zzcount_2(itime,inc,icarriage)

!           set load increment
!           Fmag says where the load is applied
!           Done this way in case distributed load applied
            call load_interp(tforce1,maxforce,k1old,time,tf1)

            !Ye,debug=======
            !IF ((itime.eq.istart0).AND.(ksub.eq.1)) THEN
            !do i=1, neq
            !   write(201,*)i,fmag(i),grav(i)
            !enddo
            !ENDIF
            !===============

            !reassign the distributed force to nodal load
            sc  = 1.0d0 - exp(-ramp*time) ! load ramp factor

            do i= 1, neq  
               pload(i) =  sc*(fmag(i)*scale_p1 + grav(i)*scale_p2) 
               pload_tmp(i) = pload(i)
            enddo

            !Ye,debug
            !IF ((itime.eq.istart0).AND.(ksub.eq.1)) THEN
            !do i=1, neq
            !   write(202,*)tf1,pload_tmp(i)
            !enddo
            !ENDIF

            ! initialize the current substep
            xyzti    =xyzt
            dispfulti=dispfult
            disptmi  =disptm
            dispti   =dispt

            !Ye,use the current position to calc triElem norm===
            do i=1,nPtsBM!xuhao
               BM_x0(i)=xyzt(BM_num(i))
               BM_y0(i)=xyzt(maxnode+BM_num(i))
               BM_z0(i)=xyzt(2*maxnode+BM_num(i))
            enddo

            !Ye,test contact symmetry===
            !CASE1:L1-L2
            flag12 = 0
            if (abs(BM_x0(19)-BM_x0(3549))<5.0d0)then!0.005
               flag12 = 1
               goto 500
            else
               xavg = (BM_x0(19)+BM_x0(3549))/2.0d0
            endif
            do i=1,nPtsBM!y&z
               BM_y1(i)=BM_y0(i)
               BM_z1(i)=BM_z0(i)
            enddo
            do i=1,nPtsBM/3!L1
               BM_x1(i)=BM_x0(i)+(xavg-BM_x0(19))
            enddo
            do i=nPtsBM/3+1,2*nPtsBM/3!L2
               BM_x1(i)=BM_x0(i)+(xavg-BM_x0(3549))
            enddo
            do i=2*nPtsBM/3+1,nPtsBM!L3
               BM_x1(i)=BM_x0(i)
            enddo

 500        continue

            !CASE2:L1-L3
            flag13 = 0
            if (abs(BM_x0(19)-BM_x0(7079))<5.0d0)then
               flag13 = 1
               goto 600
            else
               xavg = (BM_x0(19)+BM_x0(7079))/2.0d0
            endif
            do i=1,nPtsBM!y&z
               BM_y2(i)=BM_y0(i)
               BM_z2(i)=BM_z0(i)
            enddo
            do i=1,nPtsBM/3!L1
               BM_x2(i)=BM_x0(i)+(xavg-BM_x0(19))
            enddo
            do i=2*nPtsBM/3+1,nPtsBM!L3
               BM_x2(i)=BM_x0(i)+(xavg-BM_x0(7079))
            enddo
            do i=nPtsBM/3+1,2*nPtsBM/3!L2
               BM_x2(i)=BM_x0(i)
            enddo

 600        continue

            !CASE3:L2-L3
            flag23 = 0
            if (abs(BM_x0(3549)-BM_x0(7079))<5.0d0)then
               flag23 = 1
               goto 700
            else
               xavg = (BM_x0(3549)+BM_x0(7079))/2.0d0
            endif
            do i=1,nPtsBM!y&z
               BM_y3(i)=BM_y0(i)
               BM_z3(i)=BM_z0(i)
            enddo
            do i=1,nPtsBM/3!L1
               BM_x3(i)=BM_x0(i)
            enddo
            do i=2*nPtsBM/3+1,nPtsBM!L3
               BM_x3(i)=BM_x0(i)+(xavg-BM_x0(7079))
            enddo
            do i=nPtsBM/3+1,2*nPtsBM/3!L2
               BM_x3(i)=BM_x0(i)+(xavg-BM_x0(3549))
            enddo

 700        continue

            !Ye,debug
            !write(fname2,'(I7.7)')  itime
            !fname2 = 'adjust13_'//trim(fname2)
            !print*,'Write to file: ',fname2
            !OPEN(UNIT=234,FILE=fname2,STATUS='UNKNOWN')
            !WRITE(234,*)'TITLE="3D TRIANGULAR SURFACE DATA"'
            !WRITE(234,'(a)')'VARIABLES="X","Y","Z"'
            !WRITE(234,*)'ZONE T="unstruc"','N=',nPtsBM,
!     &               'E=',totNumElem,'F=FEPOINT  ET=TRIANGLE'
            !DO i=1,nPtsBM
            !   write(234,'(3(1X,1PE13.5))')BM_x2(i),BM_y2(i),BM_z2(i)
            !ENDDO
            !DO j=1,totNumElem
            !WRITE(234,'(3I10)')ElemNeig(1,j),ElemNeig(2,j),ElemNeig(3,j)
            !ENDDO
            !CLOSE(234)

            do i=1,nPtsBM
               BM_x(i)=BM_x0(i)
               BM_y(i)=BM_y0(i)
               BM_z(i)=BM_z0(i)
            enddo

            !Ye, translation will not change norm, but ElemCentx/y/z
            call calculate_arclength_norm_cnt(itime,istart0,ksub)

!            if (itime.eq.6651) then
!               write(211,*)iter_fsi
!               write(211,'(3F12.7)')ElemNorm_x(8220),ElemNorm_y(8220),
!     &                              ElemNorm_z(8220)
!               write(211,'(3F12.7)')xNormBM(3445),yNormBM(3445),
!     &                              zNormBM(3445)
!               write(211,*)'===='
!            endif
!=============================
       !reassign the distributed force to nodal load
       sc  = 1.0d0 - exp(-ramp*time) ! load ramp factor

       !Ye, contact force initialization
       do i=1,nPtsBM
          cnt_fx(i)=0.0d0
          cnt_fy(i)=0.0d0
          cnt_fz(i)=0.0d0
          cnt_force(i)=0.0d0
       enddo

       !IF (FSI_on.eq.0) GOTO 102

       do i=1,cnt_n/3
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

       call system_clock(clock2, clock_rate)
       !write(*,*)'Nonstad: entering contact force.'

       if (flag12.eq.0)then
          do i=1,nPtsBM
             BM_x(i)=BM_x1(i)
             BM_y(i)=BM_y1(i)
             BM_z(i)=BM_z1(i)
          enddo
       else
          do i=1,nPtsBM
             BM_x(i)=BM_x0(i)
             BM_y(i)=BM_y0(i)
             BM_z(i)=BM_z0(i)
          enddo
       endif
       !--Centroid of the each element, replace the one calculated!
       do m=1,totNumElem
          node1 = ElemNeig(1,m)
          node2 = ElemNeig(2,m)
          node3 = ElemNeig(3,m)

          ElemCentx(m)=(BM_x(node1)+BM_x(node2)+BM_x(node3))/3.0d0
          ElemCenty(m)=(BM_y(node1)+BM_y(node2)+BM_y(node3))/3.0d0
          ElemCentz(m)=(BM_z(node1)+BM_z(node2)+BM_z(node3))/3.0d0
       enddo 

       cnt_flg12 = 0

#ifdef _OPENMP
       call    OMP_SET_DYNAMIC(.FALSE.)
       call    OMP_SET_NUM_THREADS(16)
#endif

!$OMP PARALLEL PRIVATE(nthreads, tid)
!Obtain thread number
      !IF (ksub.eq.1) THEN
      !   tid = OMP_GET_THREAD_NUM()
      !   write(*,*) 'Hello from sub: tid=', tid

!Only master thread does this
         !IF (tid .EQ. 0) THEN
         !   nthreads = OMP_GET_NUM_THREADS()
         !   write(*,*)'Number of threads = ', nthreads
         !ENDIF
      !ENDIF
!$OMP END PARALLEL

            flg = 12

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum,f_ext,ii)
!$OMP DO
            !Ye,Contact Btw LEAFLET_1 & LEAFLET_2
            !calculate cnt_d(i)
            do i=1,cnt_n/3
               ii = cnt_xh(i)
               !tid = OMP_GET_THREAD_NUM()
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+5

               !IF ( (ksub.eq.1).or.(ksub.eq.2) )THEN
               !   write(*,'(a,I5,a,I5)'),'i1=', i, ' by tid ', tid
               !ENDIF

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii,flg
     &         ,nPtsBM/3+1,2*nPtsBM/3,cnt_d,ElemNum)

               cnt_d12(ii) = cnt_d
               ElemNum12(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               !calculate total external force at node i
               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -BI_normX(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -BI_normY(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -BI_normZ(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = BI_normX(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = BI_normY(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = BI_normZ(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,central upstream bending----
               IF (cnt_d-cnt_h.le.0.0d0)THEN
               if ( (ii.eq.2113).or.(ii.eq.  17).or.(ii.eq.2355).or.
     &              (ii.eq.  18).or.(ii.eq.1938).or.(ii.eq.  19).or.
     &              (ii.eq.1932).or.(ii.eq.  20).or.(ii.eq.1283).or.
     &              (ii.eq.  21).or.(ii.eq.3360) )then
                  contact_fx = -xyz_f(ii)
               endif
               ENDIF

               !Ye,interp the contact force
               cnt_fx12(i)=contact_fx*tf1
               cnt_fy12(i)=contact_fy*tf1
               cnt_fz12(i)=contact_fz*tf1
            enddo
!$OMP END PARALLEL

            cnt_flg21 = 0

            flg = 21

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum,f_ext,ii)
!$OMP DO
            do i=cnt_n/3+1,2*cnt_n/3
               ii = cnt_xh(i)
               !tid = OMP_GET_THREAD_NUM()
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+5

               !IF ( (ksub.eq.1).or.(ksub.eq.2) )THEN
               !   write(*,'(a,I5,a,I5)'),'i2=', i, ' by tid ', tid
               !ENDIF


               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii,flg
     &              ,1,nPtsBM/3,cnt_d,ElemNum)

               cnt_d21(ii) = cnt_d
               ElemNum21(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -BI_normX(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -BI_normY(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -BI_normZ(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = BI_normX(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = BI_normY(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = BI_normZ(ii)*cnt_k*(cnt_d-cnt_h)
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

               !Ye,central upstream bending----
               IF (cnt_d-cnt_h.le.0.0d0)THEN
               if ( (ii.eq.5643).or.(ii.eq.3547).or.(ii.eq.5885).or.
     &              (ii.eq.3548).or.(ii.eq.5468).or.(ii.eq.3549).or.
     &              (ii.eq.5462).or.(ii.eq.3550).or.(ii.eq.4813).or.
     &              (ii.eq.3551).or.(ii.eq.6890) )then
                  contact_fx = -xyz_f(ii)
               endif
               ENDIF

               !Ye,interp the contact force
               cnt_fx21(i-cnt_n/3)=contact_fx*tf1
               cnt_fy21(i-cnt_n/3)=contact_fy*tf1
               cnt_fz21(i-cnt_n/3)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL

       if (flag13.eq.0)then
          do i=1,nPtsBM
             BM_x(i)=BM_x2(i)
             BM_y(i)=BM_y2(i)
             BM_z(i)=BM_z2(i)
          enddo
       else
          do i=1,nPtsBM
             BM_x(i)=BM_x0(i)
             BM_y(i)=BM_y0(i)
             BM_z(i)=BM_z0(i)
          enddo
       endif
       !--Centroid of the each element, replace the one calculated!
       do m=1,totNumElem
          node1 = ElemNeig(1,m)
          node2 = ElemNeig(2,m)
          node3 = ElemNeig(3,m)

          ElemCentx(m)=(BM_x(node1)+BM_x(node2)+BM_x(node3))/3.0d0
          ElemCenty(m)=(BM_y(node1)+BM_y(node2)+BM_y(node3))/3.0d0
          ElemCentz(m)=(BM_z(node1)+BM_z(node2)+BM_z(node3))/3.0d0
       enddo

       cnt_flg13 = 0

       flg = 13

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum,f_ext,ii)
!$OMP DO
            !Ye,Contact Btw LEAFLET_1 & LEAFLET_3
            !calculate cnt_d(i)
            do i=1,cnt_n/3
               ii = cnt_xh(i)
               !tid = OMP_GET_THREAD_NUM()
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+5

               !IF ( (ksub.eq.1).or.(ksub.eq.2) )THEN
               !   write(*,'(a,I5,a,I5)'),'i3=', i, ' by tid ', tid
               !ENDIF

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii,flg
     &         ,2*nPtsBM/3+1,nPtsBM,cnt_d,ElemNum)

               cnt_d13(ii) = cnt_d
               ElemNum13(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -BI_normX(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -BI_normY(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -BI_normZ(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = BI_normX(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = BI_normY(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = BI_normZ(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,central upstream bending----
               IF (cnt_d-cnt_h.le.0.0d0)THEN
               if ( (ii.eq.2113).or.(ii.eq.  17).or.(ii.eq.2355).or.
     &              (ii.eq.  18).or.(ii.eq.1938).or.(ii.eq.  19).or.
     &              (ii.eq.1932).or.(ii.eq.  20).or.(ii.eq.1283).or.
     &              (ii.eq.  21).or.(ii.eq.3360) )then
                  contact_fx = -xyz_f(ii)
               endif
               ENDIF

               !Ye,interp the contact force
               cnt_fx13(i)=contact_fx*tf1
               cnt_fy13(i)=contact_fy*tf1
               cnt_fz13(i)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL

            cnt_flg31 = 0

            flg = 31

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum,f_ext,ii)
!$OMP DO
            do i=2*cnt_n/3+1,cnt_n
               ii = cnt_xh(i)
               !tid = OMP_GET_THREAD_NUM()
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+5

               !IF ( (ksub.eq.1).or.(ksub.eq.2) )THEN
               !   write(*,'(a,I5,a,I5)'),'i4=', i, ' by tid ', tid
               !ENDIF

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii,flg
     &              ,1,nPtsBM/3,cnt_d,ElemNum)

               cnt_d31(ii) = cnt_d
               ElemNum31(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -BI_normX(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -BI_normY(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -BI_normZ(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = BI_normX(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = BI_normY(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = BI_normZ(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,central upstream bending----
               IF (cnt_d-cnt_h.le.0.0d0)THEN
               if ( (ii.eq.10420).or.(ii.eq.7081).or.(ii.eq.8343).or.
     &              (ii.eq. 7080).or.(ii.eq.8992).or.(ii.eq.7079).or.
     &              (ii.eq. 8998).or.(ii.eq.7078).or.(ii.eq.9415).or.
     &              (ii.eq. 7077).or.(ii.eq.9173) )then
                  contact_fx = -xyz_f(ii)
               endif
               ENDIF

               !Ye,interp the contact force
               cnt_fx31(i-2*cnt_n/3)=contact_fx*tf1
               cnt_fy31(i-2*cnt_n/3)=contact_fy*tf1
               cnt_fz31(i-2*cnt_n/3)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL

       if (flag23.eq.0)then
          do i=1,nPtsBM
             BM_x(i)=BM_x3(i)
             BM_y(i)=BM_y3(i)
             BM_z(i)=BM_z3(i)
          enddo
       else
          do i=1,nPtsBM
             BM_x(i)=BM_x0(i)
             BM_y(i)=BM_y0(i)
             BM_z(i)=BM_z0(i)
          enddo
       endif
       !--Centroid of the each element, replace the one calculated!
       do m=1,totNumElem
          node1 = ElemNeig(1,m)
          node2 = ElemNeig(2,m)
          node3 = ElemNeig(3,m)

          ElemCentx(m)=(BM_x(node1)+BM_x(node2)+BM_x(node3))/3.0d0
          ElemCenty(m)=(BM_y(node1)+BM_y(node2)+BM_y(node3))/3.0d0
          ElemCentz(m)=(BM_z(node1)+BM_z(node2)+BM_z(node3))/3.0d0
       enddo

       cnt_flg32 = 0

       flg = 32

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum,f_ext,ii)
!$OMP DO
            !Ye,Contact Btw LEAFLET_3 & LEAFLET_2
            !calculate cnt_d(i)
            do i=2*cnt_n/3+1,cnt_n
               ii = cnt_xh(i)
               !tid = OMP_GET_THREAD_NUM()
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+5

               !IF ( (ksub.eq.1).or.(ksub.eq.2) )THEN
               !   write(*,'(a,I5,a,I5)'),'i5=', i, ' by tid ', tid
               !ENDIF

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii,flg
     &         ,nPtsBM/3+1,2*nPtsBM/3,cnt_d,ElemNum)

               cnt_d32(ii) = cnt_d
               ElemNum32(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -BI_normX(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -BI_normY(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -BI_normZ(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = BI_normX(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = BI_normY(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = BI_normZ(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,central upstream bending----
               IF (cnt_d-cnt_h.le.0.0d0)THEN
               if ( (ii.eq.10420).or.(ii.eq.7081).or.(ii.eq.8343).or.
     &              (ii.eq. 7080).or.(ii.eq.8992).or.(ii.eq.7079).or.
     &              (ii.eq. 8998).or.(ii.eq.7078).or.(ii.eq.9415).or.
     &              (ii.eq. 7077).or.(ii.eq.9173) )then
                  contact_fx = -xyz_f(ii)
               endif
               ENDIF

               !Ye,interp the contact force
               cnt_fx32(i-2*cnt_n/3)=contact_fx*tf1
               cnt_fy32(i-2*cnt_n/3)=contact_fy*tf1
               cnt_fz32(i-2*cnt_n/3)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL

            cnt_flg23 = 0

            flg = 23

!$OMP PARALLEL
     & PRIVATE(i,cnt_d,contact_fx,contact_fy,contact_fz
     & ,ElemNum,f_ext,ii)
!$OMP DO
            do i=cnt_n/3+1,2*cnt_n/3
               ii = cnt_xh(i)
               !tid = OMP_GET_THREAD_NUM()
               contact_fx = 0.0d0
               contact_fy = 0.0d0
               contact_fz = 0.0d0
               cnt_d = 1.0E+5

               !IF ( (ksub.eq.1).or.(ksub.eq.2) )THEN
               !   write(*,'(a,I5,a,I5)'),'i6=', i, ' by tid ', tid
               !ENDIF

               !only consider LV side markers
               IF( (mark(ii).eq.0).OR.( (idbc(nmap(ii),1).eq.0).or.
     &         (idbc(nmap(ii),2).eq.0).or.(idbc(nmap(ii),3).eq.0)))
     &         CYCLE

               call calc_bodyIntercept_Unstruc_cnt(ii,flg
     &         ,2*nPtsBM/3+1,nPtsBM,cnt_d,ElemNum)

               cnt_d23(ii) = cnt_d
               ElemNum23(ii) = ElemNum

               IF (cnt_d.ge.cnt_c) CYCLE

               f_ext = sqrt(xyz_f(ii)**2+xyz_f(nshare+ii)**2+
     &                      xyz_f(nshare*2+ii)**2)+csi

               IF (cnt_d-cnt_h.gt.0.0d0) THEN
                  contact_fx = -BI_normX(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fy = -BI_normY(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
                  contact_fz = -BI_normZ(ii)*f_ext*exp(-cnt_k*
     &                          (cnt_d-cnt_h)/(csi*f_ext))
               ELSE
                  contact_fx = BI_normX(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(ii)
                  contact_fy = BI_normY(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare+ii)
                  contact_fz = BI_normZ(ii)*cnt_k*(cnt_d-cnt_h)
     &                        -xyz_f(nshare*2+ii)
               ENDIF

               !Ye,central upstream bending----
               IF (cnt_d-cnt_h.le.0.0d0)THEN
               if ( (ii.eq.5643).or.(ii.eq.3547).or.(ii.eq.5885).or.
     &              (ii.eq.3548).or.(ii.eq.5468).or.(ii.eq.3549).or.
     &              (ii.eq.5462).or.(ii.eq.3550).or.(ii.eq.4813).or.
     &              (ii.eq.3551).or.(ii.eq.6890) )then
                  contact_fx = -xyz_f(ii)
               endif
               ENDIF

               !Ye,interp the contact force
               cnt_fx23(i-cnt_n/3)=contact_fx*tf1
               cnt_fy23(i-cnt_n/3)=contact_fy*tf1
               cnt_fz23(i-cnt_n/3)=contact_fz*tf1
               !=====================================
            enddo
!$OMP END PARALLEL


            !Ye, add the contact forces from two other leaflets
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

            call system_clock(clock3, clock_rate)
!            WRITE(*,*) ' Total time for contact force is:',
!     &                   REAL(clock3-clock2)/REAL(clock_rate)

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

              call system_clock(clock6, clock_rate)

!             forces from body stresses
              call body_stress_3D( dispfulti, iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
!    &                    xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadt,iter,maxdim,maxdof
     &                           )

              call system_clock(clock7, clock_rate)
!            WRITE(*,*) ' Total time for body stress is:',
!     &                   REAL(clock7-clock6)/REAL(clock_rate)

              !Ye,added mass====================
              !if (itime.eq.istart0) then
              !   do i=1, neq
              !      vel_old(i) = 0.0d0
              !        vel_s(i) = 0.0d0
              !        acc_s(i) = 0.0d0
              !    b_AddMass(i) = 0.0d0
              !   enddo
              !endif

              !do i=1, neq
              !   mass2(i) = mass(i) + b_AddMass(i) * rho_f
              !enddo
              !do i=1, neq
              !   f_AddMass(i) = b_AddMass(i)*rho_f*acc_s(i)
              !enddo
              !do i=1, neq
              !   pload(i) = pload(i) + tf1*f_AddMass(i)
              !enddo
              !==================================

              do i= 1, neq  
                 aterm =  (pt2)*(2*dispti(i)-disptmi(i))
                 vterm =  (pt1)*disptmi(i)
                 dload(i) = pload(i) - floadt(i) 
     &                      + mass(i)*aterm + dmm*mass(i)*vterm 
                 coeff = pt2 + 1.0*dmm*pt1
                 disptp(i) = dload(i)/(coeff*mass(i))
              enddo
              do i= 1, neq
                 acc(i) = (disptp(i)-2*dispti(i)+disptmi(i))*pt2
                 vel(i) = (disptp(i)-disptmi(i))*pt1
              enddo
              do i= 1, neq
                 disptmi(i) = dispti(i)
                 dispti(i)  = disptp(i)
              enddo
             call update_geom_3D(dispti, dispfulti,maxdim, idbc,maxnode,
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
                 !write(1111,'(100F15.8)') time,(velout(n), n=1,nout)
!===============================================

!////////////////////////////////////////////////////////	

        xyzt    =xyzti
        dispfult=dispfulti
        disptm  =disptmi
        dispt   =dispti

           !Print out results of interest within iteration
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

      ! End of current substep
      if(ksub.lt.isub)then
        call system_clock(clock4, clock_rate)
!            WRITE(*,*) ' Total time for current substep is:',
!     &                   REAL(clock4-clock2)/REAL(clock_rate)

        go to 130           ! go back to do another substep

      else  ! ksub = isub   ! substeps finished; so is the current time step within FSI.

             !Ye,added mass====================
             !if (itime.gt.istart0 .and. iter_fsi>=1) then
             !    do i=1, neq
             !       dv(i) = vel_s(i) - vel_old(i)
             !    enddo
                 !update b_AddMass
             !    do i=1, neq
             !       b_AddMass(i) = -dp(i)*deltat/(dv(i)*rho_f)
             !       b_AddMass(i) = abs( b_AddMass(i) )
             !       if(b_AddMass(i)*rho_f/mass(i) > 100.0d0) then
             !         b_AddMass(i) = 100.0d0*mass(i)/rho_f
             !       endif
             !    enddo
             !endif

              !Ye,disable the added mass effect
              !b_AddMass = 0.0d0

              !Ye,debug added mass
              !write(*,'(a,2(E12.6, 5X))') 'Maximum added mass:',
!     &             maxval(b_AddMass)*rho_f,maxval(mass)

              !do i= 1, neq
              !   vel_s(i) = vel(i)
              !   acc_s(i) = (vel(i) - vel_old(i))/deltat
              !enddo
             !==================================
  
!================
!          acc=(disptp-2*dispti+disptmi)*pt2
!        do i=1,neq
!           vel2(i) = (disptp(i)-dispt_s(i))/deltat
!        enddo

!        do n=1,nout
!           iprnode=inod(n,1)
!           if ( inod(n,2) .eq. 1 ) then
!              velout2(n)=vel2(iprnode)
!           else
!              velout2(n)=0.0
!           endif
!        enddo
!-===============

         !=================================================
         call system_clock(clock5, clock_rate)
         WRITE(*,*) ' Total time for solid (within FSI) is:',
     &               REAL(clock5-clock0)/REAL(clock_rate)

         !===========
         ! Communicate with the flow solver and send the velocity and displacement
         !===========
          print*, 'Nonstad sending ... '             

         ! pack and send the node position
          do i=1,nshare
            xyz_s(i)=xyzti(nmap(i))
            xyz_s(nshare+i)=xyzti(maxnode+nmap(i))
            xyz_s(2*nshare+i)=xyzti(2*maxnode+nmap(i))
!            if (i.eq.4180) then
!               write(66,'(1E16.8,1X,i6,1X,3(E16.8,1X))')time,i,
!     &              xyz_s(i),xyz_s(nshare+i),xyz_s(2*nshare+i)
!             endif
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
!             do i=1,nshare
!               ii=nmap(i)
!               ij=idbc(ii,1)
!               xyz_s(i)=vel(ij)
!               ij=idbc(ii,2)
!               xyz_s(nshare+i)=vel(ij)
!               ij=idbc(ii,3)
!               xyz_s(2*nshare+i)=vel(ij)
!               if (ii.eq.4669) then
!                  write(88,'(1E16.8,1X,i6,1X,3(E16.8,1X))')time,i,
!     &                 xyz_s(i),xyz_s(nshare+i),xyz_s(2*nshare+i)
!               endif
!             enddo
 
         !Ye, filter by using beta to reduce the oscillations in velocity before it is sent
         !to flow solver
         do i=1,nshare
            ii=nmap(i)
            ij=idbc(ii,1)
            if (ij.ne.0) then
               xyz_s(i)=vel(ij)
            else
               xyz_s(i)=0.0d0
            endif
            ij=idbc(ii,2)
            if (ij.ne.0) then
               xyz_s(nshare+i)=vel(ij)
            else
               xyz_s(nshare+i)=0.0d0
            endif
            ij=idbc(ii,3)
            if (ij.ne.0) then
               xyz_s(2*nshare+i)=vel(ij)
            else
               xyz_s(2*nshare+i)=0.0d0
            endif
            !if (ii.eq.934) then
            !   write(88,'(1E16.8,1X,i6,1X,3(E16.8,1X))')itime,i,
!     &              xyz_s(i),xyz_s(nshare+i),xyz_s(2*nshare+i)
!            endif
         enddo
 
         call MPI_SEND(xyz_s(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare+1  ), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare*2+1), nshare, 
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

! debug lines
!         call write_coord(xyzti,maxelem,                ! Fangbao Tian
!     &                     npijkm,maxnode,maxxyz,itime-1)

         !write(*,666) 'xvel in solv: ',xyz_v(1:10)
         !write(*,666) 'yvel in solv: ',xyz_v(maxnode+1: maxnode+10)
         !write(*,666) 'zvel in solv: ',xyz_v(maxnode*2+1:maxnode*2+10)

         ! Receive signal from the flow solver about the FSI convergence.
         call MPI_RECV(icvg, 1, MPI_INTEGER, irankf, 1,
     &                 MPI_COMM_WORLD,istatus,ierr)

         !Ye,debug,iterate only once
         !IF (itime.eq.2505) THEN
         !   icvg=1
         !ENDIF

         ! If not converged, go to back to load update
         if(icvg .eq. 0) then
            goto 2013  ! not converged; update force and resolve
         else
            write(*,*) 'NONSTAD: FSI iteration done.'
         endif

      endif  ! end of (ksub .lt. isub)
!=====================================

            !====================================
            !Ye,send contact force to flow slover
            CALL MPI_SEND(nshare, 1, MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_flg12, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_flg21, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_flg13, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_flg31, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_flg23, nshare,
     &           MPI_INTEGER,
     &           irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            CALL MPI_SEND(cnt_flg32, nshare,
     &           MPI_INTEGER,
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

!====================================

!          Print out results of interest at each time step
           do 174 n=1,nout
              iprnode=inod(n,1)
              irate  =inod(n,2)
              if (iprnode .eq. 0) then
                 velout(n) = 0.0
                 goto 174
              endif
              if (inod(n,2) .eq. 0) then
                 velout(n)=dispt(iprnode)
              elseif (inod(n,2) .eq. 1) then
                 velout(n)=vel(iprnode)
              elseif (inod(n,2) .eq. 2) then
                 velout(n)=acc(iprnode)
              endif
 174          continue
              tf1s=tf1*scale_p1
              write(idyn,182) time,tf1s,(velout(n), n=1,nout ),tf1


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
!                write(idyn,182) time,tf1s,(velout(n), n=1,nout ),tf1

!               call write_coordHEX20(xyzt,                ! Fangbao Tian
!     &                     npijkm,maxnode,maxelem,maxxyz,itime-1)!
                call write_coord(xyzt,maxelem,                ! Fangbao Tian
     &                     npijkm,maxnode,maxxyz,itime-1)

             !Ye,write 2nd PK stress to file
             write(fname2,'(I7.7)')  itime-1
             fname2 = 'pks2_'//trim(fname2)
             OPEN(UNIT=234,FILE=fname2,STATUS='UNKNOWN')
             DO i=1, nel
                DO j=1, 27
                 WRITE(234,'(6(1X,1PE13.5))')pks2(i,j,1),
     &               pks2(i,j,2),pks2(i,j,3),pks2(i,j,4),
     &               pks2(i,j,5),pks2(i,j,6)
                ENDDO
             ENDDO
             CLOSE(234)

             !Ye,write Cauchy stress to file
             write(fname2,'(I7.7)')  itime-1
             fname2 = 'cauchy_'//trim(fname2)
             OPEN(UNIT=234,FILE=fname2,STATUS='UNKNOWN')
             DO i=1, nel
                DO j=1, 27
                 WRITE(234,'(6(1X,1PE13.5))')cauchy(i,j,1),
     &               cauchy(i,j,2),cauchy(i,j,3),cauchy(i,j,4),
     &               cauchy(i,j,5),cauchy(i,j,6)
                ENDDO
             ENDDO
             CLOSE(234)

             !Ye,write strain to file
             write(fname2,'(I7.7)')  itime-1
             fname2 = 'strn_'//trim(fname2)
             OPEN(UNIT=234,FILE=fname2,STATUS='UNKNOWN')
             DO i=1, nel
                DO j=1, 27
                 WRITE(234,'(6(1X,1PE13.5))')strn(i,j,1),
     &               strn(i,j,2),strn(i,j,3),strn(i,j,4),
     &               strn(i,j,5),strn(i,j,6)
                ENDDO
             ENDDO
             CLOSE(234)

             endif

             if (isnap .eq. ksnap) then
                 write(isnp) time, neq, nsnapmax 
!                 do i= 1, neq  
!                    write(isnp) dispt(i)
!                 enddo
                 ksnap=0
             endif

!  Fang-Bao Tian
! Writing the restart file
       if(mod(itime-1,irestart) .eq. 0) then
        call nond_3D_exp_restart_write(time,itime,tf1,tf1s,k1old,
     &  maxnode,maxxyz,nPtsBM,xyzt,dispt,disptm,dispfult,vel,acc,xyz_f)
       endif
!       ij1=idbc(1312,1)
!       ij2=idbc(1388,1)
!       ij3=idbc(1445,1)
!      if(itime.eq.2)then
!        open(unit=11,file=fdn//'deflection.dat')
!        write(11,*)"variables=t x1 y1 z1 x2 y2 z2 x3 y3 z3"
!        write(11,*)time,xyzt(1312),xyzt(maxnode+1312),xyzt(2*maxnode+1312)
!     &               ,xyzt(1388),xyzt(maxnode+1388),xyzt(2*maxnode+1388)
!     &               ,xyzt(1445),xyzt(maxnode+1445),xyzt(2*maxnode+1445)
!      close(11)
!      open(unit=11,file=fdn//'vel_deflection.dat')
!	write(11,*)"variables=t v1 v2 v3"
!	write(11,*)time, vel(ij1),vel(ij2),vel(ij3)
!      close(11)
!        else
!        open(unit=11,file=fdn//'deflection.dat',position='append')
!        write(11,*)time,xyzt(1312),xyzt(maxnode+1312),xyzt(2*maxnode+1312)
!     &               ,xyzt(1388),xyzt(maxnode+1388),xyzt(2*maxnode+1388)
!     &               ,xyzt(1445),xyzt(maxnode+1445),xyzt(2*maxnode+1445)
!      close(11)
!      open(unit=11,file=fdn//'vel_deflection.dat',position='append')
!	write(11,*)time, vel(ij1),vel(ij2),vel(ij3)
!      close(11)
!      endif

       call system_clock(clock1, clock_rate)
       WRITE(*,*) ' Total time for entire step (from solid) is:',
     &            REAL(clock1-clock0)/REAL(clock_rate)

 200  continue
 210  continue
         write(*,'(a)')' @@'
c     END BIG TIME LOOP
c
 81          format(1x,a,1x,90(1x,g12.6))
 82          format(1x,90(1x,g12.6))
 83        format(1x,4(1x,g12.6),1x,a)
 84        format(1x,a,1x,2(1x,g12.6),1x,2(i5,1x))
 85          format(1x,90(1x,g12.6))
 86          format(1x, 6(1x,g12.6))
 182       format(60(g16.8,1x),/,60(g12.6,1x),/,60(g12.6,1x))
 183       format(1x,i5,4x,3(g12.6,1x))
c
      return
      end
