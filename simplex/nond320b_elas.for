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
c
          implicit real*8 (a-h,o-z)
          include 'commons.std'
!
         !--------------------------- F.-B. Tian
         include "mpif.h"
         INTEGER ::  IERR,IPES,ierror,npes,nrank,proc_m
         integer ::  istatus(MPI_STATUS_SIZE)
         integer ::  icvg
         !---------------------------
!
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=45000)
         integer nmap [allocatable] (:)
         real*8  xyz_s [allocatable] (:)
         real*8  xyz_f [allocatable] (:),contact_f [allocatable] (:)
         real*8  xyz_f_temp [allocatable] (:)
         real*8  A0 [allocatable] (:),r0 [allocatable] (:)
c
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
         real*8 dispfult(maxnode,3),dispfults(maxnode,3)
         real*8 disptm(neq),disptp(neq),dispts(neq)
         real*8 disptms(neq)

!=================for vel calculation
!	   real*8 dispt_1(neq),dispt_0(neq),dispt2(neq),dispt1(neq)
!==========F.-B. Tian
           real*8 dispti(neq),dispfulti(maxnode,3),disptmi(neq)
	   real*8 xyzti(maxxyz),xyzts(maxxyz)
!===================

           real*8 dload(neq),floadt(3*maxnode)

         maxdim=maxnode
         maxdof=3*maxdim
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          write(ilog,*)'@@  non_3D_exp ver270  November 2007 '

!        algorithm parameters
!          write(*,*) maxstiff,maxnode,neq,' ms'
!          stop
 
         kref=2



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
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         dt=deltat
!==========================================================================
! Communicate with the flow solver

         irankf = 0     ! the processor from the flow side

         call MPI_SEND(deltat, 1, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(is_res, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(irestart, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(npt, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)
         
! 
!      open(unit=1,file=fdn//'nmap.dat')
       open(unit=1,file='nmap.dat')
       read(1,*)nshare
       allocate(nmap(nshare))
       allocate(xyz_s(3*nshare))
       allocate(xyz_f(3*nshare))
       allocate(contact_f(nshare))
       allocate(xyz_f_temp(3*nshare))
       allocate(A0(nshare))
       allocate(r0(3*nshare))
       xyz_f = 0.0
       do i=1,nshare
        read(1,*)j,nmap(i)
!,A0(i),r0(i)
        xyz_s(i)=xyzt(nmap(i))
        xyz_s(nshare+i)=xyzt(maxnode+nmap(i))
        xyz_s(2*nshare+i)=xyzt(2*maxnode+nmap(i))
       enddo
       close(1)
!
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

!           form lumped mass
            ilump=1
!           write(*,*) maxstiff,maxnode,neq,' ms'
            write(ilog,*)'@@ FORMing lumped mass matrix '
!            call zztimer(ilog,'--> FORMmass')
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

!        Initialize  times, disp
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

!         call write_coordHEX20(xyzt,                ! Fangbao Tian
!     &                     npijkm,maxnode,maxelem,maxxyz,0)!
         call write_coord(xyzt,maxelem,                ! Fangbao Tian
     &                     npijkm,maxnode,maxxyz,0)



!        BIG TIME LOOP
!============================
         ramp=2.0
         kvelc=0
         do 200 itime= istart0,iend0
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
            dt0=deltat/isub
            ksub=0

!       save the current itime results
        xyzts=xyzt
        dispfults=dispfult
        disptms=disptm
        dispts=dispt

2013        continue
            ksub=ksub+1
            time = real(itime-2)*deltat+dt0*ksub
            pt1=1.0/(2.0*dt0)
            pt2=1.0/dt0**2
!           call zzcount_2(itime,inc,icarriage)

!           set load increment
!           Fmag says where the load is applied
!           Done this way in case distributed load applied
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            do i= 1, neq  
               pload(i) =  tf1*(fmag(i)*scale_p1 + grav(i)*scale_p2) 
               pload_tmp(i) = pload(i)
            enddo
!============================

 130    continue
        if (ksub .eq. isub+1) then
            xyzti=xyzts
            dispfulti=dispfults
            disptmi=disptms
            dispti=dispts
            ksub=1
        else
            xyzti=xyzt
            dispfulti=dispfult
            disptmi=disptm
            dispti=dispt
        endif
!------------------------------------------------------
! Communicate with the flow solver to receive the load
! for each node.  --F.-B. Tian

        if(ksub .eq. isub) then

         xyz_f_temp = xyz_f  ! dai
         xyz_f = 0.0

         call MPI_SEND(itime-1, 1, MPI_INTEGER, irankf, 1,  
     &                 MPI_COMM_WORLD,istatus,ierr)

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

         endif
!=============================
        !reassign the distributed force to nodal load
       sc       = 1.0d0 - exp(-ramp*time) ! load ramp factor

       !contact force
       do i=1,nshare
          contact_f(i)=0.0d0
       enddo

       do i=1,nshare
       ii=nmap(i)
       ij=idbc(ii,1)       
       pload(ij)=pload_tmp(ij)+xyz_f(i)*sc
       ij=idbc(ii,2)  
       pload(ij)=pload_tmp(ij)+xyz_f(nshare+i)*sc
       ij=idbc(ii,3) 
       pload(ij)=pload_tmp(ij)+xyz_f(nshare*2+i)*sc+contact_f(i)
       enddo
!============================


!             forces from body stresses
              call body_stress_3D( dispfulti, iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
!    &                    xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadt,iter,maxdim,maxdof
     &                           )

              do i= 1, neq  
                 aterm =  (pt2)*(2*dispti(i)-disptmi(i))
                 vterm =  (pt1)*disptmi(i)
                 dload(i) = pload(i) - floadt(i) 
     &                      + mass(i)*aterm + dmm*mass(i)*vterm 
                 coeff = pt2 + 1.0*dmm*pt1
                 disptp(i) = dload(i)/( coeff*mass(i))
              enddo
              do i= 1, neq
              acc(i)=(disptp(i)-2*dispti(i)+disptmi(i))*pt2
              vel(i)=(disptp(i)-disptmi(i))*pt1
              enddo
              do i= 1, neq
                 disptmi(i) = dispti(i)
                 dispti(i)  = disptp(i)
              enddo
             call update_geom_3D(dispti, dispfulti,maxdim, idbc,maxnode,
     &                     xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                     xyzti(1),xyzti(maxnode+1),xyzti(maxnode*2+1))

!	do i=1,nshare
!       ii=nmap(i)
!       ij=idbc(ii,3)
!       tempy=xyzti(2*maxnode+ii)
!       if(tempy.gt. (-0.05))then !siyuan
!       if (tempy .lt. -0.029 .and. tempy .gt. -0.039) then
!          xyzti(2*maxnode+ii)=-0.029
!            else if (tempy .gt. -0.049 .and. tempy .lt. -0.039) then
!               xyzti(2*maxnode+ii)=-0.049
!        xyzti(2*maxnode+ii)=-0.05
!        disptp(ij)=xyzti(2*maxnode+ii)-xyz0(2*maxnode+ii)
!        disptmi(ij) = dispti(ij)
!        dispti(ij)  = disptp(ij)
!	 endif
!	enddo               !siyuan
!////////////////////////////////////////////////////////
!-----------------------------------------------------
! Communicate with the flow solver and send the velocity and displacement
! --F.-B. Tian
         icvg = 1
  
         if( ksub.eq.isub) then
!================
!          acc=(disptp-2*dispti+disptmi)*pt2
!          vel=(disptp-disptmi)*pt1
!-===============
          print*, 'Nonstad sending ... '             
         ! send the node position
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
         endif ! (ksub.eq.isub)

         ! send the velcocity
         if( ksub.eq.isub) then
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

! debug lines
!         call write_coord(xyzti,maxelem,                ! Fangbao Tian
!     &                     npijkm,maxnode,maxxyz,itime-1)

         !write(*,666) 'xvel in solv: ',xyz_v(1:10)
         !write(*,666) 'yvel in solv: ',xyz_v(maxnode+1: maxnode+10)
         !write(*,666) 'zvel in solv: ',xyz_v(maxnode*2+1:maxnode*2+10)

         ! Receive signal from the flow solver about the FSI convergence.
         call MPI_RECV(icvg, 1, MPI_INTEGER, irankf, 1,
     &                 MPI_COMM_WORLD,istatus,ierr)

         endif !(ksub.eq.isub)

!================TEST FSI LOOP=================
                do 1722 n=1,nout
                   iprnode=inod(n,1)
                   if (inod(n,2) .eq. 0) then
                       velout(n)=dispt(iprnode)
                   elseif (inod(n,2) .eq. 1) then
                       velout(n)=vel(iprnode)
                   elseif (inod(n,2) .eq. 2) then
                       velout(n)=acc(iprnode)
                   endif
 1722            continue
                write(1111,'(15F15.8)') time,(velout(n), n=1,nout)
!===============================================

         ! If not converged, go to back to load update
         if(icvg .eq. 0) then
            ksub=ksub+1
            goto 130
         elseif (ksub .eq. isub) then
            write(*,*) 'NONSTAD: FSI iteration done.'
         endif
!////////////////////////////////////////////////////////	
      
        xyzt=xyzti
        dispfult=dispfulti
        disptm=disptmi
        dispt=dispti

      if(ksub.lt.isub)go to 2013

!          Print out results of interest at each large time step
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
!               write(idyn,182) time,tf1s,(velout(n), n=1,nout ),tf1
                write(*,*)'Itime=',itime

!               call write_coordHEX20(xyzt,                ! Fangbao Tian
!     &                     npijkm,maxnode,maxelem,maxxyz,itime-1)!
                call write_coord(xyzt,maxelem,                ! Fangbao Tian
     &                     npijkm,maxnode,maxxyz,itime-1)
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
     &  maxnode, maxxyz, xyzt,dispt,disptm, dispfult, vel,acc)
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
 182       format(1x,20(g12.6,1x),/,20(g12.6,1x),/,20(g12.6,1x))
 183       format(1x,i5,4x,3(g12.6,1x))
c
      return
      end
