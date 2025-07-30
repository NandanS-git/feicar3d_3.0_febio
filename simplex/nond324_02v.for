c
c
c4400
c     NON-liner 3D EXPlicit incremental analysis 
      subroutine non_3D_exp_tet10(wk,pload, disp, fmag,
     &                        dispt , 
     &                   mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt)
c
          implicit real*8 (a-h,o-z)
          include 'commons.std'
!===========for FSI===============! 
          include "mpif.h"
          INTEGER ::  IERR,IPES,ierror,npes,nrank,proc_m
          integer ::  istatus(MPI_STATUS_SIZE)
          integer ::  icvg
!==========end====================! 
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=45000)
c
         integer nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer inod(ireadmax,2),mnod(iconelem)
c
         real*8 tforce1(maxforce,2)
         real*8 velout(ireadmax)
         real*8 mass(maxmass),wk(neq)
     &          ,pload(neq ),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(maxxyz)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 grav(maxxyz)
         real*8 dispt(neq)
         real*8 dispfult(maxnode,3)
         real*8 disptm(neq),disptp(neq)

!===========for FSI==================!
         integer nmap [allocatable] (:)
         real*8 xyz_s [allocatable] (:)
         real*8 xyz_f [allocatable] (:), contact_f [allocatable] (:)
         real*8 xyz_f_temp [allocatable] (:)
         real*8 A0 [allocatable] (:), r0 [allocatable] (:)
         real*8 dispti(neq), dispfulti(maxnode,3), disptmi(neq)
         real*8 xyzti(maxxyz), pload_tmp(neq)
         real*8 tempy, tempd, delta
         real*8 veli(neq), acci(neq)
         real*8 zcenter, p0, deltad, dmin, dclose
         real*8 diffdamp [allocatable] (:)
         real*8 diffdmm(neq)
!===========end======================!
c
         real*8 dload(neq),floadt(3*maxnode)
         print*,'@@  non_3D_exp ver270  November 2007 '
         maxdim=maxnode
         maxdof=3*maxnode
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         write(ilog,*)'@@  non_3D_exp ver270  November 2007 '
c
         kref=2

c        Get things ready for time integration loop
         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
     &                    ' snap count | is_res | irestart | isub'
         call zzwrt('  -->  ')
         read (ikbd,*) deltat,npt,iprcnt,isnap,is_res,irestart,isub
        
         if (is_res .eq. 1) then
              call nond_3D_exp_restart_read(time, itime0,
     &                        tf1, tf1s, k1old, maxnode, maxxyz,
     &                        xyzt, dispt, disptm, dispfult, vel, acc)
         else 
c
           xn=0.0
           yn=0.0
           zn=0.0
           do i=1,nnp
c             write(*,*) i,maxnode,maxnode*2+i
              xyzt(i)=xyz0(i)
              xyzt(maxnode+i)=xyz0(maxnode+i)
              xyzt(maxnode*2+i)=xyz0(maxnode*2+i)
              xn=xn+(xyzt(          i)-xyzt(          1))**2
              yn=yn+(xyzt(maxnode  +i)-xyzt(maxnode  +1))**2
              zn=zn+(xyzt(maxnode*2+i)-xyzt(maxnode*2+1))**2
           enddo
           gsize=( (xn+yn+zn)/neq )**0.5
           write(ilog,*)'@@ gsize: ',gsize
           time=0
           itime0=1
         endif

         dt=deltat

!================ FSI ==========================!

c!===========================================================
c! Communicate with the flow solver

         irankf = 0     ! the processor from the flow side

         call MPI_SEND(deltat, 1, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(is_res, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(irestart, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(npt, 1, MPI_INTEGER,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

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
        read(1,*)j,nmap(i),A0(i),r0(i)
!        if (r0(i) .eq. 1) then
        xyz_s(i)=xyzt(nmap(i))
        xyz_s(nshare+i)=xyzt(maxnode+nmap(i))
        xyz_s(2*nshare+i)=xyzt(2*maxnode+nmap(i))
!        else
!        xyz_s(i)=0.0
!        xyz_s(nshare+i)=0.0
!        xyz_s(2*nshare+i)=0.0
!        endif
        enddo
        close(1)

!!!--------------------------------------------
!!!--- different damping coefficients ---------! Siyuan
        open(unit=1,file='diffdamp.dat')
        allocate(diffdamp(nnp))
        do i=1,nnp
        read(1,*)diffdamp(i)
        enddo
        close(1)
        
        do i=1,nnp
           do j=1,3
               ij=idbc(i,j)
               if (ij .NE. 0) then
                  diffdmm(ij) = diffdamp(i)*0.0
               endif
           enddo
        enddo
!----------------------------------------------

c         !send the initial share coordinates
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
c!=================================================================

!================ end ==========================!
         if (isnap .eq. 0) isnap=10000
ccccc
         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'
c
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
c
c           form lumped mass
            ilump=1
*           write(*,*) maxstiff,maxnode,neq,' ms'
            write(ilog,*)'@@ FORMing lumped mass matrix '
            call zztimer(ilog,'--> FORMmass')
            call formmass_chang(mass(1), npijkm,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      idbc, maxnode,maxelem,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp)
            z1=0.0
            do i=1,neq
               z1=z1+mass(i)
            enddo
            write(iout,*)'Total Mass: ',z1/3
*              write(iout,86) (mass(j),j=1,neq)
            call zztimer(ilog,'<-- FORMmass')
            call zzflush(ilog)
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
                write(iout,86) (grav(j),j=1,neq)
            else
                do i=1,neq
                   grav(i)=0.0  
                enddo
            endif
*           stop
c
c
         rewind(ilod)
         do i=1,neq
            read(ilod) fmag(i)
         enddo
         write(*   ,*) '@@ Reloaded  {P}   OK'
         write(ilog,*) '@@ Reloaded  {P}   OK'
c
c        INPUT load history from a file and interpolate 
         call getload_hist(tforce1,maxforce,npt,deltat)
*        do i=1,npt
*           write(ilog,*)'@@ force: ',tforce1(i,1),tforce1(i,2)
*        enddo
         write(*,*)'                    {P}  | {Grav} |     |     '
         write(*,*)'TYPE force scales:   P1  |  P2    |  P3 | P4  '
         call zzwrt('  -->  ')
         read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
c
c
c
c        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
         write(ilog ,*)'@@  Beginning increment analysis   '
c
c        Initialize  times, disp
         if (is_res .eq. 1) then

         else
           write(*  ,*) '@@ INITIAL disps set to zero'
           do i= 1, neq  
              dispt(i)=0.0
              vel(i)=0.0
!              acc(i)=-dmm*vel(i)           ! Siyuan, different damping coefficients
              acc(i)=-diffdmm(i)*vel(i)
              disptm(i) = dispt(i)-dt*vel(i)+0.5*dt*dt*acc(i)
           enddo

              call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))
         endif
c
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
 71      continue
         rewind(idyn)

         if (is_res .eq. 1) then
!           istart0=itime0+1
           istart0=itime0*isub+2
           iend0=itime0*isub+(npt-1)*isub+1 
         else
         
           k1old=1
           tf1s=tforce1(1,2)*scale_p1
           tf1=tforce1(1,2)
           write(idyn,182) time,tf1,(velout(n), n=1,nout ),tf1
c
           rewind(isnp)
           write(isnp) time, neq, nsnapmax 
           do i= 1, neq  
              write(isnp) dispt(i)
           enddo
           
           istart0=2
           iend0=(npt-1)*isub+1 
         endif
c
         kount = 0
         ksnap=1
         ksnap=0
         inc=1
         inc=iprcnt
         icarriage=0
         call zzwrt(' @@ ')
         
         if (is_res .eq. 0) then
                  call     write_coordTet10(xyzt,dispfult,
     &                             maxelem,maxxyz,0,npijkm,maxnode)
         endif

         ksub=0                  
         dt0=deltat/isub
         pt1=1.0/(2.0*dt0)
         pt2=1.0/dt0**2
         dt=deltat
!         ramp=2.0
         ramp=60.0         ! siyuan to accelerate the vocal fold vibration start
c
c         
c        BIG TIME LOOP
         do 200 itime= istart0, iend0          
            time = real(itime-1)*dt0
            ksub=ksub+1
!            call zzcount_2(itime,inc,icarriage)
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            do i= 1, neq  
               pload(i) =  tf1*(fmag(i)*scale_p1 + grav(i)*scale_p2)
               pload_tmp(i)= pload(i)
            enddo 

130         continue    !Return point for FSI if not converged.

            ! Initialize the incremental displacement and intermediate
            ! geometry.  Prepare the substep marching.
            xyzti=xyzt
            dispfulti=dispfult
            disptmi=disptm
            dispti=dispt

!==================== FSI ===================================!    
c------------------------------------------------------
c Communicate with the flow solver to receive the load
c for each node.  --F.-B. Tian
        
            if(ksub .eq. isub) then

              xyz_f_temp = xyz_f  
              xyz_f = 0.0

         call MPI_SEND((itime-1)/isub, 1, MPI_INTEGER, irankf, 1,
     &                 MPI_COMM_WORLD,istatus,ierr)

         call MPI_RECV(xyz_f(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_RECV(xyz_f(nshare+1  ), nshare,
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_RECV(xyz_f(nshare*2+1), nshare,
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

!============check FSI function====================! siyuan
!      write(*,666)'-->s xload received:',xyz_f(1:10)         
!      write(*,666)'-->s yload received:',xyz_f(nshare+1:nshare+10)
!      write(*,666)'-->s zload received:',xyz_f(nshare*2+1:nshare*2+10)
!=================================================

             resMax_FSIp = 0.0d0
             do i=1,nshare*3
             resP = xyz_f_temp(i)-xyz_f(i)
             xyz_f(i) = xyz_f(i)+ resP*relax_FSIp
             resMax_FSIp  = max(resMax_FSIp,   abs(resP))
             enddo

         call MPI_SEND(resMax_fsip, 1, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

            endif
!=================== end =====================================!
            
             sc   = 1.0d0 - exp (-ramp*time)
!             sc = 1.0d0
             contact_f = 0.0
             zcenter = 0.0d0
             dmin = 0.0024d0
             p0 = 1.0d0
             dclose = 0.0012d0

             do i=1,nshare
               if (r0(i) .eq. 1) then
               ii=nmap(i)

               ij=idbc(ii,1)
               pload(ij)=pload_tmp(ij)+xyz_f(i)*sc
               ij=idbc(ii,2)
               pload(ij)=pload_tmp(ij)+xyz_f(nshare+i)*sc
               ij=idbc(ii,3)
!!!++++++++++++++++++++++++++++++++++++++++  Siyuan
!!!++++++ penalty function for contact ++++

               if (abs(xyzti(2*maxnode+ii)-zcenter) < dmin 
     &       .and.   xyz0(2*maxnode+ii) < zcenter) then
                  deltad=dmin-abs(xyzti(2*maxnode+ii)-zcenter)
                  contact_f(i) = p0*(exp(1500.0d0*deltad)-1.0d0)/36.0d0
                  contact_f(i) = -contact_f(i)*A0(i)
               else if (abs(xyzti(2*maxnode+ii)-zcenter)< dmin 
     &       .and.   xyz0(2*maxnode+ii) > zcenter) then
                  deltad=dmin-abs(xyzti(2*maxnode+ii)-zcenter)
                  contact_f(i) = p0*(exp(1500.0d0*deltad)-1.0d0)/36.0d0
                  contact_f(i) = contact_f(i)*A0(i)
               endif

!!!++++++++++++++++++++++++++++++++++++++++
               pload(ij)=pload_tmp(ij)+xyz_f(nshare*2+i)*sc+contact_f(i)
               else
               endif
             enddo

c
c             forces from body stresses
              call body_stress_3D( dispfulti, iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                    xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadt,iter,maxdim,maxdof
     &                           )
c
              do i= 1, neq  
                 aterm =  pt2*(2*dispti(i)-disptmi(i))
                 vterm =  pt1*disptmi(i)
!                 dload(i) = pload(i) - floadt(i) 
!     &                      + mass(i)*aterm + dmm*mass(i)*vterm 
!                 coeff = pt2 + dmm*pt1
!---------------- Siyuan, different damping coefficients---------
                 dload(i) = pload(i) - floadt(i) 
     &                      + mass(i)*aterm + diffdmm(i)*mass(i)*vterm 
                 coeff = pt2 + diffdmm(i)*pt1
                 disptp(i) = dload(i)/( coeff*mass(i))
              enddo 
*             do i= 1, neq  
*                write(iout,86) pload(i),floadt(i),coeff   
*             enddo 
c
              do i= 1, neq  
                 veli(i)=(disptp(i)-            disptmi(i))*pt1
                 acci(i)=(disptp(i)-2*dispti(i)+disptmi(i))*pt2
              enddo 

!!!---------------------------------------------------------------
!!!-------- kinematic constraints added by Siyuan
              do i=1,nshare
                 if (r0(i) .eq. 1) then
                    ii=nmap(i)
                    ij=idbc(ii,3)
                    if (abs(xyz0(2*maxnode+ii)+disptp(ij)-zcenter)
     &          < dclose .and. xyz0(2*maxnode+ii) < zcenter) then
                     disptp(ij)=zcenter-dclose-xyz0(2*maxnode+ii)
                     veli(ij)=(disptp(ij)-             disptmi(ij))*pt1
                     acci(ij)=(disptp(ij)-2*dispti(ij)+disptmi(ij))*pt2
                    else if (abs(xyz0(2*maxnode+ii)+disptp(ij)-zcenter)
     &          < dclose .and. xyz0(2*maxnode+ii) > zcenter) then
                     disptp(ij)=zcenter+dclose-xyz0(2*maxnode+ii)
                     veli(ij)=(disptp(ij)-             disptmi(ij))*pt1
                     acci(ij)=(disptp(ij)-2*dispti(ij)+disptmi(ij))*pt2
                    endif
                 endif
              enddo
!!!----------------------------------------------------------------

c
c             interchange time subscripts and update
              do i= 1, neq  
                 disptmi(i) = dispti(i)
                 dispti(i)  = disptp(i)
              enddo 
              call update_geom_3D(dispti,dispfulti,maxdim,idbc,maxnode,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    xyzti(1),xyzti(maxnode+1),xyzti(maxnode*2+1))

!====================== FSI ===========================!              
              if (ksub .eq. isub) then
                   vel = veli
                   acc = acci
                   print*, 'Nonstad sending ... '
         ! send the node position
                   do i=1,nshare
                    ! if (r0(i) .eq. 1) then
                     xyz_s(i)=xyzti(nmap(i))
                     xyz_s(nshare+i)=xyzti(maxnode+nmap(i))
                     xyz_s(2*nshare+i)=xyzti(2*maxnode+nmap(i))
                    ! else
                    ! xyz_s(i)=0.0
                    ! xyz_s(nshare+i)=0.0
                    ! xyz_s(2*nshare+i)=0.0
                    ! endif
                   enddo

         call MPI_SEND(xyz_s(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare+1  ), nshare,
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare*2+1), nshare,
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

!         open(unit=400,file='pos.dat')
!         if( ksub.eq.isub) then
!         do i=1,nshare
!         write(400,667)xyz_s(i),xyz_s(nshare+i),xyz_s(nshare*2+i)
!         enddo
!         endif
!         close(400)

!         reassign the velocity
               do i=1,nshare
                 if (r0(i) .eq. 1) then
                 ii=nmap(i)
                 ij=idbc(ii,1)
                 xyz_s(i)=vel(ij)
                 ij=idbc(ii,2)
                 xyz_s(nshare+i)=vel(ij)
                 ij=idbc(ii,3)
                 xyz_s(2*nshare+i)=vel(ij)
                 else
                 xyz_s(i)=0.0
                 xyz_s(nshare+i)=0.0
                 xyz_s(2*nshare+i)=0.0
                 endif
                enddo
              endif
!======================= end ==========================!                   
                icvg = 1
!======================= FSI ===========================!
         ! send the velcocity
              if( ksub.eq.isub) then

         call MPI_SEND(xyz_s(1), nshare, MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare+1  ), nshare,
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

         call MPI_SEND(xyz_s(nshare*2+1), nshare,
     &        MPI_DOUBLE_PRECISION,
     &        irankf, 1, MPI_COMM_WORLD,istatus,ierr)

!         open(unit=500,file='vel.dat')
!         if( ksub.eq.isub) then
!         do i=1,nshare
!         write(500,667)xyz_s(i),xyz_s(nshare+i),xyz_s(nshare*2+i)
!         enddo
!         endif
!         close(500)
         call MPI_RECV(icvg, 1, MPI_INTEGER, irankf, 1,
     &                 MPI_COMM_WORLD,istatus,ierr)

                endif
!======================= end ============================!
! If not converged, go to back to load update
         if(icvg .eq. 0) then
            goto 130
         else if (ksub .eq. isub) then
            write(*,*) 'NONSTAD: FSI iteration done.'
         endif

         xyzt=xyzti
         dispfult=dispfulti
         disptm=disptmi
         dispt=dispti

              if (ksub .eq. isub) then
                   ksub=0
                   kount=kount+1
                   ksnap=ksnap+1
c
c           Print out results of interest for this time loop
                if (iprcnt .eq. kount) then
                  kount = 0
                  nout1=nout-1

                  call    write_coordTet10(xyzt,dispfult,
     &                     maxelem,maxxyz,(itime-1)/isub,npijkm,maxnode)

                  if (isnap .eq. ksnap) then
                     write(isnp) time, neq, nsnapmax 
                     do i= 1, neq  
                        write(isnp) dispt(i)
                     enddo
                     ksnap=0
                  endif
                endif

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
                write(idyn,182) time,tf1s,(velout(n), n=1,nout ),tf1

              endif ! if (ksub == isub)

!======= Writing the restart file================================! siyuan
      if(mod((itime-1)/isub,irestart) .eq. 0) then

        call nond_3D_exp_restart_write(time,(itime-1)/isub,
     &  tf1, tf1s, k1old,
     &  maxnode, maxxyz, xyzt,dispt,disptm, dispfult, vel,acc)

       endif
c
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
c
c
c4200
c     NON-liner DYNamic incremental analysis 
      subroutine non_TL3D_dyn_tet10(stf,geom,wk,pload, disp, fmag,
     &                        dispt , 
     &                   stfeff,mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   iprofv, iprofh, nloc,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt
     &                  ,respceb,isize8
     &         ,isize_dof,dispi,dui,grav,fmag3
     &         ,dloadi,floadi,floadm
     &         ,xyzt,xyzi
     &                      )
c
          implicit real*8 (a-h,o-z)
          include 'commons.std'
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=9000,isize_neq=3000)
         parameter( maxdim=2000,maxdof=maxdim*3)
c
         integer iprofv(neq), iprofh(neq)
     &          ,nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer nloc(neq)
         integer inod(ireadmax,2),mnod(iconelem)
c
         real*8 tforce1(maxforce,2), tforce2(9000,2),tforce3(9000,2)
         real*8 velout(ireadmax)
         real*8 stf(maxstiff),mass(maxstiff),wk(neq), geom(maxstiff)
     &          ,pload(neq ),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 stfeff(maxstiff)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(maxxyz)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 xyzi(maxxyz)
         real*8 dispi(isize_dof),dui(isize_dof),grav(isize_dof)
         real*8 fmag3(isize_dof)
         real*8 dispt(neq)
         real*8 dsumd,dsumz,dtol,pnorm,pnormi,wnorm,wnormi
c
*        real*8 dispfult(900,3), dispfuli(900,3),delu(3000)
         real*8 dispfult(maxdim,3), dispfuli(maxdim,3),delu(maxdof)
         real*8 dloadi(isize_dof),floadi(isize_dof),floadm(isize_dof)
c
c        eigenstuff
         real*8 respceb(isize8)
         real*8 wk0(isize_neq),wk1(isize_neq),wk2(isize_neq)
         real*8 wk_eig(100)
ctag1
c
c         write(ilog,*)'@@  non_TL3D_dyn ver324  August 2011 '
c
c
                 itmp5=125
                 open(unit=itmp5,file='stadyn.tm5')
                 rewind(itmp5)
c                itmp2=127
c                open(unit=itmp2,file='stadyn.tm2')
c                rewind(itmp2)
c                itmp3=128
c                open(unit=itmp3,file='stadyn.tm3')
c                rewind(itmp3)
                 itmp4=129
                 open(unit=itmp4,file='stadyn.tm4',form='unformatted')
                 rewind(itmp4)
c
c
c        algorithm parameters
           write(*,*) maxstiff,maxnode,neq,' ms'
*          stop
         beta0=1.0
         iter_max=60
         dtol=0.0001
         gamma0=1.0
         kref=2
*        iblowup=0
*        idiv=0
c
           xn=0.0
           yn=0.0
           zn=0.0
           do i=1,nnp
c             write(*,*) i,maxnode,maxnode*2+i
              xyzt(i)=xyz0(i)
              xyzt(maxnode+i)=xyz0(maxnode+i)
              xyzt(maxnode*2+i)=xyz0(maxnode*2+i)
              xn=xn+(xyzt(          i)-xyzt(          1))**2
              yn=yn+(xyzt(maxnode  +i)-xyzt(maxnode  +1))**2
              zn=zn+(xyzt(maxnode*2+i)-xyzt(maxnode*2+1))**2
           enddo
           gsize=( (xn+yn+zn)/neq )**0.5
           write(ilog,*)'@@ gsize: ',gsize
c
c
c        Get things ready for time integration loop
c        Get things ready for time integration loop
         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
     &                                        ' snap count '
         call zzwrt('  -->  ')
         read (*,*) deltat,npt,iprcnt,isnap
               if (isnap .eq. 0) isnap=10000
ccccc
         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'
c
         nsnapmax=npt
         write(*,*)' '
*        write(*,*)' suggest     1                  1        5 '
         write(*,*)'    u = u + beta*du   K = Ke + gamma Kg'
         write(*,*)'INPUT:   beta   |  gamma  |  ramp    '
         write(*,*)' suggest    1        1        5 '
         call zzwrt('  -->  ')
         read (*,*) beta0,gamma0,maxramp0  
         write(ilog,*) beta0,gamma0,maxramp0    ,'  ::b g r '
         write(*,*)' '
         write(*,*)'@@ algor 1=full N-R  2=mod N-R '
         write(*,*)'INPUT:   algor | iter max |  tolerance  (<1.0E-5)'
         call zzwrt('  -->  ')
         read (*   ,*) imodify,iter_max,dtol
         imodify=1
         write(ilog,*) imodify,iter_max,dtol,'  ::N-R it max tol '
c
c        set integration constants for alpha-method
         zalpha=-0.04
         zalpha=-0.04
c        zalpha=-0.20
c        zalpha=-0.40
         zalpha= 0.00
         zalpha=-0.04
         zalpha=-0.20
         zalpha=-0.30
         zalpha=-0.50
         zalpha=-0.96   !u
         zalpha=-0.60   !u
         zalpha=-0.55   !u
         zalpha=-0.50
         zalpha=-0.51   !u
         zbeta = 0.25*(1.0-zalpha)**2
         zgamma= 0.50*(1-2*zalpha)
         write(ilog,81)'@@ zalpha beta gamma: ',zalpha,zbeta,zgamma
c
         write(*,*)' '
         write(*,*)'CHOOSE nodal output:'
         nout=0
27       continue
         nout=nout+1
         write(*,*)'TYPE:  node# |  DoF  | rate     <0 0 0 to end>'
         call zzwrt('  -->  ')
         read(*,*) inode,idof,irate
         write(ilog,*) inode,idof,irate,' ::node  dof '
         if (inode .ne. 0) then
             inod(nout,1)=idbc(inode,idof)
             inod(nout,2)=irate
             write(ilog,*) '@@ eqn # ',inod(nout,1)
             goto 27
         endif
         nout=nout-1
         write(*,'(a)')' @@ '
c
         write(*,*)'@@ Eigen types: 0=none  1=      2=subspace '
         write(*,*)'INPUT type |  rate | # vectors '
         call zzwrt('  -->  ')
         read (*,*) ieig_typ,ieig_cnt,ieig_vec 
         ieig_typ=0
         ieig_cnt=1
         ieig_vec=1 
                     if (ieig_cnt .lt. 1) ieig_typ=0
                     eig_max=(npt*ieig_vec*1.0)/ieig_cnt
                     ieig_max=eig_max
                     write(ilog,*)'@@ max eig snp: ',ieig_max,eig_max
         write(ilog,'(1x,3(i8,1x),a)') 
     &                 ieig_typ,ieig_cnt,ieig_vec
     &                                           ,' ::typ c v' 
c
c           form lumped mass
            ilump=1
            write(ilog,*)'@@ FORMing lumped  mass matrix '
            call zztimer(ilog,'--> FORMmass')
            call formmass(mass(1), npijkm,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      idbc, maxnode,maxelem,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
            z1=0.0
            do i=1,neq
               z1=z1+mass(i)
            enddo
            write(ilog,*)'@@Total Mass: ',z1
            write(iout,*)'@@Total Mass: ',z1
c              write(iout,86) (mass(j),j=1,neq)
            call zztimer(ilog,'<-- FORMmass')
            call zzflush(ilog)
c
            if (igravity_on .eq. 1) then
                call body_force_grav( iglobal,
     &                   npijkm,idbc,maxnode,maxelem,
     &                  xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                   prop,nmprop,ippp,
     &                   wk(1)
     &                          )
                z1=0.0
                do i=1,neq
                   grav(i)=wk(i)
                   write(iout,*) i, grav(i)
                   z1=z1+grav(i)
                enddo
                write(iout,*)'Total gravity force: ',z1
*               call gravity( respce8(ir2),respce8(ir3),idbc,maxnode)
            else
                do i=1,neq
                   grav(i)=0.0  
                enddo
            endif
c
c
         rewind(ilod)
         do i=1,neq
            read(ilod) fmag(i)
         enddo
         write(*   ,*) '@@ Reloaded  {P}   OK'
         write(ilog,*) '@@ Reloaded  {P}   OK'
c
c        INPUT load history from a file and interpolate 
         call getload_hist(tforce1,maxforce,npt,deltat)
c
c        only for second load
         if (igravity_on .eq. 1) then
             write(*,*)'@@ Gravity: '
             call getload_hist(tforce2,9000,npt,deltat)
         endif
         if (iload_num .eq. 3) then
             open(unit=itmp,file='stadyn.ld3')
             rewind(itmp)
             do i=1,neq
                read(itmp,*) fmag3(i)
             enddo
             write(*   ,*) '@@ Reloaded  {P}_3   OK'
             write(ilog,*) '@@ Reloaded  {P}_3   OK'
             close(itmp)
c
c            INPUT load history from a file and interpolate 
             call getload_hist(tforce3,9000,npt,deltat)
         endif
c
         write(*,*)'                    reg  | grav | 3rd |     '
         write(*,*)'TYPE force scales:   P1  |  P2  |  P3 | P4  '
         call zzwrt('  -->  ')
         read(*,*) scale_p1,scale_p2,scale_p3,scale_p4
         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
c
c

c/
c        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
         write(ilog ,*)'@@  Beginning increment analysis   '
c
c        Initialize  times, disp, etc
c
           write(*  ,*) '@@ INITIAL disps set to zero'
           write(*  ,*) '@@ '
           do i= 1, neq  
              disp(i)=0.0
              dispt(i)=0.0
              delu(i)=0.0
              vel(i)=0.0
              acc(i)=0.0
           enddo
           do i= 1, nnp  
              do j= 1,3  
                 dispfult(i,j)=0.0
              enddo
           enddo
           force=0.0
           forceout=force
c
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
           itime=0
           time=0
           ziter=0
           k1old=1
           k2old=1
           k3old=1
           tf1=tforce1(1,2)
           tf1s=tf1*scale_p1
           rewind(isnp)
           force=0
c
         kount = 0
         ksnap=1
         ksnap=0
         inc=1
         icarriage=0
         call zzwrt(' @@ ')
c
c
c        prepare BIG INCREMENT LOOP
         itime=0
         iter=0
         tf1=0
         tf2=0
         tf3=0
         deltf1=1
 201     continue
c
c
c        BIG TIME LOOP
!            call zzcount(0,inc,icarriage)
         do 200 itime= 1,npt
!            call zzcount(itime,inc,icarriage)
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
            time = real(itime-1)*deltat
            if (itime .eq. 2) then
                maxramp=20
                maxramp=maxramp0
            else
                maxramp=maxramp0
            endif
            write(ilog,'(a)')'@@ '
            write(iout,'(a)')'@@ '
            write(iout,*)' N: ',itime
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
            tf0=tf1
            tf10=tf1
            tf20=tf2
            tf3=0
            tf30=tf3
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            if (igravity_on .eq. 1) then
                call load_interp(tforce2,9000,k2old,time,tf2)
            endif
            if (iload_num .eq. 3) then
                call load_interp(tforce3,9000,k3old,time,tf3)
            endif
            pnorm=0.0
            do i= 1, neq  
               pload(i) =  tf1*fmag(i)*scale_p1 + tf2*grav(i)*scale_p2
     &                    +tf3*fmag3(i)*scale_p3
               pnorm = pnorm + pload(i)**2
            enddo 
            pnorm = sqrt(pnorm/neq)
c
            if (itime .eq. 1) then
c               modified N-R
*               write(iout,*) ' before iter max ',iter,iter_max
                call formstif_3D(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
c               SECOND assemble the geometric stiffness
                call formgeom_3D(geom(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
                 goto 290
            endif
c           bottom modified N-R
c
c
c           Equilibrium ITERation
c           prepare iterates
            do i=1,neq
               dispi(i)=dispt(i)
c              predictors
               aap(i) = acc(i)
               vvp(i) = vel(i) + acc(i)*dt
               uup(i) = dispt(i) + vel(i)*dt + acc(i)*0.5*dt*dt
            enddo
            do i=1,nnp  
               do j= 1,3  
                  dispfuli(i,j)=dispfult(i,j)
               enddo
            enddo
            do i=1,nnp
               xyzi(i)=   xyzt(i)
               xyzi(maxnode+i)=xyzt(maxnode+i)
               xyzi(maxnode*2+i)=xyzt(maxnode*2+i)
            enddo
c
c
            write(iout,*)' START iteration'
            iter=0
 120        continue
              iter=iter+1
              write(iout,*)'   ' 
              write(iout,*)'@@ iter: ',iter 
              if (iter .gt. iter_max) then 
                  write(iout,*)' NOT CONVERGED at max iter: '
                  goto 100
              endif
c
c
c             full N-R
c
c             forces from body stresses
              call body_stress_3D( dispfuli, iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadi,iter,maxdim,maxdof
     &                           )
c
c             form load increment
              pnormi=0.0
              pnorm =0.0
              do i= 1, neq  
                 aterm = aap(i) + (  1/(zbeta*dt**2))*(dispi(i)-uup(i))
                 vterm = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
                 dpp = -mass(i)*aterm
     &                 - (1+zalpha)*dmm*mass(i)*vterm 
     &                 +    zalpha *dmm*mass(i)*vel(i)
                 ploadm = tf10*fmag(i)*scale_p1 + tf20*grav(i)*scale_p2
     &                   +tf30*fmag3(i)*scale_p3
                 dloadi(i) = (1+zalpha)*pload(i)  - zalpha*ploadm
     &                     - (1+zalpha)*floadi(i) + zalpha*floadm(i)
     &                     + dpp
                 pnormi = pnormi + dloadi(i)**2
                 pnorm  = pnorm  + pload (i)**2
              enddo 
              pnormi=sqrt(pnormi/neq)
              pnorm =sqrt(pnorm /neq)
              write(iout,*) ' force: ',tf1
*             write(iout,*) ' dP norm ',pnormi,pnorm
*             write(iout,86) (dloadi(j),j=1,neq)
c
                  call formstif_3D(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfuli,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
c                 SECOND assemble the geometric stiffness
                  call formgeom_3D(geom(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfuli,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
c                 Form effective stiff matrix by adding geometry matrix
                  do i= 1, neq  
                     iloc=nloc(i)
                     do 24 j= 1, iprofv(i)
                        ii=iloc+j-1
                        stfeff(iloc+j-1) = 
     &                 (1+zalpha)*(stf(iloc+j-1) +gamma0*geom(iloc+j-1))
 24                  continue
                     stfeff(iloc) =  stfeff(iloc)  
     &                     +(1/(zbeta*dt**2))*mass(i)
     &                     +(1+zalpha)*(zgamma/(zbeta*dt))*dmm*mass(i)
                  enddo   
c
c                 Decompose effective stiffness matrix
                  ier1=0
                  call uduCOL_D( stfeff,maxstiff,neq,ierror
     &                          ,iprofv,nloc,z0,z1)
                  write(iout,84) '@@ min/max UDU: ',z0,z1,itime,iter
                  if (ierror.eq.0) then
                      write(*,*)'ERROR: zero diagonal term'
                      return
                  endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c             solve for Delta u
              ier1=0
              call bakCOL_D(stfeff,maxstiff,dloadi, neq, dui,ierror,
     &                                         iprofv,iprofh,nloc)
              dsumd=0.0
              do i=1,neq
                 dsumd=dsumd+ abs(dui(i))
              enddo
              write(ilog,*)'@@  |dui| ',dsumd/neq,iter
c
c             increment Ui
              if (iter .le. maxramp) then
                  beta=iter/(maxramp+0.0)
                  beta=beta*beta0
                zi=2.0**(iter)
                z0=2.0**(maxramp)
                beta=zi/z0
                beta=beta*beta0
              else
                  beta=1.0
                  beta=beta0
              endif
                  do i= 1, neq  
                     dispi(i)= dispi(i) + beta*dui(i)
                  enddo
c
              call update_geom_3D( dispi, dispfuli,maxdim,idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1))
c
c
c             work norm 
              if (pnorm .lt. dtol/10.0) pnorm=dtol/10.0
              write(iout,81) ' dPi: ',pnormi,pnorm,pnormi/pnorm
              wnormi=0.0
              wnorm =0.0
              do i= 1, neq  
                 wnormi = wnormi + dloadi(i)*dui(i)*beta
                 wnorm  = wnorm  + pload (i)*dispi(i)
              enddo 
              if (wnorm .lt. dtol/10.0) wnorm=dtol/10.0
              write(iout,81) ' dWi: ',wnormi/neq,wnorm/neq,wnormi/wnorm
c
c             test for CONvergence
              dsumd=0.0
              dsumz=0.0
              do i=1,neq
                 dsumd=dsumd+ (dui(i)*beta)**2
                 dsumz=dsumz+ dispi(i)**2
              enddo
              dsumd=(dsumd/neq)**0.5
              dsumz=(dsumz/neq)**0.5 + 0*gsize
*             write(iout,*)' dUi: ',dsumd,dsumz
              if (dsumz .lt. dtol/1000.0) dsumz=dtol/1000.0
              dnorm=dsumd/dsumz
              write(iout,*)' dUi: ',dsumd,dsumz
              write(iout,*)' norm: ',dnorm,iter
              write(iout,*)' dtol: ',dtol 
              if (itime .eq. -2) then
                  ptol=0.5e-2
                  ptol=0.25e-2
                  if (pnormi/pnorm .lt. ptol) then
c                     converged, update geom etc, write results
                      write(iout,*)' CONVERGED: '
                      goto 100
                  else
c                     not converged, iterate again
                      goto 120
                  endif
              else
                  if (dnorm .lt. dtol) then
c                     converged, update geom etc, write results
                      write(iout,*)' CONVERGED: '
                      goto 100
                  else
c                     not converged, iterate again
                      goto 120
                  endif
              endif
c
c
c           UPDATE
 100        continue
            write(iout,*)' Maxes: ',iter
c
c           set values for time t+dt -> t
            do i= 1, neq  
               delu(i)=dispi(i)-dispt(i)
               dispt(i)=dispi(i)
               vel(i) = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
               acc(i) = aap(i) + (1.0/(zbeta*dt*dt))*(dispi(i)-uup(i))
            enddo    
            deltf1=tf1-tf0
            do i= 1, nnp  
               do j= 1,3  
                  dispfult(i,j)=dispfuli(i,j)
               enddo
            enddo
            do i=1,nnp
               xyzt(i)          =xyzi(i)
               xyzt(maxnode+i)  =xyzi(maxnode+i)
               xyzt(maxnode*2+i)=xyzi(maxnode*2+i)
            enddo
            do i= 1, neq  
               floadm(i) = floadi(i)
            enddo 
c
c
c           Print out results of interest for this time loop
 290        continue
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
                ziter=iter*1.0
                tf1s=tf1*scale_p1
                tf2s=tf2*scale_p2
                tf3s=tf3*scale_p3
                write(idyn,182) time,tf1s,tf2s,tf3s
     &                              ,(velout(n), n=1,nout)
     &                              ,ziter,tf1,tf2,tf3
                 call zzflush(idyn)
             endif
*            psload=itime*dforce
*            psload=dforce
             if (isnap .eq. ksnap .OR. itime .eq. 1) then
                 write(isnp) time, neq, nsnapmax 
                 do i= 1, neq  
                    write(isnp) dispt(i)
                 enddo
                 ksnap=0
*                call zzflush(idyn)
             endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           eigen analysis
            if (ieig_typ .gt. 0) then
            endif
c           bottom eigenanalysis
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
c
      return
      end
c
