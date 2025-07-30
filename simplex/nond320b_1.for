c
c
c4400
c     NON-liner 3D EXPlicit incremental analysis 
      subroutine non_3D_exp(stf,geom,wk,pload, disp, fmag,
     &                        dispt , 
     &                   stfeff,mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   iprofv, iprofh, nloc,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt)
c
          implicit real*8 (a-h,o-z)
          include 'commons.std'
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=45000)
c
         integer iprofv(neq), iprofh(neq)
     &          ,nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer nloc(neq)
         integer inod(ireadmax,2),mnod(iconelem)
c
         real*8 tforce1(maxforce,2)
         real*8 velout(ireadmax)
         real*8 stf(maxstiff),mass(maxstiff),wk(neq), geom(maxstiff)
     &          ,pload(neq ),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 stfeff(maxstiff)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(maxxyz)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 grav(maxxyz)
         real*8 dispt(neq)
         real*8 dispfult(maxnode,3)
         real*8 disptm(neq),disptp(neq)
c
         real*8 dload(neq),floadt(3*maxnode)
         maxdim=maxnode
         maxdof=3*maxdim
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
          write(ilog,*)'@@  non_3D_exp ver270  November 2007 '
c
c        algorithm parameters
!          write(*,*) maxstiff,maxnode,neq,' ms'
!          stop
         kref=2
c
c
c
c        Get things ready for time integration loop
!         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
!     &                                        ' snap count '
!         call zzwrt('  -->  ')
!         read(ikbd,*) deltat,npt,iprcnt,isnap
         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
     &                       ' snap count | is_res | irestart | isub'
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,iprcnt,isnap, is_res,irestart,isub   ! F.-B. Tian
!-------------------dt,    Nstep, nout, nscreen, res, ires===============

	if(is_res.eq.1)then
	 call nond_3D_exp_restart_read(time,itime0,tf1,tf1s,k1old,  
     & maxnode,maxxyz,nPtsBM,xyzt,dispt,disptm,dispfult,vel,acc,xyz_f)

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
!           write(*,*) maxstiff,maxnode,neq,' ms'
            write(ilog,*)'@@ FORMing lumped mass matrix '
!            call zztimer(ilog,'--> FORMmass')
            call formmass(mass(1), npijkm,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      idbc, maxnode,maxelem,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
            z1=0.0
            do i=1,neq
               z1=z1+mass(i)
            enddo
            write(iout,*)'Total Mass: ',z1/3
*              write(iout,86) (mass(j),j=1,neq)
!            call zztimer(ilog,'<-- FORMmass')
!            call zzflush(ilog)
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
*                write(iout,86) (grav(j),j=1,neq)
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
         write(*,*)'                    {P}  | {Grav} |     |     '
         write(*,*)'TYPE force scales:   P1  |  P2    |  P3 | P4  '
         call zzwrt('  -->  ')
         read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
c
c        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
         write(ilog ,*)'@@  Beginning increment analysis   '
c
c        Initialize  times, disp
	if(is_res.eq.1)then


	else
	     write(*  ,*) '@@ INITIAL disps set to zero'
           do i= 1, neq  
              dispt(i)=0.0
              vel(i)=0.0
              acc(i)=-dmm*vel(i)
              disptm(i) = dispt(i)-dt*vel(i)+0.5*dt*dt*acc(i)
           enddo
*          do i= 1, nnp  
*             do j= 1,3  
*                dispfult(i,j)=0.0
*             enddo
*          enddo
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
 71        continue
           rewind(idyn)

	if(is_res.eq.1)then
	istart0=itime0+1
	iend0=itime0+npt-1
	else
!           itime=0
!           time=0

           k1old=1
           tf1s=tforce1(1,2)*scale_p1
           tf1=tforce1(1,2)
           write(idyn,182) time,tf1,(velout(n), n=1,nout ),tf1
c
           rewind(isnp)
           write(isnp) time, neq, nsnapmax 
 !          do i= 1, neq  
 !             write(isnp) dispt(i)
 !          enddo
c
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
c
c
c        BIG TIME LOOP
         do 200 itime= istart0,iend0
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
!            time = real(itime-1)*deltat
	      dt0=deltat/isub
	      ksub=0
2013        continue
            ksub=ksub+1
         	  time = real(itime-2)*deltat+dt0*ksub
	      pt1=1.0/(2.0*dt0)
	      pt2=1.0/dt0**2
!            call zzcount_2(itime,inc,icarriage)
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            do i= 1, neq  
               pload(i) =  tf1*(fmag(i)*scale_p1 + grav(i)*scale_p2)
            enddo 
c
c             forces from body stresses
              call body_stress_3D( dispfult, iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                    xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadt,iter,maxdim,maxdof
     &                           )
c
              do i= 1, neq  
                 aterm =  (pt2)*(2*dispt(i)-disptm(i))
                 vterm =  (pt1)*disptm(i)
                 dload(i) = pload(i) - floadt(i) 
     &                      + mass(i)*aterm + dmm*mass(i)*vterm 
                 coeff = pt2 + 1.0*dmm*pt1
                 disptp(i) = dload(i)/( coeff*mass(i))
! Fangbao commanded lines by "!" to optimize the code 
!              enddo 
*             do i= 1, neq  
*                write(iout,86) pload(i),floadt(i),coeff   
*             enddo 
c
!              do i= 1, neq  
                 vel(i)=(disptp(i)-           disptm(i))*pt1
                 acc(i)=(disptp(i)-2*dispt(i)+disptm(i))*pt2
!              enddo 
c
c             interchange time subscripts and update
!              do i= 1, neq  
                 disptm(i) = dispt(i)
                 dispt(i)  = disptp(i)
              enddo 
              call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))

	if(ksub.lt.isub)go to 2013
c
c
c           Print out results of interest for this time loop
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
                write(idyn,182) time,tf1s,(velout(n), n=1,nout ),tf1
	         write(*,*)'Itime=',itime
               call write_coordHEX20(xyzt,                ! Fangbao Tian
     &                     npijkm,maxnode,maxelem,maxxyz,itime-1)!

             endif
             if (isnap .eq. ksnap) then
                 write(isnp) time, neq, nsnapmax 
!                 do i= 1, neq  
!                    write(isnp) dispt(i)
!                 enddo
                 ksnap=0
             endif
c
!  Fang-Bao Tian
! Writing the restart file
       if(mod(itime-1,irestart) .eq. 0) then
	 call nond_3D_exp_restart_write(time,itime,tf1,tf1s,k1old,
     & maxnode,maxxyz,nPtsBM,xyzt,dispt,disptm,dispfult,vel,acc,xyz_f)
       endif
             if(itime.eq.2)then
	open(unit=111,file=fdn//'deflection.dat')
	write(111,*)"variables=t y"
	write(111,*)time, xyzt(maxnode+27)
      close(111)
	       else
	open(unit=111,file=fdn//'deflection.dat',position='append')
	write(111,*)time, xyzt(maxnode+27)
      close(111)
	       endif
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

	subroutine nond_3D_exp_restart_write(time,itime,tf1,tf1s,k1old, 
     & maxnode,maxxyz,nPtsBM,xyzt,dispt,disptm,dispfult,vel,acc,xyz_f)

	implicit real*8 (a-h,o-z)
	include 'commons.std'
      parameter( maxdim=20000,maxdof=maxdim*3)
	real*8::xyzt(maxxyz),dispt(neq  ),disptm(neq  ),
     &        vel(neq),acc(neq),dispfult(maxnode,3),
     &        xyz_f(3*nPtsBM)


      ifile = 21 
      open(unit=ifile, file=fdn//'nond_restart_out.dat', 
     &     form='unformatted')

      write(*  ,*) '@@ Write restart file nond_restart_out.dat'    
      write(ifile) time,itime,tf1,tf1s,k1old,gsize
      write(ifile) xyzt(1:maxxyz),
     &             dispt(1:neq),
     &             disptm(1:neq),
     &             dispfult(1:maxnode,1:3),
     &             vel  (1:neq),
     &             acc  (1:neq)

      !Ye,write Force================
      write(ifile) xyz_f(1:3*nPtsBM)
      !END write FOrce

!      do i=1,neq
!      write(ifile)dispt(i),disptm(i),vel(i),acc(i)
!	enddo
!	do i=1,maxxyz
!      write(ifile) xyzt(i)
!	enddo
!	do i=1,maxnode
!      write(ifile)dispfult(i,1),dispfult(i,2),dispfult(i,3)
!	enddo         
      close(ifile)
      return
	end subroutine

	subroutine nond_3D_exp_restart_read(time,itime,tf1,tf1s,k1old, 
     & maxnode,maxxyz,nPtsBM,xyzt,dispt,disptm,dispfult,vel,acc,xyz_f)
	implicit real*8 (a-h,o-z)
	include 'commons.std'
      parameter( maxdim=20000,maxdof=maxdim*3)
	real*8::xyzt(maxxyz),dispt(neq  ),disptm(neq  ),
     &        vel(neq),acc(neq),dispfult(maxnode,3),
     &        xyz_f(3*nPtsBM)

!       parameter(maxnode0=300000, maxxyz0=900000)
!       real*8::xyzt_tmp(maxxyz0),disp_tmp(maxnode0,3)

      ifile = 21
      open(unit=ifile, file=fdn//'nond_restart_in.dat', 
     &     form='unformatted')

      write(*  ,*) '@@ Read restart file nond_restart_in.dat'    
      read(ifile) time,itime,tf1,tf1s,k1old,gsize
      read(ifile) xyzt(1:maxxyz),
     &             dispt(1:neq),
     &             disptm(1:neq),
     &             dispfult(1:maxnode,1:3),
     &             vel  (1:neq),
     &             acc  (1:neq)

      !Ye,read Force================
      read(ifile) xyz_f(1:3*nPtsBM)
      !END read Force===============

c------------------
!      read(ifile) xyzt_tmp(1:maxxyz0),
!     &             dispt(1:neq),
!     &             disptm(1:neq),
!     &             disp_tmp(1:maxnode0,1:3),
!     &             vel  (1:neq),
!     &             acc  (1:neq)
!       print*,'Reassign the varialbes...',maxnode,maxnode0
!       do i=1,nnp
!          xyzt(i) = xyzt_tmp(i)
!          xyzt(maxnode+i)=xyzt_tmp(maxnode0+i)
!          xyzt(maxnode*2+i)=xyzt_tmp(maxnode0*2+i)          
!       enddo
!       do i=1,nnp
!          dispfult(i,1:3) = disp_tmp(i,1:3)
!       enddo
c--------------------

      close(ifile)
      return
	end subroutine
