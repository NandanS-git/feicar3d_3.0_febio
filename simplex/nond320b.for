c
c4240
c     NON-liner DYNamic incremental analysis mult loads 
      subroutine non_TL3D_dyn_mul(stf,geom,wk,pload, disp, fmag,
     &                        dispt , 
     &                   stfeff,mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   iprofv, iprofh, nloc,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt
     &                  ,respce,isize8
     &                  ,tforce_m
     &                       )
c
          implicit real*8 (a-h,o-z)
          include 'commons.std'
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=5000,isize_neq=3000)
         parameter( maxdim=2000,maxdof=maxdim*3)
c
         integer iprofv(neq), iprofh(neq)
     &          ,nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer nloc(neq)
         integer inod(ireadmax,2),mnod(iconelem)
         integer imload(100)
c
         real*8 tforce1(maxforce,2)
         real*8 tforce_m(maxforce,40), tforce(maxforce)
         real*8 velout(ireadmax)
         real*8 stf(maxstiff),mass(maxstiff),wk(neq), geom(maxstiff)
     &          ,pload(neq ),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 stfeff(maxstiff)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(10000)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 xyzi(10000)
         real*8 dispi(3000),dui(3000),grav(3000)
         real*8 dispt(neq)
         real*8 dsumd,dsumz,dtol,pnorm,pnormi,wnorm,wnormi
c
*        real*8 dispfult(900,3), dispfuli(900,3),delu(3000)
*        real*8 dloadi(3000),floadi(3000),floadm(3000)
         real*8 dispfult(maxdim,3), dispfuli(maxdim,3),delu(maxdof)
         real*8 dloadi(maxdof),floadi(maxdof),floadm(maxdof)
c
c        eigenstuff
         real*8 respceb(isize8)
         real*8 wk0(isize_neq),wk1(isize_neq),wk2(isize_neq)
         real*8 wk_eig(100)

ctag1
c
*         write(ilog,*)'@@  non_TL3D_dyn_mul ver250  August 2006 '
          write(ilog,*)'@@  non_TL3D_dyn_mul ver262  November 2006 '
c
*                do i=1,nnp
*                   xyz0(i)=xyz999(i)
*                   xyz0(maxnode+i)=xyz999(maxnode+i)
*                   xyz0(maxnode*2+i)=xyz999(maxnode*2+i)
*                   write(iout,*) i,xyz0(i)
*                enddo
c
                 itmp5=125
                 open(unit=itmp5,file=fdn//'stadyn.tm5')
                 rewind(itmp5)
c                itmp2=127
c                open(unit=itmp2,file='stadyn.tm2')
c                rewind(itmp2)
c                itmp3=128
c                open(unit=itmp3,file='stadyn.tm3')
c                rewind(itmp3)
                 itmp4=129
              open(unit=itmp4,file=fdn//'stadyn.tm4',form='unformatted')
                 rewind(itmp4)
c
c
c        algorithm parameters
           write(*,*) maxstiff,maxnode,neq,' ms'
*          stop
*        beta0=1.0
         iter_max=60
*        dtol=0.0001
*        gamma0=1.0
         kref=2
*        iblowup=0
*        idiv=0
c
           xn=0.0
           yn=0.0
           zn=0.0
	     time=0.0
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
         read(ikbd,*) deltat,npt,iprcnt,isnap
               if (isnap .eq. 0) isnap=10000
ccccc
         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'
c
*        write(*,*)' '
*        write(*,*)'TYPE: force max  | # of incs | max inc'
*        call zzwrt('  -->  ')
*        read (ikbd,*) force_max,npt,npt_max
*        write(ilog,*) force_max,npt,npt_max,'  ::force,npt,max'
*        dforce=force_max/npt
*        dforce_max=dforce
*        dforce_min=dforce/8
*        dforce0=dforce
*        iprcnt=1
         nsnapmax=npt
         write(*,*)' '
         write(*,*)' suggest     1                  1        5 '
         write(*,*)'    u = u + beta*du   K = Ke + gamma Kg'
         write(*,*)'INPUT:   beta   |  gamma  |  ramp    '
         call zzwrt('  -->  ')
         read(ikbd,*) beta0,gamma0,maxramp0  
         write(ilog,*) beta0,gamma0,maxramp0    ,'  ::b g r '
         write(*,*)' '
         write(*,*)'@@ algor 1=full N-R  2=mod N-R '
         write(*,*)'INPUT:   algor | iter max |  tolerance  (<1.0E-5)'
         call zzwrt('  -->  ')
         read(ikbd,*) imodify,iter_max,dtol
         write(ilog,*) imodify,iter_max,dtol,'  ::N-R it max tol '
c
c        set integration constants for alpha-method
         zalpha=-0.04
         zalpha=-0.02
*        zalpha=-0.00
         zalpha= 0.00
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
         write(*,*)'@@ Eigen types: 0=none  1=      2=subspace '
         write(*,*)'INPUT type |  rate | # vectors '
         call zzwrt('  -->  ')
         read(ikbd,*) ieig_typ,ieig_cnt,ieig_vec 
*                    if (ieig_cnt .lt. 1) ieig_cnt=1
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
*           write(*,*) maxstiff,maxnode,neq,' ms'
            write(ilog,*)'@@ FORMing lumped  mass matrix '
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
            write(ilog,*)'@@Total Mass: ',z1
            write(iout,*)'@@Total Mass: ',z1
c              write(iout,86) (mass(j),j=1,neq)
!            call zztimer(ilog,'<-- FORMmass')
!            call zzflush(ilog)
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
*                   write(iout,*) i, grav(i)
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
c        call getload_hist(tforce1,maxforce,npt,deltat)
*        do i=1,npt
*           write(ilog,*)'@@ force: ',tforce1(i,1),tforce1(i,2)
*        enddo
c        write(*,*)'                    reg  | grav |     |     '
c        write(*,*)'TYPE force scales:   P1  |  P2  |  P3 | P4  '
c        call zzwrt('  -->  ')
c        read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
c        write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
c
c        INPUT load history from a file and interpolate 
         write(*,*)'TYPE:  # of loads  | where <style: 1=node 2=xyz>'
         call zzwrt('  -->  ')
         read(ikbd,*) iload,norxyz 
         write(ilog,*) iload,norxyz,' ::# loads style'
         deltat_d=deltat
         do i=1,iload
            call getload_m(time, npt,deltat_d,tforce,maxforce,
     &                       )
            if (norxyz .eq. 1) then
                write(*,*) ' @@ TYPE:   node  |  DoF | P_scale '
                call zzwrt(' @@ --> ')
                read(ikbd,*) jnode,jdof,pscale
                write(*,*)' '
                write(ilog,*) jnode,jdof,pscale,' ::node DoF P_scale'
c
            elseif (norxyz .eq. 2) then
                write(*,*) ' @@ TYPE:   x y z  |  DoF | P_scale '
                call zzwrt(' @@ --> ')
                read(ikbd,*) x1,y1,z1,jdof,pscale
                write(*,*)' '
                write(ilog,184) x1,y1,z1,jdof,pscale,' ::xyz DoF P_s'
c
c               find nearest node
                zlen0=1.0e20
                do n=1,nnp
                   xn=xyz0(          n)
                   yn=xyz0(maxnode  +n)
                   zn=xyz0(maxnode*2+n)
                   zlen = sqrt((x1-xn)**2 + (y1-yn)**2
     &                                         + (z1-zn)**2)
                   if (zlen .lt. zlen0) then
                       zlen0 = zlen
                       jnode = n
                   endif
                enddo
                write(ilog,*)'@@ nearest node: ',jnode
            endif
c
            do j=1,npt
               tforce_m(j,i)=tforce(j)*pscale
            enddo
            idof = idbc(jnode,jdof)
            imload(i)=idof
         enddo
c
c
c        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
         write(ilog ,*)'@@  Beginning increment analysis   '
c
c        Initialize  times, disp, etc
c
           write(*  ,*) '@@ INITIAL disps set to zero'
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
*          tf1=tforce1(1,2)
*          tf1s=tf1*scale_p1
*          write(idyn,182) time,tf1s,(velout(n), n=1,nout),ziter,tf1
           rewind(isnp)
           force=0
*          write(isnp) time, neq, nsnapmax 
*          do i= 1, neq  
*             write(isnp) dispt(i)
*          enddo
c
         kount = 0
         ksnap=1
         ksnap=0
         inc=1
         icarriage=0
         call zzwrt(' @@ ')
c
*          tf1=tforce1(1,2)*scale_p1
*          write(idyn,182) time,tf1,(velout(n), n=1,nout )
*    &                         ,tf2,tf3
c
c        prepare BIG INCREMENT LOOP
*        imodify=1  !full N-R
*        imodify=2  !modified N-R
*        imodify=1  
         itime=0
         iter=0
         tf1=0
         deltf1=1
 201     continue
c
c
c        BIG TIME LOOP
         do 200 itime= 1,npt
            call zzcount_2(itime,inc,icarriage)
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
*220           continue
            write(ilog,'(a)')'@@ '
*            write(iout,'(a)')'@@ '
            write(iout,*)' N: ',itime
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
*           tf0=tf1
*           call load_interp(tforce1,maxforce,k1old,time,tf1)
*           pnorm=0.0
*           do i= 1, neq  
*              pload(i) =  tf1*(fmag(i)*scale_p1 + grav(i)*scale_p2)
*              tf1 = tf1*scale_p1
*              pnorm = pnorm + pload(i)**2
*           enddo 
            do i= 1, neq  
               pload(i)=0.0
            enddo 
            do j=1, iload  
               idof=imload(j)
               pload(idof)=tforce_m(itime,j)
            enddo 
            pnorm=0.0
            do i= 1, neq  
               pnorm = pnorm + pload(i)**2
            enddo 
            pnorm = sqrt(pnorm/neq)
c
            if (imodify .eq. 2) then
c           modified N-R
*           write(iout,*) ' before iter max ',iter,iter_max
            call formstif_3D(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
c           SECOND assemble the geometric stiffness
            call formgeom_3D(geom(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
*           write(iout,*) ' after formgeom iter max ',iter,iter_max
            if (itime .eq. 1) goto 290
c
c           Form effective stiffness matrix by adding geometry matrix
*           do i= 1, neq  
*              iloc=nloc(i)
*              do j= 1, iprofv(i)
*                 stfeff(iloc+j-1) = stf(iloc+j-1)+gamma0*geom(iloc+j-1)
*              enddo   
*           enddo   
                do i= 1, neq  
                   iloc=nloc(i)
                   do 25 j= 1, iprofv(i)
                      ii=iloc+j-1
                       stfeff(iloc+j-1) = 
     &               (1+zalpha)*(stf(iloc+j-1) +gamma0*geom(iloc+j-1))
  25                continue
                    stfeff(iloc) =  stfeff(iloc)  
     &                     +(1/(zbeta*dt**2))*mass(i)
     &                     +(1+zalpha)*(zgamma/(zbeta*dt))*dmm*mass(i)
                enddo   
*           write(iout,*) ' before write STIFF iter max ',iter,iter_max
*           write(iout,*)'STIFFs: diag e+g:'
*           write(iout,86) (stfeff(nloc(j)),j=1,neq)
c
c           Decompose effective stiffness matrix
            ier1=0
            call uduCOL_D(stfeff,maxstiff,neq,ierror,iprofv,nloc,z0,z1)
            write(iout,84) '@@ min/max UDU: ',z0,z1,itime,iter
            if (ierror.eq.0) then
                write(*,*)'ERROR: zero diagonal term'
                return
            endif
            write(iout,*) ' iter max ',iter,iter_max
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if (iter .eq. 1000) then 
c                 estimate disp based on previous rate
                  fract=1.0
                  if (itime .gt. 22222) then
                                     tf0=0
                      fract=0.5*(tf1-tf0)/deltf1
                      write(iout,*)'@@ fract: ',fract,itime
                      write(iout,82) tf1,tf0,deltf1
*                      do i=1,neq
*                         dispi(i)=dispi(i)+delu(i)*fract
*                         write(iout,82) dispi(i)-delu(i)*fract
*     &                                 ,delu(i)*fract
*     &                                 ,dispi(i),fract
*                      enddo
*                   call update_geom_3D( dispi, dispfuli, idbc,maxnode,
*    &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                      xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1))
                  endif
              endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                 full N-R
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
*             if (iter .le. 2) then
*              do i= 1, neq  
*                 pload(i) =  0.5*(tf1+tf0)*fmag(i)
*              enddo 
*           else
*              do i= 1, neq  
*                 pload(i) =  tf1*fmag(i)
*              enddo 
*           endif
              pnormi=0.0
              pnorm =0.0
              do i= 1, neq  
                 aterm = aap(i) + (  1/(zbeta*dt**2))*(dispi(i)-uup(i))
                 vterm = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
                 dpp = -mass(i)*aterm
     &                 - (1+zalpha)*dmm*mass(i)*vterm 
c!!! &                - zalpha*dmm*mass(i)*vel(i)
     &                + zalpha*dmm*mass(i)*vel(i)
c!!!             ploadm = tf0*fmag(i)
*                ploadm = tf0*fmag(i)*scale_p1
                 ploadm = tf0*fmag(i)*pscale1
                 dloadi(i) = (1+zalpha)*pload(i)  - zalpha*ploadm
c!!! &                     - (1+zalpha)*floadi(i) - zalpha*floadm(i)
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
              if (imodify .eq. 1) then 
c                 full N-R
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
*                 write(iout,*)'STIFFs: diag e+g:'
*                 write(iout,86) (stfeff(nloc(j)),j=1,neq)
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
*                 write(iout,*) ' iter max ',iter,iter_max
              endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c             solve for Delta u
              ier1=0
              call bakCOL_D(stfeff,maxstiff,dloadi, neq, dui,ierror,
     &                                         iprofv,iprofh,nloc)
*             write(ilog,*)'@@  dui ',dui(neq),iter
*             write(iout,*) ' dui : '
*             write(iout,86) (dui(j), j=1,neq)
              dsumd=0.0
              do i=1,neq
                 dsumd=dsumd+ abs(dui(i))
              enddo
              write(ilog,*)'@@  |dui| ',dsumd/neq,iter
c
c             increment Ui
*             maxramp=10 
              if (iter .le. maxramp) then
*                 beta=1.0/(maxramp+1-iter)
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
*             write(iout,*) ' ui : ',beta
*             write(iout,86) (dispi(j), j=1,neq)
c
              call update_geom_3D( dispi, dispfuli,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1))
c
*                 write(iout,*) ' updated [ui] : '
*                 do i= 1, nnp  
*                  write(iout,82) (dispfuli(i,j), j=1,3) 
*                 enddo
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
              dsumz=(dsumz/neq)**0.5 + gsize
*             write(iout,*)' dUi: ',dsumd,dsumz
*             if (dsumz .lt. dtol/1000.0) dsumz=dtol/1000.0
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
*           iblowup=0
*           idiv   =0
c
c           set values for time t+dt -> t
            do i= 1, neq  
               delu(i)=dispi(i)-dispt(i)
               dispt(i)=dispi(i)
*              write(ilog,*) i,dispt(i),dispi(i)
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
*               tf1s=tf1*scale_p1
*               write(idyn,182) time,tf1s,(velout(n), n=1,nout)
                write(idyn,182) time,(tforce_m(itime,j),j=1,iload)
     &                              ,(velout(n), n=1,nout)
     &                              ,ziter,tf1
             endif
*            psload=itime*dforce
*            psload=dforce
             if (isnap .eq. ksnap .OR. itime .eq. 1) then
                 write(isnp) time, neq, nsnapmax 
                 do i= 1, neq  
                    write(isnp) dispt(i)
                 enddo
                 ksnap=0
!                 call zzflush(idyn)
                call write_coord(xyzt,                ! Fangbao Tian
     &                     npijkm,maxnode,maxxyz,itime-1)!
             endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           eigen analysis
            if (ieig_typ .gt. 0) then
                if (mod(itime+1,ieig_cnt) .eq. 0 .OR. itime .eq. 1
     &                                       .OR. itime .eq. npt) then
c
c                   monitors
c                   compute kinetic energy, etc 
                    s_energy=0.0
                    t_energy=0.0
                    do i=1,neq
                       t_energy = t_energy + 0.5*mass(i)*vel(i)**2
                    enddo
                    p_energy=0.0
                    do i=1,neq
                       p_energy = p_energy + 0.5*pload(i)*dispt(i)
                    enddo
c
                        xc=0.0
                        yc=0.0
                        zc=0.0
                        do i=1,nnp
                           xc = xc + xyzt(i)   
                           yc = yc + xyzt(maxnode+i)
                           zc = zc + xyzt(maxnode*2+i)
                        enddo
                        xc=xc/nnp
                        yc=yc/nnp
                        zc=zc/nnp
                        gsize=0.0
                        do i=1,nnp
                           gsize = gsize + (xyzt(i) - xc)**2   
     &                                 + (xyzt(maxnode+i) - yc)**2
     &                                 + (xyzt(maxnode*2+i) - zc)**2
                        enddo
                        gsize=sqrt(gsize/nnp)
                        t_disp=0.0
                        do i=1,neq
                           t_disp = t_disp + dispt(i)**2
                        enddo
                        u_aver=sqrt(t_disp/neq)
*                   write(itmp3,85) time,tf3,(sens(j),j=1,min(40,neq))
*    &                         ,size,u_aver
c
c                   reconstruct tangent stiffness
                    do i= 1, neq  
                       iloc=nloc(i)
                       do j= 1, iprofv(i)
                          stfeff(iloc+j-1) = stf(iloc+j-1)
     &                                     + gamma0*geom(iloc+j-1)
                       enddo
                    enddo
c
                    if (ieig_typ .eq. 0) then
c                       none
                    elseif (ieig_typ .eq. 1) then
c                       vector iter
                    elseif (ieig_typ .eq. 2) then
c                       subspace iteration for multiple eigenvalues
c
                        call sub_spc(stfeff,isize8,dloadi,dui, 
     &                           mass,     iprofv,iprofh,nloc,
     &                           respceb,ieig_loc,wk_eig)
c
c                       write(itmp5,85) time,(wk_eig(j),j=1,10)
c    &                                  ,tf1,tf2,tf3
c    &                             ,s_energy,t_energy,p_energy
                        zlam=wk_eig(1)
                        do k=1,ieig_vec
                           do j=1,neq
                              wk0(j) = respceb(ieig_loc+neq*(k-1) -1+j)
                           enddo
c                          normalize shape
                           wkmax0 = 0.0
                           do i= 1, neq  
                              if (abs(wk0(i)).gt. wkmax0) 
     &                                wkmax0=abs(wk0(i))
                           enddo
                           do i= 1, neq  
                              wk0(i) = wk0(i)/wkmax0
                           enddo
**                         alpha1=1
**                         alpha2=0
**                         do i= 1, neq  
**                            wk(i) = alpha1*wk1(i) + alpha2*wk(i) 
**                         enddo
**                         write(itmp4) time, neq,ieig_max 
**                         do i= 1, neq  
**                            write(itmp4) wk1(i)
**                         enddo
                           write(itmp4) time, neq, ieig_max 
                           do i= 1, neq  
                              write(itmp4) wk0(i)
                           enddo
c
c                          keep shape for ping
                           if (k .eq. 1) then
                               do i= 1, neq  
                                  wk1(i) = wk0(i) 
                               enddo
                           endif
                           if (k .eq. 2) then
                               do i= 1, neq  
                                  wk2(i) = wk0(i) 
                               enddo
                           endif
                        enddo
                    endif
c                   bottom eigen type
c
                tf1s=tf1*scale_p1
*                   write(itmp5,85) time,tf1s,(wk_eig(j),j=1,ieig_vec)
                    write(imon ,85) time,tf1s,(wk_eig(j),j=1,ieig_vec)
*    &                                  ,tf1,tf2,tf3
     &                             ,s_energy,t_energy,p_energy
     &                         ,gsize,u_aver
!                 call zzflush(imon)
ctag1
                endif
            endif
c           bottom eigenanalysis
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            go to load increment
*            goto 201
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
 184       format(1x,3(g12.6,1x),i5,1x,g12.6,1x,a)
c
      return
      end
c-------------------
c
c  GET applied transient LOAD Multiple history and interpolate
c   
c  Only reads in the frequency and phase delay of the load. --Fangbao Tian
c--------------------
      subroutine getload_m(time, npt,dt,tforce,nmax
     &                      )
         implicit real*8 (a-h,o-z)
	   include 'commons.std'
         real*8 tforce(nmax)
         real*8 time,ddt
         real*8 wk0(100)
         character*40 fyl2
c
c        read data file
	PI = 4.0d0*atan(1.0d0)

         call zzwrt(' @@ TYPE:  Force_Filename --> ')
         read(ikbd,'(1a40)') fyl2
         write(*,*)' '
         write(ilog,'(a3)') '@@ force '
         write(ilog,*) fyl2
c
         write(*,*) ' @@ TYPE:   time col  |  force col | max cols  '
         call zzwrt(' @@ --> ')
         read(ikbd,*) jcolt,jcolf,jcol
         write(*,*)' '
         write(ilog,*) jcolt,jcolf,jcol,' ::col t f m '
c
         open(unit=itmp ,file=fdn//adjustl(fyl2))
         rewind itmp 
c
c        adjust times and interpolate
c 
         read(itmp, *, end=240) iflag, thetam,freq, phase, ramp,
     &                    alpha0, alpham , phase1, phase2,ramp1
         write(ilog,*) ' @@ load function:   Freq  |  Phase  '

         do i = 1,npt
           ! flapping force
           tt       = time + (i-1)*dt
           sc       = 1.0d0 - exp(-ramp*tt)       ! load ramp factor
           tforce(i)= thetam*sin(2.0d0*PI*freq*tt + phase) * sc

           if(iflag .eq. 1) then  ! twisting force
          !   theta = thetam*sin(2.0d0*PI*freq*tt + phase) * sc
             sc1   = 1.0d0 - exp(-ramp1*tt)      ! load ramp factor
             tforce(i)=(alpha0+alpham*sin(2*pi*freq*tt+phase1))
     &                *sin(tforce(i)+phase2)*sc1
             !write(*,'(2F12.5)') tt,  tforce(i)
            endif
         enddo

         close (itmp )
         return
c 
 240     continue
         do i = 1,npt
            tforce(i)=0.0
         enddo
         write(ilog,*)'@@ NO DATA in load file'
c
      return
      end
c
c
c4420?
c     NON-liner STATIC incremental analysis 
      subroutine non_gen_exp(stresst_ip,floadt,wk,ploadt, disp, fmag,
     &                        dispt , 
     &                   mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   iprofv, iprofh, nloc,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt
     &                   )
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
         real*8 mass(maxstiff),wk(neq)
     &          ,ploadt(neq ),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(10000)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 grav(3000)
         real*8 dispt(neq)
         real*8 dispfult(900,3)
         real*8 disptm(3000),disptp(3000)
c
*        real*8 dload(3000)
         real*8 floadt(neq)
         real*8 stresst_ip(nel,6,27)
c
             do i=1,nel
             do j=1,6  
             do k=1,27  
                stresst_ip(i,j,k)=0.0
             enddo
             enddo
             enddo
c
*                do i=1,nnp
*                   xyz0(i)=xyz999(i)
*                   xyz0(maxnode+i)=xyz999(maxnode+i)
*                   xyz0(maxnode*2+i)=xyz999(maxnode*2+i)
*                   write(iout,*) i,xyz0(i)
*                enddo
c
c        algorithm parameters
           write(*,*) maxstiff,maxnode,neq,' ms'
*          stop
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
         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
     &                                        ' snap count '
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,iprcnt,isnap
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
         write(*,*)'TYPE:  node# |  DoF  | rate       <0 0 to end>'
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
c           form mass
            ilump=1
            write(ilog,*)'@@ FORM mass matrix '
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
            write(iout,*)'Total (active) Mass: ',z1
*               write(iout,86) (mass(j),j=1,neq)
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
*        do i=1,npt
*           write(ilog,*)'@@ force: ',tforce1(i,1),tforce1(i,2)
*        enddo
         write(*,*)'TYPE force scales:   P1  |  P2  |  P3 | P4  '
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
           write(*  ,*) '@@ INITIAL disps set to zero'
           do i= 1, neq  
              dispt(i)=0.0
              vel(i)=0.0
              acc(i)=-dmm*vel(i)
              disptm(i) = dispt(i)-dt*vel(i)+0.5*dt*dt*acc(i)
           enddo
*          stop
*          do i= 1, nnp  
*             do j= 1,3  
*                dispfult(i,j)=0.0
*             enddo
*          enddo
              call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))
*          stop
c
c          output initial dyn and snap
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
           k1old=1
           tf1=tforce1(1,2)*scale_p1
           write(idyn,182) time,tf1,(velout(n), n=1,nout )
c
           rewind(isnp)
           write(isnp) time, neq, nsnapmax 
           do i= 1, neq  
              write(isnp) dispt(i)
           enddo
c
c
         kount = 0
         ksnap=1
         ksnap=0
         inc=1
         inc=iprcnt
         icarriage=0
         call zzwrt(' @@ ')
*          stop
c
c
c        BIG TIME LOOP
         do 200 itime= 2,npt
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
            time = real(itime-1)*deltat
            call zzcount_2(itime,inc,icarriage)
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            do i= 1, neq  
               ploadt(i) =  tf1*(fmag(i)*scale_p1 + grav(i)*scale_p2)
            enddo 
*          stop
c
c             forces from body stresses
              call body_stress_gen( dispfult, iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadt,iter,
     &                    stresst_ip
     &                           )
*             stop
c
              do i= 1, neq  
                 aterm =  (1.0/(dt**2))*(2*dispt(i)-disptm(i))
                 vterm =  (1.0/(2*dt))*disptm(i)
                 dpp      = ploadt(i) - floadt(i) 
     &                      + mass(i)*aterm + dmm*mass(i)*vterm 
                 coeff = 1.0/(dt**2) + 1.0*dmm/(2.0*dt)
                 disptp(i) = dpp/( coeff*mass(i))
              enddo 
*             write(iout,*)' PFu ',n
*             do i= 1, neq  
*                write(iout,86) ploadt(i),floadt(i),disptp(i)
*             enddo 
c
              do i= 1, neq  
                 vel(i)=(disptp(i)-           disptm(i))/(2*dt)
                 acc(i)=(disptp(i)-2*dispt(i)+disptm(i))/(dt**2)
              enddo 
c
c             interchange time subscripts and update
              do i= 1, neq  
                 disptm(i) = dispt(i)
                 dispt(i)  = disptp(i)
              enddo 
              call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))
*             write(iout,*)' dispful ',n
*             do i= 1, nnp  
*                write(iout,86) (dispfult(i,j), j=1,3)  
*             enddo 
c
              call update_stress ( 
     &                   stresst_ip,
     &                   dispfult, 
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                   prop,nmprop,ippp
     &                )
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
                write(idyn,182) time,tf1,(velout(n), n=1,nout )
             endif
             if (isnap .eq. ksnap) then
                 write(isnp) time, neq, nsnapmax 
                 do i= 1, neq  
                    write(isnp) dispt(i)
                 enddo
                 ksnap=0
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
c     Nodal loads due to body stresses for 3D solids
      subroutine body_stress_gen( dispful, iglobal,
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xord0,yord0,zord0,
     &                                    xord,yord,zord,
     &                   prop,nmprop,ippp,
     &                   floadt,iter,stresst_ip
     &                )
c
          implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3)
         real*8  dispful(900,3)
         real*8  xord0(nnp), yord0(nnp), zord0(nnp)
         real*8  xord(nnp), yord(nnp), zord(nnp)
c
         real*8 prop(ippp,10),dd(6,6)
         integer nmprop(nel)
c
         real*8  stress(6,27), force(60),force20(3,20)
         real*8  floadt(neq)
         real*8  uvw(3,20),xyz(3,20),ss(6)
         real*8  stresst_ip(nel,6,27)
c
c
*        write(iout,*)' << in body_stress >>'
c
         do i=1,neq
            floadt(i)=0.0
         enddo
c
*        rewind(igeo)
c
c          For each element, calculate the strain, stress at centroid
           do 50 n=1,nel
              mat=nmprop(n)
              neltype=npijkm(n,1)
c
c                 plate
*                 e0=prop(mat,1)
*                 g0=prop(mat,2)
*                 t0=prop(mat,3)
*                 r0=prop(mat,4)
*                pl0=prop(mat,5)
c
*               call dmat(e0,g0,dd)
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
                  enddo
                  do i=1,6
                     do k=1,27
                        stress(i,k)=stresst_ip(n,i,k)
                     enddo
                  enddo
*          write(iout,*)' [stress 1] '
*          do i=1,6
*             write(iout,82) (stress(i,j), j=1,27)
*          enddo
c
                  ired_int=3
                  call int_SBE_Hex20(neltype,ired_int, 
     &                            xyz,uvw,
     &                              stress,
     &                              force20)
*                 do i=1,3
*                    write(iout,82) (force20(i,k), k=1,6)
*                 enddo
c
                  call assembFOR_HEX(floadt,force20,idbc,maxnode
     &                          ,npijkm,maxelem,n)
c
              endif
 50        continue
c          end of loop over elements
c
 82        format(1x,40(g12.6,1x))
 83        format(1x,6(g12.6,1x))
 86        format(1x,6(g12.6,1x))
c
 999  continue
      return
      end
c
c
c     INTegrate Stress [BE] dV  for HEX 20 elem
      subroutine int_SBE_HEX20(neltype,ired_int, 
     &                           xyz,uvw, stress,force
     &                        )
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),uvw(3,20)
         real*8 force(3,20),ff(60),ss(6),stress(6,27)
c
         real*8 BE(6,60),Bd(9,60)
         real*8 xg(4,4),wgt(4,4)
c        integer*4  icoord(8,3)
c        gauss-legendre sampling points
         data xg/ 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     &       -.5773502691896d0,  .5773502691896d0, 0.0d0, 0.0d0, 
     &       -.7745966692415d0, 0.0d0, .7745966692415d0,0.0d0, 
     &       -.8611363115941d0,
     &       -.3399810435849d0,  .3399810435849d0,  .8611363115941d0 /
c        gauss-legendre weights
         data wgt/ 2.0d0, 0.0d0, 0.0d0, 0.0d0, 
     &             1.0d0, 1.0d0, 0.0d0, 0.0d0, 
     &    .5555555555556d0, .8888888888889d0, .5555555555556d0, 0.0d0, 
     &    .3478548451375d0, .6521451548625d0,
     &    .6521451548625d0, .3478548451375d0 /
c
c        data icoord/1,2,2,1,1,2,2,1
c    &              ,1,1,2,2,1,1,2,2
c    &              ,1,1,1,1,2,2,2,2 /
c
         iout=23
         kmax=neltype-1
c
*          write(iout,*)' [x y z] '
*          do i=1,20
*             write(iout,82) (xyz(j,i), j=1,3)
*          enddo
*          write(*,*) ' '
*          write(iout,*)' [stress] '
*          do i=1,6
*             write(iout,82) (stress(i,j), j=1,27)
*          enddo
c
         do i=1,kmax*3
            ff(i) = 0.0 
         enddo
         ip_n=0
         nint=ired_int
c   
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
         ip_n = ip_n + 1
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
c
c                 deriv oper and jacobian
                  call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
c
c                 FORCES
*                 write(iout,*) ' BE ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (BE(i,j), j=1,24)
*                 enddo
c
                  do i=1,6
                     ss(i)=stress(i,ip_n)
                  enddo
*                 write(iout,*) ' ss '
*                    write(iout,82) (ss(j), j=1,6)
c                 add contib to force
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
                  do 70 j=1,kmax*3
                     sum = 0.0
                     do k=1,6
                        sum = sum + ss(k)*BE(k,j)
                     enddo
                     ff(j) = ff(j) + sum*wt
 70               continue
*                 write(iout,*) ' ff '
*                    write(iout,82) (ff(j), j=1,60)
 80      continue
c
c        assign to f[3,20]
         nn=0
         do j=1,kmax
            do i=1,3
               nn=nn+1
               force(i,j)=ff(nn)
            enddo
         enddo
*                 write(iout,*) ' force '
*                 do i=1,3
*                    write(iout,82) (force(i,j), j=1,20)
*                 enddo
c        check equil
c        fx=0
c        fy=0
c        fz=0
c        do j=1,kmax
c           fx=fx+force(1,j)
c           fy=fy+force(2,j)
c           fz=fz+force(3,j)
c        enddo
c                   write(iout,* ) ' equil ? '    
c                   write(iout,82) fx,fy,fz
c
 82        format(1x,70(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
c
c     UPDATE STRESS etc 
      subroutine update_stress( 
     &                   stresst_ip,
     &                   dispful, 
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xord0,yord0,zord0,
*    &                                    xord,yord,zord,
     &                   prop,nmprop,ippp
     &                )
c
          implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3)
         real*8  xord0(nnp), yord0(nnp), zord0(nnp)
*        real*8  xord(nnp), yord(nnp), zord(nnp)
         real*8  dispful(900,3)
c
         real*8 prop(ippp,10),dd(6,6)
         integer nmprop(nel)
c
*        real*8  stress(6,27)
         real*8  uvw(3,20),xyz(3,20)
         real*8  stresst_ip(nel,6,27)
*        real*8  duful(900,3),duvw(3,20)
ccccccccccccccc
c
         real*8 uu(60),ee(6),ud(9),ss(6)
         real*8 BE(6,60),Bd(9,60)
         real*8 xg(4,4),wgt(4,4)
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
*        write(iout,*)' << in update_plastic >>'
c
c          For each element, calculate the strain, stress at centroid
           do 50 n=1,nel
              mat=nmprop(n)
              neltype=npijkm(n,1)
c
c                 plate
                  e0=prop(mat,1)
                  g0=prop(mat,2)
                  t0=prop(mat,3)
                  r0=prop(mat,4)
                 pl0=prop(mat,5)
*               call dmat(e0,g0,dd)
c
*                 if (mat .eq. 1) then
*                     hprime0=etan1/(1.0 - etan1/e0)
*                     sy0=sy1
*                 elseif (mat .eq. 2) then
*                     hprime0=etan2/(1.0 - etan2/e0)
*                     sy0=sy2
*                 endif
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
*                    write(iout,* ) ' uu '         
*                    write(iout,82) (uu(k), k=1,60)
c
                  ired_int=3
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
c
               iconstit=2
                     if (iconstit .eq. 1) then
                         call dmat(e0,g0,dd)
                         call elast_lin(ud,dd,ss)
                     elseif (iconstit .eq. 2) then
                         a10=e0
                         a01=g0
                         call elast_rub(ud,e0,g0,ss)
                     endif
*                    write(iout,* ) ' ES '         
*                    write(iout,82) (ud(k), k=1,9)
*                    write(iout,82) (ee(k), k=1,6)
*                    write(iout,82) (ss(k), k=1,6)
c
                     do i=1,6
                        stresst_ip(n,i,ip_n)=ss(i)
                     enddo
 80             continue
c               bottom loop over IP
c
              endif
c             bottom diff elems
 50        continue
c          end of loop over elements
c
 82        format(1x,40(g12.6,1x))
 83        format(1x,6(g12.6,1x))
 86        format(1x,6(g12.6,1x))
 866       format(1x,6(g12.6,1x),a)
 182       format(1x,40(g12.6,1x))
c
 999  continue
      return
      end
c
c
      subroutine elast_lin(ud,dd,ss)
          implicit real*8 (a-h,o-z)
          real*8 ud(9),ss(6),dd(6,6),ee(6)
c
c                    strain ref to zero 
                     z8=1.0
                     ee(1)= ud(1) 
     &                    + z8*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
                     ee(2)= ud(5) 
     &                    + z8*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
                     ee(3)= ud(9) 
     &                    + z8*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
                     ee(4)= ud(2)+ud(4) 
     &                    + z8*(ud(1)*ud(2)+ud(4)*ud(5) +ud(7)*ud(8))
                     ee(5)= ud(6)+ud(8) 
     &                    + z8*(ud(2)*ud(3)+ud(5)*ud(6) +ud(8)*ud(9))
                     ee(6)= ud(3)+ud(7) 
     &                    + z8*(ud(1)*ud(3)+ud(4)*ud(6) +ud(7)*ud(9))
c
c                    assume stress is linear elastic
                     do i=1,6
                        sum = 0.0  
                        do k=1,6
                           sum = sum + dd(i,k)*ee(k)
                        enddo
                        ss(i)=sum
                     enddo
      return
      end
c
c
c
      subroutine elast_rub(ud,e0,g0,ss)
          implicit real*8 (a-h,o-z)
          real*8 ud(9),ss(6),dd(6,6),def(3,3),cc(3,3),del(3,3),ssk(3,3)
          real*8 cci(3,3)
c
c                 assume mooney-rivlin relation  
                c00=e0
                gam=g0
                a10 = c00*(1-gam)
                a01 = c00*(  gam)
                zmu=2*c00
                zk0 = 19.666667*zmu
c
                   do i=1,3
                      do j=1,3
                         del(i,j)=0.0  
                      enddo
                   enddo
                         del(1,1)=1.0    
                         del(2,2)=1.0    
                         del(3,3)=1.0    
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
*                  write(ilog,*)'@@ J: ', zjj
c
c                    form [C] and invatiants
                     do i=1,3
                        do j=1,3
                           sum = 0.0  
                           do k=1,3
                              sum = sum + def(k,i)*def(k,j)
                           enddo
                           cc(i,j)=sum
                        enddo
                     enddo
c                    test inverse
*                    nn=0
*                    do i=1,3
*                       do j=1,3
*                          nn=nn+1
*                          cc(i,j)=nn  
*                       enddo
*                    enddo
*                    cc(1,1)=2
*                    cc(2,1)=cc(1,2)
*                    cc(3,1)=cc(1,3)
*                    cc(3,2)=cc(2,3)
*                    write(*,*) ' '
*                    do i=1,3
*                       write(*,*) (cc(i,j), j=1,3)
*                    enddo
                     call inv3x3(cc,cci)
*                    write(*,*) ' '
*                    do i=1,3
*                       write(*,*) (cci(i,j), j=1,3)
*                    enddo
*                    write(*,*) ' '
*                    do i=1,3
*                       do j=1,3
*                          sum = 0.0  
*                          do k=1,3
*                             sum = sum + cci(i,k)*cc(k,j)
*                          enddo
*                          write(*,*) sum,i,j
*                       enddo
*                    enddo
*                    stop
c
               zi1 = cc(1,1) + cc(2,2) + cc(3,3)
               zi2 = cc(1,1)*cc(2,2) +cc(2,2)*cc(3,3) +cc(3,3)*cc(1,1)
     &             - cc(1,2)**2      -cc(2,3)**2      -cc(1,3)**2     
               zi3 = cc(1,1)*cc(2,2)*cc(3,3) +2*cc(1,2)*cc(2,3)*cc(1,3)
     &             - cc(1,1)*cc(2,3)**2 -cc(2,2)*cc(1,3)**2
     &             - cc(3,3)*cc(1,2)**2     
               zj3 = sqrt(zi3)
*              zj1 = zi1/(zi3**(1.0/3.0))
*              zj2 = zi2/(zi3**(2.0/3.0))
*              zj3 = zi3**(1.0/2.0)
c
               zi3m3 = 1.0/(zi3**(1.0/3.0)) 
               b1 = 2*a10*zi3m3 + 2*a01*zi1*zi3m3**2
               b2 =             - 2*a01    *zi3m3**2
               b3 = -( 2*a10*zi1*zi3m3 + 4*a01*zi2*zi3m3**2 )/3.0
c
                  do i=1,3
                     do j=1,3
                        ssk(i,j) = b1*del(i,j) +b2*cc(i,j) +b3*cci(i,j)
     &                            + zk0*(zj3-1)*zj3*cci(i,j)
                     enddo
                  enddo
                  ss(1)=ssk(1,1)
                  ss(2)=ssk(2,2)
                  ss(3)=ssk(3,3)
                  ss(4)=ssk(1,2)
                  ss(5)=ssk(1,3)
                  ss(6)=ssk(2,3)
      return
      end
c
c
      subroutine inv3x3(cc,cci)
          implicit real*8 (a-h,o-z)
          real*8 cc(3,3),cci(3,3),dd(3,3)
c
                d00 = cc(1,1)*(cc(2,2)*cc(3,3)-cc(3,2)*cc(2,3))
     &               -cc(1,2)*(cc(2,1)*cc(3,3)-cc(3,1)*cc(2,3))
     &               +cc(1,3)*(cc(2,1)*cc(3,2)-cc(3,1)*cc(2,2))
*           write(*,*) ' '
*           write(*,*) d00
                dd(1,1) =  cc(2,2)*cc(3,3)-cc(3,2)*cc(2,3)
                dd(1,2) = -cc(2,1)*cc(3,3)+cc(3,1)*cc(2,3)
                dd(1,3) =  cc(2,1)*cc(3,2)-cc(3,1)*cc(2,2)
c
                dd(2,1) = -cc(1,2)*cc(3,3)+cc(3,2)*cc(1,3)
                dd(2,2) =  cc(1,1)*cc(3,3)-cc(3,1)*cc(1,3)
                dd(2,3) = -cc(1,1)*cc(3,2)+cc(3,1)*cc(1,2)
c
                dd(3,1) =  cc(1,2)*cc(2,3)-cc(2,2)*cc(1,3)
                dd(3,2) = -cc(1,1)*cc(2,3)+cc(2,1)*cc(1,3)
                dd(3,3) =  cc(1,1)*cc(2,2)-cc(2,1)*cc(1,2)
*                    write(*,*) '[d] '
*                    do i=1,3
*                       write(*,*) (dd(i,j), j=1,3)
*                    enddo
c
                do i=1,3
                   do j=1,3
                      cci(i,j) = dd(j,i)/d00
                   enddo
                enddo
      return
      end
c
c
c4200
c     NON-liner DYNamic incremental analysis 
!FB
      subroutine non_TL3D_dyn(stf,geom,wk,pload, disp, fmag,
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
         parameter( maxdim=20000,maxdof=maxdim*3)
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
*         write(ilog,*)'@@  non_TL3D_dyn ver240  february 2006 '
*         write(ilog,*)'@@  non_TL3D_dyn ver248  April 2006 '
!Fangbao          write(ilog,*)'@@  non_TL3D_dyn ver320  Ocyober 2010 '
c
*                do i=1,nnp
*                   xyz0(i)=xyz999(i)
*                   xyz0(maxnode+i)=xyz999(maxnode+i)
*                   xyz0(maxnode*2+i)=xyz999(maxnode*2+i)
*                   write(iout,*) i,xyz0(i)
*                enddo
c
                 itmp5=125
                 open(unit=itmp5,file=fdn//'stadyn.tm5')
                 rewind(itmp5)
c                itmp2=127
c                open(unit=itmp2,file='stadyn.tm2')
c                rewind(itmp2)
c                itmp3=128
c                open(unit=itmp3,file='stadyn.tm3')
c                rewind(itmp3)
                 itmp4=129
              open(unit=itmp4,file=fdn//'stadyn.tm4',form='unformatted')
                 rewind(itmp4)
c
c
c        algorithm parametersl
           write(*,*) maxstiff,maxnode,neq,' ms'
! comments by Fangbao
! maxstiff is the element of matrix
! maxnode is the maxnode allowed in the simulation
! neq is the number of freedom, if a 2D problem has np nodes, neq=2.0*np, for a 3D neq=2.0*np  
*          stop
         beta0=1.0
         iter_max=60
         dtol=0.0001
         gamma0=1.0
         kref=2
*        iblowup=0
*        idiv=0
c
c
c
c        Get things ready for time integration loop
c        Get things ready for time integration loop
         write(*,*)'TYPE:  time inc  |  # of incs  |  print count |',
     &                       ' snap count | is_res | irestart'
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,  iprcnt,isnap, is_res,irestart   ! F.-B. Tian
!-------------------dt,    Nstep, nout, nscreen, res, ires===============
               if (isnap .eq. 0) isnap=10000
ccccc
         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'
!  added by F.-B. Tian to achieve restart capability here
      if(is_res .eq. 1) then
! NOTES: to restart, one must read time,xyzt,dispfult,dispt,vel,acc 
!              disp(i)=0.0
!              dispt(i)=0.0
!              delu(i)=0.0
!              vel(i)=0.0
!              acc(i)=0.0	
	call nond_restart_read(time,itime0,tf1,tf2,tf3,maxxyz,xyzt, 
     & dispfult,disp,dispt, delu, vel, acc )
	else
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
      	time = 0.0
	    itime0=1
	endif
c
*        write(*,*)' '
*        write(*,*)'TYPE: force max  | # of incs | max inc'
*        call zzwrt('  -->  ')
*        read (ikbd,*) force_max,npt,npt_max
*        write(ilog,*) force_max,npt,npt_max,'  ::force,npt,max'
*        dforce=force_max/npt
*        dforce_max=dforce
*        dforce_min=dforce/8
*        dforce0=dforce
*        iprcnt=1
         nsnapmax=npt
         write(*,*)' '
*        write(*,*)' suggest     1                  1        5 '
         write(*,*)'    u = u + beta*du   K = Ke + gamma Kg'
         write(*,*)'INPUT:   beta   |  gamma  |  ramp    '
         write(*,*)' suggest    1        1        5 '
         call zzwrt('  -->  ')
         read(ikbd,*) beta0,gamma0,maxramp0  
!relaxtion parameters: beta0 & gamma0, maxramp0
         write(ilog,*) beta0,gamma0,maxramp0    ,'  ::b g r '
         write(*,*)' '
         write(*,*)'@@ algor 1=full N-R  2=mod N-R '
         write(*,*)'INPUT:   algor | iter max |  tolerance  (<1.0E-5)'
         call zzwrt('  -->  ')
         read(ikbd  ,*) imodify,iter_max,dtol
         write(ilog,*) imodify,iter_max,dtol,'  ::N-R it max tol '
c
c        set integration constants for alpha-method
         zalpha=-0.04
         zalpha=-0.04
c        zalpha=-0.20
c        zalpha=-0.40
         zalpha= 0.00
         zbeta = 0.25*(1.0-zalpha)**2
         zgamma= 0.50*(1-2*zalpha)
         write(ilog,81)'@@ zalpha beta gamma: ',zalpha,zbeta,zgamma
c
         write(*,*)' '
! What's the function of this section?
! choose the point and its degree and (0 displacement 1 velocity and 2 acceleration)
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
         write(*,*)'@@ Eigen types: 0=none  1=      2=subspace '
         write(*,*)'INPUT type |  rate | # vectors '
         call zzwrt('  -->  ')
         read(ikbd,*) ieig_typ,ieig_cnt,ieig_vec 
*                    if (ieig_cnt .lt. 1) ieig_cnt=1
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
*           write(*,*) maxstiff,maxnode,neq,' ms'
            write(ilog,*)'@@ FORMing lumped  mass matrix '
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
            write(ilog,*)'@@Total Mass: ',z1
            write(iout,*)'@@Total Mass: ',z1
c              write(iout,86) (mass(j),j=1,neq)
!            call zztimer(ilog,'<-- FORMmass')
!            call zzflush(ilog)
c
!Input the gravity if it is activated, igravity_on is controled by the special meterials 
!This section does not need to input data from instad and force files
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
*                   write(iout,*) i, grav(i)
                   z1=z1+grav(i)
                enddo
                write(iout,*)'Total gravity force: ',z1
*               call gravity( respce8(ir2),respce8(ir3),idbc,maxnode)
            else
                do i=1,neq
                   grav(i)=0.0  
                enddo
            endif
!==grav(i) now is the gravity term: (rho*s)*g
c
c
         rewind(ilod)
         do i=1,neq
            read(ilod) fmag(i) 
         enddo
         write(*   ,*) '@@ Reloaded  {P}   OK'
         write(ilog,*) '@@ Reloaded  {P}   OK'
c
c        INPUT load history from a file force.53
         call getload_hist(tforce1,maxforce,npt,deltat)
*        do i=1,npt
*           write(ilog,*)'@@ force: ',tforce1(i,1),tforce1(i,2)
*        enddo
c
c        only for second load
c        INPUT load history from a file force.53
         if (igravity_on .eq. 1) then
             write(*,*)'@@ Gravity: '
             call getload_hist(tforce2,9000,npt,deltat)
         endif
         if (iload_num .eq. 3) then
             open(unit=itmp,file=fdn//'stadyn.ld3')
             rewind(itmp)
             do i=1,neq
                read(itmp,*) fmag3(i)
             enddo
             write(*   ,*) '@@ Reloaded  {P}_3   OK'
             write(ilog,*) '@@ Reloaded  {P}_3   OK'
             close(itmp)
c
c            INPUT load history from a file and interpolate 
*            call getload( npt,deltat,tforce2, maxforce, ilog,ikbd,itmp)
             call getload_hist(tforce3,90000,npt,deltat)
         endif
c
         write(*,*)'                    reg  | grav | 3rd |     '
         write(*,*)'TYPE force scales:   P1  |  P2  |  P3 | P4  '
         call zzwrt('  -->  ')
         read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
! The scaled parameters of the load histories
         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
c
c
c
c        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
         write(ilog ,*)'@@  Beginning increment analysis   '
c
c        Initialize  times, disp, etc
c
! Fang-Bao Tian setup the restart.
      if(is_res .eq. 1) then
	! do nothing!

	else
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
	endif
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
           itime=1
!           time=0
           ziter=0
           k1old=1
           k2old=1
           k3old=1
!           tf1=tforce1(1,2)
!           tf1s=tf1*scale_p1
           tf_1=tforce1(1,2)
           tf1s=tf_1*scale_p1
*          write(idyn,182) time,tf1s,(velout(n), n=1,nout),ziter,tf1
           rewind(isnp)
           force=0
*          write(isnp) time, neq, nsnapmax 
*          do i= 1, neq  
*             write(isnp) dispt(i)
*          enddo
c
         kount = 0
         ksnap=1
         ksnap=0
         inc=1
         icarriage=0
         call zzwrt(' @@ ')
c
*          tf1=tforce1(1,2)*scale_p1
*          write(idyn,182) time,tf1,(velout(n), n=1,nout )
*    &                         ,tf2,tf3
c
c        prepare BIG INCREMENT LOOP
*        imodify=1  !full N-R
*        imodify=2  !modified N-R
*        imodify=1  
         itime=0
         iter=0
! Fang-Bao Tian setup the restart.
      if(is_res .eq. 1) then
	! do nothing!

	else
         tf1=0
         tf2=0
         tf3=0
	endif
         deltf1=1
 201     continue

      if(is_res .eq. 1) then
	! do nothing!
	istart0=itime0+1
	iend0=itime0+npt-1
	else
            call formstif_3D(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
c           SECOND assemble the geometric stiffness
            call formgeom_3D(geom(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
	istart0=2
	iend0=npt
	endif
c
c
!============Problem 1: restart
!============Problem 2:contact force
!============Problem 3: coupling Solid-FEM
c        BIG TIME LOOP
!            call zzcount(0,inc,icarriage)
!         do 200 itime= itime0+1,itime0+npt
	   do 200 itime= istart0, iend0
		  write(*,*)'Begin of Itime=', itime
!            call zzcount(itime,inc,icarriage)
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
!            time = real(itime-1)*deltat
		  time  = time + deltat        ! Fang-Bao Tian
            if (itime .eq. 2) then
                maxramp=20
                maxramp=maxramp0
            else
                maxramp=maxramp0
            endif
*220           continue
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
! interpolate the force from the force.53 files
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            if (igravity_on .eq. 1) then
                call load_interp(tforce2,9000,k2old,time,tf2)
            endif
            if (iload_num .eq. 3) then
                call load_interp(tforce3,9000,k3old,time,tf3)
            endif
            pnorm=0.0
            do i= 1, neq  
!=======================================================================
! tf1/tf2/tf3 time changing part, fmag/grav/fmag3 constant part, scale_p* scale of the force
! This is important part to applicate this code!
! Dr. Fangbao Tian add the comments.  
! The fluid can be applied here first, as pload(i)=pload(i)+fxyz(i)
               pload(i) =  tf1*fmag(i)*scale_p1 + tf2*grav(i)*scale_p2
     &                    +tf3*fmag3(i)*scale_p3
               pnorm = pnorm + pload(i)**2
            enddo 
            pnorm = sqrt(pnorm/neq)
!=======================================================================
c
!            if (imodify .eq. 2 .OR. itime .eq. 1) then
            if (imodify .eq. 2) then
c               modified N-R
*               write(iout,*) ' before iter max ',iter,iter_max
            call formstif_3D(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
c           SECOND assemble the geometric stiffness
            call formgeom_3D(geom(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
!                if (itime .eq. 1) goto 290 !Fang-Bao Tian
c
c              Form effective stiffness matrix by adding geometry matrix
                do i= 1, neq  
                   iloc=nloc(i)
                   do 25 j= 1, iprofv(i)
                      ii=iloc+j-1
                       stfeff(iloc+j-1) = 
     &               (1+zalpha)*(stf(iloc+j-1) +gamma0*geom(iloc+j-1))
  25                continue
                    stfeff(iloc) =  stfeff(iloc)  
     &                     +(1/(zbeta*dt**2))*mass(i)
     &                     +(1+zalpha)*(zgamma/(zbeta*dt))*dmm*mass(i)
                enddo   
c
c               Decompose effective stiffness matrix
                ier1=0
            call uduCOL_D( stfeff,maxstiff,neq,ierror,iprofv,nloc,z0,z1)
                write(iout,84) '@@ min/max UDU: ',z0,z1,itime,iter
                if (ierror.eq.0) then
                    write(*,*)'ERROR: zero diagonal term'
                    return
                endif
                write(iout,*) ' iter max ',iter,iter_max
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                 full N-R
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
!=========Fangbao Tian=========
!=This part further deal with the force.
! The fluid can be applied here second, as ploadm(i)=ploadm(i)+fxyz0(i). What is fxyz0? 
! This part is different from that in shell element code. Is it equivalent to that?
              pnormi=0.0
              pnorm =0.0
              do i= 1, neq  
                 aterm = aap(i) + (  1/(zbeta*dt**2))*(dispi(i)-uup(i))
                 vterm = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
                 dpp = -mass(i)*aterm
     &                 - (1+zalpha)*dmm*mass(i)*vterm 
     &                 +    zalpha *dmm*mass(i)*vel(i)
                 ploadm = tf0*fmag(i)*scale_p1
                 ploadm =  tf10*fmag(i)*scale_p1 + tf20*grav(i)*scale_p2
     &                    +tf30*fmag3(i)*scale_p3
               pnorm = pnorm + pload(i)**2
c!!!                      include gravity?
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
!=====================================================================
              if (imodify .eq. 1) then 
c                 full N-R
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
*                 write(iout,*)'STIFFs: diag e+g:'
*                 write(iout,86) (stfeff(nloc(j)),j=1,neq)
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
*                 write(iout,*) ' iter max ',iter,iter_max
              endif
!================================================================
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c             solve for Delta u
              ier1=0
              call bakCOL_D(stfeff,maxstiff,dloadi, neq, dui,ierror,
     &                                         iprofv,iprofh,nloc)
*             write(ilog,*)'@@  dui ',dui(neq),iter
*             write(iout,*) ' dui : '
*             write(iout,86) (dui(j), j=1,neq)
              dsumd=0.0
              do i=1,neq
                 dsumd=dsumd+ abs(dui(i))
              enddo
              write(ilog,*)'@@  |dui| ',dsumd/neq,iter
c
c             increment Ui
*             maxramp=10 
              if (iter .le. maxramp) then
*                 beta=1.0/(maxramp+1-iter)
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
*             write(iout,*) ' ui : ',beta
*             write(iout,86) (dispi(j), j=1,neq)
c
              call update_geom_3D( dispi, dispfuli,maxdim,idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1))
c
*                 write(iout,*) ' updated [ui] : '
*                 do i= 1, nnp  
*                  write(iout,82) (dispfuli(i,j), j=1,3) 
*                 enddo
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
*           iblowup=0
*           idiv   =0
c
c           set values for time t+dt -> t
            do i= 1, neq  
               delu(i)=dispi(i)-dispt(i)
               dispt(i)=dispi(i)
*              write(ilog,*) i,dispt(i),dispi(i)
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
!                 call zzflush(idyn)
               call write_coordHEX20(xyzt,                ! Fangbao Tian
     &                     npijkm,maxnode,maxelem,maxxyz,itime-1)!
             endif
*            psload=itime*dforce
*            psload=dforce
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
                if (mod(itime+1,ieig_cnt) .eq. 0 .OR. itime .eq. 1
     &                                       .OR. itime .eq. npt) then
c
c                   monitors
c                   compute kinetic energy, etc 
                    s_energy=0.0
                    t_energy=0.0
                    do i=1,neq
                       t_energy = t_energy + 0.5*mass(i)*vel(i)**2
                    enddo
                    p_energy=0.0
                    do i=1,neq
                       p_energy = p_energy + 0.5*pload(i)*dispt(i)
                    enddo
c
                        xc=0.0
                        yc=0.0
                        zc=0.0
                        do i=1,nnp
                           xc = xc + xyzt(i)   
                           yc = yc + xyzt(maxnode+i)
                           zc = zc + xyzt(maxnode*2+i)
                        enddo
                        xc=xc/nnp
                        yc=yc/nnp
                        zc=zc/nnp
                        gsize=0.0
                        do i=1,nnp
                           gsize = gsize + (xyzt(i) - xc)**2   
     &                                 + (xyzt(maxnode+i) - yc)**2
     &                                 + (xyzt(maxnode*2+i) - zc)**2
                        enddo
                        gsize=sqrt(gsize/nnp)
                        t_disp=0.0
                        do i=1,neq
                           t_disp = t_disp + dispt(i)**2
                        enddo
                        u_aver=sqrt(t_disp/neq)
*                   write(itmp3,85) time,tf3,(sens(j),j=1,min(40,neq))
*    &                         ,size,u_aver
c
c                   reconstruct tangent stiffness
                    do i= 1, neq  
                       iloc=nloc(i)
                       do j= 1, iprofv(i)
                          stfeff(iloc+j-1) = stf(iloc+j-1)
     &                                     + gamma0*geom(iloc+j-1)
                       enddo
                    enddo
c
                    if (ieig_typ .eq. 0) then
c                       none
                    elseif (ieig_typ .eq. 1) then
c                       vector iter
                    elseif (ieig_typ .eq. 2) then
c                       subspace iteration for multiple eigenvalues
c
                        call sub_spc(stfeff,isize8,dloadi,dui, 
     &                           mass,     iprofv,iprofh,nloc,
     &                           respceb,ieig_loc,wk_eig)
c
c                       write(itmp5,85) time,(wk_eig(j),j=1,10)
c    &                                  ,tf1,tf2,tf3
c    &                             ,s_energy,t_energy,p_energy
                        zlam=wk_eig(1)
                        do k=1,ieig_vec
                           do j=1,neq
                              wk0(j) = respceb(ieig_loc+neq*(k-1) -1+j)
                           enddo
c                          normalize shape
                           wkmax0 = 0.0
                           do i= 1, neq  
                              if (abs(wk0(i)).gt. wkmax0) 
     &                                wkmax0=abs(wk0(i))
                           enddo
                           do i= 1, neq  
                              wk0(i) = wk0(i)/wkmax0
                           enddo
**                         alpha1=1
**                         alpha2=0
**                         do i= 1, neq  
**                            wk(i) = alpha1*wk1(i) + alpha2*wk(i) 
**                         enddo
**                         write(itmp4) time, neq,ieig_max 
**                         do i= 1, neq  
**                            write(itmp4) wk1(i)
**                         enddo
                           write(itmp4) time, neq, ieig_max 
                           do i= 1, neq  
                              write(itmp4) wk0(i)
                           enddo
c
c                          keep shape for ping
                           if (k .eq. 1) then
                               do i= 1, neq  
                                  wk1(i) = wk0(i) 
                               enddo
                           endif
                           if (k .eq. 2) then
                               do i= 1, neq  
                                  wk2(i) = wk0(i) 
                               enddo
                           endif
                        enddo
                    endif
c                   bottom eigen type
c
                tf1s=tf1*scale_p1
*                   write(itmp5,85) time,tf1s,(wk_eig(j),j=1,ieig_vec)
                    write(imon ,85) time,tf1s,(wk_eig(j),j=1,ieig_vec)
*    &                                  ,tf1,tf2,tf3
     &                             ,s_energy,t_energy,p_energy
     &                         ,gsize,u_aver
!                 call zzflush(imon)
ctag1
                endif
            endif
c           bottom eigenanalysis
c
!  Fang-Bao Tian
! Writing the restart file
            if(mod(itime-1,irestart) .eq. 0) then
	       call nond_restart_write(time,itime,tf1,tf2,tf3,maxxyz,xyzt, 
     & dispfult,disp,dispt, delu, vel, acc )
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            go to load increment
*            goto 201
c
 200  continue
 210  continue
 !  Fang-Bao Tian
! Writing the restart file
! 	       call nond_restart_write(time,itime,tf1,tf2,tf3,maxxyz,xyzt, 
!     & dispfult,disp,dispt, delu, vel, acc )
!         write(*,'(a)')' @@'
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

!=================
!Fang-Bao Tian, 2012.01.10
!=================
	subroutine nond_restart_write(time,itime,tf1,tf2,tf3,maxxyz,xyzt, 
     & dispfult,disp,dispt, delu, vel, acc )
	implicit real*8 (a-h,o-z)
	include 'commons.std'
      parameter( maxdim=20000,maxdof=maxdim*3)
	real*8::xyzt(maxxyz),dispfult(maxdim,3),disp(neq  ),dispt(neq  ),
     &        delu(maxdof),vel(neq),acc(neq)


      ifile = 2012 
      open(unit=ifile, file=fdn//'nond_restart_out.dat', 
     &     form='unformatted')

      write(*  ,*) '@@ Write restart file nond_restart_out.dat'    
      write(ifile) time,itime,tf1,tf2,tf3
      write(ifile) xyzt(1:maxxyz),
     &             dispfult(1:maxdim,1:3),
     &             disp(1:neq),
     &             dispt(1:neq),
     &             delu(1:maxdof),
     &             vel  (1:neq),
     &             acc  (1:neq)

      close(ifile)
      return
	end subroutine

	subroutine nond_restart_read(time,itime,tf1,tf2,tf3,maxxyz,xyzt, 
     & dispfult,disp,dispt, delu, vel, acc )
	implicit real*8 (a-h,o-z)
	include 'commons.std'
      parameter( maxdim=20000,maxdof=maxdim*3)
	real*8::xyzt(maxxyz),dispfult(maxdim,3),disp(neq  ),dispt(neq  ),
     &        delu(maxdof),vel(neq),acc(neq)


      ifile = 2012 
      open(unit=ifile, file=fdn//'nond_restart_in.dat', 
     &     form='unformatted')

      write(*  ,*) '@@ Read restart file nond_restart_in.dat'    
      read(ifile) time,itime,tf1,tf2,tf3
      read(ifile) xyzt(1:maxxyz),
     &             dispfult(1:maxdim,1:3),
     &             disp(1:neq),
     &             dispt(1:neq),
     &             delu(1:maxdof),
     &             vel  (1:neq),
     &             acc  (1:neq)

      close(ifile)
      return
	end subroutine
