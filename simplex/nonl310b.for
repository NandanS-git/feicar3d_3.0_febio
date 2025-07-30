c4100
c     NON-liner STATIC incremental analysis 
      subroutine non_TL3D(stf,geom,wk,pload, disp, fmag,
     &                   dispt , 
     &                   stfeff,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   iprofv, iprofh, nloc,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt)
c
          implicit real*8 (a-h,o-z)
          include 'commons.std'
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=100)
         parameter( maxdim=2000,maxdof=maxdim*3)
c
         integer iprofv(neq), iprofh(neq)
     &          ,nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer nloc(neq)
         integer inod(ireadmax,2),mnod(iconelem),nfollow(ireadmax,4)
c
         real*8 tforce1(maxforce,2), tforce2(5000,2)
         real*8 velout(ireadmax),angle(ireadmax,3)
         real*8 stf(maxstiff),wk(neq), geom(maxstiff),pload(neq )
         real*8 stfeff(maxstiff)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(10000)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 xyzi(10000)
         real*8 dispt(neq)
         real*8 dsumd,dsumz,dtol,pnorm,pnormi,wnorm,wnormi
*        real*8 dispfult(900,3), dispfuli(900,3),delu(3000)
*        real*8 dloadi(3000),floadi(3000)
*        real*8 dispi(3000),dui(3000),grav(3000)
         real*8 dispfult(maxdim,3), dispfuli(maxdim,3),delu(maxdof)
         real*8 dloadi(maxdof),floadi(maxdof)
         real*8 dispi(maxdof),dui(maxdof),grav(maxdof)
cccccccccccccccccccccccccccccc
c
*         write(ilog,*)'@@  non_TL3D ver244  March 2006 '
*         write(ilog,*)'@@  non_TL3D ver262  November 2006 '
          write(ilog,*)'@@  non_TL3D ver264  November 2006 '
c
c
*                do i=1,nnp
*                   xyz0(i)=xyz999(i)
*                   xyz0(maxnode+i)=xyz999(maxnode+i)
*                   xyz0(maxnode*2+i)=xyz999(maxnode*2+i)
*                   write(iout,*) i,xyz0(i)
*                enddo
c
c        algorithm parameters
           write(*,*) maxstiff,' ms'
*          stop
         beta0=1.0
         iter_max=60
         dtol=0.0001
         gamma0=1.0
         kref=2
*        iblowup=0
*        idiv=0
c
           do i=1,nnp
*        write(*,*) i,maxnode,maxnode*2+i
              xyzt(i)=xyz0(i)
              xyzt(maxnode+i)=xyz0(maxnode+i)
              xyzt(maxnode*2+i)=xyz0(maxnode*2+i)
           enddo
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
*        read(ikbd,*) force_max,npt,npt_max
*        write(ilog,*) force_max,npt,npt_max,'  ::force,npt,max'
*        dforce=force_max/npt
*        dforce_max=dforce
*        dforce_min=dforce/8
*        dforce0=dforce
*        iprcnt=1
         nsnapmax=npt
         write(*,*)' '
         write(*,*)'    u = u + beta*du   K = K_E + gamma K_G'
*        write(*,*)' suggest:    1                   1         5'
         write(*,*)'INPUT:   beta   |  gamma  |  ramp    '
         write(*,*)' suggest:   1        1        5  '
         call zzwrt('  -->  ')
         read(ikbd,*) beta0,gamma0,maxramp0  
         write(ilog,*) beta0,gamma0,maxramp0    ,'  ::b g r '
         write(*,*)' '
         write(*,*)'INPUT:   iter max   |  tolerance  (<1.0E-5) '
         call zzwrt('  -->  ')
         read(ikbd   ,*) iter_max,dtol
         write(ilog,*) iter_max,dtol,'  ::it max tol '
c
         write(*,*)' '
         write(*,*)'CHOOSE nodal output:'
         nout=0
27       continue
         nout=nout+1
         write(*,*)'TYPE:  node# |  DoF               <0 0 to end>'
         call zzwrt('  -->  ')
         read(ikbd,*) inode,idof
         write(ilog,*) inode,idof,' ::node  dof '
         if (inode .ne. 0) then
             inod(nout,1)=idbc(inode,idof)
             irate=0
             inod(nout,2)=irate
             write(ilog,*) '@@ eqn # ',inod(nout,1)
             goto 27
         endif
         nout=nout-1
         write(*,'(a)')' @@ '
c
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
*           write(ilog,*)'@@ force: ',tforce1(i,2)
*        enddo
c
c        only for second load
         if (igravity_on .eq. 1) then
             write(*,*)'@@ Gravity: '
             call getload_hist(tforce2,5000,npt,deltat)
         endif
c
*        write(*,*)'TYPEs 0=regular   1=follower '
         write(*,*)'                    reg  | grav |     | foll'
         write(*,*)'TYPE force scales:   P1  |  P2  |  P3 | P4  '
         call zzwrt('  -->  ')
         read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
         ifollower = (scale_p4 + 0.1)
c
c? 
         nfoll=0
         if (ifollower .eq. 1) then
             nfoll=0
28           continue
             nfoll=nfoll+1
             write(*,*)'TYPE:  node #s  1 | 2 | 3 | loc <0 0 0 0to end>'
             call zzwrt('  -->  ')
             read(ikbd,*) (nfollow(nfoll,j), j=1,4)
             write(ilog,*) (nfollow(nfoll,j), j=1,4),' ::nodes '
             if (nfollow(nfoll,1) .ne. 0) then
                 goto 28
             endif
         endif
         nfoll=nfoll-1
         write(*,'(a)')' @@ '
c
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
*          tf1=tforce1(1,2)*scale_p1
*          write(idyn,182) time,tf1,(velout(n), n=1,nout )
*    &                             ,(angle(1,n),n=1,3),ziter
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
c        prepare BIG INCREMENT LOOP
         itime=0
         iter=0
         tf1=0
         tf2=0
         deltf1=1
 201     continue
c
            write(iout,*) ' 0000   iter max ',iter_max
c        BIG TIME LOOP
         do 200 itime= 1,npt
            call zzcount_2(itime,inc,icarriage)
            kount = kount + 1
            ksnap = ksnap + 1
            time = real(itime-1)*deltat
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
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            if (igravity_on .eq. 1) then
                call load_interp(tforce2,5000,k2old,time,tf2)
            endif
*           tf1 = tf1*scale_p1
            pnorm=0.0
            do i= 1, neq  
               pload(i) =  tf1*fmag(i)*scale_p1
     &                    +tf2*grav(i)*scale_p2
               pnorm = pnorm + pload(i)**2
            enddo 
            pnorm = sqrt(pnorm/neq)
            if (itime .eq. 1) then
                kount=iprcnt
                goto 290
            endif
c
c           modified N-R
*           write(iout,*) ' before iter max ',iter,iter_max
*           call formstif_3D(stf(1), npijkm,idbc,maxnode,maxelem,
*    &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1),
*    &                      dispfult,
*    &                      nelt,iglobal,
*    &                      wk(1),
*    &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
c           SECOND assemble the geometric stiffness
*           call formgeom_3D(geom(1), npijkm,idbc,maxnode,maxelem,
*    &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1),
*    &                      dispfult,
*    &                      nelt,iglobal,
*    &                      wk(1),
*    &                      prop,nmprop,ippp,iprofv,nloc )
*           write(iout,*) ' after formgeom iter max ',iter,iter_max
c           Form effective stiffness matrix by adding geometry matrix
*           do i= 1, neq  
*              iloc=nloc(i)
*              do j= 1, iprofv(i)
*                stfeff(iloc+j-1) = stf(iloc+j-1) +gamma0*geom(iloc+j-1)
*              enddo   
*           enddo   
*           write(iout,*) ' before write STIFF iter max ',iter,iter_max
*           write(iout,*)'STIFFs: diag e+g:'
*           write(iout,86) (stfeff(nloc(j)),j=1,neq)
c
c           Decompose effective stiffness matrix
*           ier1=0
*           call uduCOL_D( stfeff,maxstiff,neq,ierror,iprofv,nloc,z0,z1)
*           write(iout,84) '@@ min/max: ',z0,z1,itime,iter
*           if (ierror.eq.0) then
*               write(*,*)'ERROR: zero diagonal term'
*               return
*           endif
*           write(iout,*) ' iter max ',iter,iter_max
c
c
c           Equilibrium ITERation
c           prepare iterates
            do i=1,neq
               dispi(i)=dispt(i)
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
c                full N-R
c
c           if (itime .eq. 2) then
*           if (itime .le. 2000) then
c               full N-R
*               write(iout,*) ' before iter max ',iter,iter_max
*               call formstif_3D(stf(1), npijkm,idbc,maxnode,maxelem,
*    &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                      xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1),
*    &                      dispfuli,
*    &                      nelt,iglobal,
*    &                      wk(1),
*    &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
c               SECOND assemble the geometric stiffness
*               call formgeom_3D(geom(1), npijkm,idbc,maxnode,maxelem,
*    &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                      xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1),
*    &                      dispfuli,
*    &                      nelt,iglobal,
*    &                      wk(1),
*    &                      prop,nmprop,ippp,iprofv,nloc )
*               write(iout,*) ' after formgeom iter max ',iter,iter_max
c               Form effective stiff matrix by adding geometry matrix
*               do i= 1, neq  
*                  iloc=nloc(i)
*                  do j= 1, iprofv(i)
*                    stfeff(iloc+j-1) = stf(iloc+j-1) 
*    &                                 +gamma0*geom(iloc+j-1)
*                  enddo   
*               enddo   
*               write(iout,*)'STIFFs: diag e+g:'
*               write(iout,86) (stfeff(nloc(j)),j=1,neq)
c
c               Decompose effective stiffness matrix
*               ier1=0
*               call uduCOL_D( stfeff,maxstiff,neq,ierror
*    &                        ,iprofv,nloc,z0,z1)
*               write(iout,84) '@@ min/max: ',z0,z1,itime,iter
*               if (ierror.eq.0) then
*                   write(*,*)'ERROR: zero diagonal term'
*                   return
*               endif
*               write(iout,*) ' iter max ',iter,iter_max
*         endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if (iter .eq. 1000) then 
c               estimate disp based on previous rate
                fract=1.0
                if (itime .gt. 2) then
                    if ( abs(deltf1) .le. 1.0e-6 ) deltf1=1.0e-6
                    fract=0.5*(tf1-tf0)/deltf1
                    if ( abs(fract) .le. 0.1) fract=0.1
                    if ( abs(fract) .ge. 1.0) fract=1.0
                    write(iout,*)'@@ fract: ',fract,itime
                    write(iout,82) tf1,tf0,deltf1
                    do i=1,neq
                       dispi(i)=dispi(i)+delu(i)*fract
                       write(iout,82) dispi(i)-delu(i)*fract
     &                               ,delu(i)*fract
     &                               ,dispi(i),fract
                    enddo
*                   call update_geom_3D( dispi, dispfuli, idbc,maxnode,
*    &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                    xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1))
                endif
            endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                full N-R
c
c           forces from body stresses
            call body_stress_3D( dispfuli, iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                    xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadi,iter,maxdim,maxdof
     &                 )
c
c           form load increment
*           if (iter .le. 2) then
*              do i= 1, neq  
*                 pload(i) =  0.5*(tf1+tf0)*fmag(i)
*              enddo 
*           else
*              do i= 1, neq  
*                 pload(i) =  tf1*fmag(i)
*              enddo 
*           endif
ctag1
            if (ifollower .eq. 1) then
                call follow_load(nfollow,ireadmax,nfoll,maxnode,
     &                    xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1),
     &                    angle      )
                pnorm=0.0
                do i= 1, neq  
                   pload(i) = 0.0
                enddo 
                do n=1,nfoll
                   n1=nfollow(n,1)
                   ind1=idbc(n1,1)
                   ind2=idbc(n1,2)
                   ind3=idbc(n1,3)
                       idof=nfollow(n,4)
                       ind4=idbc(n1,idof)
                       pp4 =  tf1*fmag(ind4)
                       aax =     angle(n,1) 
                       aay =     angle(n,2) 
                       aaz =     angle(n,3) 
                       ppx = pp4*angle(n,1) 
                       ppy = pp4*angle(n,2) 
                       ppz = pp4*angle(n,3) 
                   if (ind1 .gt. 0) pload(ind1)=ppx
                   if (ind2 .gt. 0) pload(ind2)=ppy
                   if (ind3 .gt. 0) pload(ind3)=ppz
                enddo
                pnorm=0.0
                do i= 1, neq  
                   pnorm = pnorm + pload(i)**2
                enddo 
                pnorm = sqrt(pnorm/neq)
            endif
            pnormi=0.0
            pnorm =0.0
            do i= 1, neq  
               dloadi(i) = pload(i) - floadi(i)
               pnormi = pnormi + dloadi(i)**2
               pnorm  = pnorm  + pload (i)**2
            enddo 
            pnormi=sqrt(pnormi/neq)
            pnorm =sqrt(pnorm /neq)
            write(iout,*) ' force: ',tf1
*           write(iout,*) ' dP norm ',pnormi,pnorm
*           write(iout,86) (dloadi(j),j=1,neq)
c
                call formstif_3D(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                      xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1),
     &                      dispfuli,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
c               SECOND assemble the geometric stiffness
                call formgeom_3D(geom(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
*    &                      xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1),
     &                      dispfuli,maxdim,
     &                      nelt,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
*               write(iout,*) ' after formgeom iter max ',iter,iter_max
c               Form effective stiff matrix by adding geometry matrix
                do i= 1, neq  
                   iloc=nloc(i)
                   do j= 1, iprofv(i)
                     stfeff(iloc+j-1) = stf(iloc+j-1) 
     &                                 +gamma0*geom(iloc+j-1)
                   enddo   
                enddo   
*               write(iout,*)'STIFFs: diag e+g:'
*               write(iout,86) (stfeff(nloc(j)),j=1,neq)
c
c               Decompose effective stiffness matrix
                ier1=0
                call uduCOL_D( stfeff,maxstiff,neq,ierror
     &                        ,iprofv,nloc,z0,z1)
                write(iout,84) '@@ min/max UDU: ',z0,z1,itime,iter
                if (ierror.eq.0) then
                    write(*,*)'ERROR: zero diagonal term'
                    return
                endif
*               write(iout,*) ' iter max ',iter,iter_max
*         endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c           solve for Delta u
            ier1=0
            call bakCOL_D(stfeff,maxstiff,dloadi, neq, dui,ierror,
     &                                       iprofv,iprofh,nloc)
*           write(ilog,*)'@@  dui ',dui(neq),iter
*           write(iout,*) ' dui : '
*           write(iout,86) (dui(j), j=1,neq)
            dsumd=0.0
            do i=1,neq
               dsumd=dsumd+ abs(dui(i))
            enddo
            write(ilog,*)'@@  |dui| ',dsumd/neq,iter
c
c           increment Ui
*           maxramp=10 
            if (iter .le. maxramp) then
*               beta=1.0/(maxramp+1-iter)
                beta=iter/(maxramp+0.0)
                beta=beta*beta0
                zi=2.0**(iter-1)
                z0=2.0**(maxramp-1)
                beta=zi/z0
                beta=beta*beta0
            else
                beta=1.0
                beta=beta0
            endif
                do i= 1, neq  
                   dispi(i)= dispi(i) + beta*dui(i)
                enddo
*           write(iout,*) ' ui : ',beta
*           write(iout,86) (dispi(j), j=1,neq)
c
            call update_geom_3D( dispi, dispfuli,maxdim, idbc,maxnode,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1))
c
*               write(iout,*) ' updated [ui] : '
*               do i= 1, nnp  
*                  write(iout,82) (dispfuli(i,j), j=1,3) 
*               enddo
c
c           work norm 
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
c           test for CONvergence
            dsumd=0.0
            dsumz=0.0
            do i=1,neq
               dsumd=dsumd+ (dui(i)*beta)**2
               dsumz=dsumz+ dispi(i)**2
            enddo
            dsumd=(dsumd/neq)**0.5
            dsumz=(dsumz/neq)**0.5
            if (dsumz .lt. dtol/10.0) dsumz=dtol/10.0
            dnorm=dsumd/dsumz
            write(iout,*)' dUi: ',dsumd,dsumz
            write(iout,*)' norm: ',dnorm,iter
            write(iout,*)' dtol: ',dtol 
            if (itime .eq. -2) then
                ptol=0.5e-2
                ptol=0.25e-2
                if (pnormi/pnorm .lt. ptol) then
c                   converged, update geom etc, write results
                    write(iout,*)' CONVERGED: '
                    goto 100
                else
c                   not converged, iterate again
                    goto 120
                endif
            else
                if (dnorm .lt. dtol) then
c                   converged, update geom etc, write results
                    write(iout,*)' CONVERGED: '
                    goto 100
                else
c                   not converged, iterate again
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
c
c
c           Print out results of interest for this time loop
 290        contin ue
            if (iprcnt .eq. kount) then
                kount = 0
                nout1=nout-1
                do 172 n=1,nout
                   iprnode=inod(n,1)
                   if (iprnode .eq. 0) then
                       velout(n) = 0.0
                       goto 172
                   endif
                   if (inod(n,2) .eq. 0) then
                       velout(n)=dispt(iprnode)
                   endif
 172            continue
                ziter=iter*1.0
                tf1s=tf1*scale_p1
                tf2s=tf2*scale_p2
                write(idyn,182) time,tf1s,tf2s,(velout(n), n=1,nout )
     &                              ,ziter,tf1,tf2
*    &                              ,aax,aay,aaz,ziter
!                call zzflush(idyn)
             endif
*            psload=itime*dforce
*            psload=dforce
             if (isnap .eq. ksnap .OR. itime .eq. 1) then
                 write(isnp) time, neq, nsnapmax 
                 do i= 1, neq  
                    write(isnp) dispt(i)
                 enddo
                 ksnap=0
             endif
c
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
c
      return
      end
c
c
c     GET applied transient LOAD HISTory 
      subroutine getload_hist(tforce,maxforce,npt,deltat)
         implicit real*8 (a-h,o-z)
         include 'commons.std'
         real*8 tforce(maxforce,2),wk(100)
         character*40 fylname
c
c        read data file
         call zzwrt(' @@ TYPE:  Load_Filename --> ')
         read(ikbd,'(1a40)') fylname
         write(*,*)' '
         write(ilog,'(a3)') '@@ load filenmae '
         write(ilog,*) fylname
c
         write(*,*)'@@ INPUT: Time col | Load col | # cols '
         call zzwrt(' --> ')
         read(ikbd,*) icol,jcol,kcol
         write(ilog,*) icol,jcol,kcol,' ::t f #cols'
c
         open(unit=itmp ,file=fdn//adjustl(fylname))
         rewind itmp 
         do i=1,maxforce-1 
            read(itmp,*,end=220) (wk(k),k=1,kcol)
            tforce(i,1)=wk(icol)
            tforce(i,2)=wk(jcol)
         enddo
c        filled up to max, make space available for tail
*        i=maxforce-0
c
 220     continue
         nmax=i-1
         if (nmax .lt. 1) nmax=1
         write(ilog,*)'@@ READ load points # =  ', nmax
         write(*   ,*)'@@ READ load points # =  ', nmax
c        set beginning time at zero
         do i=1,nmax
            tforce(i,1)=tforce(i,1)-tforce(1,1)
         enddo
         tmax=tforce(nmax,1)
         if (tmax .lt. npt*deltat) then
*            tforce(nmax+0,2)=0.0         
             tforce(nmax+1,1)=(npt+2)*deltat
             tforce(nmax+1,2)=0.0         
         endif
         do i=1,nmax
            write(ilog,*)'@@ load: ',tforce(i,1),tforce(i,2)
         enddo
 86      format(1x,a,1x,20(g12.6,1x))
         close (itmp )
c 
      return
      end
c
cccccccccccccc
c
c     LOAD INTERPolation
      subroutine load_interp(tforce,maxforce,kold,time,tf)
         implicit real*8 (a-h,o-z)
         real*8 tforce(maxforce,2)
c
c        adjust times and interpolate
c 
         do 210 k = 1,maxforce
            t0=tforce(k+0,1)
            t1=tforce(k+1,1)
            f0=tforce(k+0,2)
            f1=tforce(k+1,2)
c
            if (t1 .gt. time) then
                tf =(f1-f0)/(t1-t0)*(time-t0)+f0
                kold = k 
                return    
            endif
 210     continue
c 
      return
      end
c
c
c
c     Nodal loads due to body stresses for 3D solids
      subroutine body_stress_3D( dispful, iglobal,
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xord0,yord0,zord0,
*    &                                    xord,yord,zord,
     &                   prop,nmprop,ippp,
     &                   gforce,iter,maxdim,maxdof
     &                )

          use contact

          implicit real*8 (a-h,o-z)
          include 'commons.std'
#ifdef _OPENMP
          include 'omp_lib.h'
#endif
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3)
         real*8  xord0(nnp), yord0(nnp), zord0(nnp)
*        real*8  xord(nnp), yord(nnp), zord(nnp)
c
         real*8 prop(ippp,10),dd(6,6)
         integer nmprop(nel)
c
         real*8  stress(6,27), force(60),force8(3,20)
         real*8  uvw(3,20),xyz(3,20),ss(6)
c
*        real*8  dispful(900,3)
*        real*8  gforce(3000)
         real*8  dispful(maxdim,3)
         real*8  gforce(maxdof)
c
c
         write(ilog,*)'@@ << in body_stress_3D >>'

         !Ye,initialization of 2nd PK stress & strain
         pks2 = 0.0d0
         strn = 0.0d0
         cauchy = 0.0d0

         do i=1,neq
            gforce(i)=0.0
         enddo

#ifdef _OPENMP
       call    OMP_SET_DYNAMIC(.FALSE.)
       call    OMP_SET_NUM_THREADS(16)
!$OMP PARALLEL 
     & PRIVATE(n,mat,neltype,e0,g0,t0,r0,pl0,
     & dd,k,kmax,node,xyz,uvw,force8,
     & iter)
!$OMP DO
#endif
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
c
                call dmat(e0,g0,dd)
*             write(iout,*)' [D] ',i
*             do ii=1,6
*                write(iout,82) (dd(ii,jj), jj=1,6)
*             enddo
c
!              if (neltype .eq. 5) then
c                 tetra
!                  do k=1,4 
!                     node=npijkm(n,1+k)
!                     xyz(1,k)=xord0(node)
!                     xyz(2,k)=yord0(node)
!                     xyz(3,k)=zord0(node)
c
!                     uvw(1,k)=dispful(node,1)
!                     uvw(2,k)=dispful(node,2)
!                     uvw(3,k)=dispful(node,3)
!                  enddo
!                  call stress_TET( xyz,uvw,dd,ss)
*                     write(iout,*) '@@ stress '
*                     write(iout,82) (ss(jj), jj=1,6)
!                  call bforce_TET( xyz,uvw,ss,force )
*                     write(iout,*) '@@ force '
*                     write(iout,82) (force(jj), jj=1,12)
c
c!!!              replace with HEX
*                 call assembFOR_3D(gforce,force,idbc,maxnode
*    &                          ,npijkm,maxelem,n)
c
!              elseif (neltype .eq. 9) then
!                     write(iout,*) '@@ uvw '
!                  do k=1,8
!                     node=npijkm(n,1+k)
!                     xyz(1,k)=xord0(node)
!                     xyz(2,k)=yord0(node)
!                     xyz(3,k)=zord0(node)
c
!                     uvw(1,k)=dispful(node,1)
!                     uvw(2,k)=dispful(node,2)
!                     uvw(3,k)=dispful(node,3)
*                        write(iout,82) (uvw(jj,k), jj=1,3),node*1.0
!                  enddo
!                  call stress_HEX8( xyz,uvw,dd,stress)
*                    write(iout,*) '@@ stress '
*                    do ii=1,8
*                       write(iout,82) (stress(jj,ii), jj=1,6)
*                    enddo
!                  call bforce_HEX8( xyz,uvw,stress,force8 )
!                     write(iout,*) '@@ force '
*                  do ii=1,8
*                     write(iout,82) (force8(jj,ii), jj=1,3)
*                  enddo
!                  call assembFOR_HEX(gforce,force8,idbc,maxnode
!     &                          ,npijkm,maxelem,n)
*             write(iout,*)'@@ partial gFORCE: ',n
*             write(iout,86) (gforce(k),k=1,neq)
c
!              elseif (neltype .eq. 21) then
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
c
                  ired_int=2
                  ired_int=3
                  ired_int=ri3
c                 call stress_HEX20(neltype,ired_int, xyz,uvw,dd,stress,
c    &                              iter)
*                    write(iout,*) '@@ stress ',n
*                    do ii=1,27
*                       write(iout,82) (stress(jj,ii), jj=1,6)
*                    enddo
*                    write(iout,*) ' '
c                 call bforce_HEX20(neltype,ired_int, xyz,uvw,stress,
c    &                              force8,iter )
                  call nodeforce_HEX20(n,neltype,ired_int, xyz,uvw,dd,
     &                              force8,iter )
*                    write(iout,*) '@@ force ',n
*                 do ii=1,20
*                    write(iout,82) (force8(jj,ii), jj=1,3)
*                 enddo
                  call assembFOR_HEX(gforce,force8,idbc,maxnode
     &                          ,npijkm,maxelem,n)
c
!              endif
 50        continue
!$OMP END DO nowait
!$OMP END PARALLEL
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
c
c
c     Body FORCE from stress for HEXahedron 8 nodes
      subroutine nodeforce_HEX20(n,neltype,ired_int, 
     &                           xyz, uvw,dd,force,
     &                        iter)
         use contact
      
         Implicit real*8 (a-h,o-z)
         include 'commons.std'

         real*8 xyz(3,20),uvw(3,20),u(60),ud(9),ee(6),dd(6,6)
         real*8 force(3,20),ff(60),ss(6),df(60)
         real*8 skk(3,3),def(3,3),scc(3,3),cc(3,3)
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
         nn=0
         do i=1,kmax
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
*                  write(iout,*) nn,u(nn)
            enddo
         enddo
c
         do i=1,kmax*3
            ff(i) = 0.0 
         enddo
         ip_n=0
*        nint=3
*        nint=2
*        nint=3
         nint=ired_int
*        nint=2
c        nint=1
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8   
         ip_n = ip_n + 1
*            lx = icoord(ip_n,1)
*            ly = icoord(ip_n,2)
*            lz = icoord(ip_n,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
*               write(*,84) lx,ly,lz,ip_n,(icoord(ip_n,kk),kk=1,3)
c
c                 deriv oper and jacobian
                  call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det,ud)
c
ccccccccccccccccccccccccccccccccccccc
c                 STRESS
c                 write(iout,*) ' [Bd] '
c                 do i=1,9
c                    write(iout,82) (Bd(i,j), j=1,60)
c                 enddo
                  !do i=1,9
                  !   sum=0.0
                  !   do k=1,kmax*3
                  !      sum=sum+Bd(i,k)*u(k)
                  !   enddo
                  !   ud(i)=sum
                  !enddo
*                 write(iout,*) ' u,x  ',ip_n
*                    write(iout,82) ( u(j), j=1,kmax*3)
*                    write(iout,82) ( ud(j), j=1,9)
c
*        z9=0.5
         z9=1.0
*        if (iter .eq. 1) then
*            z9=0.0
*        endif
         ee(1)= ud(1) + z9*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
         ee(2)= ud(5) + z9*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
         ee(3)= ud(9) + z9*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
         ee(4)= ud(2)+ud(4) + z9*(ud(1)*ud(2)+ud(4)*ud(5) +ud(7)*ud(8))
         ee(5)= ud(6)+ud(8) + z9*(ud(2)*ud(3)+ud(5)*ud(6) +ud(8)*ud(9))
         ee(6)= ud(3)+ud(7) + z9*(ud(1)*ud(3)+ud(4)*ud(6) +ud(7)*ud(9))

         !Ye,store the Lagrangian strain
         strn(n,ip_n,1)=1.0*ee(1)
         strn(n,ip_n,2)=1.0*ee(2)
         strn(n,ip_n,3)=1.0*ee(3)
         strn(n,ip_n,4)=0.5*ee(4)
         strn(n,ip_n,5)=0.5*ee(5)
         strn(n,ip_n,6)=0.5*ee(6)

         !Ye,debug
         !if ((n.eq.425).and.(lx.eq.1).and.(ly.eq.1).and.(lz.eq.1) )then
         !   write(81,*)nel,ee(1),ee(2),ee(3),ee(4),ee(5),ee(6)
         !   write(81,*)'--------------'
         !endif

*                 write(iout,*) ' Exx  ',ip_n
*                 write(iout,82) ( ee(j), j=1,6)
c
c        z8=0.0
c        en(1)= ud(1) + z8*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
c        en(2)= ud(5) + z8*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
c        en(3)= ud(9) + z8*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
c        en(4)= ud(2)+ud(4) + z8*(ud(1)*ud(2)+ud(4)*ud(5) +ud(7)*ud(8))
c        en(5)= ud(6)+ud(8) + z8*(ud(2)*ud(3)+ud(5)*ud(6) +ud(8)*ud(9))
c        en(6)= ud(3)+ud(7) + z8*(ud(1)*ud(3)+ud(4)*ud(6) +ud(7)*ud(9))
c
c                       write(iout,* ) ' strain N L d ',iter
c                    do j=1,6
c                       write(iout,82) ee(j),en(j),en(j)-ee(j)
c                    enddo
*                    do j=1,6
*                       sum = 0.0  
*                       do kk=1,24
*                          sum = sum + bb(j,kk)*u(kk)
*                       enddo
*                       ee0 = sum
*                       write(iout,82) ee(j),ee0 
*                    enddo
*                 write(iout,82) (ee(j)*1.0e6, j=1,6)
                  do i=1,6
                     sum = 0.0  
                     do kk=1,6
                        sum = sum + dd(i,kk)*ee(kk)
                     enddo
                     ss(i)=sum
                     !Ye,store the 2nd PK stress
                     pks2(n,ip_n,i)=ss(i)
                  enddo

         !Ye,debug
         !if ((n.eq.425).and.(lx.eq.1).and.(ly.eq.1).and.(lz.eq.1) )then
         !   write(82,*)nel,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6)
         !   write(82,*)'========='
         !endif

                 !Ye,tranfer to Cauchy stress
                 skk(1,1)=ss(1)
                 skk(2,2)=ss(2)
                 skk(3,3)=ss(3)
                 skk(1,2)=ss(4)
                 skk(1,3)=ss(6)
                 skk(2,3)=ss(5)
                 skk(2,1)=skk(1,2)
                 skk(3,1)=skk(1,3)
                 skk(3,2)=skk(2,3)

                 do j=1,3
                    def(1,j)=ud(0+j)
                    def(2,j)=ud(3+j)
                    def(3,j)=ud(6+j)
                 enddo
                 def(1,1)=def(1,1)+1.0
                 def(2,2)=def(2,2)+1.0
                 def(3,3)=def(3,3)+1.0

                 zjj =def(1,1)*(def(2,2)*def(3,3)-def(3,2)*def(2,3))
     &               -def(1,2)*(def(2,1)*def(3,3)-def(3,1)*def(2,3))
     &               +def(1,3)*(def(2,1)*def(3,2)-def(3,1)*def(2,2))

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
                       scc(i,j)=sum/zjj
                    enddo
                 enddo

                 !Ye,store Cauchy stress
                 cauchy(n,ip_n,1)=scc(1,1)
                 cauchy(n,ip_n,2)=scc(2,2)
                 cauchy(n,ip_n,3)=scc(3,3)
                 cauchy(n,ip_n,4)=scc(1,2)
                 cauchy(n,ip_n,5)=scc(2,3)
                 cauchy(n,ip_n,6)=scc(1,3)

*                 write(iout,82) (ss(j), j=1,4)
*                 s1 = ss(1)+ss(2)+ss(3)
*                 s2 = ss(1)*ss(2) +ss(2)*ss(3) +ss(3)*ss(1) 
*    &                -ss(4)**2 - ss(5)**2 - ss(6)**2
*                 svm = sqrt(s1*s1 - 3.0*s2)
c                 write(iout,82) (ss(j), j=1,6),svm
*                 write(iout,*) '@@ SvM: ',ii,svm
*          write(*,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (dk(i,j), j=1,24)
*          enddo
c            do k=1,6
c               stress(k,ip_n)=ss(k)
c               strain(k,ip_n)=ee(k)
c            enddo
ccccccccccccccccccccccccccccccccccccc
c
c                 FORCES
c                 write(iout,*) ' BE ',lx,ly,lz
c                 do i=1,6
c                    write(iout,82) (BE(i,j), j=1,24)
c                 enddo
c                 do i=1,6
c                    ss(i)=stress(i,ip_n)
c                 enddo
c
c                 add contib to force
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
*                 write(iout,*) ' wt ',wt
                  do 70 j=1,kmax*3
                     sum = 0.0
                     do k=1,6
                        sum = sum + ss(k)*BE(k,j)
                     enddo
                     ff(j) = ff(j) + sum*wt  !Ye, summation for 27 interpolation pts
                     df(j) =         sum*wt
 70               continue
*                   write(iout,* ) ' part dF ',ip_n    
*                   write(iout,82) (df(j), j=1,60)
*                   write(iout,82) (ff(j), j=1,60)
 80      continue
c
         nn=0
         do j=1,kmax
            do i=1,3
               nn=nn+1
               force(i,j)=ff(nn)
            enddo
         enddo
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
c     CONVERT displacements to stresses etc
      subroutine stress_HEX20(neltype,ired_int,
     &                            xyz,uvw,dd,
     &                   stress,inonlinear )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
*        common /gauss_int/ xg,wgt
c
         real*8  u(60),xyz(3,20),uvw(3,20)
         real*8 dd(6,6),BE(6,60),Bd(9,60)
         real*8 xg(4,4),wgt(4,4),ss(8),ee(8),ud(9),en(6)
         real*8 stress(6,27), strain(6,27)
c        integer*4 icoord(8,3)
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
c        data icoord/1,2,2,1,1,2,2,1
c    &              ,1,1,2,2,1,1,2,2
c    &              ,1,1,1,1,2,2,2,2 /
c
        iout=23
         kmax=neltype-1
         nn=0
         do i=1,kmax
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
*                  write(iout,*) nn,u(nn)
            enddo
         enddo
*                write(iout,*) ' D  '
*                do i=1,6
*                   write(iout,82) (dd(i,j), j=1,6)
*                enddo
c
*        nint=3
*        nint=2
*        nint=3
         nint=ired_int
ctag2
         ip_n=0
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8
*           irr=(-1)**lx
*           iss=(-1)**ly
*           itt=(-1)**lz
            ip_n=ip_n+1
*           write(iout,*) irr,iss,itt
*            lx = icoord(ii,1)
*            ly = icoord(ii,2)
*            lz = icoord(ii,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
c
c                 deriv oper and jacobian
                  call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
*                 write(iout,*) ' BE ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (BE(i,j), j=1,60)
*                 enddo
*                 write(iout,*) ' [Bd] '
*                 do i=1,9
*                    write(iout,82) (Bd(i,j), j=1,60)
*                 enddo
         do i=1,9
            sum=0.0
            do k=1,kmax*3
               sum=sum+Bd(i,k)*u(k)
            enddo
            ud(i)=sum
         enddo
*                 write(iout,*) ' u,x  ',ip_n
*                    write(iout,82) ( u(j), j=1,kmax*3)
*                    write(iout,82) ( ud(j), j=1,9)
         z9=1.0
*        z9=0.5
         z9=1.0
         if (inonlinear .eq. 0) z9=0
*        if (iter .eq. 1) then
*            z9=0.0
*        endif

         !Ye,Lagrangian Strain tensor
         ee(1)= ud(1) + z9*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
         ee(2)= ud(5) + z9*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
         ee(3)= ud(9) + z9*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
         ee(4)= ud(2)+ud(4) + z9*(ud(1)*ud(2)+ud(4)*ud(5) +ud(7)*ud(8))
         ee(5)= ud(6)+ud(8) + z9*(ud(2)*ud(3)+ud(5)*ud(6) +ud(8)*ud(9))
         ee(6)= ud(3)+ud(7) + z9*(ud(1)*ud(3)+ud(4)*ud(6) +ud(7)*ud(9))
c
         z8=0.0
         en(1)= ud(1) + z8*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
         en(2)= ud(5) + z8*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
         en(3)= ud(9) + z8*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
         en(4)= ud(2)+ud(4) + z8*(ud(1)*ud(2)+ud(4)*ud(5) +ud(7)*ud(8))
         en(5)= ud(6)+ud(8) + z8*(ud(2)*ud(3)+ud(5)*ud(6) +ud(8)*ud(9))
         en(6)= ud(3)+ud(7) + z8*(ud(1)*ud(3)+ud(4)*ud(6) +ud(7)*ud(9))
c
*                       write(iout,* ) ' strain N L d ',iter
*                    do j=1,6
*                       write(iout,82) ee(j),en(j),en(j)-ee(j)
*                    enddo
*                    do j=1,6
*                       sum = 0.0  
*                       do kk=1,24
*                          sum = sum + bb(j,kk)*u(kk)
*                       enddo
*                       ee0 = sum
*                       write(iout,82) ee(j),ee0 
*                    enddo
*                 write(iout,82) (ee(j)*1.0e6, j=1,6)
                  do i=1,6
                     sum = 0.0  
                     do kk=1,6
                        sum = sum + dd(i,kk)*ee(kk)
                     enddo
                     ss(i)=sum
                  enddo
*                 write(iout,82) (ss(j), j=1,6)
*                 s1 = ss(1)+ss(2)+ss(3)
*                 s2 = ss(1)*ss(2) +ss(2)*ss(3) +ss(3)*ss(1) 
*    &                -ss(4)**2 - ss(5)**2 - ss(6)**2
*                 svm = sqrt(s1*s1 - 3.0*s2)
c                 write(iout,82) (ss(j), j=1,6),svm
*                 write(iout,*) '@@ SvM: ',ii,svm
*          write(*,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (dk(i,j), j=1,24)
*          enddo
             do k=1,6
                stress(k,ip_n)=ss(k)
                strain(k,ip_n)=ee(k)
             enddo
 80      continue
*          write(iout,*) ' [stress 1] '
*          do i=1,6
*             write(iout,82) (stress(i,j), j=1,27)
*          enddo
 82      format(1x,70(g12.6,1x))
c
c
      return
      end
c
c
c
      subroutine BEmat20(neltype,xyz,uvw,r,s,t,BE,Bd,det,ud)
c
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),bb(6,60),p(3,20),xj(3,3),xji(3,3)
         real*8          aa(3,20)
         real*8 BE(6,60),BL(6,60),BN(6,60),Bd(9,60),Bu(9,60)
         real*8 uvw(3,20),u(60),ud(9)
c        real*8 bn2(6,60),pp2(3,20)
c        integer*4 indx(3)
c
         iout=23
         kmax=neltype-1
c
         rp = 1.0 + r
         sp = 1.0 + s
         tp = 1.0 + t
         rm = 1.0 - r
         sm = 1.0 - s
         tm = 1.0 - t
c
c        interp fns
*        h(1) = 0.125*rm*sm*tm*(-r-s-t-2)
*        h(2) = 0.125*rp*sm*tm*(+r-s-t-2)
*        h(3) = 0.125*rp*sp*tm*(+r+s-t-2)
*        h(4) = 0.125*rm*sp*tm*(-r+s-t-2)
*          h(5) = 0.125*rm*sm*tp*(-r-s+t-2)
*          h(6) = 0.125*rp*sm*tp*(+r-s+t-2)
*          h(7) = 0.125*rp*sp*tp*(+r+s+t-2)
*          h(8) = 0.125*rm*sp*tp*(-r+s+t-2)
*        h(9 ) = 0.25*(1-r*r)*sm*tm
*        h(11) = 0.25*(1-r*r)*sp*tm
*        h(17) = 0.25*(1-r*r)*sm*tp
*        h(19) = 0.25*(1-r*r)*sp*tp
*          h(10) = 0.25*rp*(1-s*s)*tm
*          h(12) = 0.25*rm*(1-s*s)*tm
*          h(18) = 0.25*rp*(1-s*s)*tp
*          h(20) = 0.25*rm*(1-s*s)*tp
*            h(13) = 0.25*rm*sm*(1-t*t)
*            h(14) = 0.25*rp*sm*(1-t*t)
*            h(15) = 0.25*rp*sp*(1-t*t)
*            h(16) = 0.25*rm*sp*(1-t*t)
c
c        nat coords deriv wrt r
                  do i=1,3
                     do j=1,kmax
                        p(i,j)=0.0
                     enddo
                  enddo
         p(1,1) = -0.125   *sm*tm*(-r-s-t-2)
     &            +0.125*rm*sm*tm*(-1        )
         p(1,2) = +0.125   *sm*tm*(+r-s-t-2)
     &            +0.125*rp*sm*tm*(+1        )
         p(1,3) = +0.125   *sp*tm*(+r+s-t-2)
     &            +0.125*rp*sp*tm*(+1        )
         p(1,4) = -0.125   *sp*tm*(-r+s-t-2)
     &            +0.125*rm*sp*tm*(-1        )
          p(1,5) = -0.125   *sm*tp*(-r-s+t-2)
     &             +0.125*rm*sm*tp*(-1        )
          p(1,6) = +0.125   *sm*tp*(+r-s+t-2)
     &             +0.125*rp*sm*tp*(+1        )
          p(1,7) = +0.125   *sp*tp*(+r+s+t-2)
     &             +0.125*rp*sp*tp*(+1        )
          p(1,8) = -0.125   *sp*tp*(-r+s+t-2)
     &             +0.125*rm*sp*tp*(-1        )
c        nat coords deriv wrt s
         p(2,1) = -0.125*rm   *tm*(-r-s-t-2)
     &            +0.125*rm*sm*tm*(  -1     )
         p(2,2) = -0.125*rp   *tm*(+r-s-t-2)
     &            +0.125*rp*sm*tm*(  -1     )
         p(2,3) = +0.125*rp   *tm*(+r+s-t-2)
     &            +0.125*rp*sp*tm*(  +1     )
         p(2,4) = +0.125*rm   *tm*(-r+s-t-2)
     &            +0.125*rm*sp*tm*(  +1     )
          p(2,5) = -0.125*rm   *tp*(-r-s+t-2)
     &             +0.125*rm*sm*tp*(  -1     )
          p(2,6) = -0.125*rp   *tp*(+r-s+t-2)
     &             +0.125*rp*sm*tp*(  -1     )
          p(2,7) = +0.125*rp   *tp*(+r+s+t-2)
     &             +0.125*rp*sp*tp*(  +1     )
          p(2,8) = +0.125*rm   *tp*(-r+s+t-2)
     &             +0.125*rm*sp*tp*(  +1     )
c        nat coords deriv wrt t
         p(3,1) = -0.125*rm*sm   *(-r-s-t-2)
     &            +0.125*rm*sm*tm*(    -1     )
         p(3,2) = -0.125*rp*sm   *(+r-s-t-2)
     &            +0.125*rp*sm*tm*(    -1     )
         p(3,3) = -0.125*rp*sp   *(+r+s-t-2)
     &            +0.125*rp*sp*tm*(    -1     )
         p(3,4) = -0.125*rm*sp   *(-r+s-t-2)
     &            +0.125*rm*sp*tm*(    -1     )
          p(3,5) = +0.125*rm*sm   *(-r-s+t-2)
     &             +0.125*rm*sm*tp*(    +1     )
          p(3,6) = +0.125*rp*sm   *(+r-s+t-2)
     &             +0.125*rp*sm*tp*(    +1     )
          p(3,7) = +0.125*rp*sp   *(+r+s+t-2)
     &             +0.125*rp*sp*tp*(    +1     )
          p(3,8) = +0.125*rm*sp   *(-r+s+t-2)
     &             +0.125*rm*sp*tp*(    +1     )
         p(1,9 ) = 0.25*( -2*r)*sm*tm
         p(1,11) = 0.25*( -2*r)*sp*tm
         p(1,17) = 0.25*( -2*r)*sm*tp
         p(1,19) = 0.25*( -2*r)*sp*tp
           p(1,10) = 0.25*(+1)*(1-s*s)*tm
           p(1,12) = 0.25*(-1)*(1-s*s)*tm
           p(1,18) = 0.25*(+1)*(1-s*s)*tp
           p(1,20) = 0.25*(-1)*(1-s*s)*tp
             p(1,13) = 0.25*(-1)*sm*(1-t*t)
             p(1,14) = 0.25*(+1)*sm*(1-t*t)
             p(1,15) = 0.25*(+1)*sp*(1-t*t)
             p(1,16) = 0.25*(-1)*sp*(1-t*t)
         p(2,9 ) = 0.25*(1-r*r)*(-1)*tm
         p(2,11) = 0.25*(1-r*r)*(+1)*tm
         p(2,17) = 0.25*(1-r*r)*(-1)*tp
         p(2,19) = 0.25*(1-r*r)*(+1)*tp
           p(2,10) = 0.25*rp*( -2*s)*tm
           p(2,12) = 0.25*rm*( -2*s)*tm
           p(2,18) = 0.25*rp*( -2*s)*tp
           p(2,20) = 0.25*rm*( -2*s)*tp
             p(2,13) = 0.25*rm*(-1)*(1-t*t)
             p(2,14) = 0.25*rp*(-1)*(1-t*t)
             p(2,15) = 0.25*rp*(+1)*(1-t*t)
             p(2,16) = 0.25*rm*(+1)*(1-t*t)
         p(3,9 ) = 0.25*(1-r*r)*sm*(-1)
         p(3,11) = 0.25*(1-r*r)*sp*(-1)
         p(3,17) = 0.25*(1-r*r)*sm*(+1)
         p(3,19) = 0.25*(1-r*r)*sp*(+1)
           p(3,10) = 0.25*rp*(1-s*s)*(-1)
           p(3,12) = 0.25*rm*(1-s*s)*(-1)
           p(3,18) = 0.25*rp*(1-s*s)*(+1)
           p(3,20) = 0.25*rm*(1-s*s)*(+1)
             p(3,13) = 0.25*rm*sm*( -2*t)
             p(3,14) = 0.25*rp*sm*( -2*t)
             p(3,15) = 0.25*rp*sp*( -2*t)
             p(3,16) = 0.25*rm*sp*( -2*t)
*                 write(iout,*) ' [P] '   
*                 do i=1,3
*                    write(iout,82) (p(i,j), j=1,20)
*                 enddo
c
c        interp fns
*        h(1)     =  0.125*rm  *sm  *tm  *(-r-s-t-2)
*        pp2(1,1) =  0.125*(-1)*sm  *tm  *(-r-s-t-2)
*    &            +  0.125*rm  *sm  *tm  *(-1        )
*        pp2(2,1) =  0.125*rm  *(-1)*tm  *(-r-s-t-2)
*    &            +  0.125*rm  *sm  *tm  *(  -1      )
*        pp2(3,1) =  0.125*rm  *sm  *(-1)*(-r-s-t-2)
*    &            +  0.125*rm  *sm  *tm  *(    -1    )
*        h(2) = 0.125*rp*sm*tm*(+r-s-t-2)
*        h(2)     =  0.125*rp  *sm  *tm  *(+r-s-t-2)
*        pp2(1,2) =  0.125*(+1)*sm  *tm  *(+r-s-t-2)
*    &            +  0.125*rp  *sm  *tm  *(+1        )
*        pp2(2,2) =  0.125*rp  *(-1)*tm  *(+r-s-t-2)
*    &            +  0.125*rp  *sm  *tm  *(  -1      )
*        pp2(3,2) =  0.125*rp  *sm  *(-1)*(+r-s-t-2)
*    &            +  0.125*rp  *sm  *tm  *(    -1    )
*        h(3)     =  0.125*rp  *sp  *tm  *(+r+s-t-2)
*        pp2(1,3) =  0.125*(+1)*sp  *tm  *(+r+s-t-2)
*    &             + 0.125*rp  *sp  *tm  *(+1      )
*        pp2(2,3) =  0.125*rp  *(+1)*tm  *(+r+s-t-2)
*    &             + 0.125*rp  *sp  *tm  *(  +1    )
*        pp2(3,3) =  0.125*rp  *sp  *(-1)*(+r+s-t-2)
*    &             + 0.125*rp  *sp  *tm  *(    -1  )
*        h(4)     = 0.125*rm  *sp  *tm  *(-r+s-t-2)
*        pp2(1,4) = 0.125*(-1)*sp  *tm  *(-r+s-t-2)
*    &             +0.125*rm  *sp  *tm  *(-1      )
*        pp2(2,4) = 0.125*rm  *(+1)*tm  *(-r+s-t-2)
*    &             +0.125*rm  *sp  *tm  *(  +1    )
*        pp2(3,4) = 0.125*rm  *sp  *(-1)*(-r+s-t-2)
*    &             +0.125*rm  *sp  *tm  *(    -1  )
*          h(5)     = 0.125*rm  *sm  *tp  *(-r-s+t-2)
*          pp2(1,5) = 0.125*(-1)*sm  *tp  *(-r-s+t-2)
*    &               +0.125*rm  *sm  *tp  *(-1      )
*          pp2(2,5) = 0.125*rm  *(-1)*tp  *(-r-s+t-2)
*    &               +0.125*rm  *sm  *tp  *(  -1    )
*          pp2(3,5) = 0.125*rm  *sm  *(+1)*(-r-s+t-2)
*    &               +0.125*rm  *sm  *tp  *(    +1  )
*          h(6)     = 0.125*rp  *sm  *tp  *(+r-s+t-2)
*          pp2(1,6) = 0.125*(+1)*sm  *tp  *(+r-s+t-2)
*    &               +0.125*rp  *sm  *tp  *(+1      )
*          pp2(2,6) = 0.125*rp  *(-1)*tp  *(+r-s+t-2)
*    &               +0.125*rp  *sm  *tp  *(  -1    )
*          pp2(3,6) = 0.125*rp  *sm  *(+1)*(+r-s+t-2)
*    &               +0.125*rp  *sm  *tp  *(    +1  )
*          h(7)     = 0.125*rp  *sp  *tp  *(+r+s+t-2)
*          pp2(1,7) = 0.125*(+1)*sp  *tp  *(+r+s+t-2)
*    &               +0.125*rp  *sp  *tp  *(+1      )
*          pp2(2,7) = 0.125*rp  *(+1)*tp  *(+r+s+t-2)
*    &               +0.125*rp  *sp  *tp  *(  +1    )
*          pp2(3,7) = 0.125*rp  *sp  *(+1)*(+r+s+t-2)
*    &               +0.125*rp  *sp  *tp  *(    +1  )
*          h(8)     = 0.125*rm  *sp  *tp  *(-r+s+t-2)
*          pp2(1,8) = 0.125*(-1)*sp  *tp  *(-r+s+t-2)
*    &               +0.125*rm  *sp  *tp  *(-1      )
*          pp2(2,8) = 0.125*rm  *(+1)*tp  *(-r+s+t-2)
*    &               +0.125*rm  *sp  *tp  *(  +1    )
*          pp2(3,8) = 0.125*rm  *sp  *(+1)*(-r+s+t-2)
*    &               +0.125*rm  *sp  *tp  *(    +1  )
*        h(9 )    = 0.25*(1-r*r)*sm  *tm
*        pp2(1,9) = 0.25*( -2*r)*sm  *tm
*        pp2(2,9) = 0.25*(1-r*r)*(-1)*tm
*        pp2(3,9) = 0.25*(1-r*r)*sm  *(-1)
*        h(11)    = 0.25*(1-r*r)*sp  *tm
*        pp2(1,11)= 0.25*( -2*r)*sp  *tm
*        pp2(2,11)= 0.25*(1-r*r)*(+1)*tm
*        pp2(3,11)= 0.25*(1-r*r)*sp  *(-1)
*        h(17)    = 0.25*(1-r*r)*sm  *tp
*        pp2(1,17)= 0.25*( -2*r)*sm  *tp
*        pp2(2,17)= 0.25*(1-r*r)*(-1)*tp
*        pp2(3,17)= 0.25*(1-r*r)*sm  *(+1)
*        h(19)    = 0.25*(1-r*r)*sp  *tp
*        pp2(1,19)= 0.25*( -2*r)*sp  *tp
*        pp2(2,19)= 0.25*(1-r*r)*(+1)*tp
*        pp2(3,19)= 0.25*(1-r*r)*sp  *(+1)
*          h(10)    = 0.25*rp  *(1-s*s)*tm
*          pp2(1,10)= 0.25*(+1)*(1-s*s)*tm
*          pp2(2,10)= 0.25*rp  *( -2*s)*tm
*          pp2(3,10)= 0.25*rp  *(1-s*s)*(-1)
*          h(12)    = 0.25*rm  *(1-s*s)*tm
*          pp2(1,12)= 0.25*(-1)*(1-s*s)*tm
*          pp2(2,12)= 0.25*rm  *( -2*s)*tm
*          pp2(3,12)= 0.25*rm  *(1-s*s)*(-1)
*          h(18)    = 0.25*rp  *(1-s*s)*tp
*          pp2(1,18)= 0.25*(+1)*(1-s*s)*tp
*          pp2(2,18)= 0.25*rp  *( -2*s)*tp
*          pp2(3,18)= 0.25*rp  *(1-s*s)*(+1)
*          h(20)    = 0.25*rm  *(1-s*s)*tp
*          pp2(1,20)= 0.25*(-1)*(1-s*s)*tp
*          pp2(2,20)= 0.25*rm  *( -2*s)*tp
*          pp2(3,20)= 0.25*rm  *(1-s*s)*(+1)
*            h(13)    = 0.25*rm  *sm  *(1-t*t)
*            pp2(1,13)= 0.25*(-1)*sm  *(1-t*t)
*            pp2(2,13)= 0.25*rm  *(-1)*(1-t*t)
*            pp2(3,13)= 0.25*rm  *sm  *( -2*t)
*            h(14)    = 0.25*rp  *sm  *(1-t*t)
*            pp2(1,14)= 0.25*(+1)*sm  *(1-t*t)
*            pp2(2,14)= 0.25*rp  *(-1)*(1-t*t)
*            pp2(3,14)= 0.25*rp  *sm  *( -2*t)
*            h(15)    = 0.25*rp  *sp  *(1-t*t)
*            pp2(1,15)= 0.25*(+1)*sp  *(1-t*t)
*            pp2(2,15)= 0.25*rp  *(+1)*(1-t*t)
*            pp2(3,15)= 0.25*rp  *sp  *( -2*t)
*            h(16)    = 0.25*rm  *sp  *(1-t*t)
*            pp2(1,16)= 0.25*(-1)*sp  *(1-t*t)
*            pp2(2,16)= 0.25*rm  *(+1)*(1-t*t)
*            pp2(3,16)= 0.25*rm  *sp  *( -2*t)
*
*                 write(iout,*) ' [PP2-PP] '   
*                 do i=1,3
*                    write(iout,82) (pp2(i,k)-p(i,k), k=1,20)
*                 enddo
c
c        jacob at (r,s,t)
         do 30 i=1,3
            do 30 j=1,3
               dum = 0.0
               do 20 k=1,kmax
                  dum = dum + p(i,k)*xyz(j,k)
 20            continue
               xj(i,j)=dum
 30      continue
*                 write(iout,*) ' J '   
*                 do i=1,3
*                    write(iout,82) (xj(i,j), j=1,3)
*                 enddo
         det = xj(1,1)*( xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3) )
     &        -xj(1,2)*( xj(2,1)*xj(3,3) - xj(3,1)*xj(2,3) )
     &        +xj(1,3)*( xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2) )
*                 write(iout,*) ' det ',det   
c        jacob inverse
                 xji(1,1) = (xj(2,2) * xj(3,3) - xj(2,3) * xj(3,2))/det
                 xji(1,2) = (xj(3,2) * xj(1,3) - xj(3,3) * xj(1,2))/det
                 xji(1,3) = (xj(1,2) * xj(2,3) - xj(1,3) * xj(2,2))/det
                 xji(2,1) = (xj(2,3) * xj(3,1) - xj(2,1) * xj(3,3))/det
                 xji(2,2) = (xj(1,1) * xj(3,3) - xj(1,3) * xj(3,1))/det
                 xji(2,3) = (xj(1,3) * xj(2,1) - xj(1,1) * xj(2,3))/det
                 xji(3,1) = (xj(2,1) * xj(3,2) - xj(2,2) * xj(3,1))/det
                 xji(3,2) = (xj(1,2) * xj(3,1) - xj(1,1) * xj(3,2))/det
                 xji(3,3) = (xj(1,1) * xj(2,2) - xj(1,2) * xj(2,1))/det
c
*                 write(iout,*) ' J^-1 explicit '   
*                 do i=1,3
*                    write(iout,82) (xji(i,j), j=1,3)
*                 enddo
*        if (det .gt. 0.00000001) goto 40
         if (det .gt. 1.0D-16) goto 40
             write(ilog,*)'@@ ERROR det < 0.0 '
             write(ilog,*)'@@ det = ',det
         stop
c
 40      continue
c        deriv oper B
         do i=1,3
            do j=1,kmax
               sum=0.0d0
               do k=1,3 
                  sum = sum + xji(i,k)*p(k,j)
               enddo
               aa(i,j)=sum
            enddo
         enddo
*                write(iout,*) ' [aa]  '   
*                do i=1,3
*                   write(iout,82) (aa(i,j), j=1,kmax)
*                enddo
c
         do i=1,6
            do j=1,kmax*3
               bb(i,j) = 0.0d0
            enddo
         enddo
         k3=0
         do k=1,kmax
            k3 = k3 + 3
            bb(1,k3-2) = aa(1,k)
            bb(2,k3-1) = aa(2,k)
            bb(3,k3-0) = aa(3,k)
                  bb(4,k3-2) = aa(2,k)
                  bb(4,k3-1) = aa(1,k)
              bb(5,k3-1) = aa(3,k)
              bb(5,k3  ) = aa(2,k)
                bb(6,k3-2) = aa(3,k)
                bb(6,k3  ) = aa(1,k)
         enddo
c    
         do i=1,6
            do j=1,kmax*3
               BL(i,j) = 0.0d0
            enddo
         enddo
         k3=0
         do k=1,kmax
            k3 = k3 + 3
            BL(1,k3-2) = aa(1,k)
            BL(2,k3-1) = aa(2,k)
            BL(3,k3-0) = aa(3,k)
                  BL(4,k3-2) = aa(2,k)
                  BL(4,k3-1) = aa(1,k)
              BL(5,k3-1) = aa(3,k)
              BL(5,k3  ) = aa(2,k)
                BL(6,k3-2) = aa(3,k)
                BL(6,k3  ) = aa(1,k)
         enddo
c    
         do i=1,9
            do j=1,kmax*3
               Bd(i,j) = 0.0d0
            enddo
         enddo
         k3=0
         do k=1,kmax
            k3 = k3 + 3
            Bd(1,k3-2) = aa(1,k)
            Bd(2,k3-2) = aa(2,k)
            Bd(3,k3-2) = aa(3,k)
                  Bd(4,k3-1) = aa(1,k)
                  Bd(5,k3-1) = aa(2,k)
                  Bd(6,k3-1) = aa(3,k)
              Bd(7,k3-0) = aa(1,k)
              Bd(8,k3-0) = aa(2,k)
              Bd(9,k3-0) = aa(3,k)
         enddo
*                write(iout,*) ' [Bd]  '   
*                do i=1,9
*                   write(iout,82) (Bd(i,j), j=1,kmax*3)
*                enddo
c
         do i=1,6
            do j=1,9
               Bu(i,j) = 0.0d0
            enddo
         enddo
         nn=0
         do i=1,kmax
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
            enddo
         enddo
c
*        write(iout,*)' u'
*        write(iout,82) (u(k),k=1,kmax*3)
         do i=1,9
            sum=0.0
            do k=1,kmax*3
               sum=sum+Bd(i,k)*u(k)
*                 write(iout,*) sum,i,k
            enddo
            ud(i)=sum
         enddo
c
*        write(iout,*)' ud vs sum'
*        do j=1,3
*        do i=1,3
*           sum=0.0
*           do k=1,kmax
*              k3=k*3-3+j
*              sum=sum+aa(i,k)*u(k3)
*           enddo
*           write(iout,82) ud((j-1)*3+i),sum
*        enddo
*        enddo
c
         k3=0
         do k=1,3
            k3 = k3 + 3
            Bu(1,k3-2)=ud(1+k3-3)
            Bu(2,k3-1)=ud(2+k3-3)
            Bu(3,k3-0)=ud(3+k3-3)
              Bu(4,k3-2)=ud(2+k3-3)
              Bu(4,k3-1)=ud(1+k3-3)
                Bu(5,k3-1)=ud(3+k3-3)
                Bu(5,k3-0)=ud(2+k3-3)
                  Bu(6,k3-2)=ud(3+k3-3)
                  Bu(6,k3-0)=ud(1+k3-3)
          enddo
         do i=1,6
            do j=1,kmax*3
               sum=0.0
               do k=1,9
                  sum=sum+Bu(i,k)*Bd(k,j)
               enddo
               BN(i,j)=sum
            enddo
         enddo
*           write(iout,* ) 'test BN ' 
*           write(iout,82) (bN(1,  k),ud((k-1)*3+1)*aa(1,1),k=1,3)
*           write(iout,82) (bN(1,3+k),ud((k-1)*3+1)*aa(1,2),k=1,3)
*           write(iout,82) (bN(1,6+k),ud((k-1)*3+1)*aa(1,3),k=1,3)
*           write(iout,82) (bN(2,  k),ud((k-1)*3+2)*aa(2,1),k=1,3)
*    &                    ,(bN(2,3+k),ud((k-1)*3+2)*aa(2,2),k=1,3)
*    &                    ,(bN(2,6+k),ud((k-1)*3+2)*aa(2,3),k=1,3)
*           write(iout,82) bN(2,1),ud(2)*aa(2,1)
*           write(iout,82) bN(3,1),ud(3)*aa(3,1)
*           write(iout,82) (bN(3,  k),ud((k-1)*3+3)*aa(3,1),k=1,3)
*    &                    ,(bN(3,3+k),ud((k-1)*3+3)*aa(3,2),k=1,3)
*    &                    ,(bN(3,6+k),ud((k-1)*3+3)*aa(3,3),k=1,3)
*           write(iout,82) bN(4,1),(ud(1)*aa(2,1)+ud(2)*aa(1,1))
*           write(iout,82) bN(5,1),(ud(2)*aa(3,1)+ud(3)*aa(2,1))
*       write(iout,82) ((bN(5,(j-1)*3+k),(ud((k-1)*3+2)*aa(3,j)
*    &                          +ud((k-1)*3+3)*aa(2,j)),k=1,3),j=1,3)
*           write(iout,82) bN(6,1),(ud(1)*aa(3,1)+ud(3)*aa(1,1))
*       write(iout,82) (bN(6,k),(ud((k-1)*3+1)*aa(3,1)
*    &                          +ud((k-1)*3+3)*aa(1,1)),k=1,3)
*           write(iout,* ) ' ' 
*           do k=1,20
*              k3=(k-1)*3
*              do j=1,3
*                 bn2(1,k3+j) = ud((j-1)*3+1)*aa(1,k)
*                 bn2(2,k3+j) = ud((j-1)*3+2)*aa(2,k)
*                 bn2(3,k3+j) = ud((j-1)*3+3)*aa(3,k)
*                 bN2(4,k3+j) =  ud((j-1)*3+1)*aa(2,k)
*    &                          +ud((j-1)*3+2)*aa(1,k)
*                 bN2(5,k3+j) =  ud((j-1)*3+2)*aa(3,k)
*    &                          +ud((j-1)*3+3)*aa(2,k)
*                 bN2(6,k3+j) =  ud((j-1)*3+1)*aa(3,k)
*    &                          +ud((j-1)*3+3)*aa(1,k)
*              enddo
*           enddo
c    
         z9=1.0
*        z9=0.5
         z9=1.0
         do i=1,6
            do j=1,kmax*3
               BE(i,j) = BL(i,j) + z9*BN(i,j) 
            enddo
         enddo
c
*                write(iout,*) ' [BL] '   
*                do i=1,6
*                   write(iout,82) (BL(i,j), j=1,60)
*                enddo
*                write(iout,*) ' [BN] '   
*                do i=1,6
*                   write(iout,82) (BN(i,j), j=1,60)
*                enddo
*                write(iout,*) ' [BN2-BN] '   
*                do i=1,6
*                   write(iout,82) (BN2(i,j)-bn(i,j), j=1,60)
*                enddo
*        do i=1,6
*           sum=0.0
*           do k=1,24
*              sum=sum+bb(i,k)*u(k)
*           enddo
*           write(iout,*) sum
*        enddo
c
c    
 82      format(1x,70(g12.6,1x))
      return
      end
c
c
c     CONVERT displacements to stresses etc
      subroutine stress_HEX8( 
     &                            xyz,uvw,dd,
     &                   stress )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
*        common /gauss_int/ xg,wgt
c
         real*8  u(60),xyz(3,20),uvw(3,20)
         real*8 dd(6,6),BE(6,24),Bd(9,24)
         real*8 xg(4,4),wgt(4,4),ss(8),ee(8),ud(9)
         real*8 stress(6,27), strain(6,27)
c        integer*4 icoord(8,3)
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
c        data icoord/1,2,2,1,1,2,2,1
c    &              ,1,1,2,2,1,1,2,2
c    &              ,1,1,1,1,2,2,2,2 /
c
        iout=23
         nn=0
         do i=1,8
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
*                  write(iout,*) nn,u(nn)
            enddo
         enddo
c
         nint=2
         nint=2
         ip_n=0
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8
*           irr=(-1)**lx
*           iss=(-1)**ly
*           itt=(-1)**lz
            ip_n=ip_n+1
*           write(iout,*) irr,iss,itt
*            lx = icoord(ii,1)
*            ly = icoord(ii,2)
*            lz = icoord(ii,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
c
c                 deriv oper and jacobian
                  call BEmat8(xyz,uvw,ri,si,ti,BE,Bd,det)
*                 write(iout,*) ' BE ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (BE(i,j), j=1,24)
*                 enddo
*                 write(iout,*) ' Bd '
*                 do i=1,9
*                    write(iout,82) (Bd(i,j), j=1,24)
*                 enddo
         do i=1,9
            sum=0.0
            do k=1,24
               sum=sum+Bd(i,k)*u(k)
            enddo
            ud(i)=sum
         enddo
                  write(iout,*) ' u  '
                     write(iout,82) ( u(j), j=1,24)
                     write(iout,82) ( ud(j), j=1,9)
         z9=1.0
         ee(1)= ud(1) + z9*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
         ee(2)= ud(5) + z9*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
         ee(3)= ud(9) + z9*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
         ee(4)= ud(2)+ud(4) + z9*(ud(1)*ud(2)+ud(4)*ud(5) +ud(7)*ud(8))
         ee(5)= ud(6)+ud(8) + z9*(ud(2)*ud(3)+ud(5)*ud(6) +ud(8)*ud(9))
         ee(6)= ud(3)+ud(7) + z9*(ud(1)*ud(3)+ud(4)*ud(6) +ud(7)*ud(9))
c
                        write(iout,* ) ' strain '
                     do j=1,6
                        write(iout,82) ee(j) 
                     enddo
*                    do j=1,6
*                       sum = 0.0  
*                       do kk=1,24
*                          sum = sum + bb(j,kk)*u(kk)
*                       enddo
*                       ee0 = sum
*                       write(iout,82) ee(j),ee0 
*                    enddo
*                 write(iout,82) (ee(j)*1.0e6, j=1,6)
                  do i=1,6
                     sum = 0.0  
                     do kk=1,6
                        sum = sum + dd(i,kk)*ee(kk)
                     enddo
                     ss(i)=sum
                  enddo
*                 write(iout,82) (ss(j), j=1,6)
*                 s1 = ss(1)+ss(2)+ss(3)
*                 s2 = ss(1)*ss(2) +ss(2)*ss(3) +ss(3)*ss(1) 
*    &                -ss(4)**2 - ss(5)**2 - ss(6)**2
*                 svm = sqrt(s1*s1 - 3.0*s2)
c                 write(iout,82) (ss(j), j=1,6),svm
*                 write(iout,*) '@@ SvM: ',ii,svm
*          write(*,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (dk(i,j), j=1,24)
*          enddo
             do k=1,6
                stress(k,ip_n)=ss(k)
                strain(k,ip_n)=ee(k)
             enddo
 80      continue
 82      format(1x,30(g12.6,1x))
c
c
      return
      end
c
c
c     Body FORCE from stress for HEXahedron 8 nodes
      subroutine bforce_HEX8( xyz, uvw,stress,force)
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),uvw(3,20)
         real*8 stress(6,27), force(3,20),ff(60),ss(6)
c
         real*8 BE(6,24),Bd(9,24)
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
c
*          write(*,*)' [x y z] '
*          do i=1,3
*             write(*,82) (xyz(i,j), j=1,8)
*          enddo
*          write(*,*) ' '
c
         do i=1,24
            ff(i) = 0.0 
         enddo
         ip_n=0
         nint=2
c        nint=1
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8   
         ip_n = ip_n + 1
*            lx = icoord(ip_n,1)
*            ly = icoord(ip_n,2)
*            lz = icoord(ip_n,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
*               write(*,84) lx,ly,lz,ip_n,(icoord(ip_n,kk),kk=1,3)
c
c                 deriv oper and jacobian
                  call BEmat8(xyz,uvw,ri,si,ti,BE,Bd,det)
*                 write(iout,*) ' BE ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (BE(i,j), j=1,24)
*                 enddo
                  do i=1,6
                     ss(i)=stress(i,ip_n)
                  enddo
c
c                 add contib to force
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
*                 write(iout,*) ' wt ',wt
                  do 70 j=1,24
                     sum = 0.0
                     do k=1,6
                        sum = sum + ss(k)*BE(k,j)
                     enddo
                     ff(j) = ff(j) + sum*wt
 70               continue
 80      continue
c
         nn=0
         do j=1,8
            do i=1,3
               nn=nn+1
               force(i,j)=ff(nn)
            enddo
         enddo
c
 82        format(1x,40(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
c
      subroutine BEmat8(xyz,uvw,r,s,t,BE,Bd,det)
c
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,8),bb(6,24),p(3,8),xj(3,3),xji(3,3)
         real*8          aa(3,8)
         real*8 BE(6,24),BL(6,24),BN(6,24),Bd(9,24),Bu(9,24)
         real*8 uvw(3,20),u(24),ud(9)
*        real*8 h(8)
c        integer*4 indx(3)
c
         iout=23
c
         rp = 1.0 + r
         sp = 1.0 + s
         tp = 1.0 + t
         rm = 1.0 - r
         sm = 1.0 - s
         tm = 1.0 - t
c
c        interp fns
*        h(1) = 0.125*rm*sm*tm
*        h(2) = 0.125*rp*sm*tm
*        h(3) = 0.125*rp*sp*tm
*        h(4) = 0.125*rm*sp*tm
*          h(5) = 0.125*rm*sm*tp
*          h(6) = 0.125*rp*sm*tp
*          h(7) = 0.125*rp*sp*tp
*          h(8) = 0.125*rm*sp*tp
c
c        nat coords deriv wrt r
         p(1,1) = -0.125   *sm*tm
         p(1,2) =  0.125   *sm*tm
         p(1,3) =  0.125   *sp*tm
         p(1,4) = -0.125   *sp*tm
          p(1,5) = -0.125   *sm*tp
          p(1,6) =  0.125   *sm*tp
          p(1,7) =  0.125   *sp*tp
          p(1,8) = -0.125   *sp*tp
c        nat coords deriv wrt s
         p(2,1) = -0.125*rm   *tm
         p(2,2) = -0.125*rp   *tm
         p(2,3) =  0.125*rp   *tm
         p(2,4) =  0.125*rm   *tm
          p(2,5) = -0.125*rm   *tp
          p(2,6) = -0.125*rp   *tp
          p(2,7) =  0.125*rp   *tp
          p(2,8) =  0.125*rm   *tp
c        nat coords deriv wrt t
         p(3,1) = -0.125*rm*sm
         p(3,2) = -0.125*rp*sm
         p(3,3) = -0.125*rp*sp
         p(3,4) = -0.125*rm*sp
          p(3,5) =  0.125*rm*sm
          p(3,6) =  0.125*rp*sm
          p(3,7) =  0.125*rp*sp
          p(3,8) =  0.125*rm*sp
*                 write(*,*) ' [P] ',r,s,t   
*                 do i=1,3
*                    write(*,82) (p(i,j), j=1,8)
*                 enddo
c
c        jacob at (r,s,t)
         do 30 i=1,3
            do 30 j=1,3
               dum = 0.0
               do 20 k=1,8
                  dum = dum + p(i,k)*xyz(j,k)
 20            continue
               xj(i,j)=dum
 30      continue
*                 write(*,*) ' J ',r,s,t   
*                 do i=1,3
*                    write(*,82) (xj(i,j), j=1,3)
*                 enddo
         det = xj(1,1)*( xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3) )
     &        -xj(1,2)*( xj(2,1)*xj(3,3) - xj(3,1)*xj(2,3) )
     &        +xj(1,3)*( xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2) )
*                 write(*,*) ' det ',det   
c        jacob inverse
                 xji(1,1) = (xj(2,2) * xj(3,3) - xj(2,3) * xj(3,2))/det
                 xji(1,2) = (xj(3,2) * xj(1,3) - xj(3,3) * xj(1,2))/det
                 xji(1,3) = (xj(1,2) * xj(2,3) - xj(1,3) * xj(2,2))/det
                 xji(2,1) = (xj(2,3) * xj(3,1) - xj(2,1) * xj(3,3))/det
                 xji(2,2) = (xj(1,1) * xj(3,3) - xj(1,3) * xj(3,1))/det
                 xji(2,3) = (xj(1,3) * xj(2,1) - xj(1,1) * xj(2,3))/det
                 xji(3,1) = (xj(2,1) * xj(3,2) - xj(2,2) * xj(3,1))/det
                 xji(3,2) = (xj(1,2) * xj(3,1) - xj(1,1) * xj(3,2))/det
                 xji(3,3) = (xj(1,1) * xj(2,2) - xj(1,2) * xj(2,1))/det
c
*                 write(*,*) ' J^-1 explicit '   
*                 do i=1,3
*                    write(*,82) (aj(i,j), j=1,3)
*                 enddo
         if (det .gt. 0.00000001) goto 40
         write(ilog,*)'@@ ERROR det  '
         stop
c
 40      continue
c        jacob inverse
c        deriv oper B
         do i=1,3
            do j=1,8
               sum=0.0d0
               do k=1,3 
                  sum = sum + xji(i,k)*p(k,j)
               enddo
               aa(i,j)=sum
            enddo
         enddo
         do i=1,6
            do j=1,24
               bb(i,j) = 0.0d0
            enddo
         enddo
         k3=0
         do k=1,8
            k3 = k3 + 3
            bb(1,k3-2) = aa(1,k)
            bb(2,k3-1) = aa(2,k)
            bb(3,k3-0) = aa(3,k)
              bb(5,k3-1) = aa(3,k)
              bb(5,k3  ) = aa(2,k)
                bb(6,k3-2) = aa(3,k)
                bb(6,k3  ) = aa(1,k)
                  bb(4,k3-2) = aa(2,k)
                  bb(4,k3-1) = aa(1,k)
         enddo
c    
         do i=1,6
            do j=1,24
               BL(i,j) = 0.0d0
            enddo
         enddo
         k3=0
         do k=1,8
            k3 = k3 + 3
            BL(1,k3-2) = aa(1,k)
            BL(2,k3-1) = aa(2,k)
            BL(3,k3-0) = aa(3,k)
                  BL(4,k3-2) = aa(2,k)
                  BL(4,k3-1) = aa(1,k)
              BL(5,k3-1) = aa(3,k)
              BL(5,k3  ) = aa(2,k)
                BL(6,k3-2) = aa(3,k)
                BL(6,k3  ) = aa(1,k)
         enddo
c    
         do i=1,9
            do j=1,24
               Bd(i,j) = 0.0d0
            enddo
         enddo
         k3=0
         do k=1,8
            k3 = k3 + 3
            Bd(1,k3-2) = aa(1,k)
            Bd(2,k3-2) = aa(2,k)
            Bd(3,k3-2) = aa(3,k)
                  Bd(4,k3-1) = aa(1,k)
                  Bd(5,k3-1) = aa(2,k)
                  Bd(6,k3-1) = aa(3,k)
              Bd(7,k3-0) = aa(1,k)
              Bd(8,k3-0) = aa(2,k)
              Bd(9,k3-0) = aa(3,k)
         enddo
c
         do i=1,6
            do j=1,9
               Bu(i,j) = 0.0d0
            enddo
         enddo
         nn=0
         do i=1,8
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
            enddo
         enddo
         do i=1,9
            sum=0.0
            do k=1,24
               sum=sum+Bd(i,k)*u(k)
            enddo
            ud(i)=sum
         enddo
         k3=0
         do k=1,3
            k3 = k3 + 3
            Bu(1,k3-2)=ud(1+k3-3)
            Bu(2,k3-1)=ud(2+k3-3)
            Bu(3,k3-0)=ud(3+k3-3)
              Bu(4,k3-2)=ud(2+k3-3)
              Bu(4,k3-1)=ud(1+k3-3)
                Bu(5,k3-1)=ud(3+k3-3)
                Bu(5,k3-0)=ud(2+k3-3)
                  Bu(6,k3-2)=ud(3+k3-3)
                  Bu(6,k3-0)=ud(1+k3-3)
          enddo
         do  i=1,6
            do j=1,24
               sum=0.0
               do k=1,9
                  sum=sum+Bu(i,k)*Bd(k,j)
               enddo
               BN(i,j)=sum
            enddo
         enddo
c    
         z9=1.0
         do i=1,6
            do j=1,24
               BE(i,j) = BL(i,j) + z9*BN(i,j) 
            enddo
         enddo
c
*                write(iout,*) ' [BB] '   
*                do i=1,6
*                   write(iout,82) (bb(i,j), j=1,24)
*                enddo
*        do i=1,6
*           sum=0.0
*           do k=1,24
*              sum=sum+bb(i,k)*u(k)
*           enddo
*           write(iout,*) sum
*        enddo
c
c    
 82      format(1x,40(g12.6,1x))
      return
      end
c
c
c     ELeMent STiFfness for TETrahedron
      subroutine stress_TET( xyz,uvw,dd,ss )
c
         implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),aa(3,3),bb(6,12),ss(6),dd(6,6),uvw(3,20)
         real*8 ee(6),u(12),ud(9),bd(9,12),BL(6,12)
c
         iout=23
         do i=1,6
            do j=1,12
               bb(i,j)=0.
               BL(i,j)=0.
            enddo
         enddo
         do i=1,9
            do j=1,12
               Bd(i,j)=0.
            enddo
         enddo
c
         x14 = xyz(1,1)-xyz(1,4)
         x24 = xyz(1,2)-xyz(1,4)
         x34 = xyz(1,3)-xyz(1,4)
           y14 = xyz(2,1)-xyz(2,4)
           y24 = xyz(2,2)-xyz(2,4)
           y34 = xyz(2,3)-xyz(2,4)
             z14 = xyz(3,1)-xyz(3,4)
             z24 = xyz(3,2)-xyz(3,4)
             z34 = xyz(3,3)-xyz(3,4)
*        write(*,*)' Jij'
*        write(*,*)  x14,y14,z14
*        write(*,*)  x24,y24,z24
*        write(*,*)  x34,y34,z34
         aa(1,1) = y24*z34 - y34*z24
         aa(2,1) = z24*x34 - z34*x24
         aa(3,1) = x24*y34 - x34*y24
            aa(1,2) = y34*z14 - y14*z34
            aa(2,2) = z34*x14 - z14*x34
            aa(3,2) = x34*y14 - x14*y34
               aa(1,3) = y14*z24 - y24*z14
               aa(2,3) = z14*x24 - z24*x14
               aa(3,3) = x14*y24 - x24*y14
         detj = x14*(y24*z34-y34*z24)
     &        + y14*(z24*x34-z34*x24)
     &        + z14*(x24*y34-x34*y24)
         do i=1,3
            do j=1,3
               aa(i,j) = aa(i,j)/detj
            enddo
         enddo
         vol=abs(detj)/6.0
c
*         write(iout,*)' detJ: ',detj
*        write(*,*)' AA: '
*             do i=1,3
*                write(*,82) (aa(i,j), j=1,3)
*             enddo
c
c          {e}=[B]{u}
         BL(1,1) = aa(1,1) 
         BL(4,1) = aa(2,1)
         BL(6,1) = aa(3,1)
           BL(2,2) = aa(2,1) 
           BL(4,2) = aa(1,1)
           BL(5,2) = aa(3,1)
             BL(3,3) = aa(3,1) 
             BL(5,3) = aa(2,1)
             BL(6,3) = aa(1,1)
         BL(1,4) = aa(1,2) 
         BL(4,4) = aa(2,2)
         BL(6,4) = aa(3,2)
           BL(2,5) = aa(2,2) 
           BL(4,5) = aa(1,2)
           BL(5,5) = aa(3,2)
             BL(3,6) = aa(3,2) 
             BL(5,6) = aa(2,2)
             BL(6,6) = aa(1,2)
         BL(1,7) = aa(1,3) 
         BL(4,7) = aa(2,3)
         BL(6,7) = aa(3,3)
           BL(2,8) = aa(2,3) 
           BL(4,8) = aa(1,3)
           BL(5,8) = aa(3,3)
             BL(3,9) = aa(3,3) 
             BL(5,9) = aa(2,3)
             BL(6,9) = aa(1,3)
                ab1 = aa(1,1)+aa(1,2)+aa(1,3)
                ab2 = aa(2,1)+aa(2,2)+aa(2,3)
                ab3 = aa(3,1)+aa(3,2)+aa(3,3)
         BL(1,10) = -ab1    
         BL(4,10) = -ab2   
         BL(6,10) = -ab3   
           BL(2,11) = -ab2    
           BL(4,11) = -ab1   
           BL(5,11) = -ab3   
             BL(3,12) = -ab3    
             BL(5,12) = -ab2   
             bb(6,12) = -ab1   
*        write(*,*)' BB: '
*             do i=1,6
*                write(*,82) (bb(i,j), j=4,6)
*             enddo
c
c         {e} = [B]{u}  
         nn=0
         do i=1,4
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
*                  write(*,*) nn,u(nn)
            enddo
         enddo
         do i=1,6
               sum=0.0
               do k=1,12
                  sum = sum + bL(i,k)*u(k)
               enddo
               ee(i) = sum
*                   write(iout,*) i ,ee(i),' e ori'
         enddo
c
         nn=0
         do i=1,4
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
            enddo
         enddo
         Bd(1,1) = aa(1,1) 
         Bd(1,4) = aa(1,2) 
         Bd(1,7) = aa(1,3) 
         Bd(1,10) = -ab1    
           Bd(2,1) = aa(2,1) 
           Bd(2,4) = aa(2,2) 
           Bd(2,7) = aa(2,3) 
           Bd(2,10) = -ab2    
             Bd(3,1) = aa(3,1) 
             Bd(3,4) = aa(3,2) 
             Bd(3,7) = aa(3,3) 
             Bd(3,10) = -ab3    
         Bd(4,2) = aa(1,1) 
         Bd(4,5) = aa(1,2) 
         Bd(4,8) = aa(1,3) 
         Bd(4,11) = -ab1    
           Bd(5,2) = aa(2,1) 
           Bd(5,5) = aa(2,2) 
           Bd(5,8) = aa(2,3) 
           Bd(5,11) = -ab2    
             Bd(6,2) = aa(3,1) 
             Bd(6,5) = aa(3,2) 
             Bd(6,8) = aa(3,3) 
             Bd(6,11) = -ab3    
         Bd(7,3) = aa(1,1) 
         Bd(7,6) = aa(1,2) 
         Bd(7,9) = aa(1,3) 
         Bd(7,12) = -ab1    
           Bd(8,3) = aa(2,1) 
           Bd(8,6) = aa(2,2) 
           Bd(8,9) = aa(2,3) 
           Bd(8,12) = -ab2    
             Bd(9,3) = aa(3,1) 
             Bd(9,6) = aa(3,2) 
             Bd(9,9) = aa(3,3) 
             Bd(9,12) = -ab3    
         do i=1,9
            sum=0.0
            do k=1,12
               sum=sum+Bd(i,k)*u(k)
            enddo
            ud(i)=sum
         enddo
*        do i=1,6
*              sum=0.0
*              do k=1,12
*                 sum = sum + bb(i,k)*u(k)
*              enddo
*              ee(i) = sum
*                  write(*,*) i ,ee(i),' e'
*        enddo
         ee(1) = ud(1) + 1*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
         ee(2) = ud(5) + 1*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
         ee(3) = ud(9) + 1*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
         ee(4) = ud(2)+ud(4) + 1*(ud(1)*ud(2)+ud(4)*ud(5) + ud(7)*ud(8))
         ee(5) = ud(6)+ud(8) + 1*(ud(2)*ud(3)+ud(5)*ud(6) + ud(8)*ud(9))
         ee(6) = ud(3)+ud(7) + 1*(ud(1)*ud(3)+ud(4)*ud(6) + ud(7)*ud(9))
         do i=1,6
               sum=0.0
               do k=1,6
                  sum = sum + dd(i,k)*ee(k)
               enddo
               ss(i) = sum
*                   write(iout,*) i ,ee(i),ss(i),' s'
         enddo
c
 82     format(1x,20(g12.6,1x))
      return
      end
c
c
c     ELeMent STiFfness for TETrahedron
      subroutine bforce_TET( xyz,uvw,ss,ff )
c
         implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),aa(3,3),bb(6,12),ss(6),ff(12)
         real*8 uvw(3,20),u(12),ud(9)
         real*8 BL(6,12),BN(6,12),Bu(6,9),Bd(9,12)
c
         do i=1,6
            do j=1,12
               bb(i,j)=0.
            enddo
         enddo
         do i=1,6
            do j=1,12
               BL(i,j)=0.
               BN(i,j)=0.
            enddo
         enddo
         do i=1,6
            do j=1,9
               Bu(i,j)=0.
            enddo
         enddo
         do i=1,9
            do j=1,12
               Bd(i,j)=0.
            enddo
         enddo
c
         x14 = xyz(1,1)-xyz(1,4)
         x24 = xyz(1,2)-xyz(1,4)
         x34 = xyz(1,3)-xyz(1,4)
           y14 = xyz(2,1)-xyz(2,4)
           y24 = xyz(2,2)-xyz(2,4)
           y34 = xyz(2,3)-xyz(2,4)
             z14 = xyz(3,1)-xyz(3,4)
             z24 = xyz(3,2)-xyz(3,4)
             z34 = xyz(3,3)-xyz(3,4)
*        write(*,*)' Jij'
*        write(*,*)  x14,y14,z14
*        write(*,*)  x24,y24,z24
*        write(*,*)  x34,y34,z34
         aa(1,1) = y24*z34 - y34*z24
         aa(2,1) = z24*x34 - z34*x24
         aa(3,1) = x24*y34 - x34*y24
            aa(1,2) = y34*z14 - y14*z34
            aa(2,2) = z34*x14 - z14*x34
            aa(3,2) = x34*y14 - x14*y34
               aa(1,3) = y14*z24 - y24*z14
               aa(2,3) = z14*x24 - z24*x14
               aa(3,3) = x14*y24 - x24*y14
         detj = x14*(y24*z34-y34*z24)
     &        + y14*(z24*x34-z34*x24)
     &        + z14*(x24*y34-x34*y24)
         do i=1,3
            do j=1,3
               aa(i,j) = aa(i,j)/detj
            enddo
         enddo
         vol=abs(detj)/6.0
c
*        write(*,*)' detJ: ',detj
*        write(*,*)' AA: '
*             do i=1,3
*                write(*,82) (aa(i,j), j=1,3)
*             enddo
c
c          {e}=[B]{u}
         BL(1,1) = aa(1,1) 
         BL(4,1) = aa(2,1)
         BL(6,1) = aa(3,1)
           BL(2,2) = aa(2,1) 
           BL(4,2) = aa(1,1)
           BL(5,2) = aa(3,1)
             BL(3,3) = aa(3,1) 
             BL(5,3) = aa(2,1)
             BL(6,3) = aa(1,1)
         BL(1,4) = aa(1,2) 
         BL(4,4) = aa(2,2)
         BL(6,4) = aa(3,2)
           BL(2,5) = aa(2,2) 
           BL(4,5) = aa(1,2)
           BL(5,5) = aa(3,2)
             BL(3,6) = aa(3,2) 
             BL(5,6) = aa(2,2)
             BL(6,6) = aa(1,2)
         BL(1,7) = aa(1,3) 
         BL(4,7) = aa(2,3)
         BL(6,7) = aa(3,3)
           BL(2,8) = aa(2,3) 
           BL(4,8) = aa(1,3)
           BL(5,8) = aa(3,3)
             BL(3,9) = aa(3,3) 
             BL(5,9) = aa(2,3)
             BL(6,9) = aa(1,3)
                ab1 = aa(1,1)+aa(1,2)+aa(1,3)
                ab2 = aa(2,1)+aa(2,2)+aa(2,3)
                ab3 = aa(3,1)+aa(3,2)+aa(3,3)
         BL(1,10) = -ab1    
         BL(4,10) = -ab2   
         BL(6,10) = -ab3   
           BL(2,11) = -ab2    
           BL(4,11) = -ab1   
           BL(5,11) = -ab3   
             BL(3,12) = -ab3    
             BL(5,12) = -ab2   
             BL(6,12) = -ab1   
*        write(*,*)' BB: '
*             do i=1,6
*                write(*,82) (bb(i,j), j=1,12)
*             enddo
         nn=0
         do i=1,4
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
            enddo
         enddo
         Bd(1,1) = aa(1,1) 
         Bd(1,4) = aa(1,2) 
         Bd(1,7) = aa(1,3) 
         Bd(1,10) = -ab1    
           Bd(2,1) = aa(2,1) 
           Bd(2,4) = aa(2,2) 
           Bd(2,7) = aa(2,3) 
           Bd(2,10) = -ab2    
             Bd(3,1) = aa(3,1) 
             Bd(3,4) = aa(3,2) 
             Bd(3,7) = aa(3,3) 
             Bd(3,10) = -ab3    
         Bd(4,2) = aa(1,1) 
         Bd(4,5) = aa(1,2) 
         Bd(4,8) = aa(1,3) 
         Bd(4,11) = -ab1    
           Bd(5,2) = aa(2,1) 
           Bd(5,5) = aa(2,2) 
           Bd(5,8) = aa(2,3) 
           Bd(5,11) = -ab2    
             Bd(6,2) = aa(3,1) 
             Bd(6,5) = aa(3,2) 
             Bd(6,8) = aa(3,3) 
             Bd(6,11) = -ab3    
         Bd(7,3) = aa(1,1) 
         Bd(7,6) = aa(1,2) 
         Bd(7,9) = aa(1,3) 
         Bd(7,12) = -ab1    
           Bd(8,3) = aa(2,1) 
           Bd(8,6) = aa(2,2) 
           Bd(8,9) = aa(2,3) 
           Bd(8,12) = -ab2    
             Bd(9,3) = aa(3,1) 
             Bd(9,6) = aa(3,2) 
             Bd(9,9) = aa(3,3) 
             Bd(9,12) = -ab3    
         do i=1,9
            sum=0.0
            do k=1,12
               sum=sum+Bd(i,k)*u(k)
            enddo
            ud(i)=sum
         enddo
         Bu(1,1)=ud(1)
         Bu(1,4)=ud(4)
         Bu(1,7)=ud(7)
           Bu(2,2)=ud(2)
           Bu(2,5)=ud(5)
           Bu(2,8)=ud(8)
             Bu(3,3)=ud(3)
             Bu(3,6)=ud(6)
             Bu(3,9)=ud(9)
         Bu(4,1)=ud(2)
         Bu(4,2)=ud(1)
         Bu(4,4)=ud(5)
         Bu(4,5)=ud(4)
         Bu(4,7)=ud(8)
         Bu(4,8)=ud(7)
           Bu(5,2)=ud(3)
           Bu(5,3)=ud(2)
           Bu(5,5)=ud(6)
           Bu(5,6)=ud(5)
           Bu(5,8)=ud(9)
           Bu(5,9)=ud(8)
             Bu(6,1)=ud(3)
             Bu(6,3)=ud(1)
             Bu(6,4)=ud(6)
             Bu(6,6)=ud(4)
             Bu(6,7)=ud(9)
             Bu(6,9)=ud(7)
         do i=1,6
            do j=1,12
               sum=0.0
               do k=1,9
                  sum=sum+Bu(i,k)*Bd(k,j)
               enddo
               BN(i,j)=sum
            enddo
         enddo
c
c
c         {F} = {S}[B] vol
         do i=1,12
               sum=0.0
               do k=1,6
                  sum = sum + ss(k)*(BL(k,i)+BN(k,i))
               enddo
               ff(i) = sum*vol
         enddo
c
 82     format(1x,20(g12.6,1x))
      return
      end
c
c
c     ELeMent STiFfness for TETrahedron
      subroutine stiffE_TET( xyz,uvw,dd,ss,ek )
c
         implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),aa(3,3),bb(6,12),ss(6),dd(6,6),uvw(3,20)
         real*8 ee(6),u(12),ud(9),ek(60,60),etemp(6,12)
         real*8 BL(6,12),BN(6,12),Bu(6,9),Bd(9,12)
c
         iout=23
         do i=1,12
            do j=1,12
               ek(i,j)=0.0
            enddo
         enddo
         do i=1,6
            do j=1,12
               etemp(i,j)=0.0
            enddo
         enddo
         do i=1,6
            do j=1,12
               bb(i,j)=0.
               BL(i,j)=0.
               BN(i,j)=0.
            enddo
         enddo
         do i=1,6
            do j=1,9
               Bu(i,j)=0.
            enddo
         enddo
         do i=1,9
            do j=1,12
               Bd(i,j)=0.
            enddo
         enddo
c
         x14 = xyz(1,1)-xyz(1,4)
         x24 = xyz(1,2)-xyz(1,4)
         x34 = xyz(1,3)-xyz(1,4)
           y14 = xyz(2,1)-xyz(2,4)
           y24 = xyz(2,2)-xyz(2,4)
           y34 = xyz(2,3)-xyz(2,4)
             z14 = xyz(3,1)-xyz(3,4)
             z24 = xyz(3,2)-xyz(3,4)
             z34 = xyz(3,3)-xyz(3,4)
*        write(*,*)' Jij'
*        write(*,*)  x14,y14,z14
*        write(*,*)  x24,y24,z24
*        write(*,*)  x34,y34,z34
         aa(1,1) = y24*z34 - y34*z24
         aa(2,1) = z24*x34 - z34*x24
         aa(3,1) = x24*y34 - x34*y24
            aa(1,2) = y34*z14 - y14*z34
            aa(2,2) = z34*x14 - z14*x34
            aa(3,2) = x34*y14 - x14*y34
               aa(1,3) = y14*z24 - y24*z14
               aa(2,3) = z14*x24 - z24*x14
               aa(3,3) = x14*y24 - x24*y14
         detj = x14*(y24*z34-y34*z24)
     &        + y14*(z24*x34-z34*x24)
     &        + z14*(x24*y34-x34*y24)
         do i=1,3
            do j=1,3
               aa(i,j) = aa(i,j)/detj
            enddo
         enddo
         vol=abs(detj)/6.0
c
*        write(*,*)' detJ: ',detj
*        write(*,*)' AA: '
*             do i=1,3
*                write(*,82) (aa(i,j), j=1,3)
*             enddo
c
c          {e}=[B]{u}
         bL(1,1) = aa(1,1) 
         bL(4,1) = aa(2,1)
         bL(6,1) = aa(3,1)
           bL(2,2) = aa(2,1) 
           bL(4,2) = aa(1,1)
           bL(5,2) = aa(3,1)
             bL(3,3) = aa(3,1) 
             bL(5,3) = aa(2,1)
             bL(6,3) = aa(1,1)
         bL(1,4) = aa(1,2) 
         bL(4,4) = aa(2,2)
         bL(6,4) = aa(3,2)
           bL(2,5) = aa(2,2) 
           bL(4,5) = aa(1,2)
           bL(5,5) = aa(3,2)
             bL(3,6) = aa(3,2) 
             bL(5,6) = aa(2,2)
             bL(6,6) = aa(1,2)
         bL(1,7) = aa(1,3) 
         bL(4,7) = aa(2,3)
         bL(6,7) = aa(3,3)
           bL(2,8) = aa(2,3) 
           bL(4,8) = aa(1,3)
           bL(5,8) = aa(3,3)
             bL(3,9) = aa(3,3) 
             bL(5,9) = aa(2,3)
             bL(6,9) = aa(1,3)
                ab1 = aa(1,1)+aa(1,2)+aa(1,3)
                ab2 = aa(2,1)+aa(2,2)+aa(2,3)
                ab3 = aa(3,1)+aa(3,2)+aa(3,3)
         bL(1,10) = -ab1    
         bL(4,10) = -ab2   
         bL(6,10) = -ab3   
           bL(2,11) = -ab2    
           bL(4,11) = -ab1   
           bL(5,11) = -ab3   
             bL(3,12) = -ab3    
             bL(5,12) = -ab2   
             bL(6,12) = -ab1   
*        write(*,*)' BB: '
*             do i=1,6
*                write(*,82) (bb(i,j), j=1,12)
*             enddo
         nn=0
         do i=1,4
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
            enddo
         enddo
         Bd(1,1) = aa(1,1) 
         Bd(1,4) = aa(1,2) 
         Bd(1,7) = aa(1,3) 
         Bd(1,10) = -ab1    
           Bd(2,1) = aa(2,1) 
           Bd(2,4) = aa(2,2) 
           Bd(2,7) = aa(2,3) 
           Bd(2,10) = -ab2    
             Bd(3,1) = aa(3,1) 
             Bd(3,4) = aa(3,2) 
             Bd(3,7) = aa(3,3) 
             Bd(3,10) = -ab3    
         Bd(4,2) = aa(1,1) 
         Bd(4,5) = aa(1,2) 
         Bd(4,8) = aa(1,3) 
         Bd(4,11) = -ab1    
           Bd(5,2) = aa(2,1) 
           Bd(5,5) = aa(2,2) 
           Bd(5,8) = aa(2,3) 
           Bd(5,11) = -ab2    
             Bd(6,2) = aa(3,1) 
             Bd(6,5) = aa(3,2) 
             Bd(6,8) = aa(3,3) 
             Bd(6,11) = -ab3    
         Bd(7,3) = aa(1,1) 
         Bd(7,6) = aa(1,2) 
         Bd(7,9) = aa(1,3) 
         Bd(7,12) = -ab1    
           Bd(8,3) = aa(2,1) 
           Bd(8,6) = aa(2,2) 
           Bd(8,9) = aa(2,3) 
           Bd(8,12) = -ab2    
             Bd(9,3) = aa(3,1) 
             Bd(9,6) = aa(3,2) 
             Bd(9,9) = aa(3,3) 
             Bd(9,12) = -ab3    
         do i=1,9
            sum=0.0
            do k=1,12
               sum=sum+Bd(i,k)*u(k)
            enddo
            ud(i)=sum
         enddo
         Bu(1,1)=ud(1)
         Bu(1,4)=ud(4)
         Bu(1,7)=ud(7)
           Bu(2,2)=ud(2)
           Bu(2,5)=ud(5)
           Bu(2,8)=ud(8)
             Bu(3,3)=ud(3)
             Bu(3,6)=ud(6)
             Bu(3,9)=ud(9)
         Bu(4,1)=ud(2)
         Bu(4,2)=ud(1)
         Bu(4,4)=ud(5)
         Bu(4,5)=ud(4)
         Bu(4,7)=ud(8)
         Bu(4,8)=ud(7)
           Bu(5,2)=ud(3)
           Bu(5,3)=ud(2)
           Bu(5,5)=ud(6)
           Bu(5,6)=ud(5)
           Bu(5,8)=ud(9)
           Bu(5,9)=ud(8)
             Bu(6,1)=ud(3)
             Bu(6,3)=ud(1)
             Bu(6,4)=ud(6)
             Bu(6,6)=ud(4)
             Bu(6,7)=ud(9)
             Bu(6,9)=ud(7)
         do i=1,6
            do j=1,12
               sum=0.0
               do k=1,9
                  sum=sum+Bu(i,k)*Bd(k,j)
               enddo
               BN(i,j)=sum
            enddo
         enddo
c
c         [k] = [B^T][D][B] vol
         do i=1,6
            do j=1,12
               sum=0.0
               do k=1,6
                  sum = sum + dd(i,k)*(BL(k,j)+  1.0*BN(k,j))
               enddo
               etemp(i,j)=sum
            enddo
         enddo
         do i=1,12
            do j=1,12
               sum=0.0
               do k=1,6
                  sum=sum+(BL(k,i)+  1.0*BN(k,i))*etemp(k,j)
               enddo
               ek(i,j)=sum*vol
            enddo
         enddo
*         write(iout,*)' [BL]: '
*              do i=1,6
*                 write(iout,82) (BL(i,j), j=1,12)
*              enddo
*         write(iout,*)' [BN]: '
*              do i=1,6
*                 write(iout,82) (BN(i,j), j=1,12)
*              enddo
*        write(*,*)' ek: '
*             do i=1,12
*                write(*,82) (ek(i,j), j=1,12)
*             enddo
*       stop
c
c
 82     format(1x,20(g12.6,1x))
      return
      end
c
c
c     ELeMent STiFfness for TETrahedron
      subroutine stiffG_TET( xyz,uvw,ss,ekg )
c
         implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),aa(3,3),bd(9,12),ss(6),dd(6,6),uvw(3,20)
         real*8 ee(6),u(12),ekg(60,60),etemp(9,12),sm(9,9),sm0(3,3)
c
         do i=1,12
            do j=1,12
               ekg(i,j)=0.0
            enddo
         enddo
         do i=1,6
            do j=1,12
               etemp(i,j)=0.0
            enddo
         enddo
         do i=1,9
            do j=1,12
               bd(i,j)=0.
            enddo
         enddo
c
c        construct [s]=[9x9]
         do i=1,9
            do j=1,9
               sm(i,j)=0.0
            enddo
         enddo
         sm0(1,1)=ss(1)
         sm0(2,2)=ss(2)
         sm0(3,3)=ss(3)
           sm0(1,2)=ss(4)
           sm0(2,3)=ss(5)
           sm0(1,3)=ss(6)
             sm0(2,1)=ss(4)
             sm0(3,2)=ss(5)
             sm0(3,1)=ss(6)
          do i=1,3
             do j=1,3
                sm(0+i,0+j)=sm0(i,j)
                sm(3+i,3+j)=sm0(i,j)
                sm(6+i,6+j)=sm0(i,j)
             enddo
          enddo
c
         x14 = xyz(1,1)-xyz(1,4)
         x24 = xyz(1,2)-xyz(1,4)
         x34 = xyz(1,3)-xyz(1,4)
           y14 = xyz(2,1)-xyz(2,4)
           y24 = xyz(2,2)-xyz(2,4)
           y34 = xyz(2,3)-xyz(2,4)
             z14 = xyz(3,1)-xyz(3,4)
             z24 = xyz(3,2)-xyz(3,4)
             z34 = xyz(3,3)-xyz(3,4)
*        write(*,*)' Jij'
*        write(*,*)  x14,y14,z14
*        write(*,*)  x24,y24,z24
*        write(*,*)  x34,y34,z34
         aa(1,1) = y24*z34 - y34*z24
         aa(2,1) = z24*x34 - z34*x24
         aa(3,1) = x24*y34 - x34*y24
            aa(1,2) = y34*z14 - y14*z34
            aa(2,2) = z34*x14 - z14*x34
            aa(3,2) = x34*y14 - x14*y34
               aa(1,3) = y14*z24 - y24*z14
               aa(2,3) = z14*x24 - z24*x14
               aa(3,3) = x14*y24 - x24*y14
         detj = x14*(y24*z34-y34*z24)
     &        + y14*(z24*x34-z34*x24)
     &        + z14*(x24*y34-x34*y24)
         do i=1,3
            do j=1,3
               aa(i,j) = aa(i,j)/detj
            enddo
         enddo
         vol=abs(detj)/6.0
c
*        write(*,*)' detJ: ',detj
*        write(*,*)' AA: '
*             do i=1,3
*                write(*,82) (aa(i,j), j=1,3)
*             enddo
c
c          {e}=[B]{u}
                ab1 = aa(1,1)+aa(1,2)+aa(1,3)
                ab2 = aa(2,1)+aa(2,2)+aa(2,3)
                ab3 = aa(3,1)+aa(3,2)+aa(3,3)
         bd(1,1) = aa(1,1) 
         bd(1,4) = aa(1,2) 
         bd(1,7) = aa(1,3) 
         bd(1,10) = -ab1    
           bd(2,1) = aa(2,1) 
           bd(2,4) = aa(2,2) 
           bd(2,7) = aa(2,3) 
           bd(2,10) = -ab2    
             bd(3,1) = aa(3,1) 
             bd(3,4) = aa(3,2) 
             bd(3,7) = aa(3,3) 
             bd(3,10) = -ab3    
         bd(4,2) = aa(1,1) 
         bd(4,5) = aa(1,2) 
         bd(4,8) = aa(1,3) 
         bd(4,11) = -ab1    
           bd(5,2) = aa(2,1) 
           bd(5,5) = aa(2,2) 
           bd(5,8) = aa(2,3) 
           bd(5,11) = -ab2    
             bd(6,2) = aa(3,1) 
             bd(6,5) = aa(3,2) 
             bd(6,8) = aa(3,3) 
             bd(6,11) = -ab3    
         bd(7,3) = aa(1,1) 
         bd(7,6) = aa(1,2) 
         bd(7,9) = aa(1,3) 
         bd(7,12) = -ab1    
           bd(8,3) = aa(2,1) 
           bd(8,6) = aa(2,2) 
           bd(8,9) = aa(2,3) 
           bd(8,12) = -ab2    
             bd(9,3) = aa(3,1) 
             bd(9,6) = aa(3,2) 
             bd(9,9) = aa(3,3) 
             bd(9,12) = -ab3    
*        write(*,*)' BD: '
*             do i=1,9
*                write(*,82) (bd(i,j), j=1,12)
*             enddo
c
c         [k_G] = [Bd^T][S][Bd] vol
         do i=1,9
            do j=1,12
               sum=0.0
               do k=1,9
                  sum = sum + sm(i,k)*bd(k,j)
               enddo
               etemp(i,j)=sum
            enddo
         enddo
         do i=1,12
            do j=1,12
               sum=0.0
               do k=1,9
                  sum=sum+bd(k,i)*etemp(k,j)
               enddo
               ekg(i,j)=sum*vol
            enddo
         enddo
*        write(*,*)' ekg: '
*             do i=1,12
*                write(*,82) (ekg(i,j), j=1,12)
*             enddo
*       stop
c
c
 82     format(1x,20(g12.6,1x))
      return
      end
c
c
c     ELeMent STiFfness for HEXahedron 8 nodes
c     based on iso quad from Bathe  pp295-297
      subroutine stiffE_HEX8(dd, xyz,uvw, ek)
         Implicit real*8 (a-h,o-z)
c
         real*8 xyz(3,20),uvw(3,20),ek(60,60)
         real*8 dd(6,6),BE(6,24),Bd(9,24)
c
         real*8 xg(4,4),wgt(4,4),db(6)
         real*8 dk(24,24)
         integer*4  icoord(8,3)
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
         data icoord/1,2,2,1,1,2,2,1
     &              ,1,1,2,2,1,1,2,2
     &              ,1,1,1,1,2,2,2,2 /
c
c        element stiffness
         iout=23
 20      continue
         do i=1,24
            do j=1,24
               ek(i,j)=0.0
               dk(i,j)=0.0
            enddo
         enddo
c
*          write(*,*)' [x y z] '
*          do i=1,3
*             write(*,82) (xyz(i,j), j=1,8)
*          enddo
*          write(*,*) ' '
c
         ip_n=0
         nint=2
c        nint=1
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8   
*        ip_n = ip_n + 1
*            lx = icoord(ip_n,1)
*            ly = icoord(ip_n,2)
*            lz = icoord(ip_n,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
*               write(*,84) lx,ly,lz,ip_n,(icoord(ip_n,kk),kk=1,3)
c
c                 deriv oper and jacobian
                  call BEmat8(xyz,uvw,ri,si,ti,BE,Bd,det)
*                 write(iout,*) ' BE ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (BE(i,j), j=1,24)
*                 enddo
c
c                 add contib to stiff
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
*                 write(iout,*) ' wt ',wt
                  do 70 j=1,24
                     do 40 k=1,6
                        db(k) = 0.0
                        do 40 kk=1,6
                           db(k) = db(k) + dd(k,kk)*BE(kk,j)
 40                  continue
                     do 60 i=j,24   
                        sum = 0.0
                        do kk=1,6
                           sum = sum + BE(kk,i)*db(kk)
                        enddo
                        ek(i,j) = ek(i,j) + sum*wt
                        dk(i,j) =           sum*wt
 60                  continue
 70               continue
*          write(iout,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(iout,82) (dk(i,j), j=1,24)
*          enddo
 80      continue
*          write(*,*) ' [k] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (ek(i,j), j=1,24)
*          enddo
c
c        impose symmetry
         do j=1,24
            do i=j,24
               ek(j,i)=ek(i,j)
            enddo
         enddo
c    
 82        format(1x,40(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
c
c     STIFFness Geometric for HEXahedron 8 nodes
c     based on iso quad from Bathe  pp295-297
      subroutine stiffG_HEX8( xyz,uvw, stress, ekg)
         Implicit real*8 (a-h,o-z)
c
         real*8 xyz(3,20),uvw(3,20),ekg(60,60)
         real*8 dd(6,6),BE(6,24),Bd(9,24),stress(6,27)
         real*8 sm(9,9),sm0(3,3)
c
         real*8 xg(4,4),wgt(4,4),db(9)
         real*8 dk(24,24)
         integer*4  icoord(8,3)
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
         data icoord/1,2,2,1,1,2,2,1
     &              ,1,1,2,2,1,1,2,2
     &              ,1,1,1,1,2,2,2,2 /
c
c        element stiffness
         iout=23
 20      continue
         do i=1,24
            do j=1,24
               ekg(i,j)=0.0
               dk(i,j)=0.0
            enddo
         enddo
c
*          write(*,*)' [x y z] '
*          do i=1,3
*             write(*,82) (xyz(i,j), j=1,8)
*          enddo
*          write(*,*) ' '
c
         ip_n=0
         nint=2
c        nint=1
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8   
         ip_n = ip_n + 1
*            lx = icoord(ip_n,1)
*            ly = icoord(ip_n,2)
*            lz = icoord(ip_n,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
*               write(*,84) lx,ly,lz,ip_n,(icoord(ip_n,kk),kk=1,3)
c
c                 deriv oper and jacobian
                  call BEmat8(xyz,uvw,ri,si,ti,BE,Bd,det)
                  write(iout,*) ' BE ',lx,ly,lz
                  do i=1,6
                     write(iout,82) (BE(i,j), j=1,24)
                  enddo
c
c                 construct [s]=[9x9]
                  do i=1,9
                     do j=1,9
                        sm(i,j)=0.0
                     enddo
                  enddo
                  sm0(1,1)=stress(1,ip_n)
                  sm0(2,2)=stress(2,ip_n)
                  sm0(3,3)=stress(3,ip_n)
                    sm0(1,2)=stress(4,ip_n)
                    sm0(2,3)=stress(5,ip_n)
                    sm0(1,3)=stress(6,ip_n)
                      sm0(2,1)=stress(4,ip_n)
                      sm0(3,2)=stress(5,ip_n)
                      sm0(3,1)=stress(6,ip_n)
                   do i=1,3
                      do j=1,3
                         sm(0+i,0+j)=sm0(i,j)
                         sm(3+i,3+j)=sm0(i,j)
                         sm(6+i,6+j)=sm0(i,j)
                      enddo
                   enddo
c
c                 add contib to stiff
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
                  write(iout,*) ' wt ',wt
                  do 70 j=1,24
                     do 40 k=1,9
                        db(k) = 0.0
                        do 40 kk=1,9
                           db(k) = db(k) + sm(k,kk)*Bd(kk,j)
 40                  continue
                     do 60 i=j,24   
                        sum = 0.0
                        do kk=1,9
                           sum = sum + Bd(kk,i)*db(kk)
                        enddo
                        ekg(i,j) = ekg(i,j) + sum*wt
*                       dk(i,j) =           sum*wt
 60                  continue
 70               continue
*          write(iout,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(iout,82) (dk(i,j), j=1,24)
*          enddo
 80      continue
*          write(*,*) ' [k] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (ek(i,j), j=1,24)
*          enddo
c
c        impose symmetry
         do j=1,24
            do i=j,24
               ekg(j,i)=ekg(i,j)
            enddo
         enddo
c    
 82        format(1x,40(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
c
c     ELeMent STiFfness for HEXahedron 20 nodes
c     based on iso quad from Bathe  pp295-297
      subroutine stiffE_HEX20(neltype,ired_int,dd, xyz,uvw, ek)
         Implicit real*8 (a-h,o-z)
c
         real*8 xyz(3,20),uvw(3,20),ek(60,60)
         real*8 dd(6,6),BE(6,60),Bd(9,60)
c
         real*8 xg(4,4),wgt(4,4),db(6)
*        real*8 dk(24,24)
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
c        element stiffness
         iout=23
         kmax=neltype-1
 20      continue
         do i=1,kmax*3
            do j=1,kmax*3
               ek(i,j)=0.0
*              dk(i,j)=0.0
            enddo
         enddo
c
*          write(*,*)' [x y z] '
*          do i=1,3
*             write(*,82) (xyz(i,j), j=1,8)
*          enddo
*          write(*,*) ' '
c
         ip_n=0
*        nint=3
*        nint=2
c        nint=1
*        ired_int=3
         nint=ired_int
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8   
*        ip_n = ip_n + 1
*            lx = icoord(ip_n,1)
*            ly = icoord(ip_n,2)
*            lz = icoord(ip_n,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
*               write(*,84) lx,ly,lz,ip_n,(icoord(ip_n,kk),kk=1,3)
c
c                 deriv oper and jacobian
                  call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
*                 write(iout,*) ' BE ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (BE(i,j), j=1,24)
*                 enddo
c
c                 add contib to stiff
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
*                 write(iout,*) ' wt ',wt
                  do 70 j=1,kmax*3
                     do 40 k=1,6
                        db(k) = 0.0
                        do 40 kk=1,6
                           db(k) = db(k) + dd(k,kk)*BE(kk,j)
 40                  continue
                     do 60 i=j,kmax*3   
                        sum = 0.0
                        do kk=1,6
                           sum = sum + BE(kk,i)*db(kk)
                        enddo
                        ek(i,j) = ek(i,j) + sum*wt
*                       dk(i,j) =           sum*wt
 60                  continue
 70               continue
*          write(iout,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(iout,82) (dk(i,j), j=1,24)
*          enddo
 80      continue
*          write(*,*) ' [k] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (ek(i,j), j=1,24)
*          enddo
c
c        impose symmetry
         do j=1,kmax*3
            do i=j,kmax*3
               ek(j,i)=ek(i,j)
            enddo
         enddo
c    
 82        format(1x,40(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
c
c     STIFFness Geometric for HEXahedron 8 nodes
c     based on iso quad from Bathe  pp295-297
      subroutine stiffG_HEX20(neltype,ired_int, xyz,uvw, stress, ekg)
         Implicit real*8 (a-h,o-z)
c
         real*8 xyz(3,20),uvw(3,20),ekg(60,60)
         real*8 dd(6,6),BE(6,60),Bd(9,60),stress(6,27)
         real*8 sm(9,9),sm0(3,3)
c
         real*8 xg(4,4),wgt(4,4),db(9)
*        real*8 dk(24,24)
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
c        element stiffness
         iout=23
         kmax=neltype-1
 20      continue
         do i=1,kmax*3
            do j=1,kmax*3
               ekg(i,j)=0.0
*              dk(i,j)=0.0
            enddo
         enddo
c
*          write(*,*)' [x y z] '
*          do i=1,3
*             write(*,82) (xyz(i,j), j=1,8)
*          enddo
*          write(*,*) ' '
c
         ip_n=0
*        nint=3
*        nint=2
c        nint=1
*        ired_int=3
*        ired_int=2
         nint=ired_int
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8   
         ip_n = ip_n + 1
*            lx = icoord(ip_n,1)
*            ly = icoord(ip_n,2)
*            lz = icoord(ip_n,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
*               write(*,84) lx,ly,lz,ip_n,(icoord(ip_n,kk),kk=1,3)
c
c                 deriv oper and jacobian
                  call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
*                 write(iout,*) ' BE ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (BE(i,j), j=1,kmax*3)
*                 enddo
c
c                 construct [s]=[9x9]
                  do i=1,9
                     do j=1,9
                        sm(i,j)=0.0
                     enddo
                  enddo
                  sm0(1,1)=stress(1,ip_n)
                  sm0(2,2)=stress(2,ip_n)
                  sm0(3,3)=stress(3,ip_n)
                    sm0(1,2)=stress(4,ip_n)
                    sm0(2,3)=stress(5,ip_n)
                    sm0(1,3)=stress(6,ip_n)
                      sm0(2,1)=stress(4,ip_n)
                      sm0(3,2)=stress(5,ip_n)
                      sm0(3,1)=stress(6,ip_n)
                   do i=1,3
                      do j=1,3
                         sm(0+i,0+j)=sm0(i,j)
                         sm(3+i,3+j)=sm0(i,j)
                         sm(6+i,6+j)=sm0(i,j)
                      enddo
                   enddo
c
c                 add contib to stiff
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
*                 write(iout,*) ' wt ',wt
                  do 70 j=1,kmax*3
                     do 40 k=1,9
                        db(k) = 0.0
                        do 40 kk=1,9
                           db(k) = db(k) + sm(k,kk)*Bd(kk,j)
 40                  continue
                     do 60 i=j,kmax*3   
                        sum = 0.0
                        do kk=1,9
                           sum = sum + Bd(kk,i)*db(kk)
                        enddo
                        ekg(i,j) = ekg(i,j) + sum*wt
*                       dk(i,j) =           sum*wt
 60                  continue
 70               continue
*          write(iout,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(iout,82) (dk(i,j), j=1,24)
*          enddo
 80      continue
*          write(*,*) ' [k] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (ek(i,j), j=1,24)
*          enddo
c
c        impose symmetry
         do j=1,kmax*3
            do i=j,kmax*3
               ekg(j,i)=ekg(i,j)
            enddo
         enddo
*          do i=1,kmax*3
*             write(iout,82) (ekg(i,j), j=1,kmax*3)
*          enddo
c    
 82        format(1x,90(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccc
c
c     FORM STIFfness matrix  [K]  
      subroutine formstif_3D( stf,npijkm,idbc,maxnode,maxelem,
     &                       xord0, yord0,zord0,
*    &                                     xord, yord,zord,
     &                       dispfult,maxdim,
     &                       nelt,iglobal,
     &                       qms,
     &                       prop,nmprop,ippp,iprofv,nloc,kref)
c
          implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(nel)
         integer nmprop(nel)
         real*8  prop(ippp,10),dd(6,6)
         integer iprofv(neq), nloc(neq)
c
         real*8 xord0(nnp), yord0(nnp),zord0(nnp)
*        real*8 xord(nnp), yord(nnp),zord(nnp)
         real*8 stf(maxstiff)
         real*8 ek(60,60)
         real*8 qms(neq), temp,ss(6)
c
         real*8 xyz(3,20),uvw(3,20)    
c
*        real*8  dispfult(900,3)
         real*8  dispfult(maxdim,3)
c
         write(ilog,*)'@@ << in formstiff_3D >>'
         iecho=0
c
c        initialize [K]  to zero
         do i=1,maxstiff
            stf(i)=0.0
         enddo
c
c        form each element matrix, and assemble
         do 50 n=1,nel
            mat=nmprop(n)
            neltype = npijkm(n,1)
c
                e0=prop(mat,1)
                g0=prop(mat,2)
                t0=prop(mat,3)
                r0=prop(mat,4)
               pl0=prop(mat,5)
                call dmat(e0,g0,dd)
c
              if (neltype .eq. 5) then
c                 tetra
                  do k=1,4 
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
                     uvw(1,k)=dispfult(node,1)
                     uvw(2,k)=dispfult(node,2)
                     uvw(3,k)=dispfult(node,3)
                  enddo
                  call stiffE_TET( xyz,uvw,dd,ss,ek )
*                 write(iout,*)'STIFF:'
*                 do ii=10,12
*                    write(iout,86) (ek(ii,jj),jj=10,12)
*                 enddo
c
              elseif (neltype .eq. 9) then
c                 Hex8
                  do k=1,8 
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
                     uvw(1,k)=dispfult(node,1)
                     uvw(2,k)=dispfult(node,2)
                     uvw(3,k)=dispfult(node,3)
                  enddo
                  call stiffE_HEX8(dd, xyz,uvw,ek )
*                 write(iout,*)'STIFF:'
*                 do ii=10,12
*                    write(iout,86) (ek(ii,jj),jj=10,12)
*                 enddo
              elseif (neltype .eq. 21) then
c                 Hex20
                  kmax=neltype-1
                  do k=1,kmax 
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
                     uvw(1,k)=dispfult(node,1)
                     uvw(2,k)=dispfult(node,2)
                     uvw(3,k)=dispfult(node,3)
                  enddo
                  ired_int=3
                  ired_int=ri3
                  call stiffE_HEX20(neltype,ired_int,dd, xyz,uvw,ek )
*                 write(iout,*)'STIFF:'
*                 do ii=1,60
*                    write(iout,82) (ek(ii,jj),jj=1,60)
*                 enddo
*                 stop
c
               endif
c
 82            format(1x,70(g12.6,1x))
 59            format(1x,18(g12.6))
 86            format(1x, 6(g12.6,1x))
               ielm=n
               call assembCOL_3D(stf,ek,idbc,maxnode,npijkm,maxelem,
     &                           ielm,nloc)
*                 do ii=1,maxstiff
*                    write(iout,*) ii,stf(ii)
*                 enddo
*                 stop
c
 50      continue 
c        END loop over all elements
*        write(ilog,*) '@@ FORMSTFF:   Formed  [K]  OK'
      return
      end
c
c
c     FORM GEOMetric stiffness matrices  [KGx],[KGy],[KGxy]  
      subroutine formgeom_3D(geo,npijkm,idbc,maxnode,maxelem,
     &                     xord0, yord0,zord0,
*    &                                         xord, yord,zord,
     &                       dispfult,maxdim,
     &                    nelt,
     &                    iglobal,
     &                    qms,
     &                    prop,nmprop,ippp,iprofv,nloc)
c
          implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npijkm(maxelem,21),idbc(maxnode,3)
         integer nelt(nel), iprofv(neq), nloc(neq)
         integer nmprop(nel)
         real*8  prop(ippp,10)
c
         real*8  xord0(nnp), yord0(nnp),zord0(nnp)
*        real*8  xord(nnp), yord(nnp),zord(nnp)
         real*8 ekg(60,60)
         real*8 geo(maxstiff),qms(neq)
         real*8 ss(6),dd(6,6)
         real*8 xyz(3,20),uvw(3,20),stress(6,27)
c
*        real*8  dispfult(900,3)
         real*8  dispfult(maxdim,3)
c
         write(ilog,*)'@@  << in formgeom_3D >>'
         iecho=0
c
c        initialize [G] to zero
         do i=1,maxstiff
            geo(i)=0.0
         enddo
c
c        form each element matrix, and assemble
         rewind(igeo)
         do 50 n=1,nel
                 mat=nmprop(1)
                 e0=prop(mat,1)
                 g0=prop(mat,2)
                 t0=prop(mat,3)
            neltype=npijkm(n,1)
                call dmat(e0,g0,dd)
c
              if (neltype .eq. 5) then
c                 tetra
                  do k=1,4 
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
                     uvw(1,k)=dispfult(node,1)
                     uvw(2,k)=dispfult(node,2)
                     uvw(3,k)=dispfult(node,3)
                  enddo
                  call stress_TET( xyz,uvw,dd,ss)
*                     write(iout,*) '@@ stress G ',n
*                     write(iout,82) (ss(jj), jj=1,6)
                  call stiffG_TET( xyz,uvw,ss,ekg )
*                 write(iout,*)'STIFF G:'
*                 do ii=10,12
*                    write(iout,82) (ekg(ii,jj),jj=10,12)
*                 enddo
c
              elseif (neltype .eq. 9) then
c                 Hex8 
                  kmax=neltype-1
                  do k=1,kmax 
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
                     uvw(1,k)=dispfult(node,1)
                     uvw(2,k)=dispfult(node,2)
                     uvw(3,k)=dispfult(node,3)
                  enddo
                  call stress_HEX8( xyz,uvw,dd,stress)
*                  write(iout,*) '@@ stress G ',n
*                  do ii=1,8
*                     write(iout,82) (stress(jj,ii), jj=1,6)
*                  enddo
                  call stiffG_HEX8( xyz,uvw,stress,ekg )
*                 write(iout,*)'STIFF G:'
*                 do ii=10,12
*                    write(iout,82) (ekg(ii,jj),jj=10,12)
*                 enddo
c
              elseif (neltype .eq. 21) then
c                 Hex8 
                  kmax=neltype-1
                  do k=1,kmax 
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
                     uvw(1,k)=dispfult(node,1)
                     uvw(2,k)=dispfult(node,2)
                     uvw(3,k)=dispfult(node,3)
                  enddo
                  ired_int=2
                  ired_int=3
                  ired_int=ri3
                  inonlinear=1
                  call stress_HEX20(neltype,ired_int, xyz,uvw,dd,stress,
     &                              inonlinear)
*                 write(iout,*) '@@ stress G ',n
*                 do ii=1,27
*                    write(iout,82) (stress(jj,ii), jj=1,6)
*                 enddo
                  call stiffG_HEX20(neltype,ired_int,xyz,uvw,stress,ekg)
*                 write(iout,*)'STIFF G:'
c
              endif
               ielm=n
               call assembCOL_3D(geo,ekg,idbc,maxnode,npijkm,maxelem,
     &                           ielm,nloc)
 50       continue
*       stop
 82       format(1x,60(g12.6,1x))
c
*        write(ilog,*) '@@ GEOMSTFF:   Formed  [G]  OK'
      return
      end
c
c
c     Update geometry
      subroutine update_geom_3D( disp,dispful,maxdim, idbc,maxnode,
     &                   xord0,yord0,zord0,
     &                   xord,yord,zord
     &                  )
c
          implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer idbc(maxnode,3)
         real*8  disp(neq  )
         real*8  xord0(nnp), yord0(nnp), zord0(nnp)
         real*8  xord(nnp), yord(nnp), zord(nnp)
c
*        real*8  dispful(900,3)
         real*8  dispful(maxdim,3)
c
c          reassign displacements to dispful( )
           do 20 n=1,nnp
              do k=1,3
                 ieqn1=idbc(n,k)
                 if (ieqn1 .eq. 0) then
                     dispful(n,k)=0.0
                 else
                     dispful(n,k)=disp(ieqn1)
                 endif
              enddo
 20        continue
c
           do n=1,nnp
              xord(n)=xord0(n)+dispful(n,1)
              yord(n)=yord0(n)+dispful(n,2)
              zord(n)=zord0(n)+dispful(n,3)
           enddo
c      
       return
       end
c
c
c     ASSEMBle consistent FORce 
      subroutine assembFOR_HEX(gf,ff, idbc,maxnode,npijkm,maxelem,ielm)
          implicit real*8 (a-h,o-z)
                include 'commons.std'
         integer  idbc(maxnode,3), npijkm(maxelem,21)
         real*8   ff(3,20), gf(neq  )
         integer  all_eqn(60)
c
c        Set idof to posible DoF number of each nodes
            neltype = npijkm(ielm,1)
            kmax  = neltype-1
              do kk=1,kmax
                 nodek = npijkm(ielm,1+kk)
                 all_eqn( (kk-1)*3 + 1) = idbc(nodek,1)
                 all_eqn( (kk-1)*3 + 2) = idbc(nodek,2)
                 all_eqn( (kk-1)*3 + 3) = idbc(nodek,3)
              enddo
c
c        Store the values for individual array in global array
         nn=0
         do 20 j= 1,kmax
            do i= 1,3
               nn=nn+1
               ieqn1 = all_eqn(nn)
               if (ieqn1 .gt. 0) then
!$OMP ATOMIC
                   gf(ieqn1) = gf(ieqn1) + ff(i,j)
               endif
            enddo
 20      continue
c
      return
      end
c
c
c     ASSEMBle consistent FORce 
      subroutine assembFOR_3D( aa, a, idbc,maxnode, npijkm,maxelem,ielm)
          implicit real*8 (a-h,o-z)
                include 'commons.std'
         integer  idbc(maxnode,3), npijkm(maxelem,21)
         real*8   a(60), aa(neq  )
         integer  all_eqn(60)
c
c        Set idof to posible DoF number of each nodes
            neltype = npijkm(ielm,1)
            kmax  = neltype-1
              do kk=1,kmax
                 nodek = npijkm(ielm,1+kk)
                 all_eqn( (kk-1)*3 + 1) = idbc(nodek,1)
                 all_eqn( (kk-1)*3 + 2) = idbc(nodek,2)
                 all_eqn( (kk-1)*3 + 3) = idbc(nodek,3)
              enddo
c
c        Store the values for individual array in global array
         do 20 i= 1, kmax*3
            ieqn1 = all_eqn(i)
            if (ieqn1 .gt. 0) then
                aa(ieqn1) = aa(ieqn1) + a(i)
            endif
 20      continue
c
      return
      end
c
c
c     ASSEMBle element matrices in COLumn form for HEX 20 nodes
      subroutine assembCOL_3D( aa, a, idbc,maxnode, 
     &                          npijkm,maxelem,ielm,nloc)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
*        integer  idof_node(60),idbc(maxnode,3),npijkm(maxelem,21)
         integer  idbc(maxnode,3),npijkm(maxelem,21)
         integer  nloc(neq)
         integer  all_eqn(60)
c
         real*8   a(60,60)
         real*8   aa(maxstiff)
c
c        Set idof to posible DoF number of each nodes
            neltype = npijkm(ielm,1)
            kmax  = neltype-1
              do kk=1,kmax
                 nodek = npijkm(ielm,1+kk)
                 all_eqn( (kk-1)*3 + 1) = idbc(nodek,1)
                 all_eqn( (kk-1)*3 + 2) = idbc(nodek,2)
                 all_eqn( (kk-1)*3 + 3) = idbc(nodek,3)
              enddo
*           write(iout,81) ielm,(all_eqn(k), k=1,kmax*3)
 81         format(1x,61(i5,1x))
c
c        Store the values for individual array in global array
         do 20 i= 1, kmax*3
            ieqn1 = all_eqn(i)
            if (ieqn1 .gt. 0) then
               do 30 j= i, kmax*3
                  ieqn2 = all_eqn(j)
                  if (ieqn2 .gt. 0) then
                     if (ieqn1 .gt. ieqn2) then
                        jband= (ieqn1-ieqn2)+1
                        iloc = nloc(ieqn1) 
                        aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
c                       aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
                     else
                        jband= (ieqn2-ieqn1)+1
                        iloc = nloc(ieqn2) 
                        aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
c                       aa(ieqn1,jband) = aa(ieqn1,jband) + a(i,j)
                     endif
                  endif
 30            continue
           endif
 20      continue
c
      return
      end
c
c
c     NonStad's POST ANalysis
      subroutine non_postan( disp,dispful,
     &                       idbc,npijkm,maxnode,maxelem, 
     &                       xord0,yord0,zord0,
     &                       strn,strs,forc,
     &                       e_elm,s_elm,
     &                       e_elm_ip,s_elm_ip,sc_elm_ip,wk_elm_n,
     &                       prop,nmprop,ippp,iglobal )
c
          implicit real*8 (a-h,o-z)
             include 'commons.std'
         parameter (ireadmax=101)
c
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3 )
         integer mnod(ireadmax,2)
         integer iscan(2000,3)
c
         real*8  disp(neq  ),  dispful(nnp,3)
         real*8  strs(nnp,6),strn(nnp,6),forc(nel,60)
         real*8  e_elm(nel,6),s_elm(nel,6)
         real*8  e_elm_ip(nel,27,6),s_elm_ip(nel,27,6)
         real*8  sc_elm_ip(nel,27,6)
         real*8  wk_elm_n(nel,20,6)
         real*8  xord0(nnp), yord0(nnp), zord0(nnp)
         real*8  xyz(3,20),uvw(3,20)
         real*8  ssi(27,6),ssn(20,6),se_av(6)
         real*8  strain_ip(6,27),stress_ip(6,27), force_n(3,20)
         real*8                  stressc_ip(6,27)
         real*8  prop(ippp,10),dd(6,6)
         integer nmprop(nel)
c
         rewind(idyn)
c
             write(*,*)'CHOOSE snapshot:'
             write(*,*)'       0=return   '
             write(*,*)'       #=which snap '
             write(*,*)'      -#=all   snaps '
             call zzwrt('  -->  ')
             read(ikbd,*) iwhich 
             write(ilog,*) iwhich , '   ::which snap '
             if (iwhich .lt. 0) goto 2000
c
c            ONE snap
             if (iwhich .gt. 0) then
                 rewind (isnp)
                 read(isnp) time,neqm
*                write(ilog,*)'@@max snaps: ',isnapmax,neqm
*                if (iwhich .gt. isnapmax) iwhich=isnapmax
c
c                stop at this snap and write displacements
                 rewind(isnp)
                 do iw=1,iwhich
                    read(isnp,end=110) time,neqm
                    do i=1,neqm
                       read(isnp,end=110) disp(i)
                    enddo
                    call zzwrti(iw) 
                 enddo
                 iw=iw+1
 110             continue
                 iw=iw-1
                 write(*,'(a)')' @@ '
                 gnum1=time
                 write(ilog,*)'@@time     : ',gnum1,iw
                 rewind(idis)
                 do i=1,neqm
                    write(idis) disp(i)
                 enddo
c
c
                 iloop=0
 1               continue
                 write(*,*)'  POST menu: '
                 write(*,*)'           0=return '
*                write(*,*)'               LOCAL (member coords)  '
*                write(*,*)'           2=nodal strain [ue] '
*                write(*,*)'           3=nodal stress  '
*                write(*,*)'           4=nodal force '
*                write(*,*)'           5=element stress '
*                write(*,*)'           6=element (average) stress '
                 write(*,*)'               GLOBAL '
                 write(*,*)'          11=Global displacements'
                 write(*,*)'          12=Global loads '
*                write(*,*)'               CONTOURs '
*                write(*,*)'          31=store displacement data '
*                write(*,*)'          32=store strain       data '
*                write(*,*)'          33=store stress       data '
*                write(*,*)'               <<deformed>>          '
*                write(*,*)'          41=store displacement data '
*                write(*,*)'          42=store strain       data '
*                write(*,*)'          43=store stress       data '
                 write(*,*)'               <<SHAPES>>            '
                 write(*,*)'         141=store displacement data '
                 write(*,*)'               <<Hex20>>          '
                 write(*,*)'         41x=documented displacement '
                 write(*,*)'             411=displacement     '
                 write(*,*)'         42x=documented strain '
                 write(*,*)'             421=element strain IP '
                 write(*,*)'             422=element strain N  '
                 write(*,*)'             424=nodal   strain average '
                 write(*,*)'             427=element strain special '
                 write(*,*)'             428=nodal   strain special '
                 write(*,*)'         43x=documented Kirchhoff stress '
                 write(*,*)'             431=element stress IP '
                 write(*,*)'             432=element stress N  '
                 write(*,*)'             434=nodal   stress average '
                 write(*,*)'             437=element stress special '
                 write(*,*)'             438=nodal   stress special '
                 write(*,*)'         44x=documented Cauchy stress '
                 write(*,*)'             441=element stress IP '
                 write(*,*)'             442=element stress N  '
                 write(*,*)'             444=nodal   stress average '
                 write(*,*)'             447=element stress special '
                 write(*,*)'             448=nodal   stress special '
                 write(*,*)'         450=force '
                 call zzwrt(' SELECT--> ')
                 read(ikbd,*,err=1)  ioutput
                 write(ilog,*) ioutput, ' ::1=disp'
                 if (ioutput .eq. 0) then
                       return
                 endif
c
 153             continue
                 rewind (idis)
                 do i=1,neq
                    read(idis) disp(i)
                 enddo
c                reassign displacements to dispful( )
                 do i=1,nnp
                    do j=1,3
                       ieqnum=idbc(i,j)
                       if (ieqnum .eq. 0) then
                           dispful(i,j) = 0.0 
                       else
                           dispful(i,j) = disp(ieqnum)
                       endif
                    enddo    
                 enddo    
*                  do i=1,nnp
*                     xordt(i)=xord(i) + dispful(i,1)
*                     yordt(i)=yord(i) + dispful(i,2)
*                     zordt(i)=zord(i) + dispful(i,3)
*                  enddo
c
                 if (ioutput .ge. 141 .AND. ioutput .le. 149) then
                     goto 100
                 elseif (ioutput .ge. 410 .AND. ioutput .le. 419) then
                     goto 100
                 endif
c
c                no need to loop through elements multiple times
                 if (iloop   .eq. 1) goto 100
                     iloop=1
c
c                each element, calculate strain, stress at centroid
                 write(*,*)'@@ loop over elements ',nel
                 write(*,*)' '
                 nel51=nel/50+1
                 call zzwrti(nel51)
                 call zzwrt ('>')
                 do 50 i=1,nel
                    if (mod(i,50) .eq. 0) call zzwrt('.')
                    neltype=npijkm(i,1)
c
                    if (neltype .eq. 3) then
c                       frame
                    elseif (neltype .eq. 4) then
c                       plate
                    elseif (neltype .eq. 5) then
c                       3-D solid TET
                    elseif (neltype .eq. 9) then
c                       3-D solid HEX8
c
                    elseif (neltype .eq. 21) then
c                       3-D solid HEX20
                        mat=nmprop(i)
                        e0=prop(mat,1)
                        g0=prop(mat,2)
                        r0=prop(mat,4)
                        call dmat(e0,g0,dd)
c
                        kmax=neltype-1
                        do k=1,kmax
                           node=npijkm(i,1+k)
                           xyz(1,k)=xord0(node)
                           xyz(2,k)=yord0(node)
                           xyz(3,k)=zord0(node)
c
                           uvw(1,k)=dispful(node,1)
                           uvw(2,k)=dispful(node,2)
                           uvw(3,k)=dispful(node,3)
                        enddo
c
ctag2
*                       write(iout,*) '@@ Elem: ',i
                        ired_int=3
                        call convert_HEX20_non (neltype,xyz,uvw,dd,
     &                              ired_int,
     &                        strain_ip, stress_ip,stressc_ip,force_n )
*                       call convert_HEX20_non2(neltype,xyz,uvw,dd,
*    &                              ired_int)
*    &                              strain_ip, stress_ip,force_n )
*                       write(iout,*) '@@ All: ',i
*                       do k=1,27
*                          write(iout,81) (stressc_ip(j,k), j=1,6)
*                       enddo
*                       write(iout,81) (stressc_ip(j,1), j=1,6)
c
c                       SAVE in array of element values
                        ielm=i
                        max_int=ired_int**3
                        do k=1,max_int   
                           do j=1,6   
                              e_elm_ip(ielm,k,j)=strain_ip(j,k)
                              s_elm_ip(ielm,k,j)=stress_ip(j,k)
                             sc_elm_ip(ielm,k,j)=stressc_ip(j,k)
                           enddo
                        enddo
                        nn=0
                        do k=1,kmax   
                           do j=1,3   
                              nn=nn+1
                              forc(ielm,nn)=force_n(j,k)
                           enddo
                        enddo
c                       element averages
                        do j=1,6   
                           sume=0.0
                           sums=0.0
                           do k=1,max_int   
                              sume=sume+strain_ip(j,k)
                              sums=sums+stress_ip(j,k)
                           enddo
                           e_elm(i,j)=(sume/max_int)*1.0e6
                           s_elm(i,j)=sums/max_int
                        enddo
c
                    endif
c                   bottom distinction between elements
 50              continue
*                       write(iout,*) '@@ All 2: ',ielm
*                       do k=1,27
*                          write(iout,81) (sc_elm_ip(ielm,k,j), j=1,6)
*                       enddo
c
c                STORE results in <<OUT>> file
 100             continue
                        ipmax=ired_int**3
                 write(ilog,'(a)')' @@'
*                write(iout,'(a)')' @@'
                 write(*   ,'(a)')' writing output'
                 if (ioutput .eq. 1) then
c
                 elseif (ioutput .eq. 141 .OR. ioutput .eq. 142
     &                                      .OR. ioutput .eq. 143) then
                     rewind(iout)
                     evout=gnum1 
                     call shapes_HEX_D(npijkm,maxelem,
     &                         xord0,yord0,zord0,dispful,evout,ioutput)
c
                 elseif (ioutput .eq. 411) then
                     call outdis (dispful,'Single snap')
c
c                STRAIN
                 elseif (ioutput .eq. 421) then
                     write(iout,*)'Element STRAIN int pt: '
                     write(iout,'(3a)')'  Elem    IP: ',
     &                     '  Exx         Eyy          Ezz',
     &                     '        2Exy          2Eyz         2Exz'
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        j=1
                           write(iout,83) n,j
     &                               ,(e_elm_ip(n,1,k), k=1,6)
                        do j=2,ipmax
                           write(iout,84) j
     &                               ,(e_elm_ip(n,j,k), k=1,6)
                        enddo 
                     enddo 
c
                 elseif (ioutput .eq. 422) then
                     write(iout,*)'Element STRAIN nodal: '
                     write(iout,'(3a)')'  Elem     N: ',
     &                     '  Exx         Eyy          Ezz',
     &                     '        2Exy          2Eyz         2Exz'
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=e_elm_ip(n,j,k)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        write(iout,83) n,npijkm(n,1+1),(ssn(1,k), k=1,6)
                        do j=2,kmax
                          write(iout,84) npijkm(n,1+j),(ssn(j,k), k=1,6)
                        enddo 
                     enddo 
ctag1
                 elseif (ioutput .eq. 424) then
c                    get element nodal values
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        do i=1,ipmax
                           do j=1,6
                              ssi(i,j)=e_elm_ip(n,i,j)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        do i=1,kmax
                           do j=1,6
                              wk_elm_n(n,i,j) = ssn(i,j)
                           enddo 
                        enddo 
                     enddo 
                     call node_avg_HEX20(npijkm,maxelem
     &                                  ,wk_elm_n,strn)
                     write(iout,*)'NODAL STRAIN averages: '
                     write(iout,*)' node:   ',
     &                               '  Exx        Eyy         Ezz  ',
     &                     '           2Exy       2Eyz         2Exz'
                     do i=1,nnp
                        write(iout,23) i,(strn(i,j), j=1,6)
                     enddo 
                     call out_max_TET(strn,nnp)
c
                 elseif (ioutput .eq. 427) then
                     write(iout,*)'element STRAIN specials: '
                 write(iout,*)' Elem:    E1              E2         E3',
     &                     '                  von Mises  '
                     call out_prin_TET(e_elm,nel)
c
                 elseif (ioutput .eq. 428) then
                     write(iout,*)'nodal STRAIN average specials: '
                  write(iout,*)' Node:    E1              E2        E3',
     &                     '                  von Mises '
                     call out_prin_TET(strn,nnp)
c
c                Kirchhoff STRESS
                 elseif (ioutput .eq. 431) then
                     write(iout,*)'Element Kirchhoff STRESS int pt: '
                     write(iout,'(3a)')'  Elem    IP: ',
     &                     '  Sxx         Syy          Szz',
     &                     '         Sxy           Syz          Sxz'
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        j=1
                           write(iout,83) n,j
     &                               ,(s_elm_ip(n,1,k), k=1,6)
                        do j=2,ipmax
                           write(iout,84) j
     &                               ,(s_elm_ip(n,j,k), k=1,6)
                        enddo 
                     enddo 
c
                 elseif (ioutput .eq. 432) then
                     write(iout,*)'Element Kirchhoff STRESS nodal: '
                     write(iout,'(3a)')'  Elem     N: ',
     &                     '  Sxx         Syy          Szz',
     &                     '         Sxy           Syz          Sxz'
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=s_elm_ip(n,j,k)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        write(iout,83) n,npijkm(n,1+1),(ssn(1,k), k=1,6)
                        do j=2,kmax
                          write(iout,84) npijkm(n,1+j),(ssn(j,k), k=1,6)
                        enddo 
                     enddo 
c
                 elseif (ioutput .eq. 434) then
c                    get element nodal values
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        do i=1,ipmax
                           do j=1,6
                              ssi(i,j)=s_elm_ip(n,i,j)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        do i=1,kmax
                           do j=1,6
                              wk_elm_n(n,i,j) = ssn(i,j)
                           enddo 
                        enddo 
                     enddo 
                     call node_avg_HEX20(npijkm,maxelem
     &                                  ,wk_elm_n,strs)
                     write(iout,*)'NODAL Kirchhoff STRESS averages: '
                     write(iout,*)' node:   ',
     &                               '  Sxx        Syy         Szz  ',
     &                     '            Sxy        Syz          Sxz'
                     do i=1,nnp
                        write(iout,23) i,(strs(i,j), j=1,6)
                     enddo 
                     call out_max_TET(strs,nnp)
c
                 elseif (ioutput .eq. 437) then
                     write(iout,*)'element STRESS specials: '
                 write(iout,*)' Elem:    S1              S2         S3',
     &                     '                  von Mises  '
                     call out_prin_TET(s_elm,nel)
c
                 elseif (ioutput .eq. 438) then
                     write(iout,*)'nodal STRESS average specials: '
                  write(iout,*)' Node:    S1              S2        S3',
     &                     '                  von Mises '
                     call out_prin_TET(strs,nnp)
c
c
c                Cauchy STRESS
                 elseif (ioutput .eq. 441) then
                     write(iout,*)'Element Cauchy STRESS int pt: '
                     write(iout,'(3a)')'  Elem    IP: ',
     &                     '  Sxx         Syy          Szz',
     &                     '         Sxy           Syz          Sxz'
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        j=1
                           write(iout,83) n,j
     &                               ,(sc_elm_ip(n,1,k), k=1,6)
                        do j=2,ipmax
                           write(iout,84) j
     &                               ,(sc_elm_ip(n,j,k), k=1,6)
                        enddo 
                     enddo 
ctag2
                 elseif (ioutput .eq. 442) then
                     write(iout,*)'Element Cauchy STRESS nodal: '
                     write(iout,'(3a)')'  Elem     N: ',
     &                     '  Sxx         Syy          Szz',
     &                     '         Sxy           Syz          Sxz'
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=sc_elm_ip(n,j,k)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        write(iout,83) n,npijkm(n,1+1),(ssn(1,k), k=1,6)
                        do j=2,kmax
                          write(iout,84) npijkm(n,1+j),(ssn(j,k), k=1,6)
                        enddo 
                     enddo 
c
                 elseif (ioutput .eq. 444) then
c                    get element nodal values
                     do n=1,nel
                        neltype = npijkm(n,1)
                        kmax=neltype-1
                        do i=1,ipmax
                           do j=1,6
                              ssi(i,j)=sc_elm_ip(n,i,j)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        do i=1,kmax
                           do j=1,6
                              wk_elm_n(n,i,j) = ssn(i,j)
                           enddo 
                        enddo 
                     enddo 
                     call node_avg_HEX20(npijkm,maxelem
     &                                  ,wk_elm_n,strs)
                     write(iout,*)'NODAL Cauchy STRESS averages: '
                     write(iout,*)' node:   ',
     &                               '  Sxx        Syy         Szz  ',
     &                     '            Sxy        Syz          Sxz'
                     do i=1,nnp
                        write(iout,23) i,(strs(i,j), j=1,6)
                     enddo 
                     call out_max_TET(strs,nnp)
c
                 elseif (ioutput .eq. 447) then
                     write(iout,*)'element STRESS specials: '
                 write(iout,*)' Elem:    S1              S2         S3',
     &                     '                  von Mises  '
                     call out_prin_TET(s_elm,nel)
c
                 elseif (ioutput .eq. 448) then
                     write(iout,*)'nodal STRESS average specials: '
                  write(iout,*)' Node:    S1              S2        S3',
     &                     '                  von Mises '
                     call out_prin_TET(strs,nnp)
c
                 endif
                 goto 1
             endif
c            bottom single snap
c
c
c            ALL snaps
 2000        continue
             if (iwhich .lt. 0) then
c                we look at all snaps
                 write(*,*)'CHOOSE output: '
                 write(*,*)'           0=return  '
                 write(*,*)'          51=element stress'
                 write(*,*)'                <movie>     '
                 write(*,*)'         111=node displacement'
                 write(*,*)'         141=store displacement data '
                 write(*,*)'               <<Hex20>>          '
                 write(*,*)'         41x=displacement history '
                 write(*,*)'             411=displacement     '
                 write(*,*)'         42x=Lagrange strain  history'
                 write(*,*)'             421=element strain IP '
                 write(*,*)'             422=element strain N  '
                 write(*,*)'             424=nodal   strain average '
                 write(*,*)'             427=element strain special '
                 write(*,*)'             428=nodal   strain special '
                 write(*,*)'         43x=Kirchhoff stress  history'
                 write(*,*)'             431=element stress IP '
                 write(*,*)'             432=element stress N  '
                 write(*,*)'             434=nodal   stress average '
                 write(*,*)'             437=element stress special '
                 write(*,*)'             438=nodal   stress special '
                 write(*,*)'         44x=Cauchy stress  history'
                 write(*,*)'             441=element stress IP '
                 write(*,*)'             442=element stress N  '
                 write(*,*)'             444=nodal   stress average '
                 write(*,*)'             447=element stress special '
                 write(*,*)'             448=nodal   stress special '
                 write(*,*)'         440=force  history'
                 call zzwrt('  -->  ')
                 read(ikbd,*) ioutput
                 write(ilog,*) ioutput,' ::ioutput(411=disp 42=e 43=s)'
c
                 if (ioutput .eq. 0) then
                     return
                 elseif (ioutput .eq. 51) then
                     write(*,*)'INPUT: elem # '
                     call zzwrt('  -->  ')
                     read(ikbd,*) iwhich 
                     write(ilog,*) iwhich , ' ::iwhich node '
*                           mout=1
*                           mnod(mout)=iwhich
                 elseif (ioutput .eq. 111) then
                     goto 390
c                elseif (ioutput .ne. 4) then
c                    write(*,*)'INPUT: node # '
c                    call zzwrt('  -->  ')
c                    read(ikbd,*) iwhich 
c                    write(ilog,*) iwhich , ' ::iwhich node '
c                    now find the elements connected to this node
c                    mout=0
c                    node=iwhich 
c                    do n=1,nel
*                       ni=npi(n)
*                       nj=npj(n)
*                       nk=npk(n)
*                       if (node .eq. ni .or. node .eq. nj 
*    &                                   .or. node .eq. nk) then
*                           mout=mout+1
*                           mnod(mout)=n
*                       endif
c                    enddo
c                    write(ilog,'(a)')'@@  '
c                    write(ilog,'(a,i5)')'@@ # elems: ',mout
c                    if (mout .eq. 0) then
c                        write(ilog,'(a)')'@@ no surrounding elems !!!'
c                        write(*   ,'(a)')'@@ no surrounding elems !!!'
c                        goto 999
c                    endif
c                    mlines=1+ mout/8
c                    do i=1,mlines
c                       write(ilog,29)'@@ elems: ',(mnod(j), j=1,8)
c                    enddo
c
                 elseif (ioutput .eq. 411) then
c                    nodal disp
                     write(*,*)'INPUT: node # '
                     call zzwrt('  -->  ')
                     read(ikbd,*) inode 
                     write(ilog,*) inode , ' ::node '
                 elseif (ioutput .eq. 421 .OR. ioutput .eq. 431
     &              .OR. ioutput .eq. 441                      ) then
c                    element IP value 
                     write(*,*)'INPUT: IP #  |  elem #   '
                     call zzwrt('  -->  ')
                     read(ikbd,*)  ip_n,ielem 
                     write(ilog,*)  ip_n,ielem , ' ::ip_n,elem '
                            mout=1
                            mnod(mout,1)=ielem
                 elseif (ioutput .eq. 422 .OR. ioutput .eq. 432
     &              .OR. ioutput .eq. 442                      ) then
c                    element Nodal value 
                     write(*,*)'INPUT: Node #  |  elem #   '
                     call zzwrt('  -->  ')
                     read(ikbd,*) inode,ielem
                     write(ilog,*) inode,ielem,  ' ::node,elem '
                            mout=1
                            mnod(mout,1)=ielem
                            do k=2,21
                               if (npijkm(ielem,k) .eq. inode) then
                                   klocn = k
                                   mnod(mout,2)=klocn
                                   write(ilog,*) '@@ klocn: ',klocn
                                   goto 211
                               endif
                            enddo
 211                        continue
c
                 elseif (ioutput .eq. 424 .OR. ioutput .eq. 434
     &              .OR. ioutput .eq. 444                      ) then
c                    Nodal average
                     write(*,*)'INPUT:  Node # '
                     call zzwrt('  -->  ')
                     read(ikbd,*) inode 
                     write(ilog,*) inode, ' ::node'
c                    find all connected elements
                     mout=0
                     do i=1,nel
                        do k=2,21
                           if (npijkm(i,k) .eq. inode) then
                               mout = mout + 1
                               mnod(mout,1)=i
                               mnod(mout,2)=k
                                   write(ilog,*) '@@conn elms:',i,k
                                   goto 212
                            endif
                        enddo
 212                    continue
                     enddo
                 endif
c
c                READ all snaps
 390             continue
                 write(*,*)'@@ READing all snaps'
                 write(ilog,*)'@@ READing all snaps'
                 rewind (isnp)
                 read(isnp,err=901) psload,neqm
 901             continue
                 write(ilog,*)'@@ psload neqm : ', psload,neqm
                 rewind (isnp)
                 ikont=0
                 rewind (iout)
                    write(*,'(a)')' @@ '
c
 300             continue
                 ikont=ikont+1
                 read(isnp,err=902,end=999) psload,neqm
 902             continue
                 write(ilog,*)'@@ psload  : ', psload,ikont
                 time=psload
                 if (mod(ikont,10) .eq. 0) call zzwrti(ikont)
c
                 ddz=0
                 do i=1,neqm
                    read(isnp,end=999) disp(i)
                    ddz=ddz+abs(disp(i))
                 enddo
                 gnum1=psload
                 write(ilog,*)'@@ tot disp  : ',ddz, ddz/neqm
c
c                reassign displacements to dispful( )
                   do i=1,nnp
                      do j=1,3
                         ieqnum=idbc(i,j)
                         if (ieqnum .eq. 0) then
                             dispful(i,j) = 0.0 
                         else
                             dispful(i,j) = disp(ieqnum)
                         endif
                      enddo    
                   enddo    
c
 310             continue
                 if (ioutput .le. 1 .OR. ioutput .eq. 411) then
c                    no need for strains etc
                     goto 359
                 endif
                 if (ioutput .eq. 111) then
c                    no need for strains etc
                     goto 359
                 endif
c
c                   zero matrices at this time
                    do j=1,6
                       se_av(j)=0.0
                    enddo
c
ctag3
c                For each element, get strain, stress etc
                 do 350 m=1,mout
                    i = mnod(m,1)
                    klocn=mnod(m,2)-1
                    neltype=npijkm(i,1)
                    mat=nmprop(i)
                           e0=prop(mat,1)
                           g0=prop(mat,2)
                           t0=prop(mat,3)
                           r0=prop(mat,4)
                    call dmat(e0,g0,dd)
c
                        kmax=neltype-1
                        do k=1,kmax
                           node=npijkm(i,1+k)
                           xyz(1,k)=xord0(node)
                           xyz(2,k)=yord0(node)
                           xyz(3,k)=zord0(node)
c
                           uvw(1,k)=dispful(node,1)
                           uvw(2,k)=dispful(node,2)
                           uvw(3,k)=dispful(node,3)
                        enddo
c
*                   write(iout,*) '@@ Elem: ',i
                    ired_int=3
                    ipmax=ired_int**3
                    call convert_HEX20_non (neltype,xyz,uvw,dd,
     &                              ired_int,
     &                strain_ip, stress_ip,stressc_ip,force_n )
c
                    if (ioutput .eq. 421 .OR. ioutput .eq. 431      
     &             .OR. ioutput .eq. 422 .OR. ioutput .eq. 432      
     &             .OR. ioutput .eq. 441 .OR. ioutput .eq. 442      
     &                                                        ) then
c                       no need for nodal averaging
                        goto 359
                    endif
c
c                   average nodal values
                    if (ioutput .eq. 424) then
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=strain_ip(k,j)
                           enddo 
                        enddo 
                    elseif (ioutput .eq. 434) then
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=stress_ip(k,j)
                           enddo 
                        enddo 
                    elseif (ioutput .eq. 444) then
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=stressc_ip(k,j)
                           enddo 
                        enddo 
                    endif
                    call node_HEX20(ssi,ssn)
                    do j=1,6
                       se_av(j)=se_av(j)+ssn(klocn,j)
                    enddo
c
 350             continue
c                bottom of elements
c                average nodal values
                 attach=mout
                 if (attach .lt. 0.5) attach=1
                 do j=1,6
                    se_av(j)=se_av(j)/attach
                 enddo
ctag3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                OUTput section
 359             continue
                 ired_int=3
                 ipmax=ired_int**3
                 if (ioutput .eq. 0) then
                     return
c
*                elseif (ioutput .eq. 51) then
*                        write(idyn,555) time, (strs(melm,j),j=1,6) 
c
                 elseif (ioutput .eq. 111) then
c                    movie
                     call scan_line(npijkm,maxelem,iscan,linemax)
*                    do ii=1,linemax
*                       write(ilog,*) (iscan(ii,jj), jj=1,3)
*                    enddo
                     kdiv=2
                     kline=linemax*kdiv
                     write(iout,*) kline,gnum1,ikont
                     call mov_shapes(npijkm,maxelem,
     &                             xord0,yord0,zord0,
     &                             dispful,ioutput,iscan,linemax,kdiv)
c
                 elseif (ioutput .eq. 411) then
                     write(idyn,25) time,(dispful(inode,j), j=1,3),
     &                                   ' 0 0 0 '
c
                 elseif (ioutput .eq. 421) then
                     write(idyn,24) time,(strain_ip(j,ip_n), j=1,6)
c
                 elseif (ioutput .eq. 422) then
                     do j=1,ipmax
                        do k=1,6
                           ssi(j,k)=strain_ip(k,j)
                        enddo 
                     enddo 
                     call node_HEX20(ssi,ssn)
                     klocn=mnod(1,2)-1
                     write(idyn,24) time,(ssn(klocn,j), j=1,6)
c
                 elseif (ioutput .eq. 424) then
                     write(idyn,24) time,(se_av(j), j=1,6)
c
                 elseif (ioutput .eq. 431) then
                     write(idyn,24) time,(stress_ip(j,ip_n), j=1,6)
c
                 elseif (ioutput .eq. 432) then
                     do j=1,ipmax
                        do k=1,6
                           ssi(j,k)=stress_ip(k,j)
                        enddo 
                     enddo 
                     call node_HEX20(ssi,ssn)
                     klocn=mnod(1,2)-1
                     write(idyn,24) time,(ssn(klocn,j), j=1,6)
c
                 elseif (ioutput .eq. 434) then
                     write(idyn,24) time,(se_av(j), j=1,6)
c
                 elseif (ioutput .eq. 441) then
                     write(idyn,24) time,(stressc_ip(j,ip_n), j=1,6)
c
                 elseif (ioutput .eq. 442) then
                     do j=1,ipmax
                        do k=1,6
                           ssi(j,k)=stressc_ip(k,j)
                        enddo 
                     enddo 
                     call node_HEX20(ssi,ssn)
                     klocn=mnod(1,2)-1
                     write(idyn,24) time,(ssn(klocn,j), j=1,6)
c
                 elseif (ioutput .eq. 444) then
                     write(idyn,24) time,(se_av(j), j=1,6)
ctag3
                 endif
                 goto 300
                 write(*,'(a)')' @@ '
c
             endif
c            bottom all snaps
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
 29            format(1x,a,8(i5,1x))
 23            format(1x,i5,2x,6(g12.6,1x))
 24            format(1x,7(g12.6,1x))
 25            format(1x,4(g12.6,1x),a)
 555           format(1x,13(g12.6,1x))
 83            format(1x,i5,1x,i5,1x,6(g12.6,1x))
 84            format(1x,5x,1x,i5,1x,6(g12.6,1x))
 81            format(1x,30(g12.6,1x))
c
 999  continue
      return
      end
c
c
c     CONVERT displacements to stresses etc
      subroutine convert_HEX20_non (neltype,
     &                            xyz,uvw,dd,
     &                   ired_int,
     &                   strain,stress,stressc,force )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
*        common /gauss_int/ xg,wgt
c
         real*8  u(60),xyz(3,20),uvw(3,20)
         real*8 dd(6,6),hh(20),be(6,60)
         real*8 xg(4,4),wgt(4,4),bb(6,60),ss(8),ee(8),ud(9),bd(9,60)
         real*8 stress(6,27), strain(6,27),force(3,20),ek(60,60)
         real*8 stressc(6,27)
         real*8 def(3,3),skk(3,3),scc(3,3),cc(3,3)
c        integer*4 icoord(8,3)
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
c        data icoord/1,2,2,1,1,2,2,1
c    &              ,1,1,2,2,1,1,2,2
c    &              ,1,1,1,1,2,2,2,2 /
c
        iout=23
         kmax=neltype-1
         nn=0
         do i=1,kmax
            do j=1,3
               nn=nn+1
               u(nn)=uvw(j,i)
*                  write(iout,*) nn,u(nn)
            enddo
         enddo
c
         nint=3
         nint=2
         nint=ired_int
         z9=0.0
         ip_n=0
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
*        do 80 ii=1,8
*           irr=(-1)**lx
*           iss=(-1)**ly
*           itt=(-1)**lz
            ip_n=ip_n+1
*           write(iout,*) irr,iss,itt
*            lx = icoord(ii,1)
*            ly = icoord(ii,2)
*            lz = icoord(ii,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
c
c                 deriv oper and jacobian
*                 call  stdm20(xyz,hh,bb,bd,det,ri,si,ti,nel,log)
                  call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
*                 write(iout,*) ' Bb ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (Bb(i,j), j=1,60)
*                 enddo
*                 write(iout,*) ' Bd '
*                 do i=1,9
*                    write(iout,82) (Bd(i,j), j=1,24)
*                 enddo
                  do i=1,9
                     sum=0.0
                     do k=1,kmax*3
                        sum=sum+Bd(i,k)*u(k)
                     enddo
                     ud(i)=sum
                  enddo
*                 do i=1,6
*                    sum=0.0
*                    do k=1,kmax*3
*                       sum=sum+Bb(i,k)*u(k)
*                    enddo
*                    ee(i)=sum
*                 enddo
*                 write(iout,*) ' u  '
*                    write(iout,83) ( u(j), j=1,kmax*3)
*                    write(iout,86) ( ud(j), j=1,9)
         z9=1.0
         ee(1)= ud(1) + z9*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
         ee(2)= ud(5) + z9*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
         ee(3)= ud(9) + z9*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
         ee(4)= ud(2)+ud(4) + z9*(ud(1)*ud(2)+ud(4)*ud(5) +ud(7)*ud(8))
         ee(5)= ud(6)+ud(8) + z9*(ud(2)*ud(3)+ud(5)*ud(6) +ud(8)*ud(9))
         ee(6)= ud(3)+ud(7) + z9*(ud(1)*ud(3)+ud(4)*ud(6) +ud(7)*ud(9))
cctag2
*                       write(iout,* ) ' strain Exx ',ip_n
*                    do j=1,6
*                       write(iout,82) ee(j) 
*                    enddo
*                    do j=1,6
*                       sum = 0.0  
*                       do kk=1,24
*                          sum = sum + bb(j,kk)*u(kk)
*                       enddo
*                       ee0 = sum
*                       write(iout,82) ee(j),ee0 
*                    enddo
*                 write(iout,82) (ee(j)*1.0e6, j=1,6)
                  do i=1,6
                     sum = 0.0  
                     do kk=1,6
                        sum = sum + dd(i,kk)*ee(kk)
                     enddo
                     ss(i)=sum
                  enddo
*                 write(iout,82) (ss(j), j=1,6)
*                 s1 = ss(1)+ss(2)+ss(3)
*                 s2 = ss(1)*ss(2) +ss(2)*ss(3) +ss(3)*ss(1) 
*    &                -ss(4)**2 - ss(5)**2 - ss(6)**2
*                 svm = sqrt(s1*s1 - 3.0*s2)
c                 write(iout,82) (ss(j), j=1,6),svm
*                 write(iout,*) '@@ SvM: ',ii,svm
*          write(*,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (dk(i,j), j=1,24)
*          enddo
                      do k=1,6
                         stress(k,ip_n)=ss(k)
                         strain(k,ip_n)=ee(k)
                      enddo
c
c                     Cauchy stress
                      skk(1,1)=ss(1)
                      skk(2,2)=ss(2)
                      skk(3,3)=ss(3)
                      skk(1,2)=ss(4)
                      skk(1,3)=ss(5)
                      skk(2,3)=ss(6)
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
*                  write(ilog,*)'@@ J: ', zjj
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
*                           write(iout,*) ' [ K stress] ',ip_n
*                           do i=1,3
*                              write(iout,82) (skk(i,j), j=1,3)
*                           enddo
*                           write(iout,*) ' [ du/dx   ] ',ip_n
*                           do i=1,3
*                              write(iout,82) (def(i,j), j=1,3)
*                           enddo
*                           write(iout,*) ' [ C stress] ',ip_n
*                           do i=1,3
*                              write(iout,82) (scc(i,j), j=1,3)
*                           enddo
                         stressc(1,ip_n)=scc(1,1)
                         stressc(2,ip_n)=scc(2,2)
                         stressc(3,ip_n)=scc(3,3)
                         stressc(4,ip_n)=scc(1,2)
                         stressc(5,ip_n)=scc(1,3)
                         stressc(6,ip_n)=scc(2,3)
*                           write(iout,*) ' [  stressC] ',ip_n
*                           do i=1,6
*                              write(iout,82) stressc(i,ip_n)
*                           enddo
 80      continue
*          write(iout,*) ' [stress] '
*          do i=1,27
*             write(iout,82) (stress(k,i), k=1,6)
*          enddo
c
         call elmstf_HEX20(dd,xyz,ired_int,ek)
         nn=0
         do i=1,kmax
            do j=1,3
               nn=nn+1
               sum=0.0
               do k=1,kmax*3
                  sum=sum+ek(nn,k)*u(k)
               enddo
               force(j,i)=sum
            enddo
         enddo
*          write(iout,*) ' [force] '
*          do i=1,20
*             write(iout,82) (force(j,i), j=1,3)
*          enddo
c
 82      format(1x,70(g12.6,1x))
 83      format(1x, 3(g12.6,1x))
 86      format(1x, 6(g12.6,1x))
c
c
      return
      end
c
c
c     SHAPES  nodal values
      subroutine mov_shapes(npijkm,maxelem,
     &                    xord,yord,zord,
     &                    dispful,ioutput,iscan,linemax,kdiv)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21),icon(12,3),icon8(6,3),nn(8)
         integer icon4(3,3)
         integer iscan(2000,3)
c
*        real*8  strs(nnp,6), strn(nnp ,6)
         real*8  xord(nnp), yord(nnp), zord(nnp)
         real*8  dispful(nnp,3)
*        real*8  sav(6),eav(6),x(12),y(12),z(12)
*        data icon / 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8
*    &              ,2, 3, 4, 1, 5, 6, 7, 8, 6, 7, 8, 5
*    &              ,9,10,11,12,13,14,15,16,17,18,19,20/
*        data icon8/ 1, 3, 1, 2, 3, 4  
*    &              ,2, 4, 5, 6, 7, 8
*    &              ,3, 1, 8, 5, 6, 7 /
*        data icon4/ 1, 3, 3  
*    &              ,2, 4, 1
*    &              ,3, 2, 4 /
c
ctag31
*        do i=1,nnp
*          write(iout,81) (dispful(i,kk), kk=1,3)
*        enddo
*        rewind(iout)
         do i=1,linemax
            if (kdiv .eq. 1) then
                do k=1,3
                   nodek=iscan(i,k)
                   x1  =xord(nodek)
                   y1  =yord(nodek)
                   z1  =zord(nodek)
                   write(iout,81) x1,y1,z1,
     &                           (dispful(nodek,kk), kk=1,3),' 0 0 0'
                enddo
            else
               node1=iscan(i,1)
               node2=iscan(i,3)
               node3=iscan(i,2)
               x1  =xord(node1)
               y1  =yord(node1)
               z1  =zord(node1)
                 x2  =xord(node2)
                 y2  =yord(node2)
                 z2  =zord(node2)
               x3  =xord(node3)
               y3  =yord(node3)
               z3  =zord(node3)
                    u1  =dispful(node1,1)
                    v1  =dispful(node1,2)
                    w1  =dispful(node1,3)
                      u2  =dispful(node2,1)
                      v2  =dispful(node2,2)
                      w2  =dispful(node2,3)
                    u3  =dispful(node3,1)
                    v3  =dispful(node3,2)
                    w3  =dispful(node3,3)
c
c              interpolate
               dr=2.0/kdiv
               do k=1,kdiv
                  r=-1+dr*(k-1)
                  h1 = 0.5*(1-r) - 0.5*(1-r*r)
                  h2 = 0.5*(1+r) - 0.5*(1-r*r)
                  h3 =                 (1-r*r)
                  xx1 = x1*h1 + x2*h2 + x3*h3
                  yy1 = y1*h1 + y2*h2 + y3*h3
                  zz1 = z1*h1 + z2*h2 + z3*h3
                      uu1 = u1*h1 + u2*h2 + u3*h3
                      vv1 = v1*h1 + v2*h2 + v3*h3
                      ww1 = w1*h1 + w2*h2 + w3*h3
                  r=-1+dr*(k-0.5)
                  h1 = 0.5*(1-r) - 0.5*(1-r*r)
                  h2 = 0.5*(1+r) - 0.5*(1-r*r)
                  h3 =                 (1-r*r)
                  xx2 = x1*h1 + x2*h2 + x3*h3
                  yy2 = y1*h1 + y2*h2 + y3*h3
                  zz2 = z1*h1 + z2*h2 + z3*h3
                      uu2 = u1*h1 + u2*h2 + u3*h3
                      vv2 = v1*h1 + v2*h2 + v3*h3
                      ww2 = w1*h1 + w2*h2 + w3*h3
                  r=-1+dr*(k-0)
                  h1 = 0.5*(1-r) - 0.5*(1-r*r)
                  h2 = 0.5*(1+r) - 0.5*(1-r*r)
                  h3 =                 (1-r*r)
                  xx3 = x1*h1 + x2*h2 + x3*h3
                  yy3 = y1*h1 + y2*h2 + y3*h3
                  zz3 = z1*h1 + z2*h2 + z3*h3
                      uu3 = u1*h1 + u2*h2 + u3*h3
                      vv3 = v1*h1 + v2*h2 + v3*h3
                      ww3 = w1*h1 + w2*h2 + w3*h3
                  write(iout,82) xx1,yy1,zz1,
     &                        uu1,vv1,ww1,' 0 0 0'
                  write(iout,82) xx2,yy2,zz2,
     &                        uu2,vv2,ww2,' 0 0 0'
                  write(iout,82) xx3,yy3,zz3,
     &                        uu3,vv3,ww3,' 0 0 0'
               enddo
            endif
         enddo
         icont=ioutput
c        displacements do not need averaging
*        if (icont .eq. 111) then
*            do n=1,nel
*               neltype=npijkm(n,1)
*               if (neltype .eq. 21) then
c                   20-noded hex
*                   do j=1,12
*                      do kkk=1,3
*                         if (kkk .eq. 1) k=1
*                         if (kkk .eq. 2) k=3
*                         if (kkk .eq. 3) k=2
*                         nodek=npijkm(n,1+icon(j,k))
*                         x1  =xord(nodek)
*                         y1  =yord(nodek)
*                         z1  =zord(nodek)
*                         write(iout,81) x1,y1,z1,
*    &                          (dispful(nodek,kk), kk=1,3),' 0 0 0'
*                      enddo
*                   enddo
*               elseif (neltype .eq. 9) then
c                   8-noded hex
*                   do j=1,6
*                      do k=1,3
*                         nodek=npijkm(n,1+icon8(j,k))
*                         x1  =xord(nodek)
*                         y1  =yord(nodek)
*                         z1  =zord(nodek)
*                         write(iout,81) x1,y1,z1,
*    &                          (dispful(nodek,kk), kk=1,3),' 0 0 0'
*                      enddo
*                   enddo
*               elseif (neltype .eq. 5) then
c                   4-noded tetra
*                   do j=1,3
*                      do k=1,3
*                         nodek=npijkm(n,1+icon4(j,k))
*                         x1  =xord(nodek)
*                         y1  =yord(nodek)
*                         z1  =zord(nodek)
*                         write(iout,81) x1,y1,z1,
*    &                          (dispful(nodek,kk), kk=1,3),' 0 0 0'
*                      enddo
*                   enddo
*               endif
*            enddo
*        endif
*        return
c
 81      format(1x,20(g12.6,1x))
 82      format(1x, 6(g12.6,1x),1x,a)
c
 999     continue
*        close(itmp)
      return
      end
c
c
c     scan unique plotting lines  
      subroutine scan_line(npijkm,maxelem,iscan,linemax
     &                                  )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21),icon(12,3),icon8(6,3)
         integer icon4(3,3)
         integer iscan(2000,3)
         data icon / 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8
     &              ,2, 3, 4, 1, 5, 6, 7, 8, 6, 7, 8, 5
     &              ,9,10,11,12,13,14,15,16,17,18,19,20/
         data icon8/ 1, 3, 1, 2, 3, 4  
     &              ,2, 4, 5, 6, 7, 8
     &              ,3, 1, 8, 5, 6, 7 /
         data icon4/ 1, 3, 3  
     &              ,2, 4, 1
     &              ,3, 2, 4 /
c
ctag31
             nn=0
             do n=1,nel
                neltype=npijkm(n,1)
                if (neltype .eq. 21) then
c                   20-noded hex
                    do j=1,12
                       node1 = npijkm(n,1+icon(j,1))
                       node2 = npijkm(n,1+icon(j,2))
                       nodem = npijkm(n,1+icon(j,3))
                       if (node2 .lt. node1) then
                           jj=node1
                           node1=node2
                           node2=jj
                       endif
                       do i=1,nn
                          if (node1 .eq. iscan(i,1)) then
                              if (node2 .eq. iscan(i,3)) then
c                                 line exists
                                  goto 100
                              endif
                          endif
                       enddo
c                      new line
                       nn=nn+1
                       iscan(nn,1)=node1
                       iscan(nn,2)=nodem
                       iscan(nn,3)=node2
 100                   continue
                    enddo
                endif
             enddo
             linemax=nn
c
 81      format(1x,20(g12.6,1x))
c
 999     continue
      return
      end
c
c
      subroutine node_HEX20(ssi,ssn)
         implicit real*8 (a-h,o-z)
         real*8  ssi(27,6),aai(8,8),ssn(20,6),aa20(27,20)
c
         data aai /
     &        2.54904,    -.683013,     .183013,    -.683013,   
     &              -.683013,     .183013,    -.490381E-01, .183013,   
     &       -.683013,     2.54904,    -.683013,     .183013, 
     &               .183013,    -.683013,     .183013,  -.490381E-01,
     &        .183013,    -.683013,     2.54904,    -.683013, 
     &              -.490381E-01, .183013,    -.683013,     .183013,    
     &       -.683013,     .183013,    -.683013,     2.54904, 
     &               .183013,    -.490381E-01, .183013,    -.683013,    
     &       -.683013,     .183013,    -.490381E-01, .183013,  
     &              2.54904,    -.683013,     .183013,    -.683013,    
     &        .183013,    -.683013,     .183013,    -.490381E-01,
     &              -.683013,     2.54904,    -.683013,     .183013,    
     &       -.490381E-01, .183013,    -.683013,     .183013, 
     &               .183013,    -.683013,     2.54904,    -.683013,    
     &        .183013,    -.490381E-01, .183013,    -.683013, 
     &              -.683013,     .183013,    -.683013,     2.54904 /  
c
         data aa20 /
     &  2.37499, -.31246, -.12559, -.31246, -.48824,  .31481, -.12559,
     &   .31481, -.16145, -.31246, -.48824,  .31481, -.48824, -.29630,
     &   .22898,  .31481,  .22898, -.16902, -.12559,  .31481, -.16145,
     &   .31481,  .22898, -.16902, -.16145, -.16902,  .11575,
     &  -.12559, -.31246, 2.37499,  .31481, -.48824, -.31246, -.16145,
     &   .31481, -.12559,  .31481, -.48824, -.31246,  .22898, -.29630,
     &  -.48824, -.16902,  .22898,  .31481, -.16145,  .31481, -.12559,
     &  -.16902,  .22898,  .31481,  .11575, -.16902, -.16145,
     &  -.16145,  .31481, -.12559,  .31481, -.48824, -.31246, -.12559,
     &  -.31246, 2.37499, -.16902,  .22898,  .31481,  .22898, -.29630,
     &  -.48824,  .31481, -.48824, -.31246,  .11575, -.16902, -.16145,
     &  -.16902,  .22898,  .31481, -.16145,  .31481, -.12559,
     &  -.12559,  .31481, -.16145, -.31246, -.48824,  .31481, 2.37499,
     &  -.31246, -.12559,  .31481,  .22898, -.16902, -.48824, -.29630,
     &   .22898, -.31246, -.48824,  .31481, -.16145, -.16902,  .11575,
     &   .31481,  .22898, -.16902, -.12559,  .31481, -.16145,
     &  -.12559,  .31481, -.16145,  .31481,  .22898, -.16902, -.16145,
     &  -.16902,  .11575, -.31246, -.48824,  .31481, -.48824, -.29630,
     &   .22898,  .31481,  .22898, -.16902, 2.37499, -.31246, -.12559,
     &  -.31246, -.48824,  .31481, -.12559,  .31481, -.16145,
     &  -.16145,  .31481, -.12559, -.16902,  .22898,  .31481,  .11575,
     &  -.16902, -.16145,  .31481, -.48824, -.31246,  .22898, -.29630,
     &  -.48824, -.16902,  .22898,  .31481, -.12559, -.31246, 2.37499,
     &   .31481, -.48824, -.31246, -.16145,  .31481, -.12559,
     &   .11575, -.16902, -.16145, -.16902,  .22898,  .31481, -.16145,
     &   .31481, -.12559, -.16902,  .22898,  .31481,  .22898, -.29630,
     &  -.48824,  .31481, -.48824, -.31246, -.16145,  .31481, -.12559,
     &   .31481, -.48824, -.31246, -.12559, -.31246, 2.37499,
     &  -.16145, -.16902,  .11575,  .31481,  .22898, -.16902, -.12559,
     &   .31481, -.16145,  .31481,  .22898, -.16902, -.48824, -.29630,
     &   .22898, -.31246, -.48824,  .31481, -.12559,  .31481, -.16145,
     &  -.31246, -.48824,  .31481, 2.37499, -.31246, -.12559,
     &   .32628, 1.28439,  .32628, -.27072,  .05556, -.27072,  .11111,
     &  -.19444,  .11111, -.27072,  .05556, -.27072, -.22222, -.11111,
     &  -.22222,  .15961,  .05556,  .15961,  .11111, -.19444,  .11111,
     &   .15961,  .05556,  .15961, -.10405, -.00661, -.10405,
     &   .11111, -.27072,  .32628, -.19444,  .05556, 1.28439,  .11111,
     &  -.27072,  .32628,  .15961, -.22222, -.27072,  .05556, -.11111,
     &   .05556,  .15961, -.22222, -.27072, -.10405,  .15961,  .11111,
     &  -.00661,  .05556, -.19444, -.10405,  .15961,  .11111,
     &   .11111, -.19444,  .11111, -.27072,  .05556, -.27072,  .32628,
     &  1.28439,  .32628,  .15961,  .05556,  .15961, -.22222, -.11111,
     &  -.22222, -.27072,  .05556, -.27072, -.10405, -.00661, -.10405,
     &   .15961,  .05556,  .15961,  .11111, -.19444,  .11111,
     &   .32628, -.27072,  .11111, 1.28439,  .05556, -.19444,  .32628,
     &  -.27072,  .11111, -.27072, -.22222,  .15961,  .05556, -.11111,
     &   .05556, -.27072, -.22222,  .15961,  .11111,  .15961, -.10405,
     &  -.19444,  .05556, -.00661,  .11111,  .15961, -.10405,
     &   .32628, -.27072,  .11111, -.27072, -.22222,  .15961,  .11111,
     &   .15961, -.10405, 1.28439,  .05556, -.19444,  .05556, -.11111,
     &   .05556, -.19444,  .05556, -.00661,  .32628, -.27072,  .11111,
     &  -.27072, -.22222,  .15961,  .11111,  .15961, -.10405,
     &   .11111, -.27072,  .32628,  .15961, -.22222, -.27072, -.10405,
     &   .15961,  .11111, -.19444,  .05556, 1.28439,  .05556, -.11111,
     &   .05556, -.00661,  .05556, -.19444,  .11111, -.27072,  .32628,
     &   .15961, -.22222, -.27072, -.10405,  .15961,  .11111,
     &  -.10405,  .15961,  .11111,  .15961, -.22222, -.27072,  .11111,
     &  -.27072,  .32628, -.00661,  .05556, -.19444,  .05556, -.11111,
     &   .05556, -.19444,  .05556, 1.28439, -.10405,  .15961,  .11111,
     &   .15961, -.22222, -.27072,  .11111, -.27072,  .32628,
     &   .11111,  .15961, -.10405, -.27072, -.22222,  .15961,  .32628,
     &  -.27072,  .11111, -.19444,  .05556, -.00661,  .05556, -.11111,
     &   .05556, 1.28439,  .05556, -.19444,  .11111,  .15961, -.10405,
     &  -.27072, -.22222,  .15961,  .32628, -.27072,  .11111,
     &   .11111, -.19444,  .11111,  .15961,  .05556,  .15961, -.10405,
     &  -.00661, -.10405, -.27072,  .05556, -.27072, -.22222, -.11111,
     &  -.22222,  .15961,  .05556,  .15961,  .32628, 1.28439,  .32628,
     &  -.27072,  .05556, -.27072,  .11111, -.19444,  .11111,
     &  -.10405,  .15961,  .11111, -.00661,  .05556, -.19444, -.10405,
     &   .15961,  .11111,  .15961, -.22222, -.27072,  .05556, -.11111,
     &   .05556,  .15961, -.22222, -.27072,  .11111, -.27072,  .32628,
     &  -.19444,  .05556, 1.28439,  .11111, -.27072,  .32628,
     &  -.10405, -.00661, -.10405,  .15961,  .05556,  .15961,  .11111,
     &  -.19444,  .11111,  .15961,  .05556,  .15961, -.22222, -.11111,
     &  -.22222, -.27072,  .05556, -.27072,  .11111, -.19444,  .11111,
     &  -.27072,  .05556, -.27072,  .32628, 1.28439,  .32628,
     &   .11111,  .15961, -.10405, -.19444,  .05556, -.00661,  .11111,
     &   .15961, -.10405, -.27072, -.22222,  .15961,  .05556, -.11111,
     &   .05556, -.27072, -.22222,  .15961,  .32628, -.27072,  .11111,
     &  1.28439,  .05556, -.19444,  .32628, -.27072,  .11111  /
c
ctag9
c          [A] is stored as transpose
           do j=1,6
              do i=1,20
                 sum=0.0
                 do k=1,27
                    sum = sum + aa20(k,i)*ssi(k,j)
                 enddo
                 ssn(i,j) = sum
              enddo
           enddo
*                 do j=1,8
*                    write(iout,83) n,j,(ssn(j,k), k=1,6)
*                 enddo 
c
 83      format(1x,i5,1x,i5,1x,6(g12.6,1x))
      return
      end
c
c
c     AVERAGE element values to get nodal values
      subroutine node_avg_HEX20(npijkm,maxelem,wk_elm_n,strn)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21)
         real*8  strn(nnp,6)
         real*8  wk_elm_n(nel,20,6)
         real*8  eav(6)
ctag1
c        NODAL AVERAGE
         ipmax=27
         do 60 node=1,nnp
            do k=1,6
               eav(k)=0.0
            enddo
c
c           find elements attached to this node
            attach=0.0
            do 70 jelm=1,nel
               neltype=npijkm(jelm,1)
               kmax=neltype-1
               do k=1,kmax
                  if (node .eq. npijkm(jelm,1+k)) then
                      nth=k
                      goto 72
                  endif
               enddo
               goto 70
c
 72            continue  
               attach=attach+1
               do j=1,6
                  eav(j)=eav(j)+wk_elm_n(jelm,nth,j)
               enddo
c
 70         continue
            if (attach .lt. 0.5) attach=1
            do k=1,6
               strn(node,k)=eav(k)/attach
            enddo
c
 52         format(1x,i5,1x,8(g12.5,1x))
 60      continue
c        bottom of nodal loop
c
      return
      end
c
c
      subroutine follow_load(nfollow,ireadmax,nfoll,maxnode,
     &                    xord,yord,zord,angle)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer nfollow(ireadmax,4)
         real*8  xord(maxnode),yord(maxnode),zord(maxnode)
         real*8  angle(ireadmax,3)
c
c                 plate
         do n=1,nfoll
            n1=nfollow(n,1)
            n2=nfollow(n,2)
            n3=nfollow(n,3)
            n4=nfollow(n,4)
                   x1=xord(n1)
                   x2=xord(n2)
                   x3=xord(n3)
*                  x4=xord(n4)
                   y1=yord(n1)
                   y2=yord(n2)
                   y3=yord(n3)
*                  y4=yord(n4)
                     z1=zord(n1)
                     z2=zord(n2)
                     z3=zord(n3)
*                    z4=zord(n4)
*                  dx = xord(n1) - xord(i1)
*                  dy = yord(n1) - yord(i1)
*                  dz = zord(j1) - zord(i1)
*                  xl0=sqrt(dx*dx+dy*dy+dz*dz) 
*                    xll1=dx/xl0 
*                    xmm1=dy/xl0 
*                    xnn1=dz/xl0 
*                  triad_n1(1,1,n)=xll1
*                  triad_n1(2,1,n)=xmm1
*                  triad_n1(3,1,n)=xnn1
c
c               determine vector area
                axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
                ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
                azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
                area=sqrt( axy*axy + ayz*ayz + azx*azx)
                xll3=ayz/area
                xmm3=azx/area
                xnn3=axy/area
                   angle(n,1)=xll3
                   angle(n,2)=xmm3
                   angle(n,3)=xnn3
*                  triad_n1(1,3,n)=xll3
*                  triad_n1(2,3,n)=xmm3
*                  triad_n1(3,3,n)=xnn3
c
*               xll2=xmm3*xnn1 - xnn3*xmm1
*               xmm2=xnn3*xll1 - xll3*xnn1
*               xnn2=xll3*xmm1 - xmm3*xll1
*                  triad_n1(1,2,n)=xll2
*                  triad_n1(2,2,n)=xmm2
*                  triad_n1(3,2,n)=xnn2
          enddo
c
      return
      end
c
