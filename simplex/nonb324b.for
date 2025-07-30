c
c
c 4740
c     NON-liner RUBber 3D TL incremental analysis 
      subroutine nonabt_vis_TL(stf,geom,wk,pload, disp, fmag,
     &                   dispt , 
     &                   stfeff,mass,vel,acc,uup,vvp,aap
     &                  ,dmp,wk22
     &                  ,idbc,npijkm,maxnode,maxelem
     &                  ,iprofv, iprofh, nloc
     &                  ,xyz0,maxxyz
     &                  ,iglobal, ippp
     &                  ,prop, nmprop,nelt
     &                  ,respceb,isize8
     &         ,isize_elem,stresst,straint,stressi,straini,stressci
     &         ,isize_dof,dispi,dui,grav,fmag3
     &         ,dloadi,floadi,floadm
     &         ,xyzt,xyzi
     &                      )
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=5000,isize_neq=3000)
         parameter( maxdim=2000,maxdof=maxdim*3,idload=1000 )
c
         integer iprofv(neq), iprofh(neq)
     &          ,nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer nloc(neq)
         integer inod(ireadmax,2)
*        integer inod(ireadmax,2),mnod(iconelem)
c
         real*8 tforce1(maxforce,2), tforce2(5000,2),tforce3(5000,2)
         real*8 velout(ireadmax)
         real*8 stf(maxstiff),mass(maxstiff),wk(neq), geom(maxstiff)
     &          ,pload(neq ),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 stfeff(maxstiff),dmp(maxstiff),wk22(neq)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(maxxyz)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 xyzi(maxxyz)
         real*8 dispi(isize_dof),dui(isize_dof),grav(isize_dof)
         real*8 fmag3(isize_dof)
         real*8 dispt(neq)
         real*8 dsumd,dsumz,dtol,pnorm,pnormi,wnorm,wnormi
*        real*8 dispfult(900,3), dispfuli(900,3)
         real*8 dispfult(maxdim,3), dispfuli(maxdim,3)
c
         real*8 dloadi(isize_dof),floadi(isize_dof),floadm(isize_dof)
c        real*8 dloadi(idload   ),floadi(idload   ),floadm(idload   )
         real*8 eforce(3,20)
c
c        eigenstuff
         real*8 respceb(isize8)
         real*8 wk0(isize_neq),wk1(isize_neq),wk2(isize_neq)
         real*8 wk_eig(100)
ctag1
c
cccccccccccccccccccccccccccc
c        tissue   
*        parameter( isize_elem=200)
         real*8 stresst(isize_elem,6,27),straint(isize_elem,6,27)
         real*8 stressi(isize_elem,6,27)
         real*8 straini(isize_elem,6,27)
         real*8 stressci(isize_elem,6,27)
         integer*4 npijkmv(1000,2),npiv(1000),npjv(1000),jbc(1000)
         real*8 ploadv(1000),veli(1000)
cccccccccccccccccccccccccccc
c
*        write(ilog,*)'@@  nonabt_vis_TL ver310  October 2010 '
         write(ilog,*)'@@  nonabt_vis_TL ver324  October 2010 '
c
c
                 iout3=133
                 open(unit=iout3,file=fdn//'simplex.333')
                 rewind(iout3)
c
                 itmp4=129
              open(unit=itmp4,file=fdn//'stadyn.tm4',form='unformatted')
                 rewind(itmp4)
c
        call visco_file(npijkmv,nelv,etav)
           write(*,*) '@@ '
        do i=1,nelv
           write(*,*)i, (npijkmv(i,j),j=1,2)
           npiv(i)=npijkmv(i,0+1)
           npjv(i)=npijkmv(i,0+2)
           do k=1,3
              write(*,*) k, (idbc(npijkmv(i,j),k),j=1,2)
           enddo
        enddo
*       jbc=0
*       do i=1,nnp
*          ind=i*6
*          jbc(ind-5)=idbc(i,1)
*          jbc(ind-4)=idbc(i,2)
*          jbc(ind-3)=idbc(i,3)
*       enddo
*       do i=1,nnp
*          write(iout,87) i,(jbc( (i-1)*6 + j ),j=1,6)
*       enddo
        call visco_mat( dmp,nelv,npiv, npjv 
     &                 ,xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1)
     &                 ,idbc,maxnode
     &                 ,iprof,nloc)
        dmp = etav*dmp
*       stop
c
c
c        algorithm parameters
           write(*,*) maxstiff,' ms'
*          stop
         beta0=1.0
         iter_max=60
         dtol=0.0001
         gamma0=1.0
         kref=2
c
           do i=1,nnp
              xyzt(i)=xyz0(i)
              xyzt(maxnode+i)=xyz0(maxnode+i)
              xyzt(maxnode*2+i)=xyz0(maxnode*2+i)
c              write(iout,82) xyzt(i),xyzt(maxnode+i),xyzt(maxnode*2+i)
           enddo
c
c        Get things ready for time integration loop
         write(*,*)'INPUT:  time inc  |  # of incs  |  print count |',
     &                                        ' snap count '
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,iprcnt,isnap
               if (isnap .eq. 0) isnap=10000
c
         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'
c
         nsnapmax=npt
         write(*,*)' '
         write(*,*)'    u = u + beta*du   K = Ke + gamma Kg'
         write(*,*)'INPUT:   beta   |  gamma  |  ramp    '
         write(*,*)' suggest:   1        1         5 '
         call zzwrt('  -->  ')
         read(ikbd,*) beta0,gamma0,maxramp0  
         write(ilog,*) beta0,gamma0,maxramp0    ,'  ::b g r '
         write(*,*)' '
         write(*,*)'@@ algor 1=full N-R  2=mod N-R '
         write(*,*)'INPUT:   algor | iter max |  tolerance  (<1.0E-5)'
         call zzwrt('  -->  ')
         read(ikbd   ,*) imodify,iter_max,dtol
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
c
         write(*,*)' '
         write(*,*)'CHOOSE nodal output:'
         nout=0
27       continue
         nout=nout+1
         write(*,*)'TYPE:  node# |  DoF   |  rate   <0 0 0 to end>'
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
     &                                           ,' ::eig t c v' 
c
c
c        form lumped mass
         ilump=1
!Fangbao         write(ilog,*)'@@ FORMing lumped  mass matrix '
!         call zztimer(ilog,'--> FORMmass')
         call formmass(mass(1), npijkm,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      idbc, maxnode,maxelem,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
            z1=0.0
            do i=1,neq
               z1=z1+mass(i)
            enddo
!Fangbao         write(ilog,*)'@@Total Mass: ',z1/3
         write(iout,*)'@@Total Mass: ',z1/3
c              write(iout,86) (mass(j),j=1,neq)
!         call zztimer(ilog,'<-- FORMmass')
!         call zzflush(ilog)
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
!Fangbao         write(ilog,*) '@@ Reloaded  {P}   OK'
c
c        INPUT load history from a file and interpolate later 
         call getload_hist(tforce1,maxforce,npt,deltat)
c
c        only for second load
         if (igravity_on .eq. 1) then
             write(*,*) '@@ get {Grav} history: '
             call getload_hist(tforce2,5000,npt,deltat)
         endif
         if (iload_num .eq. 3) then
             open(unit=itmp,file=fdn//'stadyn.ld3')
             rewind(itmp)
             do i=1,neq
                read(itmp,*) fmag3(i)
             enddo
             write(*   ,*) '@@ Reloaded  {P}_3   OK'
!Fangbao             write(ilog,*) '@@ Reloaded  {P}_3   OK'
             close(itmp)
c
c            INPUT load history from a file and interpolate  later
             call getload_hist(tforce3,5000,npt,deltat)
         endif
         write(*,*)'                    {P}  | {Grav} | 3rd |     '
         write(*,*)'TYPE force scales:   P1  |  P2    |  P3 | P4  '
         call zzwrt('  -->  ')
         read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
!Fangbao         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
c
c
c
c        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
         write(ilog ,*)'@@  Beginning increment analysis   '
c
c        Initialize  times, disp
c
         write(*  ,*) '@@ INITIAL disps set to zero'
         veli=0.0
         do i= 1, neq  
              disp(i)=0.0
              dispt(i)=0.0
*             delu(i)=0.0
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
         do i=1,nel
            do j=1,6
               do k=1,27
                   straint(i,j,k)=0.0
                   stresst(i,j,k)=0.0
               enddo
            enddo
         enddo
c        all elements initially not yielded
c
         do 71 n=1,nout
            iprnode=inod(n,1)
!Fangbao            write(ilog,*)'@@ iprnode ',n,iprnode
            if (iprnode .eq. 0) then
                  velout(n) = 0.0
                  goto 71
            endif
            if (inod(n,2) .eq. 0) then
                  velout(n)=dispt(iprnode)
              endif
 71      continue
         rewind(idyn)
         itime=0
         time=0
           ziter=0
           iter=0
         zlam11=0.0
         zlam00=1.0
           k1old=1
           k2old=1
           k3old=1
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
         tf1=0
         tf2=0
         tf3=0
         deltf1=1
 201     continue
c
         write(iout,*) ' 0000   iter max ',iter,iter_max
         call zzcount_2(0,inc,icarriage)
c        BIG TIME LOOP
         do 200 itime= 1,npt
            call zzcount_2(itime,inc,icarriage)
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
            time = real(itime-1)*deltat
            if (itime .eq. 2) then
                maxramp=20
                maxramp=1 
                maxramp=maxramp0
            else
                maxramp=maxramp0
            endif
            write(iout,'(a)')'@@ '
            write(iout,*)' N: ',itime
               if (zlam11 .gt. 0.0) then
                   cc3= dcc/sqrt(zlam11/zlam00)
               else
                   cc3= dcc/sqrt( 0.01 )
               endif
!Fangbao            write(ilog,'(a)')'@@ '
!Fangbao            write(ilog,*    )'@@  cc3: ',cc3
c           write(ilog,*    )'@@  cc3: ',cc3,zlam00,zlam11
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
            tf0=tf1
            tf3=0
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            if (igravity_on .eq. 1) then
                call load_interp(tforce2,5000,k2old,time,tf2)
            endif
            if (iload_num .eq. 3) then
                call load_interp(tforce3,5000,k3old,time,tf3)
            endif
            pnorm=0.0
            do i= 1, neq  
               pload(i) =  tf1*fmag(i)*scale_p1 + tf2*grav(i)*scale_p2
     &                    +tf3*fmag3(i)*scale_p3
               pnorm = pnorm + pload(i)**2
            enddo 
*           call visco_load(npijkmv,nelv,etav
*    &                ,npijkm,idbc,maxnode,maxelem
*    &                ,xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1)
*    &                ,xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1)
*    &                ,vel,ploadv
*    &                     )
            ieqn20 = idbc(20,2)
            ieqn32 = idbc(32,2)
            ieqn41 = idbc(41,2)
*           iter=0
*           write(iout3,82) time,pload(ieqn41)
*    &                    ,ploadv(ieqn1),ploadv(ieqn2)
*    &                    ,vel(ieqn1),vel(ieqn2),iter         
            do i=1,neq
*              pload(i) = pload(i) + ploadv(i)
            enddo
c
            pnorm = sqrt(pnorm/neq)
c
c
            if (imodify .eq. 2) then
c               modified N-R
*               write(iout,*) ' before iter max ',iter,iter_max
c
            call formstif_abt(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      stresst,isize_elem,
     &                      nelt,iglobal,
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
            call formgeom_abt(geom(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,stresst,isize_elem,
     &                      nelt,iglobal,
     &                      prop,nmprop,ippp,iprofv,nloc )
ctag3
                if (itime .eq. 1) goto 290
c
c           Form effective stiffness matrix by adding geometry matrix
!Fangbao            write(ilog,*)'@@ dcc cc3 ',dcc,cc3
            do i= 1, neq  
               iloc=nloc(i)
               do 25 j= 1, iprofv(i)
                      ii=iloc+j-1
                  stfeff(iloc+j-1) = 
     &               (1+zalpha)*(stf(iloc+j-1) 
     &                  +gamma0*geom(iloc+j-1))
*    &                     +     dmp(iloc+j-1)*zgamma/(zbeta*dt)
  25           continue
               damp = (dmm*mass(i)+dkk*stfeff(iloc)) + cc3
               stfeff(iloc) =  stfeff(iloc)  
     &                     +(1/(zbeta*dt**2))*mass(i)
     &                     +(1+zalpha)*(zgamma/(zbeta*dt))*damp
c    &                     +(1+zalpha)*(zgamma/(zbeta*dt))*dmm*mass(i)
                enddo   
c
c               Decompose effective stiffness matrix
                ierror=0
            call uduCOL_D( stfeff,maxstiff,neq,ierror,iprofv,nloc,z0,z1)
                write(iout,84) '@@ min/max UDU: ',z0,z1,itime,iter
                zrat=z1/z0
!Fangbao                write(ilog,*)'@@ zrat: ',zrat
                if (ierror.eq.0) then
                    write(*,*)'ERROR: zero diagonal term'
                    return
                endif
                write(iout,*) ' iter max ',iter,iter_max
            endif
c           bottom modified N-R
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
*                 dufuli(i,j)=0.0
                  dispfuli(i,j)=dispfult(i,j)
               enddo
            enddo
            do i=1,nnp
               xyzi(i)=   xyzt(i)
               xyzi(maxnode+i)=xyzt(maxnode+i)
               xyzi(maxnode*2+i)=xyzt(maxnode*2+i)
            enddo
            do n=1,nel
               do i=1,6
                  do k=1,27
                     stressi(n,i,k) = stresst(n,i,k)
                     straini(n,i,k) = straint(n,i,k)
                  enddo
               enddo
            enddo
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
            tf_load =        tf1      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c               full N-R
c           forces from body stresses
c           send current stress = disp increment
            call body_stress_abt( stressi,isize_elem,
     &                    dispfuli,maxdim,
     &                    iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadi
     &                 )
ctag1
c
c           form load increment
            write(iout,*) ' force: ',tf1,tf_load
            pnormi=0.0
            pnorm =0.0
c                          visco x vel 
            do i= 1, neq  
*              wk(i) = vel(i)*(zgamma-zbeta)/(zbeta*dt)
               vterm = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
               wk(i) = -vterm
            enddo       
            wk22=0.0
            call AxBCOL(dmp,maxstiff,wk,wk22,
     &                                   neq,iprofv,iprofh,nloc)
            do i=1,neq  
               iloc=nloc(i)
               stfeff(iloc) = stf(iloc) +gamma0*geom(iloc)
               damp = (dmm*mass(i)+dkk*stfeff(iloc))+ cc3
               aterm = aap(i) + (  1/(zbeta*dt**2))*(dispi(i)-uup(i))
               vterm = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
               dpp = -mass(i)*aterm
     &                 - (1+zalpha)*damp*vterm 
     &                + zalpha*dmm*mass(i)*vel(i)
               ploadm = tf0*fmag(i)*scale_p1
               dloadi(i) = (1+zalpha)*pload(i)  - zalpha*ploadm
     &                     - (1+zalpha)*floadi(i) + zalpha*floadm(i)
     &                     + dpp
     &                     + wk22(i)
*    &                     + ploadv(i)
               pnormi = pnormi + dloadi(i)**2
               pnorm  = pnorm  + pload (i)**2
               write(iout,82) iter,pload(i),dloadi(i),ploadv(i),veli(i)
            enddo 
            write(iout3,82) time,pload(ieqn41),dloadi(ieqn20)
     &                          ,wk22(ieqn20),wk22(ieqn32)
     &                          ,vel (ieqn20),vel (ieqn32)
     &                          ,vel (ieqn32)-vel (ieqn20)
            pnormi=sqrt(pnormi/neq)
            pnorm =sqrt(pnorm /neq)
c           write(iout,*) ' force: ',tf1,tf_load
c
c
            if (imodify .eq. 1) then 
c               full N-R
                call formstif_abt(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfuli,maxdim,
     &                      stresst,isize_elem,
     &                      nelt,iglobal,
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
*              do i=1,neq
*                 write(iout,*) i,dloadi(i),'aft stif'
*              enddo
ctag2
c
c                SECOND assemble the geometric stiffness
                 call formgeom_abt(geom(1), npijkm,idbc,maxnode,maxelem
     &                     ,xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1)
     &                     ,dispfuli,maxdim,stressi,isize_elem
     &                     ,nelt,iglobal
     &                     ,prop,nmprop,ippp,iprofv,nloc )
ctag3
c                Form effective stiff matrix by adding geometry matrix
                 do i=1,neq  
                    iloc=nloc(i)
                    do 24 j=1,iprofv(i)
                       ii=iloc+j-1
                       stfeff(iloc+j-1) = 
     &                 (1+zalpha)*(stf(iloc+j-1) 
     &                     + gamma0*geom(iloc+j-1))
     &                     +         dmp(iloc+j-1)*zgamma/(zbeta*dt)
 24                 continue
                    damp = (dmm*mass(i)+dkk*stfeff(iloc)) + cc3
                    damp =  dmm*mass(i)                            
                    stfeff(iloc) =  stfeff(iloc)  
     &                          +(1/(zbeta*dt**2))*mass(i)
     &                          +(1+zalpha)*(zgamma/(zbeta*dt))*damp
                 enddo   
c
c                Decompose effective stiffness matrix
                 ierror=0
                 call uduCOL_D( stfeff,maxstiff,neq,ierror
     &                        ,iprofv,nloc,z0,z1)
                 write(iout,84) '@@ min/max UDU: ',z0,z1,itime,iter
                 if (ierror.eq.0) then
                     write(*,*)'ERROR: zero diagonal term'
                     return
                 endif
            endif
c           bottom full N-R
c
c           solve for Delta u
            ier1=0
            call bakCOL_D(stfeff,maxstiff,dloadi, neq, dui,ierror,
     &                                       iprofv,iprofh,nloc)
            dsumd=0.0
            do i=1,neq
               dsumd=dsumd+ abs(dui(i))
            enddo
!Fangbao            write(ilog,*)'@@  |dui| ',dsumd/neq,iter
!Fangbao            write(iout,*)'@@  |dui| ',dsumd/neq,iter
c
c           increment Ui
*           maxramp=10 
            if (iter .le. maxramp) then
                beta=iter/(maxramp+0.0)
                beta=beta*beta0
                zi=2.0**(iter  )
                z0=2.0**(maxramp  )
                beta=zi/z0
                beta=beta*beta0
            else
                beta=1.0
                beta=beta0
            endif
            do i=1,neq  
               dispi(i)= dispi(i) + beta*dui(i)
            enddo
            do i= 1, neq  
               veli(i) = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
            enddo    
*           write(iout3,82) time,pload(ieqn41),veli(ieqn1),vel(ieqn2)
*    &                          ,iter
            call visco_load(npijkmv,nelv,etav
     &                ,npijkm,idbc,maxnode,maxelem
     &                ,xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1)
     &                ,xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1)
     &                ,veli,ploadv
     &                     )
*           write(iout3,82) time,pload(ieqn41)
*    &                    ,ploadv(ieqn1),ploadv(ieqn2)
*    &                    ,veli(ieqn1),vel(ieqn2),iter
*           write(iout,*) ' beta: ',beta,1/beta
c
            call update_geom_3D(dispi,dispfuli,maxdim,idbc,maxnode,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1))
c
                isizeflag=isize_elem
                nall=-1
                iextra=0
                call update_stress_abt(
     &                   stressi,straini,stressci,
     &                   isize_elem,
     &                   dispfuli,maxdim,
     &                   npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                   prop,nmprop,ippp,nall,iextra,eforce
     &                )
ctag4
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
c
c           set values for time t+dt -> t
            do i= 1, neq  
*              delu(i)=dispi(i)-dispt(i)
               dispt(i)=dispi(i)
               vel(i) = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
               acc(i) = aap(i) + (1.0/(zbeta*dt*dt))*(dispi(i)-uup(i))
            enddo    
            deltf1=tf1-tf0
            do i=1,nnp  
               do j= 1,3  
                  dispfult(i,j)=dispfuli(i,j)
               enddo
            enddo
            do i=1,nnp
               xyzt(i)          =xyzi(i)
               xyzt(maxnode+i)  =xyzi(maxnode+i)
               xyzt(maxnode*2+i)=xyzi(maxnode*2+i)
            enddo
            do i=1,neq  
               floadm(i) = floadi(i)
            enddo 
c
            do n=1,nel
               do i=1,6
                  do k=1,27
                     stresst(n,i,k) = stressi(n,i,k)
                     straint(n,i,k) = straini(n,i,k)
                  enddo
               enddo
            enddo
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
!                call zzflush(idyn)
             endif
             if (isnap .eq. ksnap .OR. itime .eq. 1) then
                 write(isnp) time, neq, nsnapmax 
                 do i= 1, neq  
                    write(isnp) dispt(i)
                 enddo
                 ksnap=0
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
                        zlam11=wk_eig(1)
                        if (itime .eq. 1) zlam00 = zlam11
                        do k=1,ieig_vec
                           zlam=wk_eig(k)
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
                           write(itmp4) zlam, neq, ieig_max 
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
                    write(imon ,85) time,tf1s,(wk_eig(j),j=1,ieig_vec)
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
c
 200     continue
*210  continue
         write(*,'(a)')' @@'
c        END BIG TIME LOOP
c
 81   format(1x,a,1x,90(1x,g12.6))
 82   format(1x,90(1x,g12.6))
 83   format(1x,4(1x,g12.6),1x,a)
 84   format(1x,a,1x,2(1x,g12.6),1x,2(i5,1x))
 85   format(1x,90(1x,g12.6))
 86   format(1x, 6(1x,g12.6))
 87   format(1x, 20(1x,i5))
 182  format(1x,20(g12.6,1x),/,20(g12.6,1x),/,20(g12.6,1x))
 183  format(1x,i5,4x,3(g12.6,1x))
 184  format(1x,30(i2,1x))
 185  format(1x,i2,1x,': ',30(i2,1x))
c
      return
      end
c
c
      subroutine visco_file(npijkmv,nelv,etav)
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         character*40 filena
         integer*4 npijkmv(1000,2)
c
         write(*,'(a)')' @@'
         call zzwrt(' INPUT:  viscous datafilename-->  ')
         read(ikbd   ,'(a40)') filena
         write(ilog,*) filena
         open(idat, file = fdn//adjustl(filena))
         rewind(idat)
c
         read(idat,*) nelv
         do i=1,nelv
            read(idat,*) nn,(npijkmv(nn,j),j=1,2)
         enddo
         read(idat,*) nmatv
         do i=1,nmatv
            read(idat,*) nn,etav
         enddo
         close(idat)
c
      return 
      end
c
c
      subroutine visco_load(npijkmv,nelv,etav
     &                ,npijkm,idbc,maxnode,maxelem
     &                ,xord0,yord0,zord0
     &                ,xord ,yord ,zord 
     &                ,vel,ploadv
     &                     )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
c
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3)
         real*8  xord0(nnp),yord0(nnp),zord0(nnp),vel(neq)
         real*8  xord (nnp),yord (nnp),zord (nnp)  
         integer*4 npijkmv(1000,2)
         real*8  ploadv(1000),triad_ee(3,3),fb(6),fg(6)
c
c        do i=1,nnp
c           write(*,82) i,xord0(i),yord0(i)
c        enddo
            ploadv=0.0
            do i=1,nelv
               node1 = npijkmv(i,1)
               node2 = npijkmv(i,2)
               x01=xord0(node1)
               y01=yord0(node1)
               z01=zord0(node1)
                 x02=xord0(node2)
                 y02=yord0(node2)
                 z02=zord0(node2)
               zl0=sqrt( (x02-x01)**2 +(y02-y01)**2 +(z02-z01)**2)
               ieqn1 = idbc(node1,2)
               ieqn2 = idbc(node2,2)
               edot = (vel(ieqn2)-vel(ieqn1))/zl0
               vf = etav*edot
*                 ploadv(ieqn1) = ploadv(ieqn1) +vf
*                 ploadv(ieqn2) = ploadv(ieqn2) -vf
               v2 = vel(ieqn2)
               v1 = vel(ieqn1)
*              goto 100
c
               call get_triadv(xord,yord,zord,nnp,node1,node2
     &                        ,triad_ee)
*              if (i .eq. 4) then
*                  write(*,*) ' triad'
*                  do jj=1,3
*                     write(*,84) (triad_ee(jj,kk),kk=1,3)
*                  enddo
*              endif
               fb=0.0
               fg=0.0
               do ii=1,3
                  ieqn1 = idbc(node1,ii)
                  ieqn2 = idbc(node2,ii)
                  fg(0+ii) =vel(ieqn1)
                  fg(3+ii) =vel(ieqn2)
               enddo
               do ii=1,3
                  do kk=1,3
                     fb(0+ii) = fb(0+ii) + triad_ee(kk,ii)*fg(0+kk)
                     fb(3+ii) = fb(3+ii) + triad_ee(kk,ii)*fg(3+kk)
                  enddo
               enddo
*              if (i .eq. 4) then
*                  write(*,*) ' vels'
*                  do jj=1,6
*                     write(*,84) fb(jj),fg(jj),v1,v2,fb(0+1),fb(3+1)
*                  enddo
*              endif
               edot = (fb(4)-fb(1))/zl0
               vf = etav*edot
               fb=0.0
               fg=0.0
               fb(1) = +vf
               fb(4) = -vf
               do ii=1,3
                  do kk=1,3
                     fg(0+ii) = fg(0+ii) + triad_ee(ii,kk)*fb(0+kk)
                     fg(3+ii) = fg(3+ii) + triad_ee(ii,kk)*fb(3+kk)
                  enddo
               enddo
               do ii=1,3
                  ieqn1 = idbc(node1,ii)
                  ieqn2 = idbc(node2,ii)
                  ploadv(ieqn1) = ploadv(ieqn1) + fg(0+ii)
                  ploadv(ieqn2) = ploadv(ieqn2) + fg(3+ii)
               enddo
*              if (i .eq. 4) then
*                  write(*,*) ' force '
*                  do jj=1,6
*                     write(*,84) fb(jj),fg(jj)
*                  enddo
*              endif
c             write(*,82) x01,x02,y01,y02
c             write(*,*) zl0,node1,node2
c
 100          continue
            enddo
c
 82   format(1x,20(g12.6,1x))
 84   format(1x,20(g16.10,1x))
c
      return 
      end
c
c
c     form visco matrix  [C]  
      subroutine visco_mat( stf,nelv,npiv, npjv, xord,yord,zord
     &                    ,idbc,maxnode
     &                    ,iprof,nloc)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npiv(nelv),npjv(nelv)
         integer idbc(maxnode,3)
         integer iprof(neq), nloc(neq)
c
         real*8 xord(nnp), yord(nnp),zord(nnp)
         real*8 stf(maxstiff)
         real*8 wk(1000)
         real*8 ek6x6(6,6),tt6x6(6,6)
c
c        form each element matrix, and assemble
         do n=1,nelv
            write(*,*) 'npi j',n,npiv(n),npjv(n)
         enddo
            write(*,*) 'maxnode 2 ',maxnode
         do 50 n=1,nelv
            i1=npiv(n)
            j1=npjv(n)
                        write(*,*) i1,j1
                dx= xord(j1) - xord(i1) 
                dy= yord(j1) - yord(i1) 
                dz= zord(j1) - zord(i1) 
                xl=sqrt(dx*dx+dy*dy+dz*dz) 
                xll=dx/xl 
                xmm=dy/xl 
                xnn=dz/xl 
c
c            Calculate the stiffness matrix and assemble 
*            emlen = 1.0/xl
*            ek6x6=0.0
*            ek6x6(1,1)= emlen
*            ek6x6(1,4)=-emlen
*            ek6x6(4,4)= emlen
*            ek6x6(4,1)=-emlen
c
*               write(iout,*) 'ek6x6 before'
*            do i=1,6
*               write(iout,82) (ek6x6(i,j),j=1,6)
*            enddo
c
c            use trans3d to transform from local to global stifness
*            call trans3d_6x6(xll,xmm,xnn,ek6x6)
c
*               write(iout,*) 'ek6x6 after '
*            do i=1,6
*               write(iout,82) (ek6x6(i,j),j=1,6)
*            enddo
             call truss_6x6(xl,xll,xmm,xnn,ek6x6)
*               write(iout,*) 'tt6x6 after '
*            do i=1,6
*               write(iout,82) (tt6x6(i,j),j=1,6)
*            enddo
             call assembCOL_TRUSS(stf,ek6x6,idbc,maxnode ,i1,j1,nloc)
*            write(iout,*) ' elem #',n
*            iwd=18
*            call storeBND_f(stf,maxstiff,iprof,nloc,iout,neq,iwd,wk)
c
 50      continue 
*        stop
         write(*,'(a)') '@@ end elems'
c        END loop over all elements
c
 22                   format(1x,6(g13.6))
 23             format(1x,g13.6,1x,i5 )
 82                   format(1x,20(g13.6))
      return
      end
c
c
      subroutine truss_6x6(xl,xll,xmm,xnn,tt6x6)
         implicit real*8 (a-h,o-z)
         real*8 ek6x6(6,6),tt6x6(6,6)
c
         tt6x6(1,1) = xll*xll/xl
         tt6x6(1,2) = xll*xmm/xl
         tt6x6(1,3) = xll*xnn/xl
           tt6x6(2,1) = tt6x6(1,2)
           tt6x6(2,2) = xmm*xmm/xl
           tt6x6(2,3) = xmm*xnn/xl
             tt6x6(3,1) = tt6x6(1,3)
             tt6x6(3,2) = tt6x6(2,3)
             tt6x6(3,3) = xnn*xnn/xl
          do i=1,3
             do j=1,3
                tt6x6(0+i,3+j) =-tt6x6(i,j)
                tt6x6(3+i,0+j) =-tt6x6(i,j)
                tt6x6(3+i,3+j) = tt6x6(i,j)
             enddo
          enddo
 82                   format(1x,20(g13.6))
      return
      end
c
c
      subroutine get_triadv(xord,yord,zord,nnp,node1,node2,triad_ee)
         implicit real*8 (a-h,o-z)
         real*8  xord (nnp),yord (nnp),zord (nnp)  
         real*8  triad_ee(3,3)
c
         dx = xord(node2) - xord(node1)
         dy = yord(node2) - yord(node1)
         dz = zord(node2) - zord(node1)
         xl = sqrt(dx*dx + dy*dy + dz*dz)
            xll1=dx/xl
            xmm1=dy/xl
            xnn1=dz/xl
            dd=sqrt(xll1*xll1 + xmm1*xmm1)
         triad_ee(1,1)=xll1
         triad_ee(2,1)=xmm1
         triad_ee(3,1)=xnn1
         if (dd .lt. 0.001) then
             triad_ee(1,2)=0.0 
             triad_ee(2,2)=1.0 
             triad_ee(3,2)=0.0 
               triad_ee(1,3)=-xnn1 
               triad_ee(2,3)=0.0 
               triad_ee(3,3)=0.0 
         else
             triad_ee(1,2)=-xmm1/dd
             triad_ee(2,2)=+xll1/dd
             triad_ee(3,2)=0.0 
               triad_ee(1,3)=-xll1*xnn1/dd
               triad_ee(2,3)=-xmm1*xnn1/dd
               triad_ee(3,3)= dd 
         endif
c
      return 
      end
c
c
c SUBROUTINE TRANS3D
c      This subroutine makes 3-D coordinate transformations.
c
       subroutine trans3d_6x6 (l,m,n,ek)
         implicit real*8 (a-h,o-z)
          real*8  ek(6,6),rt(3,3),r(3,3),ktemp(6,6)
          real*8 m,n,l
          pi=4.0*atan(1.0)
c
          beta=0.0
          sb=sin(beta*pi/180)
          cb=cos(beta*pi/180)
          d=sqrt(1.0-n**2)
       if (d .lt. 1.0e-3) then
           r(1,1)  =  0.0
           r(1,2)  =  0.0
           r(1,3)  =  n
           r(2,1)  = -n*sb
           r(2,2)  =  cb 
           r(2,3)  =  0.0
           r(3,1)  = -n*cb
           r(3,2)  = -sb    
           r(3,3)  =  0.0
       else
           r(1,1)  =  l
           r(1,2)  =  m
           r(1,3)  =  n
           if (abs(beta) .le. .01) then
              r(2,1)  =  -m/d
              r(2,2)  =  l/d
              r(2,3)  =  0.0 
              r(3,1)  =  -l*n/d
              r(3,2)  =  -m*n/d
              r(3,3)  =  d
           else
              r(2,1)  =  -(m*cb+l*n*sb)/d
              r(2,2)  =  (l*cb-m*n*sb)/d
              r(2,3)  =  d*sb
              r(3,1)  =  (m*sb-l*n*cb)/d
              r(3,2)  =  -(l*sb+m*n*cb)/d
              r(3,3)  =  d*cb
           endif
       endif
c
       do 7 in=1,3
         do 7 jn=1,3
7          rt(jn,in)=r(in,jn)
c
c  take [Rtrans][K][R] using the nature of [R] to speed computation.
c  k is sectioned off into 3x3s then multiplied [rtrans][k][r]
c
       do 22 i=0,1
        do 22 j=0,1
         do 23 k=1,3
          do 23 ii=1,3
           j1=i*3
           j2=j*3
           ktemp(j1+k,j2+ii)=0.0
           do 23 jj=1,3
            ktemp(j1+k,j2+ii)=ktemp(j1+k,j2+ii)+ek(j1+k,j2+jj)*r(jj,ii)
23       continue
         do 24 k=1,3
          do 24 ii=1,3
           ek(j1+k,j2+ii)=0.0
           do 24 jj=1,3
            ek(j1+k,j2+ii)=ek(j1+k,j2+ii)+rt(k,jj)*ktemp(j1+jj,j2+ii)
24       continue
22     continue
       return 
       end
c
c
c 4720
c     NON-liner Arruda-Boyce Tissue 3D EXPlicit incremental analysis 
      subroutine nonabt_vis_exp(stf,geom,wk,pload, disp, fmag,
     &                   dispt , 
     &                   stfeff,mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   iprofv, iprofh, nloc,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt
     &                  ,respceb,isize8
     &         ,isize_elem,stresst,straint,stressct,straini,stressci
     &         ,isize_dof,disptm,disptp,grav,fmag3
     &         ,dloadt,floadt,dut
     &         ,xyzt,xyzi
     &                      )
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=45000,isize_neq=3000)
         parameter( maxdim=2000,maxdof=maxdim*3 )
c
         integer iprofv(neq), iprofh(neq)
     &          ,nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer nloc(neq)
         integer inod(ireadmax,2)
*        integer inod(ireadmax,2),mnod(iconelem)
c
         real*8 tforce1(maxforce,2),tforce2(maxforce,2)
     &         ,tforce3(maxforce,2)
         real*8 velout(ireadmax)
         real*8 stf(maxstiff),mass(maxstiff),wk(neq), geom(maxstiff)
     &          ,pload(neq ),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 stfeff(maxstiff)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(maxxyz)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 grav(isize_dof)
         real*8 dispt(neq)
         real*8 disptm(isize_dof),disptp(isize_dof)
         real*8 dispfult(maxdim,3)
         real*8 dloadt(isize_dof),floadt(isize_dof),dut(isize_dof)
c
         real*8 fmag3(isize_dof),eforce(3,20)
c
c        eigenstuff
         real*8 respceb(isize8)
         real*8 wk0(isize_neq),wk1(isize_neq),wk2(isize_neq)
         real*8 wk_eig(100)
ctag1
c
cccccccccccccccccccccccccccc
c        tissue 
*        parameter( isize_elem=2000)
         real*8 stresst(isize_elem,6,27),straint(isize_elem,6,27)
         real*8 stressct(isize_elem,6,27)
         integer*4 npijkmv(1000,2)
         real*8 ploadv(1000)
cccccccccccccccccccccccccccc
c
!Fangbao          write(ilog,*)'@@  nonabt_vis_ver320  october 2010 '
c
c
                 itmp4=129
              open(unit=itmp4,file=fdn//'stadyn.tm4',form='unformatted')
                 rewind(itmp4)
c
                 iout3=133
                 open(unit=iout3,file=fdn//'simplex.333')
                 rewind(iout3)
c
        call visco_file(npijkmv,nelv,etav)
           write(*,*) '@@ '
        do i=1,nelv
           write(*,*)i, (npijkmv(i,j),j=1,2)
           do k=1,3
           write(*,*) k, (idbc(npijkmv(i,j),k),j=1,2)
           enddo
        enddo
*       stop
c
c
c        algorithm parameters
         beta0=1.0
         gamma0=1.0
         kref=2
c
           do i=1,nnp
              xyzt(          i)=xyz0(          i)
              xyzt(maxnode  +i)=xyz0(maxnode  +i)
              xyzt(maxnode*2+i)=xyz0(maxnode*2+i)
           enddo
c
c
c        Get things ready for time integration loop
         write(*,*)'INPUT:  time inc  |  # of incs  |  print count |',
     &                                        ' snap count '
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,iprcnt,isnap
               if (isnap .eq. 0) isnap=10000
c
         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'
c
         nsnapmax=npt
c
         write(*,*)' '
         write(*,*)'CHOOSE nodal output:'
         nout=0
27       continue
         nout=nout+1
         write(*,*)'TYPE:  node# |  DoF  | rate    <0 0 0 to end>'
         call zzwrt('  -->  ')
         read(ikbd,*) inode,idof,irate
         write(ilog,*) inode,idof,irate,' ::node  dof '
         if (inode .ne. 0) then
             inod(nout,1)=idbc(inode,idof)
             inod(nout,2)=irate
!Fangbao             write(ilog,*) '@@ eqn # ',inod(nout,1)
             goto 27
         endif
         nout=nout-1
         write(*,'(a)')' @@ '
c
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
     &                                           ,' ::eig t c v' 
c
c
c           form lumped mass
            ilump=1
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
!Fangbao            write(ilog,*)'@@Total Mass: ',z1/3
            write(iout,*)'@@Total Mass: ',z1/3
c              write(iout,86) (mass(j),j=1,neq)
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
*                  write(iout,*) i, grav(i)
                   z1=z1+grav(i)
                enddo
                write(iout,*)'Total gravity force: ',z1
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
c
c        only for second load
         if (igravity_on .eq. 1) then
             write(*,*) '@@ get {Grav} history: '
             call getload_hist(tforce2,maxforce,npt,deltat)
         endif
         if (iload_num .eq. 3) then
             open(unit=itmp,file=fdn//'stadyn.ld3')
             rewind(itmp)
             do i=1,neq
                read(itmp,*) fmag3(i)
             enddo
             write(*   ,*) '@@ Reloaded  {P}_3   OK'
!Fangbao             write(ilog,*) '@@ Reloaded  {P}_3   OK'
             close(itmp)
c
c            INPUT load history from a file and interpolate 
             call getload_hist(tforce3,maxforce,npt,deltat)
         endif
         write(*,*)'                    {P}  | {Grav} | 3rd |     '
         write(*,*)'TYPE force scales:   P1  |  P2    |  P3 | P4  '
         call zzwrt('  -->  ')
         read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
!Fangbao         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
c
c
c        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
!Fangbao         write(ilog ,*)'@@  Beginning increment analysis   '
c
c        Initialize  times, disp
c
         write(*  ,*) '@@ INITIAL disps set to zero'
         dt=deltat
         do i= 1, neq  
              disp(i)=0.0
              dispt(i)=0.0
              vel(i)=0.0
              acc(i)=-dmm*vel(i)
              disptm(i) = dispt(i)-dt*vel(i)+0.5*dt*dt*acc(i)
         enddo
         call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))
         force=0.0
         forceout=force
c
         do i=1,nel
              do j=1,6
                 do k=1,27
                    straint(i,j,k)=0.0
                    stresst(i,j,k)=0.0
                 enddo
              enddo
         enddo
c
         do 71 n=1,nout
              iprnode=inod(n,1)
!Fangbao              write(ilog,*)'@@ iprnode ',n,iprnode
              if (iprnode .eq. 0) then
                  velout(n) = 0.0
                  goto 71
              endif
              if (inod(n,2) .eq. 0) then
                  velout(n)=dispt(iprnode)
              endif
 71      continue
         rewind(idyn)
         rewind(isnp)
         itime=0
         time=0
         ziter=0
         k1old=1
         k2old=1
         k3old=1
           tf1=0
           tf2=0
           tf3=0
         force=0
c
         kount = 0
         ksnap=1
         ksnap=0
         inc=1
         inc=iprcnt
         icarriage=0
         call zzwrt(' @@ ')
c
c
c        prepare BIG INCREMENT LOOP
         deltf1=1
 201     continue
c
         call zzcount_2(0,inc,icarriage)
c        BIG TIME LOOP
         do 200 itime= 1,npt
            call zzcount_2(itime,inc,icarriage)
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
            time = real(itime-1)*deltat
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
            tf0=tf1
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            if (igravity_on .eq. 1) then
                call load_interp(tforce2,maxforce,k2old,time,tf2)
            endif
            if (iload_num .eq. 3) then
                call load_interp(tforce3,maxforce,k3old,time,tf3)
            endif
            pnorm=0.0
*                  write(iout,82) ' force' ,time,neq        
            do i= 1, neq  
               pload(i) =  tf1*fmag(i)*scale_p1 + tf2*grav(i)*scale_p2
     &                    +tf3*fmag3(i)*scale_p3
*                  write(iout,82) pload(i)
            enddo 
            call visco_load(npijkmv,nelv,etav
     &                ,npijkm,idbc,maxnode,maxelem
     &                ,xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1)
     &                ,xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1)
     &                ,vel,ploadv
     &                     )
            ieqn1  = idbc(20,2)
            ieqn2  = idbc(32,2)
            ieqn41 = idbc(41,2)
            write(iout3,82) time,pload(ieqn41)
     &                    ,ploadv(ieqn1),ploadv(ieqn2)
     &                    ,vel(ieqn1),vel(ieqn2)         
            do i=1,neq
               pload(i) = pload(i) + ploadv(i)
            enddo
c
            if (itime .eq. 1) goto 290
c
c           forces from body stresses
c           send current stress = disp increment
            call body_stress_abt( stresst,isize_elem,
     &                    dispfult,maxdim,
     &                    iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadt
     &                 )
ctag1
c
c           form load increment
            do i= 1, neq  
c
c??                                           damp=0.0
                 aterm =  (1.0/(dt**2))*(2*dispt(i)-disptm(i))
                 vterm =  (1.0/(2*dt))*disptm(i)
                 dloadt(i) = pload(i) - floadt(i) 
     &                      + mass(i)*aterm + dmm*mass(i)*vterm 
*                coeff = mass(i)/(dt**2) + damp/(2.0*dt)
                 coeff =     1.0/(dt**2) + 1.0*dmm/(2.0*dt)
                 disptp(i) = dloadt(i)/(coeff*mass(i))
            enddo 
*           write(iout,*)' loads  P F dF dis ',time,new
*           do i= 1, neq  
*                  write(iout,82) pload(i),floadt(i),dloadt(i),disptp(i)
*           enddo 
c
c           update vel and acc
            do i= 1, neq  
               vel(i)=(disptp(i)-           disptm(i))/(2*dt)
               acc(i)=(disptp(i)-2*dispt(i)+disptm(i))/(dt**2)
            enddo 
c
c           interchange time subscripts and update
            do i= 1, neq  
                 disptm(i) = dispt(i)
                 dispt(i)  = disptp(i)
            enddo 
            call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))
*               write(iout,*)' nodal disp'
*               do nn=1,nnp
*                  write(iout,82) (dispfult(nn,j), j=1,3),nn
*               enddo
            isizeflag=isize_elem
            nall=-1
            iextra=0
            call update_stress_abt(
     &                   stresst,straint,stressct,isize_elem,
     &                   dispfult,maxdim,
     &                   npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                   prop,nmprop,ippp,nall,iextra,eforce
     &                )
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
                tf1s=tf1*scale_p1
                tf2s=tf2*scale_p2
                tf3s=tf3*scale_p3
                write(idyn,182) time,tf1s,tf2s,tf3s
     &                              ,(velout(n), n=1,nout)
     &                              ,tf1,tf2,tf3
!                call zzflush(idyn)
            endif
            if (isnap .eq. ksnap .OR. itime .eq. 1) then
                 write(isnp) time, neq, nsnapmax 
                 do i= 1, neq  
                    write(isnp) dispt(i)
                 enddo
                 ksnap=0
            endif
c
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
*                   elseif (ieig_typ .eq. 1) then
c                       vector iter
                    elseif (ieig_typ .eq. 1) then
c                       subspace iteration for multiple eigenvalues
c
                        call sub_spc(stfeff,isize8,dloadt,dut, 
     &                           mass,     iprofv,iprofh,nloc,
     &                           respceb,ieig_loc,wk_eig)
c
c                       write(itmp5,85) time,(wk_eig(j),j=1,10)
c    &                                  ,tf1,tf2,tf3
c    &                             ,s_energy,t_energy,p_energy
                        zlam=wk_eig(1)
                        do k=1,ieig_vec
                           zlam=wk_eig(k)
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
                           write(itmp4) zlam, neq, ieig_max 
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
 214        continue
c
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
 184       format(1x,30(i2,1x))
 185       format(1x,i2,1x,': ',30(i2,1x))
c
      return
      end
c
c
c 4720
c     NON-liner Arruda-Boyce Tissue 3D EXPlicit incremental analysis 
      subroutine nonabt_3D_exp(stf,geom,wk,pload, disp, fmag,
     &                   dispt , 
     &                   stfeff,mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   iprofv, iprofh, nloc,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt
     &                  ,respceb,isize8
     &         ,isize_elem,stresst,straint,stressct,straini,stressci
     &         ,isize_dof,disptm,disptp,grav,fmag3
     &         ,dloadt,floadt,dut
     &         ,xyzt,xyzi
     &                      )
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=45000,isize_neq=3000)
         parameter( maxdim=2000,maxdof=maxdim*3 )
c
         integer iprofv(neq), iprofh(neq)
     &          ,nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer nloc(neq)
         integer inod(ireadmax,2)
*        integer inod(ireadmax,2),mnod(iconelem)
c
         real*8 tforce1(maxforce,2),tforce2(maxforce,2)
     &         ,tforce3(maxforce,2)
         real*8 velout(ireadmax)
         real*8 stf(maxstiff),mass(maxstiff),wk(neq), geom(maxstiff)
     &          ,pload(neq ),vel(neq),acc(neq)
     &          ,uup(neq ),vvp(neq),aap(neq)
         real*8 stfeff(maxstiff)
         real*8 xyz0(maxxyz), prop(ippp,10)
         real*8 xyzt(maxxyz)
         real*8 disp(neq  )
         real*8 fmag(neq  )
         real*8 grav(isize_dof)
         real*8 dispt(neq)
         real*8 disptm(isize_dof),disptp(isize_dof)
         real*8 dispfult(maxdim,3)
         real*8 dloadt(isize_dof),floadt(isize_dof),dut(isize_dof)
c
         real*8 fmag3(isize_dof),eforce(3,20)
c
c        eigenstuff
         real*8 respceb(isize8)
         real*8 wk0(isize_neq),wk1(isize_neq),wk2(isize_neq)
         real*8 wk_eig(100)
ctag1
c
cccccccccccccccccccccccccccc
c        tissue 
*        parameter( isize_elem=2000)
         real*8 stresst(isize_elem,6,27),straint(isize_elem,6,27)
         real*8 stressct(isize_elem,6,27)
cccccccccccccccccccccccccccc
*         maxdim=2000
*         maxdim=900
*         maxdof=maxdim*3
c
!Fangbao          write(ilog,*)'@@  nonabt_3D_exp ver310  october 2010 '
c
c
                 itmp4=129
              open(unit=itmp4,file=fdn//'stadyn.tm4',form='unformatted')
                 rewind(itmp4)
c
c
c        algorithm parameters
         beta0=1.0
         gamma0=1.0
         kref=2
c
           do i=1,nnp
              xyzt(          i)=xyz0(          i)
              xyzt(maxnode  +i)=xyz0(maxnode  +i)
              xyzt(maxnode*2+i)=xyz0(maxnode*2+i)
           enddo
c
c
c        Get things ready for time integration loop
         write(*,*)'INPUT:  time inc  |  # of incs  |  print count |',
     &                                        ' snap count '
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,iprcnt,isnap
               if (isnap .eq. 0) isnap=10000
c
         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'
c
         nsnapmax=npt
c
         write(*,*)' '
         write(*,*)'CHOOSE nodal output:'
         nout=0
27       continue
         nout=nout+1
         write(*,*)'TYPE:  node# |  DoF  | rate    <0 0 0 to end>'
         call zzwrt('  -->  ')
         read(ikbd,*) inode,idof,irate
         write(ilog,*) inode,idof,irate,' ::node  dof '
         if (inode .ne. 0) then
             inod(nout,1)=idbc(inode,idof)
             inod(nout,2)=irate
!Fangbao             write(ilog,*) '@@ eqn # ',inod(nout,1)
             goto 27
         endif
         nout=nout-1
         write(*,'(a)')' @@ '
c
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
     &                                           ,' ::eig t c v' 
c
c
c           form lumped mass
            ilump=1
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
!Fangbao            write(ilog,*)'@@Total Mass: ',z1/3
            write(iout,*)'@@Total Mass: ',z1/3
c              write(iout,86) (mass(j),j=1,neq)
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
*                  write(iout,*) i, grav(i)
                   z1=z1+grav(i)
                enddo
                write(iout,*)'Total gravity force: ',z1
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
!Fangbao         write(ilog,*) '@@ Reloaded  {P}   OK'
c
c        INPUT load history from a file and interpolate 
         call getload_hist(tforce1,maxforce,npt,deltat)
c
c        only for second load
         if (igravity_on .eq. 1) then
             write(*,*) '@@ get {Grav} history: '
             call getload_hist(tforce2,maxforce,npt,deltat)
         endif
         if (iload_num .eq. 3) then
             open(unit=itmp,file=fdn//'stadyn.ld3')
             rewind(itmp)
             do i=1,neq
                read(itmp,*) fmag3(i)
             enddo
             write(*   ,*) '@@ Reloaded  {P}_3   OK'
!Fangbao             write(ilog,*) '@@ Reloaded  {P}_3   OK'
             close(itmp)
c
c            INPUT load history from a file and interpolate 
             call getload_hist(tforce3,maxforce,npt,deltat)
         endif
         write(*,*)'                    {P}  | {Grav} | 3rd |     '
         write(*,*)'TYPE force scales:   P1  |  P2    |  P3 | P4  '
         call zzwrt('  -->  ')
         read(ikbd,*) scale_p1,scale_p2,scale_p3,scale_p4
!Fangbao         write(ilog,83) scale_p1,scale_p2,scale_p3,scale_p4,' ::Psc 123'
c
c
c        INITIALize INCREMENTAL LOOP
         write(*    ,*)'@@  Beginning increment analysis   '
!Fangbao         write(ilog ,*)'@@  Beginning increment analysis   '
c
c        Initialize  times, disp
c
         write(*  ,*) '@@ INITIAL disps set to zero'
         dt=deltat
         do i= 1, neq  
              disp(i)=0.0
              dispt(i)=0.0
              vel(i)=0.0
              acc(i)=-dmm*vel(i)
              disptm(i) = dispt(i)-dt*vel(i)+0.5*dt*dt*acc(i)
         enddo
         call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))
         force=0.0
         forceout=force
c
         do i=1,nel
              do j=1,6
                 do k=1,27
                    straint(i,j,k)=0.0
                    stresst(i,j,k)=0.0
                 enddo
              enddo
         enddo
c
         do 71 n=1,nout
              iprnode=inod(n,1)
!Fangbao              write(ilog,*)'@@ iprnode ',n,iprnode
              if (iprnode .eq. 0) then
                  velout(n) = 0.0
                  goto 71
              endif
              if (inod(n,2) .eq. 0) then
                  velout(n)=dispt(iprnode)
              endif
 71      continue
         rewind(idyn)
         rewind(isnp)
         itime=0
         time=0
         ziter=0
         k1old=1
         k2old=1
         k3old=1
           tf1=0
           tf2=0
           tf3=0
         force=0
c
         kount = 0
         ksnap=1
         ksnap=0
         inc=1
         inc=iprcnt
         icarriage=0
         call zzwrt(' @@ ')
c
c
c        prepare BIG INCREMENT LOOP
         deltf1=1
 201     continue
c
         call zzcount_2(0,inc,icarriage)
c        BIG TIME LOOP
         do 200 itime= 1,npt
            call zzcount_2(itime,inc,icarriage)
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
            time = real(itime-1)*deltat
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
            tf0=tf1
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            if (igravity_on .eq. 1) then
                call load_interp(tforce2,maxforce,k2old,time,tf2)
            endif
            if (iload_num .eq. 3) then
                call load_interp(tforce3,maxforce,k3old,time,tf3)
            endif
            pnorm=0.0
*                  write(iout,82) ' force' ,time,neq        
            do i= 1, neq  
               pload(i) =  tf1*fmag(i)*scale_p1 + tf2*grav(i)*scale_p2
     &                    +tf3*fmag3(i)*scale_p3
*                  write(iout,82) pload(i)
            enddo 
c
            if (itime .eq. 1) goto 290
c
c           forces from body stresses
c           send current stress = disp increment
            call body_stress_abt( stresst,isize_elem,
     &                    dispfult,maxdim,
     &                    iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadt
     &                 )
ctag1
c
c           form load increment
            do i= 1, neq  
c
c??                                           damp=0.0
                 aterm =  (1.0/(dt**2))*(2*dispt(i)-disptm(i))
                 vterm =  (1.0/(2*dt))*disptm(i)
                 dloadt(i) = pload(i) - floadt(i) 
     &                      + mass(i)*aterm + dmm*mass(i)*vterm 
*                coeff = mass(i)/(dt**2) + damp/(2.0*dt)
                 coeff =     1.0/(dt**2) + 1.0*dmm/(2.0*dt)
                 disptp(i) = dloadt(i)/(coeff*mass(i))
            enddo 
*           write(iout,*)' loads  P F dF dis ',time,new
*           do i= 1, neq  
*                  write(iout,82) pload(i),floadt(i),dloadt(i),disptp(i)
*           enddo 
c
c           update vel and acc
            do i= 1, neq  
               vel(i)=(disptp(i)-           disptm(i))/(2*dt)
               acc(i)=(disptp(i)-2*dispt(i)+disptm(i))/(dt**2)
            enddo 
c
c           interchange time subscripts and update
            do i= 1, neq  
                 disptm(i) = dispt(i)
                 dispt(i)  = disptp(i)
            enddo 
            call update_geom_3D( dispt, dispfult,maxdim, idbc,maxnode,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      xyzt(1),xyzt(maxnode+1),xyzt(maxnode*2+1))
*               write(iout,*)' nodal disp'
*               do nn=1,nnp
*                  write(iout,82) (dispfult(nn,j), j=1,3),nn
*               enddo
            isizeflag=isize_elem
            nall=-1
            iextra=0
            call update_stress_abt(
     &                   stresst,straint,stressct,isize_elem,
     &                   dispfult,maxdim,
     &                   npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                   prop,nmprop,ippp,nall,iextra,eforce
     &                )
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
                tf1s=tf1*scale_p1
                tf2s=tf2*scale_p2
                tf3s=tf3*scale_p3
                write(idyn,182) time,tf1s,tf2s,tf3s
     &                              ,(velout(n), n=1,nout)
     &                              ,tf1,tf2,tf3
!                call zzflush(idyn)
            endif
            if (isnap .eq. ksnap .OR. itime .eq. 1) then
                 write(isnp) time, neq, nsnapmax 
                 do i= 1, neq  
                    write(isnp) dispt(i)
                 enddo
                 ksnap=0
            endif
c
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
*                   elseif (ieig_typ .eq. 1) then
c                       vector iter
                    elseif (ieig_typ .eq. 1) then
c                       subspace iteration for multiple eigenvalues
c
                        call sub_spc(stfeff,isize8,dloadt,dut, 
     &                           mass,     iprofv,iprofh,nloc,
     &                           respceb,ieig_loc,wk_eig)
c
c                       write(itmp5,85) time,(wk_eig(j),j=1,10)
c    &                                  ,tf1,tf2,tf3
c    &                             ,s_energy,t_energy,p_energy
                        zlam=wk_eig(1)
                        do k=1,ieig_vec
                           zlam=wk_eig(k)
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
                           write(itmp4) zlam, neq, ieig_max 
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
 214        continue
c
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
 184       format(1x,30(i2,1x))
 185       format(1x,i2,1x,': ',30(i2,1x))
c
      return
      end
c
c
c 4700
c     NON-liner RUBber 3D TL incremental analysis 
      subroutine nonabt_3D_TL(stf,geom,wk,pload, disp, fmag,
     &                   dispt , 
     &                   stfeff,mass,vel,acc,uup,vvp,aap,
     &                   idbc,npijkm,maxnode,maxelem,
     &                   iprofv, iprofh, nloc,
     &                   xyz0,maxxyz,
     &                   iglobal, ippp,
     &                   prop, nmprop,nelt
     &                  ,respceb,isize8
     &         ,isize_elem,stresst,straint,stressi,straini,stressci
     &         ,isize_dof,dispi,dui,grav,fmag3
     &         ,dloadi,floadi,floadm
     &         ,xyzt,xyzi
c    &         ,idload,dloadi,floadi,floadm
     &                      )
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
         parameter( ireadmax=101, iconelem=ireadmax*10)
         parameter( isize_t=300,maxforce=5000,isize_neq=3000)
         parameter( maxdim=2000,maxdof=maxdim*3,idload=1000 )
c
         integer iprofv(neq), iprofh(neq)
     &          ,nmprop(maxelem)
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(maxelem)
         integer nloc(neq)
         integer inod(ireadmax,2)
*        integer inod(ireadmax,2),mnod(iconelem)
c
         real*8 tforce1(maxforce,2), tforce2(5000,2),tforce3(5000,2)
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
*        real*8 dispfult(900,3), dispfuli(900,3)
         real*8 dispfult(maxdim,3), dispfuli(maxdim,3)
c
         real*8 dloadi(isize_dof),floadi(isize_dof),floadm(isize_dof)
c        real*8 dloadi(idload   ),floadi(idload   ),floadm(idload   )
         real*8 eforce(3,20)
c
c        eigenstuff
         real*8 respceb(isize8)
         real*8 wk0(isize_neq),wk1(isize_neq),wk2(isize_neq)
         real*8 wk_eig(100)
ctag1
c
cccccccccccccccccccccccccccc
c        tissue   
*        parameter( isize_elem=200)
         real*8 stresst(isize_elem,6,27),straint(isize_elem,6,27)
         real*8 stressi(isize_elem,6,27)
         real*8 straini(isize_elem,6,27)
         real*8 stressci(isize_elem,6,27)
cccccccccccccccccccccccccccc
*        maxdim=2000
*        maxdim=900
*        maxdof=maxdim*3
c
!Fangbao         write(ilog,*)'@@  nonabt_3D_TL ver310  October 2010 '
c
c
*                itmp5=125
*                open(unit=itmp5,file='stadyn.tm5')
*                rewind(itmp5)
                 itmp4=129
              open(unit=itmp4,file=fdn//'stadyn.tm4',form='unformatted')
                 rewind(itmp4)
c
c
c        algorithm parameters
           write(*,*) maxstiff,' ms'
*          stop
         beta0=1.0
         iter_max=60
         dtol=0.0001
         gamma0=1.0
         kref=2
c
           do i=1,nnp
              xyzt(i)=xyz0(i)
              xyzt(maxnode+i)=xyz0(maxnode+i)
              xyzt(maxnode*2+i)=xyz0(maxnode*2+i)
c              write(iout,82) xyzt(i),xyzt(maxnode+i),xyzt(maxnode*2+i)
           enddo
c
c        Get things ready for time integration loop
         write(*,*)'INPUT:  time inc  |  # of incs  |  print count |',
     &                                        ' snap count '
         call zzwrt('  -->  ')
         read(ikbd,*) deltat,npt,iprcnt,isnap
               if (isnap .eq. 0) isnap=10000
c
         write(ilog,'(1x,g13.6,1x,3(i8,1x),a)') deltat,npt,iprcnt,isnap,
     &                                 ' ::dt # pr# sn#'
c
         nsnapmax=npt
         write(*,*)' '
         write(*,*)'    u = u + beta*du   K = Ke + gamma Kg'
         write(*,*)'INPUT:   beta   |  gamma  |  ramp    '
         write(*,*)' suggest:   1        1         5 '
         call zzwrt('  -->  ')
         read(ikbd,*) beta0,gamma0,maxramp0  
         write(ilog,*) beta0,gamma0,maxramp0    ,'  ::b g r '
         write(*,*)' '
         write(*,*)'@@ algor 1=full N-R  2=mod N-R '
         write(*,*)'INPUT:   algor | iter max |  tolerance  (<1.0E-5)'
         call zzwrt('  -->  ')
         read(ikbd   ,*) imodify,iter_max,dtol
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
c
         write(*,*)' '
         write(*,*)'CHOOSE nodal output:'
         nout=0
27       continue
         nout=nout+1
         write(*,*)'TYPE:  node# |  DoF   |  rate   <0 0 0 to end>'
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
     &                                           ,' ::eig t c v' 
c
c
c        form lumped mass
         ilump=1
         write(ilog,*)'@@ FORMing lumped  mass matrix '
!         call zztimer(ilog,'--> FORMmass')
         call formmass(mass(1), npijkm,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      idbc, maxnode,maxelem,iglobal,
     &                      wk(1),
     &                      prop,nmprop,ippp,iprofv,nloc )
            z1=0.0
            do i=1,neq
               z1=z1+mass(i)
            enddo
         write(ilog,*)'@@Total Mass: ',z1/3
         write(iout,*)'@@Total Mass: ',z1/3
c              write(iout,86) (mass(j),j=1,neq)
!         call zztimer(ilog,'<-- FORMmass')
!         call zzflush(ilog)
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
c        INPUT load history from a file and interpolate later 
         call getload_hist(tforce1,maxforce,npt,deltat)
c
c        only for second load
         if (igravity_on .eq. 1) then
             write(*,*) '@@ get {Grav} history: '
             call getload_hist(tforce2,5000,npt,deltat)
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
c            INPUT load history from a file and interpolate  later
             call getload_hist(tforce3,5000,npt,deltat)
         endif
         write(*,*)'                    {P}  | {Grav} | 3rd |     '
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
c
         write(*  ,*) '@@ INITIAL disps set to zero'
         do i= 1, neq  
              disp(i)=0.0
              dispt(i)=0.0
*             delu(i)=0.0
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
         do i=1,nel
            do j=1,6
               do k=1,27
                   straint(i,j,k)=0.0
                   stresst(i,j,k)=0.0
               enddo
            enddo
         enddo
c        all elements initially not yielded
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
         itime=0
         time=0
           ziter=0
           iter=0
         zlam11=0.0
         zlam00=1.0
           k1old=1
           k2old=1
           k3old=1
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
         tf1=0
         tf2=0
         tf3=0
         deltf1=1
 201     continue
c
         write(iout,*) ' 0000   iter max ',iter,iter_max
         call zzcount_2(0,inc,icarriage)
c        BIG TIME LOOP
         do 200 itime= 1,npt
            call zzcount_2(itime,inc,icarriage)
            kount = kount + 1
            ksnap = ksnap + 1
            dt=deltat
            time = real(itime-1)*deltat
            if (itime .eq. 2) then
                maxramp=20
                maxramp=1 
                maxramp=maxramp0
            else
                maxramp=maxramp0
            endif
            write(iout,'(a)')'@@ '
            write(iout,*)' N: ',itime
               if (zlam11 .gt. 0.0) then
                   cc3= dcc/sqrt(zlam11/zlam00)
               else
                   cc3= dcc/sqrt( 0.01 )
               endif
            write(ilog,'(a)')'@@ '
            write(ilog,*    )'@@  cc3: ',cc3
c           write(ilog,*    )'@@  cc3: ',cc3,zlam00,zlam11
c
c           set load increment
c           Fmag says where the load is applied
c           Done this way in case distributed load applied
            tf0=tf1
            tf3=0
            call load_interp(tforce1,maxforce,k1old,time,tf1)
            if (igravity_on .eq. 1) then
                call load_interp(tforce2,5000,k2old,time,tf2)
            endif
            if (iload_num .eq. 3) then
                call load_interp(tforce3,5000,k3old,time,tf3)
            endif
            pnorm=0.0
            do i= 1, neq  
               pload(i) =  tf1*fmag(i)*scale_p1 + tf2*grav(i)*scale_p2
     &                    +tf3*fmag3(i)*scale_p3
               pnorm = pnorm + pload(i)**2
            enddo 
            pnorm = sqrt(pnorm/neq)
c
c
            if (imodify .eq. 2) then
c               modified N-R
*               write(iout,*) ' before iter max ',iter,iter_max
c
            call formstif_abt(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,
     &                      stresst,isize_elem,
     &                      nelt,iglobal,
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
c
            call formgeom_abt(geom(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfult,maxdim,stresst,isize_elem,
     &                      nelt,iglobal,
     &                      prop,nmprop,ippp,iprofv,nloc )
ctag3
                if (itime .eq. 1) goto 290
c
c           Form effective stiffness matrix by adding geometry matrix
!Fangbao            write(ilog,*)'@@ dcc cc3 ',dcc,cc3
            do i= 1, neq  
               iloc=nloc(i)
               do 25 j= 1, iprofv(i)
                      ii=iloc+j-1
                  stfeff(iloc+j-1) = 
     &               (1+zalpha)*(stf(iloc+j-1) +gamma0*geom(iloc+j-1))
  25           continue
               damp = (dmm*mass(i)+dkk*stfeff(iloc)) + cc3
               stfeff(iloc) =  stfeff(iloc)  
     &                     +(1/(zbeta*dt**2))*mass(i)
     &                     +(1+zalpha)*(zgamma/(zbeta*dt))*damp
c    &                     +(1+zalpha)*(zgamma/(zbeta*dt))*dmm*mass(i)
                enddo   
c
c               Decompose effective stiffness matrix
                ierror=0
            call uduCOL_D( stfeff,maxstiff,neq,ierror,iprofv,nloc,z0,z1)
                write(iout,84) '@@ min/max UDU: ',z0,z1,itime,iter
                zrat=z1/z0
!Fangbao                write(ilog,*)'@@ zrat: ',zrat
                if (ierror.eq.0) then
                    write(*,*)'ERROR: zero diagonal term'
                    return
                endif
                write(iout,*) ' iter max ',iter,iter_max
            endif
c           bottom modified N-R
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
*                 dufuli(i,j)=0.0
                  dispfuli(i,j)=dispfult(i,j)
               enddo
            enddo
            do i=1,nnp
               xyzi(i)=   xyzt(i)
               xyzi(maxnode+i)=xyzt(maxnode+i)
               xyzi(maxnode*2+i)=xyzt(maxnode*2+i)
            enddo
            do n=1,nel
               do i=1,6
                  do k=1,27
                     stressi(n,i,k) = stresst(n,i,k)
                     straini(n,i,k) = straint(n,i,k)
                  enddo
               enddo
            enddo
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
            tf_load =        tf1      
c
c               full N-R
c           forces from body stresses
c           send current stress = disp increment
            call body_stress_abt( stressi,isize_elem,
     &                    dispfuli,maxdim,
     &                    iglobal,
     &                    npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    prop,nmprop,ippp,
     &                    floadi
     &                 )
ctag1
c
c           form load increment
            write(iout,*) ' force: ',tf1,tf_load
            pnormi=0.0
            pnorm =0.0
            do i=1,neq  
               iloc=nloc(i)
               stfeff(iloc) = stf(iloc) +gamma0*geom(iloc)
               damp = (dmm*mass(i)+dkk*stfeff(iloc))+ cc3
               aterm = aap(i) + (  1/(zbeta*dt**2))*(dispi(i)-uup(i))
               vterm = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
               dpp = -mass(i)*aterm
     &                 - (1+zalpha)*damp*vterm 
     &                + zalpha*dmm*mass(i)*vel(i)
               ploadm = tf0*fmag(i)*scale_p1
               dloadi(i) = (1+zalpha)*pload(i)  - zalpha*ploadm
     &                     - (1+zalpha)*floadi(i) + zalpha*floadm(i)
     &                     + dpp
               pnormi = pnormi + dloadi(i)**2
               pnorm  = pnorm  + pload (i)**2
*                 write(iout,*) dloadi(i)
            enddo 
            pnormi=sqrt(pnormi/neq)
            pnorm =sqrt(pnorm /neq)
c           write(iout,*) ' force: ',tf1,tf_load
c
c
            if (imodify .eq. 1) then 
c               full N-R
                call formstif_abt(stf(1), npijkm,idbc,maxnode,maxelem,
     &                      xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                      dispfuli,maxdim,
     &                      stresst,isize_elem,
     &                      nelt,iglobal,
     &                      prop,nmprop,ippp,iprofv,nloc,kref)
*              do i=1,neq
*                 write(iout,*) i,dloadi(i),'aft stif'
*              enddo
ctag2
c
c                SECOND assemble the geometric stiffness
                 call formgeom_abt(geom(1), npijkm,idbc,maxnode,maxelem
     &                     ,xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1)
     &                     ,dispfuli,maxdim,stressi,isize_elem
     &                     ,nelt,iglobal
     &                     ,prop,nmprop,ippp,iprofv,nloc )
ctag3
c                Form effective stiff matrix by adding geometry matrix
                 do i=1,neq  
                    iloc=nloc(i)
                    do 24 j=1,iprofv(i)
                       ii=iloc+j-1
                       stfeff(iloc+j-1) = 
     &                 (1+zalpha)*(stf(iloc+j-1) +gamma0*geom(iloc+j-1))
 24                 continue
                    damp = (dmm*mass(i)+dkk*stfeff(iloc)) + cc3
                    damp =  dmm*mass(i)                            
                    stfeff(iloc) =  stfeff(iloc)  
     &                          +(1/(zbeta*dt**2))*mass(i)
     &                          +(1+zalpha)*(zgamma/(zbeta*dt))*damp
                 enddo   
c
c                Decompose effective stiffness matrix
                 ierror=0
                 call uduCOL_D( stfeff,maxstiff,neq,ierror
     &                        ,iprofv,nloc,z0,z1)
                 write(iout,84) '@@ min/max UDU: ',z0,z1,itime,iter
                 if (ierror.eq.0) then
                     write(*,*)'ERROR: zero diagonal term'
                     return
                 endif
            endif
c           bottom full N-R
c
c           solve for Delta u
            ier1=0
            call bakCOL_D(stfeff,maxstiff,dloadi, neq, dui,ierror,
     &                                       iprofv,iprofh,nloc)
            dsumd=0.0
            do i=1,neq
               dsumd=dsumd+ abs(dui(i))
            enddo
!Fangbao            write(ilog,*)'@@  |dui| ',dsumd/neq,iter
            write(iout,*)'@@  |dui| ',dsumd/neq,iter
c
c           increment Ui
*           maxramp=10 
            if (iter .le. maxramp) then
                beta=iter/(maxramp+0.0)
                beta=beta*beta0
                zi=2.0**(iter  )
                z0=2.0**(maxramp  )
                beta=zi/z0
                beta=beta*beta0
            else
                beta=1.0
                beta=beta0
            endif
            do i=1,neq  
               dispi(i)= dispi(i) + beta*dui(i)
            enddo
*           write(iout,*) ' beta: ',beta,1/beta
c
            call update_geom_3D(dispi,dispfuli,maxdim,idbc,maxnode,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                    xyzi(1),xyzi(maxnode+1),xyzi(maxnode*2+1))
c
                isizeflag=isize_elem
                nall=-1
                iextra=0
                call update_stress_abt(
     &                   stressi,straini,stressci,
     &                   isize_elem,
     &                   dispfuli,maxdim,
     &                   npijkm,idbc,maxnode,maxelem,
     &                    xyz0(1),xyz0(maxnode+1),xyz0(maxnode*2+1),
     &                   prop,nmprop,ippp,nall,iextra,eforce
     &                )
ctag4
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
c
c           set values for time t+dt -> t
            do i= 1, neq  
*              delu(i)=dispi(i)-dispt(i)
               dispt(i)=dispi(i)
               vel(i) = vvp(i) + (zgamma/(zbeta*dt))*(dispi(i)-uup(i))
               acc(i) = aap(i) + (1.0/(zbeta*dt*dt))*(dispi(i)-uup(i))
            enddo    
            deltf1=tf1-tf0
            do i=1,nnp  
               do j= 1,3  
                  dispfult(i,j)=dispfuli(i,j)
               enddo
            enddo
            do i=1,nnp
               xyzt(i)          =xyzi(i)
               xyzt(maxnode+i)  =xyzi(maxnode+i)
               xyzt(maxnode*2+i)=xyzi(maxnode*2+i)
            enddo
            do i=1,neq  
               floadm(i) = floadi(i)
            enddo 
c
            do n=1,nel
               do i=1,6
                  do k=1,27
                     stresst(n,i,k) = stressi(n,i,k)
                     straint(n,i,k) = straini(n,i,k)
                  enddo
               enddo
            enddo
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
!                call zzflush(idyn)
             endif
             if (isnap .eq. ksnap .OR. itime .eq. 1) then
                 write(isnp) time, neq, nsnapmax 
                 do i= 1, neq  
                    write(isnp) dispt(i)
                 enddo
                 ksnap=0
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
                        zlam11=wk_eig(1)
                        if (itime .eq. 1) zlam00 = zlam11
                        do k=1,ieig_vec
                           zlam=wk_eig(k)
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
                           write(itmp4) zlam, neq, ieig_max 
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
                    write(imon ,85) time,tf1s,(wk_eig(j),j=1,ieig_vec)
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
c
 200     continue
*210  continue
         write(*,'(a)')' @@'
c        END BIG TIME LOOP
c
 81   format(1x,a,1x,90(1x,g12.6))
 82   format(1x,90(1x,g12.6))
 83   format(1x,4(1x,g12.6),1x,a)
 84   format(1x,a,1x,2(1x,g12.6),1x,2(i5,1x))
 85   format(1x,90(1x,g12.6))
 86   format(1x, 6(1x,g12.6))
 182  format(1x,20(g12.6,1x),/,20(g12.6,1x),/,20(g12.6,1x))
 183  format(1x,i5,4x,3(g12.6,1x))
 184  format(1x,30(i2,1x))
 185  format(1x,i2,1x,': ',30(i2,1x))
c
      return
      end
c
c
c     Nodal loads due to body stresses for 3D solids
      subroutine body_stress_abt( stresst,isize_elem,
     &                   dispful,maxdim,
     &                   iglobal,
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xord0,yord0,zord0,
     &                   prop,nmprop,ippp,
     &                   gforce
     &                )
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
         real*8  stress(6,27), force20(3,20)
         real*8  gforce(3000)
         real*8  uvw(3,20),xyz(3,20)
cccccccccccccccccc tissue 
         real*8  stresst(isize_elem,6,27)
ccccccccccccccc
c
ctag1
c
         do i=1,neq
            gforce(i)=0.0
         enddo
c
c        For each element, calculate the strain, stress at centroid
         do 50 n=1,nel
            mat=nmprop(n)
            neltype=npijkm(n,1)
c
                  e0=prop(mat,1)
                  g0=prop(mat,2)
                  t0=prop(mat,3)
                  r0=prop(mat,4)
                 pl0=prop(mat,5)
c
c
              if (neltype .eq. 5) then
c                 tetra
              elseif (neltype .eq. 9) then
c                 hex8
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
c
                  do i=1,6
                     do k=1,27
                        stress(i,k)=stresst(n,i,k)
                     enddo
                  enddo
c
                  ired_int=ri3
                  call fstress_abt_HEX20(neltype,ired_int, 
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
c     Body FORCE from stress for HEXahedron 20 nodes
      subroutine fstress_abt_HEX20(neltype,ired_int, 
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
         do i=1,kmax*3
            ff(i) = 0.0 
         enddo
         ip_n=0
         nint=ired_int
ctag1   
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
                  do i=1,6
                     ss(i)=stress(i,ip_n)
                  enddo
c
c                 add contrib to force
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
                  do 70 j=1,kmax*3
                     sum = 0.0
                     do k=1,6
                        sum = sum + ss(k)*BE(k,j)
                     enddo
                     ff(j) = ff(j) + sum*wt
 70               continue
 80      continue
c
         nn=0
         do j=1,kmax
            do i=1,3
               nn=nn+1
               force(i,j)=ff(nn)
            enddo
         enddo
c
 82        format(1x,70(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccc
c
c     FORM STIFfness matrix  [K]  
      subroutine formstif_abt( stf,npijkm,idbc,maxnode,maxelem,
     &                       xord0, yord0,zord0,
     &                       dispfult,maxdim,
     &                       stresst,isize_elem,
     &                       nelt,iglobal,
     &                       prop,nmprop,ippp,iprofv,nloc,kref)
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         integer npijkm(maxelem,21),idbc(maxnode,3),nelt(nel)
         integer nmprop(nel)
         real*8  prop(ippp,10),dd(6,6)
         real*8  dispfult(maxdim,3)
         integer iprofv(neq), nloc(neq)
c
         real*8 xord0(nnp), yord0(nnp),zord0(nnp)
         real*8 stf(maxstiff)
         real*8 ek(60,60)
*        real*8 qms(neq), temp
c
         real*8 xyz(3,20),uvw(3,20)    
cccccccccccccccccc plastic
         real*8  stresst(isize_elem,6,27)
         real*8  stress(6,27),pstrain(6,27)
*        integer ipflagi(isize_elem,27),iflags(27)
ccccccccccccccc
c
c
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
c
ctag2
              if (neltype .eq. 5) then
c                 tetra
              elseif (neltype .eq. 9) then
c                 Hex8
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
c
                  do i=1,6 
                     do k=1,27
                        stress(i,k)=stresst(n,i,k)
                     enddo
                  enddo
                  ired_int=ri3
                  call stiff_abt_HEX20(neltype,ired_int,dd, xyz,
     &                          uvw,
     &                          stress,pstrain,ek,
     &                          e0,g0 )
c
               endif
c
               ielm=n
               call assembCOL_3D(stf,ek,idbc,maxnode,npijkm,maxelem,
     &                           ielm,nloc)
c
 50      continue 
c        END loop over all elements
c
 82            format(1x,70(g12.6,1x))
 59            format(1x,18(g12.6))
 86            format(1x, 6(g12.6,1x))
      return
      end
c
c
c
      subroutine abt_mat(ud,e0,g0,ss,dd,iconstit)
         implicit real*8 (a-h,o-z)
         real*8 ud(9),ss(6),def(3,3),cc(3,3),del(3,3),ssk(3,3)
         real*8 cci(3,3)
         real*8 ddk2(3,3,3,3), dd(6,6),ee2(6),ss2(6)
         integer ijkl(6,2)
c
         iout=23
         iconstit=1
         iconstit=2
              iconstit=1
              iconstit=3
              iconstit=2
         if (iconstit .eq. 1) then
c             linearized
              e0=10e2
              g0=3.819e2
              g0=3.3889e2
              e0=1.798e2
              znu0=0.475
              znu0=0.4750016
              znu0=0.475
              zmu0=60
              e0=2*(1+znu0)*zmu0
*             g0=0.5*e0/(1.0+znu0)
              g0=zmu0
              call dmat(e0,g0,dd)
              call elast_lin(ud,dd,ss)
              return
         elseif (iconstit .eq. 3) then
c             linearized
              znu0=0.475
              zmu0=60
              ee0=2*(1+znu0)*zmu0
              gg0=zmu0
              call dmat(ee0,gg0,dd)
              call elast_lin(ud,dd,ss)
              write(*,*)' linear'
              do i=1,6
                 write(*,82) (dd(i,j),j=1,6) 
              enddo
         endif
                c00=e0
                gam=g0
                a10 = c00*(1-gam)
                a01 = c00*(  gam)
                zmu=2*c00
                zk0 = 20*zmu
                zk0 = 19.666667*zmu
         zmu0=e0
         zlamc=g0
         zk0 = 19.666667*zmu0
c        zk0 = 160*zmu0
c
c        assignment to [6x6]
         ijkl(1,1)=1
         ijkl(1,2)=1
              ijkl(2,1)=2
              ijkl(2,2)=2
                ijkl(3,1)=3
                ijkl(3,2)=3
            ijkl(4,1)=1
            ijkl(4,2)=2
              ijkl(5,1)=2
              ijkl(5,2)=3
                ijkl(6,1)=1
                ijkl(6,2)=3
c
         do i=1,3
                      do j=1,3
                         del(i,j)=0.0  
                      enddo
         enddo
         del(1,1)=1.0    
         del(2,2)=1.0    
         del(3,3)=1.0    
ctag1
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
c
c           form [C] and inverse    
            do i=1,3
               do j=1,3
                  sum = 0.0  
                  do k=1,3
                     sum = sum + def(k,i)*def(k,j)
                  enddo
                  cc(i,j)=sum
               enddo
            enddo
            call inv3x3(cc,cci)
c
c           [C] invariants
            zi1 = cc(1,1) + cc(2,2) + cc(3,3)
            zi2 = cc(1,1)*cc(2,2) +cc(2,2)*cc(3,3) +cc(3,3)*cc(1,1)
     &             - cc(1,2)**2      -cc(2,3)**2      -cc(1,3)**2     
            zi3 = cc(1,1)*cc(2,2)*cc(3,3) +2*cc(1,2)*cc(2,3)*cc(1,3)
     &             - cc(1,1)*cc(2,3)**2 -cc(2,2)*cc(1,3)**2
     &             - cc(3,3)*cc(1,2)**2     
            zj3 = sqrt(zi3)
c
            zj1=zi1/zi3**(1.0/3.0)
            a1=1/2.0
            a2=1/20.0
            a3=11/1050.0
            a4=19/7000.
            a5=519/673750.0
            dudj = zmu0*(a1 + 2*a2*zj1/zlamc**2 + (3*a3*zj1**2)/zlamc**4
     &                     + (4*a4*zj1**3)/zlamc**6
     &                     + (5*a4*zj1**4)/zlamc**8  )
            zi3m3 = 1.0/(zi3**(1.0/3.0)) 
            b0 = 2*dudj*zi3m3
            b2 = -(2.0/3)*dudj*zi1*zi3m3
c
c           stress
            do i=1,3
               do j=1,3
                  ssk(i,j) = b0*del(i,j) + b2*cci(i,j)
     &                            + zk0*(zj3-1/zj3)*zj3*cci(i,j)
               enddo
            enddo
            do i=1,6
               ii=ijkl(i,1)
               jj=ijkl(i,2)
                  ss(i)=ssk(ii,jj)
            enddo
c
c        tangent modulus
         d2udj2 = zmu0*(   + 2*a2    /zlamc**2 + (6*a3*zj1   )/zlamc**4
     &                     + (12*a4*zj1**2)/zlamc**6
     &                     + (20*a4*zj1**3)/zlamc**8  )
         b00 = 2*d2udj2*(zi3m3**2)
         b02 =-(2.0/3)*d2udj2*zi1*(zi3m3**2)
     &        -(2.0/3)*dudj*zi3m3
         b22 = (2.0/9)*d2udj2*zi1*zi1*(zi3m3**2)
     &        +(2.0/9)*dudj*zi1*zi3m3
c
            do i=1,3
               do j=1,3
                  do k=1,3
                     do l=1,3
              ddk2(i,j,k,l)
     &             = 2*b00*del(i,j)*del(k,l) + 2*b22*cci(i,j)*cci(k,l)
     &             + 2*b02*( del(i,j)*cci(k,l)+cci(i,j)*del(k,l) )
     &             -   b2 *( cci(i,k)*cci(j,l)+cci(i,l)*cci(j,k) )
     &             + 2*zk0*zj3*zj3*cci(i,j)*cci(k,l)
     &             - zk0*(zj3**2-1)*( cci(i,k)*cci(j,l)
     &                               +cci(i,l)*cci(j,k) )
*    &             + zk0*(2*zj3-1)*zj3*cci(i,j)*cci(k,l)
*    &             - zk0*(  zj3-1)*zj3*( cci(i,k)*cci(j,l)
*    &                                  +cci(i,l)*cci(j,k) )
                    enddo
                 enddo
               enddo
            enddo
c
c           assign to [6x6]
            do i=1,6
               ii=ijkl(i,1)
               jj=ijkl(i,2)
               do j=1,6
                  kk=ijkl(j,1)
                  ll=ijkl(j,2)
                  dd(i,j) = ddk2(ii,jj,kk,ll)
               enddo
            enddo
c
c                 test [Sij] vs {S}=[D]{E}
c                 ee2(1)=0.5*(cc(1,1)-1)
c                 ee2(2)=0.5*(cc(2,2)-1)
c                 ee2(3)=0.5*(cc(3,3)-1)
c                 ee2(4)=cc(1,2)
c                 ee2(5)=cc(1,3)
c                 ee2(6)=cc(2,3)
c                 do i=1,6
c                       sum=0.0
c                       do k=1,6
c                          sum=sum+dd(i,k)*ee2(k)
c                       enddo
c                       ss2(i)=sum
c                 enddo
c
      return
 82       format(1x,20(g12.6,1x))
 89       format(1x,6(4i1,1x))
      return
      end
c
c
c     element STIFFness for Arruda-Boyce Tissue using HEX20 elems
c     integration based on iso quad from Bathe  pp295-297
      subroutine stiff_abt_HEX20(neltype,ired_int,dd, xyz,
     &                    uvw,
     &                    stress,pstrain,ek,
     &                    e0,g0)
         Implicit real*8 (a-h,o-z)
c
         real*8 xyz(3,20),uvw(3,20),ek(60,60)
         real*8 dd(6,6),BE(6,60),Bd(9,60)
         real*8 ui(60),udi(9),ss(6)
ccccccccccccccc
         real*8 stress(6,27),pstrain(6,27),sd(6)
c
         real*8 xg(4,4),wgt(4,4),db(6)
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
            enddo
         enddo
         nn=0
         do i=1,kmax
            do j=1,3
               nn=nn+1
               ui(nn)=uvw(j,i)
            enddo
         enddo
c
ctag2
         ip_n=0
         nint=ired_int
               do 80 lz=1,nint
            do 80 ly=1,nint
         do 80 lx=1,nint
                  ip_n=ip_n+1
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
c
c                 deriv oper and jacobian
                  call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
c
c                 add contrib to stiff
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
c
c                 displacement gradients
                  do i=1,9
                        sumi=0.0
                        do k=1,kmax*3
                           sumi=sumi+Bd(i,k)*ui(k)
                        enddo
                        udi(i)=sumi
                  enddo
c                 iconstit=1
                  call abt_mat(udi,e0,g0,ss,dd,iconstit)
c
c                 [k]=[B][D][B]
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
 60                  continue
 70               continue
 80      continue
c        bottom IPs
         do i=1,kmax*3
              if (ek(i,i) .lt. 0.0) then
                  write(iout,*)' zero k diag: ',ek(i,i),i,ip_n
              endif
         enddo
c
c        impose symmetry
         do j=1,kmax*3
            do i=j,kmax*3
               ek(j,i)=ek(i,j)
            enddo
         enddo
c    
 82      format(1x,40(g12.6,1x))
 84      format(1x,40(i4,1x))
 182     format(1x,40(g12.6,1x))
      return
      end
c
c
c     FORM GEOMetric stiffness matrix  for Arruda-Boyce Tissue
      subroutine formgeom_abt(geo,npijkm,idbc,maxnode,maxelem,
     &                       xord0, yord0,zord0,
     &                       dispfult,maxdim,stresst,isize_elem,
     &                       nelt,
     &                       iglobal,
     &                       prop,nmprop,ippp,iprofv,nloc)
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         integer npijkm(maxelem,21),idbc(maxnode,3)
         integer nelt(nel), iprofv(neq), nloc(neq)
         integer nmprop(nel)
         real*8  prop(ippp,10)
         real*8  dispfult(maxdim,3)
c
         real*8  xord0(nnp), yord0(nnp),zord0(nnp)
*        real*8  xord(nnp), yord(nnp),zord(nnp)
         real*8 ekg(60,60)
         real*8 geo(maxstiff)
         real*8 xyz(3,20),uvw(3,20),stress(6,27)
cccccccccccccccccc tissue
         real*8  stresst(isize_elem,6,27) 
ccccccccccccccc
c
c
         iecho=0
c        initialize [G] to zero
         do i=1,maxstiff
            geo(i)=0.0
         enddo
c
ctag3
c        form each element matrix, and assemble
         do 50 n=1,nel
                 mat=nmprop(1)
                 e0=prop(mat,1)
                 g0=prop(mat,2)
                 t0=prop(mat,3)
                neltype=npijkm(n,1)
c
              if (neltype .eq. 5) then
c                 tetra
              elseif (neltype .eq. 9) then
c                 Hex8 
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
c
                  do i=1,6
                     do k=1,27
                        stress(i,k) = stresst(n,i,k)
                     enddo
                  enddo
                  ired_int=ri3
                  call stiffG_HEX20(neltype,ired_int,xyz,uvw,stress,ekg)
c
              endif
              ielm=n
              call assembCOL_3D(geo,ekg,idbc,maxnode,npijkm,maxelem
     &                         ,ielm,nloc)
 50       continue
 82       format(1x,60(g12.6,1x))
c
      return
      end
c
c
c
c     UPDATE STRESS for Arruda-Boyce Tissue  
      subroutine update_stress_abt( 
     &                   stresst_ip,
     &                   straint_ip,stressc_ip,isize_elem,
     &                   dispful,maxdim, 
     &                   npijkm,idbc,maxnode,maxelem,
     &                   xord0,yord0,zord0,
     &                   prop,nmprop,ippp,nall,iextra,force
     &                )
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
         real*8 stresst_ip(isize_elem,6,27)
         real*8 straint_ip(isize_elem,6,27),stressc_ip(isize_elem,6,27)
ccccccccccccccc
c
         real*8 uu(60),ee(6),ud(9),ss(6)
         real*8 skk(3,3),def(3,3),scc(3,3),cc(3,3)
         real*8 BE(6,60),Bd(9,60)
         real*8 xg(4,4),wgt(4,4)
         real*8 force(3,20),ff(60)
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
c        For each element, calculate the strain, stress
         if (nall .lt. 0) then
               nel1=1
               nel2=nel
         else
               nel1=nall
               nel2=nall
         endif
         do 50 n=nel1,nel2
            mat=nmprop(n)
            neltype=npijkm(n,1)
c
                  e0=prop(mat,1)
                  g0=prop(mat,2)
                  t0=prop(mat,3)
                  r0=prop(mat,4)
                 pl0=prop(mat,5)
c
            if (neltype .eq. 5) then
c                 tetra
            elseif (neltype .eq. 9) then
c                 hex8
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
c
c                 update at each integration point
                  nn=0
                  do i=1,kmax
                     do j=1,3
                        nn=nn+1
                        uu(nn)=uvw(j,i)
                     enddo
                  enddo
c
                  do i=1,kmax*3
                     ff(i) = 0.0 
                  enddo
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
                     call abt_mat(ud,e0,g0,ss,dd,iconstit)
                     do i=1,6
                        stresst_ip(n,i,ip_n)=ss(i)
                     enddo
c
c               extra stuff
                if (iextra .eq. 1) then
c                   add contrib to force
                    wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
                    do 70 j=1,kmax*3
                       sum = 0.0
                       do k=1,6
                          sum = sum + ss(k)*BE(k,j)
                       enddo
                       ff(j) = ff(j) + sum*wt
 70                 continue
c
c                   Lagrange strain
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
 80                continue
c                  bottom loop over IP
                  nn=0
                  do j=1,kmax
                     do i=1,3
                        nn=nn+1
                        force(i,j)=ff(nn)
                     enddo
                  enddo
c
            endif
c           bottom diff elems
 50      continue
c        end of loop over elements
c
 82      format(1x,40(g12.6,1x))
 83      format(1x,6(g12.6,1x))
 86      format(1x,6(g12.6,1x))
 866     format(1x,6(g12.6,1x),a)
 182     format(1x,90(g12.6,1x))
 183     format(1x,i5,1x,90(g12.6,1x))
c
 999  continue
      return
      end
c
c
c
c     NONlinear RUBber POST ANalysis
      subroutine nonabt_postan( disp,dispful,
     &                       idbc,npijkm,maxnode,maxelem, 
     &                       xord0,yord0,zord0,
     &                       strn,strs,forc,
     &                       e_elm,s_elm,
     &                       e_elm_ip,s_elm_ip,ep_elm_ip,wk_elm_n,
     &                       ei_elm_ip,si_elm_ip,epi_elm_ip,
     &                       sc_elm_ip,
     &                       prop,nmprop,ippp,iglobal )
c
          implicit real*8 (a-h,o-z)
             include 'commons.std'
         parameter (ireadmax=101,maxdim=2000)
c
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3 )
         integer mnod(ireadmax,2)
         integer iscan(2000,3)
c
         real*8  disp(neq  ),  dispful(nnp,3)
*        real*8                dispful9(900,3)
         real*8                dispful9(maxdim,3)
         real*8  strs(nnp,6),strn(nnp,6),forc(nel,60)
         real*8  e_elm(nel,6),s_elm(nel,6)
         real*8  e_elm_ip(nel,6,27),s_elm_ip(nel,6,27)
         real*8  ep_elm_ip(nel,6,27)
         real*8  wk_elm_n(nel,20,6)
         real*8  ei_elm_ip(nel,6,27),si_elm_ip(nel,6,27)
     &         ,epi_elm_ip(nel,6,27),sc_elm_ip(nel,6,27)
         real*8  xord0(nnp), yord0(nnp), zord0(nnp)
         real*8  xyz(3,20),uvw(3,20)
         real*8  ssi(27,6),ssn(20,6),se_av(6),average(100)
         real*8  strain_ip(6,27),stress_ip(6,27), force_n(3,20)
         real*8                  stressc_ip(6,27)
         real*8                  strainp_ip(6,27)
         real*8                  eforce(3,20)
         real*8  prop(ippp,10),dd(6,6)
         integer nmprop(nel)
c
         parameter( isize_elem=200)
c
         rewind(idyn)
         ired_int=3
c
                do i=1,nel
                    do j=1,6  
                       do k=1,27
                          s_elm_ip(i,j,k) = 0.0
                          si_elm_ip(i,j,k) = 0.0
                          e_elm_ip(i,j,k) = 0.0
                          ei_elm_ip(i,j,k) = 0.0
*                         ep_elm_ip(i,j,k) = 0.0
*                         epi_elm_ip(i,j,k) = 0.0
                          sc_elm_ip(i,j,k) = 0.0
                       enddo
                    enddo
                enddo
c
             write(*,*)'CHOOSE snapshot:'
             write(*,*)'       0=return   '
             write(*,*)'       #=which snap '
             write(*,*)'      -#=all   snaps '
             call zzwrt('  -->  ')
             read(ikbd,*) iwhich 
!Fangbao             write(ilog,*) iwhich , '   ::which snap '
             if (iwhich .lt. 0) goto 2000
c
c            ONE snap
             if (iwhich .gt. 0) then
c
                 iloop=0
 1               continue
                 write(*,*)'  POST menu: '
                 write(*,*)'           0=return '
*                write(*,*)'               GLOBAL '
*                write(*,*)'          11=Global displacements'
*                write(*,*)'          12=Global loads '
                 write(*,*)'               <<SHAPES>>            '
                 write(*,*)'             141=store displacement data '
                 write(*,*)'               <<Hex20>>          '
                 write(*,*)'         41x=documented displacement '
                 write(*,*)'             411=displacement     '
                 write(*,*)'         42x=documented Lagrange strain '
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
c
c                READ all snaps
 290             continue
                 write(*,*)'@@ READing all snaps'
                 rewind (isnp)
                 read(isnp) psload,neqm
!Fangbao                 write(ilog,*)'@@ psload neqm : ', psload,neqm
                 rewind (isnp)
                 ikont=0
                 rewind (iout)
                    write(*,'(a)')' @@ '
c
c                BIG loop over snaps, must do updating also
                 do 200 iw=1,iwhich
                    ikont=ikont+1
                    read(isnp,end=999) psload,neqm
!Fangbao                    write(ilog,*)'@@ psload  : ', psload,ikont
                    time=psload
                    if (mod(ikont,10) .eq. 0) call zzwrti(ikont)
c
                    do i=1,neqm
                       read(isnp,end=999) disp(i)
                    enddo
                    gnum1=psload
c
c                  reassign displacements to dispful( )
                   do i=1,nnp
                      do j=1,3
                         ieqnum=idbc(i,j)
                         if (ieqnum .eq. 0) then
                             dispful(i,j) = 0.0 
                             dispful9(i,j) = 0.0 
                         else
                             dispful(i,j) = disp(ieqnum)
                             dispful9(i,j) = disp(ieqnum)
                         endif
                      enddo    
                   enddo    
c
                 if (ioutput .ge. 141 .AND. ioutput .le. 149) then
                     goto 200
                 elseif (ioutput .ge. 410 .AND. ioutput .le. 419) then
                     goto 200
                 endif
c
c
c                each element, calculate strain, stress 
c                       3-D solid HEX20
                        ipmax=27
c
                    do i=1,nel
                       do j=1,6  
                          do k=1,ipmax
                             s_elm_ip(i,j,k) = si_elm_ip(i,j,k)
                             e_elm_ip(i,j,k) = ei_elm_ip(i,j,k)
*                            ep_elm_ip(i,j,k) = epi_elm_ip(i,j,k)
                          enddo
                       enddo
                    enddo
                    nall=-1
                    iextra=1
                    call update_stress_abt(
     &                   si_elm_ip,ei_elm_ip,
     &                   sc_elm_ip,nel,
     &                   dispful9,maxdim,
     &                   npijkm,idbc,maxnode,maxelem,
     &                    xord0(1),yord0(1),zord0(1),
     &                   prop,nmprop,ippp,
     &                   nall,iextra,eforce
     &                                   )
c
 200             continue
c                   final reassignment
                    do i=1,nel
                       do j=1,6  
                          do k=1,ipmax
                             e_elm_ip(i,j,k) = ei_elm_ip(i,j,k)
                             s_elm_ip(i,j,k) = si_elm_ip(i,j,k)
                          enddo
                       enddo
                    enddo
                 goto 100
c
c
c                STORE results in <<OUT>> file
 100             continue
*                       ipmax=ired_int**3
!Fangbao                 write(ilog,'(a)')' @@'
                 write(iout,'(a)')' @@'
                 write(*   ,'(a)')' writing output'
                 if (ioutput .eq. 1) then
c
                 elseif (ioutput .eq. 141 .OR. ioutput .eq. 142
     &                                      .OR. ioutput .eq. 143) then
                     rewind(iout)
                     evout=time
                     call shapes_HEX_D(npijkm,maxelem,
     &                         xord0,yord0,zord0,dispful,evout,ioutput)
c
                 elseif (ioutput .eq. 411) then
                     call outdis (dispful,'Single snap')
c
c                Lagrange STRAIN 
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
     &                               ,(e_elm_ip(n,k,1), k=1,6)
                        do j=2,ipmax
                           write(iout,84) j
     &                               ,(e_elm_ip(n,k,j), k=1,6)
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
                              ssi(j,k)=e_elm_ip(n,k,j)
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
                        do k=1,ipmax
                           do j=1,6
                              ssi(k,j)=e_elm_ip(n,j,k)
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
     &                               ,(s_elm_ip(n,k,1), k=1,6)
                        do j=2,ipmax
                           write(iout,84) j
     &                               ,(s_elm_ip(n,k,j), k=1,6)
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
                              ssi(j,k)=s_elm_ip(n,k,j)
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
                        do k=1,ipmax
                           do j=1,6
                              ssi(k,j)=s_elm_ip(n,j,k)
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
     &                               ,(sc_elm_ip(n,k,1), k=1,6)
                        do j=2,ipmax
                           write(iout,84) j
     &                               ,(sc_elm_ip(n,k,j), k=1,6)
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
                              ssi(j,k)=sc_elm_ip(n,k,j)
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
                        do k=1,ipmax
                           do j=1,6
                              ssi(k,j)=sc_elm_ip(n,j,k)
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
                 write(*,*)'              0=return  '
*                write(*,*)'             51=element stress'
                 write(*,*)'                <movie>     '
                 write(*,*)'             111=node displacement'
                 write(*,*)'             141=store displacement data '
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
                 write(*,*)'         450=force  history'
                 write(*,*)'             469=nodal   all    special '
                 call zzwrt('  -->  ')
                 read(ikbd,*) ioutput
                 write(ilog,*) ioutput,' ::ioutput(411=disp 42=e 43=s)'
c
                 if (ioutput .eq. 0) then
                     return
                 elseif (ioutput .eq. 111) then
                     goto 390
c
                 elseif (ioutput .eq. 411) then
c                    nodal disp
                     write(*,*)'INPUT: node # '
                     call zzwrt('  -->  ')
                     read(ikbd,*) inode 
                     write(ilog,*) inode , ' ::node '
c
                 elseif (ioutput .eq. 421 .OR. ioutput .eq. 431
     &              .OR. ioutput .eq. 441                     
     &              .OR. ioutput .eq. 461                     
     &                  ) then
c                    element IP value 
                     write(*,*)'INPUT: IP #  |  elem #   '
                     call zzwrt('  -->  ')
                     read(ikbd,*)  ip_n,ielem 
                     write(ilog,*)  ip_n,ielem , ' ::ip_n,elem '
                            mout=1
                            mnod(mout,1)=ielem
c
                 elseif (ioutput .eq. 422 .OR. ioutput .eq. 432
     &              .OR. ioutput .eq. 442                      
     &              .OR. ioutput .eq. 462                      
     &              .OR. ioutput .eq. 450                      
     &                  ) then
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
     &              .OR. ioutput .eq. 444                      
     &              .OR. ioutput .eq. 464                      
     &              .OR. ioutput .eq. 469                      
     &                  ) then
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
                 rewind (isnp)
                 read(isnp) psload,neqm
                 write(ilog,*)'@@ psload neqm : ', psload,neqm
                 rewind (isnp)
                 ikont=0
                 rewind (iout)
                    write(*,'(a)')' @@ '
c
c                BIG loop over snaps
 300             continue
                 ikont=ikont+1
                 read(isnp,end=999) psload,neqm
                 write(ilog,*)'@@ psload  : ', psload,ikont
                 time=psload
                 if (mod(ikont,10) .eq. 0) call zzwrti(ikont)
c
                 do i=1,neqm
                    read(isnp,end=999) disp(i)
                 enddo
                 gnum1=psload
c
c                reassign displacements to dispful( )
                   do i=1,nnp
                      do j=1,3
                         ieqnum=idbc(i,j)
                         if (ieqnum .eq. 0) then
                             dispful(i,j) = 0.0 
                             dispful9(i,j) = 0.0 
                         else
                             dispful(i,j) = disp(ieqnum)
                             dispful9(i,j) = disp(ieqnum)
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
c                    movie, no need for strains etc
                     goto 359
                 endif
c
c                   zero matrices at this time
                    do j=1,6
                       se_av(j)=0.0
                    enddo
                    do j=1,6*4
                       average(j)=0.0
                    enddo
c
ctag3
c                For each element, get strain, stress etc
                 do 350 m=1,mout
                    i = mnod(m,1)
                    klocn=mnod(m,2)-1
                    neltype=npijkm(i,1)
                    ipmax=27
c
                    do j=1,6  
                       do k=1,ipmax
                          s_elm_ip(i,j,k) = si_elm_ip(i,j,k)
                          e_elm_ip(i,j,k) = ei_elm_ip(i,j,k)
*                         ep_elm_ip(i,j,k) = epi_elm_ip(i,j,k)
                       enddo
                    enddo
c
                    nall=i
                    iextra=1
                    call update_stress_abt(
     &                   si_elm_ip,ei_elm_ip,
     &                   sc_elm_ip,nel,
     &                   dispful9,maxdim,
     &                   npijkm,idbc,maxnode,maxelem,
     &                    xord0(1),yord0(1),zord0(1),
     &                   prop,nmprop,ippp,
     &                   nall,iextra,eforce
     &                                   )
                    do j=1,6  
                       do k=1,ipmax
                          strain_ip(j,k) = ei_elm_ip(i,j,k)
                          stress_ip(j,k) = si_elm_ip(i,j,k)
*                         strainp_ip(j,k) = epi_elm_ip(i,j,k)
                          stressc_ip(j,k) = sc_elm_ip(i,j,k)
                       enddo
                    enddo
c
                    if (ioutput .eq. 421 .OR. ioutput .eq. 422      
     &             .OR. ioutput .eq. 431 .OR. ioutput .eq. 432      
     &             .OR. ioutput .eq. 441 .OR. ioutput .eq. 442      
     &             .OR. ioutput .eq. 461 .OR. ioutput .eq. 462      
     &             .OR. ioutput .eq. 450                            
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
c
                    elseif (ioutput .eq. 469) then
                        ii6=0       !Lagrange strain
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=strain_ip(k,j)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        do j=1,6
                           average(ii6+j)=average(ii6+j)+ssn(klocn,j)
                        enddo
                        ii6=6       !Kirchhoff stress
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=stress_ip(k,j)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        do j=1,6
                           average(ii6+j)=average(ii6+j)+ssn(klocn,j)
                        enddo
                        ii6=12       !Cauchy stress  
                        do j=1,ipmax
                           do k=1,6
                              ssi(j,k)=stressc_ip(k,j)
                           enddo 
                        enddo 
                        call node_HEX20(ssi,ssn)
                        do j=1,6
                           average(ii6+j)=average(ii6+j)+ssn(klocn,j)
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
                 if (ioutput .eq. 469) then
                     do j=1,6*3
                        average(j)=average(j)/attach
                     enddo
                 endif
ctag3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                OUTput section
 359             continue
                        ipmax=ired_int**3
                 if (ioutput .eq. 0) then
                     return
c
                 elseif (ioutput .eq. 111) then
c                    movie
                     call scan_line(npijkm,maxelem,iscan,linemax)
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
c                Lagrangian strain
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
c                Kirchhoff stress
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
c                Cauchy stress
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
c
                 elseif (ioutput .eq. 450) then
                     klocn=mnod(1,2)-1
                     write(idyn,24) time,(eforce(j,klocn), j=1,3)
ctag3
c
                 elseif (ioutput .eq. 469) then
                     write(idyn,82) time,(average(j), j=1,6*3)
c
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
 81            format(1x,30(g12.6,1x))
 82            format(1x,30(g12.6,1x))
 83            format(1x,i5,1x,i5,1x,6(g12.6,1x))
 84            format(1x,5x,1x,i5,1x,6(g12.6,1x))
c
 999  continue
      return
      end
c
