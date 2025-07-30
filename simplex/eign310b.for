c
c     SUBspace iterations eigensolver for VIBrations
      subroutine sub_vibs(stf,mass,load,wk,modemax,maxroot,
     &                   ar,br,eigval,eigvec,dj,r,d,ivib,
     &                   iprof,iprof2,nloc,istiff )
c        Finds first m eigenvalues and vectors using Subspace
c        iteration, based on Bathe pp 685-689.
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         character*12 str12
         character*35 str35
         character*4  str4 
         integer iprof(neq),iprof2(neq),nloc(neq)
         integer icluster(100)
         real*8 stf(maxstiff ),  mass(maxmass ) ,load(neq)
         real*8 wk(neq),r(neq,maxroot),dj(maxroot),d(maxroot)
         real*8 ar(maxroot,maxroot),br(maxroot,maxroot)
         real*8 eigval(maxroot),eigvec(maxroot,maxroot)
c
         real*8 zlam,rtol,rtolj, rt,rtolv, dzz,dif
         real*8 znorm,wnorm,vnorm,enorm
         real*8 sum,sumk,summ, eig,eigup,eignex, ev,sqrm
c
         zlam=0.0
         iconv=0
         rtol=1.0D-8
         rtolj=1.0D-12
         do i=1,maxroot
            d(i)=0.0
         enddo
c
             ibandm=1
c
ctag1
c        establish norm based on diagonal
         sumk=0.0
         summ=0.0
         do 70 i=1,neq
            iloc=nloc(i)
            sumk=sumk+abs(stf(iloc))
            summ=summ+abs(mass(i))
            load(i)=mass(i)/stf(iloc)
 70      continue
         sumk=sumk/neq
         summ=summ/neq
         znorm=sumk/summ
        write(* ,'(1x,a,3(g12.6))')'@@ NORMs s mg : ',sumk,summ,znorm
!Fangbao        write(ilog,'(1x,a,3(g12.6))')'@@ NORMs s mg: ',sumk,summ,znorm
cX
c              NORMALIZE stiff and mass
               do i=1,maxstiff
                  stf(i)=stf(i)/sumk
               enddo
               do i=1,maxmass
                  mass(i)=mass(i)/summ
               enddo
               dzz= (znorm/neq)*1.0e-12
*              dzz= (1.0  /neq)*1.0e3
c              dzz=0
               dzz= (stf(1)/sumk)*1.0e-10
               dzz= 1.0e-10
c
c       make copy of stiffness
*       open(unit=itmp,form='unformatted',status='scratch')
        open(unit=itmp,file=fdn//'stadyn.tmp',form='unformatted')
        rewind itmp
        do i=1,maxstiff
           write(itmp) stf(i)
        enddo
        write(*,*)'@@ made copy of normalized [K]'
!Fangbao        write(ilog,*)'@@ made copy of normalized [K]'
c
c       DECOMPOSE stiffness matrix
        ierror=0
        call  uduCOL_D(stf,maxstiff,neq,ierror,iprof,nloc,z0,z1)
!Fangbao        write(ilog,83)'@@ Diag min max: ',z0,z1,z1/z0,z0/z1
!Fangbao        write(ilog,* )'@@ imult/ierror: ',ierror
*       ierror=1
        if ( abs(z0/z1) .le. 1.0e-12) ierror=0
        if (ierror .eq. 0) then
            write(*,*)' '
            write(*,*)'@@ ZERO diagonal term: try  shift'
!Fangbao            write(ilog,*)'@@ ZERO diagonal term: try shift'
            rewind itmp
            do i=1,maxstiff
               read(itmp) stf(i)
            enddo
ccc         call zerodiag   (stf,          mass,
ccc                          neq,iband,ibandm,zlam,dzz,ilog,ierror2)
            call zerodiagCOL(stf,maxstiff, mass,maxmass,iprof,nloc,
     &                       neq,ibandm,zlam, dzz,ilog,ierror2)
            if (ierror2.eq.0) return
        endif
c
c       CHECK if near RIGID body modes
        do n=1,neq
           iloc=nloc(n)
           wnorm=stf(iloc)
           if (abs(wnorm) .le. 10e-10) then
               write(*,*)' '
!Fangbao               write(ilog,*)'@@ test for rigid modes', wnorm,n
               write(*   ,*)'@@ test for rigid modes', wnorm,n
               write(*,*)' '
               write(*   ,*)'@@ !!! possible RIGID modes !!!'
!Fangbao               write(ilog,*)'@@ !!! possible RIGID modes !!!'
               write(*,*)' '
*              write(*,*)'CHOOSE:  0=return 1=ignore  2=shift'
*              call zzwrt('  -->  ')
*              read(ikbd,*) ignore
               ignore=2
               ignore=1
               write(ilog,*) ignore,'  ::0=ignore'
               if (ignore .eq. 0) then
                   return
               elseif (ignore .eq. 1) then
c                  do nothing  
               elseif (ignore .eq. 2) then
                   dzz=abs(z0)*4
                   dzz=abs(z1)/neq
                   dzz=abs(z0)*40
                   dzz=dzz+1.0e-10
                   rewind itmp
                   do i=1,maxstiff
                      read(itmp) stf(i)
                   enddo
***                call zerodiag(stf,mass,neq,iband,ibandm,
***  >                           zlam,dzz,ilog,ierror2)
                 call zerodiagCOL(stf,maxstiff, mass,maxmass,iprof,nloc,
     &                          neq,ibandm,zlam, dzz,ilog,ierror2)
                   if (ierror2.eq.0) return
               endif
               goto 19
           endif
        enddo
        write(*,*)'@@ rigid body modes checked: OK'
 19     continue
c
c       INITIAL load vector: first column is mass diag
        sum=0.0
        do 20 i=1,neq  
           iloc=nloc(i)
           sum=sum+mass(i)*mass(i)
 20     continue
        sqmas=sqrt(sum)
        do 201 i=1,neq  
           iloc=nloc(i)
           r(i,1)=mass(i)/sqmas
 201    continue
        do 22 i=1,neq
           do 24 j=2,maxroot
              r(i,j)=0.0
 24        continue
 22     continue
c
c       find largest of M/K
        nd=neq/maxroot
        l=neq
c       helps distribution for equal values of M/K
        rt_large=-1.0d12
        do 30 j=2,maxroot
           l=l-nd
           rt=rt_large
           do 32 i=1,l  
              if (load(i) .ge. rt) then   
                  rt=load(i)
                  iloc=i
              endif
 32        continue
           do 34 i=l,neq
              if (load(i) .gt. rt) then 
                  rt=load(i)
                  iloc=i
              endif
 34        continue
           r(iloc,j)=1.0
           load(iloc)=rt_large
 30     continue
c
c       RANDOM vector 
        call ranvec(wk,neq,10)
        do 37 i=1,neq
           r(i,maxroot)=wk(i)
 37     continue
ctag2
*          write(*,*)'Rmatrix '
*          write(*,*)'     '
*          do i=1,neq
*             write(*,182) (r(i,j), j=1,maxroot)
*          enddo
 182       format(1x,30(g10.4,1x))
c
        write(*,*)'     '
c       ITERATIONs begin
        iter=0
 100    continue
        iter=iter+1
        write(*,'(1x,a,i4)') '@@ ITERATION: ',iter
c tag2
c
c       SOLVE {u}k+1   &    simultaneously reduce K 
        do 110 j=1,maxroot
           call  bakCOL_D(stf,maxstiff,r(1,j) ,neq,wk,ierror,
     &                                           iprof,iprof2,nloc)
           do 130 i=j,maxroot
              sum=0.0
c             since [K]{u}k+1 = [r]k then only pre-mult
              do 140 k=1,neq
                 sum=sum+wk(k)*r(k,i)
 140          continue
              ar(i,j)=sum
              ar(j,i)=sum
 130       continue
c
c          store the solution vectors for later use with mass
           do 150 i=1,neq
              r(i,j)=wk(i)
 150       continue
 110    continue
c
c       MULTIPLY  {u}[M]{u} for reduced M
        do 160 j=1,maxroot
*          if (ibandm .eq. 1) then
               do i=1,neq
                  wk(i)= mass(i)*r(i,j)
               enddo
*          else
*              do i=1,neq
*                 wk(i)= 0.0
*              enddo
*              call AxBCOL(mass,maxstiff,r(1,j),wk,
*    &                          neq,iprof,iprof2,nloc)
*          endif
c
           do 180 i=j,maxroot
              sum=0.0
              do 190 k=1,neq
                 sum=sum+r(k,i)*wk(k)
 190          continue
              br(i,j)=sum
              br(j,i)=sum
 180       continue
           if (iconv .gt. 0) goto 160
c
c          store [M]{u} as new R load vectors
           do 200 i=1,neq
              r(i,j)=wk(i)
 200       continue
 160    continue
*       do kk=1,20
*          write(*,*) ar(kk,kk),br(kk,kk)
*        enddo
*        stop
c
c       SOLVE reduced eigenvalue prob by JACOBI rotations
        nsmax=15
        n=maxroot
        rtolj=1e-12
        call jacobi(ar,br,eigval,dj,n,rtolj,nsmax,n,eigvec,n,ilog)
        write(*,'(1x,a,1x,i4)') '@@ # of sweeps = ',nsmax
        call eigsrt(eigval,eigvec,n)
*       call eigsrt_neg(eigval,eigvec,n)
*               write(*,*)'vals '
*               do i=1,maxroot
*                  write(*,'(1x,20(g11.5,1x))') eigval(i)
*               enddo
c
c       calculate R times eigenvectors for new loads
        do  420 i=1,neq
            do 422 j=1,maxroot
               dj(j)=r(i,j)
 422        continue
            do 424 j=1,maxroot
               sum=0.0
               do 430 k=1,maxroot
                  sum=sum+dj(k)*eigvec(k,j)
 430           continue
               r(i,j)=sum
 424         continue
 420    continue
c       enforce random last vector
        call ranvec(wk,neq,iter)
        do 1137 i=1,neq
           r(i,maxroot)=wk(i)
 1137   continue
*          write(*,*)'Rmatrix '
*          write(*,*)'     '
*          do i=1,neq
*             write(*,182) (r(i,j), j=1,maxroot)
*          enddo
c
        if (iconv .gt. 0) goto 500
c
c       check for convergence
        do 380 i=1,modemax
           dif=abs(eigval(i)-d(i))
           rtolv=dif/eigval(i)
           if ( abs(rtolv) .gt. rtol) then
               write(*,'(1x,a,i5)')'@@ TRIGGER rtolv: ',i
               write(*,*)'@@      ',rtolv,' vs ',rtol
               write(*,*)'@@ Eigv_1 =',eigval(1)*sumk/summ
     &                  ,  ' Eigv_2 =',eigval(2)*sumk/summ
!Fangbao               write(ilog,*)'@@ Eigv_1 =',eigval(1)*sumk/summ
!Fangbao     &                  ,  ' Eigv_2 =',eigval(2)*sumk/summ
               maxout=maxroot
               if (maxout .gt. 45) maxout=45
               write(idyn,85) iter*1.0,(eigval(j)*sumk/summ, j=1,maxout)
 85            format(1x,50(g12.6,1x))
!               call zzflush(idyn)
               goto 400
           endif
 380    continue
!Fangbao        write(ilog,*)'@@ NORMAL SUBspace convergence'
        write(*,*)'@@ NORMAL SUBspace convergence'
        write(*,*)'@@ another round'
        iconv=1
        goto 100
c
 400    continue
!Fangbao        write(ilog,*)'@@ # of sweeps = ',nsmax,
!Fangbao     &                 '    TRIGGER rtolv: ',i
c       see if max iterations have occurred
        if (iter .eq. itermax) then
            if (i .eq. 1) then
                write(*,*)'@@ NO convergence at ITERs=',itermax
            else
                write(*,*)'@@ INCOMPLETE convergence at ITERs=',itermax
            endif
            write(*,*)'@@ another round'
            iconv=2
            goto 100
        else
            do 440 i=1,maxroot
               d(i)=eigval(i)
 440        continue
            goto 100
        endif
c       ITERATIONs end
c
 500                   continue
c
c       CHECK quality
        isuspect=0
                   rewind itmp
                   do i=1,maxstiff
                      read(itmp) stf(i)
                   enddo
        write(*   ,*)'@@ checking for SUSPECTs:'
!Fangbao        write(ilog,*)'@@ checking for SUSPECTs:'
!Fangbao        write(ilog,*)'@@ zlam: ',zlam
        do 580 i=1,maxroot
           eig=eigval(i)-zlam
               do k=1,neq
                  wk(k)= 0.0
               enddo
               call AxBCOL(stf,maxstiff,r(1,i),wk,
     &                          neq,iprof,iprof2,nloc)
           vnorm=0.0
           do 590 n=1,neq
              vnorm=vnorm+wk(n)*wk(n)
 590       continue
*          if (ibandm .eq. 1) then
               do k=1,neq
                  load(k)= mass(k)*r(k,i)
               enddo
*          else
*              do k=1,neq
*                 load(k)= 0.0
*              enddo
*              call AxBCOL(mass,maxstiff,r(1,i),load,
*    &                          neq,iprof,iprof2,nloc)
*          endif
c
           wnorm=0.0
           do 600 n=1,neq
              wk(n)=wk(n)-eig*load(n)
              wnorm=wnorm+wk(n)*wk(n)
 600       continue
           vnorm=sqrt(abs(vnorm))
           wnorm=sqrt(abs(wnorm))
           if (vnorm .lt. 10e-12) vnorm=10e-12
           enorm=wnorm/vnorm
           dj(i)=0
           icluster(i)=i
           write(str35,'(1x,a,i5,1x,g13.6,1x,a)')
     &                  '@@ QUALity: ',i,enorm,' '
           call zzwrt(str35)
           if (enorm .gt. .01) then
               isuspect=isuspect+1
               icluster(i)=0
               str12=' <-- suspect'
               dj(i)=100
               write(*,'(a)') ' <-- suspect'
ccc            write(iout,'(a)') ' <-- suspect'
           else
               str12=' '
               write(*,'(a)') ' '
ccc            write(iout,'(a)') ' '
           endif
!Fangbao           write(ilog,*)'@@ NORM: ',i,enorm,str12
 580    continue
!Fangbao           write(ilog,'(a)')' '
c tag3
        write(*,*)' '
!Fangbao        write(ilog,*)'@@ SUSPECT eigenvalues ',isuspect
        write(*   ,*)'@@ SUSPECT eigenvalues ',isuspect
        write(ilog,*)' '
c
c       FLAG spurious eigens and count others
*       ibound=0
*       do i=1,maxroot
*          ic=icluster(i)
*          if (ic .gt. 0) then
*              ibound=ibound+1
*              icluster(i)=ibound
*          endif
*       enddo
c
c       STURM SEQUENCE for total number
*       if (ivib .eq. 2) then
c       Decompose effective stiffness matrix
*       eigup=1.01*eigval(modemax)
c
*       inc=0
*399    inc=inc+1
*       eignex=0.99*eigval(modemax+inc)
*       if ( abs(eigup) .ge. abs(eignex) ) goto 399
*       write(*,*)'@@ SHIFT total: ',eignex
*       maxeigs=icluster(modemax+inc-1)     
c
ccc     call shiftip (stf,         mass,    neq,iband,ibandm,eignex  )
*       call shiftCOL(stf,maxstiff,mass,maxmass,nloc,
*    &                                      neq,ibandm,eignex)
*       ierror=0
cccc    call  udu   (stf,         neq,iband,ier1)
*       call  uduCOL_D(stf,maxstiff,neq,ierror,iprof,nloc,z0,z1)
*       write(ilog,83)'@@ Diag min max: ',z0,z1,z1/z0
*       if (ierror .eq. 0) then
*           write(*,*)'@@ ZERO diagonal term AGAIN: give up'
*           return
*       endif
*       write(iout,'(a)') '  '
ccc     call countd   (stf,              neq,iband,itotal)
*       call countdCOL(stf,maxstiff,nloc,neq,itotal)
*       write(ilog,*)'@@ should be',maxeigs,' actually',itotal
*       write(*   ,*)'@@ should be',maxeigs,' actually',itotal
*       if (itotal .gt. maxeigs) then
*           write(*,*)'@@ MISSing eigenvalues'
*       elseif (itotal .eq. maxeigs) then
*           write(*,*)'@@ OK eigenvalues'
*       elseif (itotal .lt. maxeigs) then
*           write(*,*)'@@ SPURious eigenvalues'
*       endif
*       endif
c
        write(*,*)'@@ EIGENvalues & vectors in  >>StaDyn.SNP>>'
cc       write(ilog,*)'@@ snap isnp ',snap,isnp
*       rewind(isnp)
*       mdx=0
*       do 510 i=1,modemax
*          if (dj(i) .gt. 1.0) goto 510
*          mdx=mdx+1
*510    continue
*       modemax=mdx
c
c       rectify and store
        rewind(isnp)
        sgn=1.0
        do i=1,modemax
           ev=(eigval(i)-zlam)*sumk/summ
           eigval(i)=ev*sgn
           sqrm=sqrt(summ)
           do j=1,neq
              r(j,i) = r(j,i)/sqrm
           enddo
        enddo
c
c        write(isnp) modemax
c        do i=1,modemax
c*          ev=(eigval(i)+zlam)*sumk/summ
c*          sqrm=sqrt(summ)
c*          write(isnp) ev,(r(j,i)/sqrm, j=1,neq)
c           write(isnp) eigval(i),(r(j,i), j=1,neq)
c                  write(ilog,*)'@@ value: ',i,eigval(i)
c        enddo
        do i=1,modemax
           write(isnp) eigval(i),neq,modemax
           do j=1,neq
              write(isnp) r(j,i)
           enddo
!Fangbao           write(ilog,*)'@@ value: ',i,eigval(i)
        enddo
c
c            echo to OUT file
             pi2=2*4.0*atan(1.0)
             write(iout,*)'EIGENvalues: '
                 write(iout,*)'   N    omega [r/s]    freq [Hz] '
                 do i=1,modemax
***                 ev=(eigval(i)+zlam)*sumk/summ
***                 ev=sqrt(abs(ev))
                    str4='real'
                    if (eigval(i) .lt. 0.0) str4='imag'
                    ev=sqrt(abs(eigval(i)))
                    write(iout,82) i,ev,ev/pi2,str4
                 enddo   
c
 81          format(1x,i5,1x,2(g12.6,2x))
 82          format(1x,i5,1x,2(g12.6,2x),1x,a)
 83          format(1x,a,1x,6(g12.6,2x))
c 
c
      close(itmp)
      return
      end
c
c
c     SUB_SPaCe iterations for lowest eigenmodes
      subroutine sub_spc(stfeff,isize8  ,load,mass2,
     &                        mass,          iprofv,iprofh,nloc,
     &                   respce8,ieig_loc,wk_eig)
c
         implicit real*8 (a-h,o-z)
          include 'commons.std'
         integer  iprofv(neq), iprofh(neq)
         integer nloc(neq)
c
         real*8               mass(neq),load(neq),mass2(neq),wk(neq)
         real*8 stfeff(maxstiff)
         real*8 respce8(isize8),wk_eig(100),wk1(neq),wk2(neq)
c
         real*8 rtol,rtoltest,rho,rho0,sum0,sumz,sqrz
c
         iout=23
c
              do i=1,neq
                 mass2(i)=mass(i)
              enddo
              ivib=2
              maxmass=neq
              neqma=1
              ibandx=iband
              ibandm=1
              istiff=1
cc            write(*,*)'INPUT:  # of modes of interest'
cc            call zzwrt(' --> ')
cc            read(ikbd,*) modemax
              modemax=10
              neq23=(2*neq)/3
              if (modemax .gt. neq23) then
                  write(*,*)'@@ MAX possible modes =',neq,' will use it'
                  write(*,*)'@@     will use 2/3 of it'
                  modemax=neq23
              endif
!Fangbao              write(ilog,*) '@@ ',modemax,' ::max modes'
              write(*,*)' '
              i8=modemax+8
              i8=modemax+10
              i2=modemax*2
              if (i8 .ge. neq) i8=neq
              if (i2 .ge. neq) i2=neq
              if (i8 .lt. i2) then
                  maxroot=i8
              else
                  maxroot=i2
              endif
!Fangbao              write(ilog,*)'@@ maxroot = ',maxroot
              write(*,*)' '
c
c              storage pointers
               ir1=1
               ir2=ir1+maxstiff 
               ir3=ir2+maxmass
               ir4=ir3+neq
               ir5=ir4+neq
               ir6=ir5+maxroot*maxroot
               ir7=ir6+maxroot*maxroot
               ir8=ir7+maxroot
               ir9 =ir8+maxroot*maxroot
               ir10=ir9 +maxroot
               ir11=ir10+neq*maxroot
               ir12=ir11+maxroot
               call  ckSpace(ilog,isize8,        'SUBspc',ir12)
c
!               call zztimer(ilog,'--> SUBspc')
               call sub_vib( stfeff(1), mass2(1),         respce8(ir3),
     >                       respce8(ir4),
     >                       modemax,maxroot,
     >                       respce8(ir5), respce8(ir6),
     >                       respce8(ir7), respce8(ir8), respce8(ir9 ),
     >                       respce8(ir10), respce8(ir11), ivib
     >                      ,iprofv,iprofh,nloc,istiff)
!              call zztimer(ilog,'<-- SUBspc')
              do j=1,maxroot
*                write(iout,*) respce8(ir7 - 1+j)
                 wk_eig(j) =   respce8(ir7 - 1+j)
              enddo
*             do j=1,neq
*                wk1(j) =   respce8(ir10 - 1+j)
*                wk2(j) =   respce8(ir10+neq - 1+j)
*             enddo
              ieig_loc = ir10
      return
      end
c
c
c     SUBspace iteration for VIBrations 
c        Finds first m eigenvalues and vectors using Subspace
c        iteration, based on Bathe pp 685-689.
      subroutine sub_vib(stf,mass,load,wk,modemax,maxroot,
     >                   ar,br,eigval,eigvec,dj,r,d,ivib,
     &                   iprof,iprof2,nloc,istiff)
c
          implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         character*25 str25
         character*35 str35
         integer iprof(neq),iprof2(neq),nloc(neq)
         integer icluster(100)
         real*8 stf(maxstiff ),  mass(maxmass ) ,load(neq)
         real*8 wk(neq),r(neq,maxroot),dj(maxroot),d(maxroot)
         real*8 ar(maxroot,maxroot),br(maxroot,maxroot)
         real*8 eigval(maxroot),eigvec(maxroot,maxroot)
c
         real*8 zlam,rtol,rtolj, rt,rtolv, dzz,dif
         real*8 znorm,wnorm,vnorm,enorm
         real*8 sum,sumk,summ, eig,eigup,eignex, ev,sqrm
c
         zlam=0.0
         iconv=0
         rtol=1.0D-8
         rtolj=1.0D-12
         do i=1,maxroot
            d(i)=0.0
         enddo
c
c
c        establish norm based on diagonal
         sumk=0.0
         summ=0.0
         do 70 i=1,neq
            iloc=nloc(i)
            sumk=sumk+abs(stf(iloc))
            ilocm=nloc(i)
            if (ilump .eq. 1) ilocm=i
            summ=summ+abs(mass(ilocm))
            load(i)=mass(ilocm)/stf(iloc)
 70      continue
         sumk=sumk/neq
         summ=summ/neq
         znorm=sumk/summ
       write(*   ,'(1x,a,3(g12.6))')'@@NORMs s mg n: ',sumk,summ,znorm
!Fangbao       write(ilog,'(1x,a,3(g12.6))')'@@NORMs s mg n: ',sumk,summ,znorm
cX
c              NORMALIZE stiff and mass
               do i=1,maxstiff
                  stf(i)=stf(i)/sumk
               enddo
               do i=1,maxmass
                  mass(i)=mass(i)/summ
               enddo
               dzz= (znorm/neq)*1.0e-6
               dzz= (1.0  /neq)*1.0
c
c       make copy of stiffness
*       open(unit=itmp,form='unformatted',status='scratch')
*       open(unit=itmp,file='stadyn.tmp',form='unformatted')
*       rewind itmp
*       do i=1,maxstiff
*          write(itmp) stf(i)
*       enddo
*       write(*,*)'@@made copy of normalized [K]'
c
c       DECOMPOSE stiffness matrix
        ierror=0
cccc    call udu   (stf,neq,iband,ier1)
        call uduCOL_D(stf,maxstiff,neq,ierror,iprof,nloc,z0,z1)
        ierror=1
*        if (ierror .eq. 0) then
*            write(*,*)' '
*            write(*,*)'@@ ZERO diagonal term: try  shift'
*            write(ilog,*)'@@ ZERO diagonal term: try shift'
**           rewind itmp
**           do i=1,maxstiff
**              read(itmp) stf(i)
**           enddo
ccc         call zerodiag   (stf,          mass,
ccc                          neq,iband,ibandm,zlam,dzz,ilog,ierror2)
*            call zerodiagCOL(stf,maxstiff, mass,maxmass,iprof,nloc,
*     &                       neq,ibandm,zlam, dzz,ilog,ierror2)
**           if (ierror2.eq.0) return
*        endif
c
c       CHECK if near RIGID body modes
        do n=1,neq
           iloc=nloc(n)
           wnorm=stf(iloc)
           if (abs(wnorm) .le. 10e-10) then
!Fangbao               write(ilog,*)'@@ test for rigid modes', wnorm,n
               write(*   ,*)'@@ test for rigid modes', wnorm,n
               write(*,*)' '
               write(*   ,*)'@@ !!! possible RIGID modes !!!'
!Fangbao               write(ilog,*)'@@ !!! possible RIGID modes !!!'
               write(*,*)' '
*              write(*,*)'CHOOSE:  0=return 1=ignore  2=shift'
*              call zzwrt('  -->  ')
*              read(ikbd,*) ignore
               ignore=2
               ignore=1
               write(ilog,*) ignore,'  ::0=ignore'
*               if (ignore .eq. 0) then
*                   return
*               elseif (ignore .eq. 2) then
*                   rewind itmp
*                   do i=1,maxstiff
*                      read(itmp) stf(i)
*                   enddo
****                call zerodiag(stf,mass,neq,iband,ibandm,
****  >                           zlam,dzz,ilog,ierror2)
*                call zerodiagCOL(stf,maxstiff, mass,maxmass,iprof,nloc,
*     &                          neq,ibandm,zlam, dzz,ilog,ierror2)
*                   if (ierror2.eq.0) return
*               endif
                goto 19
           endif
        enddo
*       write(*,*)'rigid body modes checked  OK'
 19     continue
c
c       INITIAL load vector: first column is mass diag
        sum=0.0
        do 20 i=1,neq  
           iloc=nloc(i)
           if (ilump .eq. 1) iloc=i
           sum=sum+mass(iloc)*mass(iloc)
 20     continue
        sqmas=sqrt(sum)
        do 201 i=1,neq  
           iloc=nloc(i)
           if (ilump .eq. 1) iloc=i
           r(i,1)=mass(iloc)/sqmas
 201    continue
        do 22 i=1,neq
           do 24 j=2,maxroot
              r(i,j)=0.0
 24        continue
 22     continue
c
c       find largest of M/K
        nd=neq/maxroot
        l=neq
c       helps distribution for equal values of M/K
        rt_large=-1.0d12
        do 30 j=2,maxroot
           l=l-nd
           rt=rt_large
           do 32 i=1,l  
              if (load(i) .ge. rt) then   
                  rt=load(i)
                  iloc=i
              endif
 32        continue
           do 34 i=l,neq
              if (load(i) .gt. rt) then 
                  rt=load(i)
                  iloc=i
              endif
 34        continue
           r(iloc,j)=1.0
           load(iloc)=rt_large
 30     continue
c
c       RANDOM vector 
        call ranvec(wk,neq,10)
        do 37 i=1,neq
           r(i,maxroot)=wk(i)
 37     continue
c
*          write(*,*)'Rmatrix '
*          write(*,*)'     '
*          do i=1,neq
*             write(*,182) (r(i,j), j=1,maxroot)
*          enddo
 182       format(1x,30(g12.6,1x))
c
c       ITERATIONs begin
        iter=0
 100    continue
        iter=iter+1
        write(*,'(1x,a,i4)') '@@ ITERATION: ',iter
c
c       SOLVE {u}k+1   &    simultaneously reduce K 
        do 110 j=1,maxroot
ccc        call bak   (stf,         r(1,j) ,neq,iband,wk,ier1)
           call bakCOL_D(stf,maxstiff,r(1,j) ,neq,wk,ierror,
     &                                           iprof,iprof2,nloc)
           do 130 i=j,maxroot
              sum=0.0
c             since [K]{u}k+1 = [r]k then only pre-mult
              do 140 k=1,neq
                 sum=sum+wk(k)*r(k,i)
 140          continue
              ar(i,j)=sum
              ar(j,i)=sum
 130       continue
c
c          store the solution vectors for later use with mass
           do 150 i=1,neq
              r(i,j)=wk(i)
 150       continue
 110    continue
c
c       MULTIPLY  {u}[M]{u} for reduced M
        do 160 j=1,maxroot
ccc        call bandAxd(mass, r(1,j) ,wk,neq,ibandm)
           if (ibandm .eq. 1) then
               do i=1,neq
                  wk(i)= mass(i)*r(i,j)
               enddo
           else
               do i=1,neq
                  wk(i)= 0.0
               enddo
               call AxBCOL(mass,maxstiff,r(1,j),wk,
     &                          neq,iprof,iprof2,nloc)
           endif
c
           do 180 i=j,maxroot
              sum=0.0
              do 190 k=1,neq
                 sum=sum+r(k,i)*wk(k)
 190          continue
              br(i,j)=sum
              br(j,i)=sum
 180       continue
           if (iconv .gt. 0) goto 160
c
c          store [M]{u} as new R load vectors
           do 200 i=1,neq
              r(i,j)=wk(i)
 200       continue
 160    continue
c
c       SOLVE reduced eigenvalue prob by JACOBI rotations
        nsmax=15
        n=maxroot
        rtolj=1e-12
        call jacobi(ar,br,eigval,dj,n,rtolj,nsmax,n,eigvec,n,ilog)
        write(*,'(1x,a,1x,i4)') '@@ # of sweeps = ',nsmax
        call eigsrt_neg(eigval,eigvec,n)
*               write(*,*)'vals '
*               do i=1,maxroot
*                  write(*,'(1x,20(g11.5,1x))') eigval(i)
*               enddo
c
c       calculate R times eigenvectors for new loads
        do  420 i=1,neq
            do 422 j=1,maxroot
               dj(j)=r(i,j)
 422        continue
            do 424 j=1,maxroot
               sum=0.0
               do 430 k=1,maxroot
                  sum=sum+dj(k)*eigvec(k,j)
 430           continue
               r(i,j)=sum
 424         continue
 420    continue
c       enforce random last vector
        call ranvec(wk,neq,iter)
        do 1137 i=1,neq
           r(i,maxroot)=wk(i)
 1137   continue
*          write(*,*)'Rmatrix '
*          write(*,*)'     '
*          do i=1,neq
*             write(*,182) (r(i,j), j=1,maxroot)
*          enddo
c
        if (iconv .gt. 0) goto 500
c
c       check for convergence
        do 380 i=1,modemax
           dif=abs(eigval(i)-d(i))
           rtolv=dif/eigval(i)
           if ( abs(rtolv) .gt. rtol) then
               write(*,'(1x,a,i5)')'@@ TRIGGER rtolv: ',i
               write(*,*)'@@      ',rtolv,' vs ',rtol
               write(*,*)'@@ Eigv_1 =',eigval(1)*sumk/summ
     &                  ,  ' Eigv_2 =',eigval(2)*sumk/summ
               maxout=modemax
*              if (maxout .gt. 45) maxout=45
*              write(idyn,82)iter*1.0,(eigval(j)*sumk/summ, j=1,maxout)
*              call zzflush(idyn)
               goto 400
           endif
 380    continue
        write(ilog,*)'@@ NORMAL SUBspace convergence'
        write(*,*)'@@ NORMAL SUBspace convergence'
        write(*,*)'@@ another round'
        iconv=1
        goto 100
c
 400    continue
!Fangbao        write(ilog,*)'@@ # of sweeps = ',nsmax,
!Fangbao     >                 '    TRIGGER rtolv: ',i
c       see if max iterations have occurred
        if (iter .eq. itermax) then
            write(*,*)'@@ NO convergence at ITERns=',itermax
            write(*,*)'@@ another round'
            iconv=2
            goto 100
        else
            do 440 i=1,maxroot
               d(i)=eigval(i)
 440        continue
            goto 100
        endif
c       ITERATIONs end
c
 500                   continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
        goto 990
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
 990    continue
c       rectify for stability problems
        sgn=1.0
        if (ivib .eq. 1) sgn=-1.0
        do i=1,modemax
           ev=(eigval(i)+zlam)*sumk/summ
           eigval(i)=ev*sgn
           sqrm=sqrt(summ)
           do j=1,neq
              r(j,i) = r(j,i)/sqrm
           enddo
        enddo
c
cccccccccccccccccccccccccccccccccc
        return
ccccccccccccccccccccccccccccccccccc
c
 81          format(1x,i5,1x,20(g12.6,2x))
 82          format(1x,50(g12.6,1x))
c 
c
**    close(itmp)
      return
      end
c
c
c     EIGENsystem solver using JACobi rotations
      subroutine eigenjac(x,eigv,a,b,d,neqmas,ivib,
     &                 iprof,istiff)
c
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer iprof(neq)
         real*8  a(neq ,neq ),b(neq ,neqmas), eigv(neq ),d(neq )
         real*8  x(neq ,neq )
         character*4 str4
c
         real*8 rtol
         rtol=1.0d-12
c
             iwidth=neq
*        istiff=0
*        if (istiff .eq. 0) then
             rewind(istf)
             read(istf) neqs,ibands
!Fangbao             write(ilog,*)'@@ neq iband from STF: ',neqs,ibands
             write(*   ,*)'@@ neq iband from STF: ',neqs,ibands
             do i=1,neq
                read(istf) (x(i,j), j=1,iband)
*               write(iout,86) i*1.0,(x(i,j), j=1,iband)
             enddo
             write(*,   *)'@@ read banded      <<StaDyn.STF<< '
*        endif
c
c       reassign to full form
        do 2 i=1,neq
           jmax=min(iband,(neq-i+1))
           do 21 j=1,jmax 
              a(i,i+j-1)=x(i,j)
21            a(i+j-1,i)=x(i,j)
 2      continue
c
        if (ivib .eq. 2) then
c           vibrations
            rewind(imss)
            read(imss) neqs,ibands
!Fangbao            write(ilog,*)'@@ neq iband from MSS: ',neqs,ibands
            do i=1,neq
               jmax=iprof(i)
               if (ilump .eq. 1) jmax=1
               read(imss) (x(i,j), j=1,jmax)
            enddo
            ibandm=ibands
c           reassign to [b] matrix
            if (ibandm.eq.1) then
                do 31 i=1,neq
                   b(i,1) = x(i,1)
 31             continue
            else
                do 32 i=1,neq
                   jmax=iprof(i)
                   do j=1,jmax
                      b(i,i-j+1)=x(i,j)
                      b(i-j+1,i)=x(i,j)
                   enddo
 32             continue
            endif
            write(*,   *)'@@ EIGN: reloaded [K] [M]  OK     '
!Fangbao            write(ilog,*)'@@ EIGN: reloaded [K] [M]  OK     '
c
        elseif (ivib .eq. 1) then
c           stability
             iwidth=neq
             rewind(igeo)
             read(igeo) neqg,ibandg
!Fangbao             write(ilog,*)'@@ neq iband from GEO: ',neqg,ibandg
             do i=1,neqg
                jmax=iprof(i)
                read(igeo) (x(i,j), j=1,jmax)
             enddo
c
c           reassign to [b]
            do 4 i=1,neq
               jmax=min(iband,(neq  -i+1))
               jmax=iprof(i)
               do j=1,jmax
                  b(i,i-j+1)=x(i,j)
                  b(i-j+1,i)=x(i,j)
               enddo
 4          continue
!Fangbao            write(ilog,*)'@@ ReLOADed   [K] [G]  OK'
            write(*   ,*)'@@ ReLOADed   [K] [G]  OK'
        endif
c
        ibandm=neqmas
        nsmax=50
        write(iout,*)' [A]'
        do i=1,neq
           write(iout,86) (a(i,j),j=1,neq)
        enddo
        write(iout,*)' [B]'
        do i=1,neq
           write(iout,86) (b(i,j),j=1,1)
        enddo
        call jacobi(a,b,eigv,d,neq,rtol,nsmax,ibandm,x,neqmas,ilog)
!Fangbao        write(ilog,'(1x,a,1x,i4)') '@@ # of sweeps = ',nsmax
        write(*,'(1x,a,1x,i4)') '@@ # of sweeps = ',nsmax
        sgn=-1
        sgn=+1
        if (ivib .eq. 1) sgn=-1
        do i=1,neq
           eigv(i)=eigv(i)*sgn
        enddo
        call eigsrt(eigv,x,neq)
c
c            store to disk
             rewind isnp
             do 84 i=1,neq
                write(isnp) eigv(i),neq
                do j=1,neq
                   write(isnp) x(j,i)
                enddo
 84          continue
c
c            echo to OUT file
             pi2=2*4.0*atan(1.0)
             write(iout,*)'EIGENvalues: '
             if (ivib .eq. 1) then
                 write(iout,*)'   N    lambda [load]   '
                 do i=1,neq
                    ev = eigv(i)
                    write(iout,81) i,ev
                 enddo   
             elseif (ivib .eq. 2) then
                 write(iout,*)'   N    omega [r/s]    freq [Hz] '
                 do i=1,neq
                    str4='real'
                    if (eigv(i) .lt. 0.0) str4='imag'
                    ev=sqrt(abs(eigv(i)))
                    write(iout,82) i,ev,ev/pi2,str4
                 enddo   
             endif
 81          format(1x,i5,1x,2(g12.6,2x))
 82          format(1x,i5,1x,2(g12.6,2x),1x,a)
 86          format(1x,60(g12.6,1x))
c 
c
      return
      end
c
c
c     EIGenvalue SoRT routine
      subroutine eigsrt(eigv,x,neq)
         implicit real*8(a-h,o-z)
         dimension eigv(neq),x(neq,neq)
c
         do 13 i=1,neq-1
            k=i
            p=eigv(i)
c           search for lowest value
            do 11 j=i+1,neq
               if (abs(eigv(j)) .lt. abs(p)) then
                  k=j
                  p=eigv(j)
               endif
 11         continue
c
c           re-arrange vectors
            if (k.ne.i) then
                eigv(k)=eigv(i)
                eigv(i)=p
                do 12 j=1,neq
                   p=x(j,i)
                   x(j,i)=x(j,k)
                   x(j,k)=p
 12             continue
            endif
 13      continue
c
      return
      end
c       
c
c     EIGenvalue SoRT routine with NEGative lowest
      subroutine eigsrt_neg(eigv,x,neq)
         implicit real*8(a-h,o-z)
         dimension eigv(neq),x(neq,neq)
c
         do 13 i=1,neq-1
            k=i
            p=eigv(i)
c           search for lowest value
            do 11 j=i+1,neq
               if (   (eigv(j)) .lt.    (p)) then
                  k=j
                  p=eigv(j)
               endif
 11         continue
c
c           re-arrange vectors
            if (k.ne.i) then
                eigv(k)=eigv(i)
                eigv(i)=p
                do 12 j=1,neq
                   p=x(j,i)
                   x(j,i)=x(j,k)
                   x(j,k)=p
 12             continue
            endif
 13      continue
c
      return
      end
c       
c
c     JACOBI rotations based on Bathe pp 643-645.
      subroutine jacobi(a,b,eigv,d,n,rtol,nsmax,ibandm,x,neqmas,ilog)
         implicit real*8(a-h,o-z)
         real*8   a(n,n),b(n,neqmas),eigv(n),d(n),x(n,n)
         real*8   aiiabs,biiabs
         character*1 ch
c
         iwrite=0
         iwrite=0
         if (iwrite .eq. 1) then
              neq=n
              nqq=min(neq,10)
              write(*,*)'A mat ',nqq,neq
              do i=1,nqq
                 write(*,82)(a(i,j), j=1,nqq)
              enddo
              write(*,*)'B mat',neqmas,neq
              do i=1,nqq
                 write(*,82)(b(i,j), j=1,neqmas)
              enddo
                      stop
          endif
 82      format(1x,40(g12.6,1x))
*                 neq=n
*                 ilog=67
*           write(ilog,*)'AAAA'
*           write(ilog,*) (a(i,i), i=1,neq)
*           write(ilog,*)'BBBB'
*           write(ilog,*) (b(i,i), i=1,neq)
c
c       check for zero diagonal terms
        do 10 i=1,n
           if (ibandm.eq.1) then
              bii=b(i,1)
           else
              bii=b(i,i)
           endif
                  biiabs=abs(bii)
                  aiiabs=abs(a(i,i))
           if (aiiabs .lt. 1.0e-10) then
               a(i,i)=2.0e-10
               aiiabs=2.0e-10
           endif
***        if (a(i,i).gt.0.0 .and. biiabs .gt. 0.0) goto 4
           if (aiiabs .gt. 0.0 .AND. biiabs .gt. 0.0) goto 4
               write(*,2020)
               stop'ERROR: 1'
 4         continue
           d(i)=a(i,i)/bii
           eigv(i)=d(i)
 10     continue
c
c       initialize the modal matrix to the unit matrix
        do 30 i=1,n
        do 20 j=1,n
           x(i,j)=0.0
 20     continue
        x(i,i)=1.0
 30     continue
c
        if (n.eq.1) return
c
c       set sweep counter
        itot=0
        nsweep=0
        nr=n-1
 40     continue
        nsweep=nsweep+1
        ch=' '
c
c       sliding threshold criteria: check off-diags
        ithreshold=5
        if (nsweep .le. ithreshold) then
*           eps=(0.01**nsweep)**2
            epst=(0.01**nsweep)**2
        else
            epst=(rtol)**2
        endif
        iknt=0
c
c       run over all off-diags
        do 210 j=1,nr
ccc        if (mod(j,5) .eq. 0)  call zzwrt('.')
           jj=j+1
           do 210 k=jj,n
c             check that off-diag term exceeds the threshold
              if (ibandm.eq.1) then
c                 eptola=(a(j,k)*a(j,k))/(a(j,j)*a(k,k))
                  eptola1=a(j,k)*a(j,k)
                  eptola2=a(j,j)*a(k,k)
*                 write(*,*) a(j,k),a(j,k), a(j,j),a(k,k)
*                 write(*,*) a(j,k)*a(j,k), a(j,j)*a(k,k)
*                 write(*,*) j,k,eptola
                  if ( abs(eptola1).lt. abs(epst*eptola2)) then
c                      no need to do a rotation
c                    a(j,k)=0.0
                     goto 210
                  endif
                  akk=-b(k,1)*a(j,k)
                  ajj=-b(j,1)*a(j,k)
                  ab=a(j,j)*b(k,1)-a(k,k)*b(j,1)
                  iknt=iknt+1
              else
                  eptola=(a(j,k)*a(j,k))/(a(j,j)*a(k,k))
                  eptolb=(b(j,k)*b(j,k))/(b(j,j)*b(k,k))
                  if ( abs(eptola).lt.epst .AND. 
     &                 abs(eptolb).lt.epst) then
c                      no need to do a rotation
c                    a(j,k)=0.0
c                    b(j,k)=0.0
                     goto 210
                  endif
                  akk=a(k,k)*b(j,k)-b(k,k)*a(j,k)
                  ajj=a(j,j)*b(j,k)-b(j,j)*a(j,k)
                  ab=a(j,j)*b(k,k)-a(k,k)*b(j,j)
                  iknt=iknt+1
             endif
             radicl=(ab*ab+4.0*akk*ajj)/4.0
             if (radicl .lt. 0.0) then
                 write(*,*) akk,ajj,ab
                 stop' error 2'
*                write(*,2020)
*                write(*,*)'@@ problem case 2: rad < 0: rad=',radicl
*                write(*,*)'@@ will ignore '
                 write(ilog,*)'@@ problem case 2: rad < 0: rad=',radicl
                 write(ilog,*)'@@ will ignore '
                 goto 210
*                stop'ERROR 2'
             endif
c
           sqch=sqrt(radicl)
           d1=ab/2.0+sqch
           d2=ab/2.0-sqch
           den=d1
           if (abs(d2).gt.abs(d1)) den=d2
           if (den .eq. 0.0) then
               ca=0.0
               cg=-a(j,k)/a(k,k)
           else
               ca=akk/den
               cg=-ajj/den
           endif
c
c          do generalized rotation
           if (n.eq.2) goto 190
           jp1=j+1
           jm1=j-1
           kp1=k+1
           km1=k-1
           if (jm1.lt.0) goto 130
c          columns
           do 120 i=1,jm1
              aj=a(i,j)
              ak=a(i,k)
              a(i,j)=aj+cg*ak
              a(i,k)=ak+ca*aj
              if (ibandm.ne.1) then
                 bj=b(i,j)
                 bk=b(i,k)
                 b(i,j)=bj+cg*bk
                 b(i,k)=bk+ca*bj
              endif
 120       continue
 130       continue
           if (kp1.gt.n) goto 160
c          rows
           do 150 i=kp1,n
              aj=a(j,i)
              ak=a(k,i)
              a(j,i)=aj+cg*ak
              a(k,i)=ak+ca*aj
              if (ibandm.ne.1) then
                 bj=b(j,i)
                 bk=b(k,i)
                 b(j,i)=bj+cg*bk
                 b(k,i)=bk+ca*bj
              endif
 150       continue
 160       continue
           if (jp1.gt.km1) goto 190
c          mixture
           do 180 i=jp1,km1
              aj=a(j,i)
              ak=a(i,k)
              a(j,i)=aj+cg*ak
              a(i,k)=ak+ca*aj
              if (ibandm.ne.1) then
                 bj=b(j,i)
                 bk=b(i,k)
                 b(j,i)=bj+cg*bk
                 b(i,k)=bk+ca*bj
              endif
 180       continue
 190       continue
c          do diagonal terms
           ak=a(k,k)
           a(k,k)=ak+2.0*ca*a(j,k)+ca*ca*a(j,j)
           a(j,j)=a(j,j)+2.0*cg*a(j,k)+cg*cg*ak
c          force off-diagonal term to exact zero
           a(j,k)=0.0
           if (ibandm.eq.1) then
               bk=b(k,1)
               b(k,1)=bk+ca*ca*b(j,1)
               b(j,1)=b(j,1)+cg*cg*bk
           else
               bk=b(k,k)
               b(k,k)=bk+2.0*ca*b(j,k)+ca*ca*b(j,j)
               b(j,j)=b(j,j)+2.0*cg*b(j,k)+cg*cg*bk
c              force off-diagonal term to exact zero
               b(j,k)=0.0
           endif
c          update eigenvectors
           do 200 i=1,n
              xj=x(i,j)
              xk=x(i,k)
              x(i,j)=xj+cg*xk
              x(i,k)=xk+ca*xj
 200       continue
c
 210     continue
c
c         update after each sweep
          do 220 i=1,n
***          if (abs(a(i,i)).le.1e-12) a(i,i)=1e-12
             if (ibandm.eq.1) then
                bii=b(i,1)
             else
                bii=b(i,i)
             endif
             if (abs(bii).le.1e-20) bii   =1e-20
c            if (bii .le. 0.0) then 
c                write(*,*)' a  b',a(i,i),bii
c                write(*,*)'case 3'
c                write(*,2020)
c                stop'ERROR 3'
c            endif
             eigv(i)=a(i,i)/bii
 220      continue
c
c           check convergence
            do 240 i=1,n
c              test on eigenvalues
cccccccccccccc tol=rtol*d(i)+rtol
               tol=rtol*d(i)
               tol=abs(tol)
               diff=abs(eigv(i)-d(i))
               if (diff .gt. tol) then
                  ch = 'd'
                  goto 280
               endif
 240        continue
c
            eps=rtol**2  
            do 250 j=1,nr
c              test off-diagonals
               jj=j+1
               do 252 k=jj,n
c                 epsa=(a(j,k)*a(j,k))/(a(j,j)*a(k,k))
                  epsa1=a(j,k)*a(j,k)
                  epsa2=a(j,j)*a(k,k)
cc                epsa= a(j,k)*a(j,k)
cc                epsajk= abs(a(j,j)*a(k,k))*eps + eps
                  if (ibandm.eq.1) then
                     epsb=0.0
                  else
                     epsb=(b(j,k)*b(j,k))/(b(j,j)*b(k,k))
cc                   epsb= b(j,k)*b(j,k)
cc                   epsbjk= abs(b(j,j)*b(k,k))*eps + eps
                  endif
*                      epsa=abs(epsa)
*                      epsb=abs(epsb)
                       epsa1=abs(epsa1)
                       epsa2=abs(epsa2)
c                if ((epsa.lt.eps) .AND. (epsb.lt.eps)) goto 252
                 if ((epsa1.lt.eps*epsa2) .AND. (epsb.lt.eps)) goto 252
c                    need more iterating
                     ch='o'
                     goto 280
 252           continue
 250        continue
c
c           scale eigenvectors
            do 270 j=1,n
               if (ibandm .eq. 1) then
                   bbjj=   b(j,1)
               else
                   bbjj=   b(j,j)
               endif
                   bbjj=abs(bbjj)
                   bb=sqrt(bbjj)
               do 272 k=1,n
                  x(k,j)=x(k,j)/bb
 272           continue
 270        continue
c
c           converged  return
            nsmax=nsweep
            return
c
c         update d matrix for another round
 280      continue
          znorm=0.0 
          do 290 i=1,n
             diff = d(i)-eigv(i)
             d(i)=eigv(i)
c            znorm=amax1(znorm,abs(diff))
             znorm= max (znorm,abs(diff))
 290      continue
          itot=itot+iknt
          write(*,*)'@@D_NORM:',znorm,' Rotns='
     &                         ,iknt,itot,' ',ch,' ',nsweep
          if (nsweep.lt.nsmax) goto 40
          return
c
 2020     format( '@@ !!! matrices not positive definite')
c
          return
          end
c
c
c     EIGen OUTputs for JACobi solver
      subroutine eigoutjac(a,  eigv,  dispful,idbc,maxnode ,ivib,
     &                     xord,yord,zord,
     &                     npijkm,maxelem)
c
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer idbc(maxnode,3)
         character*50 str1
c
         real*8    a(neq,neq ), eigv(neq ), dispful(nnp,6)
         real*8    xord(nnp),yord(nnp),zord(nnp)
         integer   npijkm(maxelem,21) 
c
         rewind isnp
         do i=1,neq
            read(isnp) eigv(i),jj
            do j=1,neq
               read(isnp) a(j,i)
            enddo
         enddo
c 
 1      continue
        write(*,*)' MODAL storage:'
        write(*,*)'       0=return      |  0              <0 0 to end>'
        write(*,*)'       # of modes    |  1=modal values'
        write(*,*)'            :        |  2=mode shapes '
        write(*,*)'            :        |  3=modal vectors'
        write(*,*)'                     |'
        write(*,*)'       mode #        |  31=store contour data'
        write(*,*)'       mode #        | 141=store contour data'
         call zzwrt(' CHOOSE --> ')
         read(ikbd,*) neigen,istore
         write(*,*) neigen,istore
             if (neigen .gt. neq) then
                 write(*,*) '@@ !!! number greater than MAX, using',neq
                 neigen=neq
             endif
         write(ilog,*) neigen,istore,'   ::store 1=val 2=shape 3=vec'
         if (istore.eq.0) then
             return
c
         elseif (istore.eq.1) then
             write(iout,'(a)')' @@ '
             write(iout,*)'MODAL VALUES:'
             do nn=1,neigen
                write(iout,81) nn,eigv(nn),sqrt(eigv(nn))
 81             format(1x,i5,1x,20(g12.6,1x))
             enddo
c
         elseif (istore.eq.2) then
             do 20 nn=1,neigen
c               Fill in  the displacement vector
                do 70 i=1, nnp
                   do j=1,3 
                      ieqnum = idbc(i,j)
                      if (ieqnum .gt. 0) then
                          dispful(i,j) = a(ieqnum,nn)
                      else
                          dispful(i,j) = 0.0e0
                      endif
                   enddo
 70             continue 
c
c               Print the displacements.
                if (ivib .eq. 2) then
c                   vibrations
                    freq=abs(eigv(nn))
                    freq=sqrt(freq)
                    write(*,23)'RESONANT FREQUENCY:',freq,'rad/s'
                    write(iout,*)' '
                    write(iout,23)'RESONANT FREQUENCY:',freq,'rad/s'
                    str1= 'Mode shapes   '
                    call outdis( dispful, str1)
                elseif (ivib .eq. 1) then
c                   stability 
                    freq=eigv(nn)
                    write(*,23)'BUCKLING LOAD:',freq,' [force]'
                    write(iout,*)' '
                    write(iout,23)'BUCKLING LOAD:',freq,' [force]'
                    str1= 'Mode shapes   '
                    call outdis( dispful, str1)
                endif
 23             format(1x,a,3x,1g13.6,1x,a)
 20          continue
c
         elseif (istore.eq.3) then
             write(iout,*)' '
             write(iout,*)'MODAL VECTORS:'
             do nn=1,neigen
                write(iout,83) eigv(nn), (a(i,nn), i=1,neq)
             enddo
 83             format(1x,7(g12.6,1x))
c
*        elseif (istore.eq.4) then
*            do nn=1,neq
*84             write(isnp) eigv(nn), (a(i,nn), i=1,neq)
*            enddo
c
         elseif (istore .eq. 31) then
c            Fill in  the displacement vector
             neig=neigen
             do i = 1, nnp 
                do j = 1, 3 
                   ieqnum = idbc(i,j)
                   if (ieqnum .gt. 0) then
                       dispful(i,j) = a(ieqnum,neig)
                   else
                       dispful(i,j) = 0.0e0
                   endif
                enddo
             enddo 
*            write(iout,'(1x,i6,1x,g12.6)') nel*3, eigv(neig)
               ev=eigv(neig)
               call shapes_HEX_D(npijkm,maxelem,
     &                       xord,yord,zord,dispful,ev,istore)
*            do n=1,nel
*               x1=xord(npi(n))
*               x2=xord(npj(n))
*               x3=xord(npk(n))
*                  y1=yord(npi(n))
*                  y2=yord(npj(n))
*                  y3=yord(npk(n))
*               z1=zord(npi(n))
*               z2=zord(npj(n))
*               z3=zord(npk(n))
*                     write(iout,82) x1,y1,z1,(dispful(npi(n),j),j=1,6)
*                     write(iout,82) x2,y2,z2,(dispful(npj(n),j),j=1,6)
*                     write(iout,82) x3,y3,z3,(dispful(npk(n),j),j=1,6)
*            enddo 
c
         elseif (istore .eq. 141) then
c            Fill in  the displacement vector
             neig=neigen
             do i = 1, nnp 
                do j = 1, 3 
                   ieqnum = idbc(i,j)
                   if (ieqnum .gt. 0) then
                       dispful(i,j) = a(ieqnum,neig)
                   else
                       dispful(i,j) = 0.0e0
                   endif
                enddo
             enddo 
               rewind(iout)
               ev=eigv(neig)
               call shapes_HEX_D(npijkm,maxelem,
     &                       xord,yord,zord,dispful,ev,istore)
         endif
         goto 1
c
 82      format(1x,9(g12.6,1x))
c
      return
      end
c
c
c     EIGen OUTputs for SUBspace iteration scheme
      subroutine eigoutsub(a,   dispful,idbc,maxnode ,ivib,
     &                     xord,yord,zord,
     &                     npijkm,maxelem)
c
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer idbc(maxnode,3)
         character*50 str1
c
         real*8    dispful(nnp,6)
         real*8    xord(nnp),yord(nnp),zord(nnp)
         integer   npijkm(maxelem,21)
         real*8    a(neq ), eigv
c
 100      continue
          write(ilog,*)'@@ rewinding eigenvectors   <<StaDyn.SNP<<'
          write(*   ,*)'@@ rewinding eigenvectors   <<StaDyn.SNP<<'
          rewind isnp
          read(isnp) z1,j1,modemax
          write(ilog,*)'@@ Max modes: ',modemax
          write(*   ,*)'@@ Max modes: ',modemax
          write(*,*)' '
          rewind isnp
c 
 82      continue
         write(*,*)' MODAL storage:'
         write(*,*)'        0=return     |   0           <0 0 to end>'
         write(*,*)'        # of modes   |   1=modal values'
         write(*,*)'             ..      |   2=mode shapes '
         write(*,*)'             ..      |   3=modal vectors'
         write(*,*)'        mode #       |   6=mode shape   '
         write(*,*)'        mode #       |  31=store contour data'
         write(*,*)'        mode #       | 141=store disp file   '
         call zzwrt(' CHOOSE --> ')
         read(ikbd,*) neigen,istore
         write(*,*) neigen,istore
         if (neigen .gt. modemax) then
             neigen=modemax
         elseif (neigen .lt. 1) then
             neigen=0
         endif
             write(ilog,*) neigen,istore,'     ::# of modes'
             write(*   ,*) neigen,'     ::# of modes '
c
c
         if (istore.eq.0) then
             return
         elseif (istore .eq. 1) then
             write(iout,'(a)')' MODAL VALUES:'
             do nn=1,neigen
                read(isnp) eigv
                do i=1,neq
                   read(isnp) zjunk
                enddo
                write(iout,83) eigv, sqrt(abs(eigv))
 83                   format(1x,7(g13.6,1x))
             enddo
         elseif (istore.eq.2) then
             do 20 nn=1,neigen
c               read(isnp) eigv, (a(i), i=1,neq)
                read(isnp) eigv
                do i=1,neq
                   read(isnp) a(i) 
                enddo
c               Fill in  the displacement vector
                do 70 i=1, nnp*6 
                   node=(i+5)/6
                   jdof=i-(node-1)*6
*                  ieqnum = jbc(i)
                   if (ieqnum .gt. 0) then
                       dispful(node,jdof) = a(ieqnum)
                   else
                       dispful(node,jdof) = 0.0e0
                   endif
 70             continue 
c
c               Print the displacements.
                if (ivib .eq. 2) then
                    freq=abs(eigv)
                    freq=sqrt(freq)
                    write(*,23)'RESONANT FREQUENCY:',freq,'rad/s'
                    write(iout,*)' '
                    write(iout,23)'RESONANT FREQUENCY:',freq,'rad/s'
                    str1= 'Mode shapes   '
                    call outdis( dispful, str1)
                elseif (ivib .eq. 1) then
                    freq=eigv
                    write(*,23)'BUCKLING LOAD:',freq,' [force]'
                    write(iout,*)' '
                    write(iout,23)'BUCKLING LOAD:',freq,' [force]'
                    str1= 'Mode shapes   '
                    call outdis( dispful, str1)
                endif
 23             format(1x,a,3x,1g13.6,1x,a)
 20          continue
c
         elseif (istore.eq.3) then
             write(iout,*)' '
             write(iout,*)'MODAL VECTORS:'
             do nn=1,neigen
c               read(isnp) eigv, (a(i), i=1,neq)
                read(isnp) eigv
                do i=1,neq
                   read(isnp) a(i) 
                enddo
                write(iout,81) eigv, (a(i), i=1,neq)
             enddo
c
         elseif (istore .eq. 6) then
             do nn=1,neigen
c               read(isnp) eigv, (a(i), i=1,neq)
                read(isnp) eigv
                do i=1,neq
                   read(isnp) a(i) 
                enddo
             enddo
c               Fill in  the displacement vector
                do i=1, nnp 
                   do j = 1, 3 
                      ieqnum = idbc(i,j)
                      if (ieqnum .gt. 0) then
                          dispful(i,j ) = a(ieqnum)
                      else
                          dispful(i,j ) = 0.0e0
                      endif
                   enddo 
                enddo
c
c               Print the displacements.
                freq=abs(eigv)
                freq=sqrt(freq)
                write(str1,63)' mode:',neigen,eigv,freq
                call outdis( dispful, str1)
 63             format(1x,a,1x,i3,1x,2(g12.6,1x))
c
         elseif (istore .eq. 31) then
c            find that eigenvector
             do n=1,neigen
c               read(isnp) eigv, (a(i), i=1,neq)
                read(isnp) eigv
                do i=1,neq
                   read(isnp) a(i) 
                enddo
             enddo
c
c               Fill in  the displacement vector
                do i = 1, nnp 
*                  node=(i+5)/6
*                  jdof=i-(node-1)*6
*                  ieqnum = jbc(i)
                do j = 1, 3 
                   ieqnum = idbc(i,j)
                   if (ieqnum .gt. 0) then
                       dispful(i,j ) = a(ieqnum)
                   else
                       dispful(i,j ) = 0.0e0
                   endif
                enddo 
                enddo 
*               write(iout,'(1x,i6,1x,g12.6)') nel*3, eigv
               call shapes_HEX_D(npijkm,maxelem,
     &                       xord,yord,zord,dispful,eigv,istore)
*               do n=1,nel
*                  x1=xord(npi(n))
*                  x2=xord(npj(n))
*                  x3=xord(npk(n))
*                   y1=yord(npi(n))
*                   y2=yord(npj(n))
*                   y3=yord(npk(n))
*                  z1=zord(npi(n))
*                  z2=zord(npj(n))
*                  z3=zord(npk(n))
*                     write(iout,81) x1,y1,z1,(dispful(npi(n),j),j=1,6)
*                     write(iout,81) x2,y2,z2,(dispful(npj(n),j),j=1,6)
*                     write(iout,81) x3,y3,z3,(dispful(npk(n),j),j=1,6)
*               enddo 
         elseif (istore .eq. 141) then
c            store eigenvector as disp
             do n=1,neigen
c               read(isnp) eigv, (a(i), i=1,neq)
                read(isnp) eigv
                do i=1,neq
                   read(isnp) a(i) 
                enddo
             enddo
*            rewind(idis)
*            do i=1,neq
*               write(idis) a(i)
*            enddo
*            write(ilog,*) '@@ idis file written ' 
*            write(*   ,*) '@@ idis file written ' 
             do i = 1, nnp 
                do j = 1, 3 
                   ieqnum = idbc(i,j)
                   if (ieqnum .gt. 0) then
                       dispful(i,j) = a(ieqnum)
                   else
                       dispful(i,j) = 0.0e0
                   endif
                enddo
             enddo 
               rewind(iout)
               call shapes_HEX_D(npijkm,maxelem,
     &                       xord,yord,zord,dispful,eigv,istore)
c
         endif
         goto 100
c
 81      format(1x,17(g13.6,1x))
c
      return
      end
c
c
c     SHIFT COLumn stiffness matrix in place by zlam*mass
      subroutine shiftCOL(stf,maxstiff,mass,maxmass,nloc,
     &                                      neq,ibandm,zlam)
         implicit real*8 (a-h,o-z)
         integer nloc(neq)
         real*8 stf(maxstiff ),  mass(maxmass) 
c
c        Shift matrix by eigenvalue estimate
         if (ibandm .eq. 1) then
            do i=1,neq
               iloc=nloc(i)
               stf(iloc) = stf(iloc) + zlam*mass(i)
            enddo
         else
            do i=1,maxstiff
               stf(i) = stf(i) + zlam*mass(i)
            enddo
         endif
c
      return
      end
c
c
c     handle ZERO DIAGonal term for COLumn storage
      subroutine zerodiagCOL(stf,maxstiff, mass,maxmass,iprof,nloc,
     &                           neq,ibandm,zlam,
     >                    dzz,ilog,ierror2)
         implicit real*8 (a-h,o-z)
         real*8 stf(maxstiff ),  mass(maxmass ) 
         integer*4 iprof(neq),nloc(neq)
c
*        write(*,*)'INPUT:  shift     [recommend > ',dzz,']'
         write(ilog,*)'@@shift  [recommend > ',dzz,']'
*        call zzwrt('  -->  ')
*        read(ikbd,*) dzz
*        write(ilog,*) dzz ,'  ::shift '
         zlam=zlam+dzz
ccc      call shiftip ( stf, mass,  neq,iband,ibandm,zlam)
         call shiftCOL(stf,maxstiff,mass,maxmass,nloc,
     &                                      neq,ibandm,zlam)
c        Decompose effective stiffness matrix
         ierror2=0
ccc      call  udu   (stf,neq,iband,ier2)
         call  uduCOL_D(stf,maxstiff,neq,ierror2,iprof,nloc,z0,z1)
!Fangbao         write(ilog,83)'@@ Diag min max 2: ',z0,z1,z0/z1
!Fangbao         write(ilog,*)'@@ imult/ierror 2: ',ierror2
         if (ierror2.eq.0) then
             write(*,*)'@@ ZERO diagonal term AGAIN: give up'
             write(ilog,*)'@@ ZERO diagonal term AGAIN: give up'
             return
         endif
 83          format(1x,a,1x,6(g12.6,2x))
      return
      end
c
c
c     RANdom VECtor generator
      subroutine ranvec(wk,neq,iseed)
         implicit real*8 (a-h,o-z)
         real*8 wk(neq)
         real*4 ranval
c
cUNiX    use    ranv( )
cc       call msseed(iseed)
         init=-iseed
         iffrand=0
         call portran(init,ranval,iffrand)
         sum=0.0
         do 36 i=1,neq
            call portran(i,ranval,iffrand)
            wk(i)=ranval-0.5
            sum=sum+wk(i)*wk(i)
 36      continue
         enorm=sqrt(sum)
         do 37 i=1,neq
            wk(i)=wk(i)/enorm
 37      continue
      return
      end
c
c
c     COUNT negative udu Diagonal terms for COLumn storage
      subroutine countdCOL(a,maxstiff,nloc,neq,iless)
         implicit real*8 (a-h,o-z)
         real*8  a(maxstiff)
         integer nloc(neq)
c
         iless=0
         do 20 i=1,neq
            iloc=nloc(i)
            if (a(iloc) .lt. 0.0) iless=iless+1
 20      continue
      return
      end
c
c
      subroutine BANDAXD( matrix, vecin, vecout,neq,iband )
c     Multiplies   [ banded ]{vector] = {vector}
         implicit real*8 (a-h,o-z)
         real*8 vecout(neq), matrix(neq,iband), vecin(neq)
c
         do 10 i= 1, neq 
            sum=0.0
            do 20 j= max(1,i-iband+1),i-1
               ii = j
               jj = i-j+1
               sum = sum + matrix(ii,jj)*vecin(j)
 20         continue
            do 30 j= i, min( i+iband-1, neq)
               ii = i
               jj = j-i+1
               sum = sum + matrix(ii,jj)*vecin(j)
 30         continue
            vecout(i)=sum
 10      continue
      return
      end
c
c     subroutine msseed(iseed)
c        integer*2 i2
c        i2=iseed
c        call seed(i2)
c     return
c     end
cccc\end{verbatim}
c
c
c      PORTable RANdom number generator
       subroutine portran(idum,ranv,iffzz)
c         based on PFTV pp 196-197. RAN1
c
          save iff,ix1,ix2,ix3
          real*4 r(97)
          parameter (m1=259200, ia1=7141, ic1=54773, rm1=1.0/m1)
          parameter (m2=134456, ia2=8121, ic2=28411, rm2=1.0/m2)
          parameter (m3=243000, ia3=4561, ic3=51349)
          data iff /0/
c
          if (idum .lt. 0 .or. iff .eq. 0) then
             iff = 1
             ix1=mod(ic1-idum,m1)
             ix1=mod(ia1*ix1+ic1,m1)
             ix2=mod(ix1,m2)
             ix1=mod(ia1*ix1+ic1,m1)
             ix3=mod(ix1,m3)
             do 10 j=1,97
                ix1=mod(ia1*ix1+ic1,m1)
                ix2=mod(ia2*ix2+ic2,m2)
                r(j)=(float(ix1)+float(ix2)*rm2)*rm1
 10          continue
          endif
c
          ix1=mod(ia1*ix1+ic1,m1)
          ix2=mod(ia2*ix2+ic2,m2)
          ix3=mod(ia3*ix3+ic3,m3)
          j=1 + (97*ix3)/M3
          if (j .gt. 97 .or. j .lt. 1) then
              write(*,*)' !!! ERROR  !!!!'
              return
          endif
          ranv = r(j)
                r(j)=(float(ix1)+float(ix2)*rm2)*rm1
C
        return
        end
c
c
c
c     AxB product for COLumn storage
      subroutine AxBCOL(matrix,maxstiff,vecin,vecout,
     &                                  neq,iprof,iprof2,nloc)
         implicit real*8 (a-h,o-z)
         integer iprof(neq),iprof2(neq),nloc(neq)
         real*8 vecout(neq),matrix(maxstiff),vecin(neq)
         real*8 val,valmat
c
         do 10 i=1,neq
            jlim=max(1,(i-iprof(i)+1))
            do 20 j=jlim,i
               is=i
               io=i-j+1
               if (io .gt. iprof(is)) goto 20
               iloc = nloc(is) + io -1
               valmat=matrix(iloc)
               val = vecin(j)
               vecout(i)=vecout(i) + val*valmat
 20         continue
            jlim=min(iprof2(i),(neq-i+1))
            do 30 j=2,jlim
               is=i+j-1
               io=j
               if (io .gt. iprof(is)) goto 30
               iloc = nloc(is) + io -1
               valmat=matrix(iloc)
               val = vecin(i+j-1)
               vecout(i)=vecout(i) + val*valmat     
 30         continue
 10      continue
      return
      end
c
c
c     SUBspace iteration eigensolver for BUCKling
      subroutine sub_buck(stf,geom,load,wk,modemax,maxroot,
     &                   ar,br,eigval,eigvec,dj,r,d,ivib,
     &                   iprof,iprof2,nloc,istiff )
c        Finds first m eigenvalues and vectors using Subspace
c        iteration, based on Bathe pp 685-689.
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         character*12 str12
         character*35 str35
         character*4  str4 
         integer iprof(neq),iprof2(neq),nloc(neq)
         integer icluster(100)
         real*8 stf(maxstiff ),  geom(maxstiff) ,load(neq)
         real*8 wk(neq),r(neq,maxroot),dj(maxroot),d(maxroot)
         real*8 ar(maxroot,maxroot),br(maxroot,maxroot)
         real*8 eigval(maxroot),eigvec(maxroot,maxroot)
c
         real*8 zlam,rtol,rtolj, rt,rtolv, dzz,dif
         real*8 znorm,wnorm,vnorm,enorm
         real*8 sum,sumk,summ, eig,eigup,eignex, ev,sqrm
c
         zlam=0.0
         iconv=0
         rtol=1.0D-8
         rtolj=1.0D-12
         do i=1,maxroot
            d(i)=0.0
         enddo
c
c        STIFFNESS and MASS
*        if (istiff .eq. 0) then
*            write(*   ,*)'@@ reading     <<StaDyn.STF<<'
*            write(ilog,*)'@@ reading     <<StaDyn.STF<<'
*            rewind (istf)
*            read(istf) neqs,ibands
*            write(*   ,*)'@@ neq iband',neqs,ibands
*            write(ilog,*)'@@ neq iband',neqs,ibands
*            do 8 i=1,neq
*               iloc=nloc(i)
*               jmax=iprof(i)
*               read(istf) (stf(iloc+j-1), j=1,jmax)
*8           continue
*        endif
*        write(iout,'(a)') ' '
c
c            stability
*            write(*   ,*)'@@ reading    <<StaDyn.GEO<<'
*            write(ilog,*)'@@ reading    <<StaDyn.GEO<<'
*            rewind (igeo)
*            read(igeo) neqg,ibandg
*            write(*   ,*)'@@ neqg ibandg',neqg,ibandg
*            write(ilog,*)'@@ neqg ibandg',neqg,ibandg
*            do i=1,neq
*               iloc=nloc(i)
*               jmax=iprof(i)
*               read(igeo) (geom(iloc+j-1), j=1,jmax)
*181                   format(1x,20(g12.6,1x))
c               write(ilog,*) geom(iloc)
*            enddo
*            ibandm=ibandg
*            write(*,*)'@@ [K] [G] reloaded  OK'
c
             ibandm=iband
c
ctag1
c        establish norm based on diagonal
         sumk=0.0
         summ=0.0
         do 70 i=1,neq
            iloc=nloc(i)
            sumk=sumk+abs(stf(iloc))
            summ=summ+abs(geom(iloc))
*           ilocm=nloc(i)
*           if (ilump .eq. 1) ilocm=i
*           summ=summ+abs(geom(ilocm))
            load(i)=geom(iloc)/stf(iloc)
 70      continue
         sumk=sumk/neq
         summ=summ/neq
         znorm=sumk/summ
        write(*   ,'(a,3(g12.6))')' @@ NORMs s g: ',sumk,summ,znorm
!Fangbao        write(ilog,'(a,3(g12.6))')' @@ NORMs s g: ',sumk,summ,znorm
cX
c              NORMALIZE stiff and geom
               do i=1,maxstiff
                  stf(i)=stf(i)/sumk
               enddo
               do i=1,maxstiff
                  geom(i)=geom(i)/summ
               enddo
               dzz= (znorm/neq)*1.0e-12
*              dzz= (1.0  /neq)*1.0e3
               dzz=0
c
c       make copy of stiffness
        open(unit=itmp,form=fdn//'unformatted',status='scratch')
        rewind itmp
        do i=1,maxstiff
           write(itmp) stf(i)
        enddo
        write(*,*)'@@ made copy of normalized [K]'
!Fangbao        write(ilog,*)'@@ made copy of normalized [K]'
c
c       DECOMPOSE stiffness matrix
        ierror=0
cccc    call  udu   (stf,neq,iband,ier1)
*       write(ilog,* )'@@ ',stf(2512)
*       write(ilog,* )'@@ ',stf(nloc(2512))
        call  uduCOL_D(stf,maxstiff,neq,ierror,iprof,nloc,z0,z1)
!Fangbao        write(ilog,83)'@@ Diag min max: ',z0,z1,z1/z0
!Fangbao        write(ilog,* )'@@ imult/ierror: ',ierror
*       write(ilog,* )'@@ ',stf(nloc(2512))
c              !!!
        ierror=1
*       if (ierror .eq. 0) then
*           write(*,*)' '
*           write(*,*)'@@ ZERO diagonal term: try  shift'
*           write(ilog,*)'@@ ZERO diagonal term: try shift'
*           rewind itmp
*           do i=1,maxstiff
*              read(itmp) stf(i)
*           enddo
ccc         call zerodiag   (stf,          geom,
ccc                          neq,iband,ibandm,zlam,dzz,ilog,ierror2)
*           call zerodiagCOL(stf,maxstiff, geom,maxstiff,iprof,nloc,
*    &                       neq,ibandm,zlam, dzz,ilog,ierror2)
*           if (ierror2.eq.0) return
*       endif
c
c       CHECK if near RIGID body modes
        do n=1,neq
           iloc=nloc(n)
           wnorm=stf(iloc)
           if (abs(wnorm) .le. 10e-10) then
               write(*,*)' '
!Fangbao               write(ilog,*)'@@ test for rigid modes', wnorm,n
               write(*   ,*)'@@ test for rigid modes', wnorm,n
               write(*,*)' '
               write(*   ,*)'@@ !!! possible RIGID modes !!!'
!Fangbao               write(ilog,*)'@@ !!! possible RIGID modes !!!'
               write(*,*)' '
*              write(*,*)'CHOOSE:  0=return 1=ignore  2=shift'
*              call zzwrt('  -->  ')
*              read(ikbd,*) ignore
               ignore=2
               ignore=1
*              write(ilog,*) ignore,'  ::0=ignore'
               if (ignore .eq. 0) then
                   return
               elseif (ignore .eq. 1) then
c                  do nothing  
               elseif (ignore .eq. 2) then
                   dzz=abs(z0)*4
                   dzz=abs(z1)/neq
                   dzz=abs(z0)*40
                   dzz=dzz+1.0e-10
                   rewind itmp
                   do i=1,maxstiff
                      read(itmp) stf(i)
                   enddo
***                call zerodiag(stf,geom,neq,iband,ibandm,
***  >                           zlam,dzz,ilog,ierror2)
                 call zerodiagCOL(stf,maxstiff,geom,maxstiff,iprof,nloc,
     &                          neq,ibandm,zlam, dzz,ilog,ierror2)
                   if (ierror2.eq.0) return
               endif
               goto 19
           endif
        enddo
        write(*,*)'@@ rigid body modes checked: OK'
!Fangbao        write(ilog,*)'@@ rigid body modes checked: OK'
 19     continue
c
c       INITIAL load vector: first column is geom diag
        sum=0.0
        do 20 i=1,neq  
           iloc=nloc(i)
           if (ilump .eq. 1) iloc=i
           sum=sum+geom(iloc)*geom(iloc)
 20     continue
        sqmas=sqrt(sum)
        do 201 i=1,neq  
           iloc=nloc(i)
           if (ilump .eq. 1) iloc=i
           r(i,1)=geom(iloc)/sqmas
 201    continue
        do 22 i=1,neq
           do 24 j=2,maxroot
              r(i,j)=0.0
 24        continue
 22     continue
c
c       find largest of M/K
        nd=neq/maxroot
        l=neq
c       helps distribution for equal values of M/K
        rt_large=-1.0d12
        do 30 j=2,maxroot
           l=l-nd
           rt=rt_large
           do 32 i=1,l  
              if (load(i) .ge. rt) then   
                  rt=load(i)
                  iloc=i
              endif
 32        continue
           do 34 i=l,neq
              if (load(i) .gt. rt) then 
                  rt=load(i)
                  iloc=i
              endif
 34        continue
           r(iloc,j)=1.0
           load(iloc)=rt_large
 30     continue
c
c       RANDOM vector 
        call ranvec(wk,neq,10)
        do 37 i=1,neq
           r(i,maxroot)=wk(i)
 37     continue
ctag2
*          write(*,*)'Rmatrix '
*          write(*,*)'     '
*          do i=1,neq
*             write(*,182) (r(i,j), j=1,maxroot)
*          enddo
 182       format(1x,30(g10.4,1x))
c
        write(*,*)'     '
c       ITERATIONs begin
        iter=0
 100    continue
        iter=iter+1
        write(*,'(1x,a,i4)') '@@ ITERATION: ',iter
c tag2
c
c       SOLVE {u}k+1   &    simultaneously reduce K 
        do 110 j=1,maxroot
ccc        call  bak   (stf,         r(1,j) ,neq,iband,wk,ier1)
           call  bakCOL_D(stf,maxstiff,r(1,j) ,neq,wk,ierror,
     &                                           iprof,iprof2,nloc)
           do 130 i=j,maxroot
              sum=0.0
c             since [K]{u}k+1 = [r]k then only pre-mult
              do 140 k=1,neq
                 sum=sum+wk(k)*r(k,i)
 140          continue
              ar(i,j)=sum
              ar(j,i)=sum
 130       continue
c
c          store the solution vectors for later use with geom
           do 150 i=1,neq
              r(i,j)=wk(i)
 150       continue
 110    continue
c
c       MULTIPLY  {u}[M]{u} for reduced M
        do 160 j=1,maxroot
ccc        call bandAxd(geom, r(1,j) ,wk,neq,ibandm)
           if (ibandm .eq. 1) then
               do i=1,neq
                  wk(i)= geom(i)*r(i,j)
               enddo
           else
               do i=1,neq
                  wk(i)= 0.0
               enddo
               call AxBCOL(geom,maxstiff,r(1,j),wk,
     &                          neq,iprof,iprof2,nloc)
           endif
c
           do 180 i=j,maxroot
              sum=0.0
              do 190 k=1,neq
                 sum=sum+r(k,i)*wk(k)
 190          continue
              br(i,j)=sum
              br(j,i)=sum
 180       continue
           if (iconv .gt. 0) goto 160
c
c          store [M]{u} as new R load vectors
           do 200 i=1,neq
              r(i,j)=wk(i)
 200       continue
 160    continue
c
*       stop
c       SOLVE reduced eigenvalue prob by JACOBI rotations
        nsmax=15
        n=maxroot
        rtolj=1e-12
        call jacobi(ar,br,eigval,dj,n,rtolj,nsmax,n,eigvec,n,ilog)
        write(*,'(1x,a,1x,i4)') '@@ # of sweeps = ',nsmax
*       stop
        call eigsrt(eigval,eigvec,n)
*       call eigsrt_neg(eigval,eigvec,n)
*               write(*,*)'vals '
*               do i=1,maxroot
*                  write(*,'(1x,20(g11.5,1x))') eigval(i)
*               enddo
c
c       calculate R times eigenvectors for new loads
        do  420 i=1,neq
            do 422 j=1,maxroot
               dj(j)=r(i,j)
 422        continue
            do 424 j=1,maxroot
               sum=0.0
               do 430 k=1,maxroot
                  sum=sum+dj(k)*eigvec(k,j)
 430           continue
               r(i,j)=sum
 424         continue
 420    continue
c       enforce random last vector
        call ranvec(wk,neq,iter)
        do 1137 i=1,neq
           r(i,maxroot)=wk(i)
 1137   continue
*          write(*,*)'Rmatrix '
*          write(*,*)'     '
*          do i=1,neq
*             write(*,182) (r(i,j), j=1,maxroot)
*          enddo
*       stop
c
        if (iconv .gt. 0) goto 500
c
c       check for convergence
        do 380 i=1,modemax
           dif=abs(eigval(i)-d(i))
           rtolv=dif/eigval(i)
           if ( abs(rtolv) .gt. rtol) then
               write(*,'(1x,a,i5)')'@@ TRIGGER rtolv: ',i
               write(*,*)'@@      ',rtolv,' vs ',rtol
!Fangbao               write(*,*)'@@ Eigv_1 =',eigval(1)*sumk/summ
!Fangbao     &                  ,  ' Eigv_2 =',eigval(2)*sumk/summ
!Fangbao               write(ilog,*)'@@ Eigv_1 =',eigval(1)*sumk/summ
!Fangbao     &                  ,  ' Eigv_2 =',eigval(2)*sumk/summ
               maxout=maxroot
               if (maxout .gt. 45) maxout=45
               write(idyn,85) iter*1.0,(eigval(j)*sumk/summ, j=1,maxout)
 85            format(1x,50(g12.6,1x))
!               call zzflush(idyn)
               goto 400
           endif
 380    continue
!Fangbao        write(ilog,*)'@@ NORMAL SUBspace convergence'
        write(*,*)'@@ NORMAL SUBspace convergence'
        write(*,*)'@@ another round'
        iconv=1
        goto 100
c
 400    continue
!Fangbao        write(ilog,*)'@@ # of sweeps = ',nsmax,
!Fangbao     &                 '    TRIGGER rtolv: ',i
c       see if max iterations have occurred
*       stop
        if (iter .eq. itermax) then
            if (i .eq. 1) then
                write(*,*)'@@ NO convergence at ITERs=',itermax
            else
                write(*,*)'@@ INCOMPLETE convergence at ITERs=',itermax
            endif
            write(*,*)'@@ another round'
            iconv=2
            goto 100
        else
            do 440 i=1,maxroot
               d(i)=eigval(i)
 440        continue
            goto 100
        endif
c       ITERATIONs end
c
 500                   continue
c
c       CHECK quality
        isuspect=0
                   rewind itmp
                   do i=1,maxstiff
                      read(itmp) stf(i)
                   enddo
        write(*   ,*)'@@ checking for SUSPECTs:'
!Fangbao        write(ilog,*)'@@ checking for SUSPECTs:'
!Fangbao        write(ilog,*)'@@ zlam: ',zlam
        do 580 i=1,maxroot
           eig=eigval(i)-zlam
ccc        call bandAxd(stf, r(1,i) ,wk ,neq,iband)
               do k=1,neq
                  wk(k)= 0.0
               enddo
               call AxBCOL(stf,maxstiff,r(1,i),wk,
     &                          neq,iprof,iprof2,nloc)
           vnorm=0.0
           do 590 n=1,neq
              vnorm=vnorm+wk(n)*wk(n)
 590       continue
ccc        call bandAxd(geom, r(1,i) , load  ,neq,ibandm)
*          if (ibandm .eq. 1) then
*              do k=1,neq
*                 load(k)= geom(k)*r(k,i)
*              enddo
*          else
               do k=1,neq
                  load(k)= 0.0
               enddo
               call AxBCOL(geom,maxstiff,r(1,i),load,
     &                          neq,iprof,iprof2,nloc)
*          endif
c
           wnorm=0.0
           do 600 n=1,neq
              wk(n)=wk(n)-eig*load(n)
              wnorm=wnorm+wk(n)*wk(n)
 600       continue
           vnorm=sqrt(abs(vnorm))
           wnorm=sqrt(abs(wnorm))
           if (vnorm .lt. 10e-12) vnorm=10e-12
           enorm=wnorm/vnorm
           dj(i)=0
           icluster(i)=i
           write(str35,'(1x,a,i5,1x,g13.6,1x,a)')
     &                  '@@ QUALity: ',i,enorm,' '
           call zzwrt(str35)
           if (enorm .gt. .01) then
               isuspect=isuspect+1
               icluster(i)=0
               str12='  <--suspect'
               dj(i)=100
               write(*,'(a)') ' <-- suspect'
ccc            write(iout,'(a)') ' <-- suspect'
           else
               str12=' '
               write(*,'(a)') ' '
ccc            write(iout,'(a)') ' '
           endif
!Fangbao           write(ilog,*)'@@ NORM: ',i,enorm,str12 
 580    continue
!Fangbao           write(ilog,'(a)')' '
c tag3
        write(*,*)' '
!Fangbao        write(ilog,*)'@@ SUSPECT eigenvalues ',isuspect
        write(*   ,*)'@@ SUSPECT eigenvalues ',isuspect
        write(ilog,*)' '
c
c       FLAG spurious eigens and count others
*       ibound=0
*       do i=1,maxroot
*          ic=icluster(i)
*          if (ic .gt. 0) then
*              ibound=ibound+1
*              icluster(i)=ibound
*          endif
*       enddo
c
c       STURM SEQUENCE for total number
*       if (ivib .eq. 2) then
c           Decompose effective stiffness matrix
*           eigup=1.01*eigval(modemax)
c
*           inc=0
*399        inc=inc+1
*           eignex=0.99*eigval(modemax+inc)
*           if ( abs(eigup) .ge. abs(eignex) ) goto 399
*           write(*,*)'@@ SHIFT total: ',eignex
*           maxeigs=icluster(modemax+inc-1)     
c
ccc     call shiftip (stf,         geom,    neq,iband,ibandm,eignex  )
*       call shiftCOL(stf,maxstiff,geom,maxstiff,nloc,
*    &                                      neq,ibandm,eignex)
*       ierror=0
cccc    call  udu   (stf,         neq,iband,ier1)
*       call  uduCOL_D(stf,maxstiff,neq,ierror,iprof,nloc,z0,z1)
*       write(ilog,83)'@@ Diag min max: ',z0,z1,z1/z0
*       if (ierror .eq. 0) then
*           write(*,*)'@@ ZERO diagonal term AGAIN: give up'
*           return
*       endif
*       write(iout,'(a)') '  '
ccc     call countd   (stf,              neq,iband,itotal)
*       call countdCOL(stf,maxstiff,nloc,neq,itotal)
*       write(ilog,*)'@@ should be',maxeigs,' actually',itotal
*       write(*   ,*)'@@ should be',maxeigs,' actually',itotal
*           if (itotal .gt. maxeigs) then
*               write(*,*)'@@ MISSing eigenvalues'
*           elseif (itotal .eq. maxeigs) then
*               write(*,*)'@@ OK eigenvalues'
*           elseif (itotal .lt. maxeigs) then
*               write(*,*)'@@ SPURious eigenvalues'
*           endif
*       endif
c
*       write(*,*)'@@ EIGENvalues & vectors in  >>StaDyn.SNP>>'
cc       write(ilog,*)'@@ snap isnp ',snap,isnp
*       rewind(isnp)
*       mdx=0
*       do 510 i=1,modemax
*          if (dj(i) .gt. 1.0) goto 510
*          mdx=mdx+1
*510    continue
*       modemax=mdx
c
c       rectify for stability problems
        rewind(isnp)
        sgn=1.0
*       if (ivib .eq. 1) sgn=-1.0
        sgn=-1.0
        do i=1,modemax
           ev=(eigval(i)-zlam)*sumk/summ
           eigval(i)=ev*sgn
           sqrm=sqrt(summ)
           do j=1,neq
              r(j,i) = r(j,i)/sqrm
           enddo
        enddo
c
c        write(isnp) modemax
c        do i=1,modemax
c*          ev=(eigval(i)+zlam)*sumk/summ
c*          sqrm=sqrt(summ)
c*          write(isnp) ev,(r(j,i)/sqrm, j=1,neq)
c           write(isnp) eigval(i),(r(j,i), j=1,neq)
c                  write(ilog,*)'@@ value: ',i,eigval(i)
c        enddo
        do i=1,modemax
           write(isnp) eigval(i),neq,modemax
           do j=1,neq
              write(isnp) r(j,i)
           enddo
           write(ilog,*)'@@ value: ',i,eigval(i)
        enddo
c
c            echo to OUT file
*            pi2=2*4.0*atan(1.0)
             write(iout,*)'EIGENvalues: '
*            if (ivib .eq. 1) then
                 write(iout,*)'   N    lambda [ xload]   '
                 do i=1,modemax
***                 ev=(eigval(i)+zlam)*sumk/summ
                    write(iout,81) i,eigval(i)
                 enddo   
*            elseif (ivib .eq. 2) then
*                write(iout,*)'   N    omega [r/s]    freq [Hz] '
*                do i=1,modemax
***                 ev=(eigval(i)+zlam)*sumk/summ
***                 ev=sqrt(abs(ev))
*                   str4='real'
*                   if (eigval(i) .lt. 0.0) str4='imag'
*                   ev=sqrt(abs(eigval(i)))
*                   write(iout,82) i,ev,ev/pi2,str4
*                enddo   
*            endif
 81          format(1x,i5,1x,2(g12.6,2x))
 82          format(1x,i5,1x,2(g12.6,2x),1x,a)
 83          format(1x,a,1x,6(g12.6,2x))
c 
c
      close(itmp)
      return
      end
c
