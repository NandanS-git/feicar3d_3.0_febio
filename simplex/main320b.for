ccc
c     main for Simplex
      subroutine nonstad_main(irank)
         implicit real*8 (a-h,o-z)
         character*40 filena
c
         include 'commons.std'
c
         parameter (ippp=101,idload=1000)
         real*8   dloadi(idload),floadi(idload),floadm(idload)
         real*8   prop_D(ippp,10)
         integer   nmprop [allocatable] (:)
         integer     nloc [allocatable] (:)
c
         real*8 respce8 [allocatable] (:)
         real*8 respceb [allocatable] (:)
         real*8 respcec [allocatable] (:)
         real*8 wk [allocatable] (:)
         integer irespce [allocatable] (:)
         real*8 xyz_D [allocatable] (:)
         real*8 tforce [allocatable] (:), grav [allocatable] (:)
         real*8 tforce1 [allocatable] (:,:),tforce2 [allocatable] (:,:)
         integer idbc[allocatable] (:,:)
         integer npijkm [allocatable] (:,:), npj [allocatable] (:),
     &          nelt [allocatable] (:), npk [allocatable] (:),
     &                                  npm [allocatable] (:)
         integer iprofv [allocatable] (:), iprofh [allocatable] (:)
         character*40 str40
         character*2 st
*        STATIC respce8,respceb,wk
*        STATIC tforce,tforce1,tforce2
c
*        str40 = ' Simplex version 4.50, February 2003'
*        str40 = ' Simplex version 4.52, March 2003'
*        str40 = ' Simplex version 4.60, September 2003'
*        str40 = ' Simplex version 4.72, February 2004'
*        str40 = ' Simplex version 4.74, March 2005'
*        str40 = ' Simplex version 4.80, April 2005'
*        str40 = ' NSol3D version 2.02, April 2005'
*        str40 = ' NSol3D version 2.08, May 2005'
*        str40 = ' Simplex version 2.20, June 2005'
*        str40 = ' Simplex version 2.30, August 2005'
*        str40 = ' Simplex version 2.34, September 2005'
*        str40 = ' Simplex version 2.38, October 2005'
*        str40 = ' Simplex version 2.40, December 2005'
*        str40 = ' Simplex version 2.42, January 2006'
*        str40 = ' Simplex version 2.44, February 2006'
*        str40 = ' Simplex version 2.46, March 2006'
*        str40 = ' Simplex version 2.48, April 2006'
*        str40 = ' Simplex version 2.50, April 2006'
*        str40 = ' Simplex version 2.52, August 2006'
*        str40 = ' Simplex version 2.54, September 2006'
*        str40 = ' Simplex version 2.60, October 2006'
*        str40 = 'Simplex version 2.62, January 2007'
*        str40 = 'Simplex version 2.64, October 2007'
*        str40 = 'Simplex version 2.70, November 2007'
*        str40 = 'Simplex version 2.80, January 2008'
*        str40 = 'Simplex version 3.00, February 2009'
*        str40 = 'Simplex version 3.02, March 2009'
*        str40 = 'Simplex version 3.12, October 2010'
*        str40 = 'Simplex version 3.20, October 2010'
         str40 = 'Simplex version 3.20b, July 2011'
c
      write(*,'(///)')
      write(*,*) 'SimplexSimplexSimplexSimplexSimplexSimplexSimplex  '
      write(*,*) ' An FEM Program for Nonlinear 3-D Solids '
      write(*,*) '  ',str40 
      write(*,*) '       (c) 2003-11  ikayex SoftWare Tools'
      write(*,*) 'SimplexSimplexSimplexSimplexSimplexSimplexSimplex  '
      write(*,*) ' '
c
           ikbd=15
!	     ikbd0=15
           ierr=6
           itmp=7
           icfg=8
             idat=12
           ilog=21
           iout=23 
!           idyn=24 
           idyn=24+irank*200
           imon=228
           imon=79  
           isnp=25
           ioutps=26
           iall=27
             istf=51
             imss=52
             idis=54
             ilod=53
           icms=78
           igeo=77
c
c     OPEN files ready for work
c
        write(st,'(i1,a)') irank,'/' 
        print*,'Home directory for the solid:', st
        fdn = st
            open( istf, file=st//'stadyn.stf', form='unformatted')
            open( imss, file=st//'stadyn.mas', form='unformatted')
            open( idis, file=st//'stadyn.dis', form='unformatted')
            open( igeo, file=st//'stadyn.geo', form='unformatted')
            open( ilod, file=st//'stadyn.lod', form='unformatted')
            open( icms, file=st//'stadyn.cms', form='unformatted')
            open( isnp, file=st//'stadyn.snp', form='unformatted')
        open( ilog , file=st//'simplex.log')
        open( iout , file=st//'stadyn.out')
        open( idyn , file=st//'stadyn.dyn')
        open( imon , file=st//'stadyn.mon')
        rewind(isnp)
        rewind(ilog )
        rewind(iout )
        rewind(idyn )
        rewind(imon )

        open(ikbd, file=st//'instab', form='formatted')  !F.-B. Tian
        rewind(ikbd)                                     !F.-B. Tian
c
!Fangbao        write(ilog,*)'@@ ',str40 
!        call timday(ilog)
        call initial (ilog,isize8, fdn,
     &                    maxnode,maxelem,maxforce,maxpts,ilump,
     &                    itermax,rtol)
          isize4=maxelem*12
          isize4=maxelem*24
          isize4=maxelem*120
          maxdof=maxnode*3
          maxxyz=maxnode*3
          isizec=maxelem*6*27*5 + maxxyz*2 + maxdof*7
c
        allocate (respce8(isize8),                irespce(isize4),
     &                                                 stat=ierr01)
!Fangbao               write(ilog,*)'@@ error respce8 ',ierr01 
        allocate (respceb(isize8),
     &                                                 stat=ierr01)
!Fangbao               write(ilog,*)'@@ error respceb ',ierr01 
        allocate (respcec(isizec),      
     &                                                 stat=ierr01)
!Fangbao               write(ilog,*)'@@ error respcec ',ierr01 
        allocate (npijkm(maxelem,21), npj(maxelem), nelt(maxelem),
     &          npk(maxelem),npm(maxelem),idbc(maxnode,3),
     &          iprofv(maxdof),
     &          iprofh(maxdof),nloc(maxdof), nmprop(maxelem),
     &                                            stat=ierr02)
        allocate (tforce(maxforce),wk(maxdof), xyz_D(maxxyz),
     &            tforce1(maxforce,2), tforce2(maxforce,2), 
     &            grav(maxdof),               stat=ierr03)
c
           if (ierr01 .eq. 0 .AND. ierr02 .eq. 0
     &                       .AND. ierr03 .eq. 0 ) then
!Fangbao               write(ilog,*)'@@ ALLOCATION succeeded ' 
           else 
!Fangbao               write(ilog,*)'@@ !!! ALLOCATION failed !!!' 
!Fangbao               write(ilog,*)'@@ decrease sizes in <<simplex.CFG>>'
!Fangbao               write(ilog,*)'@@ ERRORs 1 2 3',ierr01,ierr02,ierr03
               write(*   ,*)'@@ !!!                   !!!' 
               write(*   ,*)'@@ !!! ALLOCATION failed !!!' 
               write(*   ,*)'@@ !!!                   !!!' 
               write(*   ,*)'@@ decrease sizes in <<simplex.CFG>>'
               goto 999
           endif
c
c      in case datafile not read for some operations
       do n=1,maxelem
          nmprop(n)=1
       enddo
c
 1    continue
c
      write(*,*) '    MAIN menu: '
      write(*,*) '          0: Quit '
      write(*,*) ' '
      write(*,*) '          2: Read in structure DataFile'
      write(*,*) ' '
      write(*,*) '             <<imposed large deformtion>> '
      write(*,*) '        200: Incremental 3-D  '
      write(*,*) ' '
      write(*,*) '             <<nonlinear elastic, elastic-plastic>>'
      write(*,*) '        400: Incremental 3-D  '
      write(*,*) ' '
      write(*,*) '        700: PostScript plot files'
      write(*,*) ' '
      call zzwrt(' SELECT--> ')
       read(ikbd,*,err=1) imain
       write(ilog,*) imain,' ::MAIN'
c
      if (imain .eq. 0) then
          open(unit=icfg,file=fdn//'simplex.cfg')
          rewind icfg
          write(icfg,*) isize8,'     ::MaxStorage'
          write(icfg,*) maxelem,'     ::MaxElem'
          write(icfg,*) maxnode,'     ::MaxNode'
          write(icfg,*) maxforce,maxpts,'     ::MaxForce'
          write(icfg,*) ilump,'      ::ilump'
          write(icfg,*) itermax,'      ::itermax'
          write(icfg,*) rtol ,'      ::rtol '
!          call timday(icfg)
          goto 999
      endif
      if (imain .eq. 2) goto 200
      if (imain .eq. 700) goto 7700
c
c
ccccccccccccccccccccccccccccccccc  nonlinear cccccccccccccccccccc
      if (imain .eq. 400) then
          write(*,*) '@@   '
          write(*,*) '      410: Incremental Static  E   3-D TL '
          write(*,*) '      420: Incremental Dynamic E   3-D TL '
          write(*,*) '      424:    ..         ..    E   3-D TL mult '
          write(*,*) '      440: Incremental Dynamic E   3-D EXP '
          write(*,*) '      442:    ..         ..    gen 3-D EXP '
          write(*,*) '      450: Incremental Dynamic  Rubber 3-D TL'
          write(*,*) '      452:    ..         ..       ..   3-D EXP'
          write(*,*) '      470: Incremental Dynamic  A-B tissue 3-D TL'
          write(*,*) '      472:    ..         ..       ..   3-D EXP'
          write(*,*) '      474:    ..         ..       ..   vis TL '
          write(*,*) '      476:    ..         ..       ..   vis EXP'
          call zzwrt(' SELECT--> ')
          read(ikbd,*,err=1) imain2
          write(ilog,*) imain2,' ::MAIN2'
          if (imain2 .eq. 410) goto 4100
          if (imain2 .eq. 420) goto 4200
          if (imain2 .eq. 424) goto 4240
          if (imain2 .eq. 440) goto 4400
          if (imain2 .eq. 442) goto 4420
          if (imain2 .eq. 450) goto 4500
          if (imain2 .eq. 452) goto 4520
          if (imain2 .eq. 470) goto 4700
          if (imain2 .eq. 472) goto 4720
          if (imain2 .eq. 474) goto 4740
          if (imain2 .eq. 476) goto 4760
          goto 1
      endif
c
c
c  READ structure's  data, and  bcs
 200    continue
        write(*,'(a)')' @@'
        call zzwrt(' INPUT:  structure datafilename-->  ')
        read(ikbd ,'(a40)') filena
        write(ilog,*) filena
        open(idat, file = fdn//adjustl(filena))
        rewind(idat)
c
c
!        call zztimer(ilog,'--> DATAin')
!        call zzflush(ilog)
        call readdt(ilog, 
     &             xyz_D,
     &             respce8,isize8,
     &             idbc, maxdof, maxelem, maxnode,
     &             npijkm, iglobal, iprofv,iprofh,
     &             prop_D,nmprop,ippp)
c
!        call zztimer(ilog,'<-- DATAin')
!        call zzflush(ilog)
        close(idat)
c
c       compute locations of cols in vector
        nloc(1)=1
        do i=1,neq-1
           nloc(i+1) = nloc(i) + iprofv(i)
        enddo
!       do i=1,neq
!               write(iout,*) i,nloc(i),' nloc',iprofv(i)
!       enddo
!       maxstiff=nloc(neq)+iprofv(neq)-1
        percent=(maxstiff*100.0)/(neq*iband*1.0)
        ipercent=percent
c
        write(*,212)'@@ system size = [',neq,' X ',iband,']'
        write(*,213)'@@ Prof Storage=  ',maxstiff
        write(*,213)'@@ % of NB     =  ',ipercent
 212    format(a,i5,a,i5,a)
 213    format(a,i9,a,i5,a)
        istiff=0
        ilump=1
        maxmass=neq
!       if (imain .eq. 42) then
!           ilump=2
!           maxmass=maxstiff
!       endif
        do i=1,neq
           grav(i)=0.0
        enddo
c       gravity load vector
        if (igravity_on .eq. 1) then
                call body_force_grav( iglobal,
     &                   npijkm,idbc,maxnode,maxelem,
     &                  xyz_D(1),xyz_D(maxnode+1),xyz_D(maxnode*2+1),
     &                   prop_D,nmprop,ippp,
     &                   grav
     &                )
            write(ilog,*)'@@ formed gravity vector'
            write(*   ,*)'@@ formed gravity vector'
                z1=0.0
                do i=1,neq
!                  write(iout,*) i, grav(i)
                   z1=z1+grav(i)
                enddo
                write(ilog,*)'@@ Total gravity force: ',z1
        endif
        igravity_on=1
        goto 1
c
c
cccccccc
c
cccccccccccccccccccccccccccccccccccccccc   nonlinear    cccccccccc
c
 4000   continue
 4100   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
        write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1 + maxstiff    !K_E   
           ir3 = ir2 + maxstiff    !K_G
           ir4 = ir3 + neq         !wk
           ir5 = ir4 + neq         !{P}
           ir6 = ir5 + neq         !{u}
           ir7 = ir6 + neq         !fmag
           ir8 = ir7 + neq         !{u}_t
           ir9 = ir8 + maxstiff    !K_eff
           ir10= ir9 + neq  
           call  ckSpace(ilog,isize8,'INCrement',ir10)
c
!           call zztimer(ilog,'--> INCrement')
           write(ilog,*)'@@ global ',iglobal
           call non_TL3D( respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt)
!           call zztimer(ilog,'<-- INCrement')
c
        elseif (ianal .eq. 2) then
           goto 4900
c
        elseif (ianal .eq. 3) then
*           goto 6400
        endif
        goto 4100
cccc
cccccccccccccccccccc  dynamic  cccccccccccccccccccccccccccccc
 4200   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
        write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1 + maxstiff    !K_E   
           ir3 = ir2 + maxstiff    !K_G
           ir4 = ir3 + neq         !wk
           ir5 = ir4 + neq         !{P}
           ir6 = ir5 + neq         !{u}
           ir7 = ir6 + neq         !fmag
           ir8 = ir7 + neq         !{u}_t
           ir9 = ir8 + maxstiff    !K_eff
           ir10= ir9 + maxstiff    !mass
           ir11= ir10+ neq         !vel
           ir12= ir11+ neq         !acc
           ir13= ir12+ neq         !uup
           ir14= ir13+ neq         !vvp
           ir15= ir14+ neq         !aap
           call  ckSpace(ilog,isize8,'INCrement_dyn',ir15)
c
           isize_elem=nel
           ic1 = 1
           ic2 = ic1  + nel*6*27    !stresst
           ic3 = ic2  + nel*6*27    !straint
           ic4 = ic3  + nel*6*27    !stressi
           ic5 = ic4  + nel*6*27    !straini
           ic6 = ic5  + nel*6*27    !stressci
           isize_dof=neq
           ic7 = ic6  + neq         !ui        
           ic8 = ic7  + neq         !dui        
           ic9 = ic8  + neq         !grav      
           ic10= ic9  + neq         !fmag3     
           ic11= ic10 + neq         !dloadi    
           ic12= ic11 + neq         !floadi    
           ic13= ic12 + neq         !floadm    
           ic14= ic13 + maxxyz      !xyzt      
           ic15= ic14 + maxxyz      !xyzi      
c
c
!           call zztimer(ilog,'--> INCrement_dyn')
           write(ilog,*)'@@ global ',iglobal
*          write(*,*) maxstiff,maxnode,neq,' xxxx'
           call non_TL3D_dyn( respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13), respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt
     &                  ,respceb,isize8
     &          ,isize_dof,respcec(ic6)
     &          ,respcec(ic7 ),respcec(ic8 ),respcec(ic9) 
     &          ,respcec(ic10),respcec(ic11)
     &          ,respcec(ic12),respcec(ic13),respcec(ic14) 
     &                       )
!           call zztimer(ilog,'<-- INCrement_dyn')
c
        elseif (ianal .eq. 2) then
           goto 4900
c
        elseif (ianal .eq. 3) then
*           goto 6400
        endif
        goto 4200
cccc
cccccccccccccccccccc  dynamic mult loads  ccccccccccccccccccccccccccc
 4240   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
        write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1 + maxstiff    !K_E   
           ir3 = ir2 + maxstiff    !K_G
           ir4 = ir3 + neq         !wk
           ir5 = ir4 + neq         !{P}
           ir6 = ir5 + neq         !{u}
           ir7 = ir6 + neq         !fmag
           ir8 = ir7 + neq         !{u}_t
           ir9 = ir8 + maxstiff    !K_eff
           ir10= ir9 + maxstiff    !mass
           ir11= ir10+ neq         !vel
           ir12= ir11+ neq         !acc
           ir13= ir12+ neq         !uup
           ir14= ir13+ neq         !vvp
           ir15= ir14+ neq         !aap
           ir16= ir15+ maxforce*40 !mul loads
           call  ckSpace(ilog,isize8,'INCrement_dyn',ir15)
c
!           call zztimer(ilog,'--> INCrement_dyn')
           write(ilog,*)'@@ global ',iglobal
*          write(*,*) maxstiff,maxnode,neq,' xxxx'
           call non_TL3D_dyn_mul( respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13), respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt
     &                  ,respceb,isize8
     &                  ,respce8(ir15) 
     &                       )
!           call zztimer(ilog,'<-- INCrement_dyn')
c
        elseif (ianal .eq. 2) then
           goto 4900
c
        elseif (ianal .eq. 3) then
*           goto 6400
        endif
        goto 1
cccccccccccccccccccc  dynamic explicit  cccccccccccccccccccccccccccccc
 4400   continue
!===============updated by Tian
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
        write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1    !K_E   
           ir3 = ir2    !K_G
           ir4 = ir3 + neq         !wk
           ir5 = ir4 + neq         !{P}
           ir6 = ir5 + neq         !{u}
           ir7 = ir6 + neq         !fmag
           ir8 = ir7 + neq         !{u}_t
           ir9 = ir8               !K_eff
           ir10= ir9 + maxmass    !mass
           ir11= ir10+ neq         !vel
           ir12= ir11+ neq         !acc
           ir13= ir12+ neq         !uup
           ir14= ir13+ neq         !vvp
           ir15= ir14+ neq         !aap
           call  ckSpace(ilog,isize8,'INCrement_dyn',ir15)
c
!c           call zztimer(ilog,'--> INCrement_dyn')
           write(ilog,*)'@@ global ',iglobal
           write(*,*) maxstiff,maxnode,neq,' xxxx'
           call non_3D_exp_Tian(
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13), respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt)
!c           call zztimer(ilog,'<-- INCrement_dyn')
c

        elseif (ianal .eq. 2) then
            goto 4900
c
        elseif (ianal .eq. 3) then
*           goto 6400
        endif
        goto 1
!==============================
c
!        write(*,'(a)')' @@'
!        write(*,*) 'CHOOSE:  0=return  '
!        write(*,*) '         1=analysis  '
!        write(*,*) '         2=post_analysis full snapshot'
!        write(*,*) '         3=post_analysis partial snapshot'
!        call zzwrt(' --> ')
!        read(ikbd,*) ianal
!        write(ilog,*) ianal,'  :: 1=analysis'
c
!        ibandm=iband
!        if (ianal .eq. 0) then
!            goto 1
!        elseif (ianal .eq. 1) then
c          storage pointers
!           ir1 = 1
!           ir2 = ir1 + maxstiff    !K_E   
!           ir3 = ir2 + maxstiff    !K_G
!           ir4 = ir3 + neq         !wk
!           ir5 = ir4 + neq         !{P}
!           ir6 = ir5 + neq         !{u}
!           ir7 = ir6 + neq         !fmag
!           ir8 = ir7 + neq         !{u}_t
!           ir9 = ir8 + maxstiff    !K_eff
!           ir10= ir9 + maxstiff    !mass
!           ir11= ir10+ neq         !vel
!           ir12= ir11+ neq         !acc
!           ir13= ir12+ neq         !uup
!           ir14= ir13+ neq         !vvp
!           ir15= ir14+ neq         !aap
!           call  ckSpace(ilog,isize8,'INCrement_dyn',ir15)
c
!           call zztimer(ilog,'--> INCrement_dyn')
!           write(ilog,*)'@@ global ',iglobal
*          write(*,*) maxstiff,maxnode,neq,' xxxx'
!           call non_3D_exp( respce8(ir1), respce8(ir2),
!     &                   respce8(ir3),
!     &                   respce8(ir4), respce8(ir5), 
!     &                   respce8(ir6), respce8(ir7),
!     &                   respce8(ir8), respce8(ir9), 
!     &                   respce8(ir10), respce8(ir11), 
!     &                   respce8(ir12), respce8(ir13), respce8(ir14), 
!     &                   idbc,npijkm,maxnode,maxelem, 
!     &                   iprofv, iprofh,nloc ,
!     &                   xyz_D ,maxxyz,
!     &                   iglobal, ippp,
!     &                   prop_D, nmprop,nelt)
!           call zztimer(ilog,'<-- INCrement_dyn')
c
!
!        elseif (ianal .eq. 2) then
!            goto 4900
c
!        elseif (ianal .eq. 3) then
*           goto 6400
!        endif
!        goto 1
c
cccccccccccccccccccc  dynamic explicit general  cccccccccccccccccccccc
 4420   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
        write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1 + nel*6*27    ![stress_ip]
           ir3 = ir2 + neq         !{F}
           ir4 = ir3 + neq         !wk
           ir5 = ir4 + neq         !{P}
           ir6 = ir5 + neq         !{u}
           ir7 = ir6 + neq         !fmag
           ir8 = ir7 + neq         !{u}_t
c3          ir9 = ir8 + maxstiff    !K_eff
           ir9 = ir8 + 1                      
           ir10= ir9 + maxstiff    !mass
           ir11= ir10+ neq         !vel
           ir12= ir11+ neq         !acc
           ir13= ir12+ neq         !uup
           ir14= ir13+ neq         !vvp
           ir15= ir14+ neq         !aap
           call  ckSpace(ilog,isize8,'INCrement_dyn',ir15)
c
!           call zztimer(ilog,'--> INCrement_dyn')
           write(ilog,*)'@@ global ',iglobal
*          write(*,*) maxstiff,maxnode,neq,' xxxx'
           call non_gen_exp( respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
*    &                   respce8(ir8), respce8(ir9), 
     &                                 respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13), respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt)
!           call zztimer(ilog,'<-- INCrement_dyn')
c
        elseif (ianal .eq. 2) then
            goto 4900
        elseif (ianal .eq. 3) then
*           goto 6400
        endif
        goto 1
cccc
c
cccccccccccccccccc rubber & rubberlike material         
 4500   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
*       write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1  + maxstiff    !K_E   
           ir3 = ir2  + maxstiff    !K_G
           ir4 = ir3  + neq         !wk
           ir5 = ir4  + neq         !{P}
           ir6 = ir5  + neq         !{u}
           ir7 = ir6  + neq         !fmag
           ir8 = ir7  + neq         !{u}_t
           ir9 = ir8  + maxstiff    !K_eff
           ir10= ir9  + maxstiff    !mass 
           ir11= ir10 + neq         !vel   
           ir12= ir11 + neq         !acc   
           ir13= ir12 + neq         !uup  
           ir14= ir13 + neq         !vvp  
           ir15= ir14 + neq         !aap   
           call  ckSpace(ilog,isize8,'INCrement',ir15)
c
!           call zztimer(ilog,'--> INCrement')
           write(ilog,*)'@@ global ',iglobal
           call nonrub_3D_TL( respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13),respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt
     &                  ,respceb,isize8)
!           call zztimer(ilog,'<-- INCrement')
c
        elseif (ianal .eq. 2) then
*          goto 4900
c          storage pointers
           ir1 = 1
           ir2 = ir1 + neq         !{u}   
           ir3 = ir2 + nnp*3       ![u]
            ir4 = ir3 + nnp*6       ![n strn]
            ir5 = ir4 + nnp*6       ![n strs]
            ir6 = ir5 + nel*20*3    ![en forc]
            ir7 = ir6 + nel*6       ![e strain]
            ir8 = ir7 + nel*6       ![e stress]
            ir9 = ir8 + nel*6*27    ![strain_ip]
            ir10= ir9 + nel*6*27    ![stress_ip]
            ir11= ir10+ nel*6*27    ![strainp_ip]
            ir12= ir11+ nel*6*20    ![wk_elm_node]
            ir13= ir12+ nel*6*27    ![straini_ip]
            ir14= ir13+ nel*6*27    ![stressi_ip]
            ir15= ir14+ nel*6*27    ![strainpi_ip]
            ir16= ir15+ nel*6*27    ![stressc_ip]
           call  ckSpace(ilog,isize8,'non_postan',ir16)
c
!           call zztimer(ilog,'--> non_POSTan')
           write(ilog,*)'@@ global ',iglobal
           call nonrub_postan( respce8(ir1), respce8(ir2),
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz_D(1),xyz_D(maxnode+1),xyz_D(maxnode*2+1),
     &                   respce8(ir3),respce8(ir4),respce8(ir5),
     &                   respce8(ir6),respce8(ir7),
     &                   respce8(ir8),respce8(ir9),
     &                   respce8(ir10),respce8(ir11),
     &                   respce8(ir12),respce8(ir13),respce8(ir14),
     &                   respce8(ir15),
     &                   prop_D, nmprop,ippp,iglobal)
!           call zztimer(ilog,'<-- non_POSTan')
c
        elseif (ianal .eq. 3) then
*           goto 6400
        endif
        goto 4500
cccc
c       rubber explicit
 4520   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
*       write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1  + maxstiff    !K_E   
           ir3 = ir2  + maxstiff    !K_G
           ir4 = ir3  + neq         !wk
           ir5 = ir4  + neq         !{P}
           ir6 = ir5  + neq         !{u}
           ir7 = ir6  + neq         !fmag
           ir8 = ir7  + neq         !{u}_t
           ir9 = ir8  + maxstiff    !K_eff
           ir10= ir9  + maxstiff    !mass 
           ir11= ir10 + neq         !vel   
           ir12= ir11 + neq         !acc   
           ir13= ir12 + neq         !uup  
           ir14= ir13 + neq         !vvp  
           ir15= ir14 + neq         !aap   
           call  ckSpace(ilog,isize8,'INCrement',ir15)
c
!           call zztimer(ilog,'--> INCrement')
           write(ilog,*)'@@ global ',iglobal
           call nonrub_3D_exp( respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13),respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt
     &                  ,respceb,isize8)
!           call zztimer(ilog,'<-- INCrement')
c
        elseif (ianal .eq. 2) then
*          goto 4900
c          storage pointers
           ir1 = 1
           ir2 = ir1 + neq         !{u}   
           ir3 = ir2 + nnp*3       ![u]
            ir4 = ir3 + nnp*6       ![n strn]
            ir5 = ir4 + nnp*6       ![n strs]
            ir6 = ir5 + nel*20*3    ![en forc]
            ir7 = ir6 + nel*6       ![e strain]
            ir8 = ir7 + nel*6       ![e stress]
            ir9 = ir8 + nel*6*27    ![strain_ip]
            ir10= ir9 + nel*6*27    ![stress_ip]
            ir11= ir10+ nel*6*27    ![strainp_ip]
            ir12= ir11+ nel*6*20    ![wk_elm_node]
            ir13= ir12+ nel*6*27    ![straini_ip]
            ir14= ir13+ nel*6*27    ![stressi_ip]
            ir15= ir14+ nel*6*27    ![strainpi_ip]
            ir16= ir15+ nel*6*27    ![stressc_ip]
           call  ckSpace(ilog,isize8,'non_postan',ir16)
c
!           call zztimer(ilog,'--> non_POSTan')
           write(ilog,*)'@@ global ',iglobal
           call nonrub_postan( respce8(ir1), respce8(ir2),
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz_D(1),xyz_D(maxnode+1),xyz_D(maxnode*2+1),
     &                   respce8(ir3),respce8(ir4),respce8(ir5),
     &                   respce8(ir6),respce8(ir7),
     &                   respce8(ir8),respce8(ir9),
     &                   respce8(ir10),respce8(ir11),
     &                   respce8(ir12),respce8(ir13),respce8(ir14),
     &                   respce8(ir15),
     &                   prop_D, nmprop,ippp,iglobal)
!           call zztimer(ilog,'<-- non_POSTan')
c
        elseif (ianal .eq. 3) then
*           goto 6400
        endif
        goto 4520
cccc
c
cccccccccccccccccc Arruda-Boyce tissue
 4700   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
*       write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1  + maxstiff    !K_E   
           ir3 = ir2  + maxstiff    !K_G
           ir4 = ir3  + neq         !wk
           ir5 = ir4  + neq         !{P}
           ir6 = ir5  + neq         !{u}
           ir7 = ir6  + neq         !fmag
           ir8 = ir7  + neq         !{u}_t
           ir9 = ir8  + maxstiff    !K_eff
           ir10= ir9  + maxstiff    !mass 
           ir11= ir10 + neq         !vel   
           ir12= ir11 + neq         !acc   
           ir13= ir12 + neq         !uup  
           ir14= ir13 + neq         !vvp  
           ir15= ir14 + neq         !aap   
           call  ckSpace(ilog,isize8,'INCrement',ir15)
c
           isize_elem=nel
           ic1 = 1
           ic2 = ic1  + nel*6*27    !stresst
           ic3 = ic2  + nel*6*27    !straint
           ic4 = ic3  + nel*6*27    !stressi
           ic5 = ic4  + nel*6*27    !straini
           ic6 = ic5  + nel*6*27    !stressci
           isize_dof=neq
           ic7 = ic6  + neq         !ui        
           ic8 = ic7  + neq         !dui        
           ic9 = ic8  + neq         !grav      
           ic10= ic9  + neq         !fmag3     
           ic11= ic10 + neq         !dloadi    
           ic12= ic11 + neq         !floadi    
           ic13= ic12 + neq         !floadm    
           ic14= ic13 + maxxyz      !xyzt      
           ic15= ic14 + maxxyz      !xyzi      
c
!           call zztimer(ilog,'--> INCrement')
           write(ilog,*)'@@ global ',iglobal
           call nonabt_3D_TL(
     &                   respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13),respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt
     &                  ,respceb,isize8
     &          ,isize_elem,respcec(ic1)
     &          ,respcec(ic2), respcec(ic3),respcec(ic4), respcec(ic5) 
     &          ,isize_dof,respcec(ic6)
     &          ,respcec(ic7 ),respcec(ic8 ),respcec(ic9) 
     &          ,respcec(ic10),respcec(ic11)
     &          ,respcec(ic12),respcec(ic13),respcec(ic14) 
c    &          ,idload,dloadi,floadi,floadm
     &                     )  
!           call zztimer(ilog,'<-- INCrement')
c
        elseif (ianal .eq. 2) then
c          storage pointers
           ir1 = 1
           ir2 = ir1 + neq         !{u}   
           ir3 = ir2 + nnp*3       ![u]
            ir4 = ir3 + nnp*6       ![n strn]
            ir5 = ir4 + nnp*6       ![n strs]
            ir6 = ir5 + nel*20*3    ![en forc]
            ir7 = ir6 + nel*6       ![e strain]
            ir8 = ir7 + nel*6       ![e stress]
            ir9 = ir8 + nel*6*27    ![strain_ip]
            ir10= ir9 + nel*6*27    ![stress_ip]
            ir11= ir10+ nel*6*27    ![strainp_ip]
            ir12= ir11+ nel*6*20    ![wk_elm_node]
            ir13= ir12+ nel*6*27    ![straini_ip]
            ir14= ir13+ nel*6*27    ![stressi_ip]
            ir15= ir14+ nel*6*27    ![strainpi_ip]
            ir16= ir15+ nel*6*27    ![stressc_ip]
           call  ckSpace(ilog,isize8,'non_postan',ir16)
c
!           call zztimer(ilog,'--> non_POSTan')
           write(ilog,*)'@@ global ',iglobal
           call nonabt_postan( respce8(ir1), respce8(ir2),
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz_D(1),xyz_D(maxnode+1),xyz_D(maxnode*2+1),
     &                   respce8(ir3),respce8(ir4),respce8(ir5),
     &                   respce8(ir6),respce8(ir7),
     &                   respce8(ir8),respce8(ir9),
     &                   respce8(ir10),respce8(ir11),
     &                   respce8(ir12),respce8(ir13),respce8(ir14),
     &                   respce8(ir15),
     &                   prop_D, nmprop,ippp,iglobal)
!           call zztimer(ilog,'<-- non_POSTan')
c
        elseif (ianal .eq. 3) then
        endif
        goto 4700
cccc
c       rubber explicit
 4720   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1  + maxstiff    !K_E   
           ir3 = ir2  + maxstiff    !K_G
           ir4 = ir3  + neq         !wk
           ir5 = ir4  + neq         !{P}
           ir6 = ir5  + neq         !{u}
           ir7 = ir6  + neq         !fmag
           ir8 = ir7  + neq         !{u}_t
           ir9 = ir8  + maxstiff    !K_eff
           ir10= ir9  + maxstiff    !mass 
           ir11= ir10 + neq         !vel   
           ir12= ir11 + neq         !acc   
           ir13= ir12 + neq         !uup  
           ir14= ir13 + neq         !vvp  
           ir15= ir14 + neq         !aap   
           call  ckSpace(ilog,isize8,'INCrement',ir15)
c
           isize_elem=nel
           ic1 = 1
           ic2 = ic1  + nel*6*27    !stresst
           ic3 = ic2  + nel*6*27    !straint
           ic4 = ic3  + nel*6*27    !stressi
           ic5 = ic4  + nel*6*27    !straini
           ic6 = ic5  + nel*6*27    !stressci
           isize_dof=neq
           ic7 = ic6  + neq         !ui        
           ic8 = ic7  + neq         !dui        
           ic9 = ic8  + neq         !grav      
           ic10= ic9  + neq         !fmag3     
           ic11= ic10 + neq         !dloadt    
           ic12= ic11 + neq         !floadt    
           ic13= ic12 + neq         !dut       
           ic14= ic13 + maxxyz      !xyzt      
           ic15= ic14 + maxxyz      !xyzi      
c
c
!           call zztimer(ilog,'--> INCrement')
           write(ilog,*)'@@ global ',iglobal
           call nonabt_3D_exp( respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13),respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt
     &                  ,respceb,isize8
     &          ,isize_elem,respcec(ic1)
     &          ,respcec(ic2), respcec(ic3),respcec(ic4), respcec(ic5) 
     &          ,isize_dof,respcec(ic6)
     &          ,respcec(ic7 ),respcec(ic8 ),respcec(ic9) 
     &          ,respcec(ic10),respcec(ic11)
     &          ,respcec(ic12),respcec(ic13),respcec(ic14) 
     &                     )  
!           call zztimer(ilog,'<-- INCrement')
c
        elseif (ianal .eq. 2) then
c          storage pointers
           ir1 = 1
           ir2 = ir1 + neq         !{u}   
           ir3 = ir2 + nnp*3       ![u]
            ir4 = ir3 + nnp*6       ![n strn]
            ir5 = ir4 + nnp*6       ![n strs]
            ir6 = ir5 + nel*20*3    ![en forc]
            ir7 = ir6 + nel*6       ![e strain]
            ir8 = ir7 + nel*6       ![e stress]
            ir9 = ir8 + nel*6*27    ![strain_ip]
            ir10= ir9 + nel*6*27    ![stress_ip]
            ir11= ir10+ nel*6*27    ![strainp_ip]
            ir12= ir11+ nel*6*20    ![wk_elm_node]
            ir13= ir12+ nel*6*27    ![straini_ip]
            ir14= ir13+ nel*6*27    ![stressi_ip]
            ir15= ir14+ nel*6*27    ![strainpi_ip]
            ir16= ir15+ nel*6*27    ![stressc_ip]
           call  ckSpace(ilog,isize8,'non_postan',ir16)
c
!           call zztimer(ilog,'--> non_POSTan')
           write(ilog,*)'@@ global ',iglobal
           call nonabt_postan( respce8(ir1), respce8(ir2),
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz_D(1),xyz_D(maxnode+1),xyz_D(maxnode*2+1),
     &                   respce8(ir3),respce8(ir4),respce8(ir5),
     &                   respce8(ir6),respce8(ir7),
     &                   respce8(ir8),respce8(ir9),
     &                   respce8(ir10),respce8(ir11),
     &                   respce8(ir12),respce8(ir13),respce8(ir14),
     &                   respce8(ir15),
     &                   prop_D, nmprop,ippp,iglobal)
!           call zztimer(ilog,'<-- non_POSTan')
c
        elseif (ianal .eq. 3) then
        endif
        goto 4720
cccc
c       A-B with visco fibers
 4740   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
*       write(*,*) '         3=post_analysis partial snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1  + maxstiff    !K_E   
           ir3 = ir2  + maxstiff    !K_G
           ir4 = ir3  + neq         !wk
           ir5 = ir4  + neq         !{P}
           ir6 = ir5  + neq         !{u}
           ir7 = ir6  + neq         !fmag
           ir8 = ir7  + neq         !{u}_t
           ir9 = ir8  + maxstiff    !K_eff
           ir10= ir9  + maxstiff    !mass 
           ir11= ir10 + neq         !vel   
           ir12= ir11 + neq         !acc   
           ir13= ir12 + neq         !uup  
           ir14= ir13 + neq         !vvp  
           ir15= ir14 + neq         !aap   
           ir16= ir15 + maxstiff    !dmp  
           ir17= ir16 + neq         !wk2  
           call  ckSpace(ilog,isize8,'INCrement',ir17)
c
           isize_elem=nel
           ic1 = 1
           ic2 = ic1  + nel*6*27    !stresst
           ic3 = ic2  + nel*6*27    !straint
           ic4 = ic3  + nel*6*27    !stressi
           ic5 = ic4  + nel*6*27    !straini
           ic6 = ic5  + nel*6*27    !stressci
           isize_dof=neq
           ic7 = ic6  + neq         !ui        
           ic8 = ic7  + neq         !dui        
           ic9 = ic8  + neq         !grav      
           ic10= ic9  + neq         !fmag3     
           ic11= ic10 + neq         !dloadi    
           ic12= ic11 + neq         !floadi    
           ic13= ic12 + neq         !floadm    
           ic14= ic13 + maxxyz      !xyzt      
           ic15= ic14 + maxxyz      !xyzi      
c
!           call zztimer(ilog,'--> INCrement')
           write(ilog,*)'@@ global ',iglobal
           call nonabt_vis_TL(
     &                   respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13),respce8(ir14) 
     &                  ,respce8(ir15),respce8(ir16)   
     &                  ,idbc,npijkm,maxnode,maxelem
     &                  ,iprofv, iprofh,nloc 
     &                  ,xyz_D ,maxxyz
     &                  ,iglobal, ippp
     &                  ,prop_D, nmprop,nelt
     &                  ,respceb,isize8
     &          ,isize_elem,respcec(ic1)
     &          ,respcec(ic2), respcec(ic3),respcec(ic4), respcec(ic5) 
     &          ,isize_dof,respcec(ic6)
     &          ,respcec(ic7 ),respcec(ic8 ),respcec(ic9) 
     &          ,respcec(ic10),respcec(ic11)
     &          ,respcec(ic12),respcec(ic13),respcec(ic14) 
c    &          ,idload,dloadi,floadi,floadm
     &                     )  
!           call zztimer(ilog,'<-- INCrement')
c
        elseif (ianal .eq. 2) then
c          storage pointers
           ir1 = 1
           ir2 = ir1 + neq         !{u}   
           ir3 = ir2 + nnp*3       ![u]
            ir4 = ir3 + nnp*6       ![n strn]
            ir5 = ir4 + nnp*6       ![n strs]
            ir6 = ir5 + nel*20*3    ![en forc]
            ir7 = ir6 + nel*6       ![e strain]
            ir8 = ir7 + nel*6       ![e stress]
            ir9 = ir8 + nel*6*27    ![strain_ip]
            ir10= ir9 + nel*6*27    ![stress_ip]
            ir11= ir10+ nel*6*27    ![strainp_ip]
            ir12= ir11+ nel*6*20    ![wk_elm_node]
            ir13= ir12+ nel*6*27    ![straini_ip]
            ir14= ir13+ nel*6*27    ![stressi_ip]
            ir15= ir14+ nel*6*27    ![strainpi_ip]
            ir16= ir15+ nel*6*27    ![stressc_ip]
           call  ckSpace(ilog,isize8,'non_postan',ir16)
c
!           call zztimer(ilog,'--> non_POSTan')
           write(ilog,*)'@@ global ',iglobal
           call nonabt_postan( respce8(ir1), respce8(ir2),
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz_D(1),xyz_D(maxnode+1),xyz_D(maxnode*2+1),
     &                   respce8(ir3),respce8(ir4),respce8(ir5),
     &                   respce8(ir6),respce8(ir7),
     &                   respce8(ir8),respce8(ir9),
     &                   respce8(ir10),respce8(ir11),
     &                   respce8(ir12),respce8(ir13),respce8(ir14),
     &                   respce8(ir15),
     &                   prop_D, nmprop,ippp,iglobal)
!           call zztimer(ilog,'<-- non_POSTan')
c
        elseif (ianal .eq. 3) then
        endif
        goto 4740
cccc
cccc
c       visco rubber explicit
 4760   continue
c
        write(*,'(a)')' @@'
        write(*,*) 'CHOOSE:  0=return  '
        write(*,*) '         1=analysis  '
        write(*,*) '         2=post_analysis full snapshot'
        call zzwrt(' --> ')
        read(ikbd,*) ianal
        write(ilog,*) ianal,'  :: 1=analysis'
c
        ibandm=iband
        if (ianal .eq. 0) then
            goto 1
        elseif (ianal .eq. 1) then
c          storage pointers
           ir1 = 1
           ir2 = ir1  + maxstiff    !K_E   
           ir3 = ir2  + maxstiff    !K_G
           ir4 = ir3  + neq         !wk
           ir5 = ir4  + neq         !{P}
           ir6 = ir5  + neq         !{u}
           ir7 = ir6  + neq         !fmag
           ir8 = ir7  + neq         !{u}_t
           ir9 = ir8  + maxstiff    !K_eff
           ir10= ir9  + maxstiff    !mass 
           ir11= ir10 + neq         !vel   
           ir12= ir11 + neq         !acc   
           ir13= ir12 + neq         !uup  
           ir14= ir13 + neq         !vvp  
           ir15= ir14 + neq         !aap   
           call  ckSpace(ilog,isize8,'INCrement',ir15)
c
           isize_elem=nel
           ic1 = 1
           ic2 = ic1  + nel*6*27    !stresst
           ic3 = ic2  + nel*6*27    !straint
           ic4 = ic3  + nel*6*27    !stressi
           ic5 = ic4  + nel*6*27    !straini
           ic6 = ic5  + nel*6*27    !stressci
           isize_dof=neq
           ic7 = ic6  + neq         !ui        
           ic8 = ic7  + neq         !dui        
           ic9 = ic8  + neq         !grav      
           ic10= ic9  + neq         !fmag3     
           ic11= ic10 + neq         !dloadt    
           ic12= ic11 + neq         !floadt    
           ic13= ic12 + neq         !dut       
           ic14= ic13 + maxxyz      !xyzt      
           ic15= ic14 + maxxyz      !xyzi      
c
c
!           call zztimer(ilog,'--> INCrement')
           write(ilog,*)'@@ global ',iglobal
           call nonabt_vis_exp( respce8(ir1), respce8(ir2),
     &                   respce8(ir3),
     &                   respce8(ir4), respce8(ir5), 
     &                   respce8(ir6), respce8(ir7),
     &                   respce8(ir8), respce8(ir9), 
     &                   respce8(ir10), respce8(ir11), 
     &                   respce8(ir12), respce8(ir13),respce8(ir14), 
     &                   idbc,npijkm,maxnode,maxelem, 
     &                   iprofv, iprofh,nloc ,
     &                   xyz_D ,maxxyz,
     &                   iglobal, ippp,
     &                   prop_D, nmprop,nelt
     &                  ,respceb,isize8
     &          ,isize_elem,respcec(ic1)
     &          ,respcec(ic2), respcec(ic3),respcec(ic4), respcec(ic5) 
     &          ,isize_dof,respcec(ic6)
     &          ,respcec(ic7 ),respcec(ic8 ),respcec(ic9) 
     &          ,respcec(ic10),respcec(ic11)
     &          ,respcec(ic12),respcec(ic13),respcec(ic14) 
     &                     )  
!           call zztimer(ilog,'<-- INCrement')
c
        elseif (ianal .eq. 2) then
c          storage pointers
           ir1 = 1
           ir2 = ir1 + neq         !{u}   
           ir3 = ir2 + nnp*3       ![u]
            ir4 = ir3 + nnp*6       ![n strn]
            ir5 = ir4 + nnp*6       ![n strs]
            ir6 = ir5 + nel*20*3    ![en forc]
            ir7 = ir6 + nel*6       ![e strain]
            ir8 = ir7 + nel*6       ![e stress]
            ir9 = ir8 + nel*6*27    ![strain_ip]
            ir10= ir9 + nel*6*27    ![stress_ip]
            ir11= ir10+ nel*6*27    ![strainp_ip]
            ir12= ir11+ nel*6*20    ![wk_elm_node]
            ir13= ir12+ nel*6*27    ![straini_ip]
            ir14= ir13+ nel*6*27    ![stressi_ip]
            ir15= ir14+ nel*6*27    ![strainpi_ip]
            ir16= ir15+ nel*6*27    ![stressc_ip]
           call  ckSpace(ilog,isize8,'non_postan',ir16)
c
!           call zztimer(ilog,'--> non_POSTan')
           write(ilog,*)'@@ global ',iglobal
           call nonabt_postan( respce8(ir1), respce8(ir2),
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz_D(1),xyz_D(maxnode+1),xyz_D(maxnode*2+1),
     &                   respce8(ir3),respce8(ir4),respce8(ir5),
     &                   respce8(ir6),respce8(ir7),
     &                   respce8(ir8),respce8(ir9),
     &                   respce8(ir10),respce8(ir11),
     &                   respce8(ir12),respce8(ir13),respce8(ir14),
     &                   respce8(ir15),
     &                   prop_D, nmprop,ippp,iglobal)
!           call zztimer(ilog,'<-- non_POSTan')
c
        elseif (ianal .eq. 3) then
        endif
        goto 4760
cccc
cccc
c
c       post process 
 4900   continue
c          storage pointers
           ir1 = 1
           ir2 = ir1 + neq         !{u}   
           ir3 = ir2 + nnp*3       ![u]
            ir4 = ir3 + nnp*6       ![n strn]
            ir5 = ir4 + nnp*6       ![n strs]
            ir6 = ir5 + nel*20*3    ![en forc]
            ir7 = ir6 + nel*6       ![e strain]
            ir8 = ir7 + nel*6       ![e stress]
            ir9 = ir8 + nel*6*27    ![strain_ip]
            ir10= ir9 + nel*6*27    ![stress_ip]
            ir11= ir10+ nel*6*27    ![stressc_ip]
            ir12= ir11+ nel*6*20    ![wk_elm_node]
           call  ckSpace(ilog,isize8,'non_postan',ir12)
c
!           call zztimer(ilog,'--> non_POSTan')
           write(ilog,*)'@@ global ',iglobal
           call non_postan( respce8(ir1), respce8(ir2),
     &                   idbc,npijkm,maxnode,maxelem,
     &                   xyz_D(1),xyz_D(maxnode+1),xyz_D(maxnode*2+1),
     &                   respce8(ir3),respce8(ir4),respce8(ir5),
     &                   respce8(ir6),respce8(ir7),
     &                   respce8(ir8),respce8(ir9),
     &                   respce8(ir10),
     &                   respce8(ir11),
     &                   prop_D, nmprop,ippp,iglobal)
!           call zztimer(ilog,'<-- non_POSTan')
c
        goto 1
c
cccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       POSTSCRIPT contour files
 7700   continue
c
        write(*,*)'SAVE new file as:  '
        write(*,*)'     0=return  1=StaDyn.PS  2=NEW Name'
        call zzwrt(' SELECT --> ')
        read(ikbd,*) iname
        write(ilog,*) iname,'   ::iname'
c
        if (iname .eq. 0) then
            goto 1
        elseif (iname .eq. 1) then
            filena='stadyn.ps'
            open(ioutps,file=fdn//adjustl(filena))
            rewind ioutps
        elseif (iname .eq. 2) then
            call zzwrt(' TYPE: new filename-->  ')
            read (ikbd,'(a40)') filena
            write(ilog,*) filena
            open(ioutps,file=fdn//adjustl(filena))
            rewind ioutps
        endif
c
c
           write(*,*) '@@ entering PS_contour'
!           call zztimer(ilog,'--> CONTour')
           call ps_contour(ioutps, nmprop,
     &                     respce8(1),isize8, irespce(1),isize4
     &                     )
!           call zztimer(ilog,'<-- CONTour')
c
        close (ioutps) 
        goto 1
c
 999    continue
        write(ilog,*)'@@ Simplex OK, ended from MAIN'
c
      end
c
c
c     INITIALizations for memory requirements
      subroutine initial (ilog,isize8, fdn,
     &                    maxnode,maxelem,maxforce,maxpts,ilump,
     &                    itermax,rtol)
         implicit real*8 (a-h,o-z)
         character*2 :: fdn

         logical exists
         isize8  =150000000 
         maxnode = 15000
         maxelem = 6000
         maxforce= 10000
         maxpts=10000
         ilump=1
         itermax=16
         rtol=1.0D-10
*        r4=4
         print*, fdn//'simplex.cfg'
          inquire(file=fdn//'simplex.cfg',exist=exists)
          if ( .not. exists) then
               write(ilog,*)'@@ !!! <<simplex.cfg>> does NOT exist, ',
     &                                             ' will create'
               write(*   ,*)'@@ !!! <<simplex.cfg>> does NOT exist, ',
     &                                             ' will create'
               call wrt_cfg(
     &           isize8,maxelem,maxnode,maxforce,maxpts,rtol )
           endif
           call rd_cfg(
     &           isize8,maxelem,maxnode,maxforce,maxpts,rtol )
c
          write(ilog,*)'@@ MAXimum storage   : ',isize8 
          write(ilog,*)'@@ MAXimum elements  : ',maxelem
          write(ilog,*)'@@ MAXimum nodes     : ',maxnode
          write(ilog,*)'@@ MAXimum force incs: ',maxforce,maxpts
          write(ilog,*)'@@ ITERmax           : ',itermax        
          write(ilog,*)'@@ rtol              : ',rtol         
c
      return
      end
c
c
      subroutine wrt_cfg(
     &           isize8,maxelem,maxnode,maxforce,maxpts,rtol )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         open(unit=icfg,file=fdn//'simplex.cfg')
         rewind icfg
         write(icfg,*) isize8,'     ::MaxStorage'
         write(icfg,*) maxelem,'     ::MaxElem'
         write(icfg,*) maxnode,'     ::MaxNode'
         write(icfg,*) maxforce,maxpts,' ::MaxForce pts'
         write(icfg,*) ilump,'      ::ilump'
         write(icfg,*) itermax,'      ::itermax'
         write(icfg,*) rtol ,'      ::rtol '
!         call timday(icfg)
         close (icfg)
c
      return
      end
c
c
      subroutine rd_cfg(
     &           isize8,maxelem,maxnode,maxforce,maxpts,rtol )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
          ierror=1
          open(unit=icfg,file=fdn//'simplex.cfg')
          rewind icfg
          read(icfg,*,end=98,err=98) isize8
          read(icfg,*,end=98,err=98) maxelem
          read(icfg,*,end=98,err=98) maxnode
          read(icfg,*,end=98,err=98) maxforce,maxpts
          read(icfg,*,end=98,err=98) ilump
          read(icfg,*,end=98,err=98) itermax
          read(icfg,*,end=98,err=98) rtol 
          ierror=0
 98       continue
          close(icfg)
          if (ierror .eq. 1) then
               write(ilog,*)'@@ !!! error opening <<simplex.CFG>>'
               write(*   ,*)'@@ !!! error opening <<simplex.CFG>>'
               call wrt_cfg(
     &           isize8,maxelem,maxnode,maxforce,maxpts,rtol )
          endif
c
      return
      end
c
c
c     READ structural DaTafile
      subroutine readdt(ilog,
     >             xyz,
     >             respce8,isize8,
     >             idbc, maxdof, maxelem, maxnode,
     &             npijkm, iglobal, iprofv,iprofh,
     &             prop,nmprop,ippp)
c
         implicit real*8 (a-h,o-z)
          integer idbc(maxnode,3)
          integer npijkm(maxelem,21),
     &            iprofv(maxdof),iprofh(maxdof)
          integer nmprop(maxelem)
          real*8  prop(ippp,10)
c
          real*8  xyz(maxnode*3)
          real*8  respce8(isize8)
c
c       storage pointers
        ir1 = 1
        ir2 = ir1 + maxnode     !xload   
        ir3 = ir2 + maxnode     !yload
        ir4 = ir3 + maxnode     !xload
        ir5 = ir4 + maxnode     !cmass
        ir6 = ir5 + maxnode*3   !heat
        ir7 = ir6 + maxnode*3   !wk
        call ChckSpce(ilog,isize8,respce8,'DATAin',ir7)
c
!        call zztimer(ilog,'--> DATAin')
!        call zzflush(ilog)


        call datain( 
     >             xyz(1),xyz(maxnode+1),xyz(maxnode*2+1),
     &             respce8(ir1),respce8(ir2),
     >             respce8(ir3),respce8(ir4),
     >             respce8(ir5),
     >             idbc, maxdof, maxelem, maxnode,
     &             npijkm, iglobal,
     &             prop,nmprop,ippp)

!        call zzflush(ilog)
c
         write(*,*)'@@ now converting'
         call datacnv( 
     >             xyz(1),xyz(maxnode+1),xyz(maxnode*2+1),
     &             respce8(ir1),respce8(ir2),
     >             respce8(ir3),respce8(ir4),respce8(ir5),
     >             respce8(ir6),
     >             idbc, maxdof, maxelem, maxnode,
     &             npijkm, iglobal,
     >             iprofv, iprofh,
     &             prop,nmprop,ippp)
!        call zztimer(ilog,'<-- DATAin')
!        call zzflush(ilog)
c
      return
      end
c
c
c     ADD GRAVITY self weight to the load vector.
      subroutine add_gravity( load,bforce, idbc,maxnode)
c
         implicit real*8 (a-h,o-z)
                 include 'commons.std'
c
         integer idbc(maxnode,3)
         real*8 bforce(neq) ,load(neq )
c
*        if (abs(gx) .gt. 1e-10 .OR. abs(gy) .gt. 1.0e-10
*    &                          .OR. abs(gz) .gt. 1.0e-10) then
*            goto 100
*        else
*            return
*        endif
c
 100     continue
c
c
          write(*,*)'@@ reading <<stadyn.lod<<'
          rewind (ilod)
          do i=1,neq
             read(ilod) load(i)
          enddo
ctag1
c           Fill in  the displacement vector
            do 70 i = 1, nnp
               do j=1,3
                  ieqnum = idbc(i,j)
                  if (ieqnum .gt. 0) then
                      idof=j
                      load(ieqnum) = load(ieqnum) 
     &                             + bforce(ieqnum)
                  endif
               enddo
 70         continue 
c
c
          rewind (ilod)
          do i=1,neq
             write(ilod) load(i)
             write(ilog,*) '@@grav load ',load(i),i
          enddo
c
      return
      end
c
c
