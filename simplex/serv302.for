c
c     PostScript CONTOURs of results
      subroutine ps_contour(ioutps,nmprop,respce8,isize8,irespce,isize4)
         implicit real*8 (a-h,o-z)
         include 'commons.std'
         integer nmprop(nel)
         real*8  xymov(100,2)
         character*50 str50
         real*8   respce8(isize8)
         integer  irespce(isize4)
c
         pi=4.0*atan(1.0)
c
         do i=1,100
            xymov(i,1)=0
            xymov(i,2)=0
         enddo
         write(ilog,*)'@@ in PS_Contour'
         str50=' StaDyn '
         str50=' Simplex '
c
         call ct_initial( )
         mode=0
         ivars=1
         factor=1
c
c
 1      continue
         write(*,*) ' '
         write(*,*) '   CONTOURs: '
         write(*,*) '          0:  Quit '
*        write(*,*) ' '
         write(*,*) '                RENDER PS file'
         write(*,*) '          4:  Deformed shape'
         write(*,*) '          5:  Contours  '
*        write(*,*) ' '
         write(*,*)'                 CHANGE defaults:  '
         write(*,*)'          11:  Position model '
         write(*,*)'          12:  Rotate model             '
         write(*,*)'          13:  Move tagged surfaces     '
         write(*,*)'          14:  Add caption      '
         write(*,*)'          15:  Toggle tags      '
         write(*,*)'          16:  Change pens      '
         write(*,*)'          17:  Background       '
         call zzwrt(' SELECT --> ')
         read(ikbd,*) icntr
         write(ilog,*) icntr,'     ::contour'
         write(*,*) ' '
c
c
         if (icntr .eq. 0) then
             call ct_cfg( )
             goto 999
         endif
         if (icntr .eq. 4) goto 400
         if (icntr .eq. 5) goto 500
         if (icntr .eq. 11) goto 8510
         if (icntr .eq. 12) goto 8520
         if (icntr .eq. 13) goto 8530
         if (icntr .eq. 14) goto 8540
         if (icntr .eq. 15) goto 8550
         if (icntr .eq. 16) goto 8560
         if (icntr .eq. 17) goto 8570
c
            goto 1
c
c
c        POSITION
 8510    continue
         write(*,*)' '
         write(*,*)'INPUT:     Xposn   |   Yposn   '
         write(*,84)'current: ',pxmove,pymove
         call zzwrt('  --> ')
         read(ikbd,*,err=8510) pxmove,pymove
         write(ilog,'(1x,3(g12.6,1x),a)') pxmove,pymove,'::xy '
         goto 1
c
c        ROTATION
 8520    continue
         write(*,*)' '
         write(*,*)'INPUT:     RotX   |   RotY   |   RotZ  '
         write(*,*)' suggest     20   |    45    |    0    '
         write(*,84)'current: ',rotx,roty,rotz
         call zzwrt('  --> ')
         read(ikbd,*,err=8520) rotx,roty,rotz
         write(ilog,'(1x,3(g12.6,1x),a)') rotx,roty,rotz,'::xyz'
         goto 1
c
c        MOVE TAGs
 8530    continue
         write(*,*)' '
 8532    continue
         write(*,*)'INPUT:  tag#  |  Xmove  |  Ymove      <0 0 0 ends>'
         call zzwrt('  --> ')
         read(ikbd,*,err=8532) itag,xmov,ymov
         write(ilog,'(1x,i5,1x,2(g12.6,1x),a)') itag,xmov,ymov,'::xy ms'
         if (itag .ne. 0) then
             xymov(itag,1)=xmov
             xymov(itag,2)=ymov
             goto 8532
         endif
         goto 1
c
c        CAPTION
 8540    continue
         write(*,'(a,a50)')'current: ',str50  
         call zzwrt(' TYPE: caption--> ')
         read(ikbd,'(a50)' ) str50 
         write(ilog,'(1x,a)') str50 
         write(*,*)'INPUT:     Xposn   |   Yposn   |   size  '
         write(*,84)'current: ',xcap,ycap,scap 
         call zzwrt(' --> ')
         read(ikbd,*,err=8540) xcap,ycap,scap
         write(ilog,'(1x,3(g12.6,1x),a)') xcap,ycap,scap,' ::xy cap'
          goto 1
c
c        SUB-structures
 8550    continue
         write(*,*)'INPUT subs:   1 | 2 | 3 | 4 | 5 | 6+   <1=on 0=off>'
         write(*,85)'current: ',isub1,isub2,isub3,isub4,isub5,isub6
         call zzwrt(' --> ')
         read(ikbd,*,err=8550) isub1,isub2,isub3,isub4,isub5,isub6
         write(ilog,'(1x,6(i3,1x),a)') isub1,isub2,isub3
     &                                ,isub4,isub5,isub6,' ::subs'
          goto 1
c
c        change PENs 
 8560    continue
         write(*,*)'INPUT pens:   mesh  |  contour    '
         write(*,84)'current: ',p_thick,c_thick 
         call zzwrt(' --> ')
         read(ikbd,*,err=8560) p_thick,c_thick 
         write(ilog,*) p_thick,c_thick,'    ::pens'
          goto 1
c
c        change BACKground
 8570    continue
         write(*,*)'INPUT   R | G |  B      <0.0 to 1.0>'
         write(*,84)'current: ',r_back,g_back,b_back
         call zzwrt(' --> ')
         read(ikbd,*,err=8560) r_back,g_back,b_back
        write(ilog,'(1x,3(g12.6,1x),a)') r_back,g_back,b_back,' ::RGB'
          goto 1
c
c
c        DEFORMed shape
 400     continue
*        write(*,*)'              <M:1=Y 0=N -1=BLK> <L:1=Y 0=N -1=BLK>'
         write(*,*)'CHOOSE:    '
        write(*,*) 'Shape:     7=vector 8=tri 18=curve'
        write(*,*) 'Mesh:      1=Yes    0=No  -1=Black'
        write(*,*) 'Legend:    1=Yes    0=No  -1=Black'
        write(*,*) ' '
         write(*,*)'      shape | scale | mesh? | legend?'
           call zzwrt(' --> ')
           read(ikbd,*,err=1) ivars,factor,ishow_m,ishow_l
           write(ilog,*) ivars,factor,ishow_m,ishow_l,' ::vars X m l '
           write(*,*) ' '
           goto 1000
c
 500    continue
*       write(*,*)'              <M:1=Y 0=N -1=BLK> <L:1=Y 0=N -1=BLK>'
        write(*,*) 'CHOOSE:    '
        write(*,*) 'DoF:       1=u    2=v    3=w                       '
        write(*,*) '        or 1=Exx  2=Eyy  3=Ezz  4=Exz  5=Eyz  6=Exy'
        write(*,*) '        or 1=Sxx  2=Syy  3=Szz  4=Sxz  5=Syz  6=Sxy'
        write(*,*) 'Mesh:      1=Yes 0=No -1=Black'
        write(*,*) 'Legend:    1=Yes 0=No -1=Black'
        write(*,*) ' '
        write(*,*) 'INPUT:  DoF  | mesh? | legend? '
        call zzwrt(' --> ')
           read(ikbd,*,err=500) ivars,ishow_m,ishow_l
           write(ilog,*        ) ivars,ishow_m,ishow_l,' ::dof m l '
           write(*,*) ' '
           goto 1000
c
c
 1000    continue
*            write(ilog,*)'@@ rots x y z  ',rotx,roty,rotz
c
             call head(ioutps,r_back,g_back,b_back)
             xmove=pxmove*10000.0/10.0
             ymove=pymove*10000.0/10.0
             write(ioutps,'(1x,g12.6,a,g12.6,a)')
     &                        xmove,' ',ymove,' translate'
c
c            draw mesh
                 write(ioutps,'(a)')'%%%%   begin mesh   %%%%%%%'
                 call head_m(ioutps)
                 if (ishow_m .eq. -1) then
                     write(ioutps,'(a)') ' 0.0 0.0 0.0 setrgbcolor'
                 else
                     write(ioutps,'(a)') ' 0.6 0.6 0.6 setrgbcolor'
                 endif
c
                 write(ioutps,'(a)') '[70  0] 0 setdash'
                 write(ioutps,'(1x,g11.5,1x,a)') p_thick,' setlinewidth'
                 lseg=isize4/3
                 mode=0
                 if (ivars .le.  8) then
                     call ps_drawmesh(respce8,lseg,mode,nmprop,ioutps
     &                               ,xymov,ishow_m )
                 elseif (ivars .eq.  9) then
                     call ps_drawmesh_HEX8(respce8,lseg,mode,nmprop
     &                                   ,ioutps,xymov,ishow_m )
                 elseif (ivars .eq. 18) then
*                    call ps_drawmesh_HEX20(respce8,lseg,mode,nmprop
*    &                                   ,ioutps,xymov,ishow_m )
                 endif
                 write(ioutps,'(a)')'%%%%   end mesh   %%%%%%%'
                 call foot_m(ioutps)
c
c            draw contour
                 write(ioutps,'(a)')'%%%%   begin contour   %%%%%%%'
                 call head_m(ioutps)
                 write(ioutps,'(1x,g11.5,1x,a)') c_thick,' setlinewidth'
                 if (ivars .lt.  9) then
                     nmax=isize4/2
                     call ps_stat_cont(respce8,irespce,nmax,nmprop
     &                          ,factor,ivars,mode,ioutps,xymov,ishow_l)
                 elseif (ivars .eq.  9) then
                     nmax=isize4/2
                     call ps_shape_HEX8(respce8,irespce,nmax,nmprop
     &                          ,factor,ivars,mode,ioutps,xymov,ishow_l)
                 elseif (ivars .eq. 18 .OR. ivars .eq. 181) then
                     mdiv=10
                     if (ivars .eq. 181) mdiv=2
                     write(ilog,*)'@@ mdiv: ',mdiv
                     nmax=isize4/2
c                    scale but no plot
                     ishow_xx=0
                     call ps_shape_HEX20(respce8,irespce,nmax,nmprop
     &                        ,factor,ivars,mode,ioutps,xymov,ishow_xx
     &                        ,mdiv)
                 write(ioutps,'(a)') ' 0.6 0.6 0.6 setrgbcolor'
                 write(ioutps,'(a)') '[70  0] 0 setdash'
                 write(ioutps,'(1x,g11.5,1x,a)') p_thick,' setlinewidth'
                     call ps_drawmesh_HEX20(respce8,lseg,mode,nmprop
     &                                   ,ioutps,xymov,ishow_m
     &                         ,mdiv )
                     write(ioutps,*) ' stroke '
                 write(ioutps,'(a)')'%%%%   end hex20 mesh   %%%%%%%'
                 write(ioutps,'(1x,g11.5,1x,a)') c_thick,' setlinewidth'
                     call ps_shape_HEX20(respce8,irespce,nmax,nmprop
     &                          ,factor,ivars,mode,ioutps,xymov,ishow_l
     &                         ,mdiv)
                 endif
                 write(ioutps,'(a)')'%%%%   end contour   %%%%%%%'
                 call foot_m(ioutps)
c
c         CAPtion
          if (ishow_l .eq. 1) then
              iicap=1000*scap
              write(ioutps,'(a,i6,a)') '/Times-Italic findfont ',iicap,
     &                          ' scalefont setfont'
              write(ioutps,'(a)') ' 0 setgray '
              xc=xcap*10000/7.0
              yc=ycap*10000/7.0
              ixc=xc
              iyc=yc
              write(ioutps,82) ixc,iyc,' mt'
              write(ioutps,'(a,a,a)') '(',str50,') show'
          endif
 82       format(1x,2(i6,1x),a)
c
         call foot(ioutps)
         goto 1
c
 84      format(1x,a,6(g12.6,1x))
 85      format(1x,a,6(i3,1x))
c
 999     continue
      return
      end
c
c
c     INITIALizations
      subroutine ct_initial( )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
         logical exists
c
         ngeom=1
         nprob=1
                ilump=1
             rotx=30
             roty=25
             rotz= 0
                p_thick=40
                c_thick=80
             pxmove=1.0
             pymove=1.0
         xcap=4.0
         ycap=7.0
         scap=1.0
           r_back=1.0
           g_back=1.0
           b_back=1.0
c
                imode=1
c
         isub1=1
         isub2=1
         isub3=1
         isub4=1
         isub5=0
         isub6=0
c
         inquire(file='simplex.ctr',exist=exists)
         if ( .not. exists) then
             write(ilog,*)'@@ !!! <<Simplex.CTR>> does NOT exist, ',
     &                                             ' will create'
c            write(idial  ,*)'@@ !!! <<ctour.cfg>> does NOT exist, ',
c    &                                             ' will create'
             call ct_cfg( )
         endif
c
         call ct_rdcfg( )
c
c
       return
       end
c
       subroutine ct_cfg( )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
             open(unit=icfg,file=fdn//'simplex.ctr')
             rewind icfg
          write (icfg,*) ngeom,nprob ,'  ::geom '
          write (icfg,82) pxmove,pymove         ,'  ::pos x y  '
          write (icfg,82) p_thick,c_thick       ,'  ::pen_thick m,c '
          write (icfg,83) rotx,roty,rotz        ,'  ::rotX Y Z '
          write (icfg,83) xcap,ycap,scap        ,'  ::cap X Y size'
          write (icfg,*) max_mode,ilump        ,'  ::modes mass'
          write (icfg,81) isub1,isub2,isub3,isub4,isub5,isub6,' ::subs'
          write (icfg,83) r_back,g_back,b_back  ,'  ::RGB back'
 81              format(1x,6(i1,1x),a)
 82              format(1x,2(g12.6,1x),a)
 83              format(1x,3(g12.6,1x),a)
                 write (icfg,*) 'end '
               close (icfg)
      return
      end
c
c
      subroutine ct_rdcfg( )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         ierror=1
         open(icfg,file=fdn//'simplex.ctr')
         rewind icfg
         read(icfg,*    ,end=98,err=98) ngeom,nprob 
         read(icfg,*    ,end=98,err=98) pxmove,pymove    
         read(icfg,*    ,end=98,err=98) p_thick,c_thick
         read(icfg,*    ,end=98,err=98) rotx,roty,rotz   
         read(icfg,*    ,end=98,err=98) xcap,ycap,scap   
         read(icfg,*,end=98,err=98) max_mode,ilump
         read(icfg,*,end=98,err=98) isub1,isub2,isub3,isub4,isub5,isub6
         read(icfg,*    ,end=98,err=98) r_back,g_back,b_back   
         ierror=0
 98      continue
         close(icfg)
         if (ierror .eq. 1) then
             write(ilog,*)'@@ !!! error opening <<simplex.CTR>>'
c            write(idial   ,*)'@@ !!! error opening <<ctour.cfg>>'
             call ct_cfg( )
         endif
c
      return
      end
c
c
c     DRAW the structure MESH as gray outline
      subroutine ps_drawmesh(x,lseg,mode,nmprop,ioutps,xymov,ishow_m )
c
         implicit real*8 (a-h,o-z)
                 include 'commons.std'
c
         real*8   x(lseg,3)
         real*8   rr(3,3),xx(3),yy(3),zz(3)
         real*8   xymov(100,2)
         integer  nmprop(nel)
         icol=3
c
*           write(ilog,*)'@@ in draw rots x y z  ',rotx,roty,rotz
            call gn_rot3d_mat(rr)
c
c        BIG D0-Loop over all triangles
c        get info from Simplex.OUT file
         rewind(iout)
             read(iout,*) nel3,eigv
             write(ilog,*)'@@nel3 eigv', nel3,eigv
c
c        have all elements participate in scaling
         ielement=0
         ii=0
         do 100 i =  1,921111
            read(iout ,*,end=99) xx(1),yy(1),zz(1)
            read(iout ,*)        xx(2),yy(2),zz(2)
            read(iout ,*)        xx(3),yy(3),zz(3)
*           ielement=ielement+1
*           ielement=(i-1)/4+1
            ielement=(i-1)/6+1
*           isub=nmprop(ielement)
            isub=1
            xxm=xymov(isub,1)
            yym=xymov(isub,2)
            ii=ii+1
c
            call gn_rot3d(rr,xx,yy,zz)
c
                       ir1 = (ii-1)*3 + 1
                       ir2 = (ii-1)*3 + 2
                       ir3 = (ii-1)*3 + 3
              x(ir1,1) = xx(1) + xxm
              x(ir1,2) = yy(1) + yym
              x(ir1,3) = 0.0  
               x(ir2,1) = xx(2) + xxm
               x(ir2,2) = yy(2) + yym
               x(ir2,3) = 0.0  
                x(ir3,1) = xx(3) + xxm  
                x(ir3,2) = yy(3) + yym
                x(ir3,3) = 0.0  
 100            continue
 99   continue
                ntri =ii
                ntri3=ii*3
                iscale=2
                call ttscaledat_3(x,lseg,ntri3,iscale,ilog)
           if (ishow_m .eq. 0) return
c
ctag1
c               only plot on subs
                nstroke=0
                do 200 n=1,ntri
**!!               isub=nmprop(n)
                   isub=1
                   if (isub1 .eq. 0 .AND. isub .eq. 1) goto 200 
                   if (isub2 .eq. 0 .AND. isub .eq. 2) goto 200 
                   if (isub3 .eq. 0 .AND. isub .eq. 3) goto 200 
                   if (isub4 .eq. 0 .AND. isub .eq. 4) goto 200 
                   if (isub5 .eq. 0 .AND. isub .eq. 5) goto 200 
                   if (isub6 .eq. 0 .AND. isub .ge. 6) goto 200 
c
                   nn=(n-1)*3
                   ix1=x(nn+1,1)
                   iy1=x(nn+1,2)
                     ix2=x(nn+2,1)
                     iy2=x(nn+2,2)
                       ix3=x(nn+3,1)
                       iy3=x(nn+3,2)
                   call glnseg_m(ix1,iy1,ix2,iy2,linet,ipen,ioutps)
                   call glnseg_m(ix2,iy2,ix3,iy3,linet,ipen,ioutps)
                   call glnseg_m(ix3,iy3,ix1,iy1,linet,ipen,ioutps)
                   nstroke=nstroke+3
                   if (nstroke .gt. 500) then
                        nstroke=0
                        write(ioutps,*) ' stroke '
                   endif
c
 200            continue
c
      return
      end
c
c
c     DRAW the structure MESH as gray outline
      subroutine ps_drawmesh_HEX20(x,lseg,mode,nmprop
     &                          ,ioutps,xymov,ishow_m,mdiv )
c
         implicit real*8 (a-h,o-z)
                 include 'commons.std'
c
         real*8   x(lseg,3)
         real*8   rr(3,3),xx(3),yy(3),zz(3)
         real*8   xymov(100,2)
         integer  nmprop(nel)
         icol=3
c
*           write(ilog,*)'@@ in draw rots x y z  ',rotx,roty,rotz
            call gn_rot3d_mat(rr)
c
c        BIG D0-Loop over all triangles
c        get info from Simplex.OUT file
         write(*,*) iout
         rewind(iout)
             read(iout,*) nel3,eigv
             write(ilog,*)'@@ nel3 eigv', nel3,eigv
c
c        have all elements participate in scaling
         ielement=0
         ii=0
         do 100 i =  1,921111
            read(iout ,*,end=99) xx(1),yy(1),zz(1)
            read(iout ,*,end=99)        xx(2),yy(2),zz(2)
            read(iout ,*,end=99)        xx(3),yy(3),zz(3)
*           ielement=ielement+1
*           ielement=(i-1)/4+1
            ielement=(i-1)/6+1
*           isub=nmprop(ielement)
            isub=1
            xxm=xymov(isub,1)
            yym=xymov(isub,2)
c
            call gn_rot3d(rr,xx,yy,zz)
c
                       ii=ii+1
                       ir1 = (ii-1)*3 + 1
                       ir2 = (ii-1)*3 + 2
                       ir3 = (ii-1)*3 + 3
*          write(*,*) ir1,ir2,ir3,lseg
              x(ir1,1) = xx(1) + xxm
              x(ir1,2) = yy(1) + yym
              x(ir1,3) = 0.0  
               x(ir2,1) = xx(2) + xxm
               x(ir2,2) = yy(2) + yym
               x(ir2,3) = 0.0  
                x(ir3,1) = xx(3) + xxm  
                x(ir3,2) = yy(3) + yym
                x(ir3,3) = 0.0  
 100     continue
 99   continue
                ntri =ii
                ntri3=ii*3
                iscale=2
                iscale=0
                call ttscaledat_3(x,lseg,ntri3,iscale,ilog)
           if (ishow_m .eq. 0) return
c
ctag1
c               only plot on subs
                nstroke=0
                do 200 n=1,ntri
**!!               isub=nmprop(n)
                   isub=1
                   if (isub1 .eq. 0 .AND. isub .eq. 1) goto 200 
                   if (isub2 .eq. 0 .AND. isub .eq. 2) goto 200 
                   if (isub3 .eq. 0 .AND. isub .eq. 3) goto 200 
                   if (isub4 .eq. 0 .AND. isub .eq. 4) goto 200 
                   if (isub5 .eq. 0 .AND. isub .eq. 5) goto 200 
                   if (isub6 .eq. 0 .AND. isub .ge. 6) goto 200 
c
                   nn=(n-1)*3
                   x1=x(nn+1,1)
                   y1=x(nn+1,2)
                     x2=x(nn+2,1)
                     y2=x(nn+2,2)
                       x3=x(nn+3,1)
                       y3=x(nn+3,2)
c
c                      curved line segment
                    kmax=10
                    kmax=mdiv
                    dr=2.0/kmax        
                    do k=1,kmax+0
                       r=-1+dr*(k-1)
                       h1 = 0.5*(1-r) - 0.5*(1-r*r)
                       h2 = 0.5*(1+r) - 0.5*(1-r*r)
                       h3 =                 (1-r*r)
                       ix1 = x1*h1 + x2*h2 + x3*h3
                       iy1 = y1*h1 + y2*h2 + y3*h3
                       r=-1+dr*(k-0)
                       h1 = 0.5*(1-r) - 0.5*(1-r*r)
                       h2 = 0.5*(1+r) - 0.5*(1-r*r)
                       h3 =                 (1-r*r)
                       ix2 = x1*h1 + x2*h2 + x3*h3
                       iy2 = y1*h1 + y2*h2 + y3*h3
c
                   call glnseg_m(ix1,iy1,ix2,iy2,linet,ipen,ioutps)
*                  call glnseg_m(ix2,iy2,ix3,iy3,linet,ipen,ioutps)
*                  call glnseg_m(ix3,iy3,ix1,iy1,linet,ipen,ioutps)
                       nstroke=nstroke+3
                       if (nstroke .gt. 500) then
                            nstroke=0
                            write(ioutps,*) ' stroke '
                       endif
                   enddo
c
 200            continue
c
      return
      end
c
c
c     DRAW the structure MESH as gray outline
      subroutine ps_drawmesh_HEX8(x,lseg,mode,nmprop
     &                          ,ioutps,xymov,ishow_m )
c
         implicit real*8 (a-h,o-z)
                 include 'commons.std'
c
         real*8   x(lseg,3)
         real*8   rr(3,3),xx(3),yy(3),zz(3)
         real*8   xymov(100,2)
         integer  nmprop(nel)
         icol=3
c
*           write(ilog,*)'@@ in draw rots x y z  ',rotx,roty,rotz
            call gn_rot3d_mat(rr)
c
c        BIG D0-Loop over all triangles
c        get info from Simplex.OUT file
         rewind(iout)
             read(iout,*) nel3,eigv
             write(ilog,*)'@@nel3 eigv', nel3,eigv
c
c        have all elements participate in scaling
         ielement=0
         ii=0
         do 100 i =  1,921111
            read(iout ,*,end=99) xx(1),yy(1),zz(1)
            read(iout ,*)        xx(2),yy(2),zz(2)
            read(iout ,*)        xx(3),yy(3),zz(3)
*           ielement=ielement+1
*           ielement=(i-1)/4+1
            ielement=(i-1)/6+1
*           isub=nmprop(ielement)
            isub=1
            xxm=xymov(isub,1)
            yym=xymov(isub,2)
c
            call gn_rot3d(rr,xx,yy,zz)
c
                       ii=ii+1
                       ir1 = (ii-1)*3 + 1
                       ir2 = (ii-1)*3 + 2
                       ir3 = (ii-1)*3 + 3
              x(ir1,1) = xx(1) + xxm
              x(ir1,2) = yy(1) + yym
              x(ir1,3) = 0.0  
               x(ir2,1) = xx(2) + xxm
               x(ir2,2) = yy(2) + yym
               x(ir2,3) = 0.0  
                x(ir3,1) = xx(3) + xxm  
                x(ir3,2) = yy(3) + yym
                x(ir3,3) = 0.0  
 100     continue
 99   continue
                ntri =ii
                ntri3=ii*3
                iscale=2
                call ttscaledat_3(x,lseg,ntri3,iscale,ilog)
           if (ishow_m .eq. 0) return
c
ctag1
c               only plot on subs
                nstroke=0
                do 200 n=1,ntri
**!!               isub=nmprop(n)
                   isub=1
                   if (isub1 .eq. 0 .AND. isub .eq. 1) goto 200 
                   if (isub2 .eq. 0 .AND. isub .eq. 2) goto 200 
                   if (isub3 .eq. 0 .AND. isub .eq. 3) goto 200 
                   if (isub4 .eq. 0 .AND. isub .eq. 4) goto 200 
                   if (isub5 .eq. 0 .AND. isub .eq. 5) goto 200 
                   if (isub6 .eq. 0 .AND. isub .ge. 6) goto 200 
c
                   nn=(n-1)*3
                   x1=x(nn+1,1)
                   y1=x(nn+1,2)
                     x2=x(nn+2,1)
                     y2=x(nn+2,2)
                       x3=x(nn+3,1)
                       y3=x(nn+3,2)
c
                       ix1 = x1
                       iy1 = y1
                         ix2 = x2
                         iy2 = y2
                           ix3 = x3
                           iy3 = y3
c
                   call glnseg_m(ix1,iy1,ix2,iy2,linet,ipen,ioutps)
                   call glnseg_m(ix2,iy2,ix3,iy3,linet,ipen,ioutps)
*                  call glnseg_m(ix3,iy3,ix1,iy1,linet,ipen,ioutps)
                       nstroke=nstroke+3
                       if (nstroke .gt. 500) then
                            nstroke=0
                            write(ioutps,*) ' stroke '
                       endif
c
 200            continue
c
      return
      end
c
c
      subroutine glnseg_2 (ix1,iy1,ix2,iy2,itag,ioutps)
c
         write(ioutps,81) ix1,iy1,' mt'
         write(ioutps,81) ix2,iy2,' lt'
 81      format(1x,2(i6,1x),a)
             call color(ioutps,itag)
c
      return
      end
c
c
      subroutine grule (ix,iy,ia,ib,itag,ioutps)
c
         write(ioutps,*) ix,iy,' mt'
         write(ioutps,*) ix+ia,iy,' lt'
             call color(ioutps,itag)
c
      return
      end
c
c
      subroutine greyrule (ix,iy,ia,ib,itag,ioutps)
c
         write(ioutps,*) ix,iy+ib/2,' mt'
         write(ioutps,*) ib,' setlinewidth'
         write(ioutps,*) ix+ia,iy+ib/2,' lt'
             call color(ioutps,itag)
c
      return
      end
c
c
      subroutine gshow (ix,iy,string,ioutps)
         character*(*) string
c
         write(ioutps,*) ix,iy,' mt (',string,') show '
c
      return
      end
c
c
       subroutine color(iout,itag )
          character*25 str25
c
          if (itag .ge. 1 .AND. itag .le. 11) then
              ipen = itag
          elseif (itag .ge. 12 .AND. itag .le. 20) then
              ipen = itag-10
          elseif (itag .ge. 21 .AND. itag .le. 30) then
              ipen = itag-20
          elseif (itag .eq. 0) then
              ipen = 0
          elseif (itag .eq. 31) then
              ipen = 31
          else
              ipen = 1
          endif
          if (ipen .eq. 0) str25=' 0.0 0.0 0.0 setrgbcolor' !black
          if (ipen .eq. 1) str25=' 0.0 0.0 1.0 setrgbcolor' !blue
          if (ipen .eq. 2) str25=' 0.5 0.0 1.0 setrgbcolor'
          if (ipen .eq. 3) str25=' 0.5 0.5 1.0 setrgbcolor'
          if (ipen .eq. 4) str25=' 0.5 1.0 1.0 setrgbcolor'
          if (ipen .eq. 5) str25=' 0.0 1.0 0.5 setrgbcolor'
          if (ipen .eq. 6) str25=' .86 .84 .82 setrgbcolor' !off-white
          if (ipen .eq. 7) str25=' 0.0 1.0 0.0 setrgbcolor' !green
          if (ipen .eq. 8) str25=' 1.0 0.5 1.0 setrgbcolor'
          if (ipen .eq. 9) str25=' 1.0 0.5 0.0 setrgbcolor'
          if (ipen .eq. 10) str25=' 0.9 0.2 0.1 setrgbcolor'
          if (ipen .eq. 11) str25=' 1.0 0.0 0.0 setrgbcolor' !red
          if (ipen .ge. 12) str25=' 0.9 0.9 0.9 setrgbcolor'
          if (ipen .lt. 0 ) str25=' 0.9 0.9 0.9 setrgbcolor'
          if (ipen .eq. 31) str25=' 0.7 0.7 0.7 setrgbcolor'
             write(iout,'(a)') str25
             write(iout,'(a)') 'stroke '
c
      return
      end
c
c
      subroutine glnseg_m (ix1,iy1,ix2,iy2,linet,ipen,ioutps)
c
         write(ioutps,81) ix1,iy1,' mt'
         write(ioutps,81) ix2,iy2,' lt'
 81      format(1x,2(i6,1x),a)
c
      return
      end
c
c
       subroutine head_m(iout )
c
           write(iout,'(a)') '4.000000 setlinewidth'
           write(iout,'(a)') '[60 10] 0 setdash'
            write(iout,'(a)') 'gsave '
            write(iout,'(a)') 'newpath 0 0 moveto'
            write(iout,'(a)') ' 0.6 0.6 0.6 setrgbcolor'
c
      return
      end
c
c
       subroutine foot_m(iout )
c
             write(iout,'(a)') 'stroke '
             write(iout,'(a)') 'grestore'
             write(iout,'(a)') '[ ] 0 setdash'
c
      return
      end
c
c
      subroutine gn_rot3d_mat(rr)
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         real*8   rr1(3,3),rr2(3,3),rr3(3,3),rrt(3,3),rr(3,3)
c
         pi=4.0*atan(1.0)
c
c          rotation matrices
           tt1=rotz*pi/180
           cc=cos(tt1)
           ss=sin(tt1)
           rr1(1,1)=cc
           rr1(1,2)=-ss
           rr1(1,3)= 0
             rr1(2,1)= ss
             rr1(2,2)= cc
             rr1(2,3)=  0
           rr1(3,1)= 0
           rr1(3,2)= 0
           rr1(3,3)= 1
c
           tt2=roty*pi/180
           cc=cos(tt2)
           ss=sin(tt2)
           rr2(1,1)= cc
           rr2(1,2)= 0
           rr2(1,3)= ss
             rr2(2,1)= 0
             rr2(2,2)= 1 
             rr2(2,3)= 0 
           rr2(3,1)= -ss
           rr2(3,2)=  0 
           rr2(3,3)=  cc
c
           tt3=rotx*pi/180
           cc=cos(tt3)
           ss=sin(tt3)
           rr3(1,1)= 1
           rr3(1,2)= 0
           rr3(1,3)= 0
             rr3(2,1)= 0
             rr3(2,2)= cc
             rr3(2,3)=-ss
           rr3(3,1)= 0
           rr3(3,2)=  ss
           rr3(3,3)=  cc
c
           do i=1,3
              do j=1,3
                 temp=0.0
                 do k=1,3
                    temp=temp+rr2(i,k)*rr1(k,j)
                 enddo
                 rrt(i,j)=temp
              enddo
           enddo
           do i=1,3
              do j=1,3
                 temp=0.0
                 do k=1,3
                    temp=temp+rr3(i,k)*rrt(k,j)
                 enddo
                 rr(i,j)=temp
              enddo
           enddo
*          write(ilog,*) 'ROTation'
*          do i=1,3
*             write(ilog,*) (rr(i,j), j=1,3)
*          enddo
c
c
      return
      end
c
c
c
ctag1
c     STATic CONTours
c     31   u  v  w   
c     32  ex ey exy  etc            
c     33  sx sy sxy  etc            
      subroutine ps_stat_cont(plt,ktags,nmax,nmprop
     &           ,factor,ivars,mode,ioutps,xymov,ishow_l)
c
         implicit real*8 (a-h,o-z)
                 include 'commons.std'
c
         integer  ktags(nmax), nmprop(nel)
         real*8   plt(nmax,2)
         real*8   params(50)
         real*8   xymov(100,2)
         character*30 str30
c
c
ctag1
                if (ivars .ge. 1 .AND. ivars .le. 6) then
                  call gn_s_contour(plt,ktags,nmax,kline,params,50,
     &                     ivars,factor,nmprop,
     &                     eigv,xymov )
                elseif (ivars .ge. 7 .AND. ivars .le. 8) then
c                   deformed shape
                    call gn_s_deform(plt,ktags,nmax,kline,params,50,
     &                     ivars,factor,nmprop,
     &                     xymov)
                else 
                endif
                write(ilog,*)'@@contour lines',kline
c
 81                          format(1x,a,i3 )
c
                 nll=kline*2
                 iscale=0
                 call ttscaledat(plt,nmax,nll,iscale,ilog)
c
c                plot the line segments
                   do n=1,kline
                      nn=n*2-1
                         x1= plt(nn+0,1)
                         y1= plt(nn+0,2)
                           x2= plt(nn+1,1)
                           y2= plt(nn+1,2)
                         ix1=x1
                         iy1=y1
                         ix2=x2
                         iy2=y2
c
                         iline=7
                         itag = ktags(n) 
                         if (ishow_l .le. -1) itag=0
                             call glnseg_2(ix1,iy1,ix2,iy2,itag,ioutps)
                   enddo
c
ctag2
c                  legend for stress colors
                   if (ishow_l .le. 0) goto 99
                   ix1= 9200
                   iy1=0 
                   ipen=8  
                   write(ioutps,'(a)')'%%%%   begin legend   %%%%%%%'
                   call head_m(ioutps)
                   sleg=0.5
                   ilegend=1000*sleg
                   write(ioutps,'(a,i6,a)') '/Helvetica findfont ',
     &                      ilegend,' scalefont setfont'
                   write(ioutps,'(a)') ' 250 setlinewidth'
                   write(ioutps,'(a)') '[70 0] 0 setdash'
c
                   write(ioutps,'(a)') ' 0 setlinecap'
                   ix1=10500
                   iy1= -1000
                   iy1= -900
                      ia=3200
                      ib=11500  
                      ipen=31
                      call greyrule(ix1,iy1,ia,ib,ipen,ioutps)
                   write(ioutps,'(a)') ' 250 setlinewidth'
                   write(ioutps,'(a)') ' 1 setlinecap'
c
                   do i=1,11
                      sig=params(i)
                      write(str30,'(e9.3)') sig
c
                      ix1= 11000 
                      iy1=    0 + (i-1)*1000 
                      ipen=i
                      ia=2200
                      ib=200  
                      call grule(ix1,iy1,ia,ib,ipen,ioutps)
c
                      write(ioutps,*) ' 0.0 0.0 0.0 setrgbcolor'
                      ixw=ix1       - 100
                      iyw=iy1 - 500 - 100
                      call gshow(ixw,iyw,str30(1:12),ioutps)
                   enddo
                   call foot_m(ioutps)
                   write(ioutps,'(a)')'%%%%   end legend   %%%%%%%'
c
 99        continue
      return
      end
c
c     contour line generator for triangular data
      subroutine gn_s_contour(xyz,ktags,nmax,kline,parmm,ipmax,
     &                        ivars,factor,nmprop,
     &                        eigv,xymov )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
ctag1
         integer  ktags(nmax), nmprop(nel)
         real*8   rr(3,3)
         real*8 x(3),y(3),z(3), zzi(12), zzj(12), zzk(12) 
         real*8 parmm(ipmax),posn(6,6)
         real*8   xyz(nmax,2), xymov(100,2)
         integer icode(50)
         pi=4.0*atan(1.0)
         kline=0
c 
         call gn_rot3d_mat(rr)
c
c          change parameter
c
         itmp2=iout
         rewind itmp2
c
              read(itmp2,*) nel3,eigv
              write(ilog,*)'@@nel3 eigv', nel3,eigv
c
         write(ilog,*)'@@ ivars ',ivars
         zmax=-1e16
         zmin=1e16
         ielement=0
c
c        for contour colors only use on subs
         do i=1,951111
            read(itmp2,*,end=98) x1,y1,z1,(zzi(j), j=1,6)
            pps = sqrt( (zzi(1)-zzi(2))**2+4*zzi(3)**2)
            ppm = sqrt( (zzi(4)-zzi(5))**2+4*zzi(6)**2)
            zzi( 7)= (zzi(1)+zzi(2))/2 + pps/2     
            zzi( 8)= (zzi(1)+zzi(2))/2 - pps/2     
            zzi( 9)= (zzi(4)+zzi(5))/2 + ppm/2     
            zzi(10)= (zzi(4)+zzi(5))/2 - ppm/2     
            zzi(11)=0
            zzi(12)=0
            zz1=zzi(ivars)
              ielement=(i-1)/3 +1
*             isub=nmprop(ielement)
              isub=1
              if (isub1 .eq. 0 .AND. isub .eq. 1) goto 222 
              if (isub2 .eq. 0 .AND. isub .eq. 2) goto 222 
              if (isub3 .eq. 0 .AND. isub .eq. 3) goto 222 
              if (isub4 .eq. 0 .AND. isub .eq. 4) goto 222 
              if (isub5 .eq. 0 .AND. isub .eq. 5) goto 222 
              if (isub6 .eq. 0 .AND. isub .ge. 6) goto 222 
            if (zz1 .gt. zmax) zmax=zz1
            if (zz1 .lt. zmin) zmin=zz1
 222        continue
         enddo
 98      continue
         write(ilog,*)'@@maxes', zmin,zmax
         if ( abs(zmax-zmin) .lt. 1.0e-10) then
c             no need to do contouring
              return
         endif
         zmax=zmax/factor
         zmin=zmin/factor
c
         write(ilog,*)'@@ivars', ivars
             nparm=11
             do n=1,nparm
                dz = abs(zmax-zmin)/(nparm-0)
                parmm(n) = zmin+dz*(n-1.0)
                icode(n) = n
                write(ilog,*)'@@ params: ', n,parmm(n), icode(n)
             enddo
*            nparm=10
         write(ilog,*)'@@ # of parameters: ',nparm
c
c
c     BIG D0-Loop over all triangles
      rewind(itmp2)
              read(itmp2,*) nel3,eigv
              write(ilog,*)'@@nel3 eigv', nel3,eigv
           write(ilog,*)'@@subs:',isub1,isub2,isub3,isub4,isub5,isub6
      ielement=0
c     only plot on subs
      do 200 i =  1,921111
           write(*,*) i
         read(itmp2 ,*,end=99) x(1),y(1),z(1),(zzi(j), j=1,6)
         read(itmp2 ,*)        x(2),y(2),z(2),(zzj(j), j=1,6)
         read(itmp2 ,*)        x(3),y(3),z(3),(zzk(j), j=1,6)
              ielement=(i-1)/4 +1
*             ielement=ielement+1
*        isub=nmprop(ielement)
         isub=1
         xxm=xymov(isub,1)
         yym=xymov(isub,2)
         if (isub1 .eq. 0 .AND. isub .eq. 1) goto 200 
         if (isub2 .eq. 0 .AND. isub .eq. 2) goto 200 
         if (isub3 .eq. 0 .AND. isub .eq. 3) goto 200 
         if (isub4 .eq. 0 .AND. isub .eq. 4) goto 200 
         if (isub5 .eq. 0 .AND. isub .eq. 5) goto 200 
         if (isub6 .eq. 0 .AND. isub .ge. 6) goto 200 
c
            pps = sqrt( (zzi(1)-zzi(2))**2+4*zzi(3)**2)
            ppm = sqrt( (zzi(4)-zzi(5))**2+4*zzi(6)**2)
            zzi( 7)= (zzi(1)+zzi(2))/2 + pps/2     
            zzi( 8)= (zzi(1)+zzi(2))/2 - pps/2     
            zzi( 9)= (zzi(4)+zzi(5))/2 + ppm/2     
            zzi(10)= (zzi(4)+zzi(5))/2 - ppm/2     
c
            pps = sqrt( (zzj(1)-zzj(2))**2+4*zzj(3)**2)
            ppm = sqrt( (zzj(4)-zzj(5))**2+4*zzj(6)**2)
            zzj( 7)= (zzj(1)+zzj(2))/2 + pps/2     
            zzj( 8)= (zzj(1)+zzj(2))/2 - pps/2     
            zzj( 9)= (zzj(4)+zzj(5))/2 + ppm/2     
            zzj(10)= (zzj(4)+zzj(5))/2 - ppm/2     
c
            pps = sqrt( (zzk(1)-zzk(2))**2+4*zzk(3)**2)
            ppm = sqrt( (zzk(4)-zzk(5))**2+4*zzk(6)**2)
            zzk( 7)= (zzk(1)+zzk(2))/2 + pps/2     
            zzk( 8)= (zzk(1)+zzk(2))/2 - pps/2     
            zzk( 9)= (zzk(4)+zzk(5))/2 + ppm/2     
            zzk(10)= (zzk(4)+zzk(5))/2 - ppm/2     
c
                call gn_rot3d(rr,x,y,z)
            z1=zzi(ivars)
            z2=zzj(ivars)
            z3=zzk(ivars)
c          contour
c
         do 31  iparm=1,nparm
            parm=parmm(iparm)
            kont=0
c
c           first side
            zz=(parm-z2)*(parm-z1)
            if (zz .lt. 0.0) then   
               kont=kont+1
               sclx = x(1) + (parm-z1)/(z2-z1)*(x(2)-x(1))
               scly = y(1) + (parm-z1)/(z2-z1)*(y(2)-y(1))
               posn(kont,1)=sclx
               posn(kont,2)=scly
            endif
c
c           second side
            zz=(parm-z3)*(parm-z2)
            if (zz .lt. 0.0) then   
               kont=kont+1
               sclx = x(2) + (parm-z2)/(z3-z2)*(x(3)-x(2))
               scly = y(2) + (parm-z2)/(z3-z2)*(y(3)-y(2))
               posn(kont,1)=sclx
               posn(kont,2)=scly
            endif
            if (kont.eq.2) goto 14
c       either moveto (kont=1) or lineto (kont=2)
c
c           third side
            zz=(parm-z1)*(parm-z3)
            if (zz .lt. 0.0) then   
               kont=kont+1
               sclx = x(3) + (parm-z3)/(z1-z3)*(x(1)-x(3))
               scly = y(3) + (parm-z3)/(z1-z3)*(y(1)-y(3))
               posn(kont,1)=sclx
               posn(kont,2)=scly
            endif
c          lineto
c      
c         record line
 14       continue
          if (kont.eq.2) then
              kline=kline+1
              if (kline .gt. nmax/2) then
                  kline=kline-1
                  write(ilog,*)'@@ exceeded  NMAX/2: ',kline,nmax
                  goto 99
              endif
              nn=kline*2-1
              xyz(nn+0,1) = posn(1,1) + xxm
              xyz(nn+0,2) = posn(1,2) + yym
              xyz(nn+1,1) = posn(2,1) + xxm
              xyz(nn+1,2) = posn(2,2) + yym
              ktags(kline)=icode(iparm)
          endif
 31     continue
c       end parameter
c
 200   continue
c      end loop over triangles
c
 99   continue
      return
      end
c
c
ctag1
c     SHAPE for HEX20 elements 
c     31   u  v  w   
      subroutine ps_shape_HEX20(plt,ktags,nmax,nmprop
     &           ,factor,ivars,mode,ioutps,xymov,ishow_l,mdiv)
c
         implicit real*8 (a-h,o-z)
                 include 'commons.std'
c
         integer  ktags(nmax), nmprop(nel)
         real*8   plt(nmax,2)
         real*8   params(50)
         real*8   xymov(100,2)
         character*30 str30
c
c
ctag1
                if (ivars .ge. 1 .AND. ivars .le. 6) then
                elseif (ivars .eq. 18) then
c                   deformed shape
                    call deform_HEX20(plt,ktags,nmax,kline,params,50,
     &                     ivars,factor,nmprop,
     &                     xymov,mdiv)
                elseif (ivars .eq. 181) then
c                   deformed shape
                    call deform_HEX20(plt,ktags,nmax,kline,params,50,
     &                     ivars,factor,nmprop,
     &                     xymov,mdiv)
                else 
                endif
                write(ilog,*)'@@ contour lines',kline
c
 81                          format(1x,a,i3 )
c
                 nll=kline*2
                 iscale=0
                 iscale=2
                 call ttscaledat(plt,nmax,nll,iscale,ilog)
                 if (ishow_l .eq. 0) return
c
c                plot the line segments
                   do n=1,kline
                      nn=n*2-1
                         x1= plt(nn+0,1)
                         y1= plt(nn+0,2)
                           x2= plt(nn+1,1)
                           y2= plt(nn+1,2)
                         ix1=x1
                         iy1=y1
                         ix2=x2
                         iy2=y2
*               write(ilog,*) n,ix1,iy1,ix2,iy2
c
                         iline=7
                         itag = ktags(n) 
*               write(ilog,*) n,itag 
                         if (ishow_l .le. -1) itag=0
                             call glnseg_2(ix1,iy1,ix2,iy2,itag,ioutps)
                   enddo
c
ctag2
c                  legend for stress colors
                   if (ishow_l .le. 0) goto 99
                   ix1= 9200
                   iy1=0 
                   ipen=8  
                   write(ioutps,'(a)')'%%%%   begin legend   %%%%%%%'
                   call head_m(ioutps)
                   sleg=0.5
                   ilegend=1000*sleg
                   write(ioutps,'(a,i6,a)') '/Helvetica findfont ',
     &                      ilegend,' scalefont setfont'
                   write(ioutps,'(a)') ' 250 setlinewidth'
                   write(ioutps,'(a)') '[70 0] 0 setdash'
c
                   write(ioutps,'(a)') ' 0 setlinecap'
                   ix1=10500
                   iy1= -1000
                   iy1= -900
                      ia=3200
                      ib=11500  
                      ipen=31
                      call greyrule(ix1,iy1,ia,ib,ipen,ioutps)
                   write(ioutps,'(a)') ' 250 setlinewidth'
                   write(ioutps,'(a)') ' 1 setlinecap'
c
                   do i=1,11
                      sig=params(i)
                      write(str30,'(e9.3)') sig
c
                      ix1= 11000 
                      iy1=    0 + (i-1)*1000 
                      ipen=i
                      ia=2200
                      ib=200  
                      call grule(ix1,iy1,ia,ib,ipen,ioutps)
c
                      write(ioutps,*) ' 0.0 0.0 0.0 setrgbcolor'
                      ixw=ix1       - 100
                      iyw=iy1 - 500 - 100
                      call gshow(ixw,iyw,str30(1:12),ioutps)
                   enddo
                   call foot_m(ioutps)
                   write(ioutps,'(a)')'%%%%   end legend   %%%%%%%'
c
 99        continue
      return
      end
c
c
c     deformation of HEX20 elements  
      subroutine deform_HEX20(xyz,ktags,nmax,kline,parmm,ipmax,
     &                   ivars,
     &                     factor,nmprop,
     &                     xymov,mdiv )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
ctag1
         integer  ktags(nmax),nmprop(nel)
         real*8 rr(3,3)
         real*8 x(3),y(3),z(3), zzi(10), zzj(10), zzk(10) 
         real*8 x0(3),y0(3),z0(3) 
         real*8 u(3),v(3),w(3)
         real*8 parmm(ipmax)
         real*8   xyz(nmax,2), xymov(100,2)
         integer icode(50)
         pi=4.0*atan(1.0)
         kline=0
c 
         call gn_rot3d_mat(rr)
c
c          change parameter
c
         itmp2=iout
         rewind itmp2
c
             read(itmp2,*) nel3,eigv
             write(ilog,*)'@@ nel3 eigv', nel3,eigv
         zmax=-1e16
         zmin=1e16
         do i=1,21111
            read(itmp2,*,end=98) x1,y1,z1,(zzi(j), j=1,3)
            zz1 = sqrt(zzi(1)*zzi(1) + zzi(2)*zzi(2) + zzi(3)*zzi(3))   
            if (zz1 .gt. zmax) zmax=zz1
            if (zz1 .lt. zmin) zmin=zz1
         enddo
 98      continue
         write(ilog,*)'@@ maxes', zmin,zmax
         if ( abs(zmax-zmin) .lt. 1.0e-10) then
              zmax=1.0
              zmin=0.0
         endif
c
c            need contour values for color documentation
             n=10
             nparm=11
             do n=1,nparm
                dz = abs(zmax-zmin)/(nparm-1)
                parmm(n) = zmin+dz*(n-1.0)
                icode(n) = n
                write(ilog,*)'@@ params: ', n,parmm(n), icode(n)
             enddo
             write(ilog,*)'@@ # of parameters: ',nparm
c
c
c     BIG D0-Loop over all 3 noded line segments
      rewind(itmp2)
             read(itmp2,*) nel3,eigv
             write(ilog,*)'@@nel3 eigv', nel3,eigv
      ielement=0
      do 200 i =  1,21111
         read(itmp2 ,*,end=99) x0(1),y0(1),z0(1),(zzi(j),j=1,3)
         read(itmp2 ,*)        x0(2),y0(2),z0(2),(zzj(j),j=1,3)
         read(itmp2 ,*)        x0(3),y0(3),z0(3),(zzk(j),j=1,3)
     &                       ,z1,z2,z3,ielement
            u(1)=zzi(1)
            u(2)=zzj(1)
            u(3)=zzk(1)
              v(1)=zzi(2)
              v(2)=zzj(2)
              v(3)=zzk(2)
            w(1)=zzi(3)
            w(2)=zzj(3)
            w(3)=zzk(3)
            uu1 = sqrt(zzi(1)*zzi(1) + zzi(2)*zzi(2) + zzi(3)*zzi(3))   
            uu2 = sqrt(zzj(1)*zzj(1) + zzj(2)*zzj(2) + zzj(3)*zzj(3))   
            uu3 = sqrt(zzk(1)*zzk(1) + zzk(2)*zzk(2) + zzk(3)*zzk(3))   
*        ielement=ielement+1
*        ielement = (i-1)/4+1
*          write(*,*) ielement
         isub=nmprop(ielement)
*        isub=0
*        xxm=xymov(isub,1)
*        yym=xymov(isub,2)
         xxm=0.0
         yym=0.0
         if (isub1 .eq. 0 .AND. isub .eq. 1) goto 200 
         if (isub2 .eq. 0 .AND. isub .eq. 2) goto 200 
         if (isub3 .eq. 0 .AND. isub .eq. 3) goto 200 
         if (isub4 .eq. 0 .AND. isub .eq. 4) goto 200 
         if (isub5 .eq. 0 .AND. isub .eq. 5) goto 200 
         if (isub6 .eq. 0 .AND. isub .ge. 6) goto 200 
c
         do k=1,3
            x(k)=x0(k)+factor*u(k)
            y(k)=y0(k)+factor*v(k)
            z(k)=z0(k)+factor*w(k)
         enddo
c
         call gn_rot3d(rr,x,y,z)
         x(1)=x(1)+xxm
         x(2)=x(2)+xxm
         x(3)=x(3)+xxm
         y(1)=y(1)+yym
         y(2)=y(2)+yym
         y(3)=y(3)+yym
c
c           curved line segment
         kmax=10
         kmax=mdiv
         dr=2.0/kmax        
         do k=1,kmax+0
            r=-1+dr*(k-1)
            h1 = 0.5*(1-r) - 0.5*(1-r*r)
            h2 = 0.5*(1+r) - 0.5*(1-r*r)
            h3 =                 (1-r*r)
            xx1 = x(1)*h1 + x(2)*h2 + x(3)*h3
            yy1 = y(1)*h1 + y(2)*h2 + y(3)*h3
c           zz = z(1)*h1 + z(2)*h2 + z(3)*h3
            r=-1+dr*(k-0)
            h1 = 0.5*(1-r) - 0.5*(1-r*r)
            h2 = 0.5*(1+r) - 0.5*(1-r*r)
            h3 =                 (1-r*r)
            xx2 = x(1)*h1 + x(2)*h2 + x(3)*h3
            yy2 = y(1)*h1 + y(2)*h2 + y(3)*h3
c           zz = z(1)*h1 + z(2)*h2 + z(3)*h3
              kline=kline+1
              if (kline .gt. nmax/2) then
                  kline=kline-1
                  goto 99
              endif
              nn=kline*2-1
              xyz(nn+0,1) = xx1
              xyz(nn+0,2) = yy1
              xyz(nn+1,1) = xx2
              xyz(nn+1,2) = yy2
c
              uu = uu1*h1 + uu2*h2 + uu3*h3
              if (uu .lt. 0.0) uu=0.0
              ztag = (uu*10/zmax) + 1
              if (ztag .gt. nparm) ztag=nparm
              ktags(kline)= ztag 
*        write(ilog,*)'@@ktags', kline,ktags(kline),uu
           enddo
c
 200  continue
c
 99   continue
*     close(itmp)
      end
c
cccccccccccccccccccccccccccccccccc
c
ctag1
c     SHAPE for HEX8 elements 
c     31   u  v  w   
      subroutine ps_shape_HEX8(plt,ktags,nmax,nmprop
     &           ,factor,ivars,mode,ioutps,xymov,ishow_l)
c
         implicit real*8 (a-h,o-z)
                 include 'commons.std'
c
         integer  ktags(nmax), nmprop(nel)
         real*8   plt(nmax,2)
         real*8   params(50)
         real*8   xymov(100,2)
         character*30 str30
c
c
ctag1
                if (ivars .ge. 1 .AND. ivars .le. 6) then
                elseif (ivars .ge.  9 .AND. ivars .le.  9) then
c                   deformed shape
                    call deform_HEX8(plt,ktags,nmax,kline,params,50,
     &                     ivars,factor,nmprop,
     &                     xymov)
                elseif (ivars .ge. 18 .AND. ivars .le. 18) then
c                   deformed shape
                    mdiv=2
                    call deform_HEX20(plt,ktags,nmax,kline,params,50,
     &                     ivars,factor,nmprop,
     &                     xymov,mdiv)
                else 
                endif
                write(ilog,*)'@@contour lines',kline
c
 81                          format(1x,a,i3 )
c
                 nll=kline*2
                 iscale=0
                 call ttscaledat(plt,nmax,nll,iscale,ilog)
c
c                plot the line segments
                   do n=1,kline
                      nn=n*2-1
                         x1= plt(nn+0,1)
                         y1= plt(nn+0,2)
                           x2= plt(nn+1,1)
                           y2= plt(nn+1,2)
                         ix1=x1
                         iy1=y1
                         ix2=x2
                         iy2=y2
c
                         iline=7
                         itag = ktags(n) 
*               write(ilog,*) n,ix1,iy1,ix2,iy2,itag
                if (itag .lt. 1) itag=1
                         if (ishow_l .le. -1) itag=0
                             call glnseg_2(ix1,iy1,ix2,iy2,itag,ioutps)
                   enddo
c
ctag2
c                  legend for stress colors
                   if (ishow_l .le. 0) goto 99
                   ix1= 9200
                   iy1=0 
                   ipen=8  
                   write(ioutps,'(a)')'%%%%   begin legend   %%%%%%%'
                   call head_m(ioutps)
                   sleg=0.5
                   ilegend=1000*sleg
                   write(ioutps,'(a,i6,a)') '/Helvetica findfont ',
     &                      ilegend,' scalefont setfont'
                   write(ioutps,'(a)') ' 250 setlinewidth'
                   write(ioutps,'(a)') '[70 0] 0 setdash'
c
                   write(ioutps,'(a)') ' 0 setlinecap'
                   ix1=10500
                   iy1= -1000
                   iy1= -900
                      ia=3200
                      ib=11500  
                      ipen=31
                      call greyrule(ix1,iy1,ia,ib,ipen,ioutps)
                   write(ioutps,'(a)') ' 250 setlinewidth'
                   write(ioutps,'(a)') ' 1 setlinecap'
c
                   do i=1,11
                      sig=params(i)
                      write(str30,'(e9.3)') sig
c
                      ix1= 11000 
                      iy1=    0 + (i-1)*1000 
                      ipen=i
                      ia=2200
                      ib=200  
                      call grule(ix1,iy1,ia,ib,ipen,ioutps)
c
                      write(ioutps,*) ' 0.0 0.0 0.0 setrgbcolor'
                      ixw=ix1       - 100
                      iyw=iy1 - 500 - 100
                      call gshow(ixw,iyw,str30(1:12),ioutps)
                   enddo
                   call foot_m(ioutps)
                   write(ioutps,'(a)')'%%%%   end legend   %%%%%%%'
c
 99        continue
      return
      end
c
c
c     deformation of HEX20 elements  
      subroutine deform_HEX8(xyz,ktags,nmax,kline,parmm,ipmax,
     &                   ivars,
     &                     factor,nmprop,
     &                     xymov )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
ctag1
         integer  ktags(nmax),nmprop(nel)
         real*8 rr(3,3)
         real*8 x(3),y(3),z(3), zzi(10), zzj(10), zzk(10) 
         real*8 x0(3),y0(3),z0(3) 
         real*8 u(3),v(3),w(3)
         real*8 parmm(ipmax)
         real*8   xyz(nmax,2), xymov(100,2)
         integer icode(50)
         pi=4.0*atan(1.0)
         kline=0
c 
         call gn_rot3d_mat(rr)
c
c          change parameter
c
         itmp2=iout
         rewind itmp2
c
             read(itmp2,*) nel3,eigv
             write(ilog,*)'@@nel3 eigv', nel3,eigv
         zmax=-1e16
         zmin=1e16
         do i=1,21111
            read(itmp2,*,end=98) x1,y1,z1,(zzi(j), j=1,3)
            zz1 = sqrt(zzi(1)*zzi(1) + zzi(2)*zzi(2) + zzi(3)*zzi(3))   
            if (zz1 .gt. zmax) zmax=zz1
            if (zz1 .lt. zmin) zmin=zz1
         enddo
 98      continue
         write(ilog,*)'@@maxes', zmin,zmax
         if ( abs(zmax-zmin) .lt. 1.0e-10) then
              zmax=1.0
              zmin=0.0
         endif
c
c            need contour values for color documentation
             n=10
             nparm=11
             do n=1,nparm
                dz = abs(zmax-zmin)/(nparm-0)
                parmm(n) = zmin+dz*(n-1.0)
                icode(n) = n
                write(ilog,*)'@@ params: ', n,parmm(n), icode(n)
             enddo
             write(ilog,*)'@@ # of parameters: ',nparm
c
c
c     BIG D0-Loop over all 3 noded line segments
      rewind(itmp2)
             read(itmp2,*) nel3,eigv
             write(ilog,*)'@@nel3 eigv', nel3,eigv
      ielement=0
      do 200 i =  1,21111
         read(itmp2 ,*,end=99) x0(1),y0(1),z0(1),(zzi(j), j=1,3)
         read(itmp2 ,*)        x0(2),y0(2),z0(2),(zzj(j), j=1,3)
         read(itmp2 ,*)        x0(3),y0(3),z0(3),(zzk(j), j=1,3)
            u(1)=zzi(1)
            u(2)=zzj(1)
            u(3)=zzk(1)
              v(1)=zzi(2)
              v(2)=zzj(2)
              v(3)=zzk(2)
            w(1)=zzi(3)
            w(2)=zzj(3)
            w(3)=zzk(3)
            uu1 = sqrt(zzi(1)*zzi(1) + zzi(2)*zzi(2) + zzi(3)*zzi(3))   
            uu2 = sqrt(zzj(1)*zzj(1) + zzj(2)*zzj(2) + zzj(3)*zzj(3))   
            uu3 = sqrt(zzk(1)*zzk(1) + zzk(2)*zzk(2) + zzk(3)*zzk(3))   
*        ielement=ielement+1
*        ielement = (i-1)/4+1
*        isub=nmprop(ielement)
         isub=1
         xxm=xymov(isub,1)
         yym=xymov(isub,2)
         if (isub1 .eq. 0 .AND. isub .eq. 1) goto 200 
         if (isub2 .eq. 0 .AND. isub .eq. 2) goto 200 
         if (isub3 .eq. 0 .AND. isub .eq. 3) goto 200 
         if (isub4 .eq. 0 .AND. isub .eq. 4) goto 200 
         if (isub5 .eq. 0 .AND. isub .eq. 5) goto 200 
         if (isub6 .eq. 0 .AND. isub .ge. 6) goto 200 
c
         do k=1,3
            x(k)=x0(k)+factor*u(k)
            y(k)=y0(k)+factor*v(k)
            z(k)=z0(k)+factor*w(k)
         enddo
c
         call gn_rot3d(rr,x,y,z)
         x(1)=x(1)+xxm
         x(2)=x(2)+xxm
         x(3)=x(3)+xxm
         y(1)=y(1)+yym
         y(2)=y(2)+yym
         y(3)=y(3)+yym
c
c           curved line segment
*        kmax=10
*        dr=2.0/kmax        
*        do k=1,kmax+0
*           r=-1+dr*(k-1)
*           h1 = 0.5*(1-r) - 0.5*(1-r*r)
*           h2 = 0.5*(1+r) - 0.5*(1-r*r)
*           h3 =                 (1-r*r)
*           xx1 = x(1)*h1 + x(2)*h2 + x(3)*h3
*           yy1 = y(1)*h1 + y(2)*h2 + y(3)*h3
c*          zz = z(1)*h1 + z(2)*h2 + z(3)*h3
*           r=-1+dr*(k-0)
*           h1 = 0.5*(1-r) - 0.5*(1-r*r)
*           h2 = 0.5*(1+r) - 0.5*(1-r*r)
*           h3 =                 (1-r*r)
*           xx2 = x(1)*h1 + x(2)*h2 + x(3)*h3
*           yy2 = y(1)*h1 + y(2)*h2 + y(3)*h3
c*          zz = z(1)*h1 + z(2)*h2 + z(3)*h3
              nn=kline*2+1
            xx1 = x(1)
            yy1 = y(1)
              xx2 = x(2)
              yy2 = y(2)
                xx3 = x(3)
                yy3 = y(3)
              xyz(nn+0,1) = xx1
              xyz(nn+0,2) = yy1
              xyz(nn+1,1) = xx2
              xyz(nn+1,2) = yy2
                xyz(nn+2,1) = xx2
                xyz(nn+2,2) = yy2
                xyz(nn+3,1) = xx3
                xyz(nn+3,2) = yy3
              kline=kline+2
              if (kline .gt. nmax/2) then
                  kline=kline-1
                  goto 99
              endif
c
              uu = (uu1 + uu2      )/2.0
              ktags(kline-1)= (uu*10/zmax) + 1
              uu = (      uu2 + uu3)/2.0
              ktags(kline)= (uu*10/zmax) + 1
*          enddo
c
 200  continue
c
 99   continue
*     close(itmp)
      end
c
cccccccccccccccccccccccccccccccccc
c
c
      subroutine gn_rot3d(rr,x,y,z)
c
         implicit real*8 (a-h,o-z)
         real*8 x(3),y(3),z(3),rr(3,3)
         pi=4.0*atan(1.0)
c 
         x1=x(1)
         y1=y(1)
         z1=z(1)
         x2=x(2)
         y2=y(2)
         z2=z(2)
         x3=x(3)
         y3=y(3)
         z3=z(3)
                x(1) = rr(1,1)*x1 + rr(1,2)*y1 + rr(1,3)*z1
                y(1) = rr(2,1)*x1 + rr(2,2)*y1 + rr(2,3)*z1
                z(1) = rr(3,1)*x1 + rr(3,2)*y1 + rr(3,3)*z1
c
                x(2) = rr(1,1)*x2 + rr(1,2)*y2 + rr(1,3)*z2
                y(2) = rr(2,1)*x2 + rr(2,2)*y2 + rr(2,3)*z2
                z(2) = rr(3,1)*x2 + rr(3,2)*y2 + rr(3,3)*z2
c
                x(3) = rr(1,1)*x3 + rr(1,2)*y3 + rr(1,3)*z3
                y(3) = rr(2,1)*x3 + rr(2,2)*y3 + rr(2,3)*z3
                z(3) = rr(3,1)*x3 + rr(3,2)*y3 + rr(3,3)*z3
c
      return
      end
c
c
c
c     Static EXTRA forms of data
      subroutine gn_s_deform(xyz,ktags,nmax,kline,parmm,ipmax,
     &                   ivars,
     &                     factor,nmprop,
     &                     xymov )
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
ctag1
         integer  ktags(nmax),nmprop(nel)
         real*8 rr(3,3)
         real*8 x(3),y(3),z(3), zzi(10), zzj(10), zzk(10) 
         real*8 u(3),v(3),w(3)
         real*8 parmm(ipmax)
         real*8   xyz(nmax,2), xymov(100,2)
         integer icode(50)
         pi=4.0*atan(1.0)
         kline=0
c 
         call gn_rot3d_mat(rr)
c
c          change parameter
c
         itmp2=iout
         rewind itmp2
c
             read(itmp2,*) nel3,eigv
             write(ilog,*)'@@nel3 eigv', nel3,eigv
         zmax=-1e16
         zmin=1e16
         do i=1,21111
            read(itmp2,*,end=98) x1,y1,z1,(zzi(j), j=1,3)
            zz1 = sqrt(zzi(1)*zzi(1) + zzi(2)*zzi(2) + zzi(3)*zzi(3))   
            if (zz1 .gt. zmax) zmax=zz1
            if (zz1 .lt. zmin) zmin=zz1
         enddo
 98      continue
         write(ilog,*)'@@maxes', zmin,zmax
         if ( abs(zmax-zmin) .lt. 1.0e-10) then
              zmax=1.0
              zmin=0.0
         endif
c
c            need contour values for color documentation
             n=10
             nparm=11
             do n=1,nparm
                dz = abs(zmax-zmin)/(nparm-0)
                parmm(n) = zmin+dz*(n-1.0)
                icode(n) = n
                write(ilog,*)'@@ params: ', n,parmm(n), icode(n)
             enddo
*            nparm=10
             write(ilog,*)'@@ # of parameters: ',nparm
c
c
c     BIG D0-Loop over all triangles
      rewind(itmp2)
             read(itmp2,*) nel3,eigv
             write(ilog,*)'@@nel3 eigv', nel3,eigv
      ielement=0
      do 200 i =  1,21111
         read(itmp2 ,*,end=99) x(1),y(1),z(1),(zzi(j), j=1,3)
         read(itmp2 ,*)        x(2),y(2),z(2),(zzj(j), j=1,3)
         read(itmp2 ,*)        x(3),y(3),z(3),(zzk(j), j=1,3)
            u(1)=zzi(1)
            u(2)=zzj(1)
            u(3)=zzk(1)
              v(1)=zzi(2)
              v(2)=zzj(2)
              v(3)=zzk(2)
            w(1)=zzi(3)
            w(2)=zzj(3)
            w(3)=zzk(3)
*        ielement=ielement+1
*        ielement = (i-1)/4+1
*        isub=nmprop(ielement)
         isub=1
         xxm=xymov(isub,1)
         yym=xymov(isub,2)
         if (isub1 .eq. 0 .AND. isub .eq. 1) goto 200 
         if (isub2 .eq. 0 .AND. isub .eq. 2) goto 200 
         if (isub3 .eq. 0 .AND. isub .eq. 3) goto 200 
         if (isub4 .eq. 0 .AND. isub .eq. 4) goto 200 
         if (isub5 .eq. 0 .AND. isub .eq. 5) goto 200 
         if (isub6 .eq. 0 .AND. isub .ge. 6) goto 200 
c
c
         call gn_rot3d(rr,x,y,z)
         x(1)=x(1)+xxm
         x(2)=x(2)+xxm
         x(3)=x(3)+xxm
         y(1)=y(1)+yym
         y(2)=y(2)+yym
         y(3)=y(3)+yym
c
        if (ivars .eq. 7) then
c           vector displcement
            xc=(x(1)+x(2)+x(3))/3.0
            yc=(y(1)+y(2)+y(3))/3.0
            call gn_rot3d(rr,u,v,w)
              uxx = (u(1)+u(2)+u(3) )/3.0 
              uyy = (v(1)+v(2)+v(3) )/3.0 
              uu  = sqrt(uxx*uxx + uyy*uyy)
c
              kline=kline+1
              if (kline .gt. nmax/2) then
                  kline=kline-1
                  goto 99
              endif
              nn=kline*2-1
              xyz(nn+0,1) = xc 
              xyz(nn+0,2) = yc 
              xyz(nn+1,1) = xc+ uxx*factor 
              xyz(nn+1,2) = yc+ uyy*factor 
              ktags(kline)= (uu*10/zmax) + 1
c
        elseif (ivars .eq. 8) then
c           deformed shape u-v-w
            call gn_rot3d(rr,u,v,w)
c
c           node 1
              kline=kline+1
              if (kline .gt. nmax/2) then
                  kline=kline-1
                  goto 99
              endif
              nn=kline*2-1
              xyz(nn+0,1) = x(1)+u(1)*factor
              xyz(nn+0,2) = y(1)+v(1)*factor  
              xyz(nn+1,1) = x(2)+u(2)*factor  
              xyz(nn+1,2) = y(2)+v(2)*factor 
              uu = sqrt(zzi(1)*zzi(1)+zzi(2)*zzi(2)+zzi(3)*zzi(3))/2
     &            +sqrt(zzj(1)*zzj(1)+zzj(2)*zzj(2)+zzj(3)*zzj(3))/2
              ktags(kline)= ( (uu-zmin)*9.5/(zmax-zmin) ) + 1
c
c           node 2
              uxx = zzj(1) 
              uyy = zzj(2) 
              kline=kline+1
              if (kline .gt. nmax/2) then
                  kline=kline-1
                  goto 99
              endif
              nn=kline*2-1
              xyz(nn+0,1) = x(2)+u(2)*factor
              xyz(nn+0,2) = y(2)+v(2)*factor  
              xyz(nn+1,1) = x(3)+u(3)*factor  
              xyz(nn+1,2) = y(3)+v(3)*factor 
              uu = sqrt(zzj(1)*zzj(1)+zzj(2)*zzj(2)+zzj(3)*zzj(3))/2
     &            +sqrt(zzk(1)*zzk(1)+zzk(2)*zzk(2)+zzk(3)*zzk(3))/2
***           ktags(kline)= (uu*10/zmax) + 1
              ktags(kline)= ( (uu-zmin)*9.5/(zmax-zmin) ) + 1
c
c           node 3
              uxx = zzk(1) 
              uyy = zzk(2) 
              kline=kline+1
              if (kline .gt. nmax/2) then
                  kline=kline-1
                  goto 99
              endif
              nn=kline*2-1
              xyz(nn+0,1) = x(3)+u(3)*factor
              xyz(nn+0,2) = y(3)+v(3)*factor  
              xyz(nn+1,1) = x(1)+u(1)*factor  
              xyz(nn+1,2) = y(1)+v(1)*factor 
              uu = sqrt(zzk(1)*zzk(1)+zzk(2)*zzk(2)+zzk(3)*zzk(3))/2
     &            +sqrt(zzi(1)*zzi(1)+zzi(2)*zzi(2)+zzi(3)*zzi(3))/2
***           ktags(kline)= (uu*10/zmax) + 1
              ktags(kline)= ( (uu-zmin)*9.5/(zmax-zmin) ) + 1
        endif
c
c
 200   continue
c      end loop over triangles
c
 99   continue
*     close(itmp)
      end
c
c
c
c     SCALE 
      subroutine ttscaledat(x,maxline,nplot,iscale,ilog)
         implicit real*8 (a-h,o-z)
         real*8 x(maxline,2)
c
         common /scalemax/ xxmin,xxmax, yymin,yymax
c
         if (iscale .ge. 1) then
c            find range
             xxmin=+1e20
             xxmax=-1e20
             yymin=+1e20
             yymax=-1e20
             do n=1,nplot
                if (x(n,1) .lt. xxmin) xxmin=x(n,1)
                if (x(n,1) .gt. xxmax) xxmax=x(n,1)
                if (x(n,2) .lt. yymin) yymin=x(n,2)
                if (x(n,2) .gt. yymax) yymax=x(n,2)
             enddo
             write(ilog,'(1x,a,6(g12.5,1x))') '@@ Orig Maxes: ',
     &                    xxmin,xxmax,yymin,yymax
c
c            when axes are the same
             if (iscale .eq. 2) then
                 dxx=xxmax-xxmin
                 dyy=yymax-yymin
                 if (dyy .gt. dxx) then
                     dxx=dyy 
                 else
                     dyy=dxx 
                 endif
                 xxmax=xxmin+dxx
                 yymax=yymin+dyy
             endif
c
c            choose bigger
*            write(ilog,*)'@@ PLOT Maxes: ', xxmin,xxmax,yymin,yymax
c
c            in case zero data
             if ( abs(xxmax-xxmin) .lt. 1.0e-20) then
                 xxmax=1.0
                 xxmin=-1.0
                 write(ilog,*)'@@ NEW xxMaxes: ', xxmin,xxmax
             endif
             if ( abs(yymax-yymin) .lt. 1.0e-20) then
                 yymax=1.0
                 yymin=-1.0
                 write(ilog,*)'@@ NEW yyMaxes: ',yymin,yymax
             endif
         endif
c
c        UPDATE PLOT
 30      continue
         xscale=10000.0/(xxmax-xxmin)
         yscale=10000.0/(yymax-yymin)
          ipx1=0
          ipx2=10000
          ipy1=0
          ipy2=10000
          xpscale=(ipx2-ipx1)/10000.
          ypscale=(ipy2-ipy1)/10000.
          xp0=ipx1
          yp0=ipy1
c
                do n=1,nplot
                   xs=x(n,1 )*xscale 
                   ys=x(n,2 )*yscale
                   xs=xs -xxmin*xscale 
                   ys=ys -yymin*yscale
                      xp=(xs   )*xpscale+xp0
                      yp=(ys   )*ypscale+yp0
                   x(n,1)=xp
                   x(n,2)=yp
                enddo
c
      return
      end
c
c
c     SCALE it
      subroutine ttscaledat_3(x,lseg,nplot,iscale,ilog)
         implicit real*8 (a-h,o-z)
         real*8 x(lseg,3)
c
         common /scalemax/ xxmin,xxmax, yymin,yymax
c
c
         if (iscale .ge. 1) then
c            find range
             xxmin=+1e20
             xxmax=-1e20
             yymin=+1e20
             yymax=-1e20
             zzmin=+1e20
             zzmax=-1e20
             do n=1,nplot
                if (x(n,1) .lt. xxmin) xxmin=x(n,1)
                if (x(n,1) .gt. xxmax) xxmax=x(n,1)
                if (x(n,2) .lt. yymin) yymin=x(n,2)
                if (x(n,2) .gt. yymax) yymax=x(n,2)
                if (x(n,3) .lt. zzmin) zzmin=x(n,3)
                if (x(n,3) .gt. zzmax) zzmax=x(n,3)
             enddo
             write(ilog,'(1x,a,6(g12.5,1x))') '@@ Orig Maxes: ',
     &                    xxmin,xxmax,yymin,yymax,zzmin,zzmax
c
c            when axes are the same
             if (iscale .eq. 2) then
                 dxx=xxmax-xxmin
                 dyy=yymax-yymin
                 if (dyy .gt. dxx) then
                     dxx=dyy 
                 else
                     dyy=dxx 
                 endif
                 xxmax=xxmin+dxx
                 yymax=yymin+dyy
             endif
*            write(ilog,*)'@@ PLOT Maxes: ', xxmin,xxmax,yymin,yymax
c
c            in case zero data
             if ( abs(xxmax-xxmin) .lt. 1.0e-20) then
                  xxmax=1.0
                  xxmin=-1.0
                  write(ilog,*)'@@ NEW xxMaxes: ',xxmin,xxmax
             endif
             if ( abs(yymax-yymin) .lt. 1.0e-20) then
                  yymax=1.0
                  yymin=-1.0
                  write(ilog,*)'@@ NEW yyMaxes: ',yymin,yymax
             endif
         endif
c
c        UPDATE PLOT
 30      continue
         ipx1=0
         ipy1=0
         ipx2=10000
         ipy2=10000
         xscale=10000.0/(xxmax-xxmin)
         yscale=10000.0/(yymax-yymin)
          xpscale=(ipx2-ipx1)/10000.
          ypscale=(ipy2-ipy1)/10000.
          xp0=ipx1
          yp0=ipy1
          zp0=0
c
c               SCALE data
                do n=1,nplot
                   xs=x(n,1 )*xscale 
                   ys=x(n,2 )*yscale 
                   zs=x(n,3 )*yscale 
                   xs=xs -xxmin*xscale 
                   ys=ys -yymin*yscale
                   zs=zs
                      xp=(xs   )*xpscale+xp0
                      yp=(ys   )*ypscale+yp0
                      zp=(zs   )*ypscale+zp0
                   x(n,1)=xp
                   x(n,2)=yp
                   x(n,3)=zp
                enddo
c
      return
      end
c
c
c
      subroutine head(iout,r_back,g_back,b_back)
         implicit real*8 (a-h,o-z)
c
          write(iout,'(a)') '%!PS-Adobe-2.0' 
          write(iout,'(a)') '%%from PlotForm Ver 1.64 ikayex, 1993'
          write(iout,'(a)') '%%BoundingBox: 0 0 0 0'
          write(iout,'(a)') '%begin(plot)'
        write(iout,'(a)') '/Times-Italic findfont 300 scalefont setfont'
           write(iout,'(a)') '/Sh {show} def'
           write(iout,'(a)') '/lt { lineto} def'
           write(iout,'(a)') '/rlt { rlineto} def'
           write(iout,'(a)') '/mt { moveto} def'
           write(iout,'(a)') '/rmt { rmoveto} def'
           write(iout,'(a)') '/in {72 mul} def'
          write(iout,'(a)') '%end(defs)'
         xmove=0.0
         ymove=0.0
*        xmove=pxmove*10000.0/10.0
*        ymove=pymove*10000.0/10.0
         xscale=100
         yscale=100
         write(iout,'(1x,g12.6,a,g12.6,a)')
     &                xmove,' in ',ymove,' in translate'
            write(iout,'(a)') '.04 .04 scale'
            write(iout,'(1x,2(g12.6,1x),a)') xscale,yscale,'  scale'
            write(iout,'(a)') '.01 .01 scale'
            write(iout,'(a)') '40 setlinewidth'
            write(iout,'(a)') '2 setlinecap'
            write(iout,'(a)') '1 setlinecap'
            write(iout,'(a)') '1 setlinejoin'
            write(iout,'(a)') 'gsave'
          write(iout,'(a)') '%%inner region'
          write(iout,'(a)') '%end(scale)'
                 write(iout,'(a)') '10.000000 setlinewidth'
           write(iout,'(a)') '[70  0] 0 setdash'
c
            write(iout,'(a)')'%%%%%%% background %%%%'
            back = r_back**2 + g_back**2 + b_back**2
            write(ilog,*)'@@ background ',back
            if (back .lt. 2.9) then
                write(iout,'(a)')'0 setlinecap'
                    write(iout,'(a)')'00       7000 mt'
                  write(iout,'(a)')'14000 setlinewidth'
                  write(iout,'(a)')'15000        7000 lt'
                write(iout,'(a)')'0.0 0.0 0.0 setrgbcolor '
                write(iout,'(1x,3(g12.6,1x),a)') r_back,g_back,b_back
     &                                      ,' setrgbcolor'
               write(iout,'(a)')'stroke '
            endif
            write(iout,'(a)')'%%%%%%%%%% end bg %%%%%'
c
            write(iout,'(a)') 'gsave '
            write(iout,'(a)') 'newpath 0 0 moveto'
c
      return
      end
c
c
       subroutine foot(iout )
c
             write(iout,'(a)') 'stroke '
             write(iout,'(a)') 'grestore'
             write(iout,'(a)') '[ ] 0 setdash'
           write(iout,'(a)') '%end(plot)'
             write(iout,'(a)') 'showpage'
c
      return
      end
c
c
