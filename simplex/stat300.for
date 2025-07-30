c
c1000 62 63
c     STATIC solution
      subroutine static_pre(stf,load,wk,
     &                      idbc,maxnode,
     &                      iprofv,iprofh,nloc,istiff,
     &                      xord,yord,zord
     &                      ,imain2,pscale1,pscale2)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer iprofv(neq), iprofh(neq)
         integer nloc(neq),idbc(maxnode,3)
         real*8 stf(maxstiff)
         real*8 load(neq), wk(neq)
         real*8 xord(nnp),yord(nnp),zord(nnp)
c
c
c        READ  stiffness matrix and load vector
         if (istiff .eq. 0) then
             write(*,*)'@@ reading <<StaDyn.STF<<'
             rewind(istf)
             read(istf) j1,j2  
             write(ilog,*)'@@ neq iband ',j1,j2  
             iloc=1
             do 8 i=1,neq  
                read(istf) (stf(iloc+j-1), j=1,iprofv(i))
                iloc = iloc + iprofv(i)
 8           continue
         endif
         if (istiff .eq. 2) then
             write(*,*)'@@ reading <<StaDyn.GEO<<'
             rewind(igeo)
             read(igeo) j1,j2  
             write(ilog,*)'@@ neq iband geo: ',j1,j2  
             iloc=1
                 do i=1,neq
                    iloc = nloc(i)
                    jmax = iprofv(i)
                    read(igeo) (wk(j), j=1,jmax)
                    do j=1,jmax
*                      write(*,*) stf(iloc+j-1) , gscale*wk(j),i,j
                       stf(iloc+j-1) = stf(iloc+j-1) + pscale2*wk(j)
                    enddo
                 enddo
         endif
c
*        write(*,*)'@@ reading <<StaDyn.LOD<<'
*        rewind(ilod)
*        do i=1,neq
*           read(ilod) load(i)
*           load(i)=pscale*load(i)
c           write(*,*) i,load(i)
*        enddo
c
c
         if (imain2 .eq. 63) then
c            use applied loads {P}
             write(*,*)'@@ reading <<StaDyn.LOD<<'
             rewind(ilod)
             do i=1,neq
                read(ilod) load(i)
                load(i) = pscale1*load(i)
*               write(*,*) i,load(i),grav(i)
             enddo
c
         elseif (imain2 .eq. 62) then
c            read applied loads {dP}
c        INPUT loads  
         do i=1,neq
            load(i)=0.0
         enddo
         write(*,*)'TYPE:  # of loads  | where <style: 1=node 2=xyz>'
         call zzwrt('  -->  ')
         read(ikbd,*) iload,norxyz 
         if (iload .gt. 40) then
             iload=40
             write(ilog,*)'@@  !! # loads > 40 !!'
         endif
         write(ilog,*) iload,norxyz,' ::# loads style'
         do i=1,iload
            write(*,*)'TYPE forces:   P_x | P_y | P_z '
            call zzwrt('  -->  ')
            read(ikbd,*) px1,py1,pz1
            write(ilog,86) px1,py1,pz1,' ::Pxyz'
            if (norxyz .eq. 1) then
c               node
                write(*,*) ' @@ TYPE:   node    '
                call zzwrt(' @@ --> ')
                read(ikbd,*) jnode
                write(*,*)' '
                write(ilog,*) jnode,' ::node  '
c
            elseif (norxyz .eq. 2) then
c               nearest
                write(*,*) ' @@ TYPE:   x y z   '
                call zzwrt(' @@ --> ')
                read(ikbd,*) x1,y1,z1
                write(*,*)' '
                write(ilog,86) x1,y1,z1,' ::xyz  '
c
c               find nearest node
                zlen0=1.0e20
                do n=1,nnp
                   zlen = sqrt((x1-xord(n))**2 + (y1-yord(n))**2
     &                                         + (z1-zord(n))**2)
                   if (zlen .lt. zlen0) then
                       zlen0 = zlen
                       jnode = n
                   endif
                enddo
                write(ilog,*)'@@ nearest node: ',jnode
            endif
c
                idof = idbc(jnode,1)
            if (idof .gt. 0) load(idof)=px1 
                idof = idbc(jnode,2)
            if (idof .gt. 0) load(idof)=py1
                idof = idbc(jnode,3)
            if (idof .gt. 0) load(idof)=pz1
         enddo
         endif
c
c
         write(*   ,*)'@@ STATIC:  Reloaded  [K]  & {P}  OK'
         write(ilog,*)'@@ STATIC:  Reloaded  [K]  & {P}  OK'
c
c         Decompose effective stiffness matrix
          write(*,'(a,i5,a,i5,a)')'@@  doing UDU decomp on [',
     &                                 neq,' X',iband,']'
          ierror=0
          call uduCOL_D(stf,maxstiff,neq,ierror,iprofv,nloc,z0,z1)
          write(ilog,83)'@@ Diag min max: ',z0,z1,z1/(abs(z0)+1.0e-10)
          if (ierror .eq. 0) then
              write(ilog,*)'&& ERROR: zero diagonal term'
              write(*   ,*)'&& ERROR: zero diagonal term'
              return
          endif
          write(*,'(a)')'@@ '
          write(*   ,*)'@@ UDU: mults & divs',ierror 
          write(ilog,*)'@@ UDU: mults & divs',ierror
*           call storeBND_f(stf,maxstiff,iprofv,nloc,iout,neq,iband,qms)
c
c         SOLVE for new displacements and save in case ...
            call bakCOL_D(stf,maxstiff,load,neq,wk,ierror,
     &                                      iprofv,iprofh,nloc)
c
          write(*,*)'@@ writing {u} to   >>StaDyn.DIS>>'
          write(ilog,*)'@@ writing {u} to   >>StaDyn.DIS>>'
          rewind(idis)
          do i=1,neq
             write(idis) wk(i)
          enddo
          rewind(isnp)
          do i=1,neq
             write(isnp) wk(i)
          enddo
 22       format(1x,7(g12.5,1x))
 81       format(1x,20(g12.6,1x))
 83       format(a,1x,20(g12.6,1x))
 86       format(1x,3(g12.6,1x), a)
      return
      end
c
c
c6000 61
c     STATIC solution
      subroutine static(stf,load,grav,wk,
     &                  iprofv,iprofh,nloc,istiff,pscale1,pscale2)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer iprofv(neq), iprofh(neq)
         integer nloc(neq)
         real*8 stf(maxstiff)
         real*8 load(neq), grav(neq),wk(neq)
c
c        READ  stiffness matrix and load vector
*        if (istiff .eq. 0) then
*           write(*,*)'@@ reading <<StaDyn.STF<<'
*           rewind(istf)
*           read(istf) j1,j2  
*           write(ilog,*)'@@ neq iband ',j1,j2  
*           iloc=1
*           do 8 i=1,neq  
*              read(istf) (stf(iloc+j-1), j=1,iprofv(i))
*              iloc = iloc + iprofv(i)
*8          continue
*        endif
*           do i=1,20  
*                  write(*,*) stf(i),i,' stat'
*           enddo
         write(*,*)'@@ reading <<StaDyn.LOD<<'
         rewind(ilod)
         do i=1,neq
            read(ilod) load(i)
            load(i) = pscale1*load(i) + pscale2*grav(i)
c           write(*,*) i,load(i)
         enddo
c
         write(*   ,*)'@@ STATIC:  Reloaded  [K]  & {P}  OK'
         write(ilog,*)'@@ STATIC:  Reloaded  [K]  & {P}  OK'
c
c         Decompose effective stiffness matrix
          write(*,'(1x,a,i5,a,i5,a)')'@@  doing UDU decomp on [',
     &                                 neq,' X',iband,']'
          ierror=0
          call uduCOL_D(stf,maxstiff,neq,ierror,iprofv,nloc,z0,z1)
          write(ilog,83)'@@ Diag min max: ',z0,z1,z1/(abs(z0)+1.0e-10)
*     call storeBND_f(stf,maxstiff,iprofv,nloc,iout,neq,iband,qms)
          if (ierror .eq. 0) then
              write(ilog,*)'@@ ERROR: zero diagonal term'
              write(*   ,*)'@@ ERROR: zero diagonal term'
              return
          endif
          write(*,'(a)')'@@ '
          write(*   ,*)'@@ UDU: mults & divs',ierror 
          write(ilog,*)'@@ UDU: mults & divs',ierror
*           write(iout,*)'[UDU]'
*           call storeBND_f(stf,maxstiff,iprofv,nloc,iout,neq,iband,qms)
c
c         SOLVE for new displacements and save in case ...
            call bakCOL_D(stf,maxstiff,load,neq,wk,ierror,
     &                                      iprofv,iprofh,nloc)
c
          write(*,*)'@@ writing {u} to   >>StaDyn.DIS>>'
          write(ilog,*)'@@ writing {u} to   >>StaDyn.DIS>>'
          rewind(idis)
          do i=1,neq
             write(idis) wk(i)
          enddo
c         create double snapshot
          rewind(isnp)
          psload=0
          wk00=0.0
          write(isnp) psload,neq
          do i=1,neq
             write(isnp) wk00 
          enddo
          psload=1
          write(isnp) psload,neq
          do i=1,neq
             write(isnp) wk(i)
          enddo
c
 22       format(1x,7(g12.5,1x))
 81       format(1x,20(g12.6,1x))
 83       format(1x,a,1x,20(g12.6,1x))
      return
      end
c
c
c
c     STATic OUTput routine for displacements and nodal loads
      subroutine statout( disp, dispful, idbc,maxnode)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer idbc(maxnode,3)
         real*8  disp(neq), dispful(nnp,3)  
         character*50 str1
c
c
c        RELOAD nodal  displacements
         rewind(idis)
         do i=1,neq
            read (idis) disp(i)
         enddo
         write(*,*)'@@ reloaded displacements'
c 
c        Fill in  the full displacement vector
         do 701 i=1, nnp
            do j=1,3
               ieqnum = idbc(i,j)
               if (ieqnum .gt. 0) then
                   dispful(i,j) = disp(ieqnum)
               else
                   dispful(i,j) = 0.0e0
               endif
            enddo
 701      continue 
          write(*,*)'@@ reassigned displacements'
c
          str1='STATIC output: Global DoF'
          call outdis (dispful,str1)
c
      return 
      end
c
c
c     OUTput DISplacements for static loads
      subroutine outdis(dispful, str1)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         character*(*) str1
         real*8 dispful(nnp,3)
c
         xmax=0.0 
         ymax=0.0 
         zmax=0.0 
c
c        Do_Loop for storing displacements. It also calculates the
c        maximum values for easy reference
c
         write(iout,1100) str1
         do 10 i= 1, nnp
            node = i
            val = dispful(i,1)
            if (abs(val) .gt. abs(xmax)) then
                xmax= val 
                nxmax= node
            endif
            val = dispful(i,2)
            if (abs(val) .gt. abs(ymax)) then
                ymax= val 
                nymax= node
            endif
            val = dispful(i,3)
            if (abs(val) .gt. abs(zmax)) then
                zmax= val 
                nzmax= node
            endif
            write(iout ,1110) node, (dispful(i,j),j=1,3)
 10      continue
c
         write(iout ,1120) 
         write(iout ,1130) nxmax,nymax,nzmax
         write(iout ,1140) xmax,ymax,zmax
c
c             FORMATS
 1100  format(1x,a50,/,
     & '  node    x-disp        y-disp        z-disp      rots' )  
 1110  format(1x,i5,2x,3(g12.6,2x),2x,' 0 0 0')
 1120  format(1x,  ' maximum displacements')
 1130  format(1x,'node',6x,6(i5, 8x)) 
 1140  format(1x,'value',2x,6(1pg12.4,1x))
c
      return
      end
c
c
c     OUTput MAXimum for TETrahedral quantities
      subroutine out_max_TET(quant,np1)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
*        character*(*) str1
         integer*4 max_node(6)
         real*8 quant(np1,6)
c
         do k=1,6
            max_node(k)=1
         enddo
c
c        calculates the  maximum values for easy reference
c
*       write(iout,1100) str1
        do i=1,np1
           do k=1,6
              node = max_node(k)
              if ( abs(quant(i,k)) .gt. abs(quant(node,k)) ) then
                   max_node(k)=i
              endif
           enddo
         enddo
c
        write(iout ,1120) 
        write(iout ,1130) (max_node(k), k=1,6)
        write(iout ,1140) (quant(max_node(k),k), k=1,6)
c
c             FORMATS
 1120  format(1x,  ' maximum values:')
 1130  format(1x,'node',6x,6(i5, 8x)) 
 1140  format(1x,'value',2x,6(1pg12.4,1x))
c
      return
      end
c
c
c     OUTput PRINcipal values for TETrahedral quantities
      subroutine out_prin_TET(quant,np1)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer*4 max_node(6)
         real*8 quant(np1,6),valmax(6),val(6)
         pi=4.0*atan(1.0)
c
           do k=1,6
              valmax(k)=0.0
           enddo
c
        do i=1,np1
           sxx = quant(i,1)
           syy = quant(i,2)
           szz = quant(i,3)
               sxy = quant(i,4)
               syz = quant(i,5)
               sxz = quant(i,6)
           zi1 = sxx + syy + szz 
           zi2 = sxx*syy + sxx*szz + syy*szz 
     &                - sxy**2 - sxz**2 - syz**2
           zi3 = sxx*syy*szz + 2*sxy*syz*sxz  
     &                 -sxx*syz**2 - syy*sxz**2 - szz*sxy**2
           svm1 = sqrt(zi1**2 - 3*zi2)
c
                  qq = (zi1*zi1 - 3*zi2)/9.0
                  rr = (-2*zi1**3 + 9*zi1*zi2 - 27*zi3)/54.0
*              write(*,*) qq,rr         
                  tt = acos(rr/sqrt(qq*qq*qq))
           zz1 = zi1/3.0 - 2*sqrt(qq)*cos( (tt + 0*pi)/3.0 )
           zz2 = zi1/3.0 - 2*sqrt(qq)*cos( (tt + 2*pi)/3.0 )
           zz3 = zi1/3.0 - 2*sqrt(qq)*cos( (tt - 2*pi)/3.0 )
                     zlam1=-20e12
                     if ( zz1 .gt. zlam1) zlam1 = zz1
                     if ( zz2 .gt. zlam1) zlam1 = zz2
                     if ( zz3 .gt. zlam1) zlam1 = zz3
                  zlam2=-20e12
                  if ( zz1 .gt. zlam2 .AND. zz1 .lt. zlam1) zlam2 = zz1
                  if ( zz2 .gt. zlam2 .AND. zz2 .lt. zlam1) zlam2 = zz2
                  if ( zz3 .gt. zlam2 .AND. zz3 .lt. zlam1) zlam2 = zz3
                    zlam3=-20e12
           if ( zz1 .gt. zlam3 .AND. zz1 .lt. zlam2) zlam3 = zz1
           if ( zz2 .gt. zlam3 .AND. zz2 .lt. zlam2) zlam3 = zz2
           if ( zz3 .gt. zlam3 .AND. zz3 .lt. zlam2) zlam3 = zz3
                  val(1)=zlam1
                  val(2)=zlam2
                  val(3)=zlam3
                  val(4)=svm1
                  val(5)=svm2
           svm2 = sqrt( (zlam1-zlam2)**2 + (zlam2-zlam3)**2
     &                 +(zlam3-zlam1)**2 )/sqrt(2.0)
                  write(iout,83) i,zlam1,zlam2,zlam3,svm1,svm2 
           do k=1,5
              if ( abs(val(k)) .gt. abs(valmax(k)) ) then
                   valmax(k)=val(k)
                   max_node(k)=i
              endif
           enddo
         enddo
c
        write(iout ,1120) 
        write(iout ,1130) (max_node(k), k=1,5)
        write(iout ,1140) (valmax(k), k=1,5)
c
c             FORMATS
 83    format(1x,i5,6x,20(g12.6, 4x)) 
 1120  format(1x,  ' maximum values:')
 1130  format(1x,'node',6x,6(i5, 8x)) 
 1140  format(1x,'value',2x,6(1pg12.4,1x))
c
      return
      end
c
c
c
c SUBROUTINE MATAbd 
c     This subroutine multiplies two matrices [A]{b}={c}
c     where A is a square matrix, b and c are vectors.
c
      subroutine matAbd( a, b, c, n) 
         implicit real*8 (a-h,o-z)
          real*8   a(n,n), b(n), c(n)
c
          do 10 i = 1,n
             c(i) = 0.0e0
             do 20 j = 1,n
                c(i) = c(i) + a(i,j) * b(j)
 20          continue
 10       continue
c
      return 
      end
c
c
c TRANS3D
c     This subroutine makes 3-D vector transformations.
c
      subroutine trns3dv (gloads,l,m,n,lloads,beta)
         real*8 gloads(12),r(3,3),lloads(12),lload3(3) ,gload3(3)
         real m,n,l,d
c
         pi=4.0*atan(1.0)
         cb=cos(beta*pi/180)
         sb=sin(beta*pi/180)
         d=sqrt(1-n**2)
       if (d.lt.1e-3)then
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
          r(2,1)  = -(m*cb+l*n*sb)/d
          r(2,2)  =  (l*cb-m*n*sb)/d
          r(2,3)  =  d*sb
          r(3,1)  =  (m*sb-l*n*cb)/d
          r(3,2)  = -(l*sb+m*n*cb)/d
          r(3,3)  =  d*cb
       endif
c
c
       do 11 i=1,12,3
          gload3(1)=gloads(i)
          gload3(2)=gloads(i+1)
          gload3(3)=gloads(i+2)
             call matABd(r,gload3,lload3,3)
          lloads(i)=lload3(1)
          lloads(i+1)=lload3(2)
          lloads(i+2)=lload3(3)
11     continue
c
      return
      end
c
c
c
c
c     POST ANalysis of displacements to get stresses etc
      subroutine postan ( disp,dispful, idbc, iglobal,
     &                   npijkm,maxnode,maxelem,xord,yord,zord,
     &                   strn,strs,forc,e_elm,s_elm,
     &                   e_elm_ip,s_elm_ip,wk_elm_n,
     &                   prop,nmprop,ippp )
c
             implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npijkm(maxelem,21), idbc(maxnode,3)
         real*8  disp(neq  ),  dispful(nnp,3)
         real*8  strs(nnp,6),strn(nnp,6), forc(nel,60)
         real*8  e_elm(nel,6),s_elm(nel,6)
         real*8  e_elm_ip(nel,27,6),s_elm_ip(nel,27,6)
         real*8  wk_elm_n(nel,20,6)
         real*8  xord(nnp), yord(nnp), zord(nnp)
*        real*8  strain(24),stress(24), force(12)
         real*8  strain(48),stress(48), force(24)
         real*8  disploc(3,6), temp(6)
         real*8  delu(24),nload(12),wk(20)
         real*8  xyz(3,20),ssi(27,6),ssn(20,6),se_av(6)
         real*8            ssi8(8,6),ssn8(8,6)
         real*8  uvw(3,20),dd(6,6)
         real*8  strain_ip(6,27),stress_ip(6,27), force_n(3,20)
c
         real*8 prop(ippp,10)
         integer nmprop(nel)
c
         pi=4.0*atan(1.0)
         ired_int=3
c
         iloop=0
 1       continue
         write(*,*)'  POST menu: '
         write(*,*)'           0=return '
*        write(*,*)'               LOCAL (member coords)  '
*        write(*,*)'           4=nodal force '
*        write(*,*)'           6=element (average) stress '
*        write(*,*)'           7=force sum '
*        write(*,*)'           8=specials  '
         write(*,*)'               GLOBAL '
         write(*,*)'          11=Global displacements'
         write(*,*)'          12=Global loads '
*        write(*,*)'          14=Global assembled DoF loads '
*        write(*,*)'          15=Global assembled nodal loads '
c
         write(*,*)'          '
         write(*,*)'          21=element strain [ue] '
         write(*,*)'          22=nodal   strain [ue] '
         write(*,*)'          25=element stress    '
         write(*,*)'          26=nodal   stress  '
         write(*,*)'          27=element stress specials  '
         write(*,*)'          28=nodal   stress specials  '
*        write(*,*)'               CONTOURs '
         write(*,*)'          31=store displacement data '
*        write(*,*)'          32=store strain       data '
         write(*,*)'          33=store stress       data '
*        write(*,*)'                       frames '
*        write(*,*)'          34=store displacement data '
*        write(*,*)'          35=store strain       data '
*        write(*,*)'          36=store stress       data '
         write(*,*)'          '
         write(*,*)'         121=element stress IP  '
         write(*,*)'         122=element stress  N  '
         write(*,*)'         124=nodal   stress  '
         write(*,*)'         127=element stress specials  '
         write(*,*)'         128=nodal   stress specials  '
         write(*,*)'               CONTOURs '
         write(*,*)'         131=store displacement data '
         write(*,*)'         132=store strain       data '
         write(*,*)'         133=store stress       data '
         write(*,*)'               SHAPES   '
         write(*,*)'         141=store displacement data '
         write(*,*)'               FORCES   '
         write(*,*)'         151=element nodal force        '
                 write(*,*)'               <<Hex20>>          '
                 write(*,*)'         41x=documented displacement '
                 write(*,*)'             411=displacement     '
                 write(*,*)'         42x=documented strain '
                 write(*,*)'             421=element strain IP '
                 write(*,*)'             422=element strain N  '
                 write(*,*)'             424=nodal   strain average '
                 write(*,*)'             427=element strain special '
                 write(*,*)'             428=nodal   strain special '
                 write(*,*)'         43x=documented stress '
                 write(*,*)'             431=element stress IP '
                 write(*,*)'             432=element stress N  '
                 write(*,*)'             434=nodal   stress average '
                 write(*,*)'             437=element stress special '
                 write(*,*)'             438=nodal   stress special '
                 write(*,*)'         440=force '
         call zzwrt(' SELECT--> ')
         read(ikbd,*,err=1)  ioutput
         write(ilog,*) ioutput, ' ::1=disp'
         if (ioutput .eq. 0) then
               return
         endif
c
 153     continue
         rewind (idis)
         do i=1,neq
            read(idis) disp(i)
         enddo
c          reassign displacements to dispful( )
           do 20 i=1,nnp
              do j=1,3
                 ieqnum=idbc(i,j)
                 if (ieqnum .eq. 0) then
                     dispful(i,j) = 0.0 
                 else
                     dispful(i,j) = disp(ieqnum)
                 endif
              enddo
 20        continue
c
           if (ioutput .eq. 11) then
               call outdis (dispful,'Static')
               goto 1
*          elseif (ioutput .ge. 12 .AND. ioutput .le. 15) then
*              call outload_D( disp,dispful, idbc, iglobal,
*    &                       npijkm,maxelem,maxnode,xord,yord,zord,
*    &                       prop,nmprop,ippp,ioutput )
*              goto 1
           elseif (ioutput .ge. 141 .AND. ioutput .le. 149) then
c              shapes, need only disps
               goto 100
           elseif (ioutput .ge. 410 .AND. ioutput .le. 419) then
               goto 100
           endif
c
c          no need to loop through elements multiple times
           if (iloop   .eq. 1) goto 100
               iloop=1
c
c          For each element, calculate the strain, stress at IP
           write(*,*)'@@ loop over elements ',nel
           write(*,*)' '
           nel51=nel/50+1
           call zzwrti(nel51)
           call zzwrt ('>')
           do 50 i=1,nel
              if (mod(i,50) .eq. 0) call zzwrt('.')
              neltype=npijkm(i,1)
              i1=npijkm(i,2)
              j1=npijkm(i,3)
              k1=npijkm(i,4)
              m1=npijkm(i,5)
                  x1=xord(i1)
                  x2=xord(j1)
                  x3=xord(k1)
                  x4=xord(m1)
                   y1=yord(i1)
                   y2=yord(j1)
                   y3=yord(k1)
                   y4=yord(m1)
                    z1=zord(i1)
                    z2=zord(j1)
                    z3=zord(k1)
                    z4=zord(m1)
              mat=nmprop(i)
                 e0=prop(mat,1)
                 g0=prop(mat,2)
                 r0=prop(mat,4)
              call dmat(e0,g0,dd)
c
              if (neltype .eq. 3) then
c                 frame
              elseif (neltype .eq. 4) then
c                 plate
              elseif (neltype .eq. 5) then
c                3-D solid TET
                 e0=prop(mat,1)
                 g0=prop(mat,2)
                 r0=prop(mat,4)
c
c                  obtain element disps
                   do kk=1,3
                      delu(kk+0)=dispful(i1,kk)
                      delu(kk+3)=dispful(j1,kk)
                      delu(kk+6)=dispful(k1,kk)
                      delu(kk+9)=dispful(m1,kk)
                   enddo
c
                 call convertTET_D(delu, 
     &                            e0,g0,
     &                            xyz,
*    &                            x1,x2,x3,x4,y1,y2,y3,y4,
*    &                            z1,z2,z3,z4,
     &                                  strain, stress,force )
c
c                 SAVE in array of element values
                  do j=1,6   
                     e_elm(i,j)=strain(j)*1.0e6
                     s_elm(i,j)=stress(j)
                  enddo
                  do j=1,12   
                     forc(i,j)=force(j)
                  enddo
c
              elseif (neltype .eq. 9) then
c                3-D solid HEX8
                 e0=prop(mat,1)
                 g0=prop(mat,2)
                 r0=prop(mat,4)
c
                 do k=1,8
                    node=npijkm(i,1+k)
                    xyz(1,k)=xord(node)
                    xyz(2,k)=yord(node)
                    xyz(3,k)=zord(node)
                 enddo
c
c                  obtain element disps
                 do j=1,8
                    do k=1,3
                       ielm=i
                       node=npijkm(ielm,1+j)
                       delu(k+(j-1)*3)=dispful(node,k)
                    enddo
                 enddo
c
*                write(iout,*) '@@ Elem: ',i
*                   write(iout,833) (delu(j), j=1,24)
*                   write(iout,*  ) ' '
 833             format(1x,3(g12.6,1x))
                 call convert_HEX(delu, 
     &                            e0,g0,
     &                            xyz,
     &                                  strain, stress,force )
c
c
c                 SAVE in array of element values
                  do j=1,8   
                     do k=1,6   
                        jk=(j-1)*6
                        e_elm_ip(i,j,k)=strain(jk+k)*1.0e6
                        s_elm_ip(i,j,k)=stress(jk+k)
                     enddo
                  enddo
                  do j=1,24   
                     forc(i,j)=force(j)
                  enddo
                  do k=1,6   
                     sume=0.0
                     sums=0.0
                     do j=1,8   
                        jk=(j-1)*6
                        sume=sume+strain(jk+k)
                        sums=sums+stress(jk+k)
                     enddo
                     e_elm(i,k)=(sume/8.0)*1.0e6
                     s_elm(i,k)=sums/8.0
                  enddo
c
              elseif (neltype .eq. 21) then
c                3-D solid HEX20
                 e0=prop(mat,1)
                 g0=prop(mat,2)
                 r0=prop(mat,4)
                 call dmat(e0,g0,dd)
c
                 kmax=neltype-1
                 do k=1,kmax
                    node=npijkm(i,1+k)
                    xyz(1,k)=xord(node)
                    xyz(2,k)=yord(node)
                    xyz(3,k)=zord(node)
c
                    uvw(1,k)=dispful(node,1)
                    uvw(2,k)=dispful(node,2)
                    uvw(3,k)=dispful(node,3)
                 enddo
c
ctag
*                write(iout,*) '@@ Elem: ',i
                 ired_int=2
                 ired_int=3
                 call convert_HEX20(neltype,xyz,uvw,dd,
     &                              ired_int,
     &                              strain_ip, stress_ip,force_n )
*                write(iout,81) (stress_ip(j,1), j=1,6)
c
c
c                 SAVE in array of element values
                  ielm=i
                  max_int=ired_int**3
                  do k=1,max_int   
                     do j=1,6   
                        e_elm_ip(ielm,k,j)=strain_ip(j,k)
                        s_elm_ip(ielm,k,j)=stress_ip(j,k)
                     enddo
                  enddo
                  nn=0
                  do k=1,kmax   
                     do j=1,3   
                        nn=nn+1
                        forc(ielm,nn)=force_n(j,k)
                     enddo
                  enddo
c                 element averages
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
c             bottom distinction between elements
 50        continue
c          end of loop over elements
*          write(iout,*)'before'
*          do i=1,nel
*             write(iout,82) i,(s_elm(i,j), j=1,6)
*          enddo
c
c          get nodal averages
           if (neltype .eq. 5) then
           call node_avg_TET_D(e_elm,s_elm, npijkm,maxelem,strn,strs)
           elseif (neltype .eq. 9) then
           call node_avg_HEX_D(e_elm,s_elm_ip,npijkm,maxelem,strn,strs)
           endif
c
c          write(iout,*)'after '
c          do i=1,nel
c             write(iout,82) i,(s_elm(i,j), j=1,6)
c          enddo
c
c
c          STORE results in <<OUT>> file
 100       continue
                  ipmax=ired_int**3
           write(ilog,'(a)')' @@'
*          write(iout,'(a)')' @@'
           write(*   ,'(a)')' writing output'
           if (ioutput .eq. 1) then
c
           elseif (ioutput .eq. 12) then
               write(iout,*)'Member NODAL FORCES: '
               write(iout,'(a,a)')' Elem   Node:    Fx         Fy    Fz'
               do i=1,nel
                  neltype=npijkm(i,1)
                  i1=npijkm(i,2)
                  j1=npijkm(i,3)
                  k1=npijkm(i,4)
                  m1=npijkm(i,5)
                  write(iout,83)   i,i1,(forc(i,j), j=1,3)
                  write(iout,84)     j1,(forc(i,j), j=4,6)
                  if (neltype .eq. 4) then
                      write(iout,84) k1,(forc(i,j), j=7,9)
                  endif
                  if (neltype .eq. 5) then
                      write(iout,84) k1,(forc(i,j), j=7,9)
                      write(iout,84) m1,(forc(i,j), j=10,12)
                  endif
               enddo
c
           elseif (ioutput .eq. 21) then
               write(iout,*)'Element STRAIN [ue]: '
               write(iout,*)' Elem:    Exx        Eyy          Ezz   ',
     &                     '          Exy        Eyz          Exz'
               do i=1,nel
                     write(iout,83) i,(e_elm(i,j), j=1,6)
               enddo 
*              do i=1,nel
*                 write(iout,23) i,(e_elm(i,j), j=1,6)
*              enddo 
               call out_max_TET(e_elm,nel)
               goto 1
c
ctag
           elseif (ioutput .eq. 22) then
               write(iout,*)'nodal STRAIN averages [ue]: '
               write(iout,*)' node:    Exx        Eyy          Ezz   ',
     &                     '          Exy        Eyz          Exz'
               do i=1,nnp
                  write(iout,23) i,(strn(i,j), j=1,6)
               enddo 
               call out_max_TET(strn,nnp)
c
           elseif (ioutput .eq. 24) then
               write(iout,*)'nodal STRAIN average specials: '
               write(iout,*)' Node:    E1        E2      E3   ',
     &                     '          '
               call out_prin_TET(strn,nnp)
c
           elseif (ioutput .eq. 25) then
               write(iout,*)'Element STRESS: '
               write(iout,'(3a)')'  Elem  : ',
     &                     '  Sxx         Syy          Szz',
     &                     '         Sxy           Syz          Sxz'
               do i=1,nel
                     write(iout,82) i,(s_elm(i,j), j=1,6)
               enddo 
*              call out_max_TET(s_elm,nel)
               goto 1
c
           elseif (ioutput .eq. 26) then
               write(iout,*)'NODAL STRESS averages: '
               write(iout,*)' node:    Sxx        Syy          Szz   ',
     &                     '          Sxy        Syz          Sxz'
               do i=1,nnp
                  write(iout,23) i,(strs(i,j), j=1,6)
               enddo 
               call out_max_TET(strs,nnp)
c
           elseif (ioutput .eq. 27) then
               write(iout,*)'element STRESS specials: '
               write(iout,*)' Elem:    S1              S2           S3',
     &                     '               von Mises  '
               call out_prin_TET(s_elm,nel)
c
           elseif (ioutput .eq. 28) then
               write(iout,*)'nodal STRESS average specials: '
               write(iout,*)' Node:    S1              S2           S3',
     &                     '              von Mises '
               call out_prin_TET(strs,nnp)
c
*          elseif (ioutput .eq. 7) then
*              call forsum_D(forc, npi,npj,npk)
c
c
           elseif (ioutput .eq. 31 .OR. ioutput .eq. 32
     &                             .OR. ioutput .eq. 33) then
ctag31
               rewind(iout)
               call contours_TET_D(strn,strs, npijkm,maxelem,
     &                       xord,yord,zord,dispful,ioutput)
c
*          elseif (ioutput .eq. 34) then
*              rewind(iout)
*              write(iout,*) nel*3,' 0.0 '
*              do i=1,nel
*                 i1=npi(i)
*                 j1=npj(i)
*                 k1=npk(i)
*        write(iout,81)xord(i1),yord(i1),zord(i1),(dispful(i1,j),j=1,6)
*        write(iout,81)xord(j1),yord(j1),zord(j1),(dispful(j1,j),j=1,6)
*        write(iout,81)xord(k1),yord(k1),zord(k1),(dispful(k1,j),j=1,6)
*              enddo
*          elseif (ioutput .eq. 35) then
*              rewind(iout)
*              write(iout,*) nel*3,' 0.0 '
*              do i=1,nel
*                 i1=npi(i)
*                 j1=npj(i)
*                 k1=npk(i)
*           write(iout,81)xord(i1),yord(i1),zord(i1),(strn(i,j),j=1,6)
*           write(iout,81)xord(j1),yord(j1),zord(j1),(strn(i,j),j=7,12)
*           write(iout,81)xord(k1),yord(k1),zord(k1),(strn(i,j),j=13,18)
*              enddo 
*          elseif (ioutput .eq. 36) then
*              rewind(iout)
*              write(iout,*) nel*3,' 0.0 '
*              do i=1,nel
*                 i1=npi(i)
*                 j1=npj(i)
*                 k1=npk(i)
*           write(iout,81)xord(i1),yord(i1),zord(i1),(strs(i,j),j=1,6)
*           write(iout,81)xord(j1),yord(j1),zord(j1),(strs(i,j),j=7,12)
*           write(iout,81)xord(k1),yord(k1),zord(k1),(strs(i,j),j=13,18)
*              enddo 
c
           elseif (ioutput .eq. 121) then
               write(iout,*)'Element STRESS int pt: '
               write(iout,'(3a)')'  Elem    IP: ',
     &                     '  Sxx         Syy          Szz',
     &                     '         Sxy           Syz          Sxz'
               do n=1,nel
                  eltype = npijkm(n,1)
                  kmax=neltype-1
                  ipmax=ired_int**3
                  j=1
                     write(iout,83) n,j
     &                               ,(s_elm_ip(n,1,k), k=1,6)
                  do j=2,ipmax
                     write(iout,84) j
     &                               ,(s_elm_ip(n,j,k), k=1,6)
                  enddo 
               enddo 
*              call out_max_TET(s_elm,nel)
               goto 1
c
           elseif (ioutput .eq. 122) then
               write(iout,*)'Element STRESS nodal: '
               write(iout,'(3a)')'  Elem     N: ',
     &                     '  Sxx         Syy          Szz',
     &                     '         Sxy           Syz          Sxz'
               do n=1,nel
                  do j=1,8
                     do k=1,6
                        ssi8(j,k)=s_elm_ip(n,j,k)
                     enddo 
                  enddo 
                  call node_HEX(ssi8,ssn8)
                     write(iout,83) n,npijkm(n,1+1),(ssn8(1,k), k=1,6)
                  do j=2,8
                     write(iout,84) npijkm(n,1+j),(ssn8(j,k), k=1,6)
                  enddo 
               enddo 
               goto 1
c
           elseif (ioutput .eq. 124) then
*          call node_avg_HEX_D(e_elm,s_elm_ip,npijkm,maxelem,strn,strs)
               write(iout,*)'NODAL STRESS averages: '
               write(iout,*)' node:    Sxx        Syy          Szz   ',
     &                     '          Sxy        Syz          Sxz'
               do i=1,nnp
                  write(iout,23) i,(strs(i,j), j=1,6)
               enddo 
               call out_max_TET(strs,nnp)
c
           elseif (ioutput .eq. 127) then
               write(iout,*)'element STRESS specials: '
               write(iout,*)' Elem:    S1              S2           S3',
     &                     '               von Mises  '
               call out_prin_TET(s_elm,nel)
c
           elseif (ioutput .eq. 128) then
               write(iout,*)'nodal STRESS average specials: '
               write(iout,*)' Node:    S1              S2           S3',
     &                     '              von Mises '
               call out_prin_TET(strs,nnp)
c
           elseif (ioutput .eq. 131 .OR. ioutput .eq. 132
     &                              .OR. ioutput .eq. 133) then
               rewind(iout)
               call contours_HEX_D(strn,strs, npijkm,maxelem,
     &                       xord,yord,zord,dispful,ioutput)
c
           elseif (ioutput .eq. 141 .OR. ioutput .eq. 142
     &                              .OR. ioutput .eq. 143) then
               rewind(iout)
               evout=0.0
*             stop
               call shapes_HEX_D(npijkm,maxelem,
     &                       xord,yord,zord,dispful,evout,ioutput)
c
           elseif (ioutput .eq. 151 .OR. ioutput .eq. 152
     &                              .OR. ioutput .eq. 153) then
c
               write(iout,*)'elem   NODAL FORCES: '
               write(iout,'(a,a)')' Elem   Node:    Fx         Fy    Fz'
               do i=1,nel
                  neltype=npijkm(i,1)
                  kmax=neltype-1
                  write(iout,83)   i,npijkm(i,2),(forc(i,j), j=1,3)
                  do k=2,kmax
                     k3=(k-1)*3
                      write(iout,84) npijkm(i,1+k),(forc(i,k3+j), j=1,3)
                  enddo
               enddo
c
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
c                STRESS
                 elseif (ioutput .eq. 431) then
                     write(iout,*)'Element STRESS int pt: '
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
                     write(iout,*)'Element STRESS nodal: '
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
                     write(iout,*)'NODAL STRESS averages: '
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
           endif
           goto 1
ctag9
c
 23        format(1x,i5,2x,6(g12.6,1x))
 81        format(1x,20(g12.6,1x))
 82        format(1x,i5,20(g11.5,1x))
 83        format(1x,i5,1x,i5,1x,6(g12.6,1x))
 84        format(1x,5x,1x,i5,1x,6(g12.6,1x))
 86        format(1x, 6(g12.6,1x))
c
 999  continue
      return
      end
c
c
c     CONVERT displacements to stresses etc
      subroutine convertCST_D(disploc, area,
     >                   e0,t0,g0,pl0,
     >                   xb1,xb2,xb3,yb1,yb2,yb3,
     >                   strain, stress,force )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         real*8  disploc(3,6)
         real*8 b(3,6),c(3,3),etemp(3,6),ekb(6,6)
         real*8 strain(18),stress(18)
         real*8 strainb(9),stressb(9),ub(6)
         real*8 force(18),forceb(6)
c
                znu0=e0/2.0/g0 -1.0
c             plane strain
                if (pl0 .lt. 0) xkap=3.-4.*znu0
c             plane strain
                if (pl0 .gt. 0) xkap=(3.-znu0)/(1.+znu0)
c
                c(1,1)=g0/(xkap-1.)*(xkap+1.)
                c(1,2)=g0/(xkap-1.)*(3.-xkap)
                c(1,3)=0.
                c(2,1)=c(1,2)
                c(2,2)=c(1,1)
                c(3,2)=0.
                c(3,1)=0.
                c(3,2)=0.
                c(3,3)=g0
c
                 area=((yb1-yb2)*(xb3-xb2) + (xb2-xb1)*(yb3-yb2))/2.
                   b(1,1)=(yb2-yb3)/2./area
                   b(1,2)=0.
                   b(1,3)=(yb3-yb1)/2./area
                   b(1,4)=0.
                   b(1,5)=(yb1-yb2)/2./area
                   b(1,6)=0.
                   b(2,1)=0.
                   b(2,2)=(xb3-xb2)/2./area
                   b(2,3)=0.
                   b(2,4)=(xb1-xb3)/2./area
                   b(2,5)=0.
                   b(2,6)=(xb2-xb1)/2./area
                   b(3,1)=b(2,2)
                   b(3,2)=b(1,1)
                   b(3,3)=b(2,4)
                   b(3,4)=b(1,3)
                   b(3,5)=b(2,6)
                   b(3,6)=b(1,5)
c
c               LOCAL displacements
                ub(1)=disploc( 1,1)
                ub(2)=disploc( 1,2)
c                w   =disploc( 1,3)
                ub(3)=disploc( 2,1)
                ub(4)=disploc( 2,2)
c                w   =disploc( 2,3)
                ub(5)=disploc( 3,1)
                ub(6)=disploc( 3,2)
c                w   =disploc( 3,3)
c
c            STRAIN
             do 30 j=1,3
                strainb(j)=0.
                do 32 k=1,6
                   strainb(j)=strainb(j)+b(j,k)*ub(k)
32              continue
30           continue
c
c            STRESS
             do 40 j=1,3
                stressb(j)=0.
                do 42 k=1,3
                   stressb(j)=stressb(j)+c(j,k)*strainb(k)
42              continue
40           continue
c
c          STIFFNESS
           do 50 i=1,3
              do 50 j=1,6
                 etemp(i,j)=0.
                 do 50 k=1,3
                      etemp(i,j)=etemp(i,j)+c(i,k)*b(k,j)
 50        continue
           do 60 i=1,6
              do 60 j=1,6
                 ekb(i,j)=0.
                 do 60 k=1,3
                    ekb(i,j)=ekb(i,j)+b(k,i)*etemp(k,j)*area*t0
 60         continue
c
c          NODAL FORCES wrt member axes
           do 70 i=1,6
              forceb(i)=0.
              do 70 k=1,6
                 forceb(i)=forceb(i)+ekb(i,k)*ub(k)
 70        continue
c          assign to large array
           force(1)=forceb(1)
           force(2)=forceb(2)
c          Fy etc = 0
           force(7)=forceb(3)
           force(8)=forceb(4)
c          Fy etc =0
           force(13)=forceb(5)
           force(14)=forceb(6)
c          Fy etc =0
c
           strain( 1)=strainb(1)
           strain( 2)=strainb(2)
           strain( 3)=strainb(3)
             strain( 7)=strainb(1)
             strain( 8)=strainb(2)
             strain( 9)=strainb(3)
               strain(13)=strainb(1)
               strain(14)=strainb(2)
               strain(15)=strainb(3)
           stress( 1)=stressb(1)
           stress( 2)=stressb(2)
           stress( 3)=stressb(3)
             stress( 7)=stressb(1)
             stress( 8)=stressb(2)
             stress( 9)=stressb(3)
               stress(13)=stressb(1)
               stress(14)=stressb(2)
               stress(15)=stressb(3)
c
c
      return
      end
c
c
ctag1
c     CONVERT displacements to stresses etc
      subroutine convertTET_D(delu, 
     &                            e0,g0,
     &                            xyz,
*    &                            x1,x2,x3,x4,y1,y2,y3,y4,
*    &                            z1,z2,z3,z4,
     >                   strain, stress,force )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         real*8  delu(12)
         real*8 bb(6,12),dd(6,6),ek12(12,12)
         real*8 strain(6),stress(6)
         real*8 force(12),xyz(3,20)
c
c
         call dmat(e0,g0,dd)
         call elmstf_TET(dd,
     &                     xyz,
     &                     ek12,bb)
c
c
c            STRAIN
             do j=1,6
                strain(j)=0.0
                do k=1,12
                   strain(j)=strain(j)+bb(j,k)*delu(k)
                enddo
             enddo
c
c            STRESS
             do j=1,6
                stress(j)=0.
                do k=1,6
                   stress(j)=stress(j)+dd(j,k)*strain(k)
                enddo
             enddo
c
c          NODAL FORCES wrt member axes
           do i=1,12
              force(i)=0.
              do k=1,12
                 force(i)=force(i)+ek12(i,k)*delu(k)
                enddo
             enddo
c
 70        continue
c
      return
      end
c
ctag1
c
c     CONVERT displacements to stresses etc
      subroutine convert_HEX20(neltype,
     &                            xyz,uvw,dd,
     &                   ired_int,
     &                   strain,stress,force )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
*        common /gauss_int/ xg,wgt
c
         real*8  u(60),xyz(3,20),uvw(3,20)
         real*8 dd(6,6),hh(20)
         real*8 xg(4,4),wgt(4,4),bb(6,60),ss(8),ee(8),ud(9),bd(9,60)
         real*8 stress(6,27), strain(6,27),force(3,20),ek(60,60)
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
                  call  stdm20(xyz,hh,bb,bd,det,ri,si,ti,nel,log)
*                 call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
*                 write(iout,*) ' Bb ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (Bb(i,j), j=1,60)
*                 enddo
*                 write(iout,*) ' Bd '
*                 do i=1,9
*                    write(iout,82) (Bd(i,j), j=1,24)
*                 enddo
         do i=1,6
            sum=0.0
            do k=1,kmax*3
               sum=sum+Bb(i,k)*u(k)
            enddo
            ee(i)=sum
         enddo
*                 write(iout,*) ' u  '
*                    write(iout,82) ( u(j), j=1,kmax*3)
*                    write(iout,82) ( ud(j), j=1,9)
*        z9=0.0
*        ee(1)= ud(1) + z9*( ud(1)**2 + ud(4)**2 + ud(7)**2 )/2.0
*        ee(2)= ud(5) + z9*( ud(2)**2 + ud(5)**2 + ud(8)**2 )/2.0
*        ee(3)= ud(9) + z9*( ud(3)**2 + ud(6)**2 + ud(9)**2 )/2.0
*        ee(4)= ud(2)+ud(4) + z9*(ud(1)*ud(2)+ud(4)*ud(5) +ud(7)*ud(8))
*        ee(5)= ud(6)+ud(8) + z9*(ud(2)*ud(3)+ud(5)*ud(6) +ud(8)*ud(9))
*        ee(6)= ud(3)+ud(7) + z9*(ud(1)*ud(3)+ud(4)*ud(6) +ud(7)*ud(9))
cctag1
*                       write(iout,* ) ' strain '
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
c
c
      return
      end
c
c
c     CONVERT displacements to stresses etc
      subroutine convert_HEX(delu, 
     &                            e0,g0,
     &                            xyz,
     &                   strain, stress,force )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         common /gauss_int/ xg,wgt
c
         real*8  delu(24),xyz(3,20)
*        real*8 ek24(24,24)
         real*8 strain(48),stress(48)
         real*8 force(12)
         real*8 dd(6,6),bb(6,24)
         real*8 xg(4,4),wgt(4,4),db(6,24),ss(8),ee(8)
         integer*4 icoord(8,3)
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
         data icoord/1,2,2,1,1,2,2,1
     &              ,1,1,2,2,1,1,2,2
     &              ,1,1,1,1,2,2,2,2 /
c
c
         call dmat(e0,g0,dd)
c
         nint=2
*        ip_n=0
c              do 80 lz=1,nint
c           do 80 ly=1,nint
c        do 80 lx=1,nint
         do 80 ii=1,8
*           irr=(-1)**lx
*           iss=(-1)**ly
*           itt=(-1)**lz
*           ip_n=ip_n+1
*           write(iout,*) irr,iss,itt
             lx = icoord(ii,1)
             ly = icoord(ii,2)
             lz = icoord(ii,3)
                  ri=xg(lx,nint)
                  si=xg(ly,nint)
                  ti=xg(lz,nint)
c
c                 deriv oper and jacabian
                  call stdm8(xyz,bb,det,ri,si,ti,nel,ilog)
*                 write(*,*) ' B ',lx,ly,lz
*                 do i=1,6
*                    write(*,82) (bb(i,j), j=1,24)
*                 enddo
c
c                 add contib to stiff
*                 wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
*                 write(*,*) ' wt ',wt
                     do j=1,6
                        sum = 0.0  
                        do kk=1,24
                           sum = sum + bb(j,kk)*delu(kk)
                        enddo
                        ee(j) = sum
                     enddo
*                 write(iout,82) (ee(j)*1.0e6, j=1,6)
                  do i=1,6
                     sum = 0.0  
                     do kk=1,6
                        sum = sum + dd(i,kk)*ee(kk)
                     enddo
                     ss(i)=sum
                  enddo
*                 write(iout,82) (ss(j), j=1,6)
*                 do 70 i=1,6
*                    do 40 j=1,24
*                       sum = 0.0  
*                       do kk=1,6
*                          sum = sum + dd(i,kk)*bb(kk,j)
*                       enddo
*                       db(i,j) = sum
*40                  continue
*70               continue
*                 do i=1,6
*                    sum = 0.0  
*                    do kk=1,24
*                       sum = sum + db(i,kk)*delu(kk)
*                    enddo
*                    ss(i)=sum
*                 enddo
                  s1 = ss(1)+ss(2)+ss(3)
                  s2 = ss(1)*ss(2) +ss(2)*ss(3) +ss(3)*ss(1) 
     &                -ss(4)**2 - ss(5)**2 - ss(6)**2
                  svm = sqrt(s1*s1 - 3.0*s2)
c                 write(iout,82) (ss(j), j=1,6),svm
                  write(iout,*) '@@ SvM: ',ii,svm
*          write(*,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (dk(i,j), j=1,24)
*          enddo
*            jk=(ip_node-1)*6
             jk=(ii     -1)*6
             do k=1,6
                stress(jk+k)=ss(k)
                strain(jk+k)=ee(k)
             enddo
 80      continue
 82      format(1x,30(g12.6,1x))
ctag9
c
         return
c
c            STRAIN
             do j=1,6
                strain(j)=0.0
                do k=1,12
                   strain(j)=strain(j)+bb(j,k)*delu(k)
                enddo
             enddo
c
c            STRESS
             do j=1,6
                stress(j)=0.
                do k=1,6
                   stress(j)=stress(j)+dd(j,k)*strain(k)
                enddo
             enddo
c
c          NODAL FORCES wrt member axes
           do i=1,12
              force(i)=0.
              do k=1,12
*                force(i)=force(i)+ek12(i,k)*delu(k)
                enddo
             enddo
c
 77        continue
c
      return
      end
c
      subroutine node_HEX(ssi,ssn)
         implicit real*8 (a-h,o-z)
         real*8  ssi(8,6),aai(8,8),ssn(8,6)
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
ctag9
*       do n=1,nel
           do j=1,6
              do i=1,8
                 sum=0.0
                 do k=1,8
                    sum = sum + aai(i,k)*ssi(k,j)
                 enddo
                 ssn(i,j) = sum
              enddo
           enddo
*                 do j=1,8
*                    write(iout,83) n,j,(ssn(j,k), k=1,6)
*                 enddo 
*       enddo
c
 83      format(1x,i5,1x,i5,1x,6(g12.6,1x))
      return
      end
c
c
c     AVERAGE element values to get nodal values
      subroutine node_avg_HEX_D(e_elm,s_elm_ip,npijkm,maxelem,strn,strs)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21)
         real*8  strs(nnp,6),strn(nnp,6),e_elm(nel,6),s_elm_ip(nel,27,6)
         real*8  sav(6),eav(6)
         real*8  ssi(8,6),ssn(8,6)
c
c        NODAL AVERAGE
         do 60 node=1,nnp
            do k=1,6
               sav(k)=0.0
               eav(k)=0.0
            enddo
c
c           find elements attached to this node
            attach=0.0
            do 70 jelm=1,nel
               neltype=npijkm(jelm,1)
               do k=1,8
                  if (node .eq. npijkm(jelm,1+k)) then
                      kth=k
                      goto 72
                  endif
               enddo
               goto 70
c
 72            continue  
                  do jj=1,8
                     do kk=1,6
                        ssi(jj,kk)=s_elm_ip(jelm,jj,kk)
                     enddo 
                  enddo 
                  call node_HEX(ssi,ssn)
                   attach=attach+1
                   do k=1,6
                      eav(k)=eav(k)+ssn(kth,k)
                      sav(k)=sav(k)+ssn(kth,k)
                   enddo
c
 70         continue
            if (attach .lt. 0.5) attach=1
            do k=1,6
               strn(node,k)=eav(k)/attach
               strs(node,k)=sav(k)/attach
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
c     AVERAGE element values to get nodal values
      subroutine node_avg_D(strs,strn,s_avg, npijkm,maxelem, ioutput)
c
ctag1
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21)
         real*8  strs(nel,18), strn(nel ,18),s_avg(nnp,6)
         real*8  sav(6),eav(6),ff(6)
c
c        NODAL AVERAGE
         do 60 i=1,nnp
            do k=1,6
               sav(k)=0.0
               eav(k)=0.0
            enddo
c
c           find elements attached to this node
            attach=0.0
            do 70 j=1,nel
               neltype=npijkm(j,1)
               if (neltype .eq. 3) goto 70
               ijk=0
c              which (1,2,3) node is it
               if (i .eq. npijkm(j,2)) ijk=1
               if (i .eq. npijkm(j,3)) ijk=2
               if (neltype .eq. 4) then
                   if (i .eq. npijkm(j,4)) ijk=3
               endif
               if (ijk .gt. 0) then
                   attach=attach+1
                   do k=1,6
                      eav(k)=eav(k)+strn(j,(ijk-1)*6+k)
                      sav(k)=sav(k)+strs(j,(ijk-1)*6+k)
                   enddo
               endif
 70         continue
            if (attach .lt. 0.5) attach=1
            do k=1,6
               sav(k)=sav(k)/attach
               eav(k)=eav(k)/attach*1e6
               s_avg(i,k)=sav(k)
            enddo
c
            if (ioutput .eq. 2) then
                write(iout,52) i,(eav(k), k=1,6)
            elseif (ioutput .eq. 3) then
                write(iout,52) i,(sav(k), k=1,6)
            elseif (ioutput .eq. 8) then
                ff(1) = sav(1) + sav(2)
                ff(2) = ( (sav(1) - sav(2))**2 + 4.0*sav(3)**2 )**0.5
                ff(3) = (ff(1) - ff(2) )/2.0
                ff(4) = (ff(1) + ff(2) )/2.0
                ff(5) = sav(1) - sav(2)
                ff(6) = sav(3) 
                write(iout,52) i,(ff(k), k=1,6)
            endif
 52               format(1x,i5,1x,8(g12.5,1x))
 60      continue
c        bottom of nodal loop
            if (ioutput .ne. 6) return
c
c        ELEMENT AVERAGE
         do i=1,nel
            ii = npijkm(i,2)
            jj = npijkm(i,3)
            kk = npijkm(i,4)
            mm = npijkm(i,5)
            do k=1,6
               sav(k)=( s_avg(ii,k)+s_avg(jj,k)+s_avg(kk,k) )/3.0
            enddo
c
            write(iout,52) i,(sav(k), k=1,6)
         enddo
c
      return
      end
c
c
c
c     AVERAGE element values to get nodal values
      subroutine node_avg_TET_D(e_elm,s_elm,npijkm,maxelem,strn,strs)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21)
         real*8  strs(nnp,6), strn(nnp,6),e_elm(nel,6),s_elm(nel,6)
         real*8  sav(6),eav(6)
c
c        NODAL AVERAGE
         do 60 i=1,nnp
            do k=1,6
               sav(k)=0.0
               eav(k)=0.0
            enddo
c
c           find elements attached to this node
            attach=0.0
            do 70 j=1,nel
               neltype=npijkm(j,1)
*              if (neltype .eq. 3) goto 70
*              ijk=0
c              which (1,2,3) node is it
               if (i .eq. npijkm(j,2)) goto 72
               if (i .eq. npijkm(j,3)) goto 72
               if (neltype .ge. 4) then
                   if (i .eq. npijkm(j,4)) goto 72
               endif
               if (neltype .eq. 5) then
                   if (i .eq. npijkm(j,5)) goto 72
               endif
               goto 70
c
 72            continue  
                   attach=attach+1
                   do k=1,6
                      eav(k)=eav(k)+e_elm(j,k)
                      sav(k)=sav(k)+s_elm(j,k)
                   enddo
c
 70         continue
            if (attach .lt. 0.5) attach=1
            do k=1,6
               strn(i,k)=eav(k)/attach
               strs(i,k)=sav(k)/attach
            enddo
c
*           if (ioutput .eq. 2) then
*               write(iout,52) i,(eav(k), k=1,6)
*           elseif (ioutput .eq. 3) then
*               write(iout,52) i,(sav(k), k=1,6)
*           elseif (ioutput .eq. 8) then
*               ff(1) = sav(1) + sav(2)
*               ff(2) = ( (sav(1) - sav(2))**2 + 4.0*sav(3)**2 )**0.5
*               ff(3) = (ff(1) - ff(2) )/2.0
*               ff(4) = (ff(1) + ff(2) )/2.0
*               ff(5) = sav(1) - sav(2)
*               ff(6) = sav(3) 
*               write(iout,52) i,(ff(k), k=1,6)
*           endif
 52               format(1x,i5,1x,8(g12.5,1x))
 60      continue
c        bottom of nodal loop
*           if (ioutput .ne. 6) return
c
c        ELEMENT AVERAGE
*        do i=1,nel
*           ii = npijkm(i,2)
*           jj = npijkm(i,3)
*           kk = npijkm(i,4)
*           mm = npijkm(i,5)
*           do k=1,6
*              sav(k)=( s_avg(ii,k)+s_avg(jj,k)+s_avg(kk,k) )/3.0
*           enddo
c
*           write(iout,52) i,(sav(k), k=1,6)
*        enddo
c
      return
      end
c
c
c     CONTOUR nodal values
      subroutine contours_TET_D(strn,strs, npijkm,maxelem,
     &                    xord,yord,zord,
     &                    dispful,ioutput              )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21)
         real*8  strs(nnp,6), strn(nnp ,6)
         real*8  xord(nnp), yord(nnp), zord(nnp)
         real*8  dispful(nnp,3)
         real*8  sav(6),eav(6)
c
ctag31
         rewind(iout)
         write(iout,*) nel,gnum1
         icont=ioutput
c        displacements do not need averaging
         if (icont .eq. 31) then
             do i=1,nel
                neltype=npijkm(i,1)
                i1=npijkm(i,2)
                j1=npijkm(i,3)
                k1=npijkm(i,4)
                m1=npijkm(i,5)
                  x1=xord(i1)
                  x2=xord(j1)
                  x3=xord(k1)
                  x4=xord(m1)
                    y1=yord(i1)
                    y2=yord(j1)
                    y3=yord(k1)
                    y4=yord(m1)
                  z1=zord(i1)
                  z2=zord(j1)
                  z3=zord(k1)
                  z4=zord(m1)
                  if ( neltype .eq. 5 ) then
c???                  only save plate nodal averages for contouring
                      write(iout,81) x1,y1,z1,(dispful(i1,j), j=1,3)
                      write(iout,81) x2,y2,z2,(dispful(j1,j), j=1,3)
                      write(iout,81) x3,y3,z3,(dispful(k1,j), j=1,3)
                        write(iout,81) x1,y1,z1,(dispful(i1,j), j=1,3)
                        write(iout,81) x2,y2,z2,(dispful(j1,j), j=1,3)
                        write(iout,81) x4,y4,z4,(dispful(m1,j), j=1,3)
                      write(iout,81) x2,y2,z2,(dispful(j1,j), j=1,3)
                      write(iout,81) x3,y3,z3,(dispful(k1,j), j=1,3)
                      write(iout,81) x4,y4,z4,(dispful(m1,j), j=1,3)
                        write(iout,81) x3,y3,z3,(dispful(k1,j), j=1,3)
                        write(iout,81) x1,y1,z1,(dispful(i1,j), j=1,3)
                        write(iout,81) x4,y4,z4,(dispful(m1,j), j=1,3)
                  endif
             enddo
             return
         endif
c
c        OTHERWISE  do stresses and strains
c
         do i=1,nel
                neltype=npijkm(i,1)
                i1=npijkm(i,2)
                j1=npijkm(i,3)
                k1=npijkm(i,4)
                m1=npijkm(i,5)
                x1=xord(i1)
                x2=xord(j1)
                x3=xord(k1)
                x4=xord(m1)
                  y1=yord(i1)
                  y2=yord(j1)
                  y3=yord(k1)
                  y4=yord(m1)
                    z1=zord(i1)
                    z2=zord(j1)
                    z3=zord(k1)
                    z4=zord(m1)
c
                      write(iout,81) x1,y1,z1,(strs(i1,j), j=1,6)
                      write(iout,81) x2,y2,z2,(strs(j1,j), j=1,6)
                      write(iout,81) x3,y3,z3,(strs(k1,j), j=1,6)
                        write(iout,81) x1,y1,z1,(strs(i1,j), j=1,6)
                        write(iout,81) x2,y2,z2,(strs(j1,j), j=1,6)
                        write(iout,81) x4,y4,z4,(strs(m1,j), j=1,6)
                      write(iout,81) x2,y2,z2,(strs(j1,j), j=1,6)
                      write(iout,81) x3,y3,z3,(strs(k1,j), j=1,6)
                      write(iout,81) x4,y4,z4,(strs(m1,j), j=1,6)
                        write(iout,81) x3,y3,z3,(strs(k1,j), j=1,6)
                        write(iout,81) x1,y1,z1,(strs(i1,j), j=1,6)
                        write(iout,81) x4,y4,z4,(strs(m1,j), j=1,6)
c
         enddo
c
 81      format(1x,20(g12.6,1x))
c
 999     continue
         close(itmp)
      return
      end
c
c
c     CONTOUR nodal values
      subroutine contours_HEX_D(strn,strs, npijkm,maxelem,
     &                    xord,yord,zord,
     &                    dispful,ioutput              )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21),icon(6,4),nn(8)
         real*8  strs(nnp,6), strn(nnp ,6)
         real*8  xord(nnp), yord(nnp), zord(nnp)
         real*8  dispful(nnp,3)
         real*8  sav(6),eav(6),x(12),y(12),z(12)
         data icon / 1,1,2,1,1,1
     &              ,2,2,7,4,8,7
     &              ,3,7,5,3,7,3
     &              ,7,5,6,8,5,8 /
c
ctag31
         rewind(iout)
         write(iout,*) nel,gnum1
         icont=ioutput
c        displacements do not need averaging
         if (icont .eq. 131) then
             do i=1,nel
                neltype=npijkm(i,1)
                if (neltype .eq. 9) then
                    idiv=1
                elseif (neltype .eq. 21) then
                    idiv=2
                endif
                do ivv=1,idiv
                   if (ivv .eq. 1) then
                       do k=1,4
                          i1=npijkm(i,1+k)
                          x(k)=xord(i1)
                          y(k)=yord(i1)
                          z(k)=zord(i1)
                          nn(k)=i1
                       enddo
                       do k=1,4
c!!!                      i1=npijkm(i,1+12+k)
                          i1=npijkm(i,1+4 +k)
                          x(4+k)=xord(i1)
                          y(4+k)=yord(i1)
                          z(4+k)=zord(i1)
                          nn(4+k)=i1
                       enddo
                   elseif (ivv .eq. 2) then
                       do k=1,4
                          i1=npijkm(i,1+12+k)
                          x(0+k)=xord(i1)
                          y(0+k)=yord(i1)
                          z(0+k)=zord(i1)
                          nn(0+k)=i1
                       enddo
                       do k=1,4
                          i1=npijkm(i,1+4+k)
                          x(4+k)=xord(i1)
                          y(4+k)=yord(i1)
                          z(4+k)=zord(i1)
                          nn(4+k)=i1
                       enddo
*                      do k=1,4
*                         i1=npijkm(i,1+12+k)
*                         x(4+k)=xord(i1)
*                         y(4+k)=yord(i1)
*                         z(4+k)=zord(i1)
*                         nn(4+k)=i1
*                      enddo
                   endif
                       do itet=1,6   
                          n1=icon(itet,1)
                          n2=icon(itet,2)
                          n3=icon(itet,3)
                          n4=icon(itet,4)
                          x1=x(n1)
                          x2=x(n2)
                          x3=x(n3)
                          x4=x(n4)
                            y1=y(n1)
                            y2=y(n2)
                            y3=y(n3)
                            y4=y(n4)
                          z1=z(n1)
                          z2=z(n2)
                          z3=z(n3)
                          z4=z(n4)
                            i1=nn(n1)
                            j1=nn(n2)
                            k1=nn(n3)
                            m1=nn(n4)
c           write(ilog,*) i1,j1,k1,m1,i,itet
              write(iout,81) x1,y1,z1,(dispful(i1,j), j=1,3),' 0 0 0'
              write(iout,81) x2,y2,z2,(dispful(j1,j), j=1,3),' 0 0 0'
              write(iout,81) x3,y3,z3,(dispful(k1,j), j=1,3),' 0 0 0'
                write(iout,81) x1,y1,z1,(dispful(i1,j), j=1,3),' 0 0 0'
                write(iout,81) x2,y2,z2,(dispful(j1,j), j=1,3),' 0 0 0'
                write(iout,81) x4,y4,z4,(dispful(m1,j), j=1,3),' 0 0 0'
              write(iout,81) x2,y2,z2,(dispful(j1,j), j=1,3),' 0 0 0'
              write(iout,81) x3,y3,z3,(dispful(k1,j), j=1,3),' 0 0 0'
              write(iout,81) x4,y4,z4,(dispful(m1,j), j=1,3),' 0 0 0'
                write(iout,81) x3,y3,z3,(dispful(k1,j), j=1,3),' 0 0 0'
                write(iout,81) x1,y1,z1,(dispful(i1,j), j=1,3),' 0 0 0'
                write(iout,81) x4,y4,z4,(dispful(m1,j), j=1,3),' 0 0 0'
                       enddo
               enddo
             enddo
             return
         endif
c
c        OTHERWISE  do stresses or strains
c
         do i=1,nel
                neltype=npijkm(i,1)
                i1=npijkm(i,2)
                x1=xord(i1)
c
                do k=1,8
                   i1=npijkm(i,1+k)
                   x(k)=xord(i1)
                   y(k)=yord(i1)
                   z(k)=zord(i1)
                   nn(k)=i1
                enddo
                do itet=1,6
                   i1=icon(itet,1)
                   j1=icon(itet,2)
                   k1=icon(itet,3)
                   m1=icon(itet,4)
                   x1=x(i1)
                   x2=x(j1)
                   x3=x(k1)
                   x4=x(m1)
                     y1=y(i1)
                     y2=y(j1)
                     y3=y(k1)
                     y4=y(m1)
                   z1=z(i1)
                   z2=z(j1)
                   z3=z(k1)
                   z4=z(m1)
                     n1=nn(i1)
                     n2=nn(j1)
                     n3=nn(k1)
                     n4=nn(m1)
                      write(iout,81) x1,y1,z1,(strs(n1,j), j=1,6)
                      write(iout,81) x2,y2,z2,(strs(n2,j), j=1,6)
                      write(iout,81) x3,y3,z3,(strs(n3,j), j=1,6)
                        write(iout,81) x1,y1,z1,(strs(n1,j), j=1,6)
                        write(iout,81) x2,y2,z2,(strs(n2,j), j=1,6)
                        write(iout,81) x4,y4,z4,(strs(n4,j), j=1,6)
                      write(iout,81) x2,y2,z2,(strs(n2,j), j=1,6)
                      write(iout,81) x3,y3,z3,(strs(n3,j), j=1,6)
                      write(iout,81) x4,y4,z4,(strs(n4,j), j=1,6)
                        write(iout,81) x3,y3,z3,(strs(n3,j), j=1,6)
                        write(iout,81) x1,y1,z1,(strs(n1,j), j=1,6)
                        write(iout,81) x4,y4,z4,(strs(n4,j), j=1,6)
                enddo
c
         enddo
c
 81      format(1x,20(g12.6,1x))
c
 999     continue
         close(itmp)
      return
      end
c
c
c     SHAPES  nodal values
      subroutine shapes_HEX_D(npijkm,maxelem,
     &                    xord,yord,zord,
     &                    dispful,evout,ioutput              )
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
         integer npijkm(maxelem,21),icon(12,3),icon8(6,3),nn(8)
         integer icon4(3,3)
*        real*8  strs(nnp,6), strn(nnp ,6)
         real*8  xord(nnp), yord(nnp), zord(nnp)
         real*8  dispful(nnp,3)
         real*8  sav(6),eav(6),x(12),y(12),z(12)
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
         rewind(iout)
         write(iout,*) nel*12,evout
         icont=ioutput
c        displacements do not need averaging
         if (icont .eq. 141 .OR. icont .eq. 31) then
             do n=1,nel
                neltype=npijkm(n,1)
                if (neltype .eq. 21) then
c                   20-noded hex
                    do j=1,12
                       do k=1,3
                          nodek=npijkm(n,1+icon(j,k))
                          x1  =xord(nodek)
                          y1  =yord(nodek)
                          z1  =zord(nodek)
                          write(iout,81) x1,y1,z1,
     &                          (dispful(nodek,kk), kk=1,3),' 0 0 0',n
                       enddo
                    enddo
                elseif (neltype .eq. 9) then
c                   8-noded hex
                    do j=1,6
                       do k=1,3
                          nodek=npijkm(n,1+icon8(j,k))
                          x1  =xord(nodek)
                          y1  =yord(nodek)
                          z1  =zord(nodek)
                          write(iout,81) x1,y1,z1,
     &                          (dispful(nodek,kk), kk=1,3),' 0 0 0',n
                       enddo
                    enddo
                elseif (neltype .eq. 5) then
c                   4-noded tetra
                    do j=1,3
                       do k=1,3
                          nodek=npijkm(n,1+icon4(j,k))
                          x1  =xord(nodek)
                          y1  =yord(nodek)
                          z1  =zord(nodek)
                          write(iout,81) x1,y1,z1,
     &                          (dispful(nodek,kk), kk=1,3),' 0 0 0',n
                       enddo
                    enddo
                endif
             enddo
         endif
         return
c
 81      format(1x,6(g12.6,1x),1x,a,1x,i5)
 82      format(1x,20(g12.6,1x))
c
 999     continue
         close(itmp)
      return
      end
c
c
c 
c
c TRANS3D
c     This subroutine makes 3-D vector transformations.
c
      subroutine trns3dv_D (gloads,l,m,n,lloads,beta)
         implicit real*8 (a-h,o-z)
         real*8 gloads(12),r(3,3),lloads(12),lload3(3) ,gload3(3)
         real*8 m,n,l,d
c
         pi=4.0*atan(1.0)
         cb=cos(beta*pi/180)
         sb=sin(beta*pi/180)
         d=sqrt(1-n**2)
       if (d.lt.1e-3)then
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
          r(2,1)  = -(m*cb+l*n*sb)/d
          r(2,2)  =  (l*cb-m*n*sb)/d
          r(2,3)  =  d*sb
          r(3,1)  =  (m*sb-l*n*cb)/d
          r(3,2)  = -(l*sb+m*n*cb)/d
          r(3,3)  =  d*cb
       endif
c
c
       do 11 i=1,12,3
          gload3(1)=gloads(i)
          gload3(2)=gloads(i+1)
          gload3(3)=gloads(i+2)
             call matABd(r,gload3,lload3,3)
          lloads(i)=lload3(1)
          lloads(i+1)=lload3(2)
          lloads(i+2)=lload3(3)
11     continue
c
      return
      end
c
c
      subroutine convertFRM_D(disploc, xl,
     &                      zix0,ziy0,ziz0,a0,e0,g0,
     &                      strain,stress,force,iglobal )
c
         implicit real*8 (a-h,o-z)
            include 'commons.std'
         real*8  disploc(3,6)
         real*8 strain(18),stress(18)
         real*8 strainb(12),stressb(12),ub(12)
         real*8 force(18),forceb(12)
c
         real*8 ekb12(12,12)
c
             call elmstfTRS_D(xl,a0,e0,g0,ek6b)
 81                   format(1x,20(g12.6,1x))
c
                   do k=1,6
                      ub(  k)=disploc(1,k)
                      ub(6+k)=disploc(2,k)
                   enddo
c
c                  multiply  {load}=[ek]{delt}
                   do i = 1,12
                      forceb(i) = 0.0e0
                      do j = 1,12
*                        forceb(i) = forceb(i) + ekb12(i,j) * ub(j)
                      enddo   
                   enddo   
c                  nodal forces wrt member axes
                   do i=1,2
                      ii=(i-1)*6
                      sgn=(-1)**i
                      force(ii+1)=forceb(ii+1)
                      force(ii+2)=forceb(ii+2)
                      force(ii+3)=forceb(ii+3)
                        force(ii+4)=forceb(ii+4)
                        force(ii+5)=forceb(ii+5)
                        force(ii+6)=forceb(ii+6)
                   enddo
                   do i=1,2
                      ii=(i-1)*6
                      sgn=(-1)**i
                      stressb(ii+1)=sgn*force(ii+1)/a0
                      stressb(ii+2)=sgn*force(ii+2)/a0
                      stressb(ii+3)=sgn*force(ii+3)/a0
                        stressb(ii+4)=sgn*force(ii+4)/zix0
                        stressb(ii+5)=sgn*force(ii+5)/ziy0
                        stressb(ii+6)=sgn*force(ii+6)/ziz0
                   enddo
                   do i=1,4
                      ii=(i-1)*3
                      strainb(ii+1)=stressb(ii+1)/e0
                      strainb(ii+2)=stressb(ii+2)/g0
                      strainb(ii+3)=stressb(ii+3)/g0
                   enddo
                   do i=1,12
                      stress(i)=stressb(i)
                      strain(i)=strainb(i)
                   enddo
                   do i=13,18
                      stress(i)=stress(i-6) 
                      strain(i)=strain(i-6) 
                      force (i)=force (i-6) 
                   enddo
c
       return
       end
c
c
      subroutine get_hh_D(xb1,xb2,xb3,yb1,yb2,yb3,hh,bb,zl,alpha,eeh)
c        [H] 
         implicit real*8 (a-h,o-z)
         real*8 x(3), y(3)
         real*8 xc(3), yc(3), xm(3),ym(3)
         real*8 gg(9,9),hh(9,9),bb(3,9),p(9,3),zl(9,3)
         real*8 bh(9,9)
         real*8 area,area2,c
         real*8 a1j, a2j, a3j, b1j, b2j,b3j
         real*8 x0, y0, xi,yi
         real*8 cj, sj,dl,dx,dy
         integer  iperm(9)
         real*8   wk(9,9),sumd,eeh(9,9)
c
c
         x(1)=xb1
         x(2)=xb2
         x(3)=xb3
         y(1)=yb1
         y(2)=yb2
         y(3)=yb3
c
         area2 = (y(2) -y(1))*(x(1)-x(3)) - (x(2)-x(1))*(y(1)-y(3))
c
         x0 = (x(1)+x(2)+x(3))/3.0
         y0 = (y(1)+y(2)+y(3))/3.0
         area=0.5*area2
         c = 1.0/sqrt(area)
c
         xc(1) = c*(x(1)-x0)
         xc(2) = c*(x(2)-x0)
         xc(3) = c*(x(3)-x0)
           yc(1) = c*(y(1)-y0)
           yc(2) = c*(y(2)-y0)
           yc(3) = c*(y(3)-y0)
         xm(1) = 0.5*(xc(2)+xc(3))
         xm(2) = 0.5*(xc(3)+xc(1))
         xm(3) = 0.5*(xc(1)+xc(2))
           ym(1) = 0.5*(yc(2)+yc(3))
           ym(2) = 0.5*(yc(3)+yc(1))
           ym(3) = 0.5*(yc(1)+yc(2))
c
c        form G^T in GT and initialize HH to 0
c
         do 2500 i=1,3
            dx = xm(i) - xc(i)
            dy = ym(i) - yc(i)
            dl = sqrt(dx*dx + dy*dy)
            ci = dx/dl
            si = dy/dl
c
c           !!!a2j b2j different than paper
            a1i = -0.5*si*ci**2
            a2i =  0.5*ci**3
            b2i = -0.5*si**3
            b3i =  0.5*si**2*ci
            a3i = -(b2i + a1i + a1i)
            b1i = -(b3i + b3i + a2i)
c
            goto 100
c
c           partition as [r] [c] [h]
 100        continue
	    gg(2*i-1,1) = 1.
	    gg(2*i-1,2) = 0.
	    gg(2*i-1,3) =  -yc(i)
	    gg(2*i-1,4)=  xc(i)
	    gg(2*i-1,5)=  0.0
	    gg(2*i-1,6)=  yc(i)
c
	    gg(2*i-0,1) = 0.
	    gg(2*i-0,2) = 1.
	    gg(2*i-0,3)=  xc(i)
	    gg(2*i-0,4)=  0.0 
	    gg(2*i-0,5)=  yc(i)
	    gg(2*i-0,6)=  xc(i)
c
	    gg(6+i-0,1)=  0.0
	    gg(6+i-0,2)=  0.0
	    gg(6+i-0,3)=   c
	    gg(6+i-0,4)=  0.0
	    gg(6+i-0,5)=  0.0
	    gg(6+i-0,6)=  0.0
c
	    do j = 1,3
               dx = xm(j) - xc(j)
               dy = ym(j) - yc(j)
               dl = sqrt(dx*dx + dy*dy)
               cj = dx/dl
               sj = dy/dl
c
c              !!!a2j b2j different than paper
               a1j = -0.5*sj*cj**2
               a2j =  0.5*cj**3
               b2j = -0.5*sj**3
               b3j =  0.5*sj**2*cj
               a3j = -(b2j + a1j + a1j)
               b1j = -(b3j + b3j + a2j)
c
	       xi =  xc(i)
	       yi =  yc(i)
	       gg(2*i-1,j+6) =   a1j*xi*xi + 2.*a2j*xi*yi + a3j*yi*yi
	       gg(2*i-0,j+6) =   b1j*xi*xi + 2.*b2j*xi*yi + b3j*yi*yi
	       gg(6+i-0,j+6) =  -c*(cj*xi+sj*yi)
            enddo
 2500    continue
C
C	Factor G' and backsolve to obtain H_h
C	       Form physical stiffness and add to incoming SM
c
        do i=1,9
           do j=1,9
              hh(i,j)=gg(i,j)
           enddo
        enddo
        call ainver(hh,9,iperm,wk)
c       inverse returned in [gt]
c
                area=((yb1-yb2)*(xb3-xb2) + (xb2-xb1)*(yb3-yb2))/2.
c
         do i=1,3
            xi =  xc(i)
            yi =  yc(i)
	    do j = 1,3
               dx = xm(j) - xc(j)
               dy = ym(j) - yc(j)
               dl = sqrt(dx*dx + dy*dy)
               cj = dx/dl
               sj = dy/dl
c
c              !!!a2j b2j different than paper
               a1j = -0.5*sj*cj**2
               a2j =  0.5*cj**3
               b2j = -0.5*sj**3
               b3j =  0.5*sj**2*cj
               a3j = -(b2j + a1j + a1j)
               b1j = -(b3j + b3j + a2j)
c
	       bh(i*3-2,j) =   c*( 2*a1j*xi + 2*a2j*yi)
	       bh(i*3-1,j) =   c*( 2*b2j*xi + 2*b3j*yi)
	       bh(i*3-0,j) =  -c*( 4*b3j*xi + 4*a1j*yi)
            enddo
         enddo
c
        do i=1,9
           do j=1,9
              sumd=0.0
              do k=1,3
                 sumd = sumd + bh(i,k)*hh(6+k,j)
              enddo
              eeh(i,j)=sumd
           enddo
        enddo
c
c
cccccccccccc  LUMPing
         x21 = x(2) - x(1)
         x32 = x(3) - x(2)
         x13 = x(1) - x(3)
         x12 = -x21          
         x23 = -x32 
         x31 = -x13 
c
         y21 = y(2) - y(1)
         y32 = y(3) - y(2)
         y13 = y(1) - y(3)
         y12 = -y21          
         y23 = -y32 
         y31 = -y13 
c
         area2 = y21*x13 - x21*y13
*        write(*,*)'@@ AREA x 2: ',area2
         p(1,1) = y23
         p(2,1) = 0.0
         p(3,1) = y31
         p(4,1) = 0.0
         p(5,1) = y12
         p(6,1) = 0.0
           p(1,2) = 0.0
           p(2,2) = x32
           p(3,2) = 0.0
           p(4,2) = x13
           p(5,2) = 0.0
           p(6,2) = x21
         p(1,3) = x32
         p(2,3) = y23
         p(3,3) = x13
         p(4,3) = y31
         p(5,3) = x21
         p(6,3) = y12
         n=6
c
         if (alpha .ne. 0.0) then
             p(7,1) = y23 * (y13-y21)*alpha/6.0
             p(7,2) = x32 * (x31-x12)*alpha/6.0
             p(7,3) =       (x31*y13-x12*y21)*alpha/3.0
               p(8,1) = y31 * (y21-y32)*alpha/6.0
               p(8,2) = x13 * (x12-x23)*alpha/6.0
               p(8,3) =       (x12*y21-x23*y32)*alpha/3.0
             p(9,1) = y12 * (y32-y13)*alpha/6.0
             p(9,2) = x21 * (x23-x31)*alpha/6.0
             p(9,3) =       (x23*y32-x31*y13)*alpha/3.0
             n=9
         endif
         do i=1,9
            do j=1,3
               zl(i,j)=p(i,j)/2.0
            enddo
         enddo
c
c
 81             format(1x,20(g11.5,1x))
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     POST ANalysis of displacements to get stresses etc
      subroutine IP_stress ( disp,dispful, idbc, iglobal,
     &                   npijkm,maxnode,maxelem,xord,yord,zord,
*    &                   strn,strs,forc,e_elm,s_elm,
*    &                                  e_elm_ip,s_elm_ip,
     &                   prop,nmprop,ippp )
c
             implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npijkm(maxelem,21), idbc(maxnode,3)
         real*8  disp(neq  ),  dispful(nnp,3)
*        real*8  strs(nnp,6),strn(nnp,6), forc(nel,60)
*        real*8  e_elm(nel,6),s_elm(nel,6)
*        real*8  e_elm_ip(nel,27,6),s_elm_ip(nel,27,6)
         real*8  xord(nnp), yord(nnp), zord(nnp)
*        real*8  strain(24),stress(24), force(12)
*        real*8  strain(48),stress(48), force(24)
*        real*8  disploc(3,6), temp(6)
*        real*8  delu(24),nload(12),wk(20)
*        real*8  xyz(3,20),ssi(8,6),ssn(8,6)
         real*8  xyz(3,20)
         real*8  uvw(3,20),dd(6,6)
         real*8  strain_ip(6,27),stress_ip(6,27), force_n(3,20)
c
         real*8 prop(ippp,10)
         integer nmprop(nel)
c
         rewind(igeo)
         rewind (idis)
c
         do i=1,neq
            read(idis) disp(i)
         enddo
c          reassign displacements to dispful( )
           do 20 i=1,nnp
              do j=1,3
                 ieqnum=idbc(i,j)
                 if (ieqnum .eq. 0) then
                     dispful(i,j) = 0.0 
                 else
                     dispful(i,j) = disp(ieqnum)
                 endif
              enddo
 20        continue
c
ctag2
c          For each element, calculate the strain, stress at IP
           write(*,*)'@@ loop over elements ',nel
           write(*,*)' '
           nel51=nel/50+1
           call zzwrti(nel51)
           call zzwrt ('>')
           do 50 i=1,nel
              if (mod(i,50) .eq. 0) call zzwrt('.')
              neltype=npijkm(i,1)
              mat=nmprop(i)
                 e0=prop(mat,1)
                 g0=prop(mat,2)
                 r0=prop(mat,4)
              call dmat(e0,g0,dd)
c
              if (neltype .eq. 3) then
c                 frame
              elseif (neltype .eq. 4) then
c                 plate
              elseif (neltype .eq. 5) then
c                3-D solid TET
              elseif (neltype .eq. 9) then
c                3-D solid HEX8
              elseif (neltype .eq. 21) then
c                3-D solid HEX20
c
                 kmax=neltype-1
                 do k=1,kmax
                    node=npijkm(i,1+k)
                    xyz(1,k)=xord(node)
                    xyz(2,k)=yord(node)
                    xyz(3,k)=zord(node)
c
                    uvw(1,k)=dispful(node,1)
                    uvw(2,k)=dispful(node,2)
                    uvw(3,k)=dispful(node,3)
                 enddo
c
ctag2
*                write(iout,*) '@@ Elem: ',i
                 ired_int=3
                 call convert_HEX20(neltype,xyz,uvw,dd,
     &                              ired_int,
     &                              strain_ip, stress_ip,force_n )
*                write(iout,81) (stress_ip(j,1), j=1,6)
                 do jj=1,27
                      write(igeo) (stress_ip(ii,jj), ii=1,6)
                 enddo
c
c
              endif
c             bottom distinction between elements
 50        continue
c          end of loop over elements
c
 81        format(1x,20(g12.6,1x))
 82        format(1x,i5,20(g11.5,1x))
c
 999  continue
      return
      end
c
c
c     Nodal loads due to BODY FORCEs from GRAVity for 3D solids
      subroutine body_force_grav( iglobal,
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
c
         real*8 prop(ippp,10),dd(6,6)
         integer nmprop(nel)
c
         real*8  force(60),force8(3,20)
         real*8  gforce(neq),grav(3)
         real*8  uvw(3,20),xyz(3,20)
c
c
*        write(iout,*)' << in body_force_grav >>'
c
         do i=1,neq
            gforce(i)=0.0
         enddo
         grav(1)=gx
         grav(2)=gy
         grav(3)=gz
*        do i=1,3
*           write(*,*) (prop(i,j), j=1,3) 
*           write(*,*) grav(i)
*        enddo
c
*        rewind(igeo)
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
                rho0=prop(mat,4)
                 pl0=prop(mat,5)
*            write(*,*) e0,g0,rho0
c
c???            call dmat(e0,g0,dd)
c
              if (neltype .eq. 5) then
c                 tetra
                  do k=1,4 
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
*                    uvw(1,k)=dispful(node,1)
*                    uvw(2,k)=dispful(node,2)
*                    uvw(3,k)=dispful(node,3)
                  enddo
                 call stress_TET( xyz,uvw,dd,ss)  !Fangbao
*                    write(iout,*) '@@ stress '
*                    write(iout,82) (ss(jj), jj=1,6)
                 call bforce_TET( xyz,uvw,ss,force ) !Fangbao
*                    write(iout,*) '@@ force '
*                    write(iout,82) (force(jj), jj=1,12)
c
c!!!              replace with HEX
                 call assembFOR_HEX(gforce,force,idbc,maxnode !Fangbao
     &                          ,npijkm,maxelem,n)
c
              elseif (neltype .eq. 9) then
*                    write(iout,*) '@@ uvw '
                  do k=1,8
                     node=npijkm(n,1+k)
                     xyz(1,k)=xord0(node)
                     xyz(2,k)=yord0(node)
                     xyz(3,k)=zord0(node)
c
*                    uvw(1,k)=dispful(node,1)
*                    uvw(2,k)=dispful(node,2)
*                    uvw(3,k)=dispful(node,3)
*                       write(iout,82) (uvw(jj,k), jj=1,3),node*1.0
                  enddo
*                 call stress_HEX8( xyz,uvw,dd,stress)
*                    write(iout,*) '@@ stress '
*                    do ii=1,8
*                       write(iout,82) (stress(jj,ii), jj=1,6)
*                    enddo
*                 call bforce_HEX8( xyz,uvw,stress,force8 )
*                    write(iout,*) '@@ force '
*                 do ii=1,8
*                    write(iout,82) (force8(jj,ii), jj=1,3)
*                 enddo
*                 call assembFOR_HEX(gforce,force8,idbc,maxnode
*    &                          ,npijkm,maxelem,n)
*             write(iout,*)'@@ partial gFORCE: ',n
*             write(iout,86) (gforce(k),k=1,neq)
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
*                    uvw(1,k)=dispful(node,1)
*                    uvw(2,k)=dispful(node,2)
*                    uvw(3,k)=dispful(node,3)
                  enddo
ctag2
                  ired_int=2
                  ired_int=3
                  ired_int=2
                  ired_int=3
                  call bforce_HEX20(rho0,grav,neltype,ired_int, xyz,
     &                              force8)
*                    write(iout,*) '@@ force '
*                 do ii=1,20
*                    write(iout,82) (force8(jj,ii), jj=1,3)
*                 enddo
                  call assembFOR_HEX(gforce,force8,idbc,maxnode
     &                          ,npijkm,maxelem,n)
c
              endif
 50        continue
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
      subroutine bforce_HEX20(rho0,grav,neltype,ired_int, 
     &                           xyz, force
     &                        )
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),uvw(3,20),u(60),ud(9),ee(6),dd(6,6)
         real*8 force(3,20),ff(60),ss(6),df(60),grav(3)
c
         real*8 BE(6,60),Bd(9,60),hh(20),bb(6,60)
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
*        iout=23
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
c        nn=0
c        do i=1,kmax
c           do j=1,3
c              nn=nn+1
c              u(nn)=uvw(j,i)
*                  write(iout,*) nn,u(nn)
c           enddo
c        enddo
c
         do i=1,kmax*3
            ff(i) = 0.0 
         enddo
         ip_n=0
         nint=ired_int
*        nint=2
ctag2    nint=1
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
c                 call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
                  call stdm20(xyz,hh,bb,bd,det,ri,si,ti,nel,ilog )
c
ccccccccccccccccccccccccccccccccccccc
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
                  do 70 j=1,kmax
                     ff(j) = ff(j) + hh(j)*rho0*wt
                     df(j) =         hh(j)*rho0*wt
 70               continue
*                   write(iout,* ) ' part dF ',ip_n    
*                   write(iout,82) (df(j), j=1,20)
*                   write(iout,82) (ff(j), j=1,20)
 80      continue
c
         nn=0
         do j=1,kmax
            do i=1,3
               force(i,j)=ff(j)*grav(i)
            enddo
         enddo
c        check equil
         fx=0
         fy=0
         fz=0
         do j=1,kmax
            fx=fx+force(1,j)
            fy=fy+force(2,j)
            fz=fz+force(3,j)
         enddo
*                   write(iout,* ) ' Element body force '    
*                   write(iout,82) fx,fy,fz
c
 82        format(1x,70(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
