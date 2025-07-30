c
c     FORM STIFfness matrix  [K]  
      subroutine formstif( stf,npijkm, xord, yord,zord,
     &                     idbc,maxnode,maxelem, iglobal,
     &                     qms,
     &                     prop,nmprop,ippp,iprofv,nloc)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         character*6 str06
         character*40 str40
         integer npijkm(maxelem,21),idbc(maxnode,3)
         integer nmprop(nel)
         real*8  prop(ippp,10)
*        integer iprofv(neq), nloc(neq),idof_node(60),all_eqn(60)
         integer iprofv(neq), nloc(neq)
c
         real*8 xord(nnp), yord(nnp),zord(nnp)
         real*8 stf(maxstiff)
         real*8 ek12(12,12)
         real*8 ek24(24,24),xyz(3,20),ek60(60,60)
         real*8         qms(neq)
         real*8 bb(6,120),dd(6,6)
c
c
         write(*,*) 'CHOOSE STIFFness storage: '
         write(*,*)'         0=DEFault              >>MEMory only>>   '
         write(*,*)'         1=BINary [prof]        >>StaDyn.STF>> '
         write(*,*)'         2=ASCii  [prof]        >>StaDyn.OUT>> '
         write(*,*)'         3=ASCii  [NxB]         >>StaDyn.OUT>> '
         write(*,*)'         4=BINary [NxB]         >>named     >> '
         write(*,*)'         5=DEFault + load       >>StaDyn.OUT>> '
         write(*,*)'        11=BINary [prof+locn]   >>named     >> '
         call zzwrt(' -->   ')
         read(ikbd,*)  iecho
         write(ilog,*) iecho,'   :: 1=BIN 2=ASC'
c
c        initialize [K]  to zero
*        do i=1,maxstiff
*           stf(i)=0.0
*        enddo
*            imaj  = iglobal/10
*            ibase = iglobal-imaj* 10
*            if (ibase .eq. 2) then
*                temp = abs(po) + abs(px) + abs(py)
*                if (temp .gt. 1.0e-10)  ipress=1
*                   rewind(ilod)
*                   do i=1,neq
*                      read(ilod) qms(i)
*                   enddo
*            endif
c
c        form each element matrix, and assemble
         write(*,'(a)')'@@ '
         write(str06,'(1x,i4,a)') (nel/50+1),'>'
         call zzwrt(str06)
         do 50 i=1,nel
            if (mod(i,50) .eq. 0) call zzwrt('.')
            neltype=npijkm(i,1)
*           i1=npijkm(i,2)
*           j1=npijkm(i,3)
*           k1=npijkm(i,4)
*           m1=npijkm(i,5)
*                 x1=xord(i1)
*                 x2=xord(j1)
*                 x3=xord(k1)
*                 x4=xord(m1)
*                  y1=yord(i1)
*                  y2=yord(j1)
*                  y3=yord(k1)
*                  y4=yord(m1)
*                   z1=zord(i1)
*                   z2=zord(j1)
*                   z3=zord(k1)
*                   z4=zord(m1)
            mat=nmprop(i)
            if (neltype .eq. 3) then
c               truss
            elseif (neltype .eq. 4) then
c               plate
c
            elseif (neltype .eq. 5) then
c               3-D tet solid
                e0=prop(mat,1)
                g0=prop(mat,2)
                r0=prop(mat,4)
                call dmat(e0,g0,dd)
c
                do k=1,8
                   node=npijkm(i,1+k)
                   xyz(1,k)=xord(node)
                   xyz(2,k)=yord(node)
                   xyz(3,k)=zord(node)
                enddo
c
c               initialize stiff
*               do j=1,12
*                  do k=1,12
*                     ek12(j,k)=0.0 
*                  enddo
*               enddo
                call elmstf_TET(dd,
     &                                xyz,
     &                                ek12,bb)
*             write(iout,*)' elem: ',i
*             do ii=1,12
*                write(iout,82) (ek12(ii,j), j=1,12)
*             enddo
c             assemble stiffness
                ielm=i
              call assembCOL_TET(stf,ek12,idbc,maxnode,
     &                            npijkm,maxelem,ielm,nloc)
c
*           write(iout,*) '  band storage: before ',i
*           call storeBND_f(stf,maxstiff,iprofv,nloc,iout,neq,iband,qms)
*           write(iout,*) '  band storage: after ',i
*           call storeBND_f(stf,maxstiff,iprofv,nloc,iout,neq,iband,qms)
c
            elseif (neltype .eq. 9) then
c               3-D solid with hex_8
                e0=prop(mat,1)
                g0=prop(mat,2)
                r0=prop(mat,4)
                call dmat(e0,g0,dd)
*             write(iout,*)' [D] ',i
*             do ii=1,6
*                write(iout,82) (dd(ii,jj), jj=1,6)
*             enddo
                do k=1,8
                   node=npijkm(i,1+k)
                   xyz(1,k)=xord(node)
                   xyz(2,k)=yord(node)
                   xyz(3,k)=zord(node)
                enddo
c
                call elmstf_HEX8(dd,xyz, ek24)
*                write(iout,82) (ek24(1,j), j=1,24)
*                write(iout,82) (ek24(j,j), j=1,24)
*                write(iout,* ) '  ' 
*             write(iout,*)' elem: ',i
*             do ii=1,24
*                write(iout,82) (ek24(ii,j), j=1,24)
*             enddo
c
c               assemble stiffness
                ielm=i
*               do ii=1,nnp
*                  write(*,*) (idbc(ii,jj), jj=1,3),' idbc'
*               enddo
                call assembCOL_HEX8(stf,ek24, idbc,maxnode, 
     &                          npijkm,maxelem,ielm,nloc)
c
            elseif (neltype .eq. 21) then
c               3-D solid with hex_20
                e0=prop(mat,1)
                g0=prop(mat,2)
                r0=prop(mat,4)
                call dmat(e0,g0,dd)
*             write(iout,*)' [D] ',i
*             do ii=1,6
*                write(iout,82) (dd(ii,jj), jj=1,6)
*             enddo
                do k=1,20
                   node=npijkm(i,1+k)
                   xyz(1,k)=xord(node)
                   xyz(2,k)=yord(node)
                   xyz(3,k)=zord(node)
                enddo
*             write(iout,*)' [X] ',i
*             do ii=1,3
*                write(iout,82) (xyz(ii,jj), jj=1,20)
*             enddo
c
                ired_int=3
                ired_int=2
                ired_int=3
                ired_int=2
                ired_int=ri3
*               write(*,*)' i ',i,ri3,ired_int
                call elmstf_HEX20(dd,xyz,ired_int, ek60)
*                write(iout,82) (ek24(1,j), j=1,24)
*                write(iout,82) (ek24(j,j), j=1,24)
*                write(iout,* ) '  ' 
ctag1
*             if (i .eq. 1) then
*                 write(iout,*)' [X] ',i,ired_int
*                 do ii=1,3
*                    write(iout,82) (xyz(ii,jj), jj=1,20)
*                 enddo
*             write(iout,*)' elem: ',i
*             do ii=1,60
*                write(iout,82) (ek60(ii,j), j=1,60)
*             enddo
*             endif
*             stop
c
c               assemble stiffness
                ielm=i
*               do ii=1,nnp
*                  write(iout,*) ii,(idbc(ii,jj), jj=1,3),'   idbc'
*               enddo
c
              kmax  = neltype-1
*             do kk=1,kmax
*                nodek = npijkm(ielm,1+kk)
*                all_eqn( (kk-1)*3 + 1) = idbc(nodek,1)
*                all_eqn( (kk-1)*3 + 2) = idbc(nodek,2)
*                all_eqn( (kk-1)*3 + 3) = idbc(nodek,3)
*             enddo
*             write(iout,*) ' all equations: ' 
*             write(iout,83) (all_eqn(jj), jj=1,60)
c
*             write(iout,*)' [ek60] reduced: ',i
*             do ii=1,kmax*3
*                ieqni=all_eqn(ii)
*                if (ieqni .gt. 0) then
*                    nn=0
*                    do jj=1,kmax*3
*                       ieqnj=all_eqn(jj)
*                       if (ieqnj .gt. 0) then
*                           nn=nn+1
*                           qms(nn)=ek60(ii,jj)
*                       endif
*                    enddo
*                    write(iout,82) (qms(jj), jj=1,nn)
*                 endif
*             enddo
                call assembCOL_HEX20(stf,ek60, idbc,maxnode, 
     &                          npijkm,maxelem,ielm,nloc)
*               do ii=1,maxstiff
*                  write(iout,*) ii,stf(ii)
*               enddo
*               stop
c
c
           endif
c
 50      continue 
         write(*,'(a)') '@@ end elems'
c        END loop over all elements
*            do i=1,20
*               write(*,*) stf(i),' in formstif' 
*            enddo
c
c        STORE stiffness  matrix on disk in case
         if (iecho .eq. 2 .OR. iecho .eq. 9) then
             write(iout,'(a)')'STIFFNESS: COLumn form '
             do 11 i=1,neq  
                iloc=nloc(i)
                write(iout,22) (stf(iloc+j-1), j=1,iprofv(i))
 22                   format(1x,6(g13.6))
 11          continue
             write(iout,*)'loads '
             do i=1,neq
                write(iout,23) qms(i),i
 23             format(1x,g13.6,1x,i5 )
             enddo
         endif
         if (iecho .eq. 5) then
             write(iout,*)'loads '
             do i=1,neq
                write(iout,23) qms(i),i
             enddo
         endif
c
*        if (ibase .eq. 2) then
*             write(*,*)'@@ writing {load} to DISK'
*             rewind(ilod)
*             do i=1,neq
*                write(ilod) qms(i)
*             enddo
*        endif
         if (iecho .eq. 1 .OR. iecho .eq. 9) then
             call storeCOL(stf,maxstiff,iprofv,nloc,istf,neq,iband)
         elseif (iecho .eq. 3) then
            write(iout,'(a)') ' STIFFNESS: band storage'
            call storeBND_f(stf,maxstiff,iprofv,nloc,iout,neq,iband,qms)
         elseif (iecho .eq. 4) then
             call zzwrt('INPUT: store_filename --> ')
             read(ikbd,'(a)') str40 
!Fangbao             write(ilog,*) str40
            open(unit=itmp,file=fdn//adjustl(str40) ,form='unformatted')
             rewind itmp
             call storeBND(stf,maxstiff,iprofv,nloc,itmp,neq,iband,qms)
             close(itmp)
         elseif (iecho .eq. 11) then
             call zzwrt('INPUT: store_filename --> ')
             read(ikbd,'(a)') str40 
!Fangbao             write(ilog,*) str40
            open(unit=itmp,file=fdn//adjustl(str40) ,form='unformatted')
             rewind itmp
             call storeCOL_n(stf,maxstiff,iprofv,nloc,itmp,neq,iband)
             close(itmp)
         endif
c
 82                   format(1x,60(g13.6))
 83                   format(1x,12(i4,1x))
c
         write(*   ,*) '@@ neq iband',neq  ,iband
         write(*   ,*) '@@ FORMSTFF:   Formed  [K]  OK'
!Fangbao         write(ilog,*) '@@ FORMSTFF:   Formed  [K]  OK'
      return
      end
c
c
c     STORE COLumn oriented matrix
      subroutine storeCOL( aa,maxstiff, iprofv,nloc,ifile,neq,iband)
         implicit real*8 (a-h,o-z)
         integer  iprofv(neq), nloc(neq)
         real*8   aa(maxstiff )
c
         write(*,*)'@@ writing to DISK'
         rewind(ifile)
         write(ifile) neq, iband   
         do i=1,neq  
            iloc = nloc(i)
            imax = iprofv(i)
            write(ifile) ( aa(iloc+j-1), j=1,imax )
         enddo
c
      return
      end
c
c     STORE COLumn oriented matrix in Named file with self location
      subroutine storeCOL_n( aa,maxstiff, iprofv,nloc,ifile,neq,iband)
         implicit real*8 (a-h,o-z)
         integer  iprofv(neq), nloc(neq)
         real*8   aa(maxstiff )
c
         write(*,*)'@@ writing to DISK'
         rewind(ifile)
*        write(*,*) maxstiff,maxnode,neq, ' zzzz' 
         write(ifile) neq, iband   
         do i=1,neq  
            iloc = nloc(i)
            imax = iprofv(i)
            write(ifile) iloc,imax
            write(ifile) ( aa(iloc+j-1), j=1,imax )
         enddo
      return
      end
c
c
c     STORE BaNDed oriented matrix
      subroutine storeBND( aa,maxstiff, iprofv,nloc,ifile,neq,iband,qms)
         implicit real*8 (a-h,o-z)
         integer  iprofv(neq), nloc(neq)
         real*8   aa(maxstiff ),qms(neq)
c
         write(*,*)'@@ writing to DISK  '
***      rewind(ifile)
         write(ifile) neq, iband   
             do i=1,neq  
                do j=1,iband 
                   ii=i+j-1
                   jj=j
                   if (j .le. neq-i+1) then
                       iloc=nloc(ii)
                       joff=j
                       jmax=iprofv(ii)
                       if (j .le. jmax) then
                           qms(j) =  aa(iloc+joff-1)
                       else
                           qms(j) =  0.0
                       endif
                   else
                       qms(j)=0.0
                   endif
                enddo
                write(ifile) (qms(j), j=1,iband)
             enddo
c
      return
      end
c
c
c     STORE BaNDed oriented matrix
      subroutine storeBND_f(aa,maxstiff,iprofv,nloc,ifile,neq,iband,qms)
         implicit real*8 (a-h,o-z)
         integer  iprofv(neq), nloc(neq)
         real*8   aa(maxstiff ),qms(neq)
c
         write(*,*)'@@ writing to DISK  '
         write(ifile,*) neq, iband   
             do i=1,neq  
                do j=1,iband 
                   ii=i+j-1
                   jj=j
                   if (j .le. neq-i+1) then
                       iloc=nloc(ii)
                       joff=j
                       jmax=iprofv(ii)
                       if (j .le. jmax) then
                           qms(j) = aa(iloc+joff-1)
                       else
                           qms(j) =  0.0
                       endif
                   else
                       qms(j)=0.0
                   endif
                enddo
                write(ifile,81) (qms(j), j=1,iband)
             enddo
c
 81   format(1x,120(g12.6,1x))
      return
      end
c
c
c     ELeMent STiFfness for Constant Strain Triangle
      subroutine elmstfCST_D(e0,g0,tt0,pl0,x1,x2,x3,y1,y2,y3,ek)
c
         implicit real*8 (a-h,o-z)
         integer inew(9)
         real*8 ek(12,12),ekb(6,6),b(3,6),c(3,3),etemp(3,6)
c
         do i=1,12
            do j=1,12
               ek(i,j)=0.
            enddo   
         enddo   
         do i=1,6
            do j=1,6
               ekb(i,j)=0.
            enddo   
         enddo   
         do i=1,3
            do j=1,6
               etemp(i,j)=0.
            enddo   
         enddo   
c
         znu0=(e0/g0/2.0)-1.0
c        plane strain
         if (pl0 .lt. 0) xkap=3.-4.*znu0
c        plane stress 
         if (pl0 .gt. 0) xkap=(3.-znu0)/(1.+znu0)
c
c          {s}=[C]{e}
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
c          {e}=[B]{u}
              a=((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
         a2=a*2.0
         b(1,1)=(y2-y3)/a2
         b(1,2)=0.
         b(1,3)=(y3-y1)/a2
         b(1,4)=0.
         b(1,5)=(y1-y2)/a2
         b(1,6)=0.
         b(2,1)=0.
         b(2,2)=(x3-x2)/a2
         b(2,3)=0.
         b(2,4)=(x1-x3)/a2
         b(2,5)=0.
         b(2,6)=(x2-x1)/a2
         b(3,1)=b(2,2)
         b(3,2)=b(1,1)
         b(3,3)=b(2,4)
         b(3,4)=b(1,3)
         b(3,5)=b(2,6)
         b(3,6)=b(1,5)
c
c         [k] = At[B][C][B]
         att0=a*tt0
         do 30 i=1,3
            do 32 j=1,6
               do 34 k=1,3
                  etemp(i,j)=etemp(i,j)+c(i,k)*b(k,j)
34             continue
32          continue
30       continue
         do 40 i=1,6
            do 42 j=1,6
               do 44 k=1,3
                  ekb(i,j)=ekb(i,j)+b(k,i)*etemp(k,j)*att0
44             continue
42          continue
40       continue
c
c        assign [ekb] to [12x12]
         inew(1)=1
         inew(2)=2
         inew(3)=4
         inew(4)=5
         inew(5)=7 
         inew(6)=8 
         do i=1,6
            ii=inew(i)
            do j=1,6
               jj=inew(j)
               ek(ii,jj) = ekb(i,j)
            enddo
         enddo
c
c        add shear contribution
         b1=(y2-y3)/a2
         b2=(y3-y1)/a2
         b3=(y1-y2)/a2
         c1=(x3-x2)/a2
         c2=(x1-x3)/a2
         c3=(x2-x1)/a2
         b(1,1) = c1
         b(1,2) = c2
         b(1,3) = c3
           b(2,1) = b1
           b(2,2) = b2
           b(2,3) = b3
c
c         [k] = G [B^T][B] Vol
         gv0=a*tt0*g0
         do i=1,3
            do j=1,3
               ekb(i,j)=0.0
               do k=1,2
                  ekb(i,j) = ekb(i,j) + b(k,j)*b(k,j)*gv0
               enddo
            enddo
         enddo
c
c        assign [ekb] to [12x12]
         inew(1)=3
         inew(2)=6
         inew(3)=9
         do i=1,3
            ii=inew(i)
            do j=1,3
               jj=inew(j)
               ek(ii,jj) = ekb(i,j)
            enddo
         enddo
c
c
      return
      end
c
c
c     ELeMent STiFfness for TETrahedron
      subroutine elmstf_TET(dd,
     &                       xyz,
     &                       ek,bb)
c
         implicit real*8 (a-h,o-z)
         real*8 ek(12,12),aa(3,3),bb(6,12),dd(6,6),etemp(6,12)
         real*8 xyz(3,20)
c
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
            enddo
         enddo
*        do i=1,6
*           do j=1,6
*              dd(i,j) = 0.0
*           enddo
*        enddo
c
*        write(*,*)          x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
c        hookes law
*        znu0=(e0/g0/2.0)-1.0
*        dd1 = (1+znu0)*(1-2*znu0)
*        dd(1,1) = 1-znu0
*        dd(2,1) =   znu0
*        dd(3,1) =   znu0
*          dd(1,2) =   znu0
*          dd(2,2) = 1-znu0
*          dd(3,2) =   znu0
*            dd(1,3) =    znu0
*            dd(2,3) =    znu0
*            dd(3,3) = 1- znu0
*        dd(4,4) = (1.0-2.0*znu0)/2.0
*        dd(5,5) = (1.0-2.0*znu0)/2.0
*        dd(6,6) = (1.0-2.0*znu0)/2.0
*        do i=1,6
*           do j=1,6
*              dd(i,j) = dd(i,j)*e0/dd1 
*           enddo
*        enddo
*        write(*,*)' DD:  nu ',znu0
*             do i=1,6
*                write(*,82) (DD(i,j), j=1,6)
*             enddo
c
*        x14 = x1-x4
*        x24 = x2-x4
*        x34 = x3-x4
*          y14 = y1-y4
*          y24 = y2-y4
*          y34 = y3-y4
*            z14 = z1-z4
*            z24 = z2-z4
*            z34 = z3-z4
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
*        write(*,*)' detJ: ',detj
*        write(*,*)' AA: '
*             do i=1,3
*                write(*,82) (aa(i,j), j=1,3)
*             enddo
c          {e}=[B]{u}
         bb(1,1) = aa(1,1) 
         bb(4,1) = aa(2,1)
         bb(6,1) = aa(3,1)
           bb(2,2) = aa(2,1) 
           bb(4,2) = aa(1,1)
           bb(5,2) = aa(3,1)
             bb(3,3) = aa(3,1) 
             bb(5,3) = aa(2,1)
             bb(6,3) = aa(1,1)
         bb(1,4) = aa(1,2) 
         bb(4,4) = aa(2,2)
         bb(6,4) = aa(3,2)
           bb(2,5) = aa(2,2) 
           bb(4,5) = aa(1,2)
           bb(5,5) = aa(3,2)
             bb(3,6) = aa(3,2) 
             bb(5,6) = aa(2,2)
             bb(6,6) = aa(1,2)
         bb(1,7) = aa(1,3) 
         bb(4,7) = aa(2,3)
         bb(6,7) = aa(3,3)
           bb(2,8) = aa(2,3) 
           bb(4,8) = aa(1,3)
           bb(5,8) = aa(3,3)
             bb(3,9) = aa(3,3) 
             bb(5,9) = aa(2,3)
             bb(6,9) = aa(1,3)
                ab1 = aa(1,1)+aa(1,2)+aa(1,3)
                ab2 = aa(2,1)+aa(2,2)+aa(2,3)
                ab3 = aa(3,1)+aa(3,2)+aa(3,3)
         bb(1,10) = -ab1    
         bb(4,10) = -ab2   
         bb(6,10) = -ab3   
           bb(2,11) = -ab2    
           bb(4,11) = -ab1   
           bb(5,11) = -ab3   
             bb(3,12) = -ab3    
             bb(5,12) = -ab2   
             bb(6,12) = -ab1   
*        write(*,*)' BB: '
*             do i=1,6
*                write(*,82) (bb(i,j), j=1,12)
*             enddo
c
c         [k] = [B^T][D][B] vol
         do i=1,6
            do j=1,12
               do k=1,6
                  etemp(i,j)=etemp(i,j)+dd(i,k)*bb(k,j)
               enddo
            enddo
         enddo
         do i=1,12
            do j=1,12
               do k=1,6
                  ek(i,j)=ek(i,j)+bb(k,i)*etemp(k,j)
               enddo
            enddo
         enddo
         vol=abs(detj)/6.0
         do i=1,12
            do j=1,12
                  ek(i,j)=ek(i,j)*vol
            enddo
         enddo
*        write(*,*)' ek: '
*             do i=1,12
*                write(*,82) (ek(i,j), j=1,12)
*             enddo
*       stop
c
 82     format(1x,20(g12.6,1x))
      return
      end
c
c
c     ASSEMBle element matrices in COLumn form
      subroutine assembCOL( aa, a, jbc, i1, j1,k1,nloc)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer  jbc(nnp*6), idof(18)
         integer  nloc(neq)
c
         real*8   a(18,18)
         real*8   aa(maxstiff)
c
c        Set idof to posible DoF number of each nodes
         do 10 i= 1, 6
            idof(i)    = (i1-1)*6 + i
            idof(i+6)  = (j1-1)*6 + i
            idof(i+12) = (k1-1)*6 + i
 10      continue
c
c        Store the values for individual array in global array
         do 20 i= 1, 18
            ieqn1 = jbc(idof(i))
            if (ieqn1 .gt. 0) then
               do 30 j= i, 18
                  ieqn2 = jbc(idof(j))
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
c     ASSEMBle element matrices in COLumn form for 3-D TRUSS
      subroutine assembCOL_TRUSS( aa, a, idbc,maxnode,  i1,j1,nloc)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer  idbc(maxnode,3),nloc(neq)
         integer  ieqn(6)
c
         real*8   a(6,6)
         real*8   aa(maxstiff)
c
c        Store the values for individual array in global array
*        write(*,*) i1,j1
*        write(*,*) ' maxnode ',maxnode   
         do 10 i=1,3
            ieqn(i+0)  = idbc(i1,i)
            ieqn(i+3)  = idbc(j1,i)
 10      continue 
         do 20 i= 1,6      
            ieqn1 = ieqn(i)
            if (ieqn1 .gt. 0) then
               do 30 j= i,6      
                  ieqn2 = ieqn(j)
*           write(*,*) ieqn1,ieqn2
                  if (ieqn2 .gt. 0) then
                     if (ieqn1 .gt. ieqn2) then
                        jband= (ieqn1-ieqn2)+1
                        iloc = nloc(ieqn1) 
                        aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
*                           write(*,*)' > ',iloc,jband,iloc+jband-1
                     else
                        jband= (ieqn2-ieqn1)+1
                        iloc = nloc(ieqn2) 
                        aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
*                           write(*,*)' <= ',iloc,jband,iloc+jband-1
                     endif
                  endif
 30            continue
           endif
 20      continue
 82      format(1x,120(g12.6,1x))
c
      return
      end
c
c
c     ASSEMBle element matrices in COLumn form for TETrahedrons
      subroutine assembCOL_TET( aa, a, idbc,maxnode, 
     &                          npijkm,maxelem,ielm,nloc)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer  idbc(maxnode,3),npijkm(maxelem,21)
         integer  nloc(neq)
         integer  all_eqn(12)
c
         real*8   a(12,12)
         real*8   aa(maxstiff)
c
c        Set idof to posible DoF number of each nodes
         neltype = npijkm(ielm,1)
         kmax = neltype - 1
            do kk=1,kmax
               nodek = npijkm(ielm,1+kk) 
               all_eqn( (kk-1)*3 + 1) = idbc(nodek,1) 
               all_eqn( (kk-1)*3 + 2) = idbc(nodek,2) 
               all_eqn( (kk-1)*3 + 3) = idbc(nodek,3) 
            enddo
c
c        Store the values for individual array in global array
         nmax=kmax*3
*        if (m1 .eq. 0) nmax=9
*        if (k1 .eq. 0) nmax=6
*           write(iout,*) 'm1 k1 nmax ',m1,k1,nmax 
         do 20 i= 1, kmax*3
*           inn = (i-1)/3 + 1
*           ii=i-(inn-1)*3
*           ieqn1 = idbc(idof(i),ii)
            ieqn1 = all_eqn(i)
            if (ieqn1 .gt. 0) then
               do 30 j= i, kmax*3
*                 jnn = (j-1)/3 + 1
*                 jj=j-(jnn-1)*3
*                 ieqn2 = idbc(idof(j),jj)
                  ieqn2 = all_eqn(j)
*           write(iout,*) i,ieqn1,j,ieqn2,idof(j)
                  if (ieqn2 .gt. 0) then
                     if (ieqn1 .gt. ieqn2) then
                        jband= (ieqn1-ieqn2)+1
                        iloc = nloc(ieqn1) 
*           write(iout,*) ieqn1,ieqn2,jband,iloc,' > '
                        aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
c                       aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
                     else
                        jband= (ieqn2-ieqn1)+1
                        iloc = nloc(ieqn2) 
*           write(iout,*) ieqn1,ieqn2,jband,iloc,' else '
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
c     ASSEMBle element matrices in COLumn form for HEX 8 nodes 
      subroutine assembCOL_HEX8( aa, a, idbc,maxnode, 
     &                          npijkm,maxelem,ielm,nloc)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer  idof_node(24),idbc(maxnode,3),npijkm(maxelem,21)
         integer  nloc(neq)
         real*8   a(24,24)
         real*8   aa(maxstiff)
c
c        Set idof to posible DoF number of each nodes
            neltype = npijkm(ielm,1)
            kmax  = neltype-1
            do i=1,3
               do k=1,kmax
                  idof_node(i+(k-1)*3)=npijkm(ielm,1+k)
               enddo
            enddo
            nmax=kmax*3
*           write(iout,81) ielm,(idof_node(k), k=1,nmax)
 81         format(1x,40(i5,1x))
c
c        Store the values for individual array in global array
*           write(iout,*) 'm1 k1 nmax ',m1,k1,nmax 
         do 20 i= 1, nmax
            inn = (i-1)/3 + 1
            ii=i-(inn-1)*3
            ieqn1 = idbc(idof_node(i),ii)
            if (ieqn1 .gt. 0) then
               do 30 j= i, nmax
                  jnn = (j-1)/3 + 1
                  jj=j-(jnn-1)*3
                  ieqn2 = idbc(idof_node(j),jj)
*           write(iout,*) i,ieqn1,j,ieqn2,idof(j)
                  if (ieqn2 .gt. 0) then
                     if (ieqn1 .gt. ieqn2) then
                        jband= (ieqn1-ieqn2)+1
                        iloc = nloc(ieqn1) 
*           write(iout,*) ieqn1,ieqn2,jband,iloc,' > '
                        aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
c                       aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
                     else
                        jband= (ieqn2-ieqn1)+1
                        iloc = nloc(ieqn2) 
*           write(iout,*) ieqn1,ieqn2,jband,iloc,' else ',i,j
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
c     ASSEMBle element matrices in COLumn form for HEX 20 nodes
      subroutine assembCOL_HEX20( aa, a, idbc,maxnode, 
     &                          npijkm,maxelem,ielm,nloc)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer  idbc(maxnode,3),npijkm(maxelem,21)
         integer  nloc(neq)
         integer  all_eqn(60)
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
*           do i=1,3
*              do k=1,kmax
*                 idof_node(i+(k-1)*3)=npijkm(ielm,1+k)
*              enddo
*           enddo
*           nmax=kmax*3
*           write(iout,81) ielm,(idof_node(k), k=1,nmax)
 81         format(1x,61(i5,1x))
c
c        Store the values for individual array in global array
*           write(iout,*) 'm1 k1 nmax ',m1,k1,nmax 
         do 20 i= 1, kmax*3
            inn = (i-1)/3 + 1
            ii=i-(inn-1)*3
c           ieqn1 = idbc(idof_node(i),ii)
c                                  ^ !
c           ieqn1 = idbc(inn,ii)
            ieqn1 = all_eqn(i)
*           write(iout,*) i,idof_node(i),ii,ieqn1,inn,'     iii' 
            if (ieqn1 .gt. 0) then
               do 30 j= i, kmax*3
                  jnn = (j-1)/3 + 1
                  jj=j-(jnn-1)*3
c                 ieqn2 = idbc(idof_node(j),jj)
c                 ieqn2 = idbc(jnn,jj)
                  ieqn2 = all_eqn(j)
*           write(iout,*) i,j,ieqn1,ieqn2,idof_node(j),jj
                  if (ieqn2 .gt. 0) then
                     if (ieqn1 .gt. ieqn2) then
                        jband= (ieqn1-ieqn2)+1
                        iloc = nloc(ieqn1) 
*           write(iout,*) ieqn1,ieqn2,jband,iloc,' > '
                        aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
c                       aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
                     else
                        jband= (ieqn2-ieqn1)+1
                        iloc = nloc(ieqn2) 
*           write(iout,*) ieqn1,ieqn2,jband,iloc,' else ',i,j
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
c     ASSEMBle LUMped element matrices
      subroutine assembLUM( aa, a, jbc, i1, j1,k1)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer  jbc(nnp*6), idof(18)
c
         real*8   a(18,18)
         real*8   aa(neq)
c
c        Set idof to posible DoF number of each nodes
         do 10 i= 1, 6
            idof(i)    = (i1-1)*6 + i
            idof(i+6)  = (j1-1)*6 + i
            idof(i+12) = (k1-1)*6 + i
 10      continue
c
c        Store the values for individual array in global array
         do 20 i= 1, 18
            ieqn1 = jbc(idof(i))
            if (ieqn1 .gt. 0) then
                aa(ieqn1) = aa(ieqn1) + a(i,i)
            endif
 20      continue
c
      return
      end
c
c
c     ASSEMBle element matrices in upper symm band form
      subroutine assemb( aa, a, jbc, i1, j1,k1)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer  jbc(nnp*6), idof(18)
         real*8   a(18,18), aa(neq,iband)
c
c        Set idof to posible DoF number of each nodes
         do 10 i= 1, 6
            idof(i)    = (i1-1)*6 + i
            idof(i+6)  = (j1-1)*6 + i
            idof(i+12) = (k1-1)*6 + i
 10      continue
c
c        Store the values for individual array in global array
         do 20 i= 1, 18
            ieqn1 = jbc(idof(i))
            if (ieqn1 .gt. 0) then
               do 30 j= i, 18
                  ieqn2 = jbc(idof(j))
                  if (ieqn2 .gt. 0) then
                     if (ieqn1 .gt. ieqn2) then
                        jband= (ieqn1-ieqn2)+1
                        aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
                     else
                        jband= (ieqn2-ieqn1)+1
                        aa(ieqn1,jband) = aa(ieqn1,jband) + a(i,j)
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
c     ASSEMBle consistent FORce for pressurized flex plate
      subroutine assembFOR( aa, a, jbc, i1, j1,k1)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
         integer  jbc(nnp*6), idof(18)
         real*8   a(18), aa(neq  )
c
c        Set idof to posible DoF number of each nodes
         do 10 i= 1, 6
            idof(i)    = (i1-1)*6 + i
            idof(i+6)  = (j1-1)*6 + i
            idof(i+12) = (k1-1)*6 + i
 10      continue
c
c        Store the values for individual array in global array
         do 20 i= 1, 18
            ieqn1 = jbc(idof(i))
            if (ieqn1 .gt. 0) then
                aa(ieqn1) = aa(ieqn1) + a(i)
            endif
 20      continue
c
      return
      end
c
c
c
c SUBROUTINE ELMstf
c      calculates the element stiffness matrices.
       subroutine elmstfTRS_D( length,  area, emod,gmod,ek )
c
         implicit real*8 (a-h,o-z)
           real*8      area, length, emod      
           real*8   ek(12,12)
c
c          initialize all ek elements to zero
           do i=1,12
              do j=1,12
                 ek(i,j)=0.0
              enddo
           enddo
c
c          STIFFNESS matrix in local coordinates
           ek(1,1)   =   emod*area/length
           ek(2,2)   =   gmod*area/length
           ek(3,3)   =   gmod*area/length
           ek(4,4)   =   ek(1,1)
           ek(5,5)   =   ek(2,2)
           ek(6,6)   =   ek(3,3)
           ek(1,4)   =  -ek(1,1)
           ek(2,5)   =  -ek(2,2)
           ek(3,6)   =  -ek(3,3)
c
c          impose symmetry
           do i=1,12
              do j=i,12
                 ek(j,i) = ek(i,j)
              enddo
           enddo
c
       return 
       end
c
c
c     ELeMent STiFfness for HEXahedron 8 nodes
c     based on iso quad from Bathe  pp295-297
      subroutine elmstf_HEX8(dd,    xyz, ek)
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,8),ek(24,24)
c
         real*8 dd(6,6),bb(6,24)
         real*8 xg(4,4),wgt(4,4),db(6)
         real*8 dk(24,24)
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
                  call stdm8(xyz,bb,det,ri,si,ti,nel,ilog)
                  write(iout,*) ' B ',lx,ly,lz
                  do i=1,6
                     write(iout,82) (bb(i,j), j=1,24)
                  enddo
c
c                 add contib to stiff
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
                  write(iout,*) ' wt ',wt
                  do 70 j=1,24
                     do 40 k=1,6
                        db(k) = 0.0
                        do 40 kk=1,6
                           db(k) = db(k) + dd(k,kk)*bb(kk,j)
 40                  continue
                     do 60 i=j,24   
                           sum = 0.0
                           do kk=1,6
                              sum = sum + bb(kk,i)*db(kk)
                           enddo
 50                        continue
                           ek(i,j) = ek(i,j) + sum*wt
                           dk(i,j) =           sum*wt
 60                  continue
 70               continue
           write(iout,*) ' [dk] ',lx,ly,lz
           do i=1,24
              write(iout,82) (dk(i,j), j=1,24)
           enddo
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
      subroutine stdm8(xyz,bb,det,r,s,t,nel,ilog)
c
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,8),bb(6,24),p(3,8),xj(3,3),xji(3,3)
         real*8          aa(3,8)
*        real*8 h(8)
c        integer*4 indx(3)
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
*        call ainver(xj,3,indx,xji)
*                 write(*,*) ' J^-1 ge '   
*                 do i=1,3
*                    write(*,82) (xji(i,j), j=1,3)
*                 enddo
         if (det .gt. 0.00000001) goto 40
!Fangbao         write(ilog,*)'@@ ERROR at elem: ',nel
         stop
c
 40      continue
c        jacob inverse
c        deriv oper B
*        k3=0
*        do 60 k=1,8
*           k3 = k3 + 3
*           bb(1,k3-2) = 0.0
*           bb(1,k3-1) = 0.0
*           bb(1,k3-0) = 0.0
*             bb(2,k3-2) = 0.0
*             bb(2,k3-1) = 0.0
*             bb(2,k3-0) = 0.0
*               bb(3,k3-2) = 0.0
*               bb(3,k3-1) = 0.0
*               bb(3,k3-0) = 0.0
*           do 50 i=1,3
*              bb(1,k3-2) = bb(1,k3-2) + xji(1,i)*p(i,k)
*              bb(2,k3-1) = bb(2,k3-1) + xji(2,i)*p(i,k)
*              bb(3,k3-0) = bb(3,k3-0) + xji(3,i)*p(i,k)
*50         continue
*           bb(4,k3  ) = bb(2,k3-1)
*           bb(4,k3-1) = bb(3,k3  )
*             bb(5,k3  ) = bb(1,k3-2)
*             bb(5,k3-2) = bb(3,k3  )
*               bb(6,k3-1) = bb(1,k3-2)
*               bb(6,k3-2) = bb(2,k3-1)
*60      continue
*                 write(*,*) ' B1 '
*                 do i=1,6
*                    write(*,82) (bb(i,j), j=1,24)
*                 enddo
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
 82      format(1x,40(g12.6,1x))
      return
      end
c
c
c     ELeMent STiFfness for HEXahedron 20 nodes
c     based on iso quad from Bathe  pp295-297
      subroutine elmstf_HEX20(dd,xyz,ired_int, ek)
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),ek(60,60)
c
         real*8 dd(6,6),bb(6,60),bd(9,60),hh(20)
         real*8 xg(4,4),wgt(4,4),db(6)
         real*8 dk(60,60)
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
 20      continue
         do i=1,60
            do j=1,60
               ek(i,j)=0.0
               dk(i,j)=0.0
            enddo
         enddo
c
*          write(iout,*)' [x y z] '
*          do j=1,20
*             write(iout,82) (xyz(i,j), i=1,3),j*1.0
*          enddo
*          write(*,*) ' '
c
c        reduced integration
         ip_n=0
*        nint=2
*        nint=3
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
*               write(*,*) lx,ly,lz,ip_n,ri,si,ti   
c
c                 deriv oper and jacobian
c?
                                                            ilog=21
                  call stdm20(xyz,hh,bb,bd,det,ri,si,ti,nel,ilog)
*           stop
*                 write(iout,*) ' B ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (bb(i,j), j=1,24)
*                 enddo
c
c                 add contib to stiff
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
*                 write(iout,*) ' wt ',wt 
*                 write(iout,*) wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)
                  do 70 j=1,60
                     do 40 k=1,6
                        db(k) = 0.0
                        do kk=1,6
                           db(k) = db(k) + dd(k,kk)*bb(kk,j)
                        enddo
 40                  continue
*          write(iout,82) (db(kk), kk=1,6)  
*          write(iout,*) ' '   
                     do 60 i=j,60   
                           sum = 0.0
                           do kk=1,6
                              sum = sum + bb(kk,i)*db(kk)
*          write(iout,82) sum,bb(kk,i)*db(kk),kk,sum*wt
                           enddo
*          write(iout,82) sum,i,j,sum*wt
 50                        continue
                           ek(i,j) = ek(i,j) + sum*wt
                           dk(i,j) =           sum*wt
 60                  continue
 70               continue
*          write(iout,*) ' [dk] ',lx,ly,lz
*          do i=1,24
*             write(iout,82) (dk(j,1), j=1,60)
*             write(iout,82) (dk(j,j), j=1,12)
*          enddo
 80      continue
*          write(*,*) ' [k] ',lx,ly,lz
*          do i=1,24
*             write(*,82) (ek(i,j), j=1,24)
*          enddo
c
c        impose symmetry
         do j=1,60
            do i=j,60
               ek(j,i)=ek(i,j)
            enddo
         enddo
c    
 82        format(1x,60(g12.6,1x))
 84        format(1x,60(i4,1x))
      return
      end
c
c
      subroutine stdm20(xyz,h,bb,bd,det,r,s,t,nel,ilog)
c
         Implicit real*8 (a-h,o-z)
         real*8 xyz(3,20),bb(6,60),p(3,20),xj(3,3),xji(3,3),bd(9,60)
         real*8          aa(3,20)
         real*8 h(20)
c        integer*4 indx(3)
c
         kmax=20
         rp = 1.0 + r
         sp = 1.0 + s
         tp = 1.0 + t
         rm = 1.0 - r
         sm = 1.0 - s
         tm = 1.0 - t
c
c        interp fns
         h(1) = 0.125*rm*sm*tm*(-r-s-t-2)
         h(2) = 0.125*rp*sm*tm*(+r-s-t-2)
         h(3) = 0.125*rp*sp*tm*(+r+s-t-2)
         h(4) = 0.125*rm*sp*tm*(-r+s-t-2)
           h(5) = 0.125*rm*sm*tp*(-r-s+t-2)
           h(6) = 0.125*rp*sm*tp*(+r-s+t-2)
           h(7) = 0.125*rp*sp*tp*(+r+s+t-2)
           h(8) = 0.125*rm*sp*tp*(-r+s+t-2)
         h(9 ) = 0.25*(1-r*r)*sm*tm
         h(11) = 0.25*(1-r*r)*sp*tm
         h(17) = 0.25*(1-r*r)*sm*tp
         h(19) = 0.25*(1-r*r)*sp*tp
           h(10) = 0.25*rp*(1-s*s)*tm
           h(12) = 0.25*rm*(1-s*s)*tm
           h(18) = 0.25*rp*(1-s*s)*tp
           h(20) = 0.25*rm*(1-s*s)*tp
             h(13) = 0.25*rm*sm*(1-t*t)
             h(14) = 0.25*rp*sm*(1-t*t)
             h(15) = 0.25*rp*sp*(1-t*t)
             h(16) = 0.25*rm*sp*(1-t*t)
c
c        nat coords deriv wrt r
                  do i=1,3
                     do j=1,20
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
c?
         iout=23
*                 write(iout,*) '@@ r,s,t: '   
*                 write(iout,*) r,s,t   
*                 write(iout,*) ' [P] '   
*                 do i=1,3
*                    write(iout,82) (p(i,j), j=1,20)
*                 enddo
c
c        jacob at (r,s,t)
         do 30 i=1,3
            do 30 j=1,3
               dum = 0.0
               do 20 k=1,20
                  dum = dum + p(i,k)*xyz(j,k)
*           if (j .eq. 1) write(iout,*) i,k,dum , p(i,k)*(xyz(j,k)-0.0)
 20            continue
               xj(i,j)=dum
 30      continue
*                 write(iout,*) '[J] '   
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
*                 write(iout,*) ' [J^-1] explicit '   
*                 do i=1,3
*                    write(iout,82) (xji(i,j), j=1,3)
*                 enddo
*        call ainver(xj,3,indx,xji)
*                 write(*,*) ' J^-1 ge '   
*                 do i=1,3
*                    write(*,82) (xji(i,j), j=1,3)
*                 enddo
*        if (det .gt. 0.00000001) goto 40
         if (det .gt. 1.0D-16    ) goto 40
!Fangbao             write(ilog,*)'@@ ERROR at elem: ',nel
!Fangbao             write(ilog,*)'@@ det = ',det
         stop
c
 40      continue
c        jacob inverse
c        deriv oper B
*        k3=0
*        do 60 k=1,8
*           k3 = k3 + 3
*           bb(1,k3-2) = 0.0
*           bb(1,k3-1) = 0.0
*           bb(1,k3-0) = 0.0
*             bb(2,k3-2) = 0.0
*             bb(2,k3-1) = 0.0
*             bb(2,k3-0) = 0.0
*               bb(3,k3-2) = 0.0
*               bb(3,k3-1) = 0.0
*               bb(3,k3-0) = 0.0
*           do 50 i=1,3
*              bb(1,k3-2) = bb(1,k3-2) + xji(1,i)*p(i,k)
*              bb(2,k3-1) = bb(2,k3-1) + xji(2,i)*p(i,k)
*              bb(3,k3-0) = bb(3,k3-0) + xji(3,i)*p(i,k)
*50         continue
*           bb(4,k3  ) = bb(2,k3-1)
*           bb(4,k3-1) = bb(3,k3  )
*             bb(5,k3  ) = bb(1,k3-2)
*             bb(5,k3-2) = bb(3,k3  )
*               bb(6,k3-1) = bb(1,k3-2)
*               bb(6,k3-2) = bb(2,k3-1)
*60      continue
*                 write(*,*) ' B1 '
*                 do i=1,6
*                    write(*,82) (bb(i,j), j=1,24)
*                 enddo
         do i=1,3
            do j=1,20
               sum=0.0d0
               do k=1,3 
                  sum = sum + xji(i,k)*p(k,j)
               enddo
               aa(i,j)=sum
            enddo
         enddo
*                write(iout,*) ' [A] '
*                do i=1,3
*                   write(iout,82) (aa(i,j), j=1,20)
*                enddo
c
c        distribut according to strain
         do i=1,6
            do j=1,60
               bb(i,j) = 0.0d0
            enddo
         enddo
         k3=0
         do k=1,20
            k3 = k3 + 3
            bb(1,k3-2) = aa(1,k)
            bb(2,k3-1) = aa(2,k)
            bb(3,k3-0) = aa(3,k)
                  bb(4,k3-2) = aa(2,k)
                  bb(4,k3-1) = aa(1,k)
              bb(5,k3-1) = aa(3,k)
              bb(5,k3-0) = aa(2,k)
                bb(6,k3-2) = aa(3,k)
                bb(6,k3-0) = aa(1,k)
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
c    
 82      format(1x,60(g12.6,1x))
      return
      end
c
c
      subroutine dmat (e0,g0,dd)  
         implicit real*8 (a-h,o-z)
         real*8 dd(6,6)
c
         do i=1,6
            do j=1,6
               dd(i,j)=0.0
            enddo
         enddo
         znu=(e0/g0/2.0)-1.0
         c1 = e0/( (1.0+znu)*(1.0-2*znu) )
         c2 = c1*(0.5-znu)
         dd(1,1) = c1*(1.0-znu)
         dd(1,2) = c1*(    znu)
         dd(1,3) = dd(1,2)
           dd(2,1) = dd(1,2)
           dd(2,2) = dd(1,1)
           dd(2,3) = dd(1,2)
             dd(3,1) = dd(1,2)
             dd(3,2) = dd(1,2)
             dd(3,3) = dd(1,1)
               dd(4,4) = c2      
               dd(5,5) = c2      
               dd(6,6) = c2      
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     FORM GEOMetric stiffness matrix  [K_G]  
      subroutine formgeom( stf,geo,npijkm, xord, yord,zord,
     &                     idbc,maxnode,maxelem, iglobal,
     &                     disp,dispful,wk,
     &                     prop,nmprop,ippp,iprofv,nloc)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         character*6 str06
         character*40 str40
         integer npijkm(maxelem,21),idbc(maxnode,3)
         integer nmprop(nel)
         real*8  prop(ippp,10)
         integer iprofv(neq), nloc(neq)
c
         real*8 xord(nnp), yord(nnp),zord(nnp)
         real*8 stf(maxstiff),geo(maxstiff)
         real*8 eke60(60,60),ekg60(60,60)
         real*8 disp(neq),dispful(maxnode,3),  wk(neq)
         real*8 dd(6,6)
         real*8 strain_ip(6,27),stress_ip(6,27),force_n(3,20)
         real*8 uvw(3,20),xyz(3,20)
c
c
         write(*,*) 'CHOOSE STIFFness storage: '
         write(*,*)'         0=DEFault[prof]        >>StaDyn.GEO>> '
         write(*,*)'         2=ASCii  [prof]        >>StaDyn.OUT>> '
         write(*,*)'         3=ASCii  [NxB]         >>StaDyn.OUT>> '
         write(*,*)'         4=BINary [NxB]         >>named     >> '
         call zzwrt(' -->   ')
         read(ikbd,*)  iecho
         write(ilog,*) iecho,'   :: 1=BIN 2=ASC'
c
c        initialize [K_G]  to zero
         do i=1,maxstiff
            stf(i)=0.0     !why recompute [K_E] ??
            geo(i)=0.0
         enddo
         rewind(igeo)
c
c          reassign displacements to dispful( )
           rewind(idis)
           do i=1,neq
              read(idis) disp(i)
           enddo
*          do i=1,20 
*             write(*,*) disp(i),' u'
*          enddo
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
c
c        form each element matrix, and assemble
         write(*,'(a)')'@@ '
         write(str06,'(1x,i4,a)') (nel/50+1),'>'
         call zzwrt(str06)
         do 50 i=1,nel
            if (mod(i,50) .eq. 0) call zzwrt('.')
            neltype=npijkm(i,1)
            mat=nmprop(i)
                e0=prop(mat,1)
                g0=prop(mat,2)
                r0=prop(mat,4)
                call dmat(e0,g0,dd)
            if (neltype .eq. 3) then
c               truss
            elseif (neltype .eq. 4) then
c               plate
c
            elseif (neltype .eq. 5) then
c               3-D tet solid
c
            elseif (neltype .eq. 9) then
c               3-D solid with hex_8
c
            elseif (neltype .eq. 21) then
c               3-D solid with Hex20
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
*               ired_int=3
                ired_int=ri3
                inonlinear=0
                call stress_HEX20(neltype,ired_int,xyz,uvw,dd,
     &                        stress_ip,inonlinear)
*               call convert_HEX20(neltype,xyz,uvw,dd,
*    &                        ired_int,strain_ip,stress_ip,force_n)
*                if (i .eq. 1) then
*                    do jj=1,27
*                         write(*,*) (stress_ip(ii,jj), ii=1,6)
*                    enddo
*                 endif
c
c               ensure linear formulation
                do k=1,kmax
                     uvw(1,k)=0.0
                     uvw(2,k)=0.0
                     uvw(3,k)=0.0
                enddo
                call stiffE_HEX20(neltype,ired_int,dd,xyz,uvw,eke60)
*               call elmgeo_HEX20(neltype,ired_int,xyz,stress_ip,ekg60)
             call stiffG_HEX20(neltype,ired_int,xyz,uvw,stress_ip,ekg60)
ctag1
*             if (i .eq. 1) then
*                 write(iout,*)' [X] ',i,ired_int
*                 do ii=1,3
*                    write(iout,82) (xyz(ii,jj), jj=1,20)
*                 enddo
*             write(iout,*)' E elem: ',i
*             do ii=1,60
*                write(iout,82) (eke60(ii,j), j=1,60)
*             enddo
*               call elmstf_HEX20(dd,xyz,ired_int, eke60)
*             write(iout,*)' E elem 2: ',i
*             do ii=1,60
*                write(iout,82) (eke60(ii,j), j=1,60)
*             enddo
*             endif
*          stop
*             write(iout,*)' G elem: ',i
*             do ii=1,60
*                write(iout,82) (ekg60(ii,j), j=1,60)
*             enddo
*             write(iout,*)' E & G diag: ',i
*             do ii=1,60
*                write(iout,82) eke60(ii,ii), ekg60(ii,ii)
*             enddo
*             stop
c
c               assemble stiffness
                ielm=i
                call assembCOL_HEX20(stf,eke60, idbc,maxnode, 
     &                          npijkm,maxelem,ielm,nloc)
                call assembCOL_HEX20(geo,ekg60, idbc,maxnode, 
     &                          npijkm,maxelem,ielm,nloc)
*               do ii=1,maxstiff
*                  write(iout,*) ii,stf(ii)
*               enddo
c
           endif
c
 50      continue 
         write(*,'(a)') '@@ end elems'
*            do i=1,20
*               write(*,*) stf(i),' in formgeom' 
*            enddo
c        END loop over all elements
c        do ii=1,maxstiff
c           stf(ii)=stf(ii)+geo(ii)
c        enddo
c???  why sum now
*              do ii=1,neq
*                 write(iout,*) ii,geo(nloc(ii))
*              enddo
c
c        STORE stiffness  matrix on disk in case
         rewind(igeo)
         write(igeo) neq,iband
         call storeCOL(geo,maxstiff,iprofv,nloc,igeo,neq,iband)
c
         if (iecho .eq. 2 .OR. iecho .eq. 9) then
             write(iout,'(a)')' GEOM_STIFFNESS: COLumn form '
             do 11 i=1,neq  
                iloc=nloc(i)
                write(iout,22) (geo(iloc+j-1), j=1,iprofv(i))
 22                   format(1x,6(g13.6))
 11          continue
*            write(iout,*)'loads '
*            do i=1,neq
*               write(iout,23) qms(i),i
 23             format(1x,g13.6,1x,i5 )
*            enddo
         endif
c
         if (iecho .eq. 1 .OR. iecho .eq. 9) then
*            call storeCOL(geo,maxstiff,iprofv,nloc,igeo,neq,iband)
c
         elseif (iecho .eq. 3) then
            write(iout,'(a)') ' E_stiffness: band storage'
            call storeBND_f(stf,maxstiff,iprofv,nloc,iout,neq,iband,wk)
            write(iout,'(a)') ' G_stiffness: band storage'
            call storeBND_f(geo,maxstiff,iprofv,nloc,iout,neq,iband,wk)
c
         elseif (iecho .eq. 4) then
             call zzwrt('INPUT: store_filename --> ')
             read(ikbd,'(a)') str40 
!Fangbao             write(ilog,*) str40
            open(unit=itmp,file=fdn//adjustl(str40) ,form='unformatted')
             rewind itmp
             call storeBND(geo,maxstiff,iprofv,nloc,itmp,neq,iband,wk)
             close(itmp)
         endif
c
 82                   format(1x,60(g13.6))
 83                   format(1x,12(i4,1x))
c
         write(*   ,*) '@@ neq iband',neq  ,iband
!Fangbao         write(*   ,*) '@@ FORMSTFF:   Formed  [K]  OK'
!Fangbao         write(ilog,*) '@@ FORMSTFF:   Formed  [K]  OK'
      return
      end
c
c
c     STIFFness Geometric for HEXahedron 8 nodes
c     based on iso quad from Bathe  pp295-297
      subroutine elmgeo_HEX20(neltype,ired_int, xyz,stress, ekg)
         Implicit real*8 (a-h,o-z)
c
         real*8 xyz(3,20),ekg(60,60)
         real*8 BB(6,60),bd(9,60),stress(6,27),hh(20)
         real*8 sm(9,9),sm0(3,3)
c
         real*8 xg(4,4),wgt(4,4),db(9)
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
c                 call BEmat20(neltype,xyz,uvw,ri,si,ti,BE,Bd,det)
                  call stdm20(xyz,hh,bb,bd,det,ri,si,ti,nel,ilog)
*                 write(iout,*) ' BB ',lx,ly,lz
*                 do i=1,6
*                    write(iout,82) (BB(i,j), j=1,kmax*3)
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
 60                  continue
 70               continue
 80      continue
c
c        impose symmetry
         do j=1,kmax*3
            do i=j,kmax*3
               ekg(j,i)=ekg(i,j)
            enddo
         enddo
c    
 82        format(1x,40(g12.6,1x))
 84        format(1x,40(i4,1x))
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     FORM GEOMetric stiffness matrix  [K_G]  
      subroutine formmass( mss,npijkm, xord, yord,zord,
     &                     idbc,maxnode,maxelem, iglobal,
     &                     wk,
     &                     prop,nmprop,ippp,iprofv,nloc)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         character*6 str06
         character*40 str40
         integer npijkm(maxelem,21),idbc(maxnode,3)
         integer nmprop(nel)
         real*8  prop(ippp,10)
         integer iprofv(neq), nloc(neq)
c
         real*8 xord(nnp), yord(nnp),zord(nnp)
         real*8 mss(maxmass )
         real*8 emm60(60,60)
         real*8 wk(neq),zmmd20(20),ff(3,20)
*        real*8 dd(6,6)
*        real*8 strain_ip(6,27),stress_ip(6,27),force_n(3,20)
         real*8 uvw(3,20),xyz(3,20)
c
c
         write(*,*) 'CHOOSE MASS storage: '
         write(*,*)'         0=DEFault[prof]        >>StaDyn.GEO>> '
         write(*,*)'         2=ASCii  [prof]        >>StaDyn.OUT>> '
         write(*,*)'         3=ASCii  [NxB]         >>StaDyn.OUT>> '
         write(*,*)'         4=BINary [NxB]         >>named     >> '
         call zzwrt(' -->   ')
         read(ikbd,*)  iecho
	   write(*,*)iecho
         write(ilog,*) iecho,'   :: 1=BIN 2=ASC'
c
c        initialize [M]  to zero
!Fangbao         write(ilog,*)'@@ in form mass, ilump: ',ilump
         do i=1,maxmass 
            mss(i)=0.0
         enddo
         rewind(imss)
c
c        form each element matrix, and assemble
         write(*,'(a)')'@@ '
         write(str06,'(1x,i4,a)') (nel/50+1),'>'
         call zzwrt(str06)
         do 50 i=1,nel
            if (mod(i,50) .eq. 0) call zzwrt('.')
            neltype=npijkm(i,1)
            mat=nmprop(i)
!               e0=prop(mat,1)
!               g0=prop(mat,2)
                r0=prop(mat,4) !Ye,solid density from qed file
!               call dmat(e0,g0,dd)
            if (neltype .eq. 3) then
c               truss
            elseif (neltype .eq. 4) then
c               plate
c
            elseif (neltype .eq. 5) then
c               3-D tet solid
c
            elseif (neltype .eq. 9) then
c               3-D solid with hex_8
c
            elseif (neltype .eq. 21) then
c               3-D solid with hex_20
                kmax=neltype-1
                do k=1,kmax
                   node=npijkm(i,1+k) !Ye, # of 20 nodes in one element
                   xyz(1,k)=xord(node)
                   xyz(2,k)=yord(node)
                   xyz(3,k)=zord(node)
                enddo
c
                ired_int=3
c               ired_int=2
                call elmmas_HEX20(neltype,ired_int,xyz,r0,emm60,
     &                            ilump,zmmd20)
*             write(iout,*)' m elem: ',i
*             do ii=1,60
*                write(iout,82) (emm60(ii,j), j=1,60)
*             enddo
*             do ii=1,20
*                write(iout,82)  zmmd20(ii)
*             enddo
*             write(iout,*)' m diag: ',i
*             do ii=1,60
*                write(iout,82) emm60(ii,ii)
*             enddo
*             stop
c
c               assemble stiffness
                if (ilump .eq. 1) then
                    do ii=1,3
                       do jj=1,kmax
                          ff(ii,jj)=zmmd20(jj)
                       enddo
                     enddo
                    ielm=i
                    call assembFOR_HEX(wk,ff, idbc,maxnode, 
     &                                   npijkm,maxelem,ielm)
                else
                    ielm=i
                    call assembCOL_HEX20(mss,emm60, idbc,maxnode, 
     &                                   npijkm,maxelem,ielm,nloc)
                endif
*               do ii=1,maxstiff
*                  write(iout,*) ii,stf(ii)
*               enddo
c
           endif
c
 50      continue 
         write(*,'(a)') '@@ end elems'
c        END loop over all elements
c
*              do ii=1,neq
*                 write(iout,*) ii,geo(nloc(ii))
*              enddo
c
c        STORE mass  matrix on disk in case
                if (ilump .eq. 1) then
                    ibandm=1
                    do ii=1,neq
                       mss(ii)=wk(ii)
                    enddo
                endif
c      modify mass matrix for concentrated masses along diagonal
        rewind(icms)
        do n=1,neq
           read(icms) temp8
           if (ilump.eq.1) then    
               mss(n) = mss(n) + temp8
           elseif (ilump.eq.2) then    
               iloc=nloc(n)
               mss(iloc) = mss(iloc) + temp8
           endif
        enddo
!Fangbao        write(*   ,*) '@@ Reloaded  {concmas}  successfully'
!Fangbao        write(ilog,*) '@@ Reloaded  {concmas}  successfully'
c
         rewind(imss)
*        write(*,*) maxstiff,maxnode,neq,' yyy'
         if (ilump .eq. 1) then
             write(imss) neq,ibandm
             do i=1,neq
                write(imss) wk(i)
             enddo
         else
             write(imss) neq,iband
             call storeCOL(mss,maxstiff,iprofv,nloc,imss,neq,iband)
         endif
c
         if (iecho .eq. 2 .OR. iecho .eq. 9) then
             write(iout,'(a)')' MASS : COLumn form '
             do 11 i=1,neq  
                iloc=nloc(i)
                write(iout,22) (mss(iloc+j-1), j=1,iprofv(i))
 22                   format(1x,6(g13.6))
 11          continue
*            write(iout,*)'loads '
*            do i=1,neq
*               write(iout,23) qms(i),i
 23             format(1x,g13.6,1x,i5 )
*            enddo
         endif
c
         if (iecho .eq. 1 .OR. iecho .eq. 9) then
*            call storeCOL(mss,maxmass,iprofv,nloc,igeo,neq,iband)
c
         elseif (iecho .eq. 3) then
            write(iout,'(a)') ' MASS: band storage'
            write(iout,* ) maxmass,neq,iband    
            call storeBND_f(mss,maxmass,iprofv,nloc,iout,neq,iband,wk)
c
         elseif (iecho .eq. 4) then
             call zzwrt('INPUT: store_filename --> ')
             read(ikbd,'(a)') str40 
!Fangbao             write(ilog,*) str40
            open(unit=itmp,file=fdn//adjustl(str40) ,form='unformatted')
             rewind itmp
             call storeBND(mss,maxmass,iprofv,nloc,itmp,neq,iband,wk)
             close(itmp)
         endif
c
 82                   format(1x,60(g13.6))
 83                   format(1x,12(i4,1x))
c
         write(*   ,*) ' @@ neq iband',neq  ,iband
         write(*   ,*) ' @@ FORMMASS:   Formed  [M]  OK'
!Fangbao         write(ilog,*) ' @@ FORMMASS:   Formed  [M]  OK'
      return
      end
c
c
!====================================================
c     FORM GEOMetric stiffness matrix  [K_G]  
      subroutine formmass_Tian( mss,npijkm, xord, yord,zord,
     &                     idbc,maxnode,maxelem, iglobal,
     &                     wk,
     &                     prop,nmprop,ippp)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         character*6 str06
         character*40 str40
         integer npijkm(maxelem,21),idbc(maxnode,3)
         integer nmprop(nel)
         real*8  prop(ippp,10)
!         integer iprofv(neq), nloc(neq)
c
         real*8 xord(nnp), yord(nnp),zord(nnp)
         real*8 mss(maxmass )
         real*8 emm60(60,60)
         real*8 wk(neq),zmmd20(20),ff(3,20)
*        real*8 dd(6,6)
*        real*8 strain_ip(6,27),stress_ip(6,27),force_n(3,20)
         real*8 uvw(3,20),xyz(3,20)
c
c
         write(*,*) 'CHOOSE MASS storage: '
         write(*,*)'         0=DEFault[prof]        >>StaDyn.GEO>> '
         write(*,*)'         2=ASCii  [prof]        >>StaDyn.OUT>> '
         write(*,*)'         3=ASCii  [NxB]         >>StaDyn.OUT>> '
         write(*,*)'         4=BINary [NxB]         >>named     >> '
         call zzwrt(' -->   ')
         read(ikbd,*)  iecho
         write(*,*)iecho
         write(ilog,*) iecho,'   :: 1=BIN 2=ASC'
c
c        initialize [M]  to zero
!Fangbao         write(ilog,*)'@@ in form mass, ilump: ',ilump
         do i=1,maxmass 
            mss(i)=0.0
         enddo
         rewind(imss)

         do i=1,neq
            wk(i)=0.0d0
         enddo
c
c        form each element matrix, and assemble
         write(*,'(a)')'@@ '
         write(str06,'(1x,i4,a)') (nel/50+1),'>'
         call zzwrt(str06)
         do 50 i=1,nel
            if (mod(i,50) .eq. 0) call zzwrt('.')
            neltype=npijkm(i,1)
            mat=nmprop(i)
*               e0=prop(mat,1)
*               g0=prop(mat,2)
                r0=prop(mat,4)
*               call dmat(e0,g0,dd)
            if (neltype .eq. 3) then
c               truss
            elseif (neltype .eq. 4) then
c               plate
c
            elseif (neltype .eq. 5) then
c               3-D tet solid
c
            elseif (neltype .eq. 9) then
c               3-D solid with hex_8
c
            elseif (neltype .eq. 21) then
c               3-D solid with hex_20
                kmax=neltype-1
                do k=1,kmax
                   node=npijkm(i,1+k)
                   xyz(1,k)=xord(node)
                   xyz(2,k)=yord(node)
                   xyz(3,k)=zord(node)
                enddo
c
                ired_int=3
c               ired_int=2
                call elmmas_HEX20(neltype,ired_int,xyz,r0,emm60,
     &                            ilump,zmmd20)
*             write(iout,*)' m elem: ',i
*             do ii=1,60
*                write(iout,82) (emm60(ii,j), j=1,60)
*             enddo
*             do ii=1,20
*                write(iout,82)  zmmd20(ii)
*             enddo
*             write(iout,*)' m diag: ',i
*             do ii=1,60
*                write(iout,82) emm60(ii,ii)
*             enddo
*             stop
c
c               assemble stiffness
                if (ilump .eq. 1) then
                    do ii=1,3
                       do jj=1,kmax
                          ff(ii,jj)=zmmd20(jj)
                       enddo
                     enddo
                    ielm=i
                    call assembFOR_HEX(wk,ff, idbc,maxnode, 
     &                                   npijkm,maxelem,ielm)
                else
!                    ielm=i
!                    call assembCOL_HEX20(mss,emm60, idbc,maxnode, 
!     &                                   npijkm,maxelem,ielm,nloc)
                  stop 'formmass_Tian wrong!'
                endif
*               do ii=1,maxstiff
*                  write(iout,*) ii,stf(ii)
*               enddo
c
           endif
c
 50      continue 
         write(*,'(a)') '@@ end elems'
c        END loop over all elements
c
*              do ii=1,neq
*                 write(iout,*) ii,geo(nloc(ii))
*              enddo
c
c        STORE mass  matrix on disk in case
                if (ilump .eq. 1) then
                    ibandm=1
                    do ii=1,neq
                       mss(ii)=wk(ii)
                    enddo
                endif
c      modify mass matrix for concentrated masses along diagonal
        rewind(icms)
        do n=1,neq
           read(icms) temp8
           if (ilump.eq.1) then    
               mss(n) = mss(n) + temp8
           elseif (ilump.eq.2) then    
!               iloc=nloc(n)
!               mss(iloc) = mss(iloc) + temp8
           endif
        enddo
!Fangbao        write(*   ,*) '@@ Reloaded  {concmas}  successfully'
!Fangbao        write(ilog,*) '@@ Reloaded  {concmas}  successfully'
c
         rewind(imss)
*        write(*,*) maxstiff,maxnode,neq,' yyy'
         if (ilump .eq. 1) then
             write(imss) neq,ibandm
             do i=1,neq
                write(imss) wk(i)
             enddo
         else
	       stop 'formmass_Tian wrong!'
!             write(imss) neq,iband
!             call storeCOL(mss,maxstiff,iprofv,nloc,imss,neq,iband)
         endif
c
         if (iecho .eq. 2 .OR. iecho .eq. 9) then
             write(iout,'(a)')' MASS : COLumn form '
             do 11 i=1,neq  
 !               iloc=nloc(i)
 !               write(iout,22) (mss(iloc+j-1), j=1,iprofv(i))
 22                   format(1x,6(g13.6))
 11          continue
*            write(iout,*)'loads '
*            do i=1,neq
*               write(iout,23) qms(i),i
 23             format(1x,g13.6,1x,i5 )
*            enddo
         endif
c
         if (iecho .eq. 1 .OR. iecho .eq. 9) then
*            call storeCOL(mss,maxmass,iprofv,nloc,igeo,neq,iband)
c
         elseif (iecho .eq. 3) then
            write(iout,'(a)') ' MASS: band storage'
            write(iout,* ) maxmass,neq,iband    
!            call storeBND_f(mss,maxmass,iprofv,nloc,iout,neq,iband,wk)
c
         elseif (iecho .eq. 4) then
             call zzwrt('INPUT: store_filename --> ')
             read(ikbd,'(a)') str40 
!Fangbao             write(ilog,*) str40
            open(unit=itmp,file=fdn//adjustl(str40) ,form='unformatted')
             rewind itmp
!             call storeBND(mss,maxmass,iprofv,nloc,itmp,neq,iband,wk)
             close(itmp)
         endif
c
 82                   format(1x,60(g13.6))
 83                   format(1x,12(i4,1x))
c
         write(*   ,*) ' @@ neq iband',neq  ,iband
         write(*   ,*) ' @@ FORMMASS:   Formed  [M]  OK'
!Fangbao         write(ilog,*) ' @@ FORMMASS:   Formed  [M]  OK'
      return
      end
!====================================================
c     STIFFness Geometric for HEXahedron 8 nodes
c     based on iso quad from Bathe  pp295-297
      subroutine elmmas_HEX20(neltype,ired_int,xyz,rho0,emm,
     &                        ilump,zmmd20)
         Implicit real*8 (a-h,o-z)
c
         real*8 xyz(3,20),emm(60,60),hh(20),emm20(20,20),zmmd20(20)
         real*8 BB(6,60),bd(9,60)
         real*8 zmmd60(60,60)
c
         real*8 xg(4,4),wgt(4,4),db(9)
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
               emm(i,j)=0.0
               zmmd60(i,j)=0.0
            enddo
         enddo
                  do i=1,kmax
                     do j=1,kmax
                        emm20(i,j)=0.0
                     enddo
                        zmmd20(i)=0.0
                  enddo
c
*          write(*,*)' [x y z] '
*          do i=1,3
*             write(*,82) (xyz(i,j), j=1,8)
*          enddo
*          write(*,*) ' '
c
*        write(iout,*)' density: ',rho0
         ip_n=0
         nint=ired_int
*       nint=2
*       nint=1
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
                  call stdm20(xyz,hh,bb,bd,det,ri,si,ti,nel,ilog)
*                 write(iout,*) ' hh ',ip_n 
*                 do i=1,6
*                    write(iout,82) (hh(j), j=1,kmax)
*                 enddo
c
c                 add contib to mass
                  wt = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)*det
*                 write(iout,86) wgt(lx,nint),wgt(ly,nint),wgt(lz,nint)
*    &                          ,ri,si,ti
                  do 70 j=1,kmax
                     do 60 i=j,kmax   
                        emm20(i,j) = emm20(i,j) + rho0*hh(i)*hh(j)*wt
 60                  continue
 70               continue
 80      continue
c
c        impose symmetry
         do j=1,kmax
            do i=j,kmax
               emm20(j,i)=emm20(i,j)
            enddo
         enddo
                  i3=0
                  do i=1,kmax
                     i3=i3+3
                     j3=0
                     do j=1,kmax   
                        j3=j3+3
                        emm(i3-2,j3-2) = emm20(i,j) 
                        emm(i3-1,j3-1) = emm20(i,j) 
                        emm(i3-0,j3-0) = emm20(i,j) 
                     enddo   
                  enddo   
*             write(iout,*)' emmm 20: '
*             do ii=1,20
*                write(iout,82) (emm20(ii,jj), jj=1,20)
*             enddo
*             write(iout,*)' emm 60: '
*             do ii=1,60
*                write(iout,86) (emm(ii,jj), jj=1,60)
*             enddo
         if (ilump .ne. 1) return
c
         zmtot=0.0
         zmd  =0.0
         do i=1,kmax
               zmd = zmd + emm20(i,i)
            do j=1,kmax
               zmtot = zmtot + emm20(i,j)
            enddo
         enddo
              do ii=1,20
                 zmmd20(ii) = emm20(ii,ii)*zmtot/zmd
              enddo
         return
c
         zm1  =0.0
         zm2  =0.0
         do i=1,8   
            zm1 = zm1 + emm20(i,i)
         enddo
         do i=9,20  
            zm2 = zm2 + emm20(i,i)
         enddo
*        write(iout,*)' total, diag mass: ',zmtot,zmd
*        write(iout,*)' part1, part2    : ',zm1,zm2   
*        write(iout,*)' part1/ part2    : ',zm1/zm2   
c
                  i3=0
                  do i=1,kmax
                     i3=i3+3
                        zmmd60(i3-2,i3-2) = zmmd20(i) 
                        zmmd60(i3-1,i3-1) = zmmd20(i) 
                        zmmd60(i3-0,i3-0) = zmmd20(i) 
                  enddo
         do i=1,kmax*3
            do j=1,kmax*3
               emm(i,j)=zmmd60(i,j)
            enddo
         enddo
              write(iout,*)' emmm 20: '
              do ii=1,20
                 write(iout,82) (emm20(ii,jj), jj=1,20)
              enddo
              write(iout,*)' zmmd 60: '
              do ii=1,60
                 write(iout,82) (zmmd60(ii,jj), jj=1,60)
              enddo
         return
c
*             write(iout,*)' m 20: '
*             do ii=1,20
*                write(iout,82) (emm20(ii,jj), jj=1,20)
*             enddo
*             write(iout,*)' lump m 20: '
*                write(iout,82) (zmmd20(ii), ii=1,20)
         alf=0.2917
         do i=1,8   
            zmmd20(i) = emm20(i,i)*(zmtot/zm1)*(alf/(1+alf))
         enddo
         do i=9,20  
            zmmd20(i) = emm20(i,i)*(zmtot/zm2)*(1.0/(1+alf))
         enddo
*             write(iout,*)' lump 2 20: '
*                write(iout,82) (zmmd20(ii), ii=1,20)
c    
 82        format(1x,40(g12.6,1x))
 84        format(1x,40(i4,1x))
 86        format(1x,60(g12.6,1x))
      return
      end
c
