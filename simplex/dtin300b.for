c
c     DATA INputs from the structural data file
      subroutine datain( 
     &                   xord, yord,zord,
     &                   xload, yload,zload,
     &                    cmass,
     &                   heat,
     &                   idbc, maxdof, maxelem, maxnode,
     &                   npijkm, iglobal,
     &                   prop,nmprop,ippp)
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3)
         integer nmprop(maxelem), nummat(100)
         integer iwk(20)
         real*8  prop(ippp,10)
         real*8  temp(20)
c
         real*8 xord(maxnode), yord(maxnode), zord(maxnode)
         real*8 xload(maxnode),yload(maxnode),zload(maxnode)
         real*8 cmass(maxnode)
         real*8 heat(maxnode,3)
         real*8 temp8(maxnode*3)
c        real*8 temp8(15000)
c
         character*50 title,endin
c
c 
c        CHECK all GROUP sizes first 

         write(*,*)'@@ reading GROUP sizes'
         call grpchk(idat,ilog)
!Fangbao         write(*   ,*)'@@ GROUP sizes OK'
!Fangbao         write(ilog,*)'@@ GROUP sizes OK'
         rewind(idat)
c
c        READ header data
         write(ilog,*)'@@ reading header lines'
         read(idat,'(1a50)',err=9999) title 
         write(*  ,*)'@@ ', title
         read(idat,*,err=9999) iglobal
         read(idat,*,err=9999) iflagC,iflagM,iflagN,iflagL
         read(idat,'(1a50)',err=9999) endin
         call chkend(endin,ichk)
         if (ichk.eq.0) then
!Fangbao             write(ilog,*)'@@ Incorrect HEADER group length'
             write(*   ,*)'@@ Incorrect HEADER group length'
             stop '!! ERROR'
         endif
         ichk=1
!Fangbao         write(ilog,1004)title,iglobal, iflagC,iflagM,iflagN,iflagL
c
c        READ elem data: el# type node#1 node#2 .... 
         write(ilog,*)'@@ reading connectivities lines'
         read(idat,*,err=9999) nel
         if (nel .lt. 1) goto 9000
         if (nel .gt. maxelem) goto 9005
         do 210 i= 1, nel
            do k=1,10
               npijkm(i,k)=1
            enddo
!           read(idat,*,err=9999) j,nelt(j),npi(j),npj(j),npk(j)
            read(idat,*,err=9999) j,neltype,(iwk(k), k=1,neltype-1)
            if (j.lt.1 .OR. j.gt.nel) goto 9025
            npijkm(j,1)=neltype
            do k=1,neltype-1
               npijkm(j,1+k)=iwk(k)
            enddo
 210     continue
         read(idat,'(1a50)') endin
         call chkend(endin,ichk)
         if (ichk.eq.0) then
             write(ilog,*)'@@ Incorrect CONNECTivities data size'
             write(*   ,*)'@@ Incorrect CONNECTivities data size'
             stop '!! ERROR'
         endif
         ichk=1
c        check the data and count the number of triangles
         ntri = 0
         ntruss= 0
         ntetra= 0
         nhex8 = 0
         nhex20= 0
         do 222 j= 1, nel
            neltype=npijkm(j,1)
            if (neltype .eq. 3) ntruss=ntruss+1
            if (neltype .eq. 4) ntri=ntri+1
            if (neltype .eq. 5) ntetra=ntetra+1
            if (neltype .eq. 9) nhex8 =nhex8 +1
            if (neltype .eq. 21) nhex20 =nhex20 +1
 222     continue 
c
c        READ elem material tags 
         write(ilog,*)'@@ reading material lines'
         do i=1,100
            nummat(i)=0
         enddo
         read(idat,*,err=9999) nmat
         if (nmat .lt. 1) goto 9007
         do 212 i= 1, nmat
            read(idat,*,err=9999) nm1,nm2,material 
            if (1.gt.nm1 .OR. nm1.gt.nel) goto 9027
            if (1.gt.nm2 .OR. nm2.gt.nel) goto 9027
            nummat(material)=1
            do 214 j=nm1,nm2
               nmprop(j)=material
 214        continue
 212     continue
         read(idat,'(1a50)',err=9999) endin
         call chkend(endin,ichk)
         if (ichk.eq.0) then
             write(ilog,*)'@@ Incorrect MATERIAL data  length'
             write(*   ,*)'@@ Incorrect MATERIAL data  length'
             stop '!! ERROR'
         endif
         ichk=1
c
         write(ilog ,1005) nel,  ntetra,nhex8,nhex20, nmat
         write(*   ,*) nel,ntetra,nhex8,nhex20,nmat,'  ::elems 123 m'
c
c        READ   nodal coordinate  data: node# x    y   
         write(ilog,*)'@@ reading coordinate lines'
         read(idat,*,err=9999)  nnp
         if (nnp .lt. 2) goto 9010
         if (nnp .gt. maxnode) goto 9015
         do 231 j= 1, nel
            kmax = npijkm(j,1) - 1
            do k=1,kmax
               if (1.gt.npijkm(j,k+1) .OR. npijkm(j,k+1).gt.nnp) then
                   goto 9030
               endif
            enddo
 231     continue
         do 230 i= 1, nnp
            read(idat,*,err=9999) j,z1,z2,z3
            if (1.gt.j .OR. j.gt.nnp) goto 9040
            xord(j)=z1
            yord(j)=z2
            zord(j)=z3
 230     continue
         read(idat,'(1a50)') endin
         call chkend(endin,ichk)
         if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect COORDinate data  length'
            write(*   ,*)'@@ Incorrect COORDinate data  length'
            stop '!! ERROR'
         endif
         ichk=1
c
c        READ  bcs      default is 1    i.e. free      
         do i=1,nnp
            do j=1,3
               idbc(i,j) = 1
            enddo
         enddo
         write(ilog,*)'@@ reading bc lines'
         read(idat,*,err=9999)  nbc
         do 240 i=1,nbc
            read(idat,*,err=9070)node,ixdof,iydof,izdof
                 if (1.gt.node .OR. node.gt.nnp) goto 9045
                 if (0.gt. ixdof .OR.  ixdof.gt.1  ) goto 9051
                 if (0.gt. iydof .OR.  iydof.gt.1  ) goto 9052
                 if (0.gt. izdof .OR.  izdof.gt.1  ) goto 9053
            idbc(node,1) =  ixdof
            idbc(node,2) =  iydof
            idbc(node,3) =  izdof
 240     continue
         read(idat,'(1a50)') endin
         call chkend(endin,ichk)
         if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect BCs data  length'
            write(*   ,*)'@@ Incorrect BCs data  length'
            stop '!! ERROR'
         endif
         ichk=1
c
c        READ  nodal load data: node#  x y z forces,    cmass
        write(ilog,*)'@@ reading load lines'
        do i=1,nnp 
           xload(i)=0.0
           yload(i)=0.0
           zload(i)=0.0
           cmass(i)=0.0
        enddo
c
        read(idat,*,err=9999) nload
        open(unit=itmp,file=fdn//'stadyn.ld3')
        rewind(itmp)
        do i=1,nnp*3
           temp8(i)=0.0
        enddo
        do 232 i=1,nload
            read(idat,*,err=9999)j,z1,z2,z3,z4,z5,z6,z7 
            if (1.gt.j .OR. j.gt.nnp) goto 9047
c           set z7=-3 for 3rd load
            if (z7 .lt. 0.0) iload_num = -z7 + 0.2
            if (z7 .lt. 0.0) then
                ii=(j-1)*3
                temp8(ii+1)=z1
                temp8(ii+2)=z2
                temp8(ii+3)=z3
                  temp8(ii+4)=z4
                  temp8(ii+4)=z5
                  temp8(ii+4)=z6
!      The three lines above were changed by Fangbao Tian to the following lines on 10/04/2011
!                  temp8(ii+4)=z4
!                  temp8(ii+5)=z5
!                  temp8(ii+6)=z6
            else
                xload(j)=z1
                yload(j)=z2
                zload(j)=z3
! Note: In the dataconv, the xyzload will be saved in a file "unit=ilod" (stadyn.lod). When used in the solver, it is loaded again.
! by Fangbao
                  heat(j,1)=z4
                  heat(j,2)=z5
                  heat(j,3)=z6
                cmass(j)=z7
            endif
 232    continue
        do i=1,nnp
           ii=(i-1)*3
           write(itmp,81) (temp8(ii+k), k=1,6)
        enddo
 81     format(1x,30(g13.6,1x))
        close(itmp)
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)
        if (ichk.eq.0) then
           write(ilog,*)'@@ Incorrect LOADS data  length'
           write(*   ,*)'@@ Incorrect LOADS data  length'
           stop '!! ERROR'
        endif
        ichk=1
c
        write(ilog ,1006) nnp, nbc,nload
c
c        READ element material properties 
         write(ilog,*)'@@ reading material props lines'
         read(idat,*,err=9999) nmatp
         if (nmatp .lt. 1) goto 9007
         do 272 i= 1, nmatp
            read(idat,*,err=9028) nmp,(temp(j), j=1,7)
            nummat(nmp)=-1
                do j=1,7
                   prop(nmp,j)=temp(j)
                enddo
 272     continue
         read(idat,'(1a50)') endin
         call chkend(endin,ichk)
         if (ichk.eq.0) then
             write(ilog,*)'@@ Incorrect MATERIAL props data length'
             write(*   ,*)'@@ Incorrect MATERIAL props data length'
             stop '!! ERROR'
         endif
         ichk=1

         do i=1,100
            if (nummat(i) .gt. 0) then
                write(ilog,*)'@@ Incomplete MATERIAL props for tag =',i
                write(*   ,*)'@@ Incomplete MATERIAL props for tag =',i
                stop '!!! ERROR !!!'
            endif
         enddo

c
c
c READ  special materials properties
        write(ilog,*)'@@ reading special material lines'
        a1=1.0
        a2=1.0
        a3=1.0
           dmm=0.0
           dkk=0.0
           dcc=0.0
        Po=0.0
        Px=0.0
        Py=0.0
           gx=0.0
           gy=0.0
           gz=0.0
        nelpress=0
        igravity_on=0
           sy1=50e3
           etan1=100e3
           pxx1=0
           sy2=50e3
           etan1=100e3
           pxx2=0
               ri1=1+0.2
               ri2=2+0.2
               ri3=3+0.2
        read(idat,*,err=9999) nspmat
        do 242 i =1,nspmat
           read(idat,*,err=9999) kode, c1,c2,c3
           if (kode .eq. 2110) then
c              damping
               do j=1,10
                  if ( nummat(j) .eq. -1) then
                        t0=prop(j,3) !thickness
                       ro0=prop(j,4) !solid density
                       goto 244
                  endif
               enddo
 244           continue
               dmm=c1/(ro0)
               dkk=c2
               dcc=c3
               write(ilog,*)'@@ Proportional damping'
               write(ilog,243)'@@ Damp c1 c2 c3: ',c1,c2,c3
               write(ilog,243)'@@      mm kk cc: ',dmm,dkk,dcc
           elseif (kode .eq. 2120) then
c              gravity
               gx=c1
               gy=c2
               gz=c3
               write(ilog,*)'@@ Gravity'
               write(ilog,243)'@@ gx gy gz :     ',c1,c2,c3
               igravity_on=0
               if (abs(gx) .gt. 1e-10 .OR. abs(gy) .gt. 1.0e-10
     &                                .OR. abs(gz) .gt. 1.0e-10) then
                   igravity_on=1
                   write(ilog,*)'@@ Gravity on: ',igravity_on
               endif
           elseif (kode .eq. 2140) then
c              anisotropic
               a1=c1
               a2=c2
               a3=c3
               write(ilog,*)'@@ Anisotropic parameters'
               write(ilog,243)'@@ a1 a2 a3 :     ',c1,c2,c3
               write(ilog,*)'@@ NOT implemented !!! '
*          elseif (kode .eq. 2150) then
*              po=c1
*              px=c2
*              py=c3
*              write(ilog,*)'@@ Pressure parameters'
*              write(ilog,243)'@@ Po Px Py :     ',c1,c2,c3
*              write(ilog,*)'@@ only Po implemented !!! '
*          elseif (kode .eq. 2152) then
*              po=c1
*              nelpress=c2
*              py=c3
*              write(ilog,*)'@@ Pressure parameters'
*              write(ilog,*)'@@ Po nel :     ',po,nelpress
*              write(ilog,*)'@@ only Po implemented !!! '
           elseif (kode .eq. 2341) then
c              plasticity parameter
               sy1 = c1 
               etan1 = c2 
               pxx1 = c3 
               write(ilog,243)'@@ plasticity 1:',sy1,etan1,pxx1
           elseif (kode .eq. 2342) then
c              plasticity parameter
               sy2 = c1 
               etan2 = c2 
               pxx2 = c3 
               write(ilog,243)'@@ plasticity 2:',sy2,etan2,pxx2
           elseif (kode .eq. 2410) then
c              reduced integration
               ri1=c1+0.2
               ri2=c2+0.2
               ri3=c3+0.2
               write(ilog,*)'@@ reduced integration'
               write(ilog,243)'@@ T H8 H20 :     ',c1,c2,c3
           endif
 242    continue
 243    format(a,1x,3(g13.6,1x))
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)
        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect SPecials data length'
            write(*   ,*)'@@ Incorrect SPecials data length'
            stop '!! ERROR'
        endif
        ichk=1
c
         goto 999
c
c                     FORMATS
c
 1004   format (' ',/,'@@ HEADER GROUP                     ',/
     >         ,'@@',3x,'Title        =',3x,a50/
     >         ,'@@',3x,'Problem type =',i4/
     >         ,'@@',3x,'element flag ='i4/
     >         ,'@@',3x,'element flag ='i4/
     >         ,'@@',3x,'node    flag ='i4/
     >         ,'@@',3x,'node    flag ='i4/)
 1005  format(' ',/,'@@ ELEMENT GROUP                 ',/
     &       ,'@@',3x,'number of elements           ='i4/
     &       ,'@@',3x,'number of tetra  elements    ='i4/
     &       ,'@@',3x,'number of hex_8  elements    ='i4/
     &       ,'@@',3x,'number of hex_20 elements    ='i4/
     &       ,'@@',3x,'number of materials          ='i4/)
 1006  format(' ',/,'@@ NODE GROUP                    ',/
     >       ,'@@',3x,'number of nodal points      ='i4/
     >       ,'@@',3x,'number of boundary nodes    ='i4/
     >       ,'@@',3x,'number of loaded nodes      ='i4/)
c
c           ERROR MESSAGES
c
 9000 write(ierr,*)'Number of elements are less than 1 !!!'
      write(ilog,*)'Number of elements are less than 1 !!!'
      write(ierr,*)'nel =',nel
      goto 998
 9005 write(ierr,*)'Number of elements greater than maximum!!!'
      write(ilog,*)'Number of elements greater than maximum!!!'
      write(ierr,*)'nel =',nel
      goto 998
 9007 write(ierr,*)'Number of element materials are less than 1 !!!'
      write(ilog,*)'Number of element materials are less than 1 !!!'
      write(ierr,*)'nmat =',nmat
      goto 998
 9010 write(ierr,*)'Number of nodes less than 2 !!!'
      write(ierr,*)'nnp =',nnp
      goto 998
 9015 write(ierr,*)'Number of nodes greater than maximum !!!'
      write(ilog,*)'Number of nodes greater than maximum !!!'
      write(ierr,*)'nnp =',nnp, maxnode
      goto 998
 9020 write(ierr,*)'Number boundary nodes greater than number nodes !!!'
      write(ilog,*)'Number boundary nodes greater than number nodes !!!'
      write(ierr,*) nbc,' > ',nnp
      goto 998
 9022 write(ierr,*)'Number load nodes greater than number nodes !!!'
      write(ilog,*)'Number load nodes greater than number nodes !!!'
      write(ierr,*) nload,' > ',nnp
      goto 998
 9025 write(ierr,*)'Element number out of range !!!'
      write(ierr,*)'data at',i
      goto 998
 9027 write(ierr,*)'Material element number out of range !!!'
      write(ierr,*)'data at',i
      goto 998
 9028 write(ierr,*)'Material data ERROR  !!!'
      write(ierr,*)'data at',i
      write(*,*) nmp
      goto 998
 9030 write(ierr,*)'Node number on element out of range !!!'
      write(ierr,*)'data at',j
      goto 998
 9035 write(ierr,*)'Element type not set to 1, 2, or 3 !!!'
      write(ierr,*)'data at',j
      goto 998
 9040 write(ierr,*)'Node number out of range !!!'
      write(ierr,*)   j,' > ',nnp
      goto 998
 9045 write(ierr,*)'Boundary node number out of range !!!'
      write(ierr,*)node,' > ',nnp
      goto 998
 9047 write(ierr,*)'Load node number out of range !!!'
      write(ierr,*)   j,' > ',nnp
      goto 998
c
 9051 write(*   ,*)'X degree of freedom not set to 0 or 1 !!!'
      write(*   ,*)' dof =',ixdof
      goto 998
 9052 write(*   ,*)'Y degree of freedom not set to 0 or 1 !!!'
      write(*   ,*)' dof =',iydof
      goto 998
 9053 write(*   ,*)'Z degree of freedom not set to 0 or 1 !!!'
      write(*   ,*)' dof =',izdof
      goto 998
*9054 write(*   ,*)'Xrot degree of freedom not set to 0 or 1 !!!'
*     write(*   ,*)' rot =',ixrot
*     goto 998
*9055 write(*   ,*)'Yrot degree of freedom not set to 0 or 1 !!!'
*     write(*   ,*)' rot =',iyrot
*     goto 998
*9056 write(*   ,*)'Zrot degree of freedom not set to 0 or 1 !!!'
*     write(*   ,*)' rot =',izrot
*     goto 998
c
 9060 write(ierr,*)'Theta degree of freedom not set to 0 or 1 !!!'
      write(ierr,*)' dof =',izdof
      goto 998
 9065 write(ierr,*)'iglobal should be set > 40  '
      goto 998
 9070 write(*   ,*)'ERROR not enough data !!! '
      goto 998
 9999 continue
      write(ilog,*)'@@ ERROR reading structure DataFile  !!!'
      write(*   ,*)'@@ ERROR reading structure DataFile  !!!'
      goto 998
c
 998     continue
         write(ilog,*)'@@ !!!! ERROR in data INPUT file  !!!!!'
         stop
c
 999     continue
         write(*,*) '@@ DATAin:  Input data read OK'
      return
      end
c
c     CHecK END of each group
      subroutine chkend(endin,ichk)
         character*50 endin
         character*3  end1,end2
         end1='end'
         end2='END'
         i1=index(endin,end1)
         i2=index(endin,end2)
         if (i1.eq.0 .and. i2.eq.0) then
             ichk = 0
         else
             ichk = 1
         endif
      return
      end
c
c     GRouP CHecKing
      subroutine grpchk(idat,ilog)
        character*50 endin
        character*1 ch
	 integer::kk(22)
c
        rewind idat
c       READ header  data
        write(ilog,*)'@@ reading HEADER group'
        read(idat,'(1a1)', end=9999,err=9999) ch    
        read(idat,'(1a1)', end=9999,err=9999) ch    
        read(idat,'(1a1)', end=9999,err=9999) ch    
        read(idat,'(1a50)',end=9999,err=9999) endin

        call chkend(endin,ichk)
        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect HEADER group length'
            write(*   ,*)'@@ Incorrect HEADER group length'
            stop'!! ERROR'
        endif
        ichk=1
c
c       READ elem GROUP sizes  
        write(ilog,*)'@@ reading connectivities GROUP'
        read(idat,*,err=9999) nel
        if (nel .lt. 1) goto 9000
        do 210 i= 1, nel
!           read(idat,'(1a1)') ch  
            read(idat,*)(kk(jj),jj=1,22)  
 210    continue
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)

        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect Connect. GROUP size'
            write(*   ,*)'@@ Incorrect Connect. GROUP size'
            stop'!! ERROR'
        endif
        ichk=1
c
c       READ elem material GROUP size
        write(ilog,*)'@@ reading material GROUP'
        read(idat,*,err=9999) nmat
        if (nmat .lt. 1) goto 9007
        do 212 i= 1, nmat
           read(idat,'(1a1)') ch    
 212    continue
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)
        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect Material GROUP size'
            write(*   ,*)'@@ Incorrect Material GROUP size'
            stop'!! ERROR'
        endif
        ichk=1
c
c       READ  nodal coordinate  GROUP size
        write(ilog,*)'@@ reading coordinate GROUP'
        read(idat,*,err=9999)  nnp
!	write(*,*)nnp
        if (nnp .lt. 2) goto 9010
        do 230 i= 1, nnp
!	      read(1,*)j
!	write(*,*)j
           read(idat,'(1a1)') ch   
 230    continue
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)
!		write(*,*)nnp, endin,ichk
!	stop
        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect Coordinate GROUP size'
            write(*   ,*)'@@ Incorrect Coordinate GROUP size'
            stop'!! ERROR'
        endif
        ichk=1
c
c       READ bcs GROUP size 
        write(ilog,*)'@@ reading bc GROUP'
        read(idat,*,err=9999)  nbc
        do 240 i=1,nbc
           read(idat,'(1a1)') ch    
 240    continue
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)
        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect BCs GROUP size'
            write(*   ,*)'@@ Incorrect BCs GROUP size'
            stop'!! ERROR'
        endif
        ichk=1
c
c       READ  nodal load GROUP size
        write(ilog,*)'@@ reading load GROUP'
        read(idat,*,err=9999) nload
        do 232 i= 1, nload
           read(idat,'(1a1)') ch    
 232    continue
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)
        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect Loads GROUP size'
            write(*   ,*)'@@ Incorrect Loads GROUP size'
            stop'!! ERROR'
        endif
        ichk=1
c
c       READ  material properties GROUP size 
        write(ilog,*)'@@ reading materials properties GROUP'
        read(idat,*,err=9999) nmatp
        do i= 1, nmatp
           read(idat,'(1a1)') ch    
        enddo
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)
        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect materials props GROUP size'
            write(*   ,*)'@@ Incorrect materials props GROUP size'
            stop'!! ERROR'
        endif
        ichk=1
c
c       READ  SPecials GROUP size 
        write(ilog,*)'@@ reading SPecial materials GROUP'
        read(idat,*,err=9999) nspmat
        do 242 i= 1, nspmat
           read(idat,'(1a1)') ch    
 242    continue
        read(idat,'(1a50)') endin
        call chkend(endin,ichk)
        if (ichk.eq.0) then
            write(ilog,*)'@@ Incorrect SPecials GROUP size'
            write(*   ,*)'@@ Incorrect SPecials GROUP size'
            stop'!! ERROR'
        endif
        ichk=1
c
        goto 999
c
c           ERROR MESSAGES
c
 9000 write(*   ,*)'Number of elements are less than 1 !!!'
      write(*   ,*)'nel =',nel
      stop
 9007 write(*   ,*)'Number of element materials are less than 1 !!!'
      write(*   ,*)'nmat =',nmat
      stop
 9010 write(*   ,*)'Number of nodes less than 2 !!!'
      write(*   ,*)'nnp =',nnp
      stop
 9999 continue
      write(ilog,*)'@@ ERROR reading structure DataFile  !!!'
      write(*   ,*)'@@ ERROR reading structure DataFile  !!!'
      stop
c
 999  continue
      write(*,*) '@@ GROUP sizes OK'
      return
      end
c
c
c     DATA CoNVerter
      subroutine datacnv( 
     &                   xord, yord, zord,
     &                   xload, yload,zload,
     &                   cmass,
     &                   heat,wk,
     &                   idbc, maxdof, maxelem, maxnode,
     &                   npijkm, iglobal,
     &                   iprofv,iprofh,
     &                   prop,nmprop,ippp)
c
         implicit real*8 (a-h,o-z)
             include 'commons.std'
c
         integer npijkm(maxelem,21),idbc(maxnode,3)
         integer iprofv(maxdof),iprofh(maxdof)
         integer nmprop(maxelem)
         real*8  prop(ippp,10)
c
         real*8  xord(maxnode), yord(maxnode), zord(maxnode)
         real*8 xload(maxnode),yload(maxnode),zload(maxnode)
         real*8 cmass(maxnode),heat(maxnode,3)
         real*8 wk(maxnode*3)
         real*8 templ(6),tempm(6)
c
         character*50 title
         character*12 eltype
c
c
c READ header data again to get flags
         rewind(idat)
         read(idat,'(1a50)') title 
         write(*    ,*)'@@ ', title
         read(idat,*,err=9999) iglobal
         read(idat,*) iflagC,iflagM,iflagN,iflagL
c
c        CHECKS
         if (iflagC .gt. 0) then
c            echo connectivities
             write(ilog,10021)
             do 260 i = 1, nel
                neltype = npijkm(i,1)
                if (neltype .eq. 3) then
                    eltype = 'truss'
                else if (neltype .eq. 4) then
                    eltype = 'triangle'
                else if (neltype .eq. 5) then
                    eltype = 'tetrahedron'
                else if (neltype .eq. 9) then
                    eltype = 'hex_8   '
                else if (neltype .eq. 21) then
                    eltype = 'hex_20'
                else
                    eltype = ' ????'
                    go to 9035
                endif
                kmax=npijkm(i,1)
                write(ilog ,10033) i,eltype,(npijkm(i,k), k=1,kmax)
 260         continue
         endif
c
c
         if (iflagM .gt. 0) then
c           echo materials
            write(ilog,'(a)')'@@  elem     mat #'
            do i=1,nel
               write(ilog ,'(2h@@,1x,3(i5,1x))') i, nmprop(i)
            enddo 
            write(ilog,1002)
            do 261 i=1,ippp 
               if (prop(i,1) .gt. 0.0) then
                   write(ilog ,1003) i, (prop(i,j), j=1,7)
               endif
 261        continue
         endif
c
c        IMPOSE conditions for special GLOBAL shapesa 
         call spclbc( maxdof,maxnode,idbc,  iglobal)
*           write(ilog,*)'@@ NUMBERs 3'
*           do i=1,nnp 
*              write(ilog ,* )i,(idbc(i,j),j=1,3) 
*           enddo   
c
c        Number the equations in JBC from 1  to  ORDER 
c        Start assigning equation numbers for non-zero dof`s
c        from 1 up; only non-zero given a number
         neq=0
         do i=1,nnp
            do j=1,3
                if (idbc(i,j) .gt. 0) then
                    neq = neq + 1
                    idbc(i,j) = neq
                else
                    idbc(i,j) = 0
                endif
            enddo
         enddo
c
         if (iflagN.gt.0) then
c           echo coordinates
            write(ilog,1028 )
            do i=1,nnp 
               write(ilog ,1027 )i, xord(i), yord(i), zord(i)
            enddo
c           echo equation numbering
            write(ilog,1008 )
            do 280 i=1,nnp 
               write(ilog ,1007 )i,(idbc(i,j),j=1,3) 
 280        continue
         endif
c
c       IDBC() contains the eqn # sequence 1 2 3 4 .....
c       and wk will contain the corresponding applied values
c       and wk will contain the corresponding concent. masses
c
         do 250 i= 1, nnp
            do k=1,3
               ieqnum=idbc(i,k)
               if (ieqnum.ne.0 )   wk(ieqnum) = cmass(i)
            enddo
 250     continue
         write(*,*)'@@ writing mass loads'
         rewind(icms)
         do i=1,neq
            write(icms) wk(i)
         enddo
c
c        compute total mass
         call compmass(npijkm,nmprop,prop,ippp, 
     &                    xord,yord,zord,cmass,maxnode,maxelem,zmass)
         write(ilog,*) '@@ total mass ',zmass
         write(*   ,*) '@@ total mass ',zmass
c
         do 252 i=1,neq  
            wk(i)=0.0
 252     continue
         do 254 i= 1, nnp
               ieqnum=idbc(i,1)
               if (ieqnum .ne. 0)   wk(ieqnum) = xload(i)
               ieqnum=idbc(i,2)
               if (ieqnum .ne. 0)   wk(ieqnum) = yload(i)
               ieqnum=idbc(i,3)
               if (ieqnum .ne. 0)   wk(ieqnum) = zload(i)
 254     continue
         write(*,*)'@@ writing loads'
         rewind(ilod)
         do i=1,neq
            write(ilod) wk(i)
         enddo
c
c        compute resultants
         xres=0.0
         yres=0.0
         zres=0.0
         do i=1,nnp
            xres=xres+xload(i) 
            yres=yres+yload(i)  
            zres=zres+zload(i)
         enddo
         write(ilog,*) '@@RESULTants:'
         write(ilog,'(a,3(g12.6,1x))') '@@force: ', xres,yres,zres
c
c
         if (iflagL.gt.0) then
c           echo loads and conc masses
            write(ilog ,1018) 
            do 284 i=1,nnp 
               do j=1,3
                  ieqnum = idbc(i,j)
                  if (ieqnum .eq. 0) then
                     templ(j)=0.0
                  else
                     templ(j)=wk(ieqnum)
                  endif
               enddo
               write(ilog ,1017) i,(templ(j),j=1,3) 
     &                            ,(heat(i,j),j=1,3),cmass(i)
 284        continue
         endif
c
!remove the commands when use other method rather than _Tian
!TIAN         call maxband(npijkm,idbc,iband,maxdof,maxelem,nel,
!TIAN     &               iprofv,maxnode)
         write(*   ,*)'@@ from datain,  half band =',iband,neq
         write(ilog,*)'@@ from datain,  half band =',iband,neq
c        compute storage
         iloc=0
         do i=1,neq
            iloc=iloc+iprofv(i)
         enddo
!TIAN         maxstiff=iloc
!TIAN         percent=(maxstiff*100.0)/(neq*iband)
!TIAN         ipercent=percent 
!TIAN         write(ilog,1022) neq,iband,neq*iband,iloc,ipercent
c
c        STORE results for future reference
*        itmp=7
*        open( itmp, file='stadyn.jbc')
*        rewind itmp
*        write(itmp,33) neq,iband, nnp,nel 
c
c         calculate other profile
          do i=1,neq
             ibandh=1
             iend=i+iband-1
             if (iend .gt. neq) iend=neq
             do j=i+1,iend
                ibandv=iprofv(j)
                ji1=j-i+1
                if (ibandv .ge. ji1) then
                    ibandh = ji1
                endif
             enddo
             iprofh(i)=ibandh
*            write(itmp,33) i,iprofv(i),ibandh
          enddo
*        write(itmp,33) (jbc(n), n=1,nnp*6)
*        write(itmp,33) (npi(n), n=1,nel)
*        write(itmp,33) (npj(n), n=1,nel)
*        write(itmp,33) (npk(n), n=1,nel)
 33      format(12(i5,1x))
 34      format(1x,10(i5,1x))
*        close (itmp)
c
c
        goto 999
c
c
c                     FORMATS
10021   format(' ',/,'@@ Connectivities',/,
     &  '@@  el.   type          #    i     j     k     m  ... ')
10031   format('@@',4i5,4x,a5)
10032   format('@@',5i5,4x,a5)
10033   format('@@',i5,4x,a6,4x,21(i4,1x))
 1002   format(' ',/,'@@ Materials',/,
     &  '@@    #     E         G      ?        density     nplane ',
     &  '     C       ? ')
 1003   format('@@',i5,1x,9(g9.3,1x)) 
*1004   format (' ',/,'@@ HEADER GROUP                     ',/
*    >         ,'@@',3x,'Title        =',3x,a50/
*    >         ,'@@',3x,'Problem type =',i4/
*    >         ,'@@',3x,'element flag ='i4/
*    >         ,'@@',3x,'element flag ='i4/
*    >         ,'@@',3x,'node    flag        ='i4/
*    >         ,'@@',3x,'node    flag        ='i4/)
*1005  format(' ',/,'@@ ELEMENT GROUP                 ',/
*    &       ,'@@',3x,'number of elements          ='i5/
c    &       ,'@@',3x,'number of truss elements    ='i5/
*    &       ,'@@',3x,'number of triangle elements ='i5/
*    &       ,'@@',3x,'number of frame elements    ='i5/
*    &       ,'@@',3x,'number of materials         ='i5/)
*1006  format(' ',/,'@@ NODE GROUP                    ',/
*    >       ,'@@',3x,'number of nodal points      ='i5/
*    >       ,'@@',3x,'number of boundary nodes    ='i5/
*    >       ,'@@',3x,'number of loaded nodes      ='i5/)
 1007   format('@@',i6,6i8)
 1008   format( /,'@@  numbering of equations'/, 
     &   '@@   node     u       v       w  ')  
 1017   format('@@',i5,6(f9.3,1x),1f10.4)
 1018   format( /,'@@  applied nodal loads and conc masses'/, 
     &   '@@   node     Px       Py        Pz  ' ,
     &           '      P_q      T_o       ?        mass' )
 1027   format('@@',i5,1x,3(f10.4,1x))
 1028   format( /,'@@  node coordinates '/, 
     &            '@@   node      x        y         z  ')
 1022  format(' ',/,'@@ SYSTEM GROUP                  ',/
     &       ,'@@',3x,'system size (order)         ='i5/
     &       ,'@@',3x,'half bandwidth              ='i5/
     &       ,'@@',3x,'NB Storage                  ='i10/
     &       ,'@@',3x,'Profile Storage             ='i10/
     &       ,'@@',3x,'% of NB                     ='i10/)
c
c
c           ERROR MESSAGES
c
 9035 write(*,*)'Element type not set to 1, 2, 3, 4, 5 or 6 !!!'
      write(*,*)'data at',i
      stop
 9999 continue
      write(ilog,*)'@@ ERROR reading structure DataFile  !!!'
      write(*   ,*)'@@ ERROR reading structure DataFile  !!!'
      stop
c
 999  continue
      write(*,*)'@@ DATAck:  Input data checks OK'
      return
      end
c
c
c     MAX BaNDwidth calculation
      subroutine maxband(npijkm,idbc,iband,maxdof,maxelem,
     &                  nel,iprofv,maxnode)
         integer npijkm(maxelem,21)
         integer idbc(maxnode,3)
         integer idof_node(60),all_eqn(60)
         integer iprofv(maxdof)
c
         do 20 n=1,maxdof
            iprofv(n)=0
 20      continue
         ihfbnd=0
         do 2108 n=1,nel
            neltype = npijkm(n,1)
            kmax  = neltype-1
              do kk=1,kmax
                 nodek = npijkm(n,1+kk)
                 all_eqn( (kk-1)*3 + 1) = idbc(nodek,1)
                 all_eqn( (kk-1)*3 + 2) = idbc(nodek,2)
                 all_eqn( (kk-1)*3 + 3) = idbc(nodek,3)
              enddo
*           do i=1,3
*              do k=1,kmax
*                 idof_node(i+(k-1)*3)=npijkm(n,1+k)
*              enddo
*           enddo
            nmax=kmax*3
            do 2104 i=1,kmax*3
*              inn = (i-1)/3 + 1
*              ii=i-(inn-1)*3
*                 ieqn1=idbc(idof_node(i),ii)
                  ieqn1=all_eqn(i)
                  if (ieqn1 .gt. 0) then
                      do 2106 j=i,kmax*3
*                        jnn = (j-1)/3 + 1
*                        jj=j-(jnn-1)*3
*                           ieqn2=idbc(idof_node(j),jj)
                         ieqn2=all_eqn(j)
*           write(*,85) n,i,ii,ieqn1,j,jj,ieqn2
                            if (ieqn2 .gt. 0) then
                                ihfbnd = max0(ihfbnd,iabs(ieqn1-ieqn2))
                                jband=abs(ieqn1-ieqn2)+1
                                ieq=max(ieqn1,ieqn2)
                                if (jband .gt. iprofv(ieq)) then
                                    iprofv(ieq)=jband
                                endif
*                               ieq2=min(ieqn1,ieqn2)
*                                if (jband .gt. iprofh(ieq2)) then
*                                    iprofh(ieq2)=jband
*                                endif
                            endif
 2106                 continue
                  endif
 2104       continue
 2108    continue
         iband=ihfbnd+1
 85      format(1x,20(i6,1x))
c
      return
      end
c
c
c     impose SPeCiaL global BCs
      subroutine spclbc( maxdof,maxnode,idbc, iglobal)
c
         implicit real*8 (a-h,o-z)
               include 'commons.std'
c
*        integer jbc(maxdof),idbc(maxnode,3)
         integer idbc(maxnode,3)
         idofx=1
         idofy=1
         idofz=1
         if (iglobal .ge. 100) then
             idofx = 1
             idofy = (iglobal - 100*idofx)/10
             idofz =  iglobal - 100*idofx - 10*idofy
         endif
         if (iglobal .lt. 100 .AND. iglobal .ge. 10) then
             idofx = 0
             idofy = (iglobal - 100*idofx)/10
             idofz =  iglobal - 100*idofx - 10*idofy
         endif
         if (iglobal .lt. 10) then
             idofx = 0
             idofy = 0
             idofz =  iglobal - 100*idofx - 10*idofy
         endif
         write(ilog,*)'@@ DoF: ',idofx,idofy,idofz
c
         do i=1,nnp
            if (idofx .eq. 0) idbc(i,1) = idofx
            if (idofy .eq. 0) idbc(i,2) = idofy
            if (idofz .eq. 0) idbc(i,3) = idofz
         enddo
      return
      end
c
c
c     COMPute total MASS
      subroutine compmass(npijkm,nmprop,prop,ippp, 
     &                    xord,yord,zord,cmass,maxnode,maxelem,zmass)
c
         implicit real*8 (a-h,o-z)
         include 'commons.std'
c
         integer nmprop(maxelem)
         integer npijkm(maxelem,21)
         real*8  prop(ippp,10)
         real*8 xord(maxnode), yord(maxnode), zord(maxnode)
         real*8 cmass(maxnode)
c
         zmass=0
         do i=1,nnp
            zmass=zmass+cmass(i)
         enddo
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
              mat=nmprop(i)
              rh0=prop(mat,4)
              if (neltype .eq. 3) then
                  a0=prop(mat,3)
                  zlen=sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                  vol=zlen*a0
c
              elseif (neltype .eq. 4) then
                  tt0=prop(mat,3)
c
c                 determine vector area
                  axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
                  ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
                  azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
                  area=sqrt( axy*axy + ayz*ayz + azx*azx)
                  vol=area*tt0
c
              elseif (neltype .eq. 5) then
                  x14 = x1-x4
                  x24 = x2-x4
                  x34 = x3-x4
                    y14 = y1-y4
                    y24 = y2-y4
                    y34 = y3-y4
                      z14 = z1-z4
                      z24 = z2-z4
                      z34 = z3-z4
                  detj = x14*(y24*z34-y34*z24)
     &                 + y14*(z24*x34-z34*x24)
     &                 + z14*(x24*y34-x34*y24)
                  vol=abs(detj)/6.0
c
              elseif (neltype .eq. 9 .OR. neltype .eq. 21) then
                   i1=npijkm(i,6)
                   j1=npijkm(i,7)
                   k1=npijkm(i,8)
                   m1=npijkm(i,9)
                         x21=xord(i1)
                         x22=xord(j1)
                         x23=xord(k1)
                         x24=xord(m1)
                    y21=yord(i1)
                    y22=yord(j1)
                    y23=yord(k1)
                    y24=yord(m1)
                         z21=zord(i1)
                         z22=zord(j1)
                         z23=zord(k1)
                         z24=zord(m1)
                  x14 = x1-x4
                  x24 = x2-x4
                  x34 = x3-x4
                    y14 = y1-y4
                    y24 = y2-y4
                    y34 = y3-y4
                      z14 = z1-z4
                      z24 = z2-z4
                      z34 = z3-z4
                  detj = x14*(y24*z34-y34*z24)
     &                 + y14*(z24*x34-z34*x24)
     &                 + z14*(x24*y34-x34*y24)
                  vol=abs(detj)/6.0
              endif
              dmass=vol*rh0
              zmass=zmass+dmass
         enddo
c
      return
      end
c
