c
c
c     UDU decomposition of COLumn profiled system
      subroutine uduCOL_D(a,maxstore,neq,imult,iprof,nloc,z0,z1)
         implicit real*8 (a-h,o-z)
         real*8    a(maxstore )
         real*8 temp,sum
         integer iprof(neq), nloc(neq)
c
          imult=0
          zmult=0
          if (a(1  ) .eq. 0.0d0) then
              imult=0
              return
          endif
c
          if (neq .eq. 1) then
              imult=1
              return
          endif
c
          loops=neq/100 + 1
*         call zzwrti(loops)
*         call zzwrt ('>')
          do 10 j=2,neq
*            if (mod(j,100) .eq. 0) call zzwrt('.')
             jm1=j-1
             j2=j-iprof(j)+1
             if (j2.lt.1) then
                 j2=1
             endif
c
c            off-diagonal terms
             if (jm1.eq.1) then
                 is=j
                 io=1
                 iloc = nloc(is) + io - 1
                 sum=a(iloc)
c                sum=a(j,1)
             else
                do 20 i=j2+1,jm1
                   im1=i-1
                   is=j
                   io=j-i+1
                   iloc = nloc(is) + io - 1
                   sum=a(iloc   )
c                  sum=a(i,j-i+1)
c
                           j3=i-iprof(i)+1
                           jcol=j3
                           if (j3 .lt. j2) then
                               jcol=j2
                           endif
c                  do 21 k=j2,im1
                   do 21 k=jcol,im1
                      is=i
                      io=i-k+1
                      iloci = nloc(is) + io - 1
                      is=j
                      io=j-k+1
                      ilocj = nloc(is) + io - 1
                      sum=sum-a(iloci  )*a(ilocj  )
c                     sum=sum-a(k,i-k+1)*a(k,j-k+1)
                      zmult=zmult+1
 21                continue
                   a(iloc   )=sum
c                  a(i,j-i+1)=sum
 20             continue
                is=j
                io=1
                iloc = nloc(is) + io - 1
                sum=a(iloc   )
c               sum=a(j,1)
             endif
c
c            diagonal terms
             do 30 k=j2,jm1
                is=j
                io=j-k+1
                ilocj = nloc(is) + io - 1
                is=k
                io=1
                iloc1 = nloc(is) + io - 1
                temp=a(ilocj  )/a(iloc1)
                sum=sum-temp*a(ilocj  )
                a(ilocj  )=temp
c               temp=a(k,j-k+1)/a(k,1)
c               sum=sum-temp*a(k,j-k+1)
c               a(k,j-k+1)=temp
                zmult=zmult+2
 30          continue
             if (sum.eq.0.0d0) then
                 imult=0
                 return
             endif
             is=j
             io=1
             iloc = nloc(is) + io - 1
             a(iloc   ) = sum
c            a(j,1)=sum
 10        continue
c
c          find largest and smallest diag term
           z0 = 1.0d20
           z1 = -1.0d20
           do i=1,neq
              zz = a(nloc(i))
              if (zz .lt. z0) z0=zz
              if (zz .gt. z1) z1=zz
           enddo
c
           if (zmult .gt. 1.0e6) then
               zmult=zmult*1.0e-6
               imult=zmult
           else
               imult = zmult
           endif
c
      return
      end
c
c
c     BAcK solver of COLumn profiled system
      subroutine bakCOL_D(a,maxstore,b,neq,wk,imult,
     &                                iprof,iprof2,nloc)
         implicit real*8 (a-h,o-z)
         real*8    a(maxstore ),b(neq),wk(neq)
         real*8 sum
         integer iprof(neq), iprof2(neq)
c
         integer  nloc(neq)
c
         imult=0
         do i=1,neq
            wk(i)=0.0
*            write(*,*) i,b(i),' in bakcol'
         enddo
c                                   write(*,*)' in bakCOL '
c
c        forward substitutions
         do 10 i=1,neq
c           j=i-iband+1
            j=i-iprof(i)+1
            if (i.le.iprof(i) ) then
                j=1
            endif
*           jb=i-iband+1
*           jbb=jb
*           if (i.le.iband    ) then
*               jbb=1
*           endif
            sum=b(i)
            km1=i-1
            if (j.gt.km1) then
                wk(i)=sum
            else
                do 11 k=j,km1
                   is=i
                   io=i-k+1
                   iloc = nloc(is) + io - 1
*                  write(*,*) sum,a(iloc   )*wk(k),i,k
                   sum=sum-a(iloc   )*wk(k)
c                  sum=sum-a(k,i-k+1)*wk(k)
                   imult=imult+1
 11             continue
                wk(i)=sum
            endif
 10       continue
c
c         middle terms
          do 30 i=1,neq
             is=i
             io=1
             iloc = nloc(is) + io - 1
             wk(i)=wk(i)/a(iloc )
c            wk(i)=wk(i)/a(i,1)
                   imult=imult+1
*            write(*,*) i,wk(i),' wk '
 30       continue
c
c         backward substitution
          do 50 i1=1,neq
             i=neq-i1+1
             j=i+iprof2(i) -1
             if (j.gt.neq) then
                 j=neq
             endif
*            jb=i+iband-1
*            jbb=jb
*            if (jb.gt.neq) then
*                jbb=neq
*            endif
             sum=wk(i)
             k2=i+1
             if (k2.gt.j) then
                wk(i)=sum
             else
                do 40 k=k2,j
                   is=k
                   io=k-i+1
                   if (io .gt. iprof(is)) goto 40
                   iloc = nloc(is) + io - 1
                   sum=sum-a(iloc   )*wk(k)
c                  sum=sum-a(i,k-i+1)*wk(k)
                   imult=imult+1
 40             continue
                wk(i)=sum
             endif
 50       continue
c
      return
      end
c
c
c     UDU decomposition of COLumn profiled system
      subroutine uduCOL(a,maxstore,neq,iband,imult,iprof,nloc)
         real*8    a(maxstore )
         real*8 temp,sum
         integer iprof(neq), nloc(neq)
c
          imult=0
          if (a(1  ) .eq. 0.0d0) then
              imult=0
              return
          endif
c
          if (neq .eq. 1) then
              imult=1
              return
          endif
c
          loops=neq/100 + 1
          call zzwrti(loops)
          call zzwrt ('>')
          do 10 j=2,neq
             if (mod(j,100) .eq. 0) call zzwrt('.')
             jm1=j-1
             j2=j-iprof(j)+1
             if (j2.lt.1) then
                 j2=1
             endif
c
c            off-diagonal terms
             if (jm1.eq.1) then
                 is=j
                 io=1
                 iloc = nloc(is) + io - 1
                 sum=a(iloc)
c                sum=a(j,1)
             else
                do 20 i=j2+1,jm1
                   im1=i-1
                   is=j
                   io=j-i+1
                   iloc = nloc(is) + io - 1
                   sum=a(iloc   )
c                  sum=a(i,j-i+1)
c
                           j3=i-iprof(i)+1
                           jcol=j3
                           if (j3 .lt. j2) then
                               jcol=j2
                           endif
c                  do 21 k=j2,im1
                   do 21 k=jcol,im1
                      is=i
                      io=i-k+1
                      iloci = nloc(is) + io - 1
                      is=j
                      io=j-k+1
                      ilocj = nloc(is) + io - 1
                      sum=sum-a(iloci  )*a(ilocj  )
c                     sum=sum-a(k,i-k+1)*a(k,j-k+1)
                      imult=imult+1
 21                continue
                   a(iloc   )=sum
c                  a(i,j-i+1)=sum
 20             continue
                is=j
                io=1
                iloc = nloc(is) + io - 1
                sum=a(iloc   )
c               sum=a(j,1)
             endif
c
c            diagonal terms
             do 30 k=j2,jm1
                is=j
                io=j-k+1
                ilocj = nloc(is) + io - 1
                is=k
                io=1
                iloc1 = nloc(is) + io - 1
                temp=a(ilocj  )/a(iloc1)
                sum=sum-temp*a(ilocj  )
                a(ilocj  )=temp
c               temp=a(k,j-k+1)/a(k,1)
c               sum=sum-temp*a(k,j-k+1)
c               a(k,j-k+1)=temp
                imult=imult+2
 30          continue
             if (sum.eq.0.0d0) then
                 imult=0
                 return
             endif
             is=j
             io=1
             iloc = nloc(is) + io - 1
             a(iloc   ) = sum
c            a(j,1)=sum
 10        continue
c
c
      return
      end
c
c
c     BAcK solver of COLumn profiled system
      subroutine bakCOL(a,maxstore,b,neq,iband,wk,imult,
     &                                iprof,iprof2,nloc)
         real*8    a(maxstore ),b(neq),wk(neq)
         real*8 sum
         integer iprof(neq), iprof2(neq)
c
         integer  nloc(neq)
c
          imult=0
c        forward substitutions
         do 10 i=1,neq
c           j=i-iband+1
            j=i-iprof(i)+1
            if (i.le.iprof(i) ) then
                j=1
            endif
            jb=i-iband+1
            jbb=jb
            if (i.le.iband    ) then
                jbb=1
            endif
            sum=b(i)
            km1=i-1
            if (j.gt.km1) then
                wk(i)=sum
            else
                do 11 k=j,km1
                   is=i
                   io=i-k+1
                   iloc = nloc(is) + io - 1
                   sum=sum-a(iloc   )*wk(k)
c                  sum=sum-a(k,i-k+1)*wk(k)
                   imult=imult+1
 11             continue
                wk(i)=sum
            endif
 10       continue
c
c         middle terms
          do 30 i=1,neq
             is=i
             io=1
             iloc = nloc(is) + io - 1
             wk(i)=wk(i)/a(iloc )
c            wk(i)=wk(i)/a(i,1)
                   imult=imult+1
 30       continue
c
c         backward substitution
          do 50 i1=1,neq
             i=neq-i1+1
             j=i+iprof2(i) -1
             if (j.gt.neq) then
                 j=neq
             endif
             jb=i+iband-1
             jbb=jb
             if (jb.gt.neq) then
                 jbb=neq
             endif
             sum=wk(i)
             k2=i+1
             if (k2.gt.j) then
                wk(i)=sum
             else
                do 40 k=k2,j
                   is=k
                   io=k-i+1
                   if (io .gt. iprof(is)) goto 40
                   iloc = nloc(is) + io - 1
                   sum=sum-a(iloc   )*wk(k)
c                  sum=sum-a(i,k-i+1)*wk(k)
                   imult=imult+1
 40             continue
                wk(i)=sum
             endif
 50       continue
c
      return
      end
c
c
c     Complex UDU decomposition of COLumn profiled system
      subroutine cuduCOL(a,maxstore,neq,iband,imult,iprof,nloc)
         implicit real*8 (a-h,o-z)
         complex*16 a(maxstore )
         complex*16 temp,sum
         integer iprof(neq), nloc(neq)
c
          imult=0
          if (a(1  ) .eq. 0.0d0) then
              imult=0
              return
          endif
c
          if (neq .eq. 1) then
              imult=1
              return
          endif
c
          loops=neq/100 + 1
          do 10 j=2,neq
             jm1=j-1
             j2=j-iprof(j)+1
             if (j2.lt.1) then
                 j2=1
             endif
c
c            off-diagonal terms
             if (jm1.eq.1) then
                 is=j
                 io=1
                 iloc = nloc(is) + io - 1
                 sum=a(iloc)
c                sum=a(j,1)
             else
                do 20 i=j2+1,jm1
                   im1=i-1
                   is=j
                   io=j-i+1
                   iloc = nloc(is) + io - 1
                   sum=a(iloc   )
c                  sum=a(i,j-i+1)
c
                           j3=i-iprof(i)+1
                           jcol=j3
                           if (j3 .lt. j2) then
                               jcol=j2
                           endif
c                  do 21 k=j2,im1
                   do 21 k=jcol,im1
                      is=i
                      io=i-k+1
                      iloci = nloc(is) + io - 1
                      is=j
                      io=j-k+1
                      ilocj = nloc(is) + io - 1
                      sum=sum-a(iloci  )*a(ilocj  )
c                     sum=sum-a(k,i-k+1)*a(k,j-k+1)
                      imult=imult+1
 21                continue
                   a(iloc   )=sum
c                  a(i,j-i+1)=sum
 20             continue
                is=j
                io=1
                iloc = nloc(is) + io - 1
                sum=a(iloc   )
c               sum=a(j,1)
             endif
c
c            diagonal terms
             do 30 k=j2,jm1
                is=j
                io=j-k+1
                ilocj = nloc(is) + io - 1
                is=k
                io=1
                iloc1 = nloc(is) + io - 1
                temp=a(ilocj  )/a(iloc1)
                sum=sum-temp*a(ilocj  )
                a(ilocj  )=temp
c               temp=a(k,j-k+1)/a(k,1)
c               sum=sum-temp*a(k,j-k+1)
c               a(k,j-k+1)=temp
                imult=imult+2
 30          continue
             if (sum.eq.0.0d0) then
                 imult=0
                 return
             endif
             is=j
             io=1
             iloc = nloc(is) + io - 1
             a(iloc   ) = sum
c            a(j,1)=sum
 10        continue
c
c
      return
      end
c
c
c     Complex BAcK solver of COLumn profiled system
      subroutine cbakCOL(a,maxstore,b,neq,iband,wk,imult,
     &                                iprof,iprof2,nloc)
         implicit real*8 (a-h,o-z)
         complex*16 a(maxstore ),b(neq),wk(neq)
         complex*16 sum
         integer iprof(neq), iprof2(neq)
c
         integer  nloc(neq)
c
          imult=0
c        forward substitutions
         do 10 i=1,neq
c           j=i-iband+1
            j=i-iprof(i)+1
            if (i.le.iprof(i) ) then
                j=1
            endif
            jb=i-iband+1
            jbb=jb
            if (i.le.iband    ) then
                jbb=1
            endif
            sum=b(i)
            km1=i-1
            if (j.gt.km1) then
                wk(i)=sum
            else
                do 11 k=j,km1
                   is=i
                   io=i-k+1
                   iloc = nloc(is) + io - 1
                   sum=sum-a(iloc   )*wk(k)
c                  sum=sum-a(k,i-k+1)*wk(k)
                   imult=imult+1
 11             continue
                wk(i)=sum
            endif
 10       continue
c
c         middle terms
          do 30 i=1,neq
             is=i
             io=1
             iloc = nloc(is) + io - 1
             wk(i)=wk(i)/a(iloc )
c            wk(i)=wk(i)/a(i,1)
                   imult=imult+1
 30       continue
c
c         backward substitution
          do 50 i1=1,neq
             i=neq-i1+1
             j=i+iprof2(i) -1
             if (j.gt.neq) then
                 j=neq
             endif
             jb=i+iband-1
             jbb=jb
             if (jb.gt.neq) then
                 jbb=neq
             endif
             sum=wk(i)
             k2=i+1
             if (k2.gt.j) then
                wk(i)=sum
             else
                do 40 k=k2,j
                   is=k
                   io=k-i+1
                   if (io .gt. iprof(is)) goto 40
                   iloc = nloc(is) + io - 1
                   sum=sum-a(iloc   )*wk(k)
c                  sum=sum-a(i,k-i+1)*wk(k)
                   imult=imult+1
 40             continue
                wk(i)=sum
             endif
 50       continue
c
      return
      end
c
c
c     ReaD COLumn oriented storage
      subroutine rdCOL( aa, iprof,nloc)
         implicit real*8 (a-h,o-z)
                include 'commons.std'
c
         integer  iprof(neq), nloc(neq)
         real*8   aa(neq*iband)
c
         write(*,*)'@@ reading <<stadyn.stf<<'
         rewind(istf)
         read(istf) j1,j2  
         write(ilog,*)'@@ neq iband ',j1,j2  
         write(iout,*)'@@ reading COL form'  
         write(iout,*)'@@ neq iband ',j1,j2  
         do 8 i=1,neq  
               iloc = nloc(i)
               imax = iprof(i)
               read (istf   ) ( aa(iloc + j-1 ), j=1,imax )
               write(iout,22) ( aa(iloc + j-1 ), j=1,imax )
 22            format(1x,6(g13.6))
 8       continue
c
      return
      end
c
c
c     UDU solution of banded system using UDU decomposition
      subroutine udu(a,neq,iband,imult)
          real*8 a(neq,iband), temp,sum
c
          imult=0
          if (a(1,1) .eq. 0.0  ) then
              imult=0
              return
          endif
          if (neq .eq. 1) then
              imult=1
              return
          endif
c
          loops=neq/100 + 1
          call zzwrti(loops)
          call zzwrt ('>')
          do 10 j=2,neq
             if (mod(j,100) .eq. 0) call zzwrt('.')
             jm1=j-1
             j2=j-iband+1
             if (j2.lt.1) then
                 j2=1
             endif
c
c            off-diagonal terms
             if (jm1.eq.1) then
                 sum=a(j,1)
             else
                do 20 i=j2+1,jm1
                   im1=i-1
                   sum=a(i,j-i+1)
                   do 21 k=j2,im1
                      sum=sum-a(k,i-k+1)*a(k,j-k+1)
                      imult=imult+1
 21                continue
                   a(i,j-i+1)=sum
 20             continue
                sum=a(j,1)
             endif
c
c            diagonal terms
             do 30 k=j2,jm1
                temp=a(k,j-k+1)/a(k,1)
                sum=sum-temp*a(k,j-k+1)
                a(k,j-k+1)=temp
                   imult=imult+2
 30          continue
             if (sum.eq.0.0  ) then
                 imult=0
                 return
             endif
             a(j,1)=sum
 10        continue
c
      return
      end
c
c
      subroutine bak(a,b,neq,iband,wk,imult)
          real*8 a(neq,iband),b(neq),wk(neq), sum
c
          imult=0
c         forward substitutions
          do 10 i=1,neq
             j=i-iband+1
             if (i.le.iband) then
                 j=1
             endif
             sum=b(i)
             km1=i-1
             if (j.gt.km1) then
                wk(i)=sum
             else
                do 11 k=j,km1
                   sum=sum-a(k,i-k+1)*wk(k)
                   imult=imult+1
 11             continue
                wk(i)=sum
             endif
 10       continue
c
c         middle terms
          do 30 i=1,neq
             wk(i)=wk(i)/a(i,1)
                   imult=imult+1
 30       continue
c
c         backward substitution
          do 50 i1=1,neq
             i=neq-i1+1
             j=i+iband-1
             if (j.gt.neq) then
                 j=neq
             endif
             sum=wk(i)
             k2=i+1
             if (k2.gt.j) then
                wk(i)=sum
             else
                do 40 k=k2,j
                   sum=sum-a(i,k-i+1)*wk(k)
                   imult=imult+1
 40             continue
                wk(i)=sum
             endif
 50        continue
c
      return
      end
c
c
c
c     UDU Complex solution of banded system using UDU decomposition
      subroutine uduc(a,neq,iband,imult)
         complex*16 a(neq,iband), temp,sum
c
          imult=0
         temp = a(1,1)
         if (temp   .eq. 0.0) then
             imult=0
             return
         endif
         if (neq .eq. 1) then
             imult=1
             return
         endif
c
         do 10 j=2,neq
            jm1=j-1
            j2=j-iband+1
            if (j2.lt.1) then
                j2=1
            endif
c
c           off-diagonal terms
            if (jm1.eq.1) then
                sum=a(j,1)
            else
                do 20 i=2,jm1
                   im1=i-1
                   if (im1.lt.j2) then
                       goto 20
                   endif
                   sum=a(i,j-i+1)
                   do 21 k=j2,im1
                      sum=sum-a(k,i-k+1)*a(k,j-k+1)
                      imult=imult+1
 21                continue
                   a(i,j-i+1)=sum
 20             continue
                sum=a(j,1)
             endif
c
c            diagonal terms
             do 30 k=j2,jm1
                temp=a(k,j-k+1)/a(k,1)
                sum=sum-temp*a(k,j-k+1)
                a(k,j-k+1)=temp
                   imult=imult+2
 30          continue
             if (sum .eq. 0.0) then
                 imult=0
                 write(*,*)'case 2'
                 return
             endif
             a(j,1)=sum
 10      continue
c
      return
      end
c
c
      subroutine bakc(a,b,neq,iband,wk,imult)
         complex*16 a(neq,iband),b(neq),wk(neq), sum
c
          imult=0
c        forward substitutions
         do 10 i=1,neq
            j=i-iband+1
            if (j .le. 1) then
                j=1
            endif
            sum=b(i)
            km1=i-1
            if (j.gt.km1) then
                wk(i)=sum
            else
                do 11 k=j,km1
                   sum=sum-a(k,i-k+1)*wk(k)
                   imult=imult+1
 11             continue
                wk(i)=sum
            endif
 10      continue
c
c        middle terms
         do 30 i=1,neq
             wk(i)=wk(i)/a(i,1)
             imult=imult+1
 30      continue
c
c        backward substitution
         do 50 i1=1,neq
            i=neq-i1+1
            j=i+iband-1
            if (j.gt.neq) then
                j=neq
            endif
            sum=wk(i)
            k2=i+1
            if (k2.gt.j) then
                wk(i)=sum
            else
                do 40 k=k2,j
                   sum=sum-a(i,k-i+1)*wk(k)
                   imult=imult+1
 40             continue
                wk(i)=sum
            endif
 50      continue
c
      return
      end
c
c
c 
      subroutine ainver(a,n,indx,yn)
         implicit real*8 (a-h,o-z)
         real*8    a(n,n), yn(n,n)
         integer  indx(n)

         do 20  i = 1,n
            do 10  j = 1,n
               yn(i,j) = 0.0
 10         continue
            yn(i,i) = 1.0
 20      continue
c
         np=n
         call ludcmp(a,n,np,indx,d)
c
         do 30  j = 1,n
            call lubksb(a,n,np,indx,yn(1,j))
 30      continue
c
         do 40  j = 1,n
            do 40  i = 1,n
               a(i,j) = yn(i,j)
 40      continue
c
      return
      end
c
c
      subroutine ludcmp(a,n,np,indx,d)
         implicit real*8 (a-h,o-z)
      parameter (nmax=100,tiny=1.0e-20)
      real*8    a(np,np), vv(nmax)
      real*8    sum
      integer   indx(n)
c
      d = 1.0
      do 12  i = 1,n
         aamax = 0.0
         do 11  j = 1,n
            if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
 11      continue
         if (aamax .eq. 0.0) pause 'Singular Matrix'
         vv(i) = 1.0/aamax
 12   continue
      do 19  j = 1,n
         do 14  i = 1,j-1
            sum = a(i,j)
            do 13  k = 1,i-1
               sum = sum - a(i,k)*a(k,j)
 13         continue
            a(i,j) = sum
 14      continue
         aamax = 0.0
         do 16  i = j,n
            sum = a(i,j)
            do 15  k = 1,j-1
               sum = sum - a(i,k)*a(k,j)
 15         continue
            a(i,j) = sum
            dum = vv(i)*abs(sum)
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            endif
 16      continue
         if (j .ne. imax) then
            do 17  k = 1,n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
 17         continue
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j) .eq. 0.0) a(j,j) = tiny
         if (j .ne. n) then 
            dum = 1.0/a(j,j)
            do 18  i = j+1,n
               a(i,j) = a(i,j)*dum
 18         continue
         endif
 19   continue
c
      return
      end 
c
c
      subroutine lubksb(a,n,np,indx,b)
         implicit real*8 (a-h,o-z)
      real*8 a(np,np), b(n)
      real*8    sum
      dimension indx(n)
c
      ii = 0
      do 12  i = 1,n
         ll = indx(i)
         sum = b(ll)
         b(ll) = b(i)
         if (ii .ne. 0) then
            do 11  j = ii , i-1
               sum = sum - a(i,j)*b(j)
 11         continue
         else if (sum .ne. 0.0) then
            ii = i
         endif
         b(i) = sum
 12   continue
      do 14 i = n,1,-1
         sum = b(i)
         if (i .lt. n) then
            do 13  j = i+1,n
               sum = sum - a(i,j)*b(j)
 13         continue
         endif
         b(i) = sum/a(i,i)
 14   continue
      return
      end
c
c
