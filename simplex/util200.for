cc \begin{verbatim}
c
c     interface to integer*2 function system [c]
c    &                        (string[reference])
c         character*1 string
c     end
c
cMS    interface to integer*2 function system [c]
cMS  >                        (string[reference])
cMS       character*1 string
cMS    end
c
c
c     CHeCKs and zeros SPaCE for real*8 array
      subroutine ChckSpce(ilog,isize,array,subname,iendp1)
         implicit real*8 (a-h,o-z)
         real*8 array(isize)
         character*(*) subname
         iend=iendp1-1
*        write(*,*) ilog,isize,'  ',subname,iendp1
         write(ilog,*)'@@ ',subname,': memory used ',iend,' of',isize
         write(*   ,*)'@@ ',subname,': memory used ',iend,' of',isize
c
         if (isize .le. iend) then
             write(ilog,*)'@@  !!!ERROR: not enough reserved memory!!!'
             write(ilog,*)'@@  Try change parameters in <<------.cfg>>'
             write(*   ,*)'@@  !!!ERROR: not enough reserved memory!!!'
             write(*   ,*)'@@  Try change parameters in <<------.cfg>>'
             stop 
         else 
            do i=1,iend
               array(i)=0.0
            enddo
            return
         endif
c
      return
      end
c
c     Integer CHeCK SPaCE with zeroing
      subroutine IChckSpce(ilog,isize,iarray,subname,iendp1)
         integer iarray(isize)
         character*(*) subname
         iend=iendp1-1
         write(ilog,*)'@@ ',subname,': memory used ',iend,' of',isize
         write(*   ,*)'@@ ',subname,': memory used ',iend,' of',isize
c
         if (isize .le. iend) then
             write(ilog,*)'@@  !!!ERROR: not enough reserved memory!!!'
             write(ilog,*)'@@  Try change parameters in <<------.cfg>>'
             write(*   ,*)'@@  !!!ERROR: not enough reserved memory!!!'
             write(*   ,*)'@@  Try change parameters in <<------.cfg>>'
             stop 
         else 
            do i=1,iend
               iarray(i)=0
            enddo
            return
         endif
c
      return
      end
c
c     ChecKs SPACE but does not zero
      subroutine ckspace(ilog,isize,subname,iend)
         implicit real*8 (a-h,o-z)
         character*(*) subname
         write(ilog,*)'@@ ',subname,': memory used ',iend,' of',isize
         write(*   ,*)'@@ ',subname,': memory used ',iend,' of',isize
c
         if (isize .le. iend) then
             write(ilog,*)'@@  !!!ERROR: not enough reserved memory!!!'
             write(ilog,*)'@@  Try change parameters in <<------.cfg>>'
             write(*   ,*)'@@  !!!ERROR: not enough reserved memory!!!'
             write(*   ,*)'@@  Try change parameters in <<------.cfg>>'
             stop 
         else 
            return
         endif
c
      return
      end
c
c
      subroutine syscom  
         integer istatus, fsystem
         integer*2 system
cMS      pause' Type: system command'
cUNiX    character*40 commnd 
cUNiX       write(*,*)' Type: system command'
cUNiX       read(ikbd,'(1a40)') commnd
cUNiX       call system(commnd)
         character*40 commnd 
c           write(*,*)'TYPE: system command >> '
c           read(ikbd,'(1a40)') commnd
cL          call system(commnd)
cWAT        istatus = fsystem(commnd)
c           istatus =  system(commnd)
           write(*,'(//)')
           write(*,*)' press <<RETURN>> to continue'
           read(ikbd,*)
      return
      end
c
c
      subroutine zzwrt(string)
         character*(*) string
         write(*,'(a\)') string 
cUNiX    write(*,'(a$)') string 
      return
      end
c
      subroutine zzwrti(number)
         write(*,'(i5\)') number 
cUNiX    write(*,'(i5$)') number 
      return
      end
c
      subroutine zzwrti2(number,string)
         character*(*) string
         write(*,'(i5,a\)') number,string 
cUNiX    write(*,'(i5,a$)') number,string 
      return
      end
c
c
!/////////Fangbao Tian commented this block
!!      subroutine zzcount(icount,inc,icarriage)
!!        use ifport
!!         use ifcore
!!         logical result
!!            if (mod(icount,inc) .eq. 0) then
!!                write(*,'(i5\)') icount
cUNiX           write(*,'(i5$)') icount
cUNiX           call flush(6)
!!                icarriage=icarriage+1
!!                if (icarriage .eq. 15) then
!!                    write(*,'(/a\)') ' @@ '
cUNiX               write(*,'(/a$)') ' @@ '
!!                    icarriage=0
!!                endif
!!                call flush(6)
!!         result = commitqq(i_unit)
!!            endif
!!      return
!!      end
c
      subroutine zzcount_2(icount,inc,icarriage)
         character*3 ch4
         save k1e5,jcount
c
         if (icount .eq. 0) then
             k1e5=0
             jcount=0
             write(ch4,'(1x,i2)') k1e5  
             return
         endif
c
         jcount=jcount+1
         if (jcount .eq. 100000) then
             jcount=0
             k1e5=k1e5+1
             write(ch4,'(1x,i2)') k1e5  
         endif
            if (mod(jcount,inc) .eq. 0) then
*                write(*,'(i5\)') jcount
cUNiX           write(*,'(i5$)') jcount
cUNiX           call flush(6)
                icarriage=icarriage+1
                if (icarriage .eq. 15) then
                    write(*,'(/a,a\)') ch4,':'   
*                   write(*,'(/a\)') ' @@ '
cUNiX               write(*,'(/a$)') ' @@ '
                    icarriage=0
                endif
            endif
      return
      end
c
c
!///////////Fangbao Tian commented this subroutine
!!      subroutine zztimer(ilog,message) 
!!         use ifport
!!         integer*2 ihr,imin,isec,i100
!!         character*(*) message
cUNiX       real tttt(2)
cUNiX       call etime(tttt)
cUNiX       write(*,*) message,tttt(1),tttt(2)
cMS         call gettim(ihr,imin,isec,i100)
cMS              it1 = 3600*ihr+60*imin+isec
cMS              write(ilog,10) message, it1,i100
cMS              write(*,10) message, it1,i100
!!            call gettim(ihr,imin,isec,i100)
!!                 it1 = 3600*ihr+60*imin+isec
!!                 write(ilog,11) message, it1,i100
**               write(*,11) message, it1,i100
!!            iticks=0
cL          call timer(iticks)
cL               write(ilog,10) message, iticks/100
cL               write(*,10) message, iticks/100
!! 11              format(4h @@ ,a,1x,i10,1h.,i2,4h sec)
!! 10              format(4h @@ ,a,1x,i10,4h sec)
!!      return
!!      end
c
c
!//////Fangbao Tian commented this sub.
!!      subroutine zzflush(i_unit) 
*        use portlib
!!         use ifport
!!         use ifcore
!!         integer i_unit 
!!         logical result
!!         result = commitqq(i_unit)
*        call flush(i_unit)
cUNiX    call flush(ifyl)
cWAT     istatus = flushunit(ifyl)
!!      return
!!      end
c
c
c
c FUNCTION INTGET
      integer function intget(min, max)
c
      do 10 i= 1, 1
           read(ikbd ,*) line
           write(*,*) line
           if (line .lt. min) then
               write(*,*)' **** Minimum answer is ',min,', using This .'
               intget = min 
               return
           elseif (line .gt. max) then
               write(*,*)' **** Maximum answer is ',max,', using This .'
               intget = max 
               return
           else
               intget = line
               return
           endif
 10   continue
      write(*   ,*)'@@ Can not get acceptable answer. Assuming 0.'
      intget = 0
c
         return
         end
c
c
c FUNCTION REALGET
      real function realget(min, max)
c
      real min, max, line
c
         do 10 i= 1, 1
            read(ikbd   ,*) line
            write(*,*) line
            if (line .lt. min) then
               write(*,*)' **** Minimum answer is ',min,', use This .'
               realget=min
            else if (line .gt. max) then
               write(*,*)' **** Maximum answer is ',max,', use This .'
               realget=max
            else
               realget = line
               return
            endif
 10      continue
         write(*,*)'@@ Can not get acceptable answer. Assuming 0.0 .'
         realget = 0.0
c
         return
         end
c
c
!////////////Fangbao commented this sub.
!!        subroutine timday(ilog)
!!          use ifport
cUNiX     character*50 timdat
cL        character*11 str11
cL        character*8 str8
!!            integer*2 ihr,imin,isec,i100,iyr,imon,iday
cMS         call getdat(iyr,imon,iday)
cMS         call gettim(ihr,imin,isec,i100)
cMS         write(ilog,10) imon,iday,iyr,ihr,imin
!!            call getdat(iyr,imon,iday)
!!            call gettim(ihr,imin,isec,i100)
!!            write(ilog,10) imon,iday,iyr,ihr,imin
!! 10      format(' @@ DATE:',1x,i2,1h-,i2,1h-,i4,6x,'TIME:',1x,i2,1h:,i2)
c
cUNiX     call system('date > time.tmp')
cUNiX     open(unit=29,file='time.tmp')
cUNiX     rewind 29
cUNiX     read(29,'(1a50)') timdat
cUNiX     write(ilog,*)'@@ ', timdat
c
cL        call date(str8)
cL        call time(str11)
cL        write(ilog,*)'@@ DATE: ',str8,'  TIME: ',str11
c
!!          return
!!          end
c
      subroutine user
cUNiX    call system('uid >> /gn5/doyle/users/stridyn.usr')
cUNiX    call system('date >> /gn5/doyle/users/stridyn.usr')
      return
      end
c
c
c
      subroutine ikayex
         character*1 ch
         write(*,*)' '
       write(*,*)'****************************************************'
       write(*,*)'                                                   '
       write(*,*)'                                                   '
       write(*,*)'               COMPANY POLICY                      '
       write(*,*)'  '
       write(*,*)' The ikayex SOFTWARE TOOLS company is dedicated to '
       write(*,*)' providing low cost, high quality software for the '
       write(*,*)' analysis of structural dynamics problems. Our pol-'
       write(*,*)' icy is to price the software (with good document- '
       write(*,*)' ation and support) in the range of the cost of a  '
       write(*,*)' good technical book.                              '
       write(*,*)'  '
       write(*,*)' For a catalog of our products, please write to or'
       write(*,*)' call:'
       write(*,*)'  '
       write(*,*)'            ikayex SOFTWARE TOOLS '
       write(*,*)'               615 Elston Road     '
       write(*,*)'           Lafayette, Indiana 47905'
       write(*,*)'                (317) 477-6103    '
       write(*,*)'  '
       write(*,*)'  '
       write(*,*)'****************************************************'
       write(*,*)'PRESS   <<return>>    to continue'
       read(ikbd,'(1a1)') ch
       write(*,*)' '
       write(*,'(///////////)')
         return
         end
c
c
       integer function lentrim(string)
           character*(*) string
c
           ihi=len(string)
           do i=ihi,1,-1
              if (string(i:i) .ne. ' ') then
                  lentrim=i
                  return
              endif
           enddo
           lentrim=0
           write(*,*)'!!!! warning  zero string !!!'
c
        return
        end
c
c
      subroutine strsize(string,ilow,ihi)
         character*(*) string
c
          ihi=len_trim(string)
          do i=1,ihi
             if (string(i:i) .ne. ' ') goto 91
          enddo
 91       continue
          ilow=i
c
      return
      end
c
c 
cccccc   WATCOM FIX
       subroutine bugfix( )
          common /ibug/ ibugfix
          data    ibugfix /-1/
       return
       end
