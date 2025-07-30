c-------------------------------------------------------
c  write_coord
c  Fangbao Tian, Oct. 2011
c-------------------------------------------------------
       subroutine write_coord(xyzt,maxelem,             ! Fangbao Tian
     &                     npijkm,maxnode,maxxyz,itime)!
!

!       subroutine write_coord(xyzt,npijkm,maxnode,maxxyz,itime)!
        implicit real*8 (a-h,o-z)
        include 'commons.std'

        real*8 :: xyzt(maxxyz)
        integer:: npijkm(maxelem,21)

        character*7  :: fname

        PRINT*,'nonstad: Writing out dump file ... '
        write(fname,'(i7.7)')itime

      OPEN(UNIT=334,FILE=fdn//'s.'//fname,STATUS='UNKNOWN')

      WRITE(334,*)' TITLE = "3D Finite-Element Data"'
      WRITE(334,*)'VARIABLES= "X","Y","Z"'

      write(334,*)'ZONE N=',nnp,' E=',nel,' DATAPACKING=POINT, 
     +               ZONETYPE=FEBRICK'
!     +               ZONETYPE=FETETRAHEDRON'



      do i=1,nnp
         write(334,'(6F15.6)')xyzt(i),xyzt(maxnode+i),xyzt(maxnode*2+i)          
      enddo
         itype=npijkm(1,1)
      do i=1,nel
       write(334,'(8I7)')(npijkm(i,j),j=2,9)  
      enddo      

      CLOSE(334)

      return
      end

c-------------------------------------------------------
c  write_coordHEX20
c  Fangbao Tian, Oct. 2011
c-------------------------------------------------------

       subroutine write_coordHEX20(xyzt,npijkm,maxnode,
     &  maxelem,maxxyz,itime)!
        implicit real*8 (a-h,o-z)
        include 'commons.std'

        real*8 :: xyzt(maxxyz)
        integer:: npijkm(maxelem,21)
        real*8 :: xtemp(nnp+7*nel),ytemp(nnp+7*nel),ztemp(nnp+7*nel)
        integer:: ep_temp(8*nel,8)

        character*7  :: fname

        PRINT*,'nonstad: Writing out dump file ... '

        !!!=========NP
        do i=1,nnp
        xtemp(i)=xyzt(i)
        ytemp(i)=xyzt(maxnode+i)
        ztemp(i)=xyzt(maxnode*2+i)
        enddo
!	J1
        do i=1, nel
        J1=nnp+(i-1)*7+1
        JJ1=npijkm(i,1+1)
        JJ2=npijkm(i,1+2)
        JJ3=npijkm(i,1+3)
        JJ4=npijkm(i,1+4)
        xtemp(J1)=(xyzt(jj1)+xyzt(jj2)+xyzt(jj3)+xyzt(jj4))/4.0
        ytemp(J1)=(xyzt(maxnode+jj1)+xyzt(maxnode+jj2)
     &                 +xyzt(maxnode+jj3)+xyzt(maxnode+jj4))/4.0
        ztemp(J1)=(xyzt(maxnode*2+jj1)+xyzt(maxnode*2+jj2)
     &                 +xyzt(maxnode*2+jj3)+xyzt(maxnode*2+jj4))/4.0

!     J2
        J2=nnp+(i-1)*7+2
        JJ1=npijkm(i,1+13)
        JJ2=npijkm(i,1+14)
        JJ3=npijkm(i,1+15)
        JJ4=npijkm(i,1+16)
        xtemp(J2)=(xyzt(jj1)+xyzt(jj2)+xyzt(jj3)+xyzt(jj4))/4.0
        ytemp(J2)=(xyzt(maxnode+jj1)+xyzt(maxnode+jj2)
     &                 +xyzt(maxnode+jj3)+xyzt(maxnode+jj4))/4.0
        ztemp(J2)=(xyzt(maxnode*2+jj1)+xyzt(maxnode*2+jj2)
     &                 +xyzt(maxnode*2+jj3)+xyzt(maxnode*2+jj4))/4.0

!     J3
        J3=nnp+(i-1)*7+3
        JJ1=npijkm(i,1+5)
        JJ2=npijkm(i,1+6)
        JJ3=npijkm(i,1+7)
        JJ4=npijkm(i,1+8)
        xtemp(J3)=(xyzt(jj1)+xyzt(jj2)+xyzt(jj3)+xyzt(jj4))/4.0
        ytemp(J3)=(xyzt(maxnode+jj1)+xyzt(maxnode+jj2)
     &                 +xyzt(maxnode+jj3)+xyzt(maxnode+jj4))/4.0
        ztemp(J3)=(xyzt(maxnode*2+jj1)+xyzt(maxnode*2+jj2)
     &                 +xyzt(maxnode*2+jj3)+xyzt(maxnode*2+jj4))/4.0

!     J4
        J4=nnp+(i-1)*7+4
        JJ1=npijkm(i,1+5)
        JJ2=npijkm(i,1+6)
        JJ3=npijkm(i,1+1)
        JJ4=npijkm(i,1+2)
        xtemp(J4)=(xyzt(jj1)+xyzt(jj2)+xyzt(jj3)+xyzt(jj4))/4.0
        ytemp(J4)=(xyzt(maxnode+jj1)+xyzt(maxnode+jj2)
     &                 +xyzt(maxnode+jj3)+xyzt(maxnode+jj4))/4.0
        ztemp(J4)=(xyzt(maxnode*2+jj1)+xyzt(maxnode*2+jj2)
     &                 +xyzt(maxnode*2+jj3)+xyzt(maxnode*2+jj4))/4.0

!     J5
        J5=nnp+(i-1)*7+5
        JJ1=npijkm(i,1+2)
        JJ2=npijkm(i,1+3)
        JJ3=npijkm(i,1+7)
        JJ4=npijkm(i,1+6)
        xtemp(J5)=(xyzt(jj1)+xyzt(jj2)+xyzt(jj3)+xyzt(jj4))/4.0
        ytemp(J5)=(xyzt(maxnode+jj1)+xyzt(maxnode+jj2)
     &                 +xyzt(maxnode+jj3)+xyzt(maxnode+jj4))/4.0
        ztemp(J5)=(xyzt(maxnode*2+jj1)+xyzt(maxnode*2+jj2)
     &                 +xyzt(maxnode*2+jj3)+xyzt(maxnode*2+jj4))/4.0

!     J6
        J6=nnp+(i-1)*7+6
        JJ1=npijkm(i,1+3)
        JJ2=npijkm(i,1+4)
        JJ3=npijkm(i,1+8)
        JJ4=npijkm(i,1+7)
        xtemp(J6)=(xyzt(jj1)+xyzt(jj2)+xyzt(jj3)+xyzt(jj4))/4.0
        ytemp(J6)=(xyzt(maxnode+jj1)+xyzt(maxnode+jj2)
     &                 +xyzt(maxnode+jj3)+xyzt(maxnode+jj4))/4.0
        ztemp(J6)=(xyzt(maxnode*2+jj1)+xyzt(maxnode*2+jj2)
     &                 +xyzt(maxnode*2+jj3)+xyzt(maxnode*2+jj4))/4.0

!     J7
        J7=nnp+(i-1)*7+7
        JJ1=npijkm(i,1+1)
        JJ2=npijkm(i,1+4)
        JJ3=npijkm(i,1+8)
        JJ4=npijkm(i,1+5)
        xtemp(J7)=(xyzt(jj1)+xyzt(jj2)+xyzt(jj3)+xyzt(jj4))/4.0
        ytemp(J7)=(xyzt(maxnode+jj1)+xyzt(maxnode+jj2)
     &                 +xyzt(maxnode+jj3)+xyzt(maxnode+jj4))/4.0
        ztemp(J7)=(xyzt(maxnode*2+jj1)+xyzt(maxnode*2+jj2)
     &                 +xyzt(maxnode*2+jj3)+xyzt(maxnode*2+jj4))/4.0
       enddo
      do i=1, nel
      k1=(i-1)*8+1
      k2=(i-1)*8+2
      k3=(i-1)*8+3
      k4=(i-1)*8+4
      k5=(i-1)*8+5
      k6=(i-1)*8+6
      k7=(i-1)*8+7
      k8=(i-1)*8+8
      J1=nnp+(i-1)*7+1
      J2=nnp+(i-1)*7+2
      J3=nnp+(i-1)*7+3
      J4=nnp+(i-1)*7+4
      J5=nnp+(i-1)*7+5
      J6=nnp+(i-1)*7+6
      J7=nnp+(i-1)*7+7
!K1
        ep_temp(k1,1)=npijkm(i,1+ 1);ep_temp(k1,2)=npijkm(i,1+ 9);
        ep_temp(k1,3)=          J1;ep_temp(k1,4)=npijkm(i,1+12);
        ep_temp(k1,5)=npijkm(i,1+13);ep_temp(k1,6)=          J4;
        ep_temp(k1,7)=          J2;ep_temp(k1,8)=          J7;
!K2
        ep_temp(k2,1)=npijkm(i,1+ 9);ep_temp(k2,2)=npijkm(i,1+ 2);
        ep_temp(k2,3)=npijkm(i,1+10);ep_temp(k2,4)=          J1;
        ep_temp(k2,5)=          J4;ep_temp(k2,6)=npijkm(i,1+14);
        ep_temp(k2,7)=          J5;ep_temp(k2,8)=          J2;
!K3
        ep_temp(k3,1)=          J1;ep_temp(k3,2)=npijkm(i,1+10);
        ep_temp(k3,3)=npijkm(i,1+ 3);ep_temp(k3,4)=npijkm(i,1+11);
        ep_temp(k3,5)=          J2;ep_temp(k3,6)=          J5;
        ep_temp(k3,7)=npijkm(i,1+15);ep_temp(k3,8)=          J6;
!K4
        ep_temp(k4,1)=npijkm(i,1+12);ep_temp(k4,2)=          J1;
        ep_temp(k4,3)=npijkm(i,1+11);ep_temp(k4,4)=npijkm(i,1+ 4);
        ep_temp(k4,5)=          J7;ep_temp(k4,6)=          J2;
        ep_temp(k4,7)=          J6;ep_temp(k4,8)=npijkm(i,1+16);
!K5
        ep_temp(k5,1)=npijkm(i,1+13);ep_temp(k5,2)=          J4;
        ep_temp(k5,3)=          J2;ep_temp(k5,4)=          J7;
        ep_temp(k5,5)=npijkm(i,1+ 5);ep_temp(k5,6)=npijkm(i,1+17);
        ep_temp(k5,7)=          J3;ep_temp(k5,8)=npijkm(i,1+20);
!K6
        ep_temp(k6,1)=          J4;ep_temp(k6,2)=npijkm(i,1+14);
        ep_temp(k6,3)=          J5;ep_temp(k6,4)=          J2;
        ep_temp(k6,5)=npijkm(i,1+17);ep_temp(k6,6)=npijkm(i,1+ 6);
        ep_temp(k6,7)=npijkm(i,1+18);ep_temp(k6,8)=          J3;
!K7
        ep_temp(k7,1)=          J2;ep_temp(k7,2)=          J5;
        ep_temp(k7,3)=npijkm(i,1+15);ep_temp(k7,4)=          J6;
        ep_temp(k7,5)=          J3;ep_temp(k7,6)=npijkm(i,1+18);
        ep_temp(k7,7)=npijkm(i,1+ 7);ep_temp(k7,8)=npijkm(i,1+19);
!K8
        ep_temp(k8,1)=          J7;ep_temp(k8,2)=          J2;
        ep_temp(k8,3)=          J6;ep_temp(k8,4)=npijkm(i,1+16);
        ep_temp(k8,5)=npijkm(i,1+20);ep_temp(k8,6)=          J3;
        ep_temp(k8,7)=npijkm(i,1+19);ep_temp(k8,8)=npijkm(i,1+ 8);
      enddo


      write(fname,'(i7.7)')itime

      OPEN(UNIT=334,FILE=fdn//'s.'//fname,STATUS='UNKNOWN')

      WRITE(334,*)' TITLE = "3D Finite-Element Data"'
      WRITE(334,*)'VARIABLES= "X","Y","Z"'

      write(334,*)'ZONE N=',nnp+7*nel,' E=',8*nel,' DATAPACKING=POINT, 
     +               ZONETYPE=FEBRICK'
!    +               ZONETYPE=FETETRAHEDRON'



      do i=1,nnp+7*nel
         write(334,'(4F15.6)')xtemp(i),ytemp(i),ztemp(i)          
      enddo
         itype=npijkm(1,1)
      do i=1,8*nel
       write(334,'(8I7)')(ep_temp(i,j),j=1,8)  
      enddo      

      CLOSE(334)

      return
      end
