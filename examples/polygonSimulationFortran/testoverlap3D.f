      subroutine testoverlap3D(p1,p2,overlap)
      
      use commondata
      
      integer ii,jj,kk,i,j,k, overlap, p1,p2
      double precision Ra, Rb, dl, L(3),d(3),norm
      
      
      overlap=1
      
      do i=1,3
      d(i)=r_cubenew(i,p1)-r_cube(i,p2)
      end do
c     first for the axis of partile p1      
      do i=1,3
          Ra=0.0
          Rb=0.0
          do j=1,3
              do k=1,3
          Ra=Ra+abs(unew(k,j,p1)*unew(k,i,p1))
          Rb=Rb+abs(u(k,j,p2)*unew(k,i,p1))
              end do
          end do
          dl=0
          do k=1,3
          dl=dl+d(k)*unew(k,i,p1)
          end do
      if(abs(dl) .ge. a*(Ra+Rb)) then
      overlap=0
      return
      end if
      end do

c     now for the particle axis p2
            do i=1,3
          Ra=0.0
          Rb=0.0
          do j=1,3
              do k=1,3
          Ra=Ra+abs(unew(k,j,p1)*u(k,i,p2))
          Rb=Rb+abs(u(k,j,p2)*u(k,i,p2))
              end do
          end do
          dl=0
          do k=1,3
          dl=dl+d(k)*u(k,i,p2)
          end do
      if(abs(dl) .ge. a*(Ra+Rb)) then
      overlap=0
      return
      end if
      end do
      
c     finally for the cross particle axis
      do i=1,3
          do ii=1,3
      L(1)=unew(2,i,p1)*u(3,ii,p2)-unew(3,i,p1)*u(2,ii,p2)
      L(2)=unew(3,i,p1)*u(1,ii,p2)-unew(1,i,p1)*u(3,ii,p2)
      L(3)=unew(1,i,p1)*u(2,ii,p2)-unew(2,i,p1)*u(1,ii,p2)
      norm=sqrt(L(1)**2+L(2)**2+L(3)**2)
      if(norm .ne. 0.0) then
      
          
          L=L/norm
      
          Ra=0.0
          Rb=0.0
          do j=1,3
              do k=1,3
          Ra=Ra+abs(unew(k,j,p1)*L(k))
          Rb=Rb+abs(u(k,j,p2)*L(k))
              end do
          end do
          dl=0
          do k=1,3
          dl=dl+d(k)*L(k)
          end do
          crit=a*(Ra+Rb)
      if(abs(dl) .ge. a*(Ra+Rb)) then
      overlap=0
      return
      end if
      end if
      end do
      end do
      
      
      end
      