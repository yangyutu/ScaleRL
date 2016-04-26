      subroutine rotate3D(axis, p1, rotv, randphi, v)
      use commondata
      double precision w, x, y, z, mat(3,3), v(3), rotv(3), randphi
      integer i,j, axis, p1
      
      w=cos(0.5*randphi)
      x=rotv(1)*sin(0.5*randphi)
      y=rotv(2)*sin(0.5*randphi)
      z=rotv(3)*sin(0.5*randphi)
      mat(1,1) = 1.d0-2.d0*(y**2+z**2)
      mat(1,2) = 2.d0*(x*y-w*z)
      mat(1,3) = 2.d0*(x*z+w*y)
      mat(2,1) = 2.d0*(x*y+w*z)
      mat(2,2) = 1.d0-2.d0*(x**2+z**2)
      mat(2,3) = 2.d0*(y*z-w*x)
      mat(3,1) = 2.d0*(x*z-w*y)
      mat(3,2) = 2.d0*(y*z+w*x)
      mat(3,3) = 1.d0-2.d0*(x**2+y**2)
      
      v=0.0
      do i=1,3
          do j=1,3
              v(i)=v(i)+mat(i,j)*u(j,axis,p1)
          end do
      end do
      
      
      end