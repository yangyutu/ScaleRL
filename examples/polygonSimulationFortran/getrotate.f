      subroutine getrotate(p1, xangle, yangle, zangle)
      use commondata
      
      integer i,j,k, p1
      
      double precision trans(3,3), xangle, yangle, zangle
      
      do i=1,3
          do j=1,3
      trans(i,j)=u(j,i,p1)
          end do
      end do
      
      
      xangle=atan2(trans(3,2), trans(3,3))
      
      yangle=atan2(-trans(3,1), sqrt(trans(3,2)**2+trans(3,3)**2))   
      
      zangle=atan2(trans(2,1), trans(1,1))
      
      xangle=xangle/pi*180.0
        yangle=yangle/pi*180.0
          zangle=zangle/pi*180.0
      
      
      
      
      
      end
      
      
      
      