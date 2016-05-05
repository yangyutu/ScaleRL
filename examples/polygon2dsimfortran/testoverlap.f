      subroutine testoverlap(p1,p2,overlap)
       use commondata
      
      integer overlap,j,k,p1,p2
      double precision b1,b2,det,det1,det2,edgetemp(2,6),phitemp,
     + dist 
      
      
      do k=1,polygon
          do j=1,2
      edgetemp(j,k)=edgenew(j,k,p1)
          end do
      end do
      
      phitemp=phinew(p1)
      

      do j=1,polygon
          do k=1,polygon
      
      
      b1=edge(1,k,p2)-edgetemp(1,j)
      b2=edge(2,k,p2)-edgetemp(2,j)
      dist=sqrt(b1**2.0+b2**2.0)
!      if(dist .ge. 2.0*edgelength) then
!          goto 555
!      end if
      
      
      det=-cos(phitemp+(j-1)*angle*pi)*sin(phi(p2)+(k-1)*angle*pi)+
     +  cos(phi(p2)+(k-1)*angle*pi)*sin(phitemp+(j-1)*angle*pi)  
      
      det1=-b1*sin(phi(p2)+(k-1)*angle*pi)+
     + b2*cos(phi(p2)+(k-1)*angle*pi)
      
      
      det2=b2*cos(phitemp+(j-1)*angle*pi)-
     + b1*sin(phitemp+(j-1)*angle*pi)
      
      if(det .ne. 0.0) then
      t1=det1/det
      t2=det2/det
!      write (*,*) p1, p2, j,k,b1,b2,dist,det,det1,det2,t1,t2
      if(abs(t1) .le. edgelength .and. abs(t2) .le. edgelength) then
          overlap=1
          return
      end if
      
      end if
555   continue      
          end do
      end do
      
      end
      
      
      
