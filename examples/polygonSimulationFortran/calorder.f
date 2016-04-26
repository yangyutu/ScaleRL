      subroutine calorder(rg,S4max,G4, rmin)
      
      use commondata
      
      
      double precision rg, psi, xmean, ymean,rxij, ryij, zmean
      double precision RP, rmin,accumpsir,accumpsii,scale,dotv
      integer i,j,k, nb(nmax)
      double precision G4(nmax), n(3), S4(3,nmax),director(3),
     + S4max 
      
      psir=0.0
	psii=0.0
      accumpsii=0
      accumpsir=0
      xmean=0.0
      ymean=0.0
      zmean=0.0
      rg=0
      
      do i=1,np
          xmean=xmean+r_cube(1,i)
          ymean=ymean+r_cube(2,i)
          zmean=zmean+r_cube(3,i)
      end do
      
      xmean=xmean/dble(np)
      ymean=ymean/dble(np)
      zmean=zmean/dble(np)
      
      do i=1,np
      rg=rg+(r_cube(1,i)-xmean)**2.0+(r_cube(2,i)-ymean)**2.0+
     + (r_cube(3,i)-zmean )**2.0
      end do
      
      rg=rg/dble(np)
      rg=sqrt(rg)
      
      do k=1,np
          do m=1,3
              S4(m,k)=0.0
              n=u(:,m,k)           
      do i=1,np
          do j=1,3
      dotv=u(1,j,i)*n(1)+u(2,j,i)*n(2)+u(3,j,i)*n(3)
      S4(m,k)=S4(m,k)+(35*dotv**4.0-30*dotv**2.0+3)
          end do
      end do
      S4(m,k)=S4(m,k)/14.0/dble(np)
          end do
      end do
      
      S4max=-1000
      
      do k=1,np
          do m=1,3
              if(S4max .le. S4(m,k) ) then
                  S4max=S4(m,k)
                  director=u(:,m,k)
              end if
          end do
      end do
      
      
              
              
              nb=0
	do i=1, NP
          G4(i)=0.0
		do j=1, NP
			if (i.ne.j) then 
				rxij=r_cube(1,j)-r_cube(1,i)		
				ryij=r_cube(2,j)-r_cube(2,i)
				rzij=r_cube(3,j)-r_cube(3,i)
				RP=sqrt(RXIJ**2+RYIJ**2+ rzij**2)
				if (RP.le.rmin) then
					nb(i)=nb(i)+1
					do k=1,3
                          do m=1,3
       dotv= u(1,k,i)*u(1,m,j)+u(2,k,i)*u(2,m,j)+u(3,k,i)*u(3,m,j)
      G4(i)=G4(i)+35.0*dotv**4.0-30.0*dotv**2.0+3.0 
                          end do
                      end do
                      
				end if
			end if
		end do
		if (nb(i).ne.0) then
			G4(i)=G4(i)/14.0/3.0/dble(nb(i))
			
		end if
	end do

       
         
      
      
      end
      