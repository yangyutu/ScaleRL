      subroutine calorder(rg, psi,psi_loc, rmin)
      
      use commondata
      
      
      double precision rg, psi, xmean, ymean,rxij, ryij
      double precision RP, rmin,accumpsir,accumpsii,scale
      integer i,j,k, nb(nmax)
      double precision psir(nmax), psii(nmax),psi_loc(nmax)
      
      psir=0.0
	psii=0.0
      accumpsii=0
      accumpsir=0
      xmean=0.0
      ymean=0.0
      rg=0

      do i=1,np
          xmean=xmean+r_cube(1,i)
          ymean=ymean+r_cube(2,i)
      end do
      
      xmean=xmean/dble(np)
      ymean=ymean/dble(np)
      
      do i=1,np
      rg=rg+(r_cube(1,i)-xmean)**2.0+(r_cube(2,i)-ymean)**2.0
      end do
      
      rg=rg/dble(np)
      rg=sqrt(rg)
      !  GET PSI6 FOR ALL PARTICLES
            if(polygon .eq. 3) then
          scale=6
      elseif(polygon .eq. 4) then
          scale=4
      elseif(polygon .eq. 6) then
          scale=3
      else
          scale=1
      end if
      nb=0
	do i=1, NP
		do j=1, NP
			if (i.ne.j) then 
				rxij=r_cube(1,j)-r_cube(1,i)		
				ryij=r_cube(2,j)-r_cube(2,i)
				
				RP=sqrt(RXIJ**2+RYIJ**2)
				if (RP.le.rmin) then
					nb(i)=nb(i)+1
					
					psir(i)=psir(i)+cos(scale*phi(j))
					psii(i)=psii(i)+sin(scale*phi(j))
				end if
			end if
		end do
		if (nb(i).ne.0) then
			psir(i)=psir(i)/dble(nb(i))
			psii(i)=psii(i)/dble(nb(i))
			
		end if
	end do

       
          psi=0
		accumpsir=0
		accumpsii=0
		do i=1,np
              psi_loc(i)=sqrt(psir(i)**2+psii(i)**2)
              
		accumpsir=accumpsir+psir(i)
		accumpsii=accumpsii+psii(i)
		end do
		accumpsir=accumpsir/np
		accumpsii=accumpsii/np
		
		psi=sqrt(accumpsir**2+accumpsii**2)
      
      
      end
      