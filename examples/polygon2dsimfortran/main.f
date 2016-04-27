        program main
      
      use commondata
      
      implicit none
      
      integer overlap,idummy,ncycle, cycleindex
      double precision drmax, rotmax, ranf, dist, move,lambda,dum,
     + rg, psi, rmin, psi_loc(nmax),drift, randdisp,
     + D0, Dr, dt 
      integer i,k,j,step,nstep,iprint,p1,p2,jj,index_x,index_y,ii,kk
      character filetemp*4

      character configfile*50
      open(1,file='run.txt')
      
      read(1,*)
      read(1,*) np
      read(1,*)
      read(1,*) polygon
      read(1,*)
      read(1,*) nstep
      read(1,*)      
      read(1,*)  iprint
      read(1,*)
      read(1,*) lambda
      read(1,*)
      read(1,*)  boxlenx, boxleny
      read(1,*)
      read(1,*) idummy
      read(1,*)
      read(1,*)  rmin
      read(1,*)
      read(1,*)  dt
      read(1,*)
      read(1,*)  D0
      read(1,*)
      read(1,*)  Dr
      read(1,*) 
      read(1,*) ncycle
      read(1,*)
      read(1,'(a50)')  configfile
      
      partwidth_x=boxlenx/(boxlenx/2.0-1.0)
      partwidth_y=boxleny/(boxleny/2.0-1.0)
      partnmax_x=boxlenx/partwidth_x
      partnmax_y=boxleny/partwidth_y
      partlimit_x=-boxlenx*0.5
      partlimit_y=-boxleny*0.5
      
       do cycleindex=1,ncycle
           if(cycleindex .lt. 10) then
               write(filetemp,'(i1)') cycleindex
           elseif(cycleindex .lt. 100) then
               write(filetemp,'(i2)') cycleindex
           elseif(cycleindex .lt. 1000) then
               write(filetemp,'(i3)') cycleindex
           else
               write(filetemp,'(i4)') cycleindex
           end if
           
            open(10,file='mc_xyz'//trim(filetemp)//'.txt')
            open(50,file='mc_order'//trim(filetemp)//'.txt')
      
            a=sin(0.5*pi-pi/polygon)
      
            open(20,file=configfile)
      
            do i=1,np
                read(20,*) dum,r_cube(1,i),r_cube(2,i),dum,phi(i)
                r_cube(1,i)=r_cube(1,i)+0.01
                r_cube(2,i)=r_cube(2,i)+0.02
            end do
            close(20)
c         generate edge center
      
            angle=2.0/polygon
            edgelength=a*tan(pi/polygon)
      
            do i=1,np

                do k=1,polygon
        edge(1,k,i)=r_cube(1,i)+a*cos(phi(i)-0.5*pi+(k-1)*angle*pi)
        edge(2,k,i)=r_cube(2,i)+a*sin(phi(i)-0.5*pi+(k-1)*angle*pi)
                end do

            end do
            edgenew=edge
            phinew=phi
      
c     first test overlapping of initial configurations
            overlap=0
            do p1=1,np-1
                do p2=p1+1,np
                    call testoverlap(p1,p2,overlap)
                    if(overlap .eq. 1) then
                        write(*,*) 'overlap in initial config'
                    end if
                end do
            end do
      
      
c         first parition the space
      

            partcount=0
            partlist=0
            do i=1,Np
        
        partindex_x=int((r_cube(1,i)-partlimit_x)/partwidth_x)+1
        partindex_y=int((r_cube(2,i)-partlimit_y)/partwidth_y)+1
        
            if(partindex_x .gt. partnmax_x) then
                partindex_x=1
            end if

            if(partindex_y .gt. partnmax_y) then
                partindex_y=1
            end if
            partcount(partindex_x,partindex_y)=
     +      partcount(partindex_x,partindex_y)+1
            partlist(partcount(partindex_x,partindex_y),
     +      partindex_x,partindex_y)=i
        
            partrecord(1,i)=partindex_x
            partrecord(2,i)=partindex_y

            end do
            open(99,file='text.txt')
        
        do i=1,partnmax_x
            do j=1,partnmax_y
                do k=1,partcount(i,j)
        write(99,*) i,j, partlist(k,i,j),partcount(i,j),
     +   partrecord(1,partlist(k,i,j)),partrecord(2,partlist(k,i,j))
                end do
            end do
        end do
        close(99)
      
      
          
      do step=1,nstep  
          
          
          
          do i=1,np
           if(step .eq. 107 .and. i .eq. 58) then
              write(*,*)
               open(99,file='text.txt')
        
        do ii=1,partnmax_x
        do jj=1,partnmax_y
        do kk=1,partcount(ii,jj)
        write(99,*) ii,jj, partlist(kk,ii,jj),partcount(ii,jj),
     + partrecord(1,partlist(kk,ii,jj)),partrecord(2,partlist(kk,ii,jj))
        end do
        end do
        end do
        close(99)
          end if            
              dist=sqrt(r_cube(1,i)**2+r_cube(2,i)**2)
              do j=1,2
                  
      drift=-(r_cube(j,i))*lambda*dt*D0
         randdisp=(ranf(idummy)-0.5)*sqrt(D0*24.0*dt)
          r_cubenew(j,i)=r_cube(j,i)+drift+randdisp
        
              end do
      phinew(i)=phi(i)+(ranf(idummy)-0.5)*sqrt(Dr*24.0*dt)
      
         do k=1,polygon
      edgenew(1,k,i)=r_cubenew(1,i)+
     + a*cos(phinew(i)-0.5*pi+(k-1)*angle*pi)
      edgenew(2,k,i)=r_cubenew(2,i)+
     + a*sin(phinew(i)-0.5*pi+(k-1)*angle*pi)
         end do          
    
      overlap=0
      
c     now find the neighbor of particle i

 
          
        do index_x=-1,1,1
            do index_y=-1,1,1
              do jj=1,partcount(index_x+partrecord(1,i)
     +             ,index_y+partrecord(2,i))
          p2=partlist(jj,index_x+partrecord(1,i)
     +     ,index_y+partrecord(2,i))
             
            if(p2 .ne. i) then
                
       dist=sqrt((r_cubenew(1,i)-r_cube(1,p2))**2.0+
     +  (r_cubenew(2,i)-r_cube(2,p2))**2.0)
       if(dist .le. 2.0*a) then
           overlap=1
           goto 199
       end if
       
         
      call testoverlap(i,p2,overlap)
         
      if(overlap .eq. 1) Then
         goto 199
      end if
             end if
         end do
       end do
      end do   
          
199      if(overlap .eq. 0 )then
          
          phi(i)=phinew(i)
               do j=1,2                 
          r_cube(j,i)=r_cubenew(j,i)
          do k=1,polygon
                  edge(j,k,i)=edgenew(j,k,i)
          end do    
               end do
               
c     first pop out particle i
      do jj=1,partcount(partrecord(1,i),partrecord(2,i))
      p2=partlist(jj,partrecord(1,i),partrecord(2,i))
      if(p2 .eq. i) then
      partlist(jj,partrecord(1,i),partrecord(2,i))=
     + partlist(partcount(partrecord(1,i),partrecord(2,i)),
     + partrecord(1,i),partrecord(2,i))
      end if
      end do
      partcount(partrecord(1,i),partrecord(2,i))=
     + partcount(partrecord(1,i),partrecord(2,i))-1
               
c         now i update the partition list
        partindex_x=int((r_cube(1,i)-partlimit_x)/partwidth_x)+1
        partindex_y=int((r_cube(2,i)-partlimit_y)/partwidth_y)+1
        
        if(partindex_x .gt. partnmax_x) partindex_x=1
        if(partindex_y .gt. partnmax_y) partindex_y=1
        
        partcount(partindex_x,partindex_y)=
     +   partcount(partindex_x,partindex_y)+1
        partlist(partcount(partindex_x,partindex_y),
     +   partindex_x,partindex_y)=i
        
        partrecord(1,i)=partindex_x
        partrecord(2,i)=partindex_y
             
                   
      end if
          
  
     
          end do
          
          
          if(mod(step,iprint) .eq. 0 )then
      call calorder(rg, psi,psi_loc,rmin) 
          write(50,'(3f14.5)') dble(step)*dt/1000.0, psi, rg
              do i=1,np
          write(10,'(i5,5f12.5)') i, r_cube(1,i),r_cube(2,i)
     +          ,0.0,phi(i), psi_loc(i)
              end do
          end if
          
      end do
      
         
      end do
      
      end
