        program main
      
      use commondata
      
      implicit none
      
      integer overlap,idummy,ncycle, cycleindex
      double precision drmax, rotmax, ranf, dist, move,lambda,dum,
     + rg, psi, rmin, psi_loc(nmax),drift, randdisp,
     + D0, Dr, dt,randphi1,randphi2,
     + rotv(3), v(3), leftcorner(3), rightcorner(3),
     + xangle, yangle, zangle,norm, S4max, G4(nmax)
      integer i,k,j,step,nstep,iprint,p1,p2,jj,index_x,index_y,ii,kk,
     + m, index_z 
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
c      read(1,*)  drmax,rotmax
c      read(1,*) 
      read(1,*)  boxlenx, boxleny, boxlenz
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
      partwidth_z=boxlenz/(boxlenz/2.0-1.0)
      partnmax_x=boxlenx/partwidth_x
      partnmax_y=boxleny/partwidth_y
      partnmax_z=boxlenz/partwidth_z
      partlimit_x=-boxlenx*0.5
      partlimit_y=-boxleny*0.5
      partlimit_z=-boxlenz*0.5
      
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
      
      
      
      open(20,file=configfile)
c      u=0.0
      do i=1,np
          read(20,*) dum,r_cube(1,i),r_cube(2,i),r_cube(3,i), dum,
     +     dum,dum,dum, u(1,1,i),u(2,1,i),u(3,1,i), u(1,2,i),
     +     u(2,2,i),u(3,2,i),
     + u(1,3,i), u(2,3,i), u(3,3,i) 
c          u(1,1,i)=1
c          u(2,2,i)=1
c          u(3,3,i)=1

      end do
      close(20)
c         generate edge center
      
      a=1/sqrt(3.0)
c      angle=2.0/polygon
      edgelength=a
      r_cubenew=r_cube
      unew=u
c     first test overlapping of initial configurations
      overlap=0
      do p1=1,np-1
          do p2=p1+1,np
              call testoverlap3D(p1,p2,overlap)
              if(overlap .eq. 1) then
                  write(*,*) 'overlap in initial config', p1, p2
              end if
          end do
      end do
      
      
c         first parition the space
      

         partcount=0
        partlist=0
        do i=1,Np
        
        partindex_x=int((r_cube(1,i)-partlimit_x)/partwidth_x)+1
        partindex_y=int((r_cube(2,i)-partlimit_y)/partwidth_y)+1
        partindex_z=int((r_cube(3,i)-partlimit_z)/partwidth_z)+1
        
        if(partindex_x .gt. partnmax_x) partindex_x=1       
        if(partindex_y .gt. partnmax_y) partindex_y=1
         if(partindex_z .gt. partnmax_z) partindex_z=1
        
        partcount(partindex_x,partindex_y, partindex_z)=
     +   partcount(partindex_x,partindex_y, partindex_z)+1
        partlist(partcount(partindex_x,partindex_y, partindex_z),
     +   partindex_x,partindex_y,partindex_z)=i
        
        partrecord(1,i)=partindex_x
        partrecord(2,i)=partindex_y
      partrecord(3,i)=partindex_z
        end do
       open(99,file='text.txt')
        
        do i=1,partnmax_x
        do j=1,partnmax_y
        do k=1,partnmax_z
        do m=1,partcount(i,j,k)
        write(99,'(9i4)') i,j,k,m, partlist(m,i,j,k),partcount(i,j,k),
     +partrecord(1,partlist(m,i,j,k)),partrecord(2,partlist(m,i,j,k)),
     + partrecord(3,partlist(m,i,j,k))   
        end do
        end do
        end do
        end do
        close(99)
      
      
          
      do step=1,nstep  
                   
          do i=1,np

              do j=1,3
      drift=-(r_cube(j,i))*lambda*dt*D0
         randdisp=(ranf(idummy)-0.5)*sqrt(D0*24.0*dt)
          r_cubenew(j,i)=r_cube(j,i)+drift+randdisp    
              end do
              
              
      randphi1=(ranf(idummy)-0.5)*sqrt(Dr*24.0*dt)
       randphi2=(ranf(idummy)-0.5)*sqrt(Dr*24.0*dt)
 
c       randphi1=0.0
c       randphi2=0.0
c      we randomly pick an axis and rot
        
           rotv(1)= ranf(idummy)
           rotv(2)=ranf(idummy)
           rotv(3)= ranf(idummy)
      norm=sqrt(rotv(1)**2.0+rotv(2)**2.0+rotv(3)**2.0)
      rotv=rotv/norm
c         aroung axis  rot axis 1, 2 and 3
           call rotate3D(1, i, rotv, randphi1, v) 
           unew(:,1,i)=v
      call rotate3D(2, i, rotv, randphi1, v)      
           unew(:,2,i)=v
c      call rotate3D(3, i, rotv, randphi1, v)     
c           unew(:,3,i)=v
c     
c      write(*,*) unew(1,1,i)*unew(1,2,i)+unew(2,1,i)*unew(2,2,i)
c     +  +unew(3,1,i)*unew(3,2,i)
c         around new axis 2 rot axis 1 and 2           
  
           
           
c           write(*,*) unew(1,1,i)*unew(1,2,i)+unew(2,1,i)*unew(2,2,i)
c     +  +unew(3,1,i)*unew(3,2,i)
c      write(*,*) unew(1,1,i)*unew(1,2,i)+unew(2,1,i)*unew(2,2,i)+
c     + unew(3,1,i)*unew(3,2,i)
      unew(1,3,i)=unew(2,1,i)*unew(3,2,i)-unew(3,1,i)*unew(2,2,i)
      unew(2,3,i)=unew(3,1,i)*unew(1,2,i)-unew(1,1,i)*unew(3,2,i)
      unew(3,3,i)=unew(1,1,i)*unew(2,2,i)-unew(2,1,i)*unew(1,2,i)
c          write(*,*) unew(1,1,i)*unew(1,2,i)+unew(2,1,i)*unew(2,2,i)
c     +  +unew(3,1,i)*unew(3,2,i)

      
c     now find the neighbor of particle i

      do index_x=-1,1,1
            do index_y=-1,1,1
                do index_z=-1,1,1
              do jj=1,partcount(index_x+partrecord(1,i)
     +        ,index_y+partrecord(2,i),index_z+partrecord(3,i))
          p2=partlist(jj,index_x+partrecord(1,i)
     +     ,index_y+partrecord(2,i),index_z+partrecord(3,i))
            if(p2 .ne. i) then
                
       dist=sqrt((r_cubenew(1,i)-r_cube(1,p2))**2.0+
     +  (r_cubenew(2,i)-r_cube(2,p2))**2.0+
     +  (r_cubenew(3,i)-r_cube(3,p2))**2.0)
       
       
       if(dist .le. 2.0*a) then
           overlap=1
           goto 199
       end if
        if(dist .ge. 2.0) then
           overlap=0
           goto 198
       end if
         
         
      call testoverlap3D(i,p2,overlap)
         
      if(overlap .eq. 1) Then
         goto 199
      end if
             end if
198         end do
                end do
            end do
      end do
      
        
199   if(overlap .eq. 0 )then
          
          
          do j=1,3                
          r_cube(j,i)=r_cubenew(j,i)
          do k=1,3
          u(k,j,i)=unew(k,j,i)
          end do    
          end do
          
  
          
               
               
c     first pop out particle i
      do jj=1,partcount(partrecord(1,i),partrecord(2,i),partrecord(3,i))
      p2=partlist(jj,partrecord(1,i),partrecord(2,i),partrecord(3,i))
      if(p2 .eq. i) then
      partlist(jj,partrecord(1,i),partrecord(2,i),partrecord(3,i))=
     + partlist(partcount(partrecord(1,i),partrecord(2,i),
     + partrecord(3,i)),
     + partrecord(1,i),partrecord(2,i),partrecord(3,i))
      end if
      end do
      partcount(partrecord(1,i),partrecord(2,i),partrecord(3,i))=
     + partcount(partrecord(1,i),partrecord(2,i),partrecord(3,i))-1
               
c         now i update the partition list
             
              partindex_x=int((r_cube(1,i)-partlimit_x)/partwidth_x)+1
        partindex_y=int((r_cube(2,i)-partlimit_y)/partwidth_y)+1
        partindex_z=int((r_cube(3,i)-partlimit_z)/partwidth_z)+1
        
        if(partindex_x .gt. partnmax_x) partindex_x=1       
        if(partindex_y .gt. partnmax_y) partindex_y=1
         if(partindex_z .gt. partnmax_z) partindex_z=1
        
        partcount(partindex_x,partindex_y, partindex_z)=
     +   partcount(partindex_x,partindex_y, partindex_z)+1
        partlist(partcount(partindex_x,partindex_y, partindex_z),
     +   partindex_x,partindex_y,partindex_z)=i
        
        partrecord(1,i)=partindex_x
        partrecord(2,i)=partindex_y
       partrecord(3,i)=partindex_z
                   
      end if
          
  
     
          end do
      
          
          if(mod(step,iprint) .eq. 0 )then
      call calorder(rg, S4max,G4,rmin) 
          write(50,'(3f14.5)') dble(step)*dt/1000.0, S4max, rg
              do i=1,np

          call getrotate(i,xangle, yangle, zangle)
c      write(*,*)(leftcorner(1)-rightcorner(1))**2+(leftcorner(2)-
c     + rightcorner(2))**2+(leftcorner(3)-rightcorner(3))**2
          write(10,'(i5,16f12.5)') i, r_cube(1,i),r_cube(2,i)
     +          ,r_cube(3,i), xangle, yangle, zangle, G4(i),
     +     u(1,1,i),u(2,1,i),u(3,1,i), u(1,2,i),u(2,2,i),u(3,2,i),
     + u(1,3,i), u(2,3,i), u(3,3,i)     
              end do
          end if
          
      end do
      
         
      end do
      
      end