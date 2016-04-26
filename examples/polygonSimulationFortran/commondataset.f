      module commondata
      
      
      integer, parameter :: nmax=1000,nmax2=50
      double precision phi(nmax), r_cube(3,nmax), phinew(nmax),
     +   edge(2,6,nmax),edgenew(2,6,nmax),r_cubenew(3,nmax),
     + u(3,3,nmax), unew(3,3,nmax) 
      double precision a, polygon, boxlen,angle,edgelength
      integer np,np3
      integer nmaxpart,partlist(nmax2,nmax2,nmax2,nmax2),
     + partindex_x,partindex_y,partlimit_x,partlimit_y,partlimit_z,
     + partcount(nmax2,nmax2,nmax2),partnmax_x,partnmax_y, partnmax_z,
     + partrecord(3,nmax),partindex_z
      double precision partwidth_y,partwidth_x,partwidth_z,
     + boxlenx, boxleny, boxlenz
      
      
      double precision, parameter:: pi=3.1415926
      
      
      
      end module commondata