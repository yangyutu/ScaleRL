      module commondata
      
      
      integer, parameter :: nmax=300,nmax2=50
      double precision phi(nmax), r_cube(2,nmax), phinew(nmax),
     +   edge(2,6,nmax),edgenew(2,6,nmax),r_cubenew(2,nmax)
      double precision a, boxlen,angle,edgelength
      integer np,np3, polygon
      integer nmaxpart,partlist(nmax2,nmax2,nmax2),
     + partindex_x,partindex_y,partlimit_x,partlimit_y,
     + partcount(nmax2,nmax2),partnmax_x,partnmax_y, 
     + partrecord(2,nmax)
      double precision partwidth_y,partwidth_x, boxlenx, boxleny
      
      
      double precision, parameter:: pi=3.1415926
      
      
      
      end module commondata
