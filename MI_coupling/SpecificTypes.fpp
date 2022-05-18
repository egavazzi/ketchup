module SpecificTypes

  type dist2v
     ! ivoffset contains the number of dvz steps the velocity vector 
     ! has been shifted.
     integer :: ivzoffset
     double precision, allocatable :: f(:,:), gamma(:,:)
  end type dist2v

  type species
   ! read from input file
     integer Nvz, Nmu
     double precision vzmin, vzmax, mumin, mumax, muexp, mass, charge, &
          n0, vz0, kTz, kTp
     double precision n0L, vz0L, kTzL, kTpL, n0R, vz0R, kTzR, kTpR
     logical vzshifting, Tfromfile, vfromfile, lossconeL, lossconeR, &
          relativistic
   ! calculated
     double precision dvz
     ! Vz is the initial velocity vector comprising Nvz points from 
     ! vz0+vzmin+dvz/2 to vz0+vzmax-dvz/2, where 
     ! dvz=(vzmax-vzmin)/Nvz
     ! mu is the magnetic moment vector comprising Nmu points. 
     ! The mu-vector is generally non-uniform, unless muexp=1.
     double precision, allocatable :: Vz(:), mu(:), dmu(:), maxf0(:)
     ! Each element of the node array contains the velocity distribution
     ! for the corresponding point in xi.
     type(dist2v), allocatable :: node(:)
  end type species

  type bufferedDistribution
     ! This data type contains one vector of offsets and one vector 
     ! of distributions. When exchanging ghost points, data from 
     ! two nodes will have to be packed into one of these structures.
     integer, allocatable :: ivzoffsets(:)
     double precision, allocatable :: fs(:)
  end type bufferedDistribution

end module SpecificTypes
