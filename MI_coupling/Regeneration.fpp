! This program reads dump files from a ketchup run, creates a new grid, 
! and saves new dump files, from which a new ketchup run can be started.
!
! This program must be run from a working directory where the new, 
! regenerated, run shall reside. This directory must also contain an 
! input file named regen_par.m, which contains the name of the directory 
! where the completed run that we regenerate from is located. It also 
! must contain information on the number of processes that are to be 
! used for the new run. Ordinary input files (inputb5.m) must be present 
! in both the to and from directories.
!------------------------------------------------------------------------

program regenerate_ketchup

  use SpecificTypes

  implicit none

! input data
  integer Niter_old, dump_period_fields_old, fields_per_file_old, &
       dump_period_distr_old, &
       dump_period_distr_1v_old, dump_start_old, &
       shift_test_period_old, Nxi_old, Nz_old, Nspecies_old, BC_Poisson_old, &
       initialiser_old, dump_period_dump_old, dump_period_distr_IonoBoundary_old
  double precision zmin_old, zmax_old, dt_old, resistance_old, voltage_old, &
       voltage_init_old, E0_old
  double precision const_a_old
  character transffilename_old*200, voltagefilename_old*200
  logical startfromdumpfile_old, exitafterdump_old


  integer Niter_new, dump_period_fields_new, fields_per_file_new, &
       dump_period_distr_new, &
       dump_period_distr_1v_new, dump_start_new, &
       shift_test_period_new, Nxi_new, Nz_new, Nspecies_new, BC_Poisson_new, &
       initialiser_new, dump_period_dump_new, dump_period_distr_IonoBoundary_new
  double precision zmin_new, zmax_new, dt_new, resistance_new, voltage_new, &
       voltage_init_new, E0_new
  double precision const_a_new
  character transffilename_new*200, voltagefilename_new*200
  logical startfromdumpfile_new, exitafterdump_new



! grid parameters
  double precision dxi_old, BBC_old(2), dBBC_old(2)
  double precision, allocatable :: Z_old(:), Zcorner_old(:), dz_old(:), &
       XI_old(:), XIcorner_old(:), gp_old(:), gpcorner_old(:), &
       vector_old(:)
  double precision Zcmp_old, Zlow_old
  ! local grid
  integer Nz_local_old
  double precision, allocatable :: Z_local_old(:), dz_local_old(:), &
       XI_local_old(:), XIcorner_local_old(:), gp_local_old(:), &
       B_local_old(:), dB_local_old(:)

  double precision dxi_new, BBC_new(2), dBBC_new(2)
  double precision, allocatable :: Z_new(:), Zcorner_new(:), dz_new(:), &
       XI_new(:), XIcorner_new(:), gp_new(:), gpcorner_new(:), &
       vector_new(:), meanVz_new(:,:)
  ! local grid
  integer Nz_local_new
  double precision, allocatable ::  Z_local_new(:), dz_local_new(:), &
       XI_local_new(:), XIcorner_local_new(:), gp_local_new(:), &
       B_local_new(:), dB_local_new(:)

! particle species
  type(species), allocatable :: particle_old(:)

  type(species), allocatable :: particle_new(:)


! processes and regenereation related variables
  integer Nprocs_old, Nprocs_new, DirNameLength, dstart, dend
  integer, allocatable :: istart_old(:), iend_old(:)
  integer, allocatable :: istart_new(:), iend_new(:)
  character FromDir*200
  double precision rho_reduce

! Profiles and reinitialisation
  double precision, allocatable :: kTzprofile_global(:,:), &
       kTpprofile_global(:,:), vzprofile_global(:,:)
  double precision, allocatable :: nprofile_local(:,:), &
       kTzprofile_local(:,:), kTpprofile_local(:,:), vzprofile_local(:,:)
  integer, allocatable :: rspecies(:)
  logical, allocatable :: reinit(:)

! other variables
  integer ierr, iter_start, Ng, ii, jj, kk, ll, mm, kkk, iz, ivz, imu, doff
  integer Nretries, attempts
  double precision this_new_Z, this_new_vz, this_new_mu, weight
  logical dovz, domu, fileexists
  ! ghost_old is used to store two interior nodes, so that we don't need 
  ! to use ghost points, except for boundary conditions. 
  type(species), allocatable :: ghost_old(:)
  integer, allocatable :: Nvz_new(:), Nmu_new(:)
  double precision, allocatable :: densities_new(:,:), rho_new(:)

  parameter(Nretries = 10)

! Get regeneration parameters
  call GetRegenerationData(DirNameLength, FromDir, Nprocs_new, rho_reduce)

! Check for files to find Nprocs_old
  call GetNprocsOld(FromDir(1:DirNameLength), Nprocs_old)


! Get old and new input data
  call GetGeneralInput(FromDir(1:DirNameLength), &
       dt_old, Niter_old, dump_period_fields_old, fields_per_file_old, &
       dump_period_distr_old, dump_period_distr_1v_old, dump_start_old, &
       shift_test_period_old,resistance_old, Nz_old, zmin_old, zmax_old, &
       Nspecies_old, const_a_old, BC_Poisson_old, voltage_old, &
       voltage_init_old, initialiser_old, E0_old, startfromdumpfile_old, &
       dump_period_dump_old, exitafterdump_old, transffilename_old, &
       voltagefilename_old, dump_period_distr_IonoBoundary_old)
  Nxi_old = Nz_old
  call Bfield(2,(/ zmin_old, zmax_old /), BBC_old, dBBC_old)

  call GetGeneralInput('', dt_new, Niter_new, dump_period_fields_new, &
       fields_per_file_new, &
       dump_period_distr_new, dump_period_distr_1v_new, dump_start_new, &
       shift_test_period_new, resistance_new, Nz_new, zmin_new, zmax_new, &
       Nspecies_new, const_a_new, BC_Poisson_new, voltage_new, &
       voltage_init_new, initialiser_new, E0_new, startfromdumpfile_new, &
       dump_period_dump_new, exitafterdump_new, transffilename_new, &
       voltagefilename_new, dump_period_distr_IonoBoundary_new)
  Nxi_new = Nz_new
  call Bfield(2,(/ zmin_new, zmax_new /), BBC_new, dBBC_new)


! Get reinitialisation parameters
  allocate(rspecies(Nspecies_old))
  allocate(reinit(Nspecies_old))
  call GetReinitialisationData(Nspecies_old, rspecies, reinit)

  allocate(istart_old(Nprocs_old))
  allocate(iend_old(Nprocs_old))
  do jj = 0, Nprocs_old-1
     call DivideRegion(Nz_old, Nprocs_old, jj, &
          istart_old(jj+1), iend_old(jj+1) )
  end do
  allocate(istart_new(Nprocs_new))
  allocate(iend_new(Nprocs_new))
  do jj = 0, Nprocs_new-1
     call DivideRegion(Nz_new, Nprocs_new, jj, &
          istart_new(jj+1), iend_new(jj+1) )
  end do

 
  allocate(XI_old(Nxi_old))
  allocate(XIcorner_old(Nxi_old+1))
  allocate(gp_old(Nxi_old))
  allocate(gpcorner_old(Nxi_old+1))
  dxi_old = 1.0d0/dble(Nxi_old)
  do ii = 1, Nxi_old
     XI_old(ii) = dxi_old*(0.5d0+dble(ii-1))
     XIcorner_old(ii) = dxi_old*dble(ii-1)
  end do
  XIcorner_old(Nxi_old+1) = 1.0d0
  allocate(Z_old(Nz_old))
  allocate(Zcorner_old(Nz_old+1))
  allocate(dz_old(Nz_old-1))
! Get z=g(xi) and gp=dg/dxi.
  call GetTransform(FromDir(1:DirNameLength), transffilename_old, &
       Nxi_old, XI_old, zmin_old, zmax_old, Z_old, gp_old)
  call GetTransform(FromDir(1:DirNameLength), transffilename_old, &
       Nxi_old+1, XIcorner_old, zmin_old, zmax_old, Zcorner_old, gpcorner_old)
! This definition of dz is only intended for interpolation purposes.
  dz_old = Z_old(2:Nz_old) - Z_old(1:Nz_old-1)

! Allocate and initialise the first old local Z-vector
  Nz_local_old = iend_old(1) - istart_old(1) + 1
  allocate(Z_local_old(-1:Nz_local_old+2))
  allocate(dz_local_old(-1:Nz_local_old+1))
  allocate(B_local_old(-1:Nz_local_old+2))
  allocate(dB_local_old(-1:Nz_local_old+2))
  if ( .not. (istart_old(1) .eq. 1) ) then
     write (*,*) 'Sanity check failed!'
     write (*,*) 'istart_old(1)=', istart_old(1), ' (should be 1)'
     stop
  end if
  Z_local_old(1:Nz_local_old)     = Z_old(istart_old(1):iend_old(1))
  Z_local_old(-1)                 = Z_local_old(1)            - 2.0d0*dz_old(1)
  Zlow_old = Z_local_old(-1)
  Z_local_old(0)                  = Z_local_old(1)            -       dz_old(1)
  dz_local_old(1:Nz_local_old-1)  = dz_old(istart_old(1):iend_old(1)-1)
  dz_local_old(-1)                = dz_old(1)
  dz_local_old(0)                 = dz_old(1)
  if ( iend_old(1)+1 <= Nz_old ) then
     Z_local_old(Nz_local_old+1)  = Z_old(iend_old(1)+1)
     dz_local_old(Nz_local_old) = dz_old(iend_old(1))
  else
     Z_local_old(Nz_local_old+1)  = Z_local_old(Nz_local_old) +dz_old(Nz_old-1)
     dz_local_old(Nz_local_old) = dz_old(Nz_old-1)
  end if
  if ( iend_old(1)+2 <= Nz_old ) then
     Z_local_old(Nz_local_old+2)  = Z_old(iend_old(1)+2)
     dz_local_old(Nz_local_old+1) = dz_old(iend_old(1)+1)
  else
     Z_local_old(Nz_local_old+2)  =Z_local_old(Nz_local_old+1)+dz_old(Nz_old-1)
     dz_local_old(Nz_local_old+1) = dz_old(Nz_old-1)
  end if

  allocate(XI_new(Nxi_new))
  allocate(XIcorner_new(Nxi_new+1))
  allocate(gp_new(Nxi_new))
  allocate(gpcorner_new(Nxi_new+1))
  dxi_new = 1.0d0/dble(Nxi_new)
  do ii = 1, Nxi_new
     XI_new(ii) = dxi_new*(0.5d0+dble(ii-1))
     XIcorner_new(ii) = dxi_new*dble(ii-1)
  end do
  XIcorner_new(Nxi_new+1) = 1.0d0
  allocate(Z_new(Nz_new))
  allocate(Zcorner_new(Nz_new+1))
  allocate(dz_new(Nz_new-1))
! Get z=g(xi) and gp=dg/dxi.
  call GetTransform('', transffilename_new, Nxi_new, XI_new, &
       zmin_new, zmax_new, Z_new, gp_new)
  call GetTransform('', transffilename_new, Nxi_new+1, XIcorner_new, &
       zmin_new, zmax_new, Zcorner_new, gpcorner_new)
! This definition of dz is only intended for interpolation purposes.
  dz_new = Z_new(2:Nz_new) - Z_new(1:Nz_new-1)

  allocate(vector_new(Nz_new))
  Nz_local_new = iend_new(1) - istart_new(1) + 1

  allocate(vector_old(Nz_old))
  allocate(kTzprofile_global(Nz_old,Nspecies_old))
  allocate(kTpprofile_global(Nz_old,Nspecies_old))
  allocate(vzprofile_global(Nz_old,Nspecies_old))
  vector_old        = 0.0d0
  kTzprofile_global = 0.0d0
  kTpprofile_global = 0.0d0
  vzprofile_global  = 0.0d0


  ! Specific parameters for each species
  allocate(particle_old(Nspecies_old))
  allocate(ghost_old(Nspecies_old))
  call GetSpecies(FromDir(1:DirNameLength), Nspecies_old, particle_old)

  allocate(particle_new(Nspecies_new))
  call GetSpecies('', Nspecies_new, particle_new)
  do ii = 1, Nspecies_new
     kTzprofile_global(:,ii) = particle_new(ii)%kTz
     kTpprofile_global(:,ii) = particle_new(ii)%kTp
     vzprofile_global(:,ii)  = particle_new(ii)%vz0
  end do

! Dirty fix to let us use particle_new for both the z and vz interpolations.
! And mu.
  allocate(Nvz_new(Nspecies_new))
  Nvz_new = particle_new%Nvz
  particle_new%Nvz = particle_old%Nvz
  allocate(Nmu_new(Nspecies_new))
  Nmu_new = particle_new%Nmu
  particle_new%Nmu = particle_old%Nmu
! end Dirty fix


  ! initialise species specific constants
  particle_old%dvz = (particle_old%vzmax-particle_old%vzmin)/particle_old%Nvz
  particle_new%dvz = particle_old%dvz
  do ii = 1, Nspecies_old
     ! make allocations of distribution quantities only if mass is finite
     if (.not. isnan(particle_old(ii)%mass)) then
  ! allocate and calculate the Vz-vector for each species
        allocate( particle_old(ii)%Vz(particle_old(ii)%Nvz) )
        do ivz = 1, particle_old(ii)%Nvz
           particle_old(ii)%Vz(ivz) = particle_old(ii)%vz0 + &
                particle_old(ii)%vzmin + &
                particle_old(ii)%dvz*(0.5d0+dble(ivz-1))
        end do
  ! allocate and calculate the mu-vector for each species
        allocate( particle_old(ii)%mu(particle_old(ii)%Nmu) )
        allocate( particle_old(ii)%dmu(particle_old(ii)%Nmu) )
        allocate( particle_old(ii)%maxf0(particle_old(ii)%Nmu) )
        particle_old(ii)%maxf0=0.0d0
        do imu = 1, particle_old(ii)%Nmu
           particle_old(ii)%mu(imu) = particle_old(ii)%mumin + &
                0.5d0*((dble(imu)**particle_old(ii)%muexp + &
                dble(imu-1)**particle_old(ii)%muexp) / &
                particle_old(ii)%Nmu**particle_old(ii)%muexp) * &
                (particle_old(ii)%mumax - particle_old(ii)%mumin)
           particle_old(ii)%dmu(imu) = ((dble(imu)**particle_old(ii)%muexp - &
                dble(imu-1)**particle_old(ii)%muexp) / &
                particle_old(ii)%Nmu**particle_old(ii)%muexp) * &
                (particle_old(ii)%mumax - particle_old(ii)%mumin)
        end do
        allocate( particle_old(ii)%node(-1:Nz_local_old+2) )
        do iz = -1, Nz_local_old+2
           allocate(particle_old(ii)%node(iz)%f(particle_old(ii)%Nvz, &
                particle_old(ii)%Nmu))
        end do
        allocate( ghost_old(ii)%node(2) )
        allocate( ghost_old(ii)%node(1)%f(particle_old(ii)%Nvz, &
             particle_old(ii)%Nmu) )
        ghost_old(ii)%node(1)%f = 0.0d0
        ghost_old(ii)%node(1)%ivzoffset = 0
        allocate( ghost_old(ii)%node(2)%f(particle_old(ii)%Nvz, &
             particle_old(ii)%Nmu) )
        ghost_old(ii)%node(2)%f = 0.0d0
        ghost_old(ii)%node(2)%ivzoffset = 0
     end if
  end do

  if (Nspecies_new .ne. Nspecies_old) then
     write (*,*) 'WARNING: You have changed the number of species.'
     write(*,*) 'This could lead to trouble without end, and it probably will.'
     write (*,*)
  end if


  do ii = 1, Nspecies_new
     ! make allocations of distribution quantities only if mass is finite
     if (.not. isnan(particle_new(ii)%mass)) then
  ! allocate and calculate the Vz-vector for each species
  ! We use the old values of Nvz, Nmu etc., because interpolation is 
  ! only in z now.
        allocate( particle_new(ii)%Vz(particle_old(ii)%Nvz) )
        do ivz = 1, particle_old(ii)%Nvz
           particle_new(ii)%Vz(ivz) = particle_old(ii)%vz0 + &
                particle_old(ii)%vzmin + &
                particle_old(ii)%dvz*(0.5d0+dble(ivz-1))
        end do
  ! allocate and calculate the mu-vector for each species
        allocate( particle_new(ii)%mu(particle_old(ii)%Nmu) )
        allocate( particle_new(ii)%dmu(particle_old(ii)%Nmu) )
        allocate( particle_new(ii)%maxf0(particle_old(ii)%Nmu) )
        particle_new(ii)%maxf0=0.0d0
        do imu = 1, particle_old(ii)%Nmu
           particle_new(ii)%mu(imu) = particle_old(ii)%mumin + &
                0.5d0*((dble(imu)**particle_old(ii)%muexp + &
                dble(imu-1)**particle_old(ii)%muexp) / &
                particle_old(ii)%Nmu**particle_old(ii)%muexp) * &
                (particle_old(ii)%mumax - particle_old(ii)%mumin)
           particle_new(ii)%dmu(imu) = ((dble(imu)**particle_old(ii)%muexp - &
                dble(imu-1)**particle_old(ii)%muexp) / &
                particle_old(ii)%Nmu**particle_old(ii)%muexp) * &
                (particle_old(ii)%mumax - particle_old(ii)%mumin)
        end do
     end if
  end do



  write (*,*) '********************************'
  write (*,*) '***   ketchup regeneration   ***'
  write (*,*) '********************************'
  write (*,*) 'Nprocs_old =', Nprocs_old
  write (*,*) 'Nz_old =', Nz_old
  write (*,*)
  write (*,*) 'Nprocs_new =', Nprocs_new
  write (*,*) 'Nz_new =', Nz_new
  write (*,*) '********************************'
  write (*,*) 'Now it starts...'
  write (*,*) 'Doing the various profiles'


  do ii = 1, Nspecies_old
     call GetSpeciesDensityProfile(FromDir(1:DirNameLength), Nz_old, &
          ii, vector_old, fileexists)
     call RegenerateVector(Nz_old, Nz_new, Z_old, &
          Z_new, vector_old, vector_new)
     if (fileexists) then
        call WriteSpeciesDensityProfile('', Nz_new, ii, vector_new)
     end if

     call InterpolateTemperatureProfile(FromDir(1:DirNameLength), '',  &
          Nz_old, Nz_new, Z_old, Z_new, ii, Nspecies_old, &
          kTzprofile_global, kTpprofile_global, particle_new)
     
     call InterpolateVelocityProfile(FromDir(1:DirNameLength), '',  &
          Nz_old, Nz_new, Z_old, Z_new, ii, Nspecies_old, &
          vzprofile_global, particle_new)
  end do
  call InterpolateVoltageProfile(FromDir(1:DirNameLength), '',  &
       Nz_old, Nz_new, Z_old, Z_new)


  write (*,*) 'Now it really starts...'


! Check if we need to interpolate in vz. We need this later
  dovz = .false.
  domu = .false.
  do ii = 1, Nspecies_new
     if ( (particle_old(ii)%Nvz /= Nvz_new(ii)) .or. &
          (particle_old(ii)%vzmin /= particle_new(ii)%vzmin) .or. &
          (particle_old(ii)%vzmax /= particle_new(ii)%vzmax) ) then
        dovz = .true.
     end if
     if ( (particle_old(ii)%Nmu /= Nmu_new(ii)) .or. &
          (particle_old(ii)%mumin /= particle_new(ii)%mumin) .or. &
          (particle_old(ii)%mumax /= particle_new(ii)%mumax) ) then
        domu = .true.
     end if
  end do


! First we interpolate in z, then we take on the velocity coordinates


! Load the first of the old files
  ll=0
  attempts = 0
  call LoadDump(FromDir(1:DirNameLength),Nspecies_old, &
       Nz_local_old, iter_start, particle_old, ll, attempts, Nretries)
  call Bfield(Nz_local_old+4,Z_local_old,B_local_old,dB_local_old)
  allocate(nprofile_local(-1:Nz_local_old+2,Nspecies_old))
  allocate(kTzprofile_local(-1:Nz_local_old+2,Nspecies_old))
  allocate(kTpprofile_local(-1:Nz_local_old+2,Nspecies_old))
  allocate(vzprofile_local(-1:Nz_local_old+2,Nspecies_old))
  if (ll==0) then 
     dstart = 0
  else
     dstart = -2
  end if
  if (ll == Nprocs_old-1) then
     dend = 0
  else
     dend = 2
  end if
  call Compute_nprofile(Nz_local_old,Nspecies_old,particle_old,nprofile_local)
  kTzprofile_local = 0.0d0
  kTzprofile_local(1+dstart:Nz_local_old+dend,:) = &
       kTzprofile_global(istart_old(ll+1)+dstart:iend_old(ll+1)+dend,:)
  kTpprofile_local = 0.0d0
  kTpprofile_local(1+dstart:Nz_local_old+dend,:) = &
       kTpprofile_global(istart_old(ll+1)+dstart:iend_old(ll+1)+dend,:)
  vzprofile_local = 0.0d0
  vzprofile_local(1+dstart:Nz_local_old+dend,:) = &
       vzprofile_global(istart_old(ll+1)+dstart:iend_old(ll+1)+dend,:)
  call Reinitialise_f(Nz_local_old, Z_local_old, B_local_old, &
       nprofile_local, kTzprofile_local, kTpprofile_local, &
       vzprofile_local, Nspecies_old, particle_old, &
       ll-1, mod(ll+2, Nprocs_old+1)-1, reinit )

  call ParticleToGhost(Nspecies_old, Nz_local_old, particle_old, ghost_old)

! Loop through the new processes
  do kk = 0, Nprocs_new-1
     write (*,*) 'New file', kk+1, 'out of', Nprocs_new
     ! Allocate and construct a local Z-vector for the present process
     Nz_local_new = iend_new(kk+1) - istart_new(kk+1) + 1
     allocate(Z_local_new(Nz_local_new))
     Z_local_new = Z_new(istart_new(kk+1):iend_new(kk+1))
     allocate(meanVz_new(Nz_local_new,Nspecies_new))
     allocate(densities_new(Nz_local_new,Nspecies_new))
     allocate(rho_new(Nz_local_new))
     allocate(B_local_new(Nz_local_new))
     allocate(dB_local_new(Nz_local_new))
     meanVz_new = 0.0d0
     densities_new = 0.0d0
     rho_new = 0.0d0
     call Bfield(Nz_local_new,Z_local_new,B_local_new,dB_local_new)

     ! Allocate the new distribution function for the present process
     do ii = 1, Nspecies_new
        if (.not. isnan(particle_new(ii)%mass)) then
           allocate( particle_new(ii)%node(-1:Nz_local_new+2) )
           do iz = -1, Nz_local_new+2
              allocate(particle_new(ii)%node(iz)%f(particle_old(ii)%Nvz, &
                   particle_old(ii)%Nmu))
              particle_new(ii)%node(iz)%f = 0.0d0
           end do
        end if
     end do

     ! Go through all the grid points in z for the new process
     do iz = -1, Nz_local_new+2
        if (iz==-1) then
           if (kk==0) then
              this_new_Z = Z_local_new(1) - 2.0d0*dz_new(1)
           else
              this_new_Z = Z_new(istart_new(kk+1)-2)
           end if
        elseif (iz==0) then
           if (kk==0) then
              this_new_Z = Z_local_new(1) - dz_new(1)
           else
              this_new_Z = Z_new(istart_new(kk+1)-1)
           end if
        elseif (iz==Nz_local_new+1) then
           if (kk==Nprocs_new-1) then
              this_new_Z = Z_local_new(Nz_local_new) + dz_new(Nz_new-1)
           else
              this_new_Z = Z_new(iend_new(kk+1)+1)
           end if
        elseif (iz==Nz_local_new+2) then 
           if (kk==Nprocs_new-1) then
              this_new_Z = Z_local_new(Nz_local_new) + 2.0d0*dz_new(Nz_new-1)
           else
              this_new_Z = Z_new(iend_new(kk+1)+2)
           end if
        else
           this_new_Z = Z_local_new(iz)
        end if
       ! Adjust to avoid errors at the boundaries
        if (this_new_Z<zmin_old) then
           this_new_Z = zmin_old - 1.5d0*dz_old(1)
        end if
        if (this_new_Z>zmax_old) then
           this_new_Z = zmax_old + 1.5d0*dz_old(Nz_old-1)
        end if

        ! Are we in the range of the already loaded old file?
        if ( (this_new_Z>Z_local_old(-1) .and. &
             this_new_Z<=Z_local_old(Nz_local_old)) .or. &
             (ll==0 .and. this_new_Z<=Z_local_old(Nz_local_old)) .or. &
             (ll==Nprocs_old-1 .and. this_new_Z>=Z_local_old(1)) ) then
           ! Yes, go right ahead
        else
           ! No, we need to load the next file
           do kkk = 1, Nprocs_old
              if (istart_old(kkk)<3) then
                 Zcmp_old = Zlow_old
              else
                 Zcmp_old =  Z_old(istart_old(kkk)-2)
              end if
              if ( this_new_Z > Zcmp_old .and. &
                   this_new_Z <= Z_old(iend_old(kkk)) ) then
                 exit
              end if
           end do
           if (kkk>Nprocs_old) then
              kkk = Nprocs_old
           end if
           ll = kkk-1
           ! Now that we know the file number, we need to 
           ! reallocate the old memory, update the old Z-vector 
           ! and load the file.
           deallocate( Z_local_old )
           deallocate( dz_local_old )
           deallocate( B_local_old )
           deallocate( dB_local_old )
           deallocate( nprofile_local )
           deallocate( kTzprofile_local )
           deallocate( kTpprofile_local )
           deallocate( vzprofile_local )
           Nz_local_old = iend_old(ll+1) - istart_old(ll+1) + 1
           allocate(Z_local_old(-1:Nz_local_old+2))
           allocate(dz_local_old(-1:Nz_local_old+1))
           allocate(B_local_old(-1:Nz_local_old+2))
           allocate(dB_local_old(-1:Nz_local_old+2))
           Z_local_old(1:Nz_local_old) = Z_old(istart_old(ll+1):iend_old(ll+1))
           dz_local_old(1:Nz_local_old-1) = &
                dz_old(istart_old(ll+1):iend_old(ll+1)-1)
           if ( istart_old(ll+1)-2 >= 1 ) then
              Z_local_old(-1)              = Z_old(istart_old(ll+1)-2) 
              dz_local_old(-1)             = dz_old(istart_old(ll+1)-2) 
           else
              Z_local_old(-1)              = Z_local_old(1)   - 2.0d0*dz_old(1)
              dz_local_old(-1)             = dz_old(1)
           end if
           if ( istart_old(ll+1)-1 >= 1 ) then
              Z_local_old(0)               = Z_old(istart_old(ll+1)-1) 
              dz_local_old(0)              = dz_old(istart_old(ll+1)-1) 
           else
              Z_local_old(0)               = Z_local_old(1)   -       dz_old(1)
              dz_local_old(0)              = dz_old(1)
           end if
           if ( iend_old(ll+1)+1 <= Nz_old ) then
              Z_local_old(Nz_local_old+1)  = Z_old(iend_old(ll+1)+1)
              dz_local_old(Nz_local_old) = dz_old(iend_old(ll+1))
           else
              Z_local_old(Nz_local_old+1)  = Z_local_old(Nz_local_old) + &
                   dz_old(Nz_old-1)
              dz_local_old(Nz_local_old) = dz_old(Nz_old-1)
           end if
           if ( iend_old(ll+1)+2 <= Nz_old ) then
              Z_local_old(Nz_local_old+2)  = Z_old(iend_old(ll+1)+2)
              dz_local_old(Nz_local_old+1) = dz_old(iend_old(ll+1)+1)
           else
              Z_local_old(Nz_local_old+2)  = Z_local_old(Nz_local_old+1) + &
                   dz_old(Nz_old-1)
              dz_local_old(Nz_local_old+1) = dz_old(Nz_old-1)
           end if

           do ii = 1, Nspecies_old
              if (.not. isnan(particle_old(ii)%mass)) then
                 deallocate( particle_old(ii)%node )
              end if
           end do
           do ii = 1, Nspecies_old
              if (.not. isnan(particle_old(ii)%mass)) then
                 allocate( particle_old(ii)%node(-1:Nz_local_old+2) )
                 do jj = -1, Nz_local_old+2
                    allocate(particle_old(ii)%node(jj)% &
                         f(particle_old(ii)%Nvz,particle_old(ii)%Nmu))
                 end do
              end if
           end do
           attempts=0
           call LoadDump(FromDir(1:DirNameLength),Nspecies_old, &
                Nz_local_old, iter_start, particle_old, ll, attempts, Nretries)
!!$           write (*,*) 'Loaded file no ', ll
           call Bfield(Nz_local_old+4,Z_local_old,B_local_old,dB_local_old)
           allocate(nprofile_local(-1:Nz_local_old+2,Nspecies_old))
           allocate(kTzprofile_local(-1:Nz_local_old+2,Nspecies_old))
           allocate(kTpprofile_local(-1:Nz_local_old+2,Nspecies_old))
           allocate(vzprofile_local(-1:Nz_local_old+2,Nspecies_old))
           if (ll==0) then 
              dstart = 0
           else
              dstart = -2
           end if
           if (ll == Nprocs_old-1) then
              dend = 0
           else
              dend = 2
           end if
           call Compute_nprofile(Nz_local_old,Nspecies_old, &
                particle_old,nprofile_local)
           kTzprofile_local = 0.0d0
           kTzprofile_local(1+dstart:Nz_local_old+dend,:) = &
                kTzprofile_global(istart_old(ll+1)+dstart:iend_old(ll+1)+ &
                dend,:)
           kTpprofile_local = 0.0d0
           kTpprofile_local(1+dstart:Nz_local_old+dend,:) = &
               kTpprofile_global(istart_old(ll+1)+dstart:iend_old(ll+1)+dend,:)
           vzprofile_local = 0.0d0
           vzprofile_local(1+dstart:Nz_local_old+dend,:) = &
                vzprofile_global(istart_old(ll+1)+dstart:iend_old(ll+1)+dend,:)
           call Reinitialise_f(Nz_local_old, Z_local_old, B_local_old, &
                nprofile_local, kTzprofile_local, kTpprofile_local, &
                vzprofile_local, Nspecies_old, particle_old, &
                ll-1, mod(ll+2, Nprocs_old+1)-1, reinit )

           call GhostToParticle(Nspecies_old, Nz_local_old, &
                particle_old, ghost_old)
           call ParticleToGhost(Nspecies_old, Nz_local_old, &
                particle_old, ghost_old)
        end if

        ! Who are our nearest neighbours, And what are their weights?
        ! (borrowing the old mustard variable name Ng)
        do mm = -1, Nz_local_old+1
           Ng = mm
           if ( (this_new_Z <= Z_local_old(Ng+1)) ) then
              exit
           end if
        end do
     
        weight = ( this_new_Z - Z_local_old(Ng) )/dz_local_old(Ng)
        ! Adjust out of range points. This should only happen if the 
        ! new Z-range is wider than the old one.
        if ( weight < 0.0d0 ) then
           weight = 0.0d0
        end if
        if ( weight > 1.0d0 ) then
           weight = 1.0d0
        end if

        ! Interpolate to find the new value!
        do ii = 1, Nspecies_new
           if (.not. isnan(particle_new(ii)%mass)) then
              ! We interpolate to find a  new parallel velocity offset
              particle_new(ii)%node(iz)%ivzoffset = &
                   nint(particle_old(ii)%node(Ng)%ivzoffset*(1.0d0-weight) + &
                   particle_old(ii)%node(Ng+1)%ivzoffset * weight)

              ! First the contribution from node Ng+1
              ! difference between offsets: doff
              doff = particle_old(ii)%node(Ng+1)%ivzoffset - &
                   particle_new(ii)%node(iz)%ivzoffset
              if (doff == 0) then
                 particle_new(ii)%node(iz)%f = &
                      particle_old(ii)%node(Ng+1)%f * weight
              elseif (doff<0) then
                 particle_new(ii)%node(iz)%f(1,:) = &
                      particle_old(ii)%node(Ng+1)%f(1-doff,:) * weight
                 do jj = 1, -doff 
                    particle_new(ii)%node(iz)%f(1,:) = &
                         particle_new(ii)%node(iz)%f(1,:) + &
                         particle_old(ii)%node(Ng+1)%f(1-doff-jj,:) * weight
                 end do
                 do jj = 2, particle_old(ii)%Nvz + doff
                    particle_new(ii)%node(iz)%f(jj,:) = &
                         particle_old(ii)%node(Ng+1)%f(jj-doff,:) * weight
                 end do
                 do jj = particle_old(ii)%Nvz + doff, particle_old(ii)%Nvz  
                    particle_new(ii)%node(iz)%f(jj,:) = 0.0d0
                 end do
              else
                 do jj = 1, doff
                    particle_new(ii)%node(iz)%f(jj,:) = 0.0d0
                 end do
                 do jj = doff+1, particle_old(ii)%Nvz
                    particle_new(ii)%node(iz)%f(jj,:) = &
                         particle_old(ii)%node(Ng+1)%f(jj-doff,:) * weight
                 end do
                 do jj = particle_old(ii)%Nvz, particle_old(ii)%Nvz + doff
                    particle_new(ii)%node(iz)%f(particle_old(ii)%Nvz,:) = &
                         particle_new(ii)%node(iz)%f(particle_old(ii)%Nvz,:)+ &
                         particle_old(ii)%node(Ng+1)%f(jj-doff,:) * weight 
                 end do
              end if

              ! Then the contribution from node Ng
              ! difference between offsets: doff
              doff = particle_old(ii)%node(Ng)%ivzoffset - &
                   particle_new(ii)%node(iz)%ivzoffset
              if (doff == 0) then
                 particle_new(ii)%node(iz)%f = particle_new(ii)%node(iz)%f + &
                      particle_old(ii)%node(Ng)%f * (1.0d0-weight)
              elseif (doff<0) then
                 particle_new(ii)%node(iz)%f(1,:) = &
                      particle_new(ii)%node(iz)%f(1,:) + &
                      particle_old(ii)%node(Ng)%f(1-doff,:) * (1.0d0-weight)
                 do jj = 1, -doff 
                    particle_new(ii)%node(iz)%f(1,:) = &
                         particle_new(ii)%node(iz)%f(1,:) + &
                         particle_old(ii)%node(Ng)%f(1-doff-jj,:) * &
                         (1.0d0-weight)
                 end do
                 do jj = 2, particle_old(ii)%Nvz + doff
                    particle_new(ii)%node(iz)%f(jj,:) = &
                         particle_new(ii)%node(iz)%f(jj,:) + &
                         particle_old(ii)%node(Ng)%f(jj-doff,:)*(1.0d0-weight)
                 end do
!                 do jj = particle_old(ii)%Nvz + doff, particle_old(ii)%Nvz  
!                    particle_new(ii)%node(iz)%f(jj,:) = &
!                         particle_new(ii)%node(iz)%f(jj,:) + 0.0d0
!                 end do
              else
!                 do jj = 1, doff
!                    particle_new(ii)%node(iz)%f(jj,:) = &
!                         particle_new(ii)%node(iz)%f(jj,:) + 0.0d0
!                 end do
                 do jj = doff+1, particle_old(ii)%Nvz
                    particle_new(ii)%node(iz)%f(jj,:) = &
                         particle_new(ii)%node(iz)%f(jj,:) + &
                         particle_old(ii)%node(Ng)%f(jj-doff,:)*(1.0d0-weight)
                 end do
                 do jj = particle_old(ii)%Nvz, particle_old(ii)%Nvz + doff
                    particle_new(ii)%node(iz)%f(particle_old(ii)%Nvz,:) =  &
                         particle_new(ii)%node(iz)%f(particle_old(ii)%Nvz,:)+ &
                         particle_old(ii)%node(Ng)%f(jj-doff,:)*(1.0d0-weight) 
                 end do
              end if

           end if
        end do
     end do

     ! If the z-range is not cut, generate a new boundary distribution 
     ! instead of the interpolated one, but if we are doing vz, this 
     ! will be done later instead.
     if ( kk==0 .and. zmin_new<=zmin_old .and. .not. dovz) then
        ! This requires a relativistic gamma for relativistic species
        do ii = 1, Nspecies_new
           if ( particle_new(ii)%relativistic ) then
              ! allocate and compute relativistic gamma variable
              allocate( particle_new(ii)%node(0)%gamma &
                   (particle_new(ii)%Nvz,particle_new(ii)%Nmu) )
              particle_new(ii)%node(0)%ivzoffset = &
                   nint((particle_new(ii)%vz0L-particle_new(ii)%vz0)/ &
                   particle_new(ii)%dvz)
              call Einstein(particle_new(ii)%Nvz, particle_new(ii)%Nmu, &
                   dble(particle_new(ii)%node(0)%ivzoffset) * &
                   particle_new(ii)%dvz+particle_new(ii)%Vz, &
                   particle_new(ii)%mu, BBC_new(1), &
                   particle_new(ii)%mass, &
                   particle_new(ii)%node(0)%gamma)
           end if
        end do
        call InitialiseBoundaryLeft(Nz_local_new, BBC_new, &
             Nspecies_new, particle_new)
     end if
     if ( kk==Nprocs_new-1 .and. zmax_new>=zmax_old .and. .not. dovz ) then
        ! This requires a relativistic gamma for relativistic species
        do ii = 1, Nspecies_new
           if ( particle_new(ii)%relativistic ) then
              allocate( particle_new(ii)%node(Nz_local_new+1)%gamma &
                   (particle_new(ii)%Nvz,particle_new(ii)%Nmu) )
              particle_new(ii)%node(Nz_local_new+1)%ivzoffset = &
                   nint((particle_new(ii)%vz0R-particle_new(ii)%vz0)/ &
                   particle_new(ii)%dvz)
              call Einstein(particle_new(ii)%Nvz, &
                   particle_new(ii)%Nmu, dble( &
                   particle_new(ii)%node(Nz_local_new+1)%ivzoffset)* &
                   particle_new(ii)%dvz+particle_new(ii)%Vz, &
                   particle_new(ii)%mu, BBC_new(2), &
                   particle_new(ii)%mass, &
                   particle_new(ii)%node(Nz_local_new+1)%gamma)
           end if
        end do
        call InitialiseBoundaryRight(Nz_local_new, BBC_new, &
             Nspecies_new, particle_new)
        ! If the left hand boundary moved, the regeneration will 
        ! create a mismatch at the right hand boundary, for which 
        ! we compensate here. It is also advisable to use the 
        ! renormalise program to get correct densities
        do ii = 1, Nspecies_new
           if (.not. isnan(particle_new(ii)%mass)) then  
              particle_new(ii)%node(Nz_local_new+1)%f = &
                   particle_new(ii)%node(Nz_local_new+1)%f * &
                   BBC_old(1) / BBC_new(1)
              particle_new(ii)%node(Nz_local_new+2)%f = &
                   particle_new(ii)%node(Nz_local_new+2)%f * &
                   BBC_old(1) / BBC_new(1)
           end if
        end do
     end if


     ! Simply copy the maxf0 vector
     do ii = 1, Nspecies_new
        if (.not. isnan(particle_new(ii)%mass)) then
           particle_new(ii)%maxf0 = particle_old(ii)%maxf0
        end if
     end do


     ! Turn shifting off, if that has been requested. The opposite operation 
     ! of turning it on shall require no action.
     do ii = 1, Nspecies_new
        if (  particle_new(ii)%vzshifting == .false. .and. &
             particle_old(ii)%vzshifting == .true. ) then
           call ShiftOff(Nz_local_new,Nspecies_new,ii,particle_new)
        end if
     end do

     ! Perform a velocity shift test, and shift if appropriate
     call ComputeMeanVz(Nz_local_new,Nspecies_new,particle_new,meanVz_new)
     call ComputeDensity(Nz_local_new,Nspecies_new,particle_new, &
          B_local_new,BBC_new,densities_new,rho_new)
     call TestAndShiftNonRelativistic(Nz_local_new,Nspecies_new, &
          particle_new,meanVz_new,densities_new)

     ! Reduce rho if requested
     if ( rho_reduce /= 1.0d0 )then
        call ReduceRho(Nz_local_new, Nspecies_new, particle_new, &
             densities_new, rho_reduce)
     end if

! Dump and deallocate
     attempts = 0
     call DumpDump(Nspecies_new, Nz_local_new, 0, particle_new, kk, &
          attempts, Nretries)
     deallocate( Z_local_new )
     deallocate( meanVz_new )
     deallocate( densities_new )
     deallocate( rho_new )
     deallocate( B_local_new )
     deallocate( dB_local_new )

     do ii = 1, Nspecies_new
        if (.not. isnan(particle_new(ii)%mass)) then
           deallocate( particle_new(ii)%node )
        end if
     end do

  end do

! Deallocate the old in preparation for what is to come, or not.
  deallocate( Z_local_old )
  deallocate( dz_local_old )
  deallocate( B_local_old )
  deallocate( dB_local_old )
  deallocate( nprofile_local )
  deallocate( kTzprofile_local )
  deallocate( kTpprofile_local )
  deallocate( vzprofile_local )
  do ii = 1, Nspecies_old
     if (.not. isnan(particle_old(ii)%mass)) then
        deallocate( particle_old(ii)%node )
        deallocate( particle_old(ii)%Vz )
     end if
  end do
  do ii = 1, Nspecies_new
     if (.not. isnan(particle_new(ii)%mass)) then
        deallocate( particle_new(ii)%Vz )
     end if
  end do

  write (*,*) 'The z-direction interpolation is over.'


! Dirty fix - part 2 - reinstates particle_new%Nvz
  particle_new%Nvz=Nvz_new
! end Dirty fix - part 2

  particle_new%dvz=(particle_new%vzmax-particle_new%vzmin)/particle_new%Nvz

  if ( dovz ) then
     write (*,*) 'Doing vz.'

     ! Loop through the new processes
     do kk = 0, Nprocs_new-1
        write (*,*) kk+1, '/', Nprocs_new

        Nz_local_new = iend_new(kk+1) - istart_new(kk+1) + 1
        allocate(meanVz_new(Nz_local_new,Nspecies_new))
        allocate(densities_new(Nz_local_new,Nspecies_new))
        allocate(rho_new(Nz_local_new))
        allocate(B_local_new(Nz_local_new))
        allocate(dB_local_new(Nz_local_new))
        allocate(Z_local_new(Nz_local_new))
        Z_local_new = Z_new(istart_new(kk+1):iend_new(kk+1))
        meanVz_new = 0.0d0
        densities_new = 0.0d0
        rho_new = 0.0d0
        call Bfield(Nz_local_new,Z_local_new,B_local_new,dB_local_new)
        ! Allocate the new distribution function for the present process
        do ii = 1, Nspecies_new
           if (.not. isnan(particle_new(ii)%mass)) then
              allocate( particle_new(ii)%node(-1:Nz_local_new+2) )
              do iz = -1, Nz_local_new+2
                 allocate(particle_new(ii)%node(iz)%f(particle_new(ii)%Nvz, &
                      particle_old(ii)%Nmu))
              end do
              allocate( particle_new(ii)%Vz(particle_new(ii)%Nvz) )
              do ivz = 1, particle_new(ii)%Nvz
                 particle_new(ii)%Vz(ivz) = particle_new(ii)%vz0 + &
                      particle_new(ii)%vzmin + &
                      particle_new(ii)%dvz*(0.5d0+dble(ivz-1))
              end do
           end if
           if ( particle_new(ii)%relativistic ) then
              ! allocate and compute relativistic gamma variable
              if (kk==0) then
                 allocate( particle_new(ii)%node(0)%gamma &
                      (particle_new(ii)%Nvz,particle_new(ii)%Nmu) )
                 particle_new(ii)%node(0)%ivzoffset = &
                      nint((particle_new(ii)%vz0L-particle_new(ii)%vz0)/ &
                      particle_new(ii)%dvz)
                 call Einstein(particle_new(ii)%Nvz, particle_new(ii)%Nmu, &
                      dble(particle_new(ii)%node(0)%ivzoffset) * &
                      particle_new(ii)%dvz+particle_new(ii)%Vz, &
                      particle_new(ii)%mu, BBC_new(1), &
                      particle_new(ii)%mass, &
                      particle_new(ii)%node(0)%gamma)
              end if
              if (kk==Nprocs_new-1) then
                 allocate( particle_new(ii)%node(Nz_local_new+1)%gamma &
                      (particle_new(ii)%Nvz,particle_new(ii)%Nmu) )
                 particle_new(ii)%node(Nz_local_new+1)%ivzoffset = &
                      nint((particle_new(ii)%vz0R-particle_new(ii)%vz0)/ &
                      particle_new(ii)%dvz)
                 call Einstein(particle_new(ii)%Nvz, &
                      particle_new(ii)%Nmu, dble( &
                      particle_new(ii)%node(Nz_local_new+1)%ivzoffset)* &
                      particle_new(ii)%dvz+particle_new(ii)%Vz, &
                      particle_new(ii)%mu, BBC_new(2), &
                      particle_new(ii)%mass, &
                      particle_new(ii)%node(Nz_local_new+1)%gamma)
              end if
           end if
        end do
        ! and the old
        do ii = 1, Nspecies_old
           if (.not. isnan(particle_old(ii)%mass)) then
              allocate( particle_old(ii)%node(-1:Nz_local_new+2) )
              do iz = -1, Nz_local_new+2
                 allocate(particle_old(ii)%node(iz)%f(particle_old(ii)%Nvz, &
                      particle_old(ii)%Nmu))
              end do
              allocate( particle_old(ii)%Vz(particle_old(ii)%Nvz) )
              do ivz = 1, particle_old(ii)%Nvz
                 particle_old(ii)%Vz(ivz) = particle_old(ii)%vz0 + &
                      particle_old(ii)%vzmin + &
                      particle_old(ii)%dvz*(0.5d0+dble(ivz-1))
              end do
           end if
        end do
        ! Load the distribution, which is now considered old.
        attempts=0
        call LoadDump('',Nspecies_old,Nz_local_new,iter_start, &
             particle_old,kk, attempts, Nretries)
        do iz = -1, Nz_local_new+2
           do ii = 1, Nspecies_new
              if (.not. isnan(particle_old(ii)%mass)) then
                 do ivz = 1, particle_new(ii)%Nvz
                    particle_new(ii)%node(iz)%ivzoffset = &
                         nint( particle_old(ii)%node(iz)%ivzoffset * &
                         particle_old(ii)%dvz / particle_new(ii)%dvz )
                    this_new_vz = particle_new(ii)%Vz(ivz) + &
                         particle_new(ii)%node(iz)%ivzoffset * &
                         particle_new(ii)%dvz
                    if ( (this_new_vz<particle_old(ii)%vzmin + &
                         particle_old(ii)%node(iz)%ivzoffset * &
                         particle_old(ii)%dvz) .or. &
                         (this_new_vz>particle_old(ii)%vzmax + &
                         particle_old(ii)%node(iz)%ivzoffset * &
                         particle_old(ii)%dvz) ) then
                       particle_new(ii)%node(iz)%f(ivz,:) = 0.0d0
                    else
                 ! Who are our nearest neighbours, And what are their weights?
                       Ng = int( (this_new_vz - &
                            particle_old(ii)%Vz(1) - &
                            particle_old(ii)%node(iz)%ivzoffset * &
                            particle_old(ii)%dvz ) / &
                            particle_old(ii)%dvz) + 1
                       weight = ( this_new_vz - &
                            particle_old(ii)%Vz(Ng) - &
                            particle_old(ii)%node(iz)%ivzoffset * &
                            particle_old(ii)%dvz )/particle_old(ii)%dvz


                 ! Adjust out of range points. This should only happen if the 
                 ! new Z-range is wider than the old one.
                       if (Ng<1) then
                          write (*,*) 'BUT THAT IS IMPOSSIBLE!'
               write (*,*) 'Ng=', Ng
               write (*,*) 'weight=', weight
               write (*,*) 'this_new_vz=',this_new_vz
               write (*,*) 'particle_old(ii)%dvz=',particle_old(ii)%dvz
               write (*,*) 'particle_old(ii)%Vz(1)=',particle_old(ii)%Vz(1)
               write (*,*) 'particle_old(ii)%node(iz)%ivzoffset=', &
                    particle_old(ii)%node(iz)%ivzoffset
                          Ng = 1
                          weight = 0.0d0
                       end if
                       if (Ng>=particle_old(ii)%Nvz) then
                          particle_new(ii)%node(iz)%f(ivz,:) = &
                               particle_old(ii)%node(iz)% &
                               f(particle_old(ii)%Nvz,:) * (1.0d0-weight)
                       else
                          particle_new(ii)%node(iz)%f(ivz,:) = &
                               particle_old(ii)%node(iz)%f(Ng,:) * &
                               (1.0d0-weight) + &
                               particle_old(ii)%node(iz)%f(Ng+1,:) * weight
                       end if
                    end if
                 end do
              end if
           end do
        end do

        ! If the z-range is not cut, generate a new boundary distribution 
        ! instead of the interpolated one.
        if ( kk==0 .and. zmin_new<=zmin_old ) then
           call InitialiseBoundaryLeft(Nz_local_new, BBC_new, &
                Nspecies_new, particle_new)
        end if
        if ( kk==Nprocs_new-1 .and. zmax_new>=zmax_old ) then
           call InitialiseBoundaryRight(Nz_local_new, BBC_new, &
                Nspecies_new, particle_new)
           ! If the left hand boundary moved, the regeneration will 
           ! create a mismatch at the right hand boundary, for which 
           ! we compensate here. It is also advisable to use the 
           ! renormalise program to get correct densities
           do ii = 1, Nspecies_new
              if (.not. isnan(particle_new(ii)%mass)) then  
                 particle_new(ii)%node(Nz_local_new+1)%f = &
                      particle_new(ii)%node(Nz_local_new+1)%f * &
                      BBC_old(1) / BBC_new(1)
                 particle_new(ii)%node(Nz_local_new+2)%f = &
                      particle_new(ii)%node(Nz_local_new+2)%f * &
                      BBC_old(1) / BBC_new(1)
              end if
           end do
        end if

        ! Simply copy the maxf0 vector
        do ii = 1, Nspecies_new
           if (.not. isnan(particle_new(ii)%mass)) then
              particle_new(ii)%maxf0 = particle_old(ii)%maxf0
           end if
        end do
        ! Perform a velocity shift test, and shift if appropriate
        call ComputeMeanVz(Nz_local_new,Nspecies_new,particle_new,meanVz_new)
        call ComputeDensity(Nz_local_new,Nspecies_new,particle_new, &
             B_local_new,BBC_new,densities_new,rho_new)
        call TestAndShiftNonRelativistic(Nz_local_new,Nspecies_new, &
             particle_new,meanVz_new,densities_new)

! Dump and deallocate
        attempts = 0
        call DumpDump(Nspecies_new, Nz_local_new, 0, particle_new, kk, &
             attempts, Nretries)
        deallocate( meanVz_new )
        deallocate( densities_new )
        deallocate( rho_new )
        deallocate( Z_local_new )
        deallocate( B_local_new )
        deallocate( dB_local_new )
        do ii = 1, Nspecies_new
           if (.not. isnan(particle_new(ii)%mass)) then
              deallocate( particle_new(ii)%node )
              deallocate( particle_new(ii)%Vz )
           end if
        end do
        do ii = 1, Nspecies_old
           if (.not. isnan(particle_old(ii)%mass)) then
              deallocate( particle_old(ii)%node )
              deallocate( particle_old(ii)%Vz )
           end if
        end do

     end do
     write (*,*) 'The vz-direction interpolation is over.'
  end if


! Dirty fix - part 3 - reinstates particle_new%Nmu
  particle_new%Nmu=Nmu_new
! end Dirty fix - part 3


  if ( domu ) then
     write (*,*) 'Doing mu.'

     ! Allocate and compute Vz and mu for all species before we start.
     particle_new%dvz=(particle_new%vzmax-particle_new%vzmin)/particle_new%Nvz
     do ii = 1, Nspecies_new
        if (.not. isnan(particle_new(ii)%mass)) then
           allocate( particle_new(ii)%Vz(particle_new(ii)%Nvz) )
           do ivz = 1, particle_new(ii)%Nvz
              particle_new(ii)%Vz(ivz) = particle_new(ii)%vz0 + &
                   particle_new(ii)%vzmin + &
                   particle_new(ii)%dvz*(0.5d0+dble(ivz-1))
           end do
           deallocate ( particle_new(ii)%mu )
           deallocate ( particle_new(ii)%dmu )
           deallocate ( particle_new(ii)%maxf0 )
           allocate( particle_new(ii)%mu(particle_new(ii)%Nmu) )
           allocate( particle_new(ii)%dmu(particle_new(ii)%Nmu) )
           allocate( particle_new(ii)%maxf0(particle_new(ii)%Nmu) )
           do imu = 1, particle_new(ii)%Nmu
              particle_new(ii)%mu(imu) = particle_new(ii)%mumin + &
                   0.5d0*((dble(imu)**particle_new(ii)%muexp + &
                   dble(imu-1)**particle_new(ii)%muexp) / &
                   particle_new(ii)%Nmu**particle_new(ii)%muexp) * &
                   (particle_new(ii)%mumax - particle_new(ii)%mumin)
              particle_new(ii)%dmu(imu) = &
                   ((dble(imu)**particle_new(ii)%muexp - &
                   dble(imu-1)**particle_new(ii)%muexp) / &
                   particle_new(ii)%Nmu**particle_new(ii)%muexp) * &
                   (particle_new(ii)%mumax - particle_new(ii)%mumin)
           end do
        end if
     end do


     ! Loop through the new processes
     do kk = 0, Nprocs_new-1
        write (*,*) kk+1, '/', Nprocs_new
        Nz_local_new = iend_new(kk+1) - istart_new(kk+1) + 1
        allocate(B_local_new(Nz_local_new))
        allocate(dB_local_new(Nz_local_new))
        allocate(Z_local_new(Nz_local_new))
        Z_local_new = Z_new(istart_new(kk+1):iend_new(kk+1))
        call Bfield(Nz_local_new,Z_local_new,B_local_new,dB_local_new)
        ! Allocate the new distribution function for the present process
        do ii = 1, Nspecies_new
           if (.not. isnan(particle_new(ii)%mass)) then
              allocate( particle_new(ii)%node(-1:Nz_local_new+2) )
              do iz = -1, Nz_local_new+2
                 allocate(particle_new(ii)%node(iz)%f(particle_new(ii)%Nvz, &
                      particle_new(ii)%Nmu))
              end do
           end if
           if ( particle_new(ii)%relativistic ) then
              ! allocate and compute relativistic gamma variable
              if (kk==0) then
                 allocate( particle_new(ii)%node(0)%gamma &
                      (particle_new(ii)%Nvz,particle_new(ii)%Nmu) )
                 particle_new(ii)%node(0)%ivzoffset = &
                      nint((particle_new(ii)%vz0L-particle_new(ii)%vz0)/ &
                      particle_new(ii)%dvz)
                 call Einstein(particle_new(ii)%Nvz, particle_new(ii)%Nmu, &
                      dble(particle_new(ii)%node(0)%ivzoffset) * &
                      particle_new(ii)%dvz+particle_new(ii)%Vz, &
                      particle_new(ii)%mu, BBC_new(1), &
                      particle_new(ii)%mass, &
                      particle_new(ii)%node(0)%gamma)
              end if
              if (kk==Nprocs_new-1) then
                 allocate( particle_new(ii)%node(Nz_local_new+1)%gamma &
                      (particle_new(ii)%Nvz,particle_new(ii)%Nmu) )
                 particle_new(ii)%node(Nz_local_new+1)%ivzoffset = &
                      nint((particle_new(ii)%vz0R-particle_new(ii)%vz0)/ &
                      particle_new(ii)%dvz)
                 call Einstein(particle_new(ii)%Nvz, &
                      particle_new(ii)%Nmu, dble( &
                      particle_new(ii)%node(Nz_local_new+1)%ivzoffset)* &
                      particle_new(ii)%dvz+particle_new(ii)%Vz, &
                      particle_new(ii)%mu, BBC_new(2), &
                      particle_new(ii)%mass, &
                      particle_new(ii)%node(Nz_local_new+1)%gamma)
              end if
           end if
        end do
        ! and the old
        do ii = 1, Nspecies_old
           if (.not. isnan(particle_old(ii)%mass)) then
              allocate( particle_old(ii)%node(-1:Nz_local_new+2) )
              do iz = -1, Nz_local_new+2
                 allocate(particle_old(ii)%node(iz)%f(particle_new(ii)%Nvz, &
                      particle_old(ii)%Nmu))
              end do
           end if
        end do
        ! Load the distribution, which is now considered old.
        attempts=0
        call LoadDump('',Nspecies_old,Nz_local_new,iter_start, &
             particle_old,kk, attempts, Nretries)
        do iz = -1, Nz_local_new+2
          do ii = 1, Nspecies_new
              if (.not. isnan(particle_old(ii)%mass)) then
                 particle_new(ii)%node(iz)%ivzoffset = &
                      particle_old(ii)%node(iz)%ivzoffset
                 do imu = 1, particle_new(ii)%Nmu
                    this_new_mu = particle_new(ii)%mu(imu)
                    if  (this_new_mu<particle_old(ii)%mu(1)) then
                       particle_new(ii)%node(iz)%f(:,imu) = &
                            particle_old(ii)%node(iz)%f(:,1)
                    elseif (this_new_mu > &
                         particle_old(ii)%mu(particle_old(ii)%Nmu)) then
                       particle_new(ii)%node(iz)%f(:,imu) = 0.0d0
                    else
                 ! Who are our nearest neighbours, And what are their weights?
                       do mm = 1, particle_old(ii)%Nmu
                          Ng = mm
                          if ( this_new_mu <= particle_old(ii)%mu(Ng+1) ) then
                             exit
                          end if
                       end do
                       weight=2.0d0*(this_new_mu - particle_old(ii)%mu(Ng))/ &
                            ( particle_old(ii)%dmu(Ng) + &
                            particle_old(ii)%dmu(Ng+1) )
                       ! Check that weight is in [0, 1]
                       if ( weight < 0.0d0 ) then
                          weight = 0.0d0
                       end if
                       if ( weight > 1.0d0 ) then
                          weight = 1.0d0
                       end if

                       if (Ng<1) then
                          write (*,*) 'BUT THAT IS IMPOSSIBLE!'
               write (*,*) 'Ng=', Ng
               write (*,*) 'weight=', weight
               write (*,*) 'this_new_vz=',this_new_vz
               write (*,*) 'particle_old(ii)%dvz=',particle_old(ii)%dvz
               write (*,*) 'particle_old(ii)%Vz(1)=',particle_old(ii)%Vz(1)
               write (*,*) 'particle_old(ii)%node(iz)%ivzoffset=', &
                    particle_old(ii)%node(iz)%ivzoffset
                          Ng = 1
                          weight = 0.0d0
                       end if
                       particle_new(ii)%node(iz)%f(:,imu) = &
                            particle_old(ii)%node(iz)%f(:,Ng) * &
                            (1.0d0-weight) + &
                            particle_old(ii)%node(iz)%f(:,Ng+1) * weight
                    end if
                 end do
              end if
           end do
        end do

        ! If the z-range is not cut, generate a new boundary distribution 
        ! instead of the interpolated one.
        if ( kk==0 .and. zmin_new<=zmin_old ) then
           call InitialiseBoundaryLeft(Nz_local_new, BBC_new, &
                Nspecies_new, particle_new)
        end if
        if ( kk==Nprocs_new-1 .and. zmax_new>=zmax_old ) then
           call InitialiseBoundaryRight(Nz_local_new, BBC_new, &
                Nspecies_new, particle_new)
           ! If the left hand boundary moved, the regeneration will 
           ! create a mismatch at the right hand boundary, for which 
           ! we compensate here. It is also advisable to use the 
           ! renormalise program to get correct densities
           do ii = 1, Nspecies_new
              if (.not. isnan(particle_new(ii)%mass)) then  
                 particle_new(ii)%node(Nz_local_new+1)%f = &
                      particle_new(ii)%node(Nz_local_new+1)%f * &
                      BBC_old(1) / BBC_new(1)
                 particle_new(ii)%node(Nz_local_new+2)%f = &
                      particle_new(ii)%node(Nz_local_new+2)%f * &
                      BBC_old(1) / BBC_new(1)
              end if
           end do
        end if

        ! Set the maxf0 vector to zero. The user must run maxf0update
        ! to get it right again.
        do ii = 1, Nspecies_new
           if (.not. isnan(particle_new(ii)%mass)) then
              particle_new(ii)%maxf0 = 0.0d0
           end if
        end do

! Dump and deallocate
        attempts = 0
        call DumpDump(Nspecies_new, Nz_local_new, 0, particle_new, kk, &
             attempts, Nretries)
        deallocate( Z_local_new )
        deallocate( B_local_new )
        deallocate( dB_local_new )
        do ii = 1, Nspecies_new
           if (.not. isnan(particle_new(ii)%mass)) then
              deallocate( particle_new(ii)%node )
           end if
        end do
        do ii = 1, Nspecies_old
           if (.not. isnan(particle_old(ii)%mass)) then
              deallocate( particle_old(ii)%node )
           end if
        end do

     end do
     write (*,*) 'The mu-direction interpolation is over.'
     write (*,*) 'Run maxf0update to get a correct maxf0 vector.'
  end if


  write (*,*) 'it ends'




end program regenerate_ketchup

!------------------------------------------------------------------------

subroutine GetRegenerationData(DirNameLength, FromDir, Nprocs_new, rho_reduce)

  implicit none

  integer DirNameLength, Nprocs_new
  character FromDir*200
  double precision rho_reduce

  character indata*300
  integer i, j, k, kk

  ! Defaults
  FromDir = '~'
  Nprocs_new = 0
  rho_reduce = 1.0d0
  open(unit=1,file='regen_par.m',status='old')


  read (1,'(a)') indata
  do while (index(indata,'%END')+index(indata,'%end') == 0)
     i=index(indata,'=')
     if (i>0) then
        do j=1,i-1
           if(indata(j:j) /= ' ') exit
        enddo
        ! If there is a semicolon on the line, ignore it and everything beyond.
        k = index(indata(i+1:300),';')
        if (k==0) then
           k = 301-i
        end if
        select case (indata(j:i-1))
        case ('FromDir')
           read (indata(i+1:i+k-1),*) FromDir
        case ('Nprocs_new')
           read (indata(i+1:i+k-1),*) Nprocs_new
        case ('rho_reduce')
           read (indata(i+1:i+k-1),*) rho_reduce
        case default
           if (index(indata(j:j+1),'%')==0) then
              write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
           end if
        end select
     end if
     read (1,'(a)') indata
  end do

  close(1)

  DirNameLength=index(FromDir,'    ')-1
  if (FromDir(DirNameLength:DirNameLength) /= '/') then
     FromDir(DirNameLength+1:DirNameLength+1) = '/'
     DirNameLength = DirNameLength+1
  end if

end subroutine GetRegenerationData

!------------------------------------------------------------------------

subroutine GetReinitialisationData(Nspecies, rspecies, reinit)

  implicit none

  integer Nspecies
  integer rspecies(Nspecies)
  logical reinit(Nspecies)
  character FromDir*200

  character indata*300
  integer i, j, k, ii, ierr

  ! Defaults
  rspecies = 0
  reinit = .false.

  open(unit=1,file='regen_par.m',status='old')


  read (1,'(a)') indata
  do while (index(indata,'%REINIT')+index(indata,'%reinit') == 0)
     read (1,'(a)') indata
  end do
  do while (index(indata,'%END')+index(indata,'%end') == 0)
     i=index(indata,'=')
     if (i>0) then
        do j=1,i-1
           if(indata(j:j) /= ' ') exit
        enddo
        ! If there is a semicolon on the line, ignore it and everything beyond.
        k = index(indata(i+1:300),';')
        if (k==0) then
           k = 301-i
        end if
        select case (indata(j:i-1))
        case ('rspecies')
           read (indata(i+1:i+k-1),*,iostat=ierr) rspecies
        case default
           if (index(indata(j:j+1),'%')==0) then
              write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
           end if
        end select
     end if
     read (1,'(a)') indata
  end do

  close(1)

  do ii = 1, Nspecies
     if ( rspecies(ii)>0 .and. rspecies(ii)<=Nspecies ) then
        reinit(rspecies(ii)) = .true.
     end if
  end do

end subroutine GetReinitialisationData

!------------------------------------------------------------------------

subroutine RegenerateVector(Nz_old, Nz_new, Z_old, Z_new, &
     vector_old, vector_new)

  implicit none

  integer Nz_old, Nz_new
  double precision Z_old(Nz_old), vector_old(Nz_old)
  double precision Z_new(Nz_new), vector_new(Nz_new)

  integer Ng, jj, kk, startsearch
  double precision weight, dz_old(Nz_old-1)

! This definition of dz is only intended for interpolation purposes.
  dz_old = Z_old(2:Nz_old) - Z_old(1:Nz_old-1)

  startsearch = 1

  do jj = 1, Nz_new

   ! Who are our nearest neighbours, And what are their weights?
     do kk = startsearch, Nz_old-1
        Ng = kk
        if ( (Z_new(jj) <= Z_old(Ng+1)) ) then
           startsearch = Ng
           exit
        end if
     end do
     
     weight = ( Z_new(jj) - Z_old(Ng) )/dz_old(Ng)
     ! Adjust out of range points. This should only happen if the 
     ! new Z-range is wider than the old one.
     if ( weight < 0.0d0 ) then
        weight = 0.0d0
     end if
     if ( weight > 1.0d0 ) then
        weight = 1.0d0
     end if
     ! Interpolate to find the new value!
     vector_new(jj) = vector_old(Ng) *(1.0d0-weight) + vector_old(Ng+1)*weight
  end do

end subroutine RegenerateVector

!------------------------------------------------------------------------

subroutine GetSpeciesDensityProfile(FromDir, Nz, species, vector, fileexists)

  implicit none

  character(len=*) :: FromDir
  integer Nz, species
  double precision vector(Nz)
  character filename*230, Nznumber*6, snumber*2
  logical fileexists

  vector = 0.0d0

! Assemble filename
  write(snumber,fmt='(i2.2)') species
  write(Nznumber,fmt='(i6.6)') Nz
  filename=FromDir//'n0'//'s'//snumber//'Nz'//Nznumber//'.inp'

! Don't read the file unless it exists
  inquire (file=filename, exist = fileexists)

  if (fileexists) then
! Read input data
     open(unit=1,file=filename,status='old')
     read(1,*) vector
     close(1)
  end if

end subroutine GetSpeciesDensityProfile

!------------------------------------------------------------------------

subroutine WriteSpeciesDensityProfile(ToDir, Nz, species, vector)

  implicit none

  character(len=*) :: ToDir
  integer Nz, species
  double precision vector(Nz)

  character filename*230, Nznumber*6, snumber*2, form*20
  integer jj

! Assemble filename
  write(snumber,fmt='(i2.2)') species
  write(Nznumber,fmt='(i6.6)') Nz
  filename=ToDir//'n0'//'s'//snumber//'Nz'//Nznumber//'.inp'

  write (form,fmt='(a,i6.6,a)') '(', Nz, '(E16.9E3,a))'


! Write input data
  open(unit=1,file=filename,status='replace')
  write(1,fmt=form) (vector(jj), ' ', jj = 1, Nz)
  close(1)

end subroutine WriteSpeciesDensityProfile

!------------------------------------------------------------------------

subroutine InterpolateTemperatureProfile(FromDir, ToDir,  &
     Nz_old, Nz_new, Z_old, Z_new, ii, Nspecies, &
     kTzprofile, kTpprofile, particle)
  
  use SpecificTypes

  implicit none

  character(len=*) :: FromDir, ToDir
  integer Nz_old, Nz_new, ii, Nspecies
  double precision Z_old(Nz_old), vector_old(Nz_old)
  double precision Z_new(Nz_new), vector_new(Nz_new), &
       kTzprofile(Nz_old,Nspecies), kTpprofile(Nz_old,Nspecies)
  type(species) particle(Nspecies)

  character filename*230, Nznumber*6, snumber*2, form*20
  integer jj
  logical fileexists

! First Tz
  vector_old = 0.0d0
  vector_new = 0.0d0

! Assemble filename
  write(snumber,fmt='(i2.2)') ii
  write(Nznumber,fmt='(i6.6)') Nz_old
  filename=FromDir//'Tzs'//snumber//'Nz'//Nznumber//'.inp'

! Only do something if the file exists
  inquire (file=filename, exist = fileexists)

  if (fileexists) then
! Read input data
     open(unit=1,file=filename,status='old')
     read(1,*) vector_old
     close(1)
     if (particle(ii)%Tfromfile) then
        kTzprofile(:,ii) = vector_old
     else
        kTzprofile(:,ii) = particle(ii)%kTz
     end if
     call RegenerateVector(Nz_old, Nz_new, Z_old, &
          Z_new, vector_old, vector_new)


! Assemble filename
     write(Nznumber,fmt='(i6.6)') Nz_new
     filename=ToDir//'Tzs'//snumber//'Nz'//Nznumber//'.inp'

     write (form,fmt='(a,i6.6,a)') '(', Nz_new, '(E16.9E3,a))'


! Write input data
     open(unit=1,file=filename,status='replace')
     write(1,fmt=form) (vector_new(jj), ' ', jj = 1, Nz_new)
     close(1)
  end if

! And then Tp
  vector_old = 0.0d0
  vector_new = 0.0d0

! Assemble filename
  write(snumber,fmt='(i2.2)') ii
  write(Nznumber,fmt='(i6.6)') Nz_old
  filename=FromDir//'Tps'//snumber//'Nz'//Nznumber//'.inp'

! Only do something if the file exists
  inquire (file=filename, exist = fileexists)

  if (fileexists) then
! Read input data
     open(unit=1,file=filename,status='old')
     read(1,*) vector_old
     close(1)
     if (particle(ii)%Tfromfile) then
        kTpprofile(:,ii) = vector_old
     else
        kTpprofile(:,ii) = particle(ii)%kTp 
     end if
     call RegenerateVector(Nz_old, Nz_new, Z_old, &
          Z_new, vector_old, vector_new)


! Assemble filename
     write(Nznumber,fmt='(i6.6)') Nz_new
     filename=ToDir//'Tps'//snumber//'Nz'//Nznumber//'.inp'

     write (form,fmt='(a,i6.6,a)') '(', Nz_new, '(E16.9E3,a))'


! Write input data
     open(unit=1,file=filename,status='replace')
     write(1,fmt=form) (vector_new(jj), ' ', jj = 1, Nz_new)
     close(1)
  end if

end subroutine InterpolateTemperatureProfile

!------------------------------------------------------------------------

subroutine InterpolateVelocityProfile(FromDir, ToDir,  &
     Nz_old, Nz_new, Z_old, Z_new, ii, Nspecies, vzprofile, particle)
  
  use SpecificTypes

  implicit none

  character(len=*) :: FromDir, ToDir
  integer Nz_old, Nz_new, ii, Nspecies
  double precision Z_old(Nz_old), vector_old(Nz_old)
  double precision Z_new(Nz_new), vector_new(Nz_new), &
       vzprofile(Nz_old,Nspecies)
  type(species) particle(Nspecies)

  character filename*230, Nznumber*6, snumber*2, form*20
  integer jj
  logical fileexists

  vector_old = 0.0d0
  vector_new = 0.0d0

! Assemble filename
  write(snumber,fmt='(i2.2)') ii
  write(Nznumber,fmt='(i6.6)') Nz_old
  filename=FromDir//'vzs'//snumber//'Nz'//Nznumber//'.inp'

! Only do something if the file exists
  inquire (file=filename, exist = fileexists)

  if (fileexists) then
! Read input data
     open(unit=1,file=filename,status='old')
     read(1,*) vector_old
     close(1)
     if (particle(ii)%vfromfile) then
        vzprofile(:,ii) = vector_old
     else
        vzprofile(:,ii) = particle(ii)%vz0
     end if
     call RegenerateVector(Nz_old, Nz_new, Z_old, &
          Z_new, vector_old, vector_new)


! Assemble filename
     write(Nznumber,fmt='(i6.6)') Nz_new
     filename=ToDir//'vzs'//snumber//'Nz'//Nznumber//'.inp'

     write (form,fmt='(a,i6.6,a)') '(', Nz_new, '(E16.9E3,a))'


! Write input data
     open(unit=1,file=filename,status='replace')
     write(1,fmt=form) (vector_new(jj), ' ', jj = 1, Nz_new)
     close(1)
  end if

end subroutine InterpolateVelocityProfile

!------------------------------------------------------------------------

subroutine InterpolateVoltageProfile(FromDir, ToDir, Nz_old, Nz_new, &
     Z_old, Z_new)
  
  implicit none

  character(len=*) :: FromDir, ToDir
  integer Nz_old, Nz_new
  double precision Z_old(Nz_old), vector_old(Nz_old)
  double precision Z_new(Nz_new), vector_new(Nz_new)
  character filename*230, Nznumber*6, form*20
  integer jj
  logical fileexists

  vector_old = 0.0d0
  vector_new = 0.0d0

! Assemble filename
  write(Nznumber,fmt='(i6.6)') Nz_old
  filename=FromDir//'VprofNz'//Nznumber//'.inp'

! Only do something if the file exists
  inquire (file=filename, exist = fileexists)

  if (fileexists) then
! Read input data
     open(unit=1,file=filename,status='old')
     read(1,*) vector_old
     close(1)

     call RegenerateVector(Nz_old, Nz_new, Z_old, &
          Z_new, vector_old, vector_new)


! Assemble filename
     write(Nznumber,fmt='(i6.6)') Nz_new
     filename=ToDir//'VprofNz'//Nznumber//'.inp'

     write (form,fmt='(a,i6.6,a)') '(', Nz_new, '(E16.9E3,a))'


! Write input data
     open(unit=1,file=filename,status='replace')
     write(1,fmt=form) (vector_new(jj), ' ', jj = 1, Nz_new)
     close(1)
  end if

end subroutine InterpolateVoltageProfile

!------------------------------------------------------------------------

subroutine GhostToParticle(Nspecies, Nz_local, particle, ghost)

  use SpecificTypes

  implicit none

  integer Nspecies, Nz_local
  type(species) particle(Nspecies), ghost(Nspecies)

  integer ii

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
     particle(ii)%node(-1)%f = ghost(ii)%node(1)%f
     particle(ii)%node(-1)%ivzoffset = ghost(ii)%node(1)%ivzoffset
     particle(ii)%node(0)%f  = ghost(ii)%node(2)%f
     particle(ii)%node(0)%ivzoffset  = ghost(ii)%node(2)%ivzoffset
     end if
  end do

end subroutine GhostToParticle

!------------------------------------------------------------------------

subroutine ParticleToGhost(Nspecies, Nz_local, particle, ghost)

  use SpecificTypes

  implicit none

  integer Nspecies, Nz_local
  type(species) particle(Nspecies), ghost(Nspecies)

  integer ii

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
     ghost(ii)%node(1)%f = particle(ii)%node(Nz_local-1)%f
     ghost(ii)%node(1)%ivzoffset = particle(ii)%node(Nz_local-1)%ivzoffset
     ghost(ii)%node(2)%f = particle(ii)%node(Nz_local)%f
     ghost(ii)%node(2)%ivzoffset = particle(ii)%node(Nz_local)%ivzoffset
     end if
  end do

end subroutine ParticleToGhost

!------------------------------------------------------------------------

! This initialisation is non-relativistic

subroutine Reinitialise_f(Nz, Z, B, nprofile, &
     kTzprofile, kTpprofile, vzprofile, Nspecies, particle, &
     neighbourLeft, neighbourRight, reinit)

  use SpecificTypes

  implicit none

  integer Nz, Nspecies
  double precision Z(-1:Nz+2), B(-1:Nz+2), nprofile(-1:Nz+2,Nspecies), &
       kTzprofile(-1:Nz+2,Nspecies), kTpprofile(-1:Nz+2,Nspecies), &
       vzprofile(-1:Nz+2,Nspecies)
  type(species) particle(Nspecies)
  integer neighbourLeft, neighbourRight
  logical reinit(Nspecies)

  integer iz, imu, ivz, ii, jj, izmin, izmax
  double precision elc, vzth(-1:Nz+2,Nspecies), vzthL(Nspecies), &
       vzthR(Nspecies), vpth(-1:Nz+2,Nspecies), vpthL(Nspecies), &
       vpthR(Nspecies), vzz, vz1, d1, d2, dz, g

  parameter (elc=1.602176487d-19)

  dz = Z(2)-Z(1)
  do ii = 1, Nspecies
     vzth(:,ii) = sqrt(elc*kTzprofile(:,ii)/particle(ii)%mass)
     vpth(:,ii) = sqrt(elc*kTpprofile(:,ii)/particle(ii)%mass)
  end do
  vzthL = sqrt(elc*particle%kTzL/particle%mass)
  vzthR = sqrt(elc*particle%kTzR/particle%mass)
  vpthL = sqrt(elc*particle%kTpL/particle%mass)
  vpthR = sqrt(elc*particle%kTpR/particle%mass)


  if (neighbourLeft == -1) then
     izmin = 1
  else
     izmin = -1
  end if
  if (neighbourRight == -1) then
     izmax = Nz
  else
     izmax = Nz + 2
  end if

  do ii = 1, Nspecies
     if ( (.not. isnan(particle(ii)%mass)) .and. reinit(ii) ) then  
        ! The interior
        do iz = 1, Nz
           ! No initial offset
           particle(ii)%node(iz)%ivzoffset = 0
           do imu = 1, particle(ii)%Nmu
              d1 = particle(ii)%mu(imu)*B(iz) / &
                   (particle(ii)%mass*vpth(iz,ii)**2.0d0)
              g = (1.0d0/(2.0d0*particle(ii)%mass*vpth(iz,ii)))*exp(-d1)
              do ivz = 1, particle(ii)%Nvz
                 vzz = particle(ii)%Vz(ivz)
                 vz1 = ( vzz - vzprofile(iz,ii) ) / vzth(iz,ii)
                 d2 = 0.5d0*(vz1*vz1)
                 particle(ii)%node(iz)%f(ivz,imu) = nprofile(iz,ii)*g*exp(-d2)
              end do
           end do
        end do
     end if
  end do

! normalise the distribution function and find its maximum local value
  do ii = 1, Nspecies
     if ( (.not. isnan(particle(ii)%mass)) .and. reinit(ii) ) then  
        do iz = 1, Nz
           particle(ii)%node(iz)%f = particle(ii)%node(iz)%f * &
                nprofile(iz,ii)/ &
                (sum(matmul(particle(ii)%node(iz)%f,particle(ii)%dmu)) * &
                particle(ii)%dvz)
        end do
     end if
  end do

end subroutine Reinitialise_f

!------------------------------------------------------------------------

subroutine Compute_nprofile(Nz, Nspecies, particle, nprofile)

  use SpecificTypes
  implicit none

  integer Nz, Nspecies
  double precision nprofile(-1:Nz+2,Nspecies)
  type(species) particle(Nspecies)

  integer ii, iz

  do ii = 1, Nspecies
     if (isnan(particle(ii)%mass)) then
! In the infinite mass case nprofile is not used, and we set it to zero.
        nprofile(iz,ii) = 0.0d0
     else
! In the case of finite mass, nprofile is normalised density.
        do iz = -1, Nz+2
           nprofile(iz,ii) = &
                sum(matmul(particle(ii)%node(iz)%f,particle(ii)%dmu)) * &
                particle(ii)%dvz
        end do
     end if
  end do


end subroutine Compute_nprofile

!------------------------------------------------------------------------
subroutine ShiftOff(Nxi,Nspecies,ii,particle)

  use SpecificTypes
  implicit none

  integer Nxi, Nspecies
  type(species) particle(Nspecies)
  double precision, allocatable :: newf(:,:)

  integer ii, ixi, ivz, imu, thisshift, neighbouroffsets(2)

  if ( .not. isnan(particle(ii)%mass)) then
     allocate(newf(particle(ii)%Nvz,particle(ii)%Nmu))
     do ixi = 1, Nxi
        ! If the current offset is not zero we shift so that it becomes zero.
        if (particle(ii)%node(ixi)%ivzoffset /= 0) then
           thisshift = -particle(ii)%node(ixi)%ivzoffset
           ! Ok, let us shift
           newf = eoshift(particle(ii)%node(ixi)%f,thisshift,0.0d0,1)
           if (thisshift>0) then
              do imu = 1, particle(ii)%Nmu
                 newf(1,imu) = newf(1,imu) + &
                      sum(particle(ii)%node(ixi)%f(1:thisshift,imu))
              end do
           end if
           if (thisshift<0) then
              do imu = 1, particle(ii)%Nmu
                 newf(particle(ii)%Nvz,imu)=newf(particle(ii)%Nvz,imu)+ &
                      sum(particle(ii)%node(ixi)%f(particle(ii)%Nvz + &
                      thisshift+1:particle(ii)%Nvz,imu))
              end do
           end if
           particle(ii)%node(ixi)%f = newf
           particle(ii)%node(ixi)%ivzoffset = 0
        end if
     end do
     deallocate(newf)
  end if

end subroutine ShiftOff

!------------------------------------------------------------------------
subroutine TestAndShiftNonRelativistic(Nxi,Nspecies,particle,meanVz,densities)

  use SpecificTypes
  implicit none

  integer Nxi, Nspecies
  double precision meanVz(Nxi,Nspecies), densities(Nxi,Nspecies), &
       vzdiff, neighbourdens(2)
  type(species) particle(Nspecies)
  double precision, allocatable :: newf(:,:)

  integer ii, ixi, ivz, imu, thisshift, neighbouroffsets(2)
  do ii = 1, Nspecies
     if ( particle(ii)%vzshifting .and. .not. isnan(particle(ii)%mass)) then
        allocate(newf(particle(ii)%Nvz,particle(ii)%Nmu))
        do ixi = 1, Nxi
           ! To shift or not to shift; that is the question.
           !     Find the shifts of the neighbour cells.
           neighbouroffsets(1)=particle(ii)%node(ixi-1)%ivzoffset
           neighbouroffsets(2)=particle(ii)%node(ixi+1)%ivzoffset
           !     Find the difference between present offset and the mean.
           vzdiff = meanVz(ixi,ii) - (particle(ii)%vz0 + &
                particle(ii)%node(ixi)%ivzoffset * particle(ii)%dvz)
           !     If we shift, we shift by thisshift.
           thisshift = nint(vzdiff/particle(ii)%dvz)
           if ( abs(vzdiff) >= 0.75d0*particle(ii)%dvz .or. &
                (maxval(abs(neighbouroffsets - &
                (particle(ii)%node(ixi)%ivzoffset+thisshift))) > &
                particle(ii)%Nvz/4) ) then
              ! If the cell to cell offset difference is greater than Nvz/4, 
              ! then thisshift will be a by density weighted average.
              if (maxval(abs(neighbouroffsets - &
                   (particle(ii)%node(ixi)%ivzoffset+thisshift))) > &
                   particle(ii)%Nvz/4) then
                 ! We can't get density information from the other processes
                 ! without message passing, so we hope this will suffice.
                 if (ixi==1) then
                    neighbourdens(1)=densities(1,ii)
                 else
                    neighbourdens(1)=densities(ixi-1,ii)
                 end if
                 if (ixi==Nxi) then
                    neighbourdens(2)=densities(Nxi,ii)
                 else
                    neighbourdens(2)=densities(ixi+1,ii)
                 end if
                 if ( (sum(neighbourdens)+densities(ixi,ii))>0 ) then
                    thisshift =nint( ( neighbouroffsets(1)*neighbourdens(1) + &
                         (particle(ii)%node(ixi)%ivzoffset + thisshift)* &
                         densities(ixi,ii) + &
                         neighbouroffsets(2)*neighbourdens(2) ) / &
                         (sum(neighbourdens)+densities(ixi,ii)) ) -  &
                         particle(ii)%node(ixi)%ivzoffset
                    else
                       thisshift = 0
                 end if
              end if
              ! Ok, let us shift
              newf = eoshift(particle(ii)%node(ixi)%f,thisshift,0.0d0,1)
              if (thisshift>0) then
                 do imu = 1, particle(ii)%Nmu
                    newf(1,imu) = newf(1,imu) + &
                         sum(particle(ii)%node(ixi)%f(1:thisshift,imu))
                 end do
              end if
              if (thisshift<0) then
                 do imu = 1, particle(ii)%Nmu
                    newf(particle(ii)%Nvz,imu)=newf(particle(ii)%Nvz,imu)+ &
                         sum(particle(ii)%node(ixi)%f(particle(ii)%Nvz + &
                         thisshift+1:particle(ii)%Nvz,imu))
                 end do
              end if
              particle(ii)%node(ixi)%f = newf
              particle(ii)%node(ixi)%ivzoffset = &
                   particle(ii)%node(ixi)%ivzoffset + thisshift
           end if
        end do
        deallocate(newf)
     end if
  end do

end subroutine TestAndShiftNonRelativistic

!------------------------------------------------------------------------

!!$subroutine InitialiseBoundaryLeft(Nxi, BBC, B, Nspecies, particle)
subroutine InitialiseBoundaryLeft(Nxi, BBC, Nspecies, particle)
 
  use SpecificTypes

  implicit none

  integer Nxi, Nspecies
  double precision BBC(2)           !!!, B(Nxi)
  type(species) particle(Nspecies)

  integer ixi, imu, ivz, ii, jj
  double precision elc, vzthL(Nspecies), vpthL(Nspecies), vzz, vz1, d1, d2, g

  parameter (elc=1.602176487d-19)
  vzthL = sqrt(elc*particle%kTzL/particle%mass)
  vpthL = sqrt(elc*particle%kTpL/particle%mass)

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then  
        particle(ii)%node(0)%ivzoffset = &
             nint((particle(ii)%vz0L-particle(ii)%vz0)/particle(ii)%dvz)
        do imu = 1, particle(ii)%Nmu
           do ivz = 1, particle(ii)%Nvz
              if ( particle(ii)%relativistic ) then
                 d1 = particle(ii)%mu(imu)*BBC(1) / &
                      (particle(ii)%node(0)%gamma(ivz,imu) * &
                      particle(ii)%mass*vpthL(ii)**2.0d0)
              else
                 d1=particle(ii)%mu(imu)*BBC(1)/ &
                      (particle(ii)%mass*vpthL(ii)**2.0d0)
              end if
              g = (1.0d0/(2.0d0*particle(ii)%mass*vpthL(ii)))*exp(-d1)
              vzz = particle(ii)%Vz(ivz) + &
                   dble(particle(ii)%node(0)%ivzoffset)*particle(ii)%dvz
              vz1 = (vzz-particle(ii)%vz0L)/vzthL(ii)
              d2 = 0.5d0*(vz1*vz1)
              particle(ii)%node(0)%f(ivz,imu) = g*exp(-d2)
              ! If this is a loss cone distribution, remove what 
              ! is in the loss cones.
              if (particle(ii)%lossconeL) then
                 if (BBC(2) >= BBC(1)) then
                    if (abs(vzz) > sqrt(2.0d0*(BBC(2)-BBC(1))* &
                         particle(ii)%mu(imu)/particle(ii)%mass) ) then
                       particle(ii)%node(0)%f(ivz,imu) = 0.0d0
                    end if
                 end if
              end if
           end do
        end do
        particle(ii)%node(0)%f = particle(ii)%node(0)%f * &
             particle(ii)%n0L/ &
             (sum(matmul(particle(ii)%node(0)%f,particle(ii)%dmu)) * &
             particle(ii)%dvz*particle(ii)%n0)
        particle(ii)%node(-1) = particle(ii)%node(0)
     end if
  end do

end subroutine InitialiseBoundaryLeft

!------------------------------------------------------------------------

subroutine InitialiseBoundaryRight(Nxi, BBC, Nspecies, particle)
 
  use SpecificTypes

  implicit none

  integer Nxi, Nspecies
  double precision BBC(2), B(Nxi)
  type(species) particle(Nspecies)

  integer ixi, imu, ivz, ii, jj
  double precision elc, vzthR(Nspecies), vpthR(Nspecies), vzz, vz1, d1, d2, g

  parameter (elc=1.602176487d-19)
  
  vzthR = sqrt(elc*particle%kTzR/particle%mass)
  vpthR = sqrt(elc*particle%kTpR/particle%mass)

     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then  
           particle(ii)%node(Nxi+1)%ivzoffset = &
                nint((particle(ii)%vz0R-particle(ii)%vz0)/particle(ii)%dvz)
           do imu = 1, particle(ii)%Nmu
              do ivz = 1, particle(ii)%Nvz        
                 if ( particle(ii)%relativistic ) then
                    d1 = particle(ii)%mu(imu)*BBC(2) / &
                         (particle(ii)%node(Nxi+1)%gamma(ivz,imu) * &
                         particle(ii)%mass*vpthR(ii)**2.0d0)
                 else
                    d1 = particle(ii)%mu(imu)*BBC(2) / &
                         (particle(ii)%mass*vpthR(ii)**2.0d0)
                 end if
                 g = (1.0d0/(2.0d0*particle(ii)%mass*vpthR(ii)))*exp(-d1)
                 vzz = particle(ii)%Vz(ivz) + &
                      particle(ii)%node(Nxi+1)%ivzoffset*particle(ii)%dvz
                 vz1 = (vzz-particle(ii)%vz0R)/vzthR(ii)
                 d2 = 0.5d0*(vz1*vz1)
                 particle(ii)%node(Nxi+1)%f(ivz,imu) = g*exp(-d2)
       ! If this is a loss cone distribution, remove what is in the loss cones.
                 if (particle(ii)%lossconeR) then
                    if (BBC(1) >= BBC(2)) then
                       if (abs(vzz) > sqrt(2.0d0*(BBC(1)-BBC(2))* &
                            particle(ii)%mu(imu)/particle(ii)%mass) ) then
                          particle(ii)%node(0)%f(ivz,imu) = 0.0d0
                       end if
                    end if
                 end if
              end do
           end do
           particle(ii)%node(Nxi+1)%f = particle(ii)%node(Nxi+1)%f * &
                particle(ii)%n0R*(BBC(1)/BBC(2))/ &
                (sum(matmul(particle(ii)%node(Nxi+1)%f,particle(ii)%dmu)) * &
                particle(ii)%dvz*particle(ii)%n0)
           particle(ii)%node(Nxi+2) = particle(ii)%node(Nxi+1)
        end if
     end do

end subroutine InitialiseBoundaryRight

!------------------------------------------------------------------------

subroutine ReduceRho(Nxi, Nspecies, particle, densities1D, rho_reduce)

  use SpecificTypes
  implicit none

  integer Nxi, Nspecies
  double precision densities1D(Nxi,Nspecies), rho_reduce
  type(species) particle(Nspecies)

  double precision maxn, adjustment, gdens(Nspecies), rho(Nxi), grho
  integer ii, kk, ixi, iimaxn, ghosts(4)

! Compute a rho vector that is internal to this function. It ignores 
! the B-factor since we are only scaling distribution functions relatively
! for each z independent of the other zs.
  rho = matmul(densities1D,particle%charge)


! Find then densest negative species, compute the adjustment factor, and 
! adjust the distribution function accordingly. 
! This shall reduce rho so that rho_new = rho_old / rho_reduce.
  do ixi = 1, Nxi
     maxn = 0.0d0
     do ii = 1, Nspecies
        if ( .not. isnan(particle(ii)%mass) ) then
           if ( particle(ii)%charge < 0.0d0 ) then
              if (densities1D(ixi,ii) > maxn )then
                 maxn = densities1D(ixi,ii)
                 iimaxn = ii
              end if
           end if
        end if
     end do
     adjustment = 1.0d0 - (1.0d0-1.0d0/rho_reduce)*rho(ixi) / &
          (densities1D(ixi,iimaxn)*particle(iimaxn)%charge)
     particle(iimaxn)%node(ixi)%f = particle(iimaxn)%node(ixi)%f * adjustment
  end do

! Now, do it for the ghost points too.
  ghosts = (/ -1, 0, Nxi+1, Nxi+2 /)
  do kk = 1, 4
     ixi= ghosts(kk)
     maxn = 0.0d0
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then
           gdens(ii)=sum(matmul(particle(ii)%node(ixi)%f,particle(ii)%dmu))* &
                particle(ii)%dvz * particle(ii)%n0
        end if
     end do
     grho = dot_product(gdens,particle%charge) ! ignoring the B factor
     do ii = 1, Nspecies
        if ( .not. isnan(particle(ii)%mass) ) then
           if ( particle(ii)%charge < 0.0d0 ) then
              if (gdens(ii) > maxn )then
                 maxn = gdens(ii)
                 iimaxn = ii
              end if
           end if
        end if
     end do
     adjustment = 1.0d0 - (1.0d0-1.0d0/rho_reduce)*grho / &
          (gdens(iimaxn)*particle(iimaxn)%charge)
     particle(iimaxn)%node(ixi)%f = particle(iimaxn)%node(ixi)%f * adjustment
  end do

end subroutine ReduceRho

!------------------------------------------------------------------------
