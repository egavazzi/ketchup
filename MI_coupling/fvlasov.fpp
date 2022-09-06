! This is ketchup version b6.
!
! This fortran code is included as online supplementary material with the 
! article ''Vlasov simulations of parallel potential drops'' by 
! Herbert Gunell, Johan De Keyser, Emmanuel Gamby, and Ingrid Mann.
!
! Multiple particle species can be specified in the input file. 
! Initial temperatures, drift velocities, and densities for these 
! distributions are specified in the input file. Initial drift velocities, 
! densities, and temperatures can be read from files if one so wishes. 
! Binary files containing the distribution function are dumped at regular 
! intervals. It is possible to start the simulation anew after 
! reading these files. It is possible to specify a time dependent total 
! voltage in a file.
!
!
! dt*dv<dz is required.
!

!------------------------------------------------------------------------

program ketchup

  use SpecificTypes
  use mpi

  implicit none


! input data
  integer Niter, dump_period_fields, fields_per_file, dump_period_distr, &
       dump_period_distr_1v, dump_start, &
       shift_test_period, Nxi, Nz, Nspecies, BC_Poisson, initialiser, &
       dump_period_dump
  double precision zmin, zmax, dt, resistance, voltage, voltage_init, E0
  character transffilename*200, voltagefilename*200
  logical startfromdumpfile, exitafterdump

! Voltage control arrays
  integer NVTable
  integer, allocatable :: VTableIter(:)
  double precision, allocatable :: VTableVoltage(:)


! grid parameters
  double precision dxi, BBC(2), dBBC(2)
  double precision, allocatable :: Z(:), Zcorner(:), dz(:), &
       XI(:), XIcorner(:), gp(:), gpcorner(:), gp_local(:), &
       E(:), B(:), dB(:), gravity(:), densities1D(:,:), rho(:), &
       meanVz(:,:), current(:), &
       nprofile(:,:), kTzprofile(:,:), kTpprofile(:,:), vzprofile(:,:), &
       Vinitial(:)
  ! local grid
  integer Nz_local, Nxi_local
  double precision, allocatable ::  Z_local(:), dz_local(:), &
       XI_local(:), XIcorner_local(:), Zcorner_local(:)

! Buffers which gather field variables that are to be dumped
  integer, allocatable :: iterationbuffer(:)
  double precision, allocatable :: Ebuffer(:,:), currentbuffer(:,:), &
       density1Dbuffer(:,:,:)
  integer fields_collected

! particle species
  type(species), allocatable :: particle(:)

! relative dielectric constant etc.
  double precision epsr, epsilon, const_a, omegap, omegap2, IntB

! mpi-related variables
  integer myid, Nprocs, comm1d, ierr, neighbourLeft, neighbourRight
  integer istart, iend, dims(1)
  logical periods(1)

! other variables
  integer ii, jj, iteration, iter_start, Nretries, attempts

! series resistor 
  double precision resistI, resistV

#ifdef _TIMING_
  double precision starttime, endtime, totalstarttime, totalendtime
#endif

! fundamental physical constants
  double precision eps0, elc
  parameter(eps0 = 8.85418781762039d-12)
  parameter(elc  = 1.602176487d-19)

  parameter(Nretries = 10)


! ionospheric boundary
  double precision Boundary_I_file(301,200,100), Boundary_I_file2(301,200,100), &
          Boundary_I_file_2D(60200,100), Boundary_I_file2_2D(60200,100), &
					Boundary_I(200,100), Boundary_I2(200,100)
  integer i1,index,imu,ivz


! Initialise mpi
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierr)
 

  if (myid==0) then
     write (*,*) '*************************************'
     write (*,*) '****           ketchup           ****'
     write (*,*) '*************************************'
	end if

!  Read the incoming flux and allocate it
  open(unit=1,file='f_response.bin',access="stream",form="unformatted")
	open(unit=2,file='f_response2.bin',access="stream",form="unformatted")
	do imu = 1,100
		do ivz=1,60200
			read(1) Boundary_I_file_2D(ivz,imu)
			read(2) Boundary_I_file2_2D(ivz,imu)
		end do
	end do
  close(1)
  close(2)

  do i1 = 1,301
	do imu = 1,100
		do ivz = 1,200
               index = (i1 - 1) * 200 + ivz
               Boundary_I_file(i1,ivz,imu) = Boundary_I_file_2D(index,imu)
               Boundary_I_file2(i1,ivz,imu) = Boundary_I_file2_2D(index,imu)
          end do
     end do
  end do


! Read input data, allocate, and initialise
  ! General parameters
  call GetGeneralInput('', dt, Niter, dump_period_fields, &
       fields_per_file, dump_period_distr, dump_period_distr_1v, dump_start, & 
       shift_test_period, resistance, Nz, zmin, zmax, Nspecies, const_a, &
       BC_Poisson, voltage, voltage_init, initialiser, E0, &
       startfromdumpfile, dump_period_dump, exitafterdump, transffilename, &
       voltagefilename)
  Nxi = Nz

#ifdef _DEBUG_
!!$  if (myid == 0) then
  write (*,*) 'Done with GetGeneralInput.  myid=', myid
!!$  end if
#endif

  ! Process 0 gets the voltage control table
  if (myid == 0) then
     call GetNVTable('', voltagefilename, NVTable)
     if (NVTable > 0) then
        allocate ( VTableIter(NVTable) )
        allocate ( VTableVoltage(NVTable) )
        call GetVTable('', voltagefilename, NVTable, VTableIter, VTableVoltage)
#ifdef _DEBUG_
        write (*,*) 'Done with GetVTable.  myid=', myid
#endif
     end if
  end if


! For simplicity we keep a xi-vector for the whole system
  allocate(XI(Nxi))
  allocate(XIcorner(Nxi+1))
  allocate(gp(Nxi))
  allocate(gpcorner(Nxi+1))
  dxi = 1.0d0/dble(Nxi)
  do ii = 1, Nxi
     XI(ii) = dxi*(0.5d0+dble(ii-1))
     XIcorner(ii) = dxi*dble(ii-1)
  end do
  XIcorner(Nxi+1) = 1.0d0
  allocate(Z(Nz))
  allocate(Zcorner(Nz+1))
  allocate(dz(Nz))
! Get z=g(xi) and gp=dg/dxi.
  call GetTransform('',transffilename, Nxi, XI, zmin, zmax, Z, gp)
  call GetTransform('',transffilename,Nxi+1,XIcorner,zmin,zmax,Zcorner,gpcorner)
  dz = Zcorner(2:Nxi+1) - Zcorner(1:Nxi)
! Process 0 saves the transformation result.
  if (myid == 0) then
     attempts = 0
     call Dumpg(Nxi,Z,gp,attempts,Nretries)
  end if

! Create the local grid. The division of the grid between 
! the different processes is provided by mpi.
  dims(1) = Nprocs
  periods(1) = .false.
  call MPI_CART_CREATE(MPI_COMM_WORLD, 1, dims, periods , .false.,comm1d,ierr)
  call MPI_CART_SHIFT(comm1d, 0, 1, neighbourLeft, neighbourRight, ierr)
  call DivideRegion(Nxi, Nprocs, myid, istart, iend)

  Nxi_local = iend - istart + 1
  Nz_local = Nxi_local
  allocate(XI_local(Nxi_local))
  allocate(XIcorner_local(Nxi_local+1))
  allocate(Z_local(Nz_local))
  allocate(Zcorner_local(Nz_local+1))
  allocate(dz_local(Nz_local))
  allocate(gp_local(Nxi_local))
  XI_local = XI(istart:iend)
  XIcorner_local = XIcorner(istart:iend+1)
  Z_local = Z(istart:iend)
  Zcorner_local = Zcorner(istart:iend+1)
  dz_local = dz(istart:iend)
  gp_local = gp(istart:iend)

  allocate(E(Nxi_local))
  allocate(nprofile(Nxi_local,Nspecies))
  allocate(vzprofile(Nxi_local,Nspecies))
  allocate(kTzprofile(Nxi_local,Nspecies))
  allocate(kTpprofile(Nxi_local,Nspecies))
  allocate(B(Nxi_local))
  allocate(dB(Nxi_local))
  allocate(gravity(Nxi_local))
  allocate(densities1D(Nxi_local,Nspecies))
  allocate(meanVz(Nxi_local,Nspecies))
  allocate(current(Nxi_local))
  allocate(rho(Nxi_local))
  allocate(Vinitial(Nxi_local))

  allocate(iterationbuffer(fields_per_file ))
  allocate(Ebuffer(Nxi_local,fields_per_file ))
  allocate(currentbuffer(Nxi_local,fields_per_file ))
  allocate(density1Dbuffer(Nxi_local,Nspecies,fields_per_file ))

! initialise what we just allocated and other variables
  nprofile = 0.0d0        ! Density profile
  vzprofile = 0.0d0       ! Parallel velocity profile
  kTzprofile = 0.0d0      ! Parallel temperature profile
  kTpprofile = 0.0d0      ! Perpendicular temperature profile
  B = 0.0d0               ! Magnetic flux density
  dB = 0.0d0              ! dB/dz
  gravity = 0.0d0         ! Gravitational acceleration
  densities1D = 0.0d0     ! Density
  meanVz = 0.0d0          ! Mean parallel velocity
  current = 0.0d0         ! Current
  rho = 0.0d0             ! Charge density
  iterationbuffer = -1    ! Buffer for the iteration numbers
  Ebuffer = 0.0d0         ! Buffer for the E-field
  currentbuffer = 0.0d0   ! Buffer for the current
  density1Dbuffer = 0.0d0 ! Buffer for the densities
  fields_collected = 0    ! Counter of iterations saved
  resistI = 0.0d0         ! Current in the series resistor
  resistV = 0.0d0         ! Voltage over the series resistor

! Initialise Vinitial and E in the way the input file specifies.
! Only the default initialiser is allowed in version b2.
  select case (initialiser)
  case default
     E = 0.0d0
     Vinitial = 0.0d0
  end select


! Get gravity information
  call GetGravity(Nz_local,Z_local,gravity)
  attempts = 0
  call DumpGravity(Nz_local,gravity,myid,attempts,Nretries)


#ifdef _DEBUG_
  write (*,*) myid, neighbourLeft, neighbourRight, istart, iend
#endif
  ! Specific parameters for each species
  allocate(particle(Nspecies))
  call GetSpecies('', Nspecies, particle)
  
#ifdef _DEBUG_
  if (myid == 0) then
     write (*,*) 'Done with GetSpecies.  myid=', myid
  end if
#endif


  ! initialise species specific constants
  particle%dvz = (particle%vzmax-particle%vzmin)/particle%Nvz
  do ii = 1, Nspecies
     ! make allocations of distribution quantities only if mass is finite
     if (.not. isnan(particle(ii)%mass)) then
  ! allocate and calculate the Vz-vector for each species
        allocate( particle(ii)%Vz(particle(ii)%Nvz) )
        do jj = 1, particle(ii)%Nvz
           particle(ii)%Vz(jj) = particle(ii)%vz0 + particle(ii)%vzmin + &
                particle(ii)%dvz*(0.5d0+dble(jj-1))
        end do
  ! allocate and calculate the mu-vector for each species
        allocate( particle(ii)%mu(particle(ii)%Nmu) )
        allocate( particle(ii)%dmu(particle(ii)%Nmu) )
        allocate( particle(ii)%maxf0(particle(ii)%Nmu) )
        particle(ii)%maxf0=0.0d0
        do jj = 1, particle(ii)%Nmu
           particle(ii)%mu(jj) = particle(ii)%mumin + &
                0.5d0*((dble(jj)**particle(ii)%muexp + &
                dble(jj-1)**particle(ii)%muexp) / &
                particle(ii)%Nmu**particle(ii)%muexp) * &
                (particle(ii)%mumax - particle(ii)%mumin)
           particle(ii)%dmu(jj) = ((dble(jj)**particle(ii)%muexp - &
                dble(jj-1)**particle(ii)%muexp) / &
                particle(ii)%Nmu**particle(ii)%muexp) * &
                (particle(ii)%mumax - particle(ii)%mumin)
        end do
  ! allocate the distribution function for each species
! The regular grid, where the distribution function is manipulated 
! by this process, runs from 1 to Nxi_local. On either side there are 
! two ghost nodes, the values of which we get from neighbouring 
! processes. If myid==0 or myid==Nprocs-1 these will be locked to 
! the boundary values.
        allocate( particle(ii)%node(-1:Nxi_local+2) )
        do jj = -1, Nxi_local+2
           allocate(particle(ii)%node(jj)%f(particle(ii)%Nvz,particle(ii)%Nmu))
        end do
        if ( particle(ii)%relativistic ) then
           do jj = 1, Nxi_local
              ! allocate relativistic gamma variable
              allocate( particle(ii)%node(jj)%gamma &
                   (particle(ii)%Nvz,particle(ii)%Nmu) )
           end do
        end if
     end if
  end do


! Set epsilon.
! For large time steps the dielectric constant is adjusted to ensure 
! stability. This decreases the plasma frequency and increases the 
! Debye length, and thus enables a coarser grid, while high frequency 
! and small spatial scale phenomena are lost. 
  omegap2 = 0.0d0
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        omegap2 = omegap2 + elc**2 * particle(ii)%n0 /(particle(ii)%mass*eps0)
     end if
  end do
  omegap = sqrt(omegap2)
  epsr = max((omegap*dt*const_a)**2,1.0d0)
  epsilon = eps0 * epsr
  
  if (myid == 0) then
     write (*,*) 'epsilon_r =', epsr     
     write (*,*) 'omegap    =', omegap
     write (*,*) 'dt        =', dt
     write (*,*) '-------------------------------------'
  end if

#ifdef _DEBUG_
!!$  if (myid == 0) then
  write (*,*) 'Done some allocations.  myid=', myid
!!$  end if
#endif

! Get ourselves a B-field
  call Bfield(Nz_local,Z_local,B,dB) ! B-vector for this process
  attempts = 0
  call DumpBfield(Nz_local,B,dB,myid,attempts,Nretries)
  call Bfield(2,(/ zmin, zmax /), BBC, dBBC)
  ! BBC(1) :  B at left hand boundary
  ! BBC(2) :  B at right hand boundary
  IntB = sum(B*dz_local)   ! Local Integral of B for 
                           ! multi-process Poisson solution.
  

#ifdef _DEBUG_
!!$  if (myid == 0) then
  write (*,*) 'Got the b-field.  myid=', myid
!!$  end if
#endif

! Initialise the distributions
  call GetDensityProfile(Nz,Nz_local,Nspecies, &
       nprofile,istart,iend)

#ifdef _DEBUG_
!!$  if (myid == 0) then
  write (*,*) 'Got the profiles.  myid=', myid
!!$  end if
#endif

! If we are starting from a dump file, load it now, otherwise 
! proceed with the initialisation of the distribution function etc.
  if (startfromdumpfile) then
     attempts = 0
     call LoadDump('', Nspecies, Nxi_local, iter_start, particle, myid, &
          attempts, Nretries)
     call ComputeInitialDensity(Nxi_local,Nspecies,particle,B,BBC,nprofile, &
          densities1D,rho)
     attempts = 0
     call DumpDistrBC(Nxi_local,Nspecies,particle,0, &
          neighbourLeft,neighbourRight,attempts,Nretries)

     ! compute gamma with just loaded ivzoffset, 
     ! or with ivzoffset=0 if we are not shifting. 
     do ii = 1, Nspecies
        if ( particle(ii)%relativistic ) then
           if ( particle(ii)%vzshifting ) then
              do jj = 1, Nxi_local
                 call Einstein(particle(ii)%Nvz, particle(ii)%Nmu, &
                      dble(particle(ii)%node(jj)%ivzoffset)*particle(ii)%dvz+ &
                      particle(ii)%Vz, particle(ii)%mu, B(jj), &
                      particle(ii)%mass, particle(ii)%node(jj)%gamma)
              end do
           else
              do jj = 1, Nxi_local
                 call Einstein(particle(ii)%Nvz, particle(ii)%Nmu, &
                      particle(ii)%Vz, particle(ii)%mu, B(jj), &
                      particle(ii)%mass, particle(ii)%node(jj)%gamma)
              end do
           end if
        end if
     end do


  else
     iter_start = 0

! Initialise the distributions
     call GetTemperatureProfile(Nz,Nz_local,Nspecies, particle, &
          kTzprofile,kTpprofile,istart,iend)
     call GetVelocityProfile(Nz,Nz_local,Nspecies, particle, &
          vzprofile,istart,iend)

     ! compute gamma for initial vzoffset=0
     do ii = 1, Nspecies
        if ( particle(ii)%relativistic ) then
           do jj = 1, Nxi_local
              call Einstein(particle(ii)%Nvz, particle(ii)%Nmu, &
                   particle(ii)%Vz, particle(ii)%mu, B(jj), &
                   particle(ii)%mass, particle(ii)%node(jj)%gamma)
           end do
           ! For the purpose of generating boundary conditions, allocate and 
           ! compute gamma for the boundaries.
           if (neighbourLeft<0) then
              allocate(particle(ii)%node(0)%gamma(particle(ii)%Nvz, &
                   particle(ii)%Nmu))
              particle(ii)%node(0)%ivzoffset = &
                   nint((particle(ii)%vz0L-particle(ii)%vz0)/particle(ii)%dvz)
              call Einstein(particle(ii)%Nvz, particle(ii)%Nmu, &
                   dble(particle(ii)%node(0)%ivzoffset) * &
                   particle(ii)%dvz+particle(ii)%Vz, particle(ii)%mu, BBC(1), &
                   particle(ii)%mass, particle(ii)%node(0)%gamma)
           end if
           if (neighbourRight<0) then
              allocate( particle(ii)%node(Nz_local+1)%gamma &
                   (particle(ii)%Nvz,particle(ii)%Nmu) )
              particle(ii)%node(Nz_local+1)%ivzoffset = &
                   nint((particle(ii)%vz0R-particle(ii)%vz0)/ &
                   particle(ii)%dvz)
              call Einstein(particle(ii)%Nvz, particle(ii)%Nmu, &
                   dble(particle(ii)%node(Nz_local+1)%ivzoffset)* &
                   particle(ii)%dvz + &
                   particle(ii)%Vz, particle(ii)%mu, BBC(2), &
                   particle(ii)%mass, particle(ii)%node(Nz_local+1)%gamma)
           end if
        end if
     end do
!!$     call DumpGamma(Nxi_local,Nspecies,particle,iter_start,myid)

     select case (initialiser)
     case(1)
        call Initialise_f(Nxi_local, XI_local, BBC, B, nprofile, &
             kTzprofile, kTpprofile, vzprofile, Nspecies, particle, &
             neighbourLeft, neighbourRight,comm1d)

     case default
        if (myid == 0) then
           write (*,*) 'initialiser =', initialiser, ' is illegal!'
        end if
        call MPI_FINALIZE(ierr)
        stop
     end select

#ifdef _DEBUG_
     write (*,*) 'Done with Initialise_f_VR.  myid=', myid
#endif

! Dump the initial distribution
     attempts = 0
     call DumpDistr(Nxi_local,Nspecies,particle,0,myid,attempts,Nretries,neighbourRight)
#ifdef _DEBUG_
     write (*,*) 'Done with DumpDistr.  myid=', myid
#endif
     attempts = 0
     call DumpDistrBC(Nxi_local,Nspecies,particle,0, &
          neighbourLeft,neighbourRight,attempts,Nretries)
#ifdef _DEBUG_
     write (*,*) 'Done with DumpDistrBC.  myid=', myid
#endif

! Compute density and field and dump the initial state. 
     call ComputeInitialDensity(Nxi_local,Nspecies,particle,B,BBC,nprofile, &
          densities1D,rho)
#ifdef _DEBUG_
!!$  if (myid == 0) then
     write (*,*) 'Done with ComputeInitialDensity.  myid=', myid
!!$  end if
#endif

! Process 0 gets the voltage by interpolating the voltage control table
     if (myid == 0) then
        if (NVTable > 0) then
           call FindVoltageNow(iter_start, voltage, &
                NVTable, VTableIter, VTableVoltage)
        end if
     end if

! Make sure all processes are in phase
     call MPI_Barrier(MPI_COMM_WORLD, ierr)

     call Poisson(Nz_local, Z_local, dz_local, rho, BC_Poisson, B, BBC, &
          IntB, voltage, E0, E, epsilon, Nprocs, myid)

#ifdef _DEBUG_
!!$  if (myid == 0) then
     write (*,*) 'Done with initial Poisson.  myid=', myid
!!$  end if
#endif
     call ComputeMeanVz(Nxi_local,Nspecies,particle,meanVz)
     current = matmul(densities1D*meanVz,particle(:)%charge)
     ! collect initial fields in buffers
     fields_collected = fields_collected+1
     iterationbuffer(fields_collected) = 0
     Ebuffer(:,fields_collected) = E
     currentbuffer(:,fields_collected) = current
     density1Dbuffer(:,:,fields_collected) = densities1D
     ! if the buffers are full, then dump
     if (fields_collected >= fields_per_file) then
        attempts = 0
        call DumpEfield(Nxi_local,fields_per_file,Nspecies, &
             iterationbuffer,fields_collected,Ebuffer, &
             myid,attempts,Nretries)
        attempts = 0
        call DumpDensity(Nxi_local,fields_per_file,Nspecies, &
             iterationbuffer,fields_collected,density1Dbuffer,B,BBC, &
             myid,attempts,Nretries)
        attempts = 0
        call DumpCurrent(Nxi_local,fields_per_file,Nspecies, &
             iterationbuffer,fields_collected,currentbuffer, &
             myid,attempts,Nretries)
        iterationbuffer = 0
        fields_collected = 0
     end if

! This is in fact a leapfrog scheme, where positions are known 
! at half time steps, and velocities are known at whole time steps.
! Therefore, start by taking half a time step backwards in vz
     call AdvectionV(Nxi_local,E,dB,Nspecies,particle,gravity,-dt/2.0d0,0)
#ifdef _DEBUG_
!!$  if (myid == 0) then
     write (*,*) 'Done with initial AdvectionV.  myid=', myid
!!$  end if
#endif
  end if

! Compute density before the loop starts
#ifdef _TIMING_
  starttime = MPI_WTIME()
#endif
  call ComputeDensity(Nxi_local,Nspecies,particle,B,BBC,densities1D,rho)
#ifdef _TIMING_
  endtime = MPI_WTIME()
  if (myid .eq. 0) then
     write(*,*) 'ComputeDensity ', endtime-starttime, 'seconds'
  end if
  totalstarttime = MPI_WTIME()
#endif

! Compute the voltage across the resistor before we start
  call ComputeMeanVz(Nxi_local,Nspecies,particle,meanVz)
  resistI = sum(densities1D(1,:)*meanVz(1,:)*particle(:)%charge)
  resistV = resistance * resistI


! Here starts the main loop
  do iteration = iter_start+1, Niter
     if ( modulo(iteration,shift_test_period)==0) then
#ifdef _TIMING_
        starttime = MPI_WTIME()
#endif
        call ComputeMeanVz(Nxi_local,Nspecies,particle,meanVz)
#ifdef _TIMING_
        endtime = MPI_WTIME()
        if (myid .eq. 0) then
           write(*,*) 'ComputeMeanVz ', endtime-starttime, 'seconds'
        end if
#endif
#ifdef _TIMING_
        starttime = MPI_WTIME()
#endif

        resistI = sum(densities1D(1,:)*meanVz(1,:)*particle(:)%charge)
        resistV = resistance * resistI

        call TestAndShift(Nxi_local,Nspecies,particle,meanVz,densities1D,B)
#ifdef _TIMING_
        endtime = MPI_WTIME()
        if (myid .eq. 0) then
           write(*,*) 'TestAndShift ', endtime-starttime, 'seconds'
        end if
#endif
     end if

! Process 0 gets the voltage by interpolating the voltage control table
     if (myid == 0) then
        if (NVTable > 0) then
           call FindVoltageNow(iteration, voltage, &
                NVTable, VTableIter, VTableVoltage)
        end if
     end if

     ! Make sure all processes are in phase
     call MPI_Barrier(MPI_COMM_WORLD, ierr)

#ifdef _TIMING_
     starttime = MPI_WTIME()
#endif
     call Poisson(Nz_local, Z_local, dz_local, rho, BC_Poisson, B, BBC, &
          IntB, voltage+resistV, E0, E, epsilon, Nprocs, myid)
#ifdef _TIMING_
     endtime = MPI_WTIME()
     if (myid .eq. 0) then
        write(*,*) 'Poisson ', endtime-starttime, 'seconds'
     end if
#endif


! Advance vz
#ifdef _TIMING_
     starttime = MPI_WTIME()
#endif
     call AdvectionV(Nxi_local,E,dB,Nspecies,particle,gravity,dt,iteration)
#ifdef _TIMING_
     endtime = MPI_WTIME()
     if (myid .eq. 0) then
        write(*,*) 'AdvectionV ', endtime-starttime, 'seconds'
     end if
#endif


! Advance positions
#ifdef _TIMING_
     starttime = MPI_WTIME()
#endif
     ! write(*,*) 'Hallo?'
     call AdvectionXI(Nxi_local,XI_local,gp_local,Nspecies,particle,dt, &
          neighbourLeft, neighbourRight, comm1d, Boundary_I_file, Boundary_I_file2, &
					Boundary_I, Boundary_I2, iteration)
#ifdef _TIMING_
     endtime = MPI_WTIME()
     if (myid .eq. 0) then
        write(*,*) 'AdvectionXI ', endtime-starttime, 'seconds'
     end if
#endif

     ! Compute density
#ifdef _TIMING_
     starttime = MPI_WTIME()
#endif
     call ComputeDensity(Nxi_local,Nspecies,particle,B,BBC,densities1D,rho)
#ifdef _TIMING_
     endtime = MPI_WTIME()
     if (myid .eq. 0) then
        write(*,*) 'ComputeDensity ', endtime-starttime, 'seconds'
     end if
#endif


! Dump the distribution if it is time for that
     ! complete distribution
     if ( ( modulo(iteration,dump_period_distr)==0 .and. &
          iteration>=dump_start ) .or. iteration == Niter ) then
        attempts = 0
        call DumpDistr(Nxi_local,Nspecies,particle,iteration,myid, &
             attempts,Nretries,neighbourRight)
     end if
     ! reduced distribution
     if ( modulo(iteration,dump_period_distr_1v)==0 .and. &
          iteration>=dump_start ) then
        attempts = 0
        call DumpDistr1v(Nxi_local,Nspecies,particle,iteration,myid, &
             attempts,Nretries)
     end if
! Dump fields if so required
     if ( (modulo(iteration,dump_period_fields)==0 .and. &
          iteration>=dump_start) .or. iteration == Niter ) then
        call ComputeMeanVz(Nxi_local,Nspecies,particle,meanVz)
        current = matmul(densities1D*meanVz,particle(:)%charge)
        ! collect fields in buffers
        fields_collected = fields_collected+1
        iterationbuffer(fields_collected) = iteration
        Ebuffer(:,fields_collected) = E
        currentbuffer(:,fields_collected) = current
        density1Dbuffer(:,:,fields_collected) = densities1D
        ! if the buffers are full, then dump
        if (fields_collected >= fields_per_file) then
           attempts = 0
           call DumpEfield(Nxi_local,fields_per_file,Nspecies, &
                iterationbuffer,fields_collected,Ebuffer, &
                myid,attempts,Nretries)
           attempts = 0
           call DumpDensity(Nxi_local,fields_per_file,Nspecies, &
                iterationbuffer,fields_collected,density1Dbuffer,B,BBC, &
                myid,attempts,Nretries)
           attempts = 0
           call DumpCurrent(Nxi_local,fields_per_file,Nspecies, &
                iterationbuffer,fields_collected,currentbuffer, &
                myid,attempts,Nretries)
           iterationbuffer = 0
           fields_collected = 0
           if (myid == 0) then
              write (*,fmt='(a,i7,a, e12.6)') 'iteration = ', iteration, &
                   '  t = ', dble(iteration)*dt
           end if
        end if
     end if

     ! At regular intervals, dump something we can start from
     if ( modulo(iteration,dump_period_dump)==0 ) then
        attempts = 0
        call DumpDumpMPI(Nspecies, Nxi_local, iteration, particle, myid, &
             attempts, Nretries)
        if (myid == 0) then
           write (*,fmt='(a,i7,a, e12.6)') 'iteration = ', iteration, &
                '  t = ', dble(iteration)*dt
        end if
        if (exitafterdump) then
           if (myid == 0) then
              write (*,*) 'Exit after dump'
           end if
           exit
        end if
     end if

  end do

  ! If there is something left in the buffers, we dump it now.
  if (fields_collected >= 1) then
     attempts = 0
     call DumpEfield(Nxi_local,fields_per_file,Nspecies, &
          iterationbuffer,fields_collected,Ebuffer, &
          myid,attempts,Nretries)
     attempts = 0
     call DumpDensity(Nxi_local,fields_per_file,Nspecies, &
          iterationbuffer,fields_collected,density1Dbuffer,B,BBC, &
          myid,attempts,Nretries)
     attempts = 0
     call DumpCurrent(Nxi_local,fields_per_file,Nspecies, &
          iterationbuffer,fields_collected,currentbuffer, &
          myid,attempts,Nretries)
     iterationbuffer = 0
     fields_collected = 0
  end if

#ifdef _TIMING_
     totalendtime = MPI_WTIME()
     if (myid .eq. 0) then
        write(*,*) 'Main loop ', totalendtime-totalstarttime, 'seconds'
     end if
#endif

  if (myid == 0) then
     write (*,*) 'This is the end.'
  end if

  call MPI_FINALIZE(ierr)

end program ketchup

!------------------------------------------------------------------------
