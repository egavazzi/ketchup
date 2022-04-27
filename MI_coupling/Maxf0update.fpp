! This program reads dump files from a ketchup run, finds the maximum 
! value of the distribution function, and updates particle%maxf0
!------------------------------------------------------------------------

program maxf0update

  use SpecificTypes

  implicit none

  type maxvector
     double precision, allocatable :: maxf0(:)
  end type maxvector
! input data

  integer Niter, dump_period_fields, fields_per_file, dump_period_distr, &
       dump_period_distr_1v, dump_start, &
       shift_test_period, Nxi, Nspecies, BC_Poisson, &
       initialiser, dump_period_dump
  double precision zmin, zmax, dt, resistance, voltage, &
       voltage_init, E0
  double precision const_a
  character transffilename*200, voltagefilename*200
  logical startfromdumpfile, exitafterdump

  ! local grid
  integer Nxi_local

! particle species
  type(species), allocatable :: particle(:)

! processes and regenereation related variables
  integer Nprocs
  integer, allocatable :: istart(:), iend(:)
  type(maxvector), allocatable :: maxf0reference(:)

! other variables
  integer ierr, iter_start, ii, jj, ixi, imu, ll
  integer Nretries, attempts


  write (*,*) '********************************'
  write (*,*) '***   ketchup maxf0 update   ***'
  write (*,*) '********************************'

! Get old and new input data
  call GetGeneralInput('', &
       dt, Niter, dump_period_fields, fields_per_file, &
       dump_period_distr, dump_period_distr_1v, dump_start, &
       shift_test_period, resistance, Nxi, zmin, zmax, Nspecies, &
       const_a, BC_Poisson, voltage, voltage_init, &
       initialiser, E0, startfromdumpfile, dump_period_dump, &
       exitafterdump, transffilename, voltagefilename)

! Check for files to find Nprocs
  call GetNprocsOld('', Nprocs)

  allocate(istart(Nprocs))
  allocate(iend(Nprocs))
  do ll = 0, Nprocs-1
     call DivideRegion(Nxi, Nprocs, ll, istart(ll+1), iend(ll+1), ierr)
   end do

  allocate(maxf0reference(Nspecies))

  ! Specific parameters for each species
  allocate(particle(Nspecies))
  call GetSpecies('', Nspecies, particle)

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        allocate( particle(ii)%maxf0(particle(ii)%Nmu) )
        allocate ( maxf0reference(ii)%maxf0(particle(ii)%Nmu) )
        maxf0reference(ii)%maxf0 = 0.0d0
     end if
  end do


! Read all files and find the maxf0 vector
  do ll = 0, Nprocs-1

     ! set Nxi_local
     Nxi_local = iend(ll+1) - istart(ll+1) + 1

     ! allocate
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then
           allocate( particle(ii)%node(-1:Nxi_local+2) )
           do jj = -1, Nxi_local+2
              allocate(particle(ii)%node(jj)% &
                   f(particle(ii)%Nvz,particle(ii)%Nmu))
           end do
        end if
     end do

     ! load
     attempts = 0
     call LoadDump('',Nspecies, Nxi_local, iter_start, &
          particle, ll, attempts, Nretries)

     ! find maxf0
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then  
           do ixi = -1, Nxi_local+2
              do imu = 1, particle(ii)%Nmu
                 maxf0reference(ii)%maxf0(imu) = &
                      max(maxf0reference(ii)%maxf0(imu), &
                      maxval(particle(ii)%node(ixi)%f(:,imu)) )
              end do
           end do
        end if
     end do

     ! deallocate
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then
           deallocate( particle(ii)%node )
        end if
     end do

  end do
 
  write (*,*) 'Found maxf0'
  
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then  
        write (*,*) 'ii=',ii
        write (*,*) 'maxf0=', maxf0reference(ii)%maxf0
        write (*,*) '---'
     end if
  end do

 

  do ll = 0, Nprocs-1

     ! set Nxi_local
     Nxi_local = iend(ll+1) - istart(ll+1) + 1

     ! allocate
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then
           allocate( particle(ii)%node(-1:Nxi_local+2) )
           do jj = -1, Nxi_local+2
              allocate(particle(ii)%node(jj)% &
                   f(particle(ii)%Nvz,particle(ii)%Nmu))
           end do
        end if
     end do

     ! load
     attempts = 0
     call LoadDump('',Nspecies, Nxi_local, iter_start, particle, &
          ll, attempts, Nretries)

     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then  
           particle(ii)%maxf0 = maxf0reference(ii)%maxf0
        end if
     end do

     attempts = 0
     call DumpDump(Nspecies, Nxi_local, iter_start, particle, ll, &
          attempts, Nretries)

     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then
           deallocate( particle(ii)%node )
        end if
     end do

  end do

  write (*,*) 'Done'

end program maxf0update
