! This program renormalises the distribution function by multiplying 
! it by a, for each species, constant factor, which is read from 
! the file renorm_par.m
!------------------------------------------------------------------------

program renormalise

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

! processes and renormalisation related variables
  integer Nprocs
  integer, allocatable :: istart(:), iend(:)
  double precision, allocatable :: renorm_factor(:)

! other variables
  integer ierr, iter_start, ii, jj, ixi, imu, ll
  integer Nretries, attempts


  write (*,*) '***********************************'
  write (*,*) '***   ketchup renormalisation   ***'
  write (*,*) '***********************************'

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

  allocate(renorm_factor(Nspecies))

  ! Specific parameters for each species
  allocate(particle(Nspecies))
  call GetSpecies('', Nspecies, particle)

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        allocate( particle(ii)%maxf0(particle(ii)%Nmu) )
     end if
  end do


  ! Renormalisation factors
  call GetRenormalisationFactors(Nspecies, renorm_factor)


  do ll = 0, Nprocs-1
     write (*,*) 'File', ll+1, 'out of', Nprocs

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
           do jj = -1, Nxi_local+2
              particle(ii)%node(jj)%f = particle(ii)%node(jj)%f * &
                   renorm_factor(ii)
           end do
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

end program renormalise

!------------------------------------------------------------------------

subroutine GetRenormalisationFactors(Nspecies, renorm_factor)

  implicit none

  integer Nspecies
  double precision renorm_factor(Nspecies)
  
  character indata*300
  integer i, j, k, ii

  renorm_factor = 1.0d0

  open(unit=1,file='renorm_par.m',status='old')

  ii = 0
  read (1,'(a)') indata
  do while ( (index(indata,'%END')+index(indata,'%end') == 0) .and. &
       ( ii < Nspecies ) )
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
        case ('factor')
           ii = ii + 1
           read (indata(i+1:i+k-1),*) renorm_factor(ii)
        case default
           if (index(indata(j:j+1),'%')==0) then
              write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
           end if
        end select
     end if
     read (1,'(a)') indata
  end do

  close(1)

end subroutine GetRenormalisationFactors
!------------------------------------------------------------------------

