subroutine GetGeneralInput(FromDir, dt, Niter, dump_period_fields, &
     fields_per_file, dump_period_distr, dump_period_distr_1v, dump_start, & 
     shift_test_period, resistance, Nz, zmin, zmax, Nspecies, const_a, &
     BC_Poisson, voltage, voltage_init, initialiser, E0, &
     startfromdumpfile, dump_period_dump, exitafterdump, transffilename, &
     voltagefilename, dump_period_distr_IonoBoundary)
  implicit none

  character(len=*) :: FromDir
  character transffilename*200, voltagefilename*200
  integer Niter, dump_period_fields, fields_per_file, dump_period_distr, &
       dump_period_distr_IonoBoundary, dump_period_distr_1v, dump_start, &
       shift_test_period, Nz, Nspecies, BC_Poisson, dump_period_dump, &
       initialiser
  double precision zmin, zmax, dt, resistance, const_a, &
       voltage, voltage_init, E0
  logical startfromdumpfile, exitafterdump

  character indata*132, filename*242
  integer i, j, k

  ! Defaults
  dt = 1.0d0
  Niter = 1
  dump_period_fields = 1
  fields_per_file = 1
  dump_period_distr = 1
  dump_period_distr_IonoBoundary = 1
  dump_period_distr_1v = 10000000
  dump_start = 1
  shift_test_period = 1
  resistance = 0.0d0
  Nz =10
  zmin = 0.0d0
  zmax = 1.0d0
  Nspecies = 0
  const_a = 10.0d0
  BC_Poisson=1
  voltage = 0.0d0
  initialiser = 1
  voltage_init = 0.0d0
  E0 = 0.0d0
  startfromdumpfile = .false.
  dump_period_dump = 1000
  exitafterdump = .false.
  transffilename = 'transfb2.dat'
  voltagefilename = 'No file name here, sorry.'

  filename=FromDir//'inputb6.m'
  open(unit=1,file=filename,status='old')

  read (1,'(a)') indata
  do while (index(indata,'%END')+index(indata,'%end') == 0)
     i=index(indata,'=')
     if (i>0) then
        do j=1,i-1
           if(indata(j:j) /= ' ') exit
        enddo
        ! If there is a semicolon on the line, ignore it and everything beyond.
        k = index(indata(i+1:132),';')
        if (k==0) then
           k = 133-i
        end if
        select case (indata(j:i-1))
        case ('dt')
           read (indata(i+1:i+k-1),*) dt
        case ('Niter')
           read (indata(i+1:i+k-1),*) Niter
        case ('dump_period_fields')
           read (indata(i+1:i+k-1),*) dump_period_fields
        case ('fields_per_file')
           read (indata(i+1:i+k-1),*) fields_per_file
        case ('dump_period_distr')
           read (indata(i+1:i+k-1),*) dump_period_distr
        case ('dump_period_distr_IonoBoundary')
           read (indata(i+1:i+k-1),*) dump_period_distr_IonoBoundary
        case ('dump_period_distr_1v')
           read (indata(i+1:i+k-1),*) dump_period_distr_1v
        case ('dump_start')
           read (indata(i+1:i+k-1),*) dump_start
        case ('shift_test_period')
           read (indata(i+1:i+k-1),*) shift_test_period
        case ('resistance')
           read (indata(i+1:i+k-1),*) resistance
        case ('Nz')
           read (indata(i+1:i+k-1),*) Nz
        case ('zmin')
           read (indata(i+1:i+k-1),*) zmin
        case ('zmax')
           read (indata(i+1:i+k-1),*) zmax
        case ('Nspecies')
           read (indata(i+1:i+k-1),*) Nspecies
        case ('const_a')
           read (indata(i+1:i+k-1),*) const_a
        case ('BC_Poisson')
           read (indata(i+1:i+k-1),*) BC_Poisson
        case ('E0')
           read (indata(i+1:i+k-1),*) E0
        case ('voltage')
           read (indata(i+1:i+k-1),*) voltage
        case ('initialiser')
           read (indata(i+1:i+k-1),*) initialiser
        case ('voltage_init')
           read (indata(i+1:i+k-1),*) voltage_init
        case ('startfromdumpfile')
           if (index(indata(i+1:i+k-1),'''yes''') > 0) then
              startfromdumpfile = .true.
           end if
        case ('dump_period_dump')
           read (indata(i+1:i+k-1),*) dump_period_dump
        case ('exit_after_dump')
           if (index(indata(i+1:i+k-1),'''yes''') > 0) then
              exitafterdump = .true.
           end if
        case ('transffilename')
           read (indata(i+1:i+k-1),*) transffilename
        case ('voltagefilename')
           read (indata(i+1:i+k-1),*) voltagefilename
        case default
           if (index(indata(j:j+1),'%')==0) then
              write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
           end if
        end select
     end if
     read (1,'(a)') indata
  end do

  close(1)

end subroutine GetGeneralInput

!------------------------------------------------------------------------

subroutine GetSpecies(FromDir,Nspecies, particle)

  use SpecificTypes

  implicit none

  character(len=*) :: FromDir
  integer Nspecies
  type(species) particle(Nspecies)

  character indata*132, filename*242
  integer ii, i, j, k

  filename=FromDir//'inputb6.m'
  open(unit=1,file=filename,status='old')

  do ii = 1, Nspecies
     ! Defaults
     particle(ii)%Nvz=2
     particle(ii)%vzmin=-1.0d0
     particle(ii)%vzmax=1.0d0
     particle(ii)%vzshifting = .false.
     particle(ii)%Nmu=2
     particle(ii)%mumin=0.0d0
     particle(ii)%mumax=1.0d0
     particle(ii)%muexp=1.0d0
     particle(ii)%mass=1.674927211d-27 ! default is neutrons with low density
     particle(ii)%charge=0.0d0
     particle(ii)%relativistic = .false.
     particle(ii)%n0=1.0d-99
     particle(ii)%vfromfile = .false.
     particle(ii)%vz0=0.0d0
     particle(ii)%Tfromfile = .false.
     particle(ii)%kTz=1.0d0
     particle(ii)%kTp=1.0d0
     particle(ii)%n0L=0.0d0
     particle(ii)%vz0L=0.0d0
     particle(ii)%kTzL=1.0d0
     particle(ii)%kTpL=1.0d0
     particle(ii)%lossconeL = .false.
     particle(ii)%n0R=0.0d0
     particle(ii)%vz0R=0.0d0
     particle(ii)%kTzR=1.0d0
     particle(ii)%kTpR=1.0d0
     particle(ii)%lossconeR = .false.

     read (1,'(a)') indata
     do while (index(indata,'%SPEC')+index(indata,'%spec') == 0)
        read (1,'(a)') indata
     end do
     do while (index(indata,'%END')+index(indata,'%end') == 0)
        i=index(indata,'=')
        if (i>0) then
           do j=1,i-1
              if(indata(j:j) /= ' ') exit
           enddo
           ! If there is a semicolon on the line, ignore it and 
           ! everything beyond.
           k = index(indata(i+1:132),';')
           if (k==0) then
              k = 133-i
           end if
           select case (indata(j:i-1))
           case ('Nvz')
              read (indata(i+1:i+k-1),*) particle(ii)%Nvz
           case ('vzmin')
              read (indata(i+1:i+k-1),*) particle(ii)%vzmin
           case ('vzmax')
              read (indata(i+1:i+k-1),*) particle(ii)%vzmax
           case ('vzshifting')
              if (index(indata(i+1:i+k-1),'''on''') > 0) then
                 particle(ii)%vzshifting = .true.
              end if
           case ('Nmu')
              read (indata(i+1:i+k-1),*) particle(ii)%Nmu
           case ('mumin')
              read (indata(i+1:i+k-1),*) particle(ii)%mumin
           case ('mumax')
              read (indata(i+1:i+k-1),*) particle(ii)%mumax
           case ('muexp')
              read (indata(i+1:i+k-1),*) particle(ii)%muexp
           case ('mass')
              read (indata(i+1:i+k-1),*) particle(ii)%mass
           case ('charge')
              read (indata(i+1:i+k-1),*) particle(ii)%charge
           case ('relativistic')
              if (index(indata(i+1:i+k-1),'''yes''') > 0) then
                 particle(ii)%relativistic = .true.
              end if
           case ('n0')
              read (indata(i+1:i+k-1),*) particle(ii)%n0
           case ('vfromfile')
              if (index(indata(i+1:i+k-1),'''yes''') > 0) then
                 particle(ii)%vfromfile = .true.
              end if
           case ('vz0')
              read (indata(i+1:i+k-1),*) particle(ii)%vz0
           case ('Tfromfile')
              if (index(indata(i+1:i+k-1),'''yes''') > 0) then
                 particle(ii)%Tfromfile = .true.
              end if
           case ('kTz')
              read (indata(i+1:i+k-1),*) particle(ii)%kTz
           case ('kTp')
              read (indata(i+1:i+k-1),*) particle(ii)%kTp

           case ('n0L')
              read (indata(i+1:i+k-1),*) particle(ii)%n0L
           case ('vz0L')
              read (indata(i+1:i+k-1),*) particle(ii)%vz0L
           case ('kTzL')
              read (indata(i+1:i+k-1),*) particle(ii)%kTzL
           case ('kTpL')
              read (indata(i+1:i+k-1),*) particle(ii)%kTpL
           case ('lossconeL')
              if (index(indata(i+1:i+k-1),'''yes''') > 0) then
                 particle(ii)%lossconeL = .true.
              end if

           case ('n0R')
              read (indata(i+1:i+k-1),*) particle(ii)%n0R
           case ('vz0R')
              read (indata(i+1:i+k-1),*) particle(ii)%vz0R
           case ('kTzR')
              read (indata(i+1:i+k-1),*) particle(ii)%kTzR
           case ('kTpR')
              read (indata(i+1:i+k-1),*) particle(ii)%kTpR
           case ('lossconeR')
              if (index(indata(i+1:i+k-1),'''yes''') > 0) then
                 particle(ii)%lossconeR = .true.
              end if

           case default
              if (index(indata(j:j+1),'%')==0) then
                 write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
              end if
           end select
        end if
        read (1,'(a)') indata
     end do
  end do
  
close(1)
 
end subroutine GetSpecies

!------------------------------------------------------------------------

subroutine GetTransform(FromDir, filename, Nxi, xi, zmin, zmax, Z, gp)

  implicit none

  character(len=*) :: FromDir, filename
  integer Nxi
  double precision xi(Nxi), zmin, zmax, Z(Nxi), gp(Nxi)

  integer over
  double precision indata(4)
  complex*16 i, one, intg, g0, UHPpole, UHPresid

  character completeFilename*242


  CompleteFilename = FromDir//filename

  i = (0.0d0, 1.0d0)
  one = (1.0d0, 0.0d0)
  intg= (0.0d0, 0.0d0)
  g0= (0.0d0, 0.0d0)
  over = 0
  Z=0.0d0
  gp=0.0d0

  open(unit=1,file=CompleteFilename,status='old')

  do while (over==0)
     read(1,*, end=99, err=99) indata
     UHPpole=indata(1) + indata(2)*i
     UHPresid=indata(3) +indata(4)*i

     g0 = g0 + UHPresid*log(-UHPpole) + conjg(UHPresid)*log(-conjg(UHPpole))
     intg = intg + UHPresid*log(one-UHPpole) + &
          conjg(UHPresid)*log(one-conjg(UHPpole))
     gp = gp + UHPresid/(xi - UHPpole) + conjg(UHPresid)/(xi-conjg(UHPpole))
     Z = Z + UHPresid*log(xi-UHPpole)+conjg(UHPresid)*log(xi-conjg(UHPpole))
  end do
99   over = 1
  close(1)
  intg = intg - g0
  gp = gp*(zmax-zmin)/intg
  Z = zmin + (-g0 + Z)*(zmax-zmin)/intg

end subroutine GetTransform

!------------------------------------------------------------------------

subroutine GetNVTable(FromDir, filename, NVTable)

  implicit none

  character(len=*), intent(in) :: FromDir, filename
  integer, intent(out):: NVTable

  integer jj, indataIter
  double precision indataVoltage

  character completeFilename*242

  CompleteFilename = FromDir//filename

  NVTable = 0

  open(unit=1,file=CompleteFilename,status='old', err=99)

  jj = 0
  ! First we count the number of data points
  do while (jj==0)
     read(1,*, end=99, err=99) indataIter, indataVoltage
     NVTable = NVTable + 1
  end do
99 close(1)

end subroutine GetNVTable

!------------------------------------------------------------------------

subroutine GetVTable(FromDir, filename, NVTable, VTableIter, VTableVoltage)

  implicit none

  character(len=*), intent(in) :: FromDir, filename
  integer, intent(in) :: NVTable
  integer, intent(out) :: VTableIter(NVTable)
  double precision, intent(out) :: VTableVoltage(NVTable)

  integer jj, previousIter

  character completeFilename*242

  CompleteFilename = FromDir//filename

  previousIter = 0
  VTableIter = 0
  VTableVoltage = 0.0d0

  open(unit=1,file=CompleteFilename,status='old',err=99)

  do jj = 1, NVTable
     read(1,*, end=99, err=99) VTableIter(jj), VTableVoltage(jj)
     if ( VTableIter(jj) <= previousIter .and. jj>1 ) then
        write (*,*) 'Warning!!!'
        write (*,*) 'The iteration numbers in the voltage ', &
             'table must be strictly increasing!'
     end if
     previousIter = VTableIter(jj)
  end do
99 close(1)


end subroutine GetVTable

!------------------------------------------------------------------------

subroutine GetDensityProfile(Nz,Nz_local,Nspecies, &
     nprofile,istart,iend)

  implicit none

  integer Nz, Nz_local, Nspecies, istart, iend

  double precision nprofile(Nz_local, Nspecies)

  double precision globalnprofile(Nz)

  character filename*30, Nznumber*6, snumber*2

  integer ii

  do ii= 1, Nspecies
! Assemble filename
     write(snumber,fmt='(i2.2)') ii
     write(Nznumber,fmt='(i6.6)') Nz
     filename='n0'//'s'//snumber//'Nz'//Nznumber//'.inp'

! If the profile is missing we initialise with 1e-20, which is to be 
! regarded as empty space
     globalnprofile = 1.0d-20

! Read input data
     open(unit=1,file=filename,status='old',err=99)
     read(1,*) globalnprofile
99   close(1)

     nprofile(:,ii) = globalnprofile(istart:iend)

  end do

end subroutine GetDensityProfile

!------------------------------------------------------------------------

subroutine GetVoltageProfile(FromDir, Nz, Nz_local, Vprofile, istart, iend)

  implicit none

  character(len=*) :: FromDir
  integer Nz, Nz_local, istart, iend

  double precision Vprofile(Nz_local), globalVprofile(Nz)

  character filename*230, Nznumber*6

  Vprofile = 0.0d0
  globalVprofile = 0.0d0

! Assemble filename
  write(Nznumber,fmt='(i6.6)') Nz
  filename=FromDir//'VprofNz'//Nznumber//'.inp'

! Read input data
  open(unit=1,file=filename,status='old')
  read(1,*) globalVprofile
  close(1)

! Pick out the part that belongs to this process
  Vprofile = globalVprofile(istart:iend)

end subroutine GetVoltageProfile

!------------------------------------------------------------------------

subroutine GetVelocityProfile(Nz,Nz_local,Nspecies, particle, &
     vzprofile,istart,iend)
  
  use SpecificTypes

  implicit none

  integer Nz, Nz_local, Nspecies, istart, iend

  double precision vzprofile(Nz_local, Nspecies)

  type(species) particle(Nspecies)

  double precision globalvzprofile(Nz)

  character vzfilename*30, Nznumber*6, snumber*2

  integer ii


  do ii= 1, Nspecies
     ! Are we reading files for this species?
     if (particle(ii)%vfromfile) then
! Assemble filenames
        write(snumber,fmt='(i2.2)') ii
        write(Nznumber,fmt='(i6.6)') Nz
        vzfilename='vzs'//snumber//'Nz'//Nznumber//'.inp'

! Read input data
        open(unit=1,file=vzfilename,status='old')
        read(1,*) globalvzprofile
        close(1)
        vzprofile(:,ii) = globalvzprofile(istart:iend)

     else
        vzprofile(:,ii) = particle(ii)%vz0
     end if
  end do

end subroutine GetVelocityProfile

!------------------------------------------------------------------------

subroutine GetTemperatureProfile(Nz,Nz_local,Nspecies, particle, &
     kTzprofile,kTpprofile,istart,iend)
  
  use SpecificTypes

  implicit none

  integer Nz, Nz_local, Nspecies, istart, iend

  double precision kTzprofile(Nz_local, Nspecies), &
       kTpprofile(Nz_local, Nspecies)

  type(species) particle(Nspecies)

  double precision globalTprofile(Nz)

  character Tzfilename*30, Tpfilename*30, Nznumber*6, snumber*2

  integer ii


  do ii= 1, Nspecies
     ! Are we reading files for this species?
     if (particle(ii)%Tfromfile) then

! Assemble filenames
        write(snumber,fmt='(i2.2)') ii
        write(Nznumber,fmt='(i6.6)') Nz
        Tzfilename='Tzs'//snumber//'Nz'//Nznumber//'.inp'
        Tpfilename='Tps'//snumber//'Nz'//Nznumber//'.inp'

! Read input data
        open(unit=1,file=Tzfilename,status='old')
        read(1,*) globalTprofile
        close(1)
        kTzprofile(:,ii) = globalTprofile(istart:iend)

        open(unit=1,file=Tpfilename,status='old')
        read(1,*) globalTprofile
        close(1)
        kTpprofile(:,ii) = globalTprofile(istart:iend)

     else
        kTzprofile(:,ii) = particle(ii)%kTz
        kTpprofile(:,ii) = particle(ii)%kTp
     end if
  end do

end subroutine GetTemperatureProfile

!------------------------------------------------------------------------

subroutine Initialise_f(Nxi, XI, BBC, B, nprofile, &
     kTzprofile, kTpprofile, vzprofile, Nspecies, particle, &
     neighbourLeft, neighbourRight, comm1d)

  use SpecificTypes

  use mpi

  implicit none
!  include 'mpif.h'

  type maxvector
     double precision, allocatable :: maxf0(:)
  end type maxvector

  integer Nxi, Nspecies
  double precision XI(Nxi), BBC(2), B(Nxi), nprofile(Nxi,Nspecies), &
       kTzprofile(Nxi,Nspecies), kTpprofile(Nxi,Nspecies), &
       vzprofile(Nxi,Nspecies)
  type(species) particle(Nspecies)
  integer neighbourLeft, neighbourRight, comm1d

  integer ixi, imu, ivz, ii, jj
  double precision elc, vzth(Nxi,Nspecies), vzthL(Nspecies), &
       vzthR(Nspecies), vpth(Nxi,Nspecies), vpthL(Nspecies), &
       vpthR(Nspecies), vzz, vz1, d1, d2, dxi, g

! Variables for use with mpi based communication
  integer (kind=MPI_ADDRESS_KIND) :: offsets(2)
  integer ierr, blockcounts(2), tag, oldtypes(2), &
       b_type_sendR, b_type_sendL, b_type_recvR, b_type_recvL, &
       istart, Nmess, rcount, fnodesize, fspecsize, maxf0size, &
       fspecshape(2), requestIndex, receiveLeftIndex, receiveRightIndex
  integer, allocatable :: req(:), status(:,:)
  
  type(bufferedDistribution) send_right, send_left, &
       receive_right, receive_left

  type(maxvector) maxf0local(Nspecies)
  double precision, allocatable :: maxf0outbuffer(:), maxf0inbuffer(:)

  parameter (elc=1.602176487d-19)

! Allocate communication buffers if necessary
  fnodesize = sum( particle(:)%Nvz * particle(:)%Nmu )
  Nmess = 0
  if (neighbourLeft>-1) then
     Nmess = Nmess + 2
     allocate( send_left%ivzoffsets(Nspecies*2) )
     allocate( send_left%fs(fnodesize*2 ))
     allocate( receive_left%ivzoffsets(Nspecies*2) )
     allocate( receive_left%fs(fnodesize*2) )
     send_left%ivzoffsets = 0
     send_left%fs = 0.0d0
     receive_left%ivzoffsets = 0
     receive_left%fs = 0.0d0
  end if
  if (neighbourRight>-1) then
     Nmess = Nmess + 2
     allocate( send_right%ivzoffsets(Nspecies*2) )
     allocate( send_right%fs(fnodesize*2) )
     allocate( receive_right%ivzoffsets(Nspecies*2) )
     allocate( receive_right%fs(fnodesize*2) )
     send_right%ivzoffsets = 0
     send_right%fs = 0.0d0
     receive_right%ivzoffsets = 0
     receive_right%fs = 0.0d0
  end if
  if (Nmess>0) then
     allocate(req(Nmess))
     allocate(status(MPI_STATUS_SIZE,Nmess))
  end if
! Buffers for exchange of maxf0 will always be necessary
  maxf0size = sum( particle(:)%Nmu )
  allocate ( maxf0outbuffer(maxf0size) )
  allocate ( maxf0inbuffer(maxf0size) )
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        allocate ( maxf0local(ii)%maxf0(particle(ii)%Nmu) )
     end if
  end do


! Build a few mpi data types for communication purposes
  oldtypes(1) = MPI_INTEGER
  blockcounts(1) = Nspecies*2
  oldtypes(2) = MPI_DOUBLE_PRECISION
  blockcounts(2) = fnodesize*2

  if (neighbourLeft>-1) then
     call MPI_GET_ADDRESS(receive_left%ivzoffsets, offsets(1), ierr)
     call MPI_GET_ADDRESS(receive_left%fs, offsets(2), ierr)
     offsets = offsets-offsets(1)
     call MPI_TYPE_CREATE_STRUCT &
          (2,blockcounts,offsets,oldtypes,b_type_recvL,ierr)
     call MPI_TYPE_COMMIT(b_type_recvL, ierr)
     
     call MPI_GET_ADDRESS(send_left%ivzoffsets, offsets(1), ierr)
     call MPI_GET_ADDRESS(send_left%fs, offsets(2), ierr)
     offsets = offsets-offsets(1)
     call MPI_TYPE_CREATE_STRUCT &
          (2,blockcounts,offsets,oldtypes,b_type_sendL,ierr)
     call MPI_TYPE_COMMIT(b_type_sendL, ierr)
  end if

  if (neighbourRight>-1) then
     call MPI_GET_ADDRESS(receive_right%ivzoffsets, offsets(1), ierr)
     call MPI_GET_ADDRESS(receive_right%fs, offsets(2), ierr)
     offsets = offsets-offsets(1)
     call MPI_TYPE_CREATE_STRUCT &
          (2,blockcounts,offsets,oldtypes,b_type_recvR,ierr)
     call MPI_TYPE_COMMIT(b_type_recvR, ierr)
     
     call MPI_GET_ADDRESS(send_right%ivzoffsets, offsets(1), ierr)
     call MPI_GET_ADDRESS(send_right%fs, offsets(2), ierr)
     offsets = offsets-offsets(1)
     call MPI_TYPE_CREATE_STRUCT &
          (2,blockcounts,offsets,oldtypes,b_type_sendR,ierr)
     call MPI_TYPE_COMMIT(b_type_sendR, ierr)
  end if

  dxi = XI(2)-XI(1)
  do ii = 1, Nspecies
     vzth(:,ii) = sqrt(elc*kTzprofile(:,ii)/particle(ii)%mass)
     vpth(:,ii) = sqrt(elc*kTpprofile(:,ii)/particle(ii)%mass)
  end do
  vzthL = sqrt(elc*particle%kTzL/particle%mass)
  vzthR = sqrt(elc*particle%kTzR/particle%mass)
  vpthL = sqrt(elc*particle%kTpL/particle%mass)
  vpthR = sqrt(elc*particle%kTpR/particle%mass)

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then  
        ! The interior
        do ixi = 1, Nxi
           ! No initial offset
           particle(ii)%node(ixi)%ivzoffset = 0
           do imu = 1, particle(ii)%Nmu
              do ivz = 1, particle(ii)%Nvz
                 if ( particle(ii)%relativistic ) then
                    d1 = particle(ii)%mu(imu)*B(ixi) / &
                         (particle(ii)%node(ixi)%gamma(ivz,imu) * &
                         particle(ii)%mass*vpth(ixi,ii)**2.0d0)
                 else
                    d1 = particle(ii)%mu(imu)*B(ixi) / &
                         (particle(ii)%mass*vpth(ixi,ii)**2.0d0)
                 end if
                 g = (1.0d0/(2.0d0*particle(ii)%mass*vpth(ixi,ii)))*exp(-d1)
                 vzz = particle(ii)%Vz(ivz)
                 vz1 = ( vzz - vzprofile(ixi,ii) ) / vzth(ixi,ii)
                 d2 = 0.5d0*(vz1*vz1)
!!$                 particle(ii)%node(ixi)%f(ivz,imu) = &
!!$                      nprofile(ixi,ii)*(BBC(1)/B(ixi))*g*exp(-d2)
                 particle(ii)%node(ixi)%f(ivz,imu) = g*exp(-d2)
              end do
           end do
        end do
     end if
  end do

! normalise the distribution function and find its maximum local value
  do ii = 1, Nspecies
     maxf0local(ii)%maxf0 = 0.0d0
     if (.not. isnan(particle(ii)%mass)) then  
        do ixi = 1, Nxi
           particle(ii)%node(ixi)%f = particle(ii)%node(ixi)%f * &
                nprofile(ixi,ii)*(BBC(1)/B(ixi))/ &
                (sum(matmul(particle(ii)%node(ixi)%f,particle(ii)%dmu)) * &
                particle(ii)%dvz)
           do imu = 1, particle(ii)%Nmu
              maxf0local(ii)%maxf0(imu) = max( maxf0local(ii)%maxf0(imu), &
                   maxval(particle(ii)%node(ixi)%f(:,imu)) )
           end do
        end do
     end if
  end do


! Exchange ghost points or generate boundary condition values

  ! Receive left
  receiveLeftIndex  = 0
  receiveRightIndex = 0
  rcount = 0
  if (neighbourLeft>-1) then
     rcount = rcount + 1
     tag = (neighbourLeft+1)*2 + 1
     call MPI_IRECV(receive_left%ivzoffsets, 1, b_type_recvL, neighbourLeft, &
          tag, comm1d, req(rcount), ierr)
     receiveLeftIndex  = rcount
  end if
  ! Receive right
  if (neighbourRight>-1) then
     rcount = rcount + 1
     tag = (neighbourRight-1)*2 + 2
     call MPI_IRECV(receive_right%ivzoffsets,1,b_type_recvR,neighbourRight, &
          tag, comm1d, req(rcount), ierr)
     receiveRightIndex  = rcount
  end if
  ! Send left
  if (neighbourLeft>-1) then
     istart = 1
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then  
           fspecsize = particle(ii)%Nvz * particle(ii)%Nmu
           send_left%ivzoffsets(ii) = particle(ii)%node(1)%ivzoffset
           send_left%fs(istart:istart+fspecsize-1) = &
                reshape(particle(ii)%node(1)%f,(/fspecsize/))
           send_left%ivzoffsets(Nspecies+ii) = particle(ii)%node(2)%ivzoffset
           send_left%fs(fnodesize+istart:fnodesize+istart+fspecsize-1) = &
                reshape(particle(ii)%node(2)%f,(/fspecsize/))
           istart= istart + fspecsize
        end if
     end do
     rcount = rcount + 1
     tag = neighbourLeft*2 + 2
     call MPI_ISEND(send_left%ivzoffsets, 1, b_type_sendL, neighbourLeft, &
          tag, comm1d, req(rcount), ierr)
  end if
  ! Send right
  if (neighbourRight>-1) then
     istart = 1
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then  
           fspecsize = particle(ii)%Nvz * particle(ii)%Nmu
           send_right%ivzoffsets(ii) = particle(ii)%node(Nxi-1)%ivzoffset
           send_right%fs(istart:istart+fspecsize-1) = &
                reshape(particle(ii)%node(Nxi-1)%f,(/fspecsize/))
           send_right%ivzoffsets(Nspecies+ii) = particle(ii)%node(Nxi)%ivzoffset
           send_right%fs(fnodesize+istart:fnodesize+istart+fspecsize-1) = &
                reshape( particle(ii)%node(Nxi)%f,(/fspecsize/))
           istart= istart + fspecsize
        end if
     end do
     rcount = rcount + 1
     tag = neighbourRight*2 + 1
     call MPI_ISEND(send_right%ivzoffsets, 1, b_type_sendR, neighbourRight, &
          tag, comm1d, req(rcount), ierr)
  end if

 
! Leftmost process generates boundary condition
  if (neighbourLeft<0) then
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
       ! If this is a loss cone distribution, remove what is in the loss cones.
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
           ! update maxf0 in case maximum is at boundary
           do imu = 1, particle(ii)%Nmu
              maxf0local(ii)%maxf0(imu) = max( maxf0local(ii)%maxf0(imu), &
                   maxval(particle(ii)%node(0)%f(:,imu)) )
           end do
        end if
     end do
  end if

 
! Rightmost process generates boundary condition
  if (neighbourRight<0) then
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
           ! update maxf0 in case maximum is at boundary
           do imu = 1, particle(ii)%Nmu
              maxf0local(ii)%maxf0(imu) = max( maxf0local(ii)%maxf0(imu), &
                   maxval(particle(ii)%node(Nxi+1)%f(:,imu)) )
           end do
        end if
     end do
  end if


! Take care of what was received
!!$  call MPI_Waitall(Nmess, req, status, ierr)
  do jj = 1, Nmess
     call MPI_Waitany(Nmess, req, requestIndex, status, ierr)

  ! Put received data in its place
     ! if we are receiving from the left
     if (requestIndex==receiveLeftIndex) then
#ifdef _DEBUG_
        if (neighbourLeft==-1) then
           write(*,*) 'Everything is wrong!!!!!!!!!!!'
           stop
        end if
#endif
        istart = 1
        do ii = 1, Nspecies
           if (.not. isnan(particle(ii)%mass)) then  
              fspecsize = particle(ii)%Nvz * particle(ii)%Nmu
              fspecshape = (/particle(ii)%Nvz, particle(ii)%Nmu/)
              particle(ii)%node(-1)%ivzoffset = receive_left%ivzoffsets(ii)
              particle(ii)%node(-1)%f = &
                   reshape(receive_left%fs(istart:istart+fspecsize-1), &
                   fspecshape)
              particle(ii)%node(0)%ivzoffset=receive_left%ivzoffsets &
                   (Nspecies+ii)
              particle(ii)%node(0)%f = &
                   reshape(receive_left%fs &
                   (fnodesize+istart:fnodesize+istart+fspecsize-1),fspecshape)
              istart= istart + fspecsize
           end if
        end do
     end if
     ! if we are receiving from the Right
     if (requestIndex==receiveRightIndex) then
#ifdef _DEBUG_
        if (neighbourRight==-1) then
           write(*,*) 'Everything is wrong!!!!!!!!!!!'
           stop
        end if
#endif
        istart = 1
        do ii = 1, Nspecies
           if (.not. isnan(particle(ii)%mass)) then  
              fspecsize = particle(ii)%Nvz * particle(ii)%Nmu
              fspecshape = (/particle(ii)%Nvz, particle(ii)%Nmu/)
              particle(ii)%node(Nxi+1)%ivzoffset = receive_right%ivzoffsets(ii)
              particle(ii)%node(Nxi+1)%f = &
                   reshape(receive_right%fs(istart:istart+fspecsize-1), &
                   fspecshape)
              particle(ii)%node(Nxi+2)%ivzoffset = &
                   receive_right%ivzoffsets(Nspecies+ii)
              particle(ii)%node(Nxi+2)%f = &
                   reshape(receive_right%fs(fnodesize+istart: &
                   fnodesize+istart+fspecsize-1),fspecshape)
              istart= istart + fspecsize
           end if
        end do
     end if
  end do


! Find global max from the other p.
  maxf0outbuffer = 0.0d0
  maxf0inbuffer = 0.0d0
  istart = 1
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        maxf0outbuffer(istart:istart+particle(ii)%Nmu-1) = maxf0local(ii)%maxf0
        istart = istart + particle(ii)%Nmu
     end if
  end do
  call MPI_Allreduce(maxf0outbuffer, maxf0inbuffer, maxf0size, &
       MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  istart = 1
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        particle(ii)%maxf0 = maxf0inbuffer(istart:istart+particle(ii)%Nmu-1)
        istart = istart + particle(ii)%Nmu
     end if
  end do


  if (neighbourLeft>-1) then
     call MPI_TYPE_FREE(b_type_recvL, ierr)
     call MPI_TYPE_FREE(b_type_sendL, ierr)
  end if
  if (neighbourRight>-1) then
     call MPI_TYPE_FREE(b_type_recvR, ierr)
     call MPI_TYPE_FREE(b_type_sendR, ierr)
  end if


end subroutine Initialise_f

!------------------------------------------------------------------------

subroutine GetNprocsOld(FromDir, Nprocs_old)

  implicit none

  character(len=*) :: FromDir
  integer Nprocs_old

  integer jj, kk, ll
  character filename*242, pnumber*4
  logical fileexists

  do jj = 0, 9999

     write (pnumber,fmt='(i4.4)') jj
     filename=FromDir//'dumps/dump_distribution.p'//pnumber// &
          '.ketchup.dump' 

     inquire (file=filename, exist = fileexists)

     if ( .not. fileexists) then
        Nprocs_old = jj
        return
     end if
  end do


end subroutine GetNprocsOld

!------------------------------------------------------------------------

subroutine DivideRegion( Nxi, Nprocs, myid, istart, iend )

  implicit none

  integer Nxi, Nprocs, myid, istart, iend
  integer Nxi_local, remainder

  Nxi_local = Nxi / Nprocs
  istart = myid * Nxi_local + 1
  remainder = mod(Nxi,Nprocs)
  istart = istart + min(myid,remainder)
  if (myid < remainder) then
     Nxi_local = Nxi_local + 1
  endif
  iend = istart + Nxi_local - 1
  if ((iend > Nxi) .or. (myid .eq. Nprocs-1)) then
     iend = Nxi
  end if
  return

end subroutine DivideRegion

!------------------------------------------------------------------------
