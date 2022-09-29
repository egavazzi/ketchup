subroutine AdvectionXI(Nxi,XI,gp,Nspecies,particle,dt, &
     neighbourLeft, neighbourRight, comm1d, Boundary_I_file, &
	   Boundary_I_file2, Boundary_I, Boundary_I2, iteration)

  use SpecificTypes
  use mpi

  implicit none
!  include 'mpif.h'

! Variables for use with mpi based communication
  integer Nxi, Nspecies
  double precision XI(Nxi), gp(Nxi), dt
  type(species) particle(Nspecies)
  integer neighbourLeft, neighbourRight, comm1d

  integer ixi, ivz, imu, ii, jj, ioffs, Nvz, Nmu, thisindex, ixistart, ixiend
  double precision alfa, dxi, ff, a0, a1, vztest(2)


  type temporaryDistribution
     double precision, allocatable :: f(:,:,:)
  end type temporaryDistribution

  type(temporaryDistribution) temp(Nspecies)

#ifdef _TIMING_
  double precision starttime, endtime
#endif


! Variables for use with mpi based communication
  integer (kind=MPI_ADDRESS_KIND) :: offsets(2)
  integer ierr, blockcounts(2), tag, oldtypes(2), &
       b_type_sendR, b_type_sendL, b_type_recvR, b_type_recvL, &
       istart, Nmess, rcount, fnodesize, fspecsize, fspecshape(2), &
       requestIndex, receiveLeftIndex, receiveRightIndex
  integer, allocatable :: req(:), status(:,:)

  type(bufferedDistribution) send_right, send_left, &
       receive_right, receive_left

  double precision Boundary_I_file(301,200,100), Boundary_I_file2(301,200,100), &
      Boundary_I(200,100), Boundary_I2(200,100) 
	integer iteration     


! Allocate communication buffers if necessary
  fnodesize = sum( particle(:)%Nvz * particle(:)%Nmu )
  Nmess = 0
  if (neighbourLeft>-1) then
     Nmess = Nmess + 2
     allocate( send_left%ivzoffsets(Nspecies*2) )
     allocate( send_left%fs(fnodesize*2 ))
     allocate( receive_left%ivzoffsets(Nspecies*2) )
     allocate( receive_left%fs(fnodesize*2) )
  end if
  if (neighbourRight>-1) then
     Nmess = Nmess + 2
     allocate( send_right%ivzoffsets(Nspecies*2) )
     allocate( send_right%fs(fnodesize*2) )
     allocate( receive_right%ivzoffsets(Nspecies*2) )
     allocate( receive_right%fs(fnodesize*2) )
  end if
  if (Nmess>0) then
     allocate(req(Nmess))
     allocate(status(MPI_STATUS_SIZE,Nmess))
  end if

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


! Exchange ghost points
  ! Receive left
  rcount = 0
  receiveLeftIndex  = 0
  receiveRightIndex = 0
  if (neighbourLeft>-1) then
     rcount = rcount + 1
     tag = (neighbourLeft+1)*2 + 1
     call MPI_IRECV(receive_left%ivzoffsets,1,b_type_recvL, neighbourLeft, &
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
           send_right%ivzoffsets(Nspecies+ii) = &
                particle(ii)%node(Nxi)%ivzoffset
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

  dxi = XI(2)-XI(1)
 
! Compute the new distribution function for all points that are not 
! affected by the other processes.
  ixistart = 3
  if (neighbourLeft < 0) then
     ixistart = 1
  end if
  ixiend = Nxi-2
  if (neighbourRight <0) then
     ixiend = Nxi
  end if

#ifdef _TIMING_
  starttime = MPI_WTIME()
#endif
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        Nvz = particle(ii)%Nvz
        Nmu = particle(ii)%Nmu
        allocate(temp(ii)%f(Nvz,Nmu,Nxi))

        ! Particle motion exceeding one spacestep in one timestep 
        ! is illegal, and this test has been introduced therefor.
#ifdef _SAFETY_
        vztest(1) = maxval((particle(ii)%Vz(Nvz) + &
             dble(particle(ii)%node(1:Nxi)%ivzoffset)*particle(ii)%dvz)/gp)
        vztest(2) = minval((particle(ii)%Vz(1) + &
             dble(particle(ii)%node(1:Nxi)%ivzoffset)*particle(ii)%dvz)/gp)
        if (maxval(abs(vztest))*dt/dxi > 1.0d0) then
           write (*,*) 'advectionXI, process ', neighbourLeft+1, &
                'vztest*dt/dxi/gp =', vztest*dt/dxi, 'species', ii, &
                'ivzoffset', particle(ii)%node(1:Nxi)%ivzoffset
           stop
        end if
#endif
        if ( particle(ii)%vzshifting ) then
           ! Special treatment for the interior points closest the 
           ! boundaries to take proper care of particles leaving the system.
           if (neighbourLeft < 0) then
              call UpdateXILeftmost(Nspecies, ii, Nvz, Nxi, Nmu, &
                   dxi, dt, gp, particle, temp)
              call UpdateXI(Nspecies, ii, Nvz, Nxi, Nmu, &
                   ixistart+1, ixiend, dxi, dt, gp, particle, temp)
           elseif (neighbourRight <0) then
              call UpdateXI(Nspecies, ii, Nvz, Nxi, Nmu, &
                   ixistart, ixiend-1, dxi, dt, gp, particle, temp)
              call UpdateXIRightmost(Nspecies, ii, Nvz, Nxi, Nmu, &
                   dxi, dt, gp, particle, temp)
           else
              call UpdateXI(Nspecies, ii, Nvz, Nxi, Nmu, &
                   ixistart, ixiend, dxi, dt, gp, particle, temp)
           end if
        else
           ! Special treatment for the interior points closest the 
           ! boundaries to take proper care of particles leaving the 
           ! system, also for the unshifted species, but these have 
           ! simpler subroutines.
           if (neighbourLeft < 0) then
              call UpdateXInoShiftLeftmost(Nspecies, ii, Nvz, Nxi, Nmu, &
                   dxi, dt, gp, particle, temp)
              call UpdateXInoShift(Nspecies, ii, Nvz, Nxi, Nmu, &
                   ixistart+1, ixiend, dxi, dt, gp, particle, temp)             
           elseif (neighbourRight <0) then
              call UpdateXInoShift(Nspecies, ii, Nvz, Nxi, Nmu, &
                   ixistart, ixiend-1, dxi, dt, gp, particle, temp)
              call UpdateXInoShiftRightmost(Nspecies, ii, Nvz, Nxi, Nmu, &
                   dxi, dt, gp, particle, temp, Boundary_I_file, &
                   Boundary_I_file2, Boundary_I, Boundary_I2, iteration)
           else
              call UpdateXInoShift(Nspecies, ii, Nvz, Nxi, Nmu, &
                   ixistart, ixiend, dxi, dt, gp, particle, temp)
           end if
        end if
    end if
  end do
#ifdef _TIMING_
  endtime = MPI_WTIME()
  if (neighbourLeft+1 .eq. 0) then
     write(*,*) 'AdvectionXI: main loop ', endtime-starttime, 'seconds'
  end if
#endif


! Take care of what was received
  do jj = 1, Nmess
     call MPI_Waitany(Nmess, req, requestIndex, status, ierr)

     ! Put received data in its place and compute for near boundary points
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
              particle(ii)%node(0)%f = reshape(receive_left%fs &
                   (fnodesize+istart:fnodesize+istart+fspecsize-1),fspecshape)
              istart= istart + fspecsize
           end if
        end do
        do ii = 1, Nspecies
           if (.not. isnan(particle(ii)%mass)) then
              Nvz = particle(ii)%Nvz
              Nmu = particle(ii)%Nmu
              ixistart = 1
              ixiend   = 2
              if ( particle(ii)%vzshifting ) then
                 call UpdateXI(Nspecies, ii, Nvz, Nxi, Nmu, &
                      ixistart, ixiend, dxi, dt, gp, particle, temp)
              else
                 call UpdateXInoShift(Nspecies, ii, Nvz, Nxi, Nmu, &
                      ixistart, ixiend, dxi, dt, gp, particle, temp)
              end if
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
              particle(ii)%node(Nxi+1)%ivzoffset = &
                   receive_right%ivzoffsets(ii)
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
        do ii = 1, Nspecies
           if (.not. isnan(particle(ii)%mass)) then
              Nmu = particle(ii)%Nmu
              Nvz = particle(ii)%Nvz
              ixistart = Nxi-1
              ixiend   = Nxi
              if ( particle(ii)%vzshifting ) then
                 call UpdateXI(Nspecies, ii, Nvz, Nxi, Nmu, &
                      ixistart, ixiend, dxi, dt, gp, particle, temp)
              else
                 call UpdateXInoShift(Nspecies, ii, Nvz, Nxi, Nmu, &
                      ixistart, ixiend, dxi, dt, gp, particle, temp)
              end if
           end if
        end do
     end if
  end do

  ! update f with the new approximation
#ifdef _TIMING_
  starttime = MPI_WTIME()
#endif
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        do imu = 1, particle(ii)%Nmu
           do ixi = 1, Nxi 
              particle(ii)%node(ixi)%f(:,imu) = temp(ii)%f(:,imu,ixi)
           end do
        end do
        deallocate(temp(ii)%f)
     end if
  end do
#ifdef _TIMING_
  endtime = MPI_WTIME()
  if (neighbourLeft+1 .eq. 0) then
     write(*,*) 'AdvectionXI: update ', endtime-starttime, 'seconds'
  end if
#endif

  if (neighbourLeft>-1) then
     call MPI_TYPE_FREE(b_type_recvL, ierr)
     call MPI_TYPE_FREE(b_type_sendL, ierr)
  end if
  if (neighbourRight>-1) then
     call MPI_TYPE_FREE(b_type_recvR, ierr)
     call MPI_TYPE_FREE(b_type_sendR, ierr)
  end if

end subroutine AdvectionXI

!------------------------------------------------------------------------

subroutine AdvectionV(Nxi,E,dB,Nspecies,particle,ag,dt,iteration)

  use SpecificTypes

  implicit none

  interface
     double precision function Flux(df0,middle,df1,alfa,a0,a1,maxf0)  
       double precision, intent(in) :: df0,middle,df1,alfa,a0,a1,maxf0
     end function Flux
  end interface

  integer, intent(in) :: Nxi, Nspecies, iteration
  double precision, intent(in) :: E(Nxi), dB(Nxi), ag(Nxi), dt
  type(species) particle(Nspecies)

  integer ixi, j, imu, jj, ii
  double precision alfa, tmp, dvz, df0, df1, ff, &
       epsiplus, epsiminus, a0, a1
  double precision, allocatable :: phi(:)

#ifdef _DEBUG_
  write (*,*) 'AdvectionV 1'
  do ixi = 1,Nxi
     if (isnan(E(ixi))) then
        write (*,*) 'AdvectionV ixi=',ixi,'E', E(ixi)
     end if
  end do
  write (*,*) 'AdvectionV 2'
#endif

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        if ( particle(ii)%relativistic ) then
           call UpdateVrelativistic(Nspecies, ii, Nxi, dt, particle, &
                E, dB, ag, iteration)
        else
           call UpdateV(Nspecies, ii, Nxi, dt, particle, E, dB, ag, iteration)
        end if

     end if
  end do

end subroutine AdvectionV

!------------------------------------------------------------------------

double precision function Flux(df0,middle,df1,alfa,a0,a1,maxf0)

  implicit none
  
  double precision, intent(in) :: df0,middle,df1,alfa,a0,a1,maxf0

  double precision epsiplus, epsiminus

  ! use limiters to ensure the positivity

  ! compute epsiminus for positivity
  if ( df0 < 0.0d0) then
     epsiminus = min( 1.0d0, -2.0d0*middle/df0 )
  elseif ( df0 > 0.0d0 ) then
     epsiminus =  min( 1.0d0, 2.0d0*(maxf0-middle)/df0 )
  else
     epsiminus = 1.0d0  
  end if
  ! compute epsiplus for positivity
  if ( df1 > 0.0d0) then
     epsiplus = min( 1.0d0, 2.0d0*middle/df1 )
  elseif ( df1 < 0.0d0 ) then
     epsiplus = min( 1.0d0, -2.0d0*(maxf0-middle)/df1 )
  else
     epsiplus = 1.0d0
  end if
  ! compute the flux
  Flux = alfa*( middle + a0*epsiplus*df1 + a1*epsiminus*df0 )

end function Flux

!------------------------------------------------------------------------

subroutine UpdateXI(Nspecies, ii, Nvz, Nxi, Nmu, ixistart, ixiend, &
     dxi, dt, gp, particle, temp)

  use SpecificTypes
  implicit none

  integer Nspecies, ii, Nvz, Nxi, Nmu, ixistart, ixiend
  double precision dxi, dt, gp(Nxi)
  type(species) particle(Nspecies)
  type temporaryDistribution
     double precision, allocatable :: f(:,:,:)
  end type temporaryDistribution

  type(temporaryDistribution) temp(Nspecies)

  integer ixi, ioffs, ivz, thisindex
  double precision alfa, a0, a1
  double precision df0v(Nmu), df1v(Nmu), df2v(Nmu), tmpv(Nmu), &
       phi0v(Nmu), phi1v(Nmu), term1v(Nmu), term2v(Nmu)

  do ixi = ixistart, ixiend
     ioffs = particle(ii)%node(ixi)%ivzoffset
     do ivz = 1, Nvz
        alfa = ( (particle(ii)%Vz(ivz) + &
             dble(ioffs)*particle(ii)%dvz)/gp(ixi) ) * dt / dxi
        if (alfa > 0.0d0) then
           ! Positive velocity
           a0 =  (1.0d0-alfa)*(2.0d0-alfa)/6.0d0
           a1 =  (1.0d0-alfa)*(1.0d0+alfa)/6.0d0

     ! There will be a lot of testing below, to see whether indeces 
     ! are out of range or not. The distribution function is assumed 
     ! to be zero for out of range values.

           thisindex = ivz+ioffs-particle(ii)%node(ixi-1)%ivzoffset
           if (thisindex<1 .or. thisindex>Nvz) then
              tmpv = 0.0d0
           else
              tmpv = particle(ii)%node(ixi-1)%f(thisindex,:) 
           end if
           thisindex = ivz+ioffs-particle(ii)%node(ixi-2)%ivzoffset
           if (thisindex<1 .or. thisindex>Nvz) then
              term2v = 0.0d0
           else
              term2v = particle(ii)%node(ixi-2)%f(thisindex,:)
           end if
           df0v = tmpv - term2v
           term1v = particle(ii)%node(ixi)%f(ivz,:)
           df1v = term1v - tmpv
           call Fluxv(Nmu,df0v,tmpv,df1v,alfa,a0,a1, &
                particle(ii)%maxf0,phi0v)
           tmpv = term1v
           thisindex = ivz+ioffs-particle(ii)%node(ixi+1)%ivzoffset
           if (thisindex<1 .or. thisindex>Nvz) then
              term1v = 0.0d0
           else
              term1v = particle(ii)%node(ixi+1)%f(thisindex,:)
           end if
           df2v = term1v - tmpv
           call Fluxv(Nmu,df1v,tmpv,df2v,alfa,a0,a1, &
                particle(ii)%maxf0,phi1v)
           
           ! compute the new approximation
           temp(ii)%f(ivz,:,ixi) = tmpv - (phi1v - phi0v)

        else
           ! Negative velocity
           a0 = -(1.0d0+alfa)*(1.0d0-alfa)/6.0d0
           a1 = -(1.0d0+alfa)*(2.0d0+alfa)/6.0d0
           
           tmpv = particle(ii)%node(ixi)%f(ivz,:) 
           thisindex = ivz+ioffs-particle(ii)%node(ixi-1)%ivzoffset
           if (thisindex<1 .or. thisindex>Nvz) then
              term2v = 0.0d0
           else
              term2v = particle(ii)%node(ixi-1)%f(thisindex,:)
           end if
           df0v = tmpv - term2v
           thisindex = ivz+ioffs-particle(ii)%node(ixi+1)%ivzoffset
           if (thisindex<1 .or. thisindex>Nvz) then
              term1v = 0.0d0
           else
              term1v = particle(ii)%node(ixi+1)%f(thisindex,:)
           end if
           df1v = term1v - tmpv
           call Fluxv(Nmu,df0v,tmpv,df1v,alfa,a0,a1, &
                particle(ii)%maxf0,phi0v)
           tmpv = term1v
           thisindex = ivz+ioffs-particle(ii)%node(ixi+2)%ivzoffset
           if (thisindex<1 .or. thisindex>Nvz) then
              term1v = 0.0d0
           else
              term1v = particle(ii)%node(ixi+2)%f(thisindex,:)
           end if
           df2v = term1v - tmpv
           call Fluxv(Nmu,df1v,tmpv,df2v,alfa,a0,a1, &
                particle(ii)%maxf0,phi1v)
 
           ! compute the new approximation
           temp(ii)%f(ivz,:,ixi) = particle(ii)%node(ixi)%f(ivz,:) - &
                (phi1v - phi0v)
        end if
     end do
  end do

end subroutine UpdateXI

!------------------------------------------------------------------------

subroutine Fluxv(Nmu,df0,middle,df1,alfa,a0,a1,maxf0,Flux)

  implicit none
  
  integer, intent(in) :: Nmu
  double precision, intent(in) :: df0(Nmu),middle(Nmu),df1(Nmu), &
       alfa,a0,a1,maxf0(Nmu)
  double precision :: Flux(Nmu)

  integer imu
  double precision epsiplus(Nmu), epsiminus(Nmu)

  ! use limiters to ensure the positivity

  ! compute epsiminus for positivity
  do imu = 1, Nmu
     if ( df0(imu) < 0.0d0) then
        epsiminus(imu) = min( 1.0d0, -2.0d0*middle(imu)/df0(imu) )
     elseif ( df0(imu) > 0.0d0 ) then
        epsiminus(imu) =  min( 1.0d0, 2.0d0*(maxf0(imu)-middle(imu))/df0(imu) )
     else
        epsiminus(imu) = 1.0d0
     end if
  end do

  ! compute epsiplus for positivity
  do imu = 1, Nmu
     if ( df1(imu) > 0.0d0) then
        epsiplus(imu) = min( 1.0d0, 2.0d0*middle(imu)/df1(imu) )
     elseif ( df1(imu) < 0.0d0 ) then
        epsiplus(imu) = min( 1.0d0, -2.0d0*(maxf0(imu)-middle(imu))/df1(imu) )
     else
        epsiplus(imu) = 1.0d0
     end if
  end do

  ! compute the flux
  Flux = alfa*( middle + a0*epsiplus*df1 + a1*epsiminus*df0 )

end subroutine Fluxv

!------------------------------------------------------------------------

subroutine UpdateXInoShift(Nspecies, ii, Nvz, Nxi, Nmu, ixistart, ixiend, &
     dxi, dt, gp, particle, temp)

  use SpecificTypes
  implicit none

  interface
     double precision function Flux(df0,middle,df1,alfa,a0,a1,maxf0)  
       double precision, intent(in) :: df0,middle,df1,alfa,a0,a1,maxf0
     end function Flux
  end interface

  integer Nspecies, ii, Nvz, Nxi, Nmu, ixistart, ixiend
  double precision dxi, dt, gp(Nxi)
  type(species) particle(Nspecies)
  type temporaryDistribution
     double precision, allocatable :: f(:,:,:)
  end type temporaryDistribution

  type(temporaryDistribution) temp(Nspecies)

  integer ixi, ivz, imu
  double precision alfa, a0, a1
  double precision df1, phi0, phi1

  do ixi = ixistart, ixiend
     do imu = 1, Nmu
        do ivz = 1, Nvz
           alfa = ( particle(ii)%Vz(ivz) / gp(ixi) ) * dt / dxi
           if (alfa > 0.0d0) then
              ! Positive velocities
              a0 =  (1.0d0-alfa)*(2.0d0-alfa)/6.0d0
              a1 =  (1.0d0-alfa)*(1.0d0+alfa)/6.0d0

              df1 = particle(ii)%node(ixi)%f(ivz,imu) - &
                   particle(ii)%node(ixi-1)%f(ivz,imu)
              phi0 = Flux(particle(ii)%node(ixi-1)%f(ivz,imu) - &
                   particle(ii)%node(ixi-2)%f(ivz,imu), &
                   particle(ii)%node(ixi-1)%f(ivz,imu),df1, &
                   alfa,a0,a1,particle(ii)%maxf0(imu))
              phi1 = Flux(df1,particle(ii)%node(ixi)%f(ivz,imu), &
                   particle(ii)%node(ixi+1)%f(ivz,imu) - &
                   particle(ii)%node(ixi)%f(ivz,imu), &
                   alfa,a0,a1,particle(ii)%maxf0(imu))

              ! compute the new approximation
              temp(ii)%f(ivz,imu,ixi) = &
                   particle(ii)%node(ixi)%f(ivz,imu) - (phi1 - phi0)
           else
              ! Negative velocities
              a0 = -(1.0d0+alfa)*(1.0d0-alfa)/6.0d0
              a1 = -(1.0d0+alfa)*(2.0d0+alfa)/6.0d0

              df1 = particle(ii)%node(ixi+1)%f(ivz,imu) - &
                   particle(ii)%node(ixi)%f(ivz,imu) 
              phi0 = Flux(particle(ii)%node(ixi)%f(ivz,imu) - &
                   particle(ii)%node(ixi-1)%f(ivz,imu), &
                   particle(ii)%node(ixi)%f(ivz,imu),df1, &
                   alfa,a0,a1,particle(ii)%maxf0(imu))
              phi1 = Flux(df1,particle(ii)%node(ixi+1)%f(ivz,imu), &
                   particle(ii)%node(ixi+2)%f(ivz,imu) - &
                   particle(ii)%node(ixi+1)%f(ivz,imu), &
                   alfa,a0,a1,particle(ii)%maxf0(imu))

              ! compute the new approximation
              temp(ii)%f(ivz,imu,ixi) = &
                   particle(ii)%node(ixi)%f(ivz,imu) - (phi1 - phi0)
           end if
        end do
     end do
  end do

end subroutine UpdateXInoShift

!------------------------------------------------------------------------

subroutine UpdateXInoShiftLeftmost(Nspecies, ii, Nvz, Nxi, Nmu, &
     dxi, dt, gp, particle, temp)

  use SpecificTypes
  implicit none

  interface
     double precision function Flux(df0,middle,df1,alfa,a0,a1,maxf0)  
       double precision, intent(in) :: df0,middle,df1,alfa,a0,a1,maxf0
     end function Flux
  end interface

  integer Nspecies, ii, Nvz, Nxi, Nmu
  double precision dxi, dt, gp(Nxi)
  type(species) particle(Nspecies)
  type temporaryDistribution
     double precision, allocatable :: f(:,:,:)
  end type temporaryDistribution

  type(temporaryDistribution) temp(Nspecies)

  integer ixi, ivz, imu
  double precision alfa, a0, a1
  double precision df1, phi0, phi1

  ixi = 1
  do imu = 1, Nmu
     do ivz = 1, Nvz
        alfa = ( particle(ii)%Vz(ivz) / gp(ixi) ) * dt / dxi
        if (alfa > 0.0d0) then
           ! Positive velocities
           a0 =  (1.0d0-alfa)*(2.0d0-alfa)/6.0d0
           a1 =  (1.0d0-alfa)*(1.0d0+alfa)/6.0d0

           df1 = particle(ii)%node(ixi)%f(ivz,imu) - &
                particle(ii)%node(ixi-1)%f(ivz,imu)
           phi0 = Flux(particle(ii)%node(ixi-1)%f(ivz,imu) - &
                particle(ii)%node(ixi-2)%f(ivz,imu), &
                particle(ii)%node(ixi-1)%f(ivz,imu),df1, &
                alfa,a0,a1,particle(ii)%maxf0(imu))
           phi1 = Flux(df1,particle(ii)%node(ixi)%f(ivz,imu), &
                particle(ii)%node(ixi+1)%f(ivz,imu) - &
                particle(ii)%node(ixi)%f(ivz,imu), &
                alfa,a0,a1,particle(ii)%maxf0(imu))

           ! compute the new approximation
           temp(ii)%f(ivz,imu,ixi) = &
                particle(ii)%node(ixi)%f(ivz,imu) - (phi1 - phi0)
        else
           ! Negative velocities
           a0 = -(1.0d0+alfa)*(1.0d0-alfa)/6.0d0
           a1 = -(1.0d0+alfa)*(2.0d0+alfa)/6.0d0

           df1 = particle(ii)%node(ixi+1)%f(ivz,imu) - &
                particle(ii)%node(ixi)%f(ivz,imu) 
           phi0 = Flux(0.0d0, particle(ii)%node(ixi)%f(ivz,imu),df1, &
                alfa,a0,a1,particle(ii)%maxf0(imu))
!!$           phi0 = particle(ii)%node(ixi)%f(ivz,imu) * alfa
           phi1 = Flux(df1,particle(ii)%node(ixi+1)%f(ivz,imu), &
                particle(ii)%node(ixi+2)%f(ivz,imu) - &
                particle(ii)%node(ixi+1)%f(ivz,imu), &
                alfa,a0,a1,particle(ii)%maxf0(imu))

           ! compute the new approximation
           temp(ii)%f(ivz,imu,ixi) = &
                particle(ii)%node(ixi)%f(ivz,imu) - (phi1 - phi0)
        end if
     end do
  end do


end subroutine UpdateXInoShiftLeftmost

!------------------------------------------------------------------------
subroutine UpdateXInoShiftRightmost(Nspecies, ii, Nvz, Nxi, Nmu, &
     dxi, dt, gp, particle, temp, Boundary_I_file, Boundary_I_file2, &
     Boundary_I, Boundary_I2, iteration)

  use SpecificTypes
  implicit none

  interface
     double precision function Flux(df0,middle,df1,alfa,a0,a1,maxf0)  
       double precision, intent(in) :: df0,middle,df1,alfa,a0,a1,maxf0
     end function Flux
  end interface

  integer Nspecies, ii, Nvz, Nxi, Nmu
  double precision dxi, dt, gp(Nxi)
  type(species) particle(Nspecies)
  type temporaryDistribution
     double precision, allocatable :: f(:,:,:)
  end type temporaryDistribution

  type(temporaryDistribution) temp(Nspecies)

  integer ixi, ivz, imu
  double precision alfa, a0, a1
  double precision df1, phi0, phi1

  double precision Boundary_I_file(301,200,100), Boundary_I_file2(301,200,100), &
            Boundary_I(200,100), Boundary_I2(200,100)
  double precision time, remainder
	logical updateIonosphere
	integer iteration, index4update

	! Test if it is time for ionospheric update
	time = dt*iteration
	updateIonosphere = .false.
	if (ii == 1) then
      if (iteration <= 1.2d4) then
         remainder = mod(iteration,40)
         if (remainder == 0.0d0) then
				index4update = 1 + iteration / 40
				updateIonosphere = .true.
				write(*,*) 'time to update! t =', time
            write(*,*) 'index4update =', index4update
			end if
      end if
	end if

	if (updateIonosphere) then
		do imu = 1,Nmu
			do ivz = 1, Nvz
				Boundary_I(ivz,imu) = Boundary_I_file(index4update,ivz,imu)
        Boundary_I2(ivz,imu) = Boundary_I_file2(index4update,ivz,imu)
			end do
		end do
	end if


  ixi = Nxi
  do imu = 1, Nmu
     do ivz = 1, Nvz
        alfa = ( particle(ii)%Vz(ivz) / gp(ixi) ) * dt / dxi
        if (alfa > 0.0d0) then
           ! Positive velocities
           a0 =  (1.0d0-alfa)*(2.0d0-alfa)/6.0d0
           a1 =  (1.0d0-alfa)*(1.0d0+alfa)/6.0d0

           df1 = particle(ii)%node(ixi)%f(ivz,imu) - &
                particle(ii)%node(ixi-1)%f(ivz,imu)
           phi0 = Flux(particle(ii)%node(ixi-1)%f(ivz,imu) - &
                particle(ii)%node(ixi-2)%f(ivz,imu), &
                particle(ii)%node(ixi-1)%f(ivz,imu),df1, &
                alfa,a0,a1,particle(ii)%maxf0(imu))
           phi1 = Flux(df1,particle(ii)%node(ixi)%f(ivz,imu), &
                0.0d0, alfa,a0,a1,particle(ii)%maxf0(imu))

           ! compute the new approximation
           temp(ii)%f(ivz,imu,ixi) = &
                particle(ii)%node(ixi)%f(ivz,imu) - (phi1 - phi0)
        else
           ! Negative velocities
           a0 = -(1.0d0+alfa)*(1.0d0-alfa)/6.0d0
           a1 = -(1.0d0+alfa)*(2.0d0+alfa)/6.0d0

           df1 = particle(ii)%node(ixi+1)%f(ivz,imu) - &
                particle(ii)%node(ixi)%f(ivz,imu) 
           phi0 = Flux(particle(ii)%node(ixi)%f(ivz,imu) - &
                particle(ii)%node(ixi-1)%f(ivz,imu), &
                particle(ii)%node(ixi)%f(ivz,imu),df1, &
                alfa,a0,a1,particle(ii)%maxf0(imu))
           phi1 = Flux(df1,particle(ii)%node(ixi+1)%f(ivz,imu), &
                particle(ii)%node(ixi+2)%f(ivz,imu) - &
                particle(ii)%node(ixi+1)%f(ivz,imu), &
                alfa,a0,a1,particle(ii)%maxf0(imu))

           ! compute the new approximation
           temp(ii)%f(ivz,imu,ixi) = &
                particle(ii)%node(ixi)%f(ivz,imu) - (phi1 - phi0)
          ! add the ionospheric response to the magnetospheric electrons		
            if (ii == 1) then
              temp(ii)%f(ivz,imu,ixi)   = Boundary_I2(ivz,imu)  !402km
              temp(ii)%f(ivz,imu,ixi-1) = Boundary_I2(ivz,imu)  !405km
              temp(ii)%f(ivz,imu,ixi-2) = Boundary_I2(ivz,imu)  !408km

              temp(ii)%f(ivz,imu,ixi-3) = Boundary_I(ivz,imu) !412km
              temp(ii)%f(ivz,imu,ixi-4) = Boundary_I(ivz,imu) !419km
              temp(ii)%f(ivz,imu,ixi-5) = Boundary_I(ivz,imu) !422km
              temp(ii)%f(ivz,imu,ixi-6) = Boundary_I(ivz,imu) !425km
            end if	
        end if

     end do
  end do


end subroutine UpdateXInoShiftRightmost

!------------------------------------------------------------------------

subroutine UpdateXILeftmost(Nspecies, ii, Nvz, Nxi, Nmu, &
     dxi, dt, gp, particle, temp)

  use SpecificTypes
  implicit none

  integer Nspecies, ii, Nvz, Nxi, Nmu
  double precision dxi, dt, gp(Nxi)
  type(species) particle(Nspecies)
  type temporaryDistribution
     double precision, allocatable :: f(:,:,:)
  end type temporaryDistribution

  type(temporaryDistribution) temp(Nspecies)

  integer ixi, ioffs, ivz, thisindex
  double precision alfa, a0, a1
  double precision df0v(Nmu), df1v(Nmu), df2v(Nmu), tmpv(Nmu), &
       phi0v(Nmu), phi1v(Nmu), term1v(Nmu), term2v(Nmu)

  ixi = 1
  ioffs = particle(ii)%node(ixi)%ivzoffset
  do ivz = 1, Nvz
     alfa = ( (particle(ii)%Vz(ivz) + &
          dble(ioffs)*particle(ii)%dvz)/gp(ixi) ) * dt / dxi
     if (alfa > 0.0d0) then
        ! Positive velocity
        a0 =  (1.0d0-alfa)*(2.0d0-alfa)/6.0d0
        a1 =  (1.0d0-alfa)*(1.0d0+alfa)/6.0d0

     ! There will be a lot of testing below, to see whether indeces 
     ! are out of range or not. The distribution function is assumed 
     ! to be zero for out of range values.

        thisindex = ivz+ioffs-particle(ii)%node(ixi-1)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           tmpv = 0.0d0
        else
           tmpv = particle(ii)%node(ixi-1)%f(thisindex,:) 
        end if
        thisindex = ivz+ioffs-particle(ii)%node(ixi-2)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           term2v = 0.0d0
        else
           term2v = particle(ii)%node(ixi-2)%f(thisindex,:)
        end if
        df0v = tmpv - term2v
        term1v = particle(ii)%node(ixi)%f(ivz,:)
        df1v = term1v - tmpv
        call Fluxv(Nmu,df0v,tmpv,df1v,alfa,a0,a1, &
             particle(ii)%maxf0,phi0v)
        tmpv = term1v
        thisindex = ivz+ioffs-particle(ii)%node(ixi+1)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           term1v = 0.0d0
        else
           term1v = particle(ii)%node(ixi+1)%f(thisindex,:)
        end if
        df2v = term1v - tmpv
        call Fluxv(Nmu,df1v,tmpv,df2v,alfa,a0,a1, &
             particle(ii)%maxf0,phi1v)
           
        ! compute the new approximation
        temp(ii)%f(ivz,:,ixi) = tmpv - (phi1v - phi0v)

     else
        ! Negative velocity
        a0 = -(1.0d0+alfa)*(1.0d0-alfa)/6.0d0
        a1 = -(1.0d0+alfa)*(2.0d0+alfa)/6.0d0
           
        tmpv = particle(ii)%node(ixi)%f(ivz,:) 
        df0v = 0.0d0
        thisindex = ivz+ioffs-particle(ii)%node(ixi+1)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           term1v = 0.0d0
        else
           term1v = particle(ii)%node(ixi+1)%f(thisindex,:)
        end if
        df1v = term1v - tmpv
        call Fluxv(Nmu,df0v,tmpv,df1v,alfa,a0,a1, &
             particle(ii)%maxf0,phi0v)
!!$        phi0v = tmpv * alfa
        tmpv = term1v
        thisindex = ivz+ioffs-particle(ii)%node(ixi+2)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           term1v = 0.0d0
        else
           term1v = particle(ii)%node(ixi+2)%f(thisindex,:)
        end if
        df2v = term1v - tmpv
        call Fluxv(Nmu,df1v,tmpv,df2v,alfa,a0,a1, &
             particle(ii)%maxf0,phi1v)

        ! compute the new approximation
        temp(ii)%f(ivz,:,ixi) = particle(ii)%node(ixi)%f(ivz,:) - &
             (phi1v - phi0v)
     end if
  end do
  
end subroutine UpdateXILeftmost

!------------------------------------------------------------------------

subroutine UpdateXIRightmost(Nspecies, ii, Nvz, Nxi, Nmu, &
     dxi, dt, gp, particle, temp)

  use SpecificTypes
  implicit none

  integer Nspecies, ii, Nvz, Nxi, Nmu
  double precision dxi, dt, gp(Nxi)
  type(species) particle(Nspecies)
  type temporaryDistribution
     double precision, allocatable :: f(:,:,:)
  end type temporaryDistribution

  type(temporaryDistribution) temp(Nspecies)

  integer ixi, ioffs, ivz, thisindex
  double precision alfa, a0, a1
  double precision df0v(Nmu), df1v(Nmu), df2v(Nmu), tmpv(Nmu), &
       phi0v(Nmu), phi1v(Nmu), term1v(Nmu), term2v(Nmu)

  ixi = Nxi
  ioffs = particle(ii)%node(ixi)%ivzoffset
  do ivz = 1, Nvz
     alfa = ( (particle(ii)%Vz(ivz) + &
          dble(ioffs)*particle(ii)%dvz)/gp(ixi) ) * dt / dxi
     if (alfa > 0.0d0) then
        ! Positive velocity
        a0 =  (1.0d0-alfa)*(2.0d0-alfa)/6.0d0
        a1 =  (1.0d0-alfa)*(1.0d0+alfa)/6.0d0

     ! There will be a lot of testing below, to see whether indeces 
     ! are out of range or not. The distribution function is assumed 
     ! to be zero for out of range values.

        thisindex = ivz+ioffs-particle(ii)%node(ixi-1)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           tmpv = 0.0d0
        else
           tmpv = particle(ii)%node(ixi-1)%f(thisindex,:) 
        end if
        thisindex = ivz+ioffs-particle(ii)%node(ixi-2)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           term2v = 0.0d0
        else
           term2v = particle(ii)%node(ixi-2)%f(thisindex,:)
        end if
        df0v = tmpv - term2v
        term1v = particle(ii)%node(ixi)%f(ivz,:)
        df1v = term1v - tmpv
        call Fluxv(Nmu,df0v,tmpv,df1v,alfa,a0,a1, &
             particle(ii)%maxf0,phi0v)
        tmpv = term1v
        df2v = 0.0d0
        call Fluxv(Nmu,df1v,tmpv,df2v,alfa,a0,a1, &
             particle(ii)%maxf0,phi1v)
!!$        phi1v = term1v * alfa
           
        ! compute the new approximation
        temp(ii)%f(ivz,:,ixi) = tmpv - (phi1v - phi0v)

     else
        ! Negative velocity
        a0 = -(1.0d0+alfa)*(1.0d0-alfa)/6.0d0
        a1 = -(1.0d0+alfa)*(2.0d0+alfa)/6.0d0

        tmpv = particle(ii)%node(ixi)%f(ivz,:) 
        thisindex = ivz+ioffs-particle(ii)%node(ixi-1)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           term2v = 0.0d0
        else
           term2v = particle(ii)%node(ixi-1)%f(thisindex,:)
        end if
        df0v = tmpv - term2v
        thisindex = ivz+ioffs-particle(ii)%node(ixi+1)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           term1v = 0.0d0
        else
           term1v = particle(ii)%node(ixi+1)%f(thisindex,:)
        end if
        df1v = term1v - tmpv
        call Fluxv(Nmu,df0v,tmpv,df1v,alfa,a0,a1, &
             particle(ii)%maxf0,phi0v)
        tmpv = term1v
        thisindex = ivz+ioffs-particle(ii)%node(ixi+2)%ivzoffset
        if (thisindex<1 .or. thisindex>Nvz) then
           term1v = 0.0d0
        else
           term1v = particle(ii)%node(ixi+2)%f(thisindex,:)
        end if
        df2v = term1v - tmpv
        call Fluxv(Nmu,df1v,tmpv,df2v,alfa,a0,a1, &
             particle(ii)%maxf0,phi1v)
 
        ! compute the new approximation
        temp(ii)%f(ivz,:,ixi) = particle(ii)%node(ixi)%f(ivz,:) - &
             (phi1v - phi0v)
     end if
  end do

end subroutine UpdateXIRightmost

!------------------------------------------------------------------------
subroutine UpdateV(Nspecies, ii, Nxi, dt, particle, E, dB, ag, iteration)

  use SpecificTypes

  implicit none

  interface
     double precision function Flux(df0,middle,df1,alfa,a0,a1,maxf0)  
       double precision, intent(in) :: df0,middle,df1,alfa,a0,a1,maxf0
     end function Flux
  end interface

  integer, intent(in) :: Nxi, Nspecies, iteration
  double precision, intent(in) :: E(Nxi), dB(Nxi), ag(Nxi), dt
  type(species) particle(Nspecies)

  integer ixi, j, imu, jj, ii
  double precision alfa, tmp, dvz, df0, df1, ff, &
       epsiplus, epsiminus, a0, a1
  double precision, allocatable :: phi(:)

#ifdef _SAFETY_
  integer maxjump
  logical warning


! initialise the warning system
  maxjump = 0
  warning = .false.
#endif

  allocate(phi(particle(ii)%Nvz+1))
  do imu = 1, particle(ii)%Nmu
     do ixi = 1, Nxi
        ff = ( particle(ii)%charge*E(ixi) - &
             particle(ii)%mu(imu)*dB(ixi) )*dt/particle(ii)%mass + &
             ag(ixi)*dt

        if ( ff > 0.0d0 ) then
! For positive accelerations:
           ! f(t,xi,vz) = 0 if |vz| -> Inf
           phi(particle(ii)%Nvz+1) = 0.0d0

           ! computation of the index for the shift:
           jj = int(ff/particle(ii)%dvz)
           alfa = (ff - dble(jj)*particle(ii)%dvz)/particle(ii)%dvz
#ifdef _SAFETY_
           if (jj>=1) then
              warning = .true.
              maxjump = max(maxjump, jj)
              if (jj>particle(ii)%Nvz) then
                 write (*,*) 'jj=', jj, 'Nvz=', particle(ii)%Nvz, &
                      'iteration =', iteration
                 write (*,*) &
                      'The simulation is out of control - will stop.'
                 stop
              end if
           end if
#endif
      ! coefficients of the polynomial of degree 2 
           a0 =  (1.0d0-alfa)*(2.0d0-alfa)/6.0d0
           a1 =  (1.0d0-alfa)*(1.0d0+alfa)/6.0d0

      ! computation of gradients for the third order reconstruction 	      

      ! case 1 : the flux is equal to zero for large |velocity|
      ! j = 1,...,jj+1
      ! f(t,xi,vz) = 0 if |vz| -> Inf
           do j = 1, jj+1
              phi(j) = 0.0d0
           end do

      ! case 2 : the rest
      ! initialise what will become df0=f(i,0) in the first step, 
      ! i.e. when j=jj+2
           df1 = particle(ii)%node(ixi)%f(1,imu)

           do j = jj+2, particle(ii)%Nvz
              tmp = particle(ii)%node(ixi)%f(j-jj-1,imu) 
              df0 = df1
              df1 = particle(ii)%node(ixi)%f(j-jj,imu)-tmp
              phi(j)=Flux(df0,tmp,df1,alfa,a0,a1,particle(ii)%maxf0(imu))
           end do

      ! computation of the new approximation	  	  
           do j = jj+1, particle(ii)%Nvz
              phi(j)=particle(ii)%node(ixi)%f(j-jj,imu)-(phi(j+1)-phi(j))
           end do

        elseif ( ff <= 0.0d0 ) then 
! For negative accelerations:
           ! f(t,xi,vz) = 0 if |vz| -> +oo
           phi(1) = 0.0d0
           phi(particle(ii)%Nvz+1) = 0.0d0
           ! computation of the index shift:   
           jj = int(-ff/particle(ii)%dvz)
#ifdef _SAFETY_
           if (jj>=1) then
              warning = .true.
              maxjump = max(maxjump, jj)
              if (jj>particle(ii)%Nvz) then
                 write (*,*) 'jj=', jj, 'Nvz=', particle(ii)%Nvz, &
                      'iteration =', iteration
                 write (*,*) &
                      'The simulation is out of control - will stop.'
                 stop
              end if
           end if
#endif
           alfa = (ff + dble(jj)*particle(ii)%dvz)/particle(ii)%dvz
           ! coefficients of the polynomial of degree 2 
           a0 = -(1.0d0+alfa)*(1.0d0-alfa)/6.0d0
           a1 = -(1.0d0+alfa)*(2.0d0+alfa)/6.0d0

	! computation of gradients for the third order reconstruction

	! case 1 : the flux is equal to zero for large velocity |vz|
	! j = Nvz-jj,...,Nvz-1
	! f(t,xi,vz) = 0 if |vz| -> +oo
           do j = particle(ii)%Nvz-jj+1, particle(ii)%Nvz+1
              phi(j) = 0.0d0
           end do

	! case 2 : the rest j = 2, Nvz-jj
	! Initialise what will become df0 in the 
	! first step, i.e. when j=Nvz-jj
           df1 = particle(ii)%node(ixi)%f(2+jj,imu) - &
                particle(ii)%node(ixi)%f(1+jj,imu)
           do j = 2, particle(ii)%Nvz-jj-1
              tmp = particle(ii)%node(ixi)%f(j+jj,imu)
              df0 = df1
              df1 = particle(ii)%node(ixi)%f(j+jj+1,imu)-tmp
              phi(j)=Flux(df0,tmp,df1,alfa,a0,a1,particle(ii)%maxf0(imu))
           end do
           tmp = particle(ii)%node(ixi)%f(particle(ii)%Nvz,imu)
           df0 = df1
           df1 = -tmp
           phi(particle(ii)%Nvz-jj) = &
                Flux(df0,tmp,df1,alfa,a0,a1,particle(ii)%maxf0(imu))

           ! compute the new approximation 
           do j = 1, particle(ii)%Nvz-jj
              phi(j)=particle(ii)%node(ixi)%f(j+jj,imu)-(phi(j+1)-phi(j))
           end do

        end if
        ! update f from phi
        particle(ii)%node(ixi)%f(:,imu) =  phi(1:particle(ii)%Nvz)
     end do
  end do

  deallocate(phi)

#ifdef _SAFETY_
  if ( warning ) then
     write (*,*) 'WARNING: AdvectionV, max(jj) =', maxjump, &
          'species', ii, 'iteration =', iteration
  end if
#endif


end subroutine UpdateV

!------------------------------------------------------------------------

subroutine UpdateVrelativistic(Nspecies, ii, Nxi, dt, particle, &
     E, dB, ag, iteration)

  use SpecificTypes

  implicit none

  interface
     double precision function Flux(df0,middle,df1,alfa,a0,a1,maxf0)  
       double precision, intent(in) :: df0,middle,df1,alfa,a0,a1,maxf0
     end function Flux
  end interface

  integer, intent(in) :: Nxi, Nspecies, iteration
  double precision, intent(in) :: E(Nxi), dB(Nxi), ag(Nxi), dt
  type(species) particle(Nspecies)

  integer ixi, ivz, imu, jj, ii
  double precision alfa, tmp, dvz, df0, df1, df2, ff, &
       epsiplus, epsiminus, a0, a1, phi0, phi1
  double precision, allocatable :: phi(:)

#ifdef _SAFETY_
  integer maxjump
  logical warning


! initialise the warning system
  maxjump = 0
  warning = .false.
#endif

  allocate(phi(particle(ii)%Nvz))
  do ixi = 1, Nxi
     do imu = 1, particle(ii)%Nmu
        do ivz = 1, particle(ii)%Nvz
           ff = ( particle(ii)%charge*E(ixi) - particle(ii)%mu(imu)*dB(ixi)/ &
                particle(ii)%node(ixi)%gamma(ivz,imu) )*dt/ &
                (particle(ii)%node(ixi)%gamma(ivz,imu)*particle(ii)%mass) + &
                ag(ixi)*dt

           if ( ff > 0.0d0 ) then
              ! For positive accelerations:

              ! computation of the index for the shift:
              jj = int(ff/particle(ii)%dvz)
              alfa = (ff - dble(jj)*particle(ii)%dvz)/particle(ii)%dvz
#ifdef _SAFETY_
              if (jj>=1) then
                 warning = .true.
                 maxjump = max(maxjump, jj)
                 if (jj>particle(ii)%Nvz) then
                    write (*,*) 'jj=', jj, 'Nvz=', particle(ii)%Nvz, &
                         'iteration =', iteration
                    write (*,*) &
                         'The simulation is out of control - will stop.'
                    stop
                 end if
              end if
#endif
      ! coefficients of the polynomial of degree 2 
              a0 =  (1.0d0-alfa)*(2.0d0-alfa)/6.0d0
              a1 =  (1.0d0-alfa)*(1.0d0+alfa)/6.0d0

              if ( ivz < jj+1 ) then
                 phi(ivz) = 0.0d0
              else
                 if ( ivz == jj+1 ) then
                    df1 = particle(ii)%node(ixi)%f(1,imu)
                    phi0 = 0.0d0
                 else
                    if ( ivz == jj+2 ) then
                       df0 = particle(ii)%node(ixi)%f(1,imu)
                    else
                       df0 = particle(ii)%node(ixi)%f(ivz-jj-1,imu) - &
                            particle(ii)%node(ixi)%f(ivz-jj-2,imu)
                    end if
                    df1 = particle(ii)%node(ixi)%f(ivz-jj,imu) - &
                         particle(ii)%node(ixi)%f(ivz-jj-1,imu)
                    phi0 = Flux(df0,particle(ii)%node(ixi)%f(ivz-jj-1,imu), &
                         df1,alfa,a0,a1,particle(ii)%maxf0(imu))
                 end if
                 if ( ivz < particle(ii)%Nvz ) then
                    phi1 = Flux(df1,particle(ii)%node(ixi)%f(ivz-jj,imu), &
                         particle(ii)%node(ixi)%f(ivz-jj+1,imu) - &
                         particle(ii)%node(ixi)%f(ivz-jj,imu), &
                         alfa,a0,a1,particle(ii)%maxf0(imu))
                 else
                    phi1 = 0.0d0
                 end if
                 phi(ivz)=particle(ii)%node(ixi)%f(ivz-jj,imu) - (phi1 - phi0)
              end if
           else 
! For negative accelerations:
              ! computation of the index shift:   
              jj = int(-ff/particle(ii)%dvz)
#ifdef _SAFETY_
              if (jj>=1) then
                 warning = .true.
                 maxjump = max(maxjump, jj)
                 if (jj>particle(ii)%Nvz) then
                    write (*,*) 'jj=', jj, 'Nvz=', particle(ii)%Nvz, &
                         'iteration =', iteration
                    write (*,*) &
                         'The simulation is out of control - will stop.'
                    stop
                 end if
              end if
#endif
              alfa = (ff + dble(jj)*particle(ii)%dvz)/particle(ii)%dvz
              a0 = -(1.0d0+alfa)*(1.0d0-alfa)/6.0d0
              a1 = -(1.0d0+alfa)*(2.0d0+alfa)/6.0d0

              if ( ivz > particle(ii)%Nvz-jj ) then
                 phi(ivz) = 0.0d0
              else
                 if ( ivz+jj == 1 ) then
                    df0 = particle(ii)%node(ixi)%f(1,imu)
                 else
                    df0 = particle(ii)%node(ixi)%f(ivz+jj,imu) - &
                         particle(ii)%node(ixi)%f(ivz+jj-1,imu)
                 end if
                 if ( ivz+jj ==  particle(ii)%Nvz ) then
                    df1 = -particle(ii)%node(ixi)%f(particle(ii)%Nvz,imu)
                 else
                    df1 = particle(ii)%node(ixi)%f(ivz+jj+1,imu) - &
                         particle(ii)%node(ixi)%f(ivz+jj,imu)
                 end if
                 if ( ivz == 1 ) then
                    phi0 = 0.0d0
                 else
                    phi0=Flux(df0,particle(ii)%node(ixi)%f(ivz+jj,imu), &
                         df1,alfa,a0,a1,particle(ii)%maxf0(imu))
                 end if
                 if ( ivz == particle(ii)%Nvz ) then
                    phi1 = 0.0d0
                 else
                    if ( ivz+jj ==  particle(ii)%Nvz ) then
                       tmp = 0.0d0
                       df2 = 0.0d0
                    elseif ( ivz+jj+1 ==  particle(ii)%Nvz ) then
                       tmp = particle(ii)%node(ixi)%f(particle(ii)%Nvz,imu)
                       df2 = -tmp
                    else
                       tmp = particle(ii)%node(ixi)%f(ivz+jj+1,imu)
                       df2 = particle(ii)%node(ixi)%f(ivz+jj+2,imu)-tmp
                    end if
                    phi1=Flux(df1,tmp,df2,alfa,a0,a1,particle(ii)%maxf0(imu))
                 end if

                 phi(ivz)=particle(ii)%node(ixi)%f(ivz+jj,imu) - (phi1 - phi0)
              end if
           end if
        end do
        ! update f from phi
        particle(ii)%node(ixi)%f(:,imu) =  phi
     end do
  end do

  deallocate(phi)

#ifdef _SAFETY_
  if ( warning ) then
     write (*,*) 'WARNING: AdvectionV, max(jj) =', maxjump, &
          'species', ii, 'iteration =', iteration
  end if
#endif

end subroutine UpdateVrelativistic

!------------------------------------------------------------------------
