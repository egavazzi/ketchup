subroutine ComputeDensity(Nxi,Nspecies,particle,B,BBC,densities1D,rho)

  use SpecificTypes
  implicit none

  integer Nxi, Nspecies
  double precision  B(Nxi), BBC(2), densities1D(Nxi,Nspecies), rho(Nxi)
  type(species) particle(Nspecies)

  integer ii, ixi, imu

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        do ixi = 1, Nxi
           densities1D(ixi,ii) = &
                sum(matmul(particle(ii)%node(ixi)%f,particle(ii)%dmu)) * &
                particle(ii)%dvz * particle(ii)%n0
        end do
     end if
  end do

  rho = matmul(densities1D,particle%charge)

end subroutine ComputeDensity

!------------------------------------------------------------------------

subroutine ComputeInitialDensity(Nxi,Nspecies,particle,B,BBC,nprofile, &
     densities1D,rho)

  use SpecificTypes
  implicit none

  integer Nxi, Nspecies
  double precision nprofile(Nxi,Nspecies), densities1D(Nxi,Nspecies), &
       rho(Nxi), B(Nxi), BBC(2)
  type(species) particle(Nspecies)

  integer ii, ixi, imu

  do ii = 1, Nspecies
     if (isnan(particle(ii)%mass)) then
        do ixi = 1, Nxi
           densities1D(ixi,ii) = nprofile(ixi,ii)*particle(ii)%n0*BBC(1)/B(ixi)
        end do
     else
        do ixi = 1, Nxi
           densities1D(ixi,ii) = &
                sum(matmul(particle(ii)%node(ixi)%f,particle(ii)%dmu)) * &
                particle(ii)%dvz * particle(ii)%n0
        end do
     end if
  end do

  rho = matmul(densities1D,particle%charge)

end subroutine ComputeInitialDensity

!------------------------------------------------------------------------

subroutine ComputeMeanVz(Nxi,Nspecies,particle,meanVz)

  use SpecificTypes
  implicit none

  integer Nxi, Nspecies
  type(species) particle(Nspecies)
  double precision meanVz(Nxi,Nspecies)

  integer ii, jj, kk
  double precision deltavz, Vsum
  double precision, allocatable :: Vzf(:,:)

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        allocate (Vzf(particle(ii)%Nvz,particle(ii)%Nmu))
        do jj = 1, Nxi
           deltavz = particle(ii)%node(jj)%ivzoffset*particle(ii)%dvz
           do kk = 1, particle(ii)%Nmu
              Vzf(:,kk)=(particle(ii)%Vz+deltavz) * &
                   particle(ii)%node(jj)%f(:,kk)*particle(ii)%dmu(kk)
           end do
           Vsum = sum(matmul(particle(ii)%node(jj)%f,particle(ii)%dmu))
           if ( .not. Vsum .eq. 0.0d0 ) then
              meanVz(jj,ii) = sum(Vzf) / Vsum
           else
              meanVz(jj,ii) = 0.0d0
           end if
        end do
        deallocate (Vzf)
     end if
  end do

end subroutine ComputeMeanVz

!------------------------------------------------------------------------

subroutine TestAndShift(Nxi,Nspecies,particle,meanVz,densities1D,B)

  use SpecificTypes
  implicit none

  integer Nxi, Nspecies
  double precision meanVz(Nxi,Nspecies), densities1D(Nxi,Nspecies), B(Nxi), &
       vzdiff, neighbourdens(2), nsum
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
                    neighbourdens(1)=densities1D(1,ii)
                 else
                    neighbourdens(1)=densities1D(ixi-1,ii)
                 end if
                 if (ixi==Nxi) then
                    neighbourdens(2)=densities1D(Nxi,ii)
                 else
                    neighbourdens(2)=densities1D(ixi+1,ii)
                 end if
                 nsum = sum(neighbourdens) + densities1D(ixi,ii)
                 if ( .not. nsum .eq. 0.0d0 ) then
                    thisshift =nint( ( neighbouroffsets(1)*neighbourdens(1) + &
                         (particle(ii)%node(ixi)%ivzoffset + thisshift)* &
                         densities1D(ixi,ii) + &
                         neighbouroffsets(2)*neighbourdens(2) ) / nsum ) -  &
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
              ! We did shift, and hence we must update gamma for species
              ! that are relativistically modelled.
              if ( particle(ii)%relativistic ) then
                 call Einstein(particle(ii)%Nvz, particle(ii)%Nmu, &
                      dble(particle(ii)%node(ixi)%ivzoffset) * &
                      particle(ii)%dvz + particle(ii)%Vz, &
                      particle(ii)%mu, B(ixi), particle(ii)%mass, &
                      particle(ii)%node(ixi)%gamma)
              end if
           end if
        end do
        deallocate(newf)
     end if
  end do

end subroutine TestAndShift

!------------------------------------------------------------------------

subroutine FindVoltageNow(iteration,voltage,N,VTableIter,VTableVoltage)

  implicit none

  integer iteration, N, VTableIter(N)
  double precision voltage, VTableVoltage(N)

  integer Ng
  double precision di

  if ( iteration < VTableIter(1) ) then
     return
  end if

  if ( iteration >= VTableIter(N) ) then
     voltage = VTableVoltage(N)
     return
  end if

  do Ng = 1, N-1
     if ( iteration < VTableIter(Ng+1) ) then
        di=dble(iteration-VTableIter(Ng))/dble(VTableIter(Ng+1)-VTableIter(Ng))
        voltage = VTableVoltage(Ng) + di*(VTableVoltage(Ng+1)-VTableVoltage(Ng))
        return
     end if
  end do

end subroutine FindVoltageNow

!------------------------------------------------------------------------
subroutine Poisson(Nz, Z, dz, rho, BC_Poisson, B, BBC, &
     IntB, voltage, Ebc, E, epsilon, Nprocs, myid)

! This subroutine computes the electric field by integrating the 
! charge density. Each process computes its local accumulated 
! charge integral (aci), the integral of B (IntB), and the integral 
! of Baci, i.e. the product of B and aci. The result is sent to 
! process 0. 
! Having received these numbers, process 0 computes the global 
! accumulated charge integral, its value at the left hand border 
! of the z range of each process, and the constant E0, which is 
! equal to the electric field at the left hand boundary of the 
! whole system. 
! Process 0 then sends every other process these two numbers. 
! Process 0 is the only process that uses the voltage variable.

  use mpi
  implicit none

  integer Nz, myid, Nprocs, BC_Poisson
  double precision Z(Nz), dz(Nz), &
       rho(Nz), B(Nz), BBC(2), IntB, voltage, E(Nz), epsilon
  double precision, intent(in) :: Ebc

  integer ii
  double precision E0, aci(Nz), Baci(Nz), aci0, &
       acisends(Nprocs-1), acisNew(Nprocs-1), &
       Baciint(Nprocs-1), Bints(Nprocs-1)

  integer pp, tag, ierr, Nmess, messindex
  integer, allocatable :: req(:), status(:,:)
  double precision send(3), receive(2), send0(2,Nprocs-1), receive0(3,Nprocs-1)


  if (myid == 0) then
     Nmess = (Nprocs-1)
  else
     Nmess = 2
  end if
  allocate(req(Nmess))
  allocate(status(MPI_STATUS_SIZE,Nmess))

! Start by issuing the receiving commands for process 0
  if (myid==0) then
     ! for process 0
     do pp = 1, Nprocs-1
        tag = pp
        call MPI_IRECV(receive0(:,pp), 3, MPI_DOUBLE_PRECISION, &
             pp, tag, MPI_COMM_WORLD, req(pp), ierr)
     end do
  end if

! Then we perform the local accumulated charge integration
  aci(1) = dz(1) * rho(1)
  Baci(1) = 0.5d0*aci(1) * B(1)
  do ii = 2, Nz
     aci(ii) = aci(ii-1) + dz(ii)*rho(ii)
     Baci(ii) = (aci(ii-1) + 0.5d0*dz(ii)*rho(ii))*B(ii)
  end do

! Now we send aci(Nz) and more to process 0
  if (myid .ne. 0) then
     send(1) = aci(Nz)
     send(2) = sum(Baci*dz)
     send(3) = IntB
     tag = myid
     call MPI_ISEND(send, 3, MPI_DOUBLE_PRECISION, &
          0, tag, MPI_COMM_WORLD, req(1), ierr)
     ! Now we are ready to receive what 0 sends below
     tag = myid + Nprocs
     call MPI_IRECV(receive, 2, MPI_DOUBLE_PRECISION, &
          0, tag, MPI_COMM_WORLD, req(2), ierr)
  end if

! Process 0 computes E0 and aci(Nz) for each process and 
! distributes those values. The other processes receive said 
! values and update their aci vectors.
  if (myid==0) then
     if (Nprocs>1) then
        ! take care of received data
        call MPI_Waitall(Nmess, req, status, ierr)

        acisends  = receive0(1,:)
        Baciint   = receive0(2,:)
        Bints     = receive0(3,:)

        acisNew(1)  = aci(Nz)
        do pp = 2, Nprocs-1
           acisNew(pp)  = aci(Nz)  + sum(acisends(1:pp-1))
        end do

        select case (BC_Poisson)
        case (1)
           E0 = -(voltage*BBC(1) + ( sum(Baci*dz) + sum(Baciint) + &
                sum(acisNew*Bints) )/epsilon)/(sum(Bints)+IntB)
        case (2)
           E0 = Ebc
        case (3)
           E0 = Ebc*BBC(1)/BBC(2) - (aci(Nz) + sum(acisends))/epsilon
        case default
           write(*,*) 'Your BC_Poisson value is', BC_Poisson, &
                'We cannot handle that. Sorry.'
           stop
        end select
        aci0 = 0
        ! send E0 and acis
        do pp = 1, Nprocs-1
           send0(1,pp) = E0
           send0(2,pp) = acisNew(pp)
           tag = pp + Nprocs
           call MPI_ISEND(send0(:,pp), 2, MPI_DOUBLE_PRECISION, &
                pp, tag, MPI_COMM_WORLD, req(pp), ierr)
        end do
     else

        select case (BC_Poisson)
        case (1)
           E0 = -(voltage*BBC(1) + sum(Baci*dz)/epsilon)/IntB
        case (2)
           E0 = Ebc
        case (3)
           E0 = Ebc*BBC(1)/BBC(2) - aci(Nz)/epsilon
        case default
           write(*,*) 'Your BC_Poisson value is', BC_Poisson, &
                'We cannot handle that. Sorry.'
           stop
        end select
        aci0 = 0
     end if
  else
     ! wait for what process 0 has sent us, and put it in its place
     call MPI_Waitall(Nmess, req, status, ierr)
     E0  = receive(1)
     aci0 = receive(2)
     aci = aci + aci0
  end if

! Compute E 
  E(1) = (E0 + (aci0 + 0.5d0*dz(1)*rho(1))/epsilon)*B(1)/BBC(1)

  do ii = 2, Nz
     E(ii) = (E0 + (aci(ii-1) + 0.5d0*dz(ii)*rho(ii))/epsilon)*B(ii)/BBC(1)
  end do

! Process 0 waits for its sends to complete before leaving this subroutine
  if (myid==0) then
     call MPI_Waitall(Nmess, req, status, ierr)
  end if

end subroutine Poisson

!------------------------------------------------------------------------

subroutine Einstein(Nvz, Nmu, vz, mu, B, m0, gamma)

  implicit none

  integer Nvz, Nmu
  double precision vz(Nvz), mu(Nmu), B, m0, gamma(Nvz,Nmu)

  integer ivz, imu
  double precision c0

  parameter(c0 = 299792458.0d0)

  do imu = 1, Nmu
     do ivz = 1, Nvz
        gamma(ivz, imu) = sqrt((1.0d0+2.0d0*B*mu(imu)/((c0**2.0d0)*m0)) / &
             (1.0d0-(vz(ivz)/c0)**2.0d0) )
     end do
  end do

end subroutine Einstein

!------------------------------------------------------------------------

