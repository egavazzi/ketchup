recursive subroutine  DumpDump(Nspecies, Nxi, iteration, particle, myid, &
     attemptno,Nretries)

  use SpecificTypes
  use ifport

  implicit none

  integer Nspecies, Nxi, iteration, myid, attemptno, Nretries
  type(species) particle(Nspecies)

  integer ii, ixi, ierr
  character filename*42, tmpfile*44, pnumber*4

  attemptno = attemptno + 1

! Process 0 writes the iteration number
  if (myid==0) then
     tmpfile = 'dumps/XXdump_iteration.ketchup.dump'
     filename = 'dumps/dump_iteration.ketchup.dump'
     open(unit=1,file=tmpfile,status='replace',err=97)
     write (1,fmt='(I7)',err=98) iteration
     close(1,err=99)
     ierr = rename(tmpfile,filename)
  end if

! Everybody writes the distribution function
  write(pnumber,fmt='(i4.4)') myid
  tmpfile='dumps/XXdump_distribution.p'//pnumber//'.ketchup.dump' 
  filename='dumps/dump_distribution.p'//pnumber//'.ketchup.dump' 
  open(unit=1,form='unformatted',file=tmpfile,status='replace',err=97)
  write (1,err=98) (particle(ii)%maxf0, ii=1, Nspecies)
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        write (1,err=98) (particle(ii)%node(ixi)%f, &
             particle(ii)%node(ixi)%ivzoffset, ixi=-1, Nxi+2)
     end if
  end do
  close (1,err=99)

  ierr = rename(tmpfile,filename)
  return

! Error handling section
97 write (*,*) 'DumpDump: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpDump: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpDump: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpDump(Nspecies, Nxi, iteration, particle, myid, &
          attemptno,Nretries)
  end if

end subroutine DumpDump

!------------------------------------------------------------------------

recursive subroutine  DumpDumpMPI(Nspecies, Nxi, iteration, particle, myid, &
     attemptno,Nretries)

  use SpecificTypes
  use mpi
  use ifport

  implicit none

  integer Nspecies, Nxi, iteration, myid, attemptno, Nretries
  type(species) particle(Nspecies)

  integer ii, ixi, ierr
  character filename*42, tmpfile*44, pnumber*4

  attemptno = attemptno + 1

! Process 0 writes the iteration number
  if (myid==0) then
     tmpfile = 'dumps/XXdump_iteration.ketchup.dump'
     filename = 'dumps/dump_iteration.ketchup.dump'
     open(unit=1,file=tmpfile,status='replace',err=97)
     write (1,fmt='(I7)',err=98) iteration
     close(1,err=99)
     ierr = rename(tmpfile,filename)

  end if

! Everybody writes the distribution function
  write(pnumber,fmt='(i4.4)') myid
  tmpfile='dumps/XXdump_distribution.p'//pnumber//'.ketchup.dump' 
  filename='dumps/dump_distribution.p'//pnumber//'.ketchup.dump' 
  open(unit=1,form='unformatted',file=tmpfile,status='replace',err=97)
  write (1,err=98) (particle(ii)%maxf0, ii=1, Nspecies)
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        write (1,err=98) (particle(ii)%node(ixi)%f, &
             particle(ii)%node(ixi)%ivzoffset, ixi=-1, Nxi+2)
     end if
  end do
  close (1,err=99)

  ! Make sure all processes have written their temporary files 
  ! and then rename them.
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ierr = rename(tmpfile,filename)

  return

! Error handling section
97 write (*,*) 'DumpDump: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpDump: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpDump: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpDumpMPI(Nspecies, Nxi, iteration, particle, myid, &
          attemptno,Nretries)
  end if

end subroutine DumpDumpMPI

!------------------------------------------------------------------------

recursive subroutine  LoadDump(FromDir,Nspecies,Nxi,iteration,particle,myid, &
     attemptno,Nretries)
  
  use SpecificTypes

  implicit none

  character(len=*) :: FromDir
  integer Nspecies, Nxi, iteration, myid, attemptno, Nretries
  type(species) particle(Nspecies)

  integer ii, ixi
  character filename*242, pnumber*4

  attemptno = attemptno + 1

! Read the iteration number
  filename = FromDir//'dumps/dump_iteration.ketchup.dump'
  open(unit=1,file=filename,status='old',err=97)
  read (1,fmt='(I7)',err=98) iteration
  close(1,err=99)

! Read the distribution function
  write(pnumber,fmt='(i4.4)') myid
  filename= FromDir//'dumps/dump_distribution.p'//pnumber//'.ketchup.dump' 
  open(unit=1,form='unformatted',file=filename,status='old',err=97)
  read (1,err=98) (particle(ii)%maxf0, ii=1,Nspecies)
  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        read (1,err=98) (particle(ii)%node(ixi)%f, &
             particle(ii)%node(ixi)%ivzoffset, ixi=-1, Nxi+2)
     end if
  end do

  close(1,err=99)
  return

! Error handling section
97 write (*,*) 'LoadDump: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'LoadDump: error in read statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'LoadDump: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call LoadDump(FromDir,Nspecies,Nxi,iteration,particle,myid, &
          attemptno,Nretries)
  else
     stop
  end if

end subroutine LoadDump

!------------------------------------------------------------------------

recursive subroutine DumpDistr(Nxi,Nspecies,particle,iteration,myid, &
     attemptno,Nretries,neighbourRight)

  use SpecificTypes
  use ifport

  implicit none

  integer Nxi, Nspecies, iteration, myid, attemptno, Nretries, neighbourRight
  type(species) particle(Nspecies)

  character filename*60, tmpfile*62, fnumber*7, snumber*2, pnumber*4, form*20
  integer ixi, ivz, ii, imu, ierr

  attemptno = attemptno + 1

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then

        ! Assemble filename
        write(fnumber,fmt='(i7.7)') iteration
        write(snumber,fmt='(i2.2)') ii
        write(pnumber,fmt='(i4.4)') myid
        tmpfile='outp/datfiles/fzvzmu/XXfzvzmu'// &
             fnumber//'s'//snumber//'p'//pnumber//'.ketchup.dat' 
        filename='outp/datfiles/fzvzmu/fzvzmu'// &
             fnumber//'s'//snumber//'p'//pnumber//'.ketchup.dat' 

        write (form,fmt='(a,i5.5,a)') '(', particle(ii)%Nmu, '(E13.6E3,a))'
  
        
        ! open file and write f(xi,vz,mu)
        ! Format is 
        ! ivzoffset = offset
        ! blank line
        ! f(vz,mu)
        ! blank line
        ! and then repeat for the next node, and the next, and ...
        open(unit=1,file=tmpfile,status='replace',err=97)
        !if (neighbourRight<0) then
        !  do ixi = 1, Nxi
        !    write (1,*,err=98) 'ivzoffset = ', particle(ii)%node(ixi)%ivzoffset
        !    write (1,*,err=98)
        !    do ivz = 1, particle(ii)%Nvz
        !        write (1,fmt=form,err=98)(particle(ii)%node(ixi)%f(ivz,imu), &
        !            ' ', imu=1,particle(ii)%Nmu)
        !    end do
        !    write (1,*,err=98)
        !  end do
        !else
          do ixi = 1, Nxi
            write (1,*,err=98) 'ivzoffset = ', particle(ii)%node(ixi)%ivzoffset
            write (1,*,err=98)
            do ivz = 1, particle(ii)%Nvz
                write (1,fmt=form,err=98)(particle(ii)%node(ixi)%f(ivz,imu), &
                    ' ', imu=1,particle(ii)%Nmu)
            end do
            write (1,*,err=98)
          end do
        !end if
        close(1,err=99)
        ierr = rename(tmpfile,filename)

     end if
  end do

  return

! Error handling section
97 write (*,*) 'DumpDistr: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpDistr: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpDistr: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpDistr(Nxi,Nspecies,particle,iteration,myid,attemptno,Nretries,neighbourRight)
  end if

end subroutine DumpDistr

!------------------------------------------------------------------------

recursive subroutine DumpDistr1v(Nxi,Nspecies,particle,iteration,myid, &
     attemptno,Nretries)

  use SpecificTypes
  use ifport

  implicit none

  integer Nxi, Nspecies, iteration, myid, attemptno, Nretries
  type(species) particle(Nspecies)

  character filename*60, tmpfile*62, fnumber*7, snumber*2, pnumber*4
  integer ixi, ivz, ii, ierr
  double precision, allocatable :: fvz(:)

  attemptno = attemptno + 1

  do ii = 1, Nspecies
     if (.not. isnan(particle(ii)%mass)) then
        ! Assemble filename
        write(fnumber,fmt='(i7.7)') iteration
        write(snumber,fmt='(i2.2)') ii
        write(pnumber,fmt='(i4.4)') myid
        tmpfile='outp/datfiles/fzvz/XXfzvz'// &
             fnumber//'s'//snumber//'p'//pnumber//'.ketchup.dat' 
        filename='outp/datfiles/fzvz/fzvz'// &
             fnumber//'s'//snumber//'p'//pnumber//'.ketchup.dat' 

        allocate(fvz(particle(ii)%Nvz))

        ! open file, compute and write f(xi,vz)
        ! Format is 
        ! ivzoffset = offset
        ! blank line
        ! f(vz)
        ! blank line
        ! and then repeat for the next node, and the next, and ...
        open(unit=1,file=tmpfile,status='replace',err=97)
        do ixi = 1, Nxi
           fvz=matmul(particle(ii)%node(ixi)%f,particle(ii)%dmu)
           write (1,*,err=98) 'ivzoffset = ', particle(ii)%node(ixi)%ivzoffset
           write (1,*,err=98)
           do ivz = 1, particle(ii)%Nvz
              write (1,'(E13.6E3)',err=98) fvz(ivz)
           end do
           write (1,*,err=98)
        end do
        close(1,err=99)
        ierr = rename(tmpfile,filename)

        deallocate(fvz)

     end if
  end do

  return

! Error handling section
97 write (*,*) 'DumpDistr1v: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpDistr1v: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpDistr1v: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpDistr1v(Nxi,Nspecies,particle,iteration,myid,attemptno,Nretries)
  end if

end subroutine DumpDistr1v

!------------------------------------------------------------------------

recursive subroutine DumpDistrBC(Nxi,Nspecies,particle,iteration, &
     neighbourLeft,neighbourRight,attemptno,Nretries)

! This subroutine dumps the boundary distributions at node(-1) and node(Nxi+2)
  use SpecificTypes
  use ifport

  implicit none

  integer Nxi, Nspecies, iteration, neighbourLeft,neighbourRight, &
       attemptno, Nretries
  double precision XI(Nxi)
  type(species) particle(Nspecies)

  character filename*60, tmpfile*62, fnumber*7, snumber*2, form*20
  integer ixi, ivz, ii, imu, ierr

  attemptno = attemptno + 1

  if (neighbourLeft<0) then
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then

           ! Assemble filename
           write (fnumber,fmt='(i7.7)') iteration
           write (snumber,fmt='(i2.2)') ii
           tmpfile = 'outp/datfiles/fBC/XXfL' &
                //fnumber//'s'//snumber//'.ketchup.dat'
           filename = 'outp/datfiles/fBC/fL' &
                //fnumber//'s'//snumber//'.ketchup.dat'

           write (form,fmt='(a,i5.5,a)') '(', particle(ii)%Nmu, '(E13.6E3,a))'

           ! open file and write f(xi,vz)
           open(unit=1,file=tmpfile,status='replace',err=97)

           do ivz = 1, particle(ii)%Nvz
              write (1,fmt=form,err=98) (particle(ii)%node(-1)%f(ivz,imu), &
                   ' ', imu=1,particle(ii)%Nmu)
           end do

           close(1,err=99)
           ierr = rename(tmpfile,filename)

        end if
     end do
  end if

  if (neighbourRight<0) then
     do ii = 1, Nspecies
        if (.not. isnan(particle(ii)%mass)) then

           ! Assemble filename
           write (fnumber,fmt='(i7.7)') iteration
           write (snumber,fmt='(i2.2)') ii
           tmpfile = 'outp/datfiles/fBC/XXfR' &
                //fnumber//'s'//snumber//'.ketchup.dat'
           filename = 'outp/datfiles/fBC/fR' &
                //fnumber//'s'//snumber//'.ketchup.dat'

           write (form,fmt='(a,i5.5,a)') '(', particle(ii)%Nmu, '(E13.6E3,a))'

           ! open file and write f(xi,vz)
           open(unit=1,file=tmpfile,status='replace',err=97)

           do ivz = 1, particle(ii)%Nvz
              write (1,fmt=form,err=98) (particle(ii)%node(Nxi+2)%f(ivz,imu), &
                   ' ', imu=1,particle(ii)%Nmu)
           end do

           close(1,err=99)
           ierr = rename(tmpfile,filename)

        end if
     end do
  end if

  return

! Error handling section
97 write (*,*) 'DumpDistrBC: error in open statement, myid=', neighbourLeft+1
  goto 100
98 write (*,*) 'DumpDistrBC: error in write statement, myid=', neighbourLeft+1
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpDistrBC: error in close statement, myid=', neighbourLeft+1
100 if (attemptno<=Nretries) then
     call DumpDistrBC(Nxi,Nspecies,particle,iteration, &
          neighbourLeft,neighbourRight,attemptno,Nretries)
  end if

end subroutine DumpDistrBC

!------------------------------------------------------------------------

recursive subroutine DumpEfield(Nxi,fields_per_file,Nspecies, &
     iterationbuffer,fields_collected,Ebuffer,myid,attemptno,Nretries)

  use ifport

  implicit none

  integer Nxi, fields_per_file, Nspecies, iterationbuffer(fields_per_file), &
       fields_collected, myid, attemptno, Nretries
  double precision Ebuffer(Nxi,fields_per_file)

  character filename*69, tmpfile*62, fnumber*7, tnumber*7,  pnumber*4, form*20
  integer ifi, ixi, ierr

  attemptno = attemptno + 1

! Electric field
  ! Assemble filename
  write (fnumber,fmt='(i7.7)') iterationbuffer(1)
  write (tnumber,fmt='(i7.7)') iterationbuffer(fields_collected)
  write (pnumber,fmt='(i4.4)') myid
  tmpfile = 'outp/datfiles/Efield/XXEfield' &
       //'xxxxxxx'//'p'//pnumber//'.ketchup.dat'
  filename = 'outp/datfiles/Efield/Efield' &
       //fnumber//'To'//tnumber//'p'//pnumber//'.ketchup.dat'


  ! open file and write Ebuffer
  open(unit=1,file=tmpfile,status='replace',err=97)
  do ifi = 1, fields_collected
     write (1,*,err=98) 'timestep = ', iterationbuffer(ifi)
     write (1,*,err=98)
     do ixi = 1, Nxi
        write (1,fmt='(E13.6E3)',err=98) Ebuffer(ixi,ifi)
     end do
     write (1,*,err=98)
  end do
  close(1,err=99)
  ierr = rename(tmpfile,filename)

  return

! Error handling section
97 write (*,*) 'DumpEfield: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpEfield: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpEfield: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpEfield(Nxi,fields_per_file,Nspecies, &
          iterationbuffer,fields_collected,Ebuffer,myid,attemptno,Nretries)
  end if

end subroutine DumpEfield

!------------------------------------------------------------------------

recursive subroutine DumpDensity(Nxi,fields_per_file,Nspecies, &
     iterationbuffer,fields_collected,density1Dbuffer,B,BBC, &
     myid,attemptno,Nretries)
  use ifport

  implicit none

  integer Nxi, fields_per_file, Nspecies, iterationbuffer(fields_per_file), &
       fields_collected, myid, attemptno, Nretries
  double precision density1Dbuffer(Nxi,Nspecies,fields_per_file), B(Nxi), BBC(2)

  character filename*69, tmpfile*62, fnumber*7, tnumber*7,  pnumber*4, form*20
  integer ifi, ixi, ii, ierr
  double precision Bfactor(Nxi)

  attemptno = attemptno + 1
  
  Bfactor = B/BBC(1)

! Density
  ! Assemble filename
  write (fnumber,fmt='(i7.7)') iterationbuffer(1)
  write (tnumber,fmt='(i7.7)') iterationbuffer(fields_collected)
  write (pnumber,fmt='(i4.4)') myid
  tmpfile = 'outp/datfiles/density/XXdensity' &
       //'xxxxxxx'//'p'//pnumber//'.ketchup.dat'
  filename = 'outp/datfiles/density/density' &
       //fnumber//'To' //tnumber//'p'//pnumber//'.ketchup.dat'
  write (form,fmt='(a,i5.5,a)') '(', Nspecies, '(E13.6E3,a))'

  ! open file and write n(xi,species)
  open(unit=1,file=tmpfile,status='replace',err=97)
  do ifi = 1, fields_collected
     write (1,*,err=98) 'timestep = ', iterationbuffer(ifi)
     write (1,*,err=98)
     do ixi = 1, Nxi
        write (1,fmt=form,err=98) &
             (density1Dbuffer(ixi,ii,ifi)*Bfactor(ixi), ' ', ii=1,Nspecies)
     end do
     write (1,*,err=98)
  end do
  close(1,err=99)
  ierr = rename(tmpfile,filename)
  return

! Error handling section
97 write (*,*) 'DumpDensity: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpDensity: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpDensity: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpDensity(Nxi,fields_per_file,Nspecies,iterationbuffer, &
          fields_collected,density1Dbuffer,B,BBC,myid,attemptno,Nretries)
  end if

end subroutine DumpDensity

!------------------------------------------------------------------------

recursive subroutine DumpCurrent(Nxi,fields_per_file,Nspecies, &
     iterationbuffer,fields_collected,currentbuffer, &
     myid,attemptno,Nretries)

  use ifport

  implicit none

  integer Nxi, fields_per_file, Nspecies, iterationbuffer(fields_per_file), &
       fields_collected, myid, attemptno, Nretries
  double precision currentbuffer(Nxi,fields_per_file)

  character filename*69, tmpfile*62, fnumber*7, tnumber*7, pnumber*4, form*20
  integer ifi, ixi, ierr

  attemptno = attemptno + 1

! Current
  ! Assemble filename
  write (fnumber,fmt='(i7.7)') iterationbuffer(1)
  write (tnumber,fmt='(i7.7)') iterationbuffer(fields_collected)
  write (pnumber,fmt='(i4.4)') myid
  tmpfile = 'outp/datfiles/current/XXcurrent' &
       //'xxxxxxx'//'p'//pnumber//'.ketchup.dat'
  filename = 'outp/datfiles/current/current' &
       //fnumber//'To' //tnumber//'p'//pnumber//'.ketchup.dat'


  ! open file and write ibuffer
  open(unit=1,file=tmpfile,status='replace',err=97)
  do ifi = 1, fields_collected
     write (1,*,err=98) 'timestep = ', iterationbuffer(ifi)
     write (1,*,err=98)
     do ixi = 1, Nxi
        write (1,fmt='(E13.6E3)',err=98) currentbuffer(ixi,ifi)
     end do
     write (1,*,err=98)
  end do
  close(1,err=99)
  ierr = rename(tmpfile,filename)
  return

! Error handling section
97 write (*,*) 'DumpCurrent: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpCurrent: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpCurrent: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpCurrent(Nxi,fields_per_file,Nspecies, &
          iterationbuffer,fields_collected,currentbuffer, &
          myid,attemptno,Nretries)
  end if

end subroutine DumpCurrent

!------------------------------------------------------------------------

recursive subroutine DumpBfield(Nxi,B,dB,myid,attemptno,Nretries)

  use ifport

  implicit none

  integer Nxi, myid, attemptno, Nretries
  double precision B(Nxi), dB(Nxi)

  character filename*60, tmpfile*62, pnumber*4
  integer ixi, ierr

  attemptno = attemptno + 1

  ! Assemble filename
  write (pnumber,fmt='(i4.4)') myid
  tmpfile = 'outp/datfiles/Bfield/XXBfield' //'p'//pnumber//'.ketchup.dat'
  filename = 'outp/datfiles/Bfield/Bfield' //'p'//pnumber//'.ketchup.dat'

  ! open file and write f(vz)
  open(unit=1,file=tmpfile,status='replace',err=97)
  do ixi = 1, Nxi
     write (1,fmt='(E16.9E3,a,E16.9E3)',err=98) B(ixi), ' ', dB(ixi)
  end do
  close(1,err=99)
  ierr = rename(tmpfile,filename)
  return

! Error handling section
97 write (*,*) 'DumpBfield: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpBfield: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpBfield: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpBfield(Nxi,B,dB,myid,attemptno,Nretries)
  end if

end subroutine DumpBfield

!------------------------------------------------------------------------

recursive subroutine DumpGravity(Nxi,gravity,myid,attemptno,Nretries)

  use ifport

  implicit none

  integer Nxi, myid, attemptno, Nretries
  double precision gravity(Nxi)

  character filename*60, tmpfile*62, pnumber*4
  integer ixi, ierr

  attemptno = attemptno + 1

  ! Assemble filename
  write (pnumber,fmt='(i4.4)') myid
  tmpfile  = 'outp/datfiles/gravity/XXgravity' //'p'//pnumber//'.ketchup.dat'
  filename = 'outp/datfiles/gravity/gravity' //'p'//pnumber//'.ketchup.dat'

  ! open file and write f(vz)
  open(unit=1,file=tmpfile,status='replace',err=97)
  do ixi = 1, Nxi
     write (1,fmt='(E16.9E3)',err=98) gravity(ixi)
  end do
  close(1,err=99)
  ierr = rename(tmpfile,filename)
  return

! Error handling section
97 write (*,*) 'DumpGravity: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpGravity: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpGravity: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpGravity(Nxi,gravity,myid,attemptno,Nretries)
  end if

end subroutine DumpGravity

!------------------------------------------------------------------------

recursive subroutine Dumpg(Nxi,Z,gp,attemptno,Nretries)

  use ifport

  implicit none

  integer Nxi, attemptno, Nretries
  double precision Z(Nxi), gp(Nxi)

  character filename*60, tmpfile*62
  integer ixi, ierr

  attemptno = attemptno + 1

  ! Assemble filename
  tmpfile = 'outp/datfiles/g/XXg.ketchup.dat'
  filename = 'outp/datfiles/g/g.ketchup.dat'

  ! open file and write gp
  open(unit=1,file=tmpfile,status='replace',err=97)
  do ixi = 1, Nxi
     write (1,fmt='(E16.9E3,a,E16.9E3)',err=98) Z(ixi), ' ', gp(ixi)
  end do
  close(1,err=99)
  ierr = rename(tmpfile,filename)
  return


! Error handling section
97 write (*,*) 'Dumpg: error in open statement'
  goto 100
98 write (*,*) 'Dumpg: error in write statement'
  close(1,err=99)
  goto 100
99 write (*,*) 'Dumpg: error in close statement'
100 if (attemptno<=Nretries) then
     call Dumpg(Nxi,Z,gp,attemptno,Nretries)
  end if

end subroutine Dumpg

!------------------------------------------------------------------------

!!$subroutine DumpGamma(Nxi,Nspecies,particle,iteration, myid)
!!$
!!$  use SpecificTypes
!!$
!!$  implicit none
!!$
!!$  integer Nxi, Nspecies, iteration, myid
!!$  type(species) particle(Nspecies)
!!$
!!$  character filename*40, fnumber*7, snumber*2, pnumber*4, form*20
!!$  integer ixi, ivz, ii, imu
!!$
!!$  do ii = 1, Nspecies
!!$     if ((.not. isnan(particle(ii)%mass)) .and. &
!!$          particle(ii)%relativistic) then
!!$
!!$        ! Assemble filename
!!$        write(fnumber,fmt='(i7.7)') iteration
!!$        write(snumber,fmt='(i2.2)') ii
!!$        write(pnumber,fmt='(i4.4)') myid
!!$        filename='outp/gamma'// &
!!$             fnumber//'s'//snumber//'p'//pnumber//'.ketchup.dat' 
!!$
!!$        write (form,fmt='(a,i5.5,a)') '(', particle(ii)%Nmu, '(E13.6E3,a))'
!!$  
!!$        
!!$        ! open file and write gamma(xi,vz,mu)
!!$        ! Format is 
!!$        ! ixi = ixi (local)
!!$        ! blank line
!!$        ! gamma(vz,mu)
!!$        ! blank line
!!$        ! and then repeat for the next node, and the next, and ...
!!$        open(unit=1,file=filename,status='replace')
!!$        do ixi = 1, Nxi
!!$           write (1,*) 'ixi = ', ixi
!!$           write (1,*)
!!$           do ivz = 1, particle(ii)%Nvz
!!$              write (1,fmt=form)(particle(ii)%node(ixi)%gamma(ivz,imu), &
!!$                   ' ', imu=1,particle(ii)%Nmu)
!!$           end do
!!$           write (1,*)
!!$        end do
!!$        close(1)
!!$
!!$     end if
!!$  end do
!!$
!!$end subroutine DumpGamma

!------------------------------------------------------------------------
