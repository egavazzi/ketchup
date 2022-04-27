subroutine Bfield(Nz,Z,B,dB)

  implicit none

  integer Nz
  double precision Z(Nz), B(Nz), dB(Nz)

  double precision Bm, Bi, Lz

! initialisation
  B=1.0d-3
  dB=0.0d0

! Parameters
  Bm  = 0.0864d-6  ! Magnetospheric B at Z = 0
  Bi  = 56.0d-6    ! Ionospheric B at Z = Lz
  Lz  = 5.5d7      ! Length of the system [metres]


! Compute B
  B = Bm*exp((Z/Lz)**2.0d0 * (log(Bi/Bm) -0.6d0 -1.8d0*(Z/Lz)**2.0d0 + &
       2.4d0*(Z/Lz)**6.0d0))


! Compute dB/dz
  dB = B* (2.0d0*(Z/Lz)*(log(Bi/Bm)-0.6d0) - 7.2d0*(Z/Lz)**3.0d0 + &
       19.2d0*(Z/Lz)**7.0d0)/Lz


end subroutine Bfield

!------------------------------------------------------------------------

subroutine GetGravity(Nz,Z,gravity)

  implicit none

  integer Nz
  double precision Z(Nz), gravity(Nz)

  double precision L, G, Re, Me, TheFraction(Nz), theta_appr(Nz), r_appr(Nz)
  double precision a1, a2, a3

! Don't mourn, initialise!
  gravity = 0.0d0

  L=7.0d0         ! L-shell number
  Re=6371.2d3     ! Radius of our planet [m]
  Me=5.9736d24    ! Mass of our planet [kg]
  G=6.67428d-11   ! Newton's constant of gravitation

  a1 = 9.678850776859896d7
  a2 = 3.516930409763753d7
  a3 = 2.664304651374910d0

  theta_appr = sqrt(( ((Z-a1)/a2)**2.0d0 - 1.0d0)/a3)
  r_appr     = L*Re*(sin(theta_appr))**2.0d0


  ! The field aligned component of gravity is Fg * Br/B
  TheFraction = 2.0d0*cos(theta_appr) / &
       sqrt(4.0d0*cos(theta_appr)**2.0d0 + sin(theta_appr)**2.0d0)

  gravity=TheFraction*G*Me/r_appr**2.0d0


end subroutine GetGravity

!------------------------------------------------------------------------
