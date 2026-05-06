
!   Modified driver for nafh2d.f code obtained from POTLIB:
!   References:  A. W. Jasper, M. D. Hack, A. Chakraborty, D. G. Truhlar, 
!                and P. Piecuch, J. Chem. Phys. 115, 7945 (2001).

! This dirver is for use with the Multi-State soft code, transforms the provided grid in Jacobi coordinates
! to the bond lengths used here. The outputs will then be read by the SOFT code to run the dynamics.
! The nafh2d.f surface file will need to be downloaded seperatley from POTLIB: https://comp.chem.umn.edu/potlib/form.cgi
! Compile this driver with the nafh2d.f file in a seperate folder and proved the path to the C++ SOFT code

program driver_nafh2d
  implicit none
  integer :: i, j, nt, nsurf
  integer :: nx, ny
  double precision :: rNaH_min, rNaH_max, rHF_min, rHF_max
  double precision :: rNaH, rHF, rNaF
  double precision :: drNaH, drHF
  double precision, allocatable :: r(:,:), e(:), de(:,:)

  ! ---- This common block exists in nafh2d.f ----
  integer g_myrank, g_nprocs
  common /com_para/ g_myrank, g_nprocs

  g_myrank = 0
  g_nprocs = 1

  call prepot

  ! -------- Grid settings (Bohr) --------
  ! (You can change these ranges)

  open(unit=111,file='input',status='old',action='read')
  read(111,*) rNaH_min,rNaH_max
  read(111,*) rHF_min,rHF_max
  read(111,*) nx
  read(111,*) ny
  
  
  !rNaH_min = 2.0d0
  !rNaH_max = 8.0d0
  !rHF_min  = 1.2d0
  !rHF_max  = 5.0d0

  !nx = 161
  !ny = 161

  drNaH = (rNaH_max - rNaH_min) / dble(nx - 1)
  drHF  = (rHF_max  - rHF_min ) / dble(ny - 1)

  nt = 1
  allocate(r(nt,3), e(nt), de(3,nt))

  do nsurf = 1, 3
     call write_surface(nsurf, nx, ny, rNaH_min, rHF_min, drNaH, drHF, r, e, de, nt)
  end do

  deallocate(r, e, de)

contains

  subroutine write_surface(nsurf, nx, ny, rNaH_min, rHF_min, drNaH, drHF, r, e, de, nt)
    implicit none
    integer, intent(in) :: nsurf, nx, ny, nt
    double precision, intent(in) :: rNaH_min, rHF_min, drNaH, drHF
    double precision, intent(inout) :: r(nt,3), e(nt), de(3,nt)
    integer :: i, j, unit
    double precision :: rNaH, rHF, rNaF, yJac, xJac, mass_F, mass_Na, mass_H
    character(len=32) :: fname

    if (nsurf == 1) fname = "U11.dat"
    if (nsurf == 2) fname = "U12.dat"
    if (nsurf == 3) fname = "U22.dat"

    mass_F = 1836.15 * 19
    mass_Na = 1836.15 * 29
    mass_H =1836.15 * 1

    unit = 10 + nsurf
    open(unit=unit, file=fname, status="replace", action="write")

    do j = 0, ny-1
       yJac = rHF_min + drHF * dble(j)
       do i = 0, nx-1
          xJac = rNaH_min + drNaH * dble(i)

          rNaF = xJac - yJac*mass_H/(mass_H + mass_F)
          rNaH = xJac + yJac*mass_F/(mass_H + mass_F)
          !rNaH = xJac


          ! Collinear Na-H-F cut:
          rHF = yJac 
          !rNaF = rNaH + rHF

          r(1,1) = rNaH
          r(1,2) = rHF
          r(1,3) = rNaF

          call pot(r, e, de, nt, nsurf)

          ! write: R_NaH  R_HF  Energy(Ha)
          write(unit,'(3(1x,ES20.10))') rNaH, rHF, e(1)
       end do
       write(unit,*)  ! blank line between rows (good for some plotting tools)
    end do

    close(unit)
    write(*,*) "Wrote ", trim(fname)
  end subroutine write_surface

end program driver_nafh2d
