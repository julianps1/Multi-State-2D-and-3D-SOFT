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
  rNaH_min = 2.0d0
  rNaH_max = 8.0d0
  rHF_min  = 1.2d0
  rHF_max  = 5.0d0

  nx = 161
  ny = 161

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
    double precision :: rNaH, rHF, rNaF
    character(len=32) :: fname

    if (nsurf == 1) fname = "U11.dat"
    if (nsurf == 2) fname = "U12.dat"
    if (nsurf == 3) fname = "U22.dat"

    unit = 10 + nsurf
    open(unit=unit, file=fname, status="replace", action="write")

    do j = 0, ny-1
       rHF = rHF_min + drHF * dble(j)
       do i = 0, nx-1
          rNaH = rNaH_min + drNaH * dble(i)

          ! Collinear Na-H-F cut:
          rNaF = rNaH + rHF

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
