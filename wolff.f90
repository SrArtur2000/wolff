program wolff
  implicit none
  integer, parameter :: L = 40
  integer, parameter :: iseed = 10
  integer, parameter :: Nth = 2000, N = L**2

  integer, dimension(2) :: S
  integer :: i,j,k
  real :: r 


  call random_init(.true.,.true.)

  call startlattice(S,3)
  
contains

subroutine startlattice(S,istart)
  implicit none
  integer, dimension(2) :: S
  integer :: istart, ir

  if(istart .eq. 1) then
    do i = 1,L
      do j = 1,L
        S(i,j) = 1
      end do
    enddo
    else if(istart .eq. 2) then
     do i = 1,L
      do j = 1,L
        S(i,j) = -1
      end do
    enddo
  else
     do i = 1,L
      do j = 1,L
        call random_number(r)
        r = (r/2 - 1)*2
        S(i,j) = r
        print*, S(i,j)
      end do
    enddo
  endif
end subroutine startlattice

subroutine wolffmc(S)
  implicit none
  integer, dimension(2) :: S


end subroutine wolffmcc

subroutine measure(S)
  implicit none
  integer, dimension(2) :: S


end subroutine measure


end program wolff
