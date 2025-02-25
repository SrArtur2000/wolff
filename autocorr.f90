program autocorr
  implicit none
  integer, parameter :: Nth = 1000, Tf = 489
  integer :: i,j,k
  integer, dimension(Nth) :: tempos
  real(8), dimension(Nth) :: obs
  real(8) :: obsmean, x

  do j = 100,Tf
    obsmean = 0.0d0
    obs = 0.0d0
    do i = 1, Nth
      read(j+100,*)tempos(i), obs(i)
      obsmean = obsmean + obs(i)
    enddo

    obsmean = obsmean/Nth
    print*, obsmean
    x = 0.0d0
    do i = 1, Nth
      x = 0.0d0
      do k = 1, Nth
        x = x + obs(k)*obs(k+i)-obsmean**2
      enddo
      print*, x
      write(2000+j,*)i,x
    enddo
  enddo
end program
