program autocorr
  implicit none
  integer, parameter :: Nth = 1000, Tf = 489
  integer :: i,j,k, tempo
  real(8), dimension(Nth) :: obs
  real(8) :: obsmean, x,autox

  do j = 2,Nth
    obsmean = 0.0d0
    obs = 0.0d0
    do i = 1, Nth
      read(j+100,*)tempo, obs(i)
      obsmean = obsmean + obs(i)
    enddo

    obsmean = obsmean/Nth

    x = 0.0d0
    do i = 1, Nth/2-1
      x = 0.0d0
      x = x + (autox(obs(1),obs(i+Nth/2),obsmean)+autox(obs(Nth/2),obs(i+Nth/2),obsmean))/2
      do k = 1, Nth-2
        x = x + (autox(obs(k),obs(k+i),obsmean)+autox(obs(k+1),obs(k+1+i),obsmean))/2
      enddo
      print*, x, obsmean
      write(2000+100+j,*)i,x/1000
    enddo
  enddo
end program


function autox(m,mt,mmean)
  implicit none

  real(8) :: autox, m, mt, mmean

  autox = m*mt - mmean**2

end function
