program binning
  implicit none
  
  integer, parameter :: Nth = 10000, tau = 100, Nnew = Nth/tau

  real(8), dimension(Nth) :: Obso
  real(8), dimension(Nnew) :: Obsnew
  real(8) :: Obsmean, Obs, T, dummy
  integer :: i, j, k,contador
  do k = 102, 592
    T = real(k-100)/100
    print*, "Estamos na temperatura", T
    Obso = 0.0d0
    Obsnew = 0.0d0
    contador = 0
    do i = 1, Nth-tau+1, tau
      contador = contador + 1
      Obsmean = 0.0d0
      do j = i, tau+i-1
        read(k,*)dummy, Obso(j)
        Obsmean = Obsmean + Obso(j)
      enddo
      Obsnew(contador) = Obsmean/tau
      write(1000+k,*)contador,Obsnew(contador)
    enddo
  enddo
end program
