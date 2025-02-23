program wolff
  implicit none

  integer, parameter :: L = 10
  integer, parameter :: iseed = 10
  integer, parameter :: Nth = 3000, N = L**2

  integer, dimension(2,L,L) :: S
  integer, dimension(4,L) :: nn
  integer :: i,j,k, ix, iy,i5,j5,j6
  real(8) :: Em,Mm,beta,T,E1,M1
  real :: r

  open(40,file='MxT.dat')
  open(50,file='ExT.dat')

  call random_init(.true.,.true.)

  do i = 1, L 
    nn(1,i) = 1 + mod(i,L)
    nn(2,i) = 1 + mod(i,L)
    nn(3,i) = L - mod(L-i+1,L)
    nn(4,i) = L - mod(L-i+1,L)
  enddo
  
  call startlattice(S,3)

  T = 0.1d0
  do while(T .le. 5.0d0)
    print*, "Estamos na Temperatura", T
    beta = 1.0d0/T
    Em = 0.0d0
    Mm = 0.0d0
    do i5 = 1, Nth
      do j5 = 1, L
        do j6 = 1, L
          call random_number(r)
          ix = L*r + 1
          call random_number(r)
          iy = L*r + 1
          call wolffmc(S,ix,iy,nn)
        enddo
      enddo
      call measure(S,E1,M1,i5,nn)
      Em = Em + E1
      Mm = Mm + M1
    enddo
    write(40,*)T,Mm/Nth
    write(50,*)T,Em/Nth
    T = T + 0.01d0
  enddo
contains

subroutine startlattice(S,istart)
  implicit none
  integer, dimension(2,L,L) :: S
  integer :: istart, ir

  if(istart .eq. 1) then
    do i = 1,L
      do j = 1,L
        S(1,i,j) = 1
        S(2,i,j) = 0
      end do
    enddo
    else if(istart .eq. 2) then
     do i = 1,L
      do j = 1,L
        S(1,i,j) = -1
        S(2,i,j) = 0
      end do
    enddo
  else
     do i = 1,L
      do j = 1,L
        call random_number(r)
        ir = 2*r
        S(1,i,j) = 1 - 2*ir
        S(2,i,j) = 0
      end do
    enddo
  endif
end subroutine startlattice

subroutine wolffmc(S,ix,iy,nn)
  implicit none

  integer, dimension(2,L,L) :: S
  integer, dimension(2,N) :: c
  integer :: icluster, ix, iy, icont, ica, dir
  integer, dimension(4,L) :: nn
  integer, dimension(2) :: ipos, ipv
  real(8) :: r1
  real :: r

  c = 0
  icont = 1
  ica = 1
  icluster = 0
  c(1,ica) = ix
  c(2,ica) = iy
  S(2,c(1,ica),c(2,ica)) = 1
  ipos = (/c(1,ica),c(2,ica)/)
  do i = 1, 4 
    dir = mod(i,2) + (1+(-1)**i)
    ipv = ipos
    ipv(dir) = nn(i,ipv(dir))
    call random_number(r1)
!   print*, 1-exp(dmin1(0.0d0,-2.0d0*beta*S(1,ipos(1),ipos(2))*S(1,ipv(1),ipv(2))))
    if(r1 .le. (1 - exp(dmin1(0.0d0,-2.0d0*beta*S(1,ipos(1),ipos(2))*S(1,ipv(1),ipv(2))))) ) then
      icont = icont + 1
      c(1,icont) = ipv(1)
      c(2,icont) = ipv(2)
      S(2,c(1,icont),c(2,icont)) = 1
      icluster = 1
    endif
  enddo
  
  do while(icluster .eq. 1)
    ica = ica + 1
    icluster = 0
    ipos = (/c(1,ica),c(2,ica)/)
    do i = 1, 4
      dir = mod(i,2) + (1+(-1)**i)
      ipv = ipos
      ipv(dir) = nn(i,ipv(dir))
      if(S(2,ipv(1),ipv(2)) .eq. 0) then
        call random_number(r1)
        if(r1 .le. (1 - exp(dmin1(0.0d0,-2.0d0*beta*S(1,ipos(1),ipos(2))*S(1,ipv(1),ipv(2))))) ) then
          icont = icont + 1
          c(1,icont) = ipv(1)
          c(2,icont) = ipv(2)
          S(2,c(1,icont),c(2,icont)) = 1
          icluster = 1
        endif
      endif
    enddo
  enddo

! print*, "Esse cluster teve", icont, "Spins"
  do i = 1, icont
    S(1,c(1,i),c(2,i)) = -S(1,c(1,i),c(2,i))
    S(2,c(1,i),c(2,i)) = 0
  enddo
end subroutine wolffmc

subroutine measure(S,E,M,i,nn)
  implicit none

  integer, dimension(2,L,L) :: S
  real(8) :: E,M,Em,Mm
  integer :: i, i1
  integer, dimension(4,L) :: nn

  E = 0.0d0
  M = 0.0d0
  do i1 = 1,L
    do j = 1,L
      M = M + S(1,i1,j)
      E = E - S(1,i1,j)*(S(1,i1,nn(2,j)) + S(1,i1,nn(4,j))+ S(1,nn(1,i1),j)+S(1,nn(3,i1),j) ) 
    enddo
  enddo
  M = M/N
  E = E/N
  Em = E
  Mm = M
! write(10,*)i, M
! write(20,*)i, E
end subroutine measure

end program wolff
