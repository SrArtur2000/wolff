      program bootstrap
        implicit real *8(a-h,o-z)
        integer, parameter :: Nth=1000, Kth = 10000
        real(8) :: Oo(Nth),O1(Kth),O2(Kth),O4(Kth), sus(Kth), binder(Kth)


        open(10,file='mboot.dat')
        open(20,file='susboot.dat')
        open(30,file='binderboot.dat')

        L = 20
        L2 = L**2

        call random_init(.true.,.true.)
        
        do k = 102,592
          aT = real(k-100)/100
          Oo = 0.0d0
!Pega os dados que gerarão os dados aleatórios para o bootstrap
          do i = 1, Nth 
            read(k,*)dummy, Oo(i)
          enddo
          O1 = 0.0d0
          O2 = 0.0d0
          O4 = 0.0d0
          Okmean = 0.0d0
          Oksmean = 0.0d0
          OKssmean = 0.0d0
          do k1 = 1,Kth
            O1mean = 0.0d0
            Os = 0.0d0
            Oss = 0.0d0
            do k2 = 1,Nth
              call random_number(r)
              ip = Nth*r + 1
              O1mean = O1mean + Oo(ip)
              Os = Os + Oo(ip)**2
              Oss = Oss + Oo(ip)**4
            enddo
            O1(k1) = O1mean/Nth
            O2(k1) = Os/Nth
            O4(k1) = Oss/Nth
            Okmean = Okmean + O1(k1)
            Oksmean = Oksmean + O2(k1)
            Okssmean = Okssmean + O4(k1)
          enddo
          Okmean = Okmean/Kth
          Oksmean = Oksmean/Kth
          Okssmean = Okssmean/Kth

          sus = 0.0d0
          binder = 0.0d0
          susmean = 0.0d0
          bindermean = 0.0d0
          do k1 = 1,Kth
            sus(k1) = L2*(O2(k1) - Okmean**2)/aT
            susmean = susmean + sus(k1)
            binder(k1) = 1.0d0 - Okssmean/(3*Oksmean**2)
            bindermean = bindermean + binder(k1)
          enddo
          susmean = susmean/Kth
          bindermean = bindermean/Kth
          sigmao1 = 0.0d0
          sigmasus = 0.0d0
          sigmabi = 0.0d0
          do k1 = 1, Kth
            sigmao1 = sigmao1 + (O1(k1)-Okmean)**2
            sigmasus = sigmasus + (sus(k1) - susmean)**2
            sigmabi = sigmabi + (binder(k1) - bindermean)**2
          enddo
          sigmao1 = sigmao1/Kth
          sigmasus = sigmasus/Kth
          sigmabi = sigmabi/Kth
          write(10,*)aT, Okmean, dsqrt(sigmao1)
          write(20,*)aT, susmean, dsqrt(sigmasus)
          write(30,*)aT, bindermean, dsqrt(sigmabi)
        enddo

      end program
