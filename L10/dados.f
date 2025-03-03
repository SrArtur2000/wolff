      program jackknife
        implicit real *8(a-h,o-z)
        parameter(idmax1 = 3000000, idmax2 = 10000000)
        dimension aMo(idmax1), aM1(idmax2), aM2(idmax2), aM4(idmax2)

        L = 10 !Tamanho da rede
        L2 = L**2

        open(10,file='Mdes.dat')
        open(20,file='susdes.dat')
        open(30,file='binderdes.dat')

        it = 10 !Tempo de correlação

        do k = 102,592
          aM1o = 0.0d0
          aMmean = 0.0d0
          na = 0 !tamanho do arquivo
          aT = real(k-100)/100 !Temperatura
          do i = 1, idmax1
            read(k,*,end=2)dummy,aMo(i)
            na = na + 1
            aMmean = aMmean + aMo(i)
          enddo

2         continue

          aMmean = aMmean/na

          nb = na/it
          aM1 = 0.0d0
          aM2s = 0.0d0
          aM2mean = 0.0d0
          aM4s = 0.0d0
          aM4mean = 0.0d0
          do i = 0, nb-1
           do j = 1, it
              aM1(i+1) = aM1(i+1) + aMo(j+it*i)
            enddo
            aM1(i+1) = aM1(i+1)/it
            aM2s = aM2s + aM1(i+1)**2
            aM2mean = aM2mean + aM1(i+1)
            aM4s = aM4s + aM1(i+1)**4
          enddo
        
          aM2s = aM2s/nb
          aM2mean = aM2mean/nb
          aM4s = aM4s/nb
  
          aM20 = L2*(aM2s - aM2mean**2)/aT
          aM40 = 1.0d0 - aM4s/(3.0d0*aM2s**2)

          sigmas1 = 0.0d0

          do i = 1, nb
            sigmas1 = sigmas1 + (aM1(i) - aMmean)**2
          enddo

          write(10,*) aT, aMmean, dsqrt(sigmas1/(nb-1))

          aM2 = 0.0d0
          aM4 = 0.0d0
          sigmam2 = 0.0d0
          sigmam4 = 0.0d0
          do i = 1, nb
            a3m = 0.0d0
            a3s = 0.0d0
            a3ss = 0.0d0
            do j = 1, nb
              if(i .ne. j) then
                a3m = a3m + aM1(i)
                a3s = a3s + aM1(i)**2
                a3ss = a3ss + aM1(i)**4
              endif
            enddo
            a3m = a3m/(nb-1)
            a3s = a3s/(nb-1)
            a3ss = a3ss/(nb-1)
            
            aM2(i) = L2*(a3s - a3m**2)/aT
            sigmam2 = sigmam2 + (aM2(i)-aM20)**2
            aM4(i) = 1.0d0 - a3ss/(3.0d0*a3s**2.0d0)
            sigmam4 = sigmam4 + (aM4(i)-aM40)**2
          enddo

          write(20,*)aT, aM20, dsqrt((nb-1)*sigmam2/nb)
          write(30,*)aT, aM40, dsqrt((nb-1)*sigmam4/nb)

        enddo

      end program
