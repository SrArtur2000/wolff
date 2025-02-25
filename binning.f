	program binning
		implicit real *8(a-h,o-z)
		parameter(idmax = 1000000, L = 490)
		dimension O(idmax)
		
	       do L1 = 1,L
		    nO = 0
		    Omean = 0.0d0
		    
		    do i = 1, idmax
		      read(100 + L1,*,end=2)a, O(i)
		      Omean = Omean + O(i)
		      nO = nO + 1
		    enddo
2	          continue
		    
		    Omean = Omean/nO
		    
		    ib = 2
		    do while (ib .le. 10)
		      sigma = 0.0d0
			iO = 1
			do i = 1,nO-ib+1,ib
			  O2 = 0.0d0
			  do j = i, i+ib-1
			    O2 = O2 + O(j)
			  enddo
			  O2 = O2/ib
			  sigma = sigma + (O2 - Omean)**2
			  iO = iO + 1
		        enddo
		        
		      sigma = sigma/(iO-1)
		      write(1000 + L1,*)ib, dsqrt(sigma*ib/nO)
		      ib = ib+1
		    enddo
		enddo
	end program
