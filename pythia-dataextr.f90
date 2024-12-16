		real(kind=8):: teta, eta(100000)
		integer i,j,c,tmul,mul(200000)
		integer pid,k
		
		open(unit=1,file="output.dat")
		
		do i=1,331
		read(1,*) 
		enddo ! skipping header lines
		c=1
10		write(*,*) "Reading event:", c		
		read(1,*,iostat=iosa) c , tmul ! total multiplicity
		if(iosa.ne.0) goto 15
		read(1,*)  ! skiping line 
		j=0
		do i=1,tmul-1 
		  read(1,*,iostat=iosb) pid, teta ! temp. eta
		  if(iosb.ne.0) goto 15
		  if(pid.eq.211) then
		    j=j+1
		    eta(j)=teta 
		  endif
		enddo
		mul(c)=j ! multiplicity of particle per event
		c=c+1
		goto 10
15		close(1)		
		write(*,*)"==========================================="
		write(*,*)"No. of total event:", c-1
		write(*,*)"==========================================="

		open(unit=2,file="ex211.dat")
		k=1
		do i=1,c-1
		if(mul(i).ne.0) then
		  write(2,20) k, mul(i)
		  do j=1,mul(i)
		    write(2,25) j, eta(j)
		  enddo
		  k=k+1
		endif
		enddo
20		format(i7,4x,i4)
25		format(i4,4X,f7.3)		
		stop
		end
		
