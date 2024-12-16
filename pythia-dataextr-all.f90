		!------------------------------------------------------
		! Code by Dibakar Dhar (on 21-Jul-2024)
		! Description: ----------------------------------------
		! This program extracts eta and phi from pythia v8.3 
		! generated data remember to check all the check points 
		! before using this code. 		
		!------------------------------------------------------
		real(kind=8):: teta, eta(2000,1000), phi(2000,1000)
		real(kind=8):: mxeta,mneta,tphi,mxphi,mnphi
		integer i,j,c,tmul,mul(2000)
		integer pid,k,smul
		
		open(unit=1,file="z13tpp.dat")
		
		do i=1,330  !<- check the opuput file.
		read(1,*) 
		enddo ! skipping header lines
		c=1
!10		write(*,*) "Reading event:", c		
10		read(1,*,iostat=iosa) c , tmul ! total multiplicity
		if(iosa.eq.5010) goto 10
		if(iosa.ne.0) goto 15
		read(1,*)  ! skiping line 
		j=0
		do i=1,tmul-1 
		  read(1,*,iostat=iosb) pid, teta, tphi ! temp. eta
		  if(iosb.ne.0) goto 15
		  if(pid.eq.211) then  ! provide pid in the place of 211
		    j=j+1
		    eta(c,j)=teta 
		    phi(c,j)=tphi
		  endif
		enddo
		mul(c)=j ! multiplicity of particle per event
		c=c+1
		goto 10
15		close(1)	
		write(*,*)"==========================================="
		write(*,*)"No. of total event:", c-1
		write(*,*)"==========================================="

		open(unit=2,file="ex-13T-pp-211-all.dat")
		open(unit=3,file="ex-13T-pp-extra.dat")
		k=1			! counting needed events
		smul=0		! summation of multiplicity
		mxeta=0.	! cal. max. value of eta
		mneta=0.	! cal. min. value of eta
		mxphi=0.	! cal. max. value of phi
		mnphi=0.    ! cal. min. value of phi
		do i=1,c-1
		if(mul(i).ne.0) then
		  write(2,20) k, mul(i)
		  do j=1,mul(i)
		    write(2,25) j, eta(i,j), phi(i,j)
		    if(mxeta.lt.eta(i,j)) mxeta=eta(i,j)
		    if(mneta.gt.eta(i,j)) mneta=eta(i,j)
		    if(mxphi.lt.phi(i,j)) mxphi=phi(i,j)
		    if(mnphi.gt.phi(i,j)) mnphi=phi(i,j)   
		  enddo
		  k=k+1
		  smul=smul+mul(i)
		  ! write(*,*) k, smul, mul(i), mxeta, mneta  ! uncheck to print values on screen
		endif
		enddo
		
		write(*,*)"==========================================="
		write(*,*)"No. of needed events:", k-1
		write(*,*)"==========================================="
		
		write(3,27) "Total no. of event :", k-1 
		write(3,27) "Total multiplicity :", smul
		write(3,30) "Max. value of eta  :", mxeta
		write(3,30) "Min. value of eta  :", mneta
		write(3,30) "Max. value of phi  :", mxphi
		write(3,30) "Min. value of phi  :", mnphi
		
		write(*,*)"==========================================="
		write(*,*)"	Work Done !!! Thank you.		  "
		write(*,*)"==========================================="
		
20		format(i7,4x,i4)
25		format(i4,4X,f7.3,4X,f7.3)	
27		format(a,2X,i7)
30		format(a,2X,f7.3)	
		stop
		end
		
