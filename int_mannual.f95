!=========================================================================
!       PROG FOR INTERMITTENCY IN  PHASE SPACE (COS) IN 1-DIM.
!       This program is to calculate intermittency in 1 dimension
!		This the edited copy developed by Dibakar Dhar (**A mannual copy)
!		This program includes all manuual for all 1d intmc 
!		Please check and if change any thing please make a copy before.
!=========================================================================
! ## This program is made made to read the data in the following format ##
! ## ---------------------------------------------------------------------
! ##   event_no		multiplicity
! ##   id			rapidity
!=========================================================================
         integer(kind=8):: NS(80000),K(80000,1000)
         real(kind=8):: N(80000),PRK(80000,1000),PUF(1000),AK(1000)
         real(kind=8):: FM2(80000),FM3(80000),FM4(80000),FM5(80000)
         real(kind=8):: FM6(80000),FM7(80000),FM8(80000)
         REAL(kind=8):: MULT2,MULT3,MULT4,MULT5,MULT6,MULT7,MULT8
         integer(kind=8):: i,JCD,NEV
!=========================================================================
! ** remember the dimensions(80000) of NS,N,FM2,FM3,FM4,FM5,FM6,FM7,
! ** FM8 will be a little bit more than the total no. of events. 
! ** And other dimensions (i.e. 2000) are a little bit more
! ** than the max. mul.          
!=========================================================================         
         OPEN(2,FILE='p13tm211.dat') !<= input here the input file name
         OPEN(4,FILE='O30pQ2.dat')
         OPEN(7,FILE='O30pQ3.dat')
         OPEN(8,FILE='O30pQ4.dat')
         OPEN(9,FILE='O30pQ5.dat')
         OPEN(10,FILE='O30pQ6.dat')
         OPEN(11,FILE='O30pQ7.dat')
         open(12,file='O30pQ8.dat')

!=========================================================================
!		provide the eta range  below
!=========================================================================
          RMAX=6
          RMIN=-6.0
!=========================================================================
		NEV=76423	! <= provide total no. of events
		do j=1,NEV
			KCD=0
			READ(2,*) temp, JCD	! reading no. of particles
			NS(j)=JCD
			write(3,*) JCD
			do i=1,JCD
				READ(2,*) temp, PUF(i) ! reading eta
				KCD = KCD +1
				PRK(j,KCD)=PUF(i)
				write(3,*) PRK(j,KCD)
			enddo
			write(*,*) "Event",j,"Done", JCD, "and", KCD
		enddo
!========================================================================
!####### CALCULATION PART (** No Changes Should be made after)	#########	
!========================================================================
      WRITE(*,*)"No. of Event: ", NEV , "kill 2"
      ENEV = NEV
      WRITE(4,81)
81    FORMAT(14X,'PHI',13X,'AVF2',10X,'ERRLNF2')
      WRITE(7,82)
82    FORMAT(14X,'PHI',13X,'AVF3',10X,'ERRLNF3') 
       WRITE(8,83)
83    FORMAT(14X,'PHI',13X,'AVF4',10X,'ERRLNF4') 
       WRITE(9,84)
84    FORMAT(14X,'PHI',13X,'AVF5',10X,'ERRLNF5') 
       WRITE(10,85)
85    FORMAT(14X,'PHI',13X,'AVF6',10X,'ERRLNF6') 
       WRITE(11,86)
86    FORMAT(14X,'PHI',13X,'AVF7',10X,'ERRLNF7') 
       WRITE(12,87)
87    FORMAT(14X,'PHI',13X,'AVF8',10X,'ERRLNF8') 
      DO 101 MN=4,40
      BINN = (RMAX - RMIN)/FLOAT(MN)
      
      DO 120 II=1,NEV
      RL = RMIN
      RR = RL + BINN
      DO 110 I = 1,MN
      K(II,I) = 0
      DO 20 L = 1,NS(II)
      IF((PRK(II,L).GE.RL).AND.(PRK(II,L).LT.RR)) then 
      	K(II,I)=K(II,I)+1.
      endif

20    CONTINUE
      RL = RR
      RR = RL +BINN
110    CONTINUE
120   CONTINUE
      M=MN
      NN=0
      DO 200 I=1,M
      AK(I)=0.
      DO 210 II=1,NEV
      AK(I)=AK(I)+K(II,I)
      NN=NN+K(II,I)
210   CONTINUE
      AK(I)=AK(I)/FLOAT(NEV)
200   CONTINUE

	
      AVN = FLOAT(NN)/FLOAT(NEV)
      write(*,*) NN, "kill ", I, AVN
      SUM2=0.
      SUM3=0.
      SUM4=0.
      SUM5=0.
      SUM6=0.
      SUM7=0.
      SUM8=0.
      
      R2=0
      R3=0
      R4=0
      R5=0
      R6=0
      R7=0
      R8=0
      
      do 220 i=1,m
      SUM2=SUM2+(AK(I)/AVN)**2
      SUM3=SUM3+(AK(I)/AVN)**3
      SUM4=SUM4+(AK(I)/AVN)**4
      SUM5=SUM5+(AK(I)/AVN)**5
      SUM6=SUM6+(AK(I)/AVN)**6
      SUM7=SUM7+(AK(I)/AVN)**7
      SUM8=SUM8+(AK(I)/AVN)**8
220   CONTINUE

      R2=SUM2*M
      R3=SUM3*(M**2)
      R4=SUM4*(M**3)
      R5=SUM5*(M**4)
      R6=SUM6*(M**5)
      R7=SUM7*(M**6)
      R8= SUM8*(M**7)
      
      AVF2 = 0.
      AVF3 = 0.
      avf4=0.
      AVF5=0.
      AVF6=0.
      AVF7=0.
      AVF8=0.
      
      DO 27 II = 1,NEV
         FM2(II)=0
        FM3(II)=0
        fm4(ii)=0
        FM5(II)=0
        FM6(II)=0
        FM7(II)=0
        FM8(II)=0
        
       F2=0.
      F3=0.
      f4=0.
      F5=0.
      F6=0.
      F7=0.
      F8=0.
      
       F22=0
      F33=0
      f44=0
      F55=0
      F66=0
      F77=0
      F88=0
      
       DO 29 I = 1,M
       A=FLOAT(K(II,I))*(FLOAT(K(II,I))-1.)
       
       MULT2 = A/(AVN**2.)
       
      MULT3 = (MULT2*(FLOAT(K(II,I))-2))/AVN
      MULT4= (MULT3*(FLOAT(K(II,I))-3))/AVN
      MULT5=(MULT4*(FLOAT(K(II,I))-4))/AVN
      MULT6=(MULT5*(FLOAT(K(II,I))-5))/AVN
      MULT7=(MULT6*(FLOAT(K(II,I))-6))/AVN
      MULT8=(MULT7*(FLOAT(K(II,I))-7))/AVN
      
      
       F2 = F2 + MULT2
      F3 = F3 + MULT3
      F4= F4+ MULT4
      F5= F5+MULT5
      F6=F6+MULT6
      F7=F7+MULT7
      F8=F8+MULT8
      
29     CONTINUE
       M=MN
       AVF2=AVF2+F2*(M)
       F22=F2*M
      AVF3=AVF3+F3*(M**2)
      F33=F3*(M**2)
      AVF4=AVF4+F4*(M**3)
      F44= F4*(M**3)
      AVF5=AVF5+F5*(M**4)
      F55=F5*(M**4)
      AVF6=AVF6+F6*(M**5)
      F66=F6*(M**5)
      AVF7=AVF7+F7*(M**6)
      F77=F7*(M**6)
      AVF8=AVF8+F8*(M**7)
      F88=F8*(M**7)

      FM2(II)=F22
      FM3(II)=F33
      FM4(II)=F44
      FM5(II)=F55
      FM6(II)=F66
      FM7(II)=F77
      FM8(II)=F88

27     CONTINUE
       AVF2=AVF2/(FLOAT(NEV))
       AVGF2=AVF2
       IF(AVF2.LE.0.) GO TO 501
       AVF2=ALOG(AVF2)
501    AVF3=AVF3/(FLOAT(NEV))
       AVGF3=AVF3
       IF(AVF3.LE.0.) GO TO 502
       AVF3=ALOG(AVF3)
502    AVF4=AVF4/(FLOAT(NEV))
       AVGF4=AVF4
       IF(AVF4.LE.0.)GO TO 503
       AVF4=ALOG(AVF4)
503    AVF5=AVF5/(FLOAT(NEV))
       AVGF5=AVF5
       IF(AVF5.LE.0.)GO TO 504
       AVF5=ALOG(AVF5)
504    AVF6=AVF6/(FLOAT(NEV))
       AVGF6=AVF6
       IF(AVF6.LE.0.)GO TO 505
       AVF6=ALOG(AVF6)
505    AVF7=AVF7/(FLOAT(NEV))
       AVGF7=AVF7
       IF(AVF7.LE.0.)GO TO 5055
       AVF7=ALOG(AVF7)

5055   AVF8=AVF8/(FLOAT(NEV))
       AVGF8=AVF8
       IF(AVF8.LE.0.)GO TO 506
       AVF8=ALOG(AVF8)


506   ETA=LOG(float(MN))
      SIG2=0
      SIG3=0
      SIG4=0
      SIG5=0
      SIG6=0
      SIG7=0
      SIG8=0
       DO 15 II=1,NEV
      SIG2=SIG2+(AVGF2-FM2(II))**2
      SIG3=SIG3+(AVGF3-FM3(II))**2
      SIG4=SIG4+(AVGF4-FM4(II))**2
      SIG5=SIG5+(AVGF5-FM5(II))**2
      SIG6=SIG6+(AVGF6-FM6(II))**2
      SIG7=SIG7+(AVGF7-FM7(II))**2
      SIG8=SIG8+(AVGF8-FM8(II))**2

15     CONTINUE
      SIG2=SIG2/NEV
      SIG3=SIG3/NEV
      SIG4=SIG4/NEV
      SIG5=SIG5/NEV
      SIG6=SIG6/NEV
      SIG7=SIG7/NEV
      SIG8=SIG8/NEV

      SIGMA2=SQRT(SIG2)
      SIGMA3=SQRT(SIG3)
      SIGMA4=SQRT(SIG4)
      SIGMA5=SQRT(SIG5)
      SIGMA6=SQRT(SIG6)
      SIGMA7=SQRT(SIG7)
      SIGMA8=SQRT(SIG8)

      ERRF2=SIGMA2/SQRT(ENEV)
      ERRF3=SIGMA3/SQRT(ENEV)
      ERRF4=SIGMA4/SQRT(ENEV)
      ERRF5=SIGMA5/SQRT(ENEV)
      ERRF6=SIGMA6/SQRT(ENEV)
      ERRF7=SIGMA7/SQRT(ENEV)
      ERRF8=SIGMA8/SQRT(ENEV)

       IF (AVGF2.EQ.0) GO TO 41
       GO TO 51
51     ERRLNF2=ERRF2/AVGF2
61    IF (AVGF3.EQ.0) GO TO 42
      GO TO 52
52    ERRLNF3=ERRF3/AVGF3
62    IF (AVGF4.EQ.0.) GO TO 43
      GO TO 53
53    ERRLNF4=ERRF4/AVGF4 
63      IF (AVGF5.EQ.0) GO TO 44
       GO TO 55
55     ERRLNF5=ERRF5/AVGF5
64       IF (AVGF6.EQ.0) GO TO 45 
       GO TO 56
56     ERRLNF6=ERRF6/AVGF6
65      IF (AVGF7.EQ.0) GO TO 46 
       GO TO 57
57     ERRLNF7=ERRF7/AVGF7
      GO TO 655

655    IF (AVGF8.EQ.0) GO TO 47 
       GO TO 577
577     ERRLNF8=ERRF8/AVGF8
        GO TO 66


41     ERRLNF2=0

      GO TO 61
42     ERRLNF3=0
       GO TO 62
43     ERRLNF4=0
       GO TO 63
44    ERRLNF5=0
      GO TO 64
45     ERRLNF6=0
       GO TO 65
46    ERRLNF7=0
      GO TO 655
47    ERRLNF8=0
      GO TO 66


66     WRITE(4,105)ETA,AVF2,ERRLNF2
       WRITE(7,105)ETA,AVF3,ERRLNF3
       WRITE(8,105)ETA,AVF4,ERRLNF4
       WRITE(9,105)ETA,AVF5,ERRLNF5
       WRITE(10,105)ETA,AVF6,ERRLNF6
       WRITE(11,105)ETA,AVF7,ERRLNF7
       WRITE(12,105)ETA,AVF8,ERRLNF8

105    FORMAT(10X,F7.3,10X,F8.3,10X,F7.3)
101    CONTINUE
	   write(*,*) "=================================="
	   write(*,*) "Work Done !!! Thank you..."
	   write(*,*) "=================================="
       STOP
       END


