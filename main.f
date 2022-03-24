! Program to read spline parameters from DFTB sk files
! and plot repulsive energy.  Units from input are in atomic units,
! output will give both hartree and eV energy units
        PROGRAM MAIN
        IMPLICIT NONE
        INTEGER     ::  interval            ! second line, nInt
        REAL        ::  cutoff              ! second line, cutoff
        REAL        ::  a_1, a_2, a_3       ! third line(exp part), a1, a2, a3
        REAL,DIMENSION(:),ALLOCATABLE :: dist_0, dist_end   ! start(r0), end
        REAL,DIMENSION(:),ALLOCATABLE :: c_0, c_1, c_2, c_3 ! coefficients
        REAL                          :: c_4, c_5           ! last line
        CHARACTER(LEN=30)    :: sp_file     ! file name
        INTEGER              :: f_in        ! file integer
        CHARACTER(LEN=3)     :: dummy
        INTEGER              :: iter, jter
        REAL                 :: repulse     ! repulsive energy

! Some variables(x-axis, interparticle distance) for output
        REAL        :: dist = 1.00544
        REAL        :: range_begin = 3.6,  stride
        INTEGER     :: step = 100, sub_step = 3

        CALL get_command_argument(1,sp_file)

! First line is empty, "Spline" text
        OPEN(NEWUNIT=f_in,FILE=sp_file,STATUS='old')
        READ(f_in,*)

! Second and third line
        READ(f_in,*) interval, cutoff
        READ(f_in,*) a_1, a_2, a_3

! Fourth line onward contain spline
        ALLOCATE(dist_0(interval))
        ALLOCATE(dist_end(interval))
        ALLOCATE(c_0(interval))
        ALLOCATE(c_1(interval))
        ALLOCATE(c_2(interval))
        ALLOCATE(c_3(interval))
        DO iter = 1, interval-1
          READ(f_in,*) dist_0(iter), dist_end(iter), 
     &      c_0(iter), c_1(iter), c_2(iter), c_3(iter)
        END DO

! Last line
          READ(f_in,*) dist_0(interval), dist_end(interval), 
     &      c_0(interval), c_1(interval), 
     &      c_2(interval), c_3(interval),
     &      c_4, c_5
        CLOSE(f_in)
! End of reading spline data file



               

!=============================================================================
! front part of repulsive energy (exponential part)
! equation 5 in manual
        stride = (dist_0(1) - range_begin) / REAL(step)
        DO iter = 0, step
          dist = range_begin + iter * stride
          CALL repulse_5(a_1, a_2, a_3, dist, repulse)
          PRINT *, dist, repulse
          !PRINT *, range_begin + iter * stride
        END DO

!=============================================================================
! middle part of repulsive energy (spline)
! equation 6 in manual
        DO iter = 1 , interval-1
            stride = (dist_end(iter) - dist_0(iter)) / REAL(sub_step)
            DO jter = 0, sub_step
              dist = dist_0(iter) + jter * stride
              CALL repulse_6(dist_0(iter), c_0(iter), c_1(iter),
     &          c_2(iter), c_3(iter),
     &          dist, repulse)
              PRINT *, dist, repulse
            END DO
        END DO


!=============================================================================
! later part of repulsive energy
! equation 7 in manual
        stride = (cutoff - dist_0(interval)) / REAL(step)

        DO iter = 0, step
          dist = dist_0(interval) + iter * stride
          CALL repulse_7(dist_0(interval), 
     &      c_0(interval), c_1(interval),c_2(interval),c_3(interval),
     &      c_4, c_5, dist, repulse)
          PRINT *, dist, repulse
        END DO


!=============================================================================
! stride checking
!        stride = (range_end - range_begin) / REAL(step)
!        PRINT *, "stride is", stride
!        DO iter = 1, step
!          dist = range_begin + iter * stride
!          PRINT *, dist
!          !PRINT *, range_begin + iter * stride
!        END DO
!        PRINT *, "end of step prints"
!=============================================================================

!=============================================================================
! read values checking
!        PRINT *, ""
!        PRINT *, "Checking read values"
!        PRINT *, ""
!
!        PRINT *, interval, cutoff
!        PRINT *, a_1, a_2, a_3
!
!
!        DO iter = 1, interval-1
!           PRINT*, dist_0(iter), dist_end(iter), 
!     &       c_0(iter), c_1(iter), c_2(iter), c_3(iter)
!        END DO
!
!        PRINT*, dist_0(interval), dist_end(interval), 
!     &    c_0(interval), c_1(interval), 
!     &    c_2(interval), c_3(interval),
!     &    c_4, c_5
!
!        PRINT *, "end of spline"
!=============================================================================

        END PROGRAM MAIN

        


! equation 5 in format manual
! exponential part
        SUBROUTINE repulse_5(a1, a2, a3, dist, repulse)
        IMPLICIT NONE
        REAL,INTENT(IN)      :: a1, a2, a3
        REAL,INTENT(IN)      :: dist
        REAL,INTENT(OUT)     :: repulse

        repulse = exp(-a1*dist + a2) 
        repulse = repulse + a3
        
        END SUBROUTINE repulse_5

! equation 6 in format manual
! coefficient part, (spline)
        SUBROUTINE repulse_6(dist_0, c0, c1, c2, c3, dist, repulse)
        IMPLICIT NONE
        REAL,INTENT(IN)      :: dist_0
        REAL,INTENT(IN)      :: c0, c1, c2, c3
        REAL,INTENT(IN)      :: dist
        REAL,INTENT(OUT)     :: repulse
        REAL                 :: diff
        
        diff = dist - dist_0
        repulse = c0 + c1*diff + c2*diff**2 + c3*diff**3

        END SUBROUTINE repulse_6


! equation 7 in format manual
! last line
        SUBROUTINE repulse_7(dist_0, c0, c1, c2, c3, c4, c5, dist 
     &          , repulse)
        IMPLICIT NONE
        REAL,INTENT(IN)      :: dist_0
        REAL,INTENT(IN)      :: c0, c1, c2, c3, c4, c5
        REAL,INTENT(IN)      :: dist
        REAL,INTENT(OUT)     :: repulse
        REAL                 :: diff

        diff = dist - dist_0
        repulse = c0 + c1*diff + c2*diff**2 + c3*diff**3 + c4*diff**4 
     &      + c5*diff**5

        END SUBROUTINE repulse_7
