program main

    implicit none

    REAL(KIND(0.D0)) :: a,b
    character ( len = 255 ) :: filename
    integer :: n
    REAL(KIND(0.D0)), parameter :: pi = 3.14159265358979323846264338327950D+00
    REAL(KIND(0.D0)), allocatable, dimension ( : ) :: w, x

    REAL(KIND(0.D0)), PARAMETER :: psmall = 3.D-14
    REAL(KIND(0.D0)), PARAMETER :: pim4 = 0.7511255444649425D0 ! pi^(-1/4)
    REAL(KIND(0.D0)), parameter :: mypi = 3.14159265358979323846264338327950D+00
    INTEGER, PARAMETER :: mxit = 1000

    INTEGER :: i, j
    INTEGER :: its, m
    REAL(KIND(0.D0)) :: p1, p2, p3, pp
    REAL(KIND(0.D0)) :: z, z1
    REAL(KIND(0.D0)), DIMENSION(:), allocatable :: phi_tld

    filename='/Users/spencerbryngelson/Documents/Caltech/Research/analysis/qbmm/D/quadpts/'

    b = 0.2d0

    do n = 11,201,10
        print*, 'n = ', n
        allocate ( w(n) )
        allocate ( x(n) )
        allocate(phi_tld(n))

        w = 0d0; x = 0d0

        ! routine for Gauss-Hermite abscissas and weights (Numerical Recipe)
        ! Roots are symmetric about the origin, then find only half of them
        m = ( n+1 )/2
        DO i = 1,m

           IF ( i==1 ) THEN
              z = DSQRT( DBLE(2*n+1) ) - 1.85575D0*( DBLE(2*n+1) ) &
                **( -0.16667D0 )
           ELSE IF ( i==2 ) THEN
              z = z - 1.14D0*( DBLE(n) )**0.426D0/z
           ELSE IF ( i==3 ) THEN
              z = 1.86D0*z - 0.86D0*phi_tld(1)
           ELSE IF ( i==4 ) THEN
              z = 1.91D0*z - 0.91D0*phi_tld(2)
           ELSE
              z = 2.D0*z - phi_tld(i-2)
           END IF
           its = 1
           DO
              IF ( its>mxit.OR.DABS(z-z1)<=psmall ) THEN
                  IF (its >=2) THEN
                        print*, n,i,its, mxit, dabs(z-z1)
                        EXIT
                  END IF
              END IF
              p1 = pim4
              p2 = 0.D0
              DO j = 1,n
                 p3 = p2
                 p2 = p1
                 p1 = z*DSQRT( 2.D0/DBLE(j) )*p2 &
                    - DSQRT( DBLE(j-1)/DBLE(j) )*p3
              END DO
              pp = DSQRT( 2.D0*DBLE(n) )*p2
              z1 = z
              z = z1 - p1/pp
              its = its + 1
           END DO
           ! assign the root
           phi_tld(i) = z
           phi_tld(n+1-i) = -z
           ! assign the weight
           w(i) = 2.D0/( pp*pp )
           w(n+1-i) = w(i)
        END DO

        ! normalize the weights
        w = w/DSQRT( mypi )
        ! transform phi_tld into R0, R0(1) = R0mx, R0(n) = R0mn
        x = DEXP( DSQRT(2.D0)*b*phi_tld )

        ! sorting R0 s.t. R0(1) = R0mn, R0(n) = R0mx
        ! phi_tld is dummy
        DO i = 1,n
           phi_tld(n+1-i) = x(i)
        END DO
        x = phi_tld

        ! print*, 'x = ', x
        ! print*, 'w = ', w

        call rule_write(n,x,w,filename)

        deallocate(x,w,phi_tld)
    end do

    stop
end program main


subroutine rule_write ( n, x, w, filename )

    implicit none

    integer ( kind = 4 ) n
    character ( len = * ) filename
    character ( len = 255 ) filename_w
    character ( len = 255 ) filename_x
    real ( kind = 8 ) w(n)
    real ( kind = 8 ) x(n)
    character(len=5) :: x1 ! format descriptor
    integer :: i 

    write (x1,'(I5.5)') n ! converting integer to string using a 'internal file'
    filename_w = trim( filename ) // 'w_' // x1 // '.dat'
    filename_x = trim( filename ) // 'x_' // x1 // '.dat'

    open(unit=1,file=filename_w)
    do i = 1,n
        write(1,*) w(i)
    end do
    close(1)

    open(unit=1,file=filename_x)
    do i = 1,n
        write(1,*) x(i)
    end do
    close(1)

end subroutine
