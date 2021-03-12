program StructureFactor
    implicit none
    integer :: lx, ly
    real(8),parameter :: pi = 3.14159265359_8
    real(8), allocatable, dimension(:,:) :: theta, phi
    real(8), dimension(2) :: q,r
    real(8) :: qr
    !complex(8), dimension(0:lx-1,0:ly-1) :: Sq
    complex(8), allocatable, dimension(:,:) :: Sq
    complex(8) :: S
    integer :: ix,iy,ix1,iy1,iqx,iqy
    real(8) :: t, delta, alphay, alphaz, Bx, Bz, cp, temperature, J_spin, gamma_N, disorder ! all these variables are ignored
    logical :: magnetic_disorder
    namelist /parameters/ lx, ly, t, delta, alphay, alphaz, Bx, Bz, J_spin, gamma_N, temperature, cp, disorder, magnetic_disorder
    open(9,file='stripe.inp',status='OLD')
    read(9,nml=parameters)
    close(9)
    allocate(theta(lx,ly), phi(lx, ly), Sq(-lx/2:lx/2,-ly/2:ly/2))

    do iy = 1,ly
        do ix = 1,lx
            read *,phi(ix,iy), theta(ix,iy)
        enddo
    enddo

    do iqx = -lx/2, lx/2  ! Brilluoin zone from -pi to pi
        q(1) = 2*pi*iqx/lx
        do iqy = -lx/2, ly/2
            q(2) = 2*pi*iqy/ly
            S = cmplx(0.0_8,0.0_8)
            do iy = 1,ly
                do ix = 2,lx-1 ! edge sites neglected
                    do iy1 = 2,ly-1
                        do ix1 = 1,lx
                            r(1) = ix1 - ix
                            r(2) = iy1 - iy
                            qr = dot_product(q,r)
                            S = S + cmplx(cos(qr),sin(qr))*(sin(theta(ix,iy))*sin(theta(ix1,iy1)) &
                                   *cos(phi(ix,iy)-phi(ix1,iy1)) + cos(theta(ix,iy))*cos(theta(ix1,iy1)))
                        enddo
                    enddo
                enddo
            enddo
            Sq(iqx,iqy) = S/(lx*ly)
            print '(2f8.4,f12.6)',q,abs(Sq(iqx,iqy))
        enddo
        print *
    enddo
    end program StructureFactor
