program StructureFactor
    implicit none
    integer, parameter :: lx = 70, ly = 2
    real(8),parameter :: pi = 3.14159265359_8
    real(8), dimension(lx,ly) :: theta, phi
    real(8), dimension(2) :: q,r
    real(8) :: qr
    complex(8), dimension(0:lx-1,0:ly-1) :: Sq
    complex(8) :: S
    integer :: ix,iy,ix1,iy1,iqx,iqy

    do iy = 1,ly
        do ix = 1,lx
            read *,phi(ix,iy), theta(ix,iy)
        enddo
    enddo

    do iqx = 0, lx - 1
        q(1) = 2*pi*iqx/lx
        do iqy = 0, ly - 1
            q(2) = 2*pi*iqy/ly
            S = cmplx(0.0_8,0.0_8)
            do iy = 1,ly
                do ix = 1,lx
                    do iy1 = 1,ly
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
        enddo
        print '(3f12.6)',q(1),abs(Sq(iqx,0)),abs(Sq(iqx,1))
    enddo
    end program StructureFactor
