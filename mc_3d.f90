! obliczenia MC 
!
module stale
  integer :: n, lx, ly
  real(8) :: t, delta, alphay, alphaz, Bx, Bz, cp, temperature, J_spin, gamma_N, phi, & ! poszerzenie delt Diraca
             disorder, & ! amplituda nieporzadku 
             pi=3.141592653589793_8 

 ! real(8),parameter::phi=10.0d0,eps=1.e-3 
  real(8), parameter:: kT=1.
  real(8), parameter::eps=1.e-3 
  
  logical :: magnetic_disorder ! .false. dla ladunkowego, .true. dla magnetycznego
end module stale
  
program majorana
  use stale
  use lapack95
  implicit none
  interface
     subroutine hamiltonian(h,ang_xy,ang_z)
       complex(8), intent(out) :: h(:,:)
       real(8), intent(in) :: ang_xy(:), ang_z(:)
     end subroutine hamiltonian
     subroutine green(h,g,omega)
       complex(8),intent(in) :: h(:,:)
       complex(8),intent(out) :: g(:,:)
       real(8) :: omega
     end subroutine green
     function topo(ang_xy_tab,ang_z_tab)
       real(8) :: topo, ang_xy_tab(:), ang_z_tab(:)
     end function topo
  end interface
  real(8), allocatable :: e(:), ang_xy(:), ang_z(:), rnd(:)
  real(8) :: free_energy, r, free_en, old_free_en, en, old_en, c_z, ang_xy_old, ang_z_old
  complex(8), allocatable :: g(:,:), h(:,:), h1(:,:)
  complex(8) :: h_old(8)
  integer :: i, j, info, in, k, l, i_therm
  logical :: old_config
  namelist /parameters/ lx, ly, t, delta, alphay, alphaz, Bx, Bz, gamma_N, temperature, cp, disorder, magnetic_disorder, J_spin

  call system("echo $PPID >MC_PID")
  open(9,file='stripe.inp',status='OLD')
  read(9,nml=parameters)
  close(9)

  n = lx*ly

  allocate(e(4*n), g(4*n,4*n), h(4*n,4*n), h1(4*n,4*n), ang_xy(n), ang_z(n), rnd(n))
  i_therm = 0
  inquire(FILE="spins_3d.dat", EXIST=old_config)
  if (old_config) then
     open(10,file="spins_3d.dat")
     do i=1,n
        read(10,*) ang_xy(i), ang_z(i)
     enddo
  else
     call random_number(rnd)
     ang_xy=2*pi*rnd
     call random_number(rnd)
     ang_z=acos(2*rnd - 1)
     !ang_z=0.5*pi
  endif

  call hamiltonian(h,ang_xy,ang_z)
  h1=h
  call heevd(h1, e, info=info)
  i = 0        
  old_free_en = free_energy(e)
  old_en = sum(e(1:2*n))
  print *,old_free_en   !,topo(ang_xy,ang_z)

  main_loop: do !###############################################################
     call random_number(r)
     in = int(n*r) + 1
     ang_xy_old = ang_xy(in)
     ang_z_old = ang_z(in)
     call random_number(r)
     ang_xy(in) = 2*pi*r  
     call random_number(r)
     ang_z(in) = acos(2*r-1)
     !ang_z(in)=0.5*pi
     h_old(1) = h(in,in+n)
     h_old(2) = h(in+n,in)
     h_old(3) = h(in+2*n,in+3*n)
     h_old(4) = h(in+3*n,in+2*n)
     h_old(5) = h(in,in)
     h_old(6) = h(in+n,in+n)
     h_old(7) = h(in+2*n,in+2*n)
     h_old(8) = h(in+3*n,in+3*n)

     h(in,in+n) = 0.5d0*J_spin*sin(ang_z(in))*cmplx(cos(ang_xy(in)),-sin(ang_xy(in))) ! S.s
     h(in+n,in) = conjg(h(in,in+n))
     h(in+2*n,in+3*n) = h(in,in+n)
     h(in+3*n,in+2*n) = conjg(h(in+2*n,in+3*n))
     c_z=cos(ang_z(in))
     h(in,in)=Bz-cp+0.5d0*J_spin*c_z
     h(in+n,in+n)=-Bz-cp-0.5d0*J_spin*c_z
     h(in+2*n,in+2*n)=Bz+cp+0.5d0*J_spin*c_z
     h(in+3*n,in+3*n)=-Bz+cp-0.5d0*J_spin*c_z

     h1=h
     call heevd(h1, e, info=info)
     free_en = free_energy(e)
     en = sum(e(1:2*n))
     call random_number(r)
     if ((free_en < old_free_en) .or. (old_free_en-free_en > temperature*log(r))) then
        old_free_en = free_en
        old_en = en
     else
        ang_xy(in) = ang_xy_old
        ang_z(in) = ang_z_old

        h(in,in+n) = h_old(1)
        h(in+n,in) = h_old(2)
        h(in+2*n,in+3*n) = h_old(3)
        h(in+3*n,in+2*n) = h_old(4)
        h(in,in) = h_old(5)
        h(in+n,in+n) = h_old(6)
        h(in+2*n,in+2*n) = h_old(7)
        h(in+3*n,in+3*n) = h_old(8)

        free_en = old_free_en
        en = old_en
     endif
     i = i + 1
     if (i > i_therm .and. mod(i,100)==0) then
        open(10,file='spins_3d.dat')
        do l = 1, n
           write(10,*) ang_xy(l), ang_z(l)
        enddo
        close(10)
        print *,free_en  !,topo(ang_xy,ang_z)
     endif
  enddo main_loop !#############################################################
end program majorana

subroutine hamiltonian(h,ang_xy,ang_z)
  ! oblicza elementy hamiltonianu
  use stale
  implicit none
  complex(8), intent(out) :: h(:,:)
  real(8), intent(in):: ang_xy(:),ang_z(:)
  real(8) :: c_z
  integer :: i,j,ix,iy,idxy
  h=0._8

  do ix = 1, lx 
     do iy = 1, ly 
        i = (iy - 1) * lx + ix
        do idxy = 0, 1                 ! idxy = 0 -> przeskok w kierunku +x, idxy = 1 -> przeskok w +y
           if (idxy == 0) then
              if (ix == lx) cycle
              j = i + 1
           else
              if (iy == ly) cycle
              j = i + lx
           endif

           h(i,j) = cmplx(t,alphay)            !  C^+_up C_up
           h(i+n,j+n) = cmplx(t,-alphay)       !  C^+_dn C_dn
           h(i+2*n,j+2*n) = cmplx(-t,alphay)   ! -C^+_dn C_dn
           h(i+3*n,j+3*n) = cmplx(-t,-alphay)  ! -C^+_up C_up
           !
           h(i,j+n) = alphaz                   !  C^+_up C_dn
           h(j,i+n) = -alphaz                  !  C^+_up C_dn
           h(i+2*n,j+3*n) = -alphaz            !  C^+_up C_dn
           h(j+2*n,i+3*n) = alphaz             !  C^+_up C_dn
           !
        enddo
        c_z=cos(ang_z(i))
        h(i,i) = Bz-cp+0.5d0*J_spin*c_z
        h(i+n,i+n) = -Bz-cp-0.5d0*J_spin*c_z
        h(i+2*n,i+2*n) = Bz+cp+0.5d0*J_spin*c_z
        h(i+3*n,i+3*n) = -Bz+cp-0.5d0*J_spin*c_z
        !
        h(i,i+n) = Bx
        h(i+2*n,i+3*n) = Bx
        !
        h(i,i+2*n) = delta
        h(i+n,i+3*n) = delta
        !
        h(i,i+n) = 0.5d0*J_spin*sin(ang_z(i))*cmplx(cos(ang_xy(i)),-sin(ang_xy(i))) ! S.s
        h(i+2*n,i+3*n) = h(i,i+n)
     enddo
  enddo
  
  do i=1,4*n-1  ! elementy pod przekatna hamiltonianu
     do j=i+1,4*n
        h(j,i)=conjg(h(i,j))
     enddo
  enddo
end subroutine hamiltonian
!
!################################################################################
!
subroutine green(h,g,omega)
 ! oblicza funkcje Greena G_{ij}(omega)=[(omega+i*eps-H)^{-1}]_{ij}
  use lapack95
  use stale
  implicit none
  complex(8),intent(in) :: h(:,:)
  complex(8),intent(out) :: g(:,:)
  real(8) :: omega
  integer :: i, info
  integer, allocatable :: ipiv(:)
  allocate(ipiv(4*n))
  g=-h
  do i=1,4*n
     g(i,i)=g(i,i)+cmplx(omega,0.5*gamma_N)
  enddo
  !------------ odwracanie macierzy g -----------
  call getrf(g,ipiv,info)
  if (info /= 0) then
     print *,'GETRF',info
     stop 'GETRF'
  endif
  call getri(g,ipiv,info)
  if (info /= 0) then
     print *,'GETRI',info
     stop 'GETRI'
  endif
end subroutine green
!
!#######################################################################
!
function free_energy(ens)
use stale
implicit none
real(8),dimension(4*n) :: ens
real(8) :: free_energy,expon
integer :: i
free_energy=0.
do i=1,4*n
   !  expon=-(ens(i)-cp)/temperature
   expon=-ens(i)/temperature
   if (expon > 20) then
      free_energy=free_energy+expon
   else if (expon > -20) then
      free_energy=free_energy+log(1._8+exp(expon))
   endif
enddo
free_energy=-temperature*free_energy
end function free_energy
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
function topo(ang_xy_tab, ang_z_tab)
  use stale
  use lapack95
  implicit none

  interface
     function ainv(a)
       complex(8),dimension(4,4) :: ainv,a
     end function ainv
     function cdot(w2,w1)
       complex(8),dimension(8,8) :: cdot,w2,w1
     end function cdot
  end interface

  complex(8) :: tt(4,4),hh(4,4),uu(8,8),mbar(8,8),ww(8,8),ww1(8,8),mm(8,8),det
  real(8) :: ang_xy, ang_z, topo, s_xy, c_xy, s_z, c_z
  real(8) :: ang_xy_tab(:), ang_z_tab(:)
  integer :: i

  uu = 0
  do i=1,4
     uu(i,i) = cmplx(1,0)
     uu(i,i+4) = cmplx(1,0)
     uu(i+4,i) = cmplx(0,1)
     uu(i+4,i+4) = cmplx(0,-1)
  enddo
  uu = uu/sqrt(2._8)

  do i=1,n-1
     ang_xy = ang_xy_tab(i)
     ang_z = ang_z_tab(i)
     s_xy = sin(ang_xy_tab(i))
     c_xy = cos(ang_xy_tab(i))
     s_z = sin(ang_z_tab(i))
     c_z = cos(ang_z_tab(i))
     tt = cmplx(0._8,0._8)
!     if (i < n) then
        tt(1,1) = cmplx(t, alphay)
        tt(2,2) = cmplx(t, -alphay)
        tt(3,3) = cmplx(-t, alphay)
        tt(4,4) = cmplx(-t, -alphay)
        tt(1,2) = cmplx(alphaz, 0)
        tt(2,1) = cmplx(-alphaz, 0)
        tt(3,4) = cmplx(-alphaz, 0)
        tt(4,3) = cmplx(alphaz, 0)
!     endif
!     if (i==n/2) then
!        tt(1,1) = 0.01*cmplx(t, alphay)
!        tt(2,2) = 0.01*cmplx(t, -alphay)
!        tt(3,3) = 0.01*cmplx(-t, alphay)
!        tt(4,4) = 0.01*cmplx(-t, -alphay)
!        tt(1,2) = 0.01*cmplx(alphaz, 0)
!        tt(2,1) = 0.01*cmplx(-alphaz, 0)
!        tt(3,4) = 0.01*cmplx(-alphaz, 0)
!        tt(4,3) = 0.01*cmplx(alphaz, 0)           
!     endif  

     hh = cmplx(0._8,0._8)
     hh(1,1) = Bz - cp + 0.5d0*J_spin*c_z
     hh(2,2) = -Bz - cp - 0.5d0*J_spin*c_z
     hh(3,3) = Bz + cp + 0.5d0*J_spin*c_z
     hh(4,4) = -Bz + cp - 0.5d0*J_spin*c_z
     hh(1,3) = delta
     hh(3,1) = delta
     hh(2,4) = delta
     hh(4,2) = delta
     hh(1,2) = 0.5d0*J_spin*s_z*cmplx(c_xy,-s_xy)
     hh(2,1) = 0.5d0*J_spin*s_z*cmplx(c_xy,s_xy)
     hh(3,4) = hh(1,2)
     hh(4,3) = hh(2,1)

     mbar(1:4,1:4) = 0
     mbar(1:4,5:8) = transpose(conjg(tt))
     mbar(5:8,1:4) = -ainv(tt)
     mbar(5:8,5:8) = -matmul(ainv(tt),hh)

     mm = matmul(matmul(transpose(conjg(uu)),mbar),uu)
     
     ww1(1:4,1:4) = -matmul(ainv(mm(5:8,5:8)),mm(5:8,1:4))
     ww1(1:4,5:8) = ainv(mm(5:8,5:8))
     ww1(5:8,1:4) = mm(1:4,1:4) - matmul(matmul(mm(1:4,5:8),ainv(mm(5:8,5:8))),mm(5:8,1:4))
     ww1(5:8,5:8) = matmul(mm(1:4,5:8),ainv(mm(5:8,5:8)))

     if (i==1) then
        ww = ww1
     else
        ww = cdot(ww1,ww)
     endif
  enddo
  topo = real(det(ww(1:4,1:4)))
end function topo
!
!####################################################################################################
!
function cdot(w2,w1)
  implicit none
  interface
     function ainv(a)
       complex(8),dimension(4,4) :: ainv,a
     end function ainv
  end interface
  complex(8),dimension(8,8) :: cdot,w2,w1
  complex(8),dimension(4,4) :: u2,v2,r2,s2,u1,v1,r1,s1,ai,tmp1,tmp2
  integer :: i
  u2 = w2(1:4,1:4)
  v2 = w2(1:4,5:8)
  r2 = w2(5:8,1:4)
  s2 = w2(5:8,5:8)
  u1 = w1(1:4,1:4)
  v1 = w1(1:4,5:8)
  r1 = w1(5:8,1:4)
  s1 = w1(5:8,5:8)
  ai = 0
  do i=1,4; ai(i,i) = 1; enddo
     tmp1 = ainv(ai - matmul(u2,s1))
     cdot(1:4,1:4) = u1 + matmul(matmul(matmul(v1,tmp1),u2),r1)
     cdot(1:4,5:8) = matmul(matmul(v1,tmp1),v2)
     tmp2 = ainv(ai - matmul(s1,u2))
     cdot(5:8,1:4) = matmul(matmul(r2,tmp2),r1)
     cdot(5:8,5:8) = s2 + matmul(matmul(matmul(r2,tmp2),s1),v2)
end function cdot

function ainv(a)
  implicit none
  complex(8),dimension(4,4) :: ainv,a
  complex(8) :: work(4)
  integer(8),dimension(4) :: ipiv
  integer :: info
  ainv = a
  call zgetrf(4,4,ainv,4,ipiv,info)
  if (info /= 0) then
     print *,'GETRF',info
     stop 'GETRF'
  endif
  call zgetri(4,ainv,4,ipiv,work,4,info)
  if (info /= 0) then
     print *,'GETRI',info
     stop 'GETRI'
  endif
end function ainv

function det(mat)
  implicit none
  integer(kind=8),parameter :: N=4_8
  complex(kind=8) :: det
  complex(kind=8), intent(inout), dimension(N,N) :: mat
  integer(kind=8) :: i, info
  integer :: ipiv(N)
  real(kind=8) :: sgn 
  ipiv = 0
  call zgetrf(N, N, mat, N, ipiv, info)
  det = dcmplx(1,0)
  do i = 1, N
     det = det*mat(i, i)
  end do
  sgn = 1_8
  do i = 1, N
     if(ipiv(i) /= i) then
        sgn = -sgn
     end if
  end do
  det = sgn*det   
end function det
!
!#####################################################################################
!
