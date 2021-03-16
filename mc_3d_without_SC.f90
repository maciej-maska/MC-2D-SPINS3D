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
  real(8) :: free_energy, r, free_en, old_free_en, en, old_en, c_z, ang_xy_old, ang_z_old, cosine, cosine_min = -1._8, an, an_up, an_dn
  complex(8), allocatable :: g(:,:), h(:,:), h1(:,:), u_up(:,:), u_dn(:,:), v_up(:,:), v_dn(:,:)
  complex(8) :: h_old(8)
  integer :: i, j, info, in, k, l, i_therm
  logical :: old_config
  namelist /parameters/ lx, ly, t, delta, alphay, alphaz, Bx, Bz, gamma_N, temperature, cp, disorder, magnetic_disorder, J_spin

  call system("echo $PPID >MC_PID")
  open(9,file='stripe.inp',status='OLD')
  read(9,nml=parameters)
  close(9)

  n = lx*ly

  allocate(e(2*n), g(2*n,2*n), h(2*n,2*n), h1(2*n,2*n), ang_xy(n), ang_z(n), rnd(n)) !, u_up(n,4*n), u_dn(n,4*n), v_up(n,4*n), v_dn(n,4*n))

  call random_seed()

  i_therm = 0
  inquire(FILE="spins_3d.dat", EXIST=old_config)
  if (old_config) then
     open(10,file="spins_3d.dat")
     do i=1,n
        read(10,*) ang_xy(i), ang_z(i)
     enddo
  else
     call random_number(rnd)
     ang_xy=2*pi*rnd                          ! phi
     call random_number(rnd)
     ang_z=acos(2*rnd - 1)                    ! theta
     !ang_z=0.5*pi
  endif
  call hamiltonian(h,ang_xy,ang_z)
  h1=h
  call heevd(h1, e, jobz='V', info=info)
!  u_up=h1(1:n,:)
!  u_dn=h1(n+1:2*n,:)
!  v_dn=h1(2*n+1:3*n,:)
!  v_up=h1(3*n+1:4*n,:)
!  an_up = sum(abs(u_up(:,:2*n))**2)/n  ! brakuje funkcji F-D (T=0)
!  an_dn = sum(abs(v_dn(:,2*n+1:))**2)/n  ! brakuje funkcji F-D (T=0)

  i = 0        
  old_free_en = free_energy(e)
  old_en = sum(e(1:n))
  an = 0.0
  do i = 1, 2*n
    if (e(i)<=0) an = an + 1.0  ! T=0, FD => theta Heviside
  enddo
  print *,old_free_en, an/n !, an_up+an_dn   !,topo(ang_xy,ang_z)
 
  main_loop: do !###############################################################
     call random_number(r)
     in = int(n*r) + 1
     ang_xy_old = ang_xy(in)
     ang_z_old = ang_z(in)
!     do
        call random_number(r)
        ang_xy(in) = 2*pi*r  
        call random_number(r)
        ang_z(in) = acos(2*r-1)
!        cosine = sin(ang_z(in))*sin(ang_z_old)*cos(ang_xy(in)-ang_xy_old) + cos(ang_z(in))*cos(ang_z_old) ! cosine between the new and old angles
!        if (cosine < cosine_min) exit
!     enddo
     !ang_z(in)=0.5*pi
     h_old(1) = h(in,in+n)
     h_old(2) = h(in+n,in)
     h_old(5) = h(in,in)
     h_old(6) = h(in+n,in+n)

     h(in,in+n) = 0.5d0*J_spin*sin(ang_z(in))*cmplx(cos(ang_xy(in)),-sin(ang_xy(in))) ! S.s
     h(in+n,in) = conjg(h(in,in+n))
     c_z=cos(ang_z(in))
     h(in,in)=Bz-cp+0.5d0*J_spin*c_z
     h(in+n,in+n)=-Bz-cp-0.5d0*J_spin*c_z

     h1=h
     call heevd(h1, e, info=info)
     free_en = free_energy(e)
     en = sum(e(1:n))
     call random_number(r)
     if ((free_en < old_free_en) .or. (old_free_en-free_en > temperature*log(r))) then
        old_free_en = free_en
        old_en = en
     else
        ang_xy(in) = ang_xy_old
        ang_z(in) = ang_z_old

        h(in,in+n) = h_old(1)
        h(in+n,in) = h_old(2)
        h(in,in) = h_old(5)
        h(in+n,in+n) = h_old(6)

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
!        h1=h
!        call heevd(h1, e, jobz='V', info=info)
!        u_up=h1(1:n,:)
!        u_dn=h1(n+1:2*n,:)
!        v_dn=h1(2*n+1:3*n,:)
!        v_up=h1(3*n+1:4*n,:)
!        an_up = sum(abs(u_up(:,:2*n))**2)/n  ! brakuje funkcji F-D (T=0)
!        an_dn = sum(abs(v_dn(:,2*n+1:))**2)/n  ! brakuje funkcji F-D (T=0)

        an = 0.0
        do i = 1, 2*n
          if (e(i)<=0) an = an + 1.0 ! T=0, FD => theta Heviside
        enddo
        an=an/n
        open(10,file='n.dat')
            write(10,*) an
        close(10)
        print *,free_en, an !, an_up+an_dn  !,topo(ang_xy,ang_z)
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
           !
           h(i,j+n) = alphaz                   !  C^+_up C_dn
           h(j,i+n) = -alphaz                  !  C^+_up C_dn
           !
        enddo
        c_z=cos(ang_z(i))
        h(i,i) = Bz-cp+0.5d0*J_spin*c_z
        h(i+n,i+n) = -Bz-cp-0.5d0*J_spin*c_z
        !
        h(i,i+n) = Bx
        !
        h(i,i+n) = 0.5d0*J_spin*sin(ang_z(i))*cmplx(cos(ang_xy(i)),-sin(ang_xy(i))) ! S.s
     enddo
  enddo
  
  do i=1,2*n-1  ! elementy pod przekatna hamiltonianu
     do j=i+1,2*n
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
  allocate(ipiv(2*n))
  g=-h
  do i=1,2*n
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
real(8),dimension(2*n) :: ens
real(8) :: free_energy,expon
integer :: i
free_energy=0.
do i=1,2*n
   ! expon=-(ens(i)-cp)/temperature
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
