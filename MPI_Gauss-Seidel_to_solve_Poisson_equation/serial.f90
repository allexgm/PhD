program poisson_serial
    !
    !--> 2D Poisson equation solver (serial code)
    !    Algorithms implemented:
    !    - Standard Gauss-Seidel (point by point)
    !    - Red-Black Gauss-Seidel (two groups of points)

    implicit none
    integer, parameter :: kind = 8
    real(kind), parameter :: pi = 3.141592653589793238462643383279502884197169399_kind
    integer, parameter :: N = 601 ! points in the grid per axis
    real(kind), parameter :: L_box = 20.0_kind ! size of the square box in Bohr
    integer, parameter :: iter_max = 400000 ! max number of iterations allowed
    real(kind) :: tol = 1.0E-7_kind ! convergence threshold
    real(kind) :: src(1:N,1:N) ! input source term (0 and N+1 are the boundaries)
    real(kind) :: sigma, x, y, src_1, src_2 ! used to build the source with gaussians
    real(kind) :: tmp(1:N,1:N)
    real(kind) :: pot_gs(0:N+1,0:N+1) ! potential obtained with standarrd GS
    real(kind) :: pot_rb(0:N+1,0:N+1) ! potential obtained with Red-Black GS
    
    integer :: i, j, iter, iter_gs, iter_rb, src_type
    real(kind) :: h, diff(1:N,1:N), res, t_start, t_end, t_gs, t_rb

    h = L_box / (N-1) ! h: stepsize in the grid

    !1) Build input source term:
    src = 0.0_kind
    src_type = 1 ! 1: electron in the center, 2: I2- anion, 3: I2- anion with gaussians
    if (src_type == 1) then
        !--> (A) puntual charge (electron) in the center:
        src((N+1)/2,(N+1)/2) = 4.0_kind * pi !rho = 1 au
    else if (src_type == 2) then
        !--> (B) puntual charge representing I2- anion:
        !    (electron divided in both atom's centers)
        src((N+1)/2-int(2.51/h),(N+1)/2) = 2.0_kind * pi
        src((N+1)/2+int(2.51/h),(N+1)/2) = 2.0_kind * pi
    else if (src_type == 3) then
        !--> (C) more realistic I2- anion charge distribution:
        !    (gaussians instead of puntual charges)
        !    (where 2*sigma = van der waals radius)
        sigma = 3.74_kind/2.0_kind
        do j = 1, N
            do i = 1, N
                x = (i-1) * h
                y = (j-1) * h
                src_1 = 1.0_kind/sigma**2 * exp(- ((x - ((L_box/2.0_kind)-2.51_kind))**2 + (y - (L_box/2.0_kind))**2) &
                        / (2.0_kind * sigma**2))
                src_2 = 1.0_kind/sigma**2 * exp(- ((x - ((L_box/2.0_kind)+2.51_kind))**2 + (y - (L_box/2.0_kind))**2) &
                        / (2.0_kind * sigma**2))
                src(i,j) = src_1 + src_2 
            end do
        end do
    end if
    !
    call write_fun(src, N, h, "src.dat")
    
    !2) Standard GS algorithm:
    pot_gs = 0.0_kind
    call cpu_time(t_start)
    iter_gs = 0
    do iter = 1, iter_max
        iter_gs = iter
        tmp(1:N,1:N) = pot_gs(1:N,1:N)
        do j = 1, N
            do i = 1, N 
                pot_gs(i,j) = 0.25_kind * (pot_gs(i+1,j) + pot_gs(i-1,j) + &
                                           pot_gs(i,j+1) + pot_gs(i,j-1) + &
                                           h * h * src(i,j))     
            end do
        end do
        call update_boundary(pot_gs,N)
        !Two ways to converge results:
        !!--> (A) converge iterations:
        !diff = abs(tmp(1:N,1:N) - pot_gs(1:N,1:N))
        !if (maxval(diff) < tol) exit
        !--> (B) converge Poisson's equation:
        call compute_residual(pot_gs, src, N, h, res)
        if (res < tol) exit
    end do
    call cpu_time(t_end)
    t_gs = t_end - t_start
    write(*,*) 'iter',iter_gs 
    write(*,*) 't',t_gs 
    call write_fun(pot_gs(1:N,1:N), N, h, "pot_gs.dat")

    !3) Red-Black GS algorithm:
    pot_rb = 0.0_kind
    call cpu_time(t_start)
    iter_rb = 0
    do iter = 1, iter_max
        iter_rb = iter
        tmp(1:N,1:N) = pot_rb(1:N,1:N)    
        !Red points
        do j = 1, N
            do i = 1, N
                if (mod(i+j, 2) == 0) then
                    pot_rb(i,j) = 0.25_kind * (pot_rb(i+1,j) + pot_rb(i-1,j) + &
                                               pot_rb(i,j+1) + pot_rb(i,j-1) + &
                                               h * h * src(i,j))     
                end if
            end do
        end do
        !Black points
        do j = 1, N
            do i = 1, N
                if (mod(i+j, 2) == 1) then
                    pot_rb(i,j) = 0.25_kind * (pot_rb(i+1,j) + pot_rb(i-1,j) + &
                                               pot_rb(i,j+1) + pot_rb(i,j-1) + &
                                               h * h * src(i,j))     
                end if
            end do
        end do
        call update_boundary(pot_rb,N)
        !Two ways to converge results:
        !!--> (A) converge iterations:
        !diff = abs(tmp(1:N,1:N) - pot_rb(1:N,1:N))
        !if (maxval(diff) < tol) exit
        !--> (B) converge Poisson's equation:
        call compute_residual(pot_rb, src, N, h, res)
        if (res < tol) exit
    end do
    call cpu_time(t_end)
    t_rb = t_end - t_start
    write(*,*) 'iter',iter_rb
    write(*,*) 't',t_rb
    call write_fun(pot_rb(1:N,1:N), N, h, "pot_rb.dat")


    contains
    
    ! Write solution as ASCII: x  y  phi   (one row per interior point)
    !
    subroutine write_fun(fun, N, h, filename)
        integer,        intent(in) :: N
        real(kind),     intent(in) :: fun(1:N,1:N), h
        character(*),   intent(in) :: filename
        integer :: i, j, unit_id
        unit_id = 20
        open(unit=unit_id, file=filename, status="replace")
        write(unit_id,'(a)') "# x   y   fun"
        do j = 1, N
            do i = 1, N
                write(unit_id,'(3es16.6e3)') i*h, j*h, fun(i,j)
            end do
            write(unit_id,*)   
        end do
        close(unit_id)
    end subroutine write_fun

    !Update points outside the grid to satisfy boundary conditions
    !
    subroutine update_boundary(pot, N)
        integer,        intent(in) :: N
        real(kind),     intent(inout) :: pot(0:N+1,0:N+1)
        !Two possible boundaries:
        !--> (A) Potential at the boundary = 0: (boundary is at 0.5 and N=0.5 in both axis)
        pot(1:N,0)   = -pot(1:N,1)
        pot(1:N,N+1) = -pot(1:N,N)
        pot(0,1:N)   = -pot(1,1:N)
        pot(N+1,1:N) = -pot(N,1:N)
        !!--> (B) Derivative of the potential at the boundary = 0:
        !!pot(1:N,0)   = pot(1:N,1)
        !!pot(1:N,N+1) = pot(1:N,N)
        !!pot(0,1:N)   = pot(1,1:N)
        !!pot(N+1,1:N) = pot(N,1:N)
    end subroutine update_boundary 
    
    !L2 norm of the residual:
    !   res(i,j) = src(i,j) - [ -nabla^2 pot ]_{i,j}
    !            = src(i,j) - (4*pot(i,j) - pot(i+h,j) - pot(i-h,j) - pot(i,j+h) - pot(i,j-h))/h^2
    ! *I return  sqrt( sum(res^2) / N^2 ) 
    !
    subroutine compute_residual(pot, src, N, h, res)
        integer, intent(in) :: N
        real(kind), intent(in) :: pot(0:N+1, 0:N+1)
        real(kind), intent(in) :: src(1:N, 1:N)
        real(kind), intent(in) :: h
        real(kind), intent(out) :: res
        real(kind) :: r, lap
        integer :: i, j
        res = 0.0_kind
        do j = 1, N
            do i = 1, N
                lap = ( pot(i+1,j) + pot(i-1,j) &
                      + pot(i,j+1) + pot(i,j-1) &
                      - 4.0_kind*pot(i,j) ) / (h * h)
                r   = src(i,j) + lap          ! src - (-lap) = src + lap
                res = res + r*r
            end do
        end do
        res = sqrt(res / real(N*N, kind))
    end subroutine compute_residual

end program poisson_serial
