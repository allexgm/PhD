program poisson_parallel
    !
    !--> 2D Poisson equation solver (parallel code)
    !    Algorithms implemented:
    !    - Standard Gauss-Seidel (point by point)
    !    - Red-Black Gauss-Seidel (two groups of points)

    use mpi
    implicit none
    integer, parameter :: kind = 8
    real(kind), parameter :: pi = 3.141592653589793238462643383279502884197169399_kind
    integer, parameter :: N = 601 ! points in the grid per axis
    real(kind), parameter :: L_box = 20.0_kind ! size of the square box in Bohr
    integer, parameter :: iter_max = 400000 ! max number of iterations allowed
    real(kind) :: tol = 1.0E-7_kind ! convergence threshold
    real(kind), allocatable :: src(:,:) ! input source term (0 and N+1 are the boundaries)
    real(kind) :: sigma, x, y, src_1, src_2 ! used to build the source with gaussians
    real(kind), allocatable :: tmp(:,:)
    !real(kind) :: pot_gs(0:N+1,0:N+1), pot_rb(0:N+1,0:N+1) 
    real(kind), allocatable :: pot_gs(:,:), pot_rb(:,:) 
    
    integer :: i, j, iter, iter_gs, iter_rb, src_type
    real(kind) :: h, diff(1:N,1:N), res, t_start, t_end, t_gs, t_rb

    !new MPI variables:
    integer :: ierr, rank, numprocs
    integer :: dims(2), coords(2)
    logical :: periods(2)
    integer :: comm2d
    integer :: proc_left, proc_right, proc_down, proc_up
    
    integer :: nx, ny, i_start, i_end, j_start, j_end, p_nx, p_ny
    integer :: remainder_x, remainder_y
    
    integer :: global_i, global_j, global_i1, global_i2, global_jc
    real(kind) :: global_diff, global_res

    !0) MPI initialization and cartesian topology:
    !Basic variables:
    call MPI_Init(ierr)
    !MPI_Init: Initializes MPI enviroment. ierr = 0 if no error in all the MPI commands.
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    !MPI_Comm_rank: Returns the rank "rank" of the calling process in the communicator.
    call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
    !MPI_Comm_size: Returns the number of processes "numprocs" in the communicator. 
    !1 processor == serial code (but parallel code should also work)
    if (numprocs == 1) then
        print *, "With 1 processor the serial version of the code can be used instead of the parallel version!"
        !call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    !Create 2D grid of procesors:
    dims = [0,0]
    call MPI_Dims_create(numprocs, 2, dims, ierr)
    !MPI_Dims_create: Creates a balanced grid of processes. dims = [nproc_x, nproc_y] (number of processors in each axis)
    periods = [.false., .false.] !non periodic
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .false., comm2d, ierr)
    !MPI_Cart_create: Creates a new communicator to which topology information has been attached.
    call MPI_Cart_coords(comm2d, rank, 2, coords, ierr)
    !MPI_Cart_coords: Determines process coordinates in created topology given the rank of the process in the group of comm2d.
    !Get neigbourgs:
    call MPI_Cart_shift(comm2d, 0, 1, proc_left, proc_right, ierr)
    call MPI_Cart_shift(comm2d, 1, 1, proc_down, proc_up, ierr)
    !MPI_Cart_shift: Determines the source and destination ranks, given a shift direction (0: x-axis, 1: y-axis) and amount (1: neighbor). (proc_left, proc_right, proc_down, proc_up)

    !Calculate local variables in the subregion of the processor
    nx = N / dims(1)
    remainder_x = mod(N, dims(1))
    if (coords(1) < remainder_x) then
        nx = nx + 1
        i_start = coords(1) * nx + 1
    else
        i_start = coords(1) * nx + remainder_x + 1
    end if 
    i_end = i_start + nx - 1

    ny = N / dims(2)
    remainder_y = mod(N, dims(2))
    if (coords(2) < remainder_y) then
        ny = ny + 1
        j_start = coords(2) * ny + 1
    else
        j_start = coords(2) * ny + remainder_y + 1
    end if 
    j_end = j_start + ny - 1
    h = L_box / (N-1) ! h: stepsize in the grid

    !1) Build input source term:
    allocate(src(1:nx, 1:ny))
    src = 0.0_kind
    src_type = 1 
    !Choosing the source term (charge distribution) to solve Poisson's equation:
    if (src_type == 1) then
        !--> (A) puntual charge (electron) in the center:
        if ((N+1)/2 >= i_start .and. (N+1)/2 <= i_end .and. &
            (N+1)/2 >= j_start .and. (N+1)/2 <= j_end) then
            src((N+1)/2 - i_start + 1,(N+1)/2 - j_start + 1) = 4.0_kind * pi !rho = 1 au
            write(*,*) "Electron in the center of the box"
        end if
    else if (src_type == 2) then
        !--> (B) puntual charge representing I2- anion:
        !    (electron divided in both atom's centers)
        if ((N+1)/2-int(2.51/h) >= i_start .and. (N+1)/2-int(2.51/h) <= i_end .and. &
            (N+1)/2 >= j_start .and. (N+1)/2 <= j_end) then
            src((N+1)/2-int(2.51/h) - i_start + 1,(N+1)/2 - j_start + 1) = 2.0_kind * pi
        end if
        if ((N+1)/2+int(2.51/h) >= i_start .and. (N+1)/2+int(2.51/h) <= i_end .and. &
            (N+1)/2 >= j_start .and. (N+1)/2 <= j_end) then
            src((N+1)/2+int(2.51/h) - i_start + 1,(N+1)/2 - j_start + 1) = 2.0_kind * pi
        end if
        write(*,*) "I2- anion with 2 punctual charges"
    else if (src_type == 3) then
        !--> (C) more realistic I2- anion charge distribution:
        !    (gaussians instead of puntual charges)
        !    (where 2*sigma = van der waals radius)
        sigma = 3.74_kind/2.0_kind
        do j = j_start, j_end
            do i = i_start, i_end
                x = (i-1) * h
                y = (j-1) * h
                src_1 = 1.0_kind/sigma**2 * exp(- ((x - ((L_box/2.0_kind)-2.51_kind))**2 + (y - (L_box/2.0_kind))**2) &
                        / (2.0_kind * sigma**2))
                src_2 = 1.0_kind/sigma**2 * exp(- ((x - ((L_box/2.0_kind)+2.51_kind))**2 + (y - (L_box/2.0_kind))**2) &
                        / (2.0_kind * sigma**2))
                src(i - i_start + 1,j - j_start + 1) = src_1 + src_2 
            end do
        end do
        write(*,*) "I2- anion with 2 gaussian charges"
    end if

    call write_fun_parallel(src, nx, ny, i_start, i_end, j_start, j_end, N, h, "src.dat")

    !!!!TEST THAT UP TO HERE IT WORKS!!!!!!!!!!!!!!!!
    ! IT WORKS!!!!!!!!!!!!!!!!!!!
    ! :) :)

     !2) Standard GS algorithm:
    allocate(pot_gs(0:nx+1, 0:ny+1))
    pot_gs = 0.0_kind
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !MPI_Barrier: Blocks the calling process until all processes in the communicator have reached this routine.
    t_start = MPI_Wtime()
    iter_gs = 0
    do iter = 1, iter_max
        iter_gs = iter
        !tmp(1:nx,1:ny) = pot_gs(1:nx,1:ny)
        do j = 1, ny
            do i = 1, nx 
                pot_gs(i,j) = 0.25_kind * (pot_gs(i+1,j) + pot_gs(i-1,j) + &
                                           pot_gs(i,j+1) + pot_gs(i,j-1) + &
                                           h * h * src(i,j))     
            end do
        end do
        call exchange_boundaries(pot_gs, nx, ny, proc_left, proc_right, proc_down, proc_up, comm2d)
        call update_boundary_parallel(pot_gs, nx, ny, proc_left, proc_right, proc_down, proc_up)
        !Two ways to converge results:
        !!--> (A) converge iterations:
        !diff = abs(tmp(1:nx,1:ny) - pot_gs(1:nx,1:ny))
        !call MPI_Allreduce(diff, global_diff, nx*ny, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        !if (global_diff < tol) exit
        !--> (B) converge Poisson's equation:
        !check convergence only every 100 iterations to save time:
        !if (mod(iter, 100) == 0) then
            call compute_residual_parallel(pot_gs, src, nx, ny, h, res)
            call MPI_Allreduce(res, global_res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            !MPI_Allreduce: Combines values from all processes and distributes the result back to all processes. Here, it sums (MPI_SUM) the (MPI_DOUBLE_PRECISION) residuals across all processes.
            global_res = sqrt(global_res/real(N*N, kind))
            if (global_res < tol) exit
        !end if
    end do
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_end = MPI_Wtime()
    t_gs = t_end - t_start
    if (rank == 0) then
        write(*,*) 'iter',iter_gs 
        write(*,*) 't',t_gs 
    end if
    call write_fun_parallel(pot_gs(1:nx, 1:ny), nx, ny, i_start, i_end, j_start, j_end, N, h, "pot_gs.dat")
!     call write_fun(pot_gs(1:N,1:N), N, h, "pot_gs.dat")

    !3) Red-Black GS algorithm:
    allocate(pot_rb(0:nx+1, 0:ny+1))
    pot_rb = 0.0_kind
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_start = MPI_Wtime()
    iter_rb = 0
    do iter = 1, iter_max
        iter_rb = iter
!        tmp(1:N,1:N) = pot_rb(1:N,1:N)    
        !Red points
        do j = 1, ny
            global_j = j_start + j - 1
            do i = 1, nx
                global_i = i_start + i - 1
                if (mod(global_i+global_j, 2) == 0) then
                    pot_rb(i,j) = 0.25_kind * (pot_rb(i+1,j) + pot_rb(i-1,j) + &
                                               pot_rb(i,j+1) + pot_rb(i,j-1) + &
                                               h * h * src(i,j))     
                end if
            end do
        end do
        call exchange_boundaries(pot_rb, nx, ny, proc_left, proc_right, proc_down, proc_up, comm2d)
        call update_boundary_parallel(pot_rb, nx, ny, proc_left, proc_right, proc_down, proc_up)
        !Black points
        do j = 1, ny
            global_j = j_start + j - 1
            do i = 1, nx
                global_i = i_start + i - 1
                if (mod(global_i+global_j, 2) == 1) then
                    pot_rb(i,j) = 0.25_kind * (pot_rb(i+1,j) + pot_rb(i-1,j) + &
                                               pot_rb(i,j+1) + pot_rb(i,j-1) + &
                                               h * h * src(i,j))     
                end if
            end do
        end do
        call exchange_boundaries(pot_rb, nx, ny, proc_left, proc_right, proc_down, proc_up, comm2d)
        call update_boundary_parallel(pot_rb, nx, ny, proc_left, proc_right, proc_down, proc_up)
        !Two ways to converge results:
        !!--> (A) converge iterations:
        !diff = abs(tmp(1:N,1:N) - pot_rb(1:N,1:N))
        !call MPI_Allreduce(diff, global_diff, nx*ny, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        !if (global_diff < tol) exit
        !--> (B) converge Poisson's equation:
        !check convergence only every 100 iterations to save time:
        !if (mod(iter, 100) == 0) then
            call compute_residual_parallel(pot_rb, src, nx, ny, h, res)
            call MPI_Allreduce(res, global_res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            global_res = sqrt(global_res/real(N*N, kind))
            if (global_res < tol) exit
        !end if
    end do
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_end = MPI_Wtime()
    t_rb = t_end - t_start
    if (rank == 0) then
        write(*,*) 'iter',iter_rb 
        write(*,*) 't',t_rb 
    end if
    call write_fun_parallel(pot_rb(1:nx, 1:ny), nx, ny, i_start, i_end, j_start, j_end, N, h, "pot_rb.dat")
!     call write_fun(pot_rb(1:N,1:N), N, h, "pot_rb.dat")

    call MPI_Finalize(ierr)
    !MPI_Finalize: Terminates MPI environment. 
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

    ! Write solution as ASCII: x  y  phi   (one row per interior point)
    !
    subroutine write_fun_parallel(fun, nx, ny, is, ie, js, je, N, h, filename)
        real(kind), intent(in) :: fun(1:nx, 1:ny), h
        integer, intent(in) :: nx, ny, is, ie, js, je, N
        character(*), intent(in) :: filename
        
        real(kind), allocatable :: fun_global(:,:), buf(:,:)
        integer :: p, st(MPI_STATUS_SIZE), ierr, bounds(6)
        integer :: p_is, p_ie, p_js, p_je, p_nx, p_ny, unit_id
        !
        if (rank == 0) then
            allocate(fun_global(1:N, 1:N))
            fun_global = 0.0_kind
            fun_global(is:ie, js:je) = fun(1:nx, 1:ny)
            do p = 1, numprocs - 1
                !Recv proc limits and sizes:
                call MPI_Recv(bounds, 6, MPI_INTEGER, p, 0, MPI_COMM_WORLD, st, ierr)
                !MPI_Recv: Receives a message from a source process. 
                p_is = bounds(1)
                p_ie = bounds(2)
                p_js = bounds(3)
                p_je = bounds(4)  
                p_nx = bounds(5)
                p_ny = bounds(6)
                !Recv proc data:
                allocate(buf(1:p_nx, 1:p_ny))
                call MPI_Recv(buf, p_nx*p_ny, MPI_DOUBLE_PRECISION, p, 1, MPI_COMM_WORLD, st, ierr)
                 fun_global(p_is:p_ie, p_js:p_je) = buf
                 deallocate(buf)
            end do
            !And rank 0 will write it as the serial version:
            unit_id = 20
            open(unit=unit_id, file=filename, status="replace")
            write(unit_id,'(a)') "# x   y   fun"
            do j = 1, N
                do i = 1, N
                    write(unit_id,'(3es16.6e3)') i*h, j*h, fun_global(i,j)
                end do
                write(unit_id,*)   
            end do
            close(unit_id)    
            deallocate(fun_global)      
        else !for ranks =/ 0 (Senders)
            bounds = [is, ie, js, je, nx, ny]
            call MPI_Send(bounds, 6, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
            !MPI_Send: Sends a message to a destination process.
            allocate(buf(1:nx, 1:ny))
            buf = fun(1:nx,1:ny)
            call MPI_Send(buf, nx*ny, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, ierr)
            deallocate(buf)
        end if
    end subroutine write_fun_parallel

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
    
    !Update points outside the grid to satisfy boundary conditions (parallel version)
    !
    subroutine update_boundary_parallel(pot, nx, ny, proc_left, proc_right, proc_down, proc_up)
        integer, intent(in) :: nx, ny, proc_left, proc_right, proc_down, proc_up
        real(kind), intent(inout) :: pot(0:nx+1, 0:ny+1)
        !Update points outside the grid to satisfy boundary conditions
        !Two possible boundaries:
        !--> (A) Potential at the boundary = 0: (boundary is at 0.5 and N=0.5 in both axis)
        if(proc_left == MPI_PROC_NULL)  pot(0,1:ny)   = -pot(1,1:ny)
        if(proc_right == MPI_PROC_NULL) pot(nx+1,1:ny) = -pot(nx,1:ny)
        if(proc_down == MPI_PROC_NULL)  pot(1:nx,0)   = -pot(1:nx,1)
        if(proc_up == MPI_PROC_NULL)    pot(1:nx,ny+1) = -pot(1:nx,ny)
        !!--> (B) Derivative of the potential at the boundary = 0:
        !!if(proc_left == MPI_PROC_NULL)  pot(1:nx,0)   = pot(1:nx,1)
        !!if(proc_right == MPI_PROC_NULL) pot(1:nx,ny+1) = pot(1:nx,ny)
        !!if(proc_down == MPI_PROC_NULL)  pot(0,1:ny)   = pot(1,1:ny)
        !!if(proc_up == MPI_PROC_NULL)    pot(nx+1,1:ny) = pot(nx,1:ny)
    end subroutine update_boundary_parallel

    !Exchange boundaries between neighboring processors (parallel version)
    !
    subroutine exchange_boundaries(pot, nx, ny, proc_left, proc_right, proc_down, proc_up, comm2d)
        integer, intent(in) :: nx, ny, proc_left, proc_right, proc_down, proc_up
        integer, intent(in) :: comm2d
        real(kind), intent(inout) :: pot(0:nx+1, 0:ny+1)
        real(kind) :: sendbuf_left(ny), sendbuf_right(ny), sendbuf_down(nx), sendbuf_up(nx)
        real(kind) :: recvbuf_left(ny), recvbuf_right(ny), recvbuf_down(nx), recvbuf_up(nx)
        integer :: ierr

        sendbuf_left = pot(1, 1:ny)
        sendbuf_right = pot(nx,1:ny)
        sendbuf_down = pot(1:nx,1)  
        sendbuf_up = pot(1:nx,ny)
        !Send/Recv left and right boundaries:
        call MPI_Sendrecv(sendbuf_left, ny, MPI_DOUBLE_PRECISION, proc_left, 0, &
                          recvbuf_right, ny, MPI_DOUBLE_PRECISION, proc_right, 0, &
                          comm2d, MPI_STATUS_IGNORE, ierr)
        !MPI_Sendrecv: Performs a send and receive operation simultaneously. 
        !Here, it sends the left boundary (sendbuf_left) to the left neighbor (proc_left) 
        !and receives the right boundary (recvbuf_right) from the right neighbor (proc_right).
        call MPI_Sendrecv(sendbuf_right, ny, MPI_DOUBLE_PRECISION, proc_right, 1, &
                          recvbuf_left, ny, MPI_DOUBLE_PRECISION, proc_left, 1,&
                          comm2d,MPI_STATUS_IGNORE,ierr)
        !Send/Recv up and down boundaries:
        call MPI_Sendrecv(sendbuf_down, nx, MPI_DOUBLE_PRECISION, proc_down, 2,&
                          recvbuf_up, nx, MPI_DOUBLE_PRECISION, proc_up, 2,&
                          comm2d,MPI_STATUS_IGNORE,ierr)
        call MPI_Sendrecv(sendbuf_up, nx, MPI_DOUBLE_PRECISION, proc_up, 3,&
                          recvbuf_down, nx, MPI_DOUBLE_PRECISION, proc_down, 3,&
                          comm2d,MPI_STATUS_IGNORE,ierr)
        !Update the ghost cells with the received data:
        if (proc_left /= MPI_PROC_NULL) pot(0, 1:ny) = recvbuf_left
        if (proc_right /= MPI_PROC_NULL) pot(nx+1, 1:ny) = recvbuf_right
        if (proc_down /= MPI_PROC_NULL) pot(1:nx, 0) = recvbuf_down
        if (proc_up /= MPI_PROC_NULL) pot(1:nx, ny+1) = recvbuf_up
    end subroutine exchange_boundaries

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

    !L2 norm of the residual (parallel version):
    !
    subroutine compute_residual_parallel(pot, src, nx, ny, h, res)
        integer, intent(in) :: nx, ny
        real(kind), intent(in) :: pot(0:nx+1, 0:ny+1)
        real(kind), intent(in) :: src(1:nx, 1:ny)
        real(kind), intent(in) :: h
        real(kind), intent(out) :: res
        real(kind) :: r, lap
        integer :: i, j
        res = 0.0_kind
        do j = 1, ny
            do i = 1, nx
                lap = ( pot(i+1,j) + pot(i-1,j) &
                      + pot(i,j+1) + pot(i,j-1) &
                      - 4.0_kind*pot(i,j) ) / (h * h)
                r   = src(i,j) + lap          ! src - (-lap) = src + lap
                res = res + r*r
            end do
        end do
    end subroutine compute_residual_parallel
end program poisson_parallel
