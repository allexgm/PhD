program new_qmatrix
!
!----------------------------------------------------------------------------------+
! Code to compute the Q-matrix (primitive radial integral)                         |
! needed for the ECPs in a new way more efficient.                                 |
! See "Shaw, R. A., & Hill, J. G. (2017). The Journal of Chemical Physics, 147(7)" |        
!----------------------------------------------------------------------------------+
  
implicit none
real(kind=8) cms(3), alpha, beta, eta, p, position_A(3), position_B(3), A, B, ka, kb, nm
integer i, j, k, max_i, max_j, max_k, N, min_N, max_N, low_n, m, max_m
real(kind=8) comb, dawson_taylor, dawson_euler, dawson_cheb, x
real(kind=8), parameter :: pi = 3.14159265358979323846d0
real(kind=8), allocatable :: Xp(:), Xm(:), H(:), F(:), Ga(:), Gb(:), Q(:,:,:)

!Set the initial parameters:
!exponents of Gaussian in the basis set centered on atoms A and B respectively
alpha = 0.3d0
beta = 0.3d0
!exponent related to the ECP radial decay factor
eta = 0.5d0
!sum of exponents
p = eta + alpha + beta
!position of atoms A, B and the center of mass
position_A = (/0.d0, 0.d0, 0.d0/)
position_B = (/0.d0, 0.d0, 4.0d0/)
cms = (/0.d0, 0.d0, 2.0d0/)
A = norm2(position_A-cms) !distance of A to the center of mass
B = norm2(position_B-cms) !distnace of B to the center of mass
!ka and kb
ka = 2.d0 * alpha * A
kb = 2.d0 * beta * B
!dimensions of matrix Qijk
max_i = 5 !i: from 0 to 'max_l' for atom A
max_j = 5 !j: from 0 to 'max_l' for atom B
max_k = 6 !k: (in the paper)--> from 2 to 5/6 

!Calculate the lowest and the highest possible N fot this ranges of i, j and k
min_N = 2 - max_i - max_j !2 because is the min_k possible --> min_N possible
max_N = max_k + max_i
!If the first N is F type (= even number) then...
if (mod(min_N,2).ne.0) stop ": The lowest N (min_k-max_i-max_j) must be even!"
!...allocate the needed vectors (base integrals) and initialize them
!Ga has the same dimension as Gb, but we only need up to N=1 in Ga
allocate(Xp(max(min_N,0):max_N), Xm(max(min_N,0):max_N), H(min_N:max_N),F(min_N:max_N), Gb(min_N:max_N), Ga(min_N:1))
Xp = 0.d0
Xm = 0.d0     
H = 0.d0
F = 0.d0
Gb = 0.d0
Ga = 0.d0

!For each X(N) with N>=0 (from lowest N or 0 to highest N)...
do N=0, max_N
  !...apply eq.(36) 
  Xp(N) = ((beta*B+alpha*A)/p)**N * exp(p*((beta*B+alpha*A)/p)**2.d0) + &
          ((beta*B-alpha*A)/p)**N * exp(p*((beta*B-alpha*A)/p)**2.d0)
  Xm(N) = ((beta*B+alpha*A)/p)**N * exp(p*((beta*B+alpha*A)/p)**2.d0) - &
          ((beta*B-alpha*A)/p)**N * exp(p*((beta*B-alpha*A)/p)**2.d0)  
enddo
!For each possible N (from highest N to lowest N): [5 options]
do N=max_N, min_N, -1
  !1)-if N > or = than 2 and even...
  if (mod(N, 2) == 0.and.N.ge.2) then
    low_n = N/2
    max_m = low_n - 1
    do m=0, max_m 
      nm = real(low_n, kind=8) - real(m, kind=8) - 0.5d0
      !... apply eq.(34) 
      H(N) = H(N) + 0.25d0 * comb(2*low_n-2,2*m) * p**(-nm) * gamma(nm) * Xp(2*m)
      F(N) = F(N) + 0.25d0 * comb(2*low_n-2,2*m) * p**(-nm) * gamma(nm) * Xm(2*m)
    enddo
  !2)-if N > or = 2 and odd...
  else if (mod(N, 2) /= 0.and.N.ge.2) then
    low_n = (N-1)/2
    max_m = low_n - 1
    do m=0, max_m
      nm = (real(N-1, kind=8))/2 - real(m, kind=8) - 0.5d0
      !... apply eq.(35)
      Gb(2*low_n+1) = Gb(2*low_n+1) + 0.25d0 * comb(2*low_n-1,2*m+1) * p**(-nm) * gamma(nm) * Xm(2*m+1)
    enddo
  !3)-if N < 2 and odd...
  else if (mod(N, 2) /= 0.and.N.lt.2) then
    if (N==1) then
    !... if N=1 apply eq.(40)...
      Gb(1) = 0.5d0 * sqrt(pi) * (exp(p*((beta*B+alpha*A)/p)**2d0)*dawson_cheb(sqrt(p)*((beta*B+alpha*A)/p))- &
                                  exp(p*((beta*B-alpha*A)/p)**2d0)*dawson_cheb(sqrt(p)*((beta*B-alpha*A)/p)))               
      Ga(1) = 0.5d0 * sqrt(pi) * (exp(p*((alpha*A+beta*B)/p)**2d0)*dawson_cheb(sqrt(p)*((alpha*A+beta*B)/p))- &
                                  exp(p*((alpha*A-beta*B)/p)**2d0)*dawson_cheb(sqrt(p)*((alpha*A-beta*B)/p)))
    else
    !... otherwise apply eq.(38)
      Gb(N) = 2.d0*p/(N-1)*Gb(N+2)-ka/(N-1)*H(N+1)-kb/(N-1)*F(N+1)
      Ga(N) = 2.d0*p/(N-1)*Ga(N+2)-kb/(N-1)*H(N+1)-ka/(N-1)*F(N+1)
    endif
  endif     
  !4)-if N < 2 and even...
  if (mod(N, 2) == 0.and.N.lt.2) then 
  !... apply eq.(37) and eq.(39)
    F(N) = 2.d0*p/(N-1)*F(N+2)-ka/(N-1)*Ga(N+1)-kb/(N-1)*Gb(N+1)
    H(N) = 2.d0*p/(N-1)*H(N+2)-kb/(N-1)*Ga(N+1)-ka/(N-1)*Gb(N+1)
  endif
enddo

!Build the Q matrix
!allocate the Q matrix with the needed dimensions:
!i and j matrix indices can be "normal" (from 0 to max_i/j), 
!but for k index we need to go from lowest N to highest N 
!because we will build Q in a recursive way (from positions outside the "matrix limits")
allocate(Q(0:max_i,0:max_j,min_N:max_N))
Q(:,:,:) = 0.d0
!For Q00N
do N=min_N, max_N
  Q(0,0,N) = F(N)/(ka*kb)
enddo
!For Q01N apply eq.(33)
do N=min_N, max_N
  Q(0,1,N) = 1.d0/(ka*kb)*(1.d0*Gb(N)-1/kb*F(N-1))
enddo
!For Q0jk apply eq.(29), but j=0 and j=1 were already calculated (start in j=2)
do j=2, max_j
  do k=min_N, max_N
    if (mod(j+k,2) /= 0) cycle
    Q(0,j,k) = 1.d0*Q(0,j-2,k)-(2*j-1)/(2.d0*beta*B)*Q(0,j-1,k-1) 
  enddo
enddo
!For Qijk apply eq.(28), but i=0 were already calculated (start in i=1)
do i=1, max_i
  do j=i, max_j
    do k=min_N, max_N
      if (mod(i+j+k,2) /= 0) cycle
      Q(i,j,k) = real(2+j-i-k, kind=8)/(2.d0*alpha*A)*Q(i-1,j,k-1)-beta*B/(alpha*A)*Q(i-1,j-1,k)+p/(alpha*A)*Q(i-1,j,k+1)
    enddo
  enddo
enddo


!Plot the Q matrix
do k=2, max_k
  write(6,'(A,I2)') "Slice k =", k
  do j=max_j, 0, -1
    write(6,'(A,I2,A,100E20.10)') " j =",j,"  ",Q(:,j,k)
  enddo
    write(6,'(A)',advance='no') "                  "
  do i=0, max_i
    write(6,'(A,I2,15X)',advance='no') "i =",i
  enddo
  write(6,*) " "
  write(6,*) " "
enddo


deallocate(Xp, Xm, H, F, Gb, Ga, Q)
end program new_qmatrix

!Define the combinatorial function:
real(kind=8) function comb(a,b)
    implicit none
    integer, intent(in) :: a,b
    !factorial(x)=gamma(x+1)
    comb = gamma(real(a+1, kind=8)) / (gamma(real(b+1, kind=8))*gamma(real(a-b+1, kind=8)))
end function comb

!Define Dawson function with Chebyshev approximation:
!See https://www.netlib.org/specfun/daw
function dawson_cheb(x) result(dawson_result)
  double precision, intent(in) :: x
  double precision :: dawson_result
  double precision :: y, sum_p, sum_q, frac, w2
  integer :: i

  double precision, parameter :: zero = 0.0d0
  double precision, parameter :: half = 0.5d0
  double precision, parameter :: one = 1.0d0
  double precision, parameter :: six25 = 6.25d0
  double precision, parameter :: one225 = 12.25d0
  double precision, parameter :: two5 = 25.0d0

  double precision, parameter :: xsmall = 1.05d-08
  double precision, parameter :: xlarge = 9.49d+07
  double precision, parameter :: xmax = 2.24d+307

  double precision, parameter :: p1(10) = (/ &
    -2.69020398788704782410d-12, 4.18572065374337710778d-10, &
    -1.34848304455939419963d-08, 9.28264872583444852976d-07, &
    -1.23877783329049120592d-05, 4.07205792429155826266d-04, &
    -2.84388121441008500446d-03, 4.70139022887204722217d-02, &
    -1.38868086253931995101d-01, 1.00000000000000000004d+00 /)
  double precision, parameter :: q1(10) = (/ &
    1.71257170854690554214d-10, 1.19266846372297253797d-08, &
    4.32287827678631772231d-07, 1.03867633767414421898d-05, &
    1.78910965284246249340d-04, 2.26061077235076703171d-03, &
    2.07422774641447644725d-02, 1.32212955897210128811d-01, &
    5.27798580412734677256d-01, 1.00000000000000000000d+00 /)
  double precision, parameter :: p2(10) = (/ &
    -1.70953804700855494930d+00, -3.79258977271042880786d+01, &
    2.61935631268825992835d+01, 1.25808703738951251885d+01, &
    -2.27571829525075891337d+01, 4.56604250725163310122d+00, &
    -7.33080089896402870750d+00, 4.65842087940015295573d+01, &
    -1.73717177843672791149d+01, 5.00260183622027967838d-01 /)
  double precision, parameter :: q2(9) = (/ &
    1.82180093313514478378d+00, 1.10067081034515532891d+03, &
    -7.08465686676573000364d+00, 4.53642111102577727153d+02, &
    4.06209742218935689922d+01, 3.02890110610122663923d+02, &
    1.70641269745236227356d+02, 9.51190923960381458747d+02, &
    2.06522691539642105009d-01 /)
  double precision, parameter :: p3(10) = (/ &
    -4.55169503255094815112d+00, -1.86647123338493852582d+01, &
    -7.36315669126830526754d+00, -6.68407240337696756838d+01, &
    4.84507265081491452130d+01, 2.69790586735467649969d+01, &
    -3.35044149820592449072d+01, 7.50964459838919612289d+00, &
    -1.48432341823343965307d+00, 4.99999810924858824981d-01 /)
  double precision, parameter :: q3(9) = (/ &
    4.47820908025971749852d+01, 9.98607198039452081913d+01, &
    1.40238373126149385228d+01, 3.48817758822286353588d+03, &
    -9.18871385293215873406d+00, 1.24018500009917163023d+03, &
    -6.88024952504512254535d+01, -2.31251575385145143070d+00, &
    2.50041492369922381761d-01 /)
  double precision, parameter :: p4(10) = (/ &
    -8.11753647558432685797d+00, -3.84043882477454453430d+01, &
    -2.23787669028751886675d+01, -2.88301992467056105854d+01, &
    -5.99085540418222002197d+00, -1.13867365736066102577d+01, &
    -6.52828727526980741590d+00, -4.50002293000355585708d+00, &
    -2.50000000088955834952d+00, 5.00000000000000488400d-01 /)
  double precision, parameter :: q4(9) = (/ &
    2.69382300417238816428d+02, 5.04198958742465752861d+01, &
    6.11539671480115846173d+01, 2.08210246935564547889d+02, &
    1.97325365692316183531d+01, -1.22097010558934838708d+01, &
    -6.99732735041547247161d+00, -2.49999970104184464568d+00, &
    7.49999999999027092188d-01 /)

  y = abs(x)

  if (y > xlarge) then
    if (y <= xmax) then
      dawson_result = half / x
    else
      dawson_result = zero
    end if
  else if (y < xsmall) then
    dawson_result = x
  else
    y = x * x
    if (y < six25) then
      ! abs(x) < 2.5
      sum_p = p1(1)
      sum_q = q1(1)
      do i = 2, 10
        sum_p = sum_p * y + p1(i)
        sum_q = sum_q * y + q1(i)
      end do
      dawson_result = x * sum_p / sum_q
    else if (y < one225) then
      ! 2.5 <= abs(x) < 3.5
      frac = zero
      do i = 1, 9
        frac = q2(i) / (p2(i) + y + frac)
      end do
      dawson_result = (p2(10) + frac) / x
    else if (y < two5) then
      ! 3.5 <= abs(x) < 5.0
      frac = zero
      do i = 1, 9
        frac = q3(i) / (p3(i) + y + frac)
      end do
      dawson_result = (p3(10) + frac) / x
    else
      ! 5.0 <= abs(x) <= xlarge
      w2 = one / (x * x)
      frac = zero
      do i = 1, 9
        frac = q4(i) / (p4(i) + y + frac)
      end do
      frac = p4(10) + frac
      dawson_result = (half + half * w2 * frac) / x
    end if
  end if
end function dawson_cheb
