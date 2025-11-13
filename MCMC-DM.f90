!===============================================================================
! PROGRAM: MCMC-DM.f90
!
! DESCRIPTION:
!   This program performs parallel Markov Chain Monte Carlo (MCMC) sampling
!   using the Metropolis–Hastings algorithm to estimate the parameters of
!   the Relativistic MOND (Modified Newtonian Dynamics) model.
!
!   The goal is to fit galaxy rotation curve data and to compare the results
!   of the Relativistic MOND model with those obtained from Newtonian gravity
!   and Dark Matter halo models. The analysis focuses on the outer stellar
!   regions of galaxies where deviations from Newtonian predictions appear.
!
!   The program uses MPI (Message Passing Interface) to run multiple MCMC
!   chains concurrently, each on a separate processor. Chains are synchronized
!   periodically to monitor convergence using the Gelman–Rubin (R̂) diagnostic.
!
! FEATURES:
!   - Parallel MCMC sampling with MPI
!   - Parameter estimation for the Relativistic MOND model
!   - Comparison with Newtonian and Dark Matter rotation curve fits
!   - High-precision modeling of outer stellar rotation velocities
!   - Dynamic convergence checking via Gelman–Rubin statistic
!   - Modular and documented Fortran 90 structure
!
! INPUT:
!   - Initial parameter values (initpoint)
!   - Prior bounds (priormin, priormax)
!   - Jump sizes (jumpsize)
!   - Observational rotation curve data files
!
! OUTPUT:
!   - Accepted parameter samples (MCMC chains) compatible with GetDist
!   - Convergence diagnostics (R̂ statistics)
!   - Best-fit rotation curve predictions for MOND, Newtonian, and Dark Matter
!
! USAGE:
!   - Compile with MPI Fortran compiler, e.g.:
!         mpiifort -O2 MCMC-DM.f90 -o MCMC-DM
!   - Run with desired number of parallel chains:
!         mpirun -np <nchains> ./MCMC-DM
!
! NOTES:
!   - Each MPI process runs an independent chain; convergence information is
!     shared periodically through MPI communications.
!   - The resulting chain files can be analyzed using GetDist or similar tools.
!
! REPOSITORY:
!   https://github.com/krezazadeh/MCMC-RMOND
!
! REFERENCE:
!   https://www.arxiv.org/abs/2511.05632
!
! AUTHOR:
!   Kazem Rezazadeh
!   School of Astronomy,
!   Institute for Research in Fundamental Sciences (IPM)
!   Email: kazem.rezazadeh@ipm.ir
!   GitHub: https://github.com/krezazadeh
!
! DATE:
!   5 October 2025
!===============================================================================

program MCMC

    use mpi

    implicit none
    
    ! Constants
    real(8), parameter :: pi = acos(-1.0d0)
    real(8), parameter :: G = 6.6743015d-11 ! m^3/(kg*s)
    real(8), parameter :: c = 2.99792458d8 ! m/s
    real(8), parameter :: hbar = 1.054571817d-34 ! J*s
    real(8), parameter :: kB = 1.380649d-23 ! (m^2*kg)/(s^2*K)
    real(8), parameter :: Gt = 4.30091d-6 ! kpc (km/s)^2 Msun^-1
    real(8), parameter :: kpc = 3.0856775814913673d19 ! m
    real(8), parameter :: MSun = 1.989d30 ! kg

    ! Model parameters

    ! Constant parameters of the model

    ! Unused parameters
    real(8), parameter :: epsilon = 0.0d0
    real(8), parameter :: r0 = 0.0d0
    
    ! Number of varying parameters of the model
    integer, parameter :: nparams = 4
    
    ! Varying parameters of the model
    real(8) :: sigma0, Rd, rhoDM0, rc
    
    ! NGC3198
    integer, parameter :: ndataNGC3198 = 43
    real(8) :: dataNGC3198(ndataNGC3198, 3)

    ! MCMC
    integer, parameter :: max_try = 500
    integer, parameter :: max_iter = 1000
    real(8), parameter :: Rhat_minus_1_tol = 0.01d0
    integer :: ierr, rank, nprocess
    integer :: try
    real(8), allocatable :: points(:, :, :) ! points(i_process, i_point, i_param)
    integer :: npoints
    real(8), allocatable :: points_local(:, :, :) ! points_local(rank + 1, i_iter, i_param) 
    real(8), allocatable :: points_temp(:, :, :) ! points_local(rank + 1, i_iter, i_param) 
    integer :: i, j, k
    character(len=30) :: filename
    integer :: unit
    real(8) :: params(nparams), params_new(nparams)
    real(8), dimension(nparams) :: initpoint
    real(8), dimension(nparams) :: jumpsize
    real(8), dimension(nparams) :: priormin, priormax
    real(8) :: chi2params, chi2params_new
    real(8) :: alpha
    real(8) :: rand
    integer :: NN
    integer :: nchains
    real(8), allocatable :: mean_chain(:, :) ! (i_process, i_param)
    real(8), allocatable :: mean_param(:), B(:), W(:), Rhat(:) ! (i_param)
    real(8) :: Rhat_minus_1
    logical :: converged
    integer, parameter :: sendcount = max_iter*nparams
    integer, allocatable :: recvcounts(:), displs(:)
    real(8), allocatable :: sendbuf(:), recvbuf(:)

    real(8) :: chi2min

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocess, ierr)

    call NGC3198()

    allocate(points_local(nprocess, max_iter, nparams))
    allocate(sendbuf(sendcount))

    if (rank == 0) then
        allocate(points(nprocess, max_iter*max_try, nparams))
        allocate(points_temp(nprocess, max_iter, nparams))
        allocate(mean_chain(nprocess, nparams))
        allocate(mean_param(nparams))
        allocate(B(nparams))
        allocate(W(nparams))
        allocate(Rhat(nparams))
        allocate(recvbuf(nprocess*sendcount))
        allocate(recvcounts(nprocess))
        allocate(displs(nprocess))
    end if

    unit = 10 + rank
    write(filename, '(A,I0,A)') 'chains_', rank, '.txt'
    open(unit=unit, file=filename, status='replace')
    
    call random_seed()

    initpoint = (/ 2.82646d0, 5.11496d0, 7.568d0, 3.06157d0 /)
    priormin = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    priormax = (/ 10.0d0, 20.0d0, 10.0d0, 100.0d0 /)
    jumpsize = (/ 0.01d0, 0.01d0, 0.01d0, 0.01d0 /)

    ! start from the given initial point.
    do i = 1, nparams
        params(i) = initpoint(i)
    end do

!     ! start from an arbitrary point in the prior intervals.
!     do i = 1, nparams
!         call random_number(rand)
!         params(i) = priormin(i) + (priormax(i) - priormin(i))*rand
!     end do

    sigma0 = params(1)*1.0d8
    Rd = params(2)
    rhoDM0 = 10.0d0**params(3)
    params(4) = 3.06157d0
    rc = params(4)
    chi2params = chi2total()

!     print *, chi2params
!     stop

    chi2min = chi2params

    converged = .false.

    if (rank == 0) then
        try = 0
    end if

    do while (.not. converged)
    
    if (rank == 0) then
        try = try + 1
        print *, "Attempt: ", try
        if (try > max_try) then
            print *, "WARNING: MCMC completed all iterations but failed to converge."
            print *, "Consider increasing the number of trys or adjusting initial chain parameters."
            stop
        end if
    end if

    do i = 1, max_iter

        call propose_jump(params, jumpsize, params_new)

        do j = 1, nparams
            if (params_new(j) < priormin(j)) params_new(j) = priormin(j)
            if (params_new(j) > priormax(j)) params_new(j) = priormax(j)
        end do

        sigma0 = params_new(1)*1.0d8
        Rd = params_new(2)
        rhoDM0 = 10.0d0**params_new(3)
        params_new(4) = 3.06157d0
        rc = params_new(4)
        chi2params_new = chi2total()

        points_local(rank + 1, i, :) = params_new

        write(10 + rank, "(9e25.16)") 1.0d0, chi2params_new/2.0d0, params_new &
        , Md(dataNGC3198(ndataNGC3198, 1)), MDM(dataNGC3198(ndataNGC3198, 1)) &
        , Mtotal(dataNGC3198(ndataNGC3198, 1))

        alpha = min(1.0d0, exp(-0.5d0 * (chi2params_new - chi2params)))
        call random_number(rand)
        if (rand <= alpha) then
            params = params_new
            chi2params = chi2params_new
        end if

    end do

    ! Computation of Gelman-Rubin parameter (Rhat)
    ! The convergence step (Gelman-Rubin) is performed at rank zero.

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call pack_2D(points_local(rank+1, :, :), sendbuf)

    if (rank == 0) then
        do i = 0, nprocess - 1
            recvcounts(i+1) = sendcount
            displs(i+1) = i*sendcount
        end do
    end if
    
    call MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE_PRECISION, &
    recvbuf, recvcounts, displs, MPI_DOUBLE_PRECISION, &
    0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then

        call unpack_3D(recvbuf, points_temp)

        npoints = nprocess*max_iter*try
        nchains = nprocess
        NN = max_iter*try
        
        do i = 1, nprocess
            do j = 1, max_iter
                do k = 1, nparams
                    points(i, (try - 1)*max_iter + j, k) = points_temp(i, j, k)
                end do
            end do
        end do

        mean_chain = 0.0
        do k = 1, nparams
            do i = 1, nprocess
                do j = 1, NN
                    mean_chain(i, k) = mean_chain(i, k) + points(i, j, k)
                end do
                mean_chain(i, k) = mean_chain(i, k) / real(NN, 8)
            end do
        end do

        mean_param = 0.0
        do k = 1, nparams
            do i = 1, nprocess
                do j = 1, NN
                    mean_param(k) = mean_param(k) + points(i, j, k)
                end do
            end do
            mean_param(k) = mean_param(k) / real(nprocess*NN, 8)
        end do

        B = 0.0
        W = 0.0
        Rhat = 0.0d0
        do k = 1, nparams
            do i = 1, nprocess
                B(k) = B(k) + (mean_chain(i, k) - mean_param(k))**2
                do j = 1, NN
                    W(k) = W(k) + (points(i, j, k) - mean_chain(i, k))**2
                end do
            end do
            B(k) = B(k)/real(nprocess - 1, 8)
            W(k) = W(k) / real(nprocess*(NN - 1), 8)
            Rhat(k) = (((real(NN - 1, 8))/real(NN, 8))*W(k) + &
            (1.0d0 + 1.0d0/real(nprocess, 8))*B(k))/W(k)
        end do

        Rhat_minus_1 = maxval(abs(Rhat - 1.0d0))
        converged = (Rhat_minus_1 < Rhat_minus_1_tol)

        if (converged) then
            print *, "Converged: Rhat - 1 = ", Rhat_minus_1
        else
            print *, "Not converged: Rhat - 1 = ", Rhat_minus_1
        end if
    
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call MPI_BCAST(converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    end do ! while

    deallocate(points_local)
    deallocate(sendbuf)

    if (rank == 0) then
        deallocate(points)
        deallocate(points_temp)
        deallocate(mean_chain)
        deallocate(mean_param)
        deallocate(B)
        deallocate(W)
        deallocate(Rhat)
        deallocate(recvbuf)
        deallocate(recvcounts)
        deallocate(displs)
    end if

    close(unit)

    if (rank == 0) print *, "Done!"

    call MPI_Finalize(ierr)

contains

subroutine propose_jump(center, jumpsize, proposal)
    real(8), intent(in) :: center(nparams)
    real(8), intent(in) :: jumpsize(nparams)
    real(8), intent(out) :: proposal(nparams)
    real(8) :: randn(nparams), z
    integer :: i
    do i = 1, nparams
        call random_number(z)
        randn(i) = sqrt(-2.d0 * log(z)) * cos(2.d0 * 3.141592653d0 * z)
        proposal(i) = center(i) + randn(i) * jumpsize(i)
    end do
end subroutine propose_jump

subroutine pack_2D(array2D, array1D)
    real(8), intent(in) :: array2D(:, :)
    real(8), intent(out) :: array1D(:)
    integer :: i, j, k
    k = 1
    do i = 1, size(array2D, 1)
        do j = 1, size(array2D, 2)
            array1D(k) = array2D(i, j)
            k = k + 1
        end do
    end do
end subroutine

subroutine unpack_3D(array1D, array3D)
    real(8), intent(in) :: array1D(:)
    real(8), intent(out) :: array3D(:, :, :)
    integer :: i, j, k, idx
    idx = 1
    do i = 1, size(array3D, 1)
        do j = 1, size(array3D, 2)
            do k = 1, size(array3D, 3)
                array3D(i, j, k) = array1D(idx)
                idx = idx + 1
            end do
        end do
    end do
end subroutine

function chi2total()

    implicit none

    real(8) :: chi2total

    chi2total = chi2NGC3198()

end function chi2total

subroutine NGC3198()

    implicit none
    
    integer :: i, j
    integer :: ios
    character(len=200) :: line

    real(8) :: Rad(ndataNGC3198), Vobs(ndataNGC3198), errV(ndataNGC3198), Vgas(ndataNGC3198), Vdisk(ndataNGC3198), Vbul(ndataNGC3198), SBdisk(ndataNGC3198), SBbul(ndataNGC3198)

    open (unit = 11, file = './NGC3198_rotmod.dat', status = 'old')

    read(11,'(A)', iostat=ios) line
    read(11,'(A)', iostat=ios) line
    read(11,'(A)', iostat=ios) line
    
    do i = 1, ndataNGC3198
        read(11, *) Rad(i), Vobs(i), errV(i), Vgas(i), Vdisk(i), Vbul(i), SBdisk(i), SBbul(i)
    end do

    close(11)

    ! We use this to prevent compiler error.    
    do i = 1, ndataNGC3198
        dataNGC3198(i, 1) = Rad(i)
        dataNGC3198(i, 2) = Vobs(i)
        dataNGC3198(i, 3) = errV(i)
    end do

end subroutine NGC3198

function chi2NGC3198()

    implicit none

    real(8) :: chi2NGC3198

    integer :: i
    real(8) :: sum

    sum = 0.0d0

    do i = 1, ndataNGC3198
        sum = sum + (vModel(dataNGC3198(i, 1)) - dataNGC3198(i, 2))**2/ &
        dataNGC3198(i, 3)**2
    end do

!      do i = 1, ndataNGC3198
!          sum = sum + (vSchwartzchild(dataNGC3198(i, 1)) - dataNGC3198(i, 2))**2/ &
!          dataNGC3198(i, 3)**2
!      end do

    chi2NGC3198 = sum

end function chi2NGC3198

function Md(r)

    implicit none

    real(8) :: Md
    real(8) :: r

    Md = 2.0d0*pi*Rd*(Rd - (r + Rd)/exp(r/Rd))*sigma0

end function Md

function MDM(r)

    implicit none

    real(8) :: MDM
    real(8) :: r

    MDM = 4.0d0*pi*rc**2*rhoDM0*(r - rc*atan(r/rc))

end function MDM

function Mtotal(r)

    implicit none

    real(8) :: Mtotal
    real(8) :: r

    Mtotal = Md(r) + MDM(r)

end function Mtotal

function vND(r)

    implicit none

    real(8) :: vND
    real(8) :: r

    vND = sqrt((Gt*Md(r))/r)

end function vND

function vDM(r)

    implicit none

    real(8) :: vDM
    real(8) :: r

    vDM = sqrt((Gt*(MDM(r) + Md(r)))/r)

end function vDM

function vRMOND(r)

    implicit none

    real(8) :: vRMOND
    real(8) :: r

    if (r <= r0) then
        vRMOND = sqrt((Gt*Md(r))/r)
    else
        vRMOND = (2.0d0**(-3.0d0 - epsilon/(1.0d0 + epsilon))* &
        c**(1.0d0 - (3.0d0*epsilon)/(2.0d0*(1.0d0 + epsilon)))* &
        sqrt(C1())*sqrt(epsilon)* &
        G**(epsilon/(2.0d0 + 2.0d0*epsilon))* &
        hbar**(epsilon/(2.0d0 + 2.0d0*epsilon)))/ &
        (125.0d0*(1.0d0 + epsilon)**(epsilon/(1.0d0 + epsilon))* &
        (kpc*r)**(epsilon/(1.0d0 + epsilon)))
    end if

end function vRMOND

function C1()

    implicit none

    real(8) :: C1

    C1 = (15625.0d0*2.0d0**(6.0d0 + (2.0d0*epsilon)/(1.0d0 + epsilon))* &
    c**(-2.0d0 + (3.0d0*epsilon)/(1.0d0 + epsilon))* &
    (1.0d0 + epsilon)**((2.0d0*epsilon)/(1.0d0 + epsilon))*Gt* &
    (kpc*r0)**((2*epsilon)/(1 + epsilon))*Md(r0))/ &
    (epsilon*G**((2.0d0*epsilon)/(2.0d0 + 2.0d0*epsilon))* &
    hbar**((2.0d0*epsilon)/(2.0d0 + 2.0d0*epsilon))*r0)

end function C1

function aND(r)

    implicit none

    real(8) :: aND
    real(8) :: r

    aND = (G*MSun*Md(r))/(kpc**2*r**2)

end function aND

function a0()

    implicit none

    real(8) :: a0

    a0 = (G*MSun*Md(r0))/(kpc**2*r0**2)

end function a0

function aDebye()

    implicit none

    real(8) :: aDebye

    aDebye = 3.0d0*a0()

end function aDebye

function rDebye()

    implicit none

    real(8) :: rDebye

    rDebye = r0/sqrt(3.0d0)

end function rDebye

function vModel(r)

    implicit none

    real(8) :: vModel
    real(8) :: r

    vModel = vDM(r)

end function vModel

end program MCMC
