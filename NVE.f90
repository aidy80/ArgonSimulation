!Aidan Fike
!March 2, 2017
!Program to simulated argon atoms in a cube based on lennard jones iteractions.
!Program will simulate atoms for 1us with 100,000 .01ps timesteps under NVE conditions

program nveSim
    implicit none

    integer, parameter :: dp = kind(1.0)!Give reals double precision

    integer :: numAtoms !Number of atoms in this simulation

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :), allocatable :: pos, vel !pos: [m], vel: [m/s]

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :), allocatable :: force ![N] 
    real(dp), dimension(:), allocatable :: distBins ![]

    integer, parameter :: numDimensions = 3 !Dimension of position 
    !and velocity space
    integer, parameter :: CvvStep = 300 !Number of items in Cvv
    integer, parameter :: MSDStep = 300 !Number of items in MSD

    real(dp) :: Epot ![J] Potential of the entire system
    real(dp) :: totEnergy ![J] Total energy of the system
    real(dp) :: potential ![J] the total potential energy of the system
    real(dp) :: vSQ ![(m/s)^2] Square velocity of a given atom
    real(dp) :: vOld ![m/s] Temp variable for velocity
    real(dp) :: Fmag ![N] Used to hold the magnitude of force btwn two atoms
    real(dp) :: kineticEnergy ![J] The total kinetic energy of the system

    real :: start_time
    real :: end_time

    !Used to record the distance between two atoms
    real(dp), dimension(numDimensions) :: distance
    real(dp) :: distanceSq
    real(dp) :: distanceMag

    real(dp), dimension(CvvStep) :: Cvv
    real(dp) :: CvvOne = 0
    real(dp), dimension(MSDStep) :: MSD
    real(dp), dimension(:, :, :), allocatable :: vStore
    real(dp), dimension(:, :, :), allocatable :: pStore
    integer :: nCorr = 0
    integer :: nCorr_MSD = 0
    real(dp) :: diffusionCoeff = 0

    !Iterators
    integer :: i, j, k, m, l, n, p

    !Constants
    real(dp), parameter :: arMass = 6.633521311332984E-26 ![kg] Mass of argon
    real(dp), parameter :: Bolz = 1.38064852E-23 ![J/K] Boltzmann constant
    real(dp), parameter :: epsilon = 120.0 * Bolz![J] Minimum in the 
                                                 !Lennard Jones equation
    real(dp), parameter :: sigma = 34.0E-11 ![m]. Constant in the Lennard 
                                            !Jones equation
    real(dp), parameter :: sigmaSq = sigma**2
    real(dp) :: sigmaDistTwo![]. Used for optimizational purposes
    real(dp) :: sigmaDistSix ![]. Used for optimizational purposes
    real(dp) :: sigmaDistTwelve![]. Used for optimizational purposes
    real(dp), parameter :: dim = 3.4700E-9![m] Size of each wall of the cubic enclosure
    real,parameter :: cutoff = 2.25 * sigma ![m] Cutoff distance for 
                                            !short-range interaction
    real, parameter :: timeStep = 1.0E-14 ![s] Time between calculations
    real(dp) :: fourEps = 4.0 * epsilon ![J] Epsilon*4. Used for optimization
    real(dp) :: cutoffSq = cutoff**2 ![m^2] The cutoff radius squared
    real(dp) :: twentyFourEps = 24.0 * epsilon ![J] Epsilon*24. 
                                               !Used for optimization
    integer, parameter :: numSteps = 24000 !Number of timesteps in the program

    integer, parameter :: zeroMomentTimeStep = 100 !Number of timesteps 
                                                    !between momentum-zeroing
    integer, parameter :: numTrajSteps = 10 !Number of timesteps between 
                                            !trajectory outputs to .ttr file

    !Set up data for finding g(r)
    real(dp), parameter :: delR = 0.003e-9
    integer :: numBins = ANint((dim*1.733)/(2 * delR)) + 1
    integer :: numFullBins = ANint(dim/(2*delR))
    integer, dimension(:), allocatable :: bins
    integer :: currIndex
    real(dp) :: rho 
    real(dp), parameter :: pi = 3.1415
    real(dp) :: sphereConst
    real(dp) :: rLower 
    real(dp) :: rUpper
    real(dp) :: shellVol

    character :: nullChar
    character :: nullChar1
    integer :: nullInt
    integer :: nullInt1

    CALL cpu_time(start_time)

    !Open file containing position and velocity information 
    !about argon atoms in the cubic boundary from argon.gro
    open(unit=11, file='NVT_final.gro') 
    open(unit=91, file='NVE_totEnergy.dat')
    open(unit=92, file='NVE_potentialEnergy.dat')
    open(unit=93, file='NVE_kineticEnergy.dat')
    open(unit=94, file='NVE.gro')
    open(unit=95, file='NVE_final.gro')
    open(unit=96, file='TCF.dat')
    open(unit=97, file='gr.dat')
    open(unit=98, file='MSD.dat')

    !Read in header information from the file
    read(11, *) nullChar
    read(11, 30) numAtoms

    !Allocate space on the heap for position, velocity, and force information
    allocate(pos(numAtoms, numDimensions))
    allocate(vel(numAtoms, numDimensions))
    allocate(force(numAtoms, numDimensions))
    allocate(vStore(numAtoms, numDimensions, CvvStep))
    allocate(pStore(numAtoms, numDimensions, MSDStep))
    allocate(bins(numBins))

    rho = numAtoms / (dim ** 3)
    sphereConst = 4.0D0 * pi * rho / 3.0D0

    !Read in position and velocity information from the .gro file
    do m = 1, numAtoms
        read(11, 10) nullInt, nullChar, nullChar1, nullInt1, pos(m,1), &
            & pos(m,2), pos(m,3), vel(m,1), vel(m,2), vel(m,3)
    end do

    !Convert read-in pos/velocity information from nm and nm/ps to m and m/s
    do i = 1, numAtoms
        do k = 1, numDimensions
            pos(i,k) = pos(i,k) * 1.0E-9
            vel(i,k) = vel(i,k) * 1.0E3
        end do
    end do

    close (unit=11)

    !Adjust the velocities in the system such that the net velocity 
    !in each direction is zero. This prevents wandering ice-cube problem
    call zeroNetMomentum(vel, numAtoms, numDimensions)

    do i = 1, CvvStep
        Cvv(i) = 0.0D0
    end do

    do i = 1, MSDStep
        MSD(i) = 0.0D0
    end do

    do i = 1, numBins
        bins(i) = 0
    end do

    do p = 1, numSteps
        !Init force vectors to zero
        call initZero(force, numAtoms, numDimensions)

        !Adjust the velocities in the system such that the net velocity 
        !in each direction is zero. This prevents wandering ice-cube problem
        if (mod(p, zeroMomentTimeStep) == 0) then
            call zeroNetMomentum(vel, numAtoms, numDimensions)
        end if

        Epot = 0.0D0
        kineticEnergy = 0.0D0

        !Compute Force between each pair of atoms if their distance is below
        !the cutoff radius. Each pair is evaluated only once
        do i = 1, numAtoms - 1
            do j = i + 1, numAtoms

                !Calculate the distance between the current atoms, applying 
                !periodic boundary conditions to find the minimum possible 
                !distance for each coordinate
                do k = 1, numDimensions
                    distance(k) = pos(i, k) - pos(j, k)
                    distance(k) = distance(k) - dim * anint(distance(k) / dim)
                end do

                distanceSq = distance(1)**2 + distance(2)**2 + distance(3)**2
                
                !If the distance between the two atoms is below the cutoff, 
                !calculate the force exerted on each of them based on the 
                !lennard jones potential
                if (distanceSq < cutoffSq) then
                    sigmaDistTwo = sigmaSq / distanceSq
                    sigmaDistSix = sigmaDistTwo**3 
                    sigmaDistTwelve = sigmaDistSix**2

                    !Calc potential from lennard jones equation between the 
                    !current pair of atoms and add it to the current 
                    !potential energy sum
                    potential = fourEps * (sigmaDistTwelve - sigmaDistSix)
                    Epot = Epot + potential

                    !Calculate the resulting force on the current two atoms 
                    !based on the lennard jones potential between them. Calculated 
                    !using the negative gradient
                    Fmag = twentyFourEps * (2 * sigmaDistTwelve - sigmaDistSix)

                    do k = 1, numDimensions
                        force(i, k) = force(i, k) + &
                                            &Fmag * (distance(k) / distanceSq)
                        force(j, k) = force(j,k) - &
                                            &Fmag * (distance(k) / distanceSq)
                    end do
                end if


                currIndex = Int(sqrt(distanceSq)/delR) + 1
                bins(currIndex) = bins(currIndex) + 2
            end do
        end do

        !Use the leap-frog verlet algorithm to calculate new position and 
        !velocity vectors for all atoms based on the forces 
        !calculated between them.
        vSQ = 0.0
        do i = 1, numAtoms
            do k = 1, numDimensions 
                vOld = vel(i,k)
                vel(i, k) = vel(i, k) + (force(i, k) / arMass) * (timeStep) 
                pos(i, k) = pos(i, k) + vel(i, k) * timeStep
                vSQ = vSQ + ((vOld + vel(i,k)) * 0.5)**2
            end do
        end do

        !Call TCF to update Cvv
        call TCF(p, Cvv, CvvStep, nCorr, vStore, vel)
        call Calc_MSD(p, MSD, MSDStep, nCorr_MSD, pStore, pos)

        !Find the kinetic energy of the system and 
        !the total energy of the system
        kineticEnergy = arMass * vSQ * 0.5
        !print *, "Curr Temp:", kineticEnergy * (2.0D0 / (3 * numAtoms * Bolz))
        totEnergy = Epot + kineticEnergy

        write(91, *) (timestep * real(p)), totEnergy
        write(92, *) (timestep * real(p)), Epot
        write(93, *) (timestep * real(p)), kineticEnergy

        !Write .ttr file 
        if (mod(p, numTrajSteps) == 0) then
            write(94, *) "Trajectory file for NVE ensemble. 1ns total time"
            write(94, 30) numAtoms
            do n = 1, numAtoms 
                write(94, 10) n, "ARGON", "AR", n, (pos(n, 1) * 1E9) - 1E9 * dim &
                            &* floor(pos(n,1)/dim), (pos(n, 2) * 1E9) - 1E9 * dim &
                            &* floor(pos(n,2)/dim), pos(n, 3) * 1E9 - 1E9 * dim &
                            &* floor(pos(n,3)/dim), vel(n, 1) * 1E-3, vel(n, 2) &
                            &* 1E-3, vel(n, 3) * 1E-3

            end do
            write(94, 20) dim * 1E9, dim * 1E9, dim * 1E9
        end if

        !Write final frame in .gro format
        if (p == numAtoms) then 
            write(95, *) "Trajectory file for last frame of NVE ensemble."
            write(95, 30) numAtoms
            do n = 1, numAtoms
                write(95, 10) n, "ARGON", "AR", n, &
                            &(pos(n, 1) * 1E9), (pos(n, 2) * 1E9), (pos(n, 3) * 1E9), &
                            &vel(n, 1) * 1E-3, vel(n, 2) * 1E-3, vel(n, 3) * 1E-3
            end do
            write(95, 20) dim * 1E9, dim * 1E9, dim * 1E9
        end if
    end do

    !Output the TCF graph and calculate the diffusion coefficient from it
    CvvOne = Cvv(1)/real(nCorr)
    do m = 1, CvvStep
        if (m < 200) then
            write(96, *) (real(m) * timestep) * 1e12, &
                         &(Cvv(m)/real(nCorr)) / CvvOne
        end if
        diffusionCoeff = diffusionCoeff + (Cvv(m) / real(nCorr)) * timestep
    end do
    print *, "DiffusionCoefficient: ", diffusionCoeff * 1e4

    !Output the MSD graph
    do m = 1, MSDStep
        write(98, *) (real(m) * timestep) * 1e12, &
                                     &(MSD(m)/real(nCorr_MSD)) * 1e20
    end do

    !Output the g(r) graph
    do m = 1, numFullBins
        rLower = real(m - 1) * delR    
        rUpper = rLower + delR
        shellVol = sphereConst * (rUpper ** 3 - rLower ** 3)
        write(97, *) (rlower + delR * 0.5) / sigma, real(bins(m)) / (real(numSteps) * &
                                                    &real(numAtoms) * shellVol)
    end do

    !Free all heap memory
    deallocate(pos)
    deallocate(vel)
    deallocate(force)
    deallocate(vStore)
    deallocate(pStore)
    deallocate(bins)

    close (unit=91)
    close (unit=92)
    close (unit=93)
    close (unit=94)
    close (unit=95)
    close (unit=96)
    close (unit=97)
    close (unit=98)

    CALL cpu_time(end_time)
    print *, "Time usage:", end_time - start_time, " seconds"

    30 format(I5)
    20 format(F10.5, F10.5, F10.5)
    10 format(i5,2a5,i5,3f8.3,3f8.3)

contains
    !Helper subroutine to initialize a 2D array of dimension 
    !(arraySize, numDimensions) with all zeros. "array" is the 2D array passed 
    !and returned with all zeros
    subroutine initZero(array, arraySize, numDimensions)
        implicit none
        
        integer, intent(in) :: arraySize
        integer, intent(in) :: numDimensions
        real(dp), dimension(arraySize, numDimensions), intent(inout) :: array

        integer :: i, k

        do i = 1, arraySize 
            do k = 1, numDimensions
                array(i,k) = 0.0D0
            end do
        end do

        return

    end subroutine initZero

    !Zeros the net magnitude of the vectors passed in. Used to prevent the 
    !Wandering ice cube problem
    subroutine zeroNetMomentum(vel, numAtoms, numDimensions)
        implicit none
        
        integer, intent(in) :: numAtoms
        integer, intent(in) :: numDimensions
        real(dp), dimension(numAtoms, numDimensions), intent(inout) :: vel


        integer :: i, k
    
        real(dp), dimension(numDimensions) :: netVel
        real(dp), dimension(numDimensions) :: velScale

        !Find the current velocity of the system
        do k = 1, numDimensions
            netVel(k) = 0.0
        end do

        do i = 1, numAtoms
            do k = 1, numDimensions
                netVel(k) = netVel(k) + vel(i, k) 
            end do
        end do

        !Calc how much to scale each velocity by such that net velocity is zero
        do k = 1, numDimensions
            velScale(k) = netVel(k) / real(numAtoms)
        end do

        !Scale the velocity of each atom in the system
        do i = 1, numAtoms
            do k = 1, numDimensions
                vel(i, k) = vel(i, k) - velScale(k)
            end do
        end do
    end subroutine zeroNetMomentum
    
    !Add the the Cvv array velocity information from the current timestep
    subroutine TCF(step, Cvv, CvvStep, nCorr, vStore, vel)
        implicit none

        integer, intent(in) :: step
        integer, intent(inout) :: nCorr
        integer, intent(in) :: CvvStep
        real(dp), dimension(CvvStep), intent(inout) :: Cvv
        real(dp), dimension(numAtoms, numDimensions, CvvStep), intent(inout) :: vStore
        real(dp), dimension(numAtoms, numDimensions), intent(in) :: vel

        real(dp), dimension(CvvStep) :: sum !Due to numerical percision reasons,
                                            !keep sum such that each dot 
                                            !product is successfully added, then
                                            !this sum is added to the Cvv eventually

        integer :: index0
        integer :: currIndex
        integer :: i, k, m

        do m = 1, CvvStep
            sum = 0.0
        end do

        currIndex = mod(step - 1, CvvStep) + 1 

        !Store record of the current velocities at the current index
        do i = 1, numAtoms
            do k = 1, numDimensions
                vStore(i, k, currIndex) = vel(i,k)
            end do
        end do

        !After all indexes in vStore have been filled with velocity information, 
        !Complete the Cvv algorithm to calculate all dot products with the current 
        !index0, dependent on the step #. 
        if (step.GE.CvvStep) then
            nCorr = nCorr + 1
            index0 = mod(step, CvvStep) + 1
            do i = 1, numAtoms
                do k = 1, numDimensions
                    do m = 1, CvvStep
                        currIndex = mod(m - 1 + step, CvvStep) + 1
                        sum(m) = sum(m) + vStore(i, k, index0) * vStore(i, k, currIndex)
                    end do
                end do
            end do

            !Add recorded sums to Cvv array
            do m = 1, CvvStep
                Cvv(m) = Cvv(m) + (sum(m) / (numAtoms * numDimensions))
            end do
        end if
    end subroutine TCF

    !Add the the MSD array position information from the current timestep
    subroutine Calc_MSD(step, MSD, MSDStep, nCorr_MSD, pStore, pos)
        implicit none

        integer, intent(in) :: step
        integer, intent(inout) :: nCorr_MSD
        integer, intent(in) :: MSDStep
        real(dp), dimension(MSDStep), intent(inout) :: MSD
        real(dp), dimension(numAtoms, numDimensions, MSDStep), intent(inout) :: pStore
        real(dp), dimension(numAtoms, numDimensions), intent(in) :: pos

        integer :: index0
        integer :: currIndex
        integer :: i, k, m

        real(dp), dimension(MSDStep) :: sum !Due to numerical percision reasons,
                                            !keep sum such that each dot 
                                            !product is successfully added, then
                                            !this sum is added to the Cvv eventually

        do m = 1, MSDStep 
            sum = 0.0
        end do

        currIndex = mod(step - 1, MSDStep) + 1 

        !Store record fo the current postion at the current index
        do i = 1, numAtoms
            do k = 1, numDimensions
                pStore(i, k, currIndex) = pos(i,k)
            end do
        end do

        !After all indexes in pStore have been fileed with position information,
        !Complete the MSD algorithm to calculate all dot products with the current
        !index0, dependent on the step #.
        if (step.GE.MSDStep) then
            nCorr_MSD = nCorr_MSD + 1
            index0 = mod(step, MSDStep) + 1
            do i = 1, numAtoms
                do k = 1, numDimensions
                    do m = 1, MSDStep
                        currIndex = mod(m - 1 + step, MSDStep) + 1
                        sum(m) = sum(m) + (pStore(i, k, index0) - pStore(i, k, currIndex))**2
                    end do
                end do
            end do
        end if

        !Add sums to the MSD array
        do m = 1, MSDStep
            MSD(m) = MSD(m) + sum(m) / numAtoms
        end do
    end subroutine Calc_MSD

end program nveSim
