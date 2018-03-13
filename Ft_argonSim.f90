!Aidan Fike
!March 2, 2017
!Program to simulated argon atoms in a cube based on lennard jones iteractions.
!Program will simulate atoms for 1us with 100,000 .01ps timesteps

program argonSim
    implicit none

    integer, parameter :: dp = kind(1.0)!Give reals double precision

    integer :: numAtoms !Number of atoms in this simulation

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :), allocatable :: pos, vel !pos: [m], vel: [m/s]

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :), allocatable :: force ![N] 

    integer, parameter :: numDimensions = 3 !Dimension of position 
                                            !and velocity space

    real(dp) :: Epot ![J] Potential of the entire system
    real(dp) :: totEnergy ![J] Total energy of the system
    real(dp) :: potential ![J] the total potential energy of the system
    real(dp) :: vSQ ![(m/s)^2] Square velocity of a given atom
    real(dp) :: vOld ![m/s] Temp variable for velocity
    real(dp) :: Fmag ![N] Used to hold the magnitude of force btwn two atoms
    real(dp) :: kineticEnergy ![J] The total kinetic energy of the system

    !Used to record the distance between two atoms
    real(dp), dimension(numDimensions) :: distance
    real(dp) :: distanceSq

    !Iterators
    integer :: i
    integer :: j
    integer :: k
    integer :: m
    integer :: p

    !Constants
    real(dp), parameter :: arMass = 6.633521311332984E-26 ![kg] Mass of argon
    real(dp), parameter :: Bolz = 1.38064852E-23 ![J/K] Boltzmann constant
    real(dp) :: epsilon ![J] Minimum in the Lennard Jones equation
    real(dp), parameter :: sigma = 34E-11 ![m]. Constant in the Lennard 
                                          !Jones equation
    real(dp) :: sigmaDistTwo![]. Sigma^2. Used for optimizational purposes
    real(dp) :: sigmaDistSix ![]. Sigma^6. Used for optimizational purposes
    real(dp) :: sigmaDistTwelve![]. Sigma^12. Used for optimizational purposes
    real(dp), parameter :: dim = 3.4700E-9 ![m] Size of each wall
                                       !of the cubic enclosure
    real(dp), parameter :: temperature = 90.0 ![K] Const temperature of the system
    real :: cutoff ![m] Cutoff distance for short-range interaction
    real, parameter :: timeStep = 1.0E-14 ![s] Time between calculations
    real(dp) :: fourEps ![J^4] Epsilon**4. Used for optimization
    real(dp) :: cutoffSq ![m^2] The cutoff radius squared
    real(dp) :: twentyFourEps ![J^24] Epsilon**24. Used for optimization
    integer, parameter :: numSteps = 10000 !Number of timesteps in the program

    character :: nullChar
    character :: nullChar1
    integer :: nullInt
    integer :: nullInt1

    !Initialization of epsilon and cutoff radius
    epsilon = 120.0 * Bolz
    fourEps = 4.0 * epsilon
    twentyFourEps = 24.0 * epsilon
    cutoff = 2.25 * sigma
    cutoffSq = cutoff**2

    !Open file containing position and velocity information 
    !about argon atoms in the cubic boundary from argon.gro
    open(unit=1, file='argon.gro') 
    
    !Read in header information from the file
    read(1, *) nullChar
    read(1, 30) numAtoms

    !Allocate space on the heap for position, velocity, and force information
    allocate(pos(numAtoms, numDimensions))
    allocate(vel(numAtoms, numDimensions))
    allocate(force(numAtoms, numDimensions))

    !Read in position and velocity information from the .gro file
    do m = 1, numAtoms
        read(1, 10) nullInt, nullChar, nullChar1, nullInt1, pos(m,1), &
                            & pos(m,2), pos(m,3), vel(m,1), vel(m,2), vel(m,3)
    end do

    !Convert read-in pos/velocity information from nm and nm/ps to m and m/s
    do i = 1, numAtoms
        do k = 1, numDimensions
            pos(i,k) = pos(i,k) * 1.0E-9
            vel(i,k) = vel(i,k) * 1.0E3
        end do
    end do

    close (unit=1)

    !call testLJ(sigma, epsilon)

    do p = 1, numSteps
        !Init force vectors to zero
        call initZero(force, numAtoms, numDimensions)

        !Adjust the velocities in the system such that the net velocity 
        !in each direction is zero. This prevents wandering ice-cube problem
        call zeroNetMomentum(vel, numAtoms, numDimensions)

        Epot = 0.0
        kineticEnergy = 0.0

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

                sigmaDistTwo = (sigma**2 / distanceSq)
                sigmaDistSix = sigmaDistTwo * sigmaDistTwo * sigmaDistTwo
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

                !If the distance between the two atoms is below the cutoff, 
                !calculate the force exerted on each of them based on the 
                !lennard jones potential
                if (distanceSq < cutoffSq) then
                    do k = 1, numDimensions
                        force(i, k) = force(i, k) + &
                                            &Fmag * (distance(k) / distanceSq)
                        force(j, k) = force(j,k) - &
                                            &Fmag * (distance(k) / distanceSq)
                    end do
                end if

                !if ((i == 1).AND.(j == 2)) then
                !    print *, "FORCE", force(1,1), force(1,2), force(1, 3)
                !end if
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

        !Find the kinetic energy of the system and 
        !the total energy of the system
        kineticEnergy = arMass * vSQ * 0.5
        totEnergy = Epot + kineticEnergy

        !print*, "Timestep:", p, "Total Energy [J]: ",totEnergy, "Potential Energy [J]: ",Epot,&
        !                                 &"Kinetic Energy [J]: ", kineticEnergy

        !Scale the velocities of the system such 
        !that the temperature is at 90 degrees
        call scaleTemp(temperature, kineticEnergy, vel, &
                                                &numAtoms, numDimensions, Bolz)
    end do

    !Free all heap memory
    deallocate(pos)
    deallocate(vel)
    deallocate(force)

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

        integer :: i
        integer :: k

        do i = 1, arraySize 
            do k = 1, numDimensions
                array(i,k) = 0.0
            end do
        end do

        return

    end subroutine initZero

    !Zeros the net magnitude of the vectors passed in. Used to prevent the 
    !Wandering ice cube problem
    subroutine zeroNetMomentum(array, arraySize, numDimensions)
        implicit none
        
        integer, intent(in) :: arraySize
        integer, intent(in) :: numDimensions
        real(dp), dimension(arraySize, numDimensions), intent(inout) :: array

        integer :: i
        integer :: k
    
        real(dp), dimension(numDimensions) :: netVel
        real(dp), dimension(numDimensions) :: velScale

        !Find the current velocity of the system
        do k = 1, numDimensions
            netVel(k) = 0.0
        end do

        do i = 1, arraySize
            do k = 1, numDimensions
                netVel(k) = netVel(k) + array(i, k) 
            end do
        end do

        !Calc how much to scale each velocity by such that net velocity is zero
        do k = 1, numDimensions
            velScale(k) = netVel(k) / real(arraySize)
        end do

        !Scale the velocity of each atom in the system
        do i = 1, arraySize
            do k = 1, numDimensions
                array(i, k) = array(i, k) - velScale(k)
            end do
        end do

        !TESTING
        !do k = 1, numDimensions
        !    netVel(k) = 0.0
        !end do
        
        !do i = 1, arraySize
        !    do k = 1, numDimensions
        !        netVel(k) = netVel(k) + array(i, k)         
        !    end do
        !end do
        !print *, "NETVELS", netVel(1), netVel(2), netVel(3)
        !print *, "random velocity:", array(60, 1), array(60, 2), array(60, 3)
    end subroutine zeroNetMomentum

    !Scale the velocities in the system such that the system has a 
    !temperature equal to desiredTemperature
    subroutine scaleTemp(desiredTemp, KE, array, arraySize, numDimensions,Bolz)
        implicit none 
    
        real(dp), intent(in) :: KE
        real(dp), intent(in) :: Bolz
        real(dp), intent(in) :: desiredTemp
        integer, intent(in) :: numDimensions
        integer, intent(in) :: arraySize
        real(dp), dimension(arraySize, numDimensions), intent(inout) :: array

        real(dp) :: currTemp 
        real(dp) :: tempScale
        integer :: degreesFreedom = 3

        integer :: i
        integer :: k

        real(dp) :: nextVSq

        !Use equipartition thm to find the current temperature of the system
        currTemp = KE * (2.0D0 / (degreesFreedom * arraySize * Bolz))
        !print *,"currTemp [K]:", currTemp
        tempScale = sqrt(desiredTemp / currTemp)

        !Scale each velocity such that the kintic energy of the system can 
        !be applied to the equipartition theorem to find that the system has 
        !a temperature of "desiredTemperature"
        do i=1, arraySize
            do k=1, numDimensions
                array(i, k) = array(i, k) * tempScale
            end do
        end do

        !TESTING
        !nextVSQ = 0.0
        !do i=1,arraySize
        !    do k=1, numDimensions
        !        nextVSQ = nextVSQ + (array(i,k))**2
        !    end do
        !end do

        !currTemp = arMass * nextVSQ * 0.5 * (2.0D0 / (degreesFreedom * arraySize * Bolz))
        !print *, "CORRECTED TEMP:", currTemp

    end subroutine scaleTemp

    subroutine testLJ(sigma, epsilon)
        implicit none

        real(dp), intent(in) :: sigma
        real(dp), intent(in) :: epsilon

        integer, parameter :: numDistances = 100000000
        real, parameter :: distStep = 1.0E-14
        real(dp) :: twoSix
        real(dp) :: theoMin
        real(dp) :: distanceSq
        real(dp) :: sigmaDistTwo
        real(dp) :: sigmaDistSix
        real(dp) :: sigmaDistTwelve

        real(dp) :: minPotential = 10000
        real(dp) :: minDistance

        integer :: i

        do i = 1, numDistances
            distanceSq = (i * distStep)**2
            sigmaDistTwo = (sigma**2 / distanceSq)
            sigmaDistSix = sigmaDistTwo * sigmaDistTwo * sigmaDistTwo
            sigmaDistTwelve = sigmaDistSix**2

            potential = 4.0 * epsilon * (sigmaDistTwelve - sigmaDistSix)
            if (potential < minPotential) then
                minPotential = potential
                minDistance = i * distStep
            end if
        end do

        twoSix = 2.0**(1.0/6.0)
        theoMin = twoSix * sigma
        
        print *, "Theoretical minimum distance, sigma * 2^(1/6):", theoMin
        print *, "Theoretical minimum potential, -epsilon:", epsilon * (-1.0)
        print *, "Measured (sigma, epsilon): (", minDistance, minPotential, ")"
    end subroutine
end program argonSim
