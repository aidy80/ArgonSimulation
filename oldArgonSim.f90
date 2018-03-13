!Aidan Fike
!March 2, 2017
!Program to simulated argon atoms in a cube based on lennard jones iteractions.
!Program will simulate atoms for 1us with 100,000 .01ps timesteps

program argonSim
    implicit none

    integer, parameter :: dp = kind(1.d0)!Give reals double precision

    integer :: numAtoms !Number of atoms in this simulation

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :), allocatable :: pos, vel !pos: [m], vel: [m/s]

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :), allocatable :: force ![N] 

    integer, parameter :: numDimensions = 3 !Dimension of position 
                                            !and velocity space

    real(dp) :: Epot ![J] Potential of the entire system
    real(dp) :: totEnergy ![J] Total energy of the system
    real(dp) :: potential ![J] contains the total potential energy of the system
    real(dp) :: vSQ ![(m/s)^2] Square velocity of a given atom
    real(dp) :: vOld ![m/s] Temp variable for velocity
    real(dp) :: Fmag ![N] Used to hold the magnitude of force btwn two atoms
    real(dp) :: kineticEnergy ![J] The total kinetic energy of the system

    !Used to record the distance between two atoms
    real(dp), dimension(numDimensions) :: distance
    real(dp) :: distanceSq

    !Iterators
    integer :: x
    integer :: y
    integer :: z
    integer :: i
    integer :: j

    !Constants
    real, parameter :: arMass = 6.633521311332984E-26 ![kg] Mass of argon
    real, parameter :: Bolz = 1.38064852E-23 ![J/K] Boltzmann constant
    real :: epsilon ![J] Minimum in the Lennard Jones equation
    real,parameter :: sigma = 34E-11 ![m]. Constant in the Lennard Jones equation
    real(dp) :: sigmaDistTwo![]. Sigma^2. Used for optimizational purposes
    real(dp) :: sigmaDistSix ![]. Sigma^6. Used for optimizational purposes
    real(dp) :: sigmaDistTwelve![]. Sigma^12. Used for optimizational purposes
    real,parameter :: dim = 3.4700E-9 ![m] Size of each wall of the cubic enclosure
    real,parameter :: temperature = 90.0 ![K] Const temperature of the system
    real :: cutoff ![m] Cutoff distance for short-range interaction
    real,parameter :: timeStep = 1.0E-14 ![s] Time between calculations
    integer,parameter :: numSteps = 10000 !Number of timesteps in the program

    character :: nullChar
    character :: nullChar1
    integer :: nullInt
    integer :: nullInt1

    !Definition of epsilon and cutoff radius
    epsilon = 120.0 * Bolz
    cutoff = 2.25 * sigma

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
    do x = 1, numAtoms
        read(1, 10) nullInt, nullChar, nullChar1, nullInt1, pos(x,1), &
                            & pos(x,2), pos(x,3), vel(x,1), vel(x,2), vel(x,3)
    end do

    !Convert read-in pos/velocity information from nm and nm/ps to m and m/s
    do i = 1, numAtoms
        do j = 1, numDimensions
            pos(i,j) = pos(i,j) * 1.0E-9
            vel(i,j) = vel(i,j) * 1.0E3
        end do
    end do

    close (unit=1)

    do x = 1, numSteps
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
                do y = 1, numDimensions
                    distance(y) = pos(i, y) - pos(j, y)
                    distance(y) = distance(y) - dim * anint(distance(y) / dim)
                end do

                distanceSq = distance(1)**2 + distance(2)**2 + distance(3)**2

                sigmaDistTwo = (sigma**2 / distanceSq)
                sigmaDistSix = sigmaDistTwo *sigmaDistTwo *sigmaDistTwo
                sigmaDistTwelve = sigmaDistSix**2

                !Calc potential from lennard jones equation between the 
                !current pair of atoms and add it to the current 
                !potential energy sum
                potential = 4.0 * epsilon * (sigmaDistTwelve - sigmaDistSix)
                Epot = Epot + potential

                !Calculate the resulting force on the current two atoms 
                !based on the lennard jones potential between them. Calculated 
                !using the negative gradient
                Fmag = 4.0 * epsilon * (12 * sigmaDistTwelve - 6 * sigmaDistSix)

                !If the distance between the two atoms is below the cutoff, 
                !calculate the force exerted on each of them based on the 
                !lennard jones potential
                if (distanceSq < cutoff**2) then
                    do y = 1, numDimensions
                        force(i, y) = force(i,y) + &
                                            &Fmag * (distance(y) / distanceSq)
                        force(j, y) = force(j,y) - &
                                            &Fmag * (distance(y) / distanceSq)
                    end do
                end if
            end do
        end do


        !Use the leap-frog verlet algorithm to calculate new position and 
        !velocity vectors for all atoms based on the forces 
        !calculated between them.
        vSQ = 0
        do y = 1, numAtoms
            do z = 1, numDimensions 
                vOld = vel(y,z)
                vel(y, z) = vel(y, z) + (force(y, z) / arMass) * (timeStep) 
                pos(y, z) = pos(y, z) + vel(y, z) * timeStep
                vSQ = vSQ + ((vOld + vel(y,z)/2.0D0))**2
            end do
        end do

        !Find the kinetic energy of the system and 
        !the total energy of the system[
        kineticEnergy = arMass * vSQ / 2.0D0
        totEnergy = Epot + kineticEnergy

        print*, "Total Energy [J]: ",totEnergy, "Potential Energy [J]: ",Epot,&
                                         &"Kinetic Energy [J]: ", kineticEnergy

        !Scale the velocities of the system such 
        !that the temperature is at 90 degrees
        call scaleTemp(temperature, kineticEnergy, vel,numAtoms, numDimensions)
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

        integer :: y
        integer :: z

        do y = 1, arraySize 
            do z = 1, numDimensions
                array(y,z) = 0.0
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
        integer :: j
    
        real(dp), dimension(numDimensions) :: netVel
        real(dp), dimension(numDimensions) :: velScale

        !Find the current velocity of the system
        do i = 1, numDimensions
            netVel(i) = 0.0
        end do

        do i = 1, arraySize
            do j = 1, numDimensions
                netVel(j) = netVel(j) + array(i, j) 
            end do
        end do

        !Calc how much to scale each velocity by such that net velocity is zero
        do i = 1, numDimensions
            velScale(i) = netVel(i) / real(arraySize)
        end do

        !Scale the velocity of each atom in the system
        do i = 1, arraySize
            do j = 1, numDimensions
                array(i, j) = array(i, j) - velScale(j)
            end do
        end do

        !TESTING
        do i = 1, numDimensions
            netVel(i) = 0
        end do

        do i=1,arraySize
            do j=1,numDimensions
                netVel(j) = netVel(j) + array(i, j)
            end do
        end do

        !print *, "Netvels:", netVel(1), netVel(2), netVel(3)

    end subroutine zeroNetMomentum

    !Scale the velocities in the system such that the system has a 
    !temperature equal to desiredTemperature
    subroutine scaleTemp(desiredTemperature, KE, array, arraySize, numDimensions)
        implicit none 
    
        real(dp), intent(in) :: KE
        real, intent(in) :: desiredTemperature
        integer, intent(in) :: numDimensions
        integer, intent(in) :: arraySize
        real(dp), dimension(arraySize, numDimensions),intent(inout) :: array

        real(dp) :: currTemp 
        real(dp) :: tempScale
        integer :: degreesFreedom = 3

        !Use equipartition thm to find the current temperature of the system
        currTemp = KE * (2.0D0 / (degreesFreedom * arraySize * Bolz))
        print *,"currTemp [K]:", currTemp
        tempScale = sqrt(desiredTemperature / currTemp)

        !Scale each velocity such that the kintic energy of the system can 
        !be applied to the equipartition theorem to find that the system has 
        !a temperature of "desiredTemperature"
        do i=1,arraySize
            do j=1,numDimensions
                array(i, j) = array(i, j) * tempScale
            end do
        end do
    end subroutine scaleTemp
end program argonSim
