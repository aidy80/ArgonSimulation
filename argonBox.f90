!Written by Aidan Fike
!February 24, 2018
!
!This program outputs a gromacs-formatted-file filled with 864 evenly-space 
!argon atoms at 94.4 Kelvin. These atoms exist in a cube  of dimensions 
!3.47x3.47x3.47nm.

program argonBox 
    implicit none

    !Initial vrms of atoms in the simulation nm/ps
    real :: vrms = 0.24278

    !Number of atoms being simulated
    integer, parameter :: numAtoms = 864
    integer :: currNumAtoms = numAtoms

    !Number of atoms per side of the cube
    integer :: xNum = 10
    integer :: yNum = 10
    integer :: zNum = 9

    !Dimension of the cube in nanometers
    real :: dim = 3.47

    !Used to hold temporary variables for velocity values
    real :: a
    real :: b
    real :: c
    real :: mag !Used to normalize the randomly-created normal vector

    !Vectors to hold velocity and position information about a given vector
    real, dimension(3) :: pos, vel

    !Counters to look through all possible atoms
    integer :: x
    integer :: y
    integer :: z
    integer :: i

    !Relation between curr iterator (x,y,z) to dimension
    real :: dmdx
    real :: dmdy 
    real :: dmdz 

    !Create a random seed for random number generation. 
    !This seed will be modifed upon each iteration using scale, increment, max
    integer :: seed = 1010

    !If you want to use rand() 
    call random_seed(size = seed)
    call srand(seed)

    dmdx = dim / xNum
    dmdy = dim / yNum
    dmdz = dim / zNum

    !Open a file for writing argon atom information formatted for gromacs
    open(unit = 1, file='initArgon.gro')

    !The file's header information
    write(1, *) "MD of 864 argon atoms, t = 0.0ns"
    write(1, 30) numAtoms

    !Treating the box as a 3D grid, loop through it and create equally 
    !spaced atoms throughout.
    do x = 0, xNum - 1
        do y = 0, yNum - 1
            do z = 0, zNum - 1
                !Because this loop will go through more than numAtoms 
                !number of atoms, stop outputting atom information when you 
                !have outputted numAtoms of information
                if (currNumAtoms == 0) then
                    goto 100
                end if

                !Create a velocity vector in a random direction with 
                !magnitude vrms (nm/ps)
                do i = 1, 3
                    vel(i) = vrms * (rand() * 2 - 1)
                end do

                !Normalize the velocity vector
                mag = sqrt(vel(1)**2 + vel(2)**2 + vel(3)**2)
                vel(1) = (real(vel(1)) / mag) * vrms
                vel(2) = (real(vel(2)) / mag) * vrms
                vel(3) = (real(vel(3)) / mag) * vrms

                !Create a position vector depending on 
                !the current values of x,y, and z 
                !(the current coordinated of the atom in the grid)
                pos(1) = real(x) * dmdx
                pos(2) = real(y) * dmdy
                pos(3) = real(z) * dmdz

                !Output the position and velocity information for the 
                !most-recently-created atom
                write(1, 10) abs(numAtoms + 1 - currNumAtoms), "ARGON", "AR", &
                             & 865 - currNumAtoms, pos(1), pos(2), pos(3), vel(1), vel(2), vel(3)

                !Decrement the number of atoms left to output
                currNumAtoms = currNumAtoms - 1 
            end do
        end do
    end do


    100 continue !Continue here after the do loops

    !Output the box's dimensions
    write(1, 20) dim, dim, dim

    !Format information for the outputted file
    30 format(I5)
    20 format(F10.5, F10.5, F10.5)
    10 format(i5,2a5,i5,3f8.3,3f8.3)

    close(unit = 1)

end program argonBox
