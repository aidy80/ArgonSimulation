!Reads in a gro files and outputs its information as positions in another gro file
!Given random velocity distribution, this should result in a circle appear with randomly
!Distributed points

program drawVelDist
    integer, parameter :: numAtoms = 864
    integer, parameter :: numDimensions = 3
    real, dimension(numAtoms, numDimensions) :: vel
    
    real :: zero = 0.0
    real :: dim

    character :: nullChar
    character :: nullChar1
    integer :: nullInt
    integer :: nullInt1
    real, dimension(numDimensions) :: nullReals

    integer :: i,j,k

    open(unit=11, file='initArgon.gro')
    open(unit=91, file='velCircle.gro')

    read (11, *) nullChar
    read (11, *) nullInt

    do i=1, numAtoms
        read(11, 10) nullInt, nullChar, nullChar1, nullInt1, &
                    & (nullReals(j), j=1,3), (vel(i,k), k=1,3)
    end do

    read(11, *) dim

    write(91, *) "Velocity circle"
    write(91, 30) numAtoms

    do i=1, numAtoms
        write(91, 10) i, "ARGON", "AR", i, (vel(i,k), k=1,3), (nullReals(j), j=1,3) 
    end do

    write(91, 20) dim, dim, dim

    30 format(I5)
    20 format(F10.5, F10.5, F10.5)
    10 format(i5,2a5,i5,3f8.3,3f8.3)
end program drawVelDist
