program force1
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)

    real(dp) :: force
    real(dp) :: twelve
    real(dp) :: twelveNum
    real(dp) :: twelveOne
    real(dp) :: twelveSix
    real(dp) :: twelveDenom
    real(dp) :: twelveDenomOne
    real(dp) :: six
    real(dp) :: sixDenom
    real(dp) :: sixNum
    real(dp) :: forceCoeff

    twelveOne = 3.4E-10
    twelveSix = twelveOne**6
    twelveNum = twelveSix**2
    twelveDenomOne = 3.86E-10
    twelveDenom = twelveDenomOne**13
    twelve = 12.0 * twelveNum / twelveDenom
    print *, "Twelve", twelve
    sixDenom = 3.86E-10
    sixNum = 3.4E-10
    six = 6 * sixNum**6 / sixDenom**7
    print *, "Six", six
    forceCoeff = 4.0 * 120.0 * 1.38064E-23
    print *, "ForceCoeff", forceCoeff

    force = forceCoeff * (twelve - six)

    print *, "Force", force
    
end program force1
