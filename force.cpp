#include <iostream>
#include<math.h>

using namespace std;

int main() {
    double twelve;
    double twelveSix;
    double twelveOne;
    double twelveNum;
    double twelveDenomOne;
    double twelveDenom;
    double sixDenom;
    double sixNum;
    double six;
    double forceCoeff;
    double force;

    twelveOne = 3.4E-10;
    twelveSix = pow(twelveOne,6);
    twelveNum = pow(twelveSix,2);
    twelveDenomOne = 3.86E-10;
    twelveDenom = pow(twelveDenomOne,13);
    twelve = 12.0 * twelveNum / twelveDenom;
    cout <<  "Twelve" << twelve;
    sixDenom = 3.86E-10;
    sixNum = 3.4E-10;
    six = 6 * pow(sixNum,6) / pow(sixDenom,7);
    cout <<  "Six" << six;
    forceCoeff = 4.0 * 120.0 * 1.38064E-23;
    cout << "ForceCoeff" << forceCoeff;

    force = forceCoeff * (twelve - six);

    cout << "Force" << force;

    return 0;
} 
