#include "TVector3.h"

void make_Tvector()
{
    TVector3 v1(2, 3, 5);
    TVector3 v2(1, 1, 1);

    cout << "Fisrt component: " << v2.X() << endl;
    cout << "Second component: " << v2.Y() << endl;
    cout << "Third component: " << v2.Z() << endl;

    cout << endl << "Dot product btw v1 and v2: " << v1.Dot(v2) << endl;
    cout << "Angle btw v1 and v2 (radian): " << v1.Angle(v2) << endl << endl;


    // ### Set X, Y, Z ###
    v1.SetX(4);
    v1.SetY(-1);
    v1.SetZ(3);
    cout << "new Fisrt component of v1: " << v1.X() << endl;
    cout << "Second component: " << v1.Y() << endl;
    cout << "Third component: " << v1.Z() << endl;
}