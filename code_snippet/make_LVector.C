#include "TLorentzVector.h"

void make_LVector()
{
    // (x, y, z, t) or (px, py, pz, E)
    TLorentzVector v1(2, 3, 5, 0);
    TLorentzVector v2(1, 1, 1, 1);

    cout << "Fisrt component: " << v2.X() << endl;
    cout << "Second component: " << v2.Y() << endl;
    cout << "Third component: " << v2.Z() << endl;
    cout << "Third component: " << v2.T() << endl;

    // Actually, v2.X() = v2.Px(), v2.T() = v2.E() etc.

    cout << endl << "Dot product btw v1 and v2: " << v1.Dot(v2) << endl;
    cout << "Angle btw v1 and v2 (radian): " << v1.Angle(v2.Vect()) << endl << endl;


    // ### Set X, Y, Z ###
    v1.SetX(4);
    v1.SetY(-1);
    v1.SetZ(3);
    cout << "new Fisrt component of v1: " << v1.X() << endl;
    cout << "Second component: " << v1.Y() << endl;
    cout << "Third component: " << v1.Z() << endl << endl;

    auto v3 = v2.Vect();
    cout << "Vect() returns vecotr (x, y, z)" << endl;
    cout << "Fisrt component: " << v3.X() << endl;
    cout << "Second component: " << v3.Y() << endl;
    cout << "Third component: " << v3.Z() << endl;
}