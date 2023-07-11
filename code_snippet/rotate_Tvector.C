#include "TVector3.h"

void rotate_Tvector()
{
    TVector3 v1(1, 1, 1);
    TVector3 v2(1, 1, 1);
    TVector3 v3(1, 1, 0);

    // Rotation around axes
    cout << "v1 vector" << endl;
    cout << "Fisrt component: " << v1.X() << endl;
    cout << "Second component: " << v1.Y() << endl;
    cout << "Third component: " << v1.Z() << endl << endl;

    v1.RotateX(TMath::Pi() / 2);
    cout << "v1 rotated 90 degree around X axis" << endl;
    cout << "Fisrt component: " << v1.X() << endl;
    cout << "Second component: " << v1.Y() << endl;
    cout << "Third component: " << v1.Z() << endl << endl;


    // Rotation around other vector
    cout << "v2 vector" << endl;
    cout << "Fisrt component: " << v2.X() << endl;
    cout << "Second component: " << v2.Y() << endl;
    cout << "Third component: " << v2.Z() << endl << endl;

    v2.Rotate(TMath::Pi() / 2, v3);
    cout << "v2 rotated 90 degree around unit vector of (1, 1, 0)" << endl;
    cout << "Fisrt component: " << v2.X() << endl;
    cout << "Second component: " << v2.Y() << endl;
    cout << "Third component: " << v2.Z() << endl << endl;


}