#include "TVector3.h"
#include "TLorentzVector.h"

void rotate_Lvector()
{
    TLorentzVector v1(100, 5, 5, 100);
    TVector3 b(1, 0, 0);

    // Rotation around axes
    cout << "v1 vector" << endl;
    cout << "Fisrt component: " << v1.X() << endl;
    cout << "Second component: " << v1.Y() << endl;
    cout << "Third component: " << v1.Z() << endl;
    cout << "4th component: " << v1.T() << endl << endl;


    cout << "Rotate v1 90 degree around (1, 0, 0)" << endl;
    v1.Rotate(TMath::Pi() / 2, b); // 기준점 앞으로 움직이는 계에서 바라본 결과를 원할 때는 -를 붙인다. ROOT6의 Boost() 함수의 정의에서.
    cout << "Fisrt component: " << v1.X() << endl;
    cout << "Second component: " << v1.Y() << endl;
    cout << "Third component: " << v1.Z() << endl;
    cout << "4th component: " << v1.T() << endl << endl;
}