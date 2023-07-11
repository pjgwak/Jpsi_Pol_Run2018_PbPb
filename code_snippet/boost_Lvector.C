#include "TVector3.h"
#include "TLorentzVector.h"

void boost_Lvector()
{
    TLorentzVector v1(100, 0, 0, 100);
    TVector3 b(0.85, 0.5, 0);

    // Rotation around axes
    cout << "v1 vector" << endl;
    cout << "Fisrt component: " << v1.X() << endl;
    cout << "Second component: " << v1.Y() << endl;
    cout << "Third component: " << v1.Z() << endl;
    cout << "4th component: " << v1.T() << endl << endl;

    cout << "v1 boosted" << endl;
    //v1.Boost(b); // 기준점 뒤로 움직이는 계에서 바라본 속력을 원할 때
    // 뒤로 달려가는 입장에서는 내가 앞으로 더 빨리 다릴리는 것처럼 보인다.
    
    v1.Boost(-b); // 기준점 앞으로 움직이는 계에서 바라본 결과를 원할 때는 -를 붙인다. ROOT6의 Boost() 함수의 정의에서.
    cout << "Fisrt component: " << v1.X() << endl;
    cout << "Second component: " << v1.Y() << endl;
    cout << "Third component: " << v1.Z() << endl;
    cout << "4th component: " << v1.T() << endl << endl;

    b = v1.BoostVector();
    cout << "Fisrt component: " << b.X() << endl;
    cout << "Second component: " << b.Y() << endl;
    cout << "Third component: " << b.Z() << endl;

    // See below
    // https://cpp.hotexamples.com/examples/-/TLorentzVector/BoostVector/cpp-tlorentzvector-boostvector-method-examples.html
}