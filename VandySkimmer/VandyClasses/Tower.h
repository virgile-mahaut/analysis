#include <iostream>
#include <limits>
#include <vector>

class Tower {
    public:
        Tower() = default;
        ~Tower();

        void CopyTo(Tower* tower);


        float px() { return Px; };
        void set_px(float tmpPx) { Px = tmpPx; };

        float py() { return Py; };
        void set_py(float tmpPy) { Py = tmpPy; };

        float pz() { return Pz; };
        void set_pz(float tmpPz) { Pz = tmpPz; };

        float e() { return E; };
        void set_e(float tmpE) { E = tmpE; };
        

        int get_calo() { return calo; };
        void set_calo(int tmpCalo) { calo = tmpCalo; };

    private:
        float Px;
        float Py;
        float Pz;
        float E;
        int calo; //0=EMCal, 1=retowered EMCal, 2=IHCal, 3=OHCal

        //ClassDef(Tower,1);

};
