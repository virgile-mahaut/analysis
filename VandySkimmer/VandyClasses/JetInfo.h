#include <iostream>
#include <limits>
#include <vector>

class JetInfo {
    public:
        JetInfo() = default;
        ~JetInfo();

        float px() { return Px; };
        void set_px( float tpx ) { Px = tpx; };
        
        float py() {return Py;};
        void set_py( float tpy ) { Py = tpy; };

        float pz() { return Pz; };
        void set_pz( float tpz ) { Pz = tpz; };

        float e() { return E; };
        void set_e( float te ) { E = te; };

        float pt() { return Pt; };
        void set_pt( float pt ) {Pt = pt;};

        float pt_uncalib() { return Pt_uncalib; };
        void set_pt_uncalib( float pt_uncalib ) {Pt_uncalib = pt_uncalib;};

        float hCaloFrac() { return HCaloFrac; };
        void set_hCaloFrac( float cf ) {HCaloFrac = cf;};

        std::vector<int> get_constituents() { return constituents; };
        void set_constituents( std::vector<int> jetCons ) { constituents = jetCons; };

	std::string parent() { return Parent; };
	void set_parent( std::string par ) { Parent = par; };

	bool isQuark() { return IsQuark; }; 
	void set_isQuark( bool iQ ) { IsQuark = iQ; };

	int parentPID() { return Parent_pID; };
	void set_parentPID( int pid ) { Parent_pID = pid; };

        void CopyTo(JetInfo *jet);
    private:
        float Px;
        float Py;
        float Pz;
        float E;
        float Pt;
        float Pt_uncalib;
        float HCaloFrac;
        std::vector<int> constituents;
	//these are only for truth jets
	std::string Parent; //quark (flavor) or gluon
	bool IsQuark; //quark jet
	int Parent_pID; //useful for offline flavor, mostly reduant with quark bool
};
