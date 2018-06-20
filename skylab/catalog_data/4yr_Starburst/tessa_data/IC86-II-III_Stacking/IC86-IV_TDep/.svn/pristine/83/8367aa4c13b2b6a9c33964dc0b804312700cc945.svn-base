#include "llh/public/classes.h"
#include "llh/public/CoordEquatorialDeg.h"
#include "llh/public/SimpleEnergyProb.h"
#include "llh/public/DecBkgProb.h"
#include "fluxus/public/FluxFunction.h"
#include "llh/public/I3Analysis.h"
#include "llh/public/MultiAnalysisSet.h"
#include "llh/public/EventLoader.h"
#include "llhTimeDep/public/LocalCoordBkgProb.h"

class Ark {
    public:
    double livetime;
    double tmin;
    double tmax;
    AnalysisSet* psData;
    EquatorialDeg mySrcLocation;

    Ark(): livetime(0) { }
    virtual ~Ark() {
        if (psData) { delete psData; };
    }  // take care of psData in non-abstract derived classes

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux){};

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf){};

    // just some user-friendly alternatives

    void SetPointSource(double raDeg, double decDeg, double fluxConstant, double index) {
        PowerLawFlux flux(fluxConstant, index);
        printf("FluxConstant: %lg   PowerLaw Index: %lg   (GeV^-1 cm^-2 s^-1)\n", flux.GetSpectralIndex(), flux.GetFluxConstant() );
        SetPointSource(EquatorialDeg(raDeg, decDeg), flux);
    }

    void SetPointSource(double raDeg, double decDeg, const char* formula) {
        FormulaFlux flux(formula);
        printf("Formula: %s\n  where x is GeV, and flux is GeV^-1 cm^-2 s^-1\n", flux.GetFormula() );
        SetPointSource(EquatorialDeg(raDeg, decDeg), flux);
    }
  
    void SetPointSource(double raDeg, double decDeg, const char* formula, TimePdf * tPdf) {
        FormulaFlux flux(formula);
        printf("Formula: %s\n  where x is GeV, and flux is GeV^-1 cm^-2 s^-1\n", flux.GetFormula() );
        SetPointSource(EquatorialDeg(raDeg, decDeg), flux, tPdf);
    }

    void SetSource(SourceModule *src) {
        (dynamic_cast<I3Analysis*> (psData))->SetSource(*src);
    }
};



class MultiArk : public Ark {
    public:
    Ark* arkPtrArray[10];  // 10 is good for now!
    int n;

    MultiArk() {
        psData = new MultiAnalysisSet();
        for (int i=0; i<10; ++i) { 
        arkPtrArray[i] = NULL;
        }
        n = 0;
    }
    ~MultiArk() { 
        for (int i=0; i<10; ++i) { 
            if (arkPtrArray[i] != NULL ) {break;}
            delete arkPtrArray[i];
        }
    }
  
    void AddArk(Ark& ark) {
        arkPtrArray[n] = &ark;
        (dynamic_cast<MultiAnalysisSet*> (psData))->AddAnalysisSet(ark.psData);
        ++n;
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux) {
        mySrcLocation = srcCoord;
        int i=0;
        while (arkPtrArray[i]) {
            arkPtrArray[i]->SetPointSource(srcCoord, flux);
            ++i;
        }
    }
  
    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf) {
        mySrcLocation = srcCoord;
        int i=0;
        while (arkPtrArray[i]) {
            arkPtrArray[i]->SetPointSource(srcCoord, flux, tPdf);
            ++i;
        }
    }  

    void SetSource(SourceModule *src) {
        int i=0;
        while (arkPtrArray[i]) {
            arkPtrArray[i]->SetSource(src);
            ++i;
        }
    }  
};

class I3Ark : public Ark {
    public:
    EventLoader evLoader;
    vector<I3Event> baseEvents;
    SimpleEnergyProb* eProb;
    DecBkgProb decBkgProb;
    LocalCoordBkgProb lcBkgProb;

    double sourceZenWidthDeg;
    // This is the +/- zenith range in degrees to select MC events from,
    // and then rotate to desired source location.
    // (Smaller range is more accurate, but trade-off is lower statistics
    // for signal simulation.  You have to pick something appropriate for the 
    // statistics you have.)

    I3Ark() : Ark() ,eProb(NULL) ,  sourceZenWidthDeg(0)   { }

    ~I3Ark() {  }


    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";
    
        vector<I3Event> sourceEvents;
        evLoader.LoadSourceEvents(sourceEvents, mySrcLocation);
        I3PointGenerator i3point(sourceEvents, flux, mySrcLocation, livetime);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";
    
        vector<I3Event> sourceEvents;
        evLoader.LoadSourceEvents(sourceEvents, mySrcLocation);
        TimePdf * tPdf1 = tPdf->Clone();
        tPdf1->CheckTimeBounds(tmin,tmax);
        I3PointGenerator i3point(sourceEvents, flux, mySrcLocation, livetime, tPdf1);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }
};

