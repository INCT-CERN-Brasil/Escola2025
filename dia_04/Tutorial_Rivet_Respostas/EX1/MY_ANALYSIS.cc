// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {

  /// @brief Distribuição de para de léptons - Exercício 1
  class MY_ANALYSIS : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MY_ANALYSIS);

    /// @name Analysis methods
    /// @{

    /// Reserva os histogramas e inicializa as projeções
    void init() {
      const Cut leptonCuts = (Cuts::abspid==PID::ELECTRON || Cuts::abspid==PID::MUON) && Cuts::pT > 10*GeV && Cuts::abseta < 2.5;
      const DileptonFinder dileptons(91.2*GeV, 0.1, leptonCuts, Cuts::OPEN); // do not set mass window for now
      declare(dileptons, "dileptons");

      // Reserva os  histogramas
      vector<double> mll_bins = { 66., 74., 78., 82., 84., 86., 88., 89., 90., 91.,
                                  92., 93., 94., 96., 98., 100., 104., 108., 116. };
      book(_h["mll"], "mass_ll", mll_bins);
    }

    /// Realiza a análise evento a evento 
    void analyze(const Event& event) {
      const Particles& bosons = apply<DileptonFinder>(event, "dileptons").bosons();
      if(bosons.size()!=1) vetoEvent;
      _h["mll"]->fill(bosons[0].mass());
    }

    /// Normaliza os  histogramas
    void finalize() {
      const double scaleFactor = crossSection() / sumOfWeights();
      scale(_h, scaleFactor);
    }

    /// @}

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    /// @}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MY_ANALYSIS);
}

