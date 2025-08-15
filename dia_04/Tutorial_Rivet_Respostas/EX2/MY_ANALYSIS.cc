// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {

  /// Distribuições de pares de léptons e jatos - Exercício 2
  class MY_ANALYSIS : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MY_ANALYSIS);

    /// @name Analysis methods
    /// @{

    /// Reserva os histogramas e inicializa as projeções
    void init() {
      const Cut leptonCuts = (Cuts::abspid==PID::ELECTRON || Cuts::abspid==PID::MUON) && Cuts::pT > 10*GeV && Cuts::abseta < 2.5;
      const DileptonFinder dileptons(91.2*GeV, 0.1, leptonCuts, Cuts::massIn(80*GeV,100*GeV));
      declare(dileptons, "dileptons");

      const FinalState fs(Cuts::abseta<4.9);
      const FastJets jets(fs, JetAlg::ANTIKT, 0.4, JetMuons::NONE, JetInvisibles::NONE);
      declare(jets, "jets");

      // Book histograms
      vector<double> mll_bins = { 66., 74., 78., 82., 84., 86., 88., 89., 90., 91.,
                                  92., 93., 94., 96., 98., 100., 104., 108., 116. };
      book(_h["mll"], "mass_ll", mll_bins);
      book(_h["njets"], "n_jets", {0,1,2,3,4,5});
      book(_h["jetPtSum"], "jetPtSum", {20, 35, 50, 65, 80, 95, 110});
    }

    /// Realiza a análise evento a evento 
    void analyze(const Event& event) {
      const Particles& bosons = apply<DileptonFinder>(event, "dileptons").bosons();
      if(bosons.size()!=1) vetoEvent;
      _h["mll"]->fill(bosons[0].mass());

      const Jets& jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT>10*GeV && Cuts::absrap<4.5);
      _h["njets"]->fill(jets.size());
      _h["jetPtSum"]->fill(sum(jets, Kin::pT, 0.0));
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

