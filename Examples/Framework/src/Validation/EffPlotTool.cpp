// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/EffPlotTool.hpp"

#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

ActsExamples::EffPlotTool::EffPlotTool(
    const ActsExamples::EffPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EffPlotTool", lvl)) {}

void ActsExamples::EffPlotTool::book(
    EffPlotTool::EffPlotCache& effPlotCache) const {
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  ACTS_DEBUG("Initialize the histograms for efficiency plots");

  const int nBinsPt = 24;
  const double xBinsPt[nBinsPt + 1] = {
      0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
      0.8,  0.9,  1.,   2.,   3.,   4.,  5.,   6.,  7.,  8.,  9.,  10.};

  const int nBinsP = 13;
  const double xBinsP[nBinsP + 1] = {1., 2., 3.,  4.,  5.,  6.,  7.,
                                     8., 9., 10., 20., 30., 40., 50.};

  // efficiency vs pT
  //   effPlotCache.trackEff_vs_pT = PlotHelpers::bookEff(
  //   "trackeff_vs_pT", "Tracking efficiency;Truth pT [GeV/c];Efficiency",
  //   bPt);
  effPlotCache.trackEff_vs_pT = new TEfficiency(
      "trackeff_vs_pT",
      "Tracking efficiency;Truth #it{p}_{T} [GeV/#it{c}];Efficiency", nBinsPt,
      xBinsPt);

  // efficiency vs p
  effPlotCache.trackEff_vs_p = new TEfficiency(
      "trackeff_vs_p",
      "Tracking efficiency;Truth #it{p} [GeV/#it{c}];Efficiency", nBinsP,
      xBinsP);

  // efficiency vs eta
  effPlotCache.trackEff_vs_eta = PlotHelpers::bookEff(
      "trackeff_vs_eta", "Tracking efficiency;Truth #eta;Efficiency", bEta);
  // efficiency vs phi
  effPlotCache.trackEff_vs_phi = PlotHelpers::bookEff(
      "trackeff_vs_phi", "Tracking efficiency;Truth #phi;Efficiency", bPhi);

  // electrons
  //   effPlotCache.trackEff_vs_pT_el = PlotHelpers::bookEff(
  //   "trackeff_vs_pT_el",
  //   "Tracking efficiency (e^{#pm});Truth pT [GeV/c];Efficiency", bPt);
  effPlotCache.trackEff_vs_pT_el = new TEfficiency(
      "trackeff_vs_pT_el",
      "Tracking efficiency (e^{#pm});Truth #it{p}_{T} [GeV/#it{c}];Efficiency",
      nBinsPt, xBinsPt);

  effPlotCache.trackEff_vs_eta_el = PlotHelpers::bookEff(
      "trackeff_vs_eta_el",
      "Tracking efficiency (e^{#pm});Truth #eta;Efficiency", bEta);
  effPlotCache.trackEff_vs_phi_el = PlotHelpers::bookEff(
      "trackeff_vs_phi_el",
      "Tracking efficiency (e^{#pm});Truth #phi;Efficiency", bPhi);

  // pions
  //   effPlotCache.trackEff_vs_pT_pi = PlotHelpers::bookEff(
  //   "trackeff_vs_pT_pi",
  //   "Tracking efficiency (#pi^{#pm});Truth pT [GeV/c];Efficiency", bPt);
  effPlotCache.trackEff_vs_pT_pi =
      new TEfficiency("trackeff_vs_pT_pi",
                      "Tracking efficiency (#pi^{#pm});Truth #it{p}_{T} "
                      "[GeV/#it{c}];Efficiency",
                      nBinsPt, xBinsPt);

  effPlotCache.trackEff_vs_eta_pi = PlotHelpers::bookEff(
      "trackeff_vs_eta_pi",
      "Tracking efficiency (#pi^{#pm});Truth #eta;Efficiency", bEta);
  effPlotCache.trackEff_vs_phi_pi = PlotHelpers::bookEff(
      "trackeff_vs_phi_pi",
      "Tracking efficiency (#pi^{#pm});Truth #phi;Efficiency", bPhi);

  // kaons
  //   effPlotCache.trackEff_vs_pT_Ka = PlotHelpers::bookEff(
  //   "trackeff_vs_pT_Ka",
  //   "Tracking efficiency (K^{#pm});Truth pT [GeV/c];Efficiency", bPt);

  effPlotCache.trackEff_vs_pT_Ka = new TEfficiency(
      "trackeff_vs_pT_Ka",
      "Tracking efficiency (K^{#pm});Truth pT [GeV/c];Efficiency", nBinsPt,
      xBinsPt);
  effPlotCache.trackEff_vs_eta_Ka = PlotHelpers::bookEff(
      "trackeff_vs_eta_Ka",
      "Tracking efficiency (K^{#pm});Truth #eta;Efficiency", bEta);
  effPlotCache.trackEff_vs_phi_Ka = PlotHelpers::bookEff(
      "trackeff_vs_phi_Ka",
      "Tracking efficiency (K^{#pm});Truth #phi;Efficiency", bPhi);

  // protons
  //   effPlotCache.trackEff_vs_pT_Pr = PlotHelpers::bookEff(
  //   "trackeff_vs_pT_Pr",
  //   "Tracking efficiency (p^{#pm});Truth pT [GeV/c];Efficiency", bPt);
  effPlotCache.trackEff_vs_pT_Pr = new TEfficiency(
      "trackeff_vs_pT_Pr",
      "Tracking efficiency (p^{#pm});Truth pT [GeV/c];Efficiency", nBinsPt,
      xBinsPt);
  effPlotCache.trackEff_vs_eta_Pr = PlotHelpers::bookEff(
      "trackeff_vs_eta_Pr",
      "Tracking efficiency (p^{#pm});Truth #eta;Efficiency", bEta);
  effPlotCache.trackEff_vs_phi_Pr = PlotHelpers::bookEff(
      "trackeff_vs_phi_Pr",
      "Tracking efficiency (p^{#pm});Truth #phi;Efficiency", bPhi);
}

void ActsExamples::EffPlotTool::clear(EffPlotCache& effPlotCache) const {
  delete effPlotCache.trackEff_vs_pT;
  delete effPlotCache.trackEff_vs_p;
  delete effPlotCache.trackEff_vs_eta;
  delete effPlotCache.trackEff_vs_phi;

  delete effPlotCache.trackEff_vs_pT_el;
  delete effPlotCache.trackEff_vs_eta_el;
  delete effPlotCache.trackEff_vs_phi_el;

  delete effPlotCache.trackEff_vs_pT_pi;
  delete effPlotCache.trackEff_vs_eta_pi;
  delete effPlotCache.trackEff_vs_phi_pi;

  delete effPlotCache.trackEff_vs_pT_Ka;
  delete effPlotCache.trackEff_vs_eta_Ka;
  delete effPlotCache.trackEff_vs_phi_Ka;

  delete effPlotCache.trackEff_vs_pT_Pr;
  delete effPlotCache.trackEff_vs_eta_Pr;
  delete effPlotCache.trackEff_vs_phi_Pr;
}

void ActsExamples::EffPlotTool::write(
    const EffPlotTool::EffPlotCache& effPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  effPlotCache.trackEff_vs_pT->Write();
  effPlotCache.trackEff_vs_p->Write();
  effPlotCache.trackEff_vs_eta->Write();
  effPlotCache.trackEff_vs_phi->Write();

  effPlotCache.trackEff_vs_pT_el->Write();
  effPlotCache.trackEff_vs_eta_el->Write();
  effPlotCache.trackEff_vs_phi_el->Write();

  effPlotCache.trackEff_vs_pT_pi->Write();
  effPlotCache.trackEff_vs_eta_pi->Write();
  effPlotCache.trackEff_vs_phi_pi->Write();

  effPlotCache.trackEff_vs_pT_Ka->Write();
  effPlotCache.trackEff_vs_eta_Ka->Write();
  effPlotCache.trackEff_vs_phi_Ka->Write();

  effPlotCache.trackEff_vs_pT_Pr->Write();
  effPlotCache.trackEff_vs_eta_Pr->Write();
  effPlotCache.trackEff_vs_phi_Pr->Write();
}

void ActsExamples::EffPlotTool::fill(EffPlotTool::EffPlotCache& effPlotCache,
                                     const ActsFatras::Particle& truthParticle,
                                     bool status) const {
  const auto t_phi = phi(truthParticle.unitDirection());
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();
  const auto t_p = truthParticle.absoluteMomentum()();

  PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT, t_pT, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_p, t_p, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta, t_eta, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_phi, t_phi, status);

  if (abs(truthParticle.pdg()) == 11) {
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_el, t_pT, status);
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta_el, t_eta, status);
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_phi_el, t_phi, status);
  } else if (abs(truthParticle.pdg()) == 211) {
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_pi, t_pT, status);
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta_pi, t_eta, status);
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_phi_pi, t_phi, status);
  } else if (abs(truthParticle.pdg()) == 321) {
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_Ka, t_pT, status);
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta_Ka, t_eta, status);
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_phi_Ka, t_phi, status);
  } else if (abs(truthParticle.pdg()) == 2212) {
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_Pr, t_pT, status);
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta_Pr, t_eta, status);
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_phi_Pr, t_phi, status);
  }
}
