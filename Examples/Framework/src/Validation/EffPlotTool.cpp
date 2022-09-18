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

  effPlotCache.trackEff_vs_pT = new TEfficiency(
      "trackeff_vs_pT",
      "Tracking efficiency;Truth #it{p}_{T} [GeV/#it{c}];Efficiency", nBinsPt,
      xBinsPt);

  // efficiency vs pT
  // effPlotCache.trackEff_vs_pT = PlotHelpers::bookEff(
  //    "trackeff_vs_pT", "Tracking efficiency;Truth pT [GeV/c];Efficiency",
  //    bPt);
  // efficiency vs eta
  effPlotCache.trackEff_vs_eta = PlotHelpers::bookEff(
      "trackeff_vs_eta", "Tracking efficiency;Truth #eta;Efficiency", bEta);
  // efficiency vs phi
  effPlotCache.trackEff_vs_phi = PlotHelpers::bookEff(
      "trackeff_vs_phi", "Tracking efficiency;Truth #phi;Efficiency", bPhi);
}

void ActsExamples::EffPlotTool::clear(EffPlotCache& effPlotCache) const {
  delete effPlotCache.trackEff_vs_pT;
  delete effPlotCache.trackEff_vs_eta;
  delete effPlotCache.trackEff_vs_phi;
}

void ActsExamples::EffPlotTool::write(
    const EffPlotTool::EffPlotCache& effPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  effPlotCache.trackEff_vs_pT->Write();
  effPlotCache.trackEff_vs_eta->Write();
  effPlotCache.trackEff_vs_phi->Write();
}

void ActsExamples::EffPlotTool::fill(EffPlotTool::EffPlotCache& effPlotCache,
                                     const ActsFatras::Particle& truthParticle,
                                     bool status) const {
  const auto t_phi = phi(truthParticle.unitDirection());
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT, t_pT, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta, t_eta, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_phi, t_phi, status);
}
