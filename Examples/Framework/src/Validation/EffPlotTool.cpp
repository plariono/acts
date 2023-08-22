// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/EffPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <TEfficiency.h>

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

  const int nBinsPt = 37;
  const double xBinsPt[nBinsPt + 1] = {
      0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3,
      0.4,  0.5,  0.6,  0.7,  0.8,  0.9, 1.,   2.,   3.,   4.,   5.,  6.,   7.,
      8.,   9.,   10.,  20.,  30.,  40., 50.,  60.,  70.,  80.,  90., 100.};

  const int nBinsPtReduced = 28;
  const double xBinsPtReduced[nBinsPtReduced + 1] = {
      0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18,
      0.2,  0.25, 0.3,  0.4,  0.5,  0.6, 0.7,  0.8,  0.9,  1.,
      2.,   3.,   4.,   5.,   6.,   7.,  8.,   9.,   10.};

  // efficiency vs pT
  effPlotCache.trackEff_vs_pT = new TEfficiency(
      "trackeff_vs_pT", "Tracking efficiency;Truth pT [GeV/c];Efficiency",
      nBinsPt, xBinsPt);
  // efficiency vs eta
  effPlotCache.trackEff_vs_eta = PlotHelpers::bookEff(
      "trackeff_vs_eta", "Tracking efficiency;Truth #eta;Efficiency", bEta);
  // efficiency vs phi
  effPlotCache.trackEff_vs_phi = PlotHelpers::bookEff(
      "trackeff_vs_phi", "Tracking efficiency;Truth #phi;Efficiency", bPhi);

  // efficiency at fixed eta
  effPlotCache.trackEff_vs_pT_eta0 = new TEfficiency(
      "trackeff_vs_pT_eta0",
      "Tracking efficiency at #it{#eta} = [0.0 - 0.2];Truth #it{p}_{T} "
      "[GeV/#it{c}];Efficiency",
      nBinsPtReduced, xBinsPtReduced);

  effPlotCache.trackEff_vs_pT_eta22 = new TEfficiency(
      "trackeff_vs_pT_eta22",
      "Tracking efficiency at #it{#eta} = [2.2 - 2.4];Truth #it{p}_{T} "
      "[GeV/#it{c}];Efficiency",
      nBinsPtReduced, xBinsPtReduced);

  effPlotCache.trackEff_vs_pT_eta36 = new TEfficiency(
      "trackeff_vs_pT_eta36",
      "Tracking efficiency at #it{#eta} = [3.6 - 3.8];Truth #it{p}_{T} "
      "[GeV/#it{c}];Efficiency",
      nBinsPtReduced, xBinsPtReduced);

  // efficiency vs pT + hit in IRIS L0
  effPlotCache.trackEff_vs_pT_IrisL0hit = new TEfficiency(
      "trackeff_vs_pT_IrisL0hit",
      "Tracking efficiency;Truth pT [GeV/c];Efficiency", nBinsPt, xBinsPt);
  // efficiency vs eta + hit in IRIS L0
  effPlotCache.trackEff_vs_eta_IrisL0hit = PlotHelpers::bookEff(
      "trackeff_vs_eta_IrisL0hit",
      "Tracking efficiency + IRIS L0 hit;Truth #eta;Efficiency", bEta);
  // efficiency vs phi + hit in IRIS L0
  effPlotCache.trackEff_vs_phi_IrisL0hit = PlotHelpers::bookEff(
      "trackeff_vs_phi_IrisL0hit",
      "Tracking efficiency + IRIS L0 hit;Truth #phi;Efficiency", bPhi);

  // 3D efficiency in multiplicity, eta, pT
  effPlotCache.trackEff_vs_mult_eta_pt =
      new TEfficiency("trackeff_vs_mult_eta_pt",
                      "Tracking efficiency;Multiplicity;Truth #it{p}_{T} "
                      "[GeV/#it{c}];#it{#eta};Efficiency",
                      20, 0.5, 3.5, 40, -4., 4., 44, -2., 1.);
}

void ActsExamples::EffPlotTool::clear(EffPlotCache& effPlotCache) const {
  delete effPlotCache.trackEff_vs_pT;
  delete effPlotCache.trackEff_vs_eta;
  delete effPlotCache.trackEff_vs_phi;

  delete effPlotCache.trackEff_vs_pT_eta0;
  delete effPlotCache.trackEff_vs_pT_eta22;
  delete effPlotCache.trackEff_vs_pT_eta36;

  delete effPlotCache.trackEff_vs_pT_IrisL0hit;
  delete effPlotCache.trackEff_vs_eta_IrisL0hit;
  delete effPlotCache.trackEff_vs_phi_IrisL0hit;

  delete effPlotCache.trackEff_vs_mult_eta_pt;
}

void ActsExamples::EffPlotTool::write(
    const EffPlotTool::EffPlotCache& effPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  effPlotCache.trackEff_vs_pT->Write();
  effPlotCache.trackEff_vs_eta->Write();
  effPlotCache.trackEff_vs_phi->Write();

  effPlotCache.trackEff_vs_pT_eta0->Write();
  effPlotCache.trackEff_vs_pT_eta22->Write();
  effPlotCache.trackEff_vs_pT_eta36->Write();

  effPlotCache.trackEff_vs_pT_IrisL0hit->Write();
  effPlotCache.trackEff_vs_eta_IrisL0hit->Write();
  effPlotCache.trackEff_vs_phi_IrisL0hit->Write();

  effPlotCache.trackEff_vs_mult_eta_pt->Write();
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

  if (t_eta >= 0. && t_eta < 0.2)
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_eta0, t_pT, status);
  else if (t_eta >= 2.2 && t_eta < 2.4)
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_eta22, t_pT, status);
  else if (t_eta >= 3.6 && t_eta < 3.8)
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_eta36, t_pT, status);
}

void ActsExamples::EffPlotTool::fillEffIrisHit(
    EffPlotTool::EffPlotCache& effPlotCache,
    const ActsFatras::Particle& truthParticle, bool status) const {
  const auto t_phi = phi(truthParticle.unitDirection());
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_IrisL0hit, t_pT, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta_IrisL0hit, t_eta, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_phi_IrisL0hit, t_phi, status);
}

void ActsExamples::EffPlotTool::fill3DEff(
    EffPlotTool::EffPlotCache& effPlotCache,
    const ActsFatras::Particle& truthParticle, const float& eventMultLog,
    bool status) const {
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();
  const auto t_pT_Log = std::log10(t_pT);

  PlotHelpers::fillEff(effPlotCache.trackEff_vs_mult_eta_pt, eventMultLog,
                       t_eta, t_pT_Log, status);
}