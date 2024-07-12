// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <map>
#include <memory>
#include <string>

class TEfficiency;
namespace ActsFatras {
class Particle;
}  // namespace ActsFatras

namespace ActsExamples {

// Tools to make efficiency plots to show tracking efficiency.
// For the moment, the efficiency is taken as the fraction of successfully
// smoothed track over all tracks
class EffPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning("#eta", 160, -4., 4.)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("pT [GeV/c]", 40, 0, 100)}};
  };

  /// @brief Nested Cache struct
  struct EffPlotCache {
    TEfficiency* trackEff_vs_pT{nullptr};   ///< Tracking efficiency vs pT
    TEfficiency* trackEff_vs_eta{nullptr};  ///< Tracking efficiency vs eta
    TEfficiency* trackEff_vs_phi{nullptr};  ///< Tracking efficiency vs phi

    TEfficiency* trackEff_vs_pT_eta0{
        nullptr};  ///< Tracking efficiencies vs pT at fixed eta = 0. - 0.2
    TEfficiency* trackEff_vs_pT_eta22{
        nullptr};  ///< Tracking efficiencies vs pT at fixed eta = 2.2 - 2.4
    TEfficiency* trackEff_vs_pT_eta36{
        nullptr};  ///< Tracking efficiencies vs pT at fixed eta = 3.6 - 3.8
    TEfficiency* trackEff_vs_pT_eta38{
        nullptr};  ///< Tracking efficiencies vs pT at fixed eta = 3.8 - 4.

    TEfficiency* trackEff_vs_pT_IrisL0hit{
        nullptr};  ///< Tracking efficiency vs pT + hit in the IRIS L0
    TEfficiency* trackEff_vs_eta_IrisL0hit{
        nullptr};  ///< Tracking efficiency vs eta + hit in the IRIS L0
    TEfficiency* trackEff_vs_phi_IrisL0hit{
        nullptr};  ///< Tracking efficiency vs phi + hit in the IRIS L0

    TEfficiency* trackEff_vs_mult_eta_pt{
        nullptr};  ///< Tracking efficiency vs mult, eta, pt

    TEfficiency* trackEff_vs_pT_el{
        nullptr};  ///< Tracking efficiency vs pT (electron)
    TEfficiency* trackEff_vs_eta_el{
        nullptr};  ///< Tracking efficiency vs eta (electron)
    TEfficiency* trackEff_vs_phi_el{
        nullptr};  ///< Tracking efficiency vs phi (electron)

    TEfficiency* trackEff_vs_pT_mu{
        nullptr};  ///< Tracking efficiency vs pT (muon)
    TEfficiency* trackEff_vs_eta_mu{
        nullptr};  ///< Tracking efficiency vs eta (muon)
    TEfficiency* trackEff_vs_phi_mu{
        nullptr};  ///< Tracking efficiency vs phi (muon)

    TEfficiency* trackEff_vs_pT_pi{
        nullptr};  ///< Tracking efficiency vs pT (pion)
    TEfficiency* trackEff_vs_eta_pi{
        nullptr};  ///< Tracking efficiency vs eta (pion)
    TEfficiency* trackEff_vs_phi_pi{
        nullptr};  ///< Tracking efficiency vs phi (pion)
    TEfficiency* trackEff_vs_eta_pi_pT_100MeV{
        nullptr};  ///< Tracking efficiency vs eta (pion), pT = 100 MeV/c
    TEfficiency* trackEff_vs_eta_pi_pT_1GeV{
        nullptr};  ///< Tracking efficiency vs eta (pion), pT = 1 GeV/c

    TEfficiency* trackEff_vs_pT_ka{
        nullptr};  ///< Tracking efficiency vs pT (kaon)
    TEfficiency* trackEff_vs_eta_ka{
        nullptr};  ///< Tracking efficiency vs eta (kaon)
    TEfficiency* trackEff_vs_phi_ka{
        nullptr};  ///< Tracking efficiency vs phi (kaon)

    TEfficiency* trackEff_vs_pT_pr{
        nullptr};  ///< Tracking efficiency vs pT (proton)
    TEfficiency* trackEff_vs_eta_pr{
        nullptr};  ///< Tracking efficiency vs eta (proton)
    TEfficiency* trackEff_vs_phi_pr{
        nullptr};  ///< Tracking efficiency vs phi (proton)
    TEfficiency* trackEff_vs_eta_pr_pT_100MeV{
        nullptr};  ///< Tracking efficiency vs eta (proton), pT = 100 MeV/c
    TEfficiency* trackEff_vs_eta_pr_pT_1GeV{
        nullptr};  ///< Tracking efficiency vs eta (proton), pT = 1 GeV/c
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  EffPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the efficiency plots
  ///
  /// @param effPlotCache the cache for efficiency plots
  void book(EffPlotCache& effPlotCache) const;

  /// @brief fill efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fill(EffPlotCache& effPlotCache,
            const ActsFatras::Particle& truthParticle, bool status) const;

  /// @brief fill efficiency plots + hit in IRIS L0
  ///
  /// @param effPlotCache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fillEffIrisHit(EffPlotCache& effPlotCache,
                      const ActsFatras::Particle& truthParticle,
                      bool status) const;

  /// @brief fill efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fill3DEff(EffPlotCache& effPlotCache,
                 const ActsFatras::Particle& truthParticle,
                 const float& eventMultLog, bool status) const;

  /// @brief write the efficiency plots to file
  ///
  /// @param effPlotCache cache object for efficiency plots
  void write(const EffPlotCache& effPlotCache) const;

  /// @brief delete the efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  void clear(EffPlotCache& effPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
