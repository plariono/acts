// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// #include "TestGsfGeneric.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
//#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Generators/MultiplicityGenerators.hpp"
#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
#include "ActsExamples/Generators/VertexGenerators.hpp"
// #include "ActsExamples/Io/Csv/CsvPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjSpacePointWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
//#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"
//#include "ActsFatras/EventData/Barcode.hpp"

#include <chrono>
#include <iostream>
#include <random>


struct PrintFinalParticleStats : ActsExamples::BareAlgorithm {
  struct Config {
    std::string inParticles;
    std::string inSimulatedHits;
  } m_cfg;

  PrintFinalParticleStats(const Config &cfg, Acts::Logging::Level lvl)
      : ActsExamples::BareAlgorithm("PrintTrueParticleStats", lvl),
        m_cfg(cfg) {}

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext &ctx) const override {
    const auto &particles =
        ctx.eventStore.get<ActsExamples::SimParticleContainer>(
            m_cfg.inParticles);
    const auto &hits = ctx.eventStore.get<ActsExamples::SimHitContainer>(
        m_cfg.inSimulatedHits);

    std::size_t i = 0;
    for (const auto &part : particles) {
      ACTS_INFO("FinalParticle #" << i++ << " - pdg: " << part.pdg()
                                  << " momentum: " << part.absoluteMomentum());

      // Find the hits
      std::vector<ActsExamples::SimHit> traj_hits;
      std::copy_if(hits.begin(), hits.end(), std::back_inserter(traj_hits),
                   [&](const auto &hit) {
                     return hit.particleId() == part.particleId();
                   });

      // Print additional information about momentum if hits are present
      if (traj_hits.size() > 0) {
        std::sort(
            traj_hits.begin(), traj_hits.end(),
            [](const auto &a, const auto &b) { return a.time() < b.time(); });

        // Estimate path lengths
        std::vector<double> approx_path_lengths;
        approx_path_lengths.push_back(
            (traj_hits.front().position() - part.position()).norm());
        for (auto it = std::next(traj_hits.cbegin()); it != traj_hits.cend();
             ++it) {
          double last_part =
              (std::prev(it)->position() - it->position()).norm();
          approx_path_lengths.push_back(approx_path_lengths.back() + last_part);
        }

        // Print the stuff
        for (auto j = 0ul; j < traj_hits.size(); ++j) {
          ACTS_INFO("  pos: "
                    << traj_hits[j].position().transpose()
                    << " pathLength: " << approx_path_lengths.at(j) << " mom: "
                    << traj_hits[j].momentum4After().segment<3>(0).norm());
        }
      }

      // For the parser
      ACTS_INFO("End FinalParticle");
    }

    return ActsExamples::ProcessCode::SUCCESS;
  }
};
