#!/usr/bin/env python3
import pathlib, acts, acts.examples, alice3

u = acts.UnitConstants
geo_dir = pathlib.Path.cwd().parent
outputDir = pathlib.Path("/Users/lrnv/alice/tools/acts/acts/bin/python/output")

detector, trackingGeometry, decorators = alice3.buildALICE3Geometry(geo_dir, False, False)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

from particle_gun import addParticleGun, MomentumConfig, EtaConfig, ParticleConfig
from fatras import addFatras
from digitization import addDigitization
from seeding import addSeeding, TruthSeedRanges
from ckf_tracks import addCKFTracks, CKFPerformanceConfig


s = acts.examples.Sequencer(events=100, numThreads=-1)
s = addParticleGun(
    s,
    MomentumConfig(1.0 * u.GeV, 1.0 * u.GeV, True),
    EtaConfig(3.0, 4.0, True),
    ParticleConfig(1, acts.PdgParticle.eMuon, True),
    rnd=rnd,
)
s = addFatras(
    s,
    trackingGeometry,
    field,
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "acts/bin/python/config/alice3-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(3.0, 4.0), nHits=(9, None)),
    geoSelectionConfigFile=geo_dir /
    "acts/bin/python/config/geoSelection-alice3-cfg1.json",
    outputDirRoot=outputDir,
)
s = addCKFTracks(
    s,
    trackingGeometry,
    field,
    CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
    outputDirRoot=outputDir,
)

s.run()
