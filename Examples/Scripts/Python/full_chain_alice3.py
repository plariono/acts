#!/usr/bin/env python3
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    SeedfinderConfigArg,
    SeedingAlgorithm,
    ParticleSmearingSigmas,
    addCKFTracks,
    CKFPerformanceConfig,
)
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)
import pathlib
import acts
import acts.examples
import alice3

u = acts.UnitConstants
geo_dir = pathlib.Path.cwd().parent
outputDir = pathlib.Path("/Users/lrnv/alice/tools/acts/acts/bin/output/python/ckf_PbPb_100ev_optuna")

detector, trackingGeometry, decorators = alice3.buildALICE3Geometry(
    geo_dir, True, False)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)


s = acts.examples.Sequencer(events=100, numThreads=-1)
s = addPythia8(
    s,
    npileup=1,
    beam=acts.PdgParticle.eLead,
    cmsEnergy=5 * acts.UnitConstants.TeV,
    vtxGen=acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0.0125 * u.mm, 0.0125 *
                            u.mm, 55.5 * u.mm, 5.0 * u.ns),
        mean=acts.Vector4(0, 0, 0, 0),
    ),
    rnd=rnd,
    outputDirRoot=outputDir,
)
s = addFatras(
    s,
    trackingGeometry,
    field,
    ParticleSelectorConfig(
        eta=(0.0, 4.0), pt=(500 * u.MeV, None), removeNeutral=False),
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "acts/bin/alice3-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(0.5 * u.GeV, None), eta=(0, 4.0), nHits=(9, None)),
    SeedfinderConfigArg(
        r=(None, 200 * u.mm),
        deltaR=(0.44 * u.mm, 88.01 * u.mm),
        collisionRegion=(-250 * u.mm, 250 * u.mm),
        z=(-2000 * u.mm, 2000 * u.mm),
        maxSeedsPerSpM=1,
        sigmaScattering=49.795,
        radLengthPerSeed=0.09,
        minPt=500 * u.MeV,
        bFieldInZ=1.99724 * u.T,
        impactMax=23.25 * u.mm,
        cotThetaMax=27.2899,
    ),
    geoSelectionConfigFile=geo_dir /
    "acts/bin/geoSelection-alice3-cfg9.json",
    outputDirRoot=outputDir,
)
s = addCKFTracks(
    s,
    trackingGeometry,
    field,
    CKFPerformanceConfig(ptMin=500.0 * u.MeV, nMeasurementsMin=6),
    outputDirRoot=outputDir,
)

s.run()

# Config that results in high efficiencies (from seeding.py)
# SeedfinderConfigArg(
#     r=(None, 200 * u.mm),
#     deltaR=(1. * u.mm, 60 * u.mm),
#     collisionRegion=(-250 * u.mm, 250 * u.mm),
#     z=(-2000 * u.mm, 2000 * u.mm),
#     maxSeedsPerSpM=1,
#     sigmaScattering=50,
#     radLengthPerSeed=0.1,
#     minPt=500 * u.MeV,
#     bFieldInZ=1.99724 * u.T,
#     impactMax=3 * u.mm,
#     cotThetaMax=27.2899,
# )