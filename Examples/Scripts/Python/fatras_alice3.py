#!/usr/bin/env python3
from pathlib import Path

import acts
import acts.examples
from alice3 import buildALICE3Geometry
import argparse
from acts.examples.simulation import addParticleGun, addFatras, EtaConfig, MomentumConfig, ParticleConfig, addDigitization

u = acts.UnitConstants


def runFatras(trackingGeometry, field, outputDir, smearing, s: acts.examples.Sequencer = None):

    s = s or acts.examples.Sequencer(events=1, numThreads=-1)
    s.config.logLevel = acts.logging.INFO
    rnd = acts.examples.RandomNumbers()
    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 1.0 * u.GeV, transverse=True),
        EtaConfig(-4.0, 4.0, uniform=True),
        ParticleConfig(10, acts.PdgParticle.eMuon, randomizeCharge=False),
        rnd=rnd,
    )
    outputDir = Path(outputDir)

    addFatras(
        s,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    if smearing:
        s = addDigitization(
            s,
            trackingGeometry,
            field,
            outputDirCsv=outputDir / "csv",
            outputDirRoot=outputDir,
            digiConfigFile=Path.cwd().parent / "acts/bin/alice3-smearing-config.json",
            rnd=rnd,
        )

    return s


if "__main__" == __name__:

    p = argparse.ArgumentParser(
        description="Script to run Fatras with the ALICE3 geometry."
    )
    p.add_argument(
        "--geo_dir",
        default=Path.cwd().parent,
        type=Path,
        help="Input directory containing the ALICE3 standalone geometry.",
    )
    p.add_argument(
        "--output_dir",
        default=Path(
            "/Users/lrnvmbp14/alice/actsdir/acts/bin/output/fatras"),
        type=Path,
        help="Directory to write outputs to",
    )
    p.add_argument(
        "--no-material", action="store_true", help="Decorate material to the geometry"
    )
    p.add_argument(
        "--no-smearing", action="store_true", help="Run without digitization"
    )

    args = p.parse_args()
    args.output_dir.mkdir(exist_ok=True, parents=True)

    geo_example_dir = Path(args.geo_dir)
    output_dir = Path(args.output_dir)

    assert geo_example_dir.exists(), "Detector example input directory missing"

    detector, trackingGeometry, decorators = buildALICE3Geometry(
        geo_example_dir,
        material=not args.no_material,
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runFatras(trackingGeometry, field, output_dir,
              smearing=not args.no_smearing).run()
