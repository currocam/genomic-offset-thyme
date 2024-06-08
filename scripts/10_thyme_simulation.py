# This script performs SLiM simulation of the thyme model using tree-sequence recording.

import subprocess
import sys
import pyslim, msprime
import numpy as np
import tempfile

## Read arguments from command line
if len(sys.argv) != 6:
    print("Usage: python <script> <slimScript> <Rscript> <sigmaFitness> <ticksAfter> <seed>", file=sys.stderr)
    sys.exit(1)
slimScript = sys.argv[1]
Rscript = sys.argv[2]
sigmaFitness = float(sys.argv[3])
ticksAfter = int(sys.argv[4])
seed = int(sys.argv[5])
rng = np.random.default_rng(seed=int(sys.argv[5]))

if __name__ == "__main__":
    ## Starting with diversity generated by coalescent
    breaks = [0, 1e8 + 1]
    recomb_map = msprime.RateMap(position=breaks, rate=[1e-8])
    demog_model = msprime.Demography()
    demog_model.add_population(initial_size=5000)
    ots = msprime.sim_ancestry(
        samples=1000,
        demography=demog_model,
        random_seed=rng.integers(1e9),
        recombination_rate=recomb_map,
    )
    ## Add SLiM annotations
    ots = pyslim.annotate(ots, model_type="nonWF", tick=1, stage="early")
    ## Add mutations
    mut_map = msprime.RateMap(position=breaks, rate=[5e-11])
    mut_model = msprime.SLiMMutationModel(type=3)
    ots = msprime.sim_mutations(
        ots, rate=mut_map, model=mut_model, keep=True, random_seed=rng.integers(1e9)
    )
    ## Add effect sizes to the mutations from N(0, 0.05)
    tables = ots.tables
    tables.mutations.clear()
    mut_map = {}
    for m in ots.mutations():
        md_list = m.metadata["mutation_list"]
        slim_ids = m.derived_state.split(",")
        assert len(slim_ids) == len(md_list)
        for sid, md in zip(slim_ids, md_list):
            if sid not in mut_map:
                mut_map[sid] = rng.normal(loc=0.0, scale=0.05)
            md["selection_coeff"] = mut_map[sid]
        tables.mutations.append(m.replace(metadata={"mutation_list": md_list}))
    ## check we didn't mess anything up
    assert tables.mutations.num_rows == ots.num_mutations
    ## Adjust metadata
    ts_metadata = tables.metadata
    ts_metadata["SLiM"]["spatial_dimensionality"] = "xyz"
    ts_metadata["SLiM"]["spatial_periodicity"] = "xy"
    tables.metadata = ts_metadata
    ots = tables.tree_sequence()
    ## Run SLiM
    with tempfile.TemporaryDirectory() as tmpdirname:
        init_file = tmpdirname + "/init.trees"
        ots.dump(init_file)
        # Redirect stdout to stdout and stderr to stderr
        try:
            subprocess.run(["slim", "-l", "0", "-s", str(seed), "-d", f"initFile='{init_file}'", "-d", f"RScript='{Rscript}'", "-d", f"sigmaFitness={sigmaFitness}", "-d", f"ticksAfter={ticksAfter}", slimScript])
        except subprocess.CalledProcessError as e:
            print(f"Error: SLiM simulation failed with return code {e.returncode}", file=sys.stderr)
            sys.exit(1)
