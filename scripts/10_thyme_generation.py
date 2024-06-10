import tskit, pyslim, msprime
import numpy as np
import sys

infiles = sys.argv[1:]
outfile = infiles.pop()
age_mother_at_birth = list()
for infile in infiles:
    ts = tskit.load(infile)
    relationships = pyslim.individual_parents(ts)
    mother_birth = ts.individuals_time[relationships[:,0]]
    child_birth = ts.individuals_time[relationships[:,1]]
    age_at_birth = mother_birth-child_birth
    print(f"Median {np.median(age_at_birth)}")
    print(f"Mean {np.mean(age_at_birth)}")
    age_mother_at_birth.extend(age_at_birth)
age_mother_at_birth = np.array(age_mother_at_birth)
print(f"Final median {np.median(age_mother_at_birth)}")
print(f"Final mean {np.mean(age_mother_at_birth)}")
np.save(outfile, age_mother_at_birth)

