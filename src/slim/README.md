# SLim simulations

This directory contains a briefly description of every experiment and SLim simulation. 

## Impact of missidentified causal *loci* in genomic offset

- "m1.1.slim" This simulates 25 demes (in a 5x5 grid), with 100 unlinked QTL with additive effects to a trait which is directly related to an environmental gradient that goes through the x-axis. 

- "m2.1.slim" This simulation was adapted from [Gain et al](https://github.com/bcm-uga/geneticgap/blob/master/slimwork/poly_small_17.slim). More details in their supplementary materials. Briefly, it is a non-Wright-Fisher simulation with two traits which are directly related to an environmental gradient that goes through the x-axis and y-axis. There are two phases, one demographic and a second one were the local adaptation occur. Mutations do not appear at random, but we set them. There are 120 QTL uniformly distributed across the genome. 
- "m2.2.slim" The exact same setting as "m2.1.slim", but the number of QTLs goes from 10 to 100. 

## Quality control

This simulations are adapted from m2.1.slim. The idea is to contrast GO when computed on locally and non locally adapted loci. 

1. "m3.1.slim" No local adaptation.
2. "m3.2.slim" Local adaptation.
3. "m3.3.slim" Local adaptation, with asymmetric variance fitness. 
3. "m3.4.slim" Local adaptation, with asymmetric variance fitness (stronger)