# m6a-model
Gillespie model for m6a spreading by synthetic read-write modules

**Publication:**
Minhee Park, Nikit Patel, Albert J. Keung and Ahmad S. Khalil. "Engineering Epigenetic Regulation Using Synthetic Read-Write Modules." Cell (2019)


**Main functions:**
1) fun_clus_noRW.m (Runs 1000 simulations on clustered reporter without Reader-writer protein)
2) fun_clus_withRW.m (Runs 1000 simulations on clustered reporter with Reader-writer protein)

**Reaction Propensities Functions:**

3) lattice_clus_noRW.m (Rxn propensities for clustered reporter without Reader-Writer protein)

4) lattice_clus_withRW.m (Rxn propensities for clustered reporter with Reader-Writer protein)

**Gillespie Function:**

5) directMethod.m

**Plotting Function:**

6) boundedline.m

**Data:**

7) Dam_properties.mat - measured data from Extended Figure 3a-c

**NOTE:** ALL SIMULATIONS OF PAPER RUN ON MATLAB 2013a

**Author:** Nikit Patel
