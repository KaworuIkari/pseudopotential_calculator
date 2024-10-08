# Bulk Convergence

The first step in generating a pseudopotential is to find the spacing
of atoms in the bulk. To do this we generate a single primitive cell
of material with the correct bond angles and a reasonable initial guess
for the bond lengths

![](figures/initial_arrangement.png)

We then perform a CASTEP geometry optimization routine in order to
find the true bond lengths of the material.

## Cutoff Energy

First we test for convergence against cutoff-energy. This parameter sets
the highest energy in the basis used for the calculation.

![](figures/energy_against_cutoff_energy.png)
![](figures/cell_length_against_cutoff_energy.png)

We find convergence at around $600 eV$, although with the default
pseudopotential and `k_points` of $10 \times 10 \times 10$ we are
still a ways off the theoretical bond length of $2.53 \AA$.
Note that this cutoff energy roughly
corresponds to the `EXTREME` setting in the `.usp` file generated by
CASTEP.

```
START COMMENT
 272 COARSE
 327 MEDIUM
 354 FINE
 566 EXTREME
```

this file can be generated by calling

```bash
castep.serial --dryrun <task_name>
```

and it is usually safe to use the `EXTREME` setting in place of full convergence testing.

## K-Points Convergence

Next we test convergence against the number of k-point in the simulation. CASTEP works with
periodic boundary conditions, and makes use of bloch's theorem to reduce the cost associated with
a calculation. The number of k-points therefore sets the effective 'size' of the simulation.

For a bulk calculation we test convergence with `k_points` equal to $n \times n \times n$
for some $n$. We repeat this calculation for a range of pseudopotentials(PBE and WC), and also test the
use of the `spin_polarized` variable

![](figures/energy_against_n_k_points.Png)
![](figures/cell_length_against_n_k_points.Png)

Notice the WC pseudopotential gives a better converged cell length compared to that given by PBE.
The `spin_polarized` variable only has a small effect at low cut off energy.

## Final Bulk Configuration

We found that the calculation converged with an energy cutoff of 600ev, using WC as the pseudopotential with spin_polarised state off.

The final configuration as shown below, can then be used to generate a slab of Copper.

![](figures/final_arrangement.png)
