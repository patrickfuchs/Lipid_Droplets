# Lipid_Droplets

This repositery contains files allowing the simulation and analysis of lipid droplet systems as described in this paper [1] :

Bacle A., Gautier R., Jackson C.L., Fuchs P.F.J., Vanni S. (2017)
Interdigitation between core triglycerides and surface phospholipids modulates the surface properties of lipid droplets.
Biophysical Journal, 112, 1417-1430. [DOI](https://doi.org/10.1016/j.bpj.2017.02.032)

## Content

The different directories contain the following:

- `pure_TO`: a box of 108 trioleins (TO or tri-C18:1).
- `LD204`: a typical system mimimicking a lipid droplet consisting of an oily phase of 204 TO sandwiched between two POPC monolayers of 100 lipids surounded by water.
- `analysis`: a script `tabouret.py` allowing the analysis of TO conformations. This script is written in Python 2 and requires the modules [Numpy](https://numpy.org/), [Scipy](https://www.scipy.org/) and [mdtraj](https://mdtraj.org).

These two systems correspond to the entries `pureTO` and `LD0` in table 1 of ref. [1].

## Force field

The POPC parameters (`popc.itp`) come from the Berger force field [2] as implemented in GROMACS by Peter Tieleman (see [here](http://wcm.ucalgary.ca/tieleman/downloads)). The TO parameters (`triolein.itp`) were adapted from Berger by Vattulainen and co-workers [3]. Typically, starting from a POPC molecule, the sn-3 was replaced by an oleoyl chain. For all oleoyl chains both in POPC and TO, a correction on the double bond [4] was applied. All the systems were simulated with GROMACS 4 [6]. 

If you use these files, please cite Bacle et al. (2017) [1] and Hall et al. (2008) [3].

We kindly thank Ilpo Vattulainen for sharing the triolein parameters.

Note: Technically, the files in this repo make use of the so called *ffgmx* force field (based on a modified GROMOS 87), which was the standard way to simulate Berger lipids as implemented by [Peter Tieleman]((http://wcm.ucalgary.ca/tieleman/downloads)). However, *ffgmx* has been deprecated in recent GROMACS versions (5.0.* and above). The files (`ffgmx.itp`, etc.) can be found on old versions 4.5.* and below within any GROMACS archive on the [GROMACS download page](https://manual.gromacs.org/). If you use lipids only (POPC and TO), as in our work [1], there is no problem, this is the classical method for simulating Berger lipids. However, if you plan to use the files in this repo with proteins, you might want to use a more recent GROMOS force field as described in [Justin Lemkul's tutorial](http://www.mdtutorials.com/gmx/membrane_protein/02_topology.html).

## References

[1] Bacle A., Gautier R., Jackson C.L., Fuchs P.F.J., Vanni S. (2017)
Interdigitation between core triglycerides and surface phospholipids modulates the surface properties of lipid droplets.
Biophysical Journal, 112, 1417-1430. [DOI](https://doi.org/10.1016/j.bpj.2017.02.032)

[2] Berger O., O. Edholm, and F. Jahnig. 1997. 
Molecular dynamics simulations of a fluid bilayer of dipalmitoylphosphatidylcholine at full hydration, constant pressure, and constant temperature.
Biophys. J., 72, 2002-2013.

[3] Hall A., J. Repakova, and I. Vattulainen. 2008. 
Modeling of the triglyceride-rich core in lipoprotein particles.
J. Phys. Chem. B., 112, 13772-13782.

[4] Bachar M., P. Brunelle, D.P. Tieleman, and A. Rauk. 2004. 
Molecular dynamics simulation of a polyunsaturated lipid bilayer susceptible to lipid peroxidation.
J. Phys. Chem. B., 108, 7170-7179.

[5] Hess, B., C. Kutzner, ..., E. Lindahl. 2008. 
GROMACS 4: algorithms for highly efficient, load-balanced, and scalable molecular simulation.
J. Chem. Theory Comput., 4, 435-447.

## Licence

All the files in this repository are under licence Creative Commons Attribution - Partage dans les Mêmes Conditions 3.0 France (CC BY-SA 3.0 FR).
