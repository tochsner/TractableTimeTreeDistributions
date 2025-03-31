# TractableTimeTreeDistributions

This repository contains the code for the julia package TractableTimeTreeDistributions. It allows to approximate the distribution of a tree given a set of MCMC posterior samples.

## 🚀 Getting Started

The easiest way to install the package is by downloading the latest binaries from the [release page](https://github.com/tochsner/TractableTimeTreeDistributions/releases). After downloading the archive, unpack it and you're ready to go:

```bash
unzip TractableTimeTreeDistributions-vX.X.X.zip
./TractableTimeTreeDistributions-vX.X.X/bin/TractableTimeTreeDistributions <arguments>
```

Alternatively, you can [install Julia](https://julialang.org/downloads/) and install the package as a Julia package:

```julia
julia> using Pkg
julia> Pkg.add("TractableTimeTreeDistributions")
julia> Pkg.add("Distributions")
```

Another option is to clone the repository and use the package directly from the source code:

```bash
git clone https://github.com/tochsner/TractableTimeTreeDistributions.git
cd TractableTimeTreeDistributions

julia
julia>]
pkg> activate .
pkg> instantiate
```

## 📖 Documentation

### Interactive Workflow with Model Selection (recommended)

The recommended workflow is the following:

1. Run the following command to fit the recommended distributions and perform model selection. Input is a file containing the MCMC posterior trees in Newick or Nexus format.

```bash
./TractableTimeTreeDistributions-vX.X.X/bin/TractableTimeTreeDistributions my_mcmc_posterior_trees.trees 
```

2. After initial processing, you will get an output like this:

```bash
--------------------
Log Data Likelihoods
--------------------
[ Info: Height (LogNormal), Ratios (Beta):                      245161.39951377752
[ Info: Height (LogNormal), Ratios (LogitNormal):               242781.9853789675
[ Info: Shortest Branch (Gamma):                                242466.63796884567
[ Info: Shortest Branch (Weibull):                              241664.00000117905
[ Info: Shortest Branch (LogNormal):                            238382.73921265977
[ Info: Last Divergence (Gamma), Branches (Gamma):              235798.51548722142
[ Info: Last Divergence (Weibull), Branches (Weibull):          235164.51171292344
[ Info: Last Divergence (LogNormal), Branches (LogNormal):      231618.75107464782
```

3. You will now be asked if you want to use the model with the highest likelihood. Alternatively, you can also specify the model you want to use:

```bash
---------------
Model Selection
---------------
We recommend to use the Height (LogNormal), Ratios (Beta) distribution.
Do you want to use it? [y/n]
y
```

4. You can now perform multiple tasks:

```bash
--------------------
Tractable Operations
--------------------
Select one or more tasks to perform using the selected distribution:
 > Create a credible-region plot
   Create a log-likelihood plot
   Create a point estimate
   Get the likelihood of a given tree
   Test if a given tree is in the smallest 95%-credible region
```

### Programmatic Usage of the Julia Package

You can use the package in your own Julia scripts. See the follwing basic example and the scripts in the `experiments` folder for more advanced usage.

```julia
using TractableTimeTreeDistributions
using Distributions

# Load MCMC trees
ref_trees = load_trees("/path/to/reference/trees")

# Fit a distribution
distribution = TractableTimeTreeDist{CCD1,HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}}(ref_trees)

# Get the log-likelihood of a tree
query_tree = load_trees("/path/to/query/tree")[begin]
@info log_density(distribution, query_tree)

# Sample a tree and write it to a file
sampled_tree = sample_tree(distribution)
write_tree("sample.trees", sampled_tree)

# Get a point estimate and write it to a file
point_estimate_tree = point_estimate(distribution)
write_tree("point_estimate.trees", point_estimate_tree)
```

##  🔗 Citation

If you use this package in your research, please cite the following paper:

```
@article{ochsner2022tractable,
  title={Tractable Time Tree Distributions},
  author={Ochsner, Tobia and Klawitter, Jonathan and Drummond, Alexei J.},
  journal={arXiv preprint arXiv:2201.12219},
  year={2022}
}
```