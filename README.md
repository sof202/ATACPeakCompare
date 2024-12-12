# Peak Compare

When using MACS3 to peak call data, you may want to know if this peak is shared
between each dataset. Due to the issues around short read data and the nature
of how they manifest in count based data this task is made difficult. A naive
approach would be to just use a set operation on the output bed files of your
two samples (intersection, for example). However, this approach has proven to
be ineffective as peaks don't align nicely between datasets due to the
dichotomous nature of peak calling (a base can be in a peak, or not in a peak).
The purpose of this repository is to provide a metric that can be used to
suggest if a peak is in both datasets (or not).

This repository contains a pipeline that peak calls two samples in the steps
outlined in the 
[official macs3 tutorial](https://macs3-project.github.io/MACS/docs/Advanced_Step-by-step_Peak_Calling.html).

## Setup

In order to run the scripts included in this repository, you will need
to complete the below required steps:

1) Download all [required software](#software-requirements)
2) Copy the configuration templates found in `Setup` and place these files in
memorable locations (recommended: next to the data you are processing). For the
most part, parameters are explained in these templates. However additional
information is given in the documentation. You will likely want configuration
files for each of your datasets so that you can remember what parameters were
used.
3) Run the `PeakCall.sh` script for each of your datasets
4) Run the `PeakCompare.sh` script.
5) Interpret the metric.

For additional information it is recommended that you read the
[documentation](https://sof202.github.io/PeakCompare/).

## Software Requirements

Please ensure that each of these pieces of software can be found on your
`PATH` environment variable.

- [bash](https://www.gnu.org/software/bash/) (>=4.2.46(2))
- [SLURM Workload Manager](https://slurm.schedmd.com/overview.html) (>=20.02.3)
- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) (>=v23.10.0)
