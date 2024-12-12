# ATAC Peak Compare

When using MACS3 to peak call ATAC-Seq, you may want to know if this peak is
shared between each dataset. Due to the issues around short read data and the
nature of how they manifest in ATAC-Seq data this task is made difficult. A
naive approach would be to just use a set operation on the output bed files
of your two samples (intersection, for example). However, this approach has
proven to be ineffective as peaks don't align nicely between datasets due to
the dichotomous nature of peak calling (a base can be in a peak, or not in a
peak). The purpose of this repository is to provide a metric that can be used
to suggest if a peak is in both datasets (or not).

This repostory contains a pipeline that peak calls two samples in the steps
outlined in the 
[official macs3 tutorial](https://macs3-project.github.io/MACS/docs/Advanced_Step-by-step_Peak_Calling.html).

Although this pipeline was designed to look at ATAC data, the only assumption
on your data is that it is peak callable (count based data). Hence ChIP-Seq
and DNase-Seq (and many more) are also applicable to be used here.

## Setup

In order to run the scripts included in this repository, you will need
to complete the below required steps:

1) Download all [required software](#software-requirements)
2) Copy the config templates found in `Setup` and place these files in
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

## Overview of pipeline

After the bias track (lambdas) have been generated, the following steps are
completed:

- Peaks are called in both datasets once without merging close peaks and one
with merging close peaks (tiny peaks are removed in both peak called datasets)
- Bases that only exist in peaks due to the merging process are identified
- Confidence intervals are calculated for the bias track and used to get 
confidence intervals on the p-values (used for the peak calling process).
- For each base, a comparison is made between the confidence intervals of each
dataset.
- For each base a Boolean value is generated where:
  - FALSE if the base's CI is lower in the comparison dataset than in the
  reference dataset.
  - TRUE if the base's CI is higher in the comparison dataset than in the
  reference dataset 
  - TRUE if the base is part of a peak in the comparison dataset and the base
  is only a peak in the reference dataset due to the merging process
- A region is then defined by the user, likely a region that is calssified as a
peak in the reference dataset. For this region, the percentage of bases that
were marked TRUE in the previous step (so called 'psuedopeaks') is calculated.
This value is the promised metric.

### How to interpret the metric

A value close to one implies that this region in the comparison dataset is a
comparable peak as to what is seen in the reference dataset. This is stronger
than "the peak is called in both data sets in this region". This is because the
metric takes into account the p-values that are used to call the peaks by MACS.
This metric will only be close to 1 if both datasets have similar levels of
significance (or the comparison dataset has even more significant evidence for
a peak). In the event that the comparison set has more significant evidence for
having peaks than the reference set, it is advised to switch the roles of the
two datasets (and rerun analysis).

A value close to zero implies that this region in the comparison dataset is not
a comparable peak as to what is seen in the reference dataset. If you were to
use MACS only, you may well find that the region may have peaks called in the
comparison dataset. However, here we require that the significance of said peak
is comparable in each dataset. In this sense the metric is more strict than
simply looking at overlaps between the datasets.

The metric in reality can vary between zero and one and it is down to the user
(and their downstream analyses) of what they want to do with this information.
