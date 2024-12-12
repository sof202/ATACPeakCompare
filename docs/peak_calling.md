# Peak Calling

Before running the pipeline that is used to calculate the metric, a few
different files are required that can be obtained via MACS. Usually, peak
calling with MACS is via the function `callpeaks`. However, this process does
not give us the outputs of intermediate steps that MACS carries out. As such
the script `PeakCall.sh` will use the fundamental functions of MACS to generate
these intemediate files:

- The coverage/pileup track
- The bias track
- The pvalues track
- The set of peaks before merge process
- The set of peaks after merge process

If you want to understand the steps behind the MACS peak calling process, it is
recommended that you read the 
[official tutorial](https://macs3-project.github.io/MACS/docs/Advanced_Step-by-step_Peak_Calling.html).

## Running

To run this script, you will need to complete the 
[setup](https://github.com/sof202/PeakCompare?tab=readme-ov-file#setup) first.
After this, you can submit the script `PeakCall.sh` to a SLURM queue with:

```bash
sbatch -p queue .../PeakCall.sh path/to/config_file.txt
```

The config file to use for this script is given in `example_config_call.txt`.
The queue to send the job to will depend on the HPC system you are using. You
can view all available partitions with:

```bash
sinfo show partitions | awk '{print $1}' | uniq
```

The script is designed to work with or without a config file and will
process a single set of reads in any format MACS can work with. It is expected
that you run this script for each sample (or merged set of samples) so it is
recommended that you have two config files (so jobs can run concurrently).

## Why merged and unmerged peaks

In the final step of the peak calling process, MACS encourages the user to 
merge any peaks that are suitably close together. The recommended minimum gap
between peaks is the read length of the dataset. The metric needs to be able to
differentiate between peaks that are called due to this merging process because
such base pair positions will have very low coverage. The metric heavily relies
on comparing pvalues and read counts, if a peak has very low pvalues and read
counts, this can be misleading if the merging process is not accounted for.

## Cutoff analysis

In order to call peaks, a line in the sand must be drawn somewhere. This often
feels arbitrary and is a difficult topic in statistics. However, MACS does
allow for cutoff analysis where different thresholds are used and some metrics
are displayed for each threshold used. These metrics are:

- The number of peaks
- The average length of the peaks
- The total length of all peaks combined

Using these there are a couple of different approaches you can take to get a
suitable value for the cutoff. 

- Elbow analysis
- Using biological knowledge
  - The expected number of peaks
  - The expected average peak length

Currently the pipeline allows you to either pick a hard cutoff value or use the
most stringent cutoff that gives you a certain average peak length. If you have
your own idea here, a file is put in your chosen output directory under
`(sample_name)_cutoff_analysis.txt`. The pipeline doesn't use elbow analysis
as this is hard to automate ("How do you define a dramatic change?").

## Time

The time MACS takes to peak call will vary with the number of reads that you
have, however the process is usually very fast. For a dataset with 11,000,000
reads the total time taken was roughly 5 minutes 
(where the cpu was Intel(R) Xeon(R) CPU E5-2640 v3 @ 2.60GHz).
