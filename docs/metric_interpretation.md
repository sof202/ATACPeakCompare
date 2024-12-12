# Metric interpretation

For a formal definition of the metric, head over to 
[this page](./peak_comparing.md#how-the-metric-is-calculated).

A value close to one implies that the region in the comparison dataset is a
comparable peak as to what is seen in the reference dataset. This is stronger
than "the peak is called in both data sets in this region". This is because the
metric takes into account the p-values that are used to call the peaks by MACS.
This metric will only be close to 1 if both datasets have similar levels of
significance (or the comparison dataset has even more significant evidence for
a peak).

A value close to zero implies that the region in the comparison dataset is
*not* a comparable peak as to what is seen in the reference dataset. If you
were to use MACS only, you may well find that the region may have peaks called
in the comparison dataset. However, here we require that the significance of
said peak is comparable in each dataset. In this sense the metric is more
strict than simply looking at overlaps between the datasets.

The metric in reality can vary between zero and one and it is down to the user
(and their downstream analyses) of what they want to do with this information.
You may choose to simply use the metric in a binary sense, where a certain
threshold means the peak is in both datasets. Alternatively you might want
to apply this methodology to get the metric for a wide range of regions or a
wide set of samples. Once obtained, correlation analysis could be used to
see how variable peaks are in certain regions of the genome for certain
samples/cell types.

## Note

In the event that the comparison set has more significant evidence for
having peaks than the reference set, it is advised to switch the roles of the
two datasets (and rerun analysis).

You may also encounter scenarios where the metric exceeds 1. This is likely to
be the case if the region selected does not house a peak in the reference
dataset. In such scenarios it is recommended to either switch the roles of the
datasets or to inspect a different region of the genome.

