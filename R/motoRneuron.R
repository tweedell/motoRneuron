#'motoRneuron: A package for computing motor unit analyses
#'
#'The motoRneuron package provides functions for time-domain analyses of motor
#'units.
#'
#'@section About: The motoRneuron package was developed at the United States
#'  Army Research Laboratory.
#'
#'@section Details: The motoRneuron package provides functions to compute
#'  time-domain synchronization between the two motor unit discharge trains.
#'  MotoRneuron's primary functions configure and assess the cross correlation
#'  histogram between motor unit discharge trains for a "peak". Peaks indicate
#'  more synchronized firing at that latency than would be expected due to
#'  chance. The magnitude of the peak is also thought to indicate the strength
#'  of common input between the motor units and is characterized by
#'  synchronization "indices", listed below. Three methods are available
#'  for assessing the peak in the histogram. Details of specific methods are
#'  documented in their respective help files. Motor unit characteristic data
#'  and interspike intervals are automatically output.
#'
#'  'Recurrence_intervals', 'bin', and 'plot_bins' are support functions used
#'  within the core functions, but can also be called separately for individual
#'  use. Recurrence_intervals() calculates recurrence intervals between two
#'  motor unit discharge trains. Bin() discretizes the recurrence intervals into
#'  user-specified bins. Plot_bins() leverages ggplot2 to plot the resulting
#'  histogram.
#'
#'@section Methods for Peak Determination: Visual method - subjective analysis
#'  of normalized cumulative sum graph to determine peak location. (Nordstrom,
#'  Fuglevand, and Enoka, 1992)
#'
#'  Z-score method - a random uniform distribution is used to calculate a mean +
#'  1.96 standard deviation threshold for the experimental cross correlation
#'  histogram. Any bins of the experimental histogram within +/- 10 ms of 0 that
#'  crossed the threshold are considered to be significantly greater than
#'  expected due to chance and subsequently used for analysis. (Defreitas, Beck,
#'  Xin, and Stock, 2013)
#'
#'  Cumulative Sum (cumsum) method - the cumulative sum of the bin count of the
#'  histogram is calculated. The peak boundaries are calculated as 10% and 90% of
#'  the range (maximum minus minimum) of the cumulative sum. The peak is
#'  considered significant if its mean bin count exceeds the sum of the mean and
#'  1.96 * standard deviation of the baseline bins (<= -60 ms and >= 60 ms). If
#'  no significant peak is detected, a default +/- 5 ms peak is used. (Keen,
#'  D.A., Chou, L., Nordstrom, M.A., Fuglevand, A.J., 2012)
#'
#'@section Synchronization Indices: CIS = frequency of synchronized discharges.
#'  (Nordstrom, Fuglevand, and Enoka, 1992)
#'
#'  k' = ratio of synchronized discharges to expected discharges in peak. (Kamen
#'  & Roy, 2000)
#'
#'  k'-1 = ratio of extra synchronized discharges to expected discharges in
#'  peak. (Kamen & Roy, 2000)
#'
#'  S = ratio of synchronized discharges to total number of discharges of both
#'  motor units. (Nordstrom, Fuglevand, and Enoka, 1992)
#'
#'  E = ratio of synchronized discharges to non-synchronized discharges. (Datta
#'  & Stephens 1990)
#'
#'  SI = ratio of synchronized discharges to reference motor unit discharges.
#'  (Defreitas, Beck, Xin, and Stock, 2013)
#'
#'@section Authors: Andrew Tweedell and Matthew Tenan, maintainer: Andrew Tweedell
#'  (andrew.j.tweedell.civ@@mail.mil)
#'
#'@section References: Datta, A.K., & Stephens, J.A. (1990) Synchronization of
#'  Motor Unit Activity During Voluntary Contraction in Man. Journal of
#'  Physiology 422: 397-419
#'
#'  DeFreitas, J.M., Beck, T.W., Xin, Y., Stock, M.S. (2013) Synchronization of
#'  Low- and High-Threshold Motor Units. Muscle & Nerve DOI 10.1002/mus.23978
#'
#'  Kamen, G. & Roy, A. (2000) Motor Unit Synchronization in Young and Elderly
#'  Adults. European Journal of Applied Physiology. 81: 403-410
#'
#'  Keen, D.A., Chou, L., Nordstrom, M.A., Fuglevand, A.J. (2012) Short-term
#'  Synchrony in Diverse Motor Nuclei Presumed to Receive Different Extents of
#'  Direct Cortical Input. Journal of Neurophysiology 108: 3264-3275
#'
#'  Nordstrom, M.A., Fuglevand, A.J., Enoka, R.M. (1992) Estimating the Strength
#'  of Common Input to Human Motoneurons from the Cross-Correlogram. Journal of
#'  Physiology 453, pp. 547-574
#'
#'@docType package
#'@name motoRneuron
NULL