% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cumsum_mu_synch.R
\name{cumsum_mu_synch}
\alias{cumsum_mu_synch}
\title{Determination of Motor Unit Synchronization from Cross Correlation Histograms
using the Cumulative Sum Method}
\usage{
cumsum_mu_synch(motor_unit_1, motor_unit_2, order = 1, binwidth =
  0.001, get_data = T, plot = F)
}
\arguments{
\item{motor_unit_1, motor_unit_2}{Numeric vectors of strictly increasing
numbers denoting sequential discharge times of a motor unit or neuron or
any strictly increasing point process.}

\item{order}{Numeric as a positive integer for the number of forward and
backward orders for calculating recurrence times. Default = 1.}

\item{binwidth}{Numeric as a positive for the bin allocation size for
computational histogram. Default = 0.001 or 1 ms.}

\item{get_data}{T/F logical for outputting motor unit data. Default is TRUE.}

\item{plot}{T/F logical for outputting the cross correlation histogram.
Default is FALSE.}
}
\value{
A list of lists containing motor unit data (the names of each
  discharge train used, number of discharges, the interspike intervals (ISI),
  mean ISI, and the recurrence times associated with each order) and
  synchronization indices. #'   CIS = frequency of synchronized discharges.
  k' = ratio of total discharges in peak to expected discharges in peak. k'-1
  = ratio of synchronized discharges to expected discharges in peak. S =
  ratio of synchronized discharges to total number of discharges of both
  motor units. E = ratio of synchronized discharges to non-synchronized
  discharges. SI = ratio of synchronized discharges to reference motor unit
  discharges.
}
\description{
Calculates the time-domain synchronization indices CIS, k',
  k'-1, S, E, SI (detailed below) between the two input motor unit discharge
  trains based on the cumulative sum method. Peak boundaries are determined
  as the bins associated with 10% and 90% of the range (maximum minus
  minimum) of the cumulative sum. The peak is considered significant if its
  mean bin count exceeds the sum of the mean and 1.96 * standard deviation of
  the baseline bins (<= -60 ms and >= 60 ms). If no significant peak is
  detected, a default +/- 5 ms peak is used.
}
\examples{
   x <- c(0.035, 0.115, 0.183, 0.250, 0.306, 0.377, 0.455, 0.512, 0.577,
  0.656, 0.739, 0.821, 0.866, 0.950, 1.014, 1.085, 1.153, 1.213, 1.279,
  1.355, 1.431, 1.482, 1.551, 1.631, 1.692, 1.749, 1.832, 1.897, 1.964,
  2.106, 2.149, 2.229, 2.302, 2.384, 2.420, 2.505, 2.592, 2.644, 2.722,
  2.801, 2.870, 2.926, 3.011, 3.098, 2.030, 3.183, 3.252, 3.319, 3.395,
  3.469, 3.560, 3.589, 3.666, 3.744, 3.828, 3.876, 3.943, 4.020, 4.104)
  x <- sort(x)
  y <- sort(jitter(x))
  y <- round(y, digits = 3)
  cumsum_mu_synch(x, y, order = 1, binwidth = 0.001, get_data = TRUE,
  plot = FALSE)
}
\references{
Keen, D.A., Chou, L., Nordstrom, M.A., Fuglevand, A.J. (2012)
  Short-term Synchrony in Diverse Motor Nuclei Presumed to Receive Different
  Extents of Direct Cortical Input. Journal of Neurophysiology 108: 3264-3275
}
\keyword{cumulative}
\keyword{motor}
\keyword{recurrence,}
\keyword{sum}
\keyword{synchronization,}
\keyword{unit,}
