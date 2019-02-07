#' Determination of Motor Unit Synchronization from Cross Correlation Histograms
#' using the Cumulative Sum Method
#'
#' @export
#' @importFrom stats "median" "sd"
#' @importFrom methods "show"
#' @keywords recurrence, motor unit, synchronization, cumulative sum
#' @description Calculates the time-domain synchronization indices CIS, k',
#'   k'-1, S, E, SI (detailed below) between the two input motor unit discharge
#'   trains based on the cumulative sum method. Peak boundaries are determined
#'   as the bins associated with 10% and 90% of the range (maximum minus
#'   minimum) of the cumulative sum. The peak is considered significant if its
#'   mean bin count exceeds the sum of the mean and 1.96 * standard deviation of
#'   the baseline bins (<= -60 ms and >= 60 ms). If no significant peak is
#'   detected, a default +/- 5 ms peak is used.
#' @usage cumsum_mu_synch(motor_unit_1, motor_unit_2, order = 1, binwidth =
#'   0.001, get_data = T, plot = F)
#' @param motor_unit_1,motor_unit_2 Numeric vectors of strictly increasing
#'   numbers denoting sequential discharge times of a motor unit or neuron or
#'   any strictly increasing point process.
#' @param order Numeric as a positive integer for the number of forward and
#'   backward orders for calculating recurrence times. Default = 1.
#' @param binwidth Numeric as a positive for the bin allocation size for
#'   computational histogram. Default = 0.001 or 1 ms.
#' @param get_data T/F logical for outputting motor unit data. Default is TRUE.
#' @param plot T/F logical for outputting the cross correlation histogram.
#'   Default is FALSE.
#' @return A list of lists containing motor unit data (the names of each
#'   discharge train used, number of discharges, the interspike intervals (ISI),
#'   mean ISI, and the recurrence times associated with each order) and
#'   synchronization indices. #'   CIS = frequency of synchronized discharges.
#'   k' = ratio of total discharges in peak to expected discharges in peak. k'-1
#'   = ratio of synchronized discharges to expected discharges in peak. S =
#'   ratio of synchronized discharges to total number of discharges of both
#'   motor units. E = ratio of synchronized discharges to non-synchronized
#'   discharges. SI = ratio of synchronized discharges to reference motor unit
#'   discharges.
#' @examples
#'    x <- c(0.035, 0.115, 0.183, 0.250, 0.306, 0.377, 0.455, 0.512, 0.577,
#'    0.656, 0.739, 0.821, 0.866, 0.950, 1.014, 1.085, 1.153, 1.213, 1.279,
#'    1.355, 1.431, 1.482, 1.551, 1.631, 1.692, 1.749, 1.832, 1.897, 1.964,
#'    2.106, 2.149, 2.229, 2.302, 2.384, 2.420, 2.505, 2.592, 2.644, 2.722,
#'    2.801, 2.870, 2.926, 3.011, 3.098, 2.030, 3.183, 3.252, 3.319, 3.395,
#'    3.469, 3.560, 3.589, 3.666, 3.744, 3.828, 3.876, 3.943, 4.020, 4.104)
#'    y <- jitter(x)
#'    cumsum_mu_synch(x, y, order = 1, binwidth = 0.001, get_data = T, plot = T)
#' @references Keen, D.A., Chou, L., Nordstrom, M.A., Fuglevand, A.J. (2012)
#'   Short-term Synchrony in Diverse Motor Nuclei Presumed to Receive Different
#'   Extents of Direct Cortical Input. Journal of Neurophysiology 108: 3264-3275

cumsum_mu_synch <- function(motor_unit_1, motor_unit_2, order = 1,
                            binwidth = 0.001, get_data = T, plot = F) {

  recurrence_intervals2 <- function(motor_unit_1, motor_unit_2, order) {


    if (!is.vector(motor_unit_1) || !is.vector(motor_unit_2)) {
      stop("'motor_unit_1' and 'motor_unit_2' must be vectors.")
    }

    if (length(motor_unit_1) <= 1 || length(motor_unit_2) <= 1) {
      stop ("'motor_unit_1' and 'motor_unit_2' must be vectors of length > 1.")
    }

    if (is.unsorted(motor_unit_1, strictly = T)
        || is.unsorted(motor_unit_2, strictly = T)) {
      stop ("'motor_unit_1' and 'motor_unit_2' must be strictly increasing.")
    }

    if (!is.numeric(order) || order%%1 != 0) {
      stop("Order must be whole number.")
    }

    # reference (ref) and event motor units (MU) assigned according to which MU
    # has more firings (length()). ISI = InterSpike Intervals
    if (length(motor_unit_1) < length(motor_unit_2)) {
      ref.name <- deparse(substitute(motor_unit_1, env = parent.frame()))
      event.name <- deparse(substitute(motor_unit_2, env = parent.frame()))
      ref.MU <- motor_unit_1
      event.MU <- motor_unit_2
      ref.MU.ISI <- diff(motor_unit_1)
      event.MU.ISI <- diff(motor_unit_2)
      mean.ref.ISI <- round(mean(ref.MU.ISI), digits = 3)
      mean.event.ISI <- round(mean(event.MU.ISI), digits = 3)
    } else {
      ref.name <- deparse(substitute(motor_unit_2, env = parent.frame()))
      event.name <- deparse(substitute(motor_unit_1, env = parent.frame()))
      ref.MU <- motor_unit_2
      event.MU <- motor_unit_1
      ref.MU.ISI <- diff(motor_unit_2)
      event.MU.ISI <- diff(motor_unit_1)
      mean.ref.ISI <- round(mean(ref.MU.ISI), digits = 3)
      mean.event.ISI <- round(mean(event.MU.ISI), digits = 3)
    }

    MU.names <- list(Reference_Unit = ref.name,
                     Number_of_Reference_Discharges = length(ref.MU),
                     Reference_ISI = ref.MU.ISI,
                     Mean_Reference_ISI = mean.ref.ISI,
                     Event_Unit = event.name,
                     Number_of_Event_Discharges = length(event.MU),
                     Event_ISI = event.MU.ISI,
                     Mean_Event_ISI = mean.event.ISI,
                     Duration = max(ref.MU, event.MU) - min(ref.MU, event.MU))

    lags <- vector('list', order)

    for (i in 1:length(ref.MU)) {

      pre_diff <- rev(event.MU[event.MU < ref.MU[i]])
      pre_diff <- pre_diff[1:order]
      pre_diff <- pre_diff - (ref.MU[i])

      post_diff <- event.MU[event.MU >= ref.MU[i]]
      post_diff <- post_diff[1:order]
      post_diff <- post_diff - (ref.MU[i])

      for (j in 1:order) {
        y <- c(pre_diff[j], post_diff[j])
        lags[[j]] <- append(lags[[j]], y)
      }

    }

    # remove NA's
    lags <- lapply(lags, Filter, f = Negate(is.na))

    names(lags) <- paste(1:order)
    lags <- append(MU.names, lags)

    return(lags)

  }

  recurrence.data <- recurrence_intervals2(motor_unit_1,
                                           motor_unit_2,
                                           order)

  mean.reference.ISI <- recurrence.data$Mean_Reference_ISI

  # Create frequency table by binning recurrence times according to specfied bin
  # width using the mean reference ISI as the positive and negative boundaries.
  frequency.data <- unlist(recurrence.data[paste(1:order)])
  frequency.data <- frequency.data[frequency.data >= -mean.reference.ISI &
                                   frequency.data <= mean.reference.ISI]
  frequency.data <- as.vector(frequency.data)
  frequency.data <- motoRneuron::bin(frequency.data, binwidth = binwidth)

  # Calculate mean frequency of baseline bins (all bins outside +/- 60 ms).
  # Keen et al 2012
  baseline.mean <- frequency.data[frequency.data$Bin <= ((min(frequency.data$Bin))
                                  + 0.060) | frequency.data$Bin >=
                                  (max(frequency.data$Bin) - 0.060), ]
  baseline.mean <- mean(as.numeric(unlist(baseline.mean["Freq"])))

  # Calculate sd frequency of baseline bins (all bins outside +/- 60 ms).
  baseline.sd <- frequency.data[frequency.data$Bin <= ((min(frequency.data$Bin))
                                                + 0.060) | frequency.data$Bin >=
                                    (max(frequency.data$Bin) - 0.060), ]
  baseline.sd <- sd(as.numeric(unlist(baseline.sd["Freq"])))

  # calculate cumulative sum as new data frame
  cumsum <- data.frame(Bin = frequency.data$Bin, Cumsum = cumsum(
                       as.numeric(frequency.data$Freq) - baseline.mean))
  cumsum <- cumsum[cumsum$Bin >= ((min(frequency.data$Bin)) + 0.060) &
                     cumsum$Bin <= (max(frequency.data$Bin) - 0.060),]

  # calculate 90th and 10th percentile of max difference in cumsum
  ninety.percent <- min(cumsum$Cumsum) + ((max(cumsum$Cumsum) - min(cumsum$Cumsum)) * 0.9)

  ten.percent <- min(cumsum$Cumsum) + ((max(cumsum$Cumsum) - min(cumsum$Cumsum)) * 0.1)

  # find upper and lower bounds (bins) associated with the 90th and 10th
  # percentile
  bounds <- vector()
  bounds[1] <- cumsum[(which(abs(cumsum$Cumsum - ten.percent)
                             == min(abs(cumsum$Cumsum - ten.percent)))), 1]

  bounds[2] <- cumsum[(which(abs(cumsum$Cumsum - ninety.percent)
                             == min(abs(cumsum$Cumsum - ninety.percent)))), 1]

  # swap upper and lower limits if needed
  if(bounds[1] > bounds[2]) {

    old.lower <- bounds[1]
    bounds[1] <- bounds[2]
    bounds[2] <- old.lower
    rm(old.lower)

  }

  # subset out peak from upper and lower bounds
  peak <- frequency.data[frequency.data$Bin >= bounds[1] &
                        frequency.data$Bin <= bounds[2],]
  peak <- as.numeric(unlist(peak["Freq"]))

  peak.mean <- mean(peak)
  peak.zscore <- (peak.mean - baseline.mean) / baseline.sd

  # if peak.zscore >= 1.96 default to +/- 5 ms and define new peak
  if(peak.zscore < 1.96) {

    bounds[1] <- -0.005
    bounds[2] <- 0.005

    peak <- frequency.data[frequency.data$Bin >= bounds[1] &
                             frequency.data$Bin <= bounds[2],]
    peak <- as.numeric(unlist(peak["Freq"]))

    peak.mean <- mean(peak)
    peak.zscore <- (peak.mean - baseline.mean) / baseline.sd

  }

  # Calculate total number of instances in peak
  total.peak <- sum(peak)

  # Calculate number of instances in peak in excess of what is
  # expected (i.e. above baseline mean)
  extra.peak <- sum((peak - baseline.mean)[which((peak - baseline.mean) > 0)])

  # Calculate total number of instances in frequency data
  total.count <- sum(as.numeric(unlist(frequency.data$Freq)))

  # Determine bins in peak below baseline mean
  q <- as.numeric(vector())
  for (m in 1:length(peak)) {
    if (peak[m] <= baseline.mean) {
      q <- c(q, peak[m])
    } else {next}
  }

  # Calculate number of instances in peak below baseline mean
  expected.peak <- baseline.mean *
                   (length(which(peak > baseline.mean))) +
                   sum(q)

  Cumsum.Synch <- list()

  if (get_data) {

    Cumsum.Synch[["Data"]] <- recurrence.data

  }

  if (plot) {

    show(plot_bins(frequency.data))

  }

  Cumsum.Synch[["Indices"]] <- list(CIS = extra.peak / recurrence.data$Duration,
                          kprime = (total.peak / expected.peak),
                          kminus1 = (extra.peak / expected.peak),
                          E = (extra.peak
                              / recurrence.data$Number_of_Reference_Discharges),
                          S = (extra.peak
                               / (recurrence.data$Number_of_Reference_Discharges
                                 + recurrence.data$Number_of_Event_Discharges)),
                          SI = (extra.peak / (total.count / 2)),
                          Peak.duration =  bounds[2] - bounds[1],
                          Peak.center = median(c(bounds[2], bounds[1])))


  return(Cumsum.Synch)

}
