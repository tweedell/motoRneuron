#' Calculate Recurrence Intervals Between Motor Units or Neuron Discharge Trains
#'
#' Take two vectors representing time points of motor unit or neuron discharges
#' and calculate multi-order recurrence times between them. This function will
#' return recurrence times in whatever unit is input but will only accept
#' numeric vectors (e.g. 0.01 sec, 25 ms, or 17.5 minutes).
#'
#' @export
#' @keywords recurrence, motor unit, action potentials, discharge trains
#' @usage recurrence_intervals(motor_unit_1, motor_unit_2, order = 1)
#' @param motor_unit_1,motor_unit_2 Numeric vectors of strictly increasing
#'   numbers denoting sequential discharge times of a motor unit or neuron or
#'   any strictly increasing point process.
#' @param order Numeric as a positive integer for the number of forward and
#'   backward orders for calculating recurrence times. Default = 1. 
#' @return A list of lists containing the names of each discharge train used,
#'   number of discharges, the interspike intervals (ISI), mean ISI, and the
#'   recurrence times associated with each order.

recurrence_intervals <- function(motor_unit_1, motor_unit_2, order = 1) {


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

  # reference (ref) and event motor units assigned according to which motor unit
  # has more firings (length()). ISI = InterSpike Intervals
  if (length(motor_unit_1) < length(motor_unit_2)) {
    ref.name <- deparse(substitute(motor_unit_1))
    event.name <- deparse(substitute(motor_unit_2))
    ref.MU <- motor_unit_1
    event.MU <- motor_unit_2
    ref.MU.ISI <- diff(motor_unit_1)
    event.MU.ISI <- diff(motor_unit_2)
    mean.ref.ISI <- round(mean(ref.MU.ISI), digits = 3)
    mean.event.ISI <- round(mean(event.MU.ISI), digits = 3)
  } else {
    ref.name <- deparse(substitute(motor_unit_2))
    event.name <- deparse(substitute(motor_unit_1))
    ref.MU <- motor_unit_2
    event.MU <- motor_unit_1
    ref.MU.ISI <- diff(motor_unit_2)
    event.MU.ISI <- diff(motor_unit_1)
    mean.ref.ISI <- round(mean(ref.MU.ISI), digits = 3)
    mean.event.ISI <- round(mean(event.MU.ISI), digits = 3)
  }

  # Initialize list of motor unit characteristics
  MU.names <- list(Reference_Unit = ref.name,
                   Number_of_Reference_Discharges = length(ref.MU),
                   Reference_ISI = ref.MU.ISI,
                   Mean_Reference_ISI = mean.ref.ISI,
                   Event_Unit = event.name,
                   Number_of_Event_Discharges = length(event.MU),
                   Event_ISI = event.MU.ISI,
                   Mean_Event_ISI = mean.event.ISI,
                   Duration = max(ref.MU, event.MU) - min(ref.MU, event.MU))
  
  # Initialize list for each number 1 to the order specified
  lags <- vector('list', order)

  # Iteratively assess each reference motor unit discharge time for discharge
  # times of the event motor unit around that time point and store recurrence
  # intervals in vector (lags).
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

#' Bin Recurrence Times
#'
#' Discretizes recurrence times vector into a data frame according to specified
#' binwidth parameter. A recurrence time of "0" indicates two events happening
#' simultaneously. The "0" bin is treated as the center.
#'
#' @export
#' @keywords bin, group, motor unit, action potentials, discharge trains
#' @usage bin(recurrences, binwidth = 0.001)
#' @param recurrences Numeric vector representing recurrence times.
#' @param binwidth Numeric. Default = 0.001 (0.001 second or 1 ms). This must
#'   have the same significant digits as the recurrences parameter.
#' @return data frame containing frequency data.

bin <- function(recurrences, binwidth = 0.001){

  if (!is.numeric(recurrences)) {
    stop("'recurrences' must be a numeric vector.")
  }

  # Subset all positive intervals
  pos <- recurrences[recurrences>=0]
  # establish break points for discretization based on binwidth
  pos_br <- seq(0, (max(pos) + binwidth), by = binwidth)
  # bin according to break points
  pos_binned <- cut(pos, breaks = pos_br, right = T, include.lowest = T) 
  pos_binned <- table(pos_binned)
  pos_binned <- as.data.frame(pos_binned)
  colnames(pos_binned) <- c("Bin", "Freq")
  
  # Subset all negative intervals
  neg <- recurrences[recurrences<0]
  # establish break points for discretization based on binwidth
  neg_br <- seq(0, (min(neg) - binwidth), by = -(binwidth))
  # bin according to break points
  neg_binned <- cut(neg, breaks = neg_br, right = F, include.lowest = F) 
  neg_binned <- table(neg_binned) 
  neg_binned <- as.data.frame(neg_binned)
  colnames(neg_binned) <- c("Bin", "Freq")

  # create and return as data frame
  count_data <- rbind(neg_binned, pos_binned)
  count_data$Bin <- c(rev(neg_br[2:length(neg_br)]), pos_br[1:(length(pos_br)-1)])

  return(count_data)
}

#' Plot binned recurrence times as a histogram
#'
#' Takes data frame containing bins labels and frequency data to produce
#' histogram.
#'
#' @export
#' @import ggplot2 
#' @keywords bin, plot, motor unit, discharge trains, histogram
#' @usage plot_bins(binned_data)
#' @param binned_data data frame containing frequency data produced by bin().

plot_bins <- function(binned_data) {

  Mean <- mean(as.numeric(binned_data$Freq))

  ggplot(binned_data) +
  geom_col(mapping = aes_(x = as.numeric(binned_data$Bin), y = as.numeric(binned_data$Freq)),
           width = 0.001) +
  theme_classic() +
  xlab("Recurrence time (sec)") +
  ylab("Frequency") +
  ggtitle("Cross Correlation Histogram") +
  geom_hline(yintercept = Mean, color = "red") +
  geom_text(aes_(min(as.numeric(binned_data$Bin)),
                Mean, label = "Mean",
                vjust = -0.5, hjust = -0.01), color = "red", size = 3)

}