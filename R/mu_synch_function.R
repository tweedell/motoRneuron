#' Determination of Motor Unit Synchronization by Various Methods of Cross
#' Correlation Histogram Examination
#'
#' @export
#' @importFrom methods "show"
#' @keywords recurrence, motor unit, synchronization, visual, Zscore, sigmax,
#'   cumulative sum
#' @description Calculates the time-domain synchronization indices CIS, k',
#'   k'-1, S, E, SI (detailed below) between the two input motor unit discharge
#'   trains. One or more methods of peak determination can be chosen ("Visual",
#'   "Zscore", or "Cumsum"). Chosen method functions are called individually and
#'   detailed descriptions are documented in their respective help files. Motor
#'   unit characteristic data and interspike intervals are automatically output
#'   with mu_synch.
#' @usage mu_synch(motor_unit_1, motor_unit_2, method = "Visual", order = 1,
#'   binwidth = 0.001, plot = F)
#' @param motor_unit_1 Numeric vectors of strictly increasing numbers denoting
#'   sequential discharge times of a motor unit or neuron or any strictly
#'   increasing point process.
#' @param motor_unit_2 Numeric vectors of strictly increasing numbers denoting
#'   sequential discharge times of a motor unit or neuron or any strictly
#'   increasing point process.
#' @param method Character vector indicating which methods of peak detection to
#'   use when quantifying synchronization. "Visual", "Zscore", and "Cumsum" are
#'   the options. Default is Visual.
#' @param order Numeric as a positive integer for the number of forward and
#'   backward orders for calculating recurrence times. Default = 1.
#' @param binwidth Numeric as a positive for the bin allocation size for
#'    histogram computation. Default = 0.001 or 1 ms.
#' @param plot T/F logical for outputting the cross correlation histogram.
#'   Default is FALSE.
#' @return A list of lists containing motor unit data (the names of each
#'   discharge train used, number of discharges, the interspike intervals (ISI),
#'   mean ISI, and the recurrence times associated with each order) and
#'   synchronization indices associated with chosen methods. 
#'   CIS = frequency of synchronized discharges. 
#'   k' = ratio of total discharges in peak to expected discharges in peak. 
#'   k'-1 = ratio of synchronized discharges to expected discharges in peak. 
#'   S = ratio of synchronized discharges to total number of discharges of both 
#'   motor units. 
#'   E = ratio of synchronized discharges to non-synchronized discharges. 
#'   SI = ratio of synchronized discharges to reference motor unit discharges.
#' @seealso visual_mu_synch, zscore_mu_synch, cumsum_mu_synch, sigmax_mu_synch
#' @references Keen, D.A., Chou, L., Nordstrom, M.A., Fuglevand, A.J. (2012)
#'   Short-term Synchrony in Diverse Motor Nuclei Presumed to Receive Different
#'   Extents of Direct Cortical Input. Journal of Neurophysiology 108: 3264-3275
#'
#'   Nordstrom, M.A., Fuglevand, A.J., Enoka, R.M. (1992) Estimating the
#'   Strength of Common Input to Human Motoneurons from the Cross-Correlogram.
#'   Journal of Physiology 453, pp. 547-574
#'
#'   DeFreitas, J.M., Beck, T.W., Xin, Y., Stock, M.S. (2013) Synchronization of
#'   Low- and High-Threshold Motor Units. Muscle & Nerve DOI 10.1002/mus.23978

mu_synch <- function(motor_unit_1, motor_unit_2, method = "Visual", order = 1,
                     binwidth = 0.001, plot = F){

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
  
  recurrence.data <- recurrence_intervals2(motor_unit_1, motor_unit_2, order)
  
  mean.reference.ISI <- recurrence.data$Mean_Reference_ISI
  
  # Create frequency table by binning recurrence times according to specfied bin
  # width using the mean reference ISI as the positive and negative boundaries.
  frequency.data <- unlist(recurrence.data[paste(1:order)])
  frequency.data <- frequency.data[frequency.data >= -mean.reference.ISI & 
                                     frequency.data <= mean.reference.ISI] 
  frequency.data <- as.vector(frequency.data)
  frequency.data <- motoRneuron::bin(frequency.data, binwidth = binwidth)
  
  synch.data <- list(Data = recurrence.data)
  
  if (plot) {
    
    show(plot_bins(frequency.data))
    
  }
  
  if ("Visual" %in% method) {

    synch.data[["Visual Indices"]] <- visual_mu_synch(motor_unit_1,
                                              motor_unit_2,
                                              order = order,
                                              binwidth = binwidth,
                                              get_data = F, 
                                              plot = F)[["Indices"]]

  }

  if ("Zscore" %in% method) {

    synch.data[["Zscore Indices"]] <- zscore_mu_synch(motor_unit_1,
                                              motor_unit_2,
                                              order = order,
                                              binwidth = binwidth,
                                              get_data = F,
                                              plot = F)[["Indices"]]

  }

  if ("Cumsum" %in% method) {

    synch.data[["Cumsum Indices"]] <- cumsum_mu_synch(motor_unit_1,
                                              motor_unit_2,
                                              order = order,
                                              binwidth = binwidth,
                                              get_data = F,
                                              plot = F)[["Indices"]]

   }

  return(synch.data)

}
