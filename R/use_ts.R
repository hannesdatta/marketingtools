# decide whether a time series is appropriate (per brand in each market)
#' @export
use_ts <- function(x) {
  # exclude series from model estimation under certain conditions (see below)
  # returns TRUE for valid series, and FALSE for invalid series

  # invalid series if all values are NA
  if (all(is.na(x))) return(F)

  tmp = x
  .tmp=table(tmp)
  .tmp=.tmp[which(.tmp==max(.tmp))][1]

  # invalid series if it only contains one unique value (no variation)
  if (length(unique(tmp))<=1) return(F)

  # invalid series if 95% are the same value
  if (.tmp/length(tmp)>.95) return(F)


  # invalid series if it has less than x months of non-zero observations
  if (length(which(!tmp==min(tmp, na.rm=T)))<12) return(F)

  return(T)
}
