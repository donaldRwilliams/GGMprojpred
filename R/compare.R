#' compare
#'
#' @param Estimate estimated graph
#' @param True the true graph
#'
#' @return results in data frame
#'
#' @export
compare <- function(Estimate, True){

  True <- as.matrix(True)
  Estimate <- as.matrix(Estimate)

  # True Negative
  TN <- ifelse(Estimate[upper.tri(Estimate)] == 0 & True[upper.tri(True)] == 0, 1, 0); TN <- sum(TN)
  # False Positive
  FP <- ifelse(Estimate[upper.tri(Estimate)] == 0 & True[upper.tri(True)] != 0, 1, 0); FP <- sum(FP)
  # True Positive
  TP <- ifelse(Estimate[upper.tri(Estimate)] != 0 & True[upper.tri(True)] != 0, 1, 0); TP <- sum(TP)
  # False Negatives
  FN <- ifelse(Estimate[upper.tri(Estimate)] != 0 & True[upper.tri(True)] == 0, 1, 0); FN <- sum(FN)

  Specificity <- TN/(TN + FP)
  Sensitivity <- TP/(TP + FN)
  Precision <- TP/(TP + FP)

  Recall <- TP / (TP + FN)

  F1_score <- 2 * ((Precision * Recall) / (Precision + Recall))

  MCC <- (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))


  results <- c(Specificity, Sensitivity, Precision, Recall,  F1_score, MCC)
  results_name <- c("Specificity", "Sensitivity", "Precision", "Recall",  "F1_score", "MCC")
  results <- cbind.data.frame(measure = results_name, score = results)
  list(results = results)


}




