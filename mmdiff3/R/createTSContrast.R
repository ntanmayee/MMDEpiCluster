#' createTSContrast
#'
#' This function creates all contrasts for a time series experiments,
#' comparing all time points tp[i] with following time point tp[i+1].
#' If condition is set, only samples with given condition are considered.
#'
#' @inheritParams getPeakReads
#' @inheritParams DBAmmd-Accessors
#' @param condition (DEFAULT: NULL)
#' @param tp time points (DEFAULT: NULL)
#' @export
#'
createTSContrast <- function(MD,condition=NULL,tp=NULL){

  samples <- Samples(MD)
  if (is.null(tp)){
    if (is.null(condition)){
      tp <- sort(as.numeric(unique(samples[['Time']])))
    } else {
      idx <- which(samples$Condition==condition)
      samples.c <- samples[idx,]
      tp <- sort(as.numeric(unique(samples.c[['Time']])))
    }
  }


  group <- rep(FALSE,nrow(samples))
  names(group) <- samples$SampleID
  contrasts <- vector(mode='list',length=length(tp)-1)
  for (t in 1:(length(tp)-1)){
    group1 <- group
    group2 <- group
    if (is.null(condition)){
      idx1 <- which(samples[['Time']]==tp[t])
      idx2 <- which(samples[['Time']]==tp[t+1])
      name1 <- tp[t]
      name2 <- tp[t+1]
    } else {
      idx1 <- which(samples[['Condition']]==condition &
                      samples[['Time']]==tp[t])
      idx2 <- which(samples[['Condition']]==condition &
                      samples[['Time']]==tp[t+1])
      name1 <- paste(condition,':T=', tp[t],sep='')
      name2 <- paste(condition,':T=', tp[t+1],sep='')
    }

    group1[idx1] <- TRUE
    group2[idx2] <- TRUE


    contrasts[[t]] <- list(group1=group1,
                           group2=group2,
                           name1=name1,
                           name2=name2)
  }
  return(contrasts)
}



#' createTSContrastBycondition
#'
#' This function creates contrasts between conditions[1] and condistions[2] for
#' all time points tp of a time series experiments.
#'
#' @inheritParams getPeakReads
#' @inheritParams DBAmmd-Accessors
#' @param conditions If NULL first two conditions are considered. (DEFAULT: NULL)
#' @param tp time points. If NULL all available time points are considred
#'  (DEFAULT: NULL)
#' @export
#'
createTSContrastBycondition <- function(MD,conditions=NULL,tp=NULL){

  samples <- Samples(MD)
  if (is.null(conditions)){
    conditions <- unique(samples[['Condition']])[1:2]
  }
  if (length(conditions)!=2){
    stop('conditions must have length 2')
  }
  idx <- which(is.element(samples$Condition,conditions))
  samples.c <- samples[idx,]



  if (is.null(tp)){
    tp <- sort(as.numeric(unique(samples.c[['Time']])))
  }


  group <- rep(FALSE,nrow(samples))
  names(group) <- samples$SampleID

  contrasts <- vector(mode='list',length=length(tp))
  for (t in 1:(length(tp))){
    group1 <- group
    group2 <- group

    idx1 <- which(samples[['Condition']]==conditions[1] &
                    samples[['Time']]==tp[t])
    idx2 <- which(samples[['Condition']]==conditions[2] &
                    samples[['Time']]==tp[t])
    name1 <- paste(conditions[1],':T=', tp[t],sep='')
    name2 <- paste(conditions[2],':T=', tp[t],sep='')

    group1[idx1] <- TRUE
    group2[idx2] <- TRUE

    contrasts[[t]] <- list(group1=group1,
                           group2=group2,
                           name1=name1,
                           name2=name2)
  }
  return(contrasts)
}



