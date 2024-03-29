#' @title H2x2Factorial Table
#'
#' @description  The function \code{table.H2x2Factorial} outputs a data frame that summarizes the required number of clusters and the predicted
#' power based on a constellation of design parameters. This function is useful when the user wants a series of table-format predictions
#' based on varying design parameters including mean cluster size (m_bar), intraclass correlation coefficient (rho), and coefficient of variation of the cluster sizes (CV).
#'
#' @usage
#' table.H2x2Factorial(power=0.8, alpha=0.05,
#'                     pi_x=0.5, pi_z=0.5,
#'                     delta_x, delta_z, delta_xz, sigma2_y=1,
#'                     m_bar, CV, rho,
#'                     estimand="controlled", test="cluster", correction=FALSE,
#'                     max_n=1e8, seed_mix=NULL, size_mix=1e4,
#'                     verbose=TRUE)
#'
#' @param power a numeric value between 0 and 1 as the desired power level for sample size estimation. Default is \code{0.8}.
#' @param alpha a numeric value between 0 and 1 as the type I error rate. Default is \code{0.05}.
#' @param pi_x a numeric value between 0 and 1 as the proportion of clusters randomized to the cluster-level treatment. Default is \code{0.5}, representing a balanced allocation.
#' @param pi_z a numeric value between 0 and 1 as the proportion of individuals randomized to the individual-level treatment within each cluster. Default is \code{0.5}, representing a balanced allocation.
#' @param delta_x a nonzero numeric value for the (unstandardized) effect size of the marginal cluster-level treatment effect. Default is \code{0.25}, which is the hypothetical value for the example in the referenced paper.
#' @param delta_z a nonzero numeric value for the (unstandardized) effect size of the marginal individual-level treatment effect. Default is \code{0.33}, which is the hypothetical value for the example in the referenced paper.
#' @param delta_xz a nonzero numeric value for the (unstandardized) effect size of the interaction effect of the two treatments. Default is \code{0.3}, which is the hypothetical value for the example in the referenced paper.
#' @param sigma2_y a positive numeric value for the total variance of the continuous outcome. Default is \code{1}.
#' @param m_bar a vector of numeric values larger than 2 for a series of mean cluster sizes.
#' @param CV a vector of positive numeric values for a series of coefficients of variation of the cluster sizes.
#' @param rho a vector of numeric values between 0 and 1 for a series of intraclass correlation coefficients.
#' @param estimand a character argument indicating the type of treatment effect estimand. Supported values include \code{"controlled"}
#' (controlled or main effect estimand) and \code{"natural"} (natural or marginal effect estimand). Default is \code{"controlled"}.
#' @param test a character argument indicating the type of hypothesis test of interest. Supported values include
#' \code{"cluster"} (test for marginal cluster-level treatment effect), \code{"individual"} (test for marginal individual-level treatment effect),
#' \code{"interaction"} (interaction test for the two treatments), \code{"joint"} (joint test for the two marginal treatment effects),
#' \code{"I-U"} (intersection-union test for the two marginal effects). Default is \code{"cluster"}.
#' @param correction a logical argument indicating whether a finite sample correction should be used. Default is \code{FALSE}.
#' @param max_n an optional setting of a maximum number of clusters, which is only functional under \code{test="cluster"}, \code{"joint"}, or \code{"I-U"}. Default is \code{1e8}.
#' @param seed_mix an optional setting of a seed for conducting the simulation-based testing under a mixed distribution, which is only functional under \code{test="joint"}. Default is \code{NULL}.
#' @param size_mix a pre-specified size for the mixed distribution in the simulation-based procedure, which is only needed under \code{test="joint"}. Default is \code{1e4}.
#' @param verbose a logical argument indicating whether the parameter reiterations and supplementary messages should be presented or suppressed. Default is \code{TRUE}.
#'
#' @details
#' If the user further requires a vector of \code{power} or other parameters like \code{pi_x}, which invokes the need for multiple tables,
#' an external loop could be easily written using this function to produce multiple data frames.
#'
#' @return \code{table.H2x2Factorial} returns a data frame with inputs of \code{m_bar}, \code{rho}, and \code{CV} varied in a factorial setting, the predicted number of clusters \code{n} under the power requirement,
#' and the actual power \code{predicted.power} the estimated sample size can help to achieve, with some suppressible messages.
#'
#'
#' @export
#'
#' @examples
#' #Make a result table by providing three mean cluster sizes, three CV, and three ICC
#' table.cluster <- table.H2x2Factorial(delta_x=0.2, delta_z=0.1,
#'                                      m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1),
#'                                      estimand="controlled", test="cluster", verbose=FALSE)
#' table.cluster
#'
#' @import mvtnorm
#' @importFrom stats qnorm pnorm qt pt qchisq pchisq rchisq qf pf rf quantile
#'
table.H2x2Factorial <- function(power=0.8, alpha=0.05,
                                pi_x=0.5, pi_z=0.5,
                                delta_x=0.25, delta_z=0.33, delta_xz=0.3,
                                sigma2_y=1,
                                m_bar, CV, rho,
                                estimand="controlled",
                                test="cluster",
                                correction=FALSE,
                                max_n=1e8,
                                seed_mix=NULL,
                                size_mix=1e4,
                                verbose=TRUE) {

  if (!is.numeric(power) || power <= 0 || power >= 1 || length(power)!=1)
    stop('Target power must be a single number in (0,1)')

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha)!=1)
    stop('Type I error rate must be a single number in (0,1)')

  if (!is.numeric(sigma2_y) || sigma2_y <= 0 || length(alpha)!=1)
    stop('Total variance of the outcome must be a single positive number')

  for (i in 1:length(m_bar)){
    if (!is.numeric(m_bar[i]) || m_bar[i] <= 2)
      stop('Each mean cluster size provided must be a number larger than 2')
  }
  if (length(unique(m_bar)) != length(m_bar)){
    m_bar <- unique(m_bar)
    m_bar.list <- m_bar[1]
    for (i in 2:length(m_bar)){
      m_bar.list <- paste0(m_bar.list, ",", m_bar[i])
    }
    warning(paste0("Duplicated elements of the m_bar input is deleted. m_bar input is changed to:\n(", m_bar.list, ")\n"))
  }

  for (i in 1:length(CV)){
    if (!is.numeric(CV[i]) || CV[i] < 0)
      stop('Each coefficient of variation of the cluster sizes provided must be a positive number')
  }
  if (length(unique(CV)) != length(CV)){
    CV <- unique(CV)
    CV.list <- CV[1]
    for (i in 2:length(CV)){
      CV.list <- paste0(CV.list, ",", CV[i])
    }
    warning(paste0("Duplicated elements of the CV input is deleted. CV input is changed to:\n(", CV.list, ")\n"))
  }

  for (i in 1:length(rho)){
    if (!is.numeric(rho[i]) || rho[i] < 0 || rho[i] >= 1)
      stop('Each intraclass correlation coefficient provided must be numeric in [0,1)')
  }
  if (length(unique(rho)) != length(rho)){
    rho <- unique(rho)
    rho.list <- rho[1]
    for (i in 2:length(rho)){
      rho.list <- paste0(rho.list, ",", rho[i])
    }
    warning(paste0("Duplicated elements of the rho input is deleted. rho input is changed to:\n(", rho.list, ")\n"))
  }

  if ( !(estimand %in% c("controlled", "natural")) || length(estimand)!=1)
    stop('Type of treatment effect estimand should be a single choice from "controlled" and "natural"')

  if ( !(test %in% c("cluster", "individual", "interaction", "joint", "I-U")) || length(test)!=1)
    stop('Type of hypothesis tests should be a single choice from "cluster", "individual", "interaction", "joint", and "I-U"')

  if (test=="cluster"){
    if (!is.numeric(delta_x) || delta_x == 0 || length(delta_x)!=1)
      stop('Effect size of the marginal cluster-level treatment effect must be a nonzero number')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1 || length(pi_x)!=1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment arm must be a single number in (0,1)')

  } else if (test=="individual"){
    if (!is.numeric(delta_z) || delta_z == 0 || length(delta_z)!=1)
      stop('Effect size of the marginal individual-level treatment effect must be a nonzero number')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1 || length(pi_z)!=1)
      stop('Proportion of individuals that are randomized to the individual-level treatment arm must be a single number in (0,1)')
    if (correction==TRUE)
      message('No finite-sample correction will be done for the test for marginal individual-level treatment effect due to adequate degrees of freedom')

  } else if (test=="interaction"){
    if (!is.numeric(delta_xz) || delta_xz == 0 || length(delta_xz)!=1)
      stop('Effect size of the interaction effect must be a nonzero number')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1 || length(pi_x)!=1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment arm must be a single number in (0,1)')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1 || length(pi_z)!=1)
      stop('Proportion of individuals that are randomized to the individual-level treatment arm must be a single number in (0,1)')
    if (correction==TRUE)
      message('No finite-sample correction will be done for the interaction test due to adequate degrees of freedom')

  } else if (test=="joint"){
    if (!is.numeric(delta_x) || delta_x == 0 || length(delta_x)!=1)
      stop('Effect size of the marginal cluster-level treatment effect must be a nonzero number')
    if (!is.numeric(delta_z) || delta_z == 0 || length(delta_z)!=1)
      stop('Effect size of the marginal individual-level treatment effect must be a nonzero number')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1 || length(pi_x)!=1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment arm must be a single number in (0,1)')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1 || length(pi_z)!=1)
      stop('Proportion of individuals that are randomized to the individual-level treatment arm must be a single number in (0,1)')

  } else if (test=="I-U"){
    if (!is.numeric(delta_x) || delta_x == 0 || length(delta_x)!=1)
      stop('Effect size of the marginal cluster-level treatment effect must be a nonzero number')
    if (!is.numeric(delta_z) || delta_z == 0 || length(delta_z)!=1)
      stop('Effect size of the marginal individual-level treatment effect must be a nonzero number')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1 || length(pi_x)!=1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment arm must be a single number in (0,1)')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1 || length(pi_z)!=1)
      stop('Proportion of individuals that are randomized to the individual-level treatment arm must be a single number in (0,1)')
  }

  if (!is.logical(correction))
    stop('Finite sample correction indicator should be a logical argument')


  if (!is.numeric(max_n) || max_n <= 0 || length(max_n)!=1)
    stop('Maximum number of clusters must be a positive number')

  if (!is.null(seed_mix)){
    if (!is.numeric(seed_mix) || length(seed_mix)!=1)
      stop('User-defined seed under the finite-sample corrected joint test must be numeric')
  }

  if (!is.numeric(size_mix) || size_mix <= 0 || length(size_mix)!=1)
    stop('Sample size for simulating the mix distribution under the finite-sample corrected joint test must be a positive number')

  if (!is.logical(verbose))
    stop('Message presentation indicator should be a logical argument')


  #Effect sizes might be negative
  delta_x <- abs(delta_x)
  delta_z <- abs(delta_z)
  delta_xz <- abs(delta_xz)

  m_bar.vector <- m_bar
  CV.vector <- CV
  rho.vector <- rho

  table <- NULL
  for (m_bar.i in 1:length(m_bar.vector)){
    for (rho.i in 1:length(rho.vector)){
      for (CV.i in 1:length(CV.vector)){
        table <- rbind(table, c(m_bar[m_bar.i], rho[rho.i], CV[CV.i]))
      }
    }
  }
  m_bar <- table[,1]
  rho <- table[,2]
  CV <- table[,3]


  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)

  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  omega_xz <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_x*(1-pi_x)*pi_z*(1-pi_z))

  omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  omega_3 <- sigma2_y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  omega_23 <- pi_z*omega_3

  #Check the omega to be positive
  for (i in length(omega_x)){
    if (omega_x[i]<=0 || omega_z[i]<=0 || omega_xz[i]<=0 || omega_2[i]<=0 || omega_3[i]<=0)
      stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")
  }


  #Re-iterate the given effect sizes and the chosen test
  if (verbose==TRUE){
    if (estimand=="controlled"){
      cat('Type of treatment effect estimand:\nContolled (main) effect')

      if (test=="cluster"){
        cat('\nType of hypothesis test:\nSeparate test for the cluster-level treatment effect')
        cat(paste0('\nEffect size:\n', delta_x, " for the controlled effect of the cluster-level treatment"))

      } else if (test=="individual"){
        cat('\nType of hypothesis test:\nSeparate test for the individual-level treatment effect')
        cat(paste0('\nEffect size:\n', delta_z, " for the controlled effect of the individual-level treatment"))

      } else if (test=="interaction"){
        cat('\nType of hypothesis test:\nInteraction test')
        cat(paste0('\nEffect size:\n', delta_xz, " for the interaction effect"))

      } else if (test=="joint"){
        cat('\nType of hypothesis test:\nJoint test')
        cat(paste0('\nEffect sizes:\n', delta_x, " for the controlled effect of the cluster-level treatmen\n", delta_z, " for the controlled effect of the individual-level treatment"))

      } else if (test=="I-U"){
        cat('\nType of hypothesis test:\nIntersection-union test')
        cat(paste0('\nEffect sizes:\n', delta_x, " for the controlled effect of the cluster-level treatment\n", delta_z, " for the controlled effect of the individual-level treatment"))

      }

    } else if (estimand=="natural"){
      cat('Type of treatment effect estimand:\nNatural (marginal) effect')

      if (test=="cluster"){
        cat('\nType of hypothesis test:\nSeparate test for the cluster-level treatment effect')
        cat(paste0('\nEffect size:\n', delta_x, " for the natural effect of the cluster-level treatment"))

      } else if (test=="individual"){
        cat('\nType of hypothesis test:\nSeparate test for the individual-level treatment effect')
        cat(paste0('\nEffect size:\n', delta_z, " for the natural effect of the individual-level treatment"))

      } else if (test=="interaction"){
        cat('\nType of hypothesis test:\nInteraction test')
        cat(paste0('\nEffect size:\n', delta_xz, " for the interaction effect"))

      } else if (test=="joint"){
        cat('\nType of hypothesis test:\nJoint test')
        cat(paste0('\nEffect sizes:\n', delta_x, " for the natural effect of the cluster-level treatmen\n", delta_z, " for the natural effect of the individual-level treatment"))

      } else if (test=="I-U"){
        cat('\nType of hypothesis test:\nIntersection-union test')
        cat(paste0('\nEffect sizes:\n', delta_x, " for the natural effect of the cluster-level treatment\n", delta_z, " for the natural effect of the individual-level treatment"))

      }
    }

  }


  if (estimand=="controlled"){

    if (test=="cluster"){

      if(correction==FALSE){

        if (verbose==TRUE){
          cat("\nA Wald z-test is used without finite-sample correction\n")
        }
        n <- (z_a+z_b)^2*omega_2/delta_x^2
        n.final <- ceiling(n)
        pred.power <- pnorm(sqrt(n.final*delta_x^2/omega_2)-z_a)

      } else if (correction==TRUE){

        if (verbose==TRUE){
          cat("\nA t-test is used for finite-sample correction\n")
        }
        cluster.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
            ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
          n <- 2
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            try.power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_2/n), lower.tail = F)
            + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_2/n))

          }
          return(c(n, try.power))
        }

        cluster.pred <- NULL
        for (i in 1:nrow(table)){
          cluster.pred <- rbind(cluster.pred, cluster.n(parameter=unlist(table[i,])))
        }
        n.final <- cluster.pred[,1]
        pred.power <- cluster.pred[,2]
      }
    }


    if (test=="individual"){
      if (verbose==TRUE){
        cat("\nA Wald z-test is automatically used\n")
      }
      n <- (z_a+z_b)^2*omega_2/delta_z^2
      n.final <- ceiling(n)
      pred.power <- pnorm(sqrt(n.final*delta_z^2/omega_2)-z_a)
    }


    if (test=="interaction"){
      if (verbose==TRUE){
        cat("\nA Wald z-test is automatically used\n")
      }
      n <- (z_a+z_b)^2*omega_xz/delta_xz^2
      n.final <- ceiling(n)
      pred.power <- pnorm(sqrt(n.final*delta_xz^2/omega_xz)-z_a)
    }


    if (test=="joint"){

      if(correction==FALSE){

        if (verbose==TRUE){
          cat("\nA Chi-square test is used without finite-sample correction\n")
        }
        joint.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
            ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
          omega_3 <- sigma2_y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
          omega_23 <- pi_z*omega_3
          n <- 1
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            beta23 <- c(delta_x, delta_z)
            omega <- matrix(c(omega_2, omega_23, omega_23, omega_3), nrow=2, byrow=T)
            theta <- n*t(beta23) %*% solve(omega) %*% beta23
            try.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)

          }
          return(c(n, try.power))
        }

        joint.pred <- NULL
        for (i in 1:nrow(table)){
          joint.pred <- rbind(joint.pred, joint.n(parameter=unlist(table[i,])))
        }

        n.final <- joint.pred[,1]
        pred.power <- joint.pred[,2]

      } else if (correction==TRUE){

        if (verbose==TRUE){
          cat("\nA F-test is used for finite-sample correction\n")
        }
        joint.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
            ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
          omega_3 <- sigma2_y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
          omega_23 <- pi_z*omega_3
          n <- 2
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            beta23 <- c(delta_x, delta_z)
            omega <- matrix(c(omega_2, omega_23, omega_23, omega_3), nrow=2, byrow=T)
            theta <- n*t(beta23) %*% solve(omega) %*% beta23
            try.power <- pf(qf(1-a, 2, n-2), 2, n-2, ncp = theta, lower.tail = F)

          }
          return(c(n, try.power))
        }

        joint.pred <- NULL
        for (i in 1:nrow(table)){
          joint.pred <- rbind(joint.pred, joint.n(parameter=unlist(table[i,])))
        }
        n.final <- joint.pred[,1]
        pred.power <- joint.pred[,2]
      }
    }


    if (test=="I-U"){

      if(correction==FALSE){
        if (verbose==TRUE){
          cat("\nA z-based intersection-union test is used without finite-sample correction\n")
        }
        IU.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
            ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
          omega_3 <- sigma2_y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
          omega_23 <- pi_z*omega_3
          n <- 1
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            mean_W <- c(delta_x/sqrt(omega_2/n), delta_z/sqrt(omega_3/n))
            cov_W <- (omega_23/n) / (sqrt(omega_2/n)*sqrt(omega_3/n))
            sigma_W <- matrix(c(1, cov_W, cov_W, 1), nrow=2, byrow=T)

            try.power <- pmvnorm(lower=rep(qnorm(1-a/2),2), upper=rep(Inf,2), mean=mean_W, sigma=sigma_W) +
              pmvnorm(lower=c(qnorm(1-a/2),-Inf), upper=c(Inf,qnorm(a/2)), mean=mean_W, sigma=sigma_W) +
              pmvnorm(lower=rep(-Inf,2), upper=rep(qnorm(a/2),2), mean=mean_W, sigma=sigma_W) +
              pmvnorm(lower=c(-Inf,qnorm(a/2)), upper=c(qnorm(a/2),Inf), mean=mean_W, sigma=sigma_W)
          }
          return(c(n, try.power))
        }

        IU.pred <- NULL
        for (i in 1:nrow(table)){
          IU.pred <- rbind(IU.pred, IU.n(parameter=unlist(table[i,])))
        }
        n.final <- IU.pred[,1]
        pred.power <- IU.pred[,2]

      } else if (correction==TRUE){
        if (verbose==TRUE){
          cat("\nA t-based intersection-union test is used for finite-sample correction\n")
        }
        IU.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
            ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
          omega_3 <- sigma2_y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
          omega_23 <- pi_z*omega_3
          n <- 2
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            mean_W <- c(delta_x/sqrt(omega_2/n), delta_z/sqrt(omega_3/n))
            cov_W <- (omega_23/n) / (sqrt(omega_2/n)*sqrt(omega_3/n))
            sigma_W <- matrix(c(1, cov_W, cov_W, 1), nrow=2, byrow=T)

            try.power <- pmvt(df=n-2, lower=c(qt(1-a/2, n-2), qt(1-a/2, n-2)), upper=rep(Inf,2), delta=mean_W, sigma=sigma_W) +
              pmvt(df=n-2, lower=c(qt(1-a/2, n-2),-Inf), upper=c(Inf,qt(a/2, n-2)), delta=mean_W, sigma=sigma_W) +
              pmvt(df=n-2, lower=rep(-Inf,2), upper=c(qt(a/2, n-2), qt(a/2, n-2)), delta=mean_W, sigma=sigma_W) +
              pmvt(df=n-2, lower=c(-Inf,qt(a/2, n-2)), upper=c(qt(a/2, n-2),Inf), delta=mean_W, sigma=sigma_W)
          }
          return(c(n, try.power))
        }

        IU.pred <- NULL
        for (i in 1:nrow(table)){
          IU.pred <- rbind(IU.pred, IU.n(parameter=unlist(table[i,])))
        }

        n.final <- IU.pred[,1]
        pred.power <- IU.pred[,2]
      }
    }

  } else if (estimand=="natural"){

    if (test=="cluster"){

      if(correction==FALSE){

        if (verbose==TRUE){
          cat("\nA Wald z-test is used without finite-sample correction\n")
        }
        n <- (z_a+z_b)^2*omega_x/delta_x^2
        n.final <- ceiling(n)
        pred.power <- pnorm(sqrt(n.final*delta_x^2/omega_x)-z_a)

      } else if (correction==TRUE){

        if (verbose==TRUE){
          cat("\nA t-test is used for finite-sample correction\n")
        }
        cluster.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
          omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
          n <- 2
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            try.power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n), lower.tail = F) + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n))

          }
          return(c(n, try.power))
        }

        cluster.pred <- NULL
        for (i in 1:nrow(table)){
          cluster.pred <- rbind(cluster.pred, cluster.n(parameter=unlist(table[i,])))
        }
        n.final <- cluster.pred[,1]
        pred.power <- cluster.pred[,2]
      }
    }


    if (test=="individual"){
      if (verbose==TRUE){
        cat("\nA Wald z-test is automatically used\n")
      }
      n <- (z_a+z_b)^2*omega_z/delta_z^2
      n.final <- ceiling(n)
      pred.power <- pnorm(sqrt(n.final*delta_z^2/omega_z)-z_a)
    }


    if (test=="interaction"){
      if (verbose==TRUE){
        cat("\nA Wald z-test is automatically used\n")
      }
      n <- (z_a+z_b)^2*omega_xz/delta_xz^2
      n.final <- ceiling(n)
      pred.power <- pnorm(sqrt(n.final*delta_xz^2/omega_xz)-z_a)
    }


    if (test=="joint"){

      if(correction==FALSE){

        if (verbose==TRUE){
          cat("\nA Chi-square test is used without finite-sample correction\n")
        }
        joint.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
          eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
          omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
          omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
          n <- 1
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            theta <- n*(delta_x^2/omega_x + delta_z^2/omega_z)
            try.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)

          }
          return(c(n, try.power))
        }

        joint.pred <- NULL
        for (i in 1:nrow(table)){
          joint.pred <- rbind(joint.pred, joint.n(parameter=unlist(table[i,])))
        }

        n.final <- joint.pred[,1]
        pred.power <- joint.pred[,2]

      } else if (correction==TRUE){

        if (verbose==TRUE){
          cat("\nA simulation-based mixed F-Chi-square test is used for finite-sample correction\n")
        }
        joint.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
          eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
          omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
          omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
          n <- 2
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            set.seed(seed_mix)
            n <- n+1
            f.distn <- rf(size_mix, 1, n-2)
            chisq.distn <- rchisq(size_mix, 1)
            mix.distn <- f.distn + chisq.distn
            crt.value <- quantile(mix.distn, 1-a)

            nc.f.distn <- rf(size_mix, 1, n-2, ncp = n*delta_x^2/omega_x)
            nc.chisq.distn <- rchisq(size_mix, 1, ncp = n*delta_z^2/omega_z)
            nc.mix.distn <- nc.f.distn + nc.chisq.distn
            try.power <- mean(nc.mix.distn>crt.value)

          }
          return(c(n, try.power))
        }

        joint.pred <- NULL
        for (i in 1:nrow(table)){
          joint.pred <- rbind(joint.pred, joint.n(parameter=unlist(table[i,])))
        }
        n.final <- joint.pred[,1]
        pred.power <- joint.pred[,2]
      }
    }


    if (test=="I-U"){

      if(correction==FALSE){
        if (verbose==TRUE){
          cat("\nA z-based intersection-union test is used without finite-sample correction\n")
        }
        IU.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
          eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
          omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
          omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
          n <- 1
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            wmean.c <- sqrt(n)*delta_x/sqrt(omega_x)
            wmean.i <- sqrt(n)*delta_z/sqrt(omega_z)
            try.power <- pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) + pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(a/2), mean = wmean.i) + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(a/2), mean = wmean.i)

          }
          return(c(n, try.power))
        }

        IU.pred <- NULL
        for (i in 1:nrow(table)){
          IU.pred <- rbind(IU.pred, IU.n(parameter=unlist(table[i,])))
        }
        n.final <- IU.pred[,1]
        pred.power <- IU.pred[,2]

      } else if (correction==TRUE){
        if (verbose==TRUE){
          cat("\nA mixed t- and z-based intersection-union test is used for finite-sample correction\n")
        }
        IU.n <- function(parameter){
          m_bar <- parameter[1]
          rho <- parameter[2]
          CV <- parameter[3]
          eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
          eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
          omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
          omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
          n <- 2
          try.power <- 0
          while ((try.power < power) & (n < max_n)){
            n <- n+1
            c.ncp <- sqrt(n)*delta_x/sqrt(omega_x)
            i.mean <- sqrt(n)*delta_z/sqrt(omega_z)
            try.power <- pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) + pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(a/2), mean = i.mean) + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(a/2), mean = i.mean)

          }
          return(c(n, try.power))
        }

        IU.pred <- NULL
        for (i in 1:nrow(table)){
          IU.pred <- rbind(IU.pred, IU.n(parameter=unlist(table[i,])))
        }

        n.final <- IU.pred[,1]
        pred.power <- IU.pred[,2]
      }
    }

  }


  table <- data.frame(cbind(m_bar, rho, CV, n.final, pred.power))
  names(table) <- c("m_bar", "rho", "CV", "n", "predicted power")
  return(table)
}



