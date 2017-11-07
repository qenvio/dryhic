fitcnv <- function (x) {

   # Constants.
   MAXITER = 500
   TOLERANCE = 1e-3
   MCHNEPS = .Machine$double.eps

###############################################
#              OPTION PROCESSING              #
###############################################

   stopifnot(is.vector(x))
   x[x <= quantile(x, .05, na.rm=TRUE)/2] = NA

   # Valid copy numbers. The last is open-ened (8 or more).
   CN <- c(1, 2, 3, 4, 5, 6, 7)

   # Dimensions of the problem.
   m <- length(CN) + 1
   n <- length(x)

   # Transition parameters.
   Q <- matrix(0.1/(m-1), ncol = m, nrow = m);
   diag(Q) <- 0.9;

   # Emission parameters (initial values).
   sz <- median(x, na.rm=TRUE) / 2
   mu <- sz * c(CN, m)
   s2 <- var(x, na.rm=TRUE)
   nu <- 6

   # Initial probabilities.
   initial.prob <- rep(1/m, m)



###############################################
#            VARIABLE DEFINITIONS             #
###############################################


    old.loglik <- -Inf;
    iter <- 0;

    # Updated in E.step (Q as well).
    emission.prob <- matrix(NA_real_, nrow = n, ncol = m)
    phi <- matrix(NA_real_, nrow = n, ncol = m)
    phi.weights <- matrix(NA_real_, nrow = n, ncol = m)
    loglik <- NA_real_


###############################################
#           FUNCTION DEFINITIONS              #
###############################################


   # Perform the E-step of the modified Baum-Welch algorithm.
   # Note: the double arrow assignment <<- looks for the 
   # variable on the enclosing scopes. The design prevents
   # passing parameters by value for objects that can be
   # very large, such as 'emission.prob' and 'phi.weights'.
   # The M-step can yield invalid variance estimates because
   # of numeric underflow (see update.S() for example). In that
   # case no update is performed for the given state and the
   # algorithm 'waits' for valid parameter values.

   E.step <- function () {

      weights <- matrix(1, nrow = n, ncol = m)

      # Emission probabilities are computed up to constant terms.
      # The following piece of code is stolen from the function
      # mahalanobis.
      for (i in 1:m) {
         emission.prob[,i] <<- dt((x-mu[i])/sqrt(s2), df=nu)

         # Use Mahalanobis distances to compute the weights.
         if (!is.infinite(nu)) {
            mahalanobis <- sum(abs(x-mu[i]), na.rm=TRUE)/sqrt(s2)
            weights[,i] <- (nu+1) / (nu + mahalanobis)
         }
      }

      # Emission probabilities of NAs are set to 1 for every states.
      emission.prob[is.na(emission.prob)] <<- 1;

      forwardback <- .Fortran("fwdb",
          as.integer(m),
          as.integer(n),
          initial.prob,
          emission.prob,
          Q,
          double(n),
          matrix(double(n * m), nrow = n),
          matrix(double(m^2), nrow = m),
          double(1),
          PACKAGE = "dryhic")

      loglik <<- forwardback[[9]]
      phi <<- forwardback[[7]]
      transitions <- matrix(forwardback[[8]], nrow=m)

      # Prevent underflow.
      transitions[transitions < MCHNEPS] <- MCHNEPS

      Q <<- transitions / rowSums(transitions)
      phi.weights <<- phi * weights

  }

  # Update the means. Here we maintain the constraints on the
  # lower states, but not on the higest that is open-ended.
  update.mu <- function() {
     sumPhi.weights.x <- colSums(phi.weights*x, na.rm = TRUE)
     sumPhi.weights <- colSums(phi.weights, na.rm = TRUE)
     sz <- sum(CN*sumPhi.weights.x[-m]) / sum(CN^2*sumPhi.weights[-m])
     opnd <- sumPhi.weights.x[m] / sumPhi.weights[m]
     mu <- sz * CN
     mu[m] <- max(sz*m, opnd)
     return (mu)
  }

  # Update the variances-covariances.
  update.s2 <- function() {
     x.minus.mu <- sweep(matrix(x, nrow=n, ncol=m), MARGIN=2, STATS=mu)
     sumPhi.weights.x.minus.mu.2 <-
        colSums(phi.weights*x.minus.mu^2, na.rm = TRUE)
     sumPhi.weights <- colSums(phi.weights, na.rm = TRUE)
     s2 <- sum(sumPhi.weights.x.minus.mu.2) / sum(sumPhi.weights)
     return (s2)
  }

  update.nu <- function() {
     weights = rowSums(phi.weights)
     RHS <- 1 + mean((log(weights) - weights), na.rm = TRUE)
     # Find an upper bound.
     no.upper.bound <- TRUE
     new.nu <- nu
     lower.bound <- 0
     while (no.upper.bound) {
         if (digamma(new.nu/2) - digamma((new.nu + 1)/2) -
             log(new.nu/2) + log((new.nu + 1)/2) < RHS) {
                 lower.bound <- new.nu
                 new.nu <- 2 * new.nu
         }
         else {
             no.upper.bound <- FALSE
             upper.bound <- new.nu
         }
     }
     # The degree of freedom new.nu is estimated with a
     # precision of 0.01
     while (upper.bound - lower.bound > 0.01) {
         new.nu <- (upper.bound + lower.bound) / 2
         if (digamma(new.nu/2) - digamma((new.nu + 1)/2) -
             log(new.nu/2) + log((new.nu + 1)/2) < RHS)
                 lower.bound <- new.nu
         else {
             upper.bound <- new.nu
         }
     }

     # Note: rounding prevents oscillation of the estimates.
     return (round(new.nu, 1))
  }



###############################################
#                 MAIN LOOP                   #
###############################################

    for (iter in 1:MAXITER) {

       #browser()
       cat(paste("\riteration:", iter))

       # Update Q, phi.weights and emission.prob.
       E.step()

       # CM-step 1: updating mu and S.
       sumPhi <- colSums(phi, na.rm = TRUE);
       sumPhi.weights <- colSums(phi.weights, na.rm = TRUE);

       mu <- update.mu();
       s2 <- update.s2();

       if (!is.infinite(nu)) {

           # Update Q, phi.weights and emission.prob.
           E.step()

           # CM-step 2: updating nu.
           nu <- update.nu()
           if (nu > 250) {
               nu <- Inf;
           }

       } # if (!is.infinite(nu))


       if (abs(loglik - old.loglik) < TOLERANCE)
           break

       old.loglik <- loglik

    } # for (iter in 1:maxiter)

    cat("\n")

    ViterbiPath <- Viterbi(Q, initial.prob, emission.prob)

    # Get the copy number and correct the open-ended state.
    cn <- ViterbiPath
    cn[ViterbiPath == m] <-  mu[m] / mu[1]

    return(list(mu[1], cn))

}
