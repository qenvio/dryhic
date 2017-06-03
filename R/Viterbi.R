Viterbi <- function (Q, initialProb, emissionProb) {

   n <- nrow(emissionProb)
   m <- as.integer(nrow(Q))
   ViterbiPath <- integer(n)
   series.length <- n

   emissionProb <- log(emissionProb)
   Q <- log(Q)

   # Prevent underflow.
   emissionProb[is.infinite(emissionProb)] <- -325
   Q[is.infinite(Q)] <- -325

   viterbi <- .Fortran("vit",
      n,
      m,
      log(initialProb),
      emissionProb,
      Q,
      integer(n),
      matrix(double(n * m), nrow = n),
      PACKAGE = "dryhic"
   )

   return(viterbi[[6]])

}
