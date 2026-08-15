#' Inference with REMLF90
#' 
#' Fits a Linear Mixed Model by Restricted Maximum Likelihood
#' 
#' @param fixed an object of class \link{formula} (or one that can be coerced to
#'   that class): a symbolic description of the fixed effects of the model to be
#'   fitted. The details of model specification are given under 'Details'.
#' @param random if not \code{NULL}, an object of class \link{formula} with the 
#'   unstructured random effects.
#' @param genetic if not \code{NULL}, a list with relevant parameters for an 
#'   additive genetic effect; see 'Details'.
#' @param spatial if not \code{NULL}, a list with relevant parameters for a 
#'   spatial random effect; see 'Details'.
#' @param generic if not \code{NULL}, a named list with an \code{incidence}
#'   matrix and either a \code{covariance} or a \code{precision} matrix; see
#'   'Details'.
#' @param genomic if not \code{NULL}, a list with genomic selection parameters
#'   for single-step GBLUP. Requires a \code{genetic} component. Elements:
#'   \describe{
#'     \item{snp_file}{(required) path to genotype file in PREGSF90 format
#'       (animal ID + fixed-width genotype string with 0/1/2/5)}
#'     \item{map_file}{path to SNP map file (header: SNP_ID, CHR, POS)}
#'     \item{whichG}{G matrix method: 1 = VanRaden 2008 (default), 2 = Amin
#'       2007, 3 = Yang 2010}
#'     \item{tunedG}{G scaling: 0 = none, 1 = diag=1, 2 = match A22 (default),
#'       3 = mean(G)=mean(A22), 4 = Fst}
#'     \item{AlphaBeta}{numeric(2). Blending weights for G and A22
#'       (default c(0.95, 0.05))}
#'     \item{minfreq}{MAF threshold (default 0.05)}
#'     \item{callrate}{SNP call rate threshold (default 0.90)}
#'     \item{callrateAnim}{animal call rate threshold (default 0.90)}
#'     \item{verify_parentage}{0-3. Mendelian conflict handling (default 3)}
#'     \item{saveG}{logical. Save G matrix to file (default FALSE)}
#'     \item{saveA22}{logical. Save A22 matrix to file (default FALSE)}
#'     \item{save_ginverse}{logical. Save the genomic G-inverse (default
#'       FALSE). Required to run a later GWAS with \code{\link{postgsf90}};
#'       set \code{TRUE} at fit time to make the model GWAS-ready}
#'     \item{extra_options}{character vector of additional PREGSF90 OPTION
#'       lines passed through directly}
#'   }
#' @param data a data frame with variables and observations
#' @param var.ini if not \code{NULL}, a named list with one item for each random
#'   component in the model. See 'Details'.
#' @param method either 'ai' or 'em' for Average-Information or 
#'   Expectation-Maximization REML respectively
#' @param breedR.bin character. The local directory where the package binaries 
#'   are stored, or any of 'remote' or 'submit' for remote computing. See 
#'   'Details'.
#' @param progsf90.options character. Passed directly to OPTIONS field in 
#'   PROGSF90. See available options for 
#'   \href{http://nce.ads.uga.edu/wiki/doku.php?id=readme.reml#options}{REMLF90}
#'   and for 
#'   \href{http://nce.ads.uga.edu/wiki/doku.php?id=readme.aireml#options}{AIREMLF90}.
#'    Option \code{sol se} is passed always and cannot be removed. No checks are
#'   performed, handle with care.
#' @param weights numeric. A vector of weights for the residual variance.
#' @param parallel logical or integer. If \code{TRUE}, the AR rho grid search
#'   runs in parallel using all available cores. If an integer, uses that many
#'   cores. Requires the \code{parallel} package (default FALSE).
#' @param progress_file character. Path of a file where the backend's output is
#'   written as it is produced, so that a long fit can be followed with
#'   \code{tail -F} or inspected with \code{\link{reml_checkpoint}} from another
#'   R session. A relative path is resolved against the working directory at the
#'   time of the call. The backend's error stream goes to a companion file with
#'   the suffix \code{.err}. Local fits only. Default \code{NULL}, which keeps
#'   the output silent until the fit returns.
#' @param cont logical. Resume a previous fit, taking its initial variances from
#'   the last usable round of \code{progress_file}. Cannot be combined with an
#'   explicit \code{var.ini}. The previous log is renamed with a timestamp
#'   before the new one is written. Default \code{FALSE}. Named after the
#'   equivalent argument of \code{\link{gibbsf90}}, though the mechanism
#'   differs: for Gibbs it is a backend \code{OPTION}, whereas AIREMLF90 has no
#'   such option and the log itself is the checkpoint.
#' @param debug logical. If \code{TRUE}, the input files for blupf90 programs
#'   and their output are shown, but results are not parsed.
#'
#' @details If either \code{genetic} or \code{spatial} are not \code{NULL}, the 
#'   model residuals are assumed to have an additive genetic effects or a 
#'   spatially structured random effect, respectively. In those cases, 
#'   \code{genetic} and \code{spatial} must be lists with named relevant 
#'   parameters.
#'   
#'   The \code{generic} component implements random effects with custom 
#'   incidence and covariance (or precision) matrices. There can be any number 
#'   of them, stored in a named list with custom unique names that will be used 
#'   to identify and label the results. Each effect in the list must have an 
#'   \code{incidence} argument and either a \code{covariance} or a 
#'   \code{precision} matrix with conforming and suitable dimensions. 
#'   Optionally, an initial variance for the REML algorithm can be specified in 
#'   a third argument \code{var.ini}.
#'   
#'   \subsection{Genetic effect}{ The available models for the genetic effect 
#'   are \code{add_animal} and \code{competition}. \code{add_animal} stands for 
#'   an additive animal model with a given pedigree. \code{competition} includes
#'   the \emph{direct} additive genetic effect and also a \emph{competition} 
#'   additive genetic effect possibly correlated with the former.
#'   
#'   The minimum elements in the list of the \code{genetic} component are: 
#'   \itemize{\item \code{model} a string, either \code{add_animal} or 
#'   \code{competition} \item \code{pedigree} either an object of class 
#'   \code{pedigree} (see \code{\link{build_pedigree}}) or a \emph{data.frame} 
#'   with exactly three columns identifying the individual, his father and his 
#'   mother respectively \item \code{id} either a vector of codes or the name of
#'   the variable with the individual identifier in the \code{data}}
#'   
#'   Optional common components are:
#'   
#'   \itemize{ \item \code{var.ini} a positive initial value for the variance 
#'   component(s). For a \code{competition} model, the same initial value is 
#'   used for both effects, with a negative initial correlation between them of 
#'   half this value.}
#'   
#'   Finally, for model \code{competition} there are further mandatory and 
#'   optional elements:
#'   
#'   \itemize{ \item mandatory elements \itemize{ \item \code{coord} a matrix, 
#'   list or data.frame with two columns for the rows and columns of the 
#'   observations, respectively. This element is necessary even if duplicated in
#'   a \code{spatial} component. \item \code{competition_decay} a positive 
#'   number. The intensity of competition is weighted by the distance according 
#'   to \eqn{1/d^\alpha}. This element specifies the exponent \eqn{\alpha} to be
#'   used. Typically \eqn{1} or \eqn{2}. }
#'   
#'   \item optional elements \itemize{ \item \code{pec} Permanent Environmental 
#'   Competition effect. If present, this must be a named list with elements 
#'   \code{present} which is either \code{TRUE} or \code{FALSE} and (optionally) 
#'   \code{var.ini} specifying the initial variance for this effect. } }
#'   
#'   The Permanent Environmental Competition (\code{pec}) effect is actually 
#'   non-genetic in nature. However, it was included as an option to the 
#'   (genetic) competition effect as it is usually used in conjunction with it. 
#'   }
#'   
#'   
#'   \subsection{Spatial effects}{ The available models for the spatial effect 
#'   are \code{splines} and \code{AR1}. \code{splines} uses a  two-dimensional 
#'   tensor product of B-splines to represent the smooth spatially structured 
#'   effect (Cappa and Cantet, 2007). \code{AR1} uses a kronecker product of 
#'   autoregressive models for the rows and columns (Dutkowski et al., 2002).
#'   
#'   In both cases, the minimum necessary components in the list are \itemize{ 
#'   \item \code{model} a string, either \code{splines} or \code{AR} \item 
#'   \code{coord} a matrix, list or data.frame with two columns for the rows and
#'   columns respectively. }
#'   
#'   Optional common components are \itemize{ \item \code{var.ini} a positive 
#'   initial value for the variance component }
#'   
#'   Finally, optional model-dependent components are \itemize{ \item For model 
#'   \code{splines} \itemize{ \item \code{n.knots} a vector of two integers with
#'   the number of \emph{internal} knots for the rows and columns of the spline 
#'   design } \item For model \code{AR} \itemize{ \item \code{rho} a vector of 
#'   two numbers strictly between -1 and 1 with the autoregressive parameters 
#'   for the rows and columns respectively. Alternatively, a matrix or 
#'   data.frame with two columns where every row contain a combination of 
#'   autoregressive parameters to be tried. } }
#'   
#'   The \emph{internal} knots cover the region with observations at regular 
#'   intervals (in each dimension). For the splines design, three additional 
#'   knots are automatically added before the first internal knot, and other 
#'   tree after the last one, in each dimension. As a result, if \code{n.knots =
#'   c(n1, n2)}, then the final number of parameters of the splines model is 
#'   \eqn{(n1 + 2)(n2 + 2)}.
#'   
#'   If \code{n.knots} is omitted, a sensible number of knots for each dimension
#'   is computed based on the number of observations in the experiment. See the 
#'   internal function \code{\link{determine.n.knots}}. This is the default 
#'   function that computes the default number of knots, but you can provide an 
#'   alternative default function through \code{\link{breedR.setOption}}.
#'   
#'   Due to limitations of the REML backend, we can only fit models with 
#'   \emph{fixed} autoregressive parameters. \code{remlf90()} will fit as many 
#'   models as rows in \code{rho}, and return the results of the most likely. It
#'   will also return the list of log-likelihoods for each tried combination of 
#'   autoregressive parameters, in the component \code{rho} of the \code{reml} 
#'   object. This is useful for visualization, and further refinement of the 
#'   search for appropriate parameters.
#'   
#'   If any of the values in either column of \code{rho} is \code{NA}, then a 
#'   default set of values for the corresponding dimension will be set. See 
#'   \code{breedR.getOption('ar.eval')}, for the current defaults. You can set 
#'   your own defaults with \code{\link{breedR.setOption}}. Each will be 
#'   combined with every other value in the other column.
#'   
#'   Omitting the specification of \code{rho} is equivalent to \code{rho = c(NA,
#'   NA)}. }
#'   
#'   \subsection{Intercept}{ An intercept is automatically introduced in the 
#'   model provided the user doesn't explicitly prevents it by using \code{0} or
#'   \code{-1} in the \code{fixed} formula (as conventional in \code{R}), 
#'   \emph{and} there are no other categorical covariates in \code{fixed}. The 
#'   latter condition is actually a limitation of (ai)remlf90 backends, which 
#'   would in any case return an estimate for each level of the categorical 
#'   covariates while returning 0 for the intercept. It does not allow 
#'   alternative parameterizations. }
#'   
#'   
#'   \subsection{Initial variances}{ Initial variance components can also be 
#'   specified through an additional argument \code{var.ini}. You can either use
#'   default initial values for the variance components (see 
#'   \code{?breedR.options}) or specify custom values for \emph{each} and 
#'   \emph{all} variance components in the model. In this case, \code{var.ini} 
#'   must be a named list with one element for each term in \code{random} with 
#'   matching names, plus one last element named \code{residual} for the initial
#'   residual variance. Furthermore if there are \code{genetic} or 
#'   \code{spatial} effects, they must as well include a numeric element 
#'   \code{var.ini} with the initial variance component specification for the 
#'   corresponding effect. }
#'   
#'   
#'   \subsection{Inference method}{ AI-REML is usually faster than EM-REML, and 
#'   it provides more results. Namely, standard errors of the variance 
#'   components estimates, and covariances as well. On the other hand, is less 
#'   robust than EM-REML and it usually gives extreme results when used with the
#'   splines spatial model (as in \code{spatial = list(model ='splines', ...)}).
#'   
#'   Even when an effect accounts for no variance at all, EM-REML will always 
#'   estimate a positive variance which will be determined by the starting 
#'   value. If AI-REML does not converge but EM-REML does with the same dataset 
#'   and model, re-run EM-REML with a small starting value for the effect. If 
#'   the estimate does not change, it is likely that there is no variance. }
#'   
#'   \subsection{Heterogeneous residual variances}{
#'   BLUPF90+ supports heterogeneous residual variances modeled as a
#'   log-linear function of covariates: var(e) = exp(a0 + a1*X1 + a2*X2).
#'   Use the \code{\link{hetres_options}} helper to generate the correct
#'   option strings. For example, to model residual variance as a function
#'   of spatial coordinates x (col 4) and y (col 5):
#'   \code{progsf90.options = hetres_options(covariate_cols = c(4, 5),
#'   initial = c(log(10), 0.01, 0.01))}.
#'   See \code{\link{hetres_options}} for details.}
#'
#'   \subsection{Remote computing}{ If \code{breedR.bin = 'remote'}, the REML 
#'   program will be run remotely and the results will be automatically 
#'   transferred back automatically. While if \code{breedR.bin = 'submit'} the 
#'   job will be submitted to the server, and the job-id and other relevant 
#'   information about the model will be returned instantly. The returned object
#'   can be used to retrieve the results or check the status of the job. Several
#'   jobs can be submitted in parallel. See \code{?remote} to learn how to
#'   configure breedR for remote computing, and how to manage submitted jobs.}
#'
#'   \subsection{Following and resuming a long fit}{ A REML fit returns
#'   everything or nothing: on a large model it can run for days, and a crash
#'   at the last round leaves no partial result and nothing to restart from.
#'   \code{progress_file} and \code{cont} address both halves of that.
#'
#'   The backend prints the residual and genetic (co)variance matrices at
#'   \emph{every} round, and those matrices are the whole state of the
#'   optimiser. A streamed log is therefore a complete per-round checkpoint,
#'   and resuming is simply a matter of reading the last usable round back into
#'   the initial variances -- which typically converges in a couple of rounds.
#'   The recovered components are distributed to the right places
#'   automatically, so unlike \code{\link{reml_checkpoint}} there is nothing to
#'   assemble by hand.
#'
#'   Because \code{remlf90} blocks while the backend runs, watching a fit in
#'   progress means reading its log from a \emph{second} R session, with
#'   \code{tail -F} or \code{\link{reml_checkpoint}}. Note \code{-F} rather than
#'   \code{-f}: \code{cont = TRUE} renames the previous log with a timestamp
#'   before starting, so a follower attached by name would otherwise be left on
#'   the old file. Those backups are never cleaned up.
#'
#'   Interrupting a streamed fit returns control immediately but leaves the
#'   backend running to completion in the background; breedR does not wait for
#'   it, because closing the connection would block until the fit finished. The
#'   orphaned process keeps files open in \code{tempdir()}, which can disturb a
#'   further fit in the same R session, so start a fresh session after
#'   interrupting one.
#'
#'   Resuming may start slightly earlier than the last round printed. The final
#'   round of a log is deliberately never used, because a log truncated by a
#'   crash ends mid-line and a number cut inside its exponent parses cleanly
#'   while being wrong by orders of magnitude; and a round whose matrices are
#'   not positive definite is skipped in favour of an earlier one. The round
#'   actually used is reported and stored in \code{res$reml$resumed_from}.
#'
#'   The most valuable case is not a crash but a fit that exhausted its
#'   iterations: that returns \code{NA} variance components with a warning,
#'   while the log holds every round intact, so \code{cont = TRUE} picks up
#'   exactly where an otherwise useless result object left off.
#'
#'   Both arguments are for local fits. They are rejected for
#'   \code{breedR.bin = 'remote'} or \code{'submit'}, and for an AR \code{rho}
#'   grid search, which fits one model per \code{rho}. With \code{genomic},
#'   they cover the REML phase only: PREGSF90 runs first and is not streamed,
#'   nor is it skipped when resuming.}
#'
#' @return An object of class 'remlf90' that can be further questioned by
#'   \code{\link{fixef}}, \code{\link{ranef}}, \code{\link{fitted}}, etc.
#'
#'   Among its components, \code{funvars} holds the functions of (co)variances
#'   requested through \code{OPTION se_covar_function} (see
#'   \code{\link{var_functions}}), and is a matrix with one column per requested
#'   function and three rows:
#'
#'   \describe{
#'   \item{\code{mean}}{the plug-in value of the function evaluated at the REML
#'   solution. This is the estimate, and is consistent with the variance
#'   components reported in the same object.}
#'   \item{\code{sample mean}}{the mean of the draws AIREMLF90 takes from the
#'   asymptotic distribution of the variance components, evaluating the function
#'   at each draw.}
#'   \item{\code{sample sd}}{the standard deviation of those draws, used as the
#'   standard error of the estimate.}
#'   }
#'
#'   \code{summary()} prints all three. A gap between the first two means the
#'   sampling distribution of the function is skewed, either by the
#'   non-linearity of the function itself or by the boundary of the parameter
#'   space, and therefore that the variance components are poorly determined.
#'   When no variance functions are requested (e.g. \code{method = 'em'}, or a
#'   model with no genetic effect), \code{funvars} is an empty list.
#'
#'   When \code{progress_file} is given, \code{res$reml$progress_file} records
#'   the path used; when \code{cont = TRUE}, \code{res$reml$resumed_from}
#'   records the round the initial variances were taken from.
#' @seealso \code{\link{reml_checkpoint}} to read the variance components of a
#'   log, including one still being written.
#'
#'   \code{\link[pedigreemm]{pedigree}}
#' @references progsf90 wiki page: \url{http://nce.ads.uga.edu/wiki/doku.php}
#'   
#'   E. P. Cappa and R. J. C. Cantet (2007). Bayesian estimation of a surface to
#'   account for a spatial trend using penalized splines in an individual-tree 
#'   mixed model. \href{http://dx.doi.org/10.1139/x07-116}{\emph{Canadian 
#'   Journal of Forest Research} \strong{37}(12):2677-2688}.
#'   
#'   G. W. Dutkowski, J. Costa e Silva, A. R. Gilmour, G. A. López (2002). 
#'   Spatial analysis methods for forest genetic trials. 
#'   \href{http://dx.doi.org/10.1139/x02-111}{\emph{Canadian Journal of Forest 
#'   Research} \strong{32}(12):2201-2214}.
#'   
#' @examples
#' \dontrun{
#' ## These examples fit models with the BLUPF90+ backend, which is
#' ## downloaded at install time and unavailable during R CMD check.
#' ## Linear model
#' n <- 1e3
#' dat <- transform(data.frame(x = runif(n)),
#'                  y = 1 + 2*x + rnorm(n, sd = sqrt(3)))
#' res.lm <- remlf90(fixed = y ~ x, data = dat)
#' summary(res.lm)
#' 
#' ## Linear Mixed model
#' f3 = factor(sample(letters[1:3], n, replace = TRUE))
#' dat <- transform(dat,
#'                  f3 = f3,
#'                  y = y + (-1:1)[f3])
#' res.lmm <- remlf90(fixed  = y ~ x,
#'                    random = ~ f3,
#'                    data   = dat)
#' 
#' ## Following a long fit, and picking it up again
#' ## While this runs, `tail -F run.log` in a shell, or from a second R
#' ## session: reml_checkpoint("run.log")
#' res.slow <- remlf90(fixed  = y ~ x,
#'                     random = ~ f3,
#'                     data   = dat,
#'                     progress_file = "run.log")
#'
#' ## Had it crashed, or exhausted its iterations, resume from that log:
#' res.again <- remlf90(fixed  = y ~ x,
#'                      random = ~ f3,
#'                      data   = dat,
#'                      progress_file = "run.log",
#'                      cont   = TRUE)
#' res.again$reml$resumed_from
#'
#' ## Generic model (used to manually fit the previous model)
#' inc.mat <- model.matrix(~ 0 + f3, dat)
#' cov.mat <- diag(3)
#' res.lmm2 <- remlf90(fixed  = y ~ x,
#'                     generic = list(f3 = list(inc.mat,
#'                                              cov.mat)),
#'                     data   = dat)
#' all.equal(res.lmm, res.lmm2, check.attributes = FALSE)  # TRUE
#'                
#' ## Animal model
#' ped <- build_pedigree(c('self', 'dad', 'mum'),
#'                       data = as.data.frame(m1))
#' res.am <- remlf90(fixed   = phe_X ~ sex,
#'                   genetic = list(model    = 'add_animal',
#'                                  pedigree = ped,
#'                                  id       = 'self'),
#'                   data    = as.data.frame(m1))
#' }
#'
#' \dontrun{
#' ## Same model with specification of initial variances
#' res.am <- remlf90(fixed   = phe_X ~ sex,
#'                   genetic = list(model    = 'add_animal',
#'                                  pedigree = ped,
#'                                  id       = 'self',
#'                                  var.ini  = 1),
#'                   data    = as.data.frame(m1),
#'                   var.ini = list(resid = 1))
#'  
#' ## Animal-spatial models
#' gen.globulus <- list(model    = 'add_animal',
#'                      pedigree = globulus[, 1:3],
#'                      id       = 'self')
#' res.bm <- remlf90(fixed   = phe_X ~ gg,
#'                   genetic = gen.globulus,
#'                   spatial = list(model = 'blocks', 
#'                                  coord = globulus[, c('x','y')],
#'                                  id    = 'bl'),
#'                   data    = globulus)
#'                   
#' res.am <- remlf90(fixed   = phe_X ~ gg,
#'                   genetic = gen.globulus,
#'                   spatial = list(model = 'AR', 
#'                                  coord = globulus[, c('x','y')],
#'                                  rho   = c(.85, .8)),
#'                   data    = globulus)
#' 
#' res.sm <- remlf90(fixed   = phe_X ~ gg,
#'                   genetic = gen.globulus,
#'                   spatial = list(model = 'splines', 
#'                                  coord = globulus[, c('x','y')],
#'                                  n.knots = c(5, 5)),
#'                   data    = globulus,
#'                   method  = 'em')   # Necessary for splines models!!!
#' 
#' 
#' ## Competition models
#' 
#' # This may take some minutes...
#' # and need to be fitted with 'em'
#' res.cm <- remlf90(fixed   = phe_X ~ 1,
#'                  genetic = list(model = 'competition',
#'                                 pedigree = globulus[, 1:3],
#'                                 id = 'self',
#'                                 coord = globulus[, c('x','y')],
#'                                 competition_decay = 1,
#'                                 pec = list(present = TRUE)),
#'                  method = 'em',
#'                  data = globulus)
#' }
#' 
#' @export
#' @importFrom stats terms model.response logLik runif
remlf90 <- function(fixed,
                    random = NULL,
                    genetic = NULL,
                    spatial = NULL,
                    generic = NULL,
                    genomic = NULL,
                    data,
                    var.ini = NULL,
                    method = c('ai', 'em'),
                    breedR.bin = breedR.getOption("breedR.bin"),
                    progsf90.options = NULL,
                    weights = NULL,
                    parallel = FALSE,
                    progress_file = NULL,
                    cont = FALSE,
                    debug = FALSE) {
  
  ## Assumptions:
  ## Only 1 pedigree
  ## Solution file header: trait effect level solution
  ## No intercept
  ## (not generalized) Linear Mixed Model

  ## TODO: 
  # Allow for multiple responses
  # Allow for generalized mixed models


  ###  Call
  mc <- mcout <- match.call()

  # Builds model frame by joining the fixed and random terms
  # and translating the intercept (if appropriate) to a fake covariate
  # Add an additional 'term.types' attribute within 'terms'
  # indicating whether the term is 'fixed' or 'random'
  # progsf90 does not allow for custom model parameterizations
  # and they don't use intercepts
  mf <- build.mf(mc)
  mt <- attr(mf, 'terms')
  
  ### Checks
  if ( missing(fixed) | missing(data) ) { 
    stop("Usage: remlf90(fixed, data, ...); see ?remlf90\n")
  }
  if ( !inherits(fixed, "formula") ) { 
    stop("'fixed' should be a formula\n") 
  }
  if ( !inherits(data, "data.frame") ) { 
    stop("'data' should be a data.frame\n") 
  }
  if ( attr(terms(fixed), 'response') != 1L ) {
    stop("There is no response in the 'fixed' argument\n")
  }
  if ( !is.null(random) ) {
    check.random = TRUE
    if ( !inherits(random, "formula") )
      check.random = FALSE
    else if( attr(terms(random), 'response') != 0L )
      check.random = FALSE
    if( !check.random )
      stop("random should be a response-less formula\n")
  }
  ## Checked before the binaries so that these guards are reachable without
  ## the backend installed.
  progress_file <- check_progress_args(progress_file, cont, breedR.bin)

  if (!check_progsf90(quiet = debug | !interactive())) {
    stop('Binary dependencies missing. See ?install_progsf90')
  }

  ### Parse arguments
  method <- tolower(method)
  method <- match.arg(method)
  

  ### Number of traits
  # ntraits <- ncol(as.matrix(model.response(mf)))
  
  ### Response matrix (ntraits = ncols)
  responsem <- as.matrix(model.response(mf))
  
  ## Genetic specification
  if (!is.null(genetic)) {
    ## TODO: Ideally, I should pass the model frame only, containing the
    ## necessary variables (also for special effects and response)
    genetic <- do.call('check_genetic',
                       c(genetic, list(data = data, 
                                       response = responsem)))
  }
  
  
  ## Spatial specification
  if (!is.null(spatial)) {
    ## TODO: Ideally, I should pass the model frame only, containing the
    ## necessary variables (also for special effects and response)
    spatial <- do.call('check_spatial', 
                       c(spatial, list(data = data,
                                       response = responsem)))
    
    # If AR model without rho specified
    # we need to fit it with several fixed rho's
    # and return the most likely
    if( spatial$model == 'AR' ) {

      if (!is.null(nrow(spatial$rho))) {
        ## grid case
        ## Each rho combination is independent — fit with tryCatch so one
        ## failure doesn't kill the entire search.

        ## A grid is N fits, not one: a single log would be overwritten once
        ## per rho, and a single checkpoint would seed every rho alike. Keyed
        ## on the grid rather than on `parallel`, which is inert outside here.
        if (!is.null(progress_file) || isTRUE(cont))
          stop("'progress_file' and 'cont' are not supported for an AR rho ",
               "grid search, which fits one model per rho.\n",
               " Fix rho to a single pair and re-run.", call. = FALSE)

        n_rho <- nrow(spatial$rho)

        # Determine number of cores for parallel execution
        n_cores <- 1L
        if (!identical(parallel, FALSE)) {
          if (isTRUE(parallel)) {
            n_cores <- max(1L,
              parallel::detectCores(logical = FALSE) - 1L)
          } else if (is.numeric(parallel) && parallel >= 2) {
            n_cores <- as.integer(parallel)
          }
          n_cores <- min(n_cores, n_rho)
        }

        if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
          ## Parallel grid search (works on all platforms including Windows).
          ##
          ## Strategy:
          ##   1. Run each rho through the full R pipeline sequentially
          ##      (including binary) — but copy the generated files to
          ##      separate directories before the next rho overwrites them
          ##   2. Run BLUPF90+ binaries in parallel via parLapply
          ##      (workers only call system2(), no R packages needed)
          ##   3. Extract log-likelihoods from parallel binary output
          ##   4. Full fit of winning rho for complete R result
          ##
          ## The sequential file generation (Step 1) is dominated by the
          ## binary execution time. However, Step 1 also produces the
          ## correct files. We skip Step 1's binary results and re-run
          ## only the binaries in parallel in Step 2.
          ##
          ## Optimization: run Step 1 with each rho to get the files,
          ## but actually we only need the PARAMETER files, not the
          ## solutions. So we use eval(mc) which runs the binary too.
          ## Net cost: 16 sequential binary runs + 16/n_cores parallel
          ## binary runs + 1 final. This is WORSE unless we can skip
          ## the binary in Step 1.
          ##
          ## Better optimization: since all rhos share the same data,
          ## pedigree, and fixed effects, run the first rho fully to get
          ## the base files, then for subsequent rhos only regenerate
          ## the precision matrix file (which is the only thing that
          ## changes with rho).
          message("AR grid search: ", n_rho, " combinations on ",
                  n_cores, " cores")

          caller_env <- parent.frame()
          blupf90_bin <- file.path(breedR.bin,
                                    progsf90_files(breedR.os.type()))

          # Run the first rho fully to establish the base model files
          mc_first <- mc
          mc_first$spatial$rho <- as.numeric(spatial$rho[1, ])
          mc_first$parallel <- FALSE
          first_result <- tryCatch(
            suppressMessages(eval(mc_first, envir = caller_env)),
            error = function(e) NULL)

          # Copy base files from tempdir to each rho-specific directory
          base_files <- list.files(tempdir(), full.names = TRUE)
          base_files <- base_files[!file.info(base_files)$isdir]

          rho_dirs <- character(n_rho)
          rho_valid <- logical(n_rho)

          for (i in seq_len(n_rho)) {
            rd <- file.path(tempdir(), paste0("rho_", i))
            dir.create(rd, showWarnings = FALSE, recursive = TRUE)
            rho_dirs[i] <- rd

            # Copy all base files
            file.copy(base_files, rd, overwrite = TRUE)

            # Now regenerate ONLY the AR precision matrix for this rho
            rho_i <- as.numeric(spatial$rho[i, ])
            tryCatch({
              ar_obj <- breedr_ar(
                coordinates = spatial$coordinates,
                rho = rho_i,
                autofill = spatial$autofill)
              # Write the updated precision matrix file
              ar_pf90 <- renderpf90(ar_obj)
              sm <- as.triplet(vcov(ar_obj))
              utils::write.table(sm, file = file.path(rd, "ar_spatial"),
                row.names = FALSE, col.names = FALSE, na = "0")
              rho_valid[i] <- TRUE
            }, error = function(e) {
              rho_valid[i] <<- FALSE
            })
          }

          n_valid <- sum(rho_valid)
          if (n_valid == 0)
            stop("Could not generate files for any rho.", call. = FALSE)

          # Run BLUPF90+ binaries in parallel.
          # Workers only call system2() — no R packages needed.
          cl <- parallel::makeCluster(min(n_cores, n_valid))
          # Guarantee cluster cleanup even if parLapply fails
          on.exit(parallel::stopCluster(cl), add = TRUE)

          binary_results <- tryCatch(
            parallel::parLapply(
              cl, rho_dirs[rho_valid], function(rho_dir) {
                bin_name <- basename(blupf90_bin)
                bin_local <- file.path(rho_dir, bin_name)
                file.copy(blupf90_bin, bin_local, overwrite = TRUE)
                cwd <- setwd(rho_dir)
                on.exit({
                  setwd(cwd)
                  unlink(bin_local)
                })
                system2(file.path(".", bin_name), input = "parameters",
                        stdout = TRUE, stderr = TRUE)
              }),
            error = function(e) {
              warning("Parallel binary execution failed: ",
                      conditionMessage(e),
                      ". Falling back to sequential.", call. = FALSE)
              NULL
            }
          )

          parallel::stopCluster(cl)
          on.exit(NULL)  # remove the on.exit since we stopped manually

          # If parallel failed, fall back to sequential
          if (is.null(binary_results)) {
            caller_env <- parent.frame()
            eval.rho.fallback <- function(rho, mc, envir) {
              mc$spatial$rho <- rho
              mc$parallel <- FALSE
              tryCatch(suppressMessages(eval(mc, envir = envir)),
                       error = function(e) NULL)
            }
            ans.rho <- lapply(seq_len(n_rho), function(i) {
              eval.rho.fallback(as.numeric(spatial$rho[i, ]),
                                mc, envir = caller_env)
            })
            succeeded <- !vapply(ans.rho, is.null, logical(1))
            if (!any(succeeded))
              stop("All rho combinations failed.", call. = FALSE)
            loglik.rho <- transform(spatial$rho,
              loglik = suppressWarnings(vapply(ans.rho, function(x) {
                if (is.null(x)) NA_real_ else as.numeric(logLik(x))
              }, numeric(1))))
            rho.idx <- which.max(loglik.rho$loglik)
            ans <- ans.rho[[rho.idx]]
            unlink(rho_dirs, recursive = TRUE)
            ans$rho <- loglik.rho
            return(ans)
          }

          # Extract log-likelihoods from binary output.
          # BLUPF90+ format: "-2logL =     5619.45354541821     : AIC = ..."
          loglik_vals <- rep(NA_real_, n_rho)
          j <- 1
          for (i in seq_len(n_rho)) {
            if (!rho_valid[i]) next
            bout <- binary_results[[j]]; j <- j + 1
            if (!is.null(attr(bout, 'status'))) next
            ll_line <- grep("^-2logL", bout, value = TRUE)
            if (length(ll_line) > 0) {
              ll_str <- tail(ll_line, 1)
              # Extract the number after "-2logL ="
              ll_match <- regmatches(ll_str,
                regexpr("-2logL\\s*=\\s*([0-9.]+)", ll_str))
              if (length(ll_match) > 0) {
                ll_num <- as.numeric(sub(".*=\\s*", "", ll_match))
                if (!is.na(ll_num)) loglik_vals[i] <- -ll_num / 2
              }
            }
          }

          # Clean up temp directories
          unlink(rho_dirs, recursive = TRUE)

          loglik.rho <- transform(spatial$rho, loglik = loglik_vals)

          if (all(is.na(loglik_vals)))
            stop("All rho combinations failed during AR grid search. ",
                 "Try specifying rho manually.", call. = FALSE)

          # Full fit of winning rho for complete R result object
          rho.idx <- which.max(loglik.rho$loglik)
          mc_best <- mc
          mc_best$spatial$rho <- as.numeric(spatial$rho[rho.idx, ])
          mc_best$parallel <- FALSE
          ans <- suppressMessages(eval(mc_best, envir = caller_env))

        } else {
          ## Sequential grid search (default)
          eval.rho <- function(rho, mc, envir) {
            mc$spatial$rho <- rho
            tryCatch(
              suppressMessages(eval(mc, envir = envir)),
              error = function(e) {
                warning("rho = (", paste(round(rho, 4), collapse = ", "),
                        ") failed: ", conditionMessage(e), call. = FALSE)
                NULL
              }
            )
          }
          caller_env <- parent.frame()
          ans.rho <- lapply(seq_len(n_rho), function(i) {
            eval.rho(as.numeric(spatial$rho[i, ]), mc, envir = caller_env)
          })

          # Check that at least one combination succeeded
          succeeded <- !vapply(ans.rho, is.null, logical(1))
          if (!any(succeeded)) {
            stop("All rho combinations failed during AR grid search. ",
                 "Try specifying rho manually.", call. = FALSE)
          }

          loglik.rho <- transform(spatial$rho,
                                  loglik = suppressWarnings(
                                    vapply(ans.rho, function(x) {
                                      if (is.null(x)) NA_real_
                                      else as.numeric(logLik(x))
                                    }, numeric(1))
                                  ))
          rho.idx <- which.max(loglik.rho$loglik)
          ans <- ans.rho[[rho.idx]]
        }

        # Include estimation information
        ans$rho <- loglik.rho

        return(ans)
      }
    }
  }

  ## Generic specification
  if (!is.null(generic)) {
    ## TODO: multitrait case shoud check initial variance conformity
    ## and return a sensible default
    generic <- check_generic(generic, response = responsem)
  }
  
  ## Genomic specification
  if (!is.null(genomic)) {
    if (is.null(genetic))
      stop("'genomic' requires a 'genetic' component (e.g. add_animal model).",
           call. = FALSE)
    genomic <- check_genomic(genomic)
  }

  ## Initial variances specification
  ## We check even the NULL case, where the function returns the
  ## default initial variances for all random effects + residuals
  var.ini <- check_var.ini(var.ini, random, responsem)
  
  ## Whether the initial variances for each component are defaults
  has_var.ini <- 
    function(x) {
      if (eval(call('is.null', as.symbol(x)))) return(NA)
      else eval(call('attr', as.symbol(x), 'var.ini.default'))
    }
  var.ini.checks <- vapply(c('genetic', 'spatial', 'generic', 'var.ini'),
                           has_var.ini,
                           TRUE)

  ## Either all initial variances specified, or no specification at all
  if (any(var.ini.checks, na.rm = TRUE) && any(!var.ini.checks, na.rm = TRUE))
    stop(paste('Some initial variances missing.\n',
               'Please specify either all or none.'))

  ## Resuming supplies every initial variance, so an explicit specification
  ## would be silently discarded. Say so instead of quietly ignoring it.
  if (isTRUE(cont) && any(!var.ini.checks, na.rm = TRUE))
    stop("'cont = TRUE' takes the initial variances from ", progress_file,
         ".\n Drop the explicit var.ini specifications, or set cont = FALSE.",
         call. = FALSE)

  ## Issue a warning in the case of no specification
  if (all(var.ini.checks, na.rm = TRUE) && !isTRUE(cont)) {
    message(paste0('Using default initial variances given by ',
                  breedR.getOption('default.initial.variance'), '()\n',
                  'See ?breedR.getOption.\n'))
  }
  
  
  # Build a list of parameters and information for each effect
  effects <- build.effects(mf, genetic, spatial, generic, var.ini)

  ## Resume from a previous run's log.
  ##
  ## Done here rather than on the arguments, because the initial variances of
  ## the genetic, spatial and generic components live inside their own lists
  ## and have already been consumed by check_genetic()/check_spatial()/
  ## check_generic() further up. Only after build.effects() are the components
  ## named the way the log keys them -- the very same vector that
  ## parse_results() will use to label the results.
  resumed_round <- NULL
  if (isTRUE(cont)) {
    ckpt <- reml_checkpoint_var.ini(progress_file, effects,
                                    ntraits = ncol(responsem))
    resumed_round <- ckpt$round
    message('Resuming from round ', resumed_round, ' of ', progress_file)

    ## Assigning cov.ini directly bypasses effect_group(), where the variance
    ## is validated; reml_checkpoint_var.ini() has already run those checks.
    for (nm in setdiff(names(ckpt$var), 'residuals'))
      effects[[nm]]$cov.ini <- ckpt$var[[nm]]
    var.ini$residuals <- ckpt$var[['residuals']]
  }

  # Generate progsf90 parameters
  # TODO: Memory efficiency. At this point there are three copies of the 
  # dataset. One in data, one in mf (only needed variables)
  # and yet one more in pf90. This is a potential problem with large datasets.
  ## Build genomic OPTION lines if genomic analysis requested
  genomic_opts <- if (!is.null(genomic)) build_genomic_options(genomic) else NULL

  # Build the progsf90 object.
  # BLUPF90+ uses OPTION keywords to select between AI-REML and EM-REML.
  # PREGSF90 and BLUPF90+ need different OPTION sets:
  # - PREGSF90 needs SNP_file, QC options, G matrix options
  # - BLUPF90+ needs sol se, method VCE, heritability, readGimA22i
  method_opts <- 'method VCE'
  if (method == 'em') method_opts <- c(method_opts, 'EM-REML')

  pf90 <- progsf90(mf,
                   weights = weights,
                   effects,
                   opt = union(c('sol se', method_opts), progsf90.options),
                   res.var.ini = var.ini$residuals)

  if (!is.null(genetic) && method == 'ai') {
    ## Compute default heritability if possible
    trait_names <- colnames(model.response(mf))  # NULL for 1 trait
    pf90$parameter$options <-
      c(pf90$parameter$options,
        pf90_default_heritability(pf90$parameter$rangroup, trait_names))
  }

  # Temporary dir
  tmpdir <- tempdir()

  ## --- Genomic pre-processing with PREGSF90 ---
  pregs_out <- NULL
  if (!is.null(genomic)) {
    ped_for_xref <- get_pedigree.breedr_modelframe(effects)
    pipeline <- run_pregsf90_pipeline(genomic, genomic_opts, pf90,
                                      tmpdir, breedR.bin,
                                      pedigree = ped_for_xref)
    pf90 <- pipeline$pf90
    pregs_out <- pipeline$pregs_out
  }

  # Write the REML parameter file (or overwrite the PREGSF90 one)
  write.progsf90(pf90, dir = tmpdir)

  # Change to temporal directory to avoid specification of long paths
  # Avoids Issue #1
  cdir <- setwd(tmpdir)
  on.exit(setwd(cdir))
  
  ## write breedR_model for recovery
  ## Under `cont` the saved `effects` carry the recovered variances rather than
  ## the ones the call asks for, so record where they came from; otherwise the
  ## bundle silently claims a starting point the user never chose. `mcout` is
  ## left as typed -- it is the call, and falsifying it to avoid a replay
  ## re-triggering the resume would be the worse trade.
  save(
    effects, mf, method, mcout, resumed_round,
    file = file.path("breedR_model.RData")
  )

  ## Determine the breedR program to use (either local or remote)
  remote = FALSE
  submit = FALSE
  submit.id = ""
  if ( tolower(breedR.bin) == "remote" || tolower(breedR.bin) == "submit" ) {
    remote = TRUE
    submit.id = paste(gsub("[ :]", "-", date()), "---", as.integer(runif(1,min=1E8,max=1E9-1)), sep="")
    remote.bin = breedR.getOption('remote.bin')
    if( remote.bin == "path_to/breedR/bin/linux" ) {
      stop('breedR is not configured for remote computing. See ?breedR.options')
    }
    
    if( breedR.os('windows') ) {
      if( !breedR.cygwin.check() ) {
        stop(paste("Cannot find the CYGWIN installation:", breedR.getOption("cygwin")))
      }
      
      # Make sure the binaries are accesible
      breedR.cygwin.setPATH()
    }
    
    breedR.call = file.path(remote.bin, progsf90_files(breedR.os.type()))

    # Run either breedR.remote or breedR.submit
    if ( tolower(breedR.bin) == "remote" ) {
      ldir <- breedR.remote(submit.id, breedR.call)

      if( !identical(normalizePath(ldir), normalizePath(tmpdir)) ) stop('This should not happen')
      reml.out <- readLines(file.path(ldir, 'LOG'))
    }
    
    if ( tolower(breedR.bin) == "submit" ) {
      submit = TRUE
      reml.out <- breedR.submit(submit.id, breedR.call)
    }
  } else {
    breedR.call = file.path(breedR.bin, progsf90_files(breedR.os.type()))

    if (is.null(progress_file)) {
      reml.out <- system2(breedR.call,
                          input  = 'parameters',
                          stdout = ifelse(debug, '', TRUE))
    } else {
      ## Read the backend through a pipe and write the log ourselves.
      ##
      ## system2(stdout = <file>) would be the obvious way, but on Windows the
      ## redirect holds an exclusive lock: readLines(), file.copy() and even
      ## `type` fail with permission denied until the process exits. That would
      ## give a progress file that cannot be followed while it matters.
      ##
      ## stderr is deliberately not merged. The plain call above leaves it at
      ## its default of "", so it reaches the console and never enters
      ## reml.out; folding it in would change what parse_results() sees, and a
      ## runtime message landing between two rows of a `new G` block breaks
      ## extract_block() outright. It goes to a companion file instead.
      ##
      ## Only the binary is quoted; both redirections use relative names in the
      ## working directory. cmd.exe misparses a command line carrying more than
      ## one quoted absolute path around a redirection, so the error stream is
      ## collected locally and moved to its final name afterwards.
      ## Set the previous log aside, as late as possible. Everything that can
      ## reject the run -- the binaries check, the component checks, the AR-grid
      ## guard, the checkpoint parse, PREGSF90 -- has already happened, so a
      ## failed resume leaves the log where the user left it, under the name
      ## their retry will use.
      bak <- NULL
      if (isTRUE(cont)) {
        bak <- paste0(progress_file, '.', format(Sys.time(), '%Y%m%d-%H%M%S'))
        if (file.rename(progress_file, bak)) {
          message('Previous log kept as ', bak)
        } else {
          bak <- NULL
          warning('Could not preserve the previous log; it will be overwritten.',
                  call. = FALSE)
        }
      }

      ## Open the log *before* starting the backend. This is the one fallible
      ## step here, and doing it first means a failure leaves no child process
      ## to orphan -- an orphaned backend holds files in tempdir() and breaks
      ## the next fit in the same session.
      ##
      ## It is also the last thing that can fail after the rename above, so put
      ## the previous log back if it does; otherwise a resume that cannot open
      ## its log would consume the very file it was resuming from.
      lcon <- tryCatch(
        file(progress_file, 'w'),
        error = function(e) {
          if (!is.null(bak)) file.rename(bak, progress_file)
          stop(e)
        })
      ## Idempotent: the connection is closed as soon as the run ends, and
      ## this is only the safety net for an error part-way through.
      on.exit(try(close(lcon), silent = TRUE), add = TRUE)

      ## `con` is deliberately not registered with on.exit: close() on a pipe
      ## waits for the child to finish, so on an interrupted multi-day fit it
      ## would hang R instead of returning. The read loop below is therefore
      ## still a window in which an error leaves the backend running.
      writeLines('parameters', 'pf90_stdin')
      con <- pipe(paste0(shQuote(breedR.call),
                         ' < pf90_stdin 2> pf90_stderr'), 'r')

      message('Streaming REML progress to ', progress_file)

      reml.out <- character(0)
      repeat {
        l <- readLines(con, n = 1L, warn = FALSE)
        if (!length(l)) break
        reml.out <- c(reml.out, l)
        writeLines(l, lcon)
        flush(lcon)
        if (debug) cat(l, '\n')
      }
      status <- close(con)
      close(lcon)

      ## Keep the error stream beside the log, but only when there is one.
      err_file <- paste0(progress_file, '.err')
      unlink(err_file)
      if (file.exists('pf90_stderr') && file.info('pf90_stderr')$size > 0)
        file.copy('pf90_stderr', err_file)

      ## Restore the contract of stdout = TRUE, which the code below relies on.
      if (!identical(as.integer(status), 0L))
        attr(reml.out, 'status') <- as.integer(status)
    }
  }
  
  
  
  ## A streamed run still reports a failing backend under debug, where the
  ## check below is skipped along with the parsing.
  if (debug && !is.null(attr(reml.out, 'status')))
    stop_progsf90_failure(reml.out, attr(reml.out, 'status'))

  if( !debug ) {
    # Check whether the PROGSF90 binary exited with an error.
    # Include the tail of the REML output so users can see the actual diagnostic.
    if (!is.null(attr(reml.out, 'status')))
      stop_progsf90_failure(reml.out, attr(reml.out, 'status'))

    if( !submit ) {
      ## Save reml.out for recovery
      writeLines(reml.out, file.path(tmpdir, "reml.out"))
      
      # Parse solutions
      ans <- parse_results(file.path(tmpdir, 'solutions'), effects, mf, reml.out, method, mcout)

      ## Recorded here rather than inside parse_results(), whose signature and
      ## output existing tests pin down. The cost is that these do not survive
      ## a breedR.qget() recovery, which rebuilds the object by calling
      ## parse_results() directly -- immaterial while remote fits are refused.
      if (!is.null(progress_file)) ans$reml$progress_file <- progress_file
      if (!is.null(resumed_round)) ans$reml$resumed_from <- resumed_round
    } else {
      # Submitted job (solutions are parsed later with breedR.qget)
      ans <- list(id = submit.id,
                  effects = effects,
                  mf = mf,
                  method = method,
                  mcout = mcout)
    }
    
    class(ans) <- c('breedR', 'remlf90')  # Update to merMod in newest version of lme4 (?)

    ## Attach genomic QC results if PREGSF90 was run
    if (!is.null(genomic) && !is.null(pregs_out)) {
      ans$genomic <- parse_pregsf90_qc(tmpdir)
      ans$genomic$pregsf90_output <- pregs_out
      # Record whether the plain G-inverse was saved, so postgsf90() can tell
      # whether this fit is GWAS-ready.
      ans$genomic$save_ginverse <- isTRUE(genomic$save_ginverse)
    }
  } else {
    file.show('parameters')
    ans = NULL
  }

  return(ans)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Internal methods  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%#



#%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Interface methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @export
coef.remlf90 <- function(object, ...) { 
  ans <- 
    rbind(
      splat(rbind)(lapply(object$fixed, splat(rbind))),
      splat(rbind)(lapply(object$ranef, splat(rbind)))
    )

  return(structure(ans[, "value"], names = rownames(ans)))
}

#' @export
extractAIC.remlf90 <- function(fit, scale, k, ...) {
  return(fit$fit$AIC)
}
  
#' @method fitted remlf90
#' @export
fitted.remlf90 <- function (object, ...) {

  ef.types <- lapply(object$effects, effect_type)
  
  mml <- model.matrix(object)
  
  vall <- c(fixef(object), ranef(object))
  
  ## Match order
  stopifnot(setequal(names(mml), names(vall)))
  vall <- vall[names(mml)]
  
  silent.matmult.drop <- function(x, y) {
    suppressMessages(drop(as.matrix(x %*% y)))
  }

  ## Multiply component-wise
  ## dimensions: observation; trait (if ntraits > 1); effect
  comp.mat <- mapply(silent.matmult.drop, mml, vall, SIMPLIFY = 'array')
  
  
  # Linear Predictor / Fitted Values
  ndim <- length(dim(comp.mat))
  eta <- rowSums(comp.mat, dims = ndim - 1)

  return(eta)
}



#' Extract the fixed-effects estimates
#'
#' Extract the estimates of the fixed-effects parameters from a fitted model.
#' @name fixef
#' @title Extract fixed-effects estimates
#' @aliases fixef fixed.effects fixef.remlf90 fixef.breedR
#' @param object any fitted model object from which fixed effects estimates can
#' be extracted.
#' @param \dots not used
#' @return a named list of dataframes with the estimated coefficients and their
#'   standard error
#' @keywords models
#' @examples
#' \dontrun{
#'     res <- remlf90(phe_X ~ gg + bl, data = globulus)
#'     fixef(res)
#' }
#' @importFrom nlme fixef
#' @export fixef
#' @export
fixef.remlf90 <- function (object, ...) {
  ans <- get_estimates(object$fixed)
  class(ans) <- 'breedR_estimates'
  return(ans)
}


#' @method logLik remlf90
#' @export
logLik.remlf90 <- function (object, REML = TRUE, ...) {
  # TODO: Revise this, N parameters, df, N obs.
  # I set up things such that the AIC gives what REMLF90 says
  # But I am not sure if it is the right way.
  reml.out <- object$reml$output

  ## Number of (estimated) parameters (a.k.a. degrees of freedom)
  npar.idx <- grep('# parameters=', reml.out)
  npar <- ifelse(identical(length(npar.idx), 1L),
                 as.numeric(strsplit(reml.out[npar.idx],
                                     split=' # parameters= +')[[1]][2]),
                 'unknown')
  
  ## Rank
  ## From stats:::logLik.lm it turns out that 
  ## df = rank + 1
  rank <- ifelse( npar != 'unknown', npar - 1, 'unknown')
  
  ## Number of obervations
  N  <- nobs(object)
  N0 <- N
  if( REML & rank != 'unknown' ) {
    N <- N - rank
  }
  
  ans = -object$fit[['-2logL']]/2
  attr(ans, 'nall') <- N0
  attr(ans, 'nobs') <- N
  attr(ans, 'df') <- npar
  class(ans) ='logLik'
  ans
}

#' @method model.frame remlf90
#' @export
model.frame.remlf90 <- function (formula, ...) {
  formula$mf
}



#' @importFrom stats nobs model.response
#' @method nobs remlf90
#' @export
nobs.remlf90 <- function (object, ...) {
  if(!is.null(w <- object$weights))
    sum(w != 0)
  else
    NROW(model.response(object$mf))
}



#' Spatial plot of a model's fit components
#' 
#' Plots the predicted values of the component effects of the phenotype.
#' 
#' @param x A \code{breedR} object
#' @param type Character. Which component is to be represented in the map.
#' @param z Optional. A numeric vector to be plotted with respect to the spatial
#'   coordinates. Overrides \code{type}.
#' @param ... Further layers passed to \code{\link[ggplot2]{ggplot}}.
#'   
#' @method plot remlf90
#' @importFrom stats model.response fitted residuals
#' @export
plot.remlf90 <- function (x, type = c('phenotype', 'fitted', 'spatial', 'fullspatial', 'residuals'), z = NULL, ...) {
  
  type = match.arg(type)
  
  if( length(coord <- coordinates(x)) == 0) {
    stop(paste('Missing spatial structure. Use coordinates(',
               deparse(substitute(x)),
               ') <- coord', sep = ''))
  }

  ## Otherwise, this should be a matrix of coordinates
  coord <- as.matrix(coord)

  # Argument z is used for plotting a custom spatial variable
  if( !is.null(z) ) {
    z <- as.vector(z)
    if( !is.numeric(z) | length(z) != nrow(coord) )
      stop(paste("'z' must be a vector of length", nrow(coord), ".\n"))
    
    # Determine the type of scale, depending on whether the vector includes 0
    if( min(z) < 0 & max(z) > 0 ) sc = 'div'
    else sc = 'seq'
    
    p <- spatial.plot(data.frame(coord, z = z), scale = sc)
    
  } else {
    
    if(type == 'phenotype') {
      
      resp <- model.response(x$mf)
      
      p <- spatial.plot(data.frame(coord, z = resp), scale = 'seq')
    }
    
    
    if(type == 'fitted') {
      
      p <- spatial.plot(data.frame(coord, z = fitted(x)), scale = 'seq')
    }
    
    
    if(type == 'spatial' | type == 'fullspatial'){
      if( x$components$spatial ) {
        
        if( type == 'spatial' ) {
          value <- as.vector(model.matrix(x)$spatial %*% ranef(x)$spatial)
        } else {
          ## fullspatial case
          ## assume only one spatial effect in the effect_group spatial
          inc.mat <- model.matrix(x$effects$spatial, fullgrid = TRUE)[[1]]
          coord <- attr(inc.mat, 'coordinates')
          stopifnot(!is.null(coord))
          value <- as.vector(inc.mat %*% ranef(x)$spatial)
        }
        spdat <- data.frame(coord,
                            z     = value,
                            model = x$call$spatial$model)
        
        
        p <- spatial.plot(spdat, scale = 'div') + facet_wrap(~ model)
      } else stop('This model has no spatial effect')
    }
    
    if(type == 'residuals') {
      p <- spatial.plot(data.frame(coord, z = residuals(x)), scale = 'div')
    }
  }
  
  # Further args
  if( !missing(...) ) {
    p <- p + ...
  }
  
  p
}



#' @method plot ranef.breedR
# @describeIn ranef.breedR
#' @export
plot.ranef.breedR <- function(x, y, ...) {
  ## dotplot for each random effect
  ## only makes sense for random effects with a few levels
  ## thus we exclude from the plot genetic, spatial or other 
  ## random effects with many levels
  max.nl <- 30
  x <- x[sapply(x, length) < 30]
  
  ranef2df <- function(x) {
    if( is.null(nm <- names(x)) ) nm <- seq.int(x)
    data.frame(level= nm,
               BLUP = as.vector(x),
               ymin = as.vector(x) - 1.96*attr(x, 'se'),
               ymax = as.vector(x) + 1.96*attr(x, 'se'))
  }
  
  if( length(x) ) {
    pl <- lapply(seq.int(x), function(i) data.frame(effect = names(x)[i],
                                                    ranef2df(x[[i]])))
    pd <- do.call(rbind, pl)
    
    ggplot(pd, aes_string(x = "level", y = "BLUP", ymin = "ymin", ymax = "ymax")) + 
      geom_pointrange() + 
      coord_flip()
  } else message('No suitable random effects to plot')
}

#' @method print breedR_estimates
# @describeIn ranef
#' @export
print.breedR_estimates <- function(x, ...) {
  attr2df <- function(x) {
    data.frame(value = x, `s.e.` = attr(x, 'se'))
  }
  ans <- lapply(x, attr2df)
  print(ans, ...)
  invisible(x)
}


# predict.remlf90 <- function (object, ...) {
#   
# }
# 
# print.remlf90 <- function (object, ...) {
#   
# }


#' Extract the modes of the random effects
#' 
#' Extract the conditional modes of the random effects from a fitted model 
#' object.  For linear mixed models the conditional modes of the random effects
#' are also the conditional means.
#' 
#' This method is modeled a bit after \code{\link[lme4]{ranef}}. However, it is
#' independent and does not inherit from it. In particular, it always returns
#' the conditional variance (argument condVar in lme4).
#' 
#' @name ranef
#' @aliases ranef ranef.remlf90 ranef.breedR
#' @param object a fitted models with random effects of class 
#'   \code{\link{remlf90}}.
#' @param ... not used
#' @return An object of class \code{ranef.breedR} composed of a list of vectors 
#'   or matrices (multitrait case), one for each random effect. The length of
#'   the vectors are the number of levels of the corresponding random effect.
#'   
#'   Each random effect has an attribute called \code{"se"} which is a vector 
#'   with the standard errors.
#'   
#'   Additionally, depending of the nature of the random effect, there may be 
#'   further attributes. The pedigree will be given for genetic random effects 
#'   and the spatial prediction grid for the spatial random effects
#'   
#' @note To produce a (list of) \dQuote{caterpillar plots} of the random effects
#'   apply \code{\link{plot}} to the result of a call to \code{ranef}.
#' @examples
#' \dontrun{
#' res <- remlf90(phe_X ~ bl,
#'                genetic = list(model = 'add_animal',
#'                               pedigree = globulus[, 1:3],
#'                               id = 'self'),
#'                data = globulus)
#' str(rr <- ranef(res))
#' plot(rr)
#' }
#' @importFrom nlme ranef
#' @export ranef
#' @export
ranef.remlf90 <- function (object, ...) {

  ## ranef() will provide the model's random effects
  ## and further methods will let the user compute their 'projection'
  ## onto observed individuals (fit) or predict over unobserved individuals (pred)
  
  ans <- get_estimates(object$ranef)
  
  ## Additional attributes
  
  ## Genetic component: names of individuals
  if( object$components$pedigree ){
    
    # Indices (in ranef) of genetic-related effects (direct and/or competition)
    gen.idx <- grep('genetic', names(ans))
    nm <- get_pedigree(object)@label
    
    for (k in gen.idx) attr(ans[[k]], 'names') <- nm
    
  }
  
  ## Other random effects with no particular treatment
  idx <- grep('genetic|spatial', names(ans), invert = TRUE)
  
  for(x in names(ans[idx])) {
    if("effect_group" %in% class(object$effects[[x]])){
      if("generic" %in% class(object$effects[[x]]$effects[[1]]) &
         ! is.null(rownames(object$effects[[x]]$effects[[1]]$structure.matrix)))
        attr(ans[[x]], 'names') <- 
          rownames(object$effects[[x]]$effects[[1]]$structure.matrix)
    } else
      attr(ans[[x]], 'names') <- 
        colnames(attr(model.matrix(object)$random[[x]], 'contrasts'))
  }
  
  class(ans) <- c('ranef.breedR', 'breedR_estimates')
  return(ans)
}


#' Covariance matrix of a fitted remlf90 object
#' 
#' Returns the variance-covariance matrix of the specified random effect.
#' 
#' @param object a fitted model of class \code{remlf90}
#' @param effect the structured random effect of interest
#' @param ... Not used.
#' 
#' @method vcov remlf90
#' @export
vcov.remlf90 <- function (object,
                          effect = c('spatial',
                                     'genetic',
                                     'genetic_direct',
                                     'genetic_competition',
                                     'pec'),
                          ...) {
  
  effect <- match.arg(effect)
  
  ## Check that the effect exists
  if (!effect %in% get_efnames(object$effects)) {
    stop(paste('There is no', effect, 'effect in this object.\n'))
  }


  efgroup <- vapply(object$effects, function(x) effect %in% names(x$effects), TRUE)
  if (sum(efgroup) == 1) {
    ef <- object$effects[[which(efgroup)]]$effects[[effect]]
  } else {
    ## genetic or spatial cases. TODO: Improve this.
      ef <- object$effects[[effect]]$effects[[1]]
  }
  ## TODO:
  ## Maybe there should be a method 'remlf90' for 'get_structure'
  ## with analogous behaviour as model.matrix or ranef
  ## So we can simply say:
  ## U <- get_structure(object)[[effect]]
  
  ## Underlying covariance matrix
  U <- get_structure(ef)

  ## need to invert?
  if(attr(U, 'type') == 'precision')
    U <- solve(U)
  
  ## Incidence matrix
  B <- model.matrix(object)[[effect]]
  
  ## Scaling parameter
  ## TODO: Fix this
  if (any(efgroup)) {
    sigma2 <- object$var[[names(which(efgroup))]][effect, effect]
  } else {
    sigma2 <- object$var[effect, 1]
  }
  
  V <- suppressMessages(sigma2 * B %*% U %*% Matrix::t(B))
  return(V)
}


#' @method residuals remlf90
#' @importFrom stats model.response
#' @export
residuals.remlf90 <- function (object, ...) {
  # TODO: to be used when na.action is included in remlf90
  #   naresid(object$na.action, res)
  # TODO: add a scale parameter to allow studentization
  #     first need to determine the right sigma for each residual
  model.response(object$mf) - fitted(object)
}

#' @method summary remlf90
#' @importFrom stats logLik AIC BIC
#' @export
summary.remlf90 <- function(object, ...) {
  
  # If this is a submitted job, return the corresponding qstat object
  # instead of a summary.remlf90 object
  if( exists('id', object) ) {
    # Submitted Job
    return(breedR.qstat(object))
  } 
  
  # Literal description of the model
  effects <- paste(names(object$components), sep=' and ')
  title <- paste('Linear Mixed Model with', 
                 paste(effects, collapse = ' and '), 
                 ifelse(length(effects)==1, 'effect', 'effects'), 
                 'fit by', object$reml$version)
  
  # Formula
  fml.spec <- names(object$components)[unlist(object$components)]
  fml <- paste(c(deparse(attr(object$mf, 'terms')), fml.spec),
               collapse = ' + ')
  
  # Coefficients
  # TODO: How to avoid showing the unused levels
  coef <- splat(rbind)(lapply(object$fixed, splat(rbind)))
#   colnames(coef) <- c('Estimate')
  
  # Model fit measures
  # AIC and BIC might fail if logLik fails to retrieve
  # appropriate df or nobs attributes
  llik <- logLik(object)
  AICframe <- data.frame(AIC = tryCatch(AIC(llik), 
                                        error = function(e) object$fit$AIC), 
                         BIC = tryCatch(BIC(llik),
                                        error = function(e) 'unknown'),
                         logLik = as.vector(llik),
#TODO                          deviance = dev[["ML"]],
#                          REMLdev = dev[["REML"]],
                         row.names = "")
  
  ## Variance components
  ## Either a numeric matrix (1 trait)  - !is.list
  ## or a matrix of matrices (>1 trait) - is.list && is.matrix
  ## or a list of matrices (>1 trait, em: no SE) - !is.matrix
  if (is.list(object$var) && is.matrix(object$var)) {
    ## multiple-trait case: a 2-col (est; SE) matrix of covariance matrices

    ## transform the list of est and se symmetric matrices into a data.frame
    ## with values from lower triangule only
    varnm_df <- function(x, nm) lmat2df(object$var[nm, ], nm)
    
    ## list of data-frames with variance estimates and SE
    var_ldf <- lapply(rownames(object$var), varnm_df, x = object$var)
    
    object$var <- splat(rbind)(var_ldf)
  }
  
  ans <- c(object, 
           model.description = title, 
           formula = fml,
           model.fit = list(AICframe),
           coefficients = list(coef)
           )
  class(ans) <- 'summary.remlf90'
  ans
}



## This is modeled a bit after  print.summary.lm :
#' @method print remlf90
#' @export
print.remlf90 <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  # If this is a submitted job, return the corresponding qstat object
  # instead of a summary.remlf90 object
  if( exists('id', x) ) {
    # Submitted Job
    return(print(breedR.qstat(x)))
  } 
  
  cat(x$model.description, '\n')
  if(!is.null(x$call$formula))
    cat("Formula:", x$formula,"\n")
  if(!is.null(x$call$data))
    cat("   Data:", deparse(x$call$data),"\n")
  if(!is.null(x$call$subset))
    cat(" Subset:", x$call$subset,"\n")
  print(x$model.fit, digits = digits)
  
  if( x$components$spatial & !is.null(x$spatial$name)) {
    switch(x$spatial$name,
           AR = cat(paste("\nAutoregressive parameters for rows and columns: (",
                          paste(x$spatial$model$param, collapse = ', '),
                          ")\n", sep = '')),
           splines = cat(paste("\nNumber of inner knots for rows and columns: (",
                               paste(x$spatial$model$param, collapse =', '),
                               ")\n", sep = ''))
    )
  }
  
  invisible(x)
}



## AIREMLF90 reports three numbers for each OPTION se_covar_function: the
## plug-in value of the function at the REML solution ('mean'), and the mean
## and standard deviation of the draws it takes from the asymptotic
## distribution of the variance components ('sample mean', 'sample sd').
## The plug-in value is the estimate; the sampling SD is its standard error.

## Divergence beyond this fraction of the sampling SD is worth flagging.
## The Monte Carlo standard error of the sample mean is sd/sqrt(5000), about
## 0.014 sd, so a tenth of a SD is well beyond mere sampling noise.
funvars_divergence_tol <- 0.1

## Reshape the 3 x n matrix of parsed functions into the printed table
funvars_table <- function(funvars) {
  ans <- t(funvars)
  colnames(ans) <- c("Estimate", "Sample Mean", "S.E.")
  ans
}

## Which functions have a sample mean far away from the plug-in estimate?
## Compare multiplicatively rather than scaling the gap by the SD, which
## would give NaN whenever the sampling SD is zero.
funvars_divergent <- function(funvars, tol = funvars_divergence_tol) {
  if (!length(funvars)) return(logical(0))
  ans <- abs(funvars['mean', ] - funvars['sample mean', ]) >
    tol * funvars['sample sd', ]
  ans[is.na(ans)] <- FALSE
  ## indexing a single-column matrix by row drops the column name
  stats::setNames(ans, colnames(funvars))
}


## This is modeled a bit after  print.summary.lm :
#' @method print summary.remlf90
#' @importFrom stats printCoefmat
#' @export
print.summary.remlf90 <- function(x, digits = max(3, getOption("digits") - 3),
                                  correlation = TRUE, symbolic.cor = FALSE,
                                  signif.stars = getOption("show.signif.stars"), ...) {
  
  # cat(x$model.description, '\n')
  if(!is.null(x$formula))
    cat("Formula:", x$formula,"\n")
  if(!is.null(x$call$data))
    cat("   Data:", deparse(x$call$data),"\n")
  if(!is.null(x$call$subset))
    cat(" Subset:", x$call$subset,"\n")
  print(x$model.fit, digits = digits)

  parl <- get_param.remlf90(x)
  if (!is.null(parl)){
    cat("\nParameters of special components:\n")  
    for (i in seq_along(parl)) {
      for (j in seq_along(parl[[i]])) {
        comp.nm <- ifelse(j==1,
                          paste0(names(parl)[i], ':'),
                          paste(character(nchar(names(parl)[i]) + 1),
                                collapse = ''))
        model.nm <- paste0(names(parl[[i]])[j], ':')
        cat(comp.nm, model.nm, parl[[i]][[j]])
      }
    }
  }
  cat("\n")

  cat("\nVariance components:\n")
  print(x$var, quote = FALSE, digits = digits, ...)
  
  if (length(x$funvars)) {
    cat("\nFunctions of variance components:\n")
    print(funvars_table(x$funvars), quote = FALSE, digits = digits, ...)
    div <- funvars_divergent(x$funvars)
    if (any(div)) {
      cat("\nNote: for ", paste(names(div)[div], collapse = ", "),
          " the Sample Mean departs from the Estimate.\n",
          "  The S.E. is obtained by drawing the variance components from their\n",
          "  asymptotic distribution at the REML solution and evaluating the\n",
          "  function at each draw. Because the function is non-linear, the mean\n",
          "  of those draws need not equal the function of the means, and draws\n",
          "  falling outside the parameter space are pulled back to its boundary.\n",
          "  Both effects grow with the S.E., so a gap here means the variance\n",
          "  components are poorly determined: read the S.E. as indicative only.\n",
          sep = "")
    }
  }

  cat('\nFixed effects:\n')
  printCoefmat(x$coefficients)
  invisible(x)
}
