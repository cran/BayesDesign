\name{optimal_TwoStage}
\alias{optimal_TwoStage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain design settings for two-stage Bayesian Single-Arm Phase II Trial with Time-to-Event Endpoints
}
\description{
Obtain design parameters, type I error, power and operating characteristics of the Bayesian Single-Arm Phase II Trial Designs with Time-to-Event Endpoints (Wu et al. 2021). The exponential distribution is assumed for the survival time. The gamma prior is used here
}
\usage{
optimal_TwoStage(alphacutoff, powercutoff, S0, x,
                 ta, tf, a = 2, delta, frac = .5,
                 ntrial, complete = "partial", seed = 8232)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphacutoff}{
the desired type I error to be controlled
  }
  \item{powercutoff}{
the desired power to be achieved
  }
  \item{S0}{
the survival probability at timepoint x
  }
  \item{x}{
the survival probability S0 at timepoint x
  }
  \item{ta}{
accrual duration
  }
  \item{tf}{
follow-up duration
  }
  \item{a}{
shape parameter of prior distribution. The default value is a = 2
  }
  \item{delta}{
hazard ratio
  }
  \item{frac}{
a information fraction for interim analysis. The fefault value is
frac = 0.5
  }
  \item{ntrial}{
the number of simulated trials
  }
  \item{complete}{
whether output the full or partial information. The default value is complete = "partial". If want to show full results, it would be complete = "complete"
  }
  \item{seed}{
the seed. The default value is seed = 8232
  }
}

\value{
\code{optimal()} depending on the argument "complete", it returns a vector of partial information/complete information which includes:

partial information:
(1) m1: number of events at stage 1
(2) n1: number of patients at stage 1
(3) k1: total observation time at stage 1
(4) m: number of events of the whole design
(5) n: number of patients of the whole design
(6) k: total observation time of the whole design
(7) typeI: type I error of the whole design
(8) power: power of the whole design
(9) PET1: early stopping probabilites under alternative hypothesis
(10) ES1: expected sample size under alternative hypothesis
(11) PET0: early stopping probabilites under null hypothesis
(12) ES0: expected sample size under null hypothesis

full information:
(1) eta: cutoff point of "Go" at final stage of analysis
(2) xi: cutoff point of "no-Go" at final stage of analysis
(3) m1: number of events at stage 1
(4) n1: number of patients at stage 1
(5) k1: total observation time at stage 1
(6) m: number of events of the whole design
(7) n: number of patients of the whole design
(8) k: total observation time of the whole design
(9) typeI: type I error of the whole design
(10) power: power of the whole design
(11) PET1: early stopping probabilites under alternative hypothesis
(12) ES1: expected sample size under alternative hypothesis
(13) PET0: early stopping probabilites under null hypothesis
(14) ES0: expected sample size under null hypothesis

}

\references{
Jianrong Wu, Haitao Pan, Chia-Wei Hsu (2021). "Bayesian Single-Arm Phase II Trial Designs with Time-to-Event Endpoints." Pharmaceutical Statistics. Accepted
}
\author{
Chia-Wei Hsu, Haitao Pan, Jianrong Wu
}

\examples{

### Design 1

# H0 vs. H1: 17\% vs. 40\% (4-month PFS)
# that is, S0 = 0.17, and hazard ratio, e.g., delta = 0.517
# x = 4

optimal_TwoStage(alphacutoff = 0.1, powercutoff = 0.8, S0 = 0.17,
                 x = 4, ta = 6, tf = 6, delta = 0.517, ntrial = 10)

### Design 2

# H0 vs. H1: 17\% vs. 30\% (4-month PFS)
# that is, S0 = 0.17, and hazard ratio, e.g., delta = 0.679
# x = 4

optimal_TwoStage(alphacutoff = 0.1, powercutoff = 0.8, S0 = 0.17,
                 x = 4, ta = 6, tf = 6, delta = 0.679, ntrial = 10)
}
