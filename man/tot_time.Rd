\name{tot_time}
\alias{tot_time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sum up transformed observation time for each patient
}
\description{
Sum up transformed observation time for each patient to get
U in order to determine the trial: (1) goes to second stage
(2) stops for futility (3) declares the treatment is promising and
warrants for further study in a large scale phase III trial
(4) declares the treatment is unpromising and
is not worth for further study.
}

\usage{
tot_time(obs_time, S0, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs_time}{
a vector. Each element represents an observation time of the patient
  }
  \item{S0}{
the survival probability at timepoint x
  }
  \item{x}{
the survival probability S0 at timepoint x
  }
}
\value{
the function returns the total transformed observation time for all patients
}
\references{
Jianrong Wu, Haitao Pan, Chia-Wei Hsu (2021). "Bayesian Single-Arm Phase II Trial Designs with Time-to-Event Endpoints." Pharmaceutical Statistics. Accepted
}
\author{
Chia-Wei Hsu, Haitao Pan, Jianrong Wu
}

\examples{
obs_time <- c(3.003, 11.987, 4.306, 2.561, 1.575, 0.329, 1.940,
              0.869, 7.481, 1.861, 7.279, 0.007, 6.485, 1.981,
              4.257, 0.967, 2.619, 0.040, 0.426, 4.628)
S0 <- 0.17
x <- 4
tot_time(obs_time = obs_time, S0 = S0, x = x)
}
