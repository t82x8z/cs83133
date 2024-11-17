java c
STAT 314 Assignment 5 2024
Question One: Optimising a Metropolis algorithm with a random walk jumping distribution (8 marks)
Computer game developers are trialling a new game. For this project, they recruited 80 participants, 40 female (sex = 1) and 40 male (sex = 0), each of whom were graded of their level of gaming experience (level: 1-4 representing minimal to high). The participants were provided with the game, and the number of times y they played the game was recorded. You are asked by the researchers to fit generalised linear models to the data, provided by you on LEARN as gaming.csv, with sex as a categorical predictor and gaming as a continuous predictor.
Note, in generalised linear models, rather than estimating effects from the response data directly, we model through a link function, η(θ), and assume η(θ)i = x′iβ. The link function can be determined by re-arranging the likelihood of interest into the exponential family format,



As the response variable is a count, the game developers suggesting fitting a Poisson regression. The Poisson pmf is
The Poisson probability mass function can be re-arranged into the exponential family format as follows,

which indicates the link function is η(λ) = log(λ). Hence when fitting the Bayesian Poisson regression, assume that log(λi) = x′iβ.
Fit a Bayesian Poisson regression using Metropolis sampling. Assume the prior for β is N (0, 100Ip), i.e. independent with zero mean and variance = 100.
To get certain quantities for the proposal distribution and predictor, fit the glm shown below:
gaming= read.csv("gaming.csv",header=TRUE)
mod<-glm(y~as.factor(sex)+level,data=gaming,family='poisson')
X=model.matrix(mod)
Sigma=vcov(mod);Sigma
From the glm, extract the design matrix X. For the proposal distribution, use a Normal distribution with mean θ(t−1) and variance-covariance matrix c2ˆΣ where Σ is the variance-covariance matrix from the glm fit. Consider three candidates for c, 1.6/√p, 2.4/√p, 3.2/√p, where p is the number of parameters estimated. Run the Metropolis algorithm for 10000 iterations, and discard the first 5000. Report the following:
• Check each resulting chain converges to the same distribution. To do this, you may find installing the R package coda helpful, and providing graphical summaries helpful. (2 marks)
• The proportion of candidate draws that were accepted (1 mark).
• The effective sample size for each chain. (1 mark)
• Based on your results to the previous bullet points, what do you think is the best choice for c. Does it approximately match results stated in class on efficiency and optimal acceptance rate? (2 marks)
There are 2 marks available for correctly written the code for implementing Metropolis sampling in this problem.
Question Two: Sequential analysis (7 marks)
The Bayesian way of thinking readily lends itself to combining information from sequential experiments. To demonstrate, consider the following data extracted from the HealthIron study.
Serum ferritin levels were measured for two samples of women, one of C282Y homozygotes (n = 88) and the other of women with neither of the key mutations (C282Y and H63D) in the HFE gene, so-called HFE ‘wildtypes’(n = 242). The information available is
• idnum: Participant id.
• homc282y: Indicator whether individual is Homozygote (1) or Wildtype (0).
• time: Time since onset of menopause, measured in years.
• logsf: The natural logarithm of the serum ferritin in µg/L.
The data required to answer this question is Hiron2021.csv, which can be downloaded from LEARN.
a) Fit a standard linear regression,
E(logsf) = β0 + β1time
with responses restricted to only those who are homozygote (homc282y = 1). This can be done using the lm function in R. Report the estimated coefficients ˆβ, estimated error variance, ˆσe2 and (X′X)−1. State what priors would need to be assumed for β and τ = (σe2)−1代 写STAT 314 Assignment 5 2024
代做程序编程语言such that estimates reported here are the parameters of the posterior distribution of β, τ . [1 mark]
b) The sequential analysis. For this, you will fit a Bayesian regression using a Gibbs sampler to only the wildtype (homc282y=0) data. However, you must use priors for β and τ such that your results will be equivalent to fitting a linear model to all the data available. These priors, and resulting conditional posteriors are shown below.
The conditional posteriors implied by using priors taken from the output in a) are very similar to the case in the later part of the Week 11-12 lectures when the prior for β is N (β0, σ2βK). To remind you of this, this is slides 16-20 in the Week 11-12 lecture, which are roughly re-produced below:
Assume the general prior for β, β ∼ N (β0, σ2βK), where K is an arbitrary variance-covariance matrix. Then we can deduce the result for various special cases.
• To make derivations easier, we will work with inverse variances (τ = (σ2)−1) rather than variances, and assume inverse variance components a priori are distributed Ga(αi, γi). Further, we will assume that K is known and does not need to be estimated.
• The likelihood p(y|β, τe, X) is

• The priors are

• The joint distribution p(y, β, τe, X, τβ) = p(y|β, τe, X)p(β|, τβ, K)p(τβ|αβ, γβ)× p(τe|αe, γe) is

• As the posterior distribution is proportional to the joint distribution, the task of determining conditional posteriors is equivalent to determining the distribution kernel for the parameter(s) of interest.
• The component of the joint distribution that is a function of τe is,

which corresponds to a gamma kernel, such that
p(τe|y, β, X) = Ga(αe + n/2, γe + (y − Xβ)′(y − Xβ)/2).
• The component of the joint distribution that is a function of τβ is,

which corresponds to a gamma kernel, such that
p(τβ|y, β, K) = Ga(αβ + p/2, γβ + (β − β0)′K−1(β − β0)/2).
• The component of the joint distribution that is a function of β is,

which corresponds to a normal kernel, such that p(β|y, X, β0, K, τe, τβ) is multivariate-normal with mean =(τeX′X + τβK−1)−1(τeX′y + τβK−1β0) and variance-covariance matrix (τeX′X + τβK−1)−1.
Note: If x is multivariate normal N(µ, Σ), the kernel is e − 2/(x−µ)′Σ −1(x−µ)∝ e− 2/µ′Σ−1µeµ′Σ−1x
In the specific case,of analysing the HealthIron study sub-group by sub-group, β0 = βˆfrom a, σ2βfixed and K = σ2(X′1X1)−1, which in turn means τe = τβ. Note by writing the prior in this form, we have made β conditional on σ2 or equivalently τ = (σ2)−1. Note y1, X1, ˆβ1, ˆσ2 1refers to results from part a), the same quantities from the subset where homc282y=0 is not sub-scripted and τ = (σ2)−1

The conditional posterior for τ is a little more complicated than that shown in the lecture notes, as we have made the prior for β conditional on τ , however the prior for τ is still Ga(α, γ) with α = (n1 − p)/2 and β = (n1 − p)ˆσ2 1/2. To determine the posterior for τ , consider the modified joint distribution,

from which we can deduce the conditional posterior for τ is,

• For the Gibbs sampler, run three chains for 10000 iterations. Discard the first 1000 iterations as burn-in and then remove every second remaining iteration to reduce auto-correlation. Check that chains have converged. When storing results, convert τ back to σ2. Report posterior means, standard deviations and 95 % central credible intervals (including interpretation) for β0, β1, σ2combining results for the three chains. If you want to check your coding is correct, compare your results to fitting a linear model to all data in the HealthIron study.
Marks are as follows.
• 2 marks for coding the Gibbs sampler correctly.
• 1 mark for posterior means, standard deviation, credible intervals.
• 1 mark for interpreting credible intervals correctly.
• 2 marks for checking convergence.







         
加QQ：99515681  WX：codinghelp  Email: 99515681@qq.com
