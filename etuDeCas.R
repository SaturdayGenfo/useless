#' ---
#' title: "Etude de cas"
#' author: Leello Tadesse Dadi
#' ---
#'\newcommand{\esp}{\mathbb{E}}
#'\def\prob{\mathbb{P}}

source("functions_ETC2018.R")

alpha = 1.7
beta = 5.1

#' # Génération des données et fonctions préliminaires
#' 
#' ### Question 1
densite = function (x)
{
  return(dbeta(x, alpha, beta))
}
#' ### Question 2
#' D'apres la loi des grands nombres, la limite presque sure est
#' $$
#' \begin{aligned}
#' \hat{\phi}^{\infty}_{1.7, 5.1}(x) &= 2^8 \text{ }\prob_{B_{1.7, 5.1}}([\text{ }2^{-8} k(x), 2^{-8} (k(x) + 1) \text{ }])\\
#' &= 2^8 \int_{2^{-8} k(x)}^{2^{-8} (k(x) + 1)}\phi_{1.7, 5.1}(t)dt
#' \end{aligned}
#' $$
#' 
#' ### Question 3
#' 
#' 
#' # Interpolation polynomiale
#' ## Interpolation à partir de l’expression exacte de $\phi_{1.7, 5.1}$
#' 
#' ### Question 1
a = 2^(-8)
b = 1 - a
FUN = densite
neval = 1000
absc = a + (b-a)/(neval-1) * (0:(neval-1))
colors <- rainbow(60)
plot(absc, densite(absc), type='n')
lines(absc, densite(absc), col="red")
for (i in seq(2, 60, 3))
{
  lines(absc, interpolLagrange(i, a, b, 1000, "equi", FUN, FALSE), col=colors[i], lty="dotted", lwd=0.5)
}

