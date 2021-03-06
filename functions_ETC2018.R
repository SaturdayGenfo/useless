
dividif=function(x,y){
  ##  Newton's Divided differences
  ## @param x: a vector containing the interpolation nodes 
  ## @param y: a vector of same size as x:
  ##           values of the interpolated function at the nodes
  ## @return : a vector of same size as x:
  ##          the divided differences
  ##          \eqn{f_[x_0, ... x_k]} of order 'length(x) -1'. 
  
  n = length(x) - 1 ## degree of the Lagrange polynomial
  d  = y 
  for (j in 2:(n+1) ) {
    # d[j : (n+1) ] = (d[j : (n+1)] - d[(j-1):n])/(-x[1: (n-j + 2)] + x[j:(n+1)])
    D=d
    for(i in j:(n+1))
    {
      
      d[i] = (D[i] - D[i-1])/(x[i] -x[i-j+1])
    }
  }
  return(d)    
}

X = c(0.1,0.2,0.3)
omun=function(x){x-X[1]}
omdeux=function(x){(x-X[1])*(x-X[2])}
a=c(11,4,3)
mypoly=function(x){a[1]+ a[2]*omun(x)+a[3]*omdeux(x)}
y=sapply(X,mypoly)
dividif(X,y)


hornerNewton = function(a,x,z){
  ## Horner's method: Evaluates  a polynom P at points z, given
  ## nodes x and the coefficients a of P in Newton's basis
  ##a + (b-a)/(n-1) * (0:(n-1))
  ## @param a : vector: the  coefficients of the polynomial in
  ##           Newton's basis
  ## @param x : the interpolation nodes. 
  ## @param z : vector of points where the polynom needs to be
  ##            evaluated. 
  ## @return  : a vector of same size as z: the value of the
  ##            polynomial at points z.
  ## 
  n = length(x);
  f  = 0
  for( i in 0:(n-1)){
    f = a[n-i] + f*(z-x[n-i])
  }
  return(f)
}


interpolDividif=function(x,y,z){
  ## Efficient Lagrange interpolation using Horner's method with  
  ## Newton basis for evaluation
  ## @param x : vector containing the interpolation nodes 
  ## @param y : vector of same size as x: values of the interpolated
  ##            function at the nodes
  ## @param z : vector of points where the  interpolating polynomial
  ##            needs to be evaluated. 
  ## @return  : vector of same size as z: the value of the
  ##            interpolating polynomial at points z.
  d = dividif(x, y);
  vals = hornerNewton(d,x,z);
  return(vals);
  ## Complete the code 
}




interpolLagrange =function(n, a, b, neval, nodes = 'equi', FUN, Plot){
  ## Generic Lagrange a + (b-a)/(n-1) * (0:(n-1))interpolation, with equidistant or Chebyshev nodes. 
  ## @param n : the degree of the interpolating polynomial on each
  ## subinterval
  ## @param a : left end-point of the interval
  ## @param b : right end-point of the interval
  ## @param neval :number of evaluation points (a regular grid will be
  ## used on [a,b]
  ## @param nodes :string, either "equi" (default) for equidistant
  ## Lagrange interpolation (on each subinterval) or "cheby" for
  ## using Chebyshev nodes.
  ## @param FUN: the function to be interpolated 
  ## @param Plot : logical. Setting 'Plot' to TRUE produces a plot
  ## showing the graph of
  ## the true functions and its interpolation.  
  ## @return : vector of size 'neval': the values of the Lagrange
  ## polynomial on an equi-distant grid.
  if (nodes == "equi"){
    x = a + (b-a)/n * (0:n)
  }
  else if (nodes == "cheby"){
    x = (a+b)/2 + (b-a)/2 * cos((1/2 + (0:n))* 3.14159265/(n+1))
  }
  else{stop("the nodes must be either 'equi' or 'cheby'") }
  
  ##
  ## Complete the code: compute a vector 'f' containing
  ## the interpolated  values on an equidistant
  ## evaluation grid 'z'. 
  ##
  ##
  
  z = a + (b-a)/(neval-1) * (0:(neval-1))
  
  
  
  y = sapply(x, FUN)
  
  
  f = interpolDividif(x, y, z)
  
  if( Plot ){
    if (nodes == "equi"){ methodName = " equidistant "}
    else {   methodName = " Chebyshev "}
    
    plot(z, sapply(z,FUN), type="l", ylim=range(c(y,f)) )
    title(main = paste("Lagrange interpolation with ",
                       toString(n+1), methodName,
                       " nodes", sep=""))
    lines(z,f, col = 'blue') 
    
    legend('topright', legend=c('function','interpolation'),
           col = c('black','red'), lwd=1)
    
  }
  
  return(f)              
}

piecewiseInterpol=function(n,nInt,a,b,neval, nodes = "equi", FUN, Plot){
  ## @param n : the degree of the interpolating polynomial on each
  ## subinterval
  ## @param nInt :  the number of sub-intervals
  ## @param a, b : endpoints of the interval
  ## @param neval : the number of points on the interpolating grid (on
  ## each subinterval)
  ## @param nodes : string, either "equi" (default) for equidistant
  ## Lagrange interpolation (on each subinterval) or "cheby" for
  ## chebyshev nodes.
  ## @param FUN the function to be interpolated
  ## @param Plot : logical. Should the result be plotted ?
  ## @return : a matrix with 2 rows and neval * nInt -neval + 1:
  ## values of the interpolated funtion on a regular grid (first row)
  ## and the corresponding abscissas (second row).
  
  intEndPoints = seq(a,b,length.out = nInt+1)
  f = c()
  z = c()
  for (m in 1:nInt){
    A = intEndPoints[m]; B = intEndPoints[m+1] 
    
    fm = interpolLagrange(n ,A, B, neval, nodes, FUN, FALSE) 
    zm = A + (B-A)/(neval-1) * (0:(neval-1))
    
    if( m >= 2 && nodes == "equi"){
      fm = fm[-1]
      zm = zm[-1]
      
    }
    z = c(z,zm)
    f = c(f,fm)
  }
  
  if (Plot == 1){
    if (nodes == "equi") {methodName = " equidistant "}
    else  {methodName = " Chebyshev "}
    
    
    plot(z, sapply(z,FUN),type="l")
    title(main = paste("Piecewise  Lagrange  interpolation with ",
                       toString(n+1), methodName, " nodes  on ",
                       toString(nInt), " Intervals", sep=""))
    lines(z,f, col='red', lwd=2)
    legend('topright', legend = c('function','interpolation'),
           lwd=c(1,2), col=c('black','red'))
  }
  return(rbind(f,z) )
}

##################################
######## Exercice 1 ##############
######## Méthode des trapezes ####

trapezeInt =function(FUN,a,b,M){
  ##' TRAPEZOIDAL INTEGRATION RULE (COMPOSITE)
  ##' @param FUN : the function to be integrated
  ##' @param a, b : interval end points 
  ##' @param M : number of intervals (each of size (b-a)/M)
  ##' @return: the value of the composite trapezoidal quadrature. 
  x = seq(a,b, length.out= M+1)
  y = sapply(x, FUN)
  h = (b-a)/M
  q = sum(1/2*h*(y[1:M]+y[2:(M+1)]))
  return(q)
}

##################################
######## Exercice 2 ##############
######## Division du pas  des trapezes par deux ####


refineTrapeze=function(FUN,a,b,M,q){
  ##' refinement of the subdivision step: incremental method
  ##' @param FUN : the function to be integrated
  ##' @param a, b : interval end points 
  ##' @param M : initial number of intervals (each of size (b-a)/M)
  ##'  having been used to compute q
  ##' @param  q : the value of the trapezoidal  quadrature method
  ##'  of stepsize (b-a)/M
  ##' @return : the value of the quadrature for a stepsize h' = h/2
  h = (b-a)/M
  x =   seq(a + h/2,b -h/2, length.out= M)
  ##  x : a vector of size M :
  ##     the additional abscissas where 'fun' must be evaluated.
  y = sapply(x, FUN)
  Q = 1/2*q + 1/2*sum(h*y) 
  return(Q)
}






##################################
######## Exercice 4 ##############
######## Des trapèzes à Simpson ##
##################################


simpsonInt = function(FUN,a,b,M){
  ##' Simpson integration via trapeze rule
  ##' uses the fact that 
  ##' simpson(h) = 4/3(trapeze(h/2) - 1/4 trapeze(h))
  h = (b-a)/M;
  qtrapeze = trapezeInt(FUN, a, b, M)
  qrefined = refineTrapeze(FUN, a, b, M, qtrapeze)
  q =  4/3*(qrefined - 1/4 * qtrapeze)
  return(q)
}




##################################################
############# Exercice 5: evaluation de l'erreur a posteriori
evalErrSimpson=function(FUN,a,b,M){
  ## Computes an approximation E of the error 
  ## for the composite Simpson rule of step h=(b-a)/(2M). 
  ##This requires computing I_M and I_{2M}. 
  ##The value  q = I_{2M} is also returned. 
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h =   4/3*(qth2 - 1/4* qth ) 
  simps_h2 =  4/3*(qth4 - 1/4 * qth2)
  q = simps_h2
  N = 3
  I = 1/(2**(N+1) - 1) * (2**(N+1) *q - simps_h)
  E = (simps_h - q)/15
  return(c(E,q))
}






########## Exercice 6: Richardson
############################

richardson = function(FUN,n,t,delta){
  ## Calcule le tableau des differences  divisees en 0 du 
  ## polynome d'interpolation en t,delta t, ... delta^n t
  ## renvoie un vecteur de taille n+1:
  ## le vecteur des A_{k,k}, k= 0 .. n 
  ## (pas la matrice).   
  ## La meilleure approximation est le dernier element A[n+1].
  ##
  lx = log(t)  +  log(delta) *(0:n)
  x = exp(lx) 
  A = sapply(x,FUN) 
  for( j in 2:(n+1)){
    A[j : (n+1) ] = 1/(1-delta**(j-1)) * (A[j:(n+1)] - delta**(j-1) * A[(j-1):n]) ## Completer le code 
  }
  return(A)
}


##########################
######## exercice 7:  Romberg


romberg =function(FUN,n,a,b,M){## methode de Romberg avec n etapes
  ## appliquee sur la fonction FUN sur l'intervalle (a,b), avec un
  ## pas initial h = (b-a)/M
  h= (b-a)/M 
  A = rep(0, n+1)
  A[1] = trapezeInt(FUN,a,b,M);
  Mc = M
  ## initialisation des differences divisees
  for( i in 2:(n+1)){
    A[i] = refineTrapeze( FUN,a,b, Mc, q= A[i-1])
    Mc = 2*Mc 
  }
  delta = 1/4;
  for (j in 2:(n+1)){
    A[j : (n+1) ] = 1/(1-delta**(j-1)) * (A[j:(n+1)] - delta**(j-1) * A[(j-1):n])
  }
  return(A)
}

