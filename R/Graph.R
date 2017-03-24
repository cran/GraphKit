#' @title Graph Invar
#' @name graph_invar
#' @description  Constructs confidence intervals for graph invariants 
#' of underlying graphical models
#' 
#' @param x A sample of data
#' @param sigmaHat A sample covariance matrix for \code{x}
#' @param thetaHat A sample precision matrix estimated from \code{sigmaHat}
#' @param invar The monotone graph invariant to examine
#' @param numB The number of bootstrap samples to take Default: 1000
#' @param alpha The significance level for the confidence interval Default: 0.05
#' @return A confidence intveral for the value of \code{invar} with sig. level \code{alpha}
#' @examples
#' 
#' data(Xs,cov.hat,t.hat)
#' graph_invar(x=Xs,sigmaHat=cov.hat,thetaHat=t.hat,invar="conn")
#' 
#' @export 
graph_invar <- function(x, sigmaHat, thetaHat, invar = c("connectivity", "longest_chain", "max_degree", "largest_clique", "chromatic_number", "num_singletons", 
    "girth"), numB = 1000, alpha = 0.05) {
    invar <- match.arg(invar)
    deOmega <- debias(sigmaHat, thetaHat)
    boot <- bootstrap(thetaHat, x, numB)

    if (invar == "connectivity") {
        return(skipDownConnCI(deOmega, boot, alpha))
    }
    
    if (invar == "longest_chain") {
        warning("Invariant not yet supported", call. = FALSE)
        return(skipDownChainCI(deOmega, boot, alpha))
    }
    
    if (invar == "max_degree") {
        return(skipDownDegCI(deOmega, boot, alpha))
    }
    
    if (invar == "largest_clique") {
        warning("Invariant not yet supported", call. = FALSE)
        return(skipDownCliqueCI(deOmega, boot, alpha))
    }
    
    if (invar == "chromatic_number") {
        warning("Invariant not yet supported", call. = FALSE)
        return(skipDownChromCI(deOmega, boot, alpha))
    }
    
    if (invar == "num_singletons") {
        return(skipDownSingleCI(deOmega, boot, alpha))
    }
    
    if (invar == "girth") {
        warning("Invariant not yet supported", call. = FALSE)
        return(skipDownGirthCI(deOmega, boot, alpha))
    }
}

#' @title Graph Prop
#' @name graph_prop
#' @description Conducts hypothesis testing on the graph properties 
#' of underlying graphical models
#' 
#' @param x A sample of data
#' @param sigmaHat A sample covariance matrix for \code{x}
#' @param thetaHat A sample precision matrix estimated from \code{sigmaHat}
#' @param prop The monotone graph property to examine
#' @param numB The number of bootstrap samples to take Default: 1000
#' @param alpha The significance level for the property test Default: 0.05
#' @param k The value to threshold numeric tests at
#' @return The value of \code{prop} with sig. level \code{alpha} 
#' @examples
#' 
#' data(Xs,cov.hat,t.hat)
#' graph_prop(x=Xs,sigmaHat=cov.hat,thetaHat=t.hat,prop="conn",k=5)
#' 
#' @export
graph_prop <- function(x, sigmaHat, thetaHat, prop = c("connectivity", "longest_chain", "max_degree", "largest_clique", "chromatic_number", "num_singletons", 
    "girth", "matching", "planarity", "bipartite", "acyclic"), k = 1, numB = 1000, alpha = 0.05) {
    prop <- match.arg(prop)
    deOmega <- debias(sigmaHat, thetaHat)
    boot <- bootstrap(thetaHat, x, numB)
    
    if (prop == "connectivity") {
        return(skipDownConn(deOmega, boot, alpha, k))
    }
    
    if (prop == "longest_chain") {
        warning("Property not yet supported", call. = FALSE)
        return(skipDownChain(deOmega, boot, alpha, k))
    }
    
    if (prop == "max_degree") {
        return(skipDownDeg(deOmega, boot, alpha, k))
    }
    
    if (prop == "largest_clique") {
        warning("Property not yet supported", call. = FALSE)
        return(skipDownClique(deOmega, boot, alpha, k))
    }
    
    if (prop == "chromatic_number") {
        warning("Property not yet supported", call. = FALSE)
        return(skipDownChrom(deOmega, boot, alpha, k))
    }
    
    if (prop == "num_singletons") {
        return(skipDownSingle(deOmega, boot, alpha, k))
    }
    
    if (prop == "girth") {
        warning("Property not yet supported", call. = FALSE)
        return(skipDownGirth(deOmega, boot, alpha, k))
    }
    
    if (prop == "matching") {
        warning("Property not yet supported", call. = FALSE)
        return(skipDownMatch(deOmega, boot, alpha, k))
    }
    
    if (prop == "planarity") {
        warning("Property not yet supported", call. = FALSE)
        return(skipDownPlan(deOmega, boot, alpha, k))
    }
    
    if (prop == "bipartite") {
        warning("Property not yet supported", call. = FALSE)
        return(skipDownBipar(deOmega, boot, alpha, k))
    }
    
    if (prop == "acyclic") {
        return(skipDownCycle(deOmega, boot, alpha))
    }
}
