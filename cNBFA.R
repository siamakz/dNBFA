cNBFA <- function(counts, X, K, Burnin = 1000L, Collections = 1000L, PGTruncation = 10L, randtry = 2017)
{
    set.seed(randtry)
    y <- as.matrix(counts)
    idx.nz <- rowSums(y)!=0
    y <- y[idx.nz,]
    V <- dim(y)[1]       # V
    J <- dim(y)[2]     # J
    # hyperparameters
    a0 <- b0 <- 1
    e0 <- f0 <- 0.01
    
    P <- dim(X)[1]
    Beta <- matrix(0, P, K)
    Phi <- matrix(1/V, V, K)
    r <- rep(1/K, K)
    c0 <- 1
    gamma0 <- 1
    alpha <- rep(1, P)
    Theta <- matrix(rgamma(K*J,shape = r, rate = exp(t(Beta) %*% X)),K,J)
    eta <- 0.01
    realmin <- .Machine$double.xmin
    logLike <- -Inf
    out <- list(Phi=Phi,Beta=Beta,r=r)
    iterMax <- Burnin+Collections
    for (iter in 1:iterMax)
    {
        cat(iter, '\n')
        
        # Sample ell
        ter <- CRT_MultR(y,Phi,Theta)
        lkj <- ter$lkj
        lvk <- ter$lvk
        
        # Sample Phi
        Phi <- dirrnd(eta+lvk)
        
        # sample pj
        pj <- rbeta(J,a0+colSums(y),b0+colSums(Theta))
        qj <- -log(pmax(1-pj,realmin))
        Psi <- t(Beta) %*% X + matrix(log(pmax(qj,realmin)),nrow = K,ncol = J,byrow = T)
        
        # sample Theta
        Theta <- matrix(rgamma(K*J,r+lkj,t(qj*(1+exp(-t(Psi))))),K,J)
        
        # sample r_k
        lk_tild <- CRT_sum_matrix(lkj,r)
        temp <- rowSums(logOnePlusExp(Psi))
        r <- rgamma(K,gamma0/K+lk_tild,c0+temp)
        
        # Sample alpha
        alpha <- rgamma(P, 1 + K/2, rate = 1 + 0.5*rowSums(Beta^2))
        
        # Sample Beta
        for (k in 1:K){
            omega <- PolyaGamRnd_Gam(lkj[k,]+r[k], Psi[k,], Truncation = PGTruncation)
            sigmak <- X %*% diag(omega) %*% t(X)
            diag(sigmak) <- diag(sigmak) + pmax(alpha,1e-3)
            errMSG <- sum(eigen(sigmak)$values<=0)
            if (errMSG == 0)
                invchol <- chol(sigmak)
            else
            {
                count_0 <- 0
                while (errMSG != 0)
                {
                    diag(sigmak) <- diag(sigmak) + 10^(count_0-6)
                    errMSG <- sum(eigen(sigmak)$values<=0)
                    count_0 <- count_0 + 1
                }
                invchol <- chol(sigmak)
            }
            invchol <- solve(invchol)
            muk <- X %*% (0.5*(lkj[k,]-r[k])-omega*log(pmax(qj,realmin)))
            Beta[,k] <- invchol %*% (rnorm(P)+t(invchol) %*% muk)
        }
        Psi <- t(Beta) %*% X + matrix(log(pmax(qj,realmin)),nrow = K,ncol = J,byrow = T)
        
        # sample gamma0
        temp <- rowSums(logOnePlusExp(Psi))
        ell <- CRT_sum(lk_tild,gamma0/K)
        one_minus_p_tilde_k = c0/(c0+temp)
        gamma0 <- rgamma(1,e0+ell,f0-sum(log(pmax(one_minus_p_tilde_k,realmin)))/K)
        
        # sample c0
        c0 <- rgamma(1,1+gamma0,1+sum(r))
        
        # sample eta
        q_k <- rbeta(K,colSums(lvk),eta*V)
        Lv <- CRT_sum(lvk[lvk!=0],eta)
        eta <- rgamma(1,0.01+Lv,0.01-V*sum(log(pmax(1-q_k,realmin))))
        
        PhiTheta <- Phi %*% Theta
        temp_logLike <- sum(lgamma(y+PhiTheta)-lgamma(PhiTheta)) +
            sum(colSums(y)*log(pmax(realmin,pj)) + colSums(Theta)*log(pmax(realmin,1-pj)))
        if (iter>Burnin && temp_logLike>logLike)
        {
            logLike <- temp_logLike
            out$Phi <- Phi
            out$Beta <- Beta
            out$r <- r
        }
    }
    
    return(out)
}
