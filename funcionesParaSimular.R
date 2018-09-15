simular <- function(nsim = 5000,
                    
                    # Diseño
                    nivelesA = 2, 
                    nivelesB = 3, 
                    replicas = 3, 
                    
                    # Valores verdaderos de los coeficientes relativos a la variancia
                    mu = 0,
                    efectosTamA = rep(0, nivelesA), # efectosTamA = c(efectoTam, 0)
                    efectosTamB = rep(0, nivelesB), # efectosTamB = c(efectoTam, efectoTam, 0)
                    efectosTamAB = rep(0, nivelesA * nivelesB), # efectosTamAB = c(efectoTam, 0, efectoTam, 0, efectoTam, 0)
                    
                    # Distribución y variancia de los errores
                    distrError = "normal", 
                    deError = 1, # deError = c(1, 1, 1, 2, 2, 2)
                    
                    # Outliers
                    propOutliers = 0,
                    tamOutlier = 8
) {
    
    # Diseño
    n <- nivelesA * nivelesB * replicas
    nparam <- 1 + (nivelesA-1) + (nivelesB-1) + (nivelesA-1)*(nivelesB-1)
    d <- data.frame(A = gl(nivelesA, nivelesB * replicas, n), B = gl(nivelesB, replicas, n))
    
    # Si es un único valor de variancia, repetirlo para generalizar el código
    if (length(deError) == 1) deError <- rep(deError, nivelesA * nivelesB)
    
    ## Generar promedios verdaderos de celda en función del efecto tamaño y de las variancias
    # Obtener el efecto tamaño para cada celda combinando efectosTamA, efectosTamB y efectosTamAB:
    efectoTamCelda <- rep(efectosTamA, each = nivelesB) + rep(efectosTamB, nivelesA) + efectosTamAB
    # Generar los promedios por celda
    promCelda <- mu + efectoTamCelda * deError
    
    # Correr simulaciones en paralelo
    res <- foreach(i = seq(nsim), .combine = rbind, .export = c("metodoMC", "metodoART", "metodoPS", "metodoVDW", "metodoRobusto", "np.anova"), .options.snow = opts) %dopar% {
    # for (i in seq(nsim)) {
        
        # cat("Iteracion ", i, "\n")
        
        # Simular los errores aleatorios
        if (distrError == "normal") {           # Distribución normal para los errores
            errores <- sapply(deError, function(x) rnorm(replicas, mean = 0, sd = x))
            # Si es variancia homogenea, deError tiene 6 veces el mismo valor, es lo mismo que rnorm(n, mean = 0, sd = deError) 
            # Si las var son heterog, ya chequee que da bien (pero se necesitan mas replicas para verlo, por ej 100): apply(errores, 2, sd)
            # errores es una matriz de dim replicas*celdas, cada columna fue generada con cada valor en deError, la primera columna tiene los errores para los datos de la celda A1B1, luego A1B2, etc, etc
        } else if (distrError == "lognormal") {          # Otra distribución para los errores
            errores <- sapply(deError, function(x) rlnorm(replicas, meanlog = 0, sdlog = x))
        }
        
        # Generar las observaciones
        y <- rep(promCelda, each = replicas) + as.vector(errores)
        
        # promCelda
        # tapply(y, paste0(d$A, d$B), mean)
        # boxplot(y ~ paste0(d$A, d$B))
        # stripchart(y ~ paste0(d$A, d$B), vertical = TRUE, method = "jitter", pch = 19, add = TRUE)
        # qqnorm(y)
        
        # Añadir outliers
        if (propOutliers > 0) {
            idx <- sample(1:n, propOutliers * n)
            y[idx] <- y[idx] + sign(y[idx] - rep(promCelda, each = replicas)[idx]) * tamOutlier * rep(deError, each = replicas)[idx]
        }
        
        # Técnicas de análisis bifactorial
        MC <- metodoMC(y, d)
        ART <- metodoART(y, d)
        PS <- metodoPS(y, d)
        VDW <- metodoVDW(y, d)
        Rob <- metodoRobusto(y, d)
        
        # Elementos a devolver
        c(MC, ART, PS, VDW, Rob)
    }
    row.names(res) <- NULL
    metodos <- c("MC", "ART", "PS", "VDW", "Rob")
    colnames(res) <- paste0(rep(c("pA.", "pB.", "pAB."), length(metodos)), rep(metodos, each = 3))
    
    # Porcentaje de casos significativos
    propRechazos <- sapply(c(0.01, 0.05, 0.1), function(alpha) {
        apply(res, 2, function(col) mean(col < alpha))
    })
    propRechazos <- as.vector(propRechazos)
    names(propRechazos) <- paste0(rep(colnames(res), 3), rep(paste0("<", c(0.01, 0.05, 0.1)), each = ncol(res)))
    
    # Elementos a devolver
    output <- list(nSim = nsim, nivelesA = nivelesA, nivelesB = nivelesB, replicas = replicas, mu = mu, 
                   efectosTamA = efectosTamA, efectosTamB = efectosTamB, efectosTamAB = efectosTamAB, 
                   distrError = distrError, deError = deError, propOutliers = propOutliers, tamOutlier = tamOutlier,
                   resultados = res, propRechazos = propRechazos)
    
    return(output)
}

metodoMC <- function(y, d) {
    modelo <- lm(y ~ A * B, d)
    pvalues <- anova(modelo)$`Pr(>F)`[1:3]
    return(pvalues)
}


metodoART <- function(y, d) {
    modelo <- ARTool::art(y ~ A * B, d)
    pvalues <- anova(modelo)$`Pr(>F)`[1:3]
    return(pvalues)
}

metodoPS <- function(y, d) {
    d$y <- y
    modelo <- np.anova(y ~ A * B, d, 0)
    pvalues <- modelo$` Pr(>Chi)`[1:3]
    return(pvalues)
}

metodoVDW <- function(y, d) {
    d$y <- y
    modelo <- np.anova(y ~ A * B, d, 1)
    pvalues <- modelo$` Pr(>Chi)`[1:3]
    return(pvalues)
}

metodoRobusto <- function(y, d) {
    d$y <- y
    modelo <- try(WRS2::pbad2way(y ~ A * B, d, est = "mom", pro.dis = FALSE), silent = T)
    if (inherits(modelo, "try-error")) {
        modelo <- WRS2::pbad2way(y ~ A * B, d, est = "mom", pro.dis = TRUE)
    }
    pvalues <- c(modelo$A.p.value, modelo$B.p.value, modelo$AB.p.value)
    return(pvalues)
}

np.anova <- function(formel, data, method = 0, compact = T) {   
    # data  :  dataframe mit 3 Gruppierungsfaktoren und
    # formel:  anova-Modell
    # method:  0: generalized Kruskal-Wallis-Friedman (Puri & Sen) tests
    #             including Iman & Davenport F-tests
    #          1: generalized van der Waerden tests
    # compact: T: all tests in one table
    #          F: seperate tables for each error term (repeated measures only)
    
    # check formula
    if (mode(formel)!="call") stop("invalid formula")
    # perform aov and check model
    aov.1 <- aov(formel,data)
    if (class(aov.1)[1]=="aov")     rep=F
    if (class(aov.1)[1]=="aovlist") rep=T
    
    # extract varnames
    vars  <- all.vars(formel)
    vdep  <- vars[1]
    nfact <- length(vars)-1
    if (rep) nfact <- length(vars)-2
    nfgrp <- nfact
    # ranking procedure
    # repeated measures
    if (rep) {
        
        vpn      <- vars[length(vars)]
        nerrterm <- length(names(aov.1))-2
        nfrep    <- log2(length(names(aov.1))-1)
        nfgrp    <- nfact-nfrep
        cfrep    <- unlist(strsplit(names(aov.1),":")[3:(nfrep+2)])[seq(2,2*nfrep,2)]
        nrep     <-1
        for (i in 1:nfrep) nrep <- nrep*length(levels(data[,cfrep[i]]))
        # rank Vpn
        ysum <- ave(data[,vdep],data[,vpn],FUN=sum)
        rsum     <-(rank(ysum)+(nrep-1)/2)/nrep
        #   normal scores for van der Waerden tests
        if (method==1)
        {  ncases <- sum(summary(aov.1[[2]])[[1]]$"Df")+1
        rsum   <- qnorm(rsum/(ncases+1)) }
        #   rank repeated
        zzz_ry <- ave(data[,vdep], data[,vpn], FUN=rank) 
        if (method==1) zzz_ry <- qnorm(zzz_ry/(nrep+1))
        zzz_rx <- zzz_ry
        
        if (nfgrp>0)  
        {if (method==0)      zzz_rx <- (rsum-1)*nrep + zzz_ry
        else if (method==1) zzz_rx <- rsum+zzz_ry }
        #   get grouping factors
        if (nfgrp>1)
        {ig <- 0 ; ifgrp <- 0
        for (i in 1:nfact)
        {if (length(grep(vars[i+1],names(aov.1)))==0)
        { ig <- ig+1 ; ifgrp[ig] <- vars[i+1] }
        }
        }
    }
    else
    { ncases <- sum(anova(aov.1)[,1])+1
    zzz_rx     <- rank(data[,vdep])
    # normal scores for van der Waerden tests
    if (method==1) zzz_rx <- qnorm(zzz_rx/(ncases+1))
    }
    # perform npar anova
    cformel <- as.character(formel)
    formel1 <- as.formula(paste("zzz_rx ~",cformel[3]))
    aov.2   <- aov(formel1,data)
    # perform puri & sen-tests
    # repeated measures
    if (rep) {
        # vpn error part (including grouping effects)
        aov.21 <- summary(aov.2[[2]])[[1]]
        df     <- aov.21$"Df"
        ncases <- sum(df)+1
        nerg   <- 0
        ngeffects <- 0
        if (nfgrp>0)                                                                                          
        {ngeffects <- dim(aov.21)[1]-1
        mstotal   <- sum(aov.21$"Sum Sq")/sum(df)
        # also grouping effects
        if (nfgrp>1)
        {                                        # type III ssq necessary
            # separate anova for grouping effects
            formel2 <- paste(c("zzz_rx",paste(ifgrp,collapse="*")),collapse="~")
            formel3 <- paste(c("zzz_rx",paste(c(vpn,ifgrp),collapse="*")),collapse="~")
            data1   <- aggregate(as.formula(formel3),data=data,FUN=mean)
            aov.31  <- drop1(aov(as.formula(formel2),data1), ~. , test="F")
            aov.21$"Sum Sq"[1:(ngeffects-1)] <- nrep*aov.31$"Sum of Sq"[2:ngeffects]
        }
        # chisquare tests
        names(aov.21)[5]   <- " "
        names(aov.21)[3:4] <-c(" Chi Sq"," Pr(>Chi)")
        chisq       <- aov.21$"Sum Sq"/mstotal
        aov.21[[3]] <- chisq
        aov.21[[4]] <- 1-pchisq(chisq,df)
        aov.21[ngeffects+1,3:4] <- NA
        aov.21[,5]  <- NA
        
        erglist     <- list(aov.21)
        names(erglist)[1] <- names(aov.1)[2]
        nerg        <- 1
        }
        for (errterm in 3:(2+nerrterm))
        {                                        # for each trial error term...
            aov.22   <- summary(aov.2[[errterm]])[[1]]
            neffects <- dim(aov.22)[1]-1
            if (method==1) names(aov.22)[5] <- c(" ")
            # cols for Iman-Davenport F
            if (method==0)
            {aov.22[[6]] <- aov.22[[5]]
            aov.22[[5]] <- aov.22[[4]]
            names(aov.22)[5:6] <- names(aov.22)[4:5] }
            names(aov.22)[3:4] <- c("Chisq","Pr(>Chi)")
            # chisquare tests
            df      <- aov.22$"Df"
            sumdf   <- sum(df)
            mstotal <- sum(aov.22$"Sum Sq")/sumdf
            chisq   <-aov.22$"Sum Sq"/mstotal
            aov.22[[3]] <- chisq
            aov.22[[4]] <- 1-pchisq(chisq,df)
            aov.22[neffects+1,3:4] <- NA
            aov.22[,5]  <- NA
            # Iman-Davenport F
            if (method==0)
            {aov.22[[5]] <- (ncases-1)*chisq/(sumdf-chisq)
            aov.22[[6]] <- 1-pf(aov.22[[5]],df,df[neffects+1]) 
            aov.22[neffects+1,5:6] <- NA }
            # type of tables
            if (compact)
            {if (errterm==3) 
            {aov.3     <- aov.22
            nteffects <- neffects+1
            row.names(aov.3)[nteffects] <- paste("Residuals",row.names(aov.3)[nteffects-ngeffects-1])
            if (nfgrp>0)
            {nteffects <- nteffects+ngeffects+1
            aov.3[(ngeffects+2):nteffects,] <- aov.22[1:(neffects+1),]
            aov.3[1:(ngeffects+1),]         <- aov.21[1:(ngeffects+1),]
            if (method==0) aov.3[1:(ngeffects+1),6] <- NA
            nd<- dim(aov.21)[1]
            row.names(aov.21)[nd]<-"Residuals Btw.Vpn"        
            row.names(aov.3)[1:(ngeffects+1)]         <- row.names(aov.21)
            row.names(aov.3)[(ngeffects+2):nteffects] <- row.names(aov.22) 
            row.names(aov.3)[nteffects] <- paste("Residuals",row.names(aov.3)[nteffects-ngeffects-1]) }
            }
                
                else
                {aov.3[(nteffects+1):(nteffects+neffects+1),] <- aov.22[1:(neffects+1),]
                nteffects <- nteffects+neffects+1
                row.names(aov.3)[nteffects] <- paste("Residuals",row.names(aov.3)[nteffects-1]) }
            }
            else
            {nerg <- nerg+1
            if (nerg==1) erglist         <- list(aov.22)
            if (nerg >1) erglist[[nerg]] <- aov.22 }
        }
        if (compact) 
            erglist <- aov.3
        
        else
        {if (nfgrp==0) names(erglist) <- names(aov.1)[3:(nerrterm+2)]
        if (nfgrp >0) names(erglist) <- names(aov.1)[2:(nerrterm+2)] }
        
    }
    else {
        # only grouping effects
        # SS Type III
        aov.3 <- drop1(aov.2, ~. , test="F")
        # as dataframe
        aov.2 <- anova(aov.2) 
        names(aov.2)[3:5] <- c(" Chi Sq"," Pr(>Chi)"," ")
        # Chisquare tests
        neffects <- dim(aov.2)[1]-1
        mstotal  <- sum(aov.2[,2])/sum(aov.2[,1])
        aov.2[1:(neffects),2] <- aov.3[2:(neffects+1),2]
        aov.2[1:neffects,3]     <- aov.3[2:(neffects+1),2]/mstotal
        aov.2[,4]               <- 1-pchisq(aov.2[,3],aov.2[,1])
        aov.2[neffects+1,3:4]   <- NA
        aov.2[,5]               <- NA
        erglist  <- aov.2
    }
    
    if (method==0) attr(erglist,"heading") <- "generalized Kruskal-Wallis/Friedman (Puri & Sen) tests including Iman & Davenport F-tests"
    if (method==1) attr(erglist,"heading") <- "generalized van der Waerden tests"
    erglist
}

definirEfectosTam <- function(modelo, efectoTam) {
    switch(modelo,
           mod1 = list(efectosTamA = rep(0, 2), efectosTamB = rep(0, 3), efectosTamAB = rep(0, 6)),
           mod2 = list(efectosTamA = c(efectoTam, 0), efectosTamB = rep(0, 3), efectosTamAB = rep(0, 6)),
           mod3 = list(efectosTamA = c(efectoTam, 0), efectosTamB = c(efectoTam, efectoTam, 0), efectosTamAB = rep(0, 6)),
           mod4 = list(efectosTamA = rep(0, 2), efectosTamB = rep(0, 3), efectosTamAB = rep(c(efectoTam, 0), 3)),
           mod5 = list(efectosTamA = c(efectoTam, 0), efectosTamB = rep(0, 3), efectosTamAB = rep(c(efectoTam, 0), 3)),
           mod6 = list(efectosTamA = c(efectoTam, 0), efectosTamB = c(efectoTam, efectoTam, 0), efectosTamAB = rep(c(efectoTam, 0), 3))
    )
} 


simularReturnDatos <- function(
                    
                    # Diseño
                    nivelesA = 2, 
                    nivelesB = 3, 
                    replicas = 3, 
                    
                    # Valores verdaderos de los coeficientes relativos a la variancia
                    mu = 0,
                    efectosTamA = rep(0, nivelesA), # efectosTamA = c(efectoTam, 0)
                    efectosTamB = rep(0, nivelesB), # efectosTamB = c(efectoTam, efectoTam, 0)
                    efectosTamAB = rep(0, nivelesA * nivelesB), # efectosTamAB = c(efectoTam, 0, efectoTam, 0, efectoTam, 0)
                    
                    # Distribución y variancia de los errores
                    distrError = "normal", 
                    deError = 1, # deError = c(1, 1, 1, 2, 2, 2)
                    
                    # Outliers
                    propOutliers = 0,
                    tamOutlier = 8
) {
  
  # Diseño
  n <- nivelesA * nivelesB * replicas
  nparam <- 1 + (nivelesA-1) + (nivelesB-1) + (nivelesA-1)*(nivelesB-1)
  d <- data.frame(A = gl(nivelesA, nivelesB * replicas, n), B = gl(nivelesB, replicas, n))
  
  # Si es un único valor de variancia, repetirlo para generalizar el código
  if (length(deError) == 1) deError <- rep(deError, nivelesA * nivelesB)
  
  ## Generar promedios verdaderos de celda en función del efecto tamaño y de las variancias
  # Obtener el efecto tamaño para cada celda combinando efectosTamA, efectosTamB y efectosTamAB:
  efectoTamCelda <- rep(efectosTamA, each = nivelesB) + rep(efectosTamB, nivelesA) + efectosTamAB
  # Generar los promedios por celda
  promCelda <- mu + efectoTamCelda * deError
  
  # Simular solo una
  
    # Simular los errores aleatorios
    if (distrError == "normal") {           # Distribución normal para los errores
      errores <- sapply(deError, function(x) rnorm(replicas, mean = 0, sd = x))
      # Si es variancia homogenea, deError tiene 6 veces el mismo valor, es lo mismo que rnorm(n, mean = 0, sd = deError) 
      # Si las var son heterog, ya chequee que da bien (pero se necesitan mas replicas para verlo, por ej 100): apply(errores, 2, sd)
      # errores es una matriz de dim replicas*celdas, cada columna fue generada con cada valor en deError, la primera columna tiene los errores para los datos de la celda A1B1, luego A1B2, etc, etc
    } else if (distrError == "lognormal") {          # Otra distribución para los errores
      errores <- sapply(deError, function(x) rlnorm(replicas, meanlog = 0, sdlog = x))
    }
    
    # Generar las observaciones
    d$y <- rep(promCelda, each = replicas) + as.vector(errores)
    return(d)
}
