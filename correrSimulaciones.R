# Correr simulaciones
###################################################

# Leer opciones para la corrida
Args <- commandArgs()
print("Nro de replicas")
(repl <- as.numeric(strsplit(Args[6], "-")[[1]]))
print("Efecto tamano")
(efTam <- as.numeric(strsplit(Args[7], "-")[[1]]))
print("Nro de desvios distintos")
(largoDesv <- as.numeric(strsplit(Args[8], "-")[[1]]))
print("Prop de outliers")
(propOut <- as.numeric(strsplit(Args[9], "-")[[1]]))
print("Distr de errores")
(disErr <- strsplit(Args[10], "-")[[1]])
print("Carpeta de rtdos")
(pathOutput <- Args[11])
print("Nro de corrida")
(corrida <- as.numeric(Args[12]))
print("Nro de cores")
(cores <- as.numeric(Args[13]))
print("Nro de sims")
(nsim <- as.numeric(Args[14]))

# repl <- 9
# efTam <- 0.8
# largoDesv <- "1"
# propOut <- 0.3
# disErr <- "normal"
# pathOutput <- "/home/estadistica/Documents/SimulacionesART/"
# corrida <- 13

# Cargar paquetes
# library(doParallel)
library(doSNOW)
library(ARTool)
library(tidyr)
# library(ggplot2)
library(WRS2)
library(tictoc)

# Arrancar reloj!
tic()

# Paralelo
cl <- makeCluster(cores)
registerDoSNOW(cl)
progressBar <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(progressBar, n)
opts <- list(progress = progress)

# Cargar funciones
source(paste0(pathOutput, "/funcionesParaSimular.R"))

# Direccion para guardar los resultados
# pathOutput <- "~/Documents/SimulacionesART/"

# Cargar todos los escenarios creados
load(paste0(pathOutput, "escenarios.Rdata"))

# Elegir que escenarios correr
idx <- with(escenarios, which(replicas %in% repl &
                                  efectoTams %in% efTam &
                                  distrErrores %in% disErr &
                                  sapply(desvios, length) %in% largoDesv &
                                  propOutliers %in% propOut)) # todos los modelos para esta combinacion
escenarios <- escenarios[idx, ]
escenarios$nrOrden <- 1:nrow(escenarios)

# Correr simulaciones
salida <- apply(escenarios, 1, function(escenario) {

    cat("\nNumero de orden del escenario:", escenario$nrOrden, " / ", nrow(escenarios), "\n")
    
    # Definir la estructura de efectos en base al modelo y al efecto tamaño
    efectos <- definirEfectosTam(escenario$modelos, escenario$efectoTams)

    # Realizar simulacion
    resSimulacion <- simular(nsim = nsim, efectosTamA = efectos[[1]], efectosTamB = efectos[[2]], efectosTamAB = efectos[[3]],
                             replicas = escenario$replicas, distrError = escenario$distrErrores,
                             deError = escenario$desvios, propOutliers = escenario$propOutliers)

    # Guardar en un archivo los pvalues de cada iteracion
    write.table(resSimulacion$resultados, file = paste0(pathOutput, escenario$id, ".txt"), quote = F, sep = "\t", row.names = F)

    # Devolver los porcentajes de rechazo de H0 para anexar a la info del escenario y guardar como salida principal
    resSimulacion$propRechazos
})
salida <- cbind(escenarios, t(salida))

# Guardar los resultados de esta corrida. Incluye los id de los escenarios evaluados
salida$desvios <- sapply(salida$desvios, function(x) paste(x, collapse = "-"))# convertir lista a caracter
write.table(salida, file = paste0(pathOutput, "corrida", corrida, ".txt"), quote = F, sep = "\t", row.names = F)

close(progressBar)
stopCluster(cl)

toc()