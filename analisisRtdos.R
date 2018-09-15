# Setear working directory
setwd("~/Documents/SimulacionesART/Resultados/")
setwd("~/ART/Resultados/")

######################################
# Cargar paquetes
######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

######################################
# Cargar funciones
######################################

source("../analisisRtdos_funciones.R")

######################################
# Leer los datos
######################################

# Cargar la primera corrida - unica vez
# corridas <- read.table("corrida1.txt", header = T, sep = "\t", stringsAsFactors = F)
# corridas <- mutate(corridas, nrOrden = NULL, corrida = 1)
# write.table(corridas, "corridas.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# Agregar corridas al archivo general
# agregarCorrida(c(2:36), path = "./")
# agregarCorrida(c(37:48), path = "./")

# Leer archivo de corridas
corridas <- read.table("corridas.txt", header = T, sep = "\t", stringsAsFactors = F)

# Armar un archivo ficitio de corridas para ir preparando los graficos
# colnames(corridas)
# load("../escenarios.Rdata")
# colnames(escenarios)
# pvalues <- matrix(runif(nrow(escenarios) * 45), ncol = 45)
# colnames(pvalues) <- colnames(corridas)[9:53]
# corridas <- cbind(escenarios, pvalues)
# corridas$corrida <- 1:nrow(corridas)
# corridas$desvios <- sapply(corridas$desvios, function(x) paste(x, collapse = "-"))

# Pasar a formato largo
corridasLong <- 
    corridas %>%
    select(-id, -corrida) %>%
    gather(key = pvalueKey, value = prob, pA.MC.0.01:pAB.Rob.0.1) %>%
    extract(pvalueKey, c("efecto", "metodo", "alfa"), "p(.+?)\\.(.+?)\\.(.+)") %>%
    mutate(efecto = factor(efecto, levels = c("A", "B", "AB"), labels = paste("Efecto", c("A", "B", "AB"))),
           modelos = factor(modelos, levels = paste0("mod", 1:6), labels = paste("Modelo", 1:6)),
           desvios = as.character(sapply(strsplit(desvios, "-"), length))) # paso desvios a indicar cuantas variancias distintas hay
head(corridasLong)

#########################################################################################
# Graficar Potencias y P(E1) estimadas
#########################################################################################

# Ejemplo un solo grafico
p <- graficarProbs(corridasLong, efTam = 0.8, disErr = "normal", desv = "1", propOut = 0, logaritmo = F, 
                   nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
plot(p)

# Realizar todos los graficos que interesen, ambos en escala logaritmica o comun
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.8, disErr = "normal", desv = "1", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.8, disErr = "lognormal", desv = "1", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.8, disErr = "normal", desv = "2", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.8, disErr = "normal", desv = "1", propOut = 0.3, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.8, disErr = "normal", desv = "6", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.2, disErr = "normal", desv = "1", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.2, disErr = "lognormal", desv = "1", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 1.5, disErr = "normal", desv = "1", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 1.5, disErr = "lognormal", desv = "1", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.8, disErr = "lognormal", desv = "1", propOut = 0, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.8, disErr = "normal", desv = "1", propOut = 0.1, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
for (logs in c(F, T)) graficarProbs(corridasLong, efTam = 0.8, disErr = "normal", desv = "1", propOut = 0.3, logaritmo = logs, 
                                    nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")

#########################################################################################
# Graficar Potencias y P(E1) estimadas para un solo metodo, comparando otros criterios
#########################################################################################

mets <- c("MC", "ART", "PS", "VDW", "Rob")

# Distribucion de los errores y variancias
for (met in mets) {
    for (logs in c(F, T)) {
        graficarProbs_PorMet_DistErr_Desv(corridasLong, met, efTam = 0.8, propOut = 0, logaritmo = logs, 
                                          nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
    }
}

# Distribucion de los errores y outliers
for (met in mets) {
    for (logs in c(F, T)) {
        graficarProbs_PorMet_DistErr_PropOut(corridasLong, met, efTam = 0.8, desv = "2", logaritmo = logs, 
                                             nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
    }
}

# Outliers y variancias
for (met in mets) {
    for (logs in c(F, T)) {
        graficarProbs_PorMet_PropOut_Var(corridasLong, met, efTam = 0.8, distr = "normal", logaritmo = logs, 
                                          nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
    }
}

# Distribucion de los errores y efecto tamaño
for (met in mets) {
    for (logs in c(F, T)) {
        graficarProbs_PorMet_DistErr_EfTam(corridasLong, met, propOut = 0, desv = "1", logaritmo = logs, 
                                             nivelSig = "0.05", titulo = "default", guardarGraf = T, path = "./Graficos/")
    }
}

#########################################################################################
# Graficar Potencias y P(E1) estimadas para un solo metodo, comparando otros criterios
#########################################################################################

p <- graficarPvalues(n = 10, efTam = 0.8, disErr = "normal", desv = 1, propOut = 0, pathArchivos = "./",
                     logaritmo = F, titulo = "default", guardarGraf = F, fileName = "default", path = "./Graficos/")
p

for (logs in c(F, T)) {
    for (n in c(3:10, 20)) {
        graficarPvalues(n, efTam = 0.8, disErr = "normal", desv = 1, propOut = 0, pathArchivos = "./",
                        logaritmo = logs, titulo = "default", guardarGraf = T, fileName = "default", path = "./Graficos/")
    }
}

