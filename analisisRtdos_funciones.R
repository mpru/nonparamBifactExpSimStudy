######################################
# Funciones para el analisis de los resultados
######################################


agregarCorrida <- function(nroCorrida, path) {
    for (i in nroCorrida) {
        corrida <- read.table(paste0(path, "corrida", i, ".txt"), header = T, sep = "\t", stringsAsFactors = F)
        corrida <- mutate(corrida, nrOrden = NULL, corrida = i)
        write.table(corrida, paste0(path, "corridas.txt"), append = T, col.names = F, quote = F, sep = "\t", row.names = F)
    }
}

graficarProbs <- function(corridasLong, 
                          efTam, 
                          disErr, 
                          desv, # "1" "2" "6"
                          propOut, 
                          logaritmo = FALSE, 
                          nivelSig = "0.05", 
                          titulo = NULL, # o "default" que usa lo mismo que fileName 
                          guardarGraf = F, 
                          fileName = "default", 
                          path = "") {
    
    # Pasar como argumento un unico valor para cada parametro efTam, disErr, desv, propOut
    # Graficamos todos los n en el eje horizontal, las probs (potencias o prob E1) en el vertical
    # Los paneles estan separados por modelo y efecto, y una linea para cada metodo
    
    # Primero quedarse con los casos de interes
    corridasLongSub <- filter(corridasLong, 
                              efectoTams %in% efTam,
                              distrErrores %in% disErr, 
                              desvios %in% desv,
                              propOutliers %in% propOut,
                              alfa %in% nivelSig)
    
    # Puede que alguna prob sea cero, lo cual impide usar el log, sumarle algo
    if (logaritmo) corridasLongSub$prob[corridasLongSub$prob == 0] <- 0.0008
    
    # Data para saber en qué paneles debo graficar la linea de la P(E1)
    modelosEfectos <- data.frame(modelos = rep(1:6, 3),
                                 efecto = rep(c("A", "B", "AB"), each = 6),
                                 dibujarAlfa = c(T, F, F, T, F, F,
                                                 T, T, F, T, T, F,
                                                 T, T, T, F, F, F),
                                 paraLeyenda = "P(E1)obj")
    modelosEfectos <- mutate(modelosEfectos,
                             modelos = factor(modelos, levels = 1:6, labels = paste("Modelo", 1:6)),
                             efecto = factor(efecto, levels = c("A", "B", "AB"), labels = paste("Efecto", c("A", "B", "AB"))))
    
    # Graficar    
    p <- ggplot(corridasLongSub, aes(x = replicas, y = prob, group = metodo, color = metodo))
    p <- p + geom_line(lwd = 1)
    # p <- p + facet_grid(efecto ~ modelos)
    p <- p + scale_color_discrete("Método")
    p <- p + theme(panel.spacing = unit(.05, "lines"),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.25), 
                   strip.background = element_rect(color = "black", size = 0.25))
    p <- p + xlab("Número de réplicas") + ylab("Probabilidad Estimada")
    # p <- p + scale_x_continuous(breaks = c(3:10, 20))
    
    if (logaritmo) {
        alfaNum <- as.numeric(nivelSig)
        p <- p + facet_grid(efecto ~ modelos, scales = "free_y") # facets con eje vertical libre
        p <- p + coord_trans(y = "log")
        p <- p + scale_y_continuous(breaks = c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1))
        p <- p + geom_hline(data = modelosEfectos[modelosEfectos$dibujarAlfa, ], aes(yintercept = dibujarAlfa * alfaNum, linetype = paraLeyenda))
        p <- p + scale_linetype_manual("", labels = "P(E1)\nobjetivo", values = "dashed")
    } else {
        p <- p + facet_grid(efecto ~ modelos)
        p <- p + scale_y_continuous(breaks = seq(0.25, 1, 0.25))
    }
    # Agregar titulo
    plotTag <- paste0("efTam", efTam, "_", disErr, "_desv", desv, "_out", propOut, "_alfa", nivelSig)
    if (logaritmo) plotTag <- paste0(plotTag, "_log")
    if (!is.null(titulo)) {
        if (titulo == "default") titulo <- plotTag
        p <- p + ggtitle(titulo)
    }
    
    if (guardarGraf) {
        if (fileName == "default") fileName <- paste0(plotTag, ".png")
        ggsave(fileName, plot = p, device = "png", path = path, width = 15, height = 7.5, units = "in", dpi = 200)
    }
    invisible(p)
}


graficarProbs_PorMet_DistErr_Desv <- function(corridasLong, 
                                              met, 
                                              efTam, 
                                              propOut, 
                                              logaritmo = FALSE, 
                                              nivelSig = "0.05", 
                                              titulo = NULL, # o "default"
                                              guardarGraf = F, 
                                              fileName = "default", 
                                              path = "") {
    
    # Pasar como argumento un unico valor para cada parametro efTam, propOut
    
    # Primero quedarse con los casos de interes
    corridasLongSub <- filter(corridasLong, 
                              metodo == met,
                              efectoTams %in% efTam,
                              propOutliers %in% propOut,
                              alfa %in% nivelSig)
    
    # Puede que alguna prob sea cero, lo cual impide usar el log, sumarle algo
    if (logaritmo) corridasLongSub$prob[corridasLongSub$prob == 0] <- 0.0008
    
    # Data para saber en qué paneles debo graficar la linea de la P(E1)
    modelosEfectos <- data.frame(modelos = rep(1:6, 3),
                                 efecto = rep(c("A", "B", "AB"), each = 6),
                                 dibujarAlfa = c(T, F, F, T, F, F,
                                                 T, T, F, T, T, F,
                                                 T, T, T, F, F, F),
                                 paraLeyenda = "P(E1)obj")
    modelosEfectos <- mutate(modelosEfectos,
                             modelos = factor(modelos, levels = 1:6, labels = paste("Modelo", 1:6)),
                             efecto = factor(efecto, levels = c("A", "B", "AB"), labels = paste("Efecto", c("A", "B", "AB"))))
    
    # Graficar    
    p <- ggplot(corridasLongSub, aes(x = replicas, y = prob, group = paste0(distrErrores, desvios), linetype = distrErrores, color = desvios))
    p <- p + geom_line(lwd = 1)
    p <- p + theme(panel.spacing = unit(.05, "lines"),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.25), 
                   strip.background = element_rect(color = "black", size = 0.25),
                   legend.key.width = unit(2, "cm"))
    p <- p + xlab("Número de réplicas") + ylab("Probabilidad Estimada")
    p <- p + scale_color_discrete("Variancias",  labels = c("1" = "Homog", "2" = "Heterog\n2 distintas", "6" = "Heterog\n6 distintas"))
    
    if (logaritmo) {
        alfaNum <- as.numeric(nivelSig)
        p <- p + facet_grid(efecto ~ modelos, scales = "free_y") # facets con eje vertical libre
        p <- p + coord_trans(y = "log")
        p <- p + scale_y_continuous(breaks = c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1))
        p <- p + geom_hline(data = modelosEfectos[modelosEfectos$dibujarAlfa, ], aes(yintercept = dibujarAlfa * alfaNum, linetype = paraLeyenda), lwd = 0.5)
        p <- p + scale_linetype_manual("Distr. de\nlos errores",
                                       values = c("normal" = 1, "lognormal" = 4, "P(E1)obj" = 2),
                                       labels = c("normal" = "Normal", "lognormal" = "Lognormal", "P(E1)obj" = "P(E1)\nobjetivo"))
        p <- p + guides(linetype = guide_legend(override.aes = list(size = 1))) # definir grosor de linea en leyenda independiete del grafico
        
    } else {
        p <- p + scale_linetype_manual("Distr. de\nlos errores", 
                                       values = c("normal" = 1, "lognormal" = 4),
                                       labels = c("normal" = "Normal", "lognormal" = "Lognormal"))
        p <- p + facet_grid(efecto ~ modelos)
        p <- p + scale_y_continuous(breaks = seq(0.25, 1, 0.25))
    }
    
    # Agregar titulo
    plotTag <- paste0(met, "_DistErrYDesv_", "efTam", efTam, "_out", propOut, "_alfa", nivelSig)
    if (logaritmo) plotTag <- paste0(plotTag, "_log")
    if (!is.null(titulo)) {
        if (titulo == "default") titulo <- plotTag
        p <- p + ggtitle(titulo)
    }
    
    if (guardarGraf) {
        if (fileName == "default") fileName <- paste0(plotTag, ".png")
        ggsave(fileName, plot = p, device = "png", path = path, width = 15, height = 7.5, units = "in", dpi = 200)
    }
    invisible(p)
}


graficarProbs_PorMet_DistErr_PropOut <- function(corridasLong, 
                                                 met, 
                                                 efTam, 
                                                 desv, # "1" "2" "6"
                                                 logaritmo = FALSE, 
                                                 nivelSig = "0.05", 
                                                 titulo = NULL, # o "default"
                                                 guardarGraf = F, 
                                                 fileName = "default", 
                                                 path = "") {
    
    # Primero quedarse con los casos de interes
    corridasLongSub <- filter(corridasLong, 
                              metodo == met,
                              efectoTams %in% efTam,
                              desvios == desv, 
                              alfa %in% nivelSig)
    
    # Puede que alguna prob sea cero, lo cual impide usar el log, sumarle algo
    if (logaritmo) corridasLongSub$prob[corridasLongSub$prob == 0] <- 0.0001
    
    # Data para saber en qué paneles debo graficar la linea de la P(E1)
    modelosEfectos <- data.frame(modelos = rep(1:6, 3),
                                 efecto = rep(c("A", "B", "AB"), each = 6),
                                 dibujarAlfa = c(T, F, F, T, F, F,
                                                 T, T, F, T, T, F,
                                                 T, T, T, F, F, F),
                                 paraLeyenda = "P(E1)obj")
    modelosEfectos <- mutate(modelosEfectos,
                             modelos = factor(modelos, levels = 1:6, labels = paste("Modelo", 1:6)),
                             efecto = factor(efecto, levels = c("A", "B", "AB"), labels = paste("Efecto", c("A", "B", "AB"))))
    
    # Graficar    
    p <- ggplot(corridasLongSub, aes(x = replicas, y = prob, group = paste0(distrErrores, propOutliers), linetype = distrErrores, color = factor(propOutliers)))
    p <- p + geom_line(lwd = 1)
    p <- p + theme(panel.spacing = unit(.05, "lines"),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.25), 
                   strip.background = element_rect(color = "black", size = 0.25),
                   legend.key.width = unit(2, "cm"))
    p <- p + xlab("Número de réplicas") + ylab("Probabilidad Estimada")
    p <- p + scale_color_discrete("Prop. de Outliers")
    
    if (logaritmo) {
        alfaNum <- as.numeric(nivelSig)
        p <- p + facet_grid(efecto ~ modelos, scales = "free_y") # facets con eje vertical libre
        p <- p + coord_trans(y = "log")
        p <- p + scale_y_continuous(breaks = c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1))
        p <- p + geom_hline(data = modelosEfectos[modelosEfectos$dibujarAlfa, ], aes(yintercept = dibujarAlfa * alfaNum, linetype = paraLeyenda), lwd = 0.5)
        p <- p + scale_linetype_manual("Distr. de\nlos errores",
                                       values = c("normal" = 1, "lognormal" = 4, "P(E1)obj" = 2),
                                       labels = c("normal" = "Normal", "lognormal" = "Lognormal", "P(E1)obj" = "P(E1)\nobjetivo"))
        p <- p + guides(linetype = guide_legend(override.aes = list(size = 1))) # definir grosor de linea en leyenda independiete del grafico
    } else {
        p <- p + scale_linetype_manual("Distr. de\nlos errores", 
                                       values = c("normal" = 1, "lognormal" = 4),
                                       labels = c("normal" = "Normal", "lognormal" = "Lognormal"))
        p <- p + facet_grid(efecto ~ modelos)
        p <- p + scale_y_continuous(breaks = seq(0.25, 1, 0.25))
    }
    
    # Agregar titulo
    plotTag <- paste0(met, "_DistErrYPropOut_", "efTam", efTam, "_desv", desv, "_alfa", nivelSig)
    if (logaritmo) plotTag <- paste0(plotTag, "_log")
    if (!is.null(titulo)) {
        if (titulo == "default") titulo <- plotTag
        p <- p + ggtitle(titulo)
    }
    
    if (guardarGraf) {
        if (fileName == "default") fileName <- paste0(plotTag, ".png")
        ggsave(fileName, plot = p, device = "png", path = path, width = 15, height = 7.5, units = "in", dpi = 200)
    }
    invisible(p)
}

graficarProbs_PorMet_PropOut_Var <- function(corridasLong, 
                                             met, 
                                             efTam, 
                                             distr, 
                                             logaritmo = FALSE, 
                                             nivelSig = "0.05", 
                                             titulo = NULL, # o "default"
                                             guardarGraf = F, 
                                             fileName = "default", 
                                             path = "") {
    
    # Primero quedarse con los casos de interes
    corridasLongSub <- filter(corridasLong, 
                              metodo == met,
                              efectoTams %in% efTam,
                              distrErrores == distr,
                              alfa %in% nivelSig)
    
    # Puede que alguna prob sea cero, lo cual impide usar el log, sumarle algo
    if (logaritmo) corridasLongSub$prob[corridasLongSub$prob == 0] <- 0.0008
    
    # Data para saber en qué paneles debo graficar la linea de la P(E1)
    modelosEfectos <- data.frame(modelos = rep(1:6, 3),
                                 efecto = rep(c("A", "B", "AB"), each = 6),
                                 dibujarAlfa = c(T, F, F, T, F, F,
                                                 T, T, F, T, T, F,
                                                 T, T, T, F, F, F),
                                 paraLeyenda = "P(E1)obj")
    modelosEfectos <- mutate(modelosEfectos,
                             modelos = factor(modelos, levels = 1:6, labels = paste("Modelo", 1:6)),
                             efecto = factor(efecto, levels = c("A", "B", "AB"), labels = paste("Efecto", c("A", "B", "AB"))))
    
    # Graficar    
    p <- ggplot(corridasLongSub, aes(x = replicas, y = prob, group = paste0(propOutliers, desvios), linetype = factor(propOutliers), color = desvios))
    p <- p + geom_line(lwd = 1)
    p <- p + theme(panel.spacing = unit(.05, "lines"),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.25), 
                   strip.background = element_rect(color = "black", size = 0.25),
                   legend.key.width = unit(2, "cm"))
    p <- p + xlab("Número de réplicas") + ylab("Probabilidad Estimada")
    p <- p + scale_color_discrete("Variancias",  labels = c("1" = "Homog", "2" = "Heterog\n2 distintas", "6" = "Heterog\n6 distintas"))
    if (logaritmo) {
        alfaNum <- as.numeric(nivelSig)
        p <- p + facet_grid(efecto ~ modelos, scales = "free_y") # facets con eje vertical libre
        p <- p + coord_trans(y = "log")
        p <- p + scale_y_continuous(breaks = c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1))
        p <- p + geom_hline(data = modelosEfectos[modelosEfectos$dibujarAlfa, ], aes(yintercept = dibujarAlfa * alfaNum, linetype = paraLeyenda), lwd = 0.5)
        p <- p + scale_linetype_manual("Prop. de\noutliers",
                                       values = c("0" = 1, "0.1" = 3, "0.3" = 4, "P(E1)obj" = 2),
                                       labels = c("0" = "0", "0.1" = "0.1", "0.3" = "0.3", "P(E1)obj" = "P(E1)\nobjetivo"))
        p <- p + guides(linetype = guide_legend(override.aes = list(size = 1))) # definir grosor de linea en leyenda independiete del grafico
        
    } else {
        p <- p + scale_linetype_manual("Prop. de\noutliers",
                                       values = c("0" = 1, "0.1" = 3, "0.3" = 4),
                                       labels = c("0" = "0", "0.1" = "0.1", "0.3" = "0.3"))
        p <- p + facet_grid(efecto ~ modelos)
        p <- p + scale_y_continuous(breaks = seq(0.25, 1, 0.25))
    }
    
    # Agregar titulo
    plotTag <- paste0(met, "_PropOutYDesv_", "efTam", efTam, "_", distr, "_alfa", nivelSig)
    if (logaritmo) plotTag <- paste0(plotTag, "_log")
    if (!is.null(titulo)) {
        if (titulo == "default") titulo <- plotTag
        p <- p + ggtitle(titulo)
    }
    
    if (guardarGraf) {
        if (fileName == "default") fileName <- paste0(plotTag, ".png")
        ggsave(fileName, plot = p, device = "png", path = path, width = 15, height = 7.5, units = "in", dpi = 200)
    }
    invisible(p)
}

graficarPvalues <- function(n = 10, 
                            efTam = 0.8, 
                            disErr = "normal", 
                            desv = "1", 
                            propOut = 0, 
                            pathArchivos = "",
                            logaritmo = F,
                            titulo = NULL, # o "default"
                            guardarGraf = F,
                            fileName = "default",
                            path = "") {
    
    desv <- as.numeric(desv)
    
    # Crear nombres de archivos para cargar
    buscarEstos <- paste0("n", n, "_ef", efTam, "_", disErr, "_desv", desv, "_out", propOut, "_mod", 1:6, collapse = "|")
    # Listar todos los archivos 
    archivos <- list.files(pathArchivos)
    
    df <- data.frame()
    for (arch in grep(buscarEstos, archivos, value = T)) {
        mod <- str_extract(arch, "(mod[1-6])")
        datos <- read.delim(paste0(pathArchivos, arch))
        datos$modelo <- mod
        df <- rbind(df, datos)
    }
    
    # Pasar a formato largo
    dfLong <- 
        df %>%
        gather(key = pvalueKey, value = pvalue, -modelo) %>%
        extract(pvalueKey, c("efecto", "metodo"), "p(.+)\\.(.+)") %>%
        mutate(modelo = factor(modelo, levels = paste0("mod", 1:6), labels = paste("Modelo", 1:6)),
               efecto = factor(efecto, levels = c("A", "B", "AB"), labels = paste("Efecto", c("A", "B", "AB"))))
    
    if (logaritmo) dfLong$pvalue[dfLong$pvalue == 0] <- 0.008
    
    # Hacer grafico
    p <- ggplot(dfLong, aes(x = metodo, y = pvalue, fill = metodo))
    p <- p + geom_boxplot()
    p <- p + theme(panel.spacing = unit(.05, "lines"),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.25), 
                   strip.background = element_rect(color = "black", size = 0.25))
    p <- p + guides(fill = FALSE)
    p <- p + xlab("Método") + ylab("P-value")
    if (logaritmo) {
        p <- p + facet_grid(efecto ~ modelo, scales = "free_y") 
        p <- p + coord_trans(y = "log")
        p <- p + scale_y_continuous(breaks = c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1))
    } else {
        p <- p + facet_grid(efecto ~ modelo)
    }
    
    plotTag <- paste0("Pvalues_", "n", n, "_efTam", efTam, "_", disErr, "_out", propOut)
    if (logaritmo) plotTag <- paste0(plotTag, "_log")
    if (!is.null(titulo)) {
        if (titulo == "default") titulo <- plotTag
        p <- p + ggtitle(titulo)
    }
    
    if (guardarGraf) {
        if (fileName == "default") fileName <- paste0(plotTag, ".png")
        ggsave(fileName, plot = p, device = "png", path = path, width = 15, height = 7.5, units = "in", dpi = 200)
    }
    invisible(p)
}

graficarProbs_PorMet_DistErr_EfTam <- function(corridasLong, 
                                                 met, 
                                                 propOut, 
                                                 desv, # "1" "2" "6"
                                                 logaritmo = FALSE, 
                                                 nivelSig = "0.05", 
                                                 titulo = NULL, # o "default"
                                                 guardarGraf = F, 
                                                 fileName = "default", 
                                                 path = "") {
    
    # Primero quedarse con los casos de interes
    corridasLongSub <- filter(corridasLong, 
                              metodo == met,
                              propOutliers == propOut,
                              desvios == desv, 
                              alfa %in% nivelSig)
    
    # Puede que alguna prob sea cero, lo cual impide usar el log, sumarle algo
    if (logaritmo) corridasLongSub$prob[corridasLongSub$prob == 0] <- 0.0001
    
    # Data para saber en qué paneles debo graficar la linea de la P(E1)
    modelosEfectos <- data.frame(modelos = rep(1:6, 3),
                                 efecto = rep(c("A", "B", "AB"), each = 6),
                                 dibujarAlfa = c(T, F, F, T, F, F,
                                                 T, T, F, T, T, F,
                                                 T, T, T, F, F, F),
                                 paraLeyenda = "P(E1)obj")
    modelosEfectos <- mutate(modelosEfectos,
                             modelos = factor(modelos, levels = 1:6, labels = paste("Modelo", 1:6)),
                             efecto = factor(efecto, levels = c("A", "B", "AB"), labels = paste("Efecto", c("A", "B", "AB"))))
    
    # Graficar    
    p <- ggplot(corridasLongSub, aes(x = replicas, y = prob, group = paste0(distrErrores, efectoTams), linetype = distrErrores, color = factor(efectoTams)))
    p <- p + geom_line(lwd = 1)
    p <- p + theme(panel.spacing = unit(.05, "lines"),
                   panel.border = element_rect(color = "black", fill = NA, size = 0.25), 
                   strip.background = element_rect(color = "black", size = 0.25),
                   legend.key.width = unit(2, "cm"))
    p <- p + xlab("Número de réplicas") + ylab("Probabilidad Estimada")
    p <- p + scale_color_discrete("Efecto tamaño")
    
    if (logaritmo) {
        alfaNum <- as.numeric(nivelSig)
        p <- p + facet_grid(efecto ~ modelos, scales = "free_y") # facets con eje vertical libre
        p <- p + coord_trans(y = "log")
        p <- p + scale_y_continuous(breaks = c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1))
        p <- p + geom_hline(data = modelosEfectos[modelosEfectos$dibujarAlfa, ], aes(yintercept = dibujarAlfa * alfaNum, linetype = paraLeyenda), lwd = 0.5)
        p <- p + scale_linetype_manual("Distr. de\nlos errores",
                                       values = c("normal" = 1, "lognormal" = 4, "P(E1)obj" = 2),
                                       labels = c("normal" = "Normal", "lognormal" = "Lognormal", "P(E1)obj" = "P(E1)\nobjetivo"))
        p <- p + guides(linetype = guide_legend(override.aes = list(size = 1))) # definir grosor de linea en leyenda independiete del grafico
    } else {
        p <- p + scale_linetype_manual("Distr. de\nlos errores", 
                                       values = c("normal" = 1, "lognormal" = 4),
                                       labels = c("normal" = "Normal", "lognormal" = "Lognormal"))
        p <- p + facet_grid(efecto ~ modelos)
        p <- p + scale_y_continuous(breaks = seq(0.25, 1, 0.25))
    }
    
    # Agregar titulo
    plotTag <- paste0(met, "_DistErrYEfTam_", "propOut", propOut, "_desv", desv, "_alfa", nivelSig)
    if (logaritmo) plotTag <- paste0(plotTag, "_log")
    if (!is.null(titulo)) {
        if (titulo == "default") titulo <- plotTag
        p <- p + ggtitle(titulo)
    }
    
    if (guardarGraf) {
        if (fileName == "default") fileName <- paste0(plotTag, ".png")
        ggsave(fileName, plot = p, device = "png", path = path, width = 15, height = 7.5, units = "in", dpi = 200)
    }
    invisible(p)
}
