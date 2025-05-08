## Paquetes de apoyo necesarios ##
# install.packages("metafor") # Si no se tuviera el paquete instalado
library(metafor) # Se usa para calcular la varianza de la estimación combinada

## -power.ma-       Función para el cálculo de la potencia para un contraste unilateral derecho ##

# Argumentos #
# -ES-          Columna de un data frame con los tamaños del efecto de cada estudio (Para metafor)
# -var-         Columna de un data frame con las varianzas de cada estudio (Para metafor)
# -data-    Data frame con los datos (TEs y varianzas) (Para metafor)
# -method-      Método de estimación de tau^2 (Para metafor). Los códigos más comunes son (ejecute ?rma.uni para todos los códigos):
  # -"FE"-        (Efecto fijo. No se estima tau^2)
  # -"REML"-      (Máxima Verosimilitud Restringida). Por defecto
# -theta-   TE paramétrico alternativo. Por defecto se calcula la potencia para los puntos de corte clásicos propuestos por Cohen (0.2, 0.5, 0.8) y siempre se ofrecen estos valores a mayores.
# -alpha-   Nivel de significación. Por defecto = 0.05

# Programación de la función #
power.ma <- function(ES, var, data, method="REML", theta=NULL, alpha=0.05) {
  
  # Vector de tamaños del efecto alternativos ordenados
  theta <- sort(c(theta, 0.2, 0.5, 0.8))
  theta <- theta[!duplicated(theta)] # Se eliminan los duplicados si los hubiere
  
  # Cálculo del error típico
  ma <- rma.uni(yi=es, vi=var, data=data, method=method, weighted=T)
  SE <- ma$se # Valor del error típico
  
  # Tamaño del efecto mínimo para alcanzar significación estadística
  ES.alpha <- qnorm(1 - alpha) * SE
  
  # Cálculo final de la potencia
  power <- 1 - pnorm((ES.alpha - theta)/SE)
  
  # Retorno de la función
  print(data.frame(ES = theta,
                   Power = power))
} 


## -range.CI-       Función para el cálculo de la amplitud del IC de la estimación combinada del TE ##

# Argumentos #
# -ES-          Columna de un data frame con los tamaños del efecto de cada estudio (Para metafor)
# -var-         Columna de un data frame con las varianzas de cada estudio (Para metafor)
# -data-    Data frame con los datos (TEs y varianzas) (Para metafor)
# -method-      Método de estimación de tau^2 (Para metafor). Los códigos más comunes son (ejecute ?rma.uni para todos los códigos):
  # -"FE"-        (Efecto fijo. No se estima tau^2)
  # -"REML"-      (Máxima Verosimilitud Restringida). Por defecto
# -alpha-   Nivel de significación. Por defecto = 0.05

range.CI <- function(ES, var, data, method="REML", alpha=0.05) {
  
  # Cálculo del error típico
  ma <- rma.uni(yi=ES, vi=var, data=data, method=method, weighted=T)
  SE <- ma$se # Valor del error típico
  
  # Cálculo de la amplitud
  range <- 2 * SE * abs(qnorm(1 - alpha/2))
  
  # Retorno de la función
  print(range)
}


## -studies.power-  Función para el cálculo prospectivo del número de estudios/sujetos necesarios para alcanzar una determinada potencia ##

# Argumentos #
# -ES-          Columna de un data frame con los tamaños del efecto de cada estudio (Para metafor)
# -var-         Columna de un data frame con las varianzas de cada estudio (Para metafor)
# -data-    Data frame con los datos (TEs y varianzas) (Para metafor)
# -method-      Método de estimación de tau^2 (Para metafor). Los códigos más comunes son (ejecute ?rma.uni para todos los códigos):
  # -"FE"-        (Efecto fijo. No se estima tau^2)
  # -"REML"-      (Máxima Verosimilitud Restringida). Por defecto
# -theta-   TE paramétrico alternativo. Por defecto se calcula la potencia para los puntos de corte clásicos propuestos por Cohen (0.2, 0.5, 0.8) y siempre se ofrecen estos valores a mayores.
# -alpha-   Nivel de significación. Por defecto = 0.05
# -power-   Valor de potencia que se desea alcanzar.
# -plot.limit-  Total de estudios que se desea estimar para el gráfico. Por defecto = 50

studies.power <- function(ES, var, data, method="REML", theta=NULL, alpha=0.05, power, plot.limit=50) {
  
  # Vector de tamaños del efecto alternativos ordenados
  theta <- sort(c(theta, 0.2, 0.5, 0.8))
  theta <- theta[!duplicated(theta)] # Se eliminan los duplicados si los hubiere
  
  # Cálculo de los datos del estudio original
  ma <- rma.uni(yi=es, vi=var, data=data, method=method, weighted=T)
  SE_k <- ma$se # Valor del error típico original
  
  W_k <- 1/SE_k^2 # Peso total de los estudios originales
  
  w.mean_k <- W_k/length(data$es) # Peso medio de los estudios originales
  
  # Cálculo del error típico de los k+j estudios
  SE_kj <- theta/(qnorm(1 - alpha) - qnorm(1 - power))
  
  # Cálculo del peso de los k+j estudios
  W_kj <- 1/SE_kj^2
  
  # Cálculo del peso de los j estudios
  W_j <- W_kj - W_k
  
  # Cálculo de la varianza de los j estudios
  var_j <- 1/W_j
  
  # Cálculo del tamaño muestral estimado (sólo si method = EF)
  if (method == "FE") {
    Subjets <- 1:length(theta)
    for (i in 1:length(Subjets)) {
      Subjets[i] <- ceiling((4/var_j[i]) + ((2 * theta[i]^2)/(4 * var_j[i])))
    }
  }
  
  # Cálculo de los estudios estimados
  Studies <- ceiling(W_j/w.mean_k)
  
  # Cálculo de la potence estudio a estudio (para el gráfico)
  # Data frame para contener los datos
  data.plot <- data.frame(km = 0:(plot.limit-1))
  
  for (i in 1:length(theta)) {
    # Vector con las potencias para un TE determinado
    c <- 1 - pnorm(((qnorm(1 - alpha) * SE_k) - theta[i])/SE_k)
    
    # Mecanismo de control del bucle
    meter <- 0
    
    # Bucle
    repeat {
      # Actualización del contador
      meter <- meter + 1
      
      # Re-cálculo del error típico con los k+j estudios
      Wkj <- (1/SE_k^2) + (meter * w.mean_k)
      SEkj <- sqrt(1/Wkj)
      
      # Re-cálculo de la potencia
      powerkj <- 1 - pnorm(((qnorm(1 - alpha) * SEkj) - theta[i])/SEkj)
      
      # Asignación de la potencia al vector
      c[meter] <- powerkj
      
      # Cierre del bucle
      if (meter == plot.limit)
        break
    }
    
    # Unión del nuevo vector al data frame
    data.plot <- cbind(data.plot, c)
  }
  
  # Renombrar columnas del data frame
  cnames <- c("km", paste("ES=", theta, sep=""))
  colnames(data.plot) <- cnames
  
  # Creación del gráfico (Primera línea)
  plot(x=data.plot[,1], y=data.plot[,2], 
       type="l", lwd = 2,
       xlab = "Número de estudios", ylab = "Potencia",
       ylim=c(min(data.plot[,2]),max(data.plot[,ncol(data.plot)])), col=2)
  # Siguiente líneas
  for (i in 3:ncol(data.plot)) {
    lines(x=data.plot[,1], y=data.plot[,i], col=i, lwd=2, type="l")
  }
  # Líneas verticales para el número de estudios que alcanzan la potencia introducida
  for (i in 1:length(cnames)) {
    abline(v=Studies[i], col=(i+1))
  }
  
  # Leyenda
  legend("bottomright",
         inset=0.03,
         legend=cnames[-1],
         col=c(2:ncol(data.plot)),
         lwd=2)
  
  # Retorno de la función (para method = FE)
  if (method == "FE"){
    results <- list(Subjets = data.frame(ES = theta,
                                         Subjets = Subjets),
                    Studies = data.frame(ES = theta,
                                         Studies = Studies), 
                    Data.plot = data.plot)
  } else {
    results <- list(Studies = data.frame(ES = theta,
                                         Studies = Studies), 
                    Data.plot = data.plot)
  }
  
  # Se coercionan los valores negativos y NAs a 0
  for (i in 1:length(names(results))){
    results[[i]][results[[i]] < 0] <- 0
    results[[i]][is.na(results[[i]])] <- 0
  }
  
  print(results)
}


## -studies.range- Función para el cálculo prospectivo del número de estudios/sujetos necesarios para alcanzar una determinada amplitud del IC ##

# Argumentos #
# -ES-          Columna de un data frame con los tamaños del efecto de cada estudio (Para metafor)
# -var-         Columna de un data frame con las varianzas de cada estudio (Para metafor)
# -data-    Data frame con los datos (TEs y varianzas) (Para metafor)
# -method-      Método de estimación de tau^2 (Para metafor). Los códigos más comunes son (ejecute ?rma.uni para todos los códigos):
  # -"FE"-        (Efecto fijo. No se estima tau^2)
  # -"REML"-      (Máxima Verosimilitud Restringida). Por defecto
# -alpha-   Nivel de significación. Por defecto = 0.05
# -range-   Valor de amplitud del IC que se desea alcanzar.
# -theta-   TE paramétrico alternativo (Solo necesario si method="FE" y metric="SMD". Por defecto se calcula la potencia para los puntos de corte clásicos propuestos por Cohen (0.2, 0.5, 0.8) y siempre se ofrecen estos valores a mayores.
# -plot.limit-  Total de estudios que se desea estimar para el gráfico. Por defecto = 50

studies.range <- function(ES, var, data, method="REML", alpha=0.05, range, theta=NULL, plot.limit=50) {
  
  # Vector de tamaños del efecto alternativos ordenados (Solo se usa si method="FE" y metric="SMD")
  es <- sort(c(theta, 0.2, 0.5, 0.8))
  es <- es[!duplicated(es)] # Se eliminan los duplicados si los hubiere
  
  # Cálculo de los datos del estudio original
  ma <- rma.uni(yi=ES, vi=var, data=data, method=method, weighted=T)
  SE_k <- ma$se # Valor del error típico original
  
  W_k <- 1/SE_k^2 # Peso total de los estudios originales
  
  w.mean_k <- W_k/length(data$ES) # Peso medio de los estudios originales
  
  # Cálculo del error típico de los k+j estudios
  SE_kj <- range/(2 * abs(qnorm(1 - alpha/2)))
  
  # Cálculo del peso de los k+j estudios
  W_kj <- 1/SE_kj^2
  
  # Cálculo del peso de los j estudios
  W_j <- W_kj - W_k
  
  # Cálculo de la varianza de los j estudios
  var_j <- 1/W_j
  
  # Cálculo del tamaño muestral estimado (solo para method="FE")
  if (method == "FE") {
    Subjets <- 1:length(es)
    for (i in 1:length(Subjets)) {
      Subjets[i] <- ceiling((4/var_j) + ((2 * es[i]^2)/(4 * var_j)))
    }
  }
  
  # Cálculo de los estudios estimados
  Studies <- ceiling(W_j/w.mean_k)
  
  # Cálculo de la amplitud estudio a estudio (para el gráfico)
  # Data frame que recoge la información necesaria para el gráfico (Nº de estudios añadidos y la amplitud correspondiente)
  data.plot <- data.frame(km = 0,
                          range = range <- 2 * SE_k * abs(qnorm(1 - alpha/2))
  )
  
  # Mecanismo de control del bucle
  meter <- 0
  
  # Bucle
  repeat {
    # Actualización del contador
    meter <- meter + 1
    
    # Re-cálculo del error típico con los k+j estudios
    Wkj <- (1/SE_k^2) + (meter * w.mean_k)
    SEkj <- sqrt(1/Wkj)
    
    # Re-cálculo de la amplitud
    rangekj <- 2 * SEkj * abs(qnorm(1 - alpha/2))
    
    # Asignación de los nuevos datos al data frame
    data.plot <- rbind(data.plot, c(meter, rangekj))
    
    # Cierre del bucle
    if (meter == (plot.limit - 1))
      break
  }
  
  # Gráfico
  plot(data.plot, type="l", xlab = "Número de estudios", ylab = "Amplitud", lwd = 2)
  abline(v = Studies)
  
  
  # Retorno de la función (para method = FE)
  if (method == "FE"){
    if (metric == "SMD"){
      results <- list(Subjets = data.frame(ES = es,
                                           Subjets = Subjets),
                      Studies = Studies, 
                      Data.plot = data.plot)
    } else if (metric == "PC"){
      results <- list(Subjets = Subjets,
                      Studies = Studies, 
                      Data.plot = data.plot)
    }
    
  } else {
    results <- list(Studies = Studies, 
                    Data.plot = data.plot)
  }
  
  # Se coercionan los valores negativos y NAs a 0
  for (i in 1:length(names(results))){
    results[[i]][results[[i]] < 0] <- 0
    results[[i]][is.na(results[[i]])] <- 0
  }
  
  print(results)
}
