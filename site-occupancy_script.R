# Máster David Ferrer 15/04/2020

################################################################################

# Funcionalidades básicas con unmarked
# Este script está sacado de: https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1
# Los datos pueden adquirirse en esa dirección


# Instalamos y cargamos las librerías necesarias
install.packages("unmarked")
library(unmarked)

# Seleccionamos nuestro directorio de trabajo
setwd("/home/javifl/IREC/master_david")

# Cargamos los datos
data <- read.csv("blgr.csv")
# Vemos las primeras filas para explorar su estructura
head(data)

# La estructura de los datos es la siguiente:
# Columna 1: ID asignado a cada uno de los sitios muestreados
# Columnas 2-4: representan las veces que hemos repetido los muestreos (réplicas), en este caso 3. 0= ausencia, 1= presencia
# Columnas 5-9: variables predictoras medidas en cada sitio (covariables)
# Columnas 10-12: número de animales detectados en cada réplica del muestreo. No los vamos a usar en principio

# Vamos a guardar en "y" los datos sobre la detección de la especie en cada uno 
# de las repeticiones de los muestreos realizados y en "n" el número de sitios visitados
y <- data[,2:4]
n <- nrow(data)

# Guardamos las covariables en el objeto "blgr.site"
blgr.site <- data[,5:9]

# Vamos a crear un factor que indique el momento del muestreo (1, 2 o 3). Esto sirve
# únicamente para ordenar de otra manera nuestros datos. Quedará claro más adelante
time <- as.factor(rep(c(1,2,3),n))
blgr.obs <- data.frame(time)

# Por último ponemos en un único objeto el histórico de las detecciones y las covariables
# en un objeto especial que unmarked entiende, lo vamos a llamar "blgr"

blgr <- unmarkedFrameOccu(y = y, siteCovs = blgr.site,obsCovs=blgr.obs)

# Miramos qué pinta tiene y le pedimos un resumen:
head(blgr)
summary(blgr)

# Ya tenemos los datos listos

################################################################################

# Ahora tenemos que correr los modelos. El pseudocódigo (sintaxis) del modelo es:
# model<-occu(~detection_formula ~occupancy_formula, dataframe)
# Donde:
# model es el nombre del objeto donde guardaremos el modelo
# occu es la función que vamos a usar para correr el modelo
# ~ detection_formula es la fórmula que describe la probabilidad de detección
# ~ occupancy_formula es la fórmula que descrbe la probabilidad de ocupación (relación entre la presencia/ausencia de la especie y las covariables)
# dataframe es el objeto creado con la función unmarkedFrameOccu()

# Si queremos ajustar un modelo en el que ni hay efecto de la probabilidad de detección (probabilidad detección constante)
# ni de las covariables (modelado del intercepto únicamente), es decir, un modelo nulo, podríamos ejecutar:

fm1<-occu(~1 ~1,blgr)
fm1

# Aquí tenemos los Estimates. Ojo, porque unmarked nos los da transformados (logit en este caso)
# Para obtenerlos sin transformar podremos usar la función backTransform() indicando el estimate
# que queremos "destransformar"

# Para la probabilidad de detección:
backTransform(fm1,'det')

# Para la probabilidad de ocupación 
backTransform(fm1,'state')

# Finalmente, podemos calcular la probabilidad de presencia de la especie en cada
# uno de nuestros sitios y sus intervalos de confianza
s <- nrow(y)
re <- ranef(fm1)
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
print("95 % EB interval on number sites occupied")
rbind(s_occup = c(Estimate = sum(EBUP), colSums(CI)) )
print("95 % EB interval on proportion of sites occupied")
rbind(PAO = c(Estimate = sum(EBUP), colSums(CI))/s )

################################################################################

# Ahora crearemos un set de modelos más complejos, en los que incluiremos la probabilidad
# de detección y las covariables al modelo de ocupación:

# Probabilidad de deteccion dependiente de la réplica, sin efecto de las covariables 
fm2<-occu(~time ~1, blgr)
# Probabilidad de deteccion constante, pero la ocupación viene predicha por la covariable bqi
fm3<-occu(~1 ~bqi, blgr)
# Probabilidad de deteccion dependiente de la réplica, ocupación predicha por la covariable bqi
fm4<-occu(~time ~bqi, blgr)
# Probabilidad de deteccion constante, ocupación predicha por el logaritmo del tamaño log(field size)
fm5<-occu(~1 ~log(field.size), blgr)
# Probabilidad de deteccion constante, ocupación predicha por el logaritmo del tamaño log(field size) y por bqi
fm6<-occu(~1 ~log(field.size) + bqi, blgr)

# Ahora podemos ponerlos todos juntos y hacer una selección de modelos en base al AIC
fmlist<-fitList(primero=fm1, segundo=fm2, tercero=fm3, cuarto=fm4, quinto=fm5, sexto=fm6)
modSel(fmlist)

# Esta selección nos dice que el modelo con mayor soporte es el nulo, el primero que hicimos
# A partir de aquí podríamos hacer un multi-model average, para quedarnos con una "media" de todos,
# quedarnos con alguno en concreto, etc.
# Para continuar: https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1

################################################################################
# 16/04/2020

# Mapear la probabilidad de ocupación (predicciones a raster con unmarked)
# Este script está sacado de https://cran.r-project.org/web/packages/unmarked/vignettes/spp-dist.pdf

library(unmarked)

# Cargamos los datos y preparamos un dataset que pueda leer unmarked
data(crossbill)
umf <- unmarkedFrameOccu(
  y=as.matrix(crossbill[,c("det991", "det992", "det993")]),
  siteCovs=crossbill[,c("ele", "forest")],
  obsCovs=list(date=crossbill[,c("date991", "date992", "date993")]))

# Estandarizamos las covariables
sc <- scale(siteCovs(umf))
siteCovs(umf) <- sc

# Vemos la pinta de la base de datos:
# Se ha muestreado la presencia/ausencia de piquituertos en 267 puntos tres veces a lo largo del año 1999
# Se han tomado dos covariables predictoras para cada uno de los sitios (elevación y cobertura forestal)
# Columna 1-3: Presencia/Ausencia de piquiruertos en cada uno de los puntos de muestreo en las tres réplicas 
# Columnas 4 y 5: Covariables predictoras
# Columnas 6-8: Fecha de cada una de las réplicas (en días julianos)

head(umf)

# Ajustamos un modelo de forma similar a lo que ya hemos hecho antes:
# Creemos que la probabilidad de detección depende del día del año que hemos hecho la réplica (~ date)
# Creemos que la probabilidad de ocupación depende de la elevación y su término cuadrático y de la cobertura de bosque
# Corremos el modelo

fm.occu <- occu(~date ~ele + I(ele^2) + forest, umf)
fm.occu

# Aquí podemos ver bien las dos partes del modelo, la parte de ocupación y la de detección,
# con los estimates y p-values. Una vez ajustado el modelo, dado que la ocupación depende de 
# la elevación y la cobertura foresta, si tenemos los valores de esas covariables para cada una
# de las celdas de un mapa, podemos "resolver la ecuación" para todas las celdas, es decir, realizar
# una predicción de la probabilidad de ocupación en cada una de las celdas. Vamos a preparar el mapa (de Suiza en este caso):

data(Switzerland)
print(levelplot(elevation ~ x + y, Switzerland, aspect="iso",
                xlab="Easting (m)", ylab="Northing (m)",
                col.regions=terrain.colors(100)))

library(raster)
elevation <- rasterFromXYZ(Switzerland[,c("x","y","elevation")],
                           crs="+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")
forest <- rasterFromXYZ(Switzerland[,c("x","y","forest")],
                        crs="+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")

# Recuerda que hemos utilizado las covariables estandarizadas, ahora tenemos que hacer lo mismo 
# con los mapas para hacer las predicciones
attr(sc, "scaled:center")
attr(sc, "scaled:scale")
ele.s <- (elevation-1189)/640
forest.s <- (forest-34.7)/27.7
ef <- stack(ele.s, forest.s)
names(ef) <- c("ele", "forest")
plot(ef, col=terrain.colors(100))

# Ahora podemos hacer las predicciones (tarda un rato porque hay bastantes celdas en el mapa)
E.psi <- predict(fm.occu, type="state", newdata=ef)
plot(E.psi, axes=FALSE, col=terrain.colors(100))







