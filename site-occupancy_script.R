# Máster David Ferrer 

# Algunos links útiles sobre site-occupancy models:


#https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1

#https://cran.r-project.org/web/packages/unmarked/vignettes/spp-dist.pdf

#https://a4688ce6-a-62cb3a1a-s-sites.googlegroups.com/site/asrworkshop/home/schedule/r-occupancy-1/unmarked.pdf?attachauth=ANoY7coxpbxl7F0uTBsRdUnVHzhJnhV9Nu7xyt6XeWvI9375NJ3dDVGEspvUqMp4arNi1pFJ_X8VggpqrsLertaWJ5xvZbIVnpxa1zG2CmihMDFAgg0VvCzTMLbFl7NCPnpLIU0Puqf3xWCnYmAbdvcFiyzCRs14JwlaRCd4f3apbefxBBXKsPEIJpGkRf8cMW5kiZ9_HO5D0ijPbt5J3iNnZkRv8rpam_20-u_Ij8m5uhoc6l3X9xOj2p1VrXP3C_ywFJHWwKX0&attredirects=0

#https://cornelllabofornithology.github.io/ebird-best-practices/occupancy.html

#https://cals.arizona.edu/classes/wfsc578/PointCountsAbundance_ELOW.R

################################################################################
# 15/04/2020
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
# que queremos "detransformar"

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


################################################################################
#23/04/2020

# Preparación de los datos con R

# Para añadir las covariables ambientales a nuestro dataset a partir de cargografía digital
# podemos usar un GIS o lo podemos hacer diréctamente con R.

# Instalamos los paquetes si fuese necesario
#install.packages("rgdal")
#install.packages("raster")

library(rgdal)
library(raster)

# Vamos a cargar la capa shape con las variables ambientales
covar_shape <- readOGR("/home/javifl/IREC/master_david/variables/variables_rx.shp")

# Vamos a echar un vistazo a la tabla de atributos de la capa:
head(covar_shape@data)

# Vamos a ver cuantas cuadrículas tiene esta capa (es decir, el número de filas de la tabla de atributos)
nrow(covar_shape@data)
# Esta capa tiene 29532 cuadrículas

# Vemos que tiene 14 columnas: SP_ID, ID, X, Y, HA, FINCA_N...
# Ahora queremos los valores de esas columnas para cada uno de nuestras cámaras

# Cargamos las localizaciones de las cámaras. Para agilizar el proceso hemos creado
# un archivo .csv con tres columnas: el ID, la coordenada X y la coordenada Y
# OJO! Las coordenadas de las cámaras están en ETRS89/UTM-30 mientras que la capa
# de covariables está en ETRS89/UTM-29. Se ha cambiado en este nuevo fichero.

xy_2015 <- read.csv("/home/javifl/IREC/master_david/xy_2015.csv")

# Echamos un vistazo a esta capa de puntos:
head(xy_2015)

# Vamos a comprobar cuantas cámaras tenemos (es decir, número de filas de "xy_2015")
nrow(xy_2015)
# 38 cámaras


# Podemos graficar ambos objetos, (la capa de covariables ve muy negra por la densidad de cuadrículas)
plot(covar_shape)
points(xy_2015[,2:3], col = "red")


# Ahora podemos extraer los valores de la capa "covar_shape" para cada uno de los puntos de la capa "xy_2015"
# Esto se suele hacer en QGIS con la herramienta "Point Sampling Tools" (el pincho moruno)
# Para ello utilizamos la función "extract" del paquete "raster".
# OJO! La primera columna de la capa "xy_2015" es el ID y no se debe utilizar. Por eso usamos "xy_2015[,2:3]"
covar_2015 <- extract(covar_shape, xy_2015[,2:3])

# Vamos a ver el número de filas (debería coincidir con el número de cámaras)
nrow(covar_2015)
# 38 filas

# Echamos un vistazo al nuevo objeto:
head(covar_2015)

# Finalmente, vemos que las dos primeras columnas nos marcan los "ID" tanto de la capa de cámaras ("point.ID")
# como de la capa de covariables ("poly.ID"). A continuación tenemos todos los campos (columnas) de la capa
# de covariables para cada una de las localizaciones de las cámaras. Podemos exportar como .csv este objeto
# para poder abrirlo con excel:
write.csv(covar_2015, "miruta/covar_2015.csv")


################################################################################
# 08/05/2020
# Calibración del N-mixture model Royle (2004)
# setwd("/home/javifl/IREC/master_david/pcount")

data <- read.csv("datos_2015limpio.csv")

# Seleccionamos una especie
dataciervo <- data[data$sp == "ciervo",]

# Nos quedamos solo con 4 ocasiones, las que mas datos tiene. Habrá que probar 
# con todas las ocasiones y ver qué pasa con los NA
dataciervo <- dataciervo[,c(1:2,9:12, 21:28)]

# De momento vamos a omitir los NA
dataciervo <- na.omit(dataciervo)

# Guardamos los conteos de cienca y el número de sitios
y <- dataciervo[,3:6]
n <- nrow(dataciervo)

# Guardamos las covariables. Ojo, vamos a estandarizarlas con la función scale
d2015.site <- data.frame(scale(dataciervo[,7:14]))

# Guardamos las ocasiones en time
time <- as.factor(rep(c(1:4),n))
d2015.obs <- data.frame(time)

# Lo juntamos todo 
d2015c <- unmarkedFramePCount(y = y, siteCovs = d2015.site, obsCovs = d2015.obs)

# Calibramos los modelos
fm1 <- pcount(~1 ~dvera, d2015c, K = 150)

fm2 <- pcount(~time ~dvera, d2015c, K = 150)

fm3 <- pcount(~1 ~dwat, d2015c, K = 150)

fm4 <- pcount(~time ~dwat, d2015c, K = 150)

fm5 <- pcount(~1 ~v1+v2+v3+v4+v5+v6, d2015c, K = 150)

fm6 <- pcount(~time ~v1+v2+v3+v4+v5+v6, d2015c, K = 150)

fm7 <- pcount(~1 ~dvera + dwat, d2015c, K = 150)

fm8 <- pcount(~time ~dvera + dwat, d2015c, K = 150)

fm9 <- pcount(~1 ~v1+v2+v3+v4+v5+v6+dvera+dwat, d2015c, K = 150)

fm10 <- pcount(~time ~v1+v2+v3+v4+v5+v6+dvera+dwat, d2015c, K = 150)


################################################################################
# 19/05/2020
# Proyección de las predicciones de N-mixture model Royle (2004)
# setwd("/home/javifl/IREC/master_david/pcount")
# library(unmarked)
# library(raster)

# Vamos a correr un modelo para ciervo usando la función pcount de igual forma 
# que la vez anterior

# Preparamos los datos
data <- read.csv("datos_2015limpio.csv")
dataciervo <- data[data$sp == "ciervo",]
dataciervo <- dataciervo[,c(1:2,9:12, 21:28)]
dataciervo <- na.omit(dataciervo)
y <- dataciervo[,3:6]
n <- nrow(dataciervo)
d2015.site <- data.frame(scale(dataciervo[,7:14]))
time <- as.factor(rep(c(1:4),n))
d2015.obs <- data.frame(time)
d2015c <- unmarkedFramePCount(y = y, siteCovs = d2015.site, obsCovs = d2015.obs)

# Corremos un modelo nulo. Si no ponemos nada, pcount() asume que la abundancia sigue una Poisson
null <- pcount(~1 ~1, data = d2015c, K = 150) 
null

# Convertimos los estimates de link a la escla original. Ojo, esta función solo sirve
# cuando los procesos de detección / abundancia (state) no tienen covariables
backTransform(null, type="state")	# Abundancia
backTransform(null, type="det")		# Probabilidad de detección

# Una forma de hacer lo mismo sería:
exp(coef(null, type="state"))		# Abundancia
plogis(coef(null, type="det"))		# Probabilidad de detección

# Vamos a ajustar un modelo más complejo:
Time.VeraWat <- pcount(~time ~dvera+dwat, d2015c, K = 150)
Time.VeraWatLC <- pcount(~time ~dvera+dwat+v1+v2+v3+v4+v5+v6, d2015c, K = 150)

# Lo inspeccionamos
summary(Time.VeraWatLC)

# Aquí podemos ver las dos partes de las que se compone el modelo. Por un lado la abundancia,
# que viene determinada por las covariables dvera y dwat, y por otro lado la probabilidad de
# detección, que viene determinada por las ocasiones en las que hemos muestreado. Podemos
# hacernos algunas preguntas:

# Según este modelo
# ¿Qué variable es la que más influye en la abundancia de ciervos? ¿Por qué?

# Examinamos los estimates (backtransformados)
coef(Time.VeraWat, type="state")      # link scale (= log)
exp(coef(Time.VeraWatLC, type="state")) # backtransformado (= elevado sobre la base e)
coef(Time.VeraWat, type="det")         # link scale (= logit)
plogis(coef(Time.VeraWat, type="det")) # backtransformado (= inversa de logit)

# Vamos a relacionar la probabilidad de detección con respecto al time
barplot(plogis(coef(Time.VeraWat, type="det")), xlab="time", ylab= "probabilidad de detección", ylim =c(0,1))

# ¿En qué ocasión es más probable que detectemos un ciervo? ¿Qué podría significar?

# Ahora vamos a cargar las variables en raster con el pequete raster
library(raster)
variables <- stack(list.files(path="/home/javifl/IREC/master_david/variables_raster",pattern='*.tif', full.names=TRUE))
plot(variables)

# Como hemos corrido el modelo con las variables transformadas, para realizar predicciones
# sobre estas variables también debemos estandarizarlas de la misma manera
# Vemos como están centradas y aplicamos la fórmula
scale(dataciervo[,7:14])

variables$dvera <- (variables$dvera - 2053) / 2086
variables$dwat <- (variables$dwat - 647.01) / 420.4
variables$v1 <- (variables$v1 - 37.9) / 53.1
variables$v2 <- (variables$v2 - 58) / 62.6
variables$v3 <- (variables$v3 - 31.8) / 54.4
variables$v4 <- (variables$v4 - 3.1) / 15.9
variables$v5 <- (variables$v5 - 7.4) / 14.3
variables$v6 <- (variables$v6 - 16.5) / 38.7

# Una vez transformadas, podemos usar estas variables para proyectar las predicciones de abundancia sobre
# el mapa de nuestras covariables. Esto sería "resolver la ecuación" para cada una de las celdillas
# de nuestro mapa, dado que conocemos el valor de todas las covariables ambientales para todo el
# área de estudio.
# Esto puede tardar un poco...
Time.VeraWat_pred <- predict(Time.VeraWat, newdata=variables, type = "state")
Time.VeraWatLC_pred <- predict(Time.VeraWatLC, newdata=variables, type = "state")

# Como resultado tenemos otro stack de rasters en el que se nos muestran las predicciones de abundancias, el
# error estandar y los máximos y mínimos
plot(Time.VeraWat_pred$Predicted, axes=FALSE, col=terrain.colors(100))
plot(Time.VeraWatLC_pred$Predicted, axes=FALSE, col=terrain.colors(100))

# Al hacerlo con las variabels de usos de suelo los resultados salen bastante mal. 
# Mi impresión es que tiene que ver con el tipo de varable, que asumimos como continuas cuando creo que
# debería ser más bien un factor...


################################################################################
# 15/06/2020
# Análisis de componentes principales para las variables de land-cover

library(raster)
library(ade4)
library(factoextra)
variables <- stack(list.files(path="/home/javifl/IREC/master_david/variables_raster",pattern='*.tif', full.names=TRUE))
cover <- variables[[3:8]]
help(dudi.pca)

cover_pca <- dudi.pca(data.frame(na.omit(cover[[1]][]),na.omit(cover[[2]][]),na.omit(cover[[3]][]),
                                 na.omit(cover[[4]][]),na.omit(cover[[5]][]),na.omit(cover[[6]][])), nf =6, scannf = FALSE)


fviz_pca_var(cover_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


cover_pca$c1

pc1 <- cover[[1]]
pc1[!is.na(pc1[])] <- cover_pca$li$Axis1

pc2 <- cover[[1]]
pc2[!is.na(pc2[])] <- cover_pca$li$Axis2

pc3 <- cover[[1]]
pc3[!is.na(pc3[])] <- cover_pca$li$Axis3

pc4 <- cover[[1]]
pc4[!is.na(pc4[])] <- cover_pca$li$Axis4

pc5 <- cover[[1]]
pc5[!is.na(pc5[])] <- cover_pca$li$Axis5

pc6 <- cover[[1]]
pc6[!is.na(pc6[])] <- cover_pca$li$Axis6

writeRaster(pc1, "/home/javifl/IREC/master_david/variables_raster/pc1.tif")
writeRaster(pc2, "/home/javifl/IREC/master_david/variables_raster/pc2.tif")
writeRaster(pc3, "/home/javifl/IREC/master_david/variables_raster/pc3.tif")
writeRaster(pc4, "/home/javifl/IREC/master_david/variables_raster/pc4.tif")
writeRaster(pc5, "/home/javifl/IREC/master_david/variables_raster/pc5.tif")
writeRaster(pc6, "/home/javifl/IREC/master_david/variables_raster/pc6.tif")

################################################################################
# 18/06/2020

library(lattice)
library(parallel)
library(Rcpp)
library(unmarked)

data <- read.csv("/home/javifl/IREC/master_david/david_17062020/datos_2015limpio5.csv")

datagamo <- data[data$sp == "gamo",]
datagamo <- na.omit(datagamo)#eliminanos los valores nulos
head(datagamo)

y <- datagamo[,5:14]
n <- nrow(datagamo)

d2015.site <- data.frame(scale(datagamo[,c(15,16,23,24,25,26,27)]))

time <- as.factor(rep(c(1:10),n))
d2015.obs <- data.frame(time)

d2015 <- unmarkedFrameOccu(y = y, siteCovs = d2015.site,obsCovs=d2015.obs)

head(d2015)
summary(d2015)

fm1<-occu(~1 ~dvera, d2015)
fm2<-occu(~time ~dvera, d2015)
fm3<-occu(~1 ~dwat, d2015)
fm4<-occu(~time ~dwat, d2015)
fm5<-occu(~1 ~pc1+pc2+pc3+pc4+pc5, d2015)
fm6<-occu(~time ~pc1+pc2+pc3+pc4+pc5, d2015)
fm7<-occu(~1 ~dvera + dwat, d2015)
fm8<-occu(~time ~dvera + dwat, d2015)
fm9<-occu(~1 ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015)
fm10<-occu(~time ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015)

fm1;fm2;fm3;fm4;fm5;fm6;fm7;fm8;fm9;fm10
fmlist<-fitList(primero=fm1, segundo=fm2, tercero=fm3, cuarto=fm4, quinto=fm5, sexto=fm6, septimo=fm7, octavo=fm8, noveno=fm9, decimo=fm10)
modSel(fmlist)

summary(fm9)
