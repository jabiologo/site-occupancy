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

# No estoy seguro de si sería necesario estandarizar las variables de los PCAs
d2015.site <- data.frame(scale(datagamo[,c(15,16,23,24,25,26,27)])) 

time <- as.factor(rep(c(1:10),n))
d2015.obs <- data.frame(time)

d2015 <- unmarkedFrameOccu(y = y, siteCovs = d2015.site, obsCovs=d2015.obs)

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

summary(fm9) # Salen NaNs y estimates muy rarunos

d2015c <- unmarkedFramePCount(y = y, siteCovs = d2015.site,obsCovs=d2015.obs)

fm1<-pcount(~1 ~dvera, d2015c, K = 150)
fm2<-pcount(~time ~dvera, d2015c, K = 150)
fm3<-pcount(~1 ~dwat, d2015c, K = 150)
fm4<-pcount(~time ~dwat, d2015c, K = 150)
fm5<-pcount(~1 ~pc1+pc2+pc3+pc4+pc5, d2015c, K = 150)
fm6<-pcount(~time ~pc1+pc2+pc3+pc4+pc5, d2015c, K = 150)
fm7<-pcount(~1 ~dvera + dwat, d2015c, K = 150)
fm8<-pcount(~time ~dvera + dwat, d2015c, K= 150)
fm9<-pcount(~1 ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm10<-pcount(~time ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)

fm1;fm2;fm3;fm4;fm5;fm6;fm7;fm8;fm9;fm10
fmlist<-fitList(primero=fm1, segundo=fm2, tercero=fm3, cuarto=fm4, quinto=fm5, sexto=fm6, septimo=fm7, octavo=fm8, noveno=fm9, decimo=fm10)
modSel(fmlist)

summary(fm9) # Sale warining y estimates demasiado altos


################################################################################

dataciervo <- data[data$sp == "ciervo",]
dataciervo <- na.omit(dataciervo)
head(dataciervo)


y <- dataciervo[,3:14]
n <- nrow(dataciervo)


d2015.site <- data.frame(scale(dataciervo[,c(15,16,23,24,25,26,27)]))


time <- as.factor(rep(c(1:12),n))
d2015.obs <- data.frame(time)
d2015c <- unmarkedFramePCount(y = y, siteCovs = d2015.site,obsCovs=d2015.obs)

head(d2015c)
summary(d2015c)

fm1<-pcount(~1 ~dvera, d2015c, K = 150)
fm2<-pcount(~time ~dvera, d2015c, K = 150)
fm3<-pcount(~1 ~dwat, d2015c, K = 150)
fm4<-pcount(~time ~dwat, d2015c, K = 150)
fm5<-pcount(~1 ~pc1+pc2+pc3+pc4+pc5, d2015c, K = 150)
fm6<-pcount(~time ~pc1+pc2+pc3+pc4+pc5, d2015c, K = 150)
fm7<-pcount(~1 ~dvera + dwat, d2015c, K = 150)
fm8<-pcount(~time ~dvera + dwat, d2015c, K= 150)
fm9<-pcount(~1 ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm10<-pcount(~time ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)

fm1;fm2;fm3;fm4;fm5;fm6;fm7;fm8;fm9;fm10
fmlist<-fitList(primero=fm1, segundo=fm2, tercero=fm3, cuarto=fm4, quinto=fm5, sexto=fm6, septimo=fm7, octavo=fm8, noveno=fm9, decimo=fm10)
modSel(fmlist)

summary(fm10) #(este es el modelo que miramos por ser el m?s ajustado, y del que obtenemos los estimates)

# Vemos como están centradas y aplicamos la fórmula
scale(dataciervo[,c(15,16,23,24,25,26,27)])

library(raster)
variables <- stack(list.files(path="/home/javifl/IREC/master_david/variables_raster",pattern='*.tif', full.names=TRUE))

variables$dvera <- (variables$dvera - 2113.52868026) / 2182.1243541
variables$dwat <- (variables$dwat - 737.95705705) / 451.7157191
variables$pc1 <- (variables$pc1 - (-0.36966496)) / 1.1953509
variables$pc2 <- (variables$pc2 - 0.65164198) / 0.7663623
variables$pc3 <- (variables$pc3 - 0.47492335) / 1.2962557
variables$pc4 <- (variables$pc4 - (-0.08217354)) / 0.8082859
variables$pc5 <- (variables$pc5 - 0.06771092) / 1.3015840

fm10_pred <- predict(fm10, newdata=variables, type = "state")

plot(fm10_pred$Predicted, axes=FALSE)

################################################################################

data <- read.csv("/home/javifl/IREC/master_david/david_19062020/datos_2015limpio7.csv")

dataciervo <- data[data$sp == "ciervo",]
head(dataciervo)

y <- dataciervo[,3:20]
n <- nrow(dataciervo)

# No estoy seguro de si sería necesario estandarizar las variables de los PCAs
d2015.site <- data.frame(scale(dataciervo[,c(21,22,29:33)])) 

d2015.obs <- data.frame(as.factor(as.vector(t(dataciervo[,35:52]))))
names(d2015.obs) <- "time"

d2015c <- unmarkedFramePCount(y = y, siteCovs = d2015.site, obsCovs=d2015.obs)

fm1<-pcount(~1 ~dvera, d2015c, K = 150)
fm2<-pcount(~time ~dvera, d2015c, K = 150)
fm3<-pcount(~1 ~dwat, d2015c, K = 150)
fm4<-pcount(~time ~dwat, d2015c, K = 150)
fm5<-pcount(~1 ~pc1+pc2+pc3+pc4+pc5, d2015c, K = 150)
fm6<-pcount(~time ~pc1+pc2+pc3+pc4+pc5, d2015c, K = 150)
fm7<-pcount(~1 ~dvera + dwat, d2015c, K = 150)
fm8<-pcount(~time ~dvera + dwat, d2015c, K= 150)
fm9<-pcount(~1 ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm10<-pcount(~time ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)
#fm11<-pcount(~time ~v1+v2+v3+v4+v5+v6+dvera+dwat, d2015c, K = 150)



scale(dataciervo[,c(21,22,29:33)])

library(raster)
variables <- stack(list.files(path="/home/javifl/IREC/master_david/variables_raster",pattern='*.tif', full.names=TRUE))

variables$dvera <- (variables$dvera - 1869.998) / 2008.5392896
variables$dwat <- (variables$dwat - 617.3486) / 408.5496179
variables$pc1 <- (variables$pc1 - 0.1626645) / 1.2750652
variables$pc2 <- (variables$pc2 - 0.5680776) / 0.6207696
variables$pc3 <- (variables$pc3 - 0.4596181) / 1.1003653
variables$pc4 <- (variables$pc4 - (-0.1316229)) / 0.6793139
variables$pc5 <- (variables$pc5 - 0.003784568) / 0.9570490

fm10_pred <- predict(fm10, newdata=variables, type = "state")

plot(fm10_pred$Predicted, axes=FALSE)

# Cn las V normales

y <- dataciervo[,3:20]
n <- nrow(dataciervo)

# No estoy seguro de si sería necesario estandarizar las variables de los PCAs
d2015.site <- data.frame(scale(dataciervo[,c(21:28)])) 

d2015.obs <- data.frame(as.factor(as.vector(t(dataciervo[,35:52]))))
names(d2015.obs) <- "time"

d2015c <- unmarkedFramePCount(y = y, siteCovs = d2015.site, obsCovs=d2015.obs)

fm11<-pcount(~time ~v1+v2+v3+v4+v5+v6+dvera+dwat, d2015c, K = 150)

scale(dataciervo[,c(21:28)])

variables <- stack(list.files(path="/home/javifl/IREC/master_david/variables_raster",pattern='*.tif', full.names=TRUE))

variables$dvera <- (variables$dvera - 1869.998) / 2008.5392896
variables$dwat <- (variables$dwat - 617.3486) / 408.5496179
variables$v1 <- (variables$v1 - 35.235294) / 50.53235
variables$v2 <- (variables$v2 - 58.647059) / 62.74425
variables$v3 <- (variables$v3 - 35.000000) / 53.19717
variables$v4 <- (variables$v4 - 2.382353) / 13.89139
variables$v5 <- (variables$v5 - 5.676471) / 12.85047
variables$v6 <- (variables$v6 - 18.088235) / 37.56398

fm11_pred <- predict(fm11, newdata=variables, type = "state")

plot(fm11_pred$Predicted, axes=FALSE)

plot(fm11_pred$Predicted[fm11_pred$Predicted<1000])

ss <- fm11_pred$Predicted

#mal
# ~pc1+pc2+pc4+pc5+dvera+dwat
# ~pc1+pc2+pc3+pc4+dvera+dwat


fm10<-pcount(~time ~pc1+pc2+pc3+pc4+dvera+dwat, d2015c, K = 150)
fm27_pred <- predict(fm19, newdata=variables, type = "state")
plot(fm27_pred$Predicted, axes=FALSE)


fm11<-pcount(~time ~pc1+dvera+dwat, d2015c, K = 150)
fm12<-pcount(~time ~pc2+dvera+dwat, d2015c, K = 150)
fm13<-pcount(~time ~pc3+dvera+dwat, d2015c, K = 150)
fm14<-pcount(~time ~pc4+dvera+dwat, d2015c, K = 150)
fm15<-pcount(~time ~pc5+dvera+dwat, d2015c, K = 150)
fm16<-pcount(~time ~pc1+pc2+dvera+dwat, d2015c, K = 150)
fm17<-pcount(~time ~pc1+pc3+dvera+dwat, d2015c, K = 150)
fm18<-pcount(~time ~pc1+pc4+dvera+dwat, d2015c, K = 150)
fm19<-pcount(~time ~pc1+pc5+dvera+dwat, d2015c, K = 150)
fm20<-pcount(~time ~pc2+pc3+dvera+dwat, d2015c, K = 150)
fm21<-pcount(~time ~pc2+pc4+dvera+dwat, d2015c, K = 150)
fm22<-pcount(~time ~pc2+pc5+dvera+dwat, d2015c, K = 150)
fm23<-pcount(~time ~pc3+pc4+dvera+dwat, d2015c, K = 150)
fm24<-pcount(~time ~pc3+pc5+dvera+dwat, d2015c, K = 150)
fm25<-pcount(~time ~pc4+pc5+dvera+dwat, d2015c, K = 150)
fm26<-pcount(~time ~pc1+pc2+pc3+dvera+dwat, d2015c, K = 150)
fm27<-pcount(~time ~pc1+pc2+pc4+dvera+dwat, d2015c, K = 150)
fm28<-pcount(~time ~pc1+pc2+pc5+dvera+dwat, d2015c, K = 150)
fm29<-pcount(~time ~pc1+pc3+pc4+dvera+dwat, d2015c, K = 150)
fm30<-pcount(~time ~pc1+pc3+pc5+dvera+dwat, d2015c, K = 150)
fm31<-pcount(~time ~pc1+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm32<-pcount(~time ~pc2+pc3+pc4+dvera+dwat, d2015c, K = 150)
fm33<-pcount(~time ~pc2+pc3+pc5+dvera+dwat, d2015c, K = 150)
fm34<-pcount(~time ~pc2+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm35<-pcount(~time ~pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm36<-pcount(~time ~pc1+pc2+pc3+pc4+dvera+dwat, d2015c, K = 150)
fm37<-pcount(~time ~pc1+pc2+pc3+pc5+dvera+dwat, d2015c, K = 150)
fm38<-pcount(~time ~pc1+pc2+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm39<-pcount(~time ~pc1+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm40<-pcount(~time ~pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm41<-pcount(~time ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)

fm42<-pcount(~time ~pc1+dvera, d2015c, K = 150)
fm43<-pcount(~time ~pc2+dvera, d2015c, K = 150)
fm44<-pcount(~time ~pc3+dvera, d2015c, K = 150)
fm45<-pcount(~time ~pc4+dvera, d2015c, K = 150)
fm46<-pcount(~time ~pc5+dvera, d2015c, K = 150)
fm47<-pcount(~time ~pc1+pc2+dvera, d2015c, K = 150)
fm48<-pcount(~time ~pc1+pc3+dvera, d2015c, K = 150)
fm49<-pcount(~time ~pc1+pc4+dvera, d2015c, K = 150)
fm50<-pcount(~time ~pc1+pc5+dvera, d2015c, K = 150)
fm51<-pcount(~time ~pc2+pc3+dvera, d2015c, K = 150)
fm52<-pcount(~time ~pc2+pc4+dvera, d2015c, K = 150)
fm53<-pcount(~time ~pc2+pc5+dvera, d2015c, K = 150)
fm54<-pcount(~time ~pc3+pc4+dvera, d2015c, K = 150)
fm55<-pcount(~time ~pc3+pc5+dvera, d2015c, K = 150)
fm56<-pcount(~time ~pc4+pc5+dvera, d2015c, K = 150)
fm57<-pcount(~time ~pc1+pc2+pc3+dvera, d2015c, K = 150)
fm58<-pcount(~time ~pc1+pc2+pc4+dvera, d2015c, K = 150)
fm59<-pcount(~time ~pc1+pc2+pc5+dvera, d2015c, K = 150)
fm60<-pcount(~time ~pc1+pc3+pc4+dvera, d2015c, K = 150)
fm61<-pcount(~time ~pc1+pc3+pc5+dvera, d2015c, K = 150)
fm62<-pcount(~time ~pc1+pc4+pc5+dvera, d2015c, K = 150)
fm63<-pcount(~time ~pc2+pc3+pc4+dvera, d2015c, K = 150)
fm64<-pcount(~time ~pc2+pc3+pc5+dvera, d2015c, K = 150)
fm65<-pcount(~time ~pc2+pc4+pc5+dvera, d2015c, K = 150)
fm66<-pcount(~time ~pc3+pc4+pc5+dvera, d2015c, K = 150)
fm67<-pcount(~time ~pc1+pc2+pc3+pc4+dvera, d2015c, K = 150)
fm68<-pcount(~time ~pc1+pc2+pc3+pc5+dvera, d2015c, K = 150)
fm69<-pcount(~time ~pc1+pc2+pc4+pc5+dvera, d2015c, K = 150)
fm70<-pcount(~time ~pc1+pc3+pc4+pc5+dvera, d2015c, K = 150)
fm71<-pcount(~time ~pc2+pc3+pc4+pc5+dvera, d2015c, K = 150)
fm72<-pcount(~time ~pc1+pc2+pc3+pc4+pc5+dvera, d2015c, K = 150)

fm73<-pcount(~time ~pc1+dwat, d2015c, K = 150)
fm74<-pcount(~time ~pc2+dwat, d2015c, K = 150)
fm75<-pcount(~time ~pc3+dwat, d2015c, K = 150)
fm76<-pcount(~time ~pc4+dwat, d2015c, K = 150)
fm77<-pcount(~time ~pc5+dwat, d2015c, K = 150)
fm78<-pcount(~time ~pc1+pc2+dwat, d2015c, K = 150)
fm79<-pcount(~time ~pc1+pc3+dwat, d2015c, K = 150)
fm80<-pcount(~time ~pc1+pc4+dwat, d2015c, K = 150)
fm81<-pcount(~time ~pc1+pc5+dwat, d2015c, K = 150)
fm82<-pcount(~time ~pc2+pc3+dwat, d2015c, K = 150)
fm83<-pcount(~time ~pc2+pc4+dwat, d2015c, K = 150)
fm84<-pcount(~time ~pc2+pc5+dwat, d2015c, K = 150)
fm85<-pcount(~time ~pc3+pc4+dwat, d2015c, K = 150)
fm86<-pcount(~time ~pc3+pc5+dwat, d2015c, K = 150)
fm87<-pcount(~time ~pc4+pc5+dwat, d2015c, K = 150)
fm88<-pcount(~time ~pc1+pc2+pc3+dwat, d2015c, K = 150)
fm89<-pcount(~time ~pc1+pc2+pc4+dwat, d2015c, K = 150)
fm90<-pcount(~time ~pc1+pc2+pc5+dwat, d2015c, K = 150)
fm91<-pcount(~time ~pc1+pc3+pc4+dwat, d2015c, K = 150)
fm92<-pcount(~time ~pc1+pc3+pc5+dwat, d2015c, K = 150)
fm93<-pcount(~time ~pc1+pc4+pc5+dwat, d2015c, K = 150)
fm94<-pcount(~time ~pc2+pc3+pc4+dwat, d2015c, K = 150)
fm95<-pcount(~time ~pc2+pc3+pc5+dwat, d2015c, K = 150)
fm96<-pcount(~time ~pc2+pc4+pc5+dwat, d2015c, K = 150)
fm97<-pcount(~time ~pc3+pc4+pc5+dwat, d2015c, K = 150)
fm98<-pcount(~time ~pc1+pc2+pc3+pc4+dwat, d2015c, K = 150)
fm99<-pcount(~time ~pc1+pc2+pc3+pc5+dwat, d2015c, K = 150)
fm100<-pcount(~time ~pc1+pc2+pc4+pc5+dwat, d2015c, K = 150)
fm101<-pcount(~time ~pc1+pc3+pc4+pc5+dwat, d2015c, K = 150)
fm102<-pcount(~time ~pc2+pc3+pc4+pc5+dwat, d2015c, K = 150)
fm103<-pcount(~time ~pc1+pc2+pc3+pc4+pc5+dwat, d2015c, K = 150)


fmlist<-fitList(m11=fm11, m12=fm12, m13=fm13, m14=fm14, m15=fm15, 
                m16=fm16, m17=fm17, m18=fm18, m19=fm19, m20=fm20,
                m21=fm21, m22=fm22, m23=fm23, m24=fm24, m25=fm25,
                m26=fm26, m27=fm27, m28=fm28, m29=fm29, m30=fm30,
                m31=fm31, m32=fm32, m33=fm33, m34=fm34, m35=fm35,
                m36=fm36, m37=fm37, m38=fm38, m39=fm39, m40=fm40,
                m41=fm41, m42=fm42, m43=fm43, m44=fm44, m45=fm45,
                m46=fm46, m47=fm47, m48=fm48, m49=fm49, m50=fm50,
                m51=fm51, m52=fm52, m53=fm53, m54=fm54, m55=fm55,
                m56=fm56, m57=fm57, m58=fm58, m59=fm59, m60=fm60,
                m61=fm61, m62=fm62, m63=fm63, m64=fm64, m65=fm65,
                m66=fm66, m67=fm67, m68=fm68, m69=fm69, m70=fm70,
                m71=fm71, m72=fm72, m73=fm73, m74=fm74, m75=fm75,
                m76=fm76, m77=fm77, m78=fm78, m79=fm79, m80=fm80,
                m81=fm81, m82=fm82, m83=fm83, m84=fm84, m85=fm85,
                m86=fm86, m87=fm87, m88=fm88, m89=fm89, m90=fm90,
                m91=fm91, m92=fm92, m93=fm93, m94=fm94, m95=fm95,
                m96=fm96, m97=fm97, m98=fm98, m99=fm99, m100=fm100,
                m101=fm101, m102=fm102, m103=fm103)








modSel(fmlist)


fm27<-pcount(~time ~pc1+pc2+    pc4+    dvera+dwat, d2015c, K = 150)
fm19<-pcount(~time ~pc1+            pc5+dvera+dwat, d2015c, K = 150)
fm29<-pcount(~time ~pc1+    pc3+pc4+    dvera+dwat, d2015c, K = 150)
fm11<-pcount(~time ~pc1+                dvera+dwat, d2015c, K = 150)
fm36<-pcount(~time ~pc1+pc2+pc3+pc4+    dvera+dwat, d2015c, K = 150)
fm41<-pcount(~time ~pc1+pc2+pc3+pc4+pc5+dvera+dwat, d2015c, K = 150)
fm16<-pcount(~time ~pc1+pc2+            dvera+dwat, d2015c, K = 150)
fm38<-pcount(~time ~pc1+pc2+    pc4+pc5+dvera+dwat, d2015c, K = 150)


library(corrplot)
vv <- data.frame(cbind(na.omit(variables$dvera[]),na.omit(variables$dwat[]),na.omit(variables$pc1[]), 
            na.omit(variables$pc2[]),na.omit(variables$pc3[]),na.omit(variables$pc4[]),na.omit(variables$pc5[])))
names(vv) <- c("dvera", "dwat", "pc1", "pc2", "pc3", "pc4", "pc5")
corrplot(cor(vv), "number")
corrplot.mixed((cor(d2015.site)), upper = "ellipse")
corrplot.mixed((cor(vv)), upper = "ellipse")



################################################################################
# 30/06/2020
# Reclasificación de un raster por quantiles y evaluación mediante Kappa de Cohen

library(raster)
library(psych)

# Cargamos los rasters
foto <- raster("/home/javifl/IREC/master_david/kappa/ciervofoto_pred.tif")
tele <- raster("/home/javifl/IREC/master_david/kappa/ciervotele_pred.tif")

# Calculamos los cuantiles y hacemos una matriz con las categorizaciones según los cuantiles
quantile(foto)
quan <- matrix(c(0.0,  1.940705299, 1,
          1.940705299,  6.142425299, 2,
          6.142425299, 17.560887337, 3,
          17.560887337, 77.89667, 4), 4,3, byrow = TRUE)

# Reclasificamos el mapa, lo categorizamos en 4 y lo ploteamos para ver que pinta tiene:
foto4 <- reclassify(foto, quan)
plot(foto4)

# Hacemos lo mismo con el siguiente mapa
quantile(tele)
quan <- matrix(c(0.0,  0.4137109, 1,
                 0.4137109,  0.4874967, 2,
                 0.4874967, 0.5519628, 3,
                 0.5519628, 0.71, 4), 4,3, byrow = TRUE)

tele4 <- reclassify(tele, quan)
plot(tele4)

# Calculamos el kappa de Cohen con el paquete psych. Nos quedamos con el weighted
# Los valores van desde -1 (concordancia inversa) hasta 1 (concordancia total). Un valor cercano
# a cero indicaría una categorización al azar
cohen.kappa(x=cbind(na.omit(foto4[]),na.omit(tele4[])))
