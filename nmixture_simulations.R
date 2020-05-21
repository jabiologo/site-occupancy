library(raster)
library(dismo)

sarea <- raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10, ymn = 0, ymx = 10)  # study area
lambda <- 4
sarea[] <- rpois(100, lambda)  
plot(sarea)
sum(sarea[])

site_ID <- sample(1:100, 15)

plot(sarea)
points(xyFromCell(sarea,site_ID))

dataset <- data.frame(site_ID)
dataset$trueN <- extract(sarea,site_ID)
dataset$O1 <- NA
dataset$O2 <- NA
dataset$O3 <- NA
dataset$O4 <- NA

for (j in 1:4){
  for (i in 1:length(site_ID)){
    dataset[i,j+2] <- rbinom(1, extract(sarea,site_ID[i]), 0.4)
  }
}

dataset
dataset$maxN <- apply(dataset[,3:6],1, max)

plot(dataset$trueN, dataset$maxN)
abline(0,1)

dataUM <- unmarkedFramePCount(y = dataset[,3:6])

m1 <- pcount(~1 ~ 1, data=dataUM, K=50) 
summary(m1)

# lambda estimada
exp(coef(m1)[1])

# p estimada (LOG ODDS SCALE = LOGIT)
exp(coef(m1)[2])/(1+exp(coef(m1)[2]))
plogis(coef(m1)[2]) # lo mismo que arriba

# predicciones de N
N_m1<- bup(ranef(m1)) # s4 class
plot(dataset$trueN,N_m1, xlab="True density", ylab="Predicted density")
abline(0,1)# a 1:1 line





# Crear una variable predictora
dwat <- scale(distanceFromPoints(sarea, c(3.5,3.5))) + 2.28668

# Ahora, la lambda de la Poisson viene determinada por esta variable en la forma:
lambda <- exp(3 - 0.9*(dwat))

sarea <- raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10, ymn = 0, ymx = 10)  # study area
sarea[1:ncell(sarea)] <- rpois(1, lambda[1:ncell(sarea)])

for (i in 1:ncell(sarea)){
  sarea[i] <- rpois(1, lambda[i])
}
# N total
sum(sarea[])
# vemos la relación entre la abundancia y la distancia al punto de agua
plot(dwat[], sarea[], ylab="Abundancia", xlab="dWat",las=1)

# Ahora simulamos 15 muestreos durante 4 ocasiones. Vamos a consderar que la probabilidad
# de detección no depende de la ocasión, sino que también depende de la distancia al agua.
# Por ejemplo, cuanta mayor distancia al agua, mayor detectabilidad porque hay menos vegetación que 
# nos impide detectar la especie en cuestión... esta relación tiene la forma:

# logit(p) = gamma_0 + gamma_1 * (dwat)
gamma_0 <- 0.2
gamma_1 <- 0.8

# Para despejar logit() de la izquierda, pasamos a la derecha todo este churro
p <- exp(gamma_0 + gamma_1 * dwat) / (1 + exp(gamma_0 + gamma_1 * dwat))

# Vemos la probabilidad de detección en nuestro área de estudio
plot(p)

# Y la relación entre probabilidad de detección y distancia al agua
plot(dwat[], p[])

# Ya estamos listos para simular un muestreo teninendo en cuenta la detectabilidad
site_ID <- sample(1:100, 15)

plot(sarea)
points(xyFromCell(sarea,site_ID))

dataset <- data.frame(site_ID)
dataset$trueN <- extract(sarea,site_ID)
dataset$O1 <- NA
dataset$O2 <- NA
dataset$O3 <- NA
dataset$O4 <- NA

for (j in 1:4){
  for (i in 1:length(site_ID)){
    dataset[i,j+2] <- rbinom(1, extract(sarea,site_ID[i]), extract(p,site_ID[i]))
  }
}

dataset
dataset$maxN <- apply(dataset[,3:6],1, max)

plot(dataset$trueN, dataset$maxN)
abline(0,1)

# Ahora preparamos los datos para ajustar el modelo. Como vamos a usar como predictor
# la distancia al agua, incluimos en nuestro dataset los valores de distancia al agua
# para cada uno de nuestros 15 puntos de muestreo

dataset$dwat <- extract(dwat, site_ID)
dataUM <- unmarkedFramePCount(y = dataset[,3:6], siteCovs=data.frame(dwat=dataset$dwat))

m2 <- pcount(~dwat ~dwat, data=dataUM, K=150)
m2

# predicciones de N
N_m2<- bup(ranef(m2))
plot(dataset$trueN, N_m2, xlab="True density", ylab="Predicted density")
abline(0,1)# a 1:1 line

# lambda estimada
exp(coef(m2)[1:2])

# p estimada (LOG ODDS SCALE = LOGIT)
exp(coef(m2)[3:4])/(1+exp(coef(m2)[3:4]))
plogis(coef(m2)[3:4]) # lo mismo que arriba

sdwat <- stack(dwat)
names(sdwat) <- "dwat"

pred <- predict(m2, newdata=sdwat, type = "state")




