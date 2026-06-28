# Journal of Animal Ecology

his R scripts come with the following research article:

\*\*“Temperature and salinity jointly shape predation by the invasive blue crab \*Callinectes sapidus\* in the Berre Lagoon (France)”\*\* (\*Chkili et al., submitted to Journal of Animal Ecology\*).



It includes (i) regression analyses to estimate temperature-dependent predation responses across salinity treatments, (ii) Generalized Additive Models (GAMs) and polynomial models to predict predation rates in three-dimensional environmental space, (iii) kriging procedures to estimate spatial and seasonal predation intensity in the Berre Lagoon, and (iv) niche analyses identifying environmental no-predation conditions.



\---

## Experimental datasets



[Predation number of clams](Sal_Temp_Nb-clams%20-%203D.txt) <br>



\[Predation biomass percentage](data/Sal\_Temp\_pourcentage\_3D.txt) <br>


[Shell breakage dataset](coquilles.txt) <br>



## Spatial datasets



[Berre Lagoon predation estimation](estimation\_predation\_Berre.txt) <br>



[Berre Lagoon shapefile](couche\_Berre.shp)





## Load required libraries

```R

library(tidyverse)

library(minpack.lm)

library(AICcmodavg)

library(mgcv)

library(plot3D)

library(plotly)

library(readr)

library(dplyr)

library(lubridate)

library(sf)

library(gstat)

library(sp)

library(ggplot2)

library(viridis)

library(fields)

```



\---

## Calculation of temperature-dependent predation regressions



In our study, independent regressions were fitted for each salinity treatment to characterize thermal responses of predation and identify threshold effects.



```R

\# set working directory

setwd("C:/Users/...")



\# import data

data <- read.table("Sal\_Temp\_Nb-clams - 3D.txt", header = TRUE)



data <- data %>%

&#x20; mutate(

&#x20;   Temperature = as.numeric(gsub(",", ".", Temperature)),

&#x20;   Predation = as.numeric(gsub(",", ".", Predation)),

&#x20;   Salinity = as.factor(Salinity)

&#x20; ) %>%

&#x20; drop\_na()



\# summarise by salinity × temperature

summary\_data <- data %>%

&#x20; group\_by(Salinity, Temperature) %>%

&#x20; summarise(

&#x20;   n = n(),

&#x20;   mean = mean(Predation),

&#x20;   se = sd(Predation)/sqrt(n)

&#x20; )



\# fit regressions

fit\_model\_by\_salinity(summary\_data, "5")

fit\_model\_by\_salinity(summary\_data, "18")

fit\_model\_by\_salinity(summary\_data, "25")

fit\_model\_by\_salinity(summary\_data, "35")

fit\_model\_by\_salinity(summary\_data, "45")

```



![Figure3](figures/Figure.3.png)



Figure 3. Temperature-dependent predation responses under different salinity treatments.



\---

## Construction of the 3D predation surface



A GAM model was used to estimate predation intensity as a function of salinity and temperature.



```R

data <- read.delim("Sal\_Temp\_Nb-clams - 3D.txt")



mod\_gam <- gam(

&#x20; Predation \~ s(Salinity, Temperature, k = 10),

&#x20; data = data,

&#x20; method = "REML"

)



sal\_grid <- seq(min(data$Salinity), max(data$Salinity), length.out = 50)

temp\_grid <- seq(min(data$Temperature), max(data$Temperature), length.out = 50)



grid <- expand.grid(

&#x20; Salinity = sal\_grid,

&#x20; Temperature = temp\_grid

)



grid$Predation\_pred <- predict(mod\_gam, newdata = grid)

```



!\[Figure4A](figures/Figure4A.png)



Figure 4A. Three-dimensional predictive surface of predation rate.



\---



\## Construction of biomass consumption surface



Predation expressed as percentage of crab body mass consumed.



```R

data <- read.delim("Sal\_Temp\_pourcentage - 3D.txt")



mod\_gam <- gam(

&#x20; Predation \~ s(Salinity, Temperature, k = 10),

&#x20; data = data,

&#x20; method = "REML"

)



grid$Predation\_pred <- predict(mod\_gam, newdata = grid)

```



!\[Figure4B](figures/Figure4B.png)



Figure 4B. Three-dimensional predictive surface of biomass consumption.



\---



\## Spatial kriging analyses



Spatial interpolation was performed to estimate monthly predation pressure across the Berre Lagoon.



```R

df <- read.delim("etimation\_predation\_Berre.txt")



lagune <- st\_read("couche\_Berre.shp")



krig <- krige(

&#x20; pred \~ 1,

&#x20; locations = points\_sp,

&#x20; newdata = grille\_sp,

&#x20; model = vgm\_model

)

```



!\[Figure5](figures/Figure5.png)



Figure 5. Monthly kriging maps showing spatial predation hotspots.



\---



\## No-predation niche



Kriging interpolation was used to identify environmental conditions where predation is absent.



```R

df <- read\_delim("coquilles.txt")



vgm\_emp <- variogram(blue \~ 1, df)



krig <- krige(

&#x20; blue \~ 1,

&#x20; locations = df,

&#x20; newdata = grid\_df,

&#x20; model = vgm\_fit

)

```



!\[Figure6](figures/Figure6.png)



Figure 6. No-predation niche under low temperature and varying salinity.



\---



\## Contact



\*\*Guillaume Marchessaux\*\*

Institut Méditerranéen d’Océanologie (MIO), Marseille, France

Corresponding author: \[guillaume.marchessaux@ird.fr](mailto:guillaume.marchessaux@ird.fr)



