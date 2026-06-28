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
# set working directory
setwd("C:/Users/...")

# import data
data <- read.table("Sal_Temp_Nb-clams - 3D.txt", header = TRUE)

# data cleaning
data <- data %>%
  mutate(
    Temperature = as.numeric(gsub(",", ".", Temperature)),
    Predation = as.numeric(gsub(",", ".", Predation)),
    Salinity = as.factor(Salinity)
  ) %>%
  drop_na() %>%
  filter(Temperature >= 14)

# summary statistics
summary_data <- data %>%
  group_by(Salinity, Temperature) %>%
  summarise(
    n = n(),
    mean = mean(Predation),
    se = sd(Predation) / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  )

# R² calculation
calc_r2 <- function(model, df) {
  pred <- predict(model, newdata = df)
  1 - sum((df$mean - pred)^2) / sum((df$mean - mean(df$mean))^2)
}

# model fitting depending on salinity
fit_model_by_salinity <- function(df, s) {

  if(s == "5") {

    model <- nlsLM(
      mean ~ A * (1 - exp(-k * (Temperature - Tcrit))),
      data = df
    )

  } else if(s == "18") {

    model <- nlsLM(
      mean ~ y0 + a * exp(b * Temperature),
      data = df %>% filter(Temperature != 24)
    )

  } else if(s == "25") {

    model <- nlsLM(
      mean ~ y0 + a * (1 - exp(-b * Temperature)),
      data = df
    )

  } else if(s %in% c("35", "45")) {

    model <- lm(
      mean ~ Temperature + I(Temperature^2),
      data = df
    )
  }

  list(
    model = model,
    AICc = AICc(model),
    R2 = calc_r2(model, df)
  )
}

# plotting
plots <- list()

for(s in levels(summary_data$Salinity)) {

  df <- summary_data %>% filter(Salinity == s)

  fit <- fit_model_by_salinity(df, s)

  newdata <- data.frame(
    Temperature = seq(16, 32, length.out = 2000)
  )

  newdata$pred <- predict(fit$model, newdata = newdata)

  p <- ggplot(df, aes(x = Temperature, y = mean)) +

    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "lightblue",
      alpha = 0.6
    ) +

    geom_line(color = "blue", linewidth = 1.2) +
    geom_point(size = 2.6, color = "blue") +

    geom_line(
      data = newdata,
      aes(x = Temperature, y = pred),
      color = "black",
      linewidth = 1.4
    ) +

    labs(
      title = paste("Salinity =", s),
      x = "Temperature (°C)",
      y = "Predation"
    ) +

    coord_cartesian(ylim = c(0, 16)) +

    theme_classic(base_size = 14)

  plots[[s]] <- p
}

# final figure
final_plot <- wrap_plots(plots, ncol = 2)

final_plot
```

![Figure4](Figure4.png)



Figure 4. Effects of temperature and salinity on the predation activity of the blue crab Callinectes sapidus expressed as the number of clams (Ruditapes philippinarum) consumed per individual per day across (A) thermal conditions for each salinity treatment (5, 18, 25, 35 and 45 psu) and (B) salinity conditions for each temperature treatment (16, 20, 24, 28 and 32 °C). Black lines represent observed mean values, grey areas indicate confidence intervals, and black/red lines correspond to fitted model trends. Predation generally increased with temperature under low to intermediate salinities, whereas high salinity conditions reduced feeding activity, particularly at elevated temperatures. Regression parameters are globalized in the Supplementary Table 1.



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



![Figure5](Figure5.png)


Figure 5. Three-dimensional predicted response surface showing the combined effects of temperature (°C) and salinity on the number of clams (Ruditapes philippinarum) consumed per individual blue crab (Callinectes sapidus) per day (clams eaten ind⁻¹ d⁻¹). Colors indicate feeding intensity (number of clams eaten ind-1 day-1), ranging from low values (blue) to high values (yellow/orange). The black dotted lines represent the temperatures and salinities experimentally measured in our study.



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



