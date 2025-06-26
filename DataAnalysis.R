
# Data analysis radial growth Diaporthe amygdali
# Created by Luca Rossini on 11 December 2024
# Last update 25 June 2025
# E-mail: luca.rossini@ulb.be


# ANALYSIS: mycelium extension rate of the three strains at different 
# temperatures


# Acquisition of the data - File 'DatasetMicelio.csv'

library(readxl)

filename_mycelium <- file.choose()

data_mycelium = read_excel(filename_mycelium, sheet = 'Elaboration', 
                           col_names = T)

head(data_mycelium)

isolate <- as.factor(data_mycelium$Isolate)
tempCode <- as.factor(data_mycelium$TempCode)
tempValue <- as.numeric(data_mycelium$Temperature)
replicates <- as.factor(data_mycelium$Replicates)
length <- as.numeric(data_mycelium$Length_1_Day)


# Check the levels

levels(isolate)
levels(tempCode)
levels(replicates)


# GLM: 'lesion' variable and twig as random effect, negative binomial as 
# distribution.

library(lme4)

GenLin_Isolates <- glmer.nb(length ~ isolate + (1 | replicates), 
                            data=data_mycelium)

summary(GenLin_Isolates)

# Check the distribution of the residuals

qqnorm(residuals(GenLin_Isolates))
qqline(residuals(GenLin_Isolates))

shapiro.test(residuals(GenLin_Isolates))

# Checks with DHARMa

library(DHARMa)

testDispersion(GenLin_Isolates)
simulationOutput_Isolates <- simulateResiduals(fittedModel = GenLin_Isolates, plot = T)
residuals(simulationOutput_Isolates)

# Pairwise comparison - Isolate + Substrate:

library(emmeans)

marginal_isolates = emmeans(GenLin_Isolates, ~ isolate)
pairs(marginal_isolates, adjust="bonferroni")

# Letters of significance:

library(multcomp)

lettere_isolates <- cld(marginal_isolates, alpha=0.05, Letters=letters, 
                        adjust="bonferroni")
lettere_isolates


# PLOTS

library(ggplot2)

  # Overall dataset

boxPlot_fullMycelium <- ggplot(data_mycelium, aes(x=tempCode, y=length,
                                        fill=isolate)) + 
                               geom_boxplot(width=0.5) + 
                               xlab("Temperature (°C)") + 
                               ylab("Mycelium extension (mm/day)") + 
                               ggtitle("Mycelium extension over temperature") +
                               theme(plot.title = element_text(hjust=0.5), 
                                     text = element_text(size=21)) +
                               labs(fill = "Isolate") +
                               scale_fill_manual(
                                 values = c("DA1" = "blue", "DA16" = "red",
                                            "DA5" = "orange"),
                                 labels = c("DA-1      b", "DA-16    a", "DA-5      a")) +
                               scale_x_discrete(guide = guide_axis(angle = 0),
                                                labels = c("5","10", "15",
                                                           "20", "25", "30",
                                                           "35", "40"))

boxPlot_fullMycelium


# Plots only isolates

boxPlot_Mycelium <- ggplot(data_mycelium, aes(x=isolate, y=length,
                                                  fill=isolate)) + 
                          geom_boxplot(width=0.5) + 
                          xlab("Temperature (°C)") + 
                          ylab("Mycelium extension (mm/day)") + 
                          ggtitle("Mycelium extension over the isolates") +
                          theme(plot.title = element_text(hjust=0.5), 
                              text = element_text(size=21), legend.position = "none") +
                              labs(fill = "Isolate") +
                          scale_x_discrete(guide = guide_axis(angle = 0),
                              labels = c("DA-1", "DA-16", "DA-5"))


boxPlot_Mycelium


# Plots of the best fit Brière equations

  # DA1 Isolate
    
    # Extract the subdataset

DA1_Dataset <- data_mycelium[data_mycelium$Isolate == 'DA1',][-1]

    # Best fit parameters from python curvefit

a_DA1 <- 1.08 * 10^(-4)
TL_DA1 <- 7.0501885
TM_DA1 <- 40.0
m_DA1 <- 1 / 0.6357534

    # Boxplot

boxPlot_Fit_DA1 <- ggplot(DA1_Dataset, aes(x=Temperature, y=Length_1_Day,
                                                  fill=as.factor(Temperature))) + 
                          geom_boxplot(width=0.5, aes(group = Temperature, fill = "Experimental data")) + 
                          xlab("Temperature (°C)") + 
                          ylab("Mycelium extension (mm/day)") + 
                          ggtitle("Isolate DA-1") +
                          theme(plot.title = element_text(hjust=0.5), 
                                text = element_text(size=21), legend.position = "left") +
                          stat_function(
                                fun = function(x) a_DA1 * x * (x - TL_DA1) * (TM_DA1 - x) ^ (m_DA1),
                                aes(color = "Best fit function"),
                                size = 0.7) +
                          scale_fill_manual(name = " ", values = c("Experimental data" = "skyblue")) +
                          scale_colour_manual(name = " ", values = c("Best fit function" = "red"))

boxPlot_Fit_DA1



  # DA5 Isolate

    # Extract the subdataset

DA5_Dataset <- data_mycelium[data_mycelium$Isolate == 'DA5',][-1]

    # Best fit parameters from python curvefit

a_DA5 <- 1.116 * 10^(-4)
TL_DA5 <- 6.4524548
TM_DA5 <- 40.0
m_DA5 <- 1 / 0.707223

    # Boxplot

boxPlot_Fit_DA5 <- ggplot(DA5_Dataset, aes(x=Temperature, y=Length_1_Day,
                                           fill=as.factor(Temperature))) + 
                          geom_boxplot(width=0.5, aes(group = Temperature, fill = "Experimental data")) + 
                          xlab("Temperature (°C)") + 
                          ylab("Mycelium extension (mm/day)") + 
                          ggtitle("Isolate DA-5") +
                          theme(plot.title = element_text(hjust=0.5), 
                          text = element_text(size=21), legend.position = 'left') +
                          stat_function(
                                        fun = function(x) a_DA5 * x * (x - TL_DA5) * (TM_DA5 - x) ^ (m_DA5),
                                        aes(color = 'Best fit function'),
                                        size = 0.7) +
                          scale_fill_manual(name = " ", values = c("Experimental data" = "skyblue")) +
                          scale_colour_manual(name = " ", values = c("Best fit function" = "orange"))

boxPlot_Fit_DA5


  # DA16 Isolate

    # Extract the subdataset

DA16_Dataset <- data_mycelium[data_mycelium$Isolate == 'DA16',][-1]

    # Best fit parameters from python curvefit

a_DA16 <- 9.31 * 10^(-5)
TL_DA16 <- 5.2668765
TM_DA16 <- 40.0
m_DA16 <- 1 / 0.7175657

    # Boxplot

boxPlot_Fit_DA16 <- ggplot(DA16_Dataset, aes(x=Temperature, y=Length_1_Day,
                                           fill=as.factor(Temperature))) + 
                            geom_boxplot(width=0.5, aes(group = Temperature, fill = "Experimental data")) + 
                            xlab("Temperature (°C)") + 
                            ylab("Mycelium extension (mm/day)") + 
                            ggtitle("Isolate DA-16") +
                            theme(plot.title = element_text(hjust=0.5), 
                            text = element_text(size=21), legend.position = 'left') +
                            stat_function(
                                          fun = function(x) a_DA16 * x * (x - TL_DA16) * (TM_DA16 - x) ^ (m_DA16),
                                          aes(color = 'Best fit function'),
                                          size = 1) +
                            scale_fill_manual(name = " ", values = c("Experimental data" = "skyblue")) +
                            scale_colour_manual(name = " ", values = c("Best fit function" = "forestgreen"))

boxPlot_Fit_DA16
                     
                     
                     
                     
                     


