library(lavaan)
library(ggplot2)
library(plyr)

#compute cumulative response categories
cumulative_probs <- function(lavFit, Item){
    #simulate values for latent variable
    Xi <- rnorm(200, 0, 2)

    #dataframe with model estimates 
    DF1 <- parameterEstimates(lavFit)
    
    #extract discrimination parameter a and thresholds b for item from lavaan-object
    a <- DF1[DF1$rhs==Item & DF1$op=="=~",names(DF1)=="est"]
    b <- get_thresholds(lavFit, Item)

    #compute exponents for all categories
    all_z <- sapply(b, linexp, slope=a, Xi=Xi)

    #probability for answering in category c or higher is given by: pnorm(z)
    #compute probabilities for all categories
    P_c <- apply(all_z, 2, pnorm)

    DF2 <- data.frame(Xi, P_c)
    names(DF2) <- c("Xi", "Thresh1", "Thresh2", "Thresh3", "Thresh4")
    DF2
}

#plot cumulative response curves
cumulative_response <- function(lavFit, Item){

    Probs <- cumulative_probs(lavFit, Item)

    Plot <- ggplot(Probs, aes(Xi, Thresh1)) + geom_line() +
    geom_line(aes(Xi, Thresh2)) +
    geom_line(aes(Xi, Thresh3)) +
    geom_line(aes(Xi, Thresh4)) +
    scale_x_continuous(limits=c(-4, 4)) +
    xlab("Xi") + ylab("Score Category Probability") 

    Plot
}

#compute score category probabilities from cumulative response probabilities 
#NOTE: can only be done for items with 5 categories as of yet
category_probs <- function(LatentScore, Prob){
    P1 <- 1 - Prob[1]
    P2 <- Prob[1] - Prob[2]
    P3 <- Prob[2] - Prob[3]
    P4 <- Prob[3] - Prob[4]
    P5 <- Prob[4] - 0
    Ps <- c(P1, P2, P3, P4, P5)

    DF <- data.frame(LatentScore, t(Ps))
    names(DF) <- c("Xi", "P_Kat1", "P_Kat2", "P_Kat3", "P_Kat4", "P_Kat5")
    DF
}


#compute category response probabilities for an item and a value of the latent variable 
category_response <- function(lavFit, Item, LatentScore){
    #extract discrimination parameter a and thresholds b for item from lavaan-object
    a <- get_slope(lavFit, Item)
    b <- get_thresholds(lavFit, Item)

    #compute exponents for all categories
    all_z <- sapply(b, linexp, slope=a, Xi=LatentScore)

    #probability for answering in category c or higher is given by: pnorm(z)
    #compute probabilities for all categories
    P_c <- sapply(all_z, pnorm)
    DF <- category_probs(LatentScore, P_c)
    DF
}

#linear exponent for category c: z_c
linexp <- function(thresholds, slope, Xi) {
    res <- slope*(Xi - thresholds)
    res
}

#extract thresholds for an item from lavaan object
#takes fitted object and item name as string
get_thresholds <- function(lavFit, Item){
    Ests <- parameterEstimates(lavFit)
    thresh <- Ests[Ests$lhs==Item & Ests$op=="|", names(Ests)=="est"]

    thresh
}

#extract item slope parameter
#takes fitted object and item name as string
get_slope <- function(lavFit, Item){
    Ests <- parameterEstimates(lavFit)
    slope <- Ests[Ests$rhs==Item & Ests$op=="=~",names(Ests)=="est"]
    slope
}

#compute item information
item_information <- function(CumResp, lavFit, Item){
    #matrix that will hold all information values
    InfoMat <- matrix(nrow=nrow(CumResp), ncol=(ncol(CumResp)-1), byrow=TRUE)

    a <- get_slope(lavFit, Item)

    numRow <- nrow(CumResp)
    numCol <- ncol(CumResp) - 1

    #computing item information after formula as given by Reckase (2009)
    for (ro in 1:numRow){
        for (co in 1:numCol){
            InfoMat[ro,co] <- (a*CumResp[ro,co] * (1 - CumResp[ro,co]) - a*CumResp[ro,co+1] * (1 - CumResp[ro,co+1]))**2 / (CumResp[ro,co] - CumResp[ro,co+1])
        }
    }
    InfoMat
}


#plot item information 
information_plot <- function(lavFit, Item, from=-3, to=3, se=FALSE, values=FALSE){
    #compute slope and threshold for Item
    a <- get_slope(lavFit, Item)
    b <- get_thresholds(lavFit, Item)

    #simulate values for latent variable
    LatentVar <- runif(100, min=from, max=to)

    #compute all linear exponents for cumulative distribution function
    Expns <- aaply(.data=LatentVar, 1, .fun=linexp, thresholds=b, slope=a)

    #compute cumulative category responses 
    CumResp <- cbind(1, pnorm(Expns), 0)

    #compute item information
    InfoMat <- item_information(CumResp, lavFit, Item)

    #dataframe that holds information values
    InfoDF <- adply(.data=InfoMat, 1, .fun=sum)
    InfoDF <- cbind(InfoDF, LatentVar)

    #plot item information
    Plot <- ggplot(InfoDF, aes(LatentVar, V1)) + geom_line()

    #plot standard error of item information
    if (se == TRUE){
        seDF <- cbind(InfoDF, SE=1/sqrt(InfoDF$V1))

        Plot <- Plot + geom_line(data=seDF, aes(LatentVar, SE))
        Plot
    }

    #print actual item information as dataframe
    if(values == TRUE) {
        print(InfoDF)
    }
    Plot
}


#ICC-curves for an item
icc <- function(lavFit, Item){
    #simulate values from latent variable
    Xi <- rnorm(200, 0, 2)

    #dataframe with model estimates
    DF1 <- parameterEstimates(lavFit)

    #extract discrimination parameter a and thresholds b for item from lavaan-object
    a <- get_slope(lavFit, Item)
    b <- get_thresholds(lavFit, Item)

    #compute exponents for all categories
    all_z <- sapply(b, linexp, slope=a, Xi=Xi)

    #probability for answering in category c or higher is given by: pnorm(z)
    #compute probabilities for all categories
    P_c <- apply(all_z, 2, pnorm)

    #Kategorienwahrscheinlichkeiten aus Schwellenwahrscheinlichkeiten berechnen
    P.1to5 <- function(Prob){
        P1 <- 1 - Prob[,1]
        P2 <- Prob[,1] - Prob[,2]
        P3 <- Prob[,2] - Prob[,3]
        P4 <- Prob[,3] - Prob[,4]
        P5 <- Prob[,4] - 0
        Ps <- c(P1, P2, P3, P4, P5)
        Frame <- data.frame(Xi, P1, P2, P3, P4, P5)
        names(Frame) <- c("Xi", "P_Kat1", "P_Kat2", "P_Kat3", "P_Kat4", "P_Kat5")
        Frame
    }
    #CumProbs <- cumulative_probs(lavFit, Item)
    #DF2 <- sapply(Xi, category_response, lavFit, Item) 
    #Dataframe mit Kategorienwahrscheinlichkeiten fuer das Item ueber alle Werte von Xi
    DF2 <- P.1to5(P_c)

    Plot <- ggplot(DF2, aes(Xi, P_Kat1)) + geom_line(color="#654747") +
    geom_line(aes(Xi, P_Kat2), color="#11F315") +
    geom_line(aes(Xi, P_Kat3), color="#1149F3") +
    geom_line(aes(Xi, P_Kat4), color="#F311D5") +
    geom_line(aes(Xi, P_Kat5), color="#0A0109") +
    xlab("Xi") + ylab("Kategorienwahrscheinlichkeit") +
    coord_cartesian(xlim=c(-4, 4)) +
    theme(axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) + annotate("text", label=Itemname, x=0, y=0.95, size=4)
    ##if you want to print out the thresholds, uncomment the next line
    #print(b)
    Plot
}
