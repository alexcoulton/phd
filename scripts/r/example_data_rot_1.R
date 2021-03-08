#NOTE --- CHANGE CHROMOSOME DISTANCE TO AVERAGE CHROMOSOME DISTANCE?

#seed for random number generation
set.seed(20)
?rnorm
library(xlsx) #import excel files package
library(ggplot2)
data = read.xlsx("~/Google Drive/University/PhD/Rotation 1/example data.xlsx", sheetIndex = 1)

#generate some random data based on normal dist. w/ mean 1.3 and s.d. 0.7
#data points 1-93 correspond to temperature of 10˚C
data$chromosome_distance[1:93] = rnorm(93, 1.3, 0.7)
#initiate dependent variable w/ e (epsilon) as the random term
e = rnorm(93, 1, 0.3)
data$num_recombinations[1:93] = 2 + 0.8 * data$chromosome_distance[1:93] + data$chromosome_distance[1:93]^2 + e #quadratic relationship


#generate data for additional temperatures
#data points 94-186 correspond to temperature of 15˚C
data$chromosome_distance[94:186] = rnorm(93, 1.3, 0.7)
e2 = rnorm(93, 1, 0.3)
data$num_recombinations[94:186] = 2 + 0.8 * data$chromosome_distance[94:186] + e2 #linear relationship

data$chromosome_distance[187:279] = rnorm(93, 1.3, 0.7)
e3 = rnorm(93, 1, 0.3)
data$num_recombinations[187:279] = 1.3 + 1.4 * data$chromosome_distance[187:279] + e3 #linear relationship

data$chromosome_distance[280:(280+92)] = rnorm(93, 1.3, 0.7)
e4 = rnorm(93, 1, 0.3)
data$num_recombinations[280:(280+92)] = 10 + (-1)*(1.1 * data$chromosome_distance[280:(280+92)] + 0.2*data$chromosome_distance[280:(280+92)]^2) + e4 #quad. relationship


#scatter plot with temperature as category
ggplot(data=data, aes(data$chromosome_distance, data$num_recombinations, colour=data$temperature)) + geom_point()


#adding another 93 rows, this time with temperature 20˚C
d=data.frame(rep(20, 93)) #initialize new data frame, repeat "20" 93 times
d$plant_number=NA #add new columns to new data frame
d$chromosome_distance=NA
d$num_recombinations=NA
names(d)=c("temperature", "plant_number","num_recombinations","chromosome_distance") #change names of columns in dataframe d to match those of "data" exactly
l=rbind(data,d) #combine both dataframes. 
data=l #set data equal to l. new rows have now been added.

#adding another 93 rows, this time with temperature 25˚C
d=data.frame(rep(25, 93)) #initialize new data frame, repeat "20" 93 times
d$plant_number=NA #add new columns to new data frame
d$chromosome_distance=NA
d$num_recombinations=NA
names(d)=c("temperature", "plant_number","num_recombinations","chromosome_distance") #change names of columns in dataframe d to match those of "data" exactly
l=rbind(data,d) #combine both dataframes. 
data=l #set data equal to l. new rows have now been added.

data$temperature[187:250] = c(20)
?rbind
#plot data
#?plot #info on plot package
plot(data$chromosome_distance, data$num_recombinations)


#trying to make any negative numbers positive... not quite working
for(i in data$num_recombinations) {
    if(data$num_recombinations[i] < 1) {
        data$num_recombinations[i] = data$num_recombinations[i]*(-1)
    }
}
