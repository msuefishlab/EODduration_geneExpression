### code to 1) import and plot the EOD duration data, 2) test for duration differences, and 3) plot the days tested for duration differences ###

library(car)
library(dplyr)
library(ggplot2)
sessionInfo <- sessionInfo()

## Part 1: plot EOD duration data

# import the file with the EOD duration data
duration <- read.csv("eod_duration.csv")


#plot
p1 <- ggplot(data = duration, aes(x = day, y = EOD, group=Treatment, color = Treatment)) +  stat_summary(fun = "mean", geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(width = 0.2, preserve = "total"), size=0.7, fatten = 0.5)
# + stat_summary(fun = "mean", geom = "line")
p2 = p1 + scale_x_continuous(name = "day",  breaks = scales::breaks_pretty(13)) + scale_y_continuous(name = "EOD duration (ms)", breaks = scales::breaks_pretty(10))
p3 =  p2 + scale_colour_brewer(palette="Dark2") + theme_classic() + theme(axis.text = element_text(size=9), legend.text = element_text(size = 8), legend.title = element_text(size = 9), legend.title.align = 1)
p4 = p3 + guides(color = guide_legend(override.aes = list(size = 0.1)))


#save
ggsave("EODbyDay.png", p4, "png")


## Part 2: test duration differences between the treatments on 1) day 0, 2) day of treatment end
#make Treatment a factor
duration$Treatment <- as.factor(duration$Treatment)

#extract data for day0
d0 <- duration %>% filter(day == 0)

#extract data for last day of treatment
dlast <- duration %>% filter((Treatment == "T1day" & day == 1) | (Treatment != "T1day" & day == 8))

#run ANOVAs, if significant, run Tukey HSD test
res.aov.d0 <- aov(EOD ~ Treatment, data = d0)
summary(res.aov.d0)

res.aov.dlast <- aov(EOD ~ Treatment, data = dlast)
summary(res.aov.dlast)
tuk.dlast <- TukeyHSD(res.aov.dlast)

#check that the ANOVA assumptions are met
#homogeneity of variances with Levene test
lev.d0 <- leveneTest(res.aov.d0)
lev.dlast <- leveneTest(res.aov.dlast)

#normality, with the Shapiro-Wilk test
resi.d0 <- residuals(object = res.aov.d0)
shp.d0 <- shapiro.test(x = resi.d0)

resi.dlast <- residuals(object = res.aov.dlast)
shp.dlast <- shapiro.test(x = resi.dlast)


## Part 3: Plot duration of EODs on day 0 and on last day of treatment 

#combine into one dataframe, so I can use facets
dummy.0 <- d0
dummy.0$facet <- 'day 0'

dummy.last <- dlast
dummy.last$facet <- 'last day of treatment'

d.plot <- rbind(dummy.0, dummy.last)

#plot
q1 <- ggplot(data = d.plot, aes(x = Treatment, y = EOD, group=Treatment, color = Treatment)) +  stat_summary(fun = "mean", geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), fatten = 1.8) + facet_wrap(~facet)
q2 = q1 + geom_point(shape = 2, position = position_jitter(width = 0.07, height = NULL, seed = 84))           
q3 = q2 + scale_y_continuous(name = "EOD duration (ms)", breaks = scales::breaks_pretty(10))
q4 = q3 + scale_colour_brewer(palette="Dark2") + theme_classic() + theme(axis.text = element_text(size=9), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA))


#save
ggsave("EODday0last.png", q4, "png")


#day 0
# q1 <- ggplot(data = d0, aes(x = Treatment, y = EOD, group=Treatment, color = Treatment)) +  stat_summary(fun = "mean", geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), fatten = 1.8)
# q2 = q1 + geom_point(shape = 2, position = position_jitter(width = 0.05, height = NULL, seed = 1))         
# q3 = q2 + scale_y_continuous(name = "EOD duration (ms)", breaks = scales::breaks_pretty(10))
# q4 = q3 + scale_colour_brewer(palette="Dark2") + theme_classic() + theme(axis.text = element_text(size=9), legend.text = element_text(size = 8), legend.title = element_text(size = 9), legend.title.align = 1)
# q5 = q4 + guides(color = guide_legend(override.aes = list(size = 0.1)))

