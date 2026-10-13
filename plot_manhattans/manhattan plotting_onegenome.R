library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggtext)
#library(normentR)

setwd("~/analysis/banana_kmerGWAS/manhattan_MartinDataset_AAv4")

input<-read.delim("AAv4ref.join.txt.windows_pval_freq.txt",header = F, sep = "\t")
names(input) <- c("chr","startbp","stopbp","pvalue","counting")


data_cum <- input |>
  group_by(chr) |>
  summarise(max_bp = max(stopbp)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(chr, bp_add)

input <- input |>
  inner_join(data_cum, by = "chr") |>
  mutate(bp_cum = stopbp + bp_add)



axis_set <- input |>
  group_by(chr) |>
  summarize(center = mean(bp_cum))

filter.input <- input  |> filter(pvalue>2.5) |> filter(counting>=10) 

ggplot(filter.input, aes(  x = bp_cum, y = pvalue,  color = as_factor(chr), size = counting)) + 
  geom_point(aes(size=counting), alpha=0.9) + scale_size(range=c(1,9), breaks=c(10,20,30,40,50,60,70,80,90)) +
  scale_x_continuous(
    label = axis_set$chr,
    breaks = axis_set$center
  )  +
  scale_color_manual(values = rep(
    c("#276FBF", "#183059"),
    unique(length(axis_set$chr))
  )) + theme_minimal()


#facets
filter.input <- input  |> filter(counting>9) 

ggplot(filter.input, aes(  x = ((startbp+stopbp)/2000000), y = pvalue,  color = pvalue, size = counting)) + 
  geom_point(alpha=0.8) + scale_size(range=c(1,10), breaks=c(5,10,15,20,25,30,25,40,45,50,500)) +
  scale_x_continuous()    + ylim(2.5,8.5) + scale_colour_gradient(low="#fdd49e",high="#990000") +
   theme_minimal() +facet_wrap(~chr,nrow=4,ncol=3, scales = "fixed")


#chr7
filter.input <- input   |> filter(chr=="chr02") |> filter(counting>9) 


ggplot(filter.input, aes(  x = ((startbp+stopbp)/2000000), y = pvalue,  color = pvalue, size = counting)) + 
  geom_point(alpha=0.8) + scale_size(range=c(1,10), breaks=c(5,10,15,20,25,30,25,40,45,50,500)) +
  scale_x_continuous()    + ylim(2.5,8.5) + scale_colour_gradient(low="#fdd49e",high="#990000") +
  theme_minimal() 
