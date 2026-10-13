library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggtext)
library(stringr)
#library(normentR)


input<-read.delim("supergenome.join.txt.windows_pval_freq.txt.gz",header = F, sep = "\t")
names(input) <- c("chr","startbp","stopbp","pvalue","counting")

#extract the chr number from the chr name and add to a new column
input$chrnum <- as.numeric( str_split_i(input$chr,"chr",2) )

# data_cum <- input |>
#   group_by(chr) |>
#   summarise(max_bp = max(stopbp)) |>
#   mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
#   select(chr, bp_add)
# 
# input <- input |>
#   inner_join(data_cum, by = "chr") |>
#   mutate(bp_cum = stopbp + bp_add)
# 
# 
# 
# axis_set <- input |>
#   group_by(chr) |>
#   summarize(center = mean(bp_cum))
# 
# filter.input <- input  |> filter(pvalue>2.5) |> filter(counting>=10) 
# 
# ggplot(filter.input, aes(  x = bp_cum, y = pvalue,  color = as_factor(chr), size = counting)) + 
#   geom_point(aes(size=counting), alpha=0.9) + scale_size(range=c(1,9), breaks=c(10,20,30,40,50,60,70,80,90)) +
#   scale_x_continuous(
#     label = axis_set$chr,
#     breaks = axis_set$center
#   )  +
#   scale_color_manual(values = rep(
#     c("#276FBF", "#183059"),
#     unique(length(axis_set$chr))
#   )) + theme_minimal()


#facets

# AAAcav
# AAAgros
# AAmalSCBG
# AAmalv4
# AAzeb
# BB
# SS

filter.input <- input  |> filter(grepl("AAAgros",chr))  |> filter(counting>14) |> arrange(as.numeric(chrnum))

# first.half <- c(1,2,3,4,5,6,
#                12,13,14,15,16,17)

first.third <- c(1,2,3,12,13,14,23,24,25)


ggplot(filter.input, aes(  x = ((startbp+stopbp)/2000000), y = pvalue,  color = pvalue, size = counting)) + 
  geom_point(alpha=0.8) + scale_size(range=c(1,10), breaks=c(5,10,15,20,25,30,25,40,45,50,500)) +
  scale_x_continuous()    + ylim(2.5,10) + scale_colour_gradient(low="#fdd49e",high="#990000") +
   theme_minimal() +facet_wrap(~chrnum,nrow=6,ncol=6, scales = "fixed")






###one chr only
filter.input <- input   |> filter(chr=="chr02") |> filter(counting>9) 


ggplot(filter.input, aes(  x = ((startbp+stopbp)/2000000), y = pvalue,  color = pvalue, size = counting)) + 
  geom_point(alpha=0.8) + scale_size(range=c(1,10), breaks=c(5,10,15,20,25,30,25,40,45,50,500)) +
  scale_x_continuous()    + ylim(2.5,8.5) + scale_colour_gradient(low="#fdd49e",high="#990000") +
  theme_minimal() 
