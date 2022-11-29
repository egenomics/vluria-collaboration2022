library(readr)
library(ggplot2)
internal_loop_lengths <- read_delim("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/internal_loop_lengths.txt","\t", escape_double = FALSE, trim_ws = TRUE)

ggplot(internal_loop_lengths, aes(loop_length, fill = group, colour = group)) +
  geom_density(alpha = 0.1) +
  xlim(0,120) +
  xlab('Loop length') +
  ylab('Density')

ggsave("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/fig3_internal_loop_lengths.png", width = 9, height = 9, units = "cm")


nrow(subset(internal_loop_lengths, group=='pfam'))
nrow(subset(internal_loop_lengths, group=='whole'))
nrow(subset(internal_loop_lengths, group=='pfam_whole'))

quantile(subset(internal_loop_lengths, group=='pfam')$loop_length, prob = c(0.75, 0.80, 0.85, 0.9, 0.95))

###### Interdomain loop distance
inter_domains_loop_lengths <- read_delim("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/inter_domains_loop_lengths.txt",    "\t", escape_double = FALSE, trim_ws = TRUE)
ggplot(inter_domains_loop_lengths, aes(loop_length)) +
  geom_density(alpha = 0.1)+
  xlim(0,100)

#Number of domains
lm_eqn = function(m) {

  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

summary_domains <- read_delim("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/summary_domains.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
summary_domains$age <- factor(summary_domains$age,levels = c("Hs_ps1-2", "Hs_ps13-17", "Hs_ps20-21", "Hs_ps31", "Hs_igen", "Hs_rand"))

library(dplyr)
# install.packages("dplyr")
mydata<-summary_domains %>% 
  group_by(age) %>% 
  dplyr::summarize(n=n(), mean=mean(num_domains), prot_length=mean(protein_length), have_domain=sum(num_domains>=1))

write.table(mydata, "/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/summary_domains_table.txt", sep="\t")

'''#Num_domains boxplot per age group
p <- ggplot(summary_domains, aes(age, num_domains))
p + geom_boxplot(aes(fill=age))  +
  ylim(0,10)+
  scale_color_brewer(palette="Set1")
'''

#Correlation plot num_domains - protein_length
p <- ggplot(data = summary_domains, aes(x = num_domains, y = protein_length)) +
  geom_jitter(aes(alpha=0.7))+
    geom_smooth(method = "lm", se=FALSE, color="orange", formula = y ~ x) +
  xlim(0,20)+
  geom_text(x = 10, y = 3500, label = lm_eqn(lm(protein_length ~ num_domains, summary_domains)), parse = TRUE) +
  xlab('Number of domains') + 
  ylab('Protein length') +
  theme(legend.position="none")
p

ggsave("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/correlation_domains_protlength.png", width = 9, height = 9, units = "cm")


########################
#Composition of domains#
########################
#stretches#

intra_domain_stretch_lengths <- read_delim("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/intra_domain_stretch_lengths.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ggplot(intra_domain_stretch_lengths, aes(stretch_length, fill = letter, colour = letter)) +
  geom_density(alpha = 0.1)+
  xlim(0,50)+
  scale_color_brewer(palette="Set1")+
  facet_wrap(~group)+
  xlab('Stretch length') + 
  ylab('Density')

ggsave("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/fig4_length_secondary_structures.png", width = 9, height = 9, units = "cm")


############
#Percentage#
############


#####################################
#Composition of secondary structures#
#####################################
library(plyr)
secondary_composition <- read_delim("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/all_proteins_composition.txt", "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)

secondary_composition$group <- factor(secondary_composition$group,levels = c("Hs_ps1-2", "Hs_ps13-17", "Hs_ps20-21", "Hs_ps31", "Hs_igen", "Hs_rand"))
data=ddply(secondary_composition, .(letter, group), function(x) sum(x$number_letter))
data$group <- factor(data$group,levels = c("Hs_ps1-2", "Hs_ps13-17", "Hs_ps20-21", "Hs_ps31", "Hs_igen", "Hs_rand"))

ggplot(data,aes(x=group,y=as.numeric(V1), fill=factor(letter))) + 
  geom_bar(position = "fill", stat='identity') +
  scale_fill_brewer(palette="Set1") +
  xlab("Age group") + 
  ylab("Percentage") +
  labs(fill="Structure") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/fig2A_all_protein_composition_secondary_structures_by_age.png", width = 16, height = 9, units = "cm")

# Load plyr
library(plyr)
intra_domain_stretch_composition <- read_delim("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/intra_domain_stretch_composition.txt",     "\t", escape_double = FALSE, trim_ws = TRUE)
head(intra_domain_stretch_composition)
# Split on the site variable, and keep all the other variables (is there an
# option to keep all variables in the final result?)
intra_domain_stretch_composition$group <- factor(intra_domain_stretch_composition$group,levels = c("Hs_ps1-2", "Hs_ps13-17", "Hs_ps20-21", "Hs_ps31", "Hs_igen"))
data=ddply(intra_domain_stretch_composition, .(letter, group), function(x) sum(x$number_letter))
data$group <- factor(data$group,levels = c("Hs_ps1-2", "Hs_ps13-17", "Hs_ps20-21", "Hs_ps31", "Hs_igen"))

ggplot(data,aes(x=group,y=as.numeric(V1), fill=factor(letter))) + 
  geom_bar(position = "fill", stat='identity') +
  scale_fill_brewer(palette="Set1")

ggsave("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/fig2B_domain_composition_secondary_structures_by_age.png", width = 16, height = 9, units = "cm")

'''
ggplot(secondary_composition, aes(number_letter, fill = letter, colour = letter)) +
  geom_density(alpha = 0.1)+
  xlim(0,50)+
  scale_color_brewer(palette="Set1")+
  facet_wrap(~group)
'''

'''
#Composition of top20 words for each age.
top20_word_count <- read_delim("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/top20_word_count.txt", " ", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
data=ddply(top20_word_count, .(age), function(x) c(L=sum(x$L), E=sum(x$E), H=sum(x$H)))
data2=ddply(top20_word_count, .(age), value=value)

data2<-group_by(data2, age) %>% mutate(percent = value/sum(value)*100,total=sum(value))

library(reshape2)
library(RColorBrewer)
colourCount = length(unique(mtcars$hp))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
data<-melt(data)

ggplot(data,aes(x=age,y=as.numeric(value), fill=factor(variable))) + 
  geom_bar(position = "fill", stat='identity') +
  scale_fill_brewer(palette="Set1")

ggplot(data2, aes(x=age, y=value, fill=factor(age))) + 
  geom_bar(position = "fill", stat='identity') + 
  scale_fill_manual(values = getPalette(colourCount))

library(ggplot2)
library(ggrepel)

data2$age <- factor(data2$age,levels = c("Hs_ps1-2", "Hs_ps13-17", "Hs_ps20-21", "Hs_ps31", "Hs_igen", "Hs_rand"))
  ggplot(data2, aes(x = 1, y = 1, size = percent, label = word, color=age)) +
  geom_text_repel(segment.size = 0, force = 100) +
  scale_size(range = c(2, 15), guide = FALSE) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = '') +
  theme_classic()+
    facet_wrap(~age)

'''
library(ggwordcloud)

word_count_short <- read_delim("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/word_count_length_1_to_3.txt", "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
word_count_short$age <- factor(word_count_short$age,levels = c("Hs_ps1-2", "Hs_ps13-17", "Hs_ps20-21", "Hs_ps31", "Hs_igen", "Hs_rand"))

word_count_short$total<-tapply(word_count_short$count, word_count_short$age, FUN=sum)[word_count_short$age]
tapply(word_count_short$count, word_count_short$age, FUN=sum)

word_count_short$percent= (word_count_short$count/word_count_short$total)*100

write.table(word_count_short, file = "/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/short_words_table.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE,fileEncoding = "")

#Word cloud for short words (1-3)
ggplot(word_count_short, aes(x = 1, y = 1, size = percent, label = word, color=age)) +
  geom_text_wordcloud() +
  scale_size(range = c(2, 15), guide = "none") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = '') +
  theme_classic() +
  facet_wrap(~age)

ggsave("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/fig5_wordcloud_short_words.png", width = 9, height = 6, units = "in")


ggplot(word_count_short,aes(x = age, y=percent, color=age)) + 
  geom_point() +
  facet_wrap(~word)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_color_brewer(palette="Set1") + 
  xlab('Age group')+
  ylab('Percentage')

ggsave("/home/jlvillanueva/Documents/analysis/vluria-collaboration2022/results/fig5_word_combinations_by_age.png", width = 9, height = 6, units = "in")


