require(this.path)
require(ggplot2)
setwd(this.path::here())
load("C1.Rdata")

truth <- read.csv("../Data/TruthC1.csv")

# Basic line plot with points

p=ggplot(data=data.frame(truth=as.numeric(truth$x),
                         pred=c(AnswerC1[,2])),
         aes(x=truth, y=c(AnswerC1[,2]))) +
  geom_point()+theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20))+
 xlab( "Truth")+ylab( "Predicted quantile")+coord_fixed()+
  geom_errorbar(aes(ymin=c(AnswerC1[,3]), ymax=c(AnswerC1[,4])), width=.2)+
  ylim(range(as.numeric(truth$x),AnswerC1[,2:4]))+
  xlim(range(as.numeric(truth$x),AnswerC1[,2:4]))+
  geom_abline()
plot(p)

ggsave(p, filename =paste0("../Figures/C1.pdf"), height = 7, width = 7, bg = "transparent" )

mean(truth>AnswerC1[,3] & truth < AnswerC1[,4])
