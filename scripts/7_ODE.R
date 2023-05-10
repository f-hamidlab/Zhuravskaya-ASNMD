## Load dependencies for this workflow
library(pracma)
library(tidyverse)
library(patchwork)

# set variables
vsr  = 1.76
t0.5 = 7.08 
kdr  = log(2)/t0.5
rt <- 0.3929945
rn0 <- 0.7817625
rnf <- 0.1187838
R0wo = vsr/kdr      # Initial steady-state mRNA abundance in early development (without rn0)
R0wi = vsr*rn0/kdr  # Initial steady-state mRNA abundance in early development (with rn0)
t0   = 0            # Start, hours
tf   = 144          # Finish, hours

Rwo = function(t, y) {
  if(tt == 0){
    vsr - kdr*y
  } else {
    ifelse(t<tt, vsr - kdr*y, vsr*rt - kdr*y)
  }
}
Rwi = function(t, y) {
  if(tt == 0){
    ifelse(t<tn, vsr*rn0 - kdr*y, vsr*rnf - kdr*y)
  } else if(tn < tt){
    ifelse(t<tn, vsr*rn0 - kdr*y, ifelse(t<tt, vsr*rnf- kdr*y, vsr*rnf*rt - kdr*y))
  } else if(tn == tt){
    ifelse(t<tt, vsr*rn0 - kdr*y, vsr*rt*rnf - kdr*y)
  } else if(tn > tt){
    ifelse(t<tt, vsr*rn0 - kdr*y, ifelse(t<tn, vsr*rt*rn0 - kdr*y, vsr*rt*rnf - kdr*y))
  }
}

dfval <- data.frame(tn = c(72, 48, 72, 96),
           tt = c(0, 72, 72, 72),
           title = c("NMD ONLY", "NMD BEFORE", "NMD SIMULTANEOUS", "NMD AFTER"))

plot_non_norm <- do.call(wrap_plots,apply(dfval, 1, function(row){
  tn <<- as.numeric(row[1])
  tt <<- as.numeric(row[2])
  
  Swo <- ode45(Rwo, t0, tf, R0wo)
  Swi <- ode45(Rwi, t0, tf, R0wi)
  Iwo <-  deval(Swo$t, Swo$y, t0:tf, idx = NULL)/R0wo
  Iwi <-  deval(Swi$t, Swi$y, t0:tf, idx = NULL)/R0wo
  
  res.df <- data.frame(time=t0:tf/24,
             with=Iwi,
             without=Iwo) %>%
    pivot_longer(with:without, names_to = "group", values_to = "I") 
  
  arrow.df <- res.df %>%
    filter((time==tn/24 & group=="with")|(time==tt/24 & group=="without")) %>%
    filter(time != 0)
  
  res.df %>% 
    ggplot(aes(x=time, y=I, group=group, colour=group)) + 
    geom_line(aes(linetype= group)) +
    scale_linetype_manual(values = c("dashed","solid")) +
    scale_colour_manual(values = c("red","black")) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0,1.2), expand = c(0,0), breaks = 0:2/2) +
    scale_x_continuous(breaks = c(0,6), minor_breaks = 1:5) +
    geom_segment(data = arrow.df, aes(x=time,xend=time, y=I+0.001,yend=I+0.002),
                 linewidth=2, lineend = "round", linejoin = "bevel",
                 arrow = arrow(length = unit(-0.1, "inches")))

})) + plot_layout(nrow = 1)


plot_norm <- do.call(wrap_plots,apply(dfval, 1, function(row){
  tn <<- as.numeric(row[1])
  tt <<- as.numeric(row[2])
  
  Swo <- ode45(Rwo, t0, tf, R0wo)
  Swi <- ode45(Rwi, t0, tf, R0wi)
  Iwo <-  deval(Swo$t, Swo$y, t0:tf, idx = NULL)/R0wo
  Iwi <-  deval(Swi$t, Swi$y, t0:tf, idx = NULL)/R0wo
  
  res.df <- data.frame(time=t0:tf/24,
                      dat=Iwo/Iwi) 
  
  arrow.df <- data.frame(time=c(tn,tt)/24,
                         group=c("with","without")) %>%
    
    left_join(res.df) %>%
    filter(time != 0) 

  
  res.df %>%
    ggplot(aes(x=time, y=dat)) + 
    geom_line(colour="black") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0,10), expand = c(0,0), breaks = 0:5*2) +
    scale_x_continuous(breaks = c(0,6), minor_breaks = 1:5) +
    geom_segment(data = arrow.df, aes(x=time,xend=time, y=dat+0.001,yend=dat+0.002,
                                      colour = group),
                 linewidth=2, lineend = "round", linejoin = "bevel",
                 arrow = arrow(length = unit(-0.1, "inches")))

})) + plot_layout(nrow = 1)

