##### Survial Plotting Script with subsetting potential from Cris Print (verion = June 2019) #####

library(survival)
library(cowplot)
library(grid)
library(survminer)

print("Do you wish the plot to contain confidence intervals? - yes=(y)?")
Fl <- readLines(con = stdin(), n = 1, ok = TRUE)

print("HJ_BioInfo_Survival_RNAseqExpectedCount_Censored.txt")
FileChoice <- file.choose("TCGA PAAD Survival Censored 750 days Filtered RSEM NORM COUNT.txt") # choose input file
a <- read.delim(FileChoice, skip = 0, sep = "\t", as.is = TRUE, header = FALSE)
print("FIle is read")


colnames(a) <- a[1, ]
a <- a[-1, ]
#a[79] <- NULL
CCC <- a[, 6:ncol(a)]
time <- as.numeric(a[, 3])
status <- as.numeric(a[, 2])

#covariate <- a[, 2]
subset <- rep(TRUE, length(CCC[,1]))
#subset<-covariate=="G3"

outt <- matrix(NA, ncol = length(CCC[1, ]), nrow = 12)
outt[1, ] <- colnames(CCC)



for (counter in 1:length(CCC[1, ])) {
  print(counter)
  
  x <- CCC[, counter]
  y <- data.frame(cbind(time, status, x))
  z <- as.matrix(y$time)
  y$time <- as.numeric(z)
  z <- as.matrix(y$status)
  y$status <- as.numeric(z)
  z <- as.matrix(y$x)
  y$x <- as.numeric(z)
  x <- as.numeric(x)
  k <- coxph(Surv(time, status) ~ x, y, subset=subset)
  d <- summary(k)
  outt[11, counter] <- d$coef[5]
  outt[12, counter] <- sign(k$coefficients)
  
  for (c2 in 1:9) {
    
    nam <- paste("pl", c2, sep = ".")
    print(nam)
    x <- CCC[, counter]
    x <- as.numeric(x)
    i <- quantile(x, c2/10, na.rm = TRUE)
    x[x <= i] <- 0
    x[x > i] <- 1
    y <- data.frame(cbind(time, status, x))
    e <- survdiff(Surv(time, status) ~ x, data = y, subset=subset)
    outt[1 + c2, counter] <- 1 - pchisq(e$chisq, 1)
    dd <- survfit(Surv(time, status) ~ x, data = y, subset=subset)
    p = round(1 - pchisq(e$chisq, 1), 8)
    mm <- textGrob(paste("--- Gene =  ", outt[1, counter], "--- \n\n Cox proportional hazards model P = ", round(d$coef[5], 8), "\n\n Correlation of scaled Schoenfeld residuals \n with time (non-significant P supports the \n proportional hazards assumption)  \n P= ", 
                         round(cox.zph(k)$table[3], 8), "\n\n number individuals tested = ",k$n), gp = gpar(fontsize = 12, col = "blue"))
    ti <- paste("                    Quantile = ", 10 * c2, "%")
    s = paste("Time  \n\n(blue <= cut n=",table(x)[1],"; red > cut n=",table(x)[2],") \n Log rank p = ", round(p, 8), "; Direction =", sign(k$coefficients), "\n\n")
    
    if (Fl != "y") {
      assign(nam, ggsurvplot(dd, risk.table = FALSE, conf.int = FALSE, ggtheme = theme_minimal(), pval = F, legend = "none", 
                             title = ti, palette = c("blue", "red"), xlab = s, font.title = "bold")$plot)
    }
    
    if (Fl == "y") {
      assign(nam, ggsurvplot(dd, risk.table = FALSE, conf.int = TRUE, ggtheme = theme_minimal(), pval = F, legend = "none", title = ti, 
                             palette = c("blue", "red"), xlab = s, font.title = "bold")$plot)
    }
  }
  
  x <- CCC[, counter]
  y <- data.frame(cbind(time, status, x))
  z <- as.matrix(y$time)
  y$time <- as.numeric(z)
  z <- as.matrix(y$status)
  y$status <- as.numeric(z)
  z <- as.matrix(y$x)
  y$x <- as.numeric(z)
  x <- as.numeric(x)
  k <- coxph(Surv(time, status) ~ x, y, subset=subset)
  d <- summary(k)
  outt[11, counter] <- d$coef[5]
  outt[12, counter] <- sign(k$coefficients)
  
  ppp <- plot_grid(mm, pl.1, pl.2, pl.3, pl.4, pl.5, pl.6, pl.7, pl.8, pl.9, ncol = 5)
  save_plot(paste("Gene =  ", outt[1, counter], "_survplots.pdf"), ppp, base_aspect_ratio = 2, limitsize = FALSE, base_height = 10)
}

out <- rbind(outt[1, ], outt[2, ], outt[3, ], outt[4, ], outt[5, ], outt[6, ], outt[7, ], outt[8, ], outt[9, ], outt[10, ], outt[11, 
], outt[12, ])

outd <- rbind(c("gene", "10% quantile ChiSq p val", "20% quantile ChiSq p val", "30% quantile ChiSq p val", "40% quantile ChiSq p val", 
                "50% quantile ChiSq p val", "60% quantile ChiSq p val", "70% quantile ChiSq p val", "80% quantile ChiSq p val", "90% quantile ChiSq p val", 
                "Cox PH p val", "association direction"), t(data.frame(out)))
write.table(outd, file = "survival.txt", sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)



