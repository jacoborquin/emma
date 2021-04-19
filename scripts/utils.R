# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# various useful variables and functions,
# ggplot theme and color specifications 


# ----------------------------------------------------------------------
# Loading packages
# ----------------------------------------------------------------------

packages <- c(
    'tidyverse', 'extrafont', 'gridExtra', 
    'xtable', 'stargazer',
    'readxl',
    'data.table',
    'lme4',
    'lmerTest',
    'metafor',  
    'DescTools',  
    'cowplot',
    'irr',  
    'clubSandwich',
    'reshape2'
)
lapply(packages, library, character.only = TRUE)


# ----------------------------------------------------------------------
# Variables 
# ----------------------------------------------------------------------

# set the root directory for the project here:
# root = "/Users/au161118/Dropbox/ASB/Admin stuff/Posters & Papers/PAPERS/EMMA/scripts/emma"
# root = "/home/hstojic/research/project/attention_meta"
root = "."

# directory paths
figsDir <- file.path(root, "figs")
tablesDir <- file.path(root, "tables")
dataDir <- file.path(root, "data")

# font setup
fontSetup <- "Arial"
fontSize <- 1.75 
pointSize <- 1
lineSize <- 0.5
themeFontSize <- 6
mt <- 0.9
inch <- 0.3937008
cm <- 2.54

# single column: 8.8cm
height <- 6*inch*mt
width <- 8.8*inch*mt
fd_1c_1x0.33 <- c(height, 0.33*width)
fd_1c_1x0.5 <- c(height, 0.5*width)
fd_1c_1x0.66 <- c(height, 0.66*width)
fd_1c_1x1 <- c(height, width)
fd_1c_1.5x2 <- c(1.5*height, 2*width)
fd_1c_2x2 <- c(2*height, 2*width)
fd_1c_2.5x2 <- c(2.5*height, 2*width)
fd_1c_3x1 <- c(3*height, 1*width)
fd_1c_3.5x1 <- c(3.5*height, 1*width)
fd_1c_3x2 <- c(3*height, 2*width)


# supplementary information
height <- 7.5*inch*mt
width <- 8.9*inch*mt
fd_SI_1x0.33 <- c(height, 0.33*width)
fd_SI_1x0.5 <- c(height, 0.5*width)
fd_SI_1x0.66 <- c(height, 0.66*width)
fd_SI_1x1 <- c(height, width)
fd_SI_1x2 <- c(height, 2*width)
fd_SI_2x2 <- c(2*height, 2*width)
fd_SI_1x1.5 <- c(height, 1.5*width)
fd_SI_1.5x1.5 <- c(1.5*height, 1.5*width)
fd_SI_1.5x2 <- c(1.5*height, 2*width)
fd_SI_2x1.5 <- c(2*height, 1.5*width)
fd_SI_3x2 <- c(3*height, 2*width)


# ----------------------------------------------------------------------
# Ggplot themes
# ----------------------------------------------------------------------

tufte <- 
    theme(
        plot.title = element_text(
            size = themeFontSize+1, hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(1.1, "lines"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(lineend = 4, linetype = 1),
        axis.ticks.y = element_line(lineend = 4, linetype = 1),
        axis.ticks = element_line (colour = "black", size = 0.3), 
        axis.text = element_text(size = themeFontSize, colour = "black"),
        axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = themeFontSize),
        axis.title.y = element_text(vjust = 1.8),
        axis.title.x = element_text(vjust = -.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.justification = c(1,0),
        legend.position = c(1,0),
        legend.text = element_text(size = themeFontSize - 1),
        legend.key = element_rect(fill = "#FFFFFF"),
        legend.key.height = unit(0.8,"line"),
        strip.text = element_text(size = themeFontSize + 1),
        strip.background = element_rect(fill = "#FFFFFF"),
        text = element_text(family = fontSetup),
        validate = TRUE
    )

mytheme <- tufte


# ----------------------------------------------------------------------
# Color palettes
# ----------------------------------------------------------------------

# Color-blind friendly color combinations

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# ----------------------------------------------------------------------
# Other utilities
# ----------------------------------------------------------------------

# simple scatter plot with linear fit
genScatter <- function(plotData, x, y, xlab, ylab) {
	x <- enquo(x)
	y <- enquo(y)
	ggplot(data = plotData, aes(!!x, !!y)) +
	geom_point(alpha = .5, size = pointSize) +
	geom_smooth(method = "lm", color = "black", size = lineSize*1) +
	scale_x_continuous(breaks = seq(0,1,by=.1)) +
	xlab(xlab) +
	ylab(ylab) +
	mytheme
}

# forest plot function 
genForest = function(rmaobj, data, IV_level, label, varname) {
	fig = 
		data %>%
		select(Study, DV = varname, Domain = IV, vi.c) %>%
		mutate(SE = sqrt(vi.c)) %>%
		select(-(vi.c))
	Study = c("                               Summary effect")
	Domain = c("")
	DV = rmaobj[1]
	SE = rmaobj$se
	temp = data.table(Study, DV, SE, Domain)
	fig = rbind(fig,temp)
	fig$DV = as.numeric(fig$DV)
	fig$UL = FisherZInv(FisherZ(fig$DV) + 1.96*fig$SE)
	fig$LL = FisherZInv(FisherZ(fig$DV) - 1.96*fig$SE)
	fig$Domain = factor(fig$Domain, 
		levels = c(IV_level,""), 
		labels = c(label,"")
	)

	plot = 
		ggplot(fig, aes(x = reorder(Study,-DV), y = DV)) +
		geom_hline(yintercept = 0, size = lineSize, color = "grey") +
		geom_point(size = pointSize) +
		geom_errorbar(
            aes(ymin = LL, ymax = UL), 
            width = .2, size = lineSize*1
        ) +
		coord_flip() +
		facet_grid(Domain ~., scales = "free_y", space = "free")+
		mytheme +
		ylim(-1.05, 1.05) +
		xlab("") +
		ylab("") + 
		scale_color_manual(values = "black")+
		theme(legend.position = "none")

	return(plot)
}

# funnel plot function
genFunnel = function(rmaobj, data, varname, label) { 
	fig = 
		data %>%
		select(DV = varname, vi.c) %>%
		mutate(SE = sqrt(vi.c)) %>%
		select(-(vi.c))  
	estimate = FisherZ(as.numeric(rmaobj[2]))  # rma estimate 
	se.seq=seq(0, sqrt(max(data$vi.c)), 0.001) # make SE vector from 0 to max SE
	
	ll95 = estimate-(1.96*se.seq) 
	ul95 = estimate+(1.96*se.seq)
	ll99 = estimate-(3.29*se.seq)
	ul99 = estimate+(3.29*se.seq)
	
	dfCI = data.frame(ll95, ul95, ll99, ul99, se.seq, estimate)

 	plot <- 
		ggplot(aes(x = SE, y = DV), data = fig) +
		geom_segment(
			aes(x=sqrt(max(data$vi.c)), xend=0, y=estimate, yend=estimate),
			size = lineSize, color = "grey"
		) +
		geom_point(shape = 16, size = pointSize) +
		xlab('Standard Error') + 
		ylab('Effect size (Fisher transformed)') +
		geom_line(
            aes(x = se.seq, y = ll95), 
            data = dfCI, size = lineSize,
            linetype = 'dashed',
            color = "grey"
        ) +
		geom_line(
            aes(x = se.seq, y = ul95), 
            data = dfCI, size = lineSize,
            linetype = 'dashed',
            color = "grey"
        ) +
 	  geom_line(aes(x = se.seq, y = ll99), linetype = 'dashed', data = dfCI) +
 	  geom_line(aes(x = se.seq, y = ul99), linetype = 'dashed', data = dfCI) +
 	  scale_x_reverse() +
		coord_flip() +
        ggtitle(label) +
		mytheme 

	return(plot)        
}

# function for saving the plots 
savePlots <- function(figure, filename, dims, grob = FALSE, eps = FALSE) {
    myprint <- if (grob) grid.draw else print

    # PDFs 
    system2('rm', c(filename))
    # pdf(filename, 
    #     height = dims[1], width = dims[2], 
    #     pointsize = themeFontSize,
    #     colormodel = "cmyk",
    #     # colormodel = "gray",
    #     # family
    #     # fonts
    #     onefile = FALSE, 
    #     useDingbats = FALSE
    # );
    cairo_pdf(filename,
        height = dims[1], width = dims[2], 
        pointsize = themeFontSize,
        onefile = FALSE, 
        family = fontSetup, 
        fallback_resolution = 1500
    )
    myprint(figure)
    dev.off()
    # loadfonts(device = "pdf")
    # embed_fonts(filename)
    # system2('pdfcrop', c(filename, filename))
    
    # eps
    if (eps) {
        cairo_ps(filename, 
        # postscript(filename, 
            height = dims[1], width = dims[2], 
            pointsize = themeFontSize,
            fallback_resolution = 1200,
            family = "Arial", 
            # paper = "special", 
            # horizontal = FALSE,
            onefile = FALSE
        )
        myprint(figure)
        dev.off()
        loadfonts(device = "postscript")
        embed_fonts(filename, options = "-dEPSCrop")
    }
}

# ----------------------------------------------------------------------
# .tex functions
# ----------------------------------------------------------------------

onecoefTex = function(model, name){
  result <- paste0(
    "$\\beta_{0}=", 
    round(model[1], 2), 
    "$, $SE=", 
    round(model[2], 2), 
    "$, ", 
    ifelse(round(model[5], 2) == 0, "$p< 0.001$", paste0("$p=", round(model[5], 2),"$"))
    )
  cat(result, file = file.path(tablesDir, name))
}

oneofmanycoefTex = function(model, name, coefnum){
  result <- paste0(
    "$\\beta=", 
    round(model[coefnum, 1], 2), 
    "$, $SE=", 
    round(model[coefnum, 2], 2), 
    "$, ", 
    ifelse(round(model[coefnum, 5], 2) == 0, "$p< 0.001$", paste0("$p=", round(model[coefnum, 5], 2),"$"))
  )
  cat(result, file = file.path(tablesDir, name))
}

twocoefTex = function(model, name){
  result <- paste0(
    "$\\beta_{0}=", 
    round(summary(model)$coef[1,1], 2), 
    "$, $SE=", 
    round(summary(model)$coef[1,2], 2), 
    "$, ", 
    ifelse(round(summary(model)$coef[1,5], 2) == 0, "$p< 0.001$", paste0("$p=", round(summary(model)$coef[1,5], 2),"$")),
    ", ",
    "$\\beta_{1}=", 
    round(summary(model)$coef[2,1], 2), 
    "$, $SE=", 
    round(summary(model)$coef[2,2], 2), 
    "$, ", 
    ifelse(round(summary(model)$coef[2,5], 2) == 0, "$p< 0.001$", paste0("$p=", round(summary(model)$coef[2,5], 2),"$")) 
  )
  cat(result, file = file.path(tablesDir, name))
}

sandwichModTest = function(model, name){
  result <- paste0(
    "$t=", round(model$tstat[2], 3),"$, $p=", round(model$p[2], 3), "$"
  )
  cat(result, file = file.path(tablesDir, name))
}
