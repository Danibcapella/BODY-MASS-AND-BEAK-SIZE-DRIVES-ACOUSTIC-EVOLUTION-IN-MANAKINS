#-------------------------------------------
## SETUP FOR PLOTS
#-------------------------------------------


library(tuneR)
library(seewave)
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)

## PLOT LABELLERS

# x label formatter
s_formatter <- function(x){
  lab <- paste0(x, " s")
}

# y label formatter
khz_formatter <- function(y){
  lab <- paste0(y, " kHz")
}

## THEMES

oscillo_theme_dark <- theme(panel.grid.major.y = element_line(color="black", linetype = "dotted"),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_rect(fill="transparent"),
                            panel.border = element_rect(linetype = "solid", fill = NA, color = "grey"),
                            axis.line = element_blank(),
                            legend.position = "none",
                            plot.background = element_rect(fill="black"),
                            plot.margin = unit(c(0,1,1,1), "lines"),
                            axis.title = element_blank(),
                            axis.text = element_text(size=14, color = "grey"),
                            axis.ticks = element_line(color="grey"))

hot_theme <- theme(panel.grid.major.y = element_line(color="black", linetype = "dotted"),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(fill="transparent"),
                   panel.border = element_rect(linetype = "solid", fill = NA, color = "grey"),
                   axis.line = element_blank(),
                   legend.position = "top",
                   legend.justification = "right",
                   legend.background = element_rect(fill="black"),
                   legend.key.width = unit(50, "native"),
                   legend.title = element_text(size=16, color="grey"),
                   legend.text = element_text(size=16, color="grey"),
                   plot.background = element_rect(fill="black"),
                   axis.title = element_blank(),
                   axis.text = element_text(size=16, color = "grey"),
                   axis.ticks = element_line(color="grey"))

hot_theme_grid <- theme(panel.grid.major.y = element_line(color="black", linetype = "dotted"),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(fill="transparent"),
                        panel.border = element_rect(linetype = "solid", fill = NA, color = "grey"),
                        axis.line = element_blank(),
                        legend.position = "top",
                        legend.justification = "right",
                        legend.background = element_rect(fill="black"),
                        legend.key.width = unit(50, "native"),
                        legend.title = element_text(size=16, color="grey"),
                        legend.text = element_text(size=16, color="grey"),
                        plot.background = element_rect(fill="black"),
                        plot.margin = margin(1,1,0,1, "lines"),
                        axis.title = element_blank(),
                        axis.text = element_text(size=16, color = "grey"),
                        axis.text.x = element_blank(),
                        axis.ticks = element_line(color="grey"))

## COLORS

hot_colors <- inferno(n=9)

#-------------------------------------------
## LOADING IN A WAV
#-------------------------------------------

# the path to .wav file



# loads a wave object from the .wav file path
wav <- readWave("N_lop_01.wav")

# builds a dataframe of the wave object values
sample <- seq(1:length(wav@left))
time <- sample/wav@samp.rate
sample.left <- as.vector(cbind(wav@left))
df <- data.frame(sample, time, sample.left)

# subsets the dataframe to a more manageable size for plotting
last.index <- tail(df$sample,1)
index <- seq(from = 1, to = last.index, by = 20)
df2 <- df[index,]


#-------------------------------------------
## GGSPECTRO PLOTS
#-------------------------------------------

# builds a spectrogram using ggspectro()
# note: no x-axis labels because the plot is designed to be aligned with the oscillogram in a grid
# for x-axis labels, use hot_theme instead of hot_theme_grid
hotplot <- ggspectro(wave = wav, f = wav@samp.rate, ovlp=90)+ 
  scale_x_continuous(labels=s_formatter, expand = c(0,0))+
  scale_y_continuous(breaks = seq(from = 5, to = 20, by=5), expand = c(0,0), labels = khz_formatter, position = "right")+
  geom_raster(aes(fill=amplitude), hjust = 0, vjust = 0, interpolate = F)+
  scale_fill_gradientn(colours = hot_colors, name = "Amplitude \n (dB)", na.value = "transparent", limits = c(-60,0))+
  hot_theme_grid

# builds an oscillogram
oscplot <- ggplot(df2)+
  geom_line(mapping = aes(x=time, y=sample.left), color="grey")+ 
  scale_x_continuous(labels=s_formatter, expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), position = "right")+
  geom_hline(yintercept = 0, color="white", linetype = "dotted")+
  oscillo_theme_dark

#-------------------------------------------
## PLOT GRID
#-------------------------------------------

gA=ggplot_gtable(ggplot_build(hotplot))
gB=ggplot_gtable(ggplot_build(oscplot))
maxWidth = grid::unit.pmax(gA$widths, gB$widths)
gA$widths <- as.list(maxWidth)
gB$widths <- as.list(maxWidth)
layo <- rbind(c(1,1,1),
              c(1,1,1),
              c(1,1,1),
              c(2,2,2))

grid.newpage()
grid.arrange(gA, gB, layout_matrix = layo)



### Espectrograma para parâmetros

song <- readWave("C_rub_06.wav")

spectro(song, osc = T, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = turbo)
spec(song, wl = 1024, flim = c(0, 10), plot = 2)


### Espectrogramas diversidade de chamados

# Antilophia galeata
#song_Antilophia <- readWave("A_gal_10.wav")

spectro(song_Antilophia, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2)


# Chiroxiphia caudata
#song_Chiro <- readWave("C_cau_17.wav")

spectro(song_Chiro, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2,collevels = seq(-40,0))


#Corapipo gutturalis
#song_Corapipo <- readWave("C_gut_05.wav")

spectro(song_Corapipo, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-50,0))


# Ceratopipra rubrocapilla
#song_Cerato <- readWave("C_rub_06.wav")

spectro(song_Cerato, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-40,0))


# Chloropipo unicolor
#song_Chloro <- readWave("C_uni_02.wav")

spectro(song_Chloro, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-40,0))


# Heterocercus linteatus
#song_Heteroc <- readWave("H_lin_03.wav")

spectro(song_Heteroc, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-50,0))


# Ilicura militaris
#song_Ilicura <- readWave("I_mil_07.wav")

spectro(song_Ilicura, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-40,0))


# Lepidothrix iris
#song_Lepido <- readWave("L_iri_01.wav")

spectro(song_Lepido, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-40,0))


# Masius chrysopterus
#song_Masius <- readWave("M_chr_02.wav")

spectro(song_Masius, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-40,0))


# Machaeropterus pyrocephalus
#song_Macha <- readWave("M_pyr_07.wav")

spectro(song_Macha, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-50,0))


# Neopelma chrysolophum
#song_Neo <- readWave("N_lop_01.wav")

spectro(song_Neo, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-50,0))


# Pipra fascicauda
#song_Pipra <- readWave("P_fas_08.wav")

spectro(song_Pipra, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-50,0))


# Pseudopipra pipra
#song_Pseudo <- readWave("P_pip_01.wav")

spectro(song_Pseudo, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2)


# Tyranneutes stolzmanni
#song_Tyra <- readWave("T_sto_07.wav")

spectro(song_Tyra, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2, collevels = seq(-50,0))


# Xenopipo uniformis
#song_Xeno <- readWave("X_uni_03.wav")

spectro(song_Xeno, scale = F, flim = c(0, 10), wn="blackman", wl = 1024, ovlp = 90, palette = reverse.gray.colors.2)

        