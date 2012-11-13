\name{arizona}
\alias{arizona}
\docType{data}
\title{Arizona vegetation dataset}
\description{
  This dataset gives the measurements of 16 vegetation communitites in the Santa Catalina Mountains, Arizona. The measurements were taken along different elevations from fir forest at high elevations, through pine forest, woodlands, and desert grassland.
}
\usage{arizona}
\format{
  A data frame with 16 observations and 8 variables. The variables refer to three latent concepts: 1) \code{ENV}=environment, 2) \code{SOIL}=soil, and 3) \code{DIV}=diversity.

  \tabular{llll}{
    \emph{Num} \tab \emph{Variable}	\tab \emph{Description}			\tab \emph{Concept}	\cr
    1   \tab env.elev 		\tab Elevation (m) 			\tab environment 	\cr
    2   \tab env.incli 		\tab Terrain inclination (degrees) 	\tab environment 	\cr
    3   \tab soil.ph 		\tab Acidity and base saturation 	\tab soil 		\cr
    4   \tab soil.orgmat 	\tab Organic matter content (perc) 	\tab soil 		\cr
    5   \tab soil.nitro 	\tab Nitrogen content (perc) 		\tab soil 		\cr    
    6   \tab div.trees 		\tab Number of species of trees 	\tab diversity 		\cr
    7   \tab div.shrubs 	\tab Numer of species of shrubs		\tab diversity 		\cr
    8   \tab div.herbs 		\tab Number of species of herbs		\tab diversity 		\cr    
  }

  The complete name of the rows are:
  1) Abies lasiocarpa, 2) Abies concolor, 3) Pseudotsuga menziesii-Abies Concolor, 4) Pseudotsuga menziesii, 
  5) Pinus ponderosa-Pinus strobiformis, 6) Pinus ponderosa, 7) Pinus ponderosa-Quercus, 8) Pinus chihuahuana,
  9) Pygmy conifer-oak scrub, 10) Open oak woodland, 11) Bouteloua curtipendula, 12) Spinose-suffrutescent,
  13) Cercidium microphyllum, 14) Larrea divaricata, 15) Cercocarpus breviflorus, 16) Populus tremuloides.
}
\source{
  Mixed data from Whittaker \emph{et al} (1968), and Whittaker and Niering (1975). See \bold{References} below.
}
\references{
  Whittaker, R. H., Buol, S. W., Niering, W. A., and Havens, Y. H. (1968) A Soil and Vegetation Pattern in the Santa Catalina Mountains, Arizona. \emph{Soil Science}, \bold{105}, pp. 440-450.

  Whittaker, R. H., and Niering, W. A. (1975) Vegetation of the Santa Catalina Mountains, Arizona. V. Biomass, Production, and Diversity Along the Elevation Gradient. \emph{Ecology}, \bold{56}, pp. 771-790.
}
\examples{
  data(arizona)
  arizona
}
\keyword{datasets}
