# Get remotes package
install.packages("remotes")
require(remotes)

# Use install_github to get RInno
install.packages("RInno")

# Require Package
require(RInno)
library(installr)
#At line 8, change [1-3] to [1-5]
trace(RInno::get_R,edit=T) 
#Repeat at line 9 
trace(RInno::code_section,edit=T)

# Use RInno to get Inno Setup
install_inno()
create_app(app_name = "Proteovis", app_dir = "App",
           pkgs =  c("ggplot2", "ggpubr", "gplots", "ggrepel", "tidyverse", "DT", "reshape2", "scales", "RColorBrewer", "colourpicker", "ggnewscale", "ggfortify", "ggrepel", "shiny", "shinyjs", "shinyWidgets", "eulerr", "webshot"))
compile_iss()

test <- read.delim("test1d.txt")
test <- test %>%
  mutate(Column = sub(".*Difference ", "", Column))
