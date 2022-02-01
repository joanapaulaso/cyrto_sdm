# Call libraries ----
library(raster)
library(ENMTML)
library(biomod2)
library(rgdal)
library(sdmpredictors)
library(maptools)
library(usdm)
library(ecospat)
library(CoordinateCleaner)
library(rgbif)
library(spocc)
library(spThin)
library(dplyr)
library(dismo)
library(gridExtra)
library(tidyr)
library(corrplot)
library(geometry)
library(sp)
library(ggfortify)
library(ggplot2)

# Set wd ----
setwd("D:/CyrtoMDE")
getwd()

# Load image ----
load("preproc_thin.RData")
# Generate pseudo-absences for bioclim2 ----
resp_occ = as.numeric("resp_occ")
str(resp_occ)
cyrtopodium_matrix_thin[,"resp_occ"] = 1
resp_occ = cyrtopodium_matrix_thin$resp_occ
is.numeric(resp_occ)
colnames(cyrtopodium_matrix_thin)
# View(cyrtopodium_matrix_thin)
resp_name <- 'cyrto'
is.numeric(cyrtopodium_matrix_thin$resp_occ)
myResp <- as.numeric(cyrtopodium_matrix_thin[,"resp_occ"])
myRespXY <- cyrtopodium_matrix_thin[which(myResp == 1), c("decimalLongitude", "decimalLatitude")]
myResp <- reclassify(subset(biostack_final, 1, drop = TRUE), c(-Inf, Inf, 0))
myResp[cellFromXY(myResp, myRespXY)] <- 1

bm_cyrto100disk = BIOMOD_FormatingData(resp.var = resp_occ,
                                          expl.var = biostack_final,
                                          resp.xy = myRespXY,
                                          resp.name = resp_name,
                                          PA.nb.rep = 5,
                                          PA.nb.absences = 510, 
                                          PA.strategy = "disk", 
                                          PA.dist.min = 50000)

bm_cyrto100disk
plot(bm_cyrto100disk)

# Model options ----

# modelo_op <- BIOMOD_ModelingOptions(
#   MAXENT.Phillips = list(
#     path_to_maxent.jar = "./maxent/maxent.jar")
# )

# Modeling ----

cyrto1model = BIOMOD_Modeling(
  data = bm_cyrto100disk,
  models = c("RF", "GLM", "SRE", "MAXENT.Phillips.2"),
  # models.options = modelo_op,
  NbRunEval = 5,
  DataSplit = 70,
  Prevalence = 0.5,
  VarImport = 5,
  models.eval.meth = c("TSS", "ROC"),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = "cyrto"
)

cyrto1model

# Summarizing varibles importance ----

var_import_cyrto1 = get_variables_importance(cyrto1model) 
var_import_cyrto1
var_import_cyrto2 = apply(var_import_cyrto1, c(0, 1, 2, 3, 4), mean)
var_import_cyrto2

dir.create("Outputs_Variables_Importance")
write.csv(var_import_cyrto1, paste0("./Outputs_Variables_Importance/", "selected_var_import_cyrto1_24-12.csv"))
write.csv(var_import_cyrto2, paste0("./Outputs_Variables_Importance/", "selected_var_import_cyrto2_24-12.csv"))

cyrto1_rf = BIOMOD_LoadModels(cyrto1model, models = 'RF')
cyrto1_glm = BIOMOD_LoadModels(cyrto1model, models = 'GLM')
cyrto1_sre = BIOMOD_LoadModels(cyrto1model, models = 'SRE')
cyrto1_maxEnt = BIOMOD_LoadModels(cyrto1model, models = 'MAXENT.Phillips.2')


resp_curve_rf <- biomod2::response.plot2(
  models = cyrto1_rf, 
  Data = get_formal_data(cyrto1model, 'expl.var'), 
  show.variables = cyrto1model@expl.var.names[c(1, 2, 3, 4, 5, 6)],
  do.bivariate = F, 
  name= "RF_curva_resposta", 
  fixed.var.metric = 'mean', 
  legend = F, 
  display_title = T, 
  data_species = get_formal_data(cyrto1model, 'resp.var') 
)

resp_curve_glm <- biomod2::response.plot2(
  models = cyrto1_glm, 
  Data = get_formal_data(cyrto1model, 'expl.var'), 
  show.variables = cyrto1model@expl.var.names[c(1, 2, 3, 4, 5, 6)], 
  do.bivariate = F, 
  name= "GLM_curva_resposta", 
  fixed.var.metric = 'mean', 
  legend = F, 
  display_title = T, 
  data_species = get_formal_data(cyrto1model, 'resp.var')
)

resp_curve_sre <- biomod2::response.plot2(
  models = cyrto1_sre, 
  Data = get_formal_data(cyrto1model, 'expl.var'), 
  show.variables = cyrto1model@expl.var.names[c(1, 2, 3, 4, 5, 6)], 
  do.bivariate = F, 
  name= "SRE_curva_resposta", 
  fixed.var.metric = 'mean', 
  legend = F, 
  display_title = T, 
  data_species = get_formal_data(cyrto1model, 'resp.var') 
)

resp_curve_maxEnt <- biomod2::response.plot2(
  models = cyrto1_maxEnt, 
  Data = get_formal_data(cyrto1model, 'expl.var'), 
  show.variables = cyrto1model@expl.var.names[c(1, 2, 3, 4, 5, 6)], 
  do.bivariate = F, 
  name= "MaxEnt_curva_resposta", 
  fixed.var.metric = 'mean', 
  legend = F, 
  display_title = T, 
  data_species = get_formal_data(cyrto1model, 'resp.var') 
)
