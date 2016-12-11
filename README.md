# Microsoft-R-Client

Data Science full setup from Microsoft R Client! (1000+ packages)

You need:

* Rtools (appropriate version, most recently if you don't know)
* Microsoft R Client
* R / Microsoft R Open (for Laurae2/sparsity package, or any other package refusing to compile in Microsoft R Client)
* RStudio
* MinGW
* Git Bash
* A PATH setup appropriately (MinGW before Rtools)
* My "hacked" `etc/profile.site`
* Windows or Linux (instructions for xgboost are for Windows but adaptable very easily for Linux)

The following adds most of the packages you will ever need.

```r
packages <- c("abind", "acepack", "actuar", "ActuDistns", "ada", 
"adabag", "ade4", "ade4TkGUI", "adegraphics", "adehabitatLT", 
"adehabitatMA", "ADGofTest", "AER", "AGD", "agricolae", "AICcmodavg", 
"akima", "alabama", "AlgDesign", "alphahull", "alr3", "alr4", 
"amap", "Amelia", "anchors", "animation", "aod", "aods3", "ape", 
"aplpack", "argparse", "arm", "arules", "arulesViz", "ascii", 
"assertthat", "AUC", "BaBooN", "backports", "barcode", "bartMachine", 
"bartMachineJARs", "base", "base64", "base64enc", "BatchJobs", 
"BayesFactor", "BayesX", "BayesXsrc", "BB", "BBmisc", "bbmle", 
"BCA", "bcp", "bdsmatrix", "betareg", "BH", "BHH2", "BiasedUrn", 
"bibtex", "biclust", "BIFIEsurvey", "biganalytics", "biglm", 
"bigmemory", "bigmemory.sri", "bigRR", "bigtabulate", "binda", 
"binGroup", "bisoreg", "bit", "bit64", "bitops", "blme", "BMA", 
"boot", "bootstrap", "Boruta", "BradleyTerry2", "breakpoint", 
"brew", "brglm", "brnn", "broom", "BsMD", "bst", "btergm", "C50", 
"ca", "Cairo", "cairoDevice", "CALIBERrfimpute", "calibrator", 
"candisc", "caper", "car", "CARBayes", "CARBayesdata", "care", 
"caret", "caretEnsemble", "catdata", "caTools", "cba", "ccaPP", 
"cclust", "CDM", "CDVine", "cem", "censReg", "CEoptim", "changepoint", 
"checkmate", "checkpoint", "chemometrics", "chron", "circlize", 
"CircStats", "Ckmeans.1d.dp", "class", "classInt", "clue", "cluster", 
"clusterSim", "clustvarsel", "clv", "clValid", "cmaes", "cmprsk", 
"coda", "codetools", "coin", "colorplaner", "colorspace", "colourpicker", 
"combinat", "compiler", "CompQuadForm", "compute.es", "conf.design", 
"contfrac", "contrast", "copula", "CORElearn", "corpcor", "corrgram", 
"corrplot", "covr", "CoxBoost", "coxme", "cplm", "crayon", "crossval", 
"crp.CSFP", "crrstep", "crs", "cshapes", "cubature", "Cubist", 
"curl", "cvAUC", "cvTools", "d3heatmap", "d3Network", "DAAG", 
"data.table", "DatABEL", "dataframes2xls", "datasets", "date", 
"dbConnect", "DBI", "dbscan", "Deducer", "DeducerExtras", "deepnet", 
"degreenet", "deldir", "dendextend", "dendextendRcpp", "dendroextras", 
"DendSer", "denstrip", "DEoptim", "DEoptimR", "deployrRserve", 
"depthTools", "Deriv", "descr", "DescTools", "deSolve", "Devore7", 
"devtools", "dfoptim", "diagram", "DiagrammeR", "DiagrammeRsvg", 
"DiceDesign", "DiceKriging", "DiceOptim", "dichromat", "digest", 
"diptest", "directlabels", "discretization", "DiscriMiner", "distr", 
"distrEx", "DistributionUtils", "diveMove", "dlm", "DMwR", "doBy", 
"DoE.base", "DoE.wrapper", "doMPI", "doParallel", "doRedis", 
"downloader", "dplR", "dplyr", "drat", "DT", "dtplyr", "dtw", 
"dygraphs", "dynamicTreeCut", "dynlm", "e1071", "eaf", "earth", 
"Ecdat", "Ecfun", "ecodist", "effects", "eha", "elasticnet", 
"ElemStatLearn", "ellipse", "elliptic", "elmNN", "emdbook", "emoa", 
"emulator", "energy", "ENmisc", "EnQuireR", "entropy", "EntropyExplorer", 
"Epi", "EpiModel", "epitools", "erer", "ergm", "ergm.count", 
"ergm.userterms", "eRm", "estimability", "etm", "evaluate", "evd", 
"ExplainPrediction", "expm", "extrafont", "extrafontdb", "extraTrees", 
"factoextra", "FactoMineR", "Fahrmeir", "fail", "faraway", "fArma", 
"fAsianOptions", "fAssets", "fastcluster", "fastdigest", "fastICA", 
"fastmatch", "fastR", "fBasics", "fCopulae", "fda", "fdrtool", 
"FeaLect", "feather", "FeatureHashing", "fExoticOptions", "fExtremes", 
"ff", "ffbase", "FFTrees", "fftw", "fGarch", "fields", "filehash", 
"fImport", "findpython", "fit.models", "fitdistrplus", "flare", 
"flashClust", "flexclust", "flexmix", "flexsurv", "FME", "FMStable", 
"fMultivar", "FNN", "fNonlinear", "fontcm", "fOptions", "foreach", 
"forecast", "foreign", "formatR", "formattable", "Formula", "fortunes", 
"forward", "fpc", "fPortfolio", "fracdiff", "FRB", "frbs", "fRegression", 
"FrF2", "FrF2.catlg128", "FSelector", "fTrading", "fts", "fUnitRoots", 
"futile.logger", "futile.options", "GA", "gam", "gamair", "GAMBoost", 
"gamboostLSS", "gamlss", "gamlss.data", "gamlss.dist", "gamm4", 
"gapminder", "gbm", "gclus", "gdata", "gdtools", "gee", "geeM", 
"geepack", "GenABEL", "GenABEL.data", "GeneralizedHyperbolic", 
"genetics", "GenSA", "geoR", "geoRglm", "geosphere", "GERGM", 
"getopt", "GGally", "ggcorrplot", "ggdendro", "ggfortify", "ggiraph", 
"ggm", "ggplot2", "ggplot2movies", "ggrepel", "ggthemes", "ggvis", 
"git2r", "glasso", "glmmML", "glmnet", "glmulti", "GlobalOptions", 
"gmailr", "gmm", "gmodels", "gmp", "gnm", "gof", "goftest", "googleVis", 
"gpairs", "GPArotation", "GPfit", "gplots", "graphics", "gRbase", 
"grDevices", "grid", "gridBase", "gridExtra", "grouped", "gsl", 
"gss", "gstat", "gsubfn", "gtable", "gtools", "Guerry", "gWidgets", 
"gWidgetsRGtk2", "gWidgetstcltk", "h2o", "haplo.stats", "haven", 
"hdi", "heatmaply", "heplots", "hergm", "hexbin", "hglm", "hglm.data", 
"HH", "HiClimR", "highlight", "highr", "hmeasure", "Hmisc", "hms", 
"HSAUR", "HSAUR2", "HSAUR3", "htmlTable", "htmltools", "htmlwidgets", 
"httpuv", "httr", "huge", "hunspell", "hwriter", "hypergeo", 
"ibdreg", "ic.infer", "ICS", "ICSNP", "igraph", "igraphdata", 
"ineq", "influenceR", "Information", "infotheo", "inline", "inlinedocs", 
"intergraph", "intervals", "intsvy", "iplots", "ipred", "irace", 
"irlba", "irr", "isa2", "Iso", "ISOcodes", "isotone", "ISwR", 
"iterators", "itertools", "its", "JavaGD", "JGR", "jomo", "jpeg", 
"jsonlite", "kappalab", "kdecopula", "Kendall", "kernlab", "KernSmooth", 
"KFAS", "kinship2", "kknn", "klaR", "kmi", "knitcitations", "knitr", 
"kohonen", "koRpus", "ks", "labeling", "labeltodendro", "laeken", 
"LaF", "laGP", "Lahman", "lambda.r", "largeVis", "lars", "lasso2", 
"latentnet", "lattice", "latticeExtra", "Laurae", "lava", "lava.tobit", 
"lavaan", "lavaan.survey", "lawstat", "lazyeval", "LCA", "lcopula", 
"leaps", "LearnBayes", "lfda", "lfe", "lhs", "LiblineaR", "likert", 
"linprog", "lintr", "lisrelToR", "lme4", "lmerTest", "lmodel2", 
"lmtest", "loa", "locfit", "logcondens", "LogicReg", "logistf", 
"logspline", "lokern", "longmemo", "loo", "lpSolve", "lpSolveAPI", 
"lqa", "lqmm", "lsmeans", "lubridate", "MAc", "MAd", "magrittr", 
"mail", "manipulate", "mapdata", "mapproj", "maps", "maptools", 
"maptree", "markdown", "MASS", "Matching", "MatchIt", "mathgraph", 
"matlab", "Matrix", "matrixcalc", "MatrixModels", "matrixStats", 
"maxLik", "maxlike", "MBA", "MBESS", "mboost", "mc2d", "mcgibbsit", 
"mclust", "mcmc", "MCMCglmm", "MCMCpack", "mco", "mda", "MDSGUI", 
"mediation", "memisc", "memoise", "MEMSS", "MetABEL", "metafor", 
"methods", "Metrics", "mets", "mgcv", "mi", "mice", "miceadds", 
"microbenchmark", "microplot", "MicrosoftR", "mime", "minerva", 
"miniUI", "minpack.lm", "minqa", "mirt", "mirtCAT", "misc3d", 
"miscTools", "missForest", "missMDA", "mitml", "mitools", "mix", 
"mlbench", "MLmetrics", "mlmRev", "mlogit", "mlr", "mnlogit", 
"mnormt", "modeest", "modeltools", "mondate", "monreg", "moonBook", 
"mosaic", "mosaicData", "movMF", "MplusAutomation", "mpmi", "MPV", 
"mratios", "mRMRe", "msm", "mstate", "MSwM", "muhaz", "multcomp", 
"multcompView", "multicool", "multiwayvcov", "MuMIn", "munsell", 
"mvinfluence", "mvnormtest", "mvoutlier", "mvtnorm", "NbClust", 
"ncdf4", "ncvreg", "ndtv", "network", "networkDynamic", "networkDynamicData", 
"networksis", "neuralnet", "NeuralNetTools", "NHANES", "nlme", 
"nloptr", "NLP", "NMF", "nnet", "nnls", "nodeHarvest", "nor1mix", 
"norm", "nortest", "np", "numDeriv", "nws", "nycflights13", "obliqueRF", 
"odfWeave", "OpenMx", "openssl", "optextras", "optimx", "optmatch", 
"orcutt", "ordinal", "ore", "orloca", "orloca.es", "orthopolynom", 
"outliers", "oz", "packrat", "pamr", "pan", "pander", "parallel", 
"parallelMap", "ParamHelpers", "partitions", "party", "partykit", 
"pastecs", "pbapply", "pbivnorm", "pbkrtest", "PBSmapping", "PBSmodelling", 
"pcalg", "pcaPP", "pec", "penalized", "PerformanceAnalytics", 
"permute", "pgirmess", "pixmap", "pkgKitten", "pkgmaker", "PKI", 
"PKPDmodels", "playwith", "plm", "plot3D", "plotly", "plotmo", 
"plotrix", "pls", "plyr", "PMCMR", "pmml", "pmmlTransformations", 
"png", "poistweedie", "poLCA", "polspline", "polyclip", "polycor", 
"polynom", "prabclus", "pracma", "praise", "PredictABEL", "prefmod", 
"prim", "pROC", "prodlim", "profdpm", "profileModel", "propagate", 
"proto", "proxy", "pryr", "pscl", "pso", "pspline", "psych", 
"psychotools", "psychotree", "purrr", "pvclust", "pwr", "qap", 
"qcc", "qgraph", "qrng", "quadprog", "quantmod", "quantreg", 
"qvcalc", "R.cache", "R.devices", "R.matlab", "R.methodsS3", 
"R.oo", "R.rsp", "R.utils", "R2BayesX", "R2Cuba", "R2HTML", "R2jags", 
"R2OpenBUGS", "R2PPT", "R2wd", "R2WinBUGS", "R6", "RandomFields", 
"RandomFieldsUtils", "randomForest", "randomForestSRC", "randtests", 
"randtoolbox", "ranger", "RankAggreg", "RANN", "RArcInfo", "rARPACK", 
"RaschSampler", "raster", "rasterVis", "rattle", "rbenchmark", 
"rbounds", "rbvs", "Rcgmin", "Rcmdr", "RcmdrMisc", "RcmdrPlugin.BCA", 
"RcmdrPlugin.coin", "RcmdrPlugin.depthTools", "RcmdrPlugin.DoE", 
"RcmdrPlugin.doex", "RcmdrPlugin.epack", "RcmdrPlugin.Export", 
"RcmdrPlugin.FactoMineR", "RcmdrPlugin.HH", "RcmdrPlugin.IPSUR", 
"RcmdrPlugin.KMggplot2", "RcmdrPlugin.mosaic", "RcmdrPlugin.orloca", 
"RcmdrPlugin.pointG", "RcmdrPlugin.qual", "RcmdrPlugin.SLC", 
"RcmdrPlugin.sos", "RcmdrPlugin.steepness", "RcmdrPlugin.survival", 
"RcmdrPlugin.TeachingDemos", "RcmdrPlugin.UCA", "RColorBrewer", 
"Rcpp", "RcppArmadillo", "RcppCNPy", "RcppDE", "RcppEigen", "RcppParallel", 
"RcppProgress", "RcppRoll", "Rcsdp", "RCurl", "readr", "readxl", 
"recommenderlab", "ref", "RefManageR", "registry", "relaimpo", 
"relations", "relax", "relevent", "reliaR", "relimp", "rem", 
"reportr", "reshape", "reshape2", "RevoIOQ", "RevoMods", "RevoUtils", 
"RevoUtilsMath", "rex", "rFerns", "rgdal", "rgenoud", "rgeos", 
"rggobi", "rgl", "Rglpk", "rglwidget", "RgoogleMaps", "rgp", 
"rgpui", "RGtk2", "RGtk2Extras", "RH2", "riskRegression", "RItools", 
"rjags", "rJava", "RJDBC", "rjson", "RJSONIO", "rknn", "rlecuyer", 
"rmarkdown", "rmeta", "Rmpfr", "Rmpi", "rms", "RMySQL", "rneos", 
"rngtools", "rngWELL", "robCompositions", "robust", "robustbase", 
"rockchalk", "ROCR", "RODBC", "Rook", "rootSolve", "rotationForest", 
"roxygen2", "rpanel", "rpart", "rpart.plot", "rpf", "RPostgreSQL", 
"rrcov", "rredis", "RRF", "rrlda", "RSclient", "rsconnect", "Rserve", 
"RSiena", "RSKC", "rsm", "RSNNS", "Rsolnp", "RSpectra", "RSQLite", 
"rstan", "rstanarm", "rstudioapi", "rsvg", "Rsymphony", "rtiff", 
"Rtsne", "Rttf2pt1", "rugarch", "RUnit", "Runuran", "rversions", 
"rvest", "rvg", "Rvmmin", "RWeka", "RWekajars", "Ryacas", "sampleSelection", 
"sampling", "sandwich", "scagnostics", "scales", "scalreg", "scatterplot3d", 
"scrypt", "sda", "selectr", "sem", "semiArtificial", "semPlot", 
"semTools", "sendmailR", "sendplot", "SensoMineR", "seriation", 
"setRNG", "sets", "sfsmisc", "sgeostat", "shape", "shapefiles", 
"shapes", "shiny", "shinyjs", "shinystan", "shinythemes", "signal", 
"SimComp", "SimDesign", "simecol", "simex", "simsem", "sirt", 
"SIS", "sjmisc", "sjPlot", "sjstats", "SkewHyperbolic", "skmeans", 
"slackr", "slam", "SLC", "Sleuth2", "sm", "smbinning", "sn", 
"sna", "snow", "SnowballC", "snowfall", "snowFT", "som", "soma", 
"sos", "sourcetools", "sp", "spacetime", "spam", "sparcl", "sparsediscrim", 
"SparseGrid", "sparseLDA", "SparseM", "sparsity", "spatial", 
"spatstat", "spc", "spd", "spdep", "speedglm", "sphet", "splancs", 
"splines", "splm", "spls", "sqldf", "sROC", "stabledist", "stabs", 
"StanHeaders", "startupmsg", "stashR", "statmod", "statnet", 
"statnet.common", "stats", "stats4", "steepness", "stepPlr", 
"stringdist", "stringi", "stringr", "strucchange", "subselect", 
"subsemble", "sudoku", "SuperLearner", "superpc", "SuppDists", 
"survey", "survival", "svd", "svglite", "svGUI", "svUnit", "svyPVpack", 
"SwarmSVM", "SweaveListingUtils", "systemfit", "tables", "tabplot", 
"tabplotd3", "TAM", "tcltk", "tcltk2", "tclust", "TeachingDemos", 
"tensor", "tensorA", "tergm", "testit", "testthat", "texreg", 
"tgp", "TH.data", "threejs", "tibble", "tidyr", "tikzDevice", 
"timeDate", "timereg", "timeSeries", "tis", "tkrplot", "tm", 
"tmvtnorm", "tnam", "tools", "TransferEntropy", "translations", 
"tree", "trimcluster", "tripack", "truncdist", "truncnorm", "truncreg", 
"trust", "TSA", "tseries", "tseriesEntropy", "tsna", "TSP", "TTR", 
"tufte", "tuneR", "tweedie", "ucminf", "unmarked", "urca", "utils", 
"V8", "VariABEL", "VarianceGamma", "vars", "vcd", "vcdExtra", 
"Vdgraph", "vegan", "verification", "VGAM", "VGAMdata", "VIM", 
"VIMGUI", "VineCopula", "vioplot", "viridis", "viridisLite", 
"visNetwork", "vtreat", "wavelets", "waveslim", "wbstats", "webp", 
"webshot", "WGCNA", "WhatIf", "whisker", "whoami", "withr", "woe", 
"wordcloud", "WrightMap", "WriteXLS", "wskm", "wsrf", "xergm", 
"xergm.common", "xgboost", "xkcd", "XLConnect", "XLConnectJars", 
"XML", "xml2", "xtable", "xts", "YaleToolkit", "yaml", "yarrr", 
"Zelig", "zipcode", "zoo", "ztable")

install.packages(packages)
```

You will encounter a bit of issues on these packages due to the order they must be installed:

```r
install.packages("jsonlite")
install.packages("curl")
install.packages(c("Cubist", "rasterVis", "RSQLite"))
```

Same here. Laurae2/sparsity can not compile under Microsoft R Client. You must use Microsoft R Open or R to install it (Tools > Options to change the R used). Then, copy the /library/sparsity folder and paste it under R Client/R_Server/library. Then, you can use sparsity package for SVMLight conversion =)

```r
install.packages("data.table", type = "source", repos = "http://Rdatatable.github.io/data.table")
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos=NULL, type="source")
devtools:::install_github("Laurae2/Laurae")
devtools:::install_github("Laurae2/woe")
devtools:::install_github("Laurae2/sparsity") # FAILS
devtools:::install_github("twitter/AnomalyDetection")
devtools:::install_github('cmpolis/datacomb', subdir='pkg', ref='1.1.2')
devtools:::install_github(c("ramnathv/htmlwidgets", "smartinsightsfromdata/rpivotTable"))
devtools:::install_github('rstudio/leaflet')
install.packages(c("dygraphs", "networkD3"))
devtools:::install_github('rstudio/DT')
devtools:::install_github("bwlewis/rthreejs")
install.packages(c("addinslist", "colourpicker", "ggExtra", "ggThemeAssist", "questionr", "citr", "QRAGadget", "gitgadget"))
devtools:::install_github("tjmahr/WrapRmd")
devtools:::install_github("MangoTheCat/tidyshiny")
devtools:::install_github("homerhanumat/addinplots")
devtools:::install_github("Stan125/limoaddin")
devtools:::install_github("dkilfoyle/rpivotGadget")
devtools:::install_github("csgillespie/addinmanager")
devtools:::iinstall_github("YvesCR/arimaUI")
devtools:::install_github("benmarwick/wordcountaddin")
devtools:::install_github("Stan125/GREA")
devtools:::install_github("digital-dharma/RStudioAddIns")
devtools:::install_github("ThinkRstat/littleboxes")
install.packages("radiant", repos = "https://radiant-rstats.github.io/minicran/", type = 'binary')
devtools:::install_github("sarupurisailalith/commonUtilAddins")
devtools:::install_github("haozhu233/giphyr")
devtools:::install_github("dracodoc/namebrowser")
devtools:::install_github("lorenzwalthert/strcode")
```

xgboost installation, set to the right directory:

```r
setwd("C:/xgboost/xgboost")
install('R-package')
```

Pre-compile xgboost (here aligned to [PR1855](https://github.com/dmlc/xgboost/pull/1855)) if you didn't:

```bash
cd c:/xgboost
git clone --recursive https://github.com/Laurae2/xgboost
cd xgboost
git submodule init
git submodule update
alias make='mingw32-make'
cd dmlc-core
make
cd ../rabit
make lib/librabit_empty.a
cd ..
cp make/mingw64.mk config.mk
make
```
