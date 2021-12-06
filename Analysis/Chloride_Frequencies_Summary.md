Analysis of LCWMD “Diurnal Exceedences” of Chronic (“CCC”) Chloride
Standards
================
Curtis C. Bohlen, Casco Bay Estuary Partnership
01/12/2021

-   [Introduction](#introduction)
    -   [Are Water Quality Criteria
        Met?](#are-water-quality-criteria-met)
    -   [Sources of Threshold Values](#sources-of-threshold-values)
        -   [Chloride](#chloride)
-   [Import Libraries](#import-libraries)
-   [Data Preparation](#data-preparation)
    -   [Folder References](#folder-references)
    -   [Data on Sites and Impervious
        Cover](#data-on-sites-and-impervious-cover)
    -   [Main Data](#main-data)
    -   [Data Corrections](#data-corrections)
        -   [Anomalous Depth Values](#anomalous-depth-values)
        -   [Single S06B Chloride Observation from
            2017](#single-s06b-chloride-observation-from-2017)
        -   [Anomalous Dissolved Oxygen and Chloride
            Values](#anomalous-dissolved-oxygen-and-chloride-values)
    -   [Remove Partial Data from Winter
        Months](#remove-partial-data-from-winter-months)
    -   [Add Stream Flow Index](#add-stream-flow-index)
    -   [Remove Site S06B, Trim Data](#remove-site-s06b-trim-data)
-   [Exploratory Cross Tabs](#exploratory-cross-tabs)
    -   [Utility Function](#utility-function)
    -   [Chloride Chronic](#chloride-chronic)
    -   [Chloride Acute](#chloride-acute)
-   [Exploratory Graphics](#exploratory-graphics)
-   [GAMM models](#gamm-models)
    -   [Model Sequence](#model-sequence)
        -   [Largest Model: Site, year\_f, month\_f, and
            FlowIndex](#largest-model-site-year_f-month_f-and-flowindex)
        -   [Model 2: Site, Year, and
            FlowIndex](#model-2-site-year-and-flowindex)
        -   [Model 3: Site and Year](#model-3-site-and-year)
        -   [Model 4: Site Only](#model-4-site-only)
        -   [Model 5: Month Only](#model-5-month-only)
-   [Examine Full Model](#examine-full-model)
-   [Extract Marginal Means](#extract-marginal-means)
    -   [Add Calls to GAM objects](#add-calls-to-gam-objects)
    -   [By Site](#by-site)
        -   [Graphics](#graphics)
    -   [By Year](#by-year)
        -   [Graphics](#graphics-1)
    -   [By Month](#by-month)
        -   [Graphic](#graphic)

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Introduction

The Long Creek Watershed, almost three and a half square miles in area,
is dominated by commercial land use. The Maine Mall is one of the
largest landowners in the watershed, and it is surrounded by a range of
commercial businesses, from medical offices, to car washes. About a
third of the watershed in impervious surfaces like roads, parking lots,
and rooftops.

Landowners with an acre or more of impervious area are required to get a
Clean Water Act permit for stormwater discharges from their property.
The LCWMD provides an alternative for landowners to working to receive
an individual permit. Landowners who elect to participate in the The
Long Creek Watershed Management District receive a General Permit, in
return for providing funding to the District, and facilitating the work
of the district by permitting access to their property for certain
activities.

For more information on LCWMD, see [their web
site](restorelongcreek.org).

Over the past decade, LCWMD has contracted with several consulting firms
to provide water quality monitoring services along Long Creek. This has
produced one of the most extensive and best documented data set from the
Northeastern US looking at water quality conditions in an urban stream.

GZA Geoenvironmental Incorporated (GZA) has been the primary monitoring
contractor for LCWMD for several years, and in 2019, they conducted a
thorough review of LCWMD data. These analyses are based on their summary
data sets, and recapitulate and extend their analyses.

## Are Water Quality Criteria Met?

The primary question we ask in this Notebook, is whether water quality
criteria pertaining to levels of chloride are met. In particular, we
explore various ways of modeling those probabilities

We ask whether the probability of failing to meet criteria each day is
changing. Secondarily, we examine differences among sites in the
probability of failing criteria.

In this data set a “TRUE” value consistently implies that water quality
criteria were met or exceeded, whether that is achieved by a value
higher than or lower than some numeric criteria. “TRUE” implies good
conditions. “FALSE” implies bad conditions.

## Sources of Threshold Values

### Chloride

Maine uses established thresholds for both chronic and acute exposure to
chloride. These are the “CCC and CMC” standards for chloride in
freshwater. (06-096 CMR 584). These terms are defined in a footnote as
follows:

> The Criteria Maximum Concentration (CMC) is an estimate of the highest
> concentration of a material in surface water to which an aquatic
> community can be exposed briefly without resulting in an unacceptable
> effect. The Criterion Continuous Concentration (CCC) is an estimate of
> the highest concentration of a material in surface water to which an
> aquatic community can be exposed indefinitely without resulting in an
> unacceptable effect.

The relevant thresholds are:

-   Chloride CCC = 230 mg/l
-   Chloride CMC = 860 mg/l

In practice, chloride in Long Creek are indirectly estimated based on
measurement of conductivity. The chloride-conductivity correlations is
fairly close and robust, but estimation is an additional source of
error, although generally on the level of 10% or less.

Further, exceedances of the acute chloride threshold, the “CMC” are
relatively rare, which limits ability to construct robust models. The
classification tree models in
“Classification\_Tree\_Models\_Summary.Rmd” are perhaps more successful.

# Import Libraries

``` r
library(glmmTMB)  # Generalized Linear Mixed models.
#> Warning: package 'glmmTMB' was built under R version 4.0.5
#> Warning in checkMatrixPackageVersion(): Package version inconsistency detected.
#> TMB was built with Matrix version 1.3.4
#> Current Matrix version is 1.2.18
#> Please re-install 'TMB' from source using install.packages('TMB', type = 'source') or ask CRAN for a binary version of 'TMB' matching CRAN's 'Matrix' package
library(mgcv)     # For mixed effects GAMMs
#> Warning: package 'mgcv' was built under R version 4.0.5
#> Loading required package: nlme
#> This is mgcv 1.8-38. For overview type 'help("mgcv-package")'.

library(tidyverse)# Has to load after MASS, so `select()` is not masked
#> Warning: package 'tidyverse' was built under R version 4.0.5
#> -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
#> v ggplot2 3.3.5     v purrr   0.3.4
#> v tibble  3.1.6     v dplyr   1.0.7
#> v tidyr   1.1.4     v stringr 1.4.0
#> v readr   2.1.0     v forcats 0.5.1
#> Warning: package 'ggplot2' was built under R version 4.0.5
#> Warning: package 'tidyr' was built under R version 4.0.5
#> Warning: package 'dplyr' was built under R version 4.0.5
#> Warning: package 'forcats' was built under R version 4.0.5
#> -- Conflicts ------------------------------------------ tidyverse_conflicts() --
#> x dplyr::collapse() masks nlme::collapse()
#> x dplyr::filter()   masks stats::filter()
#> x dplyr::lag()      masks stats::lag()
library(readr)

library(emmeans)  # Provides tools for calculating marginal means
#> Warning: package 'emmeans' was built under R version 4.0.5

library(CBEPgraphics)
load_cbep_fonts()
theme_set(theme_cbep())

library(LCensMeans)
```

# Data Preparation

## Folder References

``` r
sibfldnm    <- 'Data'
parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)

dir.create(file.path(getwd(), 'figures'), showWarnings = FALSE)
dir.create(file.path(getwd(), 'models'), showWarnings = FALSE)
```

## Data on Sites and Impervious Cover

These data were derived from Table 2 from a GZA report to the Long Creek
Watershed Management District, titled “Re: Long Creek Watershed Data
Analysis; Task 2: Preparation of Explanatory and Other Variables.” The
Memo is dated November 13, 2019 File No. 09.0025977.02.

Cumulative Area and IC calculations are our own, based on the GZA data
and the geometry of the stream channel.

``` r
# Read in data and drop the East Branch, where we have no data
fn <- "Site_IC_Data.csv"
fpath <- file.path(sibling, fn)

Site_IC_Data <- read_csv(fpath) %>%
  filter(Site != "--") 
#> Rows: 7 Columns: 8
#> -- Column specification --------------------------------------------------------
#> Delimiter: ","
#> chr (4): Site, Subwatershed, PctIC, CumPctIC
#> dbl (4): Area_ac, IC_ac, CumArea_ac, CumIC_ac
#> 
#> i Use `spec()` to retrieve the full column specification for this data.
#> i Specify the column types or set `show_col_types = FALSE` to quiet this message.

# Now, create a factor that preserves the order of rows (roughly upstream to downstream). 
Site_IC_Data <- Site_IC_Data %>%
  mutate(Site = factor(Site, levels = Site_IC_Data$Site))

# Finally, convert percent covers to numeric values
Site_IC_Data <- Site_IC_Data %>%
  mutate(CumPctIC = as.numeric(substr(CumPctIC, 1, nchar(CumPctIC)-1))) %>%
  mutate(PctIC = as.numeric(substr(PctIC, 1, nchar(PctIC)-1)))
Site_IC_Data
#> # A tibble: 6 x 8
#>   Site  Subwatershed      Area_ac IC_ac CumArea_ac CumIC_ac PctIC CumPctIC
#>   <fct> <chr>               <dbl> <dbl>      <dbl>    <dbl> <dbl>    <dbl>
#> 1 S07   Blanchette Brook     434.  87.7       434.     87.7  20.2     20.2
#> 2 S06B  Upper Main Stem      623.  80.2       623.     80.2  12.9     12.9
#> 3 S05   Middle Main Stem     279.  53.6      1336     222.   19.2     16.6
#> 4 S17   Lower Main Stem      105   65.1      1441     287.   62       19.9
#> 5 S03   North Branch Trib    298. 123         298.    123    41.2     41.2
#> 6 S01   South Branch Trib    427. 240.        427.    240.   56.1     56.1
```

## Main Data

We removed the 2019 data, as we don’t have a complete year’s worth of
data, which may bias annual summaries.

Note that this data does *not* include all of the predictors used in
some models looking at chlorides. In particular, it does not include
stream flow estimates.

``` r
fn <- "Exceeds_Data.csv"
exceeds = read_csv(file.path(sibling, fn), progress=FALSE) %>%
  mutate(IC=Site_IC_Data$CumPctIC[match(Site, Site_IC_Data$Site)]) %>%
  filter(Year < 2019) %>%
  mutate(Site = factor(Site, levels=levels(Site_IC_Data$Site)),
         year_f = factor(Year),
         month_f = factor(Month, levels = 1:12, labels = month.abb),
         DOY = as.numeric(format(sdate, format = '%j')),
         season = cut(Month, breaks = c(0,2,5,8,11,13),
                      labels = c('Winter', 'Spring',
                                 'Summer', 'Fall', 'Winter')),
         season = factor(season, levels = c('Winter', 'Spring', 
                                           'Summer', 'Fall'))) %>%
  mutate(lPrecip = log1p(Precip))
#> New names:
#> * `` -> ...1
#> Rows: 11422 Columns: 19
#> -- Column specification --------------------------------------------------------
#> Delimiter: ","
#> chr   (1): Site
#> dbl   (7): ...1, Year, Month, Precip, PPrecip, MaxT, D_Median
#> lgl  (10): ClassCDO, ClassBDO, ClassC_PctSat, ClassB_PctSat, ClassCBoth, Cla...
#> date  (1): sdate
#> 
#> i Use `spec()` to retrieve the full column specification for this data.
#> i Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

## Data Corrections

### Anomalous Depth Values

Several depth observations in the record appear highly unlikely. In
particular, several observations show daily median water depths over 15
meters. A few other observations show daily median depths over 4 meters,
which also looks unlikely in a stream of this size. All these events
also occurred in May or June of 2015 at site S05. Some sort of
malfunction of the pressure transducer appears likely.

We can trace these observations back to the raw QA/QC’d pressure and
sonde data submitted to LCWMD by GZA, so they are not an artifact of our
data preparation.

We remove these extreme values. The other daily medians in May and June
of 2015 appear reasonable, and we leave them in place, although given
possible instability of the pressure sensors, it might make sense to
remove them all.

Note that removing depth observations from Site S05 will remove those
dates from any model that uses the `FlowIndex` variable (see below) as a
predictor.

``` r
exceeds <- exceeds %>%
  mutate(D_Median = if_else(D_Median > 4, NA_real_, D_Median),
         lD_Median = log1p(D_Median))
```

### Single S06B Chloride Observation from 2017

The data includes just a single chloride observation from site S06B from
any year other than 2013. While we do not know if the data point is
legitimate or not, it has high leverage in several models, and we
suspect a transcription error of some sort.

We remove the Chloride value from the data.

``` r
exceeds <- exceeds %>%
  mutate(ChlCCC = if_else(Site == 'S06B' & Year > 2014,
                              NA, ChlCCC),
         ChlCMC = if_else(Site == 'S06B' & Year > 2014,
                              NA, ChlCMC))
```

### Anomalous Dissolved Oxygen and Chloride Values

We noted extreme dissolved oxygen data at Site S03, during the end of
2016. Values were both extreme and highly variable. (See discussion in
the DO Analysis workbooks). We decided we should remove both the
chloride and oxygen observations after October 15th.

``` r
exceeds <- exceeds %>% 
  mutate(ChlCCC = if_else(Year == 2016 & Site == 'S03' & DOY > 288,
                              NA, ChlCCC),
         ChlCMC = if_else(Year == 2016 & Site == 'S03' & DOY > 288,
                              NA, ChlCMC))
```

## Remove Partial Data from Winter Months

We have very limited data from several months. We have January data from
only one year, and February data from only two, and December data from
only four years, all older. Both March and November sample sizes vary.

The limited winter data generates severely unbalanced samples, which may
lead to estimation problems, especially in models with crossed or
potentially crossed factors and predictors. More fundamentally, the
potential bias introduced by showing data from those months from just a
handful of years could give a misleading impression of seasonal
patterns. We trim December, January and February data, but leave the
other months.

It is important to remember, even after trimming the data, that:  
1. 2010 is a partial year,  
2. The period of sampling in March may be biased due to spring melt
timing.

``` r
xtabs(~ year_f + month_f, data = exceeds)
#>       month_f
#> year_f Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec
#>   2010   0   0   0   0   0  69  97 103 120 124 120  35
#>   2011   0  15 101 120 124 120 124 124 120 124 120 112
#>   2012   0  39  93  90  93 113 124 124 120  39  96 124
#>   2013   9   0  46 128 155 140 124 124 120 138 150  15
#>   2014   0   0  53 102 155 150 155 155 150 155 120   0
#>   2015   0   0   8 141 186 180 186 186 180 160  30   0
#>   2016   0   0  10 170 186 180 186 186 180 186 168   0
#>   2017   0   0 186 180 186 180 186 186 180 186 102   0
#>   2018   0   0  16 180 186 180 186 186 180 186 126   0
```

``` r
exceeds <- exceeds %>%
  filter(Month >= 3 & Month <= 11)
```

## Add Stream Flow Index

We worked through many models on a site by site basis in which we
included data on water depth, but since the depth coordinate is
site-specific, a 10 cm depth at one site may be exceptional, while at
another it is commonplace. We generally want not a local measure of
stream depth, but a watershed-wide metric of high, medium, or low stream
flow.

Middle and Lower Maine Stem sites would be suitable for a general flow
indicator across the watershed. The monitoring sites in that stretch of
Long Creek include include S05 and S17, however only site S05 has been
in continuous operation throughout the period of record, so we use depth
data from S05 to construct our general stream flow indicator.

Stream flow at S05 is correlated with flow at other sites, although not
all that closely correlated to flow in the downstream tributaries (S01
and S03).

``` r
exceeds %>%
  select(sdate, Site, lD_Median) %>%
  pivot_wider(names_from = Site, values_from = lD_Median) %>%
  select( -sdate) %>%
  cor(use = 'pairwise', method = 'pearson')
#>            S01       S03       S05      S06B       S07       S17
#> S01  1.0000000 0.4499047 0.6506630 0.6310361 0.5594067 0.7290077
#> S03  0.4499047 1.0000000 0.4526392 0.7152403 0.4578906 0.6666414
#> S05  0.6506630 0.4526392 1.0000000 0.8043943 0.7042711 0.7906571
#> S06B 0.6310361 0.7152403 0.8043943 1.0000000 0.5882527 0.8778188
#> S07  0.5594067 0.4578906 0.7042711 0.5882527 1.0000000 0.7327432
#> S17  0.7290077 0.6666414 0.7906571 0.8778188 0.7327432 1.0000000
```

We use the log of the daily median flow at S05 as a general
watershed-wide stream flow indicator, which we call `FlowIndex`. We use
the log of the raw median, to lessen the effect of the highly skewed
distribution of stream depths on the metric.

``` r
depth_data <- exceeds %>%
  filter (Site == 'S05') %>%
  select(sdate, lD_Median)

exceeds <- exceeds %>%
  mutate(FlowIndex = depth_data$lD_Median[match(sdate, depth_data$sdate)])

rm(depth_data)
```

## Remove Site S06B, Trim Data

Site S06B only has chloride data from a single year, so including it in
temporal models causes problems. We remove the Site from further
analysis. We also drop variables we will not analyze further in this
Notebook.

``` r
exceeds <- exceeds %>%
  filter(Site != 'S06B') %>%
  select(-starts_with('Class')) %>%
  select(-contains('T_ex')) %>%
  mutate(Site = droplevels(Site))
```

# Exploratory Cross Tabs

## Utility Function

This function just adds a percent summary column to a cross-tab.

``` r
xt_pct <- function(.form, .dat) {
  xt <- xtabs(.form, data = .dat)
  xt <- cbind(xt, round(apply(xt, 1, function(X) X[1]/sum(X)), 3)*100)
  names(xt[3]) <- 'Percent Fail'
  return(xt)
}
```

## Chloride Chronic

``` r
xt_pct(~Year + ChlCCC, exceeds)
#>      FALSE TRUE     
#> 2010   334  247 57.5
#> 2011   533  197 73.0
#> 2012   444  158 73.8
#> 2013   529  324 62.0
#> 2014   396  543 42.2
#> 2015   563  485 53.7
#> 2016   988  160 86.1
#> 2017   714  399 64.2
#> 2018   813  374 68.5
```

But the Strong pattern is by sites.

``` r
xt_pct(~Site + ChlCCC, exceeds)
#>     FALSE TRUE     
#> S07   935  945 49.7
#> S05   394 1016 27.9
#> S17   598  260 69.7
#> S03  1544  501 75.5
#> S01  1843  165 91.8
```

## Chloride Acute

``` r
xt_pct(~Year + ChlCMC, exceeds)
#>      FALSE TRUE    
#> 2010     1  580 0.2
#> 2011    17  713 2.3
#> 2012     6  596 1.0
#> 2013     4  849 0.5
#> 2014     5  934 0.5
#> 2015    10 1038 1.0
#> 2016    77 1071 6.7
#> 2017    51 1062 4.6
#> 2018    61 1126 5.1
```

``` r
xt_pct(~Site + ChlCMC, exceeds)
#>     FALSE TRUE    
#> S07    27 1853 1.4
#> S05     2 1408 0.1
#> S17     1  857 0.1
#> S03    37 2008 1.8
#> S01   165 1843 8.2
```

Note that those probabilities are very low – under 1% – which is likely
to limit ability to fit models to these data.

# Exploratory Graphics

These are estimated as empirical relative frequencies, with error
estimated as two times the standard error of the estimate.

``` r
exceeds %>%
  group_by(Site, Year) %>%
  summarize(CCC_true = sum(ChlCCC, na.rm = TRUE),
            CCC_count = sum(! is.na(ChlCCC)),
            CCC_p = CCC_true/CCC_count,
            CCC_err = CCC_p*(1-CCC_p)/sqrt(CCC_count),
            .groups = 'drop') %>%
  ggplot(aes(Year, CCC_p, color = Site)) +
  geom_line() +
  geom_pointrange(aes(ymin = CCC_p-2 * CCC_err, ymax = CCC_p + 2 * CCC_err)) +
  ylab('Probability of Passing\nChloride CCC Standard')
#> Warning: Removed 1 rows containing missing values (geom_pointrange).
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/Site_exploratory_graphic-1.png" style="display: block; margin: auto;" />
2016 was a rough year at most sites. Sites S01 and S03 fail this
standard frequently.

Note that for some years at site S05, we never had a failure to meet
this chloride standard. This will probably further limit the models we
can fit.

``` r
exceeds  %>%
  group_by(month_f, Year) %>%
  summarize(CCC_true = sum(ChlCCC, na.rm = TRUE),
            CCC_count = sum(! is.na(ChlCCC)),
            CCC_p = CCC_true/CCC_count,
            CCC_err = CCC_p*(1-CCC_p)/sqrt(CCC_count)) %>%
  ggplot(aes(Year, CCC_p, color = month_f)) +
  geom_line() +
  geom_pointrange(aes(ymin = CCC_p-2 * CCC_err, ymax = CCC_p + 2 * CCC_err))
#> `summarise()` has grouped output by 'month_f'. You can override using the `.groups` argument.
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/month_exploritory_graphic-1.png" style="display: block; margin: auto;" />
That shows that 2011 and 2016 were fairly bad years. But the month to
month patterns are less obvious than the year to year patterns This
highlights the role of relativity slow-dynamic processes (drought?
groundwater flow?) in shaping chloride conditions in Long Creek.

# GAMM models

Fitting full GAMM models proves to be problematic, as models took a long
time to fit. The big slow down is the correlation structure. It speeds
things slightly to subdivide the sections in which we expect the
`corAR1()` to apply.

We need to include either a smoothed term or a random factor in the
model to use GAMM, even with a correlation term. That is problematic, as
none of our predictors are good candidates for random factors. Sites
could be considered a random selection from all possible sites, or Year
could be considered a random selection from all possible years, but we
are interested in estimates of year by year and site by site
differences.

We replace actual dates with sequential integers, simplifying model
fitting. Although R `Date` objects are integers under the hood, limiting
possible values to a smaller range was essential for fitting models with
`glmmTMB()`, and it appears to speed model fitting. We make the
transformation via factors (which are also integers under the hood).

``` r
first_date <- min(exceeds$sdate)
last_date <- max(exceeds$sdate)

exceeds <- exceeds %>%
  mutate(sdate_f = factor(as.numeric(sdate), 
                          levels = as.numeric(first_date):as.numeric(last_date),
         labels = as.character(seq(from = first_date,
                                       to = last_date,
                                       by = 1))))
```

## Model Sequence

### Largest Model: Site, year\_f, month\_f, and FlowIndex

We start with a fairly large model, and wait for it to complete. This
took slightly more than an hour and a half (5500 seconds) of
computational time, as shown by `system.time()`, but the code took over
two hours to actually run.

The code saves the model once it runs, but the model is not included in
the GitHub repository. The code does not detect changes in the
underlying data, so users need to delete the models by hand if they need
to be run again on new data.

``` r
if (! file.exists("models/ccc_gamm_1.rds")) {
  print(system.time(
    ccc_gamm_1<- gamm(ChlCCC ~ Site + year_f + month_f + s(FlowIndex),
                      correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                      family = 'binomial',
                      niterPQL = 20,
                      data = exceeds)
  ))
  saveRDS(ccc_gamm_1, file="models/ccc_gamm_1.rds")
} else {
  ccc_gamm_1 <- readRDS("models/ccc_gamm_1.rds")
}
```

Marginal means from that model are likely to suffer from the
Hauck-Donner effect, since many cells in the model structure never fail
water quality criteria. Estimation from this model is therefore probably
pointless, but it can act as the largest model in a hierarchical model
sequence.

### Model 2: Site, Year, and FlowIndex

A model that omits the Month term also takes some time to fit. Marginal
means from this model are potentially informative (if fit at median
watershed flows). We are almost certainly misrepresenting the behavior
of Site S03 in this model. Site S03 appears to be doing something quite
different from the other sites, so we should consider omitting
`Site == S03` from models entirely, but we do not do so here.

This model takes not quite as long to fit, but you’ll want to go do
something else for a while while this code runs. It takes a bit over an
hour to run.

``` r
if (! file.exists("models/ccc_gamm_2.rds")) {
  print(system.time(
    ccc_gamm_2<- gamm(ChlCCC ~ Site + year_f + s(FlowIndex),
                      correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                      family = 'binomial',
                      niterPQL = 20,
                      data = exceeds)
  ))
  saveRDS(ccc_gamm_2, file="models/ccc_gamm_2.rds")
} else {
  ccc_gamm_2 <- readRDS("models/ccc_gamm_2.rds")
}
```

In order to fit a model without the flow terms, we need to select one of
our other predictors to treat as a random factor. We see substantial
differences in probability of failing these standards by Site, so we try
some models that treat `year_f` as a random factor. You could fit a
random factor two ways, using the `random = ....` function parameter, or
by including a random effects smooth term via `s(..., type = 're')`.

### Model 3: Site and Year

``` r
if (! file.exists("models/ccc_gamm_3.rds")) {
  print(system.time(
    ccc_gamm_3<- gamm(ChlCCC ~ Site + Year, random = list(year_f = ~ 1),
                      correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                      family = 'binomial',
                      niterPQL = 20, verbosePQL = TRUE,
                      data = exceeds)
  ))
  saveRDS(ccc_gamm_3, file="models/ccc_gamm_3.rds")
} else {
  ccc_gamm_3 <- readRDS("models/ccc_gamm_3.rds")
}
```

### Model 4: Site Only

``` r
if (! file.exists("models/ccc_gamm_4.rds")) {
  print(system.time(
    ccc_gamm_4<- gamm(ChlCCC ~ Site, random = list(year_f = ~ 1),
                      correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                      family = 'binomial',
                      niterPQL = 20,
                      data = exceeds)
  ))
  saveRDS(ccc_gamm_4, file="models/ccc_gamm_4.rds")
} else {
  ccc_gamm_4 <- readRDS("models/ccc_gamm_4.rds")
}
```

### Model 5: Month Only

``` r
if (! file.exists("models/ccc_gamm_5.rds")) {
  print(system.time(
    ccc_gamm_5<- gamm(ChlCCC ~ Site + month_f, random = list(year_f = ~ 1),
                      correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                      family = 'binomial',
                      niterPQL = 20,
                      data = exceeds)
  ))
  saveRDS(ccc_gamm_5, file="models/ccc_gamm_5.rds")
} else {
  ccc_gamm_5 <- readRDS("models/ccc_gamm_5.rds")
}
```

# Examine Full Model

``` r
anova(ccc_gamm_1$gam)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> ChlCCC ~ Site + year_f + month_f + s(FlowIndex)
#> 
#> Parametric Terms:
#>         df     F  p-value
#> Site     4 62.12  < 2e-16
#> year_f   8 11.56  < 2e-16
#> month_f  8  6.26 4.27e-08
#> 
#> Approximate significance of smooth terms:
#>                edf Ref.df     F p-value
#> s(FlowIndex) 6.781  6.781 38.45  <2e-16
```

All model terms appear important, Both Month and Year are fit with
factors, which detects variation, but makes it difficult to draw clear
conclusions. (But see below for consideration of whether there is a
long-term trend or not).

``` r
summary(ccc_gamm_1$gam)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> ChlCCC ~ Site + year_f + month_f + s(FlowIndex)
#> 
#> Parametric coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.04083    0.53793   0.076 0.939497    
#> SiteS05      1.27590    0.25685   4.967 6.97e-07 ***
#> SiteS17     -0.31232    0.39785  -0.785 0.432472    
#> SiteS03     -1.29374    0.25089  -5.157 2.59e-07 ***
#> SiteS01     -2.20461    0.24714  -8.921  < 2e-16 ***
#> year_f2011  -1.40966    0.45357  -3.108 0.001893 ** 
#> year_f2012  -1.09927    0.51175  -2.148 0.031748 *  
#> year_f2013  -0.63869    0.56586  -1.129 0.259070    
#> year_f2014  -0.57285    0.43913  -1.304 0.192113    
#> year_f2015  -0.84586    0.49539  -1.707 0.087789 .  
#> year_f2016  -2.87552    0.43701  -6.580 5.10e-11 ***
#> year_f2017  -1.35453    0.43049  -3.146 0.001660 ** 
#> year_f2018  -1.53858    0.42966  -3.581 0.000345 ***
#> month_fApr   0.30509    0.32980   0.925 0.354968    
#> month_fMay   0.64513    0.35828   1.801 0.071812 .  
#> month_fJun   1.19600    0.36085   3.314 0.000924 ***
#> month_fJul   1.52115    0.37028   4.108 4.04e-05 ***
#> month_fAug   1.72782    0.37057   4.663 3.19e-06 ***
#> month_fSep   1.70716    0.36830   4.635 3.64e-06 ***
#> month_fOct   1.83442    0.35645   5.146 2.74e-07 ***
#> month_fNov   1.18372    0.34878   3.394 0.000693 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>                edf Ref.df     F p-value    
#> s(FlowIndex) 6.781  6.781 38.45  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.342   
#>   Scale est. = 1         n = 6160
```

Some sites and some years are just plain worse than others. Chloride
levels are less likely to meet threshold in the spring and more likely
to later in the year.

Note that the model fits a fairly complex smoothed term – with over six
effective degrees of freedom. A simpler model would capture most of the
important detail.

``` r
plot(ccc_gamm_1$gam)
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/plot_gamm_1-1.png" style="display: block; margin: auto;" />
\# Is there a long term trend? Model three included a linear year term
and treated years as a random factor, thus effectively testing to see if
there is a long term trend in probability of meeting chloride standards
if we pool results by year.

``` r
anova(ccc_gamm_3$gam)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> ChlCCC ~ Site + Year
#> 
#> Parametric Terms:
#>      df      F p-value
#> Site  4 97.174  <2e-16
#> Year  1  2.491   0.115
```

So, when you consider the years as random variables, there is no
long-term trend in probability to pass chlorides standards. Thus a
long-term trend is present (in the full model) but whether you consider
it important or not depends on model specifications.

# Extract Marginal Means

## Add Calls to GAM objects

``` r
the_call <-  quote(gamm(ChlCCC ~ Site + year_f + month_f + s(FlowIndex),
                       correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                       family = 'binomial',
                       niterPQL = 20,
                       data = exceeds))
ccc_gamm_1$gam$call <- the_call


the_call <-  quote(gamm(ChlCCC ~ Site + year_f + s(FlowIndex),
                       correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                       family = 'binomial',
                       niterPQL = 20,
                       data = exceeds))
ccc_gamm_2$gam$call <- the_call


the_call <-  quote(gamm(ChlCCC ~ Year + Site, random = list(year_f = ~ 1),
                       correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                       family = 'binomial',
                       niterPQL = 20, verbosePQL = TRUE,
                       data = exceeds))
ccc_gamm_3$gam$call <- the_call


the_call <-  quote(gamm(ChlCCC ~ Site, random = list(year_f = ~ 1),
                       correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                       family = 'binomial',
                       niterPQL = 20, verbosePQL = TRUE,
                       data = exceeds))
ccc_gamm_4$gam$call <- the_call

the_call <-  quote(gamm(ChlCCC ~ Site + month_f, random = list(year_f = ~ 1),
                       correlation = corAR1(form = ~ as.numeric(sdate_f) | Site),
                       family = 'binomial',
                       niterPQL = 20,
                       data = exceeds))
ccc_gamm_5$gam$call <- the_call
```

## By Site

``` r
my_ref_grid <- ref_grid(ccc_gamm_1,  cov.reduce = median) 
a <- summary(emmeans(my_ref_grid, ~ Site, 
                      type = 'response'))

my_ref_grid <- ref_grid(ccc_gamm_2,  cov.reduce = median)
b <-  summary(emmeans(my_ref_grid, ~ Site,
                       type = 'response'))

my_ref_grid <- ref_grid(ccc_gamm_3,  cov.reduce = median)
c <-  summary(emmeans(my_ref_grid, ~ Site,
                       type = 'response'))

my_ref_grid <- ref_grid(ccc_gamm_4,  cov.reduce = median) 
d <-  summary(emmeans(my_ref_grid, ~ Site, 
                       type = 'response'))

my_ref_grid <- ref_grid(ccc_gamm_5,  cov.reduce = median) 
e <-  summary(emmeans(my_ref_grid, ~ Site, 
                       type = 'response'))
```

``` r
observed <- exceeds %>%
  select(Site, ChlCCC) %>%
  group_by(Site) %>%
  summarize(Site = first(Site),
              n = sum(! is.na(ChlCCC)), 
              prob = mean(ChlCCC, na.rm = TRUE),
              SE = sqrt((prob * (1-prob))/n),
              lower.CL = prob - 1.96 * SE,
              upper.CL = prob + 1.96 * SE,
              .groups = 'drop')
```

``` r
z <- tibble(Site = factor(levels(exceeds$Site), levels = levels(exceeds$Site)),
            observed = observed$prob, 
            mod_1 = a$prob, mod_2 = b$prob,
            mod_3 = c$prob, mod_4 = d$prob)


z %>%
  pivot_longer(-c(Site, observed), names_to = 'Model', 
               values_to = 'Probability') %>%

  ggplot(aes(observed, Probability, color = Model)) +
  geom_point() +
  geom_line(aes(group = Model)) +
  annotate(geom = 'text', x = observed$prob,
           y = a$prob + .15, label = z$Site) +
  
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  xlab("Observed Value") +
  ylab('Model Predictions (Adjusted)') +
  xlim(0,1) +
  ylim(0,1)
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/compare_models_by_site_1-1.png" style="display: block; margin: auto;" />
So all the models do pretty good job predicting the observed values,
except at site S17, where predictions are consistently higher than what
we observed, presumably because the record for S17 is short, including
only data from more recent years.

``` r
rm(z)
```

Models 1 and 3 generally provide nearly identical predictions, but for
all practical purposes, the four models are all saying the same thing.

We prefer models 1, for access to data on year by year and month by
month patterns, and model 4, for its simplicity.

### Graphics

#### Model 1

``` r
s <- a %>% 
  mutate(fprob = 1-prob,
         fUCL = 1 - lower.CL,
         fLCL = 1 - upper.CL)

plt1 <- ggplot(s, aes(Site, fprob)) +
 
  geom_pointrange(aes(ymin = fLCL, ymax = fUCL),
                color = cbep_colors()[1], size = .75) +
  
  ylab('Probability of Exceeding\nChronic Chloride Standard') +
  xlab('Upstream      Maine Stem                Downstream      ') +

  theme_cbep(base_size = 12) +
  theme(axis.title.x = element_text(size = 9)) +
  ylim(0,1)
```

``` r
plt1 + 
  stat_summary(geom = 'point',
               mapping = aes(x = Site, y = as.numeric(ChlCCC)),
               data = exceeds, fun = function(x) 1 - mean(x, na.rm = TRUE),
               size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'point', x = 3.5, y = .3,
           size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'text', x = 3.75, y = .3, label = 'Observed', hjust = 0,
             size =3.5) +

    annotate(geom = 'point', x = 3.5, y = .25,
           size = 2, color = cbep_colors()[1]) +
    annotate(geom = 'text', x = 3.75, y = .25, label = 'Adjusted', hjust = 0,
             size = 3.5)
#> Warning: Removed 1067 rows containing non-finite values (stat_summary).
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/plot_model_1_site_observed-1.png" style="display: block; margin: auto;" />

#### Model 4

``` r
s <- d %>% 
  mutate(fprob = 1-prob,
         fUCL = 1 - lower.CL,
         fLCL = 1 - upper.CL)

plt2 <- ggplot(s, aes(Site, fprob)) +
 
  geom_pointrange(aes(ymin = fLCL, ymax = fUCL),
                color = cbep_colors()[1], size = .75) +
  
  ylab('Probability of Exceeding\nChronic Chloride Standard') +
  xlab('Upstream   Maine Stem       Downstream') +

  theme_cbep(base_size = 12) +
  theme(axis.title.x = element_text(size = 10)) +
  ylim(0,1)
```

``` r
plt2 + 
  stat_summary(geom = 'point',
               mapping = aes(x = Site, y = as.numeric(ChlCCC)),
               data = exceeds, fun = function(x) 1 - mean(x, na.rm = TRUE),
               size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'point', x = 3.5, y = .3,
           size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'text', x = 3.75, y = .3, label = 'Observed', hjust = 0,
             size = 3.5) +

    annotate(geom = 'point', x = 3.5, y = .25,
           size = 2, color = cbep_colors()[1]) +
    annotate(geom = 'text', x = 3.75, y = .25, label = 'Adjusted', hjust = 0,
             size = 3.5)
#> Warning: Removed 1067 rows containing non-finite values (stat_summary).
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/plot_model_4_site_observed-1.png" style="display: block; margin: auto;" />

## By Year

Only Model 1 and Model 2 tracked Year as a Fixed Factor, so we only have
two models to compare. (Model 3 treated year as a regressor and random
factor).

``` r
my_ref_grid <- ref_grid(ccc_gamm_1,  cov.reduce = median) 
a <- summary(emmeans(my_ref_grid, ~ year_f, type = 'response'))

my_ref_grid <- ref_grid(ccc_gamm_2,  cov.reduce = median) 
b <- summary(emmeans(my_ref_grid, ~ year_f, type = 'response'))
```

``` r
observed <- exceeds %>%
  select(year_f, ChlCCC) %>%
  group_by(year_f) %>%
  summarize(year_f = first(year_f),
              n = sum(! is.na(ChlCCC)), 
              prob = mean(ChlCCC, na.rm = TRUE),
              SE = sqrt((prob * (1-prob))/n),
              lower.CL = prob - 1.96 * SE,
              upper.CL = prob + 1.96 * SE,
              .groups = 'drop')
```

``` r
z <- tibble(Year = factor(levels(exceeds$year_f), 
                          levels = levels(exceeds$year_f)),
            observed = observed$prob, 
            mod_1 = a$prob, mod_2 = b$prob )

z %>%
  pivot_longer(-c(Year, observed), names_to = 'Model', 
               values_to = 'Probability') %>%

  ggplot(aes(observed, Probability, color = Model)) +
  geom_point() +
  geom_line(aes(group = Model)) +
  annotate(geom = 'text', x = observed$prob,
           y = a$prob + .15, label = z$Year) +
  
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  xlab("Observed Value") +
  ylab('Model Predictions (Adjusted)') +
  xlim(0,1) +
  ylim(0,1)
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/create_plotting_df-1.png" style="display: block; margin: auto;" />
Here our marginal means do not show up so well. The model predictions
are highly correlated, but neither model does very well reproducing
observed frequencies, at least for 2010 and 2012. 2010 is a year with no
data from the early part of the year, so marginal means from either
model are biased. It’s not clear why estimates from 2012 are also on the
high side.

``` r
rm(z)
```

### Graphics

#### Model 1

``` r
s <- a %>% 
  mutate(fprob = 1-prob,
         fUCL = 1 - lower.CL,
         fLCL = 1 - upper.CL)

plt1 <- ggplot(s, aes(year_f, fprob)) +
 
  geom_pointrange(aes(ymin = fLCL, ymax = fUCL),
                color = cbep_colors()[1], size = .75) +
  
  ylab('Probability of Exceeding\nChronic Chloride Standard') +
  xlab('Upstream      Maine Stem                Downstream      ') +

  theme_cbep(base_size = 12) +
  theme(axis.title.x = element_text(size = 9)) +
  ylim(0,1)
```

``` r
xanchor = 7
space = 0.3

plt1 + 
  stat_summary(geom = 'point',
               mapping = aes(x = year_f, y = as.numeric(ChlCCC)),
               data = exceeds, fun = function(x) 1 - mean(x, na.rm = TRUE),
               size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'point', x = xanchor, y = .3,
           size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'text', x = xanchor + space, y = .3, label = 'Observed', hjust = 0,
             size = 3.5) +

    annotate(geom = 'point', x = xanchor, y = .25,
           size = 2, color = cbep_colors()[1]) +
    annotate(geom = 'text', x = xanchor + space, y = .25, label = 'Adjusted', hjust = 0,
             size = 3.5)
#> Warning: Removed 1067 rows containing non-finite values (stat_summary).
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

#### Model 2

``` r
s <- b %>% 
  mutate(fprob = 1-prob,
         fUCL = 1 - lower.CL,
         fLCL = 1 - upper.CL)
  
plt2 <- ggplot(s, aes(year_f, fprob)) +
 
  geom_pointrange(aes(ymin = fLCL, ymax = fUCL),
                color = cbep_colors()[1], size = .75) +
  
  ylab('Probability of Exceeding\nChronic Chloride Standard') +
  xlab('') +

  theme_cbep(base_size = 12) +
  theme(axis.title.x = element_text(size = 9)) +
  ylim(0,1)
```

``` r
plt2 + 
  stat_summary(geom = 'point',
               mapping = aes(x = year_f, y = as.numeric(ChlCCC)),
               data = exceeds, fun = function(x) 1 - mean(x, na.rm = TRUE),
               size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'point', x = xanchor, y = .3,
           size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'text', x = xanchor + space, y = .3, label = 'Observed', hjust = 0,
             size =3.5) +

    annotate(geom = 'point', x = xanchor, y = .25,
           size = 2, color = cbep_colors()[1]) +
    annotate(geom = 'text', x = xanchor + space, y = .25, label = 'Adjusted', hjust = 0,
             size = 3.5)
#> Warning: Removed 1067 rows containing non-finite values (stat_summary).
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

## By Month

Only Model 5 is fit with the `month_f` term. When we included other
predictors, along with thE month, we ran into problems fitting the
model.

``` r
my_ref_grid <- ref_grid(ccc_gamm_1,  cov.reduce = median) 
a <- summary(emmeans(my_ref_grid, ~ month_f, type = 'response'))

my_ref_grid <- ref_grid(ccc_gamm_5,  cov.reduce = median) 
b <- summary(emmeans(my_ref_grid, ~ month_f, type = 'response'))

observed <- exceeds %>%
  select(month_f, ChlCCC) %>%
  group_by(month_f) %>%
  summarize(month_f = first(month_f),
              n = sum(! is.na(ChlCCC)), 
              prob = mean(ChlCCC, na.rm = TRUE),
              SE = sqrt((prob * (1-prob))/n),
              lower.CL = prob - 1.96 * SE,
              upper.CL = prob + 1.96 * SE,
              .groups = 'drop')
```

``` r
 z <- tibble(Month = factor(month.abb[3:11], levels = month.abb),
             observed = observed$prob, 
             mod_1 = a$prob, 
             mod_5 = b$prob )
```

``` r
z %>%
  pivot_longer(-c(Month, observed), names_to = 'Model',
               values_to = 'Probability') %>%
  arrange(Model, Month) %>%

  ggplot(aes(observed, Probability, color = Model)) +

  geom_segment(aes(xend=c(tail(observed, n= -1), NA),
                   yend=c(Probability[2:9], NA, Probability[11:18], NA)),
                   arrow=arrow(length=unit(0.075,"inches"), type = 'closed')) +
  geom_text(aes(x = observed,
                y = Probability,
                label = Month)) +

  geom_abline(slope = 1, intercept = 0, lty = 3) +
  xlab("Observed Value") +
  ylab('Model Predictions (Adjusted)') +
  xlim(0,.6) +
  ylim(0,.6)
#> Warning: Removed 2 rows containing missing values (geom_segment).
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />
There is a lot of scatter here. The pattern is revealing. Our “flow
adjusted” monthly values are lower in winter and higher in summer than
the observed values, while the simpler model without a flow term fits
the observed values more closely.

This reflects the influence of flow on the marginal means. The usual low
flow in summer means “observed” levels generally co-occur with low flow,
when we now know chloride levels tend to be slightly higher, and thus
probability of meeting the chloride standards slightly lower. As the
algorithm adjusts the estimated probabilities to median flow (which is
rare in the summer), it gives predictions for summer flows at median
stream flow – which is a rare event. The predictions may be accurate,
but they are not accurate regarding circumstances we actually observe in
the field . Thus the “adjusted” marginal means are slightly better than
what we actually observed.

``` r
rm(z)
```

### Graphic

Based on simplified model, model 5.

``` r
s <- b %>% 
  mutate(fprob = 1-prob,
         fUCL = 1 - lower.CL,
         fLCL = 1 - upper.CL)

plt1 <- ggplot(s, aes(month_f, fprob)) +
 
  geom_pointrange(aes(ymin = fLCL, ymax = fUCL),
                color = cbep_colors()[1], size = .75) +
  
  ylab('Probability of Exceeding\nChronic Chloride Standard') +

  theme_cbep(base_size = 12) +
  theme(axis.title.x = element_text(size = 9)) +
  ylim(0,1)
```

``` r
xanchor = 1
space = 0.3

plt1 + 
   stat_summary(geom = 'point',
                mapping = aes(x = as.numeric(month_f)-2, y = as.numeric(ChlCCC)),
                data = exceeds, fun = function(x) 1 - mean(x, na.rm = TRUE),
                size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'point', x = xanchor, y = .3,
           size = 2, shape = 23, fill = cbep_colors()[4]) +
    annotate(geom = 'text', x = xanchor + space, y = .3, 
             label = 'Observed', hjust = 0,
             size = 3.5) +

    annotate(geom = 'point', x = xanchor, y = .25,
           size = 2, color = cbep_colors()[1]) +
    annotate(geom = 'text', x = xanchor + space, y = .25, 
             label = 'Adjusted', hjust = 0,
             size = 3.5)
#> Warning: Removed 1067 rows containing non-finite values (stat_summary).
```

<img src="Chloride_Frequencies_Summary_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />
