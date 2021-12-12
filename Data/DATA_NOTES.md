# Data Notes
==============

<img
  src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
  style="position:absolute;top:10px;right:50px;" />
  
## Data Files in This Folder

1.  **Sonde_Data.csv**

A large CSV file, containing data collected principally by sondes.  Data was 
collected every 15 minutes or 1 hour, depending on the deployment.

Column Name | Contents                   | Units / Values
------------|----------------------------|--------------
Unnamed     | Arbitrary sequence numbers |  
DT          | Date and Time              | '%m/%d/%Y %H:%M'
Chl	        | Chloride, Calculated       | mg/l
D	          | Depth                      | meters
DO          | Dissolved Oxygen           | mg/l
PctSat      | Oxygen percent saturation  | Percent
pH          | pH                         | NBS pH Scale (Sonde-based)
Press	      | Water Pressure             | kPa
SpCond	    | Specific Conductance       | uS/cm
T	          | Temperature                | Degrees C
Precip	    | Precipitation  (Jetport)   | inches
Site        | Monitoring Location Code   |
------------|----------------------------|--------------

Precipitation data is derived from daily published rainfall totals for the 
Portland jetport, and are daily total values.

Chloride values were calculated from specific conductance based on a multi-year
specific conductance-chloride linear regression. R^2^ values are on the order of
0.9.  Detailed statistical analysis shows that the relationship is not strictly
linear, but deviations from linearity have little practical import. The presence
of outliers also suggests some care in interpreting individual extreme values.

2.  **Site_IC_Data.csv** --  Simplified data of
    direct and cumulative subwatershed area and imperviousness for each Long
    Creek subwatershed / Monitoring Location.  Entered by hand from a table in a
    report by GZA to LCWMD, with cumulative watershed imperviousness added.

Column Name | Contents                   | Units / Values  
------------|----------------------------|--------------  
Site        | Monitoring Location Code   |    
Subwatershed | Text name of location     | 
Area_ac	    | Area of subwatershed       | Acres
IC_ac	      | Impervious area in subwatershed | Acres
CumArea_ac	| Cumulative watershed area including areas upstream | Acres
CumIC_ac	  | Cumulative watershed imperviousness | Acres
PctIC	      | Percent of impervious cover this subwatershed | Percent
CumPctIC    | Cumulative percent imperviousness  | Percent
------------|----------------------------|--------------  

3.  **Weather Data.csv**  --  Selected Weather data (including average wind
    speeds, precipitation, snow, and temperature), along with lagged and
    weighted precipitation data (mostly lagged and weighted log of
    precipitation).  Data was downloaded from 
    
    CBEP uses a custom Python program to download data from NOAA's online data
    APIs.  Specifically, data were accessed through NOAA's National 
    Centers for Environmental Information. We downloaded daily (GHCND) weather 
    summaries via API v2. Information on this API is available here:
    https://www.ncdc.noaa.gov/cdo-web/webservices/v2

    Documentation on contents of specific datasets is available at
    https://www.ncdc.noaa.gov/cdo-web/datasets


Column Name | Contents                   | Units / Values  
------------|----------------------------|--------------  
Unnamed     | Arbitrary sequence numbers |  
sdate       | Date                       | '%m/%d/%Y'  
Year        | Year, as an integer        |  
Month       | As an integer 1 = January  |
DOY         | Julian Day                 | 1:366
Precip      | Daily precipitation        | mm  
MaxT        | Daily maximum air temperature | Degrees C
lPrecip     | natural log of Precip      |
wPrecip	    | Weighted value of prior 10 days precipitation  | mm
wlPrecip    | Weighted values of log of prior 10 days precipitation |  
------------|----------------------------|--------------  

The NOAA API provides rainfall in tenths of millimeters.  We converted to 
millimeters in this file.

Weighted precipitation values use exponential weighting, with a rate parameter
of $0.8$.  Thus yesterday's precipitation counts as $80%$ of today's rainfall,
while the day before counts as $0.8 /times 0.8 = 0.64$ times as much as today's
rainfall, and so on.  Weighted rainfall *does not include* the present-day's
precipitation.

4.  **Daily_Data.csv**  
A large CSV file, containing daily summaries of the data contained in 
`Sonde_Data.csv`.  The summaries include: minimum (min), maximum (max), median, 
mean, standard  deviation ('sd'), interquartile range ('IQr') and sample size 
('n') for all numeric values (except precipitation)  found in `sonde_data.csv`. 
The order of columns is 

Column Name | Contents                   | Units / Values  
------------|----------------------------|--------------  
Site        | Monitoring Location Code   |  
sdate       | Date                       | '%m/%d/%Y'  
Values      | Summaries as {parameter}_{function} | 
Year        | Year, as an integer        |  
Month       | As an integer 1 = January  |  
Precip      | Daily precipitation        | tenths of mm    
pPrecip     | Previous day precipitation | tenths of mm
MaxT        | Daily maximum air temp.    | tenths of degree C
------------|----------------------------|--------------  

We replaced the precipitation data from the `sonde_data.csv` file with
newly downloaded weather data. We added daily maximum temperature data from the 
same source.  Both are in tenths of their customary units.

5.  **Exceeds_Data.csv**  -- Data derived from *Daily_Data.csv* containing flags
    that indicate whether conditions that day exceeded acute or chronic exposure
    thresholds, even once. Thresholds for Chlorides, DO, and Percent Saturation
    are derived from Maine water quality criteria.  Temperature thresholds are
    derived from a study of brook trout habitat use in the upper midwest. 
    
    **Dissolved Oxygen**  
    Maine's Class B standards call for dissolved oxygen above 7 mg/l, with
    percent saturation above 75%.  The Class C Standards, which apply to almost
    all of Long Creek, call for dissolved oxygen above 5 mg/l, with percent
    saturation above 6.5 mg/l. In addition, the thirty day average dissolved
    oxygen must stay above 6.5 mg/l.

    **Chloride**  
    Maine uses established thresholds for both chronic and acute exposure to
    chloride. These are the "CCC and CMC" standards for chloride in freshwater.
    (06-096 CMR 584). These terms are defined in a footnote as follows:
    
    > The Criteria Maximum Concentration (CMC) is an estimate of the highest
    concentration of a material in surface water to which an aquatic community
    can be exposed briefly without resulting in an unacceptable effect. The
    Criterion Continuous Concentration (CCC) is an estimate of the highest
    concentration of a material in surface water to which an aquatic community
    can be exposed indefinitely without resulting in an unacceptable effect.

    **Temperature**  
    There are no criteria for maximum stream temperature, but we can back into
    thresholds based on research on thermal tolerance of brook trout in streams.
    A study from Michigan and Wisconsin, showed that trout are found in streams
    with daily mean water temperatures as high as 25.3°C, but only if the period
    of exceedance of that daily average temperature is short -- only one day.
    Similarly, the one day daily maximum temperature above which trout were not
    found was 27.6°C.
    
    >  Wehrly, Kevin E.; Wang, Lizhu; Mitro, Matthew (2007). "Field‐Based
    Estimates of Thermal Tolerance Limits for Trout: Incorporating Exposure Time
    and Temperature Fluctuation." Transactions of the American Fisheries Society
    136(2):365-374.

Column Name | Contents                   | Units / Values
------------|----------------------------|--------------
Unnamed     | Arbitrary sequence numbers |  
sdate       | Date                       | '%m/%d/%Y'  
Site        | Monitoring Location Code   |  
Year        | Year, as an integer        |  
Month       | As an integer 1 = January  |  
Precip      | Daily precipitation        | tenths of mm    
pPrecip     | Previous day precipitation | tenths of mm
MaxT        | Daily maximum air temp.    | degree C
D_Median	  | Median Daily depth         | meters
ClassCDO    | Pass Class C DO standard?  | TRUE/FALSE
ClassBDO    | Pass Class B DO standard?  | TRUE/FALSE	
ClassC_PctSat	| Pass Class C Saturation standard?  | TRUE/FALSE
ClassB_PctSat | Pass Class B Saturation standard?  | TRUE/FALSE	
ClassCBoth	| Pass Class C DO and saturation standards?  | TRUE/FALSE
ClassBBoth	| Pass Class b DO and saturation standards?  | TRUE/FALSE
ChlCCC	    | Better than Chloride chronic exposure level  | TRUE/FALSE 
ChlCMC	    | Better than Chloride chronic exposure level  | TRUE/FALSE 
MaxT_ex	    | Max temperatures below  short-term exposure limit |  TRUE/FALSE 
AvgT_ex     | Average temperatures below chronic exposure limit. | TRUE/FALSE
------------|----------------------------|--------------

Note that the depth data (`D_Median`) is unique to each monitoring station, as 
it depends on local stream morphology and flow dynamics.


6.  **Full_Data.csv** -- Data derived from *Daily_data.csv* containing
    weighted sums for time series analysis.  The data contains missing values
    where there are gaps in the data, so that autocorrelation-based regression
    models rely only on data where lag data are available.

Column Name | Contents                   | Units / Values
------------|----------------------------|--------------
Unnamed     | Arbitrary sequence numbers |  
Site        | Monitoring Location Code   |   
sdate       | Date                       | '%m/%d/%Y'  
Year        | Year, as an integer        |  
Month       | As an integer 1 = January  |
DOY         | Julian Day                 | 1:366
Precip      | Daily precipitation        | mm  
MaxT        | Daily maximum air temperature | Degrees C
lPrecip     | natural log of Precip      |
wPrecip	    | Weighted value of prior 10 days precipitation  | tenths of mm
wlPrecip    | Weighted values of log of prior 10 days precipitation |  
Chl_Median  | Median daily chloride      | mg/l
D_Median	  | Median Daily depth         | meters
lD_Median	  | Log of D_median            |  
DO_Median	  | Median daily dissolved oxygen | mg/l
PctSat_Median | Median daily saturation  | percent
pH_Median	  | Median daily pH            |  
spCond_Median | median daily specific conductance | uS/cm 
T_Median    | Median daily water temperature | Degrees C
 
7.  **grab_data** a folder containing several csv files.  Files contain subsets
    by subject area of data derived from the original excel data sheets. Data
    includes grab sample data only (no data logger or other "continuous" data).

    Data in this folder are divided into several files, to keep size of each
    file more reasonable.  These data were not analyzed in detail as part of 
    the State of Casco Bay report. Data may be incomplete, or sampling history
    may be uneven.  We include the data here for archival purposes.
    
    Each file has the  same sample identifier information, as follows:
    
Column Name | Contents                   | Units / Values
------------|----------------------------|--------------
site        |  Monitoring Location Code  |   
thedate     | Date                       | '%m/%d/%Y'  
type        | Type of Sample ('N' or Field Duplicate ('FD')) |  
taskcode    | Category                   | Baseflow, Stormwater, Spring Melt
------------|----------------------------|--------------
 
    a.  **Organics (Polycyclic Aromatic Hydrocarbons and PCBs)**
    It appears the only PCB ever measured was dibenzofuran.  All the rest of 
    these are PAHs.  Units are ug/l.
        *  2-Chloronaphthalene
        *  2-Methylnaphthalene
        *  Acenaphthene
        *  Acenaphthylene
        *  Anthracene
        *  Benzo(a)anthracene
        *  Benzo(a)pyrene
        *  Benzo(b)fluoranthene
        *  Benzo(ghi)perylene
        *  Benzo(k)fluoranthene
        *  Chrysene
        *  Dibenz(a,h)anthracene 
        *  Dibenzofuran 
        *  Fluoranthene 
        *  Fluorene   
        *  Indeno(1,2,3-cd)pyrene 
        *  Naphthalene 
        *  Phenanthrene 
        *  Pyrene 
        
    b. **Metals**
    Units are mg/l
        *  Aluminum
        *  Antimony
        *  Arsenic
        *  Barium
        *  Beryllium
        *  Cadmium
        *  Chromium 
        *  Cobalt
        *  Copper 
        *  Iron
        *  Lead
        *  Manganese
        *  Mercury 
        *  Nickel                                      
        *  Selenium
        *  Silver
        *  Thallium
        *  Vanadium                     
        *  Zinc
    
    c. **Major Cations, Anions, Salinity, Etc.**
    Units are mg/l except where specified otherwise.
        *  Calcium 
        *  Magnesium
        *  Potassium
        *  Sodium
        *  CHLORIDE (AS Cl)
        *  Chloride, Calculated
        *  Hardness (As CaCO3) (used for determining toxicity of metals)
        *  Salinity (from sodium) (PPT)
        *  Specific Conductivity  (uS/cm)
        *  Total Dissolved Solids (g/l)
    
    d. **Macronutrients**
    Units are mg/l
        *  Ammonia', 
        *  Nitrate as N', 
        *  Nitrite as N', 
        *  Nitrogen',
        *  Nitrogen, Kjeldahl, Total'
        *  Organic Nitrogen', 
        *  Phosphate Ion', 
        *  Phosphorus'
    
    e. **Other WQ Parameters**
        *  Cyanide
        *  Oxidation-Reduction Potential
        *  pH
        *  Turbidity
        *  Dissolved Oxygen
        *  Dissolved Oxygen Saturation
    
    

