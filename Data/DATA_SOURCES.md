# Data Sources
## Water Quality Data
All water quality data analyzed here is derived from data provided by GZA to the
Long Creek Watershed Management District (LCWMD).  GZA has been the primary
monitoring, data management, and data analysis contractor for the LCWMD for the
last several years. GZA manages a large water quality data database on behalf of
LCWMD. The data provided constitutes a "snapshot" of the database as of mid
2019.

Data analyzed here were derived from a ZIP file transmitted from GZA to LCWMD
and Curtis Bohlen of CBEP via on-line transfer. Bohlen accessed the files
on April 23 of 2020.

The original form of the data in excel was took up a lot of disk space, memory, 
and computation to access. It contained a high fraction of replicate data, and 
it used inefficient lookup code to create consolidated data sheets.

CBEP simplified and reorganized the data to produce the derived data tables
included here.Data QA/QC and data reorganization code are availalbe as part of
out detailed data analysis 
[archives](https://github.com/CBEP-SoCB-Details/LCWMD_Monitoring), in
the "Derived_Data" folder.

## Weather Data
CBEP used a custom Python program to download data from NOAA's online data
APIs.  Specifically, data were accessed through NOAA's National Centers
for Environmental Information.

Here, we have downloaded daily (GHCND) weather summaries via API v2. Information
on this API is available here: https://www.ncdc.noaa.gov/cdo-web/webservices/v2

Documentation on specific datasets is available at
https://www.ncdc.noaa.gov/cdo-web/datasets

Portland Jetport weather data was downloaded by Curtis C. Bohlen using a custom
python script, titled "noaaweatherdataGUI.py" on April 23, 2020.

### Units
Data is in SI units, except that NOAA provides some data in tenths of the
nominal units.  This is not well documented through the API, but obvious in 
context. Temperatures are reported in tenths of degrees C, and precipitation in
tenths of a millimeter.  For this analysis, we disregard trace rainfall
events.  Users shuld check the metadata for different derived data products to
determine whether data a salready been converted to more conventinal units.  

## Impervious Cover Data
Site IC Data was derived from a report by GZA to LCWMD. 

> GZA 2019.  Re: Long Creek Watershed Data Analysis Task 2: Preparation of
  Explanatory and Other Variables.  November 13, 2019. File No. 09.0025977.02,
  received as `25977.02 Task 2 Memo 111319.pdf`

Specifically, data was derived from Table 2: "Summary of Impervious Cover 
Analysis between 2009 and 2019"  Data was entered by hand into Microsoft Excel.

Each monitoring station is in a different subwatershed of the Long Creek 
Watershed (although not always at the downstream end of that subwatershed).  The
Table indicates total imperviousness in those subwatersheds. The geometry of
Long Creek is fairly simple, so it is simple to identify which subwatersheds 
drain to other subwatersheds. CBEP used that information to calculate (in Excel) 
the cumulative watershed area, including "upstream" subwatersheds, for each
monitoring station and the cumulative "upstream" impervious area. Data was 
exported to the CSV file to simplify data exchange.
