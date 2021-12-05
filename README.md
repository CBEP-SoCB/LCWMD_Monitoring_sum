# LCWMD_Monitoring

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />


Analysis of ten years of water quality monitoring data from the Long Creek
Watershed Management District in Maine.

# Introduction
This archive includes R scripts and related products analyzing approximately ten
years of water quality monitoring data from the Long Creek Watershed Management
District (LCWMD).The District is a Maine quasi-municipal non-profit corporation
that manages water quality in the Long Creek watershed on behalf of about 130
land owners.

The Long Creek Watershed, almost three and a half square miles in area, is
dominated by commercial land use. The Maine Mall is one of the largest
landowners in the watershed, and it is surrounded by a range of commercial
businesses, from medical offices, to car washes.  About a third of the watershed
in impervious surfaces like roads, parking lots, and rooftops.

Landowners with an acre or more of impervious area are required to get a Clean
Water Act permit for stormwater discharges from their property.  The LCWMD
provides an alternative for landowners who otherwise would have to get an
individual permit. Landowners who elect to participate in the The Long Creek
Watershed Management District receive a General Permit, in return for providing
funding to the District, and facilitating the work of the district by permitting
access to their property for certain activities.

For more information on LCWMD, see [their web site](https://restorelongcreek.org).

Over the past decade, LCWMD has contracted with several consulting firms to
provide  water quality monitoring services along Long Creek.  This has produced
one of the most extensive and best documented data set from the Northeastern US 
looking at water quality conditions in an urban stream.

GZA Geoenvironmental Incorporated (GZA) has been the primary monitoring
contractor for LCWMD for several years, and in 2019, they conducted a thorough
review of LCWMD data. These analyses are based on their summary data sets, and
recapitulate and extend their analyses.

# Statement of Purpose
CBEP is committed to the ideal of open science.  Our State of the Bay data
archives ensure the science underlying the 2020 State of the Bay report is
documented and reproducible by others. The purpose of these archives is to
release  data and data analysis code whenever possible to allow others to
review, critique, learn from, and build upon CBEP science.

# Archive Structure
CBEP 2020/2021 State of the Bay data analysis summaries contain a selection of 
data,  data analysis code, and visualization code as used to produce 
results shared via our most recent State of Casco Bay report. Usually, these
archives are organized into two or three folders, including the following:

- `Data`  folder.  Contains data in simplified or derived form as used in our
data  analysis.  Associated metadata is contained in related Markdown documents,
usually `DATA_SOURCES.md` and `DATA_NOTES.md`.

- Analysis.  Contains one or more R Notebooks proceeding through the principal
data analysis steps that underpin SoCB reporting. To simplify the archives,
much preliminary analysis, and many analysis "dead ends" have been omitted. 

- Graphics.  Contains R Notebooks stepping through development of graphics, and
also copies of resulting graphics, usually in \*.png and \*.pdf formats.  These
graphics may differ from graphics as they appear in final State of the Bay
graphical layouts. Again, most draft versions of graphics have been omitted for 
clarity.

# Are Water Quality Criteria Met?
In some notebooks included in this archive, we ask  whether water quality 
criteria have been  met.

In this Archive a "TRUE" value consistently implies that water quality criteria
were met or exceeded, whether that is achieved by a value higher than or lower
than some numeric criteria.  "TRUE" implies good conditions.  "FALSE" implies 
bad conditions.
    
## Sources of Threshold Values  
### Dissolved oxygen
Maine’s Class B water quality standards call for dissolved oxygen above 7 mg/l,
with percent saturation above 75%. The Class C Standards, which apply to almost
all of Long Creek, call for dissolved oxygen above 5 mg/l, with percent
saturation above 60%. In addition, for class C conditions, the thirty day
average dissolved oxygen must stay above 6.5 mg/l.

### Chloride
Maine uses established thresholds for both chronic and acute exposure to
chloride. These are the “CCC and CMC” standards for chloride in freshwater.
(06-096 CMR 584). These terms are defined in a footnote as follows:

>   The Criteria Maximum Concentration (CMC) is an estimate of the highest
    concentration of a material in surface water to which an aquatic community
    can be exposed briefly without resulting in an unacceptable effect. The
    Criterion Continuous Concentration (CCC) is an estimate of the highest
    concentration of a material in surface water to which an aquatic community
    can be exposed indefinitely without resulting in an unacceptable effect.

The relevant thresholds are:

*   Chloride CCC  = 230  mg/l
*   Chloride CMC  = 860  mg/l

In practice, chloride levels in Long Creek are indirectly estimated based on 
measurement of conductivity.  The chloride-conductivity correlations are fairly
close and robust, but estimation is an additional source of error, although 
generally on the level of 10% or less.

### Temperature
There are no legally binding Maine criteria for maximum stream temperature, but
we can back into thresholds based on research on thermal tolerance of brook
trout in streams. A study from Michigan and Wisconsin, showed that trout are
found in streams with daily mean water temperatures as high as 25.3°C, but only
if the period of exceedance of that daily average temperature is short – only
one day. Similarly, the one day daily maximum temperature above which trout were
never found was 27.6°C. That generates two temperature criteria, one for daily
averages, and one for daily maximums. 

These criteria should be taken as rough values only, as the  original study was
observational, and thus the key driver of suitability for trout could be
another stressor correlated with these temperature metrics.

>  Wehrly, Kevin E.; Wang, Lizhu; Mitro, Matthew (2007). “Field‐Based Estimates
   of Thermal Tolerance Limits for Trout: Incorporating Exposure Time and
   Temperature Fluctuation.” Transactions of the American Fisheries Society
   136(2):365-374.


