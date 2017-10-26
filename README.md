# SLATRRA <br>
Satellite Latency Assessment Tools for Real-time River Applications <br>
<br>
Contains the codes used in analysis and production of figures in "Global estimates of river flow wave travel times and implications for low-latency satellite data" by Allen et al. <br>
<br><br>
Code Descriptions: <br>
<br>
<b>waveCelerityCalculator.R</b><br>
Contains the data analysis and production of figures and tables in the Allen et al., "Global estimates of river flow wave travel times and implications for low-latency satellite data". The most basic sections of the code include: <br>
•	Calculation of river slope <br>
•	Joining points of interest (POI: cities, dams, and gauges) to the river network <br>
•	Modeling flow wave celerity and calculating travel time <br>
•	Validating the model with a gauge-based estimate of celerity <br>
•	Generating figures and tables, and calculating the various statistics presented in the paper <br>
<br>
<b>SWOTswathGenerator.R</b><br>
Generates the SWOT dual swath paths from Aviso ephemeris table. Data output is a polygon shapefile with attributes corresponding to track orbit, swath side, direction, satellite location, and time along the full orbit cycle. See the associated dataset descriptions section above for details of output.

