# MarsOrbitPredictor
Predicting the orbit of Mars, given oppositions data

**Task** 
The data in the spreadsheet (downloadable by clicking the down-arrow button,  04_cricket_1999to2011.csv) contains data on ODI matches from 1999 to 2011.
1. Using the first innings data alone in the above spreadsheet, find the best fit 'run production functions' in terms of wickets-in-hand w and overs-to-go u. Assume the model Z(u,w) = Z0(w)[1 - exp{-b(w)u}]. Use the sum of squared errors loss function, summed across overs and wickets. You will find one Z0(w) and b(w) for each w. In your report, you should provide a plot of the ten functions, report the (20) parameters associated with the (10) production functions, and the error per point. Your function should be named:
Z0, b = DuckworthLewis20Params() 

2. Now assume the model Z(u,w) = Z0(w)[1 - exp{-Lu/Z0(w)}] and use the sum of squared errors loss function, summed over overs and wickets. Note that this new regression forces all slopes to be equal at u = 0.  In your report, you should provide a plot of the ten functions, report the (11) parameters associated with the (10) production functions, and the error per point. This function should be named:
Z0, L = DuckworthLewis11Params() 

**Dataset: 01_data_mars_opposition_updated**
This file contains data on longitude/latitude of Mars under "opposition" with the Sun, in the ecliptic coordinate system.
1. Columns A/B/C/D/E are year/month/day/hour/minute of the opposition.
2. Columns F,G,H,I denote the ZodiacIndex, Degree, Minute, Second, respectively, of Mars's (heliocentric) longitude in the ecliptic coordinate system. 
3. ZodiacIndex refers to the zodiac (Aries 0, Taurus 1, ..., Pisces 11).
4. Longitude = ZodiacIndex*30 + Degree + Minute/60 + Second/3600 (degrees)
5. Columns J,K refer to degree, minute of the geocentric latitudinal position of Mars in the ecliptic coordinate system.
6. Columns L,M,N,O refer to Mars's mean longitude, with reference to Kepler's approximated equant. Instead of using this, you will find these based on your own equant.

**Report: **
Report.pdf

**Script:**
OrbitPrediction.py

**How to run code:**
python OrbitPrediction.py

Requirements:
python3.9.0
