import math

print 'it works'

day = 22
month = 6
year = 2016
h = -3

lon = 30.289775
lat = 59.938406

def calcJD(year, month, day):
    if month <= 2:
        year -= 1
        month += 12
    A = int(year/100)
    B = 2 - A + int(A/4)

    JD = int(365.25*(year + 4716)) + int(30.6001*(month+1)) + day + B - 1524.5

    return JD

def calcDayOfWeek(juld):
    A = (juld + 1.5) % 7
    return A

def isLeapYear(yr):
    return ((yr % 4 == 0 and yr % 100 != 0) or yr % 400 == 0)
    
def calcDayOfYear(mn, dy, lpyr):
    k = 1 if lpyr else 2
    doy = int((275 * mn)/9) - k * int((mn + 9)/12) + dy -30
    return doy

def calcTimeJulianCent(jd):
    T = (jd - 2451545.0)/36525.0
    return T

def calcMeanObliquityOfEcliptic(t):
    seconds = 21.448 - t*(46.8150 + t*(0.00059 - t*(0.001813)))
    e0 = 23.0 + (26.0 + (seconds/60.0))/60.0
    return e0

def degToRad(angleDeg):
    return math.pi * angleDeg / 180.0

def calcObliquityCorrection(t):
    e0 = calcMeanObliquityOfEcliptic(t)

    omega = 125.04 - 1934.136 * t
    e = e0 + 0.00256 * int(degToRad(omega))
    return e

def calcGeomMeanLongSun(t):
    L0 = 280.46646 + t * (36000.76983 + 0.0003032 * t)
    while(L0 > 360.0):
        L0 -= 360.0
    while(L0 < 0.0):
        L0 += 360.0
    return L0

def calcGeomMeanAnomalySun(t):
    M = 357.52911 + t * (35999.05029 - 0.0001537 * t)
    return M

def calcSunEqOfCenter(t):
    m = calcGeomMeanAnomalySun(t)

    mrad = degToRad(m)
    sinm = math.sin(mrad)
    sin2m = math.sin(mrad+mrad)
    sin3m = math.sin(mrad+mrad+mrad)

    C = sinm * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin2m * (0.019993 - 0.000101 * t) + sin3m * 0.000289
    return C

def calcSunTrueLong(t):
    l0 = calcGeomMeanLongSun(t)
    c = calcSunEqOfCenter(t)

    O = l0 + c
    return O

def calcSunApparentLong(t):
    o = calcSunTrueLong(t)

    omega = 125.04 - 1934.136 * t
    lamb = o - 0.00569 - 0.00478 * math.sin(degToRad(omega))
    return lamb

def radToDeg(angleRad):
    return 180.0 * angleRad / math.pi

def calcSunRtAscension(t):
    e = calcObliquityCorrection(t)
    lamb = calcSunApparentLong(t)
 
    tananum = (math.cos(degToRad(e)) * math.sin(degToRad(lamb)))
    tanadenom = (Math.cos(degToRad(lamb)))
    alpha = radToDeg(math.atan2(tananum, tanadenom))
    return alpha

def calcSunDeclination(t):
        e = calcObliquityCorrection(t)
        lamb = calcSunApparentLong(t)

        sint = math.sin(degToRad(e)) * math.sin(degToRad(lamb))
        theta = radToDeg(math.asin(sint))
        return theta

def calcEccentricityEarthOrbit(t):
    e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t)
    return e

def calcEquationOfTime(t):
    epsilon = calcObliquityCorrection(t)
    l0 = calcGeomMeanLongSun(t)
    e = calcEccentricityEarthOrbit(t)
    m = calcGeomMeanAnomalySun(t)

    y = math.tan(degToRad(epsilon)/2.0)
    y *= y

    sin2l0 = math.sin(2.0 * degToRad(l0))
    sinm   = math.sin(degToRad(m))
    cos2l0 = math.cos(2.0 * degToRad(l0))
    sin4l0 = math.sin(4.0 * degToRad(l0))
    sin2m  = math.sin(2.0 * degToRad(m))

    Etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0 - 0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m

    return radToDeg(Etime)*4.0

def calcJDFromJulianCent(t):
    JD = t * 36525.0 + 2451545.0
    return JD

def calcSolNoonUTC(t, longitude):
    # First pass uses approximate solar noon to calculate eqtime
    tnoon = calcTimeJulianCent(calcJDFromJulianCent(t) + longitude/360.0)
    eqTime = calcEquationOfTime(tnoon)
    solNoonUTC = 720 + (longitude * 4) - eqTime

    newt = calcTimeJulianCent(calcJDFromJulianCent(t) -0.5 + solNoonUTC/1440.0) 

    eqTime = calcEquationOfTime(newt)

    solNoonUTC = 720 + (longitude * 4) - eqTime
        
    return solNoonUTC

def calcHourAngleSunrise(lat, solarDec):
    latRad = degToRad(lat)
    sdRad  = degToRad(solarDec)

    HAarg = (math.cos(degToRad(90.833))/(math.cos(latRad)*math.cos(sdRad))-math.tan(latRad) * math.tan(sdRad))

    HA = (math.acos(math.cos(degToRad(90.833))/(math.cos(latRad)*math.cos(sdRad))-math.tan(latRad) * math.tan(sdRad)))

    return HA


def calcSunriseUTC(JD, latitude, longitude):
    t = calcTimeJulianCent(JD);

    # *** Find the time of solar noon at the location, and use
    #     that declination. This is better than start of the 
    #     Julian day

    noonmin = calcSolNoonUTC(t, longitude)
    tnoon = calcTimeJulianCent (JD+noonmin/1440.0)

    # *** First pass to approximate sunrise (using solar noon)

    eqTime = calcEquationOfTime(tnoon)
    solarDec = calcSunDeclination(tnoon)
    hourAngle = calcHourAngleSunrise(latitude, solarDec)

    delta = longitude - radToDeg(hourAngle)
    timeDiff = 4 * delta   # in minutes of time
    timeUTC = 720 + timeDiff - eqTime  # in minutes

    # *** Second pass includes fractional jday in gamma calc

    newt = calcTimeJulianCent(calcJDFromJulianCent(t) + timeUTC/1440.0) 
    eqTime = calcEquationOfTime(newt)
    solarDec = calcSunDeclination(newt)
    hourAngle = calcHourAngleSunrise(latitude, solarDec)
    delta = longitude - radToDeg(hourAngle)
    timeDiff = 4 * delta;
    timeUTC = 720 + timeDiff - eqTime; # in minutes

    return timeUTC

def calcHourAngleSunset(lat, solarDec):
    latRad = degToRad(lat);
    sdRad  = degToRad(solarDec)

    HAarg = (math.cos(degToRad(90.833))/(math.cos(latRad)*math.cos(sdRad))-math.tan(latRad) * math.tan(sdRad))

    HA = (math.acos(math.cos(degToRad(90.833))/(math.cos(latRad)*math.cos(sdRad))-math.tan(latRad) * math.tan(sdRad)))

    return -HA;     # in radians

def calcSunsetUTC(JD, latitude, longitude):
    t = calcTimeJulianCent(JD)

    # *** Find the time of solar noon at the location, and use
    #     that declination. This is better than start of the 
    #     Julian day

    noonmin = calcSolNoonUTC(t, longitude)
    tnoon = calcTimeJulianCent (JD+noonmin/1440.0)

    # First calculates sunrise and approx length of day

    eqTime = calcEquationOfTime(tnoon)
    solarDec = calcSunDeclination(tnoon)
    hourAngle = calcHourAngleSunset(latitude, solarDec)

    delta = longitude - radToDeg(hourAngle)
    timeDiff = 4 * delta
    timeUTC = 720 + timeDiff - eqTime

    #/ first pass used to include fractional day in gamma calc

    newt = calcTimeJulianCent(calcJDFromJulianCent(t) + timeUTC/1440.0)
    eqTime = calcEquationOfTime(newt)
    solarDec = calcSunDeclination(newt)
    hourAngle = calcHourAngleSunset(latitude, solarDec)

    delta = longitude - radToDeg(hourAngle)
    timeDiff = 4 * delta
    timeUTC = 720 + timeDiff - eqTime # in minutes

    return timeUTC

def calcDayFromJD(jd):
    z = int(jd + 0.5)
    f = (jd + 0.5) - z

    if (z < 2299161): 
        A = z
    else:
        alpha = int((z - 1867216.25)/36524.25)
        A = z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int((B - 122.1)/365.25)
    D = int(365.25 * C)
    E = int((B - D)/30.6001)

    pday = B - D - int(30.6001 * E) + f
    pmonth = E - 1 if (E < 14) else E - 13
    pyear = C - 4716 if (month > 2) else C - 4715

    return (("0" if day<10 else "") + str(pday) + str(pmonth))




def timeStringShortAMPM(minutes, JD):
    julianday = JD
    floatHour = minutes / 60.0
    hour = int(floatHour)
    floatMinute = 60.0 * (floatHour - int(floatHour))
    minute = int(floatMinute)
    floatSec = 60.0 * (floatMinute - int(floatMinute))
    second = int(floatSec + 0.5)
    PM = False

    minute += 1 if (second >= 30) else 0

    if (minute >= 60):
        minute -= 60
        hour += hour

    daychange = False
    if (hour > 23):
        hour -= 24
        daychange = True
        julianday += 1.0

    if (hour < 0):
        hour += 24
        daychange = True
        julianday -= 1.0

    if (hour > 12):
        hour -= 12
        PM = True

    if (hour == 12):
        PM = True
    
    if (hour == 0):
        PM = False
        hour = 12

    timeStr = str(hour) + ":"
    if (minute < 10):    #  i.e. only one digit
        timeStr += "0" + str(minute) + ("PM" if (PM) else "AM")
    else:
        timeStr += "" + str(minute) + ("PM" if (PM) else "AM")

    if (daychange):
        return timeStr + " " + str(calcDayFromJD(julianday))
    return timeStr

def timeStringDate(minutes, JD):
    julianday = JD
    floatHour = minutes / 60.0
    hour = math.floor(floatHour)
    floatMinute = 60.0 * (floatHour - int(floatHour))
    minute = int(floatMinute)
    floatSec = 60.0 * (floatMinute - math.floor(floatMinute))
    second = int(floatSec + 0.5)

    minute += 1 if (second >= 30) else 0

    if (minute >= 60):
        minute -= 60
        hour += 1

    daychange = False
    if (hour > 23):
        hour -= 24
        julianday += 1.0
        daychange = True

    if (hour < 0):
        hour += 24
        julianday -= 1.0
        daychange = True

    timeStr = str(hour) + ":"
    if (minute < 10):    # i.e. only one digit
        timeStr += "0" + str(minute)
    else:
        timeStr += str(minute)

    if (daychange):
        return timeStr + " " + str(calcDayFromJD(julianday))
    return timeStr

def calcSun(lat, lon, year, day, month, tz):
        latitude = lat
        longitude = lon
        indexRS = month-1
        if (latitude >= -90) and (latitude < -89):
            latitude = -89
        if (latitude <= 90) and (latitude > 89):
            latitude = 89
            
#***** Calculate the time of sunrise           

#*********************************************************************/
#****************   NEW STUFF   ******   January, 2001   ****************
#*********************************************************************/

        JD = calcJD(float(year), month, float(day))
        doy = calcDayOfYear(month, float(day), isLeapYear(year))
        T = calcTimeJulianCent(JD)
        
        #alpha = calcSunRtAscension(T);
        theta = calcSunDeclination(T)
        Etime = calcEquationOfTime(T)

#*********************************************************************/

        eqTime = Etime
        solarDec = theta

        #Calculate sunrise for this date
        #if no sunrise is found, set flag nosunrise

        nosunrise = False;

        riseTimeGMT = calcSunriseUTC(JD, latitude, longitude)
        try:
            nosunrise = int(riseTimeGMT)
            nosunrise = False
        except:
            nosunrise = True

        # Calculate sunset for this date
        # if no sunset is found, set flag nosunset

        nosunset = False
        setTimeGMT = calcSunsetUTC(JD, latitude, longitude)
        try:
            nosunset = int(setTimeGMT)
            nosunset = False
        except:
            nosunset = True

        daySavings = 0
        zone = tz

        if (not nosunrise):   
            riseTimeLST = riseTimeGMT - (60 * zone) + daySavings; #  in minutes
            riseStr = timeStringShortAMPM(riseTimeLST, JD);
            utcRiseStr = timeStringDate(riseTimeGMT, JD);

            print "sunrise " +  riseStr;
            print "utcsunrise " +  utcRiseStr;

        if (not nosunset):   
            setTimeLST = setTimeGMT - (60 * zone) + daySavings
            setStr = timeStringShortAMPM(setTimeLST, JD)
            utcSetStr = timeStringDate(setTimeGMT, JD)

            print "sunset " +  setStr;
            print "utcsunset " +  utcSetStr;

        


def jd_to_cal(jd):
    j1=0
    j2=0
    j3=0 
    j4=0
    j5=0

    intgr = int(jd)
    frac  = jd - intgr
    gregjd  = 2299161;
    if intgr >= gregjd: 
        tmp = int( ( (intgr - 1867216) - 0.25 ) / 36524.25 )
        j1 = intgr + 1 + tmp - int(0.25*tmp)
    else:
        j1 = intgr

    #orrection for half day offset
    dayfrac = frac + 0.5
    if dayfrac >= 1.0:
        dayfrac -= 1.0
        j1 += 1

    j2 = j1 + 1524
    j3 = int( 6680.0 + ( (j2 - 2439870) - 122.1 )/365.25 )
    j4 = int(j3*365.25)
    j5 = int( (j2 - j4)/30.6001 )

    d = int(j2 - j4 - int(j5*30.6001))
    m = int(j5 - 1)
    if m > 12:
        m -= 12
    y = int(j3 - 4715)
    if m > 2:
        y -= 1
    if  y <= 0:
        y -= 1

    # get time of day from day fraction
    hr  = int(dayfrac * 24.0)
    mn  = int((dayfrac*24.0 - hr)*60.0)
    f   = ((dayfrac*24.0 - hr)*60.0 - mn)*60.0
    sc  = int(f)
    f -= sc;
    if f > 0.5:
        sc += 1
    if y < 0:
        y = -y
    year  = y
    month = m
    day   = d
    hour          = hr-h
    minute        = mn
    second        = sc
    return str(year)+' '+str(month)+' '+str(day)+' '+str(hour)+' '+str(minute)+' '+str(second)


#calculations

a = int((14-month)/12)
y = year + 4800 - a
m = month + 12*a - 3

#Julian day
JDN = day + int((153*m+2)/5)+365*y+int(y/4)-int(y/100)+int(y/400)-32045

print "JDN = "+str(JDN)

#JDN is the Julian date;
#n is the Julian day since Jan 1st, 2000 12:00.
#0.0008 is the fractional Julian Day for leap seconds and terrestrial time.
#currently = 68.184 / 86400 without DUT1.
n = JDN - 2451545.0+0.0008
print "n = "+str(n)

#Mean solar noon
Jstar = lon/360+n
print "Jstar = "+str(Jstar)

#Solar mean anomaky
M = (357.5291 + 0.98560028*Jstar) % 360
print "M = "+str(M)

#Equation of the center
C = 1.9148*math.sin(math.pi*M/180)+0.0200*math.sin(math.pi*2*M/180)+0.0003*math.sin(math.pi*3*M/180)
print "C = "+str(C)

#Ecliptic longitude (102.9372 is a value for the argument of perihelion)
G = (M + C + 180 + 102.9372) % 360
print "G = "+str(G)

#Solar transit
Jtr = Jstar + 0.0053*math.sin(math.pi*M/180)-0.0069*math.sin(math.pi*2*G/180)
print "Jtr = "+str(Jtr)

#Declination of the Sun
SinS = math.sin(math.pi*G/180)*math.sin(math.pi*23.44/180)
CosS = math.sqrt(1-SinS*SinS)
TanS = SinS/CosS
print "decl = "+str(math.asin(SinS)*180/math.pi)
print "SinS = "+str(SinS)
print "CosS = "+str(CosS)

#Hour angle
cosw = (math.sin(-0.83*math.pi/180)-math.sin(math.pi*lon/180)*SinS)/(math.cos(math.pi*lon/180)*CosS)
cosw1 = -1*math.tan(math.pi*lon/180)*TanS
cosw2 = (SinS*math.sin(math.pi*lon/180)+0.0145)/(CosS*math.cos(math.pi*lon/180))
print "cosw = "+str(cosw)
print "cosw1 = "+str(cosw1)
print "cosw2 = "+str(cosw2)
acosw = math.acos(cosw)*180/math.pi
print "acosw = "+str(acosw)

risetime = Jtr-acosw/360+2451545.0-0.0008
settime = Jtr+acosw/360+2451545.0-0.0008

print "NJset = "+str(settime)
print "NJrise = "+str(risetime)
print jd_to_cal(risetime)
print jd_to_cal(settime)
calcSun(lat, lon, year, day, month, h)
JD = calcJD(float(year), month, float(day))
sr = calcSunriseUTC(JD, lat, lon)
ss = calcSunsetUTC(JD, lat, lon)
print "Sunrise: " + str(int(sr//60)-h) + ":" + str(int(sr%60))
print "Sunset: " + str(int(ss//60)-h) + ":" + str(int(ss%60))