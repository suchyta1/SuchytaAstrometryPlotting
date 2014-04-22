#############################################################################
#
#	$Id: util.py
#	SISPI Utilities
#
#	This is a collection of utility functions
#           ArgParser     - parse command line arguments
#           float2ra      - Convert RA (dec) from degrees (float) to HH:MM:SS.SS format
#           float2dec
#           ra2float      - Convert RA (dec) from HH:MM:SS.SS format to degrees as a float
#           dec2float
#           radec2J2000   - Convert current values for RA, dec to J2000
#           J20002radec   - Convert RA, dec (J2000) to current values
#           slew_angle    - Calculate the angular distance between two pointings
#           moon_angle    - Returns moon position at CTIO as (ra, dec) in degrees (either current or at a given date)
#           moon_separation Returns angle in degrees between moon and telescope position
#           ctio_observer - Setup pyephem observer for CTIO
#           azel2radec    - Convert azimuth, elevation to ra, dec for a given (pyephem) observer
#           hadec2radec   - Convert hour angle, declination to ra, dec for a given (pyephem) observer
#           str2dict      - convert a list of space separated key=value pairs to a dictionary
#           bytes2human   - convert diskspace information to human readable information
#
#############################################################################

import math
import datetime
import ephem

#############################################################################
#
# A simple argument parser derived from the ArgParser in Pyro.util
#
# Options are given in the format -<option>  <value>
# <value> is optional
# Parameter:
#    args       - list of command line arguments (sys.argv[1:])
#
# Methods:
#    parse      - parses the args list
#    hasOption  - returns True is the option is given in args
#    getValue   - returns the value given for an option, raises KeyError is invalid option
#    getListofOptions   - returns a list of options found in args For convenience (ArgParser.options is public)
#
############################################################################


class ArgParser:
	def __init__(self):
		self.options = {}
		
	def parse(self, args):
		self.options = {}
			
		if type(args) == type(''):
			args = args.split()            # got a string, split it
	    
		while args:
			arg = args[0]
			del args[0]
			if arg[0]=='-':
				if args:
					value = args[0]
					if value[0]=='-':    # got another option
						self.options[arg[1:]] = None
					else:
						self.options[arg[1:]] = value
						del args[0]
				else:
					self.options[arg[1:]] = None
                        else:
                                self.options['REST'] = arg
                                break
                        
	def hasOption(self, option):
		return self.options.has_key(option)
	
	def getValue(self, option, default=Exception()):
		try:
			return self.options[option]
		except KeyError:
			if not isinstance(default,Exception):
				return default
			raise KeyError('No such option: %s' % option)
	def getListofOptions(self):
		return [opt for opt in self.options]

##############################################################################

def float2ra(ra):
    """
    Convert degrees as float into time string, HH:MM:SS.SS
    dec has to be of type float or a string like '34.222'
    360 deg = 24 hrs, 360/24 = 15
    """
    if ra is None:
        return ra
    if type(ra) is str or type(ra) is unicode:
        if ':' in str(ra):
            # if the string is not properly formatted this will throw an exception
            b=ra2float(ra)
            return ra
        float_ra = float(ra)
    else:
        float_ra = ra
    assert type(float_ra) is float,'Invalid parameter format'
    if float_ra < 0.0:
        sign = '-'
    else:
        sign = ''
    float_ra = abs(float_ra)
    hrs = float_ra/15. 
    hours = math.trunc(hrs)
    min = abs(hrs - hours) * 60.
    minutes = math.trunc(min) 
    seconds = round((min - minutes) * 60,3)
    if seconds == 60.0:
        seconds = 0.0
        minutes += 1
        if minutes == 60:
            minutes = 0
            hours += 1
            if hours == 24:
                hours = 0
    return sign+'%02i:%02i:%06.3f' % (hours, minutes, seconds)

def ra2float(ra):
    """
    Convert ra to degress (float).
    ra can be given as a time string, HH:MM:SS.SS or as string like '25.6554'
    or (trivially) as a float.
    An exception is thrown if ra is invalid
    360 deg = 24 hrs, 360/24 = 15
    """
    if ra is None:
        return ra
    if type(ra) is float or type(ra) is int:
        return float(ra)
    if (type(ra) is str or type(ra) is unicode) and ra.find(':') == -1:
        return float(ra)     
    assert type(ra) is str,'Invalid parameter format'
    h,m,s = ra.strip().split(':')
    if h.find('-') != -1:
        h=h.replace('-','')
        sign = -1.0
    else:
        sign = 1.0
    return sign*(float(h)*15.0 + float(m)/4.0 + float(s)/240.0)

def float2dec(dec):
    """
    Convert degrees as float into degree string, DD:MM:SS.SS
    dec has to be of type float or a string like '34.222'
    """
    if dec is None:
        return dec
    if type(dec) is str or type(dec) is unicode:
        if ':' in str(ra):
            # if the string is not properly formatted this will throw an exception
            b=dec2float(dec)
            return dec
        float_dec = float(dec)
    else:
        float_dec = dec
    assert type(float_dec) is float,'Invalid parameter format'    
    if float_dec < 0.0:
        sign = '-'
        float_dec = abs(float_dec)
    else:
        sign = ''
    degrees = math.trunc(float_dec)
    min = abs(float_dec - degrees) * 60
    minutes = math.trunc(min)
    seconds = round((min - minutes) * 60,3)
    if seconds == 60.0:
        seconds = 0.0
        minutes += 1
        if minutes == 60:
            minutes = 0
            degrees += 1
    return '%s%02i:%02i:%06.3f' % (sign, degrees, minutes, seconds)

def dec2float(dec):
    """
    Convert dec to degress (float).
    dec can be given as a time string, HH:MM:SS.SS or as string like '25.6554'
    or (trivially) as a float.
    An exception is thrown if dec is invalid
    """
    if dec is None:
        return dec
    if type(dec) is float or type(dec) is int:
        return float(dec)
    if (type(dec) is str or type(dec) is unicode) and dec.find(':') == -1:
        return float(dec)     
    assert type(dec) is str,'Invalid parameter format'
    d,m,s = dec.strip().split(':')
    if d.find('-') != -1:
        d=d.replace('-','')
        sign = -1.0
    else:
        sign = 1.0
    return sign*(float(d) + float(m)/60.0 + float(s)/3600.0)

def _supplement(fixed, date):
    T = (float(fixed) - 2000.0)/100.0
    t = (float(date) - float(fixed))/100.0
    asec = (2306.218 + 1.397*T)*t +1.095*t*t
    bsec = (2004.311 - 0.853*T)*t -0.427*t*t
    csec = (2306.218 + 1.397*T)*t + 0.302*t*t
    torad = 180.0/math.pi * 3600.0
    m = _makeMatrixSupplement( asec/torad, bsec/torad, csec/torad)
    return m

def _makeMatrixSupplement(a,b,c):
    m = [0.0 for x in range(9)]
    cA = math.cos(a)
    sA = math.sin(a)
    cB = math.cos(b)
    sB = math.sin(b)
    cC = math.cos(c)
    sC = math.sin(c)
    m[0] = cA * cB * cC - sA * sC
    m[3] = -cA * cB *sC - sA * cC
    m[6] = - cA * sB
    m[1] = sA * cB * cC + cA * sC
    m[4] = - sA * cB * sC + cA * cC
    m[7] = -sA * sB
    m[2] = sB * cC
    m[5] = -sB * sC
    m[8] = cB
    return m

def str2dict(dict_string):
    """
    Convert a string of key=value pairs to a dictionary.
    Format is 'KEY=value KEY=other value KEY=and some more values'
    For example: str2dict('key1 = 1 key2 = a b c key3=23.4) returns the dictionary
    {'key1':'1' , 'key2':'a b c', 'key3':'23.4'}
    """
    string_bits = dict_string.split('=')
    keys = [string_bits[0].strip()]
    values = []
    for bits in string_bits[1:-1]:
        pieces = bits.strip().rsplit(' ', 1)
        if len(pieces) == 1:
            key = pieces[0]
            value = 'NONE'
        else:
            key = pieces[1]
            value = pieces[0]
        keys.append(key)
        values.append(value)
    values.append(string_bits[-1])
    return dict(zip(keys, values))

def radec2J2000(RA,DEC):
    toDegrees = 180.0/math.pi
    # to radians
    ra = ra2float(RA)/180.0*math.pi
    dec = dec2float(DEC)/180.0*math.pi
    year = datetime.datetime.now().year
    m = _supplement(2000.0, year)
    results =  _Transform(ra, dec, m)
    # back to degrees
    results[0] = results[0] * 180.0/math.pi
    results[1] = results[1] * 180.0/math.pi
    return results

def J20002radec(RA,DEC):
    toDegrees = 180.0/math.pi
    # to radians
    ra = ra2float(RA)/180.0*math.pi
    dec = dec2float(DEC)/180.0*math.pi
    year = datetime.datetime.now().year
    m = _supplement(year, 2000.0)
    results =  _Transform(ra, dec, m)
    # back to degrees
    results[0] = results[0] * 180.0/math.pi
    results[1] = results[1] * 180.0/math.pi
    return results

def _Transform(ra, dec, m):
    r0 = [math.cos(ra)*math.cos(dec),
          math.sin(ra)*math.cos(dec),
          math.sin(dec)]
    s0 = [r0[0]*m[0] + r0[1]*m[1] + r0[2]*m[2],
          r0[0]*m[3] + r0[1]*m[4] + r0[2]*m[5],
          r0[0]*m[6] + r0[1]*m[7] + r0[2]*m[8]]
    r = math.sqrt(s0[0]*s0[0] + s0[1]*s0[1] + s0[2]*s0[2])
    results = [None, None]
    results[1] = math.asin(s0[2]/r)
    cosaa = ((s0[0]/r)/math.cos(results[1]))
    sinaa = ((s0[1]/r)/math.cos(results[1]))
    results[0] = math.atan2(sinaa, cosaa)
    if results[0] < 0:
        results[0] += 2* math.pi
    return results

def bytes2human(n):
    # http://code.activestate.com/recipes/578019
    # >>> bytes2human(10000)
    # '9.8K'
    # >>> bytes2human(100001221)
    # '95.4M'
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i+1)*10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n

# moon_angle: returns moon position at CTIO as (ra, dec) in degrees (either current or at a given date)
# moon_separation: returns angle in degrees between moon and telescope position
def ctio_observer(date=None, lat = None, lon = None, elevation = None):
        observatory = {'TELESCOP':'CTIO 4.0-m telescope',
                       'OBSERVAT':'CTIO',
                       'OBS-LAT':-30.16606,
                       'OBS-LONG':-70.81489,
                       'OBS-ELEV':2215.0}
        ctio = ephem.Observer()
        if lon == None:
            ctio.lon = observatory['OBS-LONG']*ephem.degree
        else:
            ctio.lon = lon
        if lat == None:
            ctio.lat = observatory['OBS-LAT']*ephem.degree
        else:
            ctio.lat = lat
        if elevation == None:
            ctio.elevation = observatory['OBS-ELEV']*ephem.degree
        else:
            ctio.elevation = elevation*ephem.degree
        if date == None:
            ctio.date = datetime.datetime.utcnow()
        else:
            ctio.date = date
        ctio.epoch=ephem.J2000
        ctio.pressure = 0
        return ctio

def moon_angle(observer = None):
        if observer == None:
            observer = ctio_observer()
        moon = ephem.Moon(observer)
        return moon.ra / ephem.degree, moon.dec / ephem.degree

def moon_separation(ra, dec, observer = None):
        ra_rad = ra2float(str(ra))*ephem.degree
        dec_rad = dec2float(str(dec))*ephem.degree
        m_ra, m_dec = moon_angle(observer = observer)
        m_ra = m_ra * ephem.degree
        m_dec = m_dec * ephem.degree
        return float(ephem.separation((ra_rad,dec_rad),(m_ra,m_dec))) / ephem.degree

def azel2radec(az, el, observer = None):
        """ Convert azimuth and elevation (altitude) to ra and dec for a given observer """
        if observer == None:
            observer = ctio_observer()
        az_rad = dec2float(str(az))*ephem.degree
        el_rad = dec2float(str(el))*ephem.degree
        ra_rad, dec_rad = observer.radec_of(az_rad, el_rad)
        return ra_rad / ephem.degree, dec_rad / ephem.degree

def hadec2radec(ha, dec, observer = None):
        """ Convert hour angle and declination to ra and dec for a given observer """
        if observer == None:
            observer = ctio_observer()
        h = ra2float(str(ha))
        d = dec2float(str(dec))
        st = ra2float(str(observer.sidereal_time()))
        ra =float2ra(st - h)
        dec = d
        return ra, dec

# Haversine formula example in Python
# Author: Wayne Dyck
# origin, destination = (ra,dec) pairs in either sexigesimal or float notation
# the angle is returned in float degrees
def slew_angle(origin, destination):
    lon1 = ra2float(origin[0])
    lat1 = dec2float(origin[1])
    lon2 = ra2float(destination[0])
    lat2 = dec2float(destination[1])

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

    """
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    l1 = math.radians(lon1)
    l2 = math.radians(lon2)
    n1 = math.pow( math.cos(p2) * math.sin(dlon), 2 )
    n2 = math.pow( math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dlon), 2 )
    n = math.sqrt( n1 + n2 )
    d = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dlon)
    ang = math.atan2( n,d )
    print 'done', n, d
    return ang*180.0/math.pi
    """

    return c*180.0/math.pi
