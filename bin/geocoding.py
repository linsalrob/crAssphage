import sys
import os
import requests
import json
import location_db as ldb

# we store the google api in a key file that is not available on github.
# You can request your own key here: 
# https://developers.google.com/maps/documentation/geocoding/get-api-key

basepath = os.path.dirname(sys.argv[0])
basepath = basepath[:basepath.rindex('/')]

keyfile = None
if os.path.exists(os.path.join(basepath, "bin/googleapi.key")):
    keyfile = os.path.join(basepath, "bin/googleapi.key")
elif os.path.exists('googleapi.key'):
    keyfile = 'googleapi.key'
else:
    sys.stderr.write("ERROR: You need to obtain or create a Google API key before you can geocode\n")
    sys.exit(-1)

with open(keyfile, 'r') as f:
    googlekey = f.readline().strip()


localitydb = os.path.join(basepath, "data/localities.db")
if not os.path.exists(localitydb):
    sys.stderr.write("ERROR: We can't find the location SQLite database ")
    sys.stderr.write("(which should be in {})\n".format(localitydb))
    sys.stderr.write("Please check your installation\n")
    sys.exit(-1)


conn = ldb.get_database_connection(localitydb)

def place_to_latlon(city, country, verbose=False):
    """
    Convert and city and country to a lat lon
    """
    
    # a simple fix
    if 'USA' == country or 'US' == country:
        country = 'United States'

    results = ldb.get_by_ascii(city, country)
    if results:
        return results[0], results[1]

    results = ldb.get_by_locale(city, country)
    if results:
        sys.stderr.write("WARNING: No results for {}, {} using get_by_ascii but we got them using get_by_locale\n".format(city, country))
        return results[0], results[1]
    
    if verbose:
        sys.stderr.write("Looking up in the database failed. Resorting to API\n")


    url = "https://maps.googleapis.com/maps/api/geocode/json?key={}&".format(googlekey)
    components = []
    if city and city != "":
        components.append("locality:{}".format(city))
    if country and country != "":
        components.append("country:{}".format(country))
    if not components:
        sys.stderr.write("ERROR: You must specify a city or a country\n")
        sys.exit(-1)
    
    c = "|".join(components)
        
    url += "components={}".format(c)
    v = requests.get(url)
    if verbose:
        sys.stderr.write("{}\n".format(url))
        sys.stderr.write("{}\n".format(v))
    j = json.loads(v.text)
    if "ZERO_RESULTS" == j['status']:
        sys.stderr.write("Could not find a place for {},{}\n".format(lat, lon))
        return None, None

    try:
        latlon = j['results'][0]['geometry']['location']
    except:
        sys.stderr.write("There was an error getting the results for {}\n".format(url))
        return None, None

    ldb.save_location(latlon['lat'], latlon['lng'], city, country, city, country)
    return latlon['lat'], latlon['lng']


def latlon_to_place(lat, lon, verbose=False):

    results = ldb.get_by_latlon(lat, lon)
    if results:
        return results[4], results[5]

    if verbose:
        sys.stderr.write("Looking up in the database failed. Resorting to API\n")

    url = "https://maps.googleapis.com/maps/api/geocode/json?key={}&".format(googlekey)
    url += "latlng={},{}".format(lat, lon)
    v = requests.get(url)
    if verbose:
        sys.stderr.write("{}\n".format(v))
    j = json.loads(v.text)
    if "ZERO_RESULTS" == j['status']:
        sys.stderr.write("Could not find a place for {},{}\n".format(lat, lon))
        return None, None

    try:
        components = j['results'][0]['address_components']
    except:
        sys.stderr.write("There was an error getting the results for {}\n".format(url))
        return None, None
    country = locality = None
    for c in components:
        if "country" in c['types']:
            country = c['long_name']
        if "locality" in c['types']:
            locality = c['long_name']

    asciiloc = locality.encode('ascii', 'ignore').decode()
    if asciiloc != locality:
        asciiloc = ""
    asciicountry = country.encode('ascii', 'ignore').decode()
    if asciicountry != country:
        asciicountry = ""

    ldb.save_location(lat, lon, locality, country, asciiloc, asciicountry)
    return locality, country



if __name__ == "__main__":
    
    city="San Diego"
    country="USA"
    lat, lon = place_to_latlon(city, country, True)
    print("San Diego: {},{}".format(lat,lon))
    sys.exit(0);
    
    p=latlon_to_place(50.038062,19.321854)
    with open("temp", 'w', encoding='utf-8') as out:
        out.write("{}\t{}\n".format(p[0], p[1]))
    locality, country = latlon_to_place(33.116202,-117.321597)
    print("33.116202,-117.321597: {}, {}".format(locality, country))
