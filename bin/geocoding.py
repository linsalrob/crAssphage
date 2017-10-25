import sys
import os
import requests
import json

# we store the google api in a key file that is not available on github.
# You can request your own key here: 
# https://developers.google.com/maps/documentation/geocoding/get-api-key

basepath = os.path.dirname(sys.argv[0])
keyfile = None
if os.path.exists(os.path.join(basepath, "googleapi.key")):
    keyfile = os.path.join(basepath, "googleapi.key")
elif os.path.exists('googleapi.key'):
    keyfile = 'googleapi.key'
else:
    sys.stderr.write("ERROR: You need to obtain or create a Google API key before you can geocode\n")
    sys.exit(-1)

with open(keyfile, 'r') as f:
    googlekey = f.readline().strip()


def place_to_latlon(city, country, verbose=False):
    """
    Convert and city and country to a lat lon
    """

    # this is a simple fix that seems to work
    if 'USA' == country:
        country = 'US'

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

    # latlon = j['results'][0]['address_components']['location']
    try:
        latlon = j['results'][0]['geometry']['location']
    except:
        sys.stderr.write("There was an error getting the results for {}\n".format(url))
        return None, None
    return latlon['lat'], latlon['lng']


def latlon_to_place(lat, lon, verbose=False):
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
    return locality, country



if __name__ == "__main__":
    
    city="San Diego"
    country="USA"
    lat, lon = place_to_latlon(city, country, True)
    print("San Diego: {},{}".format(lat,lon))
    
    p=latlon_to_place(50.038062,19.321854)
    with open("temp", 'w', encoding='utf-8') as out:
        out.write("{}\t{}\n".format(p[0], p[1]))
    locality, country = latlon_to_place(33.116202,-117.321597)
    print("33.116202,-117.321597: {}, {}".format(locality, country))
