import sys
import os
import argparse
import requests
import json
import location_db as ldb

# we store the google api in a key file that is not available on github.
# You can request your own key here: 
# https://developers.google.com/maps/documentation/geocoding/get-api-key

basepath = os.path.dirname(sys.argv[0])
if basepath:
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
sys.stderr.write("Reading database: {}\n".format(localitydb))
if not os.path.exists(localitydb):
    localitydb = '../data/localities.db'
    if not os.path.exists(localitydb):
        sys.stderr.write("ERROR: We can't find the location SQLite database ")
        sys.stderr.write("(which should be in {})\n".format(localitydb))
        sys.stderr.write("Please check your installation\n")
        sys.exit(-1)


conn = ldb.get_database_connection(localitydb)

def place_to_latlon(city, country, verbose=False, force_api=False):
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


def latlon_to_place(lat, lon, verbose=False, force_api=False):

    results = ldb.get_by_latlon(lat, lon)
    if results and not force_api:
        return results[4], results[5]

    if verbose:
        sys.stderr.write("Looking up in the database failed. Resorting to API\n")

    url = "https://maps.googleapis.com/maps/api/geocode/json?key={}&".format(googlekey)
    url += "latlng={},{}".format(lat, lon)
    v = requests.get(url)
    if verbose:
        sys.stderr.write("{}\n{}\n".format(url, v))
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
    parser = argparse.ArgumentParser(description='Find geocoordinates')
    parser.add_argument('-l', help='Lat/Lon in format lat,lon (eg: -14.235004,-51.92528)')
    parser.add_argument('-p', help='Place (eg. San Diego) (requires a country)')
    parser.add_argument('-c', help='Country')
    parser.add_argument('-f', help='force a web search vs a db search', action='store_true')
    parser.add_argument('-v', help='verbose', action='store_true')
    args = parser.parse_args()

    
    if args.p and args.c:
        lat, lon = place_to_latlon(args.p, args.c, args.v, args.f)
        print("{}, {}: {},{}".format(args.p, args.c, lat,lon))
    elif args.c:
        lat, lon = place_to_latlon(None, args.c, args.v, args.f)
        print("{}: {},{}".format(args.c, lat,lon))
    elif args.l:
        (lat, lon) = args.l.split(',')
        locality, country = latlon_to_place(lat, lon, args.v, args.f)
        print("{},{} : {}, {}".format(lat, lon, locality, country))
