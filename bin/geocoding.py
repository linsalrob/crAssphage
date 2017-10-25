import sys
import requests
import json

# we store the google api in a key file that is not available on github.
# You can request your own key here: 
# https://developers.google.com/maps/documentation/geocoding/get-api-key
with open('googleapi.key', 'r') as f:
    googlekey = f.readline().strip()


def place_to_latlon(city, country, verbose=False):
    """
    Convert and city and country to a lat lon
    """
    url = "https://maps.googleapis.com/maps/api/geocode/json?key={}&".format(googlekey)
    url += "components=locality:{}|country={}".format(city, country)
    v = requests.get(url)
    if verbose:
        sys.stderr.write("{}\n".format(url))
        sys.stderr.write("{}\n".format(v))
    j = json.loads(v.text)
    # latlon = j['results'][0]['address_components']['location']
    latlon = j['results'][0]['geometry']['location']
    return latlon['lat'], latlon['lng']


def latlon_to_place(lat, lon, verbose=False):
    url = "https://maps.googleapis.com/maps/api/geocode/json?key={}&".format(googlekey)
    url += "latlng={},{}".format(lat, lon)
    v = requests.get(url)
    if verbose:
        sys.stderr.write("{}\n".format(v))
    j = json.loads(v.text)
    components = j['results'][0]['address_components']
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
