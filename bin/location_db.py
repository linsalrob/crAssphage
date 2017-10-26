"""

Create a location database for our data. This is a SQLite3 database that
has a single table with the following attributes:

    0. latitude
    1. longitude
    2. locality: the name of the city, municipality, or province in utf-8 encoding
    3. country: the name of the country in utf-8 encoding
    4. ascii_locality: the name of the city, municipality, or province in pure ascii
    5. ascii_country: the name of the country in pure ascii

NOTE:
    The latitude and longitude are stored as REALs (ie. floats) and so you
    may or may not need to convert them from str. I made this decision
    because they are floating point numbers and there are some
    GIS calculations we could potentially do using them.

"""


import os
import sys
import sqlite3

connection = None

def get_database_connection(dbfile='../data/localities.db'):
    """
    Create a connection to the database
    """

    global connection
    if not connection:
        connection = sqlite3.connect(dbfile)
    
    return connection


def get_database_cursor(conn=None):
    """
    Connect to the database and get a cursor
    
    :param conn: the optional connection. We will make it if not provided
    """

    if not conn:
        conn = get_database_connection()

    return conn.cursor()
    

def create_database(cursor=None):
    """
    Create the database tables

    :param cursor: the database cursor. If not provided we'll make one and return it
    """

    if not cursor:
        cursor = get_database_cursor()
    cursor.execute('''CREATE TABLE location (latitude real, longitude real, locality text, country text, ascii_locality text, ascii_country text)''')
    get_database_connection().commit()
    return cursor

def save_location(latitude, longitude, locality, country, ascii_locality, ascii_country, cursor=None):
    """
    Insert something into the database
    :param latitude: the latitude in decimal degrees (signed float)
    :param longitude: the longitude in decimal degrees (signed float)
    :param locality: the town or metropolitan area in utf-8 text
    :param country: the country in utf-8
    :param ascii_locality: the locality in ascii format
    :param ascii_country: the country in ascii format
    :param cursor: the database cursor. If not provided we'll make one. We return this
    """
    
    if not cursor:
        cursor = get_database_cursor()

    # first check to see if it exists
    cursor.execute("select * from location where latitude=? and longitude=?", (latitude, longitude))
    result = cursor.fetchone()
    if result:
        sys.stderr.write("We already have an entry for {}, {}: {}, {}\n".format(result[0], result[1], result[4], result[5]))
        return cursor
    cursor.execute("insert into location values (?, ?, ?, ?, ?, ?)", (latitude, longitude, locality, country, ascii_locality, ascii_country))
    get_database_connection().commit()
    return cursor

def get_by_latlon(latitude, longitude, cursor=None):
    """
    Get the location based on the latitude and longitude

    :param latitude: the latitude in decimal degrees (signed float)
    :param longitude: the longitude in decimal degrees (signed float)
    :param cursor: the db cursor. We can get this
    :return : the array of all data
    """

    if not cursor:
        cursor = get_database_cursor()
    cursor.execute("select * from location where latitude=? and longitude=?", (latitude, longitude))
    return cursor.fetchone()


def get_by_ascii(ascii_locality, ascii_country, cursor=None):
    """
    Get the lat lon based on the ascii location
    :param ascii_locality: the locality in ascii format
    :param ascii_country: the country in ascii format
    :return : the array of all the data
    """

    if not cursor:
        cursor = get_database_cursor()
    cursor.execute("select * from location where ascii_locality=? and ascii_country=?", (ascii_locality, ascii_country))
    return cursor.fetchone()

def get_by_locale(locality, country, cursor=None):
    """
    Get the lat lon based on the ascii location
    :param locality: the locality in ascii format
    :param country: the country in ascii format
    :return : the array of all the data
    """

    if not cursor:
        cursor = get_database_cursor()
    cursor.execute("select * from location where locality=? and country=?", (locality, country))
    return cursor.fetchone()

    


