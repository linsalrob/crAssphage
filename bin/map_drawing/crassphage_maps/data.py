"""

Data that we use for drawing the maps

"""


country2continent = {"Albania": "Europe", "Australia": "Oceania", "Azerbaijan": "Asia", "Belgium": "Europe",
                     "Brazil": "South America", "Bulgaria": "Europe", "Canada": "North America", "China": "Asia",
                     "Colombia": "South America", "Cote_dIvoire": "Africa", "Croatia": "Europe",
                     "Czech_Republic": "Europe", "Denmark": "Europe", "Ecuador": "South America",
                     "Ethiopia": "Africa", "Finland": "Europe", "France": "Europe", "Gambia": "Africa",
                     "Germany": "Europe", "Hong_Kong": "Asia", "Hungary": "Europe", "Iceland": "Europe",
                     "India": "Asia", "Iran": "Asia", "Ireland": "Europe", "Israel": "Asia", "Italy": "Europe",
                     "Japan": "Asia", "Jordan": "Asia", "Kazakhstan": "Asia", "Kenya": "Africa", "Latvia": "Europe",
                     "Luxembourg": "Europe", "Malaysia": "Asia", "Malta": "Europe", "Mexico": "North America",
                     "Moldova": "Europe", "Mongolia": "Asia", "Nepal": "Asia", "New_Zealand": "Oceania",
                     "Nigeria": "Africa", "Norway": "Europe", "Pakistan": "Asia", "Peru": "South America",
                     "Poland": "Europe", "Portugal": "Europe", "Russia": "Asia", "Singapore": "Asia",
                     "South_Africa": "Africa", "Spain": "Europe", "Sri_Lanka": "Asia", "Sudan": "Africa",
                     "Sweden": "Europe", "Switzerland": "Europe", "The_Netherlands": "Europe",
                     "United_States": "North America", "USA": "North America", "Zambia": "Africa",
                     "crAssphage_A": "crAssphage_A"}

# these are arrays of 50 colors from http://www.perbang.dk/rgbgradient/ between FF0000 (red) and 0000FF (blue) or
# #00FF00 (green) and #FFFF00 (yellow)
green2yellow = ["#00FF00", "#05FF00", "#0AFF00", "#0FFF00", "#14FF00", "#1AFF00", "#1FFF00", "#24FF00", "#29FF00",
                "#2EFF00", "#34FF00", "#39FF00", "#3EFF00", "#43FF00", "#48FF00", "#4EFF00", "#53FF00", "#58FF00",
                "#5DFF00", "#62FF00", "#68FF00", "#6DFF00", "#72FF00", "#77FF00", "#7CFF00", "#82FF00", "#87FF00",
                "#8CFF00", "#91FF00", "#96FF00", "#9CFF00", "#A1FF00", "#A6FF00", "#ABFF00", "#B0FF00", "#B6FF00",
                "#BBFF00", "#C0FF00", "#C5FF00", "#CAFF00", "#D0FF00", "#D5FF00", "#DAFF00", "#DFFF00", "#E4FF00",
                "#EAFF00", "#EFFF00", "#F4FF00", "#F9FF00", "#FFFF00"]
red2blue = ["#FF0000", "#F90005", "#F4000A", "#EF000F", "#EA0014", "#E4001A", "#DF001F", "#DA0024", "#D50029",
            "#D0002E", "#CA0034", "#C50039", "#C0003E", "#BB0043", "#B60048", "#B0004E", "#AB0053", "#A60058",
            "#A1005D", "#9C0062", "#960068", "#91006D", "#8C0072", "#870077", "#82007C", "#7C0082", "#770087",
            "#72008C", "#6D0091", "#680096", "#62009C", "#5D00A1", "#5800A6", "#5300AB", "#4E00B0", "#4800B6",
            "#4300BB", "#3E00C0", "#3900C5", "#3400CA", "#2E00D0", "#2900D5", "#2400DA", "#1F00DF", "#1A00E4",
            "#1400EA", "#0F00EF", "#0A00F4", "#0500F9", "#0000FF"]
green2red = ["#00FF00", "#05F900", "#0AF400", "#0FEF00", "#14EA00", "#1AE400", "#1FDF00", "#24DA00", "#29D500",
             "#2ED000", "#34CA00", "#39C500", "#3EC000", "#43BB00", "#48B600", "#4EB000", "#53AB00", "#58A600",
             "#5DA100", "#629C00", "#689600", "#6D9100", "#728C00", "#778700", "#7C8200", "#827C00", "#877700",
             "#8C7200", "#916D00", "#966800", "#9C6200", "#A15D00", "#A65800", "#AB5300", "#B04E00", "#B64800",
             "#BB4300", "#C03E00", "#C53900", "#CA3400", "#D02E00", "#D52900", "#DA2400", "#DF1F00", "#E41A00",
             "#EA1400", "#EF0F00", "#F40A00", "#F90500", "#FF0000"]