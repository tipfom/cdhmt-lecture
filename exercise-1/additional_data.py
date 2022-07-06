from math import dist, radians, sin, cos, acos
def great_circle(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return 6371 * (
        acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
    )


cities_raw = open("./german_cities.csv").readlines()

N = 100

import numpy as np
from dataclasses import dataclass

@dataclass
class City:
    name: str
    pop: int
    lat: float
    lon: float

cities = {}

total_pop = 0

for i in range(N):
    city_raw = cities_raw[i].split(";")
    cities[i] = City(city_raw[0], int(city_raw[1]), float(city_raw[2].replace(",",".")), float(city_raw[3].replace(",",".")))
    total_pop += cities[i].pop

distances = np.zeros((N, N))
for i in range(N):
    for j in range(i+1, N):
        distances[i, j] = great_circle(cities[i].lon, cities[i].lat, cities[j].lon, cities[j].lat)

distances = distances / np.average(distances)

with open("./cities_new.dat", "x") as f_cities:
    with open("./distances_new.dat", "x") as f_distances:
        for i in range(N):
            f_cities.write(f"{i}\t\"{cities[i].name}\"\t{cities[i].pop / total_pop}\t{cities[i].lat}\t{cities[i].lon}\n")
            for j in range(i+1, N):
                f_distances.write(f"{i}\t{j}\t{distances[i][j]}\n")
