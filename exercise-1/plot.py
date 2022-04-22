import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
from pathlib import Path

from adjustText import adjust_text

ref_cities = {
    "Berlin": (52.52, 13.38),
    "Hamburg": (53.55, 10.00),
    "Munich": (48.14, 11.58),
    "Cologne": (50.95, 6.97),
    "Frankfurt": (50.12, 8.68),
    "Stuttgart": (48.79, 9.19),
    "Düsseldorf": (51.24, 6.79),
    "Dortmund": (51.51, 7.48),
    "Essen": (51.47, 7.00),
    "Leipzig": (51.35, 12.40),
    "Bremen": (53.08, 8.81),
    "Dresden": (51.05, 13.74),
    "Hanover": (52.40, 9.73),
    "Nuremberg": (49.45, 11.05),
}

# fmt: off
ref_dots = [
    [52.52, 13.38], [53.55, 10.00], [48.14, 11.58], [50.73, 7.10], [51.51, 7.11],
    [50.95, 6.97], [50.12, 8.68], [48.79, 9.19], [51.24, 6.79], [51.51, 7.48],
    [51.47, 7.00], [51.35, 12.40], [53.08, 8.81], [51.05, 13.74], [52.40, 9.73],
    [49.45, 11.05], [51.43, 6.75], [51.48, 7.20], [51.26, 7.18], [52.03, 8.53], 
    [51.96, 7.62], [49.00, 8.40], [49.50, 8.47], [48.36, 10.89], [50.08, 8.23], 
    [51.20, 6.42], [52.27, 10.51], [50.83, 12.92], [54.32, 10.12], [50.77, 6.09],
    [51.48, 11.96], [52.13, 11.62], [47.99, 7.85], [51.33, 6.55], [53.87, 10.66],
    [51.47, 6.86], [50.99, 11.03], [50.00, 8.26], [54.09, 12.10], [51.32, 9.48],
]
# fmt: on

edges_file = "./data/delta_2_3.txt"
edges_file_content = open(edges_file).read()

local_path = "tmp/"
local_gdl_copy = local_path + "/DEU_adm1.shp"
url = "https://biogeo.ucdavis.edu/data/diva/adm/DEU_adm.zip"

if not Path(local_gdl_copy).is_file():
    import io
    import zipfile
    import requests

    r = requests.get(url)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall(path=local_path)

gdf = gpd.read_file(local_path + "/DEU_adm1.shp")

plt.rcParams["figure.dpi"] = 100

texts = []
for k in ref_cities.keys():
    texts.append(plt.text(ref_cities[k][1], ref_cities[k][0], k, zorder=3))

plt.scatter([d[1] for d in ref_dots], [d[0] for d in ref_dots], c="k", zorder=2)

adjust_text(texts)

for g in gdf.geometry:
    s = str(g).replace("MULTIPOLYGON (", "").replace("POLYGON (", "")

    polys = s.split("), (")

    for poly in polys:
        poly_points = poly.replace("(", "").replace(")", "").split(", ")
        xs = []
        ys = []
        for p in poly_points:
            x, y = p.split(" ")
            xs.append(float(x))
            ys.append(float(y))
        plt.plot(xs[::7], ys[::7], c="k", linewidth=0.4, zorder=0)

for line in edges_file_content.split("\n"):
    coords = line.split(" ")
    plt.plot(
        [float(coords[0]), float(coords[2])],
        [float(coords[1]), float(coords[3])],
        c="r",
        zorder=1,
    )

plt.axis('off')

# plt.show()

import tikzplotlib

tikzplotlib.save(edges_file + ".pgf")
