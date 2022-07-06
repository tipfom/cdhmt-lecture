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
(52.5167, 13.3833), 
(53.55, 10.0), 
(48.1372, 11.5755), 
(50.9422, 6.9578), 
(50.1136, 8.6797), 
(53.1153, 8.7975), 
(51.2311, 6.7724), 
(48.7761, 9.1775), 
(51.3333, 12.3833), 
(51.5139, 7.4653), 
(51.4508, 7.0131), 
(51.05, 13.74), 
(52.3744, 9.7386), 
(49.4539, 11.0775), 
(51.4322, 6.7611), 
(51.4833, 7.2167), 
(51.2667, 7.1833), 
(52.0167, 8.5333), 
(50.7339, 7.0997), 
(51.9625, 7.6256), 
(49.0167, 8.4), 
(49.4878, 8.4661), 
(48.3717, 10.8983), 
(51.3166, 9.4912), 
(50.0825, 8.24), 
(51.2, 6.4333), 
(51.5167, 7.1), 
(52.2692, 10.5211), 
(50.7762, 6.0838), 
(54.3233, 10.1394), 
(50.8333, 12.9167), 
(51.4828, 11.9697), 
(52.1278, 11.6292), 
(47.9947, 7.8497), 
(51.3333, 6.5667), 
(50.0, 8.2667), 
(53.8697, 10.6864), 
(51.4699, 6.8514), 
(54.0833, 12.1333), 
(50.9787, 11.0328), 
(51.3594, 7.475), 
(52.4, 13.0667), 
(49.2333, 7.0), 
(51.6667, 7.8167), 
(49.4811, 8.4353), 
(51.4275, 6.8825), 
(53.1439, 8.2139), 
(52.2789, 8.0431), 
(51.0333, 6.9833), 
(49.4122, 8.71), 
(51.1667, 7.0833), 
(49.8667, 8.65), 
(51.5426, 7.219), 
(51.2003, 6.6939), 
(49.0167, 12.0833), 
(51.7167, 8.7667), 
(48.7636, 11.4261), 
(49.7944, 9.9294), 
(49.4783, 10.9903), 
(48.3984, 9.9916), 
(49.1404, 9.218), 
(48.895, 8.705), 
(52.4231, 10.7872), 
(51.5339, 9.9356), 
(51.5232, 6.9253), 
(48.4833, 9.2167), 
(50.3597, 7.5978), 
(53.55, 8.5833), 
(51.6167, 7.5167), 
(49.5964, 11.0044), 
(50.9272, 11.5864), 
(51.1802, 7.1872), 
(49.7567, 6.6414), 
(52.1503, 10.3593), 
(51.4592, 6.6197), 
(50.8756, 8.0167), 
(52.15, 9.95), 
(51.7606, 14.3342), 
(51.9, 8.3833), 
(49.4447, 7.7689), 
(51.4333, 7.3333), 
(50.1328, 8.9169), 
(53.6333, 11.4167), 
(50.8782, 12.0824), 
(48.7406, 9.3108), 
(48.8975, 9.1919), 
(51.3833, 7.6667), 
(50.8, 6.4833), 
(48.52, 9.0556), 
(54.7819, 9.4367), 
(50.7189, 12.4961), 
(50.5833, 8.6667), 
(51.3, 6.85), 
(51.6167, 7.5167), 
(48.0603, 8.4586), 
(47.6633, 9.1753), 
(51.6667, 7.1167), 
(49.6319, 8.3653), 
(51.34, 7.0416), 
]

# fmt: on

edges_file = "./data/1_3.txt"
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

plt.show()

# import tikzplotlib

# tikzplotlib.save(edges_file + ".pgf")
