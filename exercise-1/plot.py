cities = {
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
    # "Duisburg":	(51.43, 6.75),
    # "Bochum":	(51.48, 7.20),
    # "Wuppertal":	(51.26, 7.18),
    # "Bielefeld":	(52.03, 8.53),
    # "Bonn":	(50.73, 7.10),
    # "Munster":	(51.96, 7.62),
    # "Karlsruhe":	(49.00, 8.40),
    # "Mannheim":	(49.50, 8.47),
    # "Augsburg":	(48.36, 10.89),
    # "Wiesbaden":	(50.08, 8.23),
    # "Gelsenkirchen":	(51.51, 7.11),
    # "Mönchengladbach":	(51.20, 6.42),
    # "Brunswick":	(52.27, 10.51),
    # "Chemnitz":	(50.83, 12.92),
    # "Kiel":	(54.32, 10.12),
    # "Aachen":	(50.77, 6.09),
    # "Halle":	(51.48, 11.96),
    # "Magdeburg":	(52.13, 11.62),
    # "Freiburg":	(47.99, 7.85),
    # "Krefeld":	(51.33, 6.55),
    # "Lübeck":	(53.87, 10.66),
    # "Oberhausen":	(51.47, 6.86),
    # "Erfurt":	(50.99, 11.03),
    # "Mainz":	(50.00, 8.26),
    # "Rostock":	(54.09, 12.10),
    # "Kassel":	(51.32, 9.48),
}

dots = [
    (52.52, 13.38),
    (53.55, 10.00),
    (48.14, 11.58),
    (50.95, 6.97),
    (50.12, 8.68),
    (48.79, 9.19),
    (51.24, 6.79),
    (51.51, 7.48),
    (51.47, 7.00),
    (51.35, 12.40),
    (53.08, 8.81),
    (51.05, 13.74),
    (52.40, 9.73),
    (49.45, 11.05),
    (51.43, 6.75),
    (51.48, 7.20),
    (51.26, 7.18),
    (52.03, 8.53),
    (50.73, 7.10),
    (51.96, 7.62),
    (49.00, 8.40),
    (49.50, 8.47),
    (48.36, 10.89),
    (50.08, 8.23),
    (51.51, 7.11),
    (51.20, 6.42),
    (52.27, 10.51),
    (50.83, 12.92),
    (54.32, 10.12),
    (50.77, 6.09),
    (51.48, 11.96),
    (52.13, 11.62),
    (47.99, 7.85),
    (51.33, 6.55),
    (53.87, 10.66),
    (51.47, 6.86),
    (50.99, 11.03),
    (50.00, 8.26),
    (54.09, 12.10),
    (51.32, 9.48),
]

ss = """13.38 52.52 12.4 51.35
13.38 52.52 11.62 52.13
13.38 52.52 10.66 53.87
10 53.55 8.81 53.08
10 53.55 10.51 52.27
10 53.55 10.66 53.87
10 53.55 12.1 54.09
11.58 48.14 11.05 49.45
11.58 48.14 10.89 48.36
6.97 50.95 6.79 51.24
6.97 50.95 7.18 51.26
6.97 50.95 7.1 50.73
6.97 50.95 6.09 50.77
8.68 50.12 11.05 49.45
8.68 50.12 8.23 50.08
8.68 50.12 11.03 50.99
8.68 50.12 8.26 50
8.68 50.12 9.48 51.32
9.19 48.79 10.89 48.36
9.19 48.79 8.26 50
6.79 51.24 6.75 51.43
6.79 51.24 7.2 51.48
6.79 51.24 6.42 51.2
6.79 51.24 6.55 51.33
7.48 51.51 7.2 51.48
7.48 51.51 7.18 51.26
7.48 51.51 8.53 52.03
7.48 51.51 7.62 51.96
7 51.47 7.11 51.51
7 51.47 6.86 51.47
12.4 51.35 13.74 51.05
12.4 51.35 11.05 49.45
12.4 51.35 11.96 51.48
8.81 53.08 7.2 51.48
8.81 53.08 10.51 52.27
9.73 52.4 8.53 52.03
9.73 52.4 10.51 52.27
11.05 49.45 10.51 52.27
6.75 51.43 6.86 51.47
7.2 51.48 7.18 51.26
7.2 51.48 7.11 51.51
8.53 52.03 7.62 51.96
8.53 52.03 9.48 51.32
7.1 50.73 8.23 50.08
8.4 49 8.23 50.08
8.47 49.5 8.26 50
8.23 50.08 8.26 50
10.51 52.27 11.62 52.13
10.51 52.27 9.48 51.32
12.92 50.83 11.96 51.48
10.12 54.32 10.66 53.87
11.96 51.48 11.62 52.13
11.96 51.48 11.03 50.99
7.85 47.99 8.26 50"""

import matplotlib.pyplot as plt
import io
import zipfile
import requests
import geopandas as gpd

local_path = "tmp/"
# url = "https://biogeo.ucdavis.edu/data/diva/adm/DEU_adm.zip"

# r = requests.get(url)
# z = zipfile.ZipFile(io.BytesIO(r.content))
# z.extractall(path=local_path)
import numpy as np

gdf = gpd.read_file(local_path + "/DEU_adm1.shp")
# gdf.plot(color='white',
#          edgecolor='black')
plt.rcParams["figure.dpi"] = 100

print(gdf.geometry.head())

c = []


# plt.gca().set_aspect(1.4)
# plt.gca().axis('equal')
# plt.gca().set_xlim([5, 15])
# plt.gca().set_ylim([46, 56])

# plt.show()

texts = []

for k in cities.keys():
    texts.append(plt.text(cities[k][1], cities[k][0], k, zorder=3))
    # plt.annotate(k, [])

for d in dots:
    plt.scatter([d[1]], d[0], c="k", zorder=2)

from adjustText import adjust_text

adjust_text(
    texts
)  # , force_points=1.0,arrowprops=dict(arrowstyle='->', color='blue', alpha=0.4))

for g in gdf.geometry:
    s = str(g)

    if s.startswith("MULTIPOLYGON ((("):
        s = s.replace("MULTIPOLYGON (((", "").replace(")))", "")
        polys = s.split(")), ((")

        polys2 = []
        for poly in polys:
            polys2.extend(poly.split("), ("))

        for poly in polys2:
            e = poly.split(", ")
            xs = []
            ys = []
            for p in e:
                x, y = p.split(" ")
                xs.append(float(x))
                ys.append(float(y))
            plt.plot(xs[::], ys[::], c="k", linewidth=0.4, zorder=0)
    elif s.startswith("POLYGON (("):
        polys = s.replace("POLYGON ((", "").replace("))", "").split("), (")

        for poly in polys:
            e = poly.split(", ")
            xs = []
            ys = []
            for p in e:
                x, y = p.split(" ")
                xs.append(float(x))
                ys.append(float(y))
            plt.plot(xs[::], ys[::], c="k", linewidth=0.4, zorder=0)


for k in cities.keys():
    # texts.append(plt.text(cities[k][1], cities[k][0], k, zorder=3))
    # plt.annotate(k, [])
    plt.scatter([cities[k][1]], [cities[k][0]], c="k", zorder=2)

for line in ss.split("\n"):
    coords = line.split(" ")
    plt.plot(
        [float(coords[0]), float(coords[2])],
        [float(coords[1]), float(coords[3])],
        c="r",
        zorder=1,
    )

plt.show()

import tikzplotlib

tikzplotlib.save("test.pgf")
