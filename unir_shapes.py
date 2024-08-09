import geopandas as gpd
import numpy as np
from shapely.geometry import Point

gpd.options.io_engine = "pyogrio"

layers = {
    "Massas_agua_NA":"NA",
    "cota_piezo":"cota_piezo",
}

datas = {
    "NA":[],
    "geometry":[],
}

# Layers normais
for layer, col_na in layers.items():
    shp:gpd.GeoDataFrame = gpd.read_file(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_ibirite.gpkg", layer=layer)
    shp.to_crs("EPSG:31983", inplace=True)

    datas["NA"].extend(shp[col_na].to_list())
    datas["geometry"].extend(shp["geometry"].to_list())

#Trechos
trechos:gpd.GeoDataFrame = gpd.read_file(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_ibirite.gpkg", layer="trechos_filtrados")
cotas:gpd.GeoDataFrame = gpd.read_file(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_ibirite.gpkg", layer="cotas")

trechos.to_crs("EPSG:31983", inplace=True)
cotas.to_crs("EPSG:31983", inplace=True)

for trecho in trechos.itertuples():
    print(f"Index: {trecho.Index} de {len(trechos.index)}")
    geom = trecho.geometry

    xy = np.array(geom.xy)

    points = gpd.points_from_xy(xy[0], xy[1])

    for point in points:
        distances = cotas.distance(point)
        min_dist_index = distances.idxmin()
        values = cotas.loc[min_dist_index]

        NA = values["ELEV"]

        datas["NA"].append(NA)
        datas["geometry"].append(point)


gpd.GeoDataFrame({"NA":datas["NA"]}, geometry=datas["geometry"], crs="EPSG:31983").to_file(
    r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_ibirite.gpkg",
    layer="cotas_all",
    engine="fiona",
)





