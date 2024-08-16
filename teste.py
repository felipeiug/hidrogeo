import geopandas as gpd
import pandas as pd
import numpy as np
gpd.options.io_engine = "pyogrio"
import rioxarray
import xarray as xr
from shapely.geometry import mapping, box, Point
import matplotlib.pyplot as plt
from skimage import measure
import rasterio

from pykrige.ok import OrdinaryKriging
from pykrige.uk import UniversalKriging

crs = "EPSG:31983"

# Lendo os arquivos
# altimetria = rioxarray.open_rasterio(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_NA\cota_terreno.tif")
area = gpd.read_file(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_NA\area_drenagem.zip")
litografia = gpd.read_file(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_NA\litoestratigrafia.zip")
litografia = litografia[litografia.ID_UNIDADE==467]

trechos = gpd.read_file(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_NA\hidrografia.zip")
pocos = gpd.read_file(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_NA\profundidade_pocos.zip")
outorgas = gpd.read_file(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_NA\outorgas_mg.zip")

# Reprojetando as geometrias
area.to_crs(crs, inplace = True)
litografia.to_crs(crs, inplace = True)
area:gpd.GeoDataFrame = area.clip(litografia)


trechos.to_crs(crs, inplace = True)
pocos.to_crs(crs, inplace = True)
outorgas.to_crs(crs, inplace = True)
# altimetria = altimetria.rio.reproject(crs)

# Recortando o raster
bounds = area.total_bounds ############ Area 1 ou Area 2
geom = box(*bounds)
# geo_df = gpd.GeoDataFrame({'geometry': [geom]}, crs=altimetria.rio.crs)
# altimetria = altimetria.rio.clip(geo_df.geometry.apply(mapping), geo_df.crs)

# Obtendo as coordenadas
# band = altimetria.sel(band=1)
# x_coords, y_coords = band['x'].values, band['y'].values
# altitude = band.values

# Altitudes piezométicas conhecidas
na = np.array([])
x_na = np.array([])
y_na = np.array([])


pocos = pocos[pocos.within(area.unary_union)]
x = pocos.geometry.x.values
y = pocos.geometry.y.values
na_ = pocos.nivel_esta.values
na_ = np.where(na_ == None, "-1", na_)
for i, na_now in enumerate(na_):
    na_now = na_now.replace(",", ".")
    na_now = float(na_now)
    na_[i] = na_now
x_na = np.append(x_na, x)
y_na = np.append(y_na, y)
na = np.append(na, na_)

# use_rios
# trechos:gpd.GeoDataFrame = trechos[trechos.nuStharler > nu_stharler_min]
# trechos:gpd.GeoDataFrame = trechos.clip(area.unary_union)
# trechos_simpliy:gpd.GeoDataFrame = trechos.simplify(tolerance, preserve_topology=True)

# xys = (np.array(i.xy) for i in trechos_simpliy.geometry.values)

# for xy in xys:
#     x = xy[0]
#     y = xy[1]

#     for point in range(x.size):
#         if point%n_points != 0:
#             continue

#         x_na = np.append(x_na, x[point])
#         y_na = np.append(y_na, y[point])
#         na = np.append(na, 0)

outorgas = outorgas[outorgas.within(area.unary_union)]
x = outorgas.geometry.x.values
y = outorgas.geometry.y.values
na_ = outorgas.neest_4.values
na_ = np.where((na_ == None) | (na_ == 'Outros'), "-1", na_)
for i, na_now in enumerate(na_):
    na_now = na_now.replace(",", ".")
    na_now = float(na_now)
    na_[i] = na_now
x_na = np.append(x_na, x)
y_na = np.append(y_na, y)
na = np.append(na, na_)


na = na.astype(float)
mask = na > 0
x_na = x_na[mask]
y_na = y_na[mask]
na = na[mask]

df = pd.DataFrame({"COTA": na, "X": x_na, "Y": y_na})
df['geometry'] = df.apply(lambda row: Point(row['X'], row['Y']), axis=1)
points = gpd.GeoDataFrame(df, geometry='geometry', crs=crs)

points.to_file("pontos.shp")

