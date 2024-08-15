import geopandas as gpd
import numpy as np
gpd.options.io_engine = "pyogrio"
import rioxarray
import xarray as xr
from shapely.geometry import mapping, box, LineString
import matplotlib.pyplot as plt
from skimage import measure
import rasterio

from pykrige.ok import OrdinaryKriging
from pykrige.uk import UniversalKriging


crs = "EPSG:31983"
use_pocos = True # Se é pra utilizar os poços para gerar a krigeagem
use_outorgas = True # Se é pra utilizar as outorgas
use_rios = True # Se é pra utilizar os trechos para fazer a krigeagem
nu_stharler_min = 3 # Número para filtrar os rios
n_points = 3 # A cada n pontos a altitude de um é considerada
tolerance = 0.5 # simplificação dos rios

# Lendo os arquivos
altimetria = rioxarray.open_rasterio(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_NA\cota_terreno.tif")
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
altimetria = altimetria.rio.reproject(crs)

# Recortando o raster
bounds = area.total_bounds ############ Area 1 ou Area 2
geom = box(*bounds)
geo_df = gpd.GeoDataFrame({'geometry': [geom]}, crs=altimetria.rio.crs)
altimetria = altimetria.rio.clip(geo_df.geometry.apply(mapping), geo_df.crs)

# Obtendo as coordenadas
band = altimetria.sel(band=1)
x_coords, y_coords = band['x'].values, band['y'].values
altitude = band.values

def get_altitude_by_points(pontos_na: np.ndarray):
    pontos_na = pontos_na.copy()

    # Cota altimétrica relativa ao mar
    x_res = float(abs(altimetria['x'][1] - altimetria['x'][0]))
    y_res = float(abs(altimetria['y'][1] - altimetria['y'][0]))

    # Calcular as coordenadas dos vértices das células
    x_start = x_coords - x_res / 2
    x_end = x_coords + x_res / 2
    y_start = y_coords - y_res / 2
    y_end = y_coords + y_res / 2

    for i, ponto_na in enumerate(pontos_na[:, 2]):
        x = pontos_na[i, 0]
        y = pontos_na[i, 1]

        m1 = (x >= x_start) & (x <= x_end)
        m2 = (y >= y_start) & (y <= y_end)

        alt = altitude[m2]
        if alt.size == 0:
            pontos_na[i, 2] = -1
            continue

        alt = alt[0][m1]
        if alt.size == 0:
            pontos_na[i, 2] = -1
            continue

        pontos_na[i, 2] = alt.mean() - ponto_na

    return pontos_na

def generate_contornos(raster, column_name = "ELEV", dist=10)->gpd.GeoDataFrame:
    raster = raster.squeeze()
    transform = raster.rio.transform()

    contour_lines = []
    contour_values = []

    for limiar in range(int(np.floor(raster.values.min())), int(np.ceil(raster.values.max())), dist):
        mask = (raster.values > limiar)
        contours = measure.find_contours(mask, level=0.5)

        for contour in contours:
            # Converter as coordenadas do raster para coordenadas geoespaciais
            coords = [rasterio.transform.xy(transform, int(point[0]), int(point[1])) for point in contour]
            # values = np.array([raster.values[int(point[0]), int(point[1])] for point in contour])

            # Criar a geometria LineString
            line = LineString(coords)
            contour_lines.append(line)
            contour_values.append(limiar)

    # Criar um GeoDataFrame com os contornos
    gdf = gpd.GeoDataFrame({column_name: contour_values}, geometry=contour_lines, crs=crs)
    return gdf


# Altitudes piezométicas conhecidas
na = np.array([])
x_na = np.array([])
y_na = np.array([])

# Dados altimétricos
if use_pocos:
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

if use_rios:
    trechos:gpd.GeoDataFrame = trechos[trechos.nuStharler > nu_stharler_min]
    trechos:gpd.GeoDataFrame = trechos.clip(area.unary_union)
    trechos_simpliy:gpd.GeoDataFrame = trechos.simplify(tolerance, preserve_topology=True)
    
    xys = (np.array(i.xy) for i in trechos_simpliy.geometry.values)

    for xy in xys:
        x = xy[0]
        y = xy[1]

        for point in range(x.size):
            if point%n_points != 0:
                continue

            x_na = np.append(x_na, x[point])
            y_na = np.append(y_na, y[point])
            na = np.append(na, 0)

if use_outorgas:
    outorgas = outorgas[outorgas.within(area.unary_union)]
    x = outorgas.geometry.x.values
    y = outorgas.geometry.y.values
    na_ = outorgas.neest_4.values
    na_ = np.where((na_ == None) | (na_ == 'Outros'), "-1", na_)

    outorgas: gpd.GeoDataFrame = outorgas[~((na_ == None) | (na_ == 'Outros'))]
    outorgas.to_file("outorgas.shp")

    for i, na_now in enumerate(na_):
        na_now = na_now.replace(",", ".")
        na_now = float(na_now)
        na_[i] = na_now

    x_na = np.append(x_na, x)
    y_na = np.append(y_na, y)
    na = np.append(na, na_)


# Pontos
na = na.astype(float)
pontos_na = np.column_stack((x_na, y_na, na))
pontos_na = pontos_na[~np.isnan(pontos_na[:, 2])]
pontos_na = pontos_na[pontos_na[:, 2] >= 0]

# Visualização dos pontos escolhidos
# area.plot()
# trechos.plot()
# plt.scatter(pontos_na[:, 0], pontos_na[:, 1])
# plt.show()

# Cálculo do valor do NA final
null_raster = altimetria.copy(deep=True)
null_raster.data[...] = -1

variance_raster = altimetria.copy(deep=True)
variance_raster.data[...] = -1

for limit_by_cota in [True, False]:
    for tipo_interpolacao in ["INVDIST", "KRIGE", "KRIGEALT", "INVALT"]:
        # Krigagem
        if tipo_interpolacao == "KRIGE":
            # Ruido para quando der erro
            ruido_0 = np.random.normal(scale=1, size=pontos_na.shape[0])
            ruido_1 = np.random.normal(scale=1, size=pontos_na.shape[0])

            pontos_na[:, 0] += ruido_0
            pontos_na[:, 1] += ruido_1

            pontos = get_altitude_by_points(pontos_na)

            z = pontos[:, 2]

            OK = UniversalKriging(
                x=pontos[:, 0], y=pontos[:, 1], z=z,
                variogram_model="linear",
            )

            for banda in range(altimetria.shape[0]):  # Para cada banda (se houver mais de uma)
                for index_y in range(altimetria.shape[1]):  # Para cada linha
                    if index_y%10 == 0:
                        print(index_y, altimetria.shape[1])

                    points_x:list[float] = []
                    points_y:list[float] = [y_coords[index_y]]*altimetria.shape[2]

                    for index_x in range(altimetria.shape[2]):  # Para cada coluna

                        px = x_coords[index_x]
                        points_x.append(px)

                    z_interp, ss = OK.execute('points', np.array([points_x]), np.array([points_y]))

                    null_raster.values[banda, index_y] = z_interp
                    variance_raster.values[banda, index_y] = ss

                if limit_by_cota:
                    variance_now = variance_raster.values[banda]
                    altitude_now = null_raster.values[banda]
                    altitude_max = altimetria.values[banda]

                    mask = altitude_max < altitude_now

                    altitude_now = np.where(mask, altitude_max, altitude_now)
                    variance_now = np.where(mask, -1, variance_now)


        # Krigagem
        if tipo_interpolacao == "KRIGEALT":
            # Ruido para quando der erro
            ruido_0 = np.random.normal(scale=1, size=pontos_na.shape[0])
            ruido_1 = np.random.normal(scale=1, size=pontos_na.shape[0])

            pontos_na[:, 0] += ruido_0
            pontos_na[:, 1] += ruido_1
            z = pontos_na[:, 2]

            OK = UniversalKriging(
                x=pontos_na[:, 0], y=pontos_na[:, 1], z=z,
                variogram_model="linear",
            )

            for banda in range(altimetria.shape[0]):  # Para cada banda (se houver mais de uma)
                for index_y in range(altimetria.shape[1]):  # Para cada linha
                    if index_y%10 == 0:
                        print(index_y, altimetria.shape[1])

                    points_x:list[float] = []
                    points_y:list[float] = [y_coords[index_y]]*altimetria.shape[2]

                    for index_x in range(altimetria.shape[2]):  # Para cada coluna
                        px = x_coords[index_x]
                        points_x.append(px)

                    z_interp, ss = OK.execute('points', np.array([points_x]), np.array([points_y]))

                    null_raster.values[banda, index_y] = z_interp
                    variance_raster.values[banda, index_y] = ss

                variance_now = variance_raster.values[banda]
                altitude_now = null_raster.values[banda]
                altitude_max = altimetria.values[banda]

                NA_val = altitude_max - altitude_now
                var_val = variance_now

                null_raster.values[banda] = NA_val
                variance_raster.values[banda] = var_val


        # Inverso da distância
        elif tipo_interpolacao == "INVDIST":
            pontos = get_altitude_by_points(pontos_na)

            for banda in range(altimetria.shape[0]):  # Para cada banda (se houver mais de uma)
                for index_y in range(altimetria.shape[1]):  # Para cada linha
                    if index_y%10 == 0:
                        print(index_y, altimetria.shape[1])
                        
                    points_x:list[float] = []
                    points_y:list[float] = [y_coords[index_y]]*altimetria.shape[2]

                    for index_x in range(altimetria.shape[2]):  # Para cada coluna
                        px = x_coords[index_x]
                        points_x.append(px)

                    points_x = np.array(points_x)
                    points_y = np.array(points_y)

                    dif_x = points_x[:, np.newaxis]-pontos[:, 0]
                    dif_y = points_y[:, np.newaxis]-pontos[:, 1]

                    distancias = np.hypot(dif_x, dif_y)

                    multi_media = 1/distancias
                    div_media = np.sum(multi_media, axis=1)
                    na_medio = np.sum(pontos[:, 2] * multi_media, axis=1)/div_media

                    null_raster.values[banda, index_y] = na_medio
                    variance_raster.values[banda, index_y] = -1

                if limit_by_cota:
                    altitude_max = altimetria.values[banda]
                    altitude_now = null_raster.values[banda]
                    mask = altitude_max < altitude_now
                    altitude_now = np.where(mask, altitude_max, altitude_now)
                    null_raster.values[banda] = altitude_now


        # Inverso da distância
        elif tipo_interpolacao == "INVALT":
            for banda in range(altimetria.shape[0]):  # Para cada banda (se houver mais de uma)
                for index_y in range(altimetria.shape[1]):  # Para cada linha
                    if index_y%10 == 0:
                        print(index_y, altimetria.shape[1])
                        
                    points_x:list[float] = []
                    points_y:list[float] = [y_coords[index_y]]*altimetria.shape[2]

                    for index_x in range(altimetria.shape[2]):  # Para cada coluna
                        px = x_coords[index_x]
                        points_x.append(px)

                    points_x = np.array(points_x)
                    points_y = np.array(points_y)

                    dif_x = points_x[:, np.newaxis]-pontos_na[:, 0]
                    dif_y = points_y[:, np.newaxis]-pontos_na[:, 1]

                    distancias = np.hypot(dif_x, dif_y)

                    multi_media = 1/distancias
                    div_media = np.sum(multi_media, axis=1)
                    na_medio = np.sum(pontos_na[:, 2] * multi_media, axis=1)/div_media

                    altitude_val = altimetria.values[banda, index_y]

                    NA_val = altitude_val-na_medio

                    null_raster.values[banda, index_y] = NA_val
                    variance_raster.values[banda, index_y] = -1

        # Recortando o raster para a área de interesse
        null_raster = null_raster.rio.clip(area.geometry.apply(mapping), area.crs)
        comple_name = ""
        if tipo_interpolacao == "KRIGE":
            comple_name += "_kri"
        elif tipo_interpolacao == "KRIGEALT":
            comple_name += "_krialt"
        elif tipo_interpolacao == "INVDIST":
            comple_name += "_invlin"
        elif tipo_interpolacao == "INVALT":
            comple_name += "_invalt"

        if use_pocos:
            comple_name += "_pocos"
        if use_rios:
            comple_name += "_trechos"
        if use_outorgas:
            comple_name += "_outorg"

        null_raster = null_raster.rio.clip(area.geometry.apply(mapping), area.crs)
        null_raster.rio.to_raster(f"F:/Mestrado/1° Semestre/Hidrogeologia/Trabalho Prático/dados_NA/{'' if limit_by_cota else 'no_'}limit_by_cota/saida{comple_name}.tif")

        variance_raster = variance_raster.rio.clip(area.geometry.apply(mapping), area.crs)
        variance_raster.rio.to_raster(f"F:/Mestrado/1° Semestre/Hidrogeologia/Trabalho Prático/dados_NA/{'' if limit_by_cota else 'no_'}limit_by_cota/variancia{comple_name}.tif")

        gdf = generate_contornos(null_raster)
        gdf.to_file(f"F:/Mestrado/1° Semestre/Hidrogeologia/Trabalho Prático/dados_NA/{'' if limit_by_cota else 'no_'}limit_by_cota/cota{comple_name}.gpkg", driver="GPKG")



# Plotando os resultados finais

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plotar a superfície
X, Y = np.meshgrid(x_coords, y_coords)
# ax.plot_surface(X, Y, altitude, cmap='cividis')

band = null_raster.sel(band=1)
x_coords, y_coords = band['x'].values, band['y'].values
piezo = band.values
X, Y = np.meshgrid(x_coords, y_coords)
ax.plot_surface(X, Y, piezo, cmap='viridis')

# Adicionar rótulos aos eixos
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Exibir o gráfico
plt.show()

