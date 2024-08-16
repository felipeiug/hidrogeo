import geopandas as gpd
import numpy as np
gpd.options.io_engine = "pyogrio"
import rioxarray
import xarray as xr
from shapely.geometry import mapping, box, LineString
import matplotlib.pyplot as plt
from skimage import measure
import rasterio
import matplotlib.pyplot as plt

from pykrige.ok import OrdinaryKriging
from pykrige.uk import UniversalKriging

import numpy as np


crs = "EPSG:31983"

# Lendo os arquivos
altimetria = rioxarray.open_rasterio(r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_NA\limit_by_cota\saida_invlin_pocos_trechos_outorg.tif")
altimetria = altimetria.rio.reproject(crs)

# Extraindo os valores
band = altimetria.sel(band=1)
x_coords, y_coords = band['x'].values, band['y'].values
altitude = band.values

# Supondo que 'altitude' seja o seu raster de altimetria carregado como um array numpy
grad_y, grad_x = np.gradient(altitude)

magnitude = np.sqrt(grad_x**2 + grad_y**2)
direction_x = grad_x / magnitude
direction_y = grad_y / magnitude

plt.figure(figsize=(10, 8))
plt.imshow(altitude, cmap='terrain', origin='lower')

# Definindo o espaçamento das setas para evitar sobrecarga visual
spacing = 20
plt.quiver(
    np.arange(0, altitude.shape[1], spacing),  # Posições x
    np.arange(0, altitude.shape[0], spacing),  # Posições y
    direction_x[::spacing, ::spacing],          # Componentes x das setas
    direction_y[::spacing, ::spacing],          # Componentes y das setas
    scale=50,                                   # Escala das setas
    color='blue'
)

plt.colorbar(label='Elevation (m)')
plt.title('Escoamento do Terreno')
plt.show()