import geopandas as gpd
gpd.options.io_engine = "pyogrio"

import numpy as np
from shapely.geometry import Point

from get_pocos import get_pocos
from get_altimetria import get_altimetria
from piezometria_pocos import get_piezometria_pocos

################ passo 1 área de análise
area = gpd.read_file(r"F:/Mestrado/1° Semestre/Hidrogeologia/Trabalho Prático/dados_ibirite.gpkg", layer="area_drenagem_sup")

################ passo 2 obtendo o registro dos poços no SIAGAS
pocos = get_pocos(
    area.unary_union,
    file_system=r"F:\Projetos\Python\PocosSIAGAS\dados_ibirite.gpkg",
    layer="pocos",
)

################ passo 3 obtendo o mapa altimétrico para a área
cotas = get_altimetria(
    area.unary_union,
    file_system=r"F:\Mestrado\1° Semestre\Hidrogeologia\Trabalho Prático\dados_ibirite.gpkg",
    layer="cotas",
)

################ passo 4 determinando a cota piezométrica para cada poço
piezo_pocos = get_piezometria_pocos(pocos=pocos, cotas=cotas)

################ passo 5 determinando a cota piezométrica para as massas d'água



