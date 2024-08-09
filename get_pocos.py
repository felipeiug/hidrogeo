import geopandas as gpd
from shapely import Polygon

def get_pocos(area: Polygon, file_system:str=None, layer:str=None)-> gpd.GeoDataFrame:
    bbox = area.bounds

    if file_system:
        pocos:gpd.GeoDataFrame = gpd.read_file(file_system, layer=layer)

    else:
        # TODO: GET POÇOS NO SIAGAS, e clip dos dados
        raise NotImplementedError("Obtenção dos poços no siagas ainda não implementada")
    
    pocos.to_crs("EPSG:31983", inplace=True)
    pocos = pocos.clip_by_rect(*bbox)

    return pocos