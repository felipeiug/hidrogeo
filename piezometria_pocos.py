import geopandas as gpd
gpd.options.io_engine = "pyogrio"

def get_piezometria_pocos(
    pocos:gpd.GeoDataFrame,
    cotas:gpd.GeoDataFrame,
    niv_estatico_col:str="nivel_estatico",
    elev_col:str="ELEV",
):
    pocos.to_crs("EPSG:31983", inplace=True)
    cotas.to_crs("EPSG:31983", inplace=True)

    new_pocos = {
        "ponto":[],
        "cota_piezo":[],
        "geometry":[],
    }

    for poco in pocos.itertuples():
        geom = poco.geometry
        is_nascente = poco.natureza == "Nascente"

        if is_nascente:
            cota_piezo = 0
        else:
            try:
                nv_estatico:gpd.GeoSeries = getattr(poco, niv_estatico_col)
                cota_piezo = float(nv_estatico.replace(",", "."))
            except:
                cota_piezo = None

        if cota_piezo is None:
            continue

        distances = cotas.distance(geom)
        min_dist_index = distances.idxmin()

        min_dist_cota_data = cotas.loc[min_dist_index]
        elevacao = min_dist_cota_data["ELEV"]

        cota_piezo = elevacao - cota_piezo

        new_pocos["ponto"].append(poco.ponto)
        new_pocos["cota_piezo"].append(cota_piezo)
        new_pocos["geometry"].append(geom)

    return gpd.GeoDataFrame(new_pocos, geometry="geometry", crs = "EPSG:31983")


