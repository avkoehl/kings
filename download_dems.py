
huc12s = [
'180101070402',
'180101070401',
'180101070201',
'180101070204',
'180101070207'
]

huc12wbd = WBD(layer='huc12')

for huc12 in huc12s:
    # get geom
    boundary = huc12wbd.byids('huc12', huc12)
    geom = boundary.loc[0, 'geometry']

    # get dem

    # watershed

    # split into subbasins of approx equal area?

    # get 1m dem for those geoms

    # save

