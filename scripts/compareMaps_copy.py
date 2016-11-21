'''
Compare Maps.
Clips maps to matching extent, generates a difference image, & generates a scatterplot.

Usage:
  compareMaps.py <mappath1> <mappath2> <outputdir> [--map1_band=<b1>] [--map2_band=<b2>] [--map1_scale=<s1>] [--map2_scale=<s2>] [--boundarymap=<bm>] [--meta=<m>]
  compareMaps.py -h | --help

Options:
  -h --help     	  Show this screen.
  --map1_band=<b1> 	  Band of map #1 [default: 1].
  --map2_band=<b2>    Band of map #2 [default: 1].
  --map1_scale=<s1>	  Multiply each pixel in map1 by a scale factor [default: 1].
  --map2_scale=<s2>   Multiply each pixel in map2 by a scale factor [default: 1].
  --boundarymap=<bm>  Map to base boundaries on [default: 1].
  --meta=<meta>  	  Additional notes for meta.txt files.
'''
import docopt, gdal, os, subprocess
from lthacks import *
from lthacks.intersectMask import *
from gdalconst import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import ogr
import time

#import pylab


def clipMap(srcMap, mskMap, srcBand, mskBand, outputDir):

	#define output path for clipped map
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
		print "\nNew directory created: " + outputDir

	srcMapName = os.path.splitext(os.path.basename(srcMap))[0]
	srcMapExt = os.path.splitext(os.path.basename(srcMap))[1]
	maskMapName = os.path.splitext(os.path.basename(mskMap))[0]

	clippedMapPath = os.path.join(outputDir, srcMapName + "_clippedto_" + maskMapName + srcMapExt)
	
	if not os.path.exists(clippedMapPath):

		cmdArgs = [srcMap, mskMap, clippedMapPath, srcBand, mskBand]
		cmdArgsSpaces = [str(i) + " " for i in cmdArgs]
		cmd = "intersectMask {0} {1} {2} --src_band={3} --msk_band={4}".format(*cmdArgsSpaces)
	
		#call intersectMask
		print "\nClipping '" + srcMap + "' to match extent of '" + mskMap + "' ..."
		print "\n"+cmd
		process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
		process.wait()
		print process.returncode
		
	else:
		print "\nClipped map already exists: " + clippedMapPath
	
	return clippedMapPath


def get_dif_map(ar1, ar2, nodata1, nodata2):
    
    mask1 = ar1 == nodata1
    mask2 = ar2 == nodata2
    nans = np.logical_or(mask1, mask2)
    
    dif = ar1 - ar2
    dif[nans] = nodata1
    
    return dif, nans


def get_overlapping_polys(src_shp, ovr_shp, out_shp):
    '''
    Return a shapefile of all features in ds_src that touch ds_ovr
    '''
    ds_src = ogr.Open(src_shp)
    if ds_src == None:
        print 'Shapefile does not exist or is not valid:\n%s' % src_shp
        return None
    lyr_src = ds_src.GetLayer()
    srs_src = lyr_src.GetSpatialRef()
    
    ds_ovr = ogr.Open(ovr_shp)
    if ds_ovr == None:
        print 'Shapefile does not exist or is not valid:\n%s' % ovr_shp
        return None
    lyr_ovr = ds_ovr.GetLayer()
    
    # Create the output dataset
    driver = ds_src.GetDriver()
    if os.path.exists(out_shp):
        os.remove(out_shp)   
    ds_out = driver.CreateDataSource(out_shp)
    lyr_out = ds_out.CreateLayer(os.path.basename(out_shp)[:-4], srs_src, geom_type=ogr.wkbMultiPolygon) 
    lyr_out_def = lyr_out.GetLayerDefn()
    
    #Get field definitions
    lyr_src_def = lyr_src.GetLayerDefn()
    for i in range(lyr_src_def.GetFieldCount()):
        field_def = lyr_src_def.GetFieldDefn(i)
        lyr_out.CreateField(field_def) 
    
    print 'Finding overlapping features...'
    t0 = time.time()
    # Loop through each feautre and check for overlap
    #for i in xrange(lyr_src.GetFeatureCount()):
    feat_src = lyr_src.GetNextFeature()
    while feat_src:
        t1 = time.time()
        geom_src = feat_src.GetGeometryRef()
        
        for j in xrange(lyr_ovr.GetFeatureCount()):
            feat_ovr = lyr_ovr.GetFeature(j)
            geom_ovr = feat_ovr.GetGeometryRef()
            
            # If there's any overlap, add the feature to the lyr_out
            if geom_ovr.Intersect(geom_src):
                feat_out = ogr.Feature(lyr_out_def)
                feat_out.SetGeometry(geom_src)
                # Get the fields from the source
                for i in range(lyr_out_def.GetFieldCount()):
                    feat_out.SetField(lyr_out_def.GetFieldDefn(i).GetNameRef(), feat_src.GetField(i))
                lyr_out.CreateFeature(feat_out)
                feat_out.Destroy()
                feat_ovr.Destroy()
                feat_ovr = lyr_ovr.GetNextFeature()
                break
            else:
                feat_ovr.Destroy()
                feat_ovr = lyr_ovr.GetNextFeature()
                
        feat_src.Destroy()
        feat_src = lyr_src.GetNextFeature()
    
    print 'Total time: %.1f\n' % (time.time() - t0)
    
    ds_src.Destroy()
    ds_ovr.Destroy()
    lyr_src = None
    lyr_ovr = None
    
    print 'Shapefile written to: ', out_shp     


def get_overlapping_pixels(ar, tx, feat):
    '''
    Return an array of just pixels overlapping a feature from a vector ds
    - get bounding coords of feature
    - calc offset from array
    - get xsize/ysize of feature
    - rasterize the feature
        -last 2 could be separate function since this will always be the same
    - subset the array
    - mask out array
    '''
    geom = feat.GetGeometryRef()
    wkt = geom.ExportToWkt()
                pts_list = wkt.replace('POLYGON','').replace('LINEARRING','').replace('(','').replace(')','').strip().split(',')
            x = [float(p.split(' ')[0]) for p in pts_list]
            y = [float(p.split(' ')[1]) for p in pts_list]
    
def zonal_stats(ar):
    '''
    -for each feature, get the overlapping pixels
    -mask out nodata
    -
    '''

def main(mapPath1, mapPath2, outputDir, band1, band2, scale1, scale2, boundaryMap, metaComment):

	## check projection info & ensure it is matching

	#open maps & read projections
	projections = []
	for i in [mapPath1, mapPath2]:
		if not os.path.exists(i):
			sys.exit("\nMap path does not exist: '", mapPath1, "'")
		else:
			ds = gdal.Open(i, GA_ReadOnly)
			projections.append(ds.GetProjection())

# 	if not projections[0] == projections[1]:
# 		print projections[0]
# 		print projections[1]
# 		msg = "\nMaps are not in the same projection. \
# 		   	   Please reproject before running compareMaps."
# 		sys.exit(msg)

	## clip maps to matching extents

	#determine source & mask maps
	if boundaryMap == 1:
		mskMap = mapPath1
		srcMap = mapPath2

		mskBand = band1
		srcBand = band2

	elif boundaryMap == 2:
		mskMap = mapPath2
		srcMap = mapPath1

		mskBand = band2
		srcBand = band1
		
	else:
		sys.exit("boundarymap argument can only be 1 or 2.")

	clippedMapPath = clipMap(srcMap, mskMap, srcBand, mskBand, outputDir)

	#define new maps
	if boundaryMap == 1:
		mapPath2 = clippedMapPath
	else:
		mapPath1 = clippedMapPath
		
		
	#check that maps are now of matching extents
	print "\nChecking that map extents line up..."
	
	map1 = gdal.Open(mapPath1, GA_ReadOnly)
	map1Band = map1.GetRasterBand(band1)
	map1Data = map1Band.ReadAsArray()
	
	map2 = gdal.Open(mapPath2, GA_ReadOnly)
	map2Band = map2.GetRasterBand(band2)
	map2Data = map2Band.ReadAsArray()
	
	if not map1Data.shape == map2Data.shape:
		
		if map1Data.size > map2Data.size:
			
			srcMap = mapPath1
			mskMap = mapPath2
			
			srcBand = band1
			mskBand = band2
		
		else:
			
			srcMap = mapPath2
			mskMap = mapPath1
			
			srcBand = band2
			mskBand = band1
		
		clippedMapPath = clipMap(srcMap, mskMap, srcBand, mskBand, outputDir)
		
		#define new maps
		if mskMap == mapPath1:
			mapPath2 = clippedMapPath
		else:
			mapPath1 = clippedMapPath
			
		del map1, map1Band, map1Data, map2, map2Band, map2Data
		
	else:
		
		print "\nMap extents now match!"

	## generate difference image
	
	#define difference map output path
	map1Name = os.path.splitext(os.path.basename(mapPath1))[0]
	map2Name = os.path.splitext(os.path.basename(mapPath2))[0]

	diffMapPath = os.path.join(outputDir, map1Name + "_minus_" + map2Name + ".bsq")

	map1 = gdal.Open(mapPath1, GA_ReadOnly)
	map1Band = map1.GetRasterBand(band1)
	map1Data = map1Band.ReadAsArray()

	map2 = gdal.Open(mapPath2, GA_ReadOnly)
	map2Band = map2.GetRasterBand(band2)
	map2Data = map2Band.ReadAsArray()
	
	#check for nan values
	map1nans = np.isnan(map1Data)
	map2nans = np.isnan(map2Data)
	nans = np.logical_or(map1nans, map2nans)
		
	if not os.path.exists(diffMapPath):

		#get difference of arrays
		diffData = (map1Data * scale1) - (map2Data * scale2)

		#get map info
		transform = map1.GetGeoTransform()
		#driver = map1.GetDriver()
		driver = gdal.GetDriverByName("ENVI")
		dt = map1Band.DataType	

		#save difference map w/ metadata
		saveArrayAsRaster(diffData, transform, projections[0], driver, diffMapPath, dt)
		desc = "Difference map of " + map1Name + " and " + map2Name + "."
		createMetadata(sys.argv, diffMapPath, description=desc)
		
	else:
		
		print "\nDifference map already exists: " + diffMapPath


	## generate a scatterplot or 2d histogram

	print "\nGenerating a 2-D histogram..."
	print mapPath1 + " vs " + mapPath2 + "\n"
	print np.sum(np.isnan(map1Data)), np.sum(np.isnan(map2Data))
	
	#flatten array data
	x = map2Data[~nans].flatten() * scale2 
	y = map1Data[~nans].flatten() * scale1

	#plot x vs. y
	#plt.scatter(x,y)
	
	plt.hist2d(x,y,bins=50, norm=LogNorm())
	plt.xlabel(os.path.basename(mapPath2))
	plt.ylabel(os.path.basename(mapPath1))
	
	#plt.show()
	
# 	scatter_outfile = os.path.splitext(os.path.basename(mapPath1))[0] + "_vs_" + \
# 					  os.path.splitext(os.path.basename(mapPath2))[0] + "_scatter.png"
	scatter_outfile = os.path.splitext(os.path.basename(mapPath1))[0] + "_vs_" + \
					  os.path.splitext(os.path.basename(mapPath2))[0] + "_2dhistogram.png"
	scatter_outpath = os.path.join(outputDir, scatter_outfile)
	
	#save("signal", ext="svg", close=True, verbose=True)
	plt.savefig(scatter_outpath)

"""if __name__ == '__main__':

	try:
		#parse arguments, use file docstring as parameter definition
		args = docopt.docopt(__doc__)

		#call main function
		main(args['<mappath1>'], args['<mappath2>'], args['<outputdir>'], int(args['--map1_band']), 
			int(args['--map2_band']), float(args['--map1_scale']), float(args['--map2_scale']),
			int(args['--boundarymap']),args['--meta'])
			
	#handle invalid options
	except docopt.DocoptExit as e:
		print e.message"""

src_shp = '/vol/v2/stem/extent_shp/hex20000.shp'
tch_shp = '/vol/v2/stem/extent_shp/orwaca_singlepart.shp'
out_shp = '/vol/v2/stem/extent_shp/hex20000_CAORWA.shp'

get_overlapping_polys(src_shp, tch_shp, out_shp)