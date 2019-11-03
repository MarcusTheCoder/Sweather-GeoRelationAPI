''' Have Pride  in your work!!!! :D
'''
from osgeo import gdal

from gdal import ogr
from gdal import osr

import math
import shutil
import os
import requests
import json
import sys
import numpy
OUTCONTOUR = "elevCont"
OUTPOLY = "elevPOLY"
OUTTIFF = "filteredRaster.tiff"
CENTROID = 0
EDGE_CENTROID = 1
MIN_RADIUS = 1000
STEP = 1000
TIFF_STEP = 50 # How many meters before adding a new contour line
TIFF_DIV = 10
MAX_RADIUS = 8047
def unflatten(arr,width,height):
  out2d = []
  for y in range(0,height):
    out2d.append([])
    for x in range(0,width):
      index = y * width + x
      out2d[y].append(arr[index])
  return out2d

def polyganizeRaster(rasterPath,outPoly=OUTPOLY,outRaster=OUTTIFF,divv=TIFF_DIV):
  '''
  Will take input as .tiff file. Outputs a vectorized raster,
  based on steps of 100 m in elevation
  '''
  raster = gdal.Open(rasterPath)
  in1 = raster.GetRasterBand(1)
  if os.path.exists(outPoly):
    shutil.rmtree(outPoly)
  
  rarr = in1.ReadAsArray()
  flarr = numpy.ndarray.flatten(rarr)
  outArr = numpy.copy(flarr)
  flarr.sort()
   
  minVal = int(flarr[0])
  maxVal = int(flarr[len(flarr) - 1] + 1)
  stepDiff = (maxVal - minVal) / divv

  for i in range(0,len(outArr)):
    s = stepDiff
    j = 0
    while j <= divv:
      if outArr[i] > s:
        s += stepDiff
        j += 1
      else:
        break
    outArr[i] = j


  # Export the altered array as raster
  x_pixels = len(rarr[0])
  y_pixels = len(rarr)
  driver = gdal.GetDriverByName("GTiff")
  dataset = driver.Create(outRaster,x_pixels,y_pixels,1,gdal.GDT_Float32)
  o2 = numpy.array(unflatten(outArr,x_pixels,y_pixels))
  #dataset.WriteArray(o2)
  dataset.GetRasterBand(1).WriteArray(o2)
  
  geoTrans = raster.GetGeoTransform()
  proj = raster.GetProjection()
  dataset.SetGeoTransform(geoTrans)
  dataset.SetProjection(proj)
  dataset.FlushCache()  
  return stepDiff
  

def processRaster(rasterPath,sep,outPolyFp=OUTPOLY,giveEntries = False,):
  raster = gdal.Open(rasterPath)
  in1 = raster.GetRasterBand(1)

  if os.path.exists(outPolyFp):
    shutil.rmtree(outPolyFp)

  ogr_ds = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(outPolyFp)
  polygLayer = ogr_ds.CreateLayer('polyg')

  field_defn = ogr.FieldDefn("ID", ogr.OFTInteger)
  polygLayer.CreateField(field_defn)
  field_defn = ogr.FieldDefn("ELEV", ogr.OFTReal)
  polygLayer.CreateField(field_defn)
  polyID_index = polygLayer.GetLayerDefn().GetFieldIndex("ID")
  polyElev_index = polygLayer.GetLayerDefn().GetFieldIndex("ELEV")
  entries = []
  gdal.Polygonize(in1,None,polygLayer,0,[],callback=None)
  # Marcus 
  if giveEntries:
    num = 0
    for p in polygLayer:      
      elev = float(p.GetFieldAsString(polyID_index)) * sep
      p.SetField("ELEV",str(elev))
      p.SetField("ID",num)      
      num += 1

    return entries





# Panozzo
def findSafePolys(x,y,polyFp,anticipatedWaterLevel,rasterGrad,minRadius=MIN_RADIUS,maxRadius=MAX_RADIUS,numToFind=1,shapeMode=CENTROID):
  polyData = ogr.Open(polyFp)
  polyLayer = polyData.GetLayer()

  iPoint = ogr.Geometry(ogr.wkbPoint)
  iPoint.SetPoint_2D(0,x,y)
  circBufferMAX = iPoint.Buffer(maxRadius)

  
  dist = minRadius
  
  foundCount = 0
  polyID_index = polyLayer.GetLayerDefn().GetFieldIndex("ID")
  polyElev_index = polyLayer.GetLayerDefn().GetFieldIndex("ELEV")
  goodEntries = {} # Will store polygon,(point,deltaSafety)
  ucHeight = -1
  while dist <= maxRadius and foundCount < numToFind:
    buff = iPoint.Buffer(dist)
    #print(buff.Centroid().GetX(),buff.Centroid().GetY())
    for poly in polyLayer:
      if poly.GetGeometryRef().Intersects(iPoint):
        elev = float(poly.GetFieldAsString(0)) * rasterGrad
        ucHeight = elev
      
      elif buff.Intersects(poly.GetGeometryRef()):
        elev = float(poly.GetFieldAsString(0)) * rasterGrad
        if elev > anticipatedWaterLevel:
          dt = elev - anticipatedWaterLevel
          if shapeMode == CENTROID:
            interPoint = poly.GetGeometryRef().Centroid()
          elif shapeMode == EDGE_CENTROID:
            interPoint = poly.GetGeometryRef().Intersection(buff).Centroid()
          goodEntries[poly] = (interPoint,dt)
          foundCount += 1
    polyLayer.ResetReading()
    if len(goodEntries.keys()) < numToFind:
      dist += STEP
  

  return [goodEntries,ucHeight]


def createTiffFromCoord(long,lat,expanse=[MAX_RADIUS,MAX_RADIUS],res=3):
  pass
def calculateRisk(anticipWater,myElev):
  if myElev <= anticipWater:
    return 5
  elif myElev - 1 <= anticipWater:
    return 4
  elif myElev - 2 <= anticipWater:
    return 3
  elif myElev - 4 <= anticipWater:
    return 2
  elif myElev - 8 <= anticipWater:
    return 1
  else :
    return 0

def sweather_main(long,lat,minR=MIN_RADIUS,maxR=MAX_RADIUS,nToSolve=1,anticipWater=20,threshDiv=TIFF_DIV):
  '''
  Returns a list of coordinates that are safer to go to, in ranking of height above the water
  Also returns an estimation of the current elevation of the sector that the user click was
  on (more accurate is returned earlier)
  '''
  anticipatedRise = 10 # in meters
  filepath = "./data/AlbanySimpleRaster.tif"
  step = polyganizeRaster(filepath,divv=threshDiv)  
  processRaster(OUTTIFF,step,OUTPOLY,True)

  x = long # Longitude
  y = lat # Latitude
  oRef = osr.SpatialReference()
  oRef.ImportFromEPSG(4326)
  targRef = osr.SpatialReference()
  targRef.ImportFromEPSG(3857)
  ctran = osr.CoordinateTransformation(oRef,targRef)
  revTran = osr.CoordinateTransformation(targRef,oRef)
  [longi,lati,z] = ctran.TransformPoint(x,y)

  [goodPolys, ucHeight] = findSafePolys(longi,lati,OUTPOLY,anticipWater,step,shapeMode=EDGE_CENTROID,minRadius=minR,maxRadius=maxR,numToFind=nToSolve)

  l = []
  goodPoints = [] # point tuple,diff
  

  for poly, tup in goodPolys.items():
    l.append((poly,tup[0],tup[1]))
  l.sort(key=lambda ele: ele[2],reverse=True)
  for ele in l:
    x = ele[1].GetX()
    y = ele[1].GetY()
    distAway = math.sqrt((math.pow(x,2)-math.pow(longi,2))/(math.pow(y,2)-math.pow(lati,2)))
    [longi2,lati2,z] = revTran.TransformPoint(x,y)
    
    goodPoints.append([[longi2,lati2],ele[2],distAway,calculateRisk(anticipWater,ele[2])])
  goodPoints.sort(key=lambda ls: ls[1]/ls[2],reverse=True)

  return [goodPoints,ucHeight]




