import arcpy
from arcpy import env
from arcpy.sa import *
import os
import os.path
import sys
import string  
  
dir = 'F:\\Training\\tif_sample\\tif_sample\\wet_land'  
  
# Import arcpy module  
import arcpy  
  
files = os.listdir(dir)  
for f in files:  
    if os.path.splitext(f)[1] == '.TIF':  
        # Script arguments...  
        Input_raster_file = dir + os.sep + f  
  
        # Local variables...
        Output_data_type = "FLOAT"  
        Raster_Format = "ASCII"  
        Output_Workspace = "F:\\Training\\tif_sample\\txt_sample\\wet_land"  
  
        # =============== file name process ======================  
        basename = os.path.splitext(f)[0];  
        Output_ASCII = Output_Workspace + os.sep + basename + ".txt" 
  
        if os.path.exists(Output_ASCII) == False:  
            print Input_raster_file  
            # Process: Raster To Other Format (multiple)...  
            arcpy.RasterToASCII_conversion(Input_raster_file,   
                        Output_ASCII )  
  
            print Output_ASCII  
