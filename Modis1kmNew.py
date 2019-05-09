from pyhdf.SD import SD, SDC
import numpy as np
import pandas as pd
import requests
from time import sleep
import re
from math import cos,degrees
from pathlib import Path


# *********************************** MODIS data ***********************************************

#******************************** Earthdata Authorization **************************************
# overriding requests.Session.rebuild_auth to mantain headers when redirected
class SessionWithHeaderRedirection(requests.Session):
    AUTH_HOST = 'urs.earthdata.nasa.gov'
    def __init__(self, username, password):
        super().__init__()
        self.auth = (username, password)

   # Overrides from the library to keep headers when redirected to or from
 
   # the NASA auth host.
    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url

        if 'Authorization' in headers:
            original_parsed = requests.utils.urlparse(response.request.url)
            redirect_parsed = requests.utils.urlparse(url)
            
            if (original_parsed.hostname != redirect_parsed.hostname) and redirect_parsed.hostname != self.AUTH_HOST and original_parsed.hostname != self.AUTH_HOST:
                del headers['Authorization']
        return


def downloadHDFFile(username,password,filename,url):
    
    # create session with the user credentials that will be used to authenticate access to the data
    session = SessionWithHeaderRedirection(username, password)
    
    try:
        # submit the request using the session
        response = session.get(url, stream=True)
        print(response.status_code)

        # raise an exception in case of http errors
        response.raise_for_status()  

        # save the file
        with open(filename, 'wb') as fd:
            for chunk in response.iter_content(chunk_size=1024*1024):
                fd.write(chunk)

    except requests.exceptions.HTTPError as e:
        print(e)
  
 

#*************************************** Processing HDF file ******************************************

def processHDF(filename):
    # Read HDF File 
    hdf = SD(filename, SDC.READ)
    
    # Select Optical Depth .55 micron wavelength
    optical_depth_055 = hdf.select('Optical_Depth_055')
    shape = optical_depth_055[:,:,:].shape

    # Obtain fillValue and Scale factor for optical depth
    OD_fillValue  = optical_depth_055.getfillvalue()
    OD_scaleFactor = optical_depth_055.attributes(full=1)["scale_factor"][0]

    # Flatten and merge data based on number of orbits
    optical_depth = np.array([], dtype=np.int16).reshape(0,1)
    for i in range(shape[0]):
        optical_depth = np.append(optical_depth,optical_depth_055[i,:,:].ravel())

    # Scale data based on the scale factor.
    optical_depth = optical_depth*OD_scaleFactor





    # Select water vopur readings
    Column_WV_n = hdf.select('Column_WV')

    # Obtain fillValue and Scale factor for optical depth
    WV_fillValue  = Column_WV_n.getfillvalue()
    WV_scaleFactor = Column_WV_n.attributes(full=1)["scale_factor"][0]

    # Flatten and merge data based on number of orbits
    Column_WV = np.array([], dtype=np.int16).reshape(0,1)
    for i in range(shape[0]):
        Column_WV = np.append(Column_WV,Column_WV_n[i,:,:].ravel())
    
    # Scale data based on the scale factor.
    Column_WV = Column_WV*WV_scaleFactor





    # Select water vopur readings
    AOD_Uncertainty_n = hdf.select('AOD_Uncertainty')

    # Obtain fillValue and Scale factor for optical depth
    UN_fillValue  = AOD_Uncertainty_n.getfillvalue()
    UN_scaleFactor = AOD_Uncertainty_n.attributes(full=1)["scale_factor"][0]

    # Flatten and merge data based on number of orbits
    AOD_Uncertainty = np.array([], dtype=np.int16).reshape(0,1)
    for i in range(shape[0]):
        AOD_Uncertainty = np.append(AOD_Uncertainty,AOD_Uncertainty_n[i,:,:].ravel())

    # Scale data based on the scale factor.
    AOD_Uncertainty = AOD_Uncertainty*UN_scaleFactor




    # Ravel and obtain longitude and latitude values as numpy array
    long = np.tile(np.arange(lon_min,lon_max,((lon_max-lon_min)/shape[1])),shape[1])
    lat = np.repeat(np.arange(lat_max,lat_min,((lat_min-lat_max)/shape[1])),shape[1])

    # Create a column stack numpy array of latitude
    data = np.column_stack([np.tile(lat,shape[0]), np.tile(long,shape[0]),optical_depth,Column_WV,AOD_Uncertainty])    
    
    # Create dataframe from numpy column stack
    data1 = pd.DataFrame(data=data,columns=["Lat","Long","optical_depth","Column_WV","AOD_Uncertainty"])
    
    # Remove rows with fill values
    data2 = data1[(data1.optical_depth != OD_fillValue) & (data1.Column_WV != WV_fillValue) & (data1.AOD_Uncertainty != UN_fillValue)]
    
    # Average values with similar Latitude and Longitudes
    data3 = data2.groupby(['Lat','Long'])[['optical_depth','Column_WV','AOD_Uncertainty']].mean().reset_index()
    return data3,hdf,AOD_Uncertainty,Column_WV, optical_depth





#******************************************** EPA data ********************************************

def downloadEPAdata(epaUser,epaPwd,url):
    
    # Extract date values from URL
    temp = url.split('/')[7].split('.')
    epaDate = temp[0]+temp[1]+temp[2]
    
    # Build EPA API URL using DMCSV to obtain the transaction ID
    epaUrl = "https://aqs.epa.gov/api/rawDataNotify?user="+epaUser+"&pw="+epaPwd+"&format=DMCSV&param=88101&bdate="+epaDate+"&edate="+epaDate+"&minlat="+str(lat_min[0])+"&maxlat="+str(lat_max[0])+"&minlon="+str(lon_min[0])+"&maxlon="+str(lon_max[0])+"&frmonly=Y"
    transactionID = requests.get(epaUrl).text
    
    # Build data retreival URL
    retrieveUrl = "https://aqs.epa.gov/api/retrieve?id="+transactionID
    x=0
    
    # Wait until the REST request is processed by the API and data is returned
    while(x==0):
        try:
            epaData = pd.read_csv(retrieveUrl)
            x=1
        except:
            sleep(10)
    
    # Filter the EPA dataset for the required date and Sample duration of 1 Hour
    epaData = epaData[(epaData['Sample Duration'] == "1 HOUR") & (epaData['Date GMT'] == (temp[0]+"-"+temp[1]+"-"+temp[2]))]
    return epaData




#******************************************* Integration **************************************

def timeMerger(hdf, epaData):
    
    # Extract the Core Metadata details of the retrieved HDF Files
    CoreMetadata = hdf.attributes('CoreMetadata.0')['CoreMetadata.0'][0]
    
    # Obtain the start and end time of data captured timestamp based on the field "RANGEBEGINNINGTIME" and "RANGEENDINGTIME"
    startTime = re.findall('RANGEBEGINNINGTIME.*\"(.*?)\".*RANGEBEGINNINGTIME', CoreMetadata, re.DOTALL)[0]
    endTime = re.findall('RANGEENDINGTIME.*\"(.*?)\".*RANGEENDINGTIME', CoreMetadata, re.DOTALL)[0]

    # Obtain the start hour and minute, End Hour and minute
    startHour = int(startTime.split(":")[0])
    startMin = int(startTime.split(":")[1])

    endHour = int(endTime.split(":")[0])
    endMin = int(endTime.split(":")[1])

    # If the start minute or end minute > 30 then increment the start and end Hour time respectively
    if startMin >= 30:
        if startHour != 23:
            startHour+=1
    if endMin >=30:
        if endHour != 23:
            endHour+=1

    # Filter the data with GMT hour based on the start hour and end hour 
    HrGMT = epaData['24 Hour GMT'].str[:2].astype(int)    
    epaDataNew = epaData[(HrGMT >= startHour) & (HrGMT <= endHour)]
    epaDataNew = epaDataNew.reset_index(drop=True)
    epaDataNew.Latitude = epaDataNew.Latitude.astype(float)
    epaDataNew.Longitude = epaDataNew.Longitude.astype(float)
    
    return epaDataNew


#******************************************* Distance Function **************************************


def distCalc(x,y):
    # The distance between pairs of latitude and longitude is calculated using "Havershine formula"
    # Set the earth radius value in Km
    R = 6373.0
    
    # Havershine Formula
    dlat=y[['Lat']]-x[1][0]
    dlon=y[['Long']]-x[1][1]
    
    a = np.add(np.power(np.sin(dlat/2),2),cos(x[1][0])*np.multiply(np.cos(y[['Lat']]),np.power(np.sin(dlon/2),2)))
    c = 2 * np.arctan2(np.sqrt(a),np.sqrt(1-a))
    
    dis = R * c
    
    # Obtain and return indices where the distance is less than 2 kms
    return np.where(dis<=2)[0].tolist()
    
def spaceMerger(epaDataNew,data3):
    
    # Convert Lat and Longitude from degree to radians
    b=np.radians(epaDataNew[['Latitude','Longitude']].drop_duplicates().astype(float))
    c=np.radians(data3[['Lat','Long']])
    
    # Create an empty dataframe to fill with space merged EPA data
    epaDataNew2 = pd.DataFrame(columns=epaDataNew.columns.values)
    index = []

    # Iterate through each unique latitude and Longitude (i.e Ground Station Locations) and calculate distance <= 2 kms
    for each in b.iterrows():
        
        # Ontain indices from distCalc function
        val = distCalc(each,c)
        
        # if there are no ground stations near the pixels from MODIS then continue further else filter the EPA data based on space 
        if len(val) == 0:
            continue
        else:
            index+=val
            epaDataNew2 = epaDataNew2.append(epaDataNew[(epaDataNew.Latitude == degrees(each[1][0])) & (epaDataNew.Longitude == degrees(each[1][1]))])
            
    return epaDataNew2,index
        

#*********************************** Prime Function (Sequential Execution) *******************************************
    
# the url of the file we wish to retrieve
url = "https://e4ftl01.cr.usgs.gov//MODV6_Dal_H/MOTA/MCD19A2.006/2018.09.25/MCD19A2.A2018268.h10v05.006.2018277011138.hdf"

# extract the filename from the url to be used when saving the file
filename = url[url.rfind('/')+1:]  

# Enter MODIS username and Password
username = ""
password = ""

# Enter EPA username and Password
epaUser = ""
epaPwd = ""

if not Path(filename).exists():
    downloadHDFFile(username,password,filename,url)
    
    
# Get "h" and "v" values from filename, see file documentation for further details
h = int(filename[18:20])
v = int(filename[21:23])

# Read table for conversion from sinusoidal projection to normal projection
projData = pd.read_csv("sinusoidalProjectionData.csv")

slicedData = projData[(projData.h == h) & (projData.v == v)]

# Find the latitude and longitude boundary values
lon_min = slicedData.lon_min.values
lon_max = slicedData.lon_max.values
lat_min = slicedData.lat_min.values
lat_max = slicedData.lat_max.values

# Call each function sequntially

# Process HDF file
data3, hdf,AOD_Uncertainty,Column_WV, optical_depth = processHDF(filename)

# Download EPA data
epaData = downloadEPAdata(epaUser,epaPwd,url)

# Merge EPA and MODIS data in time
epaDataNew = timeMerger(hdf,epaData)

# Merge EPA and MODIS data in space and final EPA dataset
epaDataNew2,index = spaceMerger(epaDataNew,data3)

# Obtain the final MODIS dataset
modisData = data3.iloc[index,:]

modisData.to_csv("modisData.csv")
epaDataNew2.to_csv("epaData.csv")

    






