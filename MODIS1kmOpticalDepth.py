from pyhdf.SD import SD, SDC
import numpy as np
import pandas as pd
import requests

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
 
  
 
# create session with the user credentials that will be used to authenticate access to the data
username = "raghaven"
password= "Swarnalatha5244@2"
session = SessionWithHeaderRedirection(username, password)

# the url of the file we wish to retrieve
url = "https://e4ftl01.cr.usgs.gov//MODV6_Dal_H/MOTA/MCD19A2.006/2018.09.25/MCD19A2.A2018268.h10v05.006.2018277011138.hdf"

# extract the filename from the url to be used when saving the file
filename = url[url.rfind('/')+1:]  

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


hdf = SD(filename, SDC.READ)

h = int(filename[18:20])
v = int(filename[21:23])

projData = pd.read_csv("sinusoidalProjectionData.csv")

slicedData = projData[(projData.h == h) & (projData.v == v)]
lon_min = slicedData.lon_min.values
lon_max = slicedData.lon_max.values
lat_min = slicedData.lat_min.values
lat_max = slicedData.lat_max.values

#### Optical_depth_055

optical_depth_055 = hdf.select('Optical_Depth_055')
OD_fillValue  = optical_depth_055.getfillvalue()
OD_scaleFactor = optical_depth_055.attributes(full=1)["scale_factor"][0]
optical_depth = optical_depth_055[0,:,:].ravel()

long = np.tile(np.arange(lon_min,lon_max,((lon_max-lon_min)/1200)),1200)
lat = np.repeat(np.arange(lat_max,lat_min,((lat_min-lat_max)/1200)),1200)
data = np.column_stack([lat, long,optical_depth])

data1 = pd.DataFrame(data=data,columns=["Lat","Long","optical_depth"])
data2 = data1[data1.optical_depth != OD_fillValue]

data2.to_csv('MODIS_1KMData.csv')


