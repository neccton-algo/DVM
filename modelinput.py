import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4
from netCDF4 import Dataset as NetCDFFile
import datetime
import warnings
import math
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from dateutil.relativedelta import relativedelta
from scipy.interpolate import griddata
import pickle
import os
from matplotlib import gridspec
warnings.filterwarnings('ignore')

'''
Some use cases:
nc = modelinput.read('NECCTON_BATS_COCCO')
modelinput.compare_pmesh("NECCTON_BATS","NECCTON_BATS_COCCO","sil",dbot=150,cmin=0,cmax=1)
modelinput.pmeshmodel(filename2, variable, first=first, last=last, dtop=dtop, dbot=dbot, cmax=cmax, cmin=cmin)
modelinput.vs_occci_chl("NECCTON_BATS_COCCO")
'''
def read(prefix, *args, **kwargs):
   # EXPORT gotmout='PATH_TO_MODEL_OUTPUT_FOLDER'
   # You can do this as a default in your bashrc file
   #
   # As an alternative, if you want to interactively plot from a different folder,
   # in python, you can change the gotmout environment variable by:
   # import os
   # os.environ["gotmout"] = 'PATH_TO_MODEL_OUTPUT_FOLDER'
   path = kwargs.get('path',None)
   if path is not None:
      ncfile = path + "/%s.nc" %(prefix)
   else: 
      ncfile = os.getenv('gotmout')+"/%s.nc" %(prefix)
      nc = NetCDFFile(ncfile)
   return nc

def getvar(nc,varib, *args, **kwargs):

    depth = nc.variables['z'][:,:,0,0]
    if varib == 'nuh':
       depth = nc.variables['zi'][:,:,0,0]
    nctime = nc['time']
    time = netCDF4.num2date(nctime[:], units=nctime.units)
    year = np.zeros(( time.shape[0] ))
    modeldates = netCDF4.num2date(nctime[:], units=nctime.units) 
    for t in range(time.shape[0]) :
       year[t] = time[t].year
    time = mdates.date2num(time)
    time = np.broadcast_to(time[:, np.newaxis], depth.shape)

    first = kwargs.get('first',None)
    first = year[0] if first is None else first

    last = kwargs.get('last',None)
    last = year[-1] if last is None else last

    yindex = np.where((year >= first) & (year <= last)) 

    if varib == 'time':
       var = np.copy(time)
    elif varib == 'date':
       var = np.copy(modeldates)
    elif varib == 'depth':
       var = np.copy(depth)
    elif varib == 'depthi':
       var = nc.variables['zi'][:,:,0,0]    
    elif varib == 'mld':
       var = nc.variables['mld_surf'][:,0,0]
    elif varib[0:7] == "ECO_sed":
       var = nc.variables[varib][:,0,0]
    elif varib == 'chl' :
       var1 = nc.variables['ECO_diachl'][:,:,0,0]
       var2 = nc.variables['ECO_flachl'][:,:,0,0]
       var = var1 + var2
       if 'ECO_coccochl' in nc.variables.keys():
         var3 = nc.variables['ECO_coccochl'][:,:,0,0]
         var = var + var3
       if 'ECO_bg' in nc.variables.keys():
         var4 = nc.variables['ECO_bgchl'][:,:,0,0]
         var = var + var4
    elif varib == 'chl_norwecom' :
       var = nc.variables['ECO_chla'][:,:,0,0]
    elif varib == 'pbiomass' :
       var1 = nc.variables['ECO_dia'][:,:,0,0]
       var2 = nc.variables['ECO_fla'][:,:,0,0]
       var = var1 + var2
       if 'ECO_cocco' in nc.variables.keys():
         var3 = nc.variables['ECO_cocco'][:,:,0,0]
         var = var + var3
       if 'ECO_bg' in nc.variables.keys():
         var4 = nc.variables['ECO_bg'][:,:,0,0]
         var = var + var4
    elif varib == 'zbiomass' :
       var1 = nc.variables['ECO_microzoo'][:,:,0,0]
       var2 = nc.variables['ECO_mesozoo'][:,:,0,0]
       var = var1 + var2 
    elif varib == "nit":
       var = nc.variables['ECO_no3'][:,:,0,0]
       var = var / 12.01 / 6.625
    elif varib == "nit_norwecom":
       var = nc.variables['ECO_nit'][:,:,0,0]
       var = var / 14.007
    elif varib == "pho":
       var = nc.variables['ECO_pho'][:,:,0,0]
       var = var / 12.01 / 106.
    elif varib == "sil":
       var = nc.variables['ECO_sil'][:,:,0,0]
       var = var / 12.01 / 6.625
    elif varib == "pp":
       var = nc.variables['ECO_primprod'][:,:,0,0]
       var = var * 60. * 60. * 24.
    elif varib == "npp":
       var = nc.variables['ECO_npp'][:,:,0,0]
       var = var * 60. * 60. * 24.   
    elif varib == "sp":
       var = nc.variables['ECO_secprod'][:,:,0,0]
       var = var * 60. * 60. * 24.
    elif varib == 'dia' :
       var = nc.variables['ECO_dia'][:,:,0,0]   
    elif varib == 'fla' :
       var = nc.variables['ECO_fla'][:,:,0,0]
    elif varib == 'cocco' :
       var = nc.variables['ECO_cocco'][:,:,0,0] 
    elif varib == 'microzoo' :
       var = nc.variables['ECO_microzoo'][:,:,0,0]
    elif varib == 'mesozoo' :
       var = nc.variables['ECO_mesozoo'][:,:,0,0]
    elif varib == 'oxy' :
       var = nc.variables['ECO_oxy'][:,:,0,0]
    elif varib == 'oxy_norwecom' :
       var = nc.variables['ECO_oxy'][:,:,0,0]
       var = var/32.0
    elif varib == 'dom' :
       var = nc.variables['ECO_dom'][:,:,0,0]
    elif varib == 'det' :
       var = nc.variables['ECO_det'][:,:,0,0]
    elif varib == 'dic' :
       var = nc.variables['CO2_dic'][:,:,0,0]
    elif varib == "speed":
       var = nc.variables['ECO_snkspd'][:,:,0,0]
    elif varib == "pocs":
       var1 = nc.variables['ECO_microzoo'][:,:,0,0]
       var2 = nc.variables['ECO_det'][:,:,0,0]
       var3 = nc.variables['ECO_dia'][:,:,0,0]
       var4 = nc.variables['ECO_fla'][:,:,0,0]
       var = var1 + var2 + var3 + var4
       if 'ECO_cocco' in nc.variables.keys():
         var5 = nc.variables['ECO_cocco'][:,:,0,0]
         var = var + var5
    elif varib == "pocs_norwecom":
       var1 = nc.variables['ECO_mic'][:,:,0,0]
       var2 = nc.variables['ECO_det'][:,:,0,0]
       var3 = nc.variables['ECO_dia'][:,:,0,0]
       var4 = nc.variables['ECO_fla'][:,:,0,0]
       var = var1 + var2 + var3 + var4
    elif varib == "poc":
       var1 = nc.variables['ECO_microzoo'][:,:,0,0]
       var2 = nc.variables['ECO_det'][:,:,0,0]
       var3 = nc.variables['ECO_dia'][:,:,0,0]
       var4 = nc.variables['ECO_fla'][:,:,0,0]
       var = var1 + var2 + var3 + var4
       if 'ECO_cocco' in nc.variables.keys():
         var5 = nc.variables['ECO_cocco'][:,:,0,0]
         var = var + var5
       var6 = nc.variables['ECO_mesozoo'][:,:,0,0]
       var = var + var6
    else :
       if nc.variables[varib].ndim == 3:
          var = nc.variables[varib][:,0,0]
       else:
          var = nc.variables[varib][:,:,0,0]


    if varib == 'date':
       var = var[yindex]
    elif len(var.shape) == 2:
       var = var[yindex,:]
       var = var[0,:,:]
    elif len(var.shape) == 1:
       var = var[yindex]

    return var

def get_title(variable):

   if variable == "chl":
      title = r"Chlorophyll $a$ (mgChl m$^{-3}$)"
   elif variable == "nit":
      title = r"Nitrate (mmolN m$^{-3}$)"
   elif variable == "pho":
      title = r"Phosphate (mmolP m$^{-3}$)"
   elif variable == "sil":
      title = r"Silicate (mmolSi m$^{-3}$)"
   elif variable == "dia":
      title = r"Diatoms (mgC m$^{-3}$)"
   elif variable == "fla":
      title = r"Flagellates (mgC m$^{-3}$)"
   elif variable == "cocco":
      title = r"Coccolithophores (mgC m$^{-3}$)"
   elif variable == "microzoo":
      title = r"Microzooplankton (mgC m$^{-3}$)"
   elif variable == "mesozoo":
      title = r"Mesozooplankton (mgC m$^{-3}$)"
   elif variable == "oxy":
      title = r"Oxygen (mmol m$^{-3}$)" 
   elif variable == "dom":
      title = r"Dissolved Organic Matter (mgC m$^{-3}$)"
   elif variable == "det":
      title = r"Detritus (mgC m$^{-3}$)"
   elif variable == "dic":
      title = r"Dissolved Inorganic Carbon (mmol m$^{-3}$)" 
   elif variable == "speed":
      title = r"Detritus sinking speed (m d$^{-1}$)"
   elif variable == "pbiomass":
      title = r"Phytoplankton biomass (mgC m$^{-3}$)"
   elif variable == "zbiomass":
      title = r"Zooplankton biomass (mgC m$^{-3}$)" 
   elif variable == "export":
      title = r"C-export (mgC m$^{-2}$ d$^{-1}$)"
   elif variable == 'npp':
      title = r'Net primary production (mgC m$^{-3}$ d$^{-1}$)'
   else:
      title = variable
   
   return title

def get_colormap(variable):
    import cmocean
    #import nclcmaps
    if variable == "chl" or variable == "dia" or variable == "fla" or variable == "microzoo" or variable == "mesozoo" or variable == "cocco"  or variable == "total_chlorophyll":
       #cmap = cmocean.cm.algae
       cmap = plt.cm.get_cmap('Spectral_r')
       #cmap = nclcmaps.cmap('WhiteBlueGreenYellowRed')
       #cmap = cmocean.cm.curl
       #cmap = cmocean.cm.speed_r
    elif variable == "nit" or variable == "pho" or variable == "sil":
       #cmap = cmocean.cm.balance
       cmap = plt.cm.get_cmap('Spectral_r')
    elif variable == "temp":
       cmap = cmocean.cm.thermal
    elif variable == "salt":
       cmap = cmocean.cm.haline
    elif variable == "speed":
       cmap = cmocean.cm.speed_r
    elif variable == "oxy":
       #cmap = cmocean.cm.curl 
       cmap = plt.cm.get_cmap('Spectral_r') 
    else: 
       cmap = plt.cm.get_cmap('Spectral_r')
   
    return cmap
        

def pmeshmodel(ncname,variable, *args, **kwargs):
    plt.style.use('seaborn-v0_8-notebook')

    nc  = read(ncname)
    depth  = getvar(nc,"depth")
    time  = getvar(nc,"time")
    var = getvar(nc,variable)

    first = kwargs.get('first',None)
    if first is not None:
       y = int(first[0:4]); m = int(first[5:7]); d = int(first[8:10])
       plotmin = datetime.datetime(y,m,d)
    else: 
       plotmin = time[0,0]

    last = kwargs.get('last',None)
    if last is not None:
       y = int(last[0:4]); m = int(last[5:7]); d = int(last[8:10])
       plotmax = datetime.datetime(y,m,d)       
    else: 
       plotmax = time[-1,0]

    dtop = kwargs.get('dtop',None)
    dtop = (depth).max() if dtop is None else -dtop
    dbot = kwargs.get('dbot',None)
    dbot = (depth).min() if dbot is None else -dbot 

    cmap = get_colormap(variable) 
    # the following takes the indexes of the variables in order to get the accurate colormap
    ind1 = np.abs(dtop - depth[0,:]).argmin() ; ind1 = min(ind1 + 1, depth.shape[1]-1)
    ind2 = np.abs(dbot - depth[0,:]).argmin() ; ind2 = max(ind2 - 1, 0)
    time = time[:,ind2:ind1+1]
    depth = depth[:,ind2:ind1+1]
    var = var[:,ind2:ind1+1]

    scale = kwargs.get('scale',None)
    if scale is not None:
       var = var * scale

    cmax = kwargs.get('cmax',None)
    cmax = var.max() if cmax is None else cmax
    cmin = kwargs.get('cmin',None)
    cmin = var.min() if cmin is None else cmin

    fig = plt.figure(figsize=(10,3.5),facecolor='w')
    ax  = fig.add_subplot(1,1,1)
    ax.set_position([0.125,0.175,0.75,0.725])
    pmesh = plt.pcolormesh(time,depth,var,cmap=cmap,shading='auto')
    
    ax.xaxis.axis_date()
    ax.yaxis.set_tick_params(labelsize='large')
    ax.xaxis.set_tick_params(labelsize='large')
    plt.xticks(rotation = 25)

    plt.clim(cmin,cmax)
    plt.clim(cmin,cmax)
    ax.set_xlim([plotmin,plotmax])
    ax.set_ylim([dbot,dtop])
    ax.grid(color='k', linestyle=':', linewidth=1)

    plt.title(get_title(variable), loc='left', fontsize=15)
    turn_off_secondary_title = kwargs.get('turn_off_secondary_title',None)
    if turn_off_secondary_title is not None:
       pass
    else:
       plt.title("filename: "+ncname,loc="right", fontsize=10)

    ax.set_ylabel("Depth (m)",fontsize=15)

    cbaxes = fig.add_axes([0.9, 0.25, 0.025, 0.5])
    cb = plt.colorbar(pmesh, cax = cbaxes)
    cbaxes.yaxis.set_tick_params(labelsize='large')

def vs_occci_chl(ncname, *args, **kwargs):

    plt.style.use('seaborn-v0_8-notebook')
    temperature = kwargs.get('temperature',None)

    nc  = read(ncname)
    time  = getvar(nc,"time")
    if temperature is not None:
       var = getvar(nc,"temp")
    else:
       var = getvar(nc,"chl")

    var1d = var[:,-1]
    # this function expects you to be in the experiment folder
    # and the cci_chl.dat present

    dates = list()
    satchl = list()

    satname = "./cci_chl.dat"
    
    if temperature is not None:
      satname = "./cci_sst.dat" 

    with open(satname , "r") as filestream:
     next(filestream)
     for line in filestream:
         currentline = line.split(" ")
         dates.append( datetime.datetime( int(currentline[0][0:4]),int(currentline[0][5:7]),int(currentline[0][8:10]) ) )
         values = currentline[1].split("\t")
         satchl.append(float(values[1]))

    if kwargs.get('first'):
       plotmin = datetime.datetime(first,1,1)
    else: 
       plotmin = time[0,0]

    if kwargs.get('last'):
       plotmax = datetime.datetime(last,1,1)
    else: 
       plotmax = time[-1,0]

    fig = plt.figure(figsize=(10,3.5),facecolor='w')
    ax  = fig.add_subplot(1,1,1)
    if temperature is not None:
       plt.plot(time[:,0],var1d)
    else:     
       plt.plot(time[:,0],np.log10(var1d))
    plt.plot(dates,satchl,".")
    ax.xaxis.axis_date()
    ax.yaxis.set_tick_params(labelsize='large')
    ax.xaxis.set_tick_params(labelsize='large')
    plt.xticks(rotation = 25)

    if temperature is not None:
       plt.title(r"Surface temperature", loc='left', fontsize=15)
    else:
       plt.title(r"Surface chlorophyll $a$", loc='left', fontsize=15)
    turn_off_secondary_title = kwargs.get('turn_off_secondary_title',None)
    if turn_off_secondary_title is not None:
       pass
    else:
       plt.title("filename: "+ncname,loc="right", fontsize=10)
    if temperature is not None:
       ax.set_ylabel(r"$^{o}$C",fontsize=15) 
    else:
       ax.set_ylabel(r"log$_{10}$(mgChl m$^{-3}$)",fontsize=15)

    ax.set_xlim([plotmin,plotmax])
    plt.tight_layout()
       
def vs_pft(ncname, *args, **kwargs):

   plt.style.use('seaborn-v0_8-notebook')

   nc  = read(ncname)
   time  = getvar(nc,"time")
   mchl = getvar(nc,"chl")
   mdia = getvar(nc,"ECO_diachl")
   mfla = getvar(nc,"ECO_flachl")
   mchl1d = mchl[:,-1]
   mdia1d = mdia[:,-1]
   mfla1d = mfla[:,-1]
   if 'ECO_coccochl' in nc.variables.keys():
      mcocco = getvar(nc,'ECO_coccochl')
      mcocco1d = mcocco[:,-1]
    # this function expects you to be in the experiment folder
    # and the cci_chl.dat present
   import pickle
   f = open("PFTs.dat",'rb')
   chl,dia,dino,hapto,green,prochlo,prokar,micro,nano,pico,dates = pickle.load(f)

   if kwargs.get('first'):
      plotmin = datetime.datetime(first,1,1)
   else: 
      plotmin = time[0,0]

   if kwargs.get('last'):
      plotmax = datetime.datetime(last,1,1)
   else: 
      plotmax = time[-1,0]

   fig = plt.figure(figsize=(10,3.5),facecolor='w')
   ax  = fig.add_subplot(1,1,1)
   ax.set_position([0.075,0.15,0.73,0.775])
   plt.plot(time[:,0],np.log10(mchl1d),label='Model total chla',color='#0f6fb9')
   plt.plot(time[:,0],np.log10(mdia1d),label='Model diatom chla',color='#f05e31')
   plt.plot(time[:,0],np.log10(mfla1d),label='Model flagellate chla',color='#0aa775')

   if 'ECO_coccochl' in nc.variables.keys():
      plt.plot(time[:,0],np.log10(mcocco1d),label='Model coccoliths chla',color='#da70aa') 
      plt.plot(dates,np.log10(hapto),".",label="satellite haptophytes chl",color='#da70aa')

   plt.plot(dates,np.log10(chl),".",label="satellite chl",color='#0f6fb9')
   plt.plot(dates,np.log10(micro),".",label="satellite diatom chl",color='#f05e31')
   plt.plot(dates,np.log10(dino+pico),".",label="satellite dino+pico chl",color='#0aa775')

   ax.xaxis.axis_date()
   ax.yaxis.set_tick_params(labelsize='large')
   ax.xaxis.set_tick_params(labelsize='large')
   plt.xticks(rotation = 25)

   plt.title(r"Surface chlorophyll $a$", loc='left', fontsize=15)
   turn_off_secondary_title = kwargs.get('turn_off_secondary_title',None)
   if turn_off_secondary_title is not None:
      pass
   else:
      plt.title("filename: "+ncname,loc="right", fontsize=10)
   ax.set_ylabel(r"log$_{10}$(mgChl m$^{-3}$)",fontsize=15)

   ax.set_xlim([plotmin,plotmax])
   ax.set_ylim([-2.,1.5])
   plt.legend(loc='upper right', bbox_to_anchor=(0.777, 0.5, 0.5, 0.5))
   #plt.legend()
   #plt.tight_layout()    

def compare_pmesh(filename1,filename2, variable, *args, **kwargs):
    first = kwargs.get('first',None) 
    last  = kwargs.get('last',None)
    dbot  = kwargs.get('dbot',None)
    cmax  = kwargs.get('cmax',None)
    cmin  = kwargs.get('cmin',None)
    dtop  = kwargs.get('dtop',None)

    pmeshmodel(filename1, variable, first=first, last=last, dtop=dtop, dbot=dbot, cmax=cmax, cmin=cmin)
    fig1 = plt.gcf()
    pmeshmodel(filename2, variable, first=first, last=last, dtop=dtop, dbot=dbot, cmax=cmax, cmin=cmin)
    fig2 = plt.gcf()

    fig = plt.figure(figsize=(10,6),facecolor='w')
    gs = gridspec.GridSpec(2,2)

    ax1 = fig1.axes[0]
    ax1.remove()
    ax1.figure = fig
    fig.add_axes(ax1)
    ax1.set_subplotspec(gs[0, 0])
    ax1.set_position([0.125,0.175,0.75,0.55])

    ax2 = fig2.axes[0]
    ax2.remove()
    ax2.figure = fig
    fig.add_axes(ax2)
    ax2.set_subplotspec(gs[0, 1])
    ax2.set_position([0.125,1.,0.75,0.55])

    cbaxes = fig.add_axes([0.9, 0.25, 0.025, 0.5])
    cb = plt.colorbar(cax = cbaxes)
    cbaxes.yaxis.set_tick_params(labelsize='large') 

    plt.close(fig1)
    plt.close(fig2)

def ensemble_pmesh(argo,experiments,variable, *args, **kwargs):
   first = kwargs.get('first',None) 
   last  = kwargs.get('last',None)
   dbot  = kwargs.get('dbot',None)
   cmax  = kwargs.get('cmax',None)
   cmin  = kwargs.get('cmin',None)
   dtop  = kwargs.get('dtop',None)

   topfig = plt.figure(figsize=(8,max([len(experiments)*2,12])),facecolor='w')
   gs = gridspec.GridSpec(1,len(experiments))

   count = 1
   for item in experiments:
      os.environ["gotmout"] = "/cluster/work/users/cagyum/ENSEMBLE_GOTM/ENSEMBLE_Argo"+str(argo)+"/"+str(item)  
      modelinput.pmeshmodel("result", variable, first=first, last=last, dtop=dtop, dbot=dbot, cmax=cmax, cmin=cmin)
      fig = plt.gcf()
      ax = fig.axes[0]
      ax.remove()
      ax.figure = topfig.add_subplot(len(experiment),1,count)
      #topfig.add_axes(ax)
      #topfig.add_subplot(len(experiment),1,count)
      #ax.set_subplotspec(gs[0, count])
      #ax  = fig.add_subplot(1,1,1)
      if count == 1:
         pass
#         ax.set_position([0.125,0.175,0.75,1./len(experiments)])
      else:
         #ax.set_position([0.125,count*0.175+1./len(experiments),0.75,1./len(experiments)])
         plt.setp( ax.get_xticklabels(), visible=False)
      
      plt.close(fig)
      count = count + 1

def distance_on_unit_sphere(lat1, long1, lat2, long2):

# Convert latitude and longitude to
# spherical coordinates in radians. Assumes the Earth is a perfect sphere.
  degrees_to_radians = math.pi/180.0
# phi = 90 - latitude
  phi1 = (90.0 - lat1)*degrees_to_radians
  phi2 = (90.0 - lat2)*degrees_to_radians
# theta = longitude
  theta1 = long1*degrees_to_radians
  theta2 = long2*degrees_to_radians

  cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
  math.cos(phi1)*math.cos(phi2))
  arc = math.acos( cos )
  arc = arc * 6378.137 # convert to km

  return arc

def get_radius_index(nc,lon,lat,slon,slat, *args, **kwargs):
   # used for getting the satellite indexes around the model/argo coordinates
   II = np.abs(lon-slon).argmin()
   JJ = np.abs(lat-slat).argmin()
   pick = kwargs.get('pick',None)
   pick = 10 if pick is None else pick
   radius = kwargs.get('radius',None)
   radius = 25 if radius is None else radius

   dummy = nc.variables['CHL'][0,JJ-pick:JJ+pick,II-pick:II+pick]
   dist = np.zeros((dummy.shape))
   for i in range(dummy.shape[0]):
      for j in range(dummy.shape[1]):
         dist[j,i] = distance_on_unit_sphere(lat,lon,slat[JJ-pick:JJ+pick][j],slon[II-pick:II+pick][i]) 
   index = dist<=radius
   return index,II,JJ,pick

def get_pft(*args, **kwargs):

   #nc = NetCDFFile("https://{}:{}@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D".format(str(os.environ['Copernicus_user']),str(os.environ['Copernicus_pass'])))
   #nc = NetCDFFile("https://{}:{}@my.cmems-du.eu/thredds/dodsC/c3s_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1M".format(str(os.environ['Copernicus_user']),str(os.environ['Copernicus_pass'])))
   nc = NetCDFFile("https://{}:{}@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D".format(str(os.environ['Copernicus_user']),str(os.environ['Copernicus_pass'])))
   nctime = nc['time']
   time = netCDF4.num2date(nctime[:], units=nctime.units)
   time = mdates.date2num(time)
   sdates = netCDF4.num2date(nctime[:], units=nctime.units)

   argo = kwargs.get('argo',None)
   if argo is not None:
      import argoinput
      Apres, Asince, Atime2d, Adates, Alon, Alonint, Alat, Alatint, Atimeint, Atemp, Atempqc = argoinput.readARGO(argo, 'TEMP_ADJUSTED',all=True) 
      Atime2d  = argoinput.construct(Atime2d,Atempqc,[1,5,8])
      Atime1d=[]
      for i in range(Atime2d.shape[0]):
          Atime1d.append(Atime2d[i,:].max())
      argonc = NetCDFFile("GL_PR_PF_"+str(argo)+".nc")
      if "TIME" in argonc.variables.keys():
         argotime = argonc["TIME"]
      else:
         argotime = argonc["JULD"] 
      argotime = netCDF4.num2date(argotime[:], units=argotime.units)
#      argotime = argotime[~Atime2d[:,0].mask]
      argotime = argotime[~np.isnan(np.array(Atime1d))]
#      Alon = Alon[~Atime2d[:,0].mask]
#      Alat = Alat[~Atime2d[:,0].mask]
      Alon = Alon[~np.isnan(np.array(Atime1d))]
      Alat = Alat[~np.isnan(np.array(Atime1d))]
      argoseconds = np.zeros((len(argotime)))
      for s in range(len(argotime)):
         argoseconds[s] = (argotime[s]-datetime.datetime(Asince,1,1)).total_seconds()
      satt1 = np.abs(argotime[0]-sdates).argmin()
      satt2 = np.abs(argotime[-1]-sdates).argmin()
      lonint = list()
      latint = list()
      d = datetime.datetime( sdates[satt1].year,sdates[satt1].month,sdates[satt1].day )
      d2 = datetime.datetime( sdates[satt2].year,sdates[satt2].month,sdates[satt2].day ) 
      while d <= d2 :
         satseconds = (d-datetime.datetime(Asince,1,1)).total_seconds()
         lonint.append( np.interp( satseconds, argoseconds, Alon  ) )
         latint.append( np.interp( satseconds, argoseconds, Alat  ) )
         d = d + relativedelta(days=1)
   else:
      lon = kwargs.get('lon')
      lat = kwargs.get('lat')

   first = kwargs.get('first',None)
   if first is not None:
      y = int(first[0:4]); m = int(first[5:7]); d = int(first[8:10])
      tfirst = np.abs(datetime.datetime(y,m,d)-sdates).argmin()
   else:
      if argo is not None:
         tfirst = satt1
      else:
         tfirst = 0
   last = kwargs.get('last',None) 
   if last is not None:
      y = int(last[0:4]); m = int(last[5:7]); d = int(last[8:10])
      tlast = np.abs(datetime.datetime(y,m,d)-sdates).argmin()
   else: 
      if argo is not None:
         tlast = satt2
      else:
         tlast = -1   

   slat = nc.variables['lat'][:]
   slon = nc.variables['lon'][:]

   if argo is None:
      index,II,JJ,pick = get_radius_index(nc,lon,lat,slon,slat) 

   chl = list()
   dia = list()
   dino = list()
   hapto = list()
   green = list()
   prochlo = list()
   prokar = list()
   micro = list()
   nano = list()
   pico = list()
   dates = list()
   for t in range(tfirst,tlast+1):
      print(sdates[t])
      dates.append(time[t])
      if argo is not None:
         index,II,JJ,pick = get_radius_index(nc,lonint[t-tfirst],latint[t-tfirst],slon,slat)

      try:
         chl.append( np.nanmean( nc.variables['CHL'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         chl.append(-999.)
      try:
         dia.append( np.nanmean( nc.variables['DIATO'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         dia.append(-999.)
      try:
         dino.append( np.nanmean( nc.variables['DINO'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         dino.append(-999.)
      try:
         hapto.append( np.nanmean( nc.variables['HAPTO'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         hapto.append(-999.)
      try:
         prokar.append( np.nanmean( nc.variables['PROKAR'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         prokar.append(-999.)
      try:
         prochlo.append( np.nanmean( nc.variables['PROCHLO'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         prochlo.append(-999.)
      try:
         micro.append( np.nanmean( nc.variables['MICRO'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         micro.append(-999.)
      try:
         nano.append( np.nanmean( nc.variables['NANO'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         nano.append(-999.)
      try:
         pico.append( np.nanmean( nc.variables['PICO'][t,JJ-pick:JJ+pick,II-pick:II+pick][index] ) )
      except Exception:
         pico.append(-999.)

   chl = np.ma.masked_where(np.array(chl)<0.,np.array(chl))
   dia = np.ma.masked_where(np.array(dia)<0.,np.array(dia))
   dino = np.ma.masked_where(np.array(dino)<0.,np.array(dino))
   hapto = np.ma.masked_where(np.array(hapto)<0.,np.array(hapto))
   green = np.ma.masked_where(np.array(green)<0.,np.array(green))
   prochlo = np.ma.masked_where(np.array(prochlo)<0.,np.array(prochlo))
   prokar = np.ma.masked_where(np.array(prokar)<0.,np.array(prokar))
   micro = np.ma.masked_where(np.array(micro)<0.,np.array(micro))
   nano = np.ma.masked_where(np.array(nano)<0.,np.array(nano))
   pico = np.ma.masked_where(np.array(pico)<0.,np.array(pico))

   if kwargs.get('write_it'):
      import pickle
      f = open('PFTs.dat','wb')
      pickle.dump([chl,dia,dino,hapto,green,prochlo,prokar,micro,nano,pico,dates],f)
      f.close()

   return chl,dia,dino,hapto,green,prochlo,prokar,micro,nano,pico,dates

def meshdata(x,y,z,start,end,dmax,modeldepth, *args, **kwargs):
   # I've used this to calculate statistisc
   # it basically averages data within a range in time and depth
   # I found it a bit buggy when your depth interval criteria does not capture any model points
   y1 = int(start[0:4]); m1 = int(start[5:7]); d1 = int(start[8:10])
   y2 = int(end[0:4]); m2 = int(end[5:7]); d2 = int(end[8:10])
   dates = list()

   month_from = kwargs.get('month_from',None) # these will exlude months before and after
   month_to   = kwargs.get('month_to',None) 
   month_index = list()

   xave = kwargs.get('xave',None)
   xave = 10 if xave is None else xave 
#   yave   = kwargs.get('yave',None)
#   yave = 5. if yave is None else yave

   interface = np.zeros((modeldepth.shape[0]+1))
   modeldepth = np.flipud(-modeldepth)
   for dd in range(1,len(interface)):
       interface[dd] = ( modeldepth[dd-1] - interface[dd-1] ) + modeldepth[dd-1]
   depth = modeldepth

   first_date = datetime.datetime(y1,m1,d1)
   while first_date <= datetime.datetime(y2,m2,d2):
      dates.append( (first_date - datetime.datetime(1950,1,1)).days * 86400.)
      if month_from is not None:
         if np.logical_and(first_date.month>=month_from , first_date.month<=month_to):
            month_index.append(False)
         else:
            month_index.append(True) 
      first_date = first_date + datetime.timedelta(days=xave)


#   depth = np.arange(0,dmax+yave,yave)  
   try:
      # True for Argo
      dummy = x.shape[1] # only to pass the try test
      x = x * 86400. # convert to seconds
#      modeldepth = kwargs.get('modeldepth',None)
#      interface = np.zeros((modeldepth.shape[0]+1))
#      modeldepth = np.flipud(-modeldepth)
#      for dd in range(1,len(interface)):
#          interface[dd] = ( modeldepth[dd-1] - interface[dd-1] ) + modeldepth[dd-1]
#      depth = modeldepth
   except Exception:
      # True for model
      delta = x - datetime.datetime(1950,1,1)
      x = np.zeros((len(delta)))
      for t in range(len(delta)):
         x[t] = (delta[t].days * 86400. + delta[t].seconds)
      x = np.broadcast_to(x[:, np.newaxis], z.shape)
      y = np.flipud(-y)
      z = np.flipud(z)
      x = np.flipud(x)

#      interface = np.zeros((y.shape[1]+1))
#      for dd in range(1,len(interface)):
#          interface[dd] = ( y[0,dd-1] - interface[dd-1] ) + y[0,dd-1]
#      depth = y[0,:]
   x = x.flatten(); y = y.flatten(); z = z.flatten()
   x = x[~y.mask]; z = z[~y.mask]; y = y[~y.mask]
   x = x[~z.mask]; y = y[~z.mask]; z = z[~z.mask]
   x = np.squeeze(x); y = np.squeeze(y); z = np.squeeze(z);

#   print(depth)
   dindex = np.abs(dmax-depth).argmin()
   interface = interface[0:dindex+1]
   intp = np.zeros((len(interface)-1,len(dates)))-99999.
#   print(interface)
   for b in range(len(dates)-1):
#      for a in range(len(depth)-1):
      for a in range(interface.shape[0]-1):
#         index = np.logical_and(np.logical_and(y>=depth[a],y<depth[a+1]),np.logical_and(x>=dates[b],x<dates[b+1]))
         index = np.logical_and(np.logical_and(y>=interface[a],y<interface[a+1]),np.logical_and(x>=dates[b],x<dates[b+1]))
         try:
            intp[a,b] = np.nanmean(z[index])
         except Exception:
            intp[a,b] = -99999.
#      print("date "+str(b+1)+" of "+str(len(dates)-1)) 
   intp = np.ma.masked_where(intp < -99990., intp)
   if month_from is not None:
      index2d = np.broadcast_to(np.array(month_index).T[np.newaxis,:], intp.shape)
      intp = np.ma.masked_where(index2d,intp)

   xm,ym = np.meshgrid(dates,depth[0:dindex])
   mask_below = kwargs.get('mask_below',None)
   if mask_below is not None:
      intp = np.ma.masked_where(intp<mask_below,intp)
   return xm,ym,intp

def plotlims(variable,values):
   if variable == "ECO_Nlim" or variable == "ECO_Plim" or \
      variable == "ECO_Slim" or variable == "ECO_Llim" :
         ymin = 0; ymax = 1
   elif variable == "chl" or variable == "Argo CHL":
      ymin = 0; ymax = 5
   elif variable == "oxy" or variable == "Argo OXY":
      ymin = 270; ymax = 340
   else :
      ymin = min(values); ymax = max(values)
   return ymin, ymax

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    dummy = (cumsum[N:] - cumsum[:-N]) / float(N)
    dummy2 = list()
    dummy2.append(x[0])
    for t in range(len(dummy)):
       dummy2.append(dummy[t])
    dummy2.append(x[-1])
    return dummy2

def plot1d(ncname,variable, d, *args, **kwargs):
   plt.style.use('seaborn-v0_8-notebook')
   from mpl_toolkits.axes_grid1 import host_subplot
   import mpl_toolkits.axisartist as AA

   nc  = read(ncname)
   depth  = getvar(nc,"depth")

   nctime = nc['time']
   time = netCDF4.num2date(nctime[:], units=nctime.units)
   dates = netCDF4.num2date(nctime[:], units=nctime.units) 
   time = mdates.date2num(time)

   first = kwargs.get('first',None)
   if kwargs.get('first'):
      y = int(first[0:4]); m = int(first[5:7]); dd = int(first[8:10])
      plotmin = datetime.datetime(y,m,dd)
   else: 
      plotmin = dates[0]

   last = kwargs.get('last',None)
   if kwargs.get('last'):
      y = int(last[0:4]); m = int(last[5:7]); dd = int(last[8:10])
      plotmax = datetime.datetime(y,m,dd)
   else: 
      plotmax = dates[-1] 

   axnames = ["ax1","ax2","ax3","ax4","ax5","ax6","ax7","ax8","ax9","ax10"]
   colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
   count=0

   fig = plt.figure(figsize=(10,3.5),facecolor='w')
   globals()[ axnames[count] ]  = host_subplot(111, axes_class=AA.Axes)
   #globals()[ axnames[count] ]  = fig.add_subplot(1,1,1)
   first_ax = globals()[ axnames[count] ]
#   first_ax.set_position([0.075,0.175,0.73,0.7725])

   results = {}
   max1d = -99999.0 # will be modified to full set max value
   h = getvar(nc,'h')
   integrate = kwargs.get('integrate', False)
   average = kwargs.get('average', False)
   fixedy = kwargs.get('fixedy', False)
   legend = kwargs.get('legend',False)
   addx = kwargs.get('addx',None)
   addy = kwargs.get('addy',None)
   addn = kwargs.get('addn',None) 
   if isinstance(variable, str):
      dummy = variable
      variable = list()
      variable.append(dummy)
   for selvar in variable:
      if count > 0:
         globals()[ axnames[count] ] = first_ax.twinx()
         new_fixed_axis = globals()[ axnames[count] ].get_grid_helper().new_fixed_axis
         #globals()[ axnames[count] ].set_position([0.075,0.175,0.73,0.77])
         #new_fixed_axis
         globals()[ axnames[count] ].axis["right"] = new_fixed_axis(loc="right",
                                        axes=globals()[ axnames[count] ],
                                        offset=(60*(count-1), 0))
      var = getvar(nc,selvar)
      var1d = []
      if var.ndim == 1:
         for t in range(len(time)):
            var1d.append(var[t])
      else:   
         for t in range(len(time)):
            if integrate:
               index = np.abs(-depth[t, :] - d).argmin()
               if -depth[t,index] < d:
                  index = index -1
               lastmass = (-depth[t,index] - d) * var[t,index]
               restmass = np.sum(h[t,(index+1):]* var[t,(index+1):])
               var1d.append( lastmass + restmass )
            elif average:
               index = np.abs(-depth - d).argmin()
               if -depth[t,index] < d:
                  index = index -1
               lastmass = (-depth[t,index] - d) * var[t,index]
               restmass = np.sum(h[t,(index+1):]* var[t,(index+1):])
               var1d.append( (lastmass + restmass)/float(d) )
            else:
               var1d.append( np.interp( -d, depth[t,:], var[t,:] ))
      #ymin,ymax = plotlims(selvar,var1d)
      if kwargs.get('runmean'):
         var1d = running_mean(var1d, 3)

      globals()[ axnames[count] ].plot(time,var1d,label=selvar,color=colors[count])
      globals()[ axnames[count] ].set_ylabel(selvar)
      globals()[ axnames[count] ].axis["right"].label.set_color(colors[count])
      globals()[ axnames[count] ].set_ylim(0.0, max(var1d)*1.1)
      plt.xticks(rotation = 25)
      count = count + 1

      results[selvar] = var1d

      max1d = max([max(var1d), max1d]) 
   
   if addx is not None:
      addy = kwargs.get('addy',None)
      addn = kwargs.get('addn',None)
      
      for addnumber in range(len(addx)):
         max1d = max([max(np.array(addy[addnumber])), max1d]) # for now only supports 1 addx 
         globals()[ axnames[count] ] = first_ax.twinx()
         new_fixed_axis = globals()[ axnames[count] ].get_grid_helper().new_fixed_axis
         globals()[ axnames[count] ].axis["right"] = new_fixed_axis(loc="right",
                                        axes=globals()[ axnames[count] ],
                                        offset=(60*(count-1), 0))
         globals()[ axnames[count] ].plot(addx[addnumber],addy[addnumber],label=addn[addnumber],color=colors[count])
         globals()[ axnames[count] ].set_ylabel(addn[addnumber])
         globals()[ axnames[count] ].axis["right"].label.set_color(colors[count])
         ymin,ymax = plotlims(addn[addnumber],addy[addnumber])
         globals()[ axnames[count] ].set_ylim(ymin, ymax*1.1) 
         plt.xticks(rotation = 25)
         count = count + 1                  

   if fixedy:
       countmax = 0
       for selvar in variable:
             globals()[ axnames[countmax] ].set_ylim(0.0, max1d*1.1) # set the range for the first y-axis
             countmax = countmax + 1
       if addx is not None:
             for addnumber in range(len(addx)):
                 globals()[ axnames[countmax+addnumber] ].set_ylim(0.0, max1d*1.1)   

   first_ax.xaxis.axis_date()
   first_ax.yaxis.set_tick_params(labelsize='large')
   first_ax.xaxis.set_tick_params(labelsize='large')
   first_ax.axis["left"].label.set_color(colors[0])
   first_ax.grid(color='k', linestyle=':', linewidth=1)


   turn_off_secondary_title = kwargs.get('turn_off_secondary_title',None)
   if turn_off_secondary_title is not None:
      pass
   else:
      plt.title("filename: "+ncname,loc="right", fontsize=10)

   if integrate:
       plt.title("0 - "+str(d)+" meters integrate",loc="left",fontsize=13)
   elif average:
       plt.title("0 - "+str(d)+" meters average",loc="left",fontsize=13)
   else: 
       plt.title("@ "+str(d)+" meters",loc="left",fontsize=13)

   first_ax.set_xlim([plotmin,plotmax])

   if legend: plt.legend()

   #plt.tight_layout()

   
   

   return time, results

def normalize(data):
   out = np.zeros((data.shape))
   for t in range(data.shape[1]):
      out[:,t] = data[:,t]/data[:,t].max()
#   out = data/np.percentile(data,99.5)
   out = np.ma.masked_where(out < -9999., out)
   return out

def compare_multi_model_variables(file1,file2,variables, *args, **kwargs):
    import cmocean
    scale_diff_factor = kwargs.get('scale_diff_factor',None)
    if scale_diff_factor is not None:
       sc = scale_diff_factor
    else:
       sc = np.ones((len(variables)))
    plt.style.use('seaborn-v0_8-notebook')

    nc1  = read(file1)
    nc2  = read(file2)
    depthm  = getvar(nc1,"depth")
    time  = getvar(nc1,"time")

    nctime = nc1['time']
    dates = netCDF4.num2date(nctime[:], units=nctime.units)

    first = kwargs.get('first',None)
    if kwargs.get('first'):
       y = int(first[0:4]); m = int(first[5:7]); d = int(first[8:10])
       plotmin = datetime.datetime(y,m,d)
    else: 
       plotmin = dates[0]

    last = kwargs.get('last',None)
    if kwargs.get('last'):
       y = int(last[0:4]); m = int(last[5:7]); d = int(last[8:10])
       plotmax = datetime.datetime(y,m,d)
    else: 
       plotmax = dates[-1]


    nvar = len(variables)
    data = {}
    for item in range(nvar):
        var1 = getvar(nc1,variables[item])
        var2 = getvar(nc2,variables[item])
        dtop = kwargs.get('dtop',None)
        dtop_s = max(depthm) if dtop is None else -dtop[item]
        dbot = kwargs.get('dbot',None)
        dbot_s = min(depthm) if dbot is None else -dbot[item]
        ind1 = np.abs(dtop_s - depthm[0,:]).argmin() ; ind1 = min(ind1 + 1, depthm.shape[1]-1)
        ind2 = np.abs(dbot_s - depthm[0,:]).argmin() ; ind2 = max(ind2 - 1, 0)
        timem = time[:,ind2:ind1+1]
        depth_s = depthm[:,ind2:ind1+1]
        fld1 = var1[:,ind2:ind1+1]
        fld2 = var2[:,ind2:ind1+1]
        cmax = kwargs.get('cmax',None)
        cmax_s = fld.max() if cmax is None else cmax[item]
        cmin = kwargs.get('cmin',None)
        cmin_s = fld.min() if cmin is None else cmin[item]

        data[item] = {'var1':fld1, 'var2':fld2, 'cmax': cmax_s, 'cmin': cmin_s, 'depth': depth_s, 'dtop': dtop_s, 'dbot': dbot_s, 'time': timem}

    fig = plt.figure(figsize=(5*nvar,7),facecolor='w')
    for i in range(nvar):
      for v in (1,2):
        ax  = fig.add_subplot(3,nvar,i+1+nvar*(v-1))
        cmap = get_colormap(variables[i])
        pmesh = plt.pcolormesh(data[i]['time'],data[i]['depth'],data[i]['var'+str(v)],cmap=cmap,shading='auto')
        ax.xaxis.axis_date()
        plt.xticks(rotation = 25)
        plt.clim(data[i]['cmin'],data[i]['cmax'])
        ax.set_xlim([plotmin,plotmax])
        ax.set_ylim([data[i]['dbot'],data[i]['dtop']])
        ax.grid(color='k', linestyle=':', linewidth=1)
        plt.title(get_title(variables[i]), loc='left', fontsize=15)
        if v == 1:
        #   plt.title("filename: "+file1,loc="right", fontsize=9)
           ax.set_ylabel(file1,fontsize=9)   
        if v == 2:
        #   plt.title("filename: "+file2,loc="right", fontsize=9)
           ax.set_ylabel(file2,fontsize=9)
        #ax.set_ylabel("Depth (m)",fontsize=15)
        plt.colorbar()

      ax  = fig.add_subplot(3,nvar,i+1+nvar*2)


      cmap = cmocean.cm.balance
      pmesh = plt.pcolormesh(data[i]['time'],data[i]['depth'],data[i]['var1']-data[i]['var2'],cmap=cmap,shading='auto')
      ax.xaxis.axis_date()
      plt.xticks(rotation = 25)
      plt.clim(-data[i]['cmax']*sc[i],data[i]['cmax']*sc[i])
      ax.set_xlim([plotmin,plotmax])
      ax.set_ylim([data[i]['dbot'],data[i]['dtop']])
      ax.grid(color='k', linestyle=':', linewidth=1)
      plt.title(get_title(variables[i]+' difference'), loc='left', fontsize=15)
      plt.title("top - middle",loc="right", fontsize=9)
      #ax.set_ylabel("Depth (m)",fontsize=15)
      plt.colorbar()

    plt.tight_layout()

def interp_to_1d(depth,var,d):
    var1d = []     
    for t in range(depth.shape[0]):
        var1d.append( np.interp( -d, depth[t,:], var[t,:] ))
    return(np.array(var1d))

def compare_multi_model_depths(files,variable,depths, *args, **kwargs):

    plt.style.use('seaborn-v0_8-notebook')
    nfile = len(files)
    if nfile > 1:
       nc1  = read(files[0])
    else:
       nc1  = read(files)

    depthm  = getvar(nc1,"depth")
    time  = getvar(nc1,"time")

    nctime = nc1['time']
    dates = netCDF4.num2date(nctime[:], units=nctime.units)

    first = kwargs.get('first',None)
    if kwargs.get('first'):
       y = int(first[0:4]); m = int(first[5:7]); d = int(first[8:10])
       plotmin = datetime.datetime(y,m,d)
    else: 
       plotmin = dates[0]

    last = kwargs.get('last',None)
    if kwargs.get('last'):
       y = int(last[0:4]); m = int(last[5:7]); d = int(last[8:10])
       plotmax = datetime.datetime(y,m,d)
    else: 
       plotmax = dates[-1]

    cmax = kwargs.get('cmax',None)
    cmin = kwargs.get('cmin',None)

    argo = kwargs.get('argo',None)
    if argo is not None:    
       import argoinput

    fig = plt.figure(figsize=(5*len(depths),3),facecolor='w')
    for d in range(len(depths)):
      ax  = fig.add_subplot(1,len(depths),d+1)
      for i in range(nfile):
          nc = read(files[i])
          if variable=='export':
             var1 = getvar(nc,'ECO_det')
             if 'ECO_snkspd' in nc.variables:
                var2 = getvar(nc,'speed')
                var = var1*var2
             else:
                speed = kwargs.get('speed',None)
                if speed is not None:
                   var = var1*speed[i]
                else:
                   var = var1*5.0
          else:
             var = getvar(nc,variable)
          v = interp_to_1d(depthm,var,depths[d])
          plt.plot(time[:,0],v,label=files[i])

      if argo is not None:
         d1,a1 = argoinput.get1D(argo,variable,float(depths[d]))
         plt.plot(d1,a1,label="BGC-Argo",color="k")

      cmax_s = v.max() if cmax is None else cmax
      cmin_s = v.min() if cmin is None else cmin
      plt.ylim(cmin_s,cmax_s)
      ax.xaxis.axis_date()
      plt.xticks(rotation = 25)
      ax.set_xlim([plotmin,plotmax])
      ax.grid(color='k', linestyle=':', linewidth=1)

      plt.title(get_title(variable), loc='left', fontsize=15)
      plt.title("@"+str(depths[d])+'m',loc="right", fontsize=9)        

      if d==len(depths)-1: plt.legend()
    plt.tight_layout()
