#!/usr/bin/env python3
#
############################################################################
#
# MODULE:    r.survey
# AUTHOR(S): Ivan Marchesini and Txomin Bornaxea
# PURPOSE:   Define solid angle, 3d distance and view angles area 
# from multiple survey locations (points). Survey points can be at any elevation
# above ground level. 3d points, representing drone/aerial positions, are allowed
# Outputs are:
# maximum solid angle each pixel is visible from survey points
# maximum view angle each pixel is visible from survey points
# minimum 3d distance each pixel is visible from survey points
# number of survey points each pixel is visible from
# identifier of the survey point having the minimum 3d distance to each given pixel
# identifier of the survey point having the maximum solid angle to each given pixel
# identifier of the survey point having the maximum view angle angle to each given pixel
# The modulel runs in parallel but if too many processes are used it can fail due to problem
# in the management of the temporary_regions. In this case please reduce the
# number of processes. 
#
#
#############################################################################

#%Module
#% description: Returns maps of visibility indexes (3d distance, view angle and solid angle) from multiple survey points
#% keywords: visibility
#% keywords: survey
#%end
#%option G_OPT_V_INPUT
#% key: points
#% description: Name of the input points map  (representing the survey location)
#% required: yes
#%end
#%option G_OPT_R_INPUT
#% key: dem
#% description: Name of the input DEM layer
#% required: yes
#%end
#%option G_OPT_R_OUTPUT
#% key: output
#% description: Prefix for the output visibility layers
#% required: yes
#%end
#%option
#% key: maxdist
#% type: double
#% description: max distance from the input points
#% required: yes
#% answer: 1000
#%end
#%option G_OPT_V_INPUT
#% key: treesmap
#% description: Name of the vector layer representing the forested areas
#% required: no
#%end
#%option G_OPT_V_INPUT
#% key: buildingsmap
#% description: Name of the vector layer representing the buildings
#% required: no
#%end
#%option
#% key: obs_heigh
#% description: observer_elevation
#% type: double
#% required: yes
#% answer: 1.75
#%end
#%option
#% key: treesheigh
#% description: field of the attribute table containing information about threes heigh
#% type: string
#% required: no
#%end
#%option
#% key: buildingsheigh
#% description: field of the attribute table containing information about buildings heigh
#% type: string
#% required: no
#%end
#%option
#% key: obsabselev
#% description: field of the attribute table containing information about observer absolute elevation above the same datum used by the DEM (e.g. geodetic elevation)
#% type: string
#% required: no
#%end
#%option
#% key: layer
#% description: layer of the attribute table of the point map (can be usefull for obsabselev option)
#% type: string
#% answer: 1
#% required: no
#%end
#%option
#% key: viewangle_threshold
#% description: cut the output layers at a given threshold value (between 90 and 180 degrees)
#% gisprompt: 90-180
#% type: double
#% required: yes
#%answer: 90.01
#%end
#%option
#% key: object_radius
#% description: radius of the surveyed object in unit of map (default is half the DEM resolution)
#% type: double
#% required: no
#%end
#%option
#% key: nprocs
#% description: Number of processes
#% answer: 1
#% type: integer
#required: yes
#%end
#%option
#% key: memory
#% description: Amount of memory to use in MB (for r.viewshed analysis)
#% answer: 500
#% type: integer
#required: no
#%end
#%flag
#% key: b
#% description: Create a simple visible-not visible raster map
#%end
#%flag
#% key: c
#% description: Consider the curvature of the earth (current ellipsoid)
#%end
#%flag
#% key: d
#% description: allow only downward view direction (can be used for drone view)
#%end


#Importing modules
import atexit
from subprocess import PIPE
from grass.script import parser, parse_key_val
from grass.pygrass.modules import Module
import multiprocessing
import threading
import sys
import grass.script as gscript
from grass.script import core as grasscore
from math import pi


#function for cleaning temporary layers
def cleanup():
    message = " *** Cleaning all temporary maps *** "
    gscript.message(message)
    Module("g.remove", type='vector', pattern="xxtemp*", quiet=True, flags="f")
    Module("g.remove", type='vector', pattern="zzpnt*", quiet=True, flags="f")
    Module("g.remove", type='raster', pattern="xx*", quiet=True, flags="f")
    Module("g.remove", type='raster', pattern="zz*", quiet=True, flags="f")
    dem=general.orig_dem
    find_dem_modified = gscript.find_file("zz"+dem+"_modified", element = 'cell')
    if find_dem_modified['name'] != "":
        Module("g.remove", type='raster', name="zz"+dem+"_modified", quiet=True, flags="f")
    find_dem_modified = gscript.find_file("zz"+dem+"_modified_full", element = 'cell')
    if find_dem_modified['name'] != "":
        Module("g.remove", type='raster', name="zz"+dem+"_modified_full", quiet=True, flags="f")
    if main.treesmap:
        find_treesmap = gscript.find_file(main.treesmap, element = 'cell')
        if find_treesmap['name'] != "":
            Module("g.remove", type='raster', name=main.treesmap, quiet=True, flags="f")
    if main.buildmap:
        find_buildmap = gscript.find_file(main.buildmap, element = 'cell')
        if find_buildmap['name'] != "":
            Module("g.remove", type='raster', name=main.buildmap, quiet=True, flags="f")
    #if there is a MASK remove it
    find_MASK = gscript.find_file("MASK", element = 'cell')
    if find_MASK['name'] != "":
        Module("r.mask", flags="r")
    
    
#starting function. needed for preparing the data
def general(pnt, dem, treesmap, buildmap, treesheigh, buildheigh, obsabselev):
    general.orig_dem=dem
    #raster points
    Module("v.to.rast", input=pnt, output="xxrastpnt", type="point", use="val", overwrite=True, quiet=True)
    #altering DEM in case there are buldings and trees map (to obtain a sort of DSM) and in case observer is not flying the dem is kept to the ground level at observer positions
    if treesmap and buildmap:
        Module("v.to.rast", input=treesmap, output=treesmap, use="attr", attribute_column=treesheigh, overwrite=True, quiet=True)
        Module("v.to.rast", input=buildmap, output=buildmap, use="attr", attribute_column=buildheigh, overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = if(isnull({B}),{C},{B})".format(A="zztreesbuildingmap",B=buildmap,C=treesmap), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{B} = if(isnull({A}),{C},{C}+{A})".format(A="zztreesbuildingmap", B="zz"+dem+"_modified_full", C=dem), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{B} = if(isnull({A}),{D},{C})".format(A="xxrastpnt", B="zz"+dem+"_modified",C=dem,D="zz"+dem+"_modified_full"), overwrite=True, quiet=True)
        if obsabselev:
            dem="zz"+dem+"_modified_full"
        else:
            dem="zz"+dem+"_modified"
    elif treesmap and not buildmap: 
        Module("v.to.rast", input=treesmap, output=treesmap, use="attr", attribute_column=treesheigh, overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{B} = if(isnull({A}),{C},{C}+{A})".format(A=treesmap, B="zz"+dem+"_modified_full", C=dem), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{B} = if(isnull({A}),{D},{C})".format(A="xxrastpnt", B="zz"+dem+"_modified",C=dem,D="zz"+dem+"_modified_full"), overwrite=True, quiet=True)
        if obsabselev:
            dem="zz"+dem+"_modified_full"
        else:
            dem="zz"+dem+"_modified"
    elif buildmap and not treesmap:
        Module("v.to.rast", input=buildmap, output=buildmap, use="attr", attribute_column=buildheigh, overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{B} = if(isnull({A}),{C},{C}+{A})".format(A=buildmap, B="zz"+dem+"_modified_full", C=dem), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{B} = if(isnull({A}),{D},{C})".format(A="xxrastpnt", B="zz"+dem+"_modified",C=dem,D="zz"+dem+"_modified_full"), overwrite=True, quiet=True)
        if obsabselev:
            dem="zz"+dem+"_modified_full"
        else:
            dem="zz"+dem+"_modified"
    #preparing the maps of the orientation of the DEM   
    Module("r.slope.aspect", elevation=dem, slope="zzslope", aspect="zzaspect", overwrite=True, quiet=True)
    #evaluation of the azimuth layer 
    Module("r.mapcalc", expression="zzazimuth = (450-zzaspect) - int( (450-zzaspect) / 360) * 360", overwrite=True, quiet=True) 
    #evaluation of the layer of the vertical component of the versor perpendicular to the terrain slope
    Module("r.mapcalc", expression="zzc_dem = cos(zzslope)", overwrite=True, quiet=True)
    #evaluation of the layer of the north component of the versor perpendicular to the terrain slope
    Module("r.mapcalc", expression="zzb_dem = sin(zzslope)*cos(zzazimuth)", overwrite=True, quiet=True)
    #evaluation of the layer of the east component of the versor perpendicular to the terrain slope
    Module("r.mapcalc", expression="zza_dem = sin(zzslope)*sin(zzazimuth)", overwrite=True, quiet=True)    
    #creating a list of the categories of the available points in the input layer
    ret = Module('v.info', flags='t', map=pnt, stdout_=PIPE)
    stats = parse_key_val(ret.outputs.stdout)
    npnt=stats['points']
    ctg2=Module("v.category", flags="g", input=pnt, option="print", type="point", stdout_=PIPE)
    ctg=ctg2.outputs.stdout.splitlines()
    #exporting variables for other functions
    general.ctg=ctg
    general.npnt=npnt
    general.dem=dem

#iterative function. needed for preparing the data for the parallel computation
def iterate(pnt, dem, obs_heigh, maxdist, hcurv, downward, oradius, ctg, nprocs, obsabselev, memory):
    #creating the list which will contain the parameters to be used to run the "compute function in parallel
    tasks = []
    #starting the main loop for each point in the point input layer
    k=0
    #starting the main loop for each point in the point input layer
    for i in ctg:
        #the following i=2 can be uncommented for testing the loop
        #i=2
        k=k+1
        #writing a message on the screen describing which point location is being put in the task list
        message = "Preparing the task for the point %s (cat = %s) of %s points"
        gscript.message(message % (k,i,general.npnt))
        #appending the parameters needed for running the code in parallel to the task list
        tasks.append( (pnt, dem, obs_heigh, maxdist, hcurv, downward, oradius, i, nprocs, obsabselev, memory) )
        #exporting variables needed by other functions
        iterate.tasks = tasks

#the following is needed to definea a lag between the start of the parallel processes. in fact th definition of the temporary region in grass can fail if manny processes set it at the same time 
def init(lock):
    global starting
    starting = lock

#function that is run iterativelly for each pint location. It is used in parallell using multiprocessing.pool 
def compute(pnt, dem, obs_heigh, maxdist, hcurv, downward, oradius, i, nprocs, obsabselev, memory):
    try:
        #the followig lines help to set a delay between a process and the others
        starting.acquire() # no other process can get it until it is released
        threading.Timer(0.1, starting.release).start() # release in a 0.1 seconds
        #using temporary regions (on the same mapset) for the different parallel computations
        gscript.use_temp_region()
        #extracting a point from the map of the locations
        Module("v.extract", input=pnt, output="zzpnt"+i, cats=i, flags="t", overwrite=True, quiet=True)
        #getting goordinates of the point location
        coords=Module("v.to.db", flags="p", map="zzpnt"+i, type="point", option="coor", separator="|", stdout_=PIPE)
        coords=coords.outputs.stdout.splitlines()[1:][0]
        x=float(coords.split("|")[1])
        y=float(coords.split("|")[2])
        z=float(coords.split("|")[3])
        coords=str(x)+","+str(y)
        #get elevation of the terrain at the point location
        querydem=Module("r.what", coordinates=coords.split(), map=dem, stdout_=PIPE)
        obselev=float(querydem.outputs.stdout.split("|")[3])
        #setting the working region around the point location
        Module("g.region", vector="zzpnt"+i)
        region = grasscore.region()
        E = region['e']
        W = region['w']
        N = region['n']
        S = region['s']
        #Module("g.region", flags="a", e=E+maxdist, w=W-maxdist, s=S-maxdist, n=N+maxdist) 
        Module("g.region", align=dem, e=E+maxdist, w=W-maxdist, s=S-maxdist, n=N+maxdist)
        #now we check if the size of the object for which we calculate solid angle in each pixel is equal to half the resolution or is set by the user
        if oradius == 0: 
            circle_radius=region['nsres']/2
        else:
            circle_radius=oradius
        #Executing viewshed analysis
        if obsabselev:
            relative_height = z - obselev
            #message1 = "* Considered elevation of dem/dsm is: %s *"
            #message2 = "* Relative height of observer above the dem/dsm is: %s *"
            #message3 = "* Absolute elevation of observer used in r.viewshed is: %s *"
            #gscript.message(message1 % (str(obselev)) )
            #gscript.message(message2 % (str(relative_height)))
            #gscript.message(message3 % (str(z)))
            if hcurv:
                Module("r.viewshed", input=dem, output="zzview"+i, coordinates=coords.split(), memory=memory, observer_elevation=relative_height, max_distance=maxdist,  flags="c", overwrite=True, quiet=True)
            else:
                Module("r.viewshed", input=dem, output="zzview"+i, coordinates=coords.split(), memory=memory, observer_elevation=relative_height, max_distance=maxdist, overwrite=True, quiet=True)                
            if downward:
                #Since UAV nor Satellite are not expected to see above their level (they are only looking to the ground) vertical angles above 90 are set to null. 
                Module("r.mapcalc", expression="zzview{I} = if(zzview{I}>90 && zzview{I}<180,null(),zzview{I})".format(I=i), overwrite=True, quiet=True)
        else:
            #message1 = "* Considered elevation of dem/dsm is: %s *"
            #message2 = "* Relative height of observer above the dem/dsm is: %s *"
            #message3 = "* Absolute elevation of observer used in r.viewshed is: %s *"
            #gscript.message(message1 % (str(obselev)) )
            #gscript.message(message2 % (str(obs_heigh)))
            #gscript.message(message3 % (str(obselev + obs_heigh)))            
            if hcurv:
                Module("r.viewshed", input=dem, output="zzview"+i, coordinates=coords.split(), memory=memory, observer_elevation=obs_heigh, max_distance=maxdist, flags="c", overwrite=True, quiet=True)
            else:
                Module("r.viewshed", input=dem, output="zzview"+i, coordinates=coords.split(), memory=memory, observer_elevation=obs_heigh, max_distance=maxdist, overwrite=True, quiet=True)
        #Since r.viewshed set the cell of the output visibility layer to 180 under the point, this cell is set to 0.01 (DOUBLE CHECK THIS PART)
        Module("r.mapcalc",expression="zzview{I} = if(zzview{I}==180,0,zzview{I})".format(I=i), overwrite=True, quiet=True)
        #estimating the layer of the horizontal angle between point and each visible cell (angle of the horizontal line of sight)         
        Module("r.mapcalc", expression="{A} = \
            if( y()>{py} && x()>{px}, atan(({px}-x())/({py}-y())),  \
            if( y()<{py} && x()>{px}, 180+atan(({px}-x())/({py}-y())),  \
            if( y()<{py} && x()<{px}, 180+atan(({px}-x())/({py}-y())),  \
            if( y()>{py} && x()<{px}, 360+atan(({px}-x())/({py}-y())), \
            if( y()=={py} && x()>{px}, 90, \
            if( y()<{py} && x()=={px}, 180, \
            if( y()=={py} && x()<{px}, 270, \
            if( y()>{py} && x()=={px}, 0 \
            ) ) ) ) ) ) ) )".format(A='zzview_angle'+i,py=y, px=x), overwrite=True, quiet=True)
        #estimating the layer of the vertical angle between point and each visible cell  (angle of the vertical line of sight) ()
        Module("r.mapcalc", expression="zzview90_{I} = zzview{I} - 90".format(I=i), overwrite=True, quiet=True)
        #evaluate the vertical component of the versor oriented along the line of sight         
        Module("r.mapcalc", expression="zzc_view{I} = sin(zzview90_{I})".format(I=i), overwrite=True, quiet=True)
        #evaluate the northern component of the versor oriented along the line of sight  
        Module("r.mapcalc", expression="zzb_view{I} = cos(zzview90_{I})*cos(zzview_angle{I})".format(I=i), overwrite=True, quiet=True)
        #evaluate the eastern component of the versor oriented along the line of sight  
        Module("r.mapcalc", expression="zza_view{I} = cos(zzview90_{I})*sin(zzview_angle{I})".format(I=i), overwrite=True, quiet=True)    
        #estimate the three-dimensional distance between the point and each visible cell
        if obsabselev:
            Module("r.mapcalc", expression="{D} = pow(pow(abs(y()-{py}),2)+pow(abs(x()-{px}),2)+pow(abs({dtm}-{Z}),2),0.5)".format(D='zzdistance'+i, dtm=dem, Z=z, py=y, px=x), overwrite=True, quiet=True)
        else:
            Module("r.mapcalc", expression="{D} = pow(pow(abs(y()-{py}),2)+pow(abs(x()-{px}),2)+pow(abs({dtm}-({obs}+{obs_h})),2),0.5)".format(D='zzdistance'+i, dtm=dem, obs=obselev, obs_h=obs_heigh, py=y, px=x), overwrite=True, quiet=True)
        
        #estimating the layer of the angle between the versor of the terrain and the line of sight
        Module("r.mapcalc", expression="zzangle{I} = acos((zza_view{I}*zza_dem+zzb_view{I}*zzb_dem+zzc_view{I}*zzc_dem)/(sqrt(zza_view{I}*zza_view{I}+zzb_view{I}*zzb_view{I}+zzc_view{I}*zzc_view{I})*sqrt(zza_dem*zza_dem+zzb_dem*zzb_dem+zzc_dem*zzc_dem)))".format(I=i), overwrite=True, quiet=True)
        #in rare cases the angles may results, erroneusly, less than 90. Setting them to 90
        Module("r.mapcalc", expression="zzangle{I} = if(zzangle{I} > 90, zzangle{I}, 90)".format(I=i), overwrite=True, quiet=True) 
        #filtering 3d distance based on angle{I} map
        Module("r.mapcalc", expression="{D} = if(isnull(zzangle{I}),null(),{D})".format(D="zzdistance"+str(i),I=i), overwrite=True, quiet=True)
        #calculating H1 and H2 that are the distances from the observer to the more distant and less distant points of the inclinded circle representing the pixel
        Module("r.mapcalc", expression="zzH1_{I} = pow(pow({r},2)+pow({d},2)-(2*{r}*{d}*cos(270-zzangle{I})),0.5)".format(r=circle_radius,d="zzdistance"+str(i),I=i), overwrite=True, quiet=True) 
        Module("r.mapcalc", expression="zzH2_{I} = pow(pow({r},2)+pow({d},2)-(2*{r}*{d}*cos(zzangle{I}-90)),0.5)".format(r=circle_radius,d="zzdistance"+str(i),I=i), overwrite=True, quiet=True) 
        #calculating B1 and B2 that are the angles between the line passing through the observer and the center of the pixel and the distant and less distant points of the inclinded circle representing the pixel
        Module("r.mapcalc", expression="zzB1_{I} = acos( (pow({r},2)-pow(zzH1_{I},2)-pow({d},2)) / (-2*zzH1_{I}*{d}) ) ".format(r=circle_radius,d="zzdistance"+str(i),I=i), overwrite=True, quiet=True) 
        Module("r.mapcalc", expression="zzB2_{I} = acos( (pow({r},2)-pow(zzH2_{I},2)-pow({d},2)) / (-2*zzH2_{I}*{d}) ) ".format(r=circle_radius,d="zzdistance"+str(i),I=i), overwrite=True, quiet=True) 
        #calculating solid angle considering that the area of an asimetric ellipse is equal to the one of an ellipse having the minor axis equal to the sum of the tqo unequal half minor axes 
        Module("r.mapcalc", expression="zzsangle{I} = ({pi}*{r}*( {d}*tan(zzB1_{I}) + {d}*tan(zzB2_{I}) )/2 )  / (pow({r},2)+pow({d},2)) ".format(r=circle_radius,d="zzdistance"+str(i),I=i,pi=pi), overwrite=True, quiet=True) 
        #approximations for calculating solid angle can create too much larger values under or very close the position of the oserver. in such a case we assume that the solid angle is half of the visible sphere (2*pi)
        #The same occur when it is used an object_radius that is larger than thepixel size. In some cases ths can produce negative values of zzB2 whit the effect of creating negative values 
        Module("r.mapcalc", expression="zzsangle{I} = if(zzsangle{I}>2*{pi} || zzB2_{I}>=90,2*{pi},zzsangle{I})".format(I=i, pi=pi), overwrite=True, quiet=True)
        #removing temporary region    
        gscript.del_temp_region()
    except:
        #cleaning termporary layers
        #cleanup()
        #message = " ******** Something went wrong: please try to reduce the number of CPU (parameter 'procs') ******* "
        #gscript.message(message)
        #sys.exit()
        f = open("error_cat_"+i+".txt", "x")
        f.write("error in category: "+i)
        f.close()


            
#the following functions is used in parallel to process the combination of the differnt product maps 
def collectresults(task,proc):
    for i in task:
        #print(i)
        #print(proc)
        #print("combining outputs for point " + str(i) + " and for cpu " + str(proc))
        message = "Combining outputs for point %s and for cpu %s "
        gscript.message(message % (str(i),str(proc)))
        #updating the output layer of the angle
        #updating the output layer of the best angle of view among all the points in the path 
        Module("r.mapcalc", expression="{A} = if(isnull({I}) ||| {I}==0,{A},max({A},{I}))".format(A='xxtemp_a_'+str(proc),I='zzangle'+i), overwrite=True, quiet=True)
        #updating the output layer of the category of the point who has the higher angles with the considered cell
        Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}<{Z},{A},{cat}))".format(A='xxtemp_c_'+str(proc),I='zzangle'+i,Z='xxtemp_a_'+str(proc),cat=i), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}<{Z},{A},{I}))".format(A='xxmaxangle_'+str(proc),I='zzangle'+i,Z='xxtemp_a_'+str(proc)), overwrite=True, quiet=True)
	   #updating the output layer of the 3d distance
        #Module("r.mapcalc", expression="{A} = if(isnull({I}),{A}, min({I},{A}))".format(A='xxtemp_e_'+str(proc),I='distance'+i),overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = if(isnull({I}),{A}, if({A} != 0,min({I},{A}),{I}))".format(A='xxtemp_e_'+str(proc),I='zzdistance'+i),overwrite=True, quiet=True)
        #updating the output layer of the category of the point who has the higher solid angles with the considered cell
        Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}>{Z},{A},{cat}))".format(A='xxtemp_h_'+str(proc),I='zzdistance'+i,Z='xxtemp_e_'+str(proc),cat=i), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}>{Z},{A},{I}))".format(A='xxmin3ddistance_'+str(proc),I='zzdistance'+i,Z='xxtemp_e_'+str(proc)), overwrite=True, quiet=True)
        #updating the output layer of the solid angle
        #updating the output layer of the best solid angle of view among all the points in the path 
        Module("r.mapcalc", expression="{A} = if(isnull({I}) ||| {I}==0,{A},max({A},{I}))".format(A='xxtemp_f_'+str(proc),I='zzsangle'+i), overwrite=True, quiet=True)
        #updating the output layer of the category of the point who has the higher solid angles with the considered cell
        Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}<{Z},{A},{cat}))".format(A='xxtemp_g_'+str(proc),I='zzsangle'+i,Z='xxtemp_f_'+str(proc),cat=i), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}<{Z},{A},{I}))".format(A='xxmaxsangle_'+str(proc),I='zzsangle'+i,Z='xxtemp_f_'+str(proc)), overwrite=True, quiet=True)
        #updating the output layer of the number of points a pixel is visible from
        #updating the output layer of the number of points from which a cell is visible
        Module("r.mapcalc", expression="{A} = if(isnull({I}) ||| {I}==0,{A},{A}+1)".format(A='xxtemp_b_'+str(proc),I='zzangle'+i), overwrite=True, quiet=True)

           


def main():
    options, flags = parser()
    ##PARAMETER TO BE UNCOMMENTED TO RUN A TEST OF THE PROGRAM FOR DEBUGGING
    #pnt = "pt50add"
    #dem = "Leintz_dem_50"
    #output = "test"
    #maxdist = 500.0
    #getting options from the command line
    pnt = options['points']
    dem = options['dem']
    output = options['output']
    maxdist = float(options['maxdist'])  
    obs_heigh = float(options['obs_heigh'])
    nprocs = int(options['nprocs'])
    treesmap = options['treesmap']
    buildmap =  options['buildingsmap']
    treesheigh = options['treesheigh']
    buildheigh = options['buildingsheigh']
    obsabselev = options['obsabselev']
    if options['object_radius']:
        oradius = float(options['object_radius'])  
    if options['layer']:
        layer = int(options['layer'])
    if options['viewangle_threshold']:
        viewangle_threshold = float(options['viewangle_threshold'])
    memory=int(options['memory'])
    binary = flags['b']
    hcurv = flags['c']
    downward = flags['d']
    #converting points map to 3d if a layer and elevation field is provided
    try:    
        if obsabselev:
            if layer:
                Module("v.to.3d", input=pnt, output="xxtemppnt3d", type="point", column=obsabselev, layer=layer , overwrite=True, quiet=True)
                pnt="xxtemppnt3d"
    except: 
        message = "There was an error converting the layer to 3d. Plese check if you have provided column and layer information. Exiting "
        gscript.message(message)
        #cleaning termporary layers
        cleanup()
        sys.exit()           
    try:
        oradius
    except:
        oradius=0
    #exporting some variables for other functions
    main.treesmap=treesmap 
    main.buildmap=buildmap
    main.nprocs=nprocs        
    try:
        #setting the starting region alignement to the grid of the DEM
        Module("g.region", align=dem)
        #running the "general" function
        general(pnt, dem,treesmap,buildmap,treesheigh,buildheigh, obsabselev)
        #perparing the container and limiting the number of CPU used
        #message = "Using the following number of CPU: %s"
        #gscript.message(message % str(nprocs))
        #pool = multiprocessing.Pool( nprocs )
        #creating the pool ofr the parallel processes. initialization of the single processes is delayed by the lock 
        pool = multiprocessing.Pool( processes=nprocs , initializer=init, initargs=[multiprocessing.Lock()])
        #Running the "iterate" function        
        iterate(pnt, general.dem, obs_heigh, maxdist, hcurv, downward, oradius, general.ctg, nprocs, obsabselev, memory)
        #running the "compute" function in parallel
        #result=[pool.starmap( compute, [(t) for t in iterate.tasks])]
        pool.starmap( compute, [(t) for t in iterate.tasks])
        #terminating the parallel computation
        pool.close()        
        #using the generated temporary maps to create the final maps
        try:
            #the following to split the categories in nprocs chunks                
            #chunks = list(zip(*[iter(general.ctg)]*int(len(general.ctg)/nprocs)))
            chunks = [general.ctg[x:x+nprocs] for x in range(0, len(general.ctg), nprocs)]
            #creating the tasks for the parallel processing
            tasks2=list(zip(chunks, range(len(chunks))))
            #creating the "zeros" map to be used for the combination of the different teporary maps
            for i in range(len(chunks)):
                #print("creating \"map zeros\" for the chunk "+str(i) + " of " + str(len(chunks)))
                message = "Creating \"maps zero\" for the chunk %s  of %s"
                gscript.message(message % (str(i),str(len(chunks))))
                #for storing maximum view angles
                Module("r.mapcalc", expression="xxtemp_a_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)
                #for storing number of points a pixel is visible from
                Module("r.mapcalc", expression="xxtemp_b_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)
                #for storing the category of the point a pixel is visible with the maximum angle
                Module("r.mapcalc", expression="xxtemp_c_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)
                #for storing the minimum 3ddistance
                Module("r.mapcalc", expression="xxtemp_e_{p} = 0".format(p=str(i)), overwrite=True, quiet=True) 
                #for storing the maximum solid angle
                Module("r.mapcalc", expression="xxtemp_f_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)  
                #for storing the category of the point a pixel is visible with the maximum solid angle
                Module("r.mapcalc", expression="xxtemp_g_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)  
                #for storing the category of the point a pixel is visible with the minimum distance
                Module("r.mapcalc", expression="xxtemp_h_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)     
                #needed for calculating the category of the point a pixel is visible with the maximum angle
                Module("r.mapcalc", expression="xxmaxangle_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)          
                #needed for calculating the category of the point a pixel is visible with the minimum distance
                Module("r.mapcalc", expression="xxmin3ddistance_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)  
                #needed for calculating the category of the point a pixel is visible with the maximum solid angle
                Module("r.mapcalc", expression="xxmaxsangle_{p} = 0".format(p=str(i)), overwrite=True, quiet=True)                             
            #creating the pool for the combination parallel processes
            pool2 = multiprocessing.Pool( nprocs )
            #running the "collectresults" function in parallel
            pool2.starmap( collectresults, [(t) for t in tasks2])
            #terminating the parallel computation
            pool2.close()    
        except:
            #cleaning termporary layers
            cleanup()
            message = " ******** Something went wrong during the entire process of collection of results: please try to reduce the number of CPU (parameter 'procs') ******* "
            gscript.message(message)
            sys.exit()
        else:
            #creating the "zeros" map to be used for the FINAL combination of the different teporary maps
            #print("creating Final \"zeros map \"")
            message = "Creating Final \"zeros map \" "
            gscript.message(message)
            Module("r.mapcalc", expression="xxtemp_a = 0", overwrite=True, quiet=True)
            Module("r.mapcalc", expression="xxtemp_b = 0", overwrite=True, quiet=True)
            Module("r.mapcalc", expression="xxtemp_c = 0", overwrite=True, quiet=True)
            Module("r.mapcalc", expression="xxtemp_e = 0", overwrite=True, quiet=True) 
            Module("r.mapcalc", expression="xxtemp_f = 0", overwrite=True, quiet=True)
            Module("r.mapcalc", expression="xxtemp_g = 0", overwrite=True, quiet=True) 
            Module("r.mapcalc", expression="xxtemp_h = 0", overwrite=True, quiet=True)             
            #combining the maps. This must be done in series since xxtemp_c depends on xxtemp_a for each given chunk
            for i in range(len(chunks)):
                #print("combining chunk " + str(i) + " of " + str(len(chunks)) )
                message = "Combining chunk %s of %s "
                gscript.message(message % (str(i),str(len(chunks))))
                #updating the map of the angles
                Module("r.mapcalc", expression="{A} = if(isnull({I}) ||| {I}==0,{A},max({A},{I}))".format(A='xxtemp_a',I="xxtemp_a_"+str(i)), overwrite=True,  quiet=True)
                #updating the output layer of the number of points from which a cell is visible
                Module("r.mapcalc", expression="{A} = if(isnull({I}) ||| {I}==0,{A},{A}+{I})".format(A='xxtemp_b',I='xxtemp_b_'+str(i)), overwrite=True, quiet=True)
                #updating the output layer of the category of the point who has the higher angle with the considered cell
                Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}<{Z},{A},{C}))".format(A='xxtemp_c',C='xxtemp_c_'+str(i),I='xxmaxangle_'+str(i),Z='xxtemp_a'), overwrite=True,  quiet=True)
     	        #updating the output layer of the 3d distance
                Module("r.mapcalc", expression="{A} = if({A} != 0 && {I} != 0,min({I},{A}), if({A} == 0 && {I} != 0,{I}, if({A} != 0 && {I} == 0,{A} )   )  )".format(A='xxtemp_e',I='xxtemp_e_'+str(i)),overwrite=True, quiet=True)                                  ### TXOMIN ###            
                #updating the map of the solid angles
                Module("r.mapcalc", expression="{A} = if(isnull({I}) ||| {I}==0,{A},max({A},{I}))".format(A='xxtemp_f',I="xxtemp_f_"+str(i)), overwrite=True,  quiet=True)
                #updating the output layer of the category of the point who has the higher solid angle with the considered cell
                Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}<{Z},{A},{C}))".format(A='xxtemp_g',C='xxtemp_g_'+str(i),I='xxmaxsangle_'+str(i),Z='xxtemp_f'), overwrite=True,  quiet=True)
                #updating the output layer of the category of the point who has the lower 3ddistance to the considered cell
                Module("r.mapcalc", expression="{A} = if({I}==0 ||| isnull({I}),{A}, if({I}>{Z},{A},{C}))".format(A='xxtemp_h',C='xxtemp_h_'+str(i),I='xxmin3ddistance_'+str(i),Z='xxtemp_e'), overwrite=True,  quiet=True)        
            #set to null the 0 values
            Module("r.null", map='xxtemp_a', setnull=0, quiet=True)
            Module("r.null", map='xxtemp_b', setnull=0, quiet=True)
            Module("r.null", map='xxtemp_c', setnull=0, quiet=True)
            Module("r.null", map='xxtemp_e', setnull=0, quiet=True)
            Module("r.null", map='xxtemp_f', setnull=0, quiet=True)
            Module("r.null", map='xxtemp_g', setnull=0, quiet=True)
            Module("r.null", map='xxtemp_h', setnull=0, quiet=True)    
        #creating the output layer 
        #if there is a threshold for the viewangles, set a MASK
        try:
            viewangle_threshold
            Module("r.mapcalc", expression="MASK=if({A}>{B},1,null())".format(A='xxtemp_a', B=viewangle_threshold), quiet=True, overwrite=True)
        except:
            pass
        #print("Creating final maps")
        message = "Creating final maps"
        gscript.message(message)
        Module("r.mapcalc", expression="{A} = {B}".format(A=output+'_maxViewAngle',B='xxtemp_a'), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = {B}".format(A=output+'_numberOfViews',B='xxtemp_b'), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = {B}".format(A=output+'_pointOfViewWithMmaxAngle',B='xxtemp_c'), overwrite=True, quiet=True)
        ###Module("r.colors", map=output+'_visibility_index', color='population_dens', quiet=True)
        Module("r.mapcalc", expression="{A} = {B}".format(A=output+'_min3dDistance',B='xxtemp_e'), overwrite=True, quiet=True)
        #add lines of code here for calculating the max solid angle the points from where we have min distance and max solid angle
        Module("r.mapcalc", expression="{A} = {B}".format(A=output+'_maxSolidAngle',B='xxtemp_f'), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = {B}".format(A=output+'_pointOfViewWithMmaxSolidAngle',B='xxtemp_g'), overwrite=True, quiet=True)
        Module("r.mapcalc", expression="{A} = {B}".format(A=output+'_pointOfViewWithMin3dDistance',B='xxtemp_h'), overwrite=True, quiet=True)
        #setting logaritmic colors to some maps
        Module("r.colors", map=output+'_maxSolidAngle', color='blues', flags="g", quiet=True)
        Module("r.colors", map=output+'_maxViewAngle', color='oranges', flags="g", quiet=True)        
        if binary:
            Module("r.mapcalc", expression="{A} = if(isnull({B}),null(),1)".format(A=output+'_binary',B='xxtemp_a'), overwrite=True, quiet=True)            
        #print("Succesful run! .... cleaning temporary files ")
        message = "Succesful run! "
        gscript.message(message)        
    #in case of CTRL-C
    except KeyboardInterrupt:
        cleanup()
                
if __name__ == "__main__":
        options, flags = parser()
        atexit.register(cleanup)
        sys.exit(main())
