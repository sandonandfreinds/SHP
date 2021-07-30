import yt
import math
import gc
import numpy as np

"""
This code outputs .png files which can then be turned into a movie of the simulation evolving over time using: 

    ffmpeg -i frame_%04d.png -vcodec libx264 -vf scale=640:-2,format=yuv420p movie.mp4

in the export_path directory. This will produce movie file called movie.mp4 
"""

x = np.array([1,0,0]) #x axis
y = np.array([0,1,0]) #y axis
z = np.array([0,0,1]) #z axis

frame = 0 # frame counter 


"""
########################## User defined values #####################
"""
Field = ('grid', 'nbody_mass') #user defined field 
Lens = 'perspective' 
#lens type for volume rendering using yt keywords: "plane-parallel", "perspective","stereo-perspective","fisheye","spherical","stereo-spherical"
Resolution = (2048,2048) #camera resolution (pixels)

snapshot_data_path = "/disk12/brs/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_*_covering_grid.h5"
export_path = "/home/awood/full_sim_movie"

#optional variables, set to 'None' if you do not want to use.
#all array-like arguments must be unitary.
Orientation = None #can be set to x,y,z or some more funky np.array
Width = None #default is just domain width
Mass_cutoff = None #multiple of mean mass in cell for a given frame i.e. 0.5*mean_mass, helps clear up images.
Focus_point = [0.5,0.5,0.5] #focus point of camera, default is the domain center
inital_position = [0,0.5,0.5] #inital position of camera. 

rotate = True #can be set to False if user wants static image
evolution_angle = (2*np.pi)/4 #how many radians user wants to rotate camera through as system evolves
z0_angle = (2*np.pi)/4 # if user wants to see rotation at z=0 to observe final structure, set to 0 if not
"""
##################################################################
"""

def take_image(sc, image, bounds=None):
    global export_path
    source = sc[0]
    if bounds is None:
        source.tfh.set_bounds()
    else:
        source.tfh.set_bounds(bounds)
    source.tfh.set_log(True)
    source.tfh.grey_opacity = False
    sc.save(export_path+"/frame_%04d.png" %image, sigma_clip=2)
    image += 1
    return image

def setup_cam(obj, field, lens, resolution=None, orientation=None, frame_width=None, focus_point=None, position=None):
    scene = yt.create_scene(obj, field, lens_type=lens)
    camera = scene.camera
    
    if position is None:
        camera.set_position(obj.arr([0,0,0],"unitary"))
    else:
        camera.set_position(position)
    if orientation is None and focus_point is None:
        camera.set_focus(obj.domain_center)
    elif orientation is None:
        camera.set_focus(focus_point)
    elif focus_point is None:
        camera.set_focus(obj.domain_center)
        camera.switch_orientation(north_vector=orientation)
    else:
        camera.set_focus(focus_point)
        camera.switch_orientation(north_vector=orientation)
    if resolution is None:
        camera.set_resolution(420,420)
    else:
        camera.set_resolution(resolution)
    if frame_width is None:
        camera.set_width(obj.domain_width)
    else:
        camera.set_width(frame_width)
    return camera, scene


"""
this is where it all comes together
"""

ts = yt.DatasetSeries(snapshot_data_path)
steps = int(310/len(ts)) #this is to create enough frames for ~30fps to make a decent movie.
###############################################################################################
for fn in ts.outputs:
    ds = yt.load(fn)
    
    if Mass_cutoff != None:
        mass = ds.r[('grid', 'nbody_mass')]
        mean_mass = mass.mean()
        maximum_mass = mass.max()
        Bounds = (Mass_cutoff*mean_mass, maximum_mass)
    
    #inital conditions are implimented here    
    if frame == 0:
        pos = ds.arr(inital_position, "unitary")
        focus = ds.arr(Focus_point, "unitary")
        cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
        Orientation = cam.north_vector
        
    #camera will rotate through a fraction of the evolution angle for each snapshot
    if rotate == True and ds.current_redshift !=0:
        cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
        for _ in cam.iter_rotate((evolution_angle/len(ts)), steps, cam.north_vector, cam.focus):
            frame = take_image(sc, frame) #taking an image at each incriment helps to lengthen the movie if there are few snapshots
        pos = cam.position
        
    #images taken at z=0 if user wants to see final structure from different angles
    #this section will take the same amount of time in the movie as the evolving section for simplicity 
    elif rotate==True and ds.current_redshift == 0 and z0_angle > 0:
        cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
        for _ in cam.iter_rotate(z0_angle, steps*len(ts), cam.north_vector, cam.focus):
            frame = take_image(sc, frame)
        pos = cam.position
    
    elif rotate==False:
        cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
        frame = take_image(sc, frame)
    
    #deleting and calling the 'garbage collect' function prevents the code from using up too much RAM
    del sc 
    del ds
    gc.collect()
#################################################################################################
ts = ts[::-1] # flip time order backward to see evolution in reverse
for fn in ts.outputs:
    ds = yt.load(fn)
    cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
    if rotate == True:
        for _ in cam.iter_rotate((evolution_angle/len(ts)), steps, cam.north_vector, cam.focus):
            frame = take_image(sc, frame)
        pos = cam.position
    
    elif rotate==False:
        cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
        frame = take_image(sc, frame)
    
    del sc
    del ds
    gc.collect()
################################################################################################
ts = ts[::-1] # flip time order forward to see evolution again from a new angle
for fn in ts.outputs:
    ds = yt.load(fn)     
    if rotate == True and ds.current_redshift !=0:
        cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
        for _ in cam.iter_rotate((evolution_angle/len(ts)), steps, cam.north_vector, cam.focus):
            frame = take_image(sc, frame)
        pos = cam.position
        
    elif rotate==True and ds.current_redshift == 0 and z0_angle > 0:
        cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
        for _ in cam.iter_rotate(z0_angle, steps*len(ts), cam.north_vector, cam.focus):
            frame = take_image(sc, frame)
        pos = cam.position
        
    elif rotate==False:
        cam, sc = setup_cam(ds, Field, Lens, resolution=Resolution, orientation=Orientation, frame_width=Width, focus_point=focus, position=pos)
        frame = take_image(sc, frame)
    
    del sc
    del ds
    gc.collect()


    
