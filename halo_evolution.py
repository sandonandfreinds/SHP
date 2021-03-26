import yt
import ytree
import math
import gc
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

"""
This code outputs .png files which can then be turned into a movie of a halo evolving over time using: 

    ffmpeg -i frame_%04d.png -vcodec libx264 -vf scale=640:-2,format=yuv420p halo_movie.mp4

in the export_path directory. This will produce movie file called halo_movie.mp4 
"""

frame = 0 #frame counter

"""
########################## User defined values #####################
"""
field = ('grid', 'nbody_mass') #user defined field 
lens = 'perspective' #lens type for volume rendering using yt keywords
res = (2048,2048) #camera resolution (pixels)

halo_data_path ="/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar_groups/rockstar_groups.h5"
snapshot_data_path = "/disk12/brs/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_*_covering_grid.h5"
export_path = "/home/awood/halo_evolution_movie/"

Orientation = None #orientation of camera, can be set to custom np.array([]) if user wants sepcific oreintation
halo_index = None #If user would like to specify a halo to film. Default is the largest mass halo as this shows the most mergers. 
"""
##################################################################
"""

a = ytree.load(halo_data_path) #load merger tree data 

if halo_index is None:  
    k = a['mass'].argmax() #select largest halo to focus on
    mytree = a[k]
else:
    mytree = a[halo_index]


#relevant values of progenitor halos
prog_pos = mytree['prog', 'position'] 
prog_mass = mytree['prog', 'mass']
prog_rad = mytree['prog', 'virial_radius']
prog_time = mytree['prog','time']
prog_scalefactor = mytree['prog', 'scale_factor']

prog_time = prog_time[::-1] #list had to be reversed to spline
prog_mass = prog_mass[::-1] #flip mass list to match
prog_radius = prog_rad[::-1] #same here

halo_pos = np.empty((len(prog_pos),0), int)#position array x,y,z(0,1,2) in columns
for i in range(3):
     ei = prog_pos[:,i]
     ei = ei[::-1]      #flipping same as above
     ei = np.vstack(ei) #make sure it's a column
     halo_pos=np.append(halo_pos, ei, axis=1)


def take_image(ds, sc, image):
    source = sc[0]
    source.tfh.set_bounds()
    source.tfh.set_log(True)
    sc.save(export_path+"frame_%04d.png" %image, sigma_clip=2)
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

#function to return array of smoothed positions at all times for the camera, t = progen halo time array, 
def smoothed_pos(r, t, t_snap):  #t_snap = times at which the snapshot is being taken 
    smooth_r = np.empty((len(t_snap),0), int)
    for i in range(3):
        smooth_ri_func =  interpolate.UnivariateSpline(t, r[:,i], k=3)
        smooth_ri = smooth_ri_func(t_snap)
        smooth_ri = np.vstack(smooth_ri)
        smooth_r = np.append(smooth_r, smooth_ri, axis=1)
    return smooth_r

#smoothed radius of halo will set with for camera
def smoothed_rad(rad, t, t_snap):
    s_rad = rad[rad.argmax()]*len(t) #smoothing factor 
    smooth_rad_func = interpolate.UnivariateSpline(t, rad, k=3, s = s_rad)
    smooth_rad = smooth_rad_func(t_snap) 
    smooth_rad = a.arr(smooth_rad, 'kpc')
    return smooth_rad

#snap times go further back than halos, need to find first halo point to set conditions for earlier frames 
def initial_conditions_index(prog_list, snap_list):
    found = False
    if snap_list[0]-prog_list[0] < 0:
        for t in range(0,len(snap_list)):
            for j in range(0,len(prog_list)):
                if snap_list[t]-prog_list[j] > 0 and found == False:
                    index = t
                    found = True
                    break
    elif snap_list[0]-prog_list[0] > 0:
        for t in range(0,len(snap_list)):
            for j in range(0,len(prog_list)):
                if snap_list[t]-prog_list[j] < 0 and found == False:
                    index = t
                    found = True
                    break
    return index


#frame widths are set at 10 times the virial radius in order to see halos merge 
def set_frame_widths(rad, t, t_snap):
    smooth_rad = smoothed_rad(rad, t, t_snap)
    index = initial_conditions_index(t, t_snap)
    width = []
    
    for j in range(0,len(t_snap)):
        if j < index:
            width.append((10*smooth_rad[index]))
        else:
            width.append((10*smooth_rad[j]))
    
    width = ds.arr(width, 'kpc')
    return width


#focus points for the camera at each snapshot  defined by the smoothed position of the halo
def set_focus_points(r, t, t_snap):
    smooth_pos = smoothed_pos(r, t, t_snap)
    index = initial_conditions_index(t, t_snap)
    camera_pos = np.zeros(shape=(len(t_snap), 3))
    for j in range(0,len(t_snap)):
        if j < index:
            camera_pos[j,:] = smooth_pos[index,:]
        else:
            camera_pos[j,:] = smooth_pos[j,:]
            
    return ds.arr(camera_pos, "unitary")

#halo will be placed in a spherical object to create a clearer image. 
#this function places the camera in this sphere and stops the edge of the camera going beyond the box boundary
#obj referes to the original data object, ds in this case
def initial_camera_position(sphere, obj, focus_point):
   #create a vector from the halo to the center of the box
    vec=focus_point-obj.arr([0.5,0.5,0.5], "unitary")
    vec = vec.to("code_length")
    mag = math.sqrt(np.dot(vec,vec))
    n = vec/mag # normalised [halo --> center] vector 
    v=obj.arr((sp.radius*n), "code_length") #scale the normalised vector to desired distance from the object. on 
    init_pos = sp.center-v #this creates the position at which camera will be placed
    
    #if the width of the frame goes past the box boundaries, then the position componenets are adjusted to avoid this 
    for i in range(3):
        if init_pos[i]+0.1*sphere.radius > obj.arr(100, "code_length"): #sphere radius = 10 x Frame width
            init_pos[i]= init_pos[i]-0.1*sphere.radius
        if init_pos[i]-0.1*sphere.radius < obj.arr(0, "code_length"):
            init_pos[i]= init_pos[i]+0.1*sphere.radius
        
    return init_pos


"""
now we load the snapshots, first step if finding the correct time stamps and set conditions for camera at each time
"""
ts = yt.DatasetSeries(snapshot_data_path)

snap_time = []# create snap shot time series to for halo widths
for ds in ts:
    snap_time.append(ds.current_time.to("Myr")) #ytree puts halo times in Myrs, same time units are necessary for spline 


frame_widths = set_frame_widths(prog_radius, prog_time, snap_time)
focus_points = set_focus_points(halo_pos, prog_time, snap_time)

"""
this is where it all comes together
"""
steps = int(310/len(ts)) # this is to create enough frame for ~30fps to make a decent movie.

i = 0
for fn in ts.outputs:
    ds = yt.load(fn)
    if frame == 0:
        sp = ds.sphere(focus_points[0], 100*frame_widths[0])
        initial_pos = initial_camera_position(sp,ds,focus_points[0])
        cam, sc = setup_cam(sp, field, lens, resolution=res, frame_width=(0.1*sp.radius), focus_point=sp.center, position=initial_pos)
        cam.roll(np.pi)
        for _ in cam.iter_rotate((np.pi/8), 1, rot_vector=cam.north_vector,rot_center=sp.center):
            pos = cam.position
        cam.roll(np.pi)
        orient = cam.north_vector
        
        cam, sc = setup_cam(sp, field, lens, resolution=res, orientation=orient, frame_width=(0.1*sp.radius), focus_point=sp.center, position=pos)
        for _ in cam.iter_rotate(((np.pi/12)/len(ts)), steps, rot_vector=orient, rot_center=sp.center):
            frame = take_image(sp,sc,frame)
        pos = cam.position
    
    else:
        sp = ds.sphere(focus_points[i], 100*frame_widths[i])
        pos = cam.position
        cam, sc = setup_cam(sp, field, lens, resolution=res, orientation=orient, frame_width=(0.1*sp.radius), focus_point=sp.center, position=pos)
        for _ in cam.iter_rotate(((np.pi/12)/len(ts)), steps, rot_vector=orient, rot_center=sp.center):
            frame=take_image(sp, sc, frame)
        pos = cam.position
    del sc
    del ds
    gc.collect()
    i+=1

i-=1
ts = ts[::-1]
for fn in ts.outputs:
    ds = yt.load(fn)
    if i == (len(ts)-1):
        ds = yt.load(fn)
        sp = ds.sphere(focus_points[i], 100*frame_widths[i])
        pos = cam.position
        cam, sc = setup_cam(sp, field, lens, resolution=res, orientation=orient, frame_width=(0.1*sp.radius), focus_point=sp.center, position=pos)
        for _ in cam.iter_rotate(((np.pi/12)/len(ts)), steps, rot_vector=orient, rot_center=sp.center):
            frame=take_image(sp, sc, frame)
        pos = cam.position
    else:
        sp = ds.sphere(focus_points[i], 100*frame_widths[i])
        cam, sc = setup_cam(sp, field, lens, resolution=res, orientation=orient, frame_width=(0.1*sp.radius), focus_point=sp.center, position=pos)
        for _ in cam.iter_rotate(((np.pi/12)/len(ts)), steps, rot_vector=orient, rot_center=sp.center):
            frame=take_image(sp, sc, frame)
        pos = cam.position
    del sc
    del ds
    gc.collect()
    i-=1

i=0
ts = ts[::-1]
for fn in ts.outputs:
    ds = yt.load(fn)
    if i == 0:
        sp = ds.sphere(focus_points[i], 100*frame_widths[i])
        cam, sc = setup_cam(sp, field, lens, resolution=res, orientation=orient, frame_width=(0.1*sp.radius), focus_point=sp.center, position=pos)
        for _ in cam.iter_rotate(((np.pi/12)/len(ts)), steps, rot_vector=orient, rot_center=sp.center):
            frame=take_image(sp, sc, frame)
        pos = cam.position
    
    else:
        sp = ds.sphere(focus_points[i], 100*frame_widths[i])
        cam, sc = setup_cam(sp, field, lens, resolution=res, orientation=orient, frame_width=(0.1*sp.radius), focus_point=sp.center, position=pos)
        for _ in cam.iter_rotate(((np.pi/12)/len(ts)), steps, rot_vector=orient, rot_center=sp.center):
            frame=take_image(sp, sc, frame)
        pos = cam.position
    del sc
    del ds
    gc.collect()
    i+=1