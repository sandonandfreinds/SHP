"""
Documentation of the functions created to make images and movies.
"""
#######################################################################################################
import yt  #created by Britton Smith and Megan Lang, https://doi.org/10.21105/joss.01881
import ytree
import math
import gc
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
#######################################################################################################
"""
setup_cam: returns yt objects 'camera' and 'scene' .
'scene' is the volume rendered by yt. 
'camera' is a virtual camera in the scene that can be manipulated by the user.
The compulsory arguments are 'obj', 'field' and 'lens'. 
'obj' is the data source for the volume rendering (normally obj = ds=yt.load(...)).
'field' is the field in the data source the user would like to visualise. i.e. density, mass, temperature etc.
'lens' is a keyword argument that requires one of the lens type keyords from yt: "plane-parallel", "perspective","stereo-perspective","fisheye","spherical",
"stereo-spherical".
"""
def setup_cam(obj, field, lens, resolution=None, orientation=None, frame_width=None, focus_point=None, position=None):
    scene = yt.create_scene(obj, field, lens_type=lens)
    camera = scene.camera
    
    if position is None: #position means the locataion of the virtual camera in the scene
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

#######################################################################################################
"""
take_image: produces an image of a 'scene' created by setup_cam() and exports it as a .png file to a user specified location.
'image' is an integer counter to keep track of images made. 
"""
def take_image(scene, image):
    global export_path
    source = scene[0]
    source.tfh.set_bounds()
    source.tfh.set_log(True)
    scene.save(export_path+"/frame_%04d.png" %image, sigma_clip=2)
    image += 1
    return image 

#######################################################################################################
"""
take_image: produces an annotated image of a 'scene' created by setup_cam() and exports it as a .png file to a user specified location.
The image will be annotated with a plot of the color transfer fucntion on the right and the time in the simulation of the scene.
'image' is an integer counter to keep track of images made. 
"""
def take_annotated_image(scene, image):
    source = scene[0]
    source.tfh.set_bounds()
    source.tfh.set_log(True)
    text_string = f"T = {float(ds.current_time.to('Gyr'))} Gyr"
    scene.save_annotated(export_path+"/annontated_frame_%04d.png", sigma_clip=6, text_annotate=[[(0.1, 0.95), text_string]])
    image += 1
    return image

#######################################################################################################
"""
smoothed_pos: returns array of smoothed x,y,z coordinates of the halo at each snapshot time. ([[x,y,z]
                                                                                                [....]]) 
These will set the camera's focus point.
r = halo posititions from the merger tree 
t = progenitor halo time times from merger tree
t_snap = times at which the snapshot are taken
"""
def smoothed_pos(r, t, t_snap):  
    smooth_r = np.empty((len(t_snap),0), int)
    for i in range(3):
        smooth_ri_func =  interpolate.UnivariateSpline(t, r[:,i], k=3)
        smooth_ri = smooth_ri_func(t_snap)
        smooth_ri = np.vstack(smooth_ri)
        smooth_r = np.append(smooth_r, smooth_ri, axis=1)
    return smooth_r
########################################################################################################
"""
smoothed_rad: returns array of smoothed virial radii at each snapshot, ([...], "kpc"). 
These will set the camera's frame width. 
rad = halo radii from the merger tree 
t = progenitor halo time times from merger tree
t_snap = times at which the snapshot are taken
"""
def smoothed_rad(rad, t, t_snap):
    s_rad = rad[rad.argmax()]*len(t) #smoothing factor 
    smooth_rad_func = interpolate.UnivariateSpline(t, rad, k=3, s = s_rad)
    smooth_rad = smooth_rad_func(t_snap) 
    smooth_rad = a.arr(smooth_rad, 'kpc') #merger tree radii are in 'kpc', making sure units remain consistent
    return smooth_rad
#######################################################################################################
"""
initial_conditions_index: returns an integer that allows initial conditions for the camera to be set
snapshots start before first halos form, need to find the snapshot that corresponds to the first halo in the tree.
Then at all times prior to the first halo forming the camera have the corect frame width and focus point. 

prog_list = progenitor time list/array
snap_list = snapshot times
"""
def initial_conditions_index(prog_list, snap_list):
    found = False
    if snap_list[0]-prog_list[0] < 0: #if first snap shot takes place before first halo progentior
        for t in range(0,len(snap_list)):
            for j in range(0,len(prog_list)):
                if snap_list[t]-prog_list[j] > 0 and found == False: # when this becomes positive, the snapshot corresponding to the first progenitor has been found
                    index = t #function returns this index to set initial conditions 
                    found = True
                    break
    elif snap_list[0]-prog_list[0] > 0: #if first snap shot takes place after first halo progentior
        for t in range(0,len(snap_list)):
            for j in range(0,len(prog_list)):
                if snap_list[t]-prog_list[j] < 0 and found == False:
                    index = t
                    found = True
                    break
    return index
#######################################################################################################
"""
set_frame_widths: returns array of frame widths for the camera for each snapshot. ([...], "kpc") 
rad = halo radii from the merger tree 
t = progenitor halo time times from merger tree
t_snap = times at which the snapshot are taken
"""
def set_frame_widths(rad, t, t_snap):
    smooth_rad = smoothed_rad(rad, t, t_snap)
    index = initial_conditions_index(t, t_snap)
    width = []
    
    for j in range(0,len(t_snap)):
        if j < index: #frame widths for snapshots taken before first progenitor forms get set to the width of the first progenitor
            width.append((10*smooth_rad[index])) 
        else:
            width.append((10*smooth_rad[j]))
    #factor of 10 is to make sure width is always greater than the difference between smoothed-real position of halo
    width = ds.arr(width, 'kpc') #merger tree radii are in 'kpc', making sure units remain consistent
    return width
#######################################################################################################
"""
set_focus_points: returns array of focus points for the camera for each snapshot. ([[x,y,z]
                                                                                    [....]], "unitary") 
rad = halo radii from the merger tree 
t = progenitor halo time times from merger tree
t_snap = times at which the snapshot are taken
"""
def set_focus_points(r, t, t_snap):
    smooth_pos = smoothed_pos(r, t, t_snap)
    index = initial_conditions_index(t, t_snap)
    camera_pos = np.zeros(shape=(len(t_snap), 3))
    for j in range(0,len(t_snap)):
        if j < index: #focus points for snapshots taken before first progenitor forms get set to the position of the first progenitor
            camera_pos[j,:] = smooth_pos[index,:]
        else:
            camera_pos[j,:] = smooth_pos[j,:]
            
    return ds.arr(camera_pos, "unitary") #merger tree positions are in 'unitary', making sure units remain consistent
#######################################################################################################
"""
initial_camera_position: returns a position array ([x,y,z], "code_length")
halo will be placed in a spherical region centered on the halo position with radius 10xframe_width to create a clearer image.
This function places the camera in this sphere and stops the edge of the camera going beyond the box boundary. As yt volume rendering does not obey periodicity.
'obj' is the data source for the volume rendering (normally obj = ds=yt.load(...)).
sphere = obj.sphere(focus_points[0], 100*frame_widths[0]) (this is done in the main code)
"""
def initial_camera_position(sphere, obj, focus_point):
   #create a vector from the halo to the center of the box
    vec=focus_point-obj.arr([0.5,0.5,0.5], "unitary")
    vec = vec.to("code_length")
    mag = math.sqrt(np.dot(vec,vec))
    n = vec/mag # normalised [halo --> center] vector 
    v=obj.arr((sphere.radius*n), "code_length") #scale the normalised vector to desired distance from the object. on 
    init_pos = sphere.center-v #this creates the position at which camera will be placed
    
    #if the width of the frame goes past the box boundaries, then the position componenets are adjusted to avoid this 
    for i in range(3):
        if init_pos[i]+0.1*sphere.radius > obj.arr(100, "code_length"): #sphere radius = 10 x Frame width
            init_pos[i]= init_pos[i]-0.1*sphere.radius
        if init_pos[i]-0.1*sphere.radius < obj.arr(0, "code_length"):
            init_pos[i]= init_pos[i]+0.1*sphere.radius
        
    return init_pos

