# -*- coding: utf-8 -*-
"""
Created on Mon May 16 16:14:12 2022

@author: s2132627
"""
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d


"""
######################
##### User Input #####
######################
"""
unit = "mm"
input1 = krigarr1 # Dataframe or numpy array
input2 = krigarr2 # Dataframe or numpy array
Geofile_number_of_points = 0 # In order for the points making up the surface to have the right index in the .geo and .gli files, they need to start from 1+ how many points are already in those files  
Geofile_number_of_lines = 0
Geofile_number_of_lineloops = 0
Geofile_number_of_planesurfaces = 0
Geofile_number_of_physicalsurfaces = 0
Geofile_number_of_volumes = 0
Geofile_number_of_physicalvolumes = 0

""" Set Variables """
count_P = Geofile_number_of_points
count_L = Geofile_number_of_lines
count_LL = Geofile_number_of_lineloops
count_PlS = Geofile_number_of_planesurfaces
count_PhS = Geofile_number_of_physicalsurfaces
count_V = Geofile_number_of_volumes
count_PhV = Geofile_number_of_physicalvolumes

""" Convert to meters """   # Convert to meters
if unit =="mm":
    arr1b = input1 * 10**(3)
    arr2b = input2 * 10**(3)

""" Centralise """  
# Centralise
arr1xavg = arr1b[:,0].mean();
arr1yavg = arr1b[:,1].mean();
arr1zavg = arr1b[:,2].mean();
arr2xavg = arr2b[:,0].mean();
arr2yavg = arr2b[:,1].mean();
arr2zavg = arr2b[:,2].mean();
print(f'arr1xavg: {arr1xavg}, arr1yavg: {arr1yavg}, arr1zavg: {arr1zavg}, arr2xavg: {arr2xavg}, arr2yavg: {arr2yavg}, arr2zavg: {arr2zavg}')
arr1b[:,0] = arr1b[:,0] - arr1xavg; arr1b[:,1] = arr1b[:,1] - arr1yavg; arr1b[:,2] = arr1b[:,2] - arr1zavg;
arr2b[:,0] = arr2b[:,0] - arr2xavg; arr2b[:,1] = arr2b[:,1] - arr2yavg; arr2b[:,2] = arr2b[:,2] - arr2zavg;

""" Convex Hull 1 """ # Identify the points that "wrap" each of the surfaces
hull1 = ConvexHull(arr1b[:,0:2]) # To get the convex hull of the boundary polygon in 2D, we can only use XY data, not Z [:,0:2] meaning all rows and columns from 0 to 2, where 2 not inclusive. That will give us the indexes of the convex hull points which we can then plot using the full XYZ data
hull1points = np.take(arr1b, hull1.vertices, 0);
hull1points = np.append(hull1points, [hull1points[0,:]], axis = 0)
arr1r = np.delete(arr1b, hull1.vertices, axis=0) # Remove the rows of points pertaining to the convex hull indices
# Plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(hull1points[:,0], hull1points[:,1], hull1points[:,2], 'r--', lw=2) # adds (duplicates) the first point of the convex hull at the end of the sequence in order to close it when plotted.
ax.scatter(arr1b[:,0], arr1b[:,1], arr1b[:,2], 'ro', s=0.01)
ax.set_title(f"Input1")
ax.azim = 0
ax.elev = 90
ax.dist = 10
plt.show()


""" Convex Hull 2 """ # Identify the points that "wrap" each of the surfaces
hull2 = ConvexHull(arr2b[:,0:2]) # To get the convex hull of the boundary polygon in 2D, we can only use XY data, not Z [:,0:2] meaning all rows and columns from 0 to 2, where 2 not inclusive. That will give us the indexes of the convex hull points which we can then plot using the full XYZ data
hull2points = np.take(arr2b, hull2.vertices, 0);
hull2points = np.append(hull2points, [hull2points[0,:]], axis = 0)
arr2r = np.delete(arr2b, hull2.vertices, axis=0) # Remove the rows of points pertaining to the convex hull indices
# Plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(hull2points[:,0], hull2points[:,1], hull2points[:,2], 'r--', lw=2) # adds (duplicates) the first point of the convex hull at the end of the sequence in order to close it when plotted.
ax.scatter(arr2b[:,0], arr2b[:,1], arr2b[:,2], 'ro', s=0.01)
ax.set_title(f"Input2")
plt.show()



"""write to .geo"""
# open other file in write mode 
output_file  = open(basewd + "OGS_Fracture_GEO" + ".geo", 'w') 
output_file.write("Mesh.MshFileVersion = 2.2;\n")
output_file.write("lc = 0.01;\n")
output_file.write("\n")

#############
# Surface 1 #
#############
surf1_lines = list([])
output_file.write( "\n// Surf1 Convex Hull Points\n")
# Create surf1 points
for p, index in enumerate(hull1.vertices):
    count_P = count_P + 1
    if p == 0:
        surf1_bound_start = count_P
    point = arr1b[index]
    output_file.write( f"Point({count_P}) = " + "{" + str(point[0]) + ", " + str(point[1]) + ", " + str(point[2]) + ", lc};\n" )
surf1_bound_stop = count_P

output_file.write( "\n// Surf1 Remaining Points\n")
for p, point in enumerate(arr1b):
    count_P = count_P + 1
    if p == 0:
        surf1_insurfp_start = count_P # Index of first point to add to surface
    output_file.write( f"Point({count_P}) = " + "{" + str(point[0]) + ", " + str(point[1]) + ", " + str(point[2]) + ", lc};\n" )
surf1_insurfp_stop = count_P    # Index of last point to add to surface

# Create surf1 list of sequence of lines point pairs (list) for the LineLoop
for first, second in zip(range(surf1_bound_start, surf1_bound_stop+1), range(surf1_bound_start+1, surf1_bound_stop+1)):
    surf1_lines.append([first, second])
surf1_lines.append([surf1_bound_stop, surf1_bound_start])

# Create surf1 lines that will compose the LineLoop
for l, line in enumerate(surf1_lines):
    count_L = count_L + 1    
    if l == 0:
        surf1_LL_start = count_L    # Saves the line number for the LineLoop
    output_file.write(f"Line ({count_L}) = " + "{" + str(line[0]) + ", " + str(line[1]) + "};\n" )
surf1_LL_stop = count_L

# Create surf1 LineLoop
count_LL = count_LL + 1
output_file.write(f"Line Loop ({count_LL}) = " + "{" + ', '.join([str(i) for i in range(surf1_LL_start, surf1_LL_stop + 1)]) + "};\n" ) 

#Create Plane Surface 1
count_PlS = count_PlS + 1
output_file.write( "Plane Surface(" + str(count_PlS) + ")={" + str(count_LL) + "};\n" ) # Plane Surface(2)={2};

# Add surf1 Remaining points to surface 1
for p in range(surf1_insurfp_start, surf1_insurfp_stop+1):
    output_file.write("Point{" + str(p) + "}" + " In Surface{" + str(count_PlS) + "};\n" ) # Point{9} In Surface{2};

#Create Physical Surface 1
count_PhS = count_PhS + 1
output_file.write( "Physical Surface(" + str(count_PhS) + ")={" + str(count_PlS) + "};\n" ) # Physical Surface(1) = {1}; 

 

###############################################
################ Separator#####################
for i in range(100):
    output_file.write("\n//\n")
###############################################

#############
# Surface 2 #
#############
surf2_lines = list([])
output_file.write( "\n// Surf2 Convex Hull Points\n")
# Create surf2 points
for p, index in enumerate(hull2.vertices):
    count_P = count_P + 1
    if p == 0:
        surf2_bound_start = count_P
    point = arr2b[index]
    output_file.write( f"Point({count_P}) = " + "{" + str(point[0]) + ", " + str(point[1]) + ", " + str(point[2]) + ", lc};\n" )
surf2_bound_stop = count_P

output_file.write( "\n// Surf2 Remaining Points\n")
for p, point in enumerate(arr2b):
    count_P = count_P + 1
    if p == 0:
        surf2_insurfp_start = count_P # Index of first point to add to surface
    output_file.write( f"Point({count_P}) = " + "{" + str(point[0]) + ", " + str(point[1]) + ", " + str(point[2]) + ", lc};\n" )
surf2_insurfp_stop = count_P    # Index of last point to add to surface

# Create surf2 list of sequence of lines point pairs (list) for the LineLoop
for first, second in zip(range(surf2_bound_start, surf2_bound_stop+1), range(surf2_bound_start+1, surf2_bound_stop+1)):
    surf2_lines.append([first, second])
surf2_lines.append([surf2_bound_stop, surf2_bound_start])

# Create surf2 lines that will compose the LineLoop
for l, line in enumerate(surf2_lines):
    count_L = count_L + 1    
    if l == 0:
        surf2_LL_start = count_L    # Saves the line number for the LineLoop
    output_file.write(f"Line ({count_L}) = " + "{" + str(line[0]) + ", " + str(line[1]) + "};\n" )
surf2_LL_stop = count_L

# Create surf2 LineLoop
count_LL = count_LL + 1
output_file.write(f"Line Loop ({count_LL}) = " + "{" + ', '.join([str(i) for i in range(surf2_LL_start, surf2_LL_stop + 1)]) + "};\n" ) 


#Create Plane Surface 2
count_PlS = count_PlS + 1
output_file.write( "Plane Surface(" + str(count_PlS) + ")={" + str(count_LL) + "};\n" ) # Plane Surface(2)={2};

# Add surf2 Remaining points to surface 2
for p in range(surf2_insurfp_start, surf2_insurfp_stop+1):
    output_file.write("Point{" + str(p) + "}" + " In Surface{" + str(count_PlS) + "};\n" ) # Point{9} In Surface{2};

#Create Physical Surface 2
count_PhS = count_PhS + 1
output_file.write( "Physical Surface(" + str(count_PhS) + ")={" + str(count_PlS) + "};\n" ) # Physical Surface(2) = {2}; 

#output_file.write("\n#STOP\n")
output_file.close() #Close output file to submit changes and open again to read


