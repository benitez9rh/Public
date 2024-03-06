# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:52:00 2022

@author: s2132627
"""
import os
import pathlib
from scipy import stats
from scipy.stats import zscore # imports the normal score method used to get rid of the outliers in the data
from scipy.stats import lognorm, norm
from numpy import *
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import math
# import geostatspy.GSLIB as GSLIB                                  # Geostatspy is always giving trouble importing and impoting numba so I simply copied the
# import geostatspy.geostats as geostats                            # functions I needed directly to this script
import geostatspynscore                                             # Copied function and function dependencies from geostatspy GitHub
import time
from time import ctime
from scipy.spatial import ConvexHull, convex_hull_plot_2d
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
extension_txt = ".txt"

def duration():
    finish = time.time()
    days = math.floor( (finish-stop0)/86400)
    hours = math.floor( (finish-stop0)/3600 - days*24 )
    minutes = math.floor( (finish-stop0)/60 - (days*24+hours)*60)
    seconds = math.floor( (finish-stop0) - ((days*24+hours)*60+minutes)*60 )
    print(f'\n**Duration:**\ndays: {days} \nhours: {hours} \nminutes: {minutes} \nseconds: {seconds}')
def isanydigit(n: str) -> bool: # Tests if is a digit because isdigit() doesn't recognise floats
    try:
        float(n)
        n.isdigit()
        return True
    except ValueError:
        return False

def convexhull(dataframe, xcol = "x", ycol = "y"):
    points =  np.array( list(map(list, zip( dataframe[xcol], dataframe[ycol] ) )) )
    hull = ConvexHull(points, qhull_options = None)
    vertices_indexes = list(hull.vertices)
    return vertices_indexes


"""Check if point is inside convex polygon # https://stackoverflow.com/questions/1119627/how-to-test-if-a-point-is-inside-of-a-convex-polygon-in-2d-integer-coordinates
"""
RIGHT = "RIGHT"
LEFT = "LEFT"
def inside_convex_polygon(point, vertices):
    previous_side = None
    n_vertices = len(vertices)
    for n in xrange(n_vertices):
        a, b = vertices[n], vertices[(n+1)%n_vertices]
        affine_segment = v_sub(b, a)
        affine_point = v_sub(point, a)
        current_side = get_side(affine_segment, affine_point)
        if current_side is None:
            return False #outside or over an edge
        elif previous_side is None: #first segment
            previous_side = current_side
        elif previous_side != current_side:
            return False
    return True
def get_side(a, b):
    x = cosine_sign(a, b)
    if x < 0:
        return LEFT
    elif x > 0: 
        return RIGHT
    else:
        return None
def v_sub(a, b):
    return (a[0]-b[0], a[1]-b[1])
def cosine_sign(a, b):
    return a[0]*b[1]-a[1]*b[0]


# =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Freiberg Gneiss\MidSquare').as_posix() #flips the backslashes to plug in open()
# inputfilename = pathlib.PureWindowsPath(r'Extracted M CNL_SideB_top').as_posix()
# bs="\\"; wd=wd+bs                               # Use this instead in Windows
# 
# vcol = 'z'
# save = False
# 
# df = pd.read_csv(wd + inputfilename + extension_txt, usecols=[0,1,2], names=['x', 'y', 'z'])                     # read a .csv file in as a DataFrame. The input file needs to be in a table format, i.e. columns of 'x', 'y', 'z' with #rows = #points. As is, it only reads columns 0,1,2 and ignores the rest.
# # df = pd.read_csv(wd + inputfilename + extension_csv)                  #This new line is because reading the .csv file with a header and an index column will create some issue
# # ResVarMap = pd.read_csv(wd + inputfilename + extension_csv, header=0, index_col=0)                              # Read csv's of the variogram maps
# 
# points =  np.array( list(map(list, zip( df['x'], df['y'] ) )) )
# hull = ConvexHull(points, qhull_options = None)
# vertices_indexes = list(hull.vertices)
# 
# plt.plot(points[:,0], points[:,1], 'o')
# for simplex in hull.simplices:
#     plt.plot(points[simplex, 0], points[simplex, 1], 'r-')
# vindex = list(hull.vertices); vindex.append(hull.vertices[0]); vindex = np.array(vindex) # add the first point to the end in order to connect the last point with the first
# plt.plot(points[vindex,0], points[vindex,1], 'r--', lw=2)
# # plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
# # plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
# plt.show()
# =============================================================================



"""
Jarvis March
package com.stablesort.convexhull;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Implementation of "Jarvis march" algorithm for solving Convex Hull problem: https://en.wikipedia.org/wiki/Gift_wrapping_algorithm
 * 
 * video tutorial at: https://youtu.be/B2AJoQSZf4M
 * 
 * Starts off by finding the lowest point on the y-axis for its 1st vertex.
 * Then, instead of doing any kind of sorting, it just loops through all of the points again in a brute force way to find the point that makes 
 * the smallest counterclockwise angle in reference to the previous vertex. It simply repeats this iteration through all of the points until 
 * all of the vertices are determined and it gets back to the starting point.
 * 
 * @author Andre Violentyev
 */
public class ConvexHullJarvisMarch {
	public List<Point> marchByAngle(Collection<? extends Point> points) {
		List<Point> hull = new ArrayList<>();
		
		Point startingPoint = GraphUtils.getMinY(points); // starting point guaranteed to be on the hull		
		hull.add(startingPoint);
				
		Point prev = startingPoint;
		float prevAngle = -1;
		
		while (true) {
			float minAngle = Float.MAX_VALUE;
			double maxDist = 0;
			Point next = null;
		
			/*
			 * iterate over every point and pick the one that creates the largest angle
			 */
			for (Point p : points) {
				if (p == prev) continue;
				
				float angle = GraphUtils.angle(prev, p);
				double dist = GraphUtils.dist(prev, p);
				int compareAngles = Float.compare(angle, minAngle);
				
				if (compareAngles <= 0 && angle > prevAngle) {
					if (compareAngles < 0 || dist > maxDist) {
						/*
						 * found a better Point. It either has a smaller angle, or if it's collinear, then it's further way
						 */
						minAngle = angle;
						maxDist = dist;
						next = p;						
					}
				}
			}
			
			if (next == startingPoint) break; // came back to the starting point, so we are done
			
			hull.add(next);
			
			prevAngle =  GraphUtils.angle(prev, next);
			prev = next;			
		}
		
		return hull;
	}
	
	
	/**
	 * 
	 * @param points
	 * @return
	 */
	public List<Point> march(Collection<? extends Point> points) {
		List<Point> hull = new ArrayList<>();
		
		Point startingPoint = GraphUtils.getMinY(points); // bottom most, left most point is guaranteed to be on the hull		
		hull.add(startingPoint);
				
		Point prevVertex = startingPoint;
		
		while (true) {
			Point candidate = null;		
			/*
			 * iterate over every point and pick the one that creates the smallest counterclockwise angle
			 * in reference to the previous vertex
			 */
			for (Point point : points) {
				if (point == prevVertex) continue;
				
				if (candidate == null) {
					candidate = point;
					continue;
				}
				
				int ccw = GraphUtils.ccw(prevVertex, candidate, point);
				
				if (ccw == 0 && GraphUtils.dist(prevVertex, candidate) < GraphUtils.dist(prevVertex, point)) {
					candidate = point; // collinear tie is decided by the distance
				} else if (ccw < 0) {
					/*
					 * if made a clockwise turn, then found a better point that makes a smaller 
					 * counterclockwise angle in reference to the previous vertex
					 */
					candidate = point; 
				}				
			}
			
			if (candidate == startingPoint) break; // came back to the starting point, so we are done
			
			hull.add(candidate);
			prevVertex = candidate;			
		}
		
		return hull;
	}
	
	public static void main(String[] args) {
		List<Point> points = new ArrayList<>();
		
		points.add(new Point(2, 2));
		points.add(new Point(1, 3));
		points.add(new Point(1, 1));		
		points.add(new Point(0, 2));
		points.add(new Point(3, 1));		
		
		ConvexHullJarvisMarch hull = new ConvexHullJarvisMarch();
		
		
		System.out.println("Jarvis March:" + hull.march(points));
	}
}
"""


"""
Graham Scan
package com.stablesort.convexhull;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.List;

/**
 * The Convex Hull of a set of N points is the smallest perimeter fence that encloses all of the points.
 * Alt: The smallest area convex polygon enlosing all of the points.
 * 
 * Implementation of "Graham scan" algorithm:
 * 1. find the point p with the smallest y-coordinate
 * 2. sort the point by polar angle from p
 * 3. iterate over the points in sorted order, discard unless a counter clockwise turn is created
 * 3.1 make sure that the sort places collinear points in expected sequence 
 * 
 * video tutorial at: https://youtu.be/B2AJoQSZf4M
 * @author Andre Violentyev
 */
public class ConvexHullGrahamScan {

	/**
	 * The comparator just compares the cross product of two vectors to see which one is on the
	 * left side and which is on the right side. Actual angles don't need to be calculated.
	 * 
	 * @param points
	 * @param ref
	 */
	private void sortByAngle(List<? extends Point> points, Point ref) {
		Collections.sort(points, (b, c) -> {
			/*
			 * the ref point should always be pushed to the beginning
			 */
			if (b == ref) return -1;
			if (c == ref) return 1;
			
			int ccw = GraphUtils.ccw(ref, b, c);
			
			if (ccw == 0) {
				/*
				 * Handle collinear points. We can just use the x coordinate and not 
				 * bother with the y since the ratio of y/x is going to be the same
				 */
				if (Float.compare(b.x, c.x) == 0) {
					/*
					 * rare case of floats matching up in a vertical line, we want 
					 * the closer point to be first
					 */
					return b.y < c.y ? -1 : 1;				
				} else {
					return b.x < c.x ? -1 : 1;				
				}				
			} else {
				return ccw * -1;
			}
		});			
	}
	
	/**
	 * The main algorithm. 
	 * 
	 * @param points
	 * @return
	 */
	public List<Point> scan(List<? extends Point> points) {
		Deque<Point> stack = new ArrayDeque<Point>();
		
		/*
		 * bottom most, left most point is guaranteed to be on the hull
		 */
		Point minYPoint = GraphUtils.getMinY(points);		
		sortByAngle(points, minYPoint); // sort by angle with respect to minYPoint
		
		stack.push(points.get(0)); // 1st point is guaranteed to be on the hull
		stack.push(points.get(1)); // don't know about this one yet
		
		for (int i = 2, size = points.size(); i < size; i++) {
			Point next = points.get(i);
			Point p = stack.pop();			
			
			while (stack.peek() != null && GraphUtils.ccw(stack.peek(), p, next) <= 0) { 
				p = stack.pop(); // delete points that create clockwise turn
			}
						
			stack.push(p);
			stack.push(points.get(i));
		}
		
		/*
		 * the very last point pushed in could have been collinear, so we check for that
		 */
		Point p = stack.pop();
		if (GraphUtils.ccw(stack.peek(), p, minYPoint) > 0) {
			stack.push(p); // put it back, everything is fine
		}
		
		return new ArrayList<>(stack);
	}
	
	public static void main(String[] args) {
		List<Point> points = new ArrayList<>();
		
		points.add(new Point(2, 2));
		points.add(new Point(-2, 3));
		points.add(new Point(1, 1));		
		
		ConvexHullGrahamScan hull = new ConvexHullGrahamScan();
		System.out.println("Graham Scan:" + hull.scan(points));
	}
}
"""