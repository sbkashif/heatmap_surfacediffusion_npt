# Salman Bin Kashif
# Sarupria Group
# Clemson University

#Date created: July 1, 2019
#Last updated: July 8, 2019 

import MDAnalysis as mda
import numpy as np
import math
import statistics
import sys
import argparse
import os

def create_parser():
	"""

	This is the function to pass the input values
	
	:type binx: int
	:param binx: No. of bins in X-direction
	
	"""

	parser = argparse.ArgumentParser(prog = 'HeatMapUsingMDAnalysis', usage = '%(prog)s [-h for help]', \
                                      description = 'Generate the heat map for a trajectory')
	parser.add_argument('-f', "--f", help = 'Input xtc file (Required).')
	parser.add_argument('-s',"--s",help='Input pdb file (required).')
	parser.add_argument('-binx', "--binx", help = "Enter number of x bins", default=int(0.5))
	parser.add_argument('-biny', "--biny", help = "Enter number of y bins", default=int(0.5))
	parser.add_argument('-bf',"--bf",help="Begin Frame (Frame ID)",default=int(0))
	parser.add_argument('-ef',"--ef",help="End Frame (Frame ID)")
	return parser


def main(args):
	"""
	This is the main function which basically processes the input parameters
	
	"""
	XTCFile=args.f
	GROFile=args.s
	binx=float(args.binx)
	biny=float(args.biny)
	bf=int(args.bf)
	ef=int(args.ef)
	getTrajectoryDistribution(binx,biny,bf,ef)


def obtainTrajectoryData(x_bin_size,y_bin_size):
        global coordinates_pa
        global coordinates_mbl
        global mbl
        global size
        global box
        global x_n_bins
        global y_n_bins
        global Lx
        global Ly
        global Lz
	
        #Read the full trajectory from the xtc file
        u=mda.Universe("./last10ns.gro","./last10ns.xtc")
        print(u)

        #Assign the polyamide atom coordinates 'pa' variable
        pa=u.select_atoms("resname MPD1 MPD2 TMC1 TMC2 TMC3")
        print(pa)
        
        #Assign the methylene blue atom coordinates to 'mbl' variable
        mbl=u.select_atoms("resname MBL")
        print("Following are the coordinates for MBL:\n")
        print(mbl)
        coordinates_pa=[]
        coordinates_mbl=[]
        box=[]
        cog_pa=[]
        cog_mbl=[]
        cog=[]
        coordinates_outside_atom=[]
        coordinates_inside_atom=[]
        outside_atom=u.select_atoms("name N3")
        inside_atom=u.select_atoms("name N1")
	
        for ts in u.trajectory[0:3]:
                coordinates_mbl.append(mbl.positions)
                coordinates_pa.append(pa.positions)
                box.append(ts.dimensions[:3])
                cog.append(mbl.center_of_geometry())
                coordinates_outside_atom.append(outside_atom.positions)
                coordinates_inside_atom.append(inside_atom.positions)
        
        coordinates_pa=np.array(coordinates_pa)
        coordinates_mbl=np.array(coordinates_mbl)
        coordinates_outside_atom=np.array(coordinates_outside_atom)
        coordinates_inside_atom=np.array(coordinates_inside_atom)
        box=np.array(box)
        cog=np.array(cog)
    	#Since MDAnalysis give output in Angstorm, the units are converted to nm for consistency
	
        #coordinates_pa=np.divide(coordinates_pa,10.0)
        #coordindates_mbl=np.divide(coordinates_mbl,10.0)
        #cog=np.divide(cog,10.0)
        #box=np.divide(box,10.0)
        
       
        print("COG from MD analysis:\n")
        print(cog)         
        for i in range(0,len(coordinates_mbl)):
            print("Frame:",i)
            print("Coordinates:\n")
            print("Outside atom:",coordinates_outside_atom[i])
            print("Inside atom:",coordinates_inside_atom[i])
            print("COG from my calculation",readFrame(i,coordinates_mbl))

        #Extracting the box size thorughout the trajectory
        Lx=box[:,0]
        Ly=box[:,1]
        Lz=box[:,2]
        Lx_max=max(Lx)
        Ly_max=max(Ly)

        size=len(coordinates_pa)
        print("No. of frames:%d"%(size))
        #The use of this part is in heatmap and it will be revoked later
        #x_n_bins=math.ceil(Lx_max/x_bin_size)
        #y_n_bins=math.ceil(Ly_max/y_bin_size)
    
        return coordinates_pa, coordinates_mbl, box

def readFrame(frame, group_coordinates):
        global x
        global y
        global z
        global p_x
        global p_y
        global p_z
        global Lx
        global Ly
        global Lz
        global Lx_max   
        global Ly_max
   	
	
	    #Obtaining the x,y and z coordinates for the given frame
        #The coordiates are extracted to treat raw or with PBC
        x=group_coordinates[frame][:,0]
        y=group_coordinates[frame][:,1]
        z=group_coordinates[frame][:,2]
        p_x=group_coordinates[frame][:,0]
        p_y=group_coordinates[frame][:,1]
        p_z=group_coordinates[frame][:,2]
 
        print("readFrame X coordinates read:\n",x)
        print("readFrame Y coordinates read:\n",y)
        Lx=box[:,0]
        Ly=box[:,1]
        Lz=box[:,2]
        
        #print("Box dimension of frame: %s\t is %s,%s,%s"%(frame,Lx[frame],Ly[frame],Lz[frame]))
        

        #Applying the pbc to p_* variables 
        for i in range(0,len(x)):
                #Uncomment next line to print the x-coordinate value before applying the PBC
                #print(x[i])
                data=p_x[i]/Lx[frame]
                fd=math.floor(data)
                p_x[i]=p_x[i]-(fd*Lx[frame])
        for i in range(0,len(y)):
                #Uncomment next line to print the y-coordinate value before the applyinf the PBC
                #print(y[i])
                data=p_y[i]/Ly[frame]
                fd=math.floor(data)
                p_y[i]=p_y[i]-(fd*Ly[frame])
        for i in range(0,len(z)):
                data=p_z[i]/Lz[frame]
                fd=math.floor(data)
                p_z[i]=p_z[i]-(fd*Lz[frame])
       
        print(x) 
        print(y) 
        print("Sum of y:", sum(y))
        print("Length of y:",len(y))
        #COG calculated for the raw trajectory
        cog_x=sum(x)/len(x)
        cog_y=sum(y)/len(y)
        #print("COG of frame:%s\t is %s\t%s:" %(frame,cog_x,cog_y))
        return cog_x, cog_y
    
    

def getTrajectoryDistribution(x_bin_size,y_bin_size,bf,ef):
        
        global coglist_pa
        global coglist_mbl	
        obtainTrajectoryData(x_bin_size,y_bin_size)	
        coglist_pa=np.full((ef,2),-100.0)
        coglist_mbl=np.full((ef,2),-100.0)
       
         
        for i in range(bf,ef):
            print ("Looping in Frame:%d\n"%(i))
        #Calculating the MBL cog for frame i
            cog_x, cog_y = readFrame(i,coordinates_mbl)
            coglist_mbl[i][0]=cog_x
            coglist_mbl[i][1]=cog_y
        #Calculating the PA cog for frame i
            readFrame(i,coordinates_pa)
            coglist_pa[i][0]=cog_x
            coglist_pa[i][1]=cog_y
        
        coglist_mbl_w=open('./cog_list_allComp.xvg',"w")
        coglist_mbl_w.write('Frame\tMBL X\tMBL Y\tPA X\tPA Y\n')
        for i in range(bf,ef):
            coglist_mbl_w.write('%d\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n'%(i,coglist_mbl[i][0],coglist_mbl[i][1],coglist_pa[i][0],coglist_pa[i][1]))
            
            
#if __name__=='__main__':
#	args = create_parser().parse_args()
#	main(args)
    
obtainTrajectoryData(0.1,0.1)
#readFrame(1)
#getTrajectoryDistribution(0.1,0.1,1,100)
