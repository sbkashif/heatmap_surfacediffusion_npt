# Salman Bin Kashif
# Sarupria Group
# Clemson University

#Date created: July 1, 2019
#Last updated: 

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
        global coordinates_mbl
        global size
        global box
        global x_n_bins
        global y_n_bins
	
        #Read the full trajectory from the xtc file
        u=mda.Universe("../test_files/prod-pacxb-npt-whole.gro","../test_files/prod-pacxb-npt-whole.xtc")
        print(u)

        #Assign the polyamide atom coordinates 'polyamide' variable
        pa=u.select_atoms("resname MPD1 MPD2 TMC1 TMC2 TMC3")
        print(pa)
        
        #Assign the methylene blue atom coordinates to 'mbl' variable
        mbl=u.select_atoms("resname MBL")
        print("Following are the coordinated for MBL:\n")
        print(mbl)
        coordinates_pa=[]
        coordinates_mbl=[]
        box=[]
        cog_pa=[]
        cog_mbl=[]
	
        for ts in u.trajectory:
                coordinates_mbl.append(mbl.positions)
                coordinates_pa.append(pa.positions)
                box.append(ts.dimensions[:3])
	
        coordinates_pa=np.array(coordinates_pa)
        coordinates_mbl=np.array(coordinates_mbl)
        box=np.array(box)
    
    	#Since MDAnalysis give output in Angstorm, the units are converted to nm for consistency
	
        coordinates_pa=np.divide(coordinates_pa,10.0)
        coordindates_mbl=np.divide(coordinates_mbl,10.0)
        box=np.divide(box,10.0)
        Lx=box[:,0]
        Ly=box[:,1]
        Lz=box[:,2]
        Lx_max=max(Lx)
        Ly_max=max(Ly)

        size=len(coordinates_mbl)
        print("No. of frames%d"%(size))
        x_n_bins=math.ceil(Lx_max/x_bin_size)
        y_n_bins=math.ceil(Ly_max/y_bin_size)

def readFrame(frame):
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
        global cog_x
        global cog_y
        global cog_z
    	
	
	#Obtaining the x,y and z coordinated from the frame
        x=coordinates_mbl[frame][:,0]
        y=coordinates_mbl[frame][:,1]
        z=coordinates_mbl[frame][:,2]
        p_x=coordinates_mbl[frame][:,0]
        p_y=coordinates_mbl[frame][:,1]
        p_z=coordinates_mbl[frame][:,2]
 
        Lx=box[:,0]
        Ly=box[:,1]
        Lz=box[:,2]
        
        print(Lx[frame],Ly[frame],Lz[frame])
        for i in range(0,len(x)):
                #Uncomment next line to print the x-coordinate value before applying the PBC
                #print(x[i])
                data=x[i]/Lx[frame]
                fd=math.floor(data)
                p_x[i]=x[i]-(fd*Lx[frame])
        for i in range(0,len(y)):
                #Uncomment next line to print the y-coordinate value before the applyinf the PBC
                #print(y[i])
                data=y[i]/Ly[frame]
                fd=math.floor(data)
                p_y[i]=y[i]-(fd*Ly[frame])
        for i in range(0,len(z)):
                data=z[i]/Lz[frame]
                fd=math.floor(data)
                p_z[i]=z[i]-(fd*Lz[frame])
        
        print(y)
        cog_x=sum(x)/len(x)
        cog_y=sum(y)/len(y)
        print(cog_x)
        print(cog_y)
    
    

def getTrajectoryDistribution(x_bin_size,y_bin_size,bf,ef):
        
        global coglist	
        obtainTrajectoryData(x_bin_size,y_bin_size)	
        coglist=np.full((ef,2),-100.0)
        
        for i in range(bf,ef):
            readFrame(i)
            print ("Frame%d\n"%(i))
            coglist[i][0]=cog_x
            coglist[i][1]=cog_y
        
        coglist_w=open('./cog_list.xvg',"w")
        for i in range(bf,ef):
            coglist_w.write('%d\t%8.3f\t%8.3f\n'%(i,coglist[i][0],coglist[i][1]))
            
            
#if __name__=='__main__':
#	args = create_parser().parse_args()
#	main(args)
    
#obtainTrajectoryData(0.1,0.1)
#readFrame(1)
getTrajectoryDistribution(0.1,0.1,1,100)
