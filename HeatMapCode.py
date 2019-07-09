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


def obtainTrajectoryData(x_bin_size,y_bin_size,bf,ef):
        global coordinates_pa
        global coordinates_mbl
        global size
        global box
        global x_n_bins
        global y_n_bins
        global Lx
        global Ly
        global Lz
        global logfile

        logfile=open("output.log","w")	
        #Read the full trajectory from the xtc file
        #u=mda.Universe("../test_files/last10ns.gro","../test_files/last10ns.xtc")
        u=mda.Universe("../test_files/prod-pacxb-npt-whole.gro","../test_files/prod-pacxb-npt-whole.xtc")
        print(u)

        #Assign the polyamide atom coordinates 'pa' variable
        pa=u.select_atoms("resname MPD1 MPD2 TMC1 TMC2 TMC3")
        logfile.write("\nRead input trajectory Polyamdie coordinates:\n%s"%(pa))
        
        #Assign the methylene blue atom coordinates to 'mbl' variable
        mbl=u.select_atoms("resname MBL")
        logfile.write("\nRead input trajectory Methylene Blue coordinates:\n%s"%(mbl))
        coordinates_mbl=[]
        coordinates_pa=[]
        box=[]

        bf=int(bf-1)
        ef=int(ef)
        print(u.trajectory)
        
	
        for ts in u.trajectory[bf:ef]:
                print(ts)
                logfile.write("\nReading timestep:%s\n"%(ts))
                coordinates_mbl.append(mbl.positions)
                coordinates_pa.append(pa.positions)
                box.append(ts.dimensions[:3])
        
        coordinates_pa=np.array(coordinates_pa)
        coordinates_mbl=np.array(coordinates_mbl)

        box=np.array(box)
    	
        #Since MDAnalysis give output in Angstorm, the units are converted to nm for consistency
        coordinates_pa=np.divide(coordinates_pa,10.0)
        
        coordinates_mbl=np.divide(coordinates_mbl,10.0)
        
        box=np.divide(box,10.0)
       
        #print("COG from MD analysis:\n")
        #print(cog)         

        #Extracting the box size thorughout the trajectory
        Lx=box[:,0]
        Ly=box[:,1]
        Lz=box[:,2]
        Lx_max=max(Lx)
        Ly_max=max(Ly)

        size=len(coordinates_mbl)
        print("\nNumber of frames read:%s\n"%(size))
        logfile.write("\nNumber of frames read:%s\n"%(size))
        #The use of this part is in heatmap and it will be revoked later
        x_n_bins=math.ceil(Lx_max/x_bin_size)
        y_n_bins=math.ceil(Ly_max/y_bin_size)
    
        return coordinates_pa, coordinates_mbl, box

def calculateFrameCOG(frame, input_group_coordinates):
	
	    #Obtaining the x,y and z coordinates for the given frame
        #The coordiates are extracted to treat raw or with PBC
        group_coordinates=input_group_coordinates.copy()
        x=input_group_coordinates[frame][:,0]
        y=input_group_coordinates[frame][:,1]
        z=input_group_coordinates[frame][:,2]
        p_x=group_coordinates[frame][:,0]
        p_y=group_coordinates[frame][:,1]
        p_z=group_coordinates[frame][:,2]
 
        #print("readFrame X coordinates read:\n",x)
        #print("readFrame Y coordinates read:\n",y)
        #Lx=box[:,0]
        #Ly=box[:,1]
        #Lz=box[:,2]
        
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
       
        if(frame==0 or (frame+1)%10==0 or (frame+1)==size):
            logfile.write("\ncalculateFrameCOG function X coordinates read:\n%s"%(x))
            logfile.write("\ncalculateFrameCOG function Y coordinates read:\n%s"%(y))

        #COG calculated for the raw trajectory
        cog_x=sum(x)/len(x)
        cog_y=sum(y)/len(y)
        
        #COG converted to PBC
        data_x=cog_x/Lx[frame]
        fd=math.floor(data_x)
        cog_x=cog_x-(fd*Lx[frame])

        data_y=cog_y/Ly[frame]
        fd=math.floor(data_y)
        cog_y=cog_y-(fd*Ly[frame])

        return cog_x, cog_y
    
    

def getTrajectoryDistribution(x_bin_size,y_bin_size,bf,ef):
        
        #global coglist_pa
        #global coglist_mbl	
        obtainTrajectoryData(x_bin_size,y_bin_size,bf,ef)	
        coglist_pa=np.full((ef,2),-100.0)
        coglist_mbl=np.full((ef,2),-100.0)
        freq_cog_pa=np.full((x_n_bins,y_n_bins),0)
        freq_cog_mbl=np.full((x_n_bins,y_n_bins),0)
        bf=bf-1
        
        for i in range(bf,ef):
            i=i-bf
            print ("Looping in frame corresponding to extracted trajectory timestep:%d\n"%(i))
            #Calculating the MBL cog for frame i
            logfile.write("\nCalculating COG for MBL for frame corresponding to extracted trajectory timestep:%s\n"%(i))
            cog_x, cog_y = calculateFrameCOG(i,coordinates_mbl)
            coglist_mbl[i+bf][0]=cog_x
            coglist_mbl[i+bf][1]=cog_y
        #Determining and appending the bin frequency
            bin_x_ID=math.floor(cog_x/x_bin_size)
            bin_y_ID=math.floor(cog_y/y_bin_size)
            freq_cog_mbl[bin_x_ID][bin_y_ID]+=1
	   
        #Calculating the PA cog for frame i
            logfile.write("\nCalculating COG for PA for frame corresponding to extracted trajectory timestep:%s\n"%(i))
            cog_x,cog_y = calculateFrameCOG(i,coordinates_pa)
            coglist_pa[i+bf][0]=cog_x
            coglist_pa[i+bf][1]=cog_y
            bin_x_ID=math.floor(cog_x/x_bin_size)
            bin_y_ID=math.floor(cog_y/y_bin_size)
            freq_cog_pa[bin_x_ID][bin_y_ID]+=1
	   
 

        coglist_w=open('./cog_list_allComp.xvg',"w")
        coglist_w.write('#Frame\tMBL X\tMBL Y\tPA X\tPA Y\n')
        freq_mbl_w=open('./freq_list_MBL.xvg',"w")
        freq_pa_w=open('./freq_list_PA.xvg',"w")
        
        for i in range(bf,ef):
            coglist_w.write('%d\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n'%((int)(i+1),coglist_mbl[i][0],coglist_mbl[i][1],coglist_pa[i][0],coglist_pa[i][1]))
        
        for j in range(0,x_n_bins):
            for k in range(0,y_n_bins):
            #if(j*x_bin_size <=Lx[i] and k*y_bin_size <= Ly[i]):
                n_freq_mbl=freq_cog_mbl[j][k]/size
                n_freq_pa=freq_cog_pa[j][k]/size
                freq_mbl_w.write('%8.3f\t%8.3f\t%8.3f\t%8.3f\n'%(j*x_bin_size+x_bin_size/2,k*y_bin_size+y_bin_size/2,freq_cog_mbl[j][k],n_freq_mbl))
                freq_pa_w.write('%8.3f\t%8.3f\t%8.3f\t%8.3f\n'%(j*x_bin_size+x_bin_size/2,k*y_bin_size+y_bin_size/2,freq_cog_pa[j][k],n_freq_pa))
            freq_pa_w.write("\n") 
            freq_mbl_w.write("\n")
                    
	   
#if __name__=='__main__':
#	args = create_parser().parse_args()
#	main(args)
    
#obtainTrajectoryData(0.1,0.1)
#readFrame(1)
getTrajectoryDistribution(0.3,0.3,19999,20001)
