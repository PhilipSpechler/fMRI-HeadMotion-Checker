#Addendum code for HCP scripts
#Run me -AFTER- generic fMRI Volume!
#################################
##  Written by Philip Spechler ##
#################################
#
##############################################################################
### Usage: Give two positional arguments, project folder name, subjectID
### Ex: python <this script> <project_folder_name> <subject_ID>
##############################################################################
# This code will calculate FD from the motion regressors output from mcflirt
# This code will also calculate head motion summary stats into a csv for QC
# This code will also generate pretty graphs for human QC
###############################################################################
# Look for my output in your subjects "MNINonLinear/Results/MotionSummary"
# I will make the Motion Summary directory ^ and populate it with lots of goodies
##############################################################################
import os,sys,shutil
import numpy as np
import matplotlib
matplotlib.use('agg')
import pylab as plt
#####################################################################
#	Edit the modality list for your fMRI task names
#	List element expected to be a comma separated string
######################################################################
modality=['stroop','mid','rsfMRI','sst']
######################################################################
#####   Set the voxelsize (vxsz), FD threshold (FDthrs)     ##########
### Also set number of Pre and Post volumes to be censored/scrubbed  #
######################################################################
#Used for plotting nice visuals, threshold boundaries, and summary stats
vxsz=2.25
FDthrsh=0.89
nPre=1
nPost=1
#############################################################
####### GENERALLY NO NEED TO EDIT BELOW THIS LINE ###########
#############################################################
homedir=os.path.expanduser('~/data/')
proj=sys.argv[1]
sub=sys.argv[2]
#############################################################
############     PATH TO MOTION DIRS    #####################
#############################################################
datadir=homedir+proj+'/'+sub+'/MNINonLinear/Results/'
summarydir=homedir+proj+'/'+sub+'/MNINonLinear/Results/MotionSummary/'
#
#
try:
	shutil.rmtree(summarydir)
except :
	print("\n\t\t****Motion Summary Directory Exists.... Now Deleting! ****")
#
try:
	os.mkdir(summarydir)
except :
	print("Made Fresh New Motion Summary Directory")

##########################################################
##### This is the function to calculate FD
###########################################################
def calculate_FD_P(in_file,fMRI):
    """
    Method to calculate Framewise Displacement (FD) calculations
    (Power et al., 2012)
    
    Parameters
    ----------
    in_file : string
        movement parameters vector file path
    
    Returns
    -------
    out_file : string
        Frame-wise displacement mat 
        file path
    
    """

    lines = open(in_file, 'r').readlines()
    rows = [[float(x) for x in line.split()] for line in lines]
    cols = np.array([list(col) for col in zip(*rows)])
    
    #translations = np.transpose(np.abs(np.diff(cols[3:6, :])))
    #rotations = np.transpose(np.abs(np.diff(cols[0:3, :])))
    
    #This is flipped in HCP mcflirt-- translations first, then rotations
    rotations = np.transpose(np.abs(np.diff(cols[3:6, :])))
    translations = np.transpose(np.abs(np.diff(cols[0:3, :])))
    
    FD_power = np.sum(translations, axis = 1) + (50*3.141/180)*np.sum(rotations, axis =1)
    
    #FD is zero for the first time point
    FD_power = np.insert(FD_power, 0, 0)
    
    #in_file=in_file.split(".par")[0]
    out_file = fMRI+'_FD.txt'
    
    np.savetxt(out_file, FD_power)
    
    return out_file
#
#############################################
## This function will make the graphs--   ###
#############################################
def make_pretty_motion_graph(fMRI,sub,mopar,FD):
	#Make Composite Graph
	plt.suptitle('{0} Task: All Head Motion\nSubject: {1}'.format(fMRI,sub))
	plt.subplot(2, 1, 1)
	plt.axhspan(-vxsz,vxsz,hold=True,facecolor='0.75')
	plt.plot(mopar,linewidth=2.5)
	plt.ylabel('Head  Displacement (mm)')
	plt.tight_layout(pad=3)
	#
	plt.subplot(2, 1, 2)
	plt.axhspan(FDthrsh,0,hold=True,facecolor='0.75')
	plt.plot(FD,'k',linewidth=2.5)
	plt.xlabel('Volume Index')
	plt.ylabel('FD (mm)')
	plt.ylim(ymin=0)
	plt.tight_layout(pad=3)
	plt.savefig(summarydir+sub+'_'+fMRI+'_MotionGraph.pdf')
	plt.close()
	return
#####################################################
## This function will calculate summary stats--   ###
##     -AND- Call up Write CSV Function           ###
#####################################################
def head_motion_summary(fMRI,sub,mopar,FD,nPre,nPost):
	########################################################
	#Set boundaries so not scrubbing non-existent volumes
	boundaries=[-1]
	boundaries.extend((len(FD),len(FD)+1,len(FD)+2))
	#########################################################
	#First set an empty array of equal size to FD array
	FDbinary=np.zeros(np.shape(FD))
	FDbinary[0:8]=1	
	if any(FD>FDthrsh):
    	#Set a list of indexes where FD passes threshold
		scrubvols=np.argwhere(FD>FDthrsh)
    	############################################################################
		##### Find vol index where FD>thrshold with nPre and nPost Scrubbed vols ###
		#####       As suggested by Power 2012 NeuroImage     ######################
		############################################################################
		if nPre==1:
			for k in np.argwhere(FD>FDthrsh)-1:
				if k not in scrubvols and k not in boundaries:
					scrubvols=np.append(scrubvols,k)
			else:
				pass
		if nPost==1 or nPost==2:
			for l in np.argwhere(FD>FDthrsh)+1:
				if l not in scrubvols and l not in boundaries:
					scrubvols=np.append(scrubvols,l)
				else:
					pass
			if nPost==2:
				for m in np.argwhere(FD>FDthrsh)+2:
					if m not in scrubvols and m not in boundaries:
						scrubvols=np.append(scrubvols,m)
					else:
						pass
    	#Save both the index value and binary version
		try:
			if sub+'_'+i+'_MotionScrub_Index.txt' not in os.listdir(summarydir):
				scrub_vol_index=np.savetxt(summarydir+sub+'_'+i+'_MotionScrub_Index.txt',scrubvols,fmt='%.f')
				#Set FDbinary array indices to 1 for each scrubvol index 
				FDbinary[scrubvols]=1
				scrub_vol_binary=np.savetxt(summarydir+sub+'_'+i+'_MotionScrub_Binary.txt',FDbinary,fmt='%.f')
			else:
				print("\t****Censor File Already Exists for: "+i)
		except :
			pass
	else:
		scrub_vol_index=[]
		scrub_vol_binary=[]
	##################################################
	#### Write FDbinary to "censor.txt" for 1st level#
	##################################################
	np.savetxt(datadir+i+'/'+'censor.txt',FDbinary,fmt='%i')
	################################################
	###  SEND TO WRITE SUMMARY STATS FUNCTION   ####
	################################################
	write_motion_stats(fMRI,sub,mopar,FD,FDbinary)
	return

#####################################################
## This function writes summary stats to CSV--    ###
#####################################################
def write_motion_stats(fmri,sub,mopar,FD,FDbinary):
	try:
		with open(summarydir+sub+'_MotionStats.csv','a') as f:
			f.write(fmri+','+
			str(mopar.min())+','+
			str(mopar.max())+','+
			str(FD.max())+','+
			str(FD.mean())+','+
			str(int(sum(FDbinary)))+','+
			str(int(len(FDbinary)-sum(FDbinary)))+','+
			str(sum(FDbinary)/len(FDbinary))+'\n')
			f.close()
		print("\n+++ Summary stats written for : "+i+" +++\n")
		del FDbinary
	except:
		print("\t****NO summary stats for: "+i)
	return f
##########################################
###  CD to proj folder and given sub  ##
##########################################
try:
	os.chdir(datadir)
except :
	print("\t****Can't change directory to: "+datadir)
##########################################
###  This loop will just calculate FD  ##
##########################################
for i in os.listdir('.'):
	for j in modality:
		if any(char.isdigit() for char in i):#To cd to run folders of task, rather than 2nd Level folder
			try:
				os.chdir(i)
				if i+'_FD.txt' not in os.listdir('.'):
					calculate_FD_P('Movement_Regressors.txt',i)
					print("\tFramewise Displacement Calculation Complete for: "+i)
				else:
					print("\t****Framewise Displacement already calcualted for: "+i)
				os.chdir('../')
			except:
				print("\t****NO FD for: "+j)
				os.chdir('../')
##########################################
#####  Setup Summary Stat CSV File  ######
##########################################
try:
	if sub+'_MotionStats.csv' not in os.listdir(summarydir):
		with open(summarydir+sub+'_MotionStats.csv','a') as f:
			f.write("Task,minDisp,maxDisp,maxFD,meanFD,nVol2Scrub,nNetVols,%Scrubbed"+'\n')
			f.close()
	else:
		print("\t****Summary Stat CSV Exists")
except :
	print("\t*****ERROR!!!! Unable to setup summary stat csv")
#########################################################
###  This loop will---
### 	1.) Make pretty motion graphs
###		2.) Calculate Head Motion Summary Stats
###		3.) Write summary stats to csv
#########################################################
os.chdir(datadir)
for i in os.listdir('.'):
	for j in modality:
		if j+'1' in i or j+'2' in i:
			try:
				os.chdir(i)
				print("Finding Motion Estimates Here: "+os.getcwd())
				######################################
				mopar=np.genfromtxt('Movement_Regressors.txt',usecols=(0,1,2,3,4,5))
				FD=np.genfromtxt(i+'_FD.txt')
				######################################
				try:
					print("Making Pretty Graph-- ")
					make_pretty_motion_graph(i,sub,mopar,FD)
				except :
					print("Unable to make pretty graphs for: "+i)
				######################################
				try:
					print("Calculate and Write Summary Stats for: {0}-- \n".format(i))
					head_motion_summary(i,sub,mopar,FD,nPre,nPost)
				except :
					print("Unable to calculate summary stats for: "+i)
				#######################################
				del mopar
				del FD
				os.chdir('../')				
			except:
				print("\t****WARNING****Unable to make graphs, summary stats, log, for: "+i)
				os.chdir('../')
##############################################################################				    
print("ALL DONE")
os.chdir(homedir)
