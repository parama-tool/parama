#!/usr/bin/python

import os
import map_plotter
import pdb_file_processing
import get_shrt_conts
import get_disallowed_phipsi_geometry
import sys


no_of_args=len(sys.argv)
if no_of_args>3 or no_of_args<3:
    sys.exit("Argument error : Incorrect number of arguments\n\nUsage : \tpython parama_results.py <pdb_file_name> <chain_id>\n\n")

pdb_file_name=sys.argv[1].strip()
chain_id=sys.argv[2].strip().upper()
if not os.path.isfile(pdb_file_name):
    sys.exit("File error : The entered pdb file does not exist\n\n")

pdb_id=pdb_file_name.split(".")[0]

logfile=open(pdb_id+".log","w")
logfile.write("***************PARAMA************************\n\n")
pdb_file_processing.checkMultipleModel(pdb_id) #in case of structures with multiple models, retain only one model and rewrite the saved pdb file
result=pdb_file_processing.calculateDihedrals(pdb_id,chain_id) #calculate phi-psi for all the residues in given chain
if result==0: #no phi-psi has been calculated, and hence invalid chain id (or) no valid two-linked peptide units
    print '\n\n Chain not found in PDB file \n\n'
    logfile.write('\n\n Chain not found in PDB file \n\nProgram terminated\n')
    logfile.close()
    sys.exit()
#phi-psi has been calculated, make the plots
map_plotter.plot_prob(pdb_id,chain_id) #calls the function to create the histogram
print "\n\nHistogram of phi-psi probabilities plotted in "+pdb_id+"_"+chain_id+"_probability_distribution.png in results folder\n"
logfile.write("Histogram of phi-psi probabilities plotted in "+pdb_id+"_"+chain_id+"_probability_distribution.png in results folder\n")
print "Probabilities associated with each phi-psi value is found in "+pdb_id+"_"+chain_id+"_phipsi_probabilities.csv in results folder\n"
logfile.write("\nProbabilities associated with each phi-psi value is found in "+pdb_id+"_"+chain_id+"_phipsi_probabilities.csv in results folder\n")
map_plotter.plot_map(pdb_id,chain_id) #calls the function to plot the ramachandran map
print "Ramachandran plot on ensemble map created; stored as "+pdb_id+"_"+chain_id+"_rplot.png in results folder\n "
logfile.write("\nRamachandran plot on ensemble map created; stored as "+pdb_id+"_"+chain_id+"_rplot.png in results folder\n")
if os.stat('results/'+pdb_id+'_'+chain_id+'_disall_points.csv').st_size == 0: #this file contains all the disallowed phi-psi; if that file is empty, no disallowed conformations found
    print 'This chain has no disallowed points according to classical Ramachandran Map \n'
    print 'Program completed\n\n'
    logfile.write('\nThis chain has no disallowed points according to classical Ramachandran Map \n')
    logfile.write('\nProgram completed\n\n')
    logfile.close()
else: #disallowed phi-psi are there, do short contact analysis
    flag1=0
    get_shrt_conts.find_expected_short_contacts(pdb_id,chain_id) #find the expected short contacts
    get_shrt_conts.find_actual_distances(pdb_id,chain_id) #find the actual distances within structure for each of the short contacts
    print "Short contact analysis performed for the disallowed phi-psi values. Results can be found at "+pdb_id+"_"+chain_id+"_short_contact_analysis.txt in results folder\n"
    logfile.write("\nShort contact analysis performed for the disallowed phi-psi values. Results can be found at "+pdb_id+"_"+chain_id+"_short_contact_analysis.txt in results folder\n")
    infile=open('results/'+pdb_id+'_'+chain_id+'_disall_points.csv','r')
    for line in infile: #for each disallowed phi-psi
        lineparts=line.split(",")
        if float(lineparts[4])<=0.25: #if probability less than 0.25
            flag1=1 #note that phi psi with probability lower than 0.25 is present
    infile.close() #end of loop of disallowed phi-psi
    if flag1==1: #there are probability values lower than 0.25, generate detailed results pertaining to expected geometry
        get_disallowed_phipsi_geometry.create_detailed_results(pdb_id,chain_id)
        print "This structure has phi-psi values with probability lower than 0.25. Details on expected geometry for such phi-psi values, obtained from small molecule peptide structures, can be seen in "+pdb_id+"_"+chain_id+"_detailed_results.txt in results folder\n"
        logfile.write("\nThis structure has phi-psi values with probability lower than 0.25. Details on expected geometry for such phi-psi values, obtained from small molecule peptide structures, can be seen in "+pdb_id+"_"+chain_id+"_detailed_results.txt in results folder\n")
    else:
        print "This structure has no phi-psi values with probability lower than 0.25.\n"
        logfile.close("\nThis structure has no phi-psi values with probability lower than 0.25.\n")
    print 'Program completed\n\n'
    logfile.write('Program completed\n\n')
    logfile.close()
