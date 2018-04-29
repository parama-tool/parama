#!/usr/bin/python
import sys
def create_detailed_results(pdb_id,chain_id):
    infile=open('results/'+pdb_id+"_"+chain_id+"_disall_points.csv","r")
    opfile=open('results/'+pdb_id+"_"+chain_id+"_detailed_results.txt","w")
    opfile.write('********************************Detailed Results for '+pdb_id+'****************************************\n\n')
    opfile.write('This file contains the bond lengths (in angstroms), angles (in degrees) and omega planarity values (in degrees) observed in 3 or lesser small molecule peptide structures, for which the correspoding phi-psi has been classified as allowed. The ideal Pauling-Corey bond parameter values are given in brackets next to each observed value.\n\n')
    opfile.write('Following is the naming convention used for the atoms in a two-linked peptide unit.\n\n')
    opfile.write('     O1      CB       H3\n')
    opfile.write('     ||       |        |\n')   
    opfile.write('CA1--C1--N2--CA2--C2--N3--CA3\n')
    opfile.write('          |   |   || \n')
    opfile.write('         H2  HA   O2 \n\n')
    opfile.write('The residue with the disallowed conformation in the protein corresponds to the CA2 atom\n\n')
    for line in infile:
        lineparts=line.split(",")
        prob=float(lineparts[-1])
        if prob!=0 and prob<0.25:
            opfile.write("\nExpected geometries for : "+lineparts[2]+lineparts[3]+"\t phi = "+lineparts[0]+"\tpsi = "+lineparts[1]+"\tprobability = "+str(prob)+"\n\n")
            disal_file=open("disallowed_phipsi_geometry/"+lineparts[0]+"_"+lineparts[1]+".dat","r")
            opfile.write(disal_file.read())
            opfile.write("\n----------------------------------------------------------------------------------------------------------\n")
            disal_file.close()
    opfile.close()

