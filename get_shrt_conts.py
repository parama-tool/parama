#!/usr/bin/python


def find_expected_short_contacts(pdb_id,chain_id):
	import cPickle
	infile=open('results/'+pdb_id+'_'+chain_id+'_disall_points.csv','r')
	opfile=open('results/'+pdb_id+'_'+chain_id+'_expected_short_contacts.csv','w')
	shrt_conts=cPickle.load(open('data/disal_short_contacts.pkl','rb')) #contains details of short contacts of every disallowed phi-psi according to classical map
	for line in infile:
		lineparts=line.split(',')
		phi=int(lineparts[0])
		psi=int(lineparts[1])
		opfile.write(lineparts[0]+','+lineparts[1]+','+lineparts[2]+','+lineparts[3]+'\n')
		contacts=shrt_conts[str(phi)+','+str(psi)] #get the short contact information and write to output
		for item in contacts:		
			opfile.write(item)
	opfile.close()
	infile.close()

def FindDistance(a1,a2):
	import math
	return round(math.sqrt(((a2[0]-a1[0])*(a2[0]-a1[0]))+((a2[1]-a1[1])*(a2[1]-a1[1]))+((a2[2]-a1[2])*(a2[2]-a1[2]))),3)

def CheckShortContacts(a1,a2,d1):
	dist=FindDistance(a1,a2)
	if dist<d1:
		return "disallowed: "+str(dist)
	else:
		return "allowed: "+str(dist)

def find_actual_distances(pdb_id,chain_id):
	import cPickle
	infile=open('results/'+pdb_id+'_'+chain_id+'_expected_short_contacts.csv','r')
	dist=cPickle.load(open('data/short_contact_distances.pkl','rb'))
	opfile=open('results/'+pdb_id+'_'+chain_id+'_short_contact_analysis.txt','w')
	opfile.write('residue_name,residue_number,phi,psi,atom_1(residue_number)...atom_2(residue_number),expected_distance,observed_distace,short_contact_limit\n')
	for line1 in infile:
		lineparts=line1.split(',')
		if '...' not in lineparts[0]: #this line contains info on residue name, number, phi-psi
			resname_initial=lineparts[2]
			residue_initial=int(lineparts[3])
			resname=lineparts[2]
			residue=int(lineparts[3])
			phi=lineparts[0]
			psi=lineparts[1]
			pdbfile=open(pdb_id+".pdb","r")
			atoms=dict()
			for line in pdbfile:
				if line[0:4]=="ATOM" and line[21]==chain_id and int(line[22:26])==residue and line[12:16].strip() in ["N","H","CA","CB","C","O","HA"]: #i th residue
					if line[16]!=" ":
						if line[16]=="A": #store only A conformation in case of multiple occupancy
							atom_type=line[12:16].strip()+"2"
							resname=line[17:20]
							x=float(line[30:38])
							y=float(line[38:46])
							z=float(line[46:54])
							atoms[atom_type]=[x,y,z]
					else:
						atom_type=line[12:16].strip()+"2"
						resname=line[17:20]
						x=float(line[30:38])
						y=float(line[38:46])
						z=float(line[46:54])
						atoms[atom_type]=[x,y,z]
				if line[0:4]=="ATOM" and line[21]==chain_id and int(line[22:26])==residue+1 and line[12:16].strip() in ["N","H","CA"]: #i+1 st residue
					if line[16]!=" ":
						if line[16]=="A":
							atom_type=line[12:16].strip()+"3"
							resname=line[17:20]
							x=float(line[30:38])
							y=float(line[38:46])
							z=float(line[46:54])
							atoms[atom_type]=[x,y,z]
					else:
						atom_type=line[12:16].strip()+"3"
						resname=line[17:20]
						x=float(line[30:38])
						y=float(line[38:46])
						z=float(line[46:54])
						atoms[atom_type]=[x,y,z]
				if line[0:4]=="ATOM" and line[21]==chain_id and int(line[22:26])==residue-1 and line[12:16].strip() in ["CA","C","O"]: #i-1 th residue
					if line[16]!=" ":
						if line[16]=="A":
							atom_type=line[12:16].strip()+"1"
							resname=line[17:20]
							x=float(line[30:38])
							y=float(line[38:46])
							z=float(line[46:54])
							atoms[atom_type]=[x,y,z]
					else:
						atom_type=line[12:16].strip()+"1"
						resname=line[17:20]
						x=float(line[30:38])
						y=float(line[38:46])
						z=float(line[46:54])
						atoms[atom_type]=[x,y,z]
			pdbfile.close()
		else: #lines with expected short contact info
			splitparts=lineparts[0].split('...')
			atom1=splitparts[0]
			atom2=splitparts[1]	
			if atom1 in ["CB","HA"]: #rename atoms according to convention
				atom1=atom1+"2"
			if atom2 in ["CB","HA"]:
				atom2=atom2+"2"
			try: #will throw exception if either atom1 or atom2 is missing
				#get residue numbers of atoms in short contact from the structure
				if int(atom1[-1])==1:
					atom1_print=atom1[0:-1]+"("+str(residue_initial-1)+")"
				elif int(atom1[-1])==2:
					atom1_print=atom1[0:-1]+"("+str(residue_initial)+")"
				elif int(atom1[-1])==3:
					atom1_print=atom1[0:-1]+"("+str(residue_initial+1)+")"
				if int(atom2[-1])==1:
					atom2_print=atom2[0:-1]+"("+str(residue_initial-1)+")"
				elif int(atom2[-1])==2:
					atom2_print=atom2[0:-1]+"("+str(residue_initial)+")"
				elif int(atom2[-1])==3:
					atom2_print=atom2[0:-1]+"("+str(residue_initial+1)+")"
				shrt=CheckShortContacts(atoms[atom1],atoms[atom2],dist[atom1+","+atom2]) #returns the short contact distance
				opfile.write(resname_initial+","+str(residue_initial)+","+phi+","+psi+","+atom1_print+"..."+atom2_print+","+lineparts[1][:-1]+","+shrt.split(":")[1].strip()+","+str(dist[atom1+","+atom2])+"\n")
			except KeyError, e: #missing atom in short contact pair
				opfile.write(resname_initial+","+str(residue_initial)+","+phi+","+psi+","+atom1_print+"..."+atom2_print+","+lineparts[1][:-1]+",Missing atom "+str(e)[1:-2]+","+str(dist[atom1+","+atom2])+"\n")
			

					

