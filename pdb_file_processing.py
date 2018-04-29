#!/usr/bin/python

import math
import numpy as np
import os
def FindDihedralAngle(A,B,C,D):
    BCi=C[0]-B[0]
    BCj=C[1]-B[1]
    BCk=C[2]-B[2]

    BAi=A[0]-B[0]
    BAj=A[1]-B[1]
    BAk=A[2]-B[2]

    CDi=D[0]-C[0]
    CDj=D[1]-C[1]
    CDk=D[2]-C[2]

    Q1i=(BCj*BAk)-(BCk*BAj)
    Q1j=(BCk*BAi)-(BCi*BAk)
    Q1k=(BCi*BAj)-(BCj*BAi)

    Q2i=(BCj*CDk)-(BCk*CDj)
    Q2j=(BCk*CDi)-(BCi*CDk)
    Q2k=(BCi*CDj)-(BCj*CDi)
    magQ1=math.sqrt((Q1i*Q1i)+(Q1j*Q1j)+(Q1k*Q1k))
    Q1i=Q1i/magQ1
    Q1j=Q1j/magQ1
    Q1k=Q1k/magQ1

    magQ2=math.sqrt((Q2i*Q2i)+(Q2j*Q2j)+(Q2k*Q2k))
    Q2i=Q2i/magQ2
    Q2j=Q2j/magQ2
    Q2k=Q2k/magQ2

    Q1dotQ2=(Q1i*Q2i)+(Q1j*Q2j)+(Q1k*Q2k)
    chi=math.acos(Q1dotQ2)
    chinew=math.degrees(chi)
    
    Q1=np.array([Q1i,Q1j,Q1k])
    Q2=np.array([Q2i,Q2j,Q2k])
    Q1crossQ2=np.cross(Q1,Q2)
    magBC=math.sqrt((BCi*BCi)+(BCj*BCj)+(BCk*BCk))
    unitBCi=BCi/magBC
    unitBCj=BCj/magBC
    unitBCk=BCk/magBC
    unitBC=np.array([unitBCi,unitBCj,unitBCk])
    anglesign=np.dot(Q1crossQ2,unitBC)
    if anglesign<0:
        chinew=chinew*-1
        
    return int(round(chinew))

def findpsi(system):
    if system[0][6]==system[1][6] and system[1][6]==system[2][6] and system[2][6]==system[3][6]-1:
        A=[system[0][1],system[0][2],system[0][3]]
        B=[system[1][1],system[1][2],system[1][3]]
        C=[system[2][1],system[2][2],system[2][3]]
        D=[system[3][1],system[3][2],system[3][3]]
        psi=FindDihedralAngle(A,B,C,D)
        return psi
    else:
        return -999

def findphi(system):
    if system[2][6]==system[3][6]-1 and system[3][6]==system[4][6] and system[4][6]==system[5][6]:
        A=[system[2][1],system[2][2],system[2][3]]
        B=[system[3][1],system[3][2],system[3][3]]
        C=[system[4][1],system[4][2],system[4][3]]
        D=[system[5][1],system[5][2],system[5][3]]
        phi=FindDihedralAngle(A,B,C,D)
        return phi
    else:
        return -999

def findomega(system):
    if system[1][6]==system[2][6] and system[2][6]==system[3][6]-1 and system[3][6]==system[4][6]:
        A=[system[1][1],system[1][2],system[1][3]]
        B=[system[2][1],system[2][2],system[2][3]]
        C=[system[3][1],system[3][2],system[3][3]]
        D=[system[4][1],system[4][2],system[4][3]]
        omega=FindDihedralAngle(A,B,C,D)
        return omega
    else:
        return -999
        

def send_to_queue(atom,filename,opfile):
    global phiqueue
    global psiqueue
    global omegaqueue
    phiqueue.append(atom)
    psiqueue.append(atom)
    omegaqueue.append(atom)
    if atom[0]=="N":
        if len(psiqueue)==4:
            psi=findpsi(psiqueue)
            if psi!=-9999:
                opfile.write(filename+","+psiqueue[1][4]+","+str(psi)+","+str(psiqueue[1][6])+",psi\n")
            psiqueue=[]
            psiqueue.append(atom)       
    elif atom[0]=="CA":
        if  len(omegaqueue)==5:
            omega=findomega(omegaqueue)
            if omega!=-9999:
                opfile.write(filename+","+omegaqueue[1][4]+","+str(omega)+","+str(omegaqueue[1][6])+",omega\n")
            omegaqueue=[]
            omegaqueue.append(psiqueue[0])
            omegaqueue.append(psiqueue[1])
    elif atom[0]=="C":
        if len(phiqueue)==6:
            phi=findphi(phiqueue)
            if phi!=-9999:
                opfile.write(filename+","+phiqueue[4][4]+","+str(phi)+","+str(phiqueue[4][6])+",phi\n")
            phiqueue=[]
            phiqueue.append(psiqueue[0])
            phiqueue.append(psiqueue[1])
            phiqueue.append(psiqueue[2])

def checkMultipleModel(pdb_id):
    infile=open(pdb_id+".pdb","r")
    towrite=""
    flag=0
    for line in infile:
        if line[0:5]=="MODEL":
            flag=1
            continue
        if line[0:4]=="ATOM" and flag==0:
            break
        if flag==1 and line[0:6]=="ENDMDL":
            break
        elif flag==1:
            towrite+=line
    infile.close()
    if flag!=0:
        opfile=open(pdb_id+".pdb","w")
        opfile.write(towrite)
        opfile.close()

def cleanPdbFile(pdb_id,chain_id,missing_res):
    infile=open(pdb_id+".pdb","r")
    cur_res=None
    residue=dict()
    towrite=[]
    for line in infile:
        if line[0:4]=="ATOM" and line[21]==chain_id:
            if cur_res==None:
                cur_res=int(line[22:26])
            elif cur_res!=int(line[22:26]) and cur_res not in missing_res: #new residue and not a missing residue, write N, CA and C atom lines, in that order (order important for phi-psi calculation)
                towrite.append(residue['N'])
                towrite.append(residue['CA'])
                towrite.append(residue['C'])
                for item in residue.keys():
                    if item not in ['N','CA','C']: #write rest of the atom lines of the residue
                        towrite.append(residue[item])
                residue=dict()
                cur_res=int(line[22:26])
            if cur_res in missing_res: #if missing residue, ignore
                residue=dict()
                cur_res=int(line[22:26])
            if cur_res==int(line[22:26]):
                if not residue.has_key(line[12:16].strip()):
                    residue[line[12:16].strip()]=line #store atom lines of non-missing residues
        elif line[0:6]=="ANISOU": #ignore anisou lines
            continue
        else:
            towrite.append(line)
    infile.close()
    opfile=open(pdb_id+".pdb","w") #rewrite existing pdb file
    for line in towrite:
        opfile.write(line)
    opfile.close()

def findMissingRes(pdb_id,chain_id):            
    infile=open(pdb_id+".pdb","r")
    cur_res=None
    missing_res=[]
    atoms=dict()
    for line in infile:
        if line[0:4] == "ATOM" and line[21]==chain_id:
            if cur_res==None:
                cur_res=int(line[22:26])
            if cur_res!=int(line[22:26]): #new residue
                if cur_res+1!=int(line[22:26]): #residue number currently read is not +1 of previously read residue
                    for i in range(cur_res+1,int(line[22:26])):
                        missing_res.append(i)
                else:
                    added_atoms=atoms.keys()
                    if "N" not in added_atoms or "C" not in added_atoms or "CA" not in added_atoms: #if any of the three backbone atoms not present in a residue, consider it as missing
                        missing_res.append(cur_res)
                cur_res=int(line[22:26])
                atoms=dict()
            if line[21]==chain_id and line[12:16].strip() in ["N","CA","C"]:
                atom_type=line[12:16].strip()
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                if not atoms.has_key(atom_type): #store first conformer (in case of multiple occupancy)
                    atoms[atom_type]=[x,y,z]
    infile.close()
    return missing_res

def calculateDihedrals(pdb_id,chain_id):
    global phiqueue
    global psiqueue
    global omegaqueue
    phiqueue=[]
    psiqueue=[]
    omegaqueue=[]
    opfile=open('results/'+pdb_id+"_"+chain_id+"_phipsi_temp.csv","w")
    missing_res=findMissingRes(pdb_id,chain_id)
    cleanPdbFile(pdb_id,chain_id,missing_res)
    infile=open(pdb_id+".pdb","r")
    for line in infile:
        if "ATOM" in line[0:4]:
            if chain_id in line[21]:
                if int(line[22:26]) in missing_res: #ignore atoms from missing residues
                    continue
                if line[16]!=" ": #alternate conformation present
                    if line[16]=="A": #use only A conformation
                        cur_atom_name=line[12:16].strip()
                        if cur_atom_name=="N" or cur_atom_name=="C" or cur_atom_name=="CA":
                            atom_to_consider=[cur_atom_name,float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],float(line[54:60]),int(line[22:26])]
                            send_to_queue(atom_to_consider,pdb_id,opfile) 
                else:
                    cur_atom_name=line[12:16].strip()
                    if cur_atom_name=="N" or cur_atom_name=="C" or cur_atom_name=="CA":
                        atom_to_consider=[cur_atom_name,float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],float(line[54:60]),int(line[22:26])]
                        send_to_queue(atom_to_consider,pdb_id,opfile) 
    infile.close()
    opfile.close()
    if os.stat('results/'+pdb_id+"_"+chain_id+"_phipsi_temp.csv").st_size!=0:
        infile=open('results/'+pdb_id+"_"+chain_id+"_phipsi_temp.csv","r")
        opfile=open('results/'+pdb_id+"_"+chain_id+"_phipsi.csv","w") #combine phi and psi into single line
        cur_res=None
        cur_res_name=None
        for line in infile:
            lineparts=line.split(",")
            if cur_res==None:
                cur_res=int(lineparts[3])
                phi=" "
                psi=" "
                cur_res_name=lineparts[1]
            elif cur_res!=int(lineparts[3]):
                if phi!="-999" and psi!="-999" and phi!=" " and psi!=" ":
                    opfile.write(str(cur_res)+","+cur_res_name+","+phi+","+psi+"\n")
                phi=" "
                psi=" "
                cur_res=int(lineparts[3])
                cur_res_name=lineparts[1]
            if cur_res==int(lineparts[3]):
                if lineparts[4] == "phi\n":
                    phi=lineparts[2]
                elif lineparts[4] == "psi\n":
                    psi=lineparts[2]
        infile.close()
        opfile.close()
        return 1
    else: #no phi-psi could be calculated
        return 0

def addchainid(pdb_id):
    pdbfile=open(pdb_id+'.pdb','r')
    towrite=""
    for line in pdbfile:
        if line[0:3]=="TER":
            towrite+=line
            break
        if line[0:4]=="ATOM":
            if line[21]==" ":
                towrite+=line[:21]+"A"+line[22:]
            else:
                towrite+=line
    pdbfile.close()
    opfile=open(pdb_id+'.pdb','w')
    opfile.write(towrite)
    opfile.close()
    return towrite[21]
    

