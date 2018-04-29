#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle
from PIL import Image
import os

def point_inside_polygon(x,y,poly):
	n = len(poly)
	inside = False
	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y
	
	return inside

def plot_map(pdb_id,chain_id):
	phivalues=[]
	psivalues=[]
	texts=[]
	dis_phi=[]
	dis_psi=[]
	dis_labels=[]
	temp_poly=[[-144,18],[-107,-19],[-53,-19],[-82,18]] #boundary of bridge region
	phipsi_apd=cPickle.load(open("data/ideal_rmap_apd.pkl","rb"))
	phipsi_prob=cPickle.load(open("data/phipsi_probability.pkl","rb"))
	phipsi_file=open('results/'+pdb_id+'_'+chain_id+'_phipsi.csv','r')
	for line in phipsi_file:
		lineparts=line.split(",")
		if lineparts[1]!='GLY':
			phi=int(lineparts[2])
			psi=int(lineparts[3])
			if phipsi_apd[str(phi)+","+str(psi)]==0: #if a phi-psi is disallowed
				if not point_inside_polygon(phi,psi,temp_poly): #if the point is not in bridge region, append it to disallowed phi-psi
					dis_phi.append(phi)
					dis_psi.append(psi)
					dis_labels.append(lineparts[1]+","+lineparts[0])
				else: #point in bridge region
					phivalues.append(phi)
					psivalues.append(psi)
			else: #point is allowed
				phivalues.append(phi)
				psivalues.append(psi)
	im=plt.imread("data/ensemble_allstructs_probability_outline.png")
	plt.imshow(im, extent=[-180, 180, -180, 180],zorder=1) #ensemble probability map outline
	plt.scatter(phivalues,psivalues,zorder=2,s=30,alpha=0.7,c="k") #plot the allowed points
	if len(dis_phi)!=0: #plot the disallowed points in different marker
		plt.scatter(dis_phi,dis_psi,zorder=2,s=30,alpha=1,marker="d",c="k",edgecolors="red")
	plt.xlim(-180, 180)
	plt.ylim(-180, 180)
	plt.xlabel('phi (in $\degree$)',fontsize=12)
	plt.ylabel('psi (in $\degree$)',fontsize=12)
	plt.savefig('results/'+pdb_id+'_'+chain_id+'_rplot_temp.png',format='png',dpi=200,bbox_inches="tight")
	plt.close()
	opfile=open('results/'+pdb_id+'_'+chain_id+'_disall_points.csv','w') #store the disallowed phi and psi points
	for i in range(0,len(dis_phi)):
		opfile.write(str(dis_phi[i])+','+str(dis_psi[i])+','+dis_labels[i]+','+str(round(phipsi_prob[str(dis_phi[i])+','+str(dis_psi[i])],3))+'\n')
	opfile.close()
	phipsi_file.close()
	images_list=['results/'+pdb_id+"_"+chain_id+"_rplot_temp.png","data/rplot_legend_color_new.png"] #combine the phi-psi plot and legend images
	imgs=[Image.open(i) for i in images_list]
	img1=imgs[0].resize([610,607],Image.ANTIALIAS)
	img2=imgs[1].resize([196,607],Image.ANTIALIAS)
	new_im=Image.new('RGB',(806,607))
	new_im.paste(img1,(0,0))
	new_im.paste(img2,(610,0))
	new_im.save('results/'+pdb_id+'_'+chain_id+'_rplot.png')

def plot_prob(pdb_id,chain_id):
	phipsi_prob=cPickle.load(open("data/phipsi_probability.pkl","rb"))
	phipsi_file=open('results/'+pdb_id+'_'+chain_id+'_phipsi.csv','r')
	opfile=open('results/'+pdb_id+'_'+chain_id+'_phipsi_probabilities.csv','w')
	opfile.write('residue_number,residue_name,phi,psi,probability\n')
	probs=[]
	for line in phipsi_file:
		lineparts=line.split(",")
		if lineparts[1]!='GLY':
			phi=int(lineparts[2])
			psi=int(lineparts[3])
			probs.append(phipsi_prob[str(phi)+","+str(psi)]) #store probabilities of all non-gly residues
			opfile.write(line[0:-1]+','+str(phipsi_prob[str(phi)+","+str(psi)])+'\n')
	arr=plt.hist(probs,bins=[0,0.25,0.5,0.75,1],color="#3399ff",edgecolor="black") #make histogram with mentioned bins
	plt.xlabel("Probability",fontsize=20)
	plt.ylabel("No. of points",fontsize=20)
	plt.xticks([0,0.25,0.5,0.75,1],fontsize=18)
	plt.yticks(fontsize=18)
	total=0
	for i in range(4):
		total=total+int(arr[0][i]) #add up number of phi-psi points
	for i in range(4):
		to_print=(int(arr[0][i])/float(total))*100 #find percent of points in each bin
		to_print=round(to_print,1)
		plt.text(arr[1][i]+0.1,arr[0][i]+2,str(to_print)+"%",fontsize=15) #annotate each bin with percent
	plt.savefig('results/'+pdb_id+'_'+chain_id+'_probability_distribution.png',dpi=200,format="png")
	plt.close()	
	phipsi_file.close()
	opfile.close()

