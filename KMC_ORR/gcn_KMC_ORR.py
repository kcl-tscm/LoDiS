import math
import numpy as np
import collections as co
import time
import matplotlib.pyplot as plt
import os
import sys

def update_progress(progress):
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()



to_read= raw_input("Insert the name of the file containing the catalyst you wish to use for this simulation: ")
type(to_read)
start = time.time()

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    
    
length=file_len(to_read)

counterline=0
countercluster=0

omission = 0
testcon = 4
os.system("mkdir Output")
with open("gcn_bridge_genome.dat", "w")as o:
	with open(to_read, "r") as f:
		while counterline <length:
			countercluster+=1
			natoms=int(f.readline())
			secondline=f.readline()
			x=[]
			y=[]
			z=[]
			numN=[]
			numN=[0]*natoms
                        neighlist=[]
                        for i in range(natoms):
                            neighlist.append([])
			lines=[]
			cutoffNN=3.88*math.sqrt(2)/2*1.2
                        cna_max=12.0 # maximum coordination number of fcc crystal for single atoms
			cnb_max=18.0 # maximum coordination number of fcc crystal for pair of atoms
			num_pair=0
			for i in range(natoms):
				lines.append(f.readline())
			counterline+=2
			for i in lines:
				counterline+=1
				x.append(float(i.split()[1]))
				y.append(float(i.split()[2]))
				z.append(float(i.split()[3]))
			if testcon == 4:
                            print len(numN)
                            testcon = 5


				
			def dist(i,j):
				return math.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2)
			
			
			
			# Coordination number of the sites
			for i in range(natoms-1):
				for j in range(i+1, natoms):			
					if dist(i,j)<cutoffNN:
						numN[i]+=1
						numN[j]+=1
                                                neighlist[i].append(j)
                                                neighlist[j].append(i)
                            # Test that the recorded possible neighbours of the site correspond to the coordination number of the site
                                for i in range(natoms):
                                    cnmapper = numN[i]-len(neighlist[i])
                                    if cnmapper != 0:
                                        print cnmapper , i
                                        print "Error: neighbour tracking not conserved"

                        # Generalised coordination number of atop sites
                        def gcn_a():
                                gcn_a=[]

                                for i in range(natoms):
                                    neighb_a=[]
                                    for k in range(natoms):
                                        if k!=i:
                                            if dist(i,k)<cutoffNN:
                                                neighb_a.append(numN[k])
                                    neighbours_a=np.asarray(neighb_a)
                                    if numN[i] <10:
                                        gcn_a.append(np.sum(neighbours_a)/cna_max)
                                    else:
                                        gcn_a.append(0)

                                return gcn_a



                        # Generalised coordination number of bridge sites
			def gcn_b():
				gcn_b=[]
			        omission = 0

				for i in range(natoms-1):
					for j in range(i+1,natoms):
						if dist(i,j)<cutoffNN:
							neighb_b=[]
							for k in range(natoms):
								if i!=k and j!=k:
									if dist(i,k)<cutoffNN or dist(j,k)<cutoffNN:
										neighb_b.append(numN[k])
							neighbours_b=np.asarray(neighb_b)			
							if numN[i] <10 and numN[j]<10:
								gcn_b.append(np.sum(neighbours_b)/cnb_max)
							else:
								gcn_b.append(0)
                                                        if len(gcn_b) == 1:
                                                            print("GCN bridge genome being constructed")
                                                        update_progress((len(gcn_b)+omission)/(0.5*float(natoms**2 -natoms)))
                                                else: # used for accurate percentages during updating of the gcn genome progress
                                                    omission += 1

			
				return gcn_b
			
			

                        GCN_A=gcn_a()
			
			GCN_B=gcn_b()	
			
			singleco=co.Counter(GCN_A)

			
			counter=co.Counter(GCN_B)
			
			
				
			def cga1():
				nl=0
				nm=0
				nh=0
				for key in singleco:
					if key!=0:
						if key<=5.5:
							nl+=singleco[key]
						if key>5.5 and key<7.5:
							nm+=singleco[key]
						if key>=7.5:
							nh+=singleco[key]
				
				cg_list=[nl, nm, nh]
				return cg_list
			
			def cgb1():
				nl=0
				nm=0
				nh=0
				for key in counter:
					if key!=0:
						if key<=5.5:
							nl+=counter[key]
						if key>5.5 and key<7.5:
							nm+=counter[key]
						if key>=7.5:
							nh+=counter[key]
				
				cg_list=[nl, nm, nh]
				return cg_list
							
				
				
			def cgb2():
				nl=0
				nm=0
				nh=0
				nsh=0
				for key in counter:
					if key!=0:
						if key<=5.5:
							nl+=counter[key]
						elif key>=8.3:
							nsh+=counter[key]
						elif (key>=7.5 and key<8.3) or (key>=6.59 and key<=6.61):
							nh+=counter[key]
						else:
							nm+=counter[key]
				
				cg_list=[nl, nm, nh, nsh]
				return cg_list
				
							
			def cgb3():
				nl=0
                                nlm=0
				nmm=0
                                nmh=0
				nh=0
				for key in counter:
                                        print key
					if key!=0:
						if key<=5.2:
							nl+=counter[key]
						if key>5.2 and key<6:
							nlm+=counter[key]
                                                if (key>6 and key<6.6) or (key>6.7 and key<7):
                                                        nmm+=counter[key]
                                                if (key>6.6 and key<6.63) or (key>7.1 and key<7.3):
                                                        nmh+=counter[key]
						if (key>6.63 and key<6.7) or (key>7.3):
							nh+=counter[key]

				cg_list=[nl, nlm, nmm, nmh, nh]
				return cg_list
	
				
			def cgb1_and_2():
				nl1=0
				nm1=0
				nh1=0
				nsh1=0
				for key in counter:
					if key!=0:
						if key<=5.5:
							nl1+=counter[key]
						if key>5.5 and key<7.5:
							nm1+=counter[key]
						if key>=7.5:
							nh1+=counter[key]
							
				nl2=0
				nm2=0
				nh2=0
				nsh2=0
				for key in counter:
					if key!=0:
						if key<=5.5:
							nl2+=counter[key]
						elif key>=8.3:
							nsh2+=counter[key]
						elif (key>=7.5 and key<8.3) or (key>=6.59 and key<=6.61):
							nh2+=counter[key]
						else:
							nm2+=counter[key]		
				
				print "Cluster #",countercluster
				print "\nLow gcn sites: ", nl1, "  ",nl2
				print "\nMedium gcn sites: ", nm1, "  ",nm2
				print "\nHigh gcn sites: ", nh1, "  ",nh2
				print "\nSuper-high gcn sites: ", nsh2
				print "\nTotal surface sites: ", nl1+nm1+nh1+nsh1
				print "\n**********************************************************"
				
				if countercluster==1:
					o.write("Coarse graining mode=1 and 2")
					
				
				o.write("\n\nCluster #"+str(countercluster))
				o.write("\n\nGCN fingerprint:\n")
				for key in counter:
				
					if key!=0:
						o.write(str( key))
						o.write(str("	"))
						o.write(str(counter[key]))
						o.write(str("\n"))
				o.write("\n\nLow gcn sites: "+ str( nl1)+", "+str(nl2))
				o.write("\n\nMedium gcn sites: "+ str( nm1)+", "+str(nm2))
				o.write("\n\nHigh gcn sites: "+ str( nh1)+", "+str(nh2))
				o.write("\n\nSuper-high gcn sites: "+ str( nsh1))
				o.write("\n\nTotal surface sites: "+ str( nl1+nm1+nh1))
				o.write("\n\n********************************************************")
							
				#y1=[nl1,nm1,nh1,nsh1]
				y1=[nl1*100.0/(nl1+nm1+nh1+nsh1),nm1*100.0/(nl1+nm1+nh1+nsh1),nh1*100.0/(nl1+nm1+nh1+nsh1),nsh1*100.0/(nl1+nm1+nh1+nsh1)]
				#y2=[nl2,nm2,nh2,nsh2]
				y2=[nl2*100.0/(nl2+nm2+nh2+nsh2),nm2*100.0/(nl2+nm2+nh2+nsh2),nh2*100.0/(nl2+nm2+nh2+nsh2),nsh2*100.0/(nl2+nm2+nh2+nsh2)]
				x=["Low gcn","Medium gcn","High gcn","Super-high gcn"]
				y_pos=np.arange(4)
				plt.bar(y_pos, y1, color="b", width=0.25 , align="edge")
				plt.bar(y_pos+0.25, y2, color="r", width=0.25, align="edge")
				plt.ylabel("Percentage of sites")
				#plt.ylabel("Number of sites")
				plt.xticks(y_pos, x)
				plt.title("GCN coarse-graining")
				#plt.text("Low gcn sites: ", nl1," ", nl2,"\nMedium gcn sites: ", nm1," ", nm2,"\nHigh gcn sites: ", nh1," ", nh2,"\nSuper-high gcn sites: ", nsh1)
				plt.savefig("snap"+str(countercluster))
				plt.close()			
			
			
			#Which coarse-graining?
                        ab=2 # atop or bridge sites (1 is atop, 2 is bridge)
			cgb=3
			
			if cgb==1 and ab==2:
				coarse_grain=cgb1()
				nl=coarse_grain[0]
				nm=coarse_grain[1]
				nh=coarse_grain[2]
				print "Cluster #",countercluster
				print "\nLow gcn sites: ", nl
				print "\nMedium gcn sites: ", nm
				print "\nHigh gcn sites: ", nh
				print "\nTotal surface sites: ", nl+nm+nh
				print "\n**********************************************************"
				
				if countercluster==1:
					o.write("Coarse graining mode=1")
				
				with open("inputgcn.dat", "w") as p:
					p.write(str(nl)+"\n")
					p.write(str(nm)+"\n")
					p.write(str(nh)+"\n")
					p.write(str(countercluster))
				
	
				o.write("\n\nCluster #"+str(countercluster))
				o.write("\n\nGCN fingerprint:\n")
				for key in counter:
				
					if key!=0:
						o.write(str( key))
						o.write(str("	"))
						o.write(str(counter[key]))
						o.write(str("\n"))
				o.write("\n\nLow gcn sites: "+ str( nl))
				o.write("\n\nMedium gcn sites: "+ str( nm))
				o.write("\n\nHigh gcn sites: "+ str( nh))
				o.write("\n\nTotal surface sites: "+ str( nl+nm+nh))
				o.write("\n\n********************************************************")
				
				
				x=["Low gcn","Medium gcn","High gcn"]
				y=[nl*100.0/(nl+nm+nh),nm*100.0/(nl+nm+nh),nh*100.0/(nl+nm+nh)]
				#y=[nl,nm,nh]
				y_pos=np.arange(len(x))
				plt.bar(y_pos, y, color="b",width=0.50,align="center")
				plt.title("GCN coarse-graining")
				plt.ylabel("Percentage of sites")
				#plt.ylabel("Number of sites")
	
				plt.xticks(y_pos, x)
				
				plt.savefig("snap"+str(countercluster))
				plt.close()
			
			
			elif cgb==2 and ab==2:
				coarse_grain=cgb2()
				nl=coarse_grain[0]
				nm=coarse_grain[1]
				nh=coarse_grain[2]
				nsh=coarse_grain[3]
				print "Cluster #",countercluster			
				print "\nLow gcn sites: ", nl
				print "\nMedium gcn sites: ", nm
				print "\nHigh gcn sites: ", nh
				print "\nSuper-high gcn sites: ", nsh
				print "\nTotal surface sites: ", nl+nm+nh+nsh
				print "\n**********************************************************"
				if countercluster==1:
					o.write("Coarse graining mode=2")
				
				with open("inputgcn.dat", "w") as p:
					p.write(str(nl)+"\n")
					p.write(str(nm)+"\n")
					p.write(str(nh)+"\n")
				
				o.write("\n\nCluster #"+str(countercluster))
				o.write("\n\nGCN fingerprint:\n")
				for key in counter:
				
					if key!=0:
						o.write(str( key))
						o.write(str("	"))
						o.write(str(counter[key]))
						o.write(str("\n"))
				o.write("\n\nLow gcn sites: "+ str( nl))
				o.write("\n\nMedium gcn sites: "+ str( nm))
				o.write("\n\nHigh gcn sites: "+ str( nh))
				o.write("\n\nSuper-high gcn sites: "+ str( nsh))
				o.write("\n\nTotal surface sites: "+ str( nl+nm+nh))
				o.write("\n\n********************************************************")
				
				
				x=["Low gcn","Medium gcn","High gcn","Super-high gcn"]
				#y=[nl*100.0/(nl+nm+nh),nm*100.0/(nl+nm+nh+nsh),nh*100.0/(nl+nm+nh+nsh),nsh*100.0/(nl+nm+nh+nsh)]
				y=[nl,nm,nh,nsh]
				y_pos=np.arange(len(x))
				plt.bar(y_pos, y, color="r",width=0.50, align="center")
				plt.title("GCN coarse-graining")
				#plt.ylabel("Percentage of sites")
				plt.ylabel("Number of sites")
				plt.xticks(y_pos, x)
				
				plt.savefig("snap"+str(countercluster))
				plt.close()
			
	    		elif cgb==3 and ab==2:
				coarse_grain=cgb3()
				nl=coarse_grain[0]
                                nlm=coarse_grain[1]
				nmm=coarse_grain[2]
				nmh=coarse_grain[3]
                                nh=coarse_grain[4]
				print "Cluster #",countercluster
				print "\nLow gcn sites: ", nl
				print "\nMedium gcn sites: ", nlm+nmm
				print "\nHigh gcn sites: ", nmh+nh
				print "\nTotal surface sites: ", nl+nlm+nmm+nmh+nh
				print "\n**********************************************************"
				
				if countercluster==1:
					o.write("Coarse graining mode=1")
				
				with open("inputgcn.dat", "w") as p:
					p.write(str(nl)+"\n")
					p.write(str(nlm+nmm)+"\n")
					p.write(str(nmh+nh)+"\n")
					p.write(str(countercluster)+"\n")
                                        p.write(str(cga1()[0]+cga1()[1]+cga1()[2]))
				
	
				o.write("\n\nCluster #"+str(countercluster))
				o.write("\n\nGCN fingerprint:\n")
				for key in counter:
				
					if key!=0:
						o.write(str( key))
						o.write(str("	"))
						o.write(str(counter[key]))
						o.write(str("\n"))
				o.write("\n\nLow gcn sites: "+ str( nl))
				o.write("\n\nMedium gcn sites: "+ str( nlm+nmm))
				o.write("\n\nHigh gcn sites: "+ str( nmh+nh))
				o.write("\n\nTotal surface sites: "+ str( nl+nlm+nmm+nmh+nh))
                                o.write("\n\nTotal number of sites using atop description:"+ str(cga1()[0]+cga1()[1]+cga1()[2]))
				o.write("\n\n********************************************************")
				
				
				x=["Low gcn","Medium gcn","High gcn"]
				y=[nl*100.0/(nl+nlm+nmm+nmh+nh),(nlm+nmm)*100.0/(nl+nlm+nmm+nmh+nh),(nmh+nh)*100.0/(nl+nlm+nmm+nmh+nh)]
				#y=[nl,nm,nh]
				y_pos=np.arange(len(x))
				plt.bar(y_pos, y, color="b",width=0.50,align="center")
				plt.title("GCN coarse-graining")
				plt.ylabel("Percentage of sites")
				#plt.ylabel("Number of sites")
	
				plt.xticks(y_pos, x)
				
				plt.savefig("snap"+str(countercluster))
				plt.close()
			
			elif cgb==12 and ab==2:
				cgb1_and_2()	
				
			os.system("python KMC_ORRnew.py")
end = time.time()
print "\n\nTime it took for the program to run"
print(end - start)
			
