import pandas as pd
import csv
import numpy as np
from datetime import datetime
import datetime
from math import radians, cos, sin, sqrt, atan, degrees
import itertools
from tqdm import tqdm
import time

#Loading data
df=pd.read_csv("01_data_mars_opposition_updated.csv")
#List to hold all the times 
times=list()
times.append(0)
#Finding times differences b/w next oppostions
for i in range(1,len(df)):
    date1=datetime.datetime(df["Year"][i-1],df["Month"][i-1],df["Day"][i-1],df["Hour"][i-1],df["Minute"][i-1])
    date2=datetime.datetime(df["Year"][i],df["Month"][i],df["Day"][i],df["Hour"][i],df["Minute"][i])
    diff=date2-date1
    numDays=diff.days+diff.seconds/86400
    times.append(numDays)
times=np.array(times)
#print (times)

#Defining a list for using parameter values obtained from last question
global paras
paras=list(np.ones(6))

#Calculating Heliocentric Longitudes
long_angles_degree=np.array(df["ZodiacIndex"]*30+df["Degree"]+df["Minute.1"]/60+df["Second"]/3600)
oppositions=np.stack((times,long_angles_degree),axis=1)
#global pos_global=list(np.ones(1))*6




def minimize(comb_list,oppositions):
    min_check=1000
    min_list=list()
    pos=0
    for it in tqdm(range(0,len(comb_list)), desc = 'tqdm() Progress Bar'):  
        c=comb_list[it][0]
        r=comb_list[it][3]
        e1=comb_list[it][4]
        e2=comb_list[it][5]
        z=comb_list[it][2]
        s=comb_list[it][1]        
        errors=list()
        time=0
        # let equant cartestian co-ordinate be (h,k)
        h=e1*cos(radians(e2+z))
        k=e1*sin(radians(e2+z))
        #Finding all the error angles
        fie1=z
        for i in range(len(oppositions)):
            time=oppositions[i][0]            
            fie=((s*time)%360)+fie1 #calulating angle of eqant spoke from aries reference
            fie=fie%360
            fie1=fie
            
            
            # calculating "l" parameter, which will aid in finding the intercept point
            h_dash=h-cos(radians(c))
            k_dash=k-sin(radians(c))
            
            b_dash=((h_dash*cos(radians(fie)))+(k_dash*sin(radians(fie))))*2
            c_dash=(h_dash**2)+(k_dash**2)-(r**2)
            try:
                sqr=sqrt((b_dash**2)-(4*1*c_dash))
            except:
                sqr=0
            root1=(-b_dash+sqr)/2
            root2=(-b_dash-sqr)/2  

            #Intercepts (points in the orbit model in cartesian system)    
            
            #Checking for Root1
            x_int_1=h+(root1*cos(radians(fie)))
            y_int_1=k+(root1*sin(radians(fie)))
            angle_eqMars_root1=atan(y_int_1/x_int_1)
            angle_eqMars_deg_root1=degrees(angle_eqMars_root1)                  
            
            #Converting radians to degrees keeping quadrants in mind
            if(y_int_1<0 and x_int_1>0):
                angle_eqMars_deg_root1=angle_eqMars_deg_root1%360
            if(y_int_1<0 and x_int_1<0):
                angle_eqMars_deg_root1=angle_eqMars_deg_root1+180
            if(y_int_1>0 and x_int_1>0):
                angle_eqMars_deg_root1=angle_eqMars_deg_root1
            if(y_int_1>0 and x_int_1<0):
                angle_eqMars_deg_root1=angle_eqMars_deg_root1+180
                
            if(y_int_1==0 and x_int_1<0):
                angle_eqMars_deg_root1=0
            if(y_int_1==0 and x_int_1>0):
                angle_eqMars_deg_root1=180
            if(y_int_1>0 and x_int_1==0):
                angle_eqMars_deg_root1=90
            if(y_int_1<0 and x_int_1==0):
                angle_eqMars_deg_root1=270
                
            ####################
            #Checking for Root2
            x_int_2=h+(root2*cos(radians(fie)))
            y_int_2=k+(root2*sin(radians(fie)))
            angle_eqMars_root2=atan(y_int_2/x_int_2)
            angle_eqMars_deg_root2=degrees(angle_eqMars_root2)
            
            #Converting radians to degrees keeping quadrants in mind
            if(y_int_2<0 and x_int_2>0):
                angle_eqMars_deg_root2=angle_eqMars_deg_root2%360
            if(y_int_2<0 and x_int_2<0):
                angle_eqMars_deg_root2=angle_eqMars_deg_root2+180
            if(y_int_2>0 and x_int_2>0):
                angle_eqMars_deg_root2=angle_eqMars_deg_root2
            if(y_int_2>0 and x_int_2<0):
                angle_eqMars_deg_root2=angle_eqMars_deg_root2+180
            if(y_int_2==0 and x_int_2<0):
                angle_eqMars_deg_root2=0
            if(y_int_2==0 and x_int_2>0):
                angle_eqMars_deg_root2=180
            if(y_int_2>0 and x_int_2==0):
                angle_eqMars_deg_root2=90
            if(y_int_2<0 and x_int_2==0):
                angle_eqMars_deg_root2=270

            # Finding the best error for both candidate root solutions
            err_root1=(angle_eqMars_deg_root1 -  oppositions[i][1])
            err_root2=(angle_eqMars_deg_root2 -  oppositions[i][1])
            if(abs(err_root1)<abs(err_root2)):
                err=abs(err_root1)
            else:
                err=abs(err_root2)            
            #print("First-", angle_eqMars_deg_root1,y_int_1,x_int_1, "Second-", angle_eqMars_deg_root2,y_int_2,x_int_2, "Longitude-", oppositions[i][1])
            
            errors.append(err)
            
        #print("\nIteration-",it,"Error=",errors,"\n")
        if(max(errors)<0.067):
            pos=it # To obtain the best parameters position
            min_list=errors
            break
        else:
            avg=sum(errors)/len(errors)
            if(avg<=min_check):
                min_list=errors
                pos=it       
                min_check=avg
    parameters=comb_list[pos]
    return min_list,parameters


# Question 1


def MarsEquantModel(c,r,e1,e2,z,s,oppositions):    
    # Using Grid Search to find the optimal values
    c_new=np.arange(c,c+1.5,0.05) 
    s_new=np.arange(s,s+0.01,0.01) 
    z_new=np.arange(z,z+1,0.05) 
    r_new=np.arange(r,r+0.5,0.05) 
    e1_new=np.arange(e1,e1+0.2,0.05) 
    e2_new=np.arange(e2,e2+1,0.05) 

    import itertools
    comb = [list(c_new),list(s_new),list(z_new),list(r_new),list(e1_new),list(e2_new)]

    #Obtaining all the combinations
    comb_list=list(itertools.product(*comb))
    len(comb_list)


    #pos contains the position of best parameters in "comb" list
    errors,parameters=minimize(comb_list,oppositions)
    maxError=max(errors)
    global paras
    paras=parameters
    return errors,maxError

# Calling function for fist question
# Inital value guessing 

from decimal import Decimal

print("Computing Question 1...\n")
# Input the inital guess of parameters
c_str = input("Inital guess of c: ")
c=float(c_str)
s_str = input("Inital guess of s: ")
s=float(s_str)
z_str = input("Inital guess of z: ")
z=float(z_str)
r_str = input("Inital guess of r: ")
r=float(r_str)
e1_str = input("Inital guess of e1: ")
e1=float(e1_str)
e2_str = input("Inital guess of e2: ")
e2=float(e2_str)
"""
c=148.5 
s=360/687 
z=55.5 
r=8.5 
e1=1.6 
e2=92 
"""


errors,maxError  = MarsEquantModel(c,r,e1,e2,z,s,oppositions)
print("Question 1---------\n","All Errors in list-\n",errors,"\nMaximum Error-\n",maxError)


# Question 2



def bestOrbitInnerParams(r,s,oppositions):

    # Using Grid Search to find the optimal values
    
    # To search around the parameters obtained from Question 1
    global paras
    c=paras[0]
    z=paras[2]
    e1=paras[4]
    e2=paras[5]
    #Fixing r and s
    s_new=list()
    s_new.append(s)
    r_new=list()
    r_new.append(r)
    
    c_new=np.arange(c-0.1,c+0.1,0.05)
    z_new=np.arange(z-0.1,z+0.1,0.05)
    e1_new=np.arange(e1-0.1,e1+0.1,0.05)
    e2_new=np.arange(e2-0.1,e2+0.1,0.05)
    
    import itertools
    comb = [list(c_new),list(s_new),list(z_new),list(r_new),list(e1_new),list(e2_new)]
    #Obtaining all the combinations
    comb_list=list(itertools.product(*comb))
    #pos contains the position of best parameters in "comb" list
    # Calling the minimize function that returns the best goodness of fit angles
    
    errors,parameters=minimize(comb_list,oppositions)
    maxError=max(errors)
    c=parameters[0]
    e1=parameters[4]
    e2=parameters[5]
    z=parameters[2]    
    
    #Saving the latest obtained parameters in the global parameter for the next question
    paras=parameters
    
    return c,e1,e2,z,errors,maxError


# Calling function for Question 2
#Fixing s and r with the parameters obtained fom question 1

print("Computing Question 2...\n")

s=paras[1]
r=paras[3]

c,e1,e2,z,errors,maxError=bestOrbitInnerParams(r,s,oppositions)
print("\nQuestion 2---------\n","All Errors in list-\n",errors,"\nMaximum Error-\n",maxError,"\nValues of c,e1,e2,z-\n",c,e1,e2,z)

# Question 3



def bestS(r,oppositions):
    # Using Grid Search to find the optimal values
    #Using s from question 2 output
    global paras
    c=paras[0]
    z=paras[2]
    e1=paras[4]
    e2=paras[5]    
    s=paras[1]

    #Fixing r, c, z, e1 and e2
    r_new=list()
    r_new.append(r)
    c_new=list()
    c_new.append(c)
    z_new=list()
    z_new.append(z)
    e1_new=list()
    e1_new.append(e1)
    e2_new=list()
    e2_new.append(e2)
    
    #Ranging s for discretized searcg
    s_new=np.arange(s-0.05,s+0.05,0.005)

    import itertools
    comb = [list(c_new),list(s_new),list(z_new),list(r_new),list(e1_new),list(e2_new)]
    #Obtaining all the combinations
    comb_list=list(itertools.product(*comb))
    #pos contains the position of best parameters in "comb" list
    # Calling the minimize function that returns the best goodness of fit angles
    errors,parameters=minimize(comb_list,oppositions)
    maxError=max(errors)
    s=parameters[1]
    paras=parameters
    return s,errors,maxError
    
# Calling function for Question 3  

print("Computing Question 3...\n")
  
r=paras[3]
s,errors,maxError=bestS(r,oppositions)
print("\nQuestion 3---------\n","All Errors in list-\n",errors,"\nMaximum Error-\n",maxError,"\nValues of s-\n",s)


# Question 4


def bestR(s,oppositions):
    # Using Grid Search to find the optimal values
    #Using s from question 2 output
    global paras
    c=paras[0]
    z=paras[2]
    e1=paras[4]
    e2=paras[5]
    r=paras[3]

    #Fixing s, c, z, e1 and e2
    s_new=list()
    s_new.append(s)
    c_new=list()
    c_new.append(c)
    z_new=list()
    z_new.append(z)
    e1_new=list()
    e1_new.append(e1)
    e2_new=list()
    e2_new.append(e2)
    
    #Ranging r for discretized searcg
    r_new=np.arange(r-0.05,r+0.05,0.005)

    import itertools
    comb = [list(c_new),list(s_new),list(z_new),list(r_new),list(e1_new),list(e2_new)]
    #Obtaining all the combinations
    comb_list=list(itertools.product(*comb))
    #pos contains the position of best parameters in "comb" list
    # Calling the minimize function that returns the best goodness of fit angles
    errors,parameters=minimize(comb_list,oppositions)
    maxError=max(errors)
    r=parameters[3]
    paras=parameters
    return r,errors,maxError
    
# Calling function for Question 4   

print("Computing Question 4...\n")

s=paras[1]
r,errors,maxError=bestR(s,oppositions)
print("\nQuestion 4---------\n","All Errors in list-\n",errors,"\nMaximum Error-\n",maxError,"\nValues of r-\n",r)


# Question 5


def bestMarsOrbitParams(oppositions):
    # Using Grid Search to find the optimal values
    # For inital guessing, taking the parameter values obtained from question 4
    global paras
    c=paras[0]
    z=paras[2]
    e1=paras[4]
    e2=paras[5]
    r=paras[3]
    s=paras[1]

    #Fixing c, z, e1 and e2
    c_new=list()
    c_new.append(c)
    z_new=list()
    z_new.append(z)
    e1_new=list()
    e1_new.append(e1)
    e2_new=list()
    e2_new.append(e2)
    
    #Ranging r and s iteratively
    r_new=np.arange(r-0.05,r+0.05,0.005)
    s_new=np.arange(s-0.05,s+0.05,0.005)
    

    import itertools
    comb = [list(c_new),list(s_new),list(z_new),list(r_new),list(e1_new),list(e2_new)]
    #Obtaining all the combinations
    comb_list=list(itertools.product(*comb))
    #pos contains the position of best parameters in "comb" list
    # Calling the minimize function that returns the best goodness of fit angles
    errors,parameters=minimize(comb_list,oppositions)
    maxError=max(errors)
    paras=parameters
    c=paras[0]
    z=paras[2]
    e1=paras[4]
    e2=paras[5]
    r=paras[3]
    s=paras[1]
    return r,s,c,e1,e2,z,errors,maxError
    
# Calling function for Question 4   

print("Computing Question 5...\n")

r,s,c,e1,e2,z,errors,maxError = bestMarsOrbitParams(oppositions)
print("\nQuestion 5---------\n","All Errors in list-\n",errors,"\nMaximum Error-\n",maxError,"\nValues of r,s,c,e1,e2,z-\n",r,s,c,e1,e2,z)

