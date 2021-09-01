# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 00:04:08 2021

@author: Stanzin Norbu
"""

import numpy as np
import matplotlib.pyplot as plt
import math
A = 1e-4
E = [500000,312000]
scale=1
dof=2
conn = np.array(([1,2],[2,3]))
xcord = [0,1,0]
ycord = [0,0,0.75226877894]

# conn = np.array(([1,2],[2,3],[3,4],[4,8],[7,8],[6,7],[5,6],[2,5],[1,6],[2,7],
#                  [3,6],[3,8],[4,7],[2,6],[7,3],[1,9],[2,9],[5,9],[6,9],[2,10],
#                  [3,10],[6,10],[7,10],[3,11],[4,11],[7,11],[8,11]))
# xcord = [0,1,2,3,0,1,2,3,0.5,1.5,2.5]
# ycord = [1,1,1,1,0,0,0,0,0.5,0.5,0.5] 
te = len(conn)
tn = len(xcord)
ndof= dof*tn
le = []
cos_el = []
sin_el = []
cons_el= []
sn_el = []
en_el = []
disp =np.ones((tn*2,1))
bcs = [0,1,4,5]
disp[bcs]=0
print('this is displacement value')
print(disp)  
force= np.zeros(len(xcord)*dof) 
force[3]=-1
print(force)
def truss():
    for i in range(te):
        a= conn[i,0]
        b= conn[i,1]
        x2 = xcord[b-1]
        x1 = xcord[a-1]
        y2 = ycord[b-1]
        y1 = ycord[a-1]
        
        lei = math.sqrt((x2-x1)**2+(y2-y1)**2)
        cos = (x2-x1)/lei
        sin = (y2-y1)/lei
        cons = A*E[i]/lei
        cos_el.append(cos)
        sin_el.append(sin)
        cons_el.append(cons)
        sn_el.append(a)
        en_el.append(b)
        le.append(lei) 
        #print(le[i])            
    ele_stiff_mat = []
    print(sn_el)
    print(en_el)
    for i in range(te):
        cc= float(cos_el[i])**2
        ss= float(sin_el[i])**2
        cs= float(sin_el[i])*float(cos_el[i])
        ele = cons_el[i]*np.array([[cc,cs,-cc,-cs],
                                   [cs,ss,-cs,-ss],
                                   [-cc,-cs,cc,cs],
                                   [-cs,-ss,cs,ss]])
        ele_stiff_mat.append(ele)
    temp_ele = []  
     
    for i in range(te):
        m= sn_el[i]*2
        n= en_el[i]*2
        vec = [m-1,m,n-1,n]
    
        mat = np.zeros((tn*2,tn*2)) 
        ele_mat = ele_stiff_mat[i]
        for j in range(4):
            for k in range(4):
                a=vec[j]-1
                b=vec[k]-1
                mat[a,b]=ele_mat[j,k]
                
        temp_ele.append(mat)       
    #print(temp_ele)    
    GSM = np.zeros((tn*2,tn*2))
    for i in range(te):
        GSM =GSM+temp_ele[i]
    #print('this is global stiffness matrix')    
    #print(GSM)    
    
    #reduce stiffnees matrix
    rrGsm= np.delete(GSM,bcs,0)
    crGsm= np.delete(rrGsm,bcs,1)
    reGsm = crGsm
    print('this is reduce stiffnees matreix')
    #print(reGsm)
    print(reGsm.shape)
    rforcemat = np.delete(force, bcs, 0)
    print(rforcemat)
    print(rforcemat.shape)
    dispresult = np.matmul(np.linalg.inv(reGsm),rforcemat)
    print('this is displacement matrix')
    print(np.around(dispresult,3))
    counter =0
    for i in range(tn*2):
        if disp[i,0]==1:
            disp[i,0]=dispresult[counter]
            counter=counter+1
    print('this is full displacement vector')
    print(disp)   

    
    ## force vector
    forceresult = np.matmul(GSM,disp)
    print('this is force vector')
    print(forceresult)
    return disp,forceresult
def reactions():
    ## _______new coordinates_____
    newXcoor = []
    newYcoor = []
    count=0
    for i in range(tn):
        s=xcord[i]+disp[count]*scale
        newXcoor.append(s)
        count=count+1
        t = ycord[i]+disp[count]*scale
        newYcoor.append(t)
        count=count+1
    print("this is new coordinates")    
    print(newXcoor)    
    print(newYcoor)
    
    ## new length of the member
    newlen= []
    
        
    for i in range(te):
        a,b = sn_el[i], en_el[i]
        x2 = float(newXcoor[b-1])
        x1 = float(newXcoor[a-1])
        y2 = float(newYcoor[b-1])
        y1 = float(newYcoor[a-1])
        
        l = np.sqrt((x2-x1)**2+(y2-y1)**2)
        newlen.append(l)
    ##strain 
    print("this is newlength value")
    print(newlen)
    strain = np.zeros((te,1))
    for i in range(te):
        strain[i,0] = (newlen[i]-le[i])/(le[i])  
        
    print('this is strain value -ve---compression and +ve--tension')
    print(strain)
    
    ##stress 
    stress = np.zeros((te,1))
    for i in range(te):
        stress[i,0] = E[i]*strain[i,0]
    print('this is stress value')
    print(stress)
    
    #member force
    memfor = np.zeros((te,1))
    for i in range(te):
        memfor[i,0] = A* stress[i,0]
    print("this is memeber force")
    print(memfor)  
    return strain,stress,memfor,newXcoor,newYcoor
    
 
def plot(conn,xcord,ycord,c,ls,lw,name):
    for i in range(te):
        a= conn[i,0]
        b= conn[i,1]
        x2 = xcord[b-1]
        x1 = xcord[a-1]
        y2 = ycord[b-1]
        y1 = ycord[a-1]
        line,= plt.plot([x1,x2],[y1,y2],color=c,linestyle=ls,linewidth=lw)  
    line.set_label(name)
    plt.legend(prop={'size':8})
    
 
if __name__=='__main__':    
    
    
    
    disp,forceresult=truss()    
    strain,stress,memfor,newXcoor,newYcoor=reactions()       
    plot(conn,xcord,ycord,'black','--',1,'undeformed')   
    plot(conn,newXcoor,newYcoor,'red','--',1.5,'deformed')       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    