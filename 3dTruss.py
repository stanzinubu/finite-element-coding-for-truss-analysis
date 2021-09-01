
import numpy as np
import matplotlib.pyplot as plt
import math
# from mpl_toolkits.mplot3d import axes3d
A = 200
E = 1e3
scale=10
dof=3


conn = np.array(([1,2],[1,4],[2,3],[1,5],[2,6],[2,5],[2,4],
                  [1,3],[6,1],[6,3],[5,4],[3,4],[6,5],[10,3],
                  [7,6],[9,4],[8,5],[3,8],[4,7],[6,9],[5,10],
                  [3,7],[4,8],[5,9],[6,10]))
xcord = [-37.5,37.5,-37.5,37.5,37.5,-37.5,-100,100,100,-100]
ycord = [0,0,37.5,37.5,-37.5,-37.5,100,100,-100,-100]
zcord = [200,200,100,100,100,100,0,0,0,0]
if len(xcord)==len(ycord)==len(zcord):
    print('coordinates satiesfied')
te = len(conn)
tn = len(xcord)
le = []
cx = []
cy = []
cz = []
cons_el= []
sn_el = []
en_el = []
disp =np.ones((tn*dof,1))
#bcs=[0,1,2,3,4,5,6,7,8]
bcs = [18,19,20,21,22,23,24,25,26,27,28,29]

disp[bcs]=0
print('this is displacement value')
print(disp)  
force= np.zeros(len(xcord)*dof) 
force[0]=-1000
force[1]=1000
force[2]=-10
force[3]=1000
force[4]=1000
force[5]=10



print(force)
def truss():
    for i in range(te):
        a= conn[i,0]
        b= conn[i,1]
        x2 = xcord[b-1]
        x1 = xcord[a-1]
        y2 = ycord[b-1]
        y1 = ycord[a-1]
        z2 = zcord[b-1]
        z1 = zcord[a-1]        
        lei = math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        cxe = (x2-x1)/lei
        cye = (y2-y1)/lei
        cze = (z2-z1)/lei
        cons = A*E/lei
        cx.append(cxe)
        cy.append(cye)
        cz.append(cze)
        cons_el.append(cons)
        sn_el.append(a)
        en_el.append(b)
        le.append(lei) 
        #print(le[i])            
    ele_stiff_mat = []
    print(sn_el)
    print(en_el)
    for i in range(te):
        cxx= float(cx[i])**2
        cyy= float(cy[i])**2
        czz= float(cz[i])**2
        cxy= float(cy[i])*float(cx[i])
        cxz= float(cz[i])*float(cx[i])
        czy= float(cy[i])*float(cz[i])        
        ele = cons_el[i]*np.array([[cxx,cxy,cxz,-cxx,-cxy,-cxz],
                                    [cxy,cyy,-czy,-cxy,-cyy,-czy],
                                    [cxz,czy,czz,-cxz,-czy,-czz],
                                    [-cxx,-cxy,-cxz,cxx,cxy,cxz],
                                    [-cxy,-cyy,-czy,cxy,cyy,czy],
                                    [-cxz,-czy,-czz,cxz,czy,czz]])                                    
        ele_stiff_mat.append(ele)
    temp_ele = []  
     
    for i in range(te):
        m= sn_el[i]*2
        n= en_el[i]*2
        vec = [m-2,m-1,m,n-2,n-1,n]
    
        mat = np.zeros((tn*dof,tn*dof)) 
        ele_mat = ele_stiff_mat[i]
        for j in range(6):
            for k in range(6):
                a=vec[j]-1
                b=vec[k]-1
                mat[a,b]=ele_mat[j,k]
                
        temp_ele.append(mat)       
    #print(temp_ele)    
    GSM = np.zeros((tn*dof,tn*dof))
    for i in range(te):
        GSM =GSM+temp_ele[i]
    #print('this is global stiffness matrix')    
    #print(GSM)    
    
    #reduce stiffnees matrix
    rrGsm= np.delete(GSM,bcs,0)
    crGsm= np.delete(rrGsm,bcs,1)
    reGsm = crGsm
    print('this is reduce stiffnees matreix')
    print(reGsm)
    print(reGsm.shape)
    rforcemat = np.delete(force, bcs, 0)
    print(rforcemat)
    print(rforcemat.shape)
    dispresult = np.matmul(np.linalg.inv(reGsm),rforcemat)
    print('this is displacement matrix')
    print(np.around(dispresult,3))
    counter =0
    for i in range(tn*3):
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
    newZcoor = []
    count=0
    for i in range(tn):
        s=float(xcord[i]+disp[count]*scale)
        newXcoor.append(s)
        count=count+1
        t = float(ycord[i]+disp[count]*scale)
        newYcoor.append(t)
        count=count+1
        n= float(zcord[i]+disp[count]*scale)
        newZcoor.append(n)
        count = count+1
    print("this is new coordinates")    
    print(newXcoor)    
    print(newYcoor)
    print(newZcoor)
    
    ## new length of the member
    newlen= []
    
        
    for i in range(te):
        a,b = sn_el[i], en_el[i]
        x2 = float(newXcoor[b-1])
        x1 = float(newXcoor[a-1])
        y2 = float(newYcoor[b-1])
        y1 = float(newYcoor[a-1])
        z2 = float(newZcoor[b-1])
        z1 = float(newZcoor[a-1])
        
        l = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
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
        stress[i,0] = E*strain[i,0]
    print('this is stress value')
    print(stress)
    
    #member force
    memfor = np.zeros((te,1))
    for i in range(te):
        memfor[i,0] = A* stress[i,0]
    print("this is memeber force")
    print(memfor) 
    print('this is zcoor')
    print(newZcoor)
    #print(zcord)
    # newZcord = newZcoor[0].tolist()
    print(type(newZcoor))
    return strain,stress,memfor,newXcoor,newYcoor,newZcoor
    
 
def plot(conn,xcord,ycord,zcord,c,ls,lw,name):
    plt.gca(projection='3d')
    for i in range(te):
        a= conn[i,0]
        b= conn[i,1]
        x2 = xcord[b-1]
        x1 = xcord[a-1]
        y2 = ycord[b-1]
        y1 = ycord[a-1]
        z2= zcord[b-1]
        z1 = zcord[a-1]
        line,= plt.plot([x1,x2],[y1,y2],[z1,z2],color=c,linestyle=ls,linewidth=lw)  
    line.set_label(name)
    plt.legend(prop={'size':8})
    
 
if __name__=='__main__':    
    
    
    
    disp,forceresult=truss()    
    strain,stress,memfor,newXcoor,newYcoor,newZcoor=reactions()       
    plot(conn,xcord,ycord,zcord,'black','--',1,'undeformed')   
    plot(conn,newXcoor,newYcoor,newZcoor,'red','--',1.5,'deformed')       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





