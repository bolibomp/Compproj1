import matplotlib.pyplot as plt
import numpy as np
#Parameters
tol=1e-4 #Delta v
Rtrue=6.91 #Length of throw
H=3.18 #Maximum height
B2=4*np.pi*0.03**2*1.2*0.5/2 #A*density*C/2
m=0.07 #Mass of ball
deltaT=1e-3 #Delta t
g=9.82 #acceleration to ground


def fforw(array): #f for forwards in time
    v=np.sqrt(array[1]**2+array[3]**2)
    return np.array([array[1], -B2*array[1]*v/m, 
                     +array[3], -(g+B2*array[3]*v/m)])



def fback(array): #f for backwards in time
    v=np.sqrt(array[1]**2+array[3]**2)
    return np.array([-array[1], +B2*array[1]*v/m, 
                     -array[3], g+B2*array[3]*v/m])


v0x=1.5 #Starting trying at v0x=1.5 m/s (We know the answer is around 1.8
Rtot=[0] #, makes the code little faster)
rbacklist=[]
rforwlist=[]
while True:
    rback=np.array([0,v0x,H,0])#r is our vector
    rforw=np.array([0,v0x,H,0])
    xforw=[]
    yforw=[]
    xback=[]
    yback=[]
    
    while rforw[2]>0: #Runga Kutta forwards in time
        k1=deltaT*fforw(rforw)
        k2=deltaT*fforw(rforw+1/2*k1)
        k3=deltaT*fforw(rforw+1/2*k2)
        k4=deltaT*fforw(rforw+1*k3)
        rforw=rforw+k1/6+k2/3+k3/3+k4/6
        xforw.append(rforw[0])
        yforw.append(rforw[2])
        
    #Interpolation of last element
    a=yforw[-2]/(yforw[-2]-yforw[-1])
    rforw[0]=xforw[-2]+(xforw[-1]-xforw[-2])*a 
         
    while rback[2]>0: #Runga Kutta backwards in time
        k1=deltaT*fback(rback)
        k2=deltaT*fback(rback+1/2*k1)
        k3=deltaT*fback(rback+1/2*k2)
        k4=deltaT*fback(rback+1*k3)
        rback=rback+k1/6+k2/3+k3/3+k4/6
        xback.append(rback[0])
        yback.append(rback[2])
        
    #Interpolation of last element
    a=yforw[-2]/(yforw[-2]-yforw[-1])
    rforw[0]=xforw[-2]+(xforw[-1]-xforw[-2])*a 
         
    rbacklist.append(abs(rback[0]))
    rforwlist.append(rforw[0])
    Rtot.append(rforw[0]+abs(rback[0]))
    if abs(Rtot[-1]-Rtrue)>abs(Rtot[-2]-Rtrue): #Finds best initial velocity
        v0x=v0x-tol
        break
    v0x+=tol #Just adds Delta v
print(v0x)

# This part is for plotting and getting initial values

def forwardtrajectory(v0x): #we start where y=h and vy=0 
    rforw=np.array([0,v0x,H,0])
    posx = [rforw[0]]
    posy = [rforw[2]]
    while rforw[2]>0:
        k1=deltaT*fforw(rforw)
        k2=deltaT*fforw(rforw+1/2*k1)
        k3=deltaT*fforw(rforw+1/2*k2)
        k4=deltaT*fforw(rforw+1*k3)
        rforw=rforw+k1/6+k2/3+k3/3+k4/6
        posx.append(rforw[0])
        posy.append(rforw[2])
    return posx,posy
        
def backtrajectory(v0x): #we start where y=h and vy=0, reverse 
    rback=np.array([0,v0x,H,0])
    posx = [rback[0]]
    posy = [rback[2]]
    velox = [rback[1]]
    veloy= [rback[3]]
    while rback[2]>0:
        k1=deltaT*fback(rback)
        k2=deltaT*fback(rback+1/2*k1)
        k3=deltaT*fback(rback+1/2*k2)
        k4=deltaT*fback(rback+1*k3)
        rback=rback+k1/6+k2/3+k3/3+k4/6
        posx.append(rback[0])
        velox.append(rback[1])
        posy.append(rback[2])
        veloy.append(rback[3])
    velox[-1]=velox[-2]+(velox[-1]-velox[-2])*(posy[-2]/(posy[-2]-posy[-1]))
    veloy[-1]=veloy[-2]+(veloy[-1]-veloy[-2])*(posy[-2]/(posy[-2]-posy[-1]))
    return posx,velox[-1],posy, veloy[-1]

xforw, yforw = forwardtrajectory(v0x)
xvelo, yvelo = backtrajectory(v0x)[1], backtrajectory(v0x)[3]
xback, yback = backtrajectory(v0x)[0], backtrajectory(v0x)[2]

print(xvelo,yvelo)
plt.title("Trajectory")
plt.plot(xforw, yforw,"g--", label="vx=1E-3")
plt.plot(xback,yback,"g--" )
plt.xlabel("distance traveled (m)")
plt.ylabel("height (m)")
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

        
        