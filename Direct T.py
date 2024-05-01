# Direct time period predictor

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Function to find the position of the COM
def CalcCOMass(x1,x2,M1,M2):
    COM = (1/(M1 + M2 + Mass_Rod)) *((Mass_Rod*Total_Length/2) + (M1*x1) + (M2*x2))
    return COM

def FindTotalICM(COM,L_CM_Small,L_CM_Large):
    Mass_DiskLarge = 1.369 # Mass of the larger disk
    Mass_DiskSmall = 1.004 # Mass of the smaller disk
    Mass_Rod= 1.582 # Mass of the entire rod

    Total_Length = 1.659 # total length of the rod
    End_Length = 0.34 # length between 
    Rod_Width = 0.01 # Width of the rod
    Disk_Radius = 0.047# Distance from to end to pivot

    # Compute the moment of inertia about centre of mass
    ICM = Mass_Rod*((1/12)*(Total_Length**2) + (Total_Length/2 - COM)**2)
    # Add component due to small disk
    ICM = ICM + ((1/2)*Mass_DiskSmall*(Disk_Radius)**2) + Mass_DiskSmall*(L_CM_Small**2) 
    # Add component due to large disk
    ICM = ICM + ((1/2)*Mass_DiskLarge*(Disk_Radius)**2) + Mass_DiskLarge*(L_CM_Large**2)
    # Return value
    return ICM

# Main code
# Valid constants
Mass_DiskLarge = 1.369 # Mass of the larger disk
Mass_DiskSmall = 1.004 # Mass of the smaller disk
Mass_Rod= 1.582 # Mass of the entire rod
Mass_Total = Mass_DiskLarge + Mass_DiskSmall + Mass_Rod

Total_Length = 1.659 # total length of the rod
End_Length = 0.34 # length between 
Rod_Width = 0.01 # Width of the rod
Disk_Radius = 0.047# Distance from to end to pivot
EffLength = 0.9939

g = 9.8112 # Test coeffecient of g

High_Position = 0.147 # Positon of the disk when fixed at high


# Configuration 1 positions (small on high)
DiskSmallPos = High_Position
DiskLargePos = np.linspace(0.1,0.9,num = 25)
TimePeriodsHigh = []
TimePeriodsLow = []

for i in range(0,len(DiskLargePos)):
    # First find the centre of mass (with ref to top of pole)
    COM = CalcCOMass(High_Position,(DiskLargePos[i]+End_Length),Mass_DiskSmall,Mass_DiskLarge)

    # Then calculate moment of inertia
    Dist_CM_Small = COM - 0.147 # Small in high position
    Dist_CM_Large = (End_Length + DiskLargePos[i]) - COM
    MOInertia = FindTotalICM(COM,Dist_CM_Small,Dist_CM_Large)
    #print("ICM = " + str(MOInertia))

    # Then use this to find the radius of gyration
    K_G = (MOInertia/(Mass_Total))**0.5

    # Find pivot -> COM distance
    L_CM = COM - End_Length
        
    # Now use this to find the time period
    T = 2*np.pi*((K_G**2 + L_CM**2)/(g*L_CM))**0.5

    # Display the time period
    TimePeriodsHigh.append(T)

    # for debugging purposes
    #print(f"At L = {DiskLargePos[i]} for high config MOI = {MOInertia} and LCM = {COM}")
    #print("  ")

# Configuration 2 positions (Small on low)
DiskSmallPos = Total_Length - 0.147
for i in range(0,len(DiskLargePos)):
    # First find the centre of mass (with ref to top of pole)
    COM = CalcCOMass(DiskSmallPos,(DiskLargePos[i]+End_Length),Mass_DiskSmall,Mass_DiskLarge)
    
    # Then calculate moment of inertia
    Dist_CM_Small = DiskSmallPos - COM # Small in high position
    Dist_CM_Large = abs((End_Length + DiskLargePos[i]) - COM)
    MOInertia = FindTotalICM(COM,Dist_CM_Small,Dist_CM_Large)
    #print("ICM = " + str(MOInertia))

    # Then use this to find the radius of gyration
    K_G = (MOInertia/(Mass_Total))**0.5

    # Find pivot -> COM distance
    L_CM = COM - End_Length
        
    # Now use this to find the time period
    T = 2*np.pi*((K_G**2 + L_CM**2)/(g*L_CM))**0.5

    # Display the time period
    TimePeriodsLow.append(T)

    #print(f"At L = {DiskLargePos[i]} for low config MOI = {MOInertia} and LCM = {COM}")
    #print("  ")


# Fit quadratic curves to the data

HighCurveModel = np.poly1d(np.polyfit(DiskLargePos,TimePeriodsHigh,4))
#plt.plot(DiskLargePos, HighCurveModel(DiskLargePos), color = "red",label = "High Config")
print(HighCurveModel)

LowCurveModel = np.poly1d(np.polyfit(DiskLargePos,TimePeriodsLow,4))
#plt.plot(DiskLargePos, LowCurveModel(DiskLargePos), color = "blue",label = "Low Config")
print(LowCurveModel)

# Solve for g
RootEquation = HighCurveModel - LowCurveModel
Coeff = RootEquation.coefficients

# Employ the NewtonRapson for root finding
def f(t,Coeff):
  return  Coeff[0]*t**4 + Coeff[1]*t**3 + Coeff[2]*t**2 + Coeff[3]*t**1 + Coeff[4]

def fDash(x,Coeff):
    h = 0.01
    D0 = (f(x + (h/2),Coeff) - f(x - (h/2),Coeff))/h
    DeltaDeriv = 1
    while DeltaDeriv >= 0.0001:
        Deriv = (f(x + (h/2),Coeff) - f(x - (h/2),Coeff))/h
        DeltaDeriv = ((D0 - Deriv)**2)**2
        D0 = Deriv
        h = h/10
    return Deriv

# Newton-Raphson method
X0 = 0.3

# Performing the cycles
Check = 0
while Check < 5:
    # Checking for stationary points
    if fDash(X0,Coeff) > -0.1 and fDash(X0,Coeff) < 0.1 :
        print("Stationary point found, enter a different X0")
        delta = 101
        
    # Performing a cycle
    X0 = X0 - (f(X0,Coeff)/fDash(X0,Coeff))
    delta = (f(X0,Coeff)/fDash(X0,Coeff))
    
    # checks if an acceptable answer is found
    if delta <= 0.000001:
        Check = Check + 1
    else:
        Check = 0
print(X0)

EffTimePeriod = HighCurveModel(X0)
Effg = 4*EffLength*(np.pi**2)/(EffTimePeriod)**2
print(f"Common Time period = {X0} => g = {Effg}")

# Plot the time periods
MarkerColours = ["red","green","orange","purple","blue"]
Markers = ["x","+","1","2","."]

for i in range(0,1):
    plt.plot(DiskLargePos, TimePeriodsHigh, Markers[i], color ='black')
    plt.plot(DiskLargePos, TimePeriodsLow, Markers[i], color = MarkerColours[i])
    plt.plot(DiskLargePos, LowCurveModel(DiskLargePos), color = MarkerColours[i],label = "Low Config", ls = "--")
    plt.plot(DiskLargePos, HighCurveModel(DiskLargePos), color = "black",label = "High Config", ls = "--")
    plt.legend()
    plt.title("Predicted Time period Curves")
    plt.xlabel("Pivot-Mass distance")
    plt.ylabel("Time Period")
    plt.savefig(f"Time period curves {i}")
    plt.show()

