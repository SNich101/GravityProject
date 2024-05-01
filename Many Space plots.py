# This is to plot and compare several trajectories
# primitive model

# libraries used in graph plotter 1.0
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Define global variables
Delta = 0.05 # This is the step taken each time
B = 3.154*10**-2 #3.154*10**-4
C = -9.81

# Find the new dervative (v)
def FindNewDeriv(CurrentDeriv,CurrentPos):
    if CurrentDeriv < 0:
            NewDeriv = CurrentDeriv + Delta*(C*np.sin(CurrentPos) + B*CurrentDeriv**2)
    else:
            NewDeriv = CurrentDeriv + Delta*(C*np.sin(CurrentPos) - B*CurrentDeriv**2)
    #print(NewDeriv)
    return NewDeriv

# Find the new position (Theta)
def FindNewPos(NewDeriv,CurrentPos):
    NewPos = CurrentPos + Delta*NewDeriv
    return NewPos

# Main program

MarkerColours = ["black","green","orange","purple","red"]
Markers = ["x","+","1","2","."]

# Define figure

plt.figure(1, figsize = (8,6))
plt.xlabel("Time (s)")
plt.ylabel("Angle (Rad)")
plt.tick_params(direction = "in", length = 7)
#plt.title("Pendulumn Trajectory (Space)")
plt.xlabel("Time [s]")
plt.ylabel("$\Theta$ [Radians]")

# This is to record the time periods found
MaximaTime = [[] , [] , [], [] , [] , [], [] , [] , [], []]
TPeriods = [[] , [] , [], [] , [] , [], [] , [] , [], []]

MaximaAve = [[] , [] , [], [] , [] , [], [] , [] , [], []]
TAverage = [[] , [] , [], [] , [] , [], [] , [] , [], []]

Diff = 0.2
# Generate the trajectories
for Plots in range(0,5):

    # Generate the positions
    CurrentPos = 0.1  + Diff*Plots # This is the starting position (in radians)
    CurrentDeriv = 0 # Starting velocity (in rad/s)
    PrevMaximaTime = 0 # Tracks the most recent maxima observed
    Theta = [CurrentPos] # Holds the angular displacement at time t
    Velocity = [0] # Holds the velocity of the pendulumn  at a given time
    Time = [0] # Holds the time for each point
    CurrentTime = 0
    RepPlots = 1 # This is used to make group plots (max 5)

    # This is for T perioc analysis
    MaximaTime[Plots] = []
    TPeriods[Plots] = []
    Points = 15000

    for i in range(0,Points):
        NewDeriv = FindNewDeriv(CurrentDeriv,CurrentPos)
        NewPos = FindNewPos(NewDeriv,CurrentPos)
        #print(f"Current position = {NewPos}, Current Deriv = {NewDeriv}")
        CurrentTime = CurrentTime + Delta
        
        # Add the new position to the list
        Theta.append(NewPos)
        Velocity.append(NewDeriv)
        Time.append(CurrentTime)

        # Check for maxima
        if (abs(NewDeriv) < 0.1) and (NewPos > 0):
            #Evaluate the time period
            TPeriod = CurrentTime - PrevMaximaTime
            PrevMaximaTime = CurrentTime
            if TPeriod > 0.5:
                #print(f"Plot number = {Plots}")
                #print("Time Period :" + str(TPeriod))
                #print("Maxima at Time, t = " + str(CurrentTime))
                #print("Current Derivative: " + str(NewDeriv))
                #print("")
                MaximaTime[Plots].append(CurrentTime)
                TPeriods[Plots].append(TPeriod)
                

        # Repeat the process
        CurrentPos = NewPos
        CurrentDeriv = NewDeriv

    if Plots in [0,4]:
        plt.plot(Time[0:200],Theta[0:200],ls = "none",marker = "x", color = MarkerColours[Plots], label = f"$\Theta_0 = ${round(0.1 + Diff*Plots,2)} Rad")

    # The aim of this section is to reduce the scattering in T via averaging#
    Seperation = 10
    for j in range(0,int(len(TPeriods[Plots])/Seperation)):
        #print(j)
        Sum = 0
        for x in range(0,Seperation):
            Sum = Sum + TPeriods[Plots][j*Seperation + x]
        # Now record the average point
        TAverage[Plots].append(Sum/Seperation)
        MaximaAve[Plots].append(MaximaTime[Plots][Seperation*j])
plt.legend()
plt.savefig("Pendulum Trajectory Space")
plt.show()

# Now average the time period over time


# Now I want to analyse the overall time period behavior
#print("Maxima recorded at : ",MaximaTime[0])
#print("Time periods found :  ",TPeriods[0])
plt.title("Pendulum decay over time")
plt.xlabel("Position of maxima [s]")
plt.ylabel("Time period recorded [s]")
for i in range(0,Plots + 1):
    plt.plot(MaximaAve[i][1:],TAverage[i][1:],ls = "none", marker = "x", label = f"$\Theta_0 = ${round(0.1 + Diff*i,3)} Rad")
plt.legend()
plt.savefig("TPeriod decay")
plt.show()



