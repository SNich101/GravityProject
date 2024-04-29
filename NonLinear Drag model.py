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

# Generate the positions
CurrentPos = 2 # This is the starting position (in radians)
CurrentDeriv = 0 # Starting velocity (in rad/s)
PrevMaximaTime = 0 # Tracks the most recent maxima observed
Theta = [CurrentPos] # Holds the angular displacement at time t
Velocity = [0] # Holds the velocity of the pendulumn  at a given time
Time = [0] # Holds the time for each point
CurrentTime = 0
RepPlots = 1 # This is used to make group plots (max 5)

for i in range(0,1000):
    NewDeriv = FindNewDeriv(CurrentDeriv,CurrentPos)
    NewPos = FindNewPos(NewDeriv,CurrentPos)
    #print(f"Current position = {NewPos}, Current Deriv = {NewDeriv}")
    CurrentTime = CurrentTime + Delta
    
    # Add the new position to the list
    Theta.append(NewPos)
    Velocity.append(NewDeriv)
    Time.append(CurrentTime)

    # Check for maxima
    #if (abs(NewDeriv) < 0.1) and (NewPos > 0):
        # Evaluate the time period
        #TPeriod = CurrentTime - PrevMaximaTime
        #print("Time Period :" + str(TPeriod))
        #print("Maxima at Time, t = " + str(CurrentTime))
        #print("Current Derivative: " + str(NewDeriv))
        #print("")
        #PrevMaximaTime = CurrentTime
        

    # Repeat the process
    CurrentPos = NewPos
    CurrentDeriv = NewDeriv


# Plot the results
MarkerColours = ["black","green","orange","purple","blue"]
Markers = ["x","+","1","2","."]

# Define figure
for i in range(0,RepPlots):
    plt.figure(1, figsize = (8,6))
    plt.xlabel("Time (s)")
    plt.ylabel("Angle (Rad)")
    plt.tick_params(direction = "in", length = 7)

    # custom options for my plots
    if i == 0:
        plt.title("Pendulumn Trajectory (Space)")
        plt.xlabel("Time [s]")
        plt.ylabel("$\Theta$ [Radians]")

    plt.plot(Time[0:200],Theta[0:200],ls = "none", marker = Markers[i], color = MarkerColours[i])
    plt.savefig(f"Space Trajectory Version {i}")
    plt.show()

    # Next plot the phase diagram
    plt.xlabel("Angle (Rad)")
    plt.ylabel("Velocity (Rad/s)")
    
    # custom options for my plots
    if i == 0:
        plt.title("Pendulumn Trajectory (Phase)")
        plt.xlabel("$\Theta$ [Radians]")
        plt.ylabel("Velocity [Rad/s]")
        
    plt.tick_params(direction = "in", length = 7)
    plt.plot(Theta,Velocity, color = MarkerColours[i])
    plt.savefig(f"Phase Trajectory Version {i}")
    plt.show()

# use constants given to plot results curves

# 











        
