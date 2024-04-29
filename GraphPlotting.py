# Generate the following plots
# Time period data and fit
# Ball drop data and fit
# Compare model actual curves

# First import the needed libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# First import the neccesary data
Data = pd.read_excel("./ExperimentalData.xlsx")
ExpHighTime1 = Data["HighTime1"].dropna()
ExpHighTime2 = Data["HighTime2"].dropna()
ExpHighTime3 = Data["HighTime3"].dropna()

ExpLowTime1 = Data["LowTime1"].dropna()
ExpLowTime2 = Data["LowTime2"].dropna()
ExpLowTime3 = Data["LowTime3"].dropna()

ExpMassPosHigh = Data["MassPosHigh"].dropna()
ExpMassPosLow = Data["MassPosLow"].dropna()



# ---------------------------- Average the data and plot Time curves ---------------

# Now take averages and uncertainties
ExpHighTimes = []
ExpHighUnc = []
ExpLowTimes = []
ExpLowUnc = []

for i in range(0,len(ExpHighTime1)):
    # First find the average time
    ExpHighTimes.append((ExpHighTime1[i] + ExpHighTime2[i] + ExpHighTime3[i])/3)

    # Now find the uncertainty in the time
    Mean = ExpHighTimes[i]
    ExpHighUnc.append((ExpHighTime1[i] - Mean)**2 + (ExpHighTime2[i] - Mean)**2 + (ExpHighTime3[i] - Mean)**2)
    ExpHighUnc[i] = (ExpHighUnc[i]/3)**0.5


for i in range(0,len(ExpHighTime1)):
    # First find the average time
    ExpLowTimes.append((ExpLowTime1[i] + ExpLowTime2[i] + ExpLowTime3[i])/3)

    # Now find the uncertainty in the time
    Mean = ExpLowTimes[i]
    ExpLowUnc.append((ExpLowTime1[i] - Mean)**2 + (ExpLowTime2[i] - Mean)**2 + (ExpLowTime3[i] - Mean)**2)
    ExpLowUnc[i] = (ExpLowUnc[i]/3)**0.5


# These are the predicted coeffecients by the model
# Curves = [0]t^4 + ... [4]
ModelHighCoeffecients = [1.153,-3.625,4.707,-2.431,2.311]
ModelLowCoeffecients = [0.07551,-0.4309,1.053,-0.7333,2.093]

# First plot the experimental time period curves

# Using the same routine as in the model, for fit curves
ExpHighCurveFit, Highcov = np.polyfit(ExpMassPosHigh,ExpHighTimes,4, cov = True)
ExpHighCurveFit = np.poly1d(ExpHighCurveFit)


ExpLowCurveFit, Lowcov = np.polyfit(ExpMassPosLow,ExpLowTimes,4, cov = True)
ExpLowCurveFit = np.poly1d(ExpLowCurveFit)

# Now plot the figure
plt.xlabel("Mass distance [cm]")
plt.ylabel("Time Period [s]")
plt.title("Experimental curves")
plt.errorbar(ExpMassPosHigh, ExpHighTimes, xerr = 0.02 , yerr = ExpHighUnc, marker = "x", ls = "none", color = "black")
plt.plot(ExpMassPosLow, ExpLowTimes, "x", ls = "none", color = "red")
plt.plot(ExpMassPosHigh, ExpHighCurveFit(ExpMassPosHigh), color = "black")
plt.plot(ExpMassPosLow, ExpLowCurveFit(ExpMassPosLow), color = "red")
plt.savefig("Experimental time curves")
plt.show()


#--------------------- Analyse Curves to find final value for g --------------------

# Now do Newton Raphson for g

# Solve for g
RootEquation = ExpHighCurveFit - ExpLowCurveFit
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

# Look for first crossover point
X0 = 0.3

# Performing the cycles
Check = 0
while Check < 5:
    # Checking for stationary points
    if fDash(X0,Coeff) > -0.1 and fDash(X0,Coeff) < 0.1 :
        #print("Stationary point found, enter a different X0")
        delta = 101
        
    # Performing a cycle
    X0 = X0 - (f(X0,Coeff)/fDash(X0,Coeff))
    delta = (f(X0,Coeff)/fDash(X0,Coeff))
    
    # checks if an acceptable answer is found
    if delta <= 0.000001:
        Check = Check + 1
    else:
        Check = 0


ExpCrossOver1 = X0 
ExpTimeCrossHighUnc = Highcov[0,0]**0.5 + Highcov[1,1]**0.5 + Highcov[2,2]**0.5 + Highcov[3,3]**0.5 + Highcov[4,4]**0.5  

# Look for second crossover point
X0 = 50

# Performing the cycles
Check = 0
while Check < 5:
    # Checking for stationary points
    if fDash(X0,Coeff) > -0.01 and fDash(X0,Coeff) < 0.01 :
        #print("Stationary point found, enter a different X0")
        delta = 101
        
    # Performing a cycle
    X0 = X0 - (f(X0,Coeff)/fDash(X0,Coeff))
    delta = (f(X0,Coeff)/fDash(X0,Coeff))
    
    # checks if an acceptable answer is found
    if delta <= 0.000001:
        Check = Check + 1
    else:
        Check = 0

ExpCrossOver2 = X0 

ExpTimeCrossLowUnc = Lowcov[0,0]**0.5 + Lowcov[1,1]**0.5 + Lowcov[2,2]**0.5 + Lowcov[3,3]**0.5 + Lowcov[4,4]**0.5

print(f"The first crossover occurs at {round(ExpCrossOver1,3)}cm")
print(f"The second crossover occurs at {round(ExpCrossOver2,3)}cm")

# Using the crossovers evaluate g and the uncertainties for g

TCrossHigh = []
TCrossHigh.append(ExpHighCurveFit(ExpCrossOver1))
TCrossHigh.append(ExpHighCurveFit(ExpCrossOver2))

TCrossLow = []
TCrossLow.append(ExpHighCurveFit(ExpCrossOver1))
TCrossLow.append(ExpHighCurveFit(ExpCrossOver2))



EffectLen = 99.39
EffectLenUnc = 0.01
TotalAv = 0
Residue = 0
gUncertainties = []
gValues = []
gFinal = 0
gTotal = 0

for i in range(0,2):
    CurrentHighT = TCrossHigh[i]
    CurrentLowT = TCrossLow[i]
    
    # Find g from time period values
    gHigh = EffectLen*(2*np.pi/CurrentHighT)**2
    gLow = EffectLen*(2*np.pi/CurrentLowT)**2
    
    # Now find uncertainties
    gHighUnc = ((EffectLenUnc * gHigh/EffectLen) + (ExpTimeCrossHighUnc * EffectLen * 8 * (np.pi)**2 / CurrentHighT))**0.5
    gLowUnc = ((EffectLenUnc * gLow/EffectLen) + (ExpTimeCrossLowUnc * EffectLen * 8 * (np.pi)**2 / CurrentLowT))**0.5
    
    print(f"From crossover point {i + 1} from high g = ( {round(gHigh/100,5)} +- {round(gHighUnc/100,5)} )m/s^2")
    print(f"From crossover point {i + 1} from Low g = ( {round(gLow/100,5)} +- {round(gLowUnc/100,5)} )m/s^2")
    
    gValues.append(gHigh)
    gValues.append(gLow)
    gUncertainties.append(gHighUnc)
    gUncertainties.append(gLowUnc)


# Finally find a weighted average using the uncertainties
for i in range(0,len(gUncertainties)):
    gFinal = gFinal + gValues[i]/(gUncertainties[i]**2)
    gTotal = gTotal + (gUncertainties[i])**-2
    
gFinal = gFinal/gTotal
gFinalUnc = (gTotal)**-0.5

print(f"Taking a weighted average, g = ( {round(gFinal/100,5)} +- {round(gFinalUnc/100,5)} )m/s^2")

# ---------------------------- Comparing with simulation ----------------------

# Now compare the model and the actual results
# The coeffecients were retrieved from a seperate program that fit a curve to the modelled data

def HighConfigModel(x):
    HighCo = [1.159,-3.643,4.724,-2.437,2.311]
    Point = 0
    for i in range(0,5):
        Point = Point + HighCo[i]*x**(4-i)
    return Point
        
def LowConfigModel(x):
    LowCo = [0.07551,-0.4309,1.053,-0.7333,2.093]
    Point = 0
    for i in range(0,5):
        Point = Point + LowCo[i]*x**(4-i)
    return Point

# Convert to metres
ExpMassPosHigh = ExpMassPosHigh/100
ExpMassPosLow = ExpMassPosLow/100

# Model curves
DiskLargePos = np.linspace(0,0.9,num = 1000)
plt.xlabel("Mass distance [m]")
plt.ylabel("Time Period [s]")
plt.plot(DiskLargePos, LowConfigModel(DiskLargePos), color = "red",label = "Low Config", ls = "--")
plt.plot(DiskLargePos, HighConfigModel(DiskLargePos), color = "black",label = "High Config", ls = "--")

# Experiental results
plt.errorbar(ExpMassPosHigh, ExpHighTimes, xerr = 0 , yerr = ExpHighUnc, marker = "x", ls = "none", color = "black")
plt.plot(ExpMassPosLow, ExpLowTimes, "x", ls = "none", color = "red")
plt.title("Model on Data")
plt.legend()   
plt.savefig("Model on data")     
plt.show()

# ---------------------------Plot the ball drop curve ---------------------

ExpBallHeights = Data["Height"].dropna()
ExpBallTime1 = Data["FallTime1"].dropna()
ExpBallTime2 = Data["FallTime2"].dropna()
ExpBallTime3 = Data["FallTime3"].dropna()


AverageTime = []
AverageTSq = []
TSqUnc = []
TimeUnc = []

# First average the fall times and prodoce a T^2 for each height
for i in range(0,len(ExpBallHeights)):
    AverageTime.append(0)
    TimeUnc.append(0)
    AverageTSq.append(0)
    TSqUnc.append(0)
    
    # First take the average of the three points
    AverageTime[i] = (ExpBallTime1[i] + ExpBallTime2[i] + ExpBallTime3[i])/3
    
    # Now find the deviation from the mean
    Mean = AverageTime[i]
    TimeUnc[i] = (ExpBallTime1[i] - Mean)**2 + (ExpBallTime2[i] - Mean)**2 + (ExpBallTime3[i] - Mean)**2
    TimeUnc[i] = (TimeUnc[i]/3)**0.5
    
    # Now convert T -> T^2 
    AverageTSq[i] = (AverageTime[i])**2
    TSqUnc[i] = 2*(TimeUnc[i]/AverageTime[i]) # Percentage uncertainty
    TSqUnc[i] = TSqUnc[i] * AverageTSq[i] # Convert to abs uncertainty
    

# Now plot the data alongside a straight line fit 
def str_line(x,m,c):
    return m*x + c # Straight line fit


popt, pcov = curve_fit(str_line, ExpBallHeights, AverageTSq, sigma=TSqUnc)
# popt = values for model, in order of entered into the function
# pcov = varience of values
Gradient = popt[0]
GradientError = pcov[0,0]**0.5

Intercept = popt[1]
InterceptError = pcov[1,1]**0.5


# Plotting the fitting curve
x_line = np.linspace(ExpBallHeights.min(),ExpBallHeights.max(),100) # evenly spaced x values to use

plt.title("Ball Drop")
plt.xlabel("Height [m]")
plt.ylabel("$T^2$ [s]")

plt.errorbar(ExpBallHeights, AverageTSq, yerr = TSqUnc, marker = "x", color = "red", ls = "none", label = "Data")
plt.plot(x_line, str_line(x_line,Gradient,Intercept),linestyle="--",color = "black", label = "straight line fit")
plt.legend()

plt.savefig("Ball Drop")

# Now find the value for g 
gBall = 2/Gradient

# Now find the uncertainty
gBallUnc = (GradientError/Gradient) * gBall

# The final results
print("")
print(f"From the ball drop experiment the valeu of g = ( {round(gBall,5)} +- {round(gBallUnc,5)} )m/s^2")

# Printing the values
#print("Gradient = ", round(Gradient,10), "+-", round(GradientError,10))
#print("Intercept = ", round(Intercept,8), "+-", round(InterceptError,2))
    
    








