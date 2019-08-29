# 2019.05.20
# Christian D. Ahrberg - BNTL - Sogang University
#
# To be run using Python 2.7 / Python 3.0 and above require a different format of sending commands to the syringe pumps
#
# Python code for automated optimization of gold@ironoxide nanoparticles
#

# Importing of required packages
import time
import serial
import os
import RPi.GPIO as GPIO
from scipy import signal
import numpy as np

# Initialize GPIO of raspberry pi
GPIO.setmode(GPIO.BCM)  # Setting pin numbering of GPIOs to board
DEBUG = 1               # Turning on debug mode

# Definition of classes used within this code

class Filter():
    # Class for a low pass filter, used to remove high frequency noise from the transmission data

    # Initializing the filter class
    def __init__(self,samplingrate,cutofrate,treshold,*kwargs):
        fs = samplingrate                           # Frequency for sampling data in Hz
        fc = cutofrate                              # Cut of frequency for filter in Hz
        self.threshold = treshold                   # Transmission threshold to differentiate between droplets and oil phase

        # Calculating further parameters from inputs
        w = fc / (fs / 2)                           # Normalizing frequencies
        self.b ,self.a = signal.butter(5, w, 'low') # Butterworth filer design of scipy package

    def Lowpass(self,data,*kwargs):
        # performing low pass filter on provided data and adding filtered
        # data, and mean of data below threshold to instance
        # data should be a list with just one column (transmission)

        # List to add droplet data to
        sample = []

        # filtering data using forward-backward filter of scipy package
        self.filterdata = signal.filtfilt(self.b,self.a,data)

        # Applying threshold to extract data of droplets
        for j in self.filterdata:   # Iterating through all data points
            if j < self.threshold:  # If transmission below threshold it is droplet
                sample.append(j)    # Adding droplet data to previously defined list

        self.average = np.mean(sample) # Calculating the mean transmission of extracted droplet data

class Simplex():
    # Class for 2D simplex algorithm trying to find the minimum of a curve
    #
    # Should be used as follows:
    # 1. measure first two conditions
    # 2. Initialise Simplex class
    # 3. Find optimum using Simplex function Simplexmeth
    
    # Initialization of Simplex class
    def __init__(self,flowrate1,flowrate2,abs1,abs2,**kwargs):
        # flowrate1 = first flowrate guess in uL/min
        # flowrate2 = second flowrate guess in uL/min
        # abs1 = transmission at first flowrate
        # abs2 = transmission at secoond flowrate
        
        self.s = 1              # Step size (between 0 and 1, 0 = no step, 1 = Full step)
        self.nextflowrate = 0   # Variable to save next flowrate to
        self.n = 0              # Variable counting the number of iterations of simplex algorithm
        self.Memory = []        # List to save past flowrates and transmission measurements to
        
        # Determining in which direction to make first step
        absdif = abs1-abs2                                  # Calculating difference in transmission between the first two guesses
        if absdif < 0:                                      # Transmission of flowrate 2 is larger -> make step in the direction of flowrate 1
            self.Memory.append([self.n,flowrate2,abs2])     # Saving data
            # Making first step
            self.Simplexmeth(flowrate1,abs1)
        else:                                               # Transmission of flowrate 1 is larger -> make step in the direction of flowrate 2
            self.Memory.append([self.n,flowrate1,abs1])     # Saving data
            # Making first step
            self.Simplexmeth(flowrate2,abs2)

    # Function performing one step of simplex algorithm 
    def Simplexmeth(self,flowratenew,absnew,*kwargs):
        self.n = self.n + 1                                 # Increment iteration counter by one
        absold = self.Memory[self.n-1][2]                   # Read last transmission measurement from memory
        flowrateold = self.Memory[self.n-1][1]              # Read last flowrate from memory
        flowratedif = flowratenew - flowrateold             # Calculating the difference in flowrates (step size)
        # finding point with lower transmission
        absdif = absnew - absold
        # Making step
        if absdif > 0:
            # old flowrate was better
            self.nextflowrate = flowratenew - flowratedif
        elif absdif < 0:
            # New flowrate is better
            self.nextflowrate = flowratenew + flowratedif
        elif absdif == 0:
            # Same, so decrease step size (done in next foor loop)
            self.nextflowrate = flowratenew + flowratedif

        # If flowrate already tested decrease step size
        for i in range(0,len(self.Memory)):                 # Checking if flowrate already in memory
            if self.nextflowrate == self.Memory[i][1]:
                print('decrease step size')                 # If flowrate tested decrease step size
                self.s = self.s * 0.5                       # Decreasing step size by half
                if self.s < 0.1:                            # Checking if step size above minimum step size
                    self.s = 0.1 
                if absdif >= 0:
                    # Old flowrate was better
                    self.nextflowrate = flowratenew - flowratedif*self.s
                elif absdif < 0:
                    # New flowrate is better
                    self.nextflowrate = flowratenew + flowratedif*self.s
                # Round
                self.nextflowrate = round(self.nextflowrate,1) # Rounds new flowrate to first decimal place
            else:
                pass                                        # Normal step if flowrate is new

        # updating Memory
        self.Memory.append([self.n,flowratenew,absnew])

class WirteToFile():
    # Class for writing measurments to file

    # Initilaizing class
    def __init__(self,prefix,number):
        writy = False                       # Variable to check if filename already exists
        while writy == False:
            # Creating filename to write to
            self.filename = str(prefix) + str(number) + '.txt'
            # Checking if it exists
            if os.path.isfile(self.filename) == False:
                # File does not exist create new file
                writy = True
                self.file = open(self.filename,'w')
                self.file.write('Time;Flowrate;Absorbtion\n')
            else:
            # File already exists, changing filename
                print(" Suggested file already exists, incrementing filename by one number \n")
                number = number + 1

    def fileappend(self,time,flowrate,absorbtion):
        # Appending a line to the file, for raw data
        self.file.write('{0:.3f};{1:.1f};{2:.0f} \n'.format(time,flowrate,absorbtion))

    def fileappendfilter(self,absorbtion):
        # Appending a line to the file, for filtered data
        # Filtered data written with same header as normal data, time and flowrate can be found in corresponding raw data file
        for i in absorbtion:
            self.file.write('{0:.0f} \n'.format(i))

    def fileappendmemory(self,memory):
        # Writing the memory of the simplex algorithm to a text file
        for i in memory:
            self.file.write(str(i[0])+';'+str(i[1])+';'+str(i[2])+'\n')
    
    def fileclose(self,*kwargs):
        # Closing File
        self.file.close()
        

class MeasureAbs():
    # Class for measuring the transmission

    def __init__(self,*kwargs):
        # Initializing the class and controls

        # Defining pin connections
        self.SPICLK = 18     # For ADC
        self.SPIMISO = 23    # For ADC
        self.SPIMOSI = 24    # For ADC
        self.SPICS = 25      # For ADC
        self.PINPWM = 27     # Pin for PWM
        self.PINFAN = 17     # Pin for active cooling with fan
          
        # Setup SPI interface pins, defining in and output pins
        GPIO.setup(self.SPIMOSI, GPIO.OUT)
        GPIO.setup(self.SPIMISO, GPIO.IN)
        GPIO.setup(self.SPICLK, GPIO.OUT)
        GPIO.setup(self.SPICS, GPIO.OUT)

    def readadc(self,adcnum):
        # Function for reading value from ADC
        # Checking for correct channel
        if((adcnum > 7) or (adcnum < 0)):
            return -1
        GPIO.output(self.SPICS,True)

        GPIO.output(self.SPICLK, False)    # Start clock low
        GPIO.output(self.SPICS,False)        # bring CS low

        commandout = adcnum
        commandout |= 0x10              # Check, startbit + single end bit
        commandout <<= 3                # Only sending 5 bits
        for i in range(5):
            if(commandout & 0x80):
                GPIO.output(self.SPIMOSI,True)
            else:
                GPIO.output(self.SPIMOSI,False)
            commandout <<= 1
            GPIO.output(self.SPICLK, True)
            GPIO.output(self.SPICLK, False)

        adcout = 0
        # read one empty bit, one null bit and 10 ADC bits
        for i in range(15):
            GPIO.output(self.SPICLK, True)
            GPIO.output(self.SPICLK, False)
            adcout <<= 1
            if (GPIO.input(self.SPIMISO)):
                adcout |= 0x1

        GPIO.output(self.SPICS, True)

        adcout >>= 1            # First bit is null, so drop it
        return adcout           # Returning output of ADC

class LED():
    # Class for controlling the LED

    # Initializing class
    def __init__(self,LEDPIN):
        self.LEDPIN = LEDPIN                # Defining GPIO pin used for LED
        GPIO.setup(self.LEDPIN,GPIO.OUT)    # Defining LED pin as output
    # Function turning on the LED
    def ledon(self,*kwargs):
        GPIO.output(self.LEDPIN,True)
        print ('LED on')
    # Function turning of the LED
    def ledoff(self,*kwargs):
        GPIO.output(self.LEDPIN,False)
        print ('LED off')
        

class Pumpcomuni():
    # Class for communicating with the pump (Code for Chemyx Syring pumps Fusion 200)

    # Initializing class
    def __init__(self,portname,diameter):
        portname = '/dev/' + str(portname)  # USB port connected to pump, can be checked in terminal with: dmesg | grep tty
        print portname
        self.ser = serial.Serial(           # Starting serial communication
            port = portname,                # Defining port for communication
            baudrate = 9600,                # Baudrate of pump
            parity = serial.PARITY_NONE,    # Defining parity bits
            stopbits = serial.STOPBITS_ONE, # Defining serial stop bits
            timeout = 2                     # Defining read timeout
            )

        self.ser.write('echo off \r')       # Turning of echo function of pump
        # setting initial diameter for syringe
        print('Setting syringe diameter to {0} mm'.format(diameter))
        command = 'set diameter ' + str(diameter) + '\r'    # Making string for command
        self.ser.write(command)                             # Sending command to pump
        x = self.ser.readline()                             # Reading reply of pump
        print(x)                                            # Printing pump reply to console

        # Setting units of pump to ul/min
        print('Setting pump units to ul/min')
        self.ser.write('set units 2 \r')                    # Sending command to pump
        x = self.ser.readline()                             # Reading reply of pump
        print(x)                                            # Printing pump reply to console


    def pumprate(self,rate):
        # Setting pumping rate, also starts the pump, rate = desired pumping rate in uL/min
        print('Changing pumping rate')
        command = 'set rate ' + str(rate) + ' \r'           # Making command to change pump rate
        self.ser.write(command)                             # Sending command to pump
        x = self.ser.readline()                             # Reading reply of pump
        print(x)                                            # Printing pump reply to console

    def stopump(self,*kwargs):
        # Stopping pump
        print('Stopping pump')                              # Sending command to pump
        self.ser.write('stop \r')                           # Reading reply of pump
        x = self.ser.readline()                             # Printing pump reply to console
        print(x)

    def startpump(self,*kwargs):
        # Stopping pump
        print('Starting pump')
        self.ser.write('pause \r')                          # Sending command to pump, making sure it is in known (pause) state
        self.ser.write('start \r')                          # Sending command to pump to start
        self.ser.write('run \r')                            # Sending command to pump to run
        x = self.ser.readline()                             # Reading reply of pump
        print(x)                                            # Printing pump reply to console

    def helppump(self,*kwargs):
        # Printing help menu of pump
        print('Pump Commands')
        self.ser.write('help \r')                           # Sending command to pump to return help proms
        i = 1
        while i < 100:
            x = self.ser.readline()                         # Reading reply of pump
            print(x)                                        # Printing pump reply to console
        

# Function for measurement of transmission at a given flowrate
def measureatflowrate(flowrate,tglobal):
    tlocal = 0.0                                                        # Defining local time
    Writy = WirteToFile('File',flowrate)                                # Establishing file to safe data to
    Pumpy.pumprate(flowrate)                                            # Setting flowrate in pump
    print('Measuring at flow rate of {0} ul/min'.format(flowrate))
    print('Waste')                                                      # Sending first part of sample to waste to avoid collecting data at intermediate conditions
    while tlocal<tint:                                                  # Collecting data for waste-period in intermediate file, not used for simplex algorithm
        # Measuring six times and taking the average and taking average
        meas1 = float(Absorb.readadc(1))
        meas2 = float(Absorb.readadc(1))
        meas3 = float(Absorb.readadc(1))
        meas4 = float(Absorb.readadc(1))
        meas5 = float(Absorb.readadc(1))
        meas6 = float(Absorb.readadc(1))
        meas = (meas1 + meas2 + meas3 + meas4 + meas5 + meas6)/6
        # Updating time
        time.sleep(twait)                                               # Waiting predefined time step
        tglobal = tglobal + twait                                       # Updating global time
        tlocal = tlocal + twait                                         # Updating local time
        # Writing collected data to files
        Writy.fileappend(tglobal,flowrate,meas)

    # Closing file and writing to new one for data used for simplex algorithm
    Writy.fileclose()                                                   # Closing old file
    Writy = WirteToFile('File-Stable',flowrate)                         # Creating new file
    data = []                                                           # List to write data to
    print("Collect Sample")
    while tlocal<2*tint and tlocal >= tint:                             # Measuring  transmission at stable condition
        # Measuring 6 times and taking average
        meas1 = float(Absorb.readadc(1))
        meas2 = float(Absorb.readadc(1))
        meas3 = float(Absorb.readadc(1))
        meas4 = float(Absorb.readadc(1))
        meas5 = float(Absorb.readadc(1))
        meas6 = float(Absorb.readadc(1))
        meas = (meas1 + meas2 + meas3 + meas4 + meas5 + meas6)/6
        data.append(meas)
        # Writing to files
        Writy.fileappend(tglobal,flowrate,meas)
        # Updating time
        time.sleep(twait)                                               # Waiting predefined time step
        tglobal = tglobal + twait                                       # Updating global time
        tlocal = tlocal + twait                                         # Updating local time

    # Processing data
    FilterInstance.Lowpass(data)                                        # Lowpass filter on collected data
    print('Absorbtion at {0} is {1}'.format(flowrate,FilterInstance.average))
    absorbtion = FilterInstance.average                                 # defining average transmission used for simplex algorithm
    Writy.fileclose()                                                   # Closing raw data file
    Writy = WirteToFile('File-Filtered-Data-stable',flowrate)           # Creating file for filtered data
    Writy.fileappendfilter(FilterInstance.filterdata)                   # Writing filtered data to file
    Writy.fileclose()                                                   # Closing file
    return absorbtion                                                   # Returning average transmission for simplex algorithm from this function
    
# Main Code

Led = LED(27)                                                           # Creating instance for LED 
Led.ledon()                                                             # Turning on led
Pumpy = Pumpcomuni('ttyUSB0',8.7)                                       # Establishing pump instance, specifying port and syringe diameter
Absorb = MeasureAbs()                                                   # creating an instance for the transmission measurement

tglobal = 0.0                                                           # Variable for global elapsed time in seconds
tint = 120.0                                                            # Time period for sampling and waiting when changing flowrate in seconds
twait = 0.05                                                            # Time period between measurements in seconds
cutofffreq = 10                                                         # Cut-off frequency for lowpass filter in Hz
threshhold = 7800                                                       # Transmission threshold for defining droplet data and oil phase
FilterInstance = Filter(1/twait,cutofffreq,threshhold)                  # Creating Instance for Lowpass Filter
flowlow = 0.0                                                           # Lower limit for flowrate, preventing negative flow rates
flowupp = 30.0                                                          # Upper limit for flowrate, preventing to high back pressure

found = False                                                           # Boolean variable to indicate if optimum flowrate is found
flowrate1 = 16                                                          # First guess for optimum flowrate in ul/min
flowrate2 = 12                                                          # Second guess for optimum flowrate in ul/min


# Measuring first flowrate
Pumpy.pumprate(flowrate1)                                               # Setting flowrate in pump
Pumpy.startpump()                                                       # Starting pump

# Measuring transmission at first two flowrates
abs1 = measureatflowrate(flowrate1,tglobal)
abs2 = measureatflowrate(flowrate2,tglobal)

SimplexInstance = Simplex(flowrate1,flowrate2,abs1,abs2)                # Creating instance for simplex algorithm class

# Start Simplex algorithm
while found == False:
    # Determining new flow rate and measuring transmission
    flowratenew = SimplexInstance.nextflowrate                          # Determining next flowrate guess
    absnew = measureatflowrate(flowratenew,tglobal)                     # Measuring transmission at new flowrate
    SimplexInstance.Simplexmeth(flowratenew,absnew)                     # Simplex as defined by class function 
    # Checking for convergence
    for i in range(0,len(SimplexInstance.Memory)):
            if SimplexInstance.nextflowrate == SimplexInstance.Memory[i][1]:    # If new flowrate already tested, then optimum flowrate found
                print "found"
                found = True
                break
    else:                                                               # Else keep iterating 
        print('Next flow rate to be measured {0} ul/min'.format(SimplexInstance.nextflowrate))
        found = False
    if SimplexInstance.nextflowrate <= flowlow or SimplexInstance.nextflowrate > flowupp:   # New flowrate out of set bounds, end algorithm
        print("Flow rate out of bounds")
        found = True
print("Steps taken: \n {0}".format(SimplexInstance.Memory))             # number of iterations, flowrates, and transmission measurements 

# Stopping pumps and resetting GPIO
Pumpy.stopump()
GPIO.cleanup()

# Writing memory to a file
Writy = WirteToFile('Memory',0)
Writy.fileappendmemory(SimplexInstance.Memory)
Writy.fileclose()
