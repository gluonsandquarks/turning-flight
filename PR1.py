import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as mt
from os import system
from tabulate import tabulate

headers = ["Variable", "Value"]

flag = 0
flag2 = 0

CLmax = 0
Nlim = 0
g = 9.81
pi = 3.141592


uPos = []
CL = []

uGraph = np.arange(0.7, 1.4, 0.01)

################################
################################
####### Type of aircraft #######
################################
################################


print("Choose the type of aircraft:")

while flag == 0:
	print("""[1] Propeller --- [2] Turbojet""")
	type = int(input())
	
	if (type == 1):
		print("You have chosen a propeller aircraft.")
		flag = 1
	elif (type == 2):
		print("You have chosen a turbojet aircraft.")
		flag = 2
	else:
		print("Please type only the digits [1] or [2].")


#########################################################
#########################################################
####### Aerodynamical and Structural Restrictions #######
#########################################################
#########################################################

print("Are there aerodynamic and structural restrictions? [y/n]")

while flag2 == 0:
	ans = input()
	
	if (ans == "y"):
		a = 0

		while a == 0:
			b = 0
			print("Enter the value of CLmax: ")
			CLmax = float(input())
			print("Enter the value of Nlim: ")
			Nlim = float(input())
			print("You have chosen a CLmax = %f and Nlim = %f. Are these values correct? [y/n]" % (CLmax, Nlim))
			
			while b == 0:
				typo = input()
				if (typo == "y"):
					a = 1
					b = 1
				elif (typo == "n"):
					b = 1
				else:
					print("Please type only the characters [y] or [n].")
		
		flag2 = 1

	elif (ans == "n"):
		print("There will be no aerodynamic and/or structural restrictions for the analysis of the aircraft.")
		flag2 = 2
	else:
		print("Please type only the characters [y] or [n].")


########################################################
########################################################
####### Equations and calculations for each case #######
########################################################
########################################################


if (flag == 1 and flag2 == 1):
	# 1 1 = propeller with restrictions

	flag3 = 0
	while flag3 == 0:
		c = 0
		system('cls')
		
		########################
		### Preliminary data ###
		########################

		print("Enter the weight of the aircraft (Newtons):")
		W = float(input())
		print("Enter the surface area of the wings (m^2):")
		S = float(input())
		print("Enter the zero lift drag coefficient (CD0):")
		CD0 = float(input())
		print("Enter the factor k:")
		k = float(input())
		print("Enter the air density (kg/m^3):")
		rho = float(input())
		print("Enter the efficiency of the motor (eta):")
		eta = float(input())
		print("Enter the power of the motor (kW):")
		power = float(input())

		table = [["W", W], ["S", S], ["CD0", CD0], ["k", k], ["rho", rho], ["eta", eta], ["power", power]]
		print(tabulate(table, headers, tablefmt="grid"))
		print("Are these values correct? [y/n]")

		while c == 0:
			secTypo = str(input())
			if secTypo == 'y':
				c = 1
				flag3 = 1
			elif secTypo == 'n':
				c = 1
			else:
				print("Please type only the characters [y] or [n].")

	####################
	### Calculations ###
	####################
	
	vRef = mt.sqrt((2*W)/(rho*S)) * mt.pow((k/CD0), 1.0/4.0)
	Em = 1/(2*mt.sqrt(CD0*k))
	p = (eta*power*1000*Em)/(vRef*W)
	
	### MSTR ###
	rootsU = np.roots([1, 0, 0, p, -1])
	for i in range(len(rootsU)):
		if str(rootsU[i]).find("0j") != -1:
			if float(rootsU[i]) > 0:
				uMSTR = float(rootsU[i])			

	NMSTR = mt.sqrt((2*p*uMSTR) - mt.pow(uMSTR, 4))
	CLMSTR = mt.sqrt(CD0/k) * (NMSTR/mt.pow(uMSTR, 2))

	if (CLMSTR > CLmax or NMSTR > Nlim):
		CLMSTR = CLmax
		uMSTR = mt.pow((2*p)/(1 + ((k/CD0) * mt.pow(CLmax, 2))), 1.0/3.0)
		NMSTR = mt.sqrt((2*p*uMSTR) - mt.pow(uMSTR, 4))
		if (NMSTR > Nlim):
			NMSTR = Nlim
			rootsU = np.roots([1, 0, 0, -2*p, mt.pow(Nlim, 2)])

			for i in range(len(rootsU)):
				if str(rootsU[i]).find("0j") != -1:
					if float(rootsU[i]) > 0:
						uPos.append(float(rootsU[i]))

			for root in range(len(uPos)):
				CL.append(mt.sqrt(CD0/k) * (NMSTR/mt.pow(uPos[root], 2)))

			if (max(CL) <= CLmax):
				CLMSTR = max(CL)
				uMSTR = min(uPos)
			else:
				CLMSTR = min(CL)
				uMSTR = max(uPos)

			RMSTR = (mt.pow(vRef, 2) * mt.pow(uMSTR, 2)) / (g*mt.sqrt(mt.pow(NMSTR, 2) - 1))
			omegaMSTR = (g*mt.sqrt(mt.pow(NMSTR, 2) - 1)) / (vRef * uMSTR)
			EMSTR = Em*NMSTR*uMSTR / p
			muMSTR = mt.atan(mt.sqrt(mt.pow(NMSTR, 2) - 1)) * 180 / pi

		else:
			RMSTR = (mt.pow(vRef, 2) * mt.pow(uMSTR, 2)) / (g*mt.sqrt(mt.pow(NMSTR, 2) - 1))
			omegaMSTR = (g*mt.sqrt(mt.pow(NMSTR, 2) - 1)) / (vRef * uMSTR)
			EMSTR = Em*NMSTR*uMSTR / p
			muMSTR = mt.atan(mt.sqrt(mt.pow(NMSTR, 2) - 1)) * 180 / pi
	else:
		RMSTR = (mt.pow(vRef, 2) * mt.pow(uMSTR, 2)) / (g*mt.sqrt(mt.pow(NMSTR, 2) - 1))
		omegaMSTR = (g*mt.sqrt(mt.pow(NMSTR, 2) - 1)) / (vRef * uMSTR)
		EMSTR = Em*NMSTR*uMSTR / p
		muMSTR = mt.atan(mt.sqrt(mt.pow(NMSTR, 2) - 1)) * 180 / pi



	# Resetting the arrays
	uPos = []
	CL = []

	### SST ###
	uSST = 2/(3*p)
	NSST = mt.sqrt((2*uSST*p) - mt.pow(uSST, 4))
	CLSST = mt.sqrt(CD0/k) * (NSST/mt.pow(uSST, 2))

	if (CLSST > CLmax or NSST > Nlim):
		CLSST = CLmax
		uSST = mt.pow((2*p)/(1 + ((k/CD0) * mt.pow(CLmax, 2))), 1.0/3.0)
		NSST = mt.sqrt((2*uSST*p) - mt.pow(uSST, 4))
		if (NMSTR > Nlim):
			NSST = Nlim
			rootsU = np.roots([1, 0, 0, -2*p, mt.pow(Nlim, 2)])

			for i in range(len(rootsU)):
				if str(rootsU[i]).find("0j") != -1:
					if float(rootsU[i]) > 0:
						uPos.append(float(rootsU[i]))

			for root in range(len(uPos)):
				CL.append(mt.sqrt(CD0/k) * (NSST/mt.pow(uPos[root], 2)))

			if (max(CL) <= CLmax):
				CLSST = max(CL)
				uSST = min(uPos)
			else:
				CLSST = min(CL)
				uSST = max(uPos)

			RSST = (mt.pow(vRef, 2) * mt.pow(uSST, 2)) / (g * mt.sqrt(mt.pow(NSST, 2) - 1))
			omegaSST = (g * mt.sqrt(mt.pow(NSST, 2) - 1)) / (vRef*uSST)
			ESST = (uSST * Em * NSST) / p
			muSST = mt.atan(mt.sqrt(mt.pow(NSST, 2) - 1)) * 180 / pi

		else:
			RSST = (mt.pow(vRef, 2) * mt.pow(uSST, 2)) / (g * mt.sqrt(mt.pow(NSST, 2) - 1))
			omegaSST = (g * mt.sqrt(mt.pow(NSST, 2) - 1)) / (vRef*uSST)
			ESST = (uSST * Em * NSST) / p
			muSST = mt.atan(mt.sqrt(mt.pow(NSST, 2) - 1)) * 180 / pi

	else:
		RSST = (mt.pow(vRef, 2) * mt.pow(uSST, 2)) / (g * mt.sqrt(mt.pow(NSST, 2) - 1))
		omegaSST = (g * mt.sqrt(mt.pow(NSST, 2) - 1)) / (vRef*uSST)
		ESST = (uSST * Em * NSST) / p
		muSST = mt.atan(mt.sqrt(mt.pow(NSST, 2) - 1)) * 180 / pi

	# Resetting the arrays
	uPos = []
	CL = []

	### Nmax ###
	uNmax = mt.pow((p/2), 1.0/3.0)
	NNmax = mt.sqrt((2*p*uNmax) - mt.pow(uNmax, 4))
	CLNmax = mt.sqrt(CD0/k) * (NNmax/mt.pow(uNmax, 2))

	if (CLNmax > CLmax or NNmax > Nlim):
		CLNmax = CLmax
		uNmax = mt.pow((2*p)/(1 + ((k/CD0) * mt.pow(CLmax, 2))), 1.0/3.0)
		NNmax = mt.sqrt((2*uNmax*p) - mt.pow(uNmax, 4))
		if (NNmax > Nlim):
			NNmax = Nlim
			rootsU = np.roots([1, 0, 0, -2*p, mt.pow(Nlim, 2)])

			for i in range(len(rootsU)):
				if str(rootsU[i]).find("0j") != -1:
					if float(rootsU[i]) > 0:
						uPos.append(float(rootsU[i]))

			for root in range(len(uPos)):
				CL.append(mt.sqrt(CD0/k) * (NSST/mt.pow(uPos[root], 2)))

			if (max(CL) <= CLmax):
				CLNmax = max(CL)
				uNmax = min(uPos)
			else:
				CLNmax = min(CL)
				uNmax = max(uPos)

			RNmax = (mt.pow(vRef, 2) * mt.pow(uNmax, 2)) / (g * mt.sqrt(mt.pow(NNmax, 2) - 1))
			omegaNmax = (g * mt.sqrt(mt.pow(NNmax, 2) - 1)) / (vRef*uNmax)
			ENmax = Em * (mt.sqrt(3)/2)
			muNmax = mt.atan(mt.sqrt(mt.pow(NNmax, 2) - 1)) * 180 / pi

		else:
			RNmax = (mt.pow(vRef, 2) * mt.pow(uNmax, 2)) / (g * mt.sqrt(mt.pow(NNmax, 2) - 1))
			omegaNmax = (g * mt.sqrt(mt.pow(NNmax, 2) - 1)) / (vRef*uNmax)
			ENmax = Em * (mt.sqrt(3)/2)
			muNmax = mt.atan(mt.sqrt(mt.pow(NNmax, 2) - 1)) * 180 / pi

	else:
		RNmax = (mt.pow(vRef, 2) * mt.pow(uNmax, 2)) / (g * mt.sqrt(mt.pow(NNmax, 2) - 1))
		omegaNmax = (g * mt.sqrt(mt.pow(NNmax, 2) - 1)) / (vRef*uNmax)
		ENmax = Em * (mt.sqrt(3)/2)
		muNmax = mt.atan(mt.sqrt(mt.pow(NNmax, 2) - 1)) * 180 / pi


	#######################
	## Arrays for graphs ##
	#######################

	CLGraph = []
	omegaGraph = []
	RGraph = []
	NGraph = []

	for val in range(len(uGraph)):
		NGraph.append(mt.sqrt((2*p*uGraph[val]) - mt.pow(uGraph[val], 4)))
		CLGraph.append(mt.sqrt(CD0/k) * (mt.sqrt((2*p*uGraph[val]) - mt.pow(uGraph[val], 4))/mt.pow(uGraph[val], 2)))
		RGraph.append(((mt.pow(vRef, 2)*mt.pow(uGraph[val], 2))/ (g * (mt.sqrt( mt.pow(mt.sqrt((2*p*uGraph[val]) - mt.pow(uGraph[val], 4)), 2) - 1))) ))
		omegaGraph.append(( (g * (mt.sqrt( mt.pow(mt.sqrt((2*p*uGraph[val]) - mt.pow(uGraph[val], 4)), 2) - 1))) / (uGraph[val] * vRef)))


elif (flag == 1 and flag2 == 2):
	# 1 2 propeller without restrictions
	
	flag3 = 0
	while flag3 == 0:
		c = 0
		system('cls')
		
		########################
		### Preliminary data ###
		########################

		print("Enter the weight of the aircraft (Newtons):")
		W = float(input())
		print("Enter the surface area of the wings (m^2):")
		S = float(input())
		print("Enter the zero lift drag coefficient (CD0):")
		CD0 = float(input())
		print("Enter the factor k:")
		k = float(input())
		print("Enter the air density (kg/m^3):")
		rho = float(input())
		print("Enter the efficiency of the motor (eta):")
		eta = float(input())
		print("Enter the power of the motor (kW):")
		power = float(input())

		table = [["W", W], ["S", S], ["CD0", CD0], ["k", k], ["rho", rho], ["eta", eta], ["power", power]]
		print(tabulate(table, headers, tablefmt="grid"))
		print("Are these values correct? [y/n]")

		while c == 0:
			secTypo = str(input())
			if secTypo == 'y':
				c = 1
				flag3 = 1
			elif secTypo == 'n':
				c = 1
			else:
				print("Please type only the characters [y] or [n].")


	####################
	### Calculations ###
	####################
	
	vRef = mt.sqrt((2*W)/(rho*S)) * mt.pow((k/CD0), 1.0/4.0)
	Em = 1/(2*mt.sqrt(CD0*k))
	p = (eta*power*1000*Em)/(vRef*W)
	
	### MSTR ###
	rootsU = np.roots([1, 0, 0, p, -1])
	for i in range(len(rootsU)):
		if str(rootsU[i]).find("0j") != -1:
			if float(rootsU[i]) > 0:
				uMSTR = float(rootsU[i])			

	NMSTR = mt.sqrt((2*p*uMSTR) - mt.pow(uMSTR, 4))
	CLMSTR = mt.sqrt(CD0/k) * (NMSTR/mt.pow(uMSTR, 2))
	RMSTR = (mt.pow(vRef, 2) * mt.pow(uMSTR, 2)) / (g*mt.sqrt(mt.pow(NMSTR, 2) - 1))
	omegaMSTR = (g*mt.sqrt(mt.pow(NMSTR, 2) - 1)) / (vRef * uMSTR)
	EMSTR = Em*NMSTR*uMSTR / p
	muMSTR = mt.atan(mt.sqrt(mt.pow(NMSTR, 2) - 1)) * 180 / pi


	### SST ###
	uSST = 2/(3*p)
	NSST = mt.sqrt((2*uSST*p) - mt.pow(uSST, 4))
	CLSST = mt.sqrt(CD0/k) * (NSST/mt.pow(uSST, 2))
	RSST = (mt.pow(vRef, 2) * mt.pow(uSST, 2)) / (g * mt.sqrt(mt.pow(NSST, 2) - 1))
	omegaSST = (g * mt.sqrt(mt.pow(NSST, 2) - 1)) / (vRef*uSST)
	ESST = (uSST * Em * NSST) / p
	muSST = mt.atan(mt.sqrt(mt.pow(NSST, 2) - 1)) * 180 / pi

	### Nmax ###
	uNmax = mt.pow((p/2), 1.0/3.0)
	NNmax = mt.sqrt((2*p*uNmax) - mt.pow(uNmax, 4))
	CLNmax = mt.sqrt(CD0/k) * (NNmax/mt.pow(uNmax, 2))
	RNmax = (mt.pow(vRef, 2) * mt.pow(uNmax, 2)) / (g * mt.sqrt(mt.pow(NNmax, 2) - 1))
	omegaNmax = (g * mt.sqrt(mt.pow(NNmax, 2) - 1)) / (vRef*uNmax)
	ENmax = Em * (mt.sqrt(3)/2)
	muNmax = mt.atan(mt.sqrt(mt.pow(NNmax, 2) - 1)) * 180 / pi

	#######################
	## Arrays for graphs ##
	#######################

	CLGraph = []
	omegaGraph = []
	RGraph = []
	NGraph = []

	for val in range(len(uGraph)):
		NGraph.append(mt.sqrt((2*p*uGraph[val]) - mt.pow(uGraph[val], 4)))
		CLGraph.append(mt.sqrt(CD0/k) * (mt.sqrt((2*p*uGraph[val]) - mt.pow(uGraph[val], 4))/mt.pow(uGraph[val], 2)))
		RGraph.append(((mt.pow(vRef, 2)*mt.pow(uGraph[val], 2))/ (g * (mt.sqrt( mt.pow(mt.sqrt((2*p*uGraph[val]) - mt.pow(uGraph[val], 4)), 2) - 1))) ))
		omegaGraph.append(( (g * (mt.sqrt( mt.pow(mt.sqrt((2*p*uGraph[val]) - mt.pow(uGraph[val], 4)), 2) - 1))) / (uGraph[val] * vRef)))

elif (flag == 2 and flag2 == 1):
	# 2 1 turbojet with restrictions

	flag3 = 0
	while flag3 == 0:
		c = 0
		system('cls')
		
		########################
		### Preliminary data ###
		########################

		print("Enter the weight of the aircraft (Newtons):")
		W = float(input())
		print("Enter the surface area of the wings (m^2):")
		S = float(input())
		print("Enter the zero lift drag coefficient (CD0):")
		CD0 = float(input())
		print("Enter the factor k:")
		k = float(input())
		print("Enter the air density (kg/m^3):")
		rho = float(input())
		print("Enter the thrust of the motor (Newtons):")
		thrust = float(input())

		table = [["W", W], ["S", S], ["CD0", CD0], ["k", k], ["rho", rho], ["thrust", thrust]]
		print(tabulate(table, headers, tablefmt="grid"))
		print("Are these values correct? [y/n]")

		while c == 0:
			secTypo = str(input())
			if secTypo == 'y':
				c = 1
				flag3 = 1
			elif secTypo == 'n':
				c = 1
			else:
				print("Please type only the characters [y] or [n].")


	####################
	### Calculations ###
	####################

	vRef = mt.sqrt((2*W)/(rho*S)) * mt.pow((k/CD0), 1.0/4.0)
	Em = 1/(2*mt.sqrt(CD0*k))
	z = (thrust*Em)/W

	# Resetting the arrays
	uPos = []
	CL = []

	### MSTR ###
	uMSTR = 1		
	NMSTR = mt.sqrt((2*z) - 1)
	CLMSTR = mt.sqrt((CD0/k) * ((2*z) - 1))

	if (CLMSTR > CLmax or NMSTR > Nlim):
		CLMSTR = CLmax
		uMSTR = mt.sqrt((2*z)/ ( 1 + ( (k * mt.pow(CLmax, 2)) /CD0 ) ) )
		NMSTR = mt.sqrt( (2 * z * mt.pow(uMSTR, 2)) - mt.pow(uMSTR, 4))
		if (NMSTR > Nlim):
			NMSTR = Nlim
			rootsU = np.roots([1, 0, (-1*(2*z)), 0, mt.pow(Nlim, 2)])

			for i in range(len(rootsU)):
				if str(rootsU[i]).find("0j") != -1:
					if float(rootsU[i]) > 0:
						uPos.append(float(rootsU[i]))

			for root in range(len(uPos)):
				CL.append(mt.sqrt((mt.pow(Nlim, 2) * CD0)/ (k * mt.pow(uPos[root], 4))))

			if (max(CL) <= CLmax):
				CLMSTR = max(CL)
				uMSTR = min(uPos)
			else:
				CLMSTR = min(CL)
				uMSTR = max(uPos)

			omegaMSTR = (g/vRef) * mt.sqrt(2 * (z - 1))
			RMSTR = (mt.pow(vRef, 2) * mt.pow(uMSTR, 2)) / (g * mt.sqrt(2*(z - 1)))
			EMSTR = (Em / z) * mt.sqrt((2*z) - 1)
			muMSTR = mt.acos(1/(mt.sqrt((2*z) - 1))) * (180 / pi)


		else:
			omegaMSTR = (g/vRef) * mt.sqrt(2 * (z - 1))
			RMSTR = (mt.pow(vRef, 2) * mt.pow(uMSTR, 2)) / (g * mt.sqrt(2*(z - 1)))
			EMSTR = (Em / z) * mt.sqrt((2*z) - 1)
			muMSTR = mt.acos(1/(mt.sqrt((2*z) - 1))) * (180 / pi)

	else:
		omegaMSTR = (g/vRef) * mt.sqrt(2 * (z - 1))
		RMSTR = (mt.pow(vRef, 2) * mt.pow(uMSTR, 2)) / (g * mt.sqrt(2*(z - 1)))
		EMSTR = (Em / z) * mt.sqrt((2*z) - 1)
		muMSTR = mt.acos(1/(mt.sqrt((2*z) - 1))) * (180 / pi)


	# Resetting the arrays
	uPos = []
	CL = []


	### SST ###
	uSST = 1 / mt.sqrt(z)
	NSST = (1/z) * mt.sqrt((2*z*z) - 1)
	CLSST = mt.sqrt((CD0/k) * ((2*z*z) - 1))

	if (CLSST > CLmax or NSST > Nlim):
		CLSST = CLmax
		uSST = mt.sqrt((2*z)/ ( 1 + ( (k * mt.pow(CLmax, 2)) /CD0 ) ) )
		NSST = mt.sqrt( ((2 * z * z) - 1) / z)
		if (NSST > Nlim):
			NSST = Nlim
			rootsU = np.roots([1, 0, (-1*(2*z)), 0, mt.pow(Nlim, 2)])

			for i in range(len(rootsU)):
				if str(rootsU[i]).find("0j") != -1:
					if float(rootsU[i]) > 0:
						uPos.append(float(rootsU[i]))

			for root in range(len(uPos)):
				CL.append(mt.sqrt((mt.pow(Nlim, 2) * CD0)/ (k * mt.pow(uPos[root], 4))))

			if (max(CL) <= CLmax):
				CLSST = max(CL)
				uSST = min(uPos)
			else:
				CLSST = min(CL)
				uSST = max(uPos)

			omegaSST = mt.sqrt(((z*z) - 1) / z) * (g/vRef)
			RSST = (mt.pow(vRef, 2)/g) * (1 / mt.sqrt((z*z) - 1))
			ESST = (Em / (z*z)) * mt.sqrt((2*z*z) - 1)
			muSST = mt.acos(z/(mt.sqrt((2*z*z) - 1))) * (180 / pi)


		else:
			omegaSST = mt.sqrt(((z*z) - 1) / z) * (g/vRef)
			RSST = (mt.pow(vRef, 2)/g) * (1 / mt.sqrt((z*z) - 1))
			ESST = (Em / (z*z)) * mt.sqrt((2*z*z) - 1)
			muSST = mt.acos(z/(mt.sqrt((2*z*z) - 1))) * (180 / pi)

	else:
		omegaSST = mt.sqrt(((z*z) - 1) / z) * (g/vRef)
		RSST = (mt.pow(vRef, 2)/g) * (1 / mt.sqrt((z*z) - 1))
		ESST = (Em / (z*z)) * mt.sqrt((2*z*z) - 1)
		muSST = mt.acos(z/(mt.sqrt((2*z*z) - 1))) * (180 / pi)

	# Resetting the arrays
	uPos = []
	CL = []

	### NMax ###
	uNmax = mt.sqrt(z)
	NNmax = z
	CLNmax = mt.sqrt(CD0/k)

	if (CLNmax > CLmax or NNmax > Nlim):
		CLNmax = CLmax
		uNmax = mt.sqrt((2*z)/ ( 1 + ( (k * mt.pow(CLmax, 2)) /CD0 ) ) )
		NNmax = mt.sqrt( (2 * z * mt.pow(uNmax, 2)) - mt.pow(uNmax, 4))
		if (NNmax > Nlim):
			NNmax = Nlim
			rootsU = np.roots([1, 0, (-1*(2*z)), 0, mt.pow(Nlim, 2)])

			for i in range(len(rootsU)):
				if str(rootsU[i]).find("0j") != -1:
					if float(rootsU[i]) > 0:
						uPos.append(float(rootsU[i]))

			for root in range(len(uPos)):
				CL.append(mt.sqrt((mt.pow(Nlim, 2) * CD0)/ (k * mt.pow(uPos[root], 4))))

			if (max(CL) <= CLmax):
				CLNmax = max(CL)
				uNmax = min(uPos)
			else:
				CLNmax = min(CL)
				uNmax = max(uPos)

			omegaNmax = mt.sqrt(((z*z) - 1) / z) * (g/vRef)
			RNmax = (mt.pow(vRef, 2) * z)/(g * mt.sqrt((z*z) - 1))
			ENmax = Em
			muNmax = mt.acos(1/z) * (180 / pi)

		else:
			omegaNmax = mt.sqrt(((z*z) - 1) / z) * (g/vRef)
			RNmax = (mt.pow(vRef, 2) * z)/(g * mt.sqrt((z*z) - 1))
			ENmax = Em
			muNmax = mt.acos(1/z) * (180 / pi)

	else:
		omegaNmax = mt.sqrt(((z*z) - 1) / z) * (g/vRef)
		RNmax = (mt.pow(vRef, 2) * z)/(g * mt.sqrt((z*z) - 1))
		ENmax = Em
		muNmax = mt.acos(1/z) * (180 / pi)


	#######################
	## Arrays for graphs ##
	#######################

	CLGraph = []
	omegaGraph = []
	RGraph = []
	NGraph = []

	for val in range(len(uGraph)):
		NGraph.append(mt.sqrt((2*z*mt.pow(uGraph[val], 2)) - mt.pow(uGraph[val], 4)))
		CLGraph.append((mt.sqrt(CD0/k) * ((mt.sqrt((2*z*mt.pow(uGraph[val], 2)) - mt.pow(uGraph[val], 4)))/mt.pow(uGraph[val], 2))))
		RGraph.append(((mt.pow(vRef, 2)*mt.pow(uGraph[val], 2))/ (g * (mt.sqrt( mt.pow(mt.sqrt((2*z*mt.pow(uGraph[val], 2)) - mt.pow(uGraph[val], 4)), 2) - 1))) ))
		omegaGraph.append(( (g * (mt.sqrt( mt.pow(mt.sqrt((2*z*mt.pow(uGraph[val], 2)) - mt.pow(uGraph[val], 4)), 2) - 1))) / (uGraph[val] * vRef)))


elif (flag == 2 and flag2 == 2):
	# 2 2 turbojet without restrictions

	flag3 = 0
	while flag3 == 0:
		c = 0
		system('cls')
		########################
		### Preliminary data ###
		########################
		print("Enter the weight of the aircraft (Newtons):")
		W = float(input())
		print("Enter the surface area of the wings (m^2):")
		S = float(input())
		print("Enter the zero lift drag coefficient (CD0):")
		CD0 = float(input())
		print("Enter the factor k:")
		k = float(input())
		print("Enter the air density (kg/m^3):")
		rho = float(input())
		print("Enter the thrust of the motor (Newtons):")
		thrust = float(input())

		table = [["W", W], ["S", S], ["CD0", CD0], ["k", k], ["rho", rho], ["thrust", thrust]]
		print(tabulate(table, headers, tablefmt="grid"))
		print("Are these values correct? [y/n]")

		while c == 0:
			secTypo = str(input())
			if secTypo == 'y':
				c = 1
				flag3 = 1
			elif secTypo == 'n':
				c = 1
			else:
				print("Please type only the characters [y] or [n].")

	####################
	### Calculations ###
	####################

	vRef = mt.sqrt((2*W)/(rho*S)) * mt.pow((k/CD0), 1.0/4.0)
	Em = 1/(2*mt.sqrt(CD0*k))
	z = (thrust*Em)/W

	### MSTR ###
	uMSTR = 1
	NMSTR = mt.sqrt((2*z) - 1)
	CLMSTR = mt.sqrt((CD0/k) * ((2*z) - 1))
	omegaMSTR = (g/vRef) * mt.sqrt(2 * (z - 1))
	RMSTR = (mt.pow(vRef, 2) * mt.pow(uMSTR, 2)) / (g * mt.sqrt(2*(z - 1)))
	EMSTR = (Em / z) * mt.sqrt((2*z) - 1)
	muMSTR = mt.acos(1/(mt.sqrt((2*z) - 1))) * (180 / pi)

	### SST ###
	uSST = 1 / mt.sqrt(z)
	CLSST = mt.sqrt((CD0/k) * ((2*z*z) - 1))
	NSST = (1/z) * mt.sqrt((2*z*z) - 1)
	omegaSST = mt.sqrt(((z*z) - 1) / z) * (g/vRef)
	RSST = (mt.pow(vRef, 2)/g) * (1 / mt.sqrt((z*z) - 1))
	ESST = (Em / (z*z)) * mt.sqrt((2*z*z) - 1)
	muSST = mt.acos(z/(mt.sqrt((2*z*z) - 1))) * (180 / pi)

	### NMax ###
	uNmax = mt.sqrt(z)
	CLNmax = mt.sqrt(CD0/k)
	NNmax = z
	omegaNmax = mt.sqrt(((z*z) - 1) / z) * (g/vRef)
	RNmax = (mt.pow(vRef, 2) * z)/(g * mt.sqrt((z*z) - 1))
	ENmax = Em
	muNmax = mt.acos(1/z) * (180 / pi)

	#######################
	## Arrays for graphs ##
	#######################

	CLGraph = []
	omegaGraph = []
	RGraph = []
	NGraph = []

	for val in range(len(uGraph)):
		NGraph.append(mt.sqrt((2*z*mt.pow(uGraph[val], 2)) - mt.pow(uGraph[val], 4)))
		CLGraph.append((mt.sqrt(CD0/k) * ((mt.sqrt((2*z*mt.pow(uGraph[val], 2)) - mt.pow(uGraph[val], 4)))/mt.pow(uGraph[val], 2))))
		RGraph.append(((mt.pow(vRef, 2)*mt.pow(uGraph[val], 2))/ (g * (mt.sqrt( mt.pow(mt.sqrt((2*z*mt.pow(uGraph[val], 2)) - mt.pow(uGraph[val], 4)), 2) - 1))) ))
		omegaGraph.append(( (g * (mt.sqrt( mt.pow(mt.sqrt((2*z*mt.pow(uGraph[val], 2)) - mt.pow(uGraph[val], 4)), 2) - 1))) / (uGraph[val] * vRef)))


###############
#### TABLE ####
###############

x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)




fig, axs = plt.subplots(2, 2)


plt.figure(1)

axs[0, 0].plot(uMSTR, CLMSTR, marker='x', color='r', label='CL MSTR')
axs[0, 0].plot(uSST, CLSST, marker='x', color='b', label='CL SST')
axs[0, 0].plot(uNmax, CLNmax, marker='x', color='g', label='CL Nmax')
axs[0, 0].plot(uGraph, CLGraph)
axs[0, 0].axhline(y=CLmax, color='c', linestyle='-', label='CLmax')
axs[0, 0].set_title("Cl vs. u")
axs[1, 0].plot(uMSTR, omegaMSTR, marker='x', color='r', label='omega MSTR')
axs[1, 0].plot(uSST, omegaSST, marker='x', color='b', label='omega SST')
axs[1, 0].plot(uNmax, omegaNmax, marker='x', color='g', label='omega Nmax')
axs[1, 0].plot(uGraph, omegaGraph)
axs[1, 0].set_title("omega vs. u")
axs[0, 1].plot(uMSTR, RMSTR, marker='x', color='r', label='R MSTR')
axs[0, 1].plot(uSST, RSST, marker='x', color='b', label='R SST')
axs[0, 1].plot(uNmax, RNmax, marker='x', color='g', label='R Nmax')
axs[0, 1].plot(uGraph, RGraph)
axs[0, 1].set_title("R vs u")
axs[1, 1].plot(uMSTR, NMSTR, marker='x', color='r', label='N MSTR')
axs[1, 1].plot(uSST, NSST, marker='x', color='b', label='N SST')
axs[1, 1].plot(uNmax, NNmax, marker='x', color='g', label='N Nmax')
axs[1, 1].plot(uGraph, NGraph)
axs[1, 1].axhline(y=Nlim, color='c', linestyle='-', label='Nlim')
axs[1, 1].set_title("N vs u")
axs[0, 0].legend()
axs[1, 0].legend()
axs[0, 1].legend()
axs[1, 1].legend()

fig.tight_layout()




fig, ax = plt.subplots(dpi=100)

# Hide axes
plt.figure(2)
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')

tableVariables=np.array(["u","omega (rad/s)","R (m)","N","Cl","E","Mu (deg)"])

mstrData=np.around(np.array([uMSTR, omegaMSTR, RMSTR, NMSTR, CLMSTR, EMSTR, muMSTR]),6)
string_magnitudes=str(mstrData)

sstData=np.around(np.array([uSST, omegaSST, RSST, NSST, CLSST, ESST, muSST]),6)
string_magnitudes=str(sstData)

nmaxData=np.around(np.array([uNmax, omegaNmax, RNmax, NNmax, CLNmax, ENmax, muNmax]),6)
string_magnitudes=str(nmaxData)


df = pd.DataFrame({'MSTR' : mstrData.tolist(), 'SST' : sstData.tolist(), 'Nmax' : nmaxData.tolist()})
ax.table(cellText=df.values, rowLabels=tableVariables, colLabels=df.columns, cellLoc='center', loc='center')
fig.tight_layout()
####
 

plt.show()