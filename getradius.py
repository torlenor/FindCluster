#!/usr/bin/python
import subprocess

# cmdstr = "./findcluster.x -s 8 -t 4 -n 2 pollev.list -f FRACTION | grep Radius | awk \'{print $6;}\' | awk \'{print substr($0, 1, length() - 2)}\'"
cmdstr = "./findcluster.x -s 40 -t 18 -n 2 pollev40.list -f FRACTION | grep Radius | awk \'{print $6;}\' | awk \'{print substr($0, 1, length() - 2)}\'"

F1 = 0.0
F2 = 0.90

F = F1

maxcnt=10

givenradius=input('Please enter the desired radius: ')
deltagoal=input('Please enter the desired accuracy: ')

delta = float('inf')

runcmdstr = str.replace(cmdstr, "FRACTION", repr(F1));
cmd = subprocess.Popen(runcmdstr, shell=True, stdout=subprocess.PIPE)
first_line = cmd.stdout.readline()
radius1=float(str.strip(first_line));

runcmdstr = str.replace(cmdstr, "FRACTION", repr(F2));
cmd = subprocess.Popen(runcmdstr, shell=True, stdout=subprocess.PIPE)
first_line = cmd.stdout.readline()
radius2=float(str.strip(first_line));

if radius1<givenradius:
	print("ERROR")
	print("Radius1 = {0:.5f}".format(radius1))
	sys.exit(0)
if radius2>givenradius:
	print("ERROR")
	print("Radius2 = {0:.5f}".format(radius2))
	sys.exit(0)

cnt=0

bestdelta=float('inf')
bestradius=0
bestF=0

while (delta > deltagoal and cnt < maxcnt) :
	Fnew = 0.5*(F1+F2)

	runcmdstr = str.replace(cmdstr, "FRACTION", repr(Fnew));
	cmd = subprocess.Popen(runcmdstr, shell=True, stdout=subprocess.PIPE)
	first_line = cmd.stdout.readline()
	radiusnew=float(str.strip(first_line));

	if radiusnew<givenradius:
		F2=Fnew
		radius2=radiusnew
	else:
		F1=Fnew
		radius1=radiusnew

	delta = abs(radiusnew - givenradius)
	if delta<bestdelta:
		bestdelta=delta
		bestradius=radiusnew
		bestF=Fnew

	print("Radius1 = {0:.8f}".format(radius1))
	print("Radius2 = {0:.8f}".format(radius2))
	print("Radiusnew = {0:.8f}".format(radiusnew))
	print("F1 = {0:.8f}".format(F1))
	print("F2 = {0:.8f}".format(F2))
	print("Fnew = {0:.8f}".format(Fnew))
	print("Delta = {0:.8f}".format(delta))
	print("")

	cnt = cnt+1

print("------------------------------------------------------------------------")
print("Final results after {0:.0f} iterations:".format(cnt))
print("")
print("Radiusnew = {0:.8f}".format(bestradius))
print("Fnew = {0:.8f}".format(bestF))
print("Delta = {0:.8f}".format(bestdelta))
print("")
print("Have a nice day!")


