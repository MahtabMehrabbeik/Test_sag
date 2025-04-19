from networkIRR import *

# an example adjacency matrix
Aij32=np.array([[0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1], 
				[1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1], 
				[1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1], 
				[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1], 
				[0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0],
				[1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1], 
				[1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1], 
				[1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1], 
				[1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1], 
				[1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0], 
				[1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0]])
				
# create a networkIRR object for the adjacency matrix
data32=NetworkIRR(Aij32)

# get clusters
print ("the orbits are:")
print (data32.get_orbits())

print ("the number of symmetries is",data32.get_automorphism_group().order())

# generate T matrix
tmat=data32.get_transformation_operator()

print ("the T matrix:")
#print T matrix
for row in tmat:
	printstring=""
	for x in row:
		wch=str(round(real(x),3))
		printstring=printstring+wch+" "*(int(7)-len(wch))

	printstring=printstring+"\n"
	print (printstring)
	
# generate B matrix
tinv=la.inv(tmat)
b=np.dot(tmat,np.dot(Aij32,tinv))

print ("the B matrix:")
for row in b:
	printstring=""
	for x in row:
		wch=str(round(real(x),3))
		printstring=printstring+wch+" "*(int(7)-len(wch))
	printstring=printstring+"\n"
	print (printstring)
