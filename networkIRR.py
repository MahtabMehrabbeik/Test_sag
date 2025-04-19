import scipy as sp
import numpy.random as rn
import numpy.linalg as la
import numpy as np
import sage.all as sg

class NetworkIRR:
	'''
	Class NetworkIRR:
	
	NetworkIRR objects contain an adjacency matrix that represents a 
	network, and related data related to the IRR coordinate system
	associated with that network.
	
	Instance data:
	
	_adjmat: Adjacency matrix in node coordinates
	
	_nnodes: integer, number of nodes
	
	_group: Sage Group Object, the automorphism group of the network
	_permutation_matricies: List of permutation matricies 
	in node space representing the the automorphism group,
	
	_conjugacy_classes: Result of Sage call
	
	_conjugacy_classes_matrices: List of list of matrices representing
	conjugacy classes, in the same order as _conjugacy classes
	
	_character_table: Result of Sage call 
	
	_IRR_degeneracies: Number of times each IRR is present in the 
	permutation matrix representation. order is the same as the order of 
	IRRs in the character table
	
	_projection_operators: List of projection operators, in the order of 
	IRRs in the character table
	_transformation_operator
   
	'''
	
	def __init__(self, adjmat=None):
		'''
		Constructor for NetworkIRR objects. 
		
		Arguments:
		adjmat: adjacency matrix, 2-dimensional numpy array 
		
		'''
		self._reset_data()
		self._set_adjacency_matrix(adjmat)
		
	
	def _reset_data(self):
		''' 
		Resets the NetworkIRR object. 
		This is called by both _set_adjacency_matrix and 
		the constructor. This is a private method, 
		use it at your own peril.
		
		'''
		self._adjmat=None
		self._nnodes=0
		self._group=None
		self._orbits=None
		self._permutation_matrices=None
		self._conjugacy_classes=None
		self._conjugacy_classes_matrices=None
		self._character_table=None
		self._IRR_degeneracies=None
		self._projection_operators=None
		self._transformation_operator=None

	def _set_adjacency_matrix(self,adjmat):
		'''
		Sets the adjacency matrix. 
		
		Note: calling this method clears all data 
		(i.e. group elements, projection operators)
		This is a private method, use it at your own peril.
		
		'''
		self._reset_data();
		self._adjmat=np.array(adjmat.copy())
		self._nnodes=len(self._adjmat)
		
		
	def get_adjacency_matrix(self):
		'''
		Returns a deep copy of the current adjacency matrix. 
		
		'''
		return self._adjmat.copy() 
		
	def get_automorphism_group(self):
		'''
		Returns a copy of the automorphism group of the network
		as a Sage group object
		
		'''
		if self._group==None:
			self._group,self._orbits = \
			(sg.Graph(sg.Matrix(self._adjmat))).automorphism_group(orbits=True)
		
		return sg.copy(self._group)
		
	def get_automorphism_group_matrices(self):
		'''
		Returns a list of _nnodesX_nnodes numpy matrices.
		
		These are permutation matrices in node space that represent 
		the group elements. Use caution in calling get_group_matrices.
		in a large network with many symetries, 
		the result could require a lot of memory.
		
		'''
		if self._group==None:
			self.get_automorphism_group()
		
		if self._permutation_matrices==None:
			self._permutation_matrices=[]
			
			for element in self._group.list():
				self._permutation_matrices\
				.append(np.array(element.matrix()))
				
		return list(self._permutation_matrices)
		
		
	def get_orbits(self):
		'''
		Returns the orbits of automorphism group as a list of lists

		'''
		print(self.get_adjacency_matrix())
		if self._orbits==None:
			self._group,self._orbits = \
			sg.Graph(sg.Matrix(self.get_adjacency_matrix())). \
			automorphism_group(orbits=True)
			
		return sg.copy(self._orbits)
	
	def get_character_table(self):
		'''
		Returns the character table of the group

		'''
		if self._character_table==None:
			self._character_table=self.get_automorphism_group().character_table()

		return sg.copy(self._character_table)
	
	def get_conjugacy_classes(self):
		'''
		Returns the conjugacy classes

		'''
		if self._conjugacy_classes==None:
			if self._group==None:
				self.get_automorphism_group()
				
			self._conjugacy_classes=self._group.conjugacy_classes()
			
		return sg.copy(self._conjugacy_classes)
		
		
		
	def get_conjugacy_classes_matrices(self):
		'''
		Returns the conjugacy classes, 
		as a list of lists of permutation matricies in node space

		'''
		if self._conjugacy_classes==None:
			self.get_conjugacy_classes()
		
		if self._conjugacy_classes_matrices==None:
			self._conjugacy_classes_matrices=[]
			
			for conjclass in self._conjugacy_classes:
				sublist=[]
				#This line makes no sense, but conjclass.list()
				#returns an error
				clist=sg.ConjugacyClass(self._group,\
				conjclass.representative()).list()
				
				for element in clist:
					sublist.append(np.array(element.matrix()))
					
				self._conjugacy_classes_matrices.append(sublist)
			
		return list(self._conjugacy_classes_matrices)

	def get_numIRRs(self):
		'''
		Returns the number of IRRs
		
		'''
		characters=self.get_character_table()
		numIRRs=len(characters[0])
		return numIRRs
		
	def get_IRR_degeneracies(self):
		'''
		Returns a list representing the number of times 
		each IRR is present in the permutation matrix representation.
		
		'''
		if self._IRR_degeneracies==None:
			characters=self.get_character_table()
			numIRRs=len(characters[0])
			group_order=self.get_automorphism_group().order()
			
			self._IRR_degeneracies=[]
			matricies=self.get_conjugacy_classes_matrices()
			for i in range(numIRRs): #Loop over IRRs
				total=0
				for j in range(numIRRs): #loop over classes
					total=total +float(len(matricies[j])/group_order)* np.conj(np.complex(characters[i][j]))* np.trace(matricies[j][0])
				
				self._IRR_degeneracies.append(round(np.real(total)))
				
		return list(self._IRR_degeneracies)
		
	def get_projection_operator(self, j):
		'''
		Returns a numpy array representing a matrix that
		projects into the subspace of an irreducible representaion.
		
		Arguments:
		j-index of the IRR (jth row of the character table). 
		
		'''
		degen=self.get_IRR_degeneracies()
		if self._projection_operators==None:
			self._projection_operators=[None]*len(degen)
		
		if self._projection_operators[j] == None:
			#the dimension of the IRR is the character of the identity.
			IRR_dimension=self.get_character_table()[j,0]
			group_order=self.get_automorphism_group().order()
			characters=self.get_character_table()[j]
			matricies=self.get_conjugacy_classes_matrices()
			result=np.zeros((self._nnodes,self._nnodes))
			
			for i in range(len(characters)):
				for mat in matricies[i]:
					result=result+mat*complex(characters[i])
					
			self._projection_operators[j]=result*np.float(IRR_dimension)/np.float(group_order)
		
		return self._projection_operators[j]

	def get_transformation_operator(self):
		'''
		Returns the tranformation operator into 
		the IRR coordinate system as a numpy array.

		'''
		epsilon=1e-12
		if self._transformation_operator==None:
			result=[]
			degens=self.get_IRR_degeneracies()
			total=0
			for j in range(len(degens)):
				IRR_dimension=self.get_character_table()[j,0]
				if degens[j]:
					P=self.get_projection_operator(j)
					U,W,VH = la.svd(P)
					
					
					#Uncomment this output if you want to see it.
					#print "Representation ", j, " with degeneracy ", degens[j]
					#print "W=",W
					#print "dimension=",IRR_dimension
					
					R1=int(IRR_dimension)*degens[j]
					R2=0
					total=total+R1 
					for i,w in enumerate(W):
						if np.abs(w-1.0)<epsilon:
							R2=R2+1
							result.append(VH[i])
						
					if R1!=R2:
						print ("Warning!")
						print ("Found, ",R2," singular vectors")
						print ("there should be", R1)
						raise Exception
						
					#print "Rank=",R1,R2
					#print "P=",np.real(P)
					
			#print "sum of dimension*degeneracy", total
			self._transformation_operator=np.array(result)
		return self._transformation_operator.copy()
	 
