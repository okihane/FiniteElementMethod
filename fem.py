#pylint:disable=E1101
#from ele4 import *
import ele3
import ele4
import solve
import numpy as np
class FEM():
	def __init__(self,filename):
		self.filename = filename
		self.npoints = 0
		self.positions = []
		self.elements = []
		self.eleshape = 0
		self.nconstraints = 0
		self.constraints = []
		self.fs = []
		self.K = []
	def model(self):
		lines = open(self.filename,encoding='utf-8').readlines()
		n = 3
		linesp = lines[n].split()
		self.npoints = int(linesp[0])
		dimpos = int(linesp[1])
		n += 1
		#self.positions = []
		for i in range(n,n+self.npoints):
			linesp = lines[i].split()
			if dimpos == 2:
				tposition = [float(linesp[1]),float(linesp[2])]
			elif dimpos == 3:
				tposition = [float(linesp[1]),float(linesp[2]),float(linesp[3])]
			self.positions.append(tposition)
		print(self.positions)
		#读取单元
		n += (self.npoints+3)
		linesp = lines[n].split()
		nelements = int(linesp[0])
		self.eleshape = int(linesp[1])
		dimele = int(linesp[2])
		n += 1
		self.elements = []
		for i in range(n,n+nelements):
			linesp = lines[i].split()
			if self.eleshape == 3:
				telement = [int(linesp[1]),int(linesp[2]),int(linesp[3])]
			elif self.eleshape == 4:
				telement = [int(linesp[1]),int(linesp[2]),int(linesp[3]),int(linesp[4])]
			self.elements.append(telement)
		print(self.elements)
		#读取约束
		n += (nelements+3)
		linesp = lines[n].split()
		nconstraints = int(linesp[0])
		self.nconstraints = nconstraints
		dimcons = int(linesp[1])
		n += 1
		#self.constraints = []
		for i in range(n,n+nconstraints):
			linesp = lines[i].split()
			if dimcons == 2:
				tconstraint = [int(linesp[0]),int(linesp[1]),int(linesp[2])]
			if dimcons == 3:
				tconstraint = [int(linesp[0]),int(linesp[1]),int(linesp[2]),int(linesp[3])]
			self.constraints.append(tconstraint)
		print(self.constraints)
		#读取荷载
		n += (nconstraints+3)
		linesp = lines[n].split()
		nfs = int(linesp[0])
		dimf = int(linesp[1])
		n += 1
		#self.fs = []
		for i in range(n,n+nfs):
			linesp = lines[i].split()
			if dimf == 2:
				tf = [int(linesp[0]),float(linesp[1]),float(linesp[2])]
			if dimf == 3:
				tf = [int(linesp[0]),float(linesp[1]),float(linesp[2]),float(linesp[3])]
			self.fs.append(tf)
		print(self.fs)
	def kkk(self):
		Elementset = []
		# 创建一个总体刚度矩阵的列表
		kk = np.zeros((self.npoints*2, self.npoints*2))
		if self.eleshape == 4:
			for i in range(len(self.elements)):
				#print(elements[i])
				tempK = ele4.Ele4(self.elements[i],self.positions)
				tempK.func_k()
				#tempK = K_element(self.elements[i],self.positions)
				Elementset.append(tempK)
				print("单元{}的刚度矩阵：{}".format(i+1,tempK.K))
			# 然后得到组合后的刚度矩阵
			ele4.matrix_assembly(kk, Elementset)
		elif self.eleshape == 3:
			for i in range(len(self.elements)):
				#print(elements[i])
				tempK = ele3.Ele3(self.elements[i],self.eleshape,self.positions)
				tempK.func_k()
				#tempK = K_element(self.elements[i],self.positions)
				Elementset.append(tempK)
				print("单元{}的刚度矩阵：{}".format(i+1,tempK.K))
			# 然后得到组合后的刚度矩阵
			ele3.matrix_assembly(kk, Elementset)
		self.K = kk
		print("组合后的总刚度矩阵：{}".format(kk))
		#print(len(kk))
	def solve(self):
		fs = self.fs
		fs2 = np.zeros((self.npoints,2))
		for fi in fs:
			for i in range(1,3):
			#for fii in fi:
				#fs2.append([fii])
				fs2[fi[0]-1][i-1] = fi[i]
		new_F = np.array(fs2).reshape(-1,1)
		constraints = self.constraints
		delindex = []
		for i in range(len(constraints)):
			coni = constraints[i]
			for j in range(1,len(coni)):
				if coni[j] == 1:
					delindex.append(2*(coni[0]-1)+j-1)
		print(delindex)
		new2_F = np.delete(new_F, delindex, axis=0)
		new_kk_row = np.delete(self.K, delindex, axis=0)
		new_kk = np.delete(new_kk_row, delindex, axis=1)# 在原总刚度下删除U=0所在的行和列，得到一个新的矩阵
		#mysolve = solve.Solve(new_kk, new2_F)
		'''
		fs2 = np.zeros((self.npoints-self.nconstraints,2))
		print(self.npoints-self.nconstraints)
		print(fs2)
		fs = self.fs
		#fs2 = []
		for fi in fs:
			for i in range(1,3):
			#for fii in fi:
				#fs2.append([fii])
				fs2[fi[0]-1][i-1] = fi[i]
		new2_F = np.array(fs2).reshape(-1,1)
		'''
		#new2_F1 = np.matrix([[0],[-1000],[-1000],[0],[0],[-1000],[1000],[0],[-1000],[0],[0],[1000],[1000],[0]])
		#print(new2_F1)
		#new2_F = np.matrix(fs2)
		print(new2_F)
		new2_U = solve.Solve(new_kk, new2_F).my_LUsolve()
		print(new2_U)
		npoints = self.npoints
		U1 = np.zeros((npoints*2, 1))
		nocons = np.arange(0,2*npoints)
		#print(nocons)
		nocons = np.delete(nocons,delindex)
		#print(nocons)
		for i in range(len(nocons)):
			U1[nocons[i]] = new2_U[i]
		print(U1)
		'''
		U2 = np.zeros((npoints*2, 1))
		U2[2, 0] = new2_U[0, 0]
		U2[3, 0] = new2_U[1, 0]
		for i in range(12):
		    U2[i+6, 0] = new2_U[i+2, 0]
		print(U2)
		'''
		print("各节点的位移U1为：")
		print(U1)
		print(np.shape(self.K))
		print(np.shape(U1))
		print("各节点的力F1为：")
		F1 = self.K*np.matrix(U1)
		print(F1)