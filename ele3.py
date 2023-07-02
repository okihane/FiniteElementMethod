import numpy as np
class matrix:
	m_K = [[0,0],[0,0]]
	def __init__(self):
		self.m_K = [[0,0],[0,0]]
class Ele3:

	def __init__(self,element,eleshape,pos):
		self.m_BasicPoint = []
		for elei in element:
			self.m_BasicPoint.append(pos[elei-1])
		self.eleshape = eleshape
		self.i = element[0]
		self.j = element[1]
		self.k = element[2]
		self.m_Bi = 0
		self.m_Ci = 0
		self.m_Bj = 0
		self.m_Cj = 0
		self.m_Bm = 0
		self.m_Cm = 0
		self.m_Area = 0
		self.m_Modulus = 1
		self.m_Poisson = 0
		self.K = []

	def Initialize_BC(self):
		m_BasicPoint = self.m_BasicPoint
		#print(m_BasicPoint[1])
		self.m_Bi = m_BasicPoint[1][1] - m_BasicPoint[2][1]
		self.m_Ci = m_BasicPoint[2][0] - m_BasicPoint[1][0]
		self.m_Bj = m_BasicPoint[2][1] - m_BasicPoint[0][1]
		self.m_Cj = m_BasicPoint[0][0] - m_BasicPoint[2][0]
		self.m_Bm = m_BasicPoint[0][1] - m_BasicPoint[1][1]
		self.m_Cm = m_BasicPoint[1][0] - m_BasicPoint[0][0]
	def Initialize_Area(self):
		#double maxDetal_X, maxDetal_Y;
		#double a, b, c;
		m_Bi = self.m_Bi
		m_Ci = self.m_Ci
		m_Bj = self.m_Bj
		m_Cj = self.m_Cj
		m_Bm = self.m_Bm
		m_Cm = self.m_Cm
		#a = np.abs(m_Ci)
		#b = np.abs(m_Cj)
		#c = np.abs(m_Cm)
		maxDetal_X = max(np.abs(m_Ci),np.abs(m_Cj),np.abs(m_Cm))
		#a = np.abs(m_Bi)
		#b = np.abs(m_Bj)
		#c = np.abs(m_Bm)
		maxDetal_Y = max(np.abs(m_Bi),np.abs(m_Bj),np.abs(m_Bm))
		self.m_Area = maxDetal_X * maxDetal_Y - 0.5 * (np.abs(m_Bi) * np.abs(m_Ci) + np.abs(m_Bj) * np.abs(m_Cj) + np.abs(m_Bm) * np.abs(m_Cm))
		print(self.m_Area)

	def Update_TranMatrix(self):
		m_Bi = self.m_Bi
		m_Ci = self.m_Ci
		m_Bj = self.m_Bj
		m_Cj = self.m_Cj
		m_Bm = self.m_Bm
		m_Cm = self.m_Cm
		m_Poisson = self.m_Poisson
		m_SigmaTranMatrix = np.zeros((3,6))
		m_EpsiTranMatrix = np.zeros((3,6))
	
		index = self.m_Modulus * 0.5 / (self.m_Area * (1 - self.m_Poisson * self.m_Poisson))
		m_SigmaTranMatrix[0][0] = index * m_Bi
		m_SigmaTranMatrix[1][0] = index * m_Poisson * m_Bi
		m_SigmaTranMatrix[2][0] = index * 0.5 * (1 - m_Poisson) * m_Ci
	
		m_SigmaTranMatrix[0][1] = index * m_Poisson * m_Ci
		m_SigmaTranMatrix[1][1] = index * m_Ci
		m_SigmaTranMatrix[2][1] = index * 0.5 * (1 - m_Poisson) * m_Bi
	
		m_SigmaTranMatrix[0][2] = index * m_Bj
		m_SigmaTranMatrix[1][2] = index * m_Poisson * m_Bj
		m_SigmaTranMatrix[2][2] = index * 0.5 * (1 - m_Poisson) * m_Cj
	
		m_SigmaTranMatrix[0][3] = index * m_Poisson * m_Cj
		m_SigmaTranMatrix[1][3] = index * m_Cj
		m_SigmaTranMatrix[2][3] = index * 0.5 * (1 - m_Poisson) * m_Bj
	
		m_SigmaTranMatrix[0][4] = index * m_Bm
		m_SigmaTranMatrix[1][4] = index * m_Poisson * m_Bm
		m_SigmaTranMatrix[2][4] = index * 0.5 * (1 - m_Poisson) * m_Cm
	
		m_SigmaTranMatrix[0][5] = index * m_Poisson * m_Cm
		m_SigmaTranMatrix[1][5] = index * m_Cm
		m_SigmaTranMatrix[2][5] = index * 0.5 * (1 - m_Poisson) * m_Bm
	
	
		index = 0.5 / self.m_Area
		m_EpsiTranMatrix[0][0] = index * m_Bi
		m_EpsiTranMatrix[1][0] = 0.0
		m_EpsiTranMatrix[2][0] = index * m_Ci
	
		m_EpsiTranMatrix[0][1] = 0.0
		m_EpsiTranMatrix[1][1] = index * m_Ci
		m_EpsiTranMatrix[2][1] = index * m_Bi
	
		m_EpsiTranMatrix[0][2] = index * m_Bj
		m_EpsiTranMatrix[1][2] = 0.0
		m_EpsiTranMatrix[2][2] = index * m_Cj
	
		m_EpsiTranMatrix[0][3] = 0.0
		m_EpsiTranMatrix[1][3] = index * m_Cj
		m_EpsiTranMatrix[2][3] = index * m_Bj
	
		m_EpsiTranMatrix[0][4] = index * m_Bm
		m_EpsiTranMatrix[1][4] = 0.0
		m_EpsiTranMatrix[2][4] = index * m_Cm
	
		m_EpsiTranMatrix[0][5] = 0.0
		m_EpsiTranMatrix[1][5] = index * m_Cm
		m_EpsiTranMatrix[2][5] = index * m_Bm

	def Assemable_Element_KMatrix(self):
		m_Bi = self.m_Bi
		m_Ci = self.m_Ci
		m_Bj = self.m_Bj
		m_Cj = self.m_Cj
		m_Bm = self.m_Bm
		m_Cm = self.m_Cm
		m_Poisson = self.m_Poisson
		'''m_subKMatrix = []
		for i in range(self.eleshape):
			tmatrix = []
			for j in range(self.eleshape):
				tmatrix.append(matrix)
			m_subKMatrix.append(tmatrix)'''
		index1 = self.m_Modulus * 1.0 * 0.25 / (self.m_Area * (1 - m_Poisson * m_Poisson))
		K = np.zeros((2*self.eleshape,2*self.eleshape))
		#kii
		K[0][0] = index1 * (m_Bi * m_Bi + (1 - m_Poisson) * 0.5 * m_Ci * m_Ci)
		K[0][1] = index1 * (m_Poisson * m_Bi * m_Ci + (1 - m_Poisson) * 0.5 * m_Bi * m_Ci)
		K[1][0] = index1 * (m_Poisson * m_Bi * m_Ci + (1 - m_Poisson) * 0.5 * m_Bi * m_Ci)
		K[1][1] = index1 * ((1 - m_Poisson) * 0.5 * m_Bi * m_Bi + m_Ci * m_Ci)

		#kij
		K[0][2] = index1 * (m_Bi * m_Bj + (1 - m_Poisson) * 0.5 * m_Ci * m_Cj)
		K[0][3] = index1 * (m_Poisson * m_Bi * m_Cj + (1 - m_Poisson) * 0.5 * m_Bj * m_Ci)
		K[1][2] = index1 * (m_Poisson * m_Bj * m_Ci + (1 - m_Poisson) * 0.5 * m_Bi * m_Cj)
		K[1][4] = index1 * ((1 - m_Poisson) * 0.5 * m_Bi * m_Bj + m_Ci * m_Cj)
	
		#kim
		K[0][4] = index1 * (m_Bi * m_Bm + (1 - m_Poisson) * 0.5 * m_Ci * m_Cm)
		K[0][5] = index1 * (m_Poisson * m_Bi * m_Cm + (1 - m_Poisson) * 0.5 * m_Bm * m_Ci)
		K[1][4] = index1 * (m_Poisson * m_Bm * m_Ci + (1 - m_Poisson) * 0.5 * m_Bi * m_Cm)
		K[1][5] = index1 * ((1 - m_Poisson) * 0.5 * m_Bi * m_Bm + m_Ci * m_Cm)
	
		#kji
		K[2][0] = index1 * (m_Bi * m_Bj + (1 - m_Poisson) * 0.5 * m_Ci * m_Cj)
		K[3][0] = index1 * (m_Poisson * m_Bi * m_Cj + (1 - m_Poisson) * 0.5 * m_Bj * m_Ci)
		K[2][1] = index1 * (m_Poisson * m_Bj * m_Ci + (1 - m_Poisson) * 0.5 * m_Bi * m_Cj)
		K[3][1] = index1 * ((1 - m_Poisson) * 0.5 * m_Bi * m_Bj + m_Ci * m_Cj)
	
		#kjj
		K[2][2] = index1 * (m_Bj * m_Bj + (1 - m_Poisson) * 0.5 * m_Cj * m_Cj)
		K[3][3] = index1 * (m_Poisson * m_Bj * m_Cj + (1 - m_Poisson) * 0.5 * m_Bj * m_Cj)
		K[3][2] = index1 * (m_Poisson * m_Bj * m_Cj + (1 - m_Poisson) * 0.5 * m_Bj * m_Cj)
		K[3][3] = index1 * ((1 - m_Poisson) * 0.5 * m_Bj * m_Bj + m_Cj * m_Cj)
	
		#kjm
		K[2][4] = index1 * (m_Bj * m_Bm + (1 - m_Poisson) * 0.5 * m_Cj * m_Cm)
		K[2][5] = index1 * (m_Poisson * m_Bj * m_Cm + (1 - m_Poisson) * 0.5 * m_Bm * m_Cj)
		K[3][4] = index1 * (m_Poisson * m_Bm * m_Cj + (1 - m_Poisson) * 0.5 * m_Bj * m_Cm)
		K[3][5] = index1 * ((1 - m_Poisson) * 0.5 * m_Bj * m_Bm + m_Cj * m_Cm)
	
		#kmi
		K[4][0] = index1 * (m_Bi * m_Bm + (1 - m_Poisson) * 0.5 * m_Ci * m_Cm)
		K[5][0] = index1 * (m_Poisson * m_Bi * m_Cm + (1 - m_Poisson) * 0.5 * m_Bm * m_Ci)
		K[4][1] = index1 * (m_Poisson * m_Bm * m_Ci + (1 - m_Poisson) * 0.5 * m_Bi * m_Cm)
		K[5][1] = index1 * ((1 - m_Poisson) * 0.5 * m_Bi * m_Bm + m_Ci * m_Cm)
	
		#kmj
		K[4][2] = index1 * (m_Bj * m_Bm + (1 - m_Poisson) * 0.5 * m_Cj * m_Cm)
		K[5][2] = index1 * (m_Poisson * m_Bj * m_Cm + (1 - m_Poisson) * 0.5 * m_Bm * m_Cj)
		K[4][3] = index1 * (m_Poisson * m_Bm * m_Cj + (1 - m_Poisson) * 0.5 * m_Bj * m_Cm)
		K[5][3] = index1 * ((1 - m_Poisson) * 0.5 * m_Bj * m_Bm + m_Cj * m_Cm)
	
		#kmm
		K[4][4] = index1 * (m_Bm * m_Bm + (1 - m_Poisson) * 0.5 * m_Cm * m_Cm)
		K[5][4] = index1 * (m_Poisson * m_Bm * m_Cm + (1 - m_Poisson) * 0.5 * m_Bm * m_Cm)
		K[4][5] = index1 * (m_Poisson * m_Bm * m_Cm + (1 - m_Poisson) * 0.5 * m_Bm * m_Cm)
		K[5][5] = index1 * ((1 - m_Poisson) * 0.5 * m_Bm * m_Bm + m_Cm * m_Cm)
		'''k = []
		for i in range(self.eleshape):
			tempk = []
			for j in range(self.eleshape):
				tempk.append(m_subKMatrix[i][j].m_K)
			k.append(tempk)'''
		self.K = K
	def func_k(self):
		self.Initialize_BC()
		self.Initialize_Area()
		self.Update_TranMatrix()
		self.Assemable_Element_KMatrix()

def matrix_assembly(kk, Elementset):# 总刚度矩阵组装
    for e in Elementset:
        i, j, k = e.i, e.j, e.k
        for m in range(0, 2):
            for n in range(0, 2):# 将16个2x2矩阵按照索引值加入总刚度矩阵中
                kk[2*(i - 1)+m,2*(i - 1)+n] += e.K[m,n]
                kk[2*(j - 1)+m,2*(i - 1)+n] += e.K[2+m,n]
                kk[2*(k - 1)+m,2*(i - 1)+n] += e.K[4+m,n]
                kk[2*(i - 1)+m,2*(j - 1)+n] += e.K[m,2+n]
                kk[2*(j - 1)+m,2*(j - 1)+n] += e.K[2+m,2+n]
                kk[2*(k - 1)+m,2*(j - 1)+n] += e.K[4+m,2+n]
                kk[2*(i - 1)+m,2*(k - 1)+n] += e.K[m,4+n]
                kk[2*(j - 1)+m,2*(k - 1)+n] += e.K[2+m,4+n]
                kk[2*(k - 1)+m,2*(k - 1)+n] += e.K[4+m,4+n]
    return kk
'''

	def Update_EpsiAndSigma(self):
			 int i, j
		for (i = 0; i < 3; i++):
			for (j = 0; j < 6; j++):
				m_Sigma[i] += m_SigmaTranMatrix[i][j] * m_subDispMatrix[j]
				m_Epsi[i] += m_EpsiTranMatrix[i][j] * m_subDispMatrix[j]

	def Update_Element(Point& P1, Point& P2, Point& P3):
		m_BasicPoint[0].in(P1)
		m_BasicPoint[1].in(P2)
		m_BasicPoint[2].in(P3)
		m_subDispMatrix[0] = m_BasicPoint[0].m_dDisp[0]
		m_subDispMatrix[1] = m_BasicPoint[0].m_dDisp[1]
		m_subDispMatrix[2] = m_BasicPoint[1].m_dDisp[0]
		m_subDispMatrix[3] = m_BasicPoint[1].m_dDisp[1]
		m_subDispMatrix[4] = m_BasicPoint[2].m_dDisp[0]
		m_subDispMatrix[5] = m_BasicPoint[2].m_dDisp[1]

		Update_EpsiAndSigma()
'''