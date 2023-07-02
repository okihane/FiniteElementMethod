import numpy as np
class Solve:
	def __init__(self,A,b):
		self.A = A
		self.b = b
	# LU解线性方程组
	def LU(self):
	    #生成值全位0的U矩阵，和单位矩阵L
	    A = self.A
	    L = np.eye(len(A))
	    U = np.zeros(np.shape(A))
	    for r in range(1, len(A)):  # 求U的第一行和L的第一列
	        U[0, r - 1] = A[0, r - 1]
	        L[r, 0] = A[r, 0] / A[0, 0]
	    U[0, -1] = A[0, -1]
	    for r in range(1, len(A)):  # 先求U再求L
	        for i in range(r, len(A)):
	            delta = 0
	            for k in range(0, r):  # 求∑(???∗???)
	                delta += L[r, k] * U[k, i]
	            U[r, i] = A[r, i] - delta
	            for i in range(r + 1, len(A)):  # 求L矩阵
	                theta = 0
	                for k in range(0, r):  # 求∑(?i?∗??r)
	                    theta += L[i, k] * U[k, r]
	                L[i, r] = (A[i, r] - theta) / U[r, r]
	    return L,U
	def my_LUsolve(self):
	    L, U = self.LU() # 得到L和U
	    #print("L={}".format(L))
	    #print("U={}".format(U))
	    # 求解线性方程LY=b
	    A = self.A
	    b = self.b
	    n = len(A)
	    y = np.zeros((n, 1))
	    b = np.array(b).reshape(n,1) # 把b列表格式变成向量格式
	    for i in range(len(A)):
	        t = 0
	        for j in range(i):
	            t += L[i][j]* y[j][0]
	        y[i][0] = b[i][0] - t
	    # print("y={}".format(y))
	    # 求解UX=Y
	    X = np.zeros((n, 1))
	    for i in range(len(A)-1,-1,-1):
	        t = 0
	        for j in range(i+1,len(A)):
	            t += U[i][j]*X[j][0]
	        t = y[i][0] - t
	        if t != 0 and U[i][i] == 0:
	            return 0
	        X[i] = t/U[i][i]
	    # print("X={}".format(X))
	    return X