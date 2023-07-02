import numpy as np
class Ele4:
	def __init__(self,element,pos):
		self.k = element[0]
		self.l = element[1]
		self.i = element[2]
		self.j = element[3]
		self.pos = pos
		self.K = []
		#self.K = func_k(self.i,self.j,self.k,self.l,self.pos)
	# 1个单元四个节点
	def shapefunction(self,r,s):#形函数
	    N1 = 1 / 4 * (1 - r) * (1 - s)
	    N2 = 1 / 4 * (1 + r) * (1 - s)
	    N3 = 1 / 4 * (1 + r) * (1 + s)
	    N4 = 1 / 4 * (1 - r) * (1 + s)
	    return N1,N2,N3,N4
	
	def diffNdr(self,r,s): # 求dNidr
	    dN1dr = 1 / 4 * (-1) * (1 - s)
	    dN2dr = 1 / 4 * (1) * (1 - s)
	    dN3dr = 1 / 4 * (1) * (1 + s)
	    dN4dr = 1 / 4 * (-1) * (1 + s)
	    dNdr = [dN1dr,dN2dr,dN3dr,dN4dr]
	    return dNdr
	
	
	def diffNds(self,r,s): # 求dNids
	    dN1ds = 1 / 4 * (1 - r) * (-1)
	    dN2ds = 1 / 4 * (1 + r) * (-1)
	    dN3ds = 1 / 4 * (1 + r) * (1)
	    dN4ds = 1 / 4 * (1 - r) * (1)
	    dNds = [dN1ds, dN2ds, dN3ds, dN4ds]
	    return dNds
	
	def jacobian(self,x,y,r,s): # 求J,Jinv,Jdet
	    dNdr = self.diffNdr(r,s)
	    dNds = self.diffNds(r,s)
	
	    dxdr = x[0]*dNdr[0]+x[1]*dNdr[1]+x[2]*dNdr[2]+x[3]*dNdr[3]
	    dxds = x[0]*dNds[0]+x[1]*dNds[1]+x[2]*dNds[2]+x[3]*dNds[3]
	    dydr = y[0]*dNdr[0]+y[1]*dNdr[1]+y[2]*dNdr[2]+y[3]*dNdr[3]
	    dyds = y[0]*dNds[0]+y[1]*dNds[1]+y[2]*dNds[2]+y[3]*dNds[3]
	
	    J = np.array([[dxdr,dxds],[dydr,dyds]])
	    Jdet = np.linalg.det(J)
	    # Jdet = J[0][0]*J[1][1]-J[0][1]*J[1][0]
	    Jinv = np.linalg.inv(J)
	    return Jinv,Jdet
	
	def Bmatrix(self,r,s,Jinv):# 求B
	    dNdr = self.diffNdr(r, s)
	    dNds = self.diffNds(r, s)
	
	    B1 = np.matrix([[1,0,0,0],[0,0,0,1],[0,1,1,0]])
	    B2 = np.zeros((4,4))
	    B2[0:2,0:2] = Jinv
	    B2[2:4,2:4] = Jinv
	    B3 = np.zeros((4,8))
	    B3[0,0] = dNdr[0]
	    B3[0, 2] = dNdr[1]
	    B3[0, 4] = dNdr[2]
	    B3[0, 6] = dNdr[3]
	    B3[1, 0] = dNds[0]
	    B3[1, 2] = dNds[1]
	    B3[1, 4] = dNds[2]
	    B3[1, 6] = dNds[3]
	    B3[2, 1] = dNdr[0]
	    B3[2, 3] = dNdr[1]
	    B3[2, 5] = dNdr[2]
	    B3[2, 7] = dNdr[3]
	    B3[3, 1] = dNds[0]
	    B3[3, 3] = dNds[1]
	    B3[3, 5] = dNds[2]
	    B3[3, 7] = dNds[3]
	
	    B = B1*B2*B3
	    return B
	
	def planeStressC(self,E,nu):# 弹性系数矩阵
	    C = np.zeros((3,3))
	    cons = E/(1+nu)
	    C[0, 0] = C[1, 1] = cons*1/(1-nu)
	    C[0, 1] = C[1, 0] = cons * nu / (1 - nu)
	    C[2, 2] = cons * 1 / 2
	    return C
	
	def func_k(self):# 根据参数得到4个单元刚度矩阵
	    intPoint = [-1 / np.sqrt(3), 1 / np.sqrt(3)]  # 二阶高斯点坐标
	    weight = [1.0, 1.0]  # 权重
	    E = 1 # 弹性模量
	    nu = 0.3 # 泊松比
	    C = self.planeStressC(E,nu)
	    x0 = []
	    y0 = []
	    for p in self.pos:
	    	x0.append(p[0])
	    	y0.append(p[1])
	    x = [x0[self.i - 1], x0[self.j - 1], x0[self.k - 1], x0[self.l - 1]]
	    y = [y0[self.i - 1], y0[self.j - 1], y0[self.k - 1], y0[self.l - 1]]
	    K = np.zeros((8,8))
	    for j in range(0, 2): # 数值积分求K
	        for i in range(0, 2):
	            r = intPoint[i]
	            s = intPoint[j]
	            Jinv, Jdet = self.jacobian(x, y, r, s)
	            B = self.Bmatrix(r, s, Jinv)
	            BT = np.transpose(B)
	            tmp = BT*C*B*Jdet
	            K = K + tmp
	    self.K = K
	    return K

def matrix_assembly(kk, Elementset):# 总刚度矩阵组装
    for e in Elementset:
        i, j, k, l= e.i, e.j, e.k, e.l
        for m in range(0, 2):
            for n in range(0, 2):# 将16个2x2矩阵按照索引值加入总刚度矩阵中
                kk[2*(i - 1)+m,2*(i - 1)+n] += e.K[m,n]
                kk[2*(j - 1)+m,2*(i - 1)+n] += e.K[2+m,n]
                kk[2*(k - 1)+m,2*(i - 1)+n] += e.K[4+m,n]
                kk[2 * (l - 1)+m,2 * (i - 1)+n] += e.K[6+m,n]
                kk[2*(i - 1)+m,2*(j - 1)+n] += e.K[m,2+n]
                kk[2*(j - 1)+m,2*(j - 1)+n] += e.K[2+m,2+n]
                kk[2*(k - 1)+m,2*(j - 1)+n] += e.K[4+m,2+n]
                kk[2 * (l - 1)+m,2*(j - 1)+n] += e.K[6+m,2+n]
                kk[2*(i - 1)+m,2*(k - 1)+n] += e.K[m,4+n]
                kk[2*(j - 1)+m,2*(k - 1)+n] += e.K[2+m,4+n]
                kk[2*(k - 1)+m,2*(k - 1)+n] += e.K[4+m,4+n]
                kk[2 * (l - 1)+m,2*(k - 1)+n] += e.K[6+m,4+n]
                kk[2*(i - 1)+m,2 * (l - 1)+n] += e.K[m,6+n]
                kk[2*(j - 1)+m,2 * (l - 1)+n] += e.K[2+m,6+n]
                kk[2*(k - 1)+m,2 * (l - 1)+n] += e.K[4+m,6+n]
                kk[2 * (l - 1)+m,2 * (l - 1)+n] += e.K[6+m,6+n]
    return kk
