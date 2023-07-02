import  fem

#femt = fem.FEM('example3.txt')
femt = fem.FEM('example4.txt')
femt.model()
femt.kkk()
femt.solve()