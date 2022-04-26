import numpy as np
import numpy.linalg as linalg #

# вычисляем матрицу S+
def MatrixSPlus(S):
    m = S.shape[0]
    n = S.shape[1]
    Splus = np.zeros((n,m))
    for i in range(n):  
        if S[i][i] < 10**(-10):
            Splus[i][i] = 0
        else:
            Splus[i][i] = 1/S[i][i]
    return Splus

# ищем ранг
def Rank(M):
    n = M.shape[0] # n = matrix.shape[1]
    r = 0
    while r < n and M[r][r]>10**(-10): 
        print(M[r][r])
        r += 1
    return r

# мой вариант(0)
def test1():
    A = np.array((
        (0.0, 1.0, -1.0, 0.0, 1.0),
        (1.0, 0.0, 1.0, 0.0, 0.0),
        (1.0, 0.0, 1.0, 1.0, 0.0),
        (0.0, 1.0, 1.0, 0.0, -1.0),
        (0.0, 0.0, 1.0, 0.0, -1.0),
        (0.0, 1.0, 1.0, 0.0, 1.0)
        ))
    b = np.array([1.0,1.0,1.0,0.0,1.0,0.0])
    return A,b

A,b = test1()
m,n = A.shape

# разложение для A
U,vector,VTrans = np.linalg.svd(A,True,True)
V = VTrans.transpose()
UTrans = U.transpose()

# строим диагональную матрицу S 
S = np.zeros((m,n))
for i in range(len(vector)):
    S[i][i] = vector[i]


print('Матрица U\n',U.round(6),'\n')
print('Матрица S\n',S.round(6),'\n')
print('Матрица V\n',V.round(6),'\n')

# проверка V и U на ортогональность
print('Матрица U*UTrans\n',np.dot(U,UTrans).round(6),'\n')
print('Матрица UTrans*U\n',np.dot(UTrans,U).round(6),'\n')

print('Матрица V*VTrans\n',np.dot(V,VTrans).round(6),'\n')
print('Матрица VTrans*V\n',np.dot(VTrans,V).round(6),'\n')

print('Матрица S\n',S.round(6),'\n')

# проверка правильности нахождения S
S_test = np.dot(A,V)
S_test = np.dot(UTrans,S_test)
print('Проверочная матрица S_test: \n',S_test.round(6),'\n')

SPlus = MatrixSPlus(S)
print('Матрица SPlus: \n',SPlus.round(6),'\n')

# вычисляем 'c'
rank = Rank(SPlus)
c = np.dot(UTrans,b) 
c[rank:] = 0
print('c: ',c)

# проверочное нахождение матрицы A
#A = np.dot(U,S)
#A = np.dot(A,VTrans)
#print('A: \n',A.round(6))

# находим матрицу A+
APlus = np.dot(SPlus,UTrans)
APlus = np.dot(V,APlus)
print('Псевдообратная матрица A+: \n',APlus.round(6),'\n')

##A = A.transpose()
AA = np.dot(A,APlus)
print(' Матрица A*A+: \n',AA.round(6),'\n')

AA2 = np.dot(APlus,A)
print('Матрица A+*A: \n',AA2.round(6),'\n')

# вычисляем y
y = np.dot(SPlus,c)
print('y: ',y)

# находим x
x = np.dot(V,y)
print('Решение с помощью SVD-разложения: ')
print('x: ',x.round(6))

x2 = np.dot(APlus,b)
print('Решение с помощью псевдообратной матрицы: ')
print('x: ',x2.round(6))

