import time
import math

# Вывод матрицы
def vivod(A):
    print('матрица:')
    for k in A:
        print(k)
    print('\n')

# генерация первоначальной квадратной матрицы A размера NxN
def generate(N,k=9):
    A = []
    K = N // 2
    for i in range(1,N+1):
        res = []
        for j in range(1,N+1):
            if i == j:
                res.append(min(N,10))
            elif abs(i-j) > K:
                res.append(0)
            else:
                res.append(1/j)
        A.append(res) 
        
    return A

# формирование вектора правых частей b
def vector_b(A,N,k=9):
    x = []
    for v in range(1,N+1): # здесь индексы начинаем с 1 для удобства(в задании они указаны с 1)
        if v == 1:
            x.append(1)
        elif v == (1+k):
            x.append(1)
        elif v == (1+5*k):
            x.append(1)
        else:
            x.append(0)
    print('x: ',x) 

    b = []
    for i in range(len(A)):
        res = 0 # значение одного элемента вектора b
        v = 0
        for j in range(len(A[i])):
            #print(A[i][j])
            res += A[i][j] * x[v]
            v += 1
        #print('res: ',res)
        b.append(res)

    return [b,x]


# Вращение Гивенса
def givens(a,b):
    if b == 0:
        c = 1
        s = 0
    else:
        if abs(b)>abs(a):
            t = -a/b
            s = 1/(math.sqrt(1+t*t))
            c = s*t
        else:
            t = -b/a
            c = 1/(math.sqrt(1+t*t))
            s = c*t

    return [c,s]


# Умножение на матрицу Гивенса
def row_rot(A,k,j,c,s): # i,k - индексы строк матрицы, которые будут изменяться
    q = len(A[0]) # кол-во столбцов матрицы
    for i in range(q):
        t1 = A[j][i] 
        t2 = A[k][i]
        A[j][i] = t1*c + t2*s
        A[k][i] = t2*c - t1*s
    vivod(A)

    return A


# приписываем к матрице A справа столбец b
def Add_b(A,b,N):
    v = 0
    A2 = []
    for st in A:
        res = st
        res.append(b[v])
        A2.append(res)
        v += 1

    return A2
    

# QR разложение с Гауссовским порядком обхода элементов(мой вариант)
def gauss(A,N):
    for k in range(0,N-1):
        for j in range(1+k,N):     
            c,s = givens(A[k][k],A[j][k])
            print('k: ',k,' j: ',j)
            print('cos: ',c)
            print('sin: ',s)
            A = row_rot(A,k,j,c,s)
           

# решение системы линейных уравнений после приведения матрицы к верхнетреугольному виду
def syst(A,N):
    b = []
    # отделяем вектор b от матрицы A
    for i in range(len(A)):
        b.append(A[i][-1])
        A[i] = A[i][:-1]

    X = [] 
    x = 0
    k = -1
    j = 1
    for st in reversed(A):
        x = b[k] 
        for i in range(j):
            # умножаем на bool(k != (-i-1)) чтобы не вычитался диагональный элемент
            if k < -1: # на первом шаге ничего не вычитаем
                #x = x - st[-i-1]*bool(k != (-i-1))*X[i]
                x = x - st[-i-1]*X[i]
        x = x / st[k]
        if k < -1:
            j += 1
        X.append(x)
        k -= 1
    X = X[::-1] # так как мы добавляли в конец то вектор x вышел перевернутым

    return X

# считаем абсолютную погрешность на основе второй нормы
def pogr(X_true,x_calcul): 
    A_true = 0
    A_calcul = 0
    print('X_true: ',X_true)
    print('x_calcul: ',x_calcul)
    for v in X_true:
        A_true += abs(v)

    for k in x_calcul:
        A_calcul += abs(k)
        
    return abs(A_true-A_calcul)

# тест
def test(N):
    t1 = time.time()
    A = generate(N) # генерируем матрицу размера N
    vivod(A)

    b,x = vector_b(A,N) # вычисляем вектор b
    print('b: ',b)
    print('Вектор x который должен получиться: ',x)
    A = Add_b(A,b,N) # столбец b припписываем к матрице A
    vivod(A) 
    gauss(A,N) # QR разложение

    X = syst(A,N) # решаем верхнетреугольную систему
    print('Посчитанное решение X: ',X)
    delta = pogr(x,X) # абсолютная погрешность
    print('абсолютная погрешность: ',delta)
    t2 = time.time()
    print('время вычисления при N=',N,': ',t2-t1)


test(4)
#test(10)
#test(50)
#test(100)











