import time
import math
import numpy as np

# 10 вариант -> k = N mod 5 = 0

def sign(x):  # основная
    if x >= 0:
        return 1
    elif x < 0:
        return -1


def norm(X):   # основная
    m = 0
    for i in X:
        m += i*i  
    m = math.sqrt(m)
    return m
 
def make_householder2(a): 
    v = a / (a[0] + norm(a)*sign(a[0]))
    v[0] = 1
    H = np.eye(a.shape[0])
    H -= (2 / np.dot(v, v)) * np.dot(v[:, None], v[None, :])
    return H


def qr(A):
    print('----------- qr ----------------------------')
    m, n = A.shape
    Q = np.eye(m)
    print('Q: \n',Q)
    print('m: ',m,' n: ',n) #
    #A2 = A.transpose()  #
    PP = np.eye(A.shape[1],A.shape[1])
    for i in range(n - (m == n)):

        # выбор ведущего столбца в A
        A2 = A.transpose() #
        res = {}
        for k in range(i,len(A2)): # по столбцам A
            stolb = 0 
            for v in range(i,len(A2[k])): # по строкам A
                stolb += A2[k][v]**2
            stolb = math.sqrt(stolb) # вторая норма столбца
            #print('stolb: ',stolb,' k: ',k)
            res[k] = stolb
            print('норма ',k,'столбца матрицы = ',stolb)
        maxx = max(res.values())
        print('норма максимального столбца: ',maxx)  #    
        if maxx < 0.000001: # останавливаем выполнение разложения
            pass 
        else:  # переставляем i-ый столбец с maxx в матрице A          
            for l,h in res.items(): # получаем l - индекс 'максимального' столбца матрицы A
                if h == maxx:
                   break
            
            # формируем матрицу перестановок и умножаем матрицу A на неё
            P = np.eye(A.shape[1],A.shape[1])
            temp = np.copy(P[l])
            P[l] = P[i]
            P[i] = temp     
            PP = np.dot(PP,P) # сохраняем матрицу перестановок
            #print('A: \n',A.round(6))

            print('Выбираем',l,'столбец. Переставляем его с',i,'столбцом')
            A = np.dot(A,P)
            print('A: \n',A.round(6))
        
        #print('A[i:, i]: ',A[i:, i]) # первый столбец очередного минора матрицы
        H = np.eye(m)
        H[i:, i:] = make_householder2(A[i:, i]) 
        #print('H: \n',H) #
        #print('H[i:, i:]: \n',H[i:, i:]) #
        #print('i: ',i) #
        Q = np.dot(Q, H)
        print('Q: \n',Q.round(6))
        A = np.dot(H, A) # после всех перемножений A будет верхнетреугольной матрицей
        print('A: \n',A.round(6))
        print('')
    print('----------- qr ----------------------------')
    return Q, A, PP
 
 

# решение системы линейных уравнений после приведения матрицы к верхнетреугольному виду
def syst(R,c): # c - вектор правых частей, R - верхнетреугольная матрица
    X = [] 
    x = 0
    k = -1
    j = 1
    for st in reversed(R):
        x = c[k] 
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


# считаем норму невязки
def pogr(A,b,X):
    A_calcul = 0
    D = np.dot(A,X) 
    print('A*X: \n',D.round(6))
    print('b: \n',b)
    r = 0
    for v in range(len(D)):
        r += (D[v]-b[v])**2
    r = math.sqrt(r)
    
    return r

def test1():
    a = np.array(((
        (0,1,-1,0,1),
        (1,0,1,0,0),
        (1,0,1,1,0),
        (0,1,1,0,-1),
        (0,0,1,0,-1),
        (0,1,1,0,1),
    )))
    b = np.array((1.0,1.0,1.0,0.0,1.0,0.0))
    return a,b

# пример с нулевой нормой невязки
def test2():
    a = np.array(((
        (12, -51,   4),
        ( 6, 167, -68),
        (-4,  24, -41),
    )))
    b = np.array((1.0,1.0,1.0))
    return a,b

a,b = test2()
Q, R, PP = qr(a) # PP - матрица перестановок
print('Q:\n', Q.round(6))
print('R:\n', R.round(6))
print('матрица перестановок PP:\n', PP.round(6))

# проверка: если получается единичная матрица, то q - ортогональная матрица, значит qr разложение находится верно
#print('Q.transpose(): \n',Q.transpose().round(6)) 
#vv = np.dot(Q,Q.transpose()) #
#print('Q*Q.transpose: \n',vv.round(2)) #

#print('r: \n',np.dot(q.transpose(),a).round(2))
b2 = np.dot(Q.transpose(),b)
print('b2: \n',b2.round(6)) # находим изменённое b

m, n = a.shape # число строк и столбцов

# находим 'c' и 'r2' для решения системы уравнений
c = b2[:n]
R2 = R[:n]
print('c: \n',c.round(6))
print('R2: \n',R2.round(6))
# решаем систему уравнений
X = syst(R2,c) 
X = np.dot(PP,X)
print('Посчитанное решение X: ',X.round(6))

#vvv = np.dot(a,X)
#print(vvv)


delta = pogr(a,b,X) # абсолютная погрешность
print('норма невязки: ',delta)

# проверяем верность равенства
print('проверяем верность равенства')
print('A*П: \n',np.dot(a,PP).round(6)) # PP - итоговая матрица перестановок П
print('Q*R: \n',np.dot(Q,R).round(6))


