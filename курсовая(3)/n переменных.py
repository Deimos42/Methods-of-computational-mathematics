import numpy as np
import matplotlib.pyplot as plot
from matplotlib import cm
# градиентный метод минимизации функции многих переменных

# TODO: погрешность линейных функций не уменьшается(или не использовать их или для них дробный шаг использовать)
# TODO: переполнение памяти при большем кол-ве итераций(оптимизировать)
# TODO: придумать интересный пример
# TODO: X задавать там же где и A,S(но X должен быть глобальным)

# Не у всех функций верно находит минимум

a = 0.01 # 0.01
global_epsilon = 0.000001
chislo = 0 # число итераций
global X
X = np.array((10.0,10.0)) # TODO: подумать как разместить X. Пока его приходится менять для каждого теста
# X - значения иксов. Изначально соотвествуют начальной точке, потом изменяется на каждом шаге
# A - коэффициенты при x.
def test1(): # основной
    #global X 
    #X = np.array((10.0,10.0))
    A = np.array((1.0,1.0))
    S = np.array((2,2)) # степени при X
    b = 5
    return A,S,b

def test11():
    #global X 
    #X = np.array((10.0,10.0))
    A = np.array((1.0,1.0))  # F = x**3 + y**2
    S = np.array((3,3)) # степени при X
    b = -50
    return A,S,b

def test2():
    #global X
    #X = np.array((10.0,10.0,10.0,10.0))
    A = np.array((1.0,2.0,-3.0,-4.0))
    S = np.array((3,2,1,1))
    b = 0
    return A,S,b

def test3():
    #global X
    #X = np.array((10.0,))
    A = np.array((10.0,))
    S = np.array((3,))
    b = 5.0 # слагаемое 0 степени
    return A,S,b

def test4():
    #F = x**2 - y**2
    #global X
    #X = np.array((10.0,10.0))
    A = np.array((1.0,-2.0))
    S = np.array((2,2))
    b = 0
    return A,S,b

# тест для квадратичной функции
def test_kvadra():
    #global X
    #X = np.array((10.0,10.0))
    #(x ** 3) - 5*x*y + y**2 
    A2 = np.array((
    (1.0, 5.0),
    (5.0, 1.0),
    ))
    b = 50
    return A2,b



def function1(X,index=-1,meaning=-1):
    X_new = X.copy()
    if index != -1:  
        X_new[index] = meaning

    if name_test == 'test1':
        A,S,b = test1()
    elif name_test == 'test11':
        A,S,b = test11()
    elif name_test == 'test2':
        A,S,b = test2()
    elif name_test == 'test3':
        A,S,b = test3()
    elif name_test == 'test4':
        A,S,b = test4()
    else:
        print('такой тест не задан')
        return

    F2 = 0
    F = A*X_new**S
    for i in F:
        F2 += i
    F = F2
    F += b # если добавлять b раньше, то считает непправильно(значение b умножается на два, так как вектор F состоит из двух элементов)

    return F


def function2(X,index=-1,meaning=-1):
    # --------------------------------------------------------------------------- 
    # для квадратичной функции
    X_new = X.copy()
    if index != -1:  
        X_new[index] = meaning

    if name_test == 'test_kvadra':
        A2,b = test_kvadra()
    else:
        print('такой тест не задан')
        return

    F = 0
    F2 = 0
    F = X_new * A2 * X_new.transpose() # квадратичная матрица
    #F = X * A2
    #F = F * X.transpose()
    for i in F:
        for j in i:
            F2 += j
    F = F2
    F += b # если добавлять b раньше, то считает непправильно(значение b умножается на два, так как вектор F состоит из двух элементов)
    # ---------------------------------------------------------------------------

    return F



# eplsilon - значение дифференцируемой переменной
# index - индекс переменной. Остальные значения переменных берутся из массива X
def derivative(epsilon,index):
    if name_test == 'test_kvadra':
        return (function2(X,index,epsilon + global_epsilon) - function2(X,index,epsilon)) / global_epsilon
    else: 
        return (function1(X,index,epsilon + global_epsilon) - function1(X,index,epsilon)) / global_epsilon

# TODO: если точность не уменьшается со временем, то вывести сообщение что метод не сходится
# или уменьшить шаг
def minimiz1(X,a,global_epsilon,chislo):   # для функции одной переменной
    if name_test == 'test_kvadra':
        F = function2(X) # TODO: должны получить значение функции с параметрами из X
    else:
        F = function1(X)
    y = 0
    fix = 0 # фиксируем значение. Если оно долго не изменяется - останавливаем алгоритм
    #while True: # добавить условие критерия останова

    gradient = 0
    for k in range(len(X)):
        gradient += derivative(X[k],k)
        #print('gradient: ',gradient)
    print('gradient: '.upper(),gradient)

    X2 = X.copy()
    for i in range(len(X)):
        X2[i] = X[i] - a*gradient
    print('X2: ',X2)

    if name_test == 'test_kvadra':
        F2 = function2(X2) # TODO: должны получить значение функции с параметрами из X
    else:
        F2 = function1(X2)

    print('F2(x2): ',F2)
    print('F1(x2): ',F)
    print('abs(F2-F): ',abs(F2-F),'\n')
    chislo += 1

    if abs(F2-F) <= global_epsilon: # критерий останова
        return (X2,F2,abs(F2-F),chislo) # решение, достигнутая точность, число итераций
    else:
        Res = minimiz1(X2,a,global_epsilon,chislo)
        return Res

    #return (x2,abs(F2-F),chislo)
    #return


# функция прорисовки графика для функции
# TODO: написать функцию прорисовки графика для квадратичной функции
def draw(point,X,N):  # N - число разбиений по каждой оси
    point_x, point_y, point_z = point
    fig = plot.figure(figsize = (X[0], X[1]))   # (figsize = (X_old[0], X_old[1])) 
    ax = fig.add_subplot(1, 1, 1, projection = '3d') # размещаем оси на сетке
    xval = np.linspace(-X[0],X[1],N)                                                
    yval = np.linspace(-X[0],X[1],N)  
    X, Y = np.meshgrid(xval, yval) #  список координатных сеток

    Resul = np.zeros(N*N).reshape(N,N) # вместо 4,4 X[0],X[1]
    X2 = np.zeros(2) 

    if name_test == 'test1':
        A,S,b = test1()
    elif name_test == 'test11':
        A,S,b = test11()
    elif name_test == 'test2':
        A,S,b = test2()
    elif name_test == 'test3':
        A,S,b = test3()
    elif name_test == 'test4':
        A,S,b = test4()
    else:
        print('такой тест не задан')
        return

    for i in range(len(X)):
        for j in range(len(X)):
            #print(i,j)
            X2[0], X2[1] = X[i][j], Y[i][j]
            #print(X2)

            F2 = 0.0
            F = A*X2**S
            for k in F:
                F2 += k
            F = F2
            F += b
            Resul[i][j] = F         
            F = 0
            print('')

    ax.scatter(point_x, point_y, point_z, color='red')
    surf = ax.plot_surface(X, Y, Resul, rstride = 5,cstride = 5, cmap = cm.plasma)
    plot.show()
    return



# функция прорисовки графика для функции
# TODO: написать функцию прорисовки графика для квадратичной функции
# для прорисовки графика функций, заданных квадратически
#def draw2(point,X,N):  # N - число разбиений по каждой оси
#    point_x, point_y, point_z = point
#    fig = plot.figure(figsize = (X[0], X[1]))   # (figsize = (X_old[0], X_old[1])) 
#    ax = fig.add_subplot(1, 1, 1, projection = '3d') # размещаем оси на сетке
#    xval = np.linspace(-X[0],X[1],N)                                                
#    yval = np.linspace(-X[0],X[1],N)  
#    X, Y = np.meshgrid(xval, yval) #  список координатных сеток
#
#    Resul = np.zeros(N*N).reshape(N,N) # вместо 4,4 X[0],X[1]
#    X2 = np.zeros(2)
#
#    if name_test == 'test_kvadra':
#        A,b = test_kvadra()
#    else:
#        print('такой тест не задан')
#        return
#
#    for i in range(len(X)):
#        for j in range(len(X)):
#            #print(i,j)
#            X2[0], X2[1] = X[i][j], Y[i][j]
#            #print(X2)
#
#            F2 = 0.0
#            F = X2 * A2 * X2.transpose() 
#            for s in F:
#                for d in s:
#                    F2 += d
#            F = F2
#            F += b 
#            Resul[i][j] = F
#            F = 0
#        print('')
#
#    ax.scatter(point_x, point_y, point_z, color='red')
#    surf = ax.plot_surface(X, Y, Resul, rstride = 5,cstride = 5, cmap = cm.plasma) 
#    plot.show()
#    return

# рисует график для двумерной функции
def draw_2d(X,N):
    fig = plot.figure(figsize = (10, 10))   # (figsize = (X_old[0], X_old[1])) 
    ax = fig.add_subplot(1, 1, 1, projection = '3d') # размещаем оси на сетке
    X = np.linspace(-X,X,N)
    print(type(X))
    vv = 0
    # TODO: не работает - исправить
    for g in X:
        print(X[vv])
        print(g[0])
        X[vv] = g[0]
        vv += 1
 
    #X, Y = np.meshgrid(xval, yval) #  список координатных сеток

    print(X)
    Resul = np.zeros(N) 
    X2 = np.zeros(1) 

    if name_test == 'test1':
        A,S,b = test1()
    elif name_test == 'test11':
        A,S,b = test11()
    elif name_test == 'test2':
        A,S,b = test2()
    elif name_test == 'test3':
        A,S,b = test3()
    elif name_test == 'test4':
        A,S,b = test4()
    else:
        print('такой тест не задан')
        return


    for i in range(len(X)):
        X2 = X[i]
        F2 = 0.0
        F = A*X2**S
        for k in F:
            F2 += k
        F = F2
        F += b
        Resul[i] = F
        F = 0
        print('')



    #surf = ax.plot_surface(X, Y, Resul, rstride = 5,cstride = 5, cmap = cm.plasma) 
    plot.plot(X,Resul)
    plot.show()
    return


name_test = input('Введите имя теста: ')

Resul = ()
Resul = minimiz1(X,a,global_epsilon,chislo)
print('-------------------------------------------------')
# решение, достигнутая точность, число итераций
print('найденная точка: ',Resul[0])
print('значение в этой точке: ',Resul[1])
print('достигнутая точность: ',Resul[2])
print('число итераций: ',Resul[3])
print('шаг: ',a)

min_x = Resul[0][0]
print(len(Resul))
point = []
for k in Resul[0]:
    point.append(k)
point.append(Resul[1])
print(point)

if len(Resul[0]) == 1:
    #minimum = (min_x,Resul[1])
    draw_2d(X,100)

if len(Resul[0]) == 2:
    min_y = Resul[0][1]
    #minimum = (min_x,min_y,Resul[1]) # может не нужно 
    draw(point,X,100)

