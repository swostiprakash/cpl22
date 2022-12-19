import math as mt
import matplotlib.pyplot as plt

# Input Functions

def read_integer(path): 
    """
        Taking inputs that are just mere numbers.
    """
    with open(path) as f:
        N = []
        for line in f:
            if line[0] != '#' and line.strip() != '':
                N.append(float(line))
    return N

def read_mat(path):
    """
        Taking matrix input.
    """
    with open(path) as f:
        A = []
        for line in f:
            if line[0] != '#' and line.strip() != '':
                A.append([float(num) for num in line.split(' ')])              
    return A
    
def read_vec(path):
    """
        Taking vector input.
    """
    with open(path) as f:
        A = []
        for line in f:
            if line[0] != '#' and line.strip() != '':
                A.append(float(line))
    return A
    
def split_mat(A):
    """
        Splitting a n x 2 matrix into two n x 1 matrices.
    """
    x = []
    y = []
    for i in range(len(A)):
        x.append(A[i][0])
        y.append(A[i][1])
    return x,y

def read_csv(path):    
    """
        Taking csv input in case of related data.
    """
    with open(path) as f:
        A = []
        for line in f:
            B = []
            for elem in line.split(','):
                try:
                    elem = float(elem)
                    B.append(elem)
                except:
                    B.append(elem)
            A.append(B)
        C, D = split_mat(A)
        A = []
        A.append(C)
        A.append(D)
        A[0].pop(0)
        A[1].pop(0)
    return A

# Random Number Generator

a = 10

def lcg():

    global a
    a = ((1103515245*a+12345)%32768)/32768 # Formula for the prng
    return a

# Sum/Product of different sets of numbers

def sum_odd(N):
    s = 0
    num = 1
    for i in range(0,N):
        s = s + num
        num = num + 2
    return s
    
def fact(N):
    """
        Returns the factorial of N.
    """
    fac = 1
    num = 1
    for i in range(0,N):
        fac = fac*num
        num = num + 1
    return fac
    
def AP(a,d,N):
    """
        Returns the sum of a given AP upto N terms.
    """
    s = 0
    i = 0
    t = a
    while i is not N:
        s = s + t
        t = t + d
        i = i + 1
    return s
    
def GP(a,r,N):
    """
        Returns the sum of a given GP upto N terms.
    """
    i = 0
    t = a
    s = 0
    while i is not N:
        s = s + t
        t = t*r
        i = i + 1
    return s
    
def HP(a,d,N):
    """
        Returns the sum of a given HP upto N terms.
    """
    i = 0
    s = 0
    t = a # I have taken the first tern to be a and not 1/a
    while i is not N:
        s = s + t
        r = 1/t
        t = 1/(r+d)
        i = i + 1        
    return s

# Matrix Operations

def mat_mult(A,B):
    """
        Returns the product of two matrices.
    """
    cA = len(A[0]) # Number of columns in A
    rB = len(B) # Number of rows in B
    if cA == rB: # Condition for valid multiplication
        mat = []
        for i in range(len(A)): # Initialising an empty matrix
            mat.append([])
        for i in range(len(A)): # Select the row of matrix A
            for k in range(len(B[0])): # Select the column of matrix B
                s = 0
                for j in range(len(B)): # Select the row of matrix B
                    s = s + A[i][j]*B[j][k]
                mat[i].append(round(s,2))
        return mat
    else:
        return 'Invalid Operation'
        
def dot(A,B):
    """
        Returns the dot product of two vectors.
        Note that it takes a matrix input and not a vector input.
    """
    s = 0
    for j in range(len(A)):
        s = s + A[j][0]*B[j][0]
    return s

def trans(L):
    """
        Returns the transpose of a matrix.
    """
    for i in range(len(L)):
        for j in range(len(L)):
            if i < j:
                L[i][j], L[j][i] = L[j][i], L[i][j]
    return L

# Complex Numbers

class myComplex:
    def __init__(self,real,img):
        self.r = real
        self.i = img
    
    def sum_complex(self,r2,i2):
        sr = self.r + r2
        si = self.i + i2
        return sr, si
    
    def modulus(self):
        mod = (self.r**2 + self.i**2)**0.5
        return mod
    
    def prod_complex(self,r2,i2):
        pr = self.r*r2 - self.i*i2
        pi = self.r*i2 + self.i*r2
        return pr,pi

# System of Simulataneous Equations Solvers
    
def dec(A):
    """
        Returns the LU decomposed matrix.
    """
    n = len(A)
    for j in range(0,n):
        for i in range(1,j+1):
            s = 0
            for k in range(0,i):
                s = s + A[i][k]*A[k][j]
            A[i][j] = A[i][j] - s
        for i in range(j+1,n):
            s = 0
            for k in range(0,j):
                s = s + A[i][k]*A[k][j]
            A[i][j] = (A[i][j] - s)/A[j][j]
            
def LU(A, B):
    """
        Returns the solution to a system of simultaneous linear equations using LU decomposition.
        Note that it takes matrix input and not a vector input.
    """
    dec(A)
    y = []
    for i in range(len(A)):
        y.append([None])
    for i in range(len(A)):
        s = 0
        for j in range(0,i):
            s = s + A[i][j]*y[j][0]
        y[i][0] = B[i][0] - s
    x = []
    for i in range(len(A)):
        x.append([None])
    for i in range(len(A)-1,-1,-1):
        s = 0
        for j in range(i+1,len(A)):
            s = s + A[i][j]*x[j][0]
        x[i][0] = (y[i][0] - s)/A[i][i]
    return x

def ch_dec(A):
    """
        Returns the Cholesky decomposed matrix.
    """
    for i in range(len(A)):
        for j in range(len(A[0])):
            s = 0
            if i == j:
                for t in range(0,i):
                    s = s + A[t][i]**2
                A[i][i] = round((A[i][i] - s)**0.5,3)
            if i < j:
                for k in range(0,i):
                    s = s + A[k][i]*A[k][j]
                A[i][j] = round((A[i][j] - s)/A[i][i],3)
    trans(A)
    for i in range(len(A)):
        for j in range(len(A)):
            if i < j:
                A[i][j] = 0
def trans(L):
    for i in range(len(L)):
        for j in range(len(L)):
            if i < j:
                L[i][j], L[j][i] = L[j][i], L[i][j]
                
def choleski(A, B):
    """
        Returns the solution to a system of simultaneous linear equations using Cholesky decomposition.
        Note that it takes matrix input and not a vector input.
    """
    ch_dec(A)
    y = []
    for i in range(len(A)):
        y.append([0])
    for i in range(len(A)):
        s = 0
        for j in range(0,i):
            s = s + A[i][j]*y[j][0]
        y[i][0] = (B[i][0] - s)/A[i][i]
    return y

def is_sym(A):
    """
        Checks if a matrix is symmetric.
    """
    for i in range(len(A)):
        for j in range(len(A[0])):
            if A[i][j] != A[j][i]:
                return False
    return True

def Jacobi(A,B):
    """
        Returns the solution to a system of simultaneous linear equations using Jacobi Method.
        Note that it takes both matrix (A) and vector (B) input.
    """
    x = []
    y = []
    t = 1
    E = 0.0001
    c = 0
    for i in range(len(A)):
        x.append(1)
        y.append(0)

    while abs(t) > E:
        c = c + 1
        for i in range(len(x)):
            s = 0
            for j in range(len(A)):
                if j != i:
                    s = s + A[i][j]*x[j]
            y[i] = (B[i] - s)/A[i][i]
            t = (y[i] - x[i])*100/x[i]
            x[i] = y[i]
    print("Steps =",c)
    return x

def Seidel(A,B):    
    """
        Returns the solution to a system of simultaneous linear equations using Gauss-Seidel Method.
        Note that it takes both matrix (A) and vector (B) input.
    """
    E = 0.0001
    t = 1
    c = 0
    x = []
    y = []
    for i in range(len(A)):
        x.append(0)
        y.append(0)
    while c < 20 and t > E:
        c = c + 1
        for i in range(len(x)):
            s1 = 0
            s2 = 0
            for j in range(0, i):
                s1= s1+ A[i][j]*x[j]
            for j in range(i+1,len(x)):
                s2 = s2 + A[i][j]*x[j]
            x[i] = (B[i] - s1 - s2)/A[i][i]
            if y[i] !=0:    
                t = (y[i] - x[i])*100/y[i]
                t = abs(t)
            y[i] = x[i]
        
    print("Steps =",c)
    return x

def pos_def(A):
    """
        Checks if a matrix is positive definite. Not very useful.
    """
    import SPSL.spsl as lib
    x = []
    for i in range(len(A)):
        x.append([lib.lcg(10,10)[i+1]])
    print(x)
    xt = [[]]
    for i in range(len(x)):
        xt[0].append(x[i][0])
    t = lib.mat_mult(A,x)
    if lib.mat_mult(xt,t)[0][0] <= 0:
        return False
    return True

def Gauss(aug):
    """
        Returns the solution to a system of simultaneous linear equations using Gauss-Jordan Method.
        Note that it takes an augmented matrix as an input.
    """
    for j in range(0,len(aug)):
        if aug[j][j] == 0:
            k = j
            max = aug[j][j]
            for l in range(j+1,len(aug)):
                if aug[l][j] > max:
                    max = aug[l][j]
                    k = l
            aug = swap(aug,j,k)
        if aug[j][j] != 1:
            k = j
            for l in range(j+1,len(aug)):
                if aug[l][j] != 0 and l != j:
                    k = l
                    break
            n = aug[j][j] - 1
            m = aug[k][j]
            t = n/m
            aug[j] = sub_row(aug[j],mult_row(aug[k],t))
        
        for i in range(0,len(aug)):
            if aug[i][j] != 0 and i != j:
                aug[i] = sub_row(aug[i],mult_row(aug[j],aug[i][j]))
    y = []
    for i in range(len(aug)):
        y.append(aug[i][len(aug)])
    return y

def swap(mat,i,j): # Swap two rows of a matrix
    dup = mat[i]
    mat[i] = mat[j]
    mat[j] = dup
    return mat

def mult_row(row, fac): # Multiply a row with a factor
    new_row = [i*fac for i in row]
    return new_row

def sub_row(row1,row2): # Subtract one row from another
    new = []
    for i in range(len(row1)):
        new.append(row1[i] - row2[i])
    return new

def diag_dom(A, B):
    """
        Returns the diagonally dominant matrix if possible.
    """
    for i in range(len(A)):
        max = A[i][i]
        m = i
        for j in range(len(A)):
            if A[i][j] > max:
                max = A[i][j]
                m = j
        if m != i:
            A = swap(A, i, m)
            B = swap(B, i, m)
            i = 0
    for i in range(len(A)):
        s = 0
        for j in range(len(A)):
            if j != i:
                s = s + A[i][j]
        if A[i][i] < s:
            print('Sorry this matrix cannot be made diagonally dominant')
            return False
def aug(x,y):
    """
        Retruns the augmented matrix.
    """
    aug = []
    for i in range(len(x)):
        aug.append(x[i][:])
    for i in range(len(y)):
        aug[i].append(y[i])
    return aug

# Data Fitting Functions

def fvalue(m,x,c):
    """
        A helping function for linear fit. Returns the functional value at a particular point.
    """
    return (m*x + c)

def linearfit(x,y,lab1,lab2):
    """
        Linearly fits the data and plots the same alongwith the original datapoints. 
        Returns the slope, y-intercept, the standard errors in slope and intercept and, the Pearson's correlation coefficient.
    """
    n = len(x)
    sx = 0
    sy = 0
    sx2 = 0
    sy2 = 0
    sxy = 0
    t = 0
    for i in range(len(x)):
        sx = sx + x[i]
        sy = sy + y[i]
        sx2 = sx2 + x[i]**2
        sy2 = sy2 + y[i]**2
        sxy = sxy + x[i]*y[i]
    xm = sx/n
    ym = sy/n
    m = (n*sxy - sx*sy)/(n*sx2 - sx**2) # Slope of the best fit
    c = ym - xm*m # Intercept of the best fit
    r = (n*sxy - sx*sy)/((n*sx2 - sx**2)*(n*sy2 - sy**2))**0.5 # Pearson's correlation ratio
    for i in range(len(x)): 
        t = t + (y[i] - m*x[i] - c)**2
    s = (t/(n-2))**0.5
    e1 = s*(n/(n*sx2 - sx**2))**0.5 # Standard error in slope
    e2 = s*(sx2/(sx2*n-sx**2))**0.5 # Standard error in intercept
    m1 = []
    c1 = []
    for i in range(len(x)):
        m1.append(m)
        c1.append(c)
    fit = list(map(fvalue,x,m1,c1))
    lab = 'Fit: y = ' + str(round(m,3)) + 'x + ' + str(round(c,3))
    plt.scatter(x,y,color='red',s=8)
    plt.plot(x,fit,label=lab)
    plt.xlabel(lab1)
    plt.ylabel(lab2)
    plt.legend()
    plt.savefig('Plot.png')
    plt.show()
    return m,c,e1,e2,r

def polyfit(data,n):
    """
        Fits the data into a polynomial of given order. 
        Returns the coefficients of the different terms of the polynomial.
    """
    x = []
    y = []
    for i in range(n+1):
        x.append([])
        for j in range(n+1):
            s = 0
            for k in range(len(data)):
                s = s + data[k][0]**(i+j)
            x[i].append(s)
    for i in range(n+1):
        s = 0
        for k in range(len(data)):
            s = s + (data[k][0]**i)*data[k][1]
        y.append(s)
    ag = aug(x,y)
    return Seidel(x,y)

# Polar Plot

def polarplt(theta,I,label):
    """
        Plots the polar plot of given data.
        The angles need to be given in degrees.
    """
    for i in range(len(theta)):
        theta[i] = theta[i]*np.pi/180
    plt.polar(theta,I,'r.')
    plt.polar(theta,I)
    plt.savefig(label)
    plt.show()

# Roots of Nonlinear Equations    

def NewRaph(f,fd,x1):
    """
        Returns the root of a non-linear equation using the Newton-Raphson method.
    """
    E = 0.000001
    n = 0
    x2 = 0
    while abs(x2 - x1) > E or n == 0:
        x1 = x2
        n = n + 1
        x2 = x1 - f(x1)/fd(x1)
    return x2, n
    
def deflate(P,r):
    """
        Synthetic division method implementation.
    """
    for i in range(1,len(P)):
        P[i] = P[i] + P[i-1]*r
    P.pop()

def der(P):
    """
        Returns the derivative of a polynomial. 
    """
    D = []
    for i in range(len(P)):
        j = len(P) - 1 - i
        D.append(j*P[i])
    if len(D)!=0:
        D.pop()
    return D

def Laguere(P, x):
    """
        Implementation of the Laguere's method.
    """
    val = 0
    for i in range(len(P)-1,-1,-1):
        j = len(P) -1 - i
        val = val + (x**j)*P[i]
    if abs(val) <= 10**(-4):
        return x
    else:
        n = len(P)
        D1 = der(P[:])
        D2 = der(D1[:])
        b2 = x
        b1 = 0
        a = 1
        while abs(a) > 10**(-4):
            d1 = 0
            d2 = 0
            val = 0
            for i in range(len(D1)-1,-1,-1):
                j = len(D1) - 1 - i
                d1 = d1 + (b2**j)*D1[i]
            for i in range(len(D2)-1,-1,-1):
                j = len(D2) - 1 - i
                d2 = d2 + (b2**j)*D2[i]
            for i in range(len(P)-1,-1,-1):
                j = len(P) - 1 - i
                val = val + (b2**j)*P[i]
            if abs(val) < 10**(-4):
                return b2
            G = d1/val
            H = G**2 - d2/val
            if G > 0:
                a = n/(G + ((n-1)*(n*H-G**2))**0.5)
            else:
                a = n/(G - ((n-1)*(n*H-G**2))**0.5)
            b1 = b2
            b2 = b1 - a
        val = 0
        for i in range(len(P)-1,-1,-1):
            j = len(P) -1 - i
            val = val + (b2**j)*P[i]
        if abs(val) < 10**(-4):
            return b2
        return False

def Lag(P,g):
    """
        Returns the root of a non-linear equation using the Laguere method.
    """
    r = []
    T = P[:]
    x = g
    for i in range(len(P)-1):
        x = Laguere(T[:],x)
        if x is not False:
            r.append(x)
            deflate(T,x)
        else:
            x = g
            i = i -1
    return r

def bisect(a,b,f):
    """
        Returns the root of a non-linear equation using the Bisection method.
    """
    E = 0.0001
    n = 0
    while abs(f(a) - f(b)) > E and abs(a-b) > E:
        n = n + 1
        c = (a + b)/2
        if f(c)*f(b) < 0:
            a = c
        elif f(c)*f(a) < 0:
            b = c
    root = (a + b)/2
    return root, n

def brac(a,b,f):
    """
        Returns the right interval in which bisection/Regula-Falsi method can be used successfully.
    """
    m = 0
    while f(a)*f(b) > 0:
        m+=1
        if abs(f(a)) > abs(f(b)):
            b = b + 0.01*(b-a)
        elif abs(f(a)) < abs(f(b)):
            a = a - 0.01*(b-a)
    return a,b,m

def bis_root(a,b,f):
    """
        Returns the root of a non-linear equation using the bisection method. 
        It incorporates both bracketing and the subsequent implementation of the bisection method.
    """
    a, b, m = brac(a,b,f)
    if f(a)*f(b) < 0:
        return bisect(a,b,f)
    else:
        return None, None
    
def regula(a,b,f):
    """
        Returns the root of a non-linear equation using the Regula-Falsi method.
    """
    E = 0.0001
    n = 0
    c = b - (b-a)*f(b)/(f(b) - f(a))
    dup = c
    while abs(dup-c) > E or (dup - c) == 0:
        dup = c
        n = n + 1
        c = b - ((b-a)*f(b))/(f(b) - f(a))
        if f(c)*f(b) < 0:
            a = c
        elif f(c)*f(a) < 0:
            b = c
    root = c
    return root, n

def falsi(a,b,f):
    """
        Returns the root of a non-linear equation using the Regula-Falsi method. 
        It incorporates both bracketing and the subsequent implementation of the Regula-Falsi method.
    """
    a, b, m = brac(a,b,f)
    if f(a)*f(b) < 0:
        return regula(a,b,f)
    else:
        return None, None
    
# Numerical Integration

def Monte(a,b,N,f):
    """
        Returns the integral of a function using the Monte Carlo method.
    """
    s1 = 0
    s2 = 0
    for i in range(0,N):
        x = a + (b-a)*lcg()
        t = f(x)
        s1 = s1 + t
        s2 = s2 + t**2
    F = (b-a)*s1/N
    sig = s2/N - (s1/N)**2
    return F, sig

def mid(a,b,N,f):
    """
        Returns the integral of a function using the Mid-point method.
    """
    h = (b - a)/N
    t = a + h/2
    s = 0
    for i in range(0, N):
        s = s + h*f(t+h*i)
    return s

def trap(a,b,N,f):
    """
        Returns the integral of a function using the Trapezoidal method.
    """
    h = (b - a)/N
    s = (h/2)*f(b) + (h/2)*f(a)
    for i in range(1,N):
        s = s + h*f(a + i*h)
    return s

def simp(a,b,N,f):
    """
        Returns the integral of a function using the Simpson method.
    """
    h = (b - a)/N
    t = a
    s = (h/3)*f(b) + (h/3)*f(a)
    for i in range(1, N):
        if i % 2 != 0:
            w = 4*h/3
        elif i % 2 == 0:
            w = 2*h/3
        s = s + w*f(a + h*i)
    return s

# ODE Solvers

def ForEuler(x0,y0,xn,N,dy):
    """
        Returns the solution to an ODE using the Forward-Euler method.
    """
    h = abs(xn-x0)/N
    y = y0
    Y = [y0]
    for i in range(N):
        y = y + h*dy(x0 + i*h)
        Y.append(y)
    return Y

def Pred(x0,y0,xn,N,dy):
    """
        Returns the solution to an ODE using the Predictor-Corrector method.
    """
    h = abs(xn-x0)/N
    yp = y0
    yc = y0
    x1 = x0
    y = [y0]
    for i in range(N):
        yp = yp + h*dy(x1)
        x2 = x1 + h
        yc = yc + h*(dy(x1) + dy(x2))/2
        x1 = x2
        y.append(yc)
    return y

def RK4(x0,v0,N,a,b,dxdt,dvdt):
    """
        Returns the solution to an ODE using the RK4 (Runge-Kutta) method.
    """
    h = (b-a)/N
    v = v0
    x = x0
    t = 0
    X = [x0]
    V = [v0]
    T = [t]
    for i in range(N):
        k1x = h*dxdt(x,v,t)
        k1v = h*dvdt(x,v,t)
        
        k2x = h*dxdt(x,v + k1v/2, t)
        k2v = h*dvdt(x + k1x/2, v, t)
        
        k3x = h*dxdt(x,v + k2v/2, t)
        k3v = h*dvdt(x + k2x/2, v, t)
        
        k4x = h*dxdt(x,v + k3v, t)
        k4v = h*dvdt(x + k3x, v, t)
        
        x = x + (k1x + 2*k2x + 2*k3x + k4x)/6
        v = v + (k1v + 2*k2v + 2*k3v + k4v)/6
        t = t + h
        
        X.append(x)
        V.append(v)
        T.append(t)
        
    return X,V,T

def heat(lx,lt):
    """
        Returns the solution to the heat equation and plots the same.
    """
    nx = 20
    nt = 1000
    v1 = []
    v2 = []
    v0 = []
    hx = lx/nx
    ht = lt/nt
    a = ht/hx**2
    x = []
    for i in range(nx+1):
        x.append(hx*i)
    for i in range(len(x)):
        if x[i] == lx/2:
            v0.append(573)
        else:
            v0.append(0)
    v1 = v0[:]
    for k in range(600):
        c = k%100
        v2 = []
        for i in range(nx+1):
            if i == 0:
                v2.append((1-2*a)*v1[i] + a*v1[i+1])
            elif i == nx:
                v2.append((1-2*a)*v1[i] + a*v1[i-1])
            else:
                v2.append((1-2*a)*v1[i] + a*(v1[i-1] + v1[i+1]))
        v1 = v2[:]
        if c == 0:
            plt.plot(x,v2,label=str(k)+' s')

    plt.xlabel('Position, x')
    plt.ylabel('Temperature, T')
    plt.legend()
    plt.savefig('Heat Plot.png')
    plt.show()
    
def pde(x0,v0,N,a,b,dxdt,dvdt):
    """
        Returns the solution to PDE using the RK4 method.
    """
    h = (b-a)/N
    v = v0
    x = x0
    t = 0
    X = [x0]
    V = [v0]
    T = [t]
    for i in range(N):
        k1x = h*dxdt(x,v,t)
        k1v = h*dvdt(x,v,t)
        
        k2x = h*dxdt(x,v + k1v/2, t)
        k2v = h*dvdt(x + k1x/2, v, t)
        
        k3x = h*dxdt(x,v + k2v/2, t)
        k3v = h*dvdt(x + k2x/2, v, t)
        
        k4x = h*dxdt(x,v + k3v, t)
        k4v = h*dvdt(x + k3x, v, t)
        
        x = x + (k1x + 2*k2x + 2*k3x + k4x)/6
        v = v + (k1v + 2*k2v + 2*k3v + k4v)/6
        t = t + h
        
        X.append(x)
        V.append(v)
        T.append(t)
        
    return X,V,T

# Finding the dominant eigenvalue and the corresponding eigenvector

def Evalue(A, x):
    """
        Returns the dominant eigenvalue and the corresponding eigenvector of a given matrix.
        It takes the matrix A and a guess vector x (matrix input).
    """
    k = 2
    xk1 = x[:]
    xk2 = mat_mult(A,xk1)
    e1 = 0
    e2 = dot(xk2,x)/dot(xk1,x)
    while abs(e1-e2) > 10**(-3) and k < 100:
        e1 = e2
        xk1 = xk2
        xk2 = mat_mult(A,xk1)
        e2 = dot(xk2,x)/dot(xk1,x)
        k = k + 1
    ev = []
    n = 0
    for i in range(len(xk2)):
        n = n + xk2[i][0]**2
    n = mt.sqrt(n)
    for i in range(len(xk2)):
        ev.append(xk2[i][0]/n)
        ev[i] = mt.trunc(ev[i]*10**3)/10**3
    e = mt.trunc(e2*10**3)/10**3
    print('Number of iterations =',k)
    return e,ev