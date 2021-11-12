
################################    A Library of Functions      ##################################

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#simple function which displays a matrix

def matrixDisplay(M):
    for i in range(len(M)):
        for j in range(len(M[i])):
            print((M[i][j]), end = " ")
        print()


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

#matrix product

def matrixProduct(L, M):
    if len(L[0]) != len(M): #ensuring the plausiblity
        print("Matrix multiplication not possible.")
    else:
        print("Multiplying the two matrices: ")
        P=[[0 for i in range(len(M[0]))] for j in range(len(L))] #initializing empty matrix
        for i in range(len(L)): #iterating rows
            for j in range(len(M[0])): #iterating columns
                for k in range(len(M)): #iterating elements and substituing them
                    P[i][j] = P[i][j] + (L[i][k] * M[k][j])
        matrixDisplay(P)

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#the gauss-jordan elimination code

def gaussj(a, b):
    n = len(b) #defining the range through which the loops will run
    for k in range(n): #loop to index pivot rows and eliminated columns
	#partial pivoting
        if abs(a[k][k]) < 1.0e-12: 
            for i in range(k+1, n): 
                if abs(a[i][k]) > abs(a[k][k]):
                    for j in range(k, n): 
                        a[k][j], a[i][j] = a[i][j], a[k][j] #swapping of rows
                    b[k], b[i] = b[i], b[k] 
                    break
	#division of the pivot row 
        pivot = a[k][k]
        if pivot == 0:
            print("There is no unique solution to this system of equations.")
            return
        for j in range(k, n): #index of columns of the pivot row
            a[k][j] /= pivot
        b[k] /= pivot
	#elimination loop
        for i in range(n): #index of subtracted rows
            if i == k or a[i][k] == 0: continue
            factor = a[i][k]
            for j in range(k, n): #index of columns for subtraction
                a[i][j] -= factor * a[k][j]
            b[i] -= factor * b[k]
    print(b)

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#calculation of determinant using gauss-jordan elimination


def determinant(a):
    n = len(a) #defining the range through which the loops will run
    if n != len(a[0]): #checking if determinant is possible to calculate
        print("The matrix must be a square matrix.")
    else:
        s = 0
	#code to obtain row echelon matrix using partial pivoting
        for k in range(n-1):
            if abs(a[k][k]) < 1.0e-12:
                for i in range(k+1, n):
                    if abs(a[i][k]) > abs(a[k][k]):
                        for j in range(k, n):
                            a[k][j], a[i][j] = a[i][j], a[k][j] #swapping
                            s = s + 1 #counting the number of swaps happened
            for i in range(k+1, n):
                if a[i][k] == 0: continue
                factor = a[i][k]/a[k][k]
                for j in range(k, n):
                    a[i][j] = a[i][j] - factor * a[k][j]
        d = 1
        for i in range(len(a)):
            d = d * a[i][i] #enlisting the diagonal elements
        d = d*(-1)**s
        print(d)

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#calculating inverse


def inverse(a):
    n = len(a) #defining the range through which loops will run
    #constructing the n X 2n augmented matrix
    P = [[0.0 for i in range(len(a))] for j in range(len(a))]
    for i in range(n):
        for j in range(n):
            P[j][j] = 1.0
    for i in range(len(a)):
        a[i].extend(P[i])
    #main loop for gaussian elimination begins here
    for k in range(n):
        if abs(a[k][k]) < 1.0e-12:
            for i in range(k+1, n):
                if abs(a[i][k]) > abs(a[k][k]):
                    for j in range(k, 2*n):
                        a[k][j], a[i][j] = a[i][j], a[k][j] #swapping of rows
                    break
        pivot = a[k][k] #defining the pivot
        if pivot == 0: #checking if matrix is invertible
            print("This matrix is not invertible.")
            return
        else:
            for j in range(k, 2*n): #index of columns of the pivot row
                a[k][j] /= pivot
            for i in range(n): #index the subtracted rows
                if i == k or a[i][k] == 0: continue
                factor = a[i][k]
                for j in range(k, 2*n): #index the columns for subtraction
                    a[i][j] -= factor * a[k][j]
    for i in range(len(a)): #displaying the matrix
        for j in range(n, len(a[0])):
            print("{:.2f}".format(a[i][j]), end = " ") #printing upto 2 places in decimal. 
        print()


#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#defining a function for forward substituion

def fwdSub(L, b):
    y = [0 for i in range(len(b))]
    for i in range(len(b)):
        sumj = 0
        for j in range(i):
            sumj = sumj + L[i][j]*y[j]
        y[i] = (b[i] -sumj)/L[i][i]
    return y

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#defining a function for backward substitution


def bwdSub(U, y):
    n = len(y)
    x = [0 for i in range(len(y))]
    for i in range(n-1, -1, -1):
        sumj = 0
        for j in range(i+1, n):
            sumj = sumj + U[i][j] * x[j]
        x[i] = (y[i] - sumj)/U[i][i]
    return x


#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#defining doolittle algorithm


def dlittle(A):
    n = len(A)
    L = [[0 for i in range(n)] for j in range(n)]
    U = [[0 for i in range(n)] for j in range(n)]
    for z in range(n):
        L[z][z] = 1
        c1 = 0
        c2 = 0
        c3 = 0
        for p in range(z):
            c1 = c1 + L[z][p]*U[p][z]
        U[z][z] = (A[z][z] - c1)
        for i in range(z+1, n):
            for p in range(z):
                c2 = c2 + L[z][p]*U[p][i]
            U[z][i] = (A[z][i] - c2)
        for k in range(z+1, n):
            for p in range(z):
                c3 = c3 + L[k][p]*U[p][z]
            L[k][z] = (A[k][z] - c3)/U[z][z]
    return (L, U)

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#defining crout algorithm


def crout(A):
    n = len(A)
    L = [[0 for i in range(n)] for j in range(n)]
    U = [[0 for i in range(n)] for j in range(n)]
    for z in range(n):
        U[z][z] = 1
        for j in range(z, n):
            tempL = A[j][z] 
            for k in range(z):
                tempL -= L[j][k]*U[k][z]
            L[j][z] = tempL
        for j in range(z+1, n):
            tempU = A[z][j]
            for k in range(z):
                tempU -= L[z][k]*U[k][j]
            U[z][j] = tempU/L[z][z]
    return (L, U)

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#defining a solver function


def solver(A, b, algo):
    L, U = algo(A)
    print("L = " + str(L) + "\n")
    print("U = " + str(U) + "\n")
    y = fwdSub(L, b)
    x = bwdSub(U, y)
    return x

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#defining the cholesky algorithm


def cholesky(a):
    n = len(a)
    L = [[0 for i in range(n)] for j in range(n)]
    for j in range(n):
        for i in range(j, n):
            if i == j:
                sumj = 0
                for k in range(j):
                    sumj = sumj + (L[i][k]**2)
                L[i][j] = (a[i][j] - sumj)**(1/2)
            else:
                sumk = 0
                for k in range(j):
                    sumk = sumk + (L[i][k]*L[j][k])
                L[i][j] = (a[i][j] - sumk)/L[j][j]
    return L

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#defining the solver function for cholesky algorithm 


def choleskySolver(L, U, b):
    n = len(L)
    y = [0 for i in range(n)]
    x = [0 for i in range(n)]
    for i in range(n):
        sumj = 0
        for j in range(i):
            sumj += L[i][j]*y[j]
        y[i] = (b[i]-sumj)/L[i][i]
    for i in range(n-1, -1, -1):
        sumj = 0
        for j in range(i+1, n):
            sumj += U[i][j]*x[j]
        x[i] = (y[i]-sumj)/U[i][i]
    return x

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
     


#defining a function to transpose the matrices

def transpose(a):
    n = len(a)
    astar = [[0 for i in range(n)] for j in range(n)]
    for i in range(len(a)):
        for j in range(len(a[0])):
            astar[j][i] = a[i][j]
    return astar

#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

#defining the LU decomposition function and related complementary functions for question 2

def inverseLU(M,I): #function to call inverse using LU
    if M[1][1] == 0 and M[0][1] != 0:
        swapRows(M, 0,1,4) #if diagonal element is 0, swaps to prevent 0 determinant
    LUDecomp(M)

    L = LUDecomp.L #making the given matrix into LU form
    U = LUDecomp.U 
    return fwdBwdSub(L,U,I)


##################################################################################################
##################################################################################################


def swapRows(Ab, old, new_r, cols):
    temp = []       #temp list to store old list

    for c in range (0, int(cols)):
        temp.append(Ab[int(old)][c])
        Ab[int(old)][c] = Ab[int(new_r)][c]     #swapping values
        Ab[int(new_r)][c] = temp[c]


##################################################################################################
##################################################################################################


def LUDecomp(M):
    PartialPivot(M, 0, len(M), len(M[0]))  #partial pivoting given matrix
    n = len(M)        
    lower = [[0 for x in range(n)]
             for y in range(n)]
    upper = [[0 for x in range(n)]
             for y in range(n)]
 
    # Decomposing matrix into Upper
    # and Lower triangular matrix
    for i in range(n):
 
        # Upper Triangular
        for j in range(i, n):
 
            # Summation of L(i, j) * U(j, k)
            sum = 0
            for k in range(i):
                sum += (lower[i][k] * upper[k][j])
 
            # Evaluating U(i, k)
            upper[i][j] =  M[i][j] - sum
 
        # Lower Triangular
        for k in range(i, n):
            if (i == k):
                lower[i][i] = 1  # Diagonal as 1
            else:
 
                # Summation of L(k, j) * U(j, i)
                sum = 0
                for j in range(i):
                    sum += (lower[k][j] * upper[j][i])
 
                # Evaluating L(k, i)
                lower[k][i] = ((M[k][i] - sum) /
                                  upper[i][i])
    print("Lower Triangular")

    # Displaying the result :
    for i in range(n):
        print(lower[i])

    print ("Upper Triangular")
    for i in range(n):
        print(upper[i])
    LUDecomp.L = lower
    LUDecomp.U = upper


##################################################################################################
##################################################################################################


def fwdBwdSub(L, U, b):
    y = [[0 for c in range(len(b[0]))] for r in range(len(b))]
    for i in range(len(b)):
        for k in range (len(b[0])): #looping over the coloumns to calculate y
            y[i][k] = b[i][k]
            for j in range(i):
                y[i][k]=y[i][k]-(L[i][j]*y[j][k]) #formula for y
    
            y[i][k] = y[i][k]/L[i][i] 

    n = len(y)


    x = [[0,0,0,0] for r in range(len(b))]
    if U[n-1][n-1] == 0: #checking if diagonal elements are zero
        raise ValueError

    for i in range(n-1, -1, -1):
        for k in range (len(b[0])): #iterating over coloumns to calculate x
            x[i][k] = y[i][k]
            for j in range(i+1,n):
                x[i][k] = x[i][k] -(U[i][j]*x[j][k]) #formula for x
            x[i][k] = x[i][k]/U[i][i]    
    
    print ("The inverse of the given matrix is " )
    for i in x:
        print (i)
    return(x)


##################################################################################################
##################################################################################################


def PartialPivot(Ab, m, rows, cols):
    global n,swapnumber          #global variable to store how many swap are done
    n = 0
    swapnumber = 0
    pivot = Ab[int(m)][int(m)]          #starting pivot of matrix
    for i in range (int(rows)):         
        if pivot < Ab[int(i)][int(m)]:  #checking with other elements of the same coloumn
            pivot = Ab[int(i)][int(m)]
            n += 1
            swapnumber = i
    if swapnumber != 0:
        swapRows(Ab, m, swapnumber, cols)    #swapping if condition satisfies
            
    if int(pivot) == 0:
        print ("No unique solution")   #if pivot is 0 at end it returns no solution
        return None



##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################



 
#following modules are required to execute some essential things asked in the question
import matplotlib.pyplot as plt
import math

##################################################################################################
##################################################################################################



#defining the bracketing function for bisection method
#this function also prints the plot and a table of f(x_i) vs x_i


def rootBisec(y, a, b):

    #this part of the code is to print the plot and the table

    if y(a)*y(b)<0: #no need for bracketing
        iterx, iterCount, xiCount = bisection(y, a, b, 1.0e-5) #execute the bisection function. the tolerance e is defined as 1e-05 as asked in the question
        x1 = iterx
        y1 = iterCount
        y11 = xiCount
        plt.plot(x1, y1)
        plt.show()

        print("i", end=" ")
        print("x_i")
        for i in range(len(x1)):
            print(x1[i], end=" ")
            print(y11[i])
    
    #this part of the code executes the bracketing algorithm 

    else: #executing bracketing with beta = 1.1
        if abs(y(a))<abs(y(b)):
            a = a - 1.1*(b-a)
            rootBisec(y, a, b) #call rootBisec again and this process continues until the bracketing is complete
        if abs(y(a))>abs(y(b)):
            b = b + 1.1*(b-a)
            rootBisec(y, a, b) #call rootBisec again and this process continues until the bracketing is complete

#note that this code is not domain sensitive. it will work with most functions but with
#functions like the natural logarithm which have very strict domain of definitions, the
#input values while performing bracketing will quickly go out of domain, resulting in a
#math domain error from Python3. so kindly choose your guesses with great deliberation.

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#the actual function which executes (or finds the root) via the bisection method


def bisection(f, a, b, e): #f is the function, a and b represent the interval and e is the tolerance
    i = 0 #this is defined so as to make it easy to construct the table later
    x = 1 #a dummy value for the initial value of root
    iterCount = [] #a list which stores the f(x_i)'s
    xiCount = []
    iterx = [] #a list which stores the x_i's

    #the main code

    condition = True

    #there are two conditions of the code related to the position of the initial guesses and the nature of curve
    #the code will change slightly for whether the function is increasing or decreasing in the chosen interval
    #these if statements take care of that

    #if function is increasing in chosen interval

    if f(a)<f(b):
        while condition:
            g = f(x)
            x = (a+b)/2
            if f(x)<0:
                a = x
            else:
                b = x
            iterCount.append(f(x))
            xiCount.append(x)
            i = i + 1
            iterx.append(i)
            condition = abs(f(x)-g)>e
    
    #if function is decreasing in chosen interval

    if f(a)>f(b):
        while condition:
            g = f(x)
            x = (a+b)/2
            if f(x)<0:
                b = x
            else:
                a = x
            iterCount.append(f(x))
            xiCount.append(x)
            i = i + 1
            iterx.append(i)
            condition = abs(f(x)-g)>e       
    print("Required root is: ", x) #printing the root
    return iterx, iterCount, xiCount #return the values for plotting and tabulating purposes


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#the code for the Regula Falsi or false position method


def regfalsi(f,a,b,e,n): #f is the function, a and be represent the interval, e is the tolerance and 
                         #n is an upper limit on the number of iterations
    xh = 0 #initial dummy value for the root
    iterCount = [] #a list which stores the f(x_i)'s
    xiCount = []

    #the main code starts here

    for fal in range(1,n+1):
        xh = b - (b-a)/(f(b)-f(a))*f(b)
        iterCount.append(f(xh))
        xiCount.append(xh)
        if abs(f(xh)) < e: break
        elif f(a)*f(xh)<0:
            b = xh
        else:
            a = xh
    print("Required root is: ", xh) #printing the root
    return iterCount, xiCount #return the values for plotting and tabulating purposes


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#this code is analogous to rootBisec. plots and tabulate the asked details as well


def rootRegFalsi(y, a, b):

    #this part of the code is to print the plot and the table

    if y(a)*y(b)<0: #no bracketing needed
        iterCount, xiCount = regfalsi(y, a, b, 1.0e-5, 100) #y is function, a and b represent the interval. 
                                                   #e is the tolerance (value specified according to the question)
                                                   #maximum iterations are taken to be 100
        x2 = list(range(1, len(iterCount)+1))
        y2 = iterCount
        y22 = xiCount
        plt.plot(x2, y2)
        plt.show()
        print("i", end=" ")
        print("x_i")
        for i in range(len(x2)):
            print(x2[i], end=" ")
            print(y22[i])

    #this part of the code executes the bracketing algorithm
    
    else:
        if abs(y(a))<abs(y(b)): 
            a = a - 1.1*(b-a) #executing bracketing with beta = 1.1
            rootRegFalsi(y, a, b) #call rootBisec again and this process continues until the bracketing is complete
        if abs(y(a))>abs(y(b)):
            b = b + 1.1*(b-a)
            rootRegFalsi(y, a, b) #call rootBisec again and this process continues until the bracketing is complete


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#the code for Newton-Raphson method to find roots

def newRaph(f, df, x0, e, n): #f is the function, df is the first derivative of the function
                              #x0 is the dummy valriable for root
                              #e is the tolerance (taken as 1.0e-5 as required in the question)
                              #n is the maximum iterations 
    iterCount = [f(x0)] #this enlists the f(x_i)'s for plotting and tabulating purposes
    xiCount = []

    #the main algorith for Newton-Rhapson here

    for i in range(n):
        xnew = x0 - f(x0)/df(x0)
        iterCount.append(f(xnew))
        xiCount.append(xnew)
        if abs(xnew - x0)<e: break
        x0 = xnew
    
    #the code for plotting and tabulating


    x3 = list(range(1, len(iterCount)+1))
    y3 = iterCount
    y33 = xiCount
    plt.plot(x3, y3)
    plt.show()
    print("i", end=" ")
    print("x_i")
    for i in range(len(x3)-1):
        print(x3[i], end=" ")
        print(y33[i])
    return xnew, i #returns the value for printing the results later


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


#this is a code to evaluate the co-efficients of the differentials of the polynomial
#the format in which polynomial p must be used is written below


#p,dp,ddp = evalPoly(a,x).
#p = a[0] + a[1]*x + a[2]*xˆ2 +...+ a[n]*xˆn
#with its derivatives dp = p’ and ddp = p’’
#at x.


def evalPoly(a, x):
    n = len(a) - 1
    p = a[n]
    dp = 0.0
    ddp = 0.0
    for i in range(1, n+1):
        ddp = ddp*x +2.0*dp
        dp = dp*x + p
        p = p*x + a[n-i]
    return p, dp, ddp

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################



#this is the main algorithm code for Laguerre's method

#roots = polyRoots(a).
#again the format of the polynomial should be: a[0] + a[1]*x + a[2]*xˆ2 +...+ a[n]*xˆn = 0.
#The roots are returned in the array ’roots’

def polyRoots(guess, a, tol=1.0e-12): #guess is the input from the user 
                                      #a is the list which contains the co-efficients of polynomial p
                                      #tol is tolerance, taken to be 1.0e-12 (arbitrarily)
    
    #main code for Laguerre's method begins here

    def laguerre(a, tol):
        x = guess #starting value is used as an input from the user
        n = len(a) - 1
        for i in range(30):
            p, dp, ddp = evalPoly(a,x) #calling the evalPoly to generate p, dp and ddp
            if abs(p) < tol: return x
            g = dp/p
            h = g*g - ddp/p
            f = math.sqrt((n-1)*(n*h - g*g))
            if abs(g+f) > abs(g-f): dx = n/(g+f)
            else: dx = n/(g-f)
            x = x - dx
            if abs(dx) < tol: return x
        print("Too many iterations") #exception in case of the guess is not enough to find the roots

    #the following code is to deflate the polynomial using the synthetic division method

    def deflPoly(a, root):
        n = len(a) - 1
        b = [0.0]*n
        b[n-1] = a[n]
        for i in range(n-2, -1, -1):
            b[i] = a[i+1] + root*b[i+1]
        return b

    #this final part of function returns the necessary roots neatly in a list

    n = len(a) - 1
    roots = [0.0 for i in range(n)]
    for i in range(n):
        x = laguerre(a, tol)
        roots[i] = x
        a = deflPoly(a, x)
    return roots


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


def midpoint(a, b, N, y): #a is the lower limit, b is the upper limit
                          #N is the number of divisions made in the limit
                          #y is the mathematical function to be integrated   
    h = (b-a)/N
    mids = []
    for i in range(1, N+1):
        mids.append((2*a + (2*i -1)*h)/2)
    I = 0
    for j in range(len(mids)):
        I = I + h*y(mids[j])
    return I


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


def trapezoidal(a, b, N, y): #a is the lower limit, b is the upper limit
                             #N is the number of divisions made in the limit
                             #y is the mathematical function to be integrated
    h = (b-a)/N
    ends = []
    for i in range(N+1):
        ends.append(a + i*h)
    T = []
    for j in range(1, len(ends)):
        T.append((h/2)*(y(ends[j-1])+y(ends[j])))
    trap = sum(T)
    return trap

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################


def simpson(a, b, N, y):#a is the lower limit, b is the upper limit
                        #N is the number of divisions made in the limit
                        #y is the mathematical function to be integrated
    
    #check if N is even or not, convert to even if N odd

    if (N % 2) != 0:
        N = N + 1

    #simpson algorithm starts here


    h = (b-a)/N
    ends = []
    for i in range(N+1):
        ends.append(a + i*h)
    S = []
    for j in range(0, N-1, 2):
        S.append((h/3)*(y(ends[j])+(4*y(ends[j+1]))+y(ends[j+2])))
    simp = sum(S)
    return simp


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

<<<<<<< HEAD

=======
>>>>>>> c6eeebd689ba4093dc42c20ee786c9c2b882205e
##################################################################################################

import random

##################################################################################################

<<<<<<< HEAD
def monteCarlo(a, b, N, s, iter, y): #a and b are the limits of the integration
                                     #N is the number of the divisions made within the limits
                                     #s is the step-size of increasing N
                                     #iter is the number of times N will be increased
                                     #y is the mathematical function to be integrated

=======
def monteCarlo(a, b, N, s, iter, y):
>>>>>>> c6eeebd689ba4093dc42c20ee786c9c2b882205e
    Flist = []
    xlist = list(range(N, iter*s, s))
    while N < (iter*s):
        vars = []
        for i in range(N):
            vars.append(random.uniform(a, b))
        f = []
        f2 = []
        for j in range(len(vars)):
            f.append(y(vars[j]))
        for j in range(len(vars)):
            f2.append((y(vars[j]))**2)
        sumf = sum(f)
        sumf2 = sum(f2)
        F = ((b-a)/N)*sumf
        Flist.append(F)
        sigmaf = math.sqrt((1/N)*sumf2 - ((1/N)*sumf)**2)
        N = N + s
    plt.plot(xlist, Flist)
<<<<<<< HEAD
    plt.ylabel("\u03C0")
    plt.xlabel("N")
=======
>>>>>>> c6eeebd689ba4093dc42c20ee786c9c2b882205e
    plt.show


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################