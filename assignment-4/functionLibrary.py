
################################    A Library of Functions      ##################################

##################################################################################################


#simple function which displays a matrix

def matrixDisplay(M):
    for i in range(len(M)):
        for j in range(len(M[i])):
            print((M[i][j]), end = " ")
        print()


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


#defining a solver function


def solver(A, b, algo):
    L, U = algo(A)
    print("L = " + str(L) + "\n")
    print("U = " + str(U) + "\n")
    y = fwdSub(L, b)
    x = bwdSub(U, y)
    return x

#################################################################################################


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


#defining a function to transpose the matrices

def transpose(a):
    n = len(a)
    astar = [[0 for i in range(n)] for j in range(n)]
    for i in range(len(a)):
        for j in range(len(a[0])):
            astar[j][i] = a[i][j]
    return astar

#################################################################################################


def inverseLU(M,I): #function to call inverse using LU
    if M[1][1] == 0 and M[0][1] != 0:
        swapRows(M, 0,1,4) #if diagonal element is 0, swaps to prevent 0 determinant
    LUDecomp(M)

    L = LUDecomp.L #making the given matrix into LU form
    U = LUDecomp.U 
    return fwdBwdSub(L,U,I)

def swapRows(Ab, old, new_r, cols):
    temp = []       #temp list to store old list

    for c in range (0, int(cols)):
        temp.append(Ab[int(old)][c])
        Ab[int(old)][c] = Ab[int(new_r)][c]     #swapping values
        Ab[int(new_r)][c] = temp[c]

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

