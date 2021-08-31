
################################    A Library of Functions      ##################################

##################################################################################################


#simple function which displays a matrix

def matrixDisplay(M):
    for i in range(len(M)):
        for j in range(len(M[i])-1):
            print((M[i][j]), end = " ")
        print()


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
    for i in range(3):
        for j in range(3):
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
            print(a[i][j], end = " ")
        print()
