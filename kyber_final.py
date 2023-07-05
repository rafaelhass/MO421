import numpy as np
import math

n = 256 # 256
q = 3329 # 3329
k = 2
n1 = 3
n2 = 2
cofdiv = np.zeros(n+1, dtype=int)
cofdiv[0] = 1
cofdiv[-1] = 1

def createRandomMatrixOfArraysofIntLimited(dimX: int, dimY: int, dimArrays: int, minValue: int, maxValue: int):
    matrix = np.empty((dimX, dimY), dtype=np.ndarray) # Create the matrix
    # Generate random values and assign them to each array
    for i in range(dimX):
        for j in range(dimY):
            random_values = np.random.randint(low=minValue, high=maxValue, size=dimArrays) 
            matrix[i, j] = random_values
    return matrix

def setPolynomialRing(A: np.ndarray):
    rows = A.shape[0]
    columns = A.shape[1]
    matrix = np.zeros((rows, columns), dtype=np.ndarray)

    for i in range(rows):
        for j in range(columns):

            div, d = np.polynomial.polynomial.polydiv(A[i][j], cofdiv)
            matrix[i][j] = d % q

    return matrix

def sum_arrays(A: np.ndarray,B: np.ndarray):

    if A.shape != B.shape:
        raise Exception("Dimensions incompatible:" + str( A.shape)+'x'+str(B.shape))
    rows = A.shape[0]
    columns = A.shape[1]
    matrix = np.zeros((rows, columns), dtype=np.ndarray)
    for i in range(rows):
        for j in range(columns):
            poly_sum = np.polynomial.polynomial.polyadd(A[i][j], B[i][j])
            matrix[i][j] = poly_sum

    return(setPolynomialRing(matrix))

def mult_matrix_polinomyals(A: np.ndarray,B: np.ndarray):

    rows = A.shape[0]
    columns = B.shape[1]
    matrix = np.zeros((rows, columns), dtype=np.ndarray)
    if A.shape[0] % B.shape[1] != 0:
        raise Exception("Dimensions incompatible:" + str( A.shape[0])+'x'+str(B.shape[1]))
    
    # Perform matrix multiplication
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            for k in range(B.shape[0]):
                resultMult = np.polynomial.polynomial.polymul(A[i][k], B[k][j]) # Line A x Column B
                matrix[i, j] = np.polynomial.polynomial.polyadd(matrix[i,j],resultMult)
    
    return(setPolynomialRing(matrix))

def generateKey():

    s = createRandomMatrixOfArraysofIntLimited(k,1,n,-n1,n1+1)
    A = createRandomMatrixOfArraysofIntLimited(k,k,n,-q,q)
    e = createRandomMatrixOfArraysofIntLimited(k,1,n,-n1,n1+1)

    multA_S= mult_matrix_polinomyals(A,s)
    t = sum_arrays(multA_S, e)
    return ((A,t),s)

def encription(arrayMessageB, SK):

    A = SK[0]
    t = SK[1]
    r = createRandomMatrixOfArraysofIntLimited(k,1,n,-n1,n1+1)
    e1 = createRandomMatrixOfArraysofIntLimited(k,1,n,-n2,n2+1)
    e2 = createRandomMatrixOfArraysofIntLimited(1,1,n,-n2,n2+1)
    m = math.ceil(q/2)*arrayMessageB

    #u = A.T + e1
    mult_At_r = mult_matrix_polinomyals(A.T,r)
    u = sum_arrays(mult_At_r, e1)

    #v = t_t*r + e2 + m

    mult_tt_r = mult_matrix_polinomyals(t.T,r)
    
    sumE2M = sum_arrays(e2, m)
    v = sum_arrays(mult_tt_r,sumE2M)

    return (u,v)

def decription(s, cm):
    u = cm[0]
    v = cm[1]

    #mn = v - s_t*u
    mult_st_u = mult_matrix_polinomyals(s.T,u)

    mn = sum_arrays(v, -mult_st_u)

    # print('mn: ' + str(mn))

    resultado = []
    for item in mn[0][0]:
        if abs(item - math.ceil(q/2)) < math.ceil(q/4):
            resultado.append(1)
        else:
            resultado.append(0)

    while len(resultado) < n:
        resultado.append(0)

    return [resultado]

def KEM_KeyGen():
    return generateKey()

def KEM_Encrypt(pk):

    m = createRandomMatrixOfArraysofIntLimited(1,1,n,0,2)
    ct = encription(m, pk)
    return ct,m

def KEM_Decrypt(sk, ct):
    return decription(sk, ct)

##################################### funcao principal  #####################################

def encriptMain():
    pk,sk = generateKey()
    m = createRandomMatrixOfArraysofIntLimited(1,1,n,0,2)
    ct = encription(m, pk)
    mb = decription(sk, ct)

    if np.polynomial.Polynomial(mb[0]) == np.polynomial.Polynomial(m[0][0]):
        return True
    
def KyberKem():
    pk,sk = KEM_KeyGen()
    ct,ss = KEM_Encrypt(pk)
    ss_1 = KEM_Decrypt(sk, ct)

    if np.polynomial.Polynomial(ss_1[0]) == np.polynomial.Polynomial(ss[0][0]):
        return True

if __name__ == "__main__":

    # KyberKem()
    for i in range(1000):
        if not(KyberKem()):
            print('deu ruim')

    # encriptMain()
    for i in range(1000):
        if not(encriptMain()):
            print('deu ruim')
