
#It's Linear Algebra time
import random as random
import math as math

class DimensionError(Exception):
    def __init__(self, reason):
        self.reason = reason
    def __str__(self):
        return repr(self.reason)
    
class Vector(object):
    def __init__(self, argv):
        for arg in argv:
            argtype = type(arg)
            if (argtype != int) and (argtype != float):
                raise TypeError("a Vector must only contain numbers!")
        self.values = list(argv)

    def __str__(self):
        return f"{tuple(self.values)}"
    
    def __getitem__(self, item):
        if (type(item) != int):
            return 0
        else:
            return self.values[item]

    def __len__(self):
        return len(self.values)

    @property
    def rlen(self):
        zeros = 0
        for entry in self.values:
            if (entry != 0):
                return len(self) - zeros
            else:
                zeros += 1

    @property
    def llen(self):
        zeros = 0
        for entry in self.values:
            if (entry != 0):
                return zeros
            else:
                zeros += 1
    
    @property
    def dim(self):
        contents = self.values
        rows = len(contents)
        if type(self) == Vector:
            columns = 1
        else:
            columns = len(contents[0].values)
        return (rows, columns)

    def dot_compatible(self, other):
        if (type(other) == type(self)):
            return (self.dim == other.dim)
        else:
            return False
            
    def __add__(self, other):
        if (type(other) != type(self)):
            raise TypeError(f"{other} not a Vector!")
        elif (self.dim != other.dim):
            raise DimensionError("Incompatible Dimensions!")
        else:
            tempv = []
            v1, u1 = self.values, other.values
            for index in range(len(v1)):
                tempv.append(v1[index] + u1[index])
            return Vector(tempv)

    def __sub__(self, other):
        return self + other*(-1)
    
    def dot(self, other):
        if not self.dot_compatible(other):
            raise TypeError("Incompatible Objects!")
        else:
            v1, u1 = self.values, other.values
            output = 0
            for index in range(len(v1)):
                output += (v1[index] * u1[index])
            return output

    def __mul__(self, other):
        if (type(other) == int) or (type(other) == float):
            if (other == 0):
                output = [0]*self.dim[0]
                return Vector(output)
            tempv, vector = [], self.values
            for number in vector:
                tempv.append(number*other)
            return Vector(tempv)
        #if i wish to define vector multiplication, insert code here

    @property
    def is_zero(self):
        for item in self:
            if (item != 0):
                return False
        return True

class Matrix(Vector):
    def __init__(self, argv):
        if (set([type(X) for X in argv])=={list}):
            self.values = [Vector(X) for X in argv]
        else:
            for arg in argv:
                argtype = type(arg)
                if (argtype != (Vector)):
                    raise TypeError("a Matrix must only contain Vectors!")
            self.values = list(argv)

    def __str__(self):
        output = ""
        for vector in self.values:
            output += " "
            output += repr(tuple(vector.values))
            output += ",\n"
        output = output[1:-2]
        output = "[" + output + "]"
        return output

    def __add__(self, other):
        if (type(other) != type(self)):
            raise TypeError(f"{other} not a Matrix!")
        elif (self.dim != other.dim):
            raise DimensionError("Incompatible Dimensions!")
        else:
            tempM = []
            m1, m2 = self.values, other.values
            for index in range(len(m1)):
                tempM.append(m1[index] + m2[index])
            return Matrix(tempM)
        
    def __sub__(self, other):
        return self + other*(-1)
    
    def transpose(self):
        other_matrix, tempOM, O_Rows, O_Cols = self.values, [], self.dim[0], self.dim[1]
        for column in range(O_Cols):
            temp_col = []
            for row in range(O_Rows):
                temp_col.append(other_matrix[row].values[column])
            tempOM.append(Vector(temp_col))
        return Matrix(tempOM)
        
    def __mul__(self, other):
        if (type(other) == int) or (type(other) == float):
            matrix = self.values
            output = Matrix([vector*other for vector in matrix])
            return output
        elif (type(other) == Vector):
            if (self.dim[1] != other.dim[0]):
                raise DimensionError("Incompatible Dimensions!")
            else:
                tempV, matrix, vector = [], self.values, other
                for row in matrix:
                    tempV.append(row.dot(vector))
                return Vector(tempV)
        elif (type(other) == Matrix):
            if (self.dim[1] != other.dim[0]):
                raise DimensionError("Incompatible Dimensions!")
            else:
                #columnize other matrix
                tempOM = other.transpose().values

                #tempOM is now a collection of column vectors instead
                left_matrix, tempM, L_Rows, L_Cols = self.values, [], self.dim[0], self.dim[1]
                for row_vector in left_matrix:
                    temp_row = []
                    counter = 0
                    for col_vector in tempOM:
                        temp_row.append(row_vector.dot(col_vector))
                    tempM.append(Vector(temp_row))
                return Matrix(tempM)
                        
    def __pow__(self, other):
        if (other == 0):
            rows, columns = self.dim
            tempM = []
            for row in range(rows):
                temp_row = [0]*columns
                temp_row[row] = 1
                tempM.append(Vector(temp_row))
            return Matrix(tempM)
        elif (type(other) == int):
            if (other > 0):
                output = self**0
                for i in range(other):
                    output *= self
                return output
    @property       
    def det(self):
        matrix = self.values
        if (self.dim == (2,2)):
            v1, v2 = matrix
            return (v1[0]*v2[1]) - (v1[1]*v2[0])
        elif (self.dim[0] < self.dim[1]):
            return 0
        else:
            
            coefficients = matrix[0].values
            cofactors = Matrix(matrix[1:]).transpose().values
            index, output = 0, 0
            for k in coefficients:
                cofactor_matrix = Matrix(cofactors[:index] + cofactors[index+1:])
                output += ((-1)**index)*k*cofactor_matrix.det
                index += 1
            return output

    def sort(self, key = lambda x: x[0], reverse = False):
        matrix = self.values
        matrix.sort(key = key, reverse = reverse)
        return Matrix(matrix)
    
    #row echelon form [ref]
    def ref(self):
        rows, cols = self.dim
        matrix = list(map(lambda x: x, self.values))
        matrix.sort(key = lambda x: x.llen)
        output_matrix = [matrix[0]]
        
        for row in range(1, rows):
            #sort at the start of the current iteration, unsure if needed
            #matrix.sort(key = lambda x: x.llen)
            
            #current vector is matrix[row]
            curr_vector = matrix[row]

            #subtract current vector by previous vectors scaled by leading coefficients of previous rows
            for prev_row in range(row):
                previous_vector = output_matrix[prev_row]
                #if the current vector has non-zero component of the previous vector's pivot entry
                if (curr_vector.llen == previous_vector.llen) and not (curr_vector.is_zero or previous_vector.is_zero):
                    #pivot coefficient is: vector[vector.llen]
                    scale_factor = curr_vector[curr_vector.llen]/previous_vector[previous_vector.llen]
                    curr_vector += previous_vector * scale_factor * (-1)
                else:
                    #no need to subtract if component already zero
                    continue
            
            output_matrix.append(curr_vector)
        return Matrix(output_matrix)
    
    #reduced row echelon form [rref] IMCOMPLETE, HAVE ONLY SCALED THINGS DOWN, NEED TO MAKE SURE ALL COLS ONLY CONTAIN AT MOST 1 NON ZERO ENTRY
    def rref(self):
        matrix = self.ref().values
        output_matrix = []
        for vector in matrix:
            #pivot coefficient is: vector[vector.llen]
            coefficient = vector[vector.llen]
            if (coefficient == 0) or (coefficient == 1):
                output_matrix.append(vector)
            else:
                output_matrix.append(vector*(1/coefficient))
        return Matrix(output_matrix)

    #this is my implementation of augmented Matracies
    def __or__(self, other):
        matrix = self.transpose().values
        if (type(other) == Vector):
            matrix.append(other)
        elif (type(other) == Matrix):
            matrix.extend(other.transpose().values)
        else:
            return "Only can augment with other matrix or vectors"
        
        return Matrix(matrix).transpose()
    
    def _(self, other):
        matrix = self.values
        if (type(other) == Vector):
            matrix.append(other)
        elif (type(other) == Matrix):
            matrix.extend(other.values)
        else:
            return "Only can augment with other matrix or vectors"
        return Matrix(matrix)

    
    def adj(self):
        rows, cols = self.dim
        matrix = self.values
        matrix = list(map(lambda x: x.values, matrix))
        #print(f"matrix is \n{matrix}\n")
        output = []
        for row in range(rows):
            new_vector = []
            for col in range(cols):
                #cofactor matrix
                cofactor = matrix[:row] + matrix[row+1:]
                cofactor = list(map(lambda x: x[:col] + x[col+1:], cofactor))
                cofactor = list(map(lambda x: Vector(x), cofactor))
                determinant = Matrix(cofactor).det*((-1)**(row+col))
                new_vector.append(determinant)
            output.append(Vector(new_vector)) 
        return Matrix(output).transpose()

    def rnd(self, digits):
        output = []
        matrix = self.values
        matrix = list(map(lambda x: x.values, matrix))
        for vector in matrix:
            new_vector = list(map(lambda x: round(x, digits), vector))
            output.append(Vector(new_vector))
        return Matrix(output)
    
    def int(self):
        output = []
        matrix = self.values
        matrix = list(map(lambda x: x.values, matrix))
        for vector in matrix:
            itype = type(1)
            new_vector = list(map(lambda x: itype(x), vector))
            output.append(Vector(new_vector))
        return Matrix(output)

    def inverse(self):
        if (self.det == 0):
            return "Not Invertible Matrix"
        elif (self.dim == (2,2)):
            u,v = self
            a,b=u
            c,d=v
            k=1/(a*d-b*c)
            return Matrix([Vector([d, -1*b]), Vector([-1*c,a])])*k
        else:
            scale = 1/self.det
            output = self.adj() * scale
            return output
            
            
            
#test objects to play with
a, b = (1, 2, 3), (-5, 3, 1)
v1, v2 = Vector(a), Vector(b)
u1, u2 = v1 + v2, v1.dot(v2)
M = Matrix((v1, v2, u1))
w1, w2, w3 = Vector((1, 2, 3)), Vector((5, 8, 9)), Vector((12, -1, 2))
W = Matrix((w1, w2, w3))
l1, l2, l3 = Vector((0, 2, 3)), Vector((5, 8, 9)), Vector((0, 0, 2))
L = Matrix((l1,l2,l3))
A=Matrix([Vector([6, 21]),
   Vector([21,122])])
B=Matrix([Vector([1,1,1,1,1]),
        Vector([1,2,4,6,8])])
Y=Vector([3,5,4,9,12])
'''
print(f"sum is {u1}")
print(f"dot product is {u2}")
print(m)
'''

def exp(x):
    output = x*0
    for i in range(100):
        output += (x**i) * (1/(math.factorial(i)))
    return output

