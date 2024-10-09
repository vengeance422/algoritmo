from fastaReader import fastaReader
import random
import numpy 
import copy
import itertools
from evaluadorBlosum import evaluadorBlosum

class bacteria():
    
    
    def __init__(self, path):
        self.matrix = fastaReader(path)
        self.blosumScore = 0
        self.fitness = 0
        self.interaction =0
        self.NFE = 0
        
    def showGenome(self):
     for seq in self.matrix.seqs:
        print(seq)

    #simplificar el codigo utilizando deepcody
    def clonar(self, path):
        newBacteria = bacteria(path)
        newBacteria.matrix.seqs = copy.deepcopy(self.matrix.seqs)
        return newBacteria


    def tumboNado(self, numGaps):
        
        self.cuadra()
        matrixCopy = copy.deepcopy(self.matrix.seqs)
        """convierto a lista para poder modificar"""
        matrixCopy = matrixCopy.tolist()
        gapRandomNumber = random.randint(0,numGaps)  #numero de gaps a insertar
        for i in range(gapRandomNumber):                    #cilco de gaps 
            seqnum = random.randint(0, len(matrixCopy)-1)   #selecciono secuencia
            pos = random.randint(0, len(matrixCopy[0]))
            part1 = matrixCopy[seqnum][:pos]
            part2 = matrixCopy[seqnum][pos:]
            temp = "-".join([part1, part2])     #inserto gap
            matrixCopy[seqnum] = temp
        matrixCopy = numpy.array(matrixCopy)   #convierto a numpy array de regreso para fijar tama�os
        self.matrix.seqs = matrixCopy
        
        self.cuadra()
        self.limpiaColumnas()
      
        


    def cuadra(self):
        """rellena con gaps las secuencias mas cortas"""
        import numpy
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))
        for i in range(len(seq)):
            if len(seq[i]) < maxLen:
                seq[i] = seq[i] + "-"*(maxLen-len(seq[i]))
        self.matrix.seqs = numpy.array(seq)
        

    """metodo para saber si alguna columna de self.matrix tiene  gap en todos los elementos"""
    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True
    


    """metodo que recorre la matriz y elimina las columnas con gaps en todos los elementos"""
    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1
        
            
        """metodo para eliminar un elemento especifico en cada secuencia"""
    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]





        """metodo para obtener una lista con los elementos de cada columna"""
    def getColumn(self, col):
        column = []
        for i in range(len(self.matrix.seqs)):
            column.append(self.matrix.seqs[i][col])
        return column
    


        """metodo para evaluar columnas"""
    def autoEvalua(self):   
        score = 0
        total_columns = len(self.matrix.seqs[0])
        for i in range(total_columns):  
            column = self.getColumn(i)
        gapCount = column.count("-")
        column = [x for x in column if x != "-"]
        pares = self.obtener_pares_unicos(column)

        for par in pares:
            score += self.evaluador.getScore(par[0], par[1])
        
        # Penalizar gaps de manera proporcional al número de secuencias
        score -= (gapCount * 2) / total_columns
    
        self.blosumScore = score
        self.NFE += 1

        
    #implementación de la libreria itertools
    def obtener_pares_unicos(self, columna):
        return list(itertools.combinations(set(columna), 2))
    
    def update_memory(self, neighbors):
        for neighbor in neighbors:
            if neighbor.fitness > self.best_solution.fitness:
                self.best_solution = neighbor.copy()

    def generate_new_solution(self):
        # ... lógica para generar una nueva solución basada en la memoria
        new_solution = self.best_solution.copy()
        # Aplicar pequeñas perturbaciones aleatorias a new_solution
        return new_solution
