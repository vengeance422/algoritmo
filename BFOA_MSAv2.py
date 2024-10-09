from bacteria import bacteria
from chemiotaxis import chemiotaxis

import numpy

poblacion = []
path = "C:\\secuenciasBFOA\\multiFasta.fasta"
numeroDeBacterias = 15  # Aumentamos el número de bacterias
numRandomBacteria = 2   # Insertamos más bacterias aleatorias
iteraciones = 15        # Reducimos ligeramente el número de iteraciones
tumbo = 2               # Incrementamos el tumbo para más variabilidad
nado = 4                # Incrementamos el nado para más exploración
chemio = chemiotaxis()
veryBest = bacteria(path)                #mejor bacteria   
tempBacteria = bacteria(path)            #bacteria temporal para validaciones
original = bacteria(path)                #bacteria original sin gaps
globalNFE = 0      #numero de evaluaciones de la funcion objetivo

# Ajustamos los parámetros de atracción y repulsión
dAttr = 0.15  # Aumentamos ligeramente la atracción
wAttr = 0.3   # Aumentamos la fuerza de atracción
hRep = 0.2    # Aumentamos la repulsión
wRep = 15     # Aumentamos la fuerza de repulsión

def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction
    
def validaSecuencias(path, veryBest):
    #clona a veryBest en tempBacteria   
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    #descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-","")
    #tempBacteria.tumboNado(1)    

    #valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return
      
        
    
    
    
    
    
    
    


for i in range(numeroDeBacterias):                                            #poblacion inicial
    poblacion.append(bacteria(path))


for _ in range(iteraciones):                                                  #numero de iteraciones  
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)
        #bacteria.tumboNado(nado)
        bacteria.autoEvalua()  
        #print("blosumScore: ",bacteria.blosumScore)
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)                 #d_attr, w_attr, h_rep, w_rep):
    globalNFE += chemio.parcialNFE 
    best = max(poblacion, key=lambda x: x.fitness)
    if (veryBest == None) or (best.fitness > veryBest.fitness):
         clonaBest(veryBest, best)
    print("interaccion: ",veryBest.interaction,"fitness: ",veryBest.fitness, " NFE:",globalNFE )
    
    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)                #inserta  bacterias aleatorias
    print("poblacion: ",len(poblacion))


veryBest.showGenome()
validaSecuencias(path, veryBest)