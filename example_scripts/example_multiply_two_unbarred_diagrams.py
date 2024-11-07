from pair_multiplication import *

if __name__ == "__main__":

    partition1 = (2,1)
    print("Let's make a Young Diagram corresponding to the partition: "+str(partition1))
    young_diagram1 = YoungDiagram(partition1)
    
    partition2 = (2,1)
    print("And another corresponding to the partition: "+str(partition2))
    young_diagram2 = YoungDiagram(partition2)
    
    yd1_times_yd2 = young_diagram1*young_diagram2 # multiplying them (this uses the Littewood-Richardson rule)
    print('Their tensor multiple is: ',yd1_times_yd2)
    print()
    
    #NOTE: yd1_times_yd2 is a DirectSum object, and inherits from dict. 
    print("The result of the multiplication is an object that you can use like a dictionary.")
    print("The keys are YoungDiagram or Pair objects: ",list(yd1_times_yd2.keys()))
    print("The multiplicities are stored as the values for each key: ",list(yd1_times_yd2.values()))
    print()
    # NOTE: the list with all the diagram/pair objects in the direct sum can also be accessed with 
    #       yd1_times_yd2.elements() and the multiplicities can be accessed with yd1_times_yd2.multiplicities()
    
    
    print("For Nc=3, the resulting young diagrams are: ",yd1_times_yd2.evaluate_for_Nc(Nc=3))
    # we can also get the dimensions of all the resulting Young Diagrams for a given Nc:
    print("For Nc=3, the dimensions of the resulting diagrams are: ",yd1_times_yd2.dimension_Nc(Nc=3))
    print("For Nc=3, the sum of the dimensions weighed by the multiplicities is: ",yd1_times_yd2.dimension_Nc(Nc=3).sum())
