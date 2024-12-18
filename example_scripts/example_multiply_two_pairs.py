from pair_multiplication import *

if __name__ == "__main__":

    partition1 = ((1),(1,1))
    print("Let's make a Diagram Pair corresponding to the partition pair: "+str(partition1))
    pair1 = Pair(partition1)
    
    partition2 = ((2),(1))
    print("And another corresponding to the partition pair: "+str(partition2))
    pair2 = Pair(partition2)
    
    p1_times_p2 = pair1*pair2 # multiplying them (this uses the Littewood-Richardson rule)
    print('Their tensor multiple is: ',p1_times_p2)
    print('\n\n')
    
    #NOTE: yd1_times_yd2 is a DirectSum object, and inherits from dict. 
    print("The result of the multiplication is an object that you can use like a dictionary.")
    print("The keys are YoungDiagram or Pair objects: ",list(p1_times_p2.keys()))
    print("The lowesr Nc under which they appear are stored as: ",list(p1_times_p2.lowest_Nc()))
    print("The multiplicities are stored as the values for each key: ",list(p1_times_p2.values()))
    print('\n\n')
    # NOTE: the list with all the diagram/pair objects in the direct sum can also be accessed with 
    #       yd1_times_yd2.elements() and the multiplicities can be accessed with yd1_times_yd2.multiplicities()
    
    Nc3_p1_times_p2 = p1_times_p2.evaluate_for_Nc(Nc=3)
    print("For Nc=3, the resulting young diagrams are: ",Nc3_p1_times_p2)
    print()
    Nc3_p1_times_Nc3_p2 = p1.evaluate_for_Nc(Nc=3)*p2.evaluate_for_Nc(Nc=3)
    print("We can first evaluate each pair for Nc=3, then multiply the Young diagrams: ",Nc3_p1_times_Nc3_p2)
    print()
    print("We can check that they're the same: ", Nc3_p1_times_Nc3_p2==Nc3_p1_times_p2)
    
    print('\n\n')
    
    # we can also get the dimensions of all the resulting Young Diagrams for a given Nc:
    print("For Nc=3, the dimensions of the resulting diagrams are: ",p1_times_p2.dimension_Nc(Nc=3))
    print("For Nc=3, the sum of the dimensions weighed by the multiplicities is: ",p1_times_p2.dimension_Nc(Nc=3).sum())
    print("When first evaluating Nc=3 on both pairs, the dimensions are: ", Nc3_p1_times_Nc3_p2.dimension_Nc(Nc=3))
    print("And the sum is: ", Nc3_p1_times_Nc3_p2.dimension_Nc(Nc=3))
    
