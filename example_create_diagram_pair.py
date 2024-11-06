from pair_multiplication import *

if __name__ == "__main__":

    partition_barred = (2,1)
    partition_unbarred = (2)
    pair = Pair((partition_barred,partition_unbarred)) # the partitions are always given
                                                       # in the order (barred,unbarred)
    print("Let's make a Young Diagram Pair!")
    print("With the barred diagram given by the partition: "+str(partition_barred))
    print("With the unbarred diagram given by the partition: "+str(partition_unbarred)+'\n')
    
    print("Its representation is: "+str(pair))
    print("The representation has the Nc where the diagram can first arise in square brackets: "+str(pair)[0:3])
    print("then the barred diagram is the first partition in the pair: "+str(pair)[4:10])
    print("and the unbarred diagram is the second partition in the pair: "+str(pair)[11:14])
    print("The first Nc where this pair can label a diagram is: "+str(pair.N0))
    print()
    
    dimension_Nc3 = pair.dimension_Nc(Nc = 3)
    print("Under Nc=3, "+str(pair)+" has dimension: ", dimension_Nc3)
    
    young_diagram_when_Nc_is_3 = pair.evaluate_for_Nc(Nc = 3) # this is a new object, an unbarred YoungDiagram
    print("Under Nc=3, "+str(pair)+" has representation: "+str(young_diagram_when_Nc_is_3))
    print("Is the new (Nc-dependent) representation equal to the pair representation?: ", 
    young_diagram_when_Nc_is_3 == pair)
    print()
    
    print("We can also create the Young diagrams first, and then combine them in a pair (see code).")
    yd_barred = YoungDiagram(partition_barred,barred = True)
    yd_unbarred = YoungDiagram(partition_unbarred)
    
    pair_from_yd = Pair((yd_barred,yd_unbarred))
    print("Do they produce the same object?: ",pair_from_yd==pair)
    print()
    
    print("We can pair one YoungDiagram object with another YoungDiagram, or a partition to produce a Pair (see code)")
    pair_from_partition_barred = yd_barred.pair_with(partition_unbarred)
    pair_from_yd_barred = yd_barred.pair_with(yd_unbarred)
    print("Do they all create the same objects?: ",pair_from_partition_barred==pair_from_yd_barred and pair_from_yd_barred==pair)
    
    
    
    
    
    
