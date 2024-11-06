from pair_multiplication import *

if __name__ == "__main__":
    # Create a barred Young diagram with partition (1)
    partition = (1)
    print("Let's make a barred Young Diagram corresponding to the partition: "+str(partition))
    young_diagram = YoungDiagram(partition,barred=True)

    # Display the Young diagram
    print("Young diagram representation:")
    print(young_diagram)
    
    dimension_Nc3 = young_diagram.dimension_Nc(Nc = 3)
    print("Under Nc=3, it has dimension: "+str(dimension_Nc3))
    
    print("The first Nc where it appears is: Nc="+str(young_diagram.N0))
    
    young_diagram_when_Nc_is_3 = young_diagram.evaluate_for_Nc(Nc = 3) # this is a new object, an unbarred YoungDiagram
    print("Under Nc=3, it has representation: "+str(young_diagram_when_Nc_is_3))
    print("Is the new (unbarred) representation equal to the barred one generally?: ", 
    young_diagram_when_Nc_is_3 == young_diagram)
    
    dimension_Nc3_unbarred = young_diagram_when_Nc_is_3.dimension_Nc(Nc = 3)
    print("We can also get the same dimension from the unbarred representation we recovered: ", dimension_Nc3_unbarred)
    
    print("Let's check if the two ways of getting dimensions give the same results: ",
    dimension_Nc3 == dimension_Nc3_unbarred)
    
    

