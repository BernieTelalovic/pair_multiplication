from pair_multiplication import *

if __name__ == "__main__":
    # Create a Young diagram with partition (2, 1)
    partition = (2,1)
    print("Let's make a Young Diagram corresponding to the partition: "+str(partition))
    young_diagram = YoungDiagram(partition)

    # Display the Young diagram
    print("Young diagram representation:")
    print(young_diagram)
    
    print("The first Nc where it appears is: Nc="+str(young_diagram.N0))
    
    dimension_Nc3 = young_diagram.dimension_Nc(Nc=3)
    print("Under Nc=3, it has dimension: "+str(dimension_Nc3))
    
    

