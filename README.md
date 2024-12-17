pair_multiplication is a Python library for working with Young Diagrams and pairs of Young diagrams.


## Features

- **NullDiagram**: Handles null cases of Young diagrams.
- **YoungDiagram** inherits from *NullDiagram*: A class to create and manipulate Young diagrams.
- **Pair** inherits from *YoungDiagram*: Represents a pair of Young diagrams.
- **DirectSum** inherits from *dict*: Creates a direct sum of diagrams.
- **DimensionDirectSum** inherits from *DirectSum*: Useful for quickly showing the dimensions and multiplicities of a direct sum

## Installation

Install the package using `pip`:

```bash
pip install git+https://github.com/BernieTelalovic/pair_multiplication.git

```


```python
from pair_multiplication import *
```

# YoungDiagram class construction

Starting from a partition tuple labelling a diagram, $(a_1,a_2,..)$ with $a_1\geq a_2\geq...$ being the number of boxes in each row labelled by the subscript, we can construct young diagrams:


```python
yd = YoungDiagram((3,2,1))
```

but not:


```python
yd = YoungDiagram((1,2,3)) # this is supposed to give an error, don't worry
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    Input In [3], in <cell line: 1>()
    ----> 1 yd = YoungDiagram((1,2,3))


    File ~/Documents/birdtracks/pair_multiplication/pair_multiplication/classes.py:283, in YoungDiagram.__init__(self, partition, Nc, barred, weight, inherited_N0)
        281 if len(partition)>1:
        282     if not all(i >= j for i, j in zip(partition, partition[1:])):
    --> 283         raise ValueError('Not a young diagram.')
        286 N0 = len(partition)
        288 if not (Nc is None):


    ValueError: Not a young diagram.


### Quark-like Young Diagrams

We will use conventional Young diagrams to label representations of quarks:


```python
yd1 = YoungDiagram((3,2,1))
```


```python
yd1
```




$$ (3, 2, 1) $$



We can access the tuple corresponding to the partition:


```python
yd1.partition
```




    (3, 2, 1)



the lowest Nc for which this digram could exist:


```python
yd1.N0
```




    3



its representation for a given Nc:


```python
Nc = 3
yd1_Nc3 = yd1.evaluate_for_Nc(Nc)
```


```python
yd1_Nc3
```




$$ (2, 1) $$



its dimension for a given Nc:


```python
yd1.dimension_Nc(Nc)
```




    8



### Antiquark-like Young Diagrams

We will use "barred" Young diagrams to label representations of quarks, arising from the adjoint representation:


```python
yd2 = YoungDiagram((3,2,1),barred = True)
```

they are not the same objects!


```python
yd1==yd2
```




    False



But when we evaluate a barred diagram for a given Nc, it becomes a conventional, unbarred Young diagram:


```python
yd2_Nc3 = yd2.evaluate_for_Nc(Nc)
```


```python
yd2_Nc3
```




$$ (2, 1) $$



in this case, its the same as the previous diagram when Nc=3:


```python
yd2_Nc3==yd1_Nc3
```




    True



# Multiplying Young diagrams - the DirectSum class

### (Un)barred diagram-(un)barred diagram multiplication

Uses Littlewood-Richardson rule for diagram multiplication


```python
yd_unbarred = YoungDiagram((2,1))
```


```python
ydubr_tensor_ydubr = yd_unbarred.LR_multiply(yd_unbarred)
```

The resulting object is a DirectSum class, which displays the direct (tensor) sum of YoungDiagram objects. The DirectSum class is a type of python dict, with keys being the YoungDiagram/Pair (more on these below) objects, and the values corresponding to the multiplicities. In the display, the constants are the multiplicities and their subscripts are the first Nc where this pair labels a young diagram.


```python
ydubr_tensor_ydubr
```




$$ 1_{2}\,(4, 2)\oplus1_{2}\,(3, 3)\oplus2_{3}\,(3, 2, 1)\oplus1_{3}\,(4, 1, 1)\oplus1_{3}\,(2, 2, 2)\oplus1_{4}\,(3, 1, 1, 1)\oplus1_{4}\,(2, 2, 1, 1) $$



the elements can be accessed using:


```python
ydubr_tensor_ydubr.keys() # this returns dict_keys
```




    dict_keys([(3, 3), (4, 1, 1), (4, 2), (2, 2, 1, 1), (3, 2, 1), (2, 2, 2), (3, 1, 1, 1)])



or:


```python
ydubr_tensor_ydubr.elements()# this gives a list
```




    [(3, 3), (4, 1, 1), (4, 2), (2, 2, 1, 1), (3, 2, 1), (2, 2, 2), (3, 1, 1, 1)]



the multiplicities can be recovered as a list in two ways as well:


```python
ydubr_tensor_ydubr.values()
```




    dict_values([1, 1, 1, 1, 2, 1, 1])



or:


```python
ydubr_tensor_ydubr.multiplicities()
```




    [1, 1, 1, 1, 2, 1, 1]



the lowest nc for each diagram can be separately recovered:


```python
ydubr_tensor_ydubr.lowest_Nc()
```




    [2, 3, 2, 4, 3, 3, 4]



We can evaluate it under a given Nc:


```python
Nc = 3
ydubr_tensor_ydubr_nc3 = ydubr_tensor_ydubr.evaluate_for_Nc(Nc=Nc)
```


```python
ydubr_tensor_ydubr_nc3
```




$$ 1_{0}\,()\oplus1_{1}\,(3)\oplus2_{2}\,(2, 1)\oplus1_{2}\,(4, 2)\oplus1_{2}\,(3, 3) $$



we can get the dimension in the same way:


```python
ydubr_tensor_ydubr_nc3.dimension_Nc() # here the direct sum already knows which Nc we used
```




$$ 2\cdot8+1\cdot1+2\cdot10+1\cdot27 $$



and the sum:


```python
ydubr_tensor_ydubr_nc3.dimension_Nc().sum()
```




    64



Similar when multiplying two barred young diagrams:


```python
yd_barred = YoungDiagram((2,1),barred = True)
ydbr_tensor_ydbr = yd_barred*yd_barred
```


```python
ydbr_tensor_ydbr
```




$$ 1_{2}\,\overline{(4, 2)}\oplus1_{2}\,\overline{(3, 3)}\oplus1_{3}\,\overline{(2, 2, 2)}\oplus1_{3}\,\overline{(4, 1, 1)}\oplus2_{3}\,\overline{(3, 2, 1)}\oplus1_{4}\,\overline{(2, 2, 1, 1)}\oplus1_{4}\,\overline{(3, 1, 1, 1)} $$



We can also evaluate it under an Nc:


```python
ydbr_tensor_ydbr_nc3 = ydbr_tensor_ydbr.evaluate_for_Nc(Nc=Nc)
```


```python
ydbr_tensor_ydbr_nc3
```




$$ 1_{0}\,()\oplus1_{1}\,(3)\oplus2_{2}\,(2, 1)\oplus1_{2}\,(4, 2)\oplus1_{2}\,(3, 3) $$



this is now the same as the DirectSum containing the first tensor multiple above:


```python
ydbr_tensor_ydbr_nc3==ydubr_tensor_ydubr_nc3
```




    True



we can also get the dimensions directly (this time specifying the Nc)


```python
ydbr_tensor_ydbr.dimension_Nc(Nc=Nc)
```

    /home/anduril/Documents/birdtracks/pair_multiplication/pair_multiplication/classes.py:857: UserWarning: Conjugate not admissible under given Nc.
      warnings.warn('Conjugate not admissible under given Nc.')





$$ 2\cdot8+1\cdot1+2\cdot10+1\cdot27 $$



# Barred diagram-unbarred diagram multiplication - the Pair class

We do the multiplication using King's Q rule for diagram multiplication


```python
yd_barred = YoungDiagram((2,1),barred = True)
yd_unbarred = YoungDiagram((2,1))

barred_tensor_unbarred = yd_barred.LR_multiply(yd_unbarred)
```

The results is a DirectSum of Pair objects, where the first partition is always the barred diagram, the second is always the unbarred diagram.


```python
barred_tensor_unbarred
```




$$ 1_{2}\,()\oplus2_{2}\,\left(\overline{(1)},(1)\right)\oplus1_{2}\,\left(\overline{(2)},(2)\right)\oplus1_{3}\,\left(\overline{(2)},(1, 1)\right)\oplus1_{3}\,\left(\overline{(1, 1)},(2)\right)\oplus1_{4}\,\left(\overline{(1, 1)},(1, 1)\right)\oplus1_{4}\,\left(\overline{(2, 1)},(2, 1)\right) $$



To construct a pair we can either give a tuple of partitions (the first one is always the barred one)


```python
pair_from_partitions = Pair(((2,1),(2,1)))
```

The multiple of 1 next to it stores the lowest Nc as its subscript


```python
pair_from_partitions
```




$$ 1_{4}\,\left(\overline{(2, 1)},(2, 1)\right) $$



we can also construct it using two Young diagrams:


```python
yd_barred = YoungDiagram((2,1),barred = True)
yd_unbarred = YoungDiagram((2,1))

pair_from_diagrams = Pair((yd_barred,yd_unbarred))
```

they're the same:


```python
pair_from_partitions==pair_from_diagrams
```




    True



another way is to pair one Young diagram with either a partition:


```python
pair_from_diag_and_partition = yd_barred.pair_with((2,1))
```

(when using this method, the given partition will create a diagram that is unbarred if yd_barred and vice-versa.)

We can pair a diagram with another diagram:


```python
pair_from_diag_and_diag = yd_barred.pair_with(yd_unbarred)
```

(in this case, one must be barred and one must be unbarred, but they can be given in either order.)

they're all the same:


```python
pair_from_partitions==pair_from_diag_and_diag and pair_from_diag_and_diag==pair_from_diag_and_partition
```




    True



we can evaluate the Young diagram resulting from a given Nc in the usual way:


```python
yd_Nc7 = pair_from_partitions.evaluate_for_Nc(Nc=7)
```


```python
yd_Nc7
```




$$ (4, 3, 2, 2, 2, 1) $$



For an Nc lower than this diagrams lowest Nc, we get a NullDiagram 


```python
pair_from_partitions.evaluate_for_Nc(Nc=3)
```

    /home/anduril/Documents/birdtracks/pair_multiplication/pair_multiplication/classes.py:857: UserWarning: Conjugate not admissible under given Nc.
      warnings.warn('Conjugate not admissible under given Nc.')





$$ None $$



# Myltiplying Young Diagram Pairs

Using the column adding rules


```python
pair1 = Pair(((2,1),(1,1)))
pair2 = Pair(((1),(1)))
```


```python
pair_multiple = pair1*pair2
pair_multiple
```




$$ 1_{4}\,\left(\overline{(1, 1)},(1)\right)\oplus1_{4}\,\left(\overline{(2)},(1)\right)\oplus1_{4}\,\left(\overline{(2, 1)},(2)\right)\oplus2_{4}\,\left(\overline{(2, 1)},(1, 1)\right)\oplus1_{4}\,\left(\overline{(3)},(1, 1)\right)\oplus1_{4}\,\left(\overline{(3, 1)},(2, 1)\right)\oplus1_{4}\,\left(\overline{(2, 2)},(2, 1)\right)\oplus1_{5}\,\left(\overline{(1, 1, 1)},(1, 1)\right)\oplus1_{5}\,\left(\overline{(2, 1)},(1, 1)\right)\oplus1_{5}\,\left(\overline{(2, 2)},(1, 1, 1)\right)\oplus1_{5}\,\left(\overline{(3, 1)},(1, 1, 1)\right)\oplus1_{5}\,\left(\overline{(2, 1, 1)},(2, 1)\right)\oplus1_{6}\,\left(\overline{(2, 1, 1)},(1, 1, 1)\right) $$




```python
lowest_nc = min(pair_multiple.lowest_Nc())
lowest_nc = 6
```


```python
pair_multiple.evaluate_for_Nc(lowest_nc) == pair1.evaluate_for_Nc(lowest_nc)*pair2.evaluate_for_Nc(lowest_nc)
```

    /home/anduril/Documents/birdtracks/pair_multiplication/pair_multiplication/classes.py:385: UserWarning: Diagram multiplication performed under specific Nc.
      warnings.warn('Diagram multiplication performed under specific Nc.')
    /home/anduril/Documents/birdtracks/pair_multiplication/pair_multiplication/classes.py:789: UserWarning: Diagram multiplication performed under specific Nc.
      warnings.warn('Diagram multiplication performed under specific Nc.')





    True




```python
pair_multiple.evaluate_for_Nc(lowest_nc) ==\
pair1.evaluate_for_Nc(lowest_nc).LR_multiply(pair2.evaluate_for_Nc(lowest_nc))
```

    /home/anduril/Documents/birdtracks/pair_multiplication/pair_multiplication/classes.py:631: UserWarning: Diagram multiplication performed under specific Nc.
      warnings.warn('Diagram multiplication performed under specific Nc.')





    True



## Coming soon:

 - better handling of diagram multiplicities
 - better documentation and testing
 - more Latexing functions!


```python

```
