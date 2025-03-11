pair_multiplication is a Python library for working with Young Diagrams and pairs of Young diagrams.


## Features

- **NullDiagram**: Handles null cases of Young diagrams.
- **YoungDiagram** inherits from *NullDiagram*: A class to create and manipulate Young diagrams.
- **Pair** inherits from *YoungDiagram*: Represents a pair of Young diagrams.
- **DirectSum**: Creates a direct sum of diagrams, indexable like a numpy array.
- **DimensionDirectSum**: Useful for quickly showing the dimensions and multiplicities of a direct sum

The multiplication methods have been successfully numerically tested for all partitions generating diagram pairs up to 8 boxes in total (4 boxes in barred and 4 boxes in unbarred diagrams).

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


    File ~/Documents/birdtracks/pair_multiplication/pair_multiplication/classes.py:314, in YoungDiagram.__init__(self, partition, Nc, barred, weight, inherited_N0)
        312     self.width = partition[0]
        313     if not all(i >= j for i, j in zip(partition, partition[1:])):
    --> 314         raise ValueError('Not a young diagram.')
        315 elif len(partition) > 0:
        316     self.width = partition[0]


    ValueError: Not a young diagram.


### Young Diagrams

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

We will use "barred" Young diagrams to label representations of antiquarks, arising from the complex-conjugate of the covarient representation:


```python
yd2 = YoungDiagram((3,2,1),barred = True)
```

they are not the same objects!


```python
yd1==yd2
```




    False



But when we evaluate a barred diagram for a given $N_c$, e.g., $N_c=3$, it becomes a conventional, unbarred Young diagram:


```python
yd2_Nc3 = yd2.evaluate_for_Nc(3)
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



but for a different $N_c$, it produces a different diagram


```python
yd2.evaluate_for_Nc(4)
```




$$ (3, 2, 1) $$



## Composite Diagram Pairs

To create an arbitrary tensor representation of $SU(N)$, with any number of covarient or contravarient indices, we create a composite representation of an antiquark-like diagram and a quark-like diagram, composing them in a diagram pair. 

They're constructed from two tuples, the first one defining the barred diagram, and the second the unbarred diagram.


```python
pair1 = Pair(((1),(1)))
pair1
```




$$ 1_{2} \left(\overline{(1)},(1)\right) $$



In the notation, the subscript denotes the lowest $N_c$ where this representation can exist. Like a barred diagram, we can evaluate it for a certain $N_c$:


```python
pair1_nc3 = pair1.evaluate_for_Nc(Nc=2)
pair1_nc3
```




$$ (2) $$




```python
pair1_nc3 = pair1.evaluate_for_Nc(Nc=3)
pair1_nc3
```




$$ (2, 1) $$



A pair can also be made by pairing a YoungDiagram with a barred Young diagram, or the other way around:


```python
yd.pair_with(yd2) == yd2.pair_with(yd)
```




    True



but not if they're both barred/unbarred:


```python
yd.pair_with(yd)
```


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    Input In [21], in <cell line: 1>()
    ----> 1 yd.pair_with(yd)


    File ~/Documents/birdtracks/pair_multiplication/pair_multiplication/classes.py:661, in YoungDiagram.pair_with(self, other, inherited_N0, Nc)
        658         return Pair((other,self),Nc=None, 
        659                     inherited_N0=max(self.N0,other.N0, inherited_N0))
        660     else:
    --> 661         raise AttributeError("One of the diagrams must be barred to form a pair")
        662 else:
        663     try:


    AttributeError: One of the diagrams must be barred to form a pair


They can also be constructed by pairing a diagram with a partition, and the partition will be treated as the barred diagram if the first object is unbarred, or vice-versa:


```python
YoungDiagram((2,1)).pair_with((1))
```




$$ 1_{3} \left(\overline{(1)},(2, 1)\right) $$




```python
YoungDiagram((1),barred=True).pair_with((2,1))
```




$$ 1_{3} \left(\overline{(1)},(2, 1)\right) $$



When a pair is constructed with a given $N_c$, it immediately evalues the corresponding YoungDiagram (unbarred) instead:


```python
Pair(((1),(1)),Nc=3)
```




$$ (2, 1) $$



## Common properties of all objects

Diagrams and pairs all have the following attributes, accessed the same way:

    1) .Nc: the carried value of $N_c$ (int or None)
    2) .N0: the lowest $N_c$ for which the representation occurs
    3) .partition: the partition representation of the diagram/ pair of partitions for a pair object
    4) .barred: True only for barred Young diagrams, otherwise False
    5) .n: the number of boxes in the partition

# Multiplying diagrams and pairs - the DirectSum class

Any two of these objects can be multiplied together, if they share the same value of $N_c$, i.e., for both objects it has an equal value or is None. 

The resulting decomposition is stored in a DirectSum object, which has attriutes elements and multiplicities:


```python
yd = YoungDiagram((2,1))

ds = yd * yd
ds
```




$$ \begin{array}{c}1_{2} (3, 3)\oplus1_{2} (4, 2)\oplus
1_{3} (2, 2, 2)\oplus2_{3} (3, 2, 1)\oplus1_{3} (4, 1, 1)\oplus
1_{4} (2, 2, 1, 1)\oplus1_{4} (3, 1, 1, 1)\end{array} $$



the list containing all elements is acessed via:


```python
ds.elements
```




    array([[2](3, 3), [2](4, 2), [3](2, 2, 2), [3](3, 2, 1), [3](4, 1, 1),
           [4](2, 2, 1, 1), [4](3, 1, 1, 1)], dtype=object)



and their index-corresponding multiplicities as:


```python
ds.multiplicities
```




    array([1, 1, 1, 2, 1, 1, 1])



DirectSums can collectivley be evaluated for a given $N_c$, if it was previously None


```python
ds_nc3 = ds.evaluate_for_Nc(3)
ds_nc3
```




$$ \begin{array}{c}1_{0}\,()\oplus
1_{1}\,(3)\oplus
2_{2}\,(2, 1)\oplus1_{2}\,(4, 2)\oplus1_{2}\,(3, 3)\end{array} $$



they can be filtered in a numpy-friendly way. The following gives all the irreps of dimension 10 in SU(3):


```python
ds_nc3[ds_nc3.dimension_Nc()==10]
```




$$ \begin{array}{c}1_{1}\,(3)\oplus
1_{2}\,(3, 3)\end{array} $$



You can use the same setup to investigate the multiplicities, elements, or dimensions:


```python
ds_nc3.multiplicities[ds_nc3.elements==YoungDiagram((2,1),Nc=3)]
```




    array([2])



### Standard Littlewood-Richardson Multiplication of Young diagrams only

Pair multiplication implements the column-wise multiplication method soon to be outlined in Ref.[work in progress].
An alternative way to multiply Young diagrams *only* is to use the built-in function .LR_multiply():


```python
yd.LR_multiply(yd)
```




$$ \begin{array}{c}1_{2}\,(3, 3)\oplus1_{2}\,(4, 2)\oplus
1_{3}\,(2, 2, 2)\oplus2_{3}\,(3, 2, 1)\oplus1_{3}\,(4, 1, 1)\oplus
1_{4}\,(2, 2, 1, 1)\oplus1_{4}\,(3, 1, 1, 1)\end{array} $$



It produces the same result as the multiplication we previously constructed, but uses a different implementation - a good cross check:


```python
yd.LR_multiply(yd) == yd * yd
```




    True



## Myltiplying Composite Pairs

The general multiplication method in this package extends to multiplying all combinations of YoungDiagram, Pair and DirectSum objects. 

We'll show an example of multiplying Pairs, and dow the general $N_c$ dependence is preserved:


```python
pair1 = Pair(((1,1),(2)))
pair2 = Pair(((1),(1)))
```

then multiply them $P_1\otimes P_2$:


```python
p1_times_p2 = pair1*pair2
p1_times_p2
```




$$ \begin{array}{c}1_{3} \left(\overline{(1)},(1)\right)\oplus1_{3} \left(\overline{(1, 1)},(2)\right)\oplus1_{3} \left(\overline{(2)},(2)\right)\oplus1_{3} \left(\overline{(2, 1)},(3)\right)\oplus
1_{4} \left(\overline{(1, 1)},(1, 1)\right)\oplus1_{4} \left(\overline{(1, 1)},(2)\right)\oplus1_{4} \left(\overline{(2, 1)},(2, 1)\right)\oplus1_{4} \left(\overline{(1, 1, 1)},(3)\right)\oplus
1_{5} \left(\overline{(1, 1, 1)},(2, 1)\right)\end{array} $$



multiplying $P_2\otimes P_1$ should give the same answer:


```python
p2_times_p1 = pair1*pair2

p1_times_p2 == p2_times_p1
```




    True



We can find the lowest $N_c$ for which each representation is admissible:


```python
lowest_nc = p1_times_p2.lowest_Nc()
lowest_nc
```




    array([3, 3, 3, 3, 4, 4, 4, 4, 5])



Let's pick an $N_c$ and test the tensor multiplication:


```python
Nc = min(lowest_nc)
```

Then we can check if the tensor multiple is the same when we multiply pairs and then evaluate the $N_c$, vs. first evaluating each pair for the given $N_c$ and then multiplying their results. 

In this case, we use the same multiplication algorithm (column LR):


```python
Nc
```




    3




```python
p1_times_p2.evaluate_for_Nc(Nc)
```




$$ \begin{array}{c}1_{1}\,(3)\oplus
1_{2}\,(2, 1)\oplus1_{2}\,(4, 2)\oplus1_{2}\,(5, 1)\end{array} $$



Now we check if the tensor multiple is the same when comparing it with the implemented LR algorithm:


```python
pair1.evaluate_for_Nc(Nc).LR_multiply(pair2.evaluate_for_Nc(Nc))
```




$$ \begin{array}{c}1_{1}\,(3)\oplus
1_{2}\,(2, 1)\oplus1_{2}\,(4, 2)\oplus1_{2}\,(5, 1)\end{array} $$




```python
p1_times_p2.evaluate_for_Nc(Nc) ==\
pair1.evaluate_for_Nc(Nc).LR_multiply(pair2.evaluate_for_Nc(Nc)).evaluate_for_Nc(Nc)
```




    True



# Printing Outputs

The package has several ways to output the results of calculations. Objects themselves can be printed on the command line and in IPython.

A pythonic print() statement gives the command-line friendly output:


```python
print(p1_times_p2)
```

    1[3]((1),(1))+1[3]((1, 1),(2))+1[3]((2),(2))+1[3]((2, 1),(3))+1[4]((1, 1),(1, 1))+1[4]((1, 1),(2))+1[4]((2, 1),(2, 1))+1[4]((1, 1, 1),(3))+1[5]((1, 1, 1),(2, 1))


Here, the multiplicities are shown as indices at the front of each pair of partitions, and the square brackets [] contain the lowest $N_c$ of the representation.

Simply outputting the object itself produces the IPython formatting:


```python
p1_times_p2
```




$$ \begin{array}{c}1_{3} \left(\overline{(1)},(1)\right)\oplus1_{3} \left(\overline{(1, 1)},(2)\right)\oplus1_{3} \left(\overline{(2)},(2)\right)\oplus1_{3} \left(\overline{(2, 1)},(3)\right)\oplus
1_{4} \left(\overline{(1, 1)},(1, 1)\right)\oplus1_{4} \left(\overline{(1, 1)},(2)\right)\oplus1_{4} \left(\overline{(2, 1)},(2, 1)\right)\oplus1_{4} \left(\overline{(1, 1, 1)},(3)\right)\oplus
1_{5} \left(\overline{(1, 1, 1)},(2, 1)\right)\end{array} $$



Toggling the \LaTeX strings on and off can be done in two ways. To simply print the result, use .print(tex=True/False):


```python
p1_times_p2.print() #default is tex=False
```

    1[3]((1),(1))+1[3]((1, 1),(2))+1[3]((2),(2))+1[3]((2, 1),(3))+1[4]((1, 1),(1, 1))+1[4]((1, 1),(2))+1[4]((2, 1),(2, 1))+1[4]((1, 1, 1),(3))+1[5]((1, 1, 1),(2, 1))



```python
p1_times_p2.print(tex=True) 
```

    $$ \begin{array}{c}1_{3} \left(\overline{(1)},(1)\right)\oplus1_{3} \left(\overline{(1, 1)},(2)\right)\oplus1_{3}\left(\overline{(2)},(2)\right)\oplus1_{3} \left(\overline{(2, 1)},(3)\right)\oplus    1_{4} \left(\overline{(1, 1)},(1, 1)\right)\oplus1_{4} \left(\overline{(1, 1)},(2)\right)\oplus1_{4} \left(\overline{(2, 1)},(2, 1)\right)\oplus1_{4} \left(\overline{(1, 1, 1)},(3)\right)\oplus    1_{5} \left(\overline{(1, 1, 1)},(2, 1)\right)\end{array} $$


Or, to get the string objects producing the outputs, you can call .to_str(tex=True/False):


```python
string = p1_times_p2.to_str() #default is tex=False
string
```




    '1[3]((1),(1))+1[3]((1, 1),(2))+1[3]((2),(2))+1[3]((2, 1),(3))+1[4]((1, 1),(1, 1))+1[4]((1, 1),(2))+1[4]((2, 1),(2, 1))+1[4]((1, 1, 1),(3))+1[5]((1, 1, 1),(2, 1))'




```python
tex_string = p1_times_p2.to_str(tex = True) 
tex_string
```




    '\\[\\begin{array}{c}1_{3} \\left(\\overline{(1)},(1)\\right)\\oplus1_{3} \\left(\\overline{(1, 1)},(2)\\right)\\oplus1_{3} \\left(\\overline{(2)},(2)\\right)\\oplus1_{3} \\left(\\overline{(2, 1)},(3)\\right)\\oplus\n\n1_{4} \\left(\\overline{(1, 1)},(1, 1)\\right)\\oplus1_{4} \\left(\\overline{(1, 1)},(2)\\right)\\oplus1_{4} \\left(\\overline{(2, 1)},(2, 1)\\right)\\oplus1_{4} \\left(\\overline{(1, 1, 1)},(3)\\right)\\oplus\n\n1_{5} \\left(\\overline{(1, 1, 1)},(2, 1)\\right)\\end{array}\\]'



## Coming soon:

 - ~ordering elements in the direct sum in a readable way~
 - ~algorithm speed-up for higher numbers of boxes~
 - ~better handling of diagram multiplicities~
 - better documentation and testing
 - ~more Latexing functions!~
