# Efficient column and row scaling for sparse matrices

`colScale` normalises a matrix by the sum of each column

`rowScale` normalises a matrix by the sum of each row

sqrt of both column and row normalisation

## Usage

``` r
colScaleM(A, na.rm = TRUE)

rowScaleM(A, na.rm = TRUE)

geomScaleM(A, na.rm = TRUE)
```

## Arguments

- A:

  a sparse (or dense) matrix

- na.rm:

  when `T`, the default, converts any NA values (usually resulting from
  divide by zero errors when a neuron has no partners) to 0

## Value

A scaled sparse matrix

## Details

Conventionally in connectivity matrices, rows are upstream input neurons
and columns are downstream output neurons. `colScaleM` will therefore
normalise by the total input onto each downstream neuron while
`rowScaleM` will normalise by the total output from each upstream
neuron.

These functions will work with dense inputs but are less efficient than
their base R cousins (e.g. `scale`). As an example, normalising a 4x4k
subset of VNC connectivity data took

- `colScaleM` sparse 3.5ms

- `colScaleM` dense 1460ms

- `scale` dense 365ms

In conclusion sparse matrices with `colScaleM`/`rowScaleM` are the way
to go for real connectivity matrices! For large matrices they will be
*much* faster and for small ones, both methods will be fast anyway.

Note also that these functions were originally called `colScale` etc but
then the [`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html) package
added functions of the same name. The functions in this package differ
in two respects from the
`Matrix::`[`colScale,rowScale`](https://rdrr.io/pkg/Matrix/man/dimScale.html)
functions. First they calculate the required row or column sums to
perform the scaling. Second, they handle the inevitable zero sum rows or
columns via the default `na.rm` argument.

## See also

`Matrix::`[`colScale`](https://rdrr.io/pkg/Matrix/man/dimScale.html) on
which these are based and, for the original version,
<https://stackoverflow.com/questions/39284774/column-rescaling-for-a-very-large-sparse-matrix-in-r>

## Examples

``` r
library(Matrix)
set.seed(42)
A <- Matrix(rbinom(100,10,0.05), nrow = 10)
rownames(A)=letters[1:10]
colnames(A)=LETTERS[1:10]

colScaleM(A)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 10 column names ‘A’, ‘B’, ‘C’ ... ]]
#>                                                                                
#> a 0.2222222 .         0.125 0.1666667 .         .         0.1428571 . .   0.125
#> b 0.2222222 0.1428571 .     0.1666667 .         .         0.2857143 . .   .    
#> c .         0.2857143 0.375 .         .         .         0.1428571 . .   .    
#> d 0.1111111 .         0.250 0.1666667 0.2222222 0.3333333 .         . 0.5 0.250
#> e 0.1111111 .         .     .         .         .         0.1428571 . 0.5 0.250
#> f .         0.2857143 .     0.1666667 0.2222222 0.3333333 .         1 .   0.125
#> g 0.1111111 0.2857143 .     .         0.1111111 0.3333333 .         . .   .    
#> h .         .         0.125 .         0.1111111 .         0.1428571 . .   .    
#> i 0.1111111 .         .     0.1666667 0.2222222 .         0.1428571 . .   0.125
#> j 0.1111111 .         0.125 0.1666667 0.1111111 .         .         . .   0.125
rowScaleM(A)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 10 column names ‘A’, ‘B’, ‘C’ ... ]]
#>                                                                              
#> a 0.3333333 .         0.1666667 0.1666667 .         .     0.1666667 .     .  
#> b 0.3333333 0.1666667 .         0.1666667 .         .     0.3333333 .     .  
#> c .         0.3333333 0.5000000 .         .         .     0.1666667 .     .  
#> d 0.1000000 .         0.2000000 0.1000000 0.2000000 0.100 .         .     0.1
#> e 0.2000000 .         .         .         .         .     0.2000000 .     0.2
#> f .         0.2500000 .         0.1250000 0.2500000 0.125 .         0.125 .  
#> g 0.2000000 0.4000000 .         .         0.2000000 0.200 .         .     .  
#> h .         .         0.3333333 .         0.3333333 .     0.3333333 .     .  
#> i 0.1666667 .         .         0.1666667 0.3333333 .     0.1666667 .     .  
#> j 0.2000000 .         0.2000000 0.2000000 0.2000000 .     .         .     .  
#>            
#> a 0.1666667
#> b .        
#> c .        
#> d 0.2000000
#> e 0.4000000
#> f 0.1250000
#> g .        
#> h .        
#> i 0.1666667
#> j 0.2000000
geomScaleM(A)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 10 column names ‘A’, ‘B’, ‘C’ ... ]]
#>                                                                        
#> a 0.2721655 .         0.1443376 0.1666667 .         .         0.1543033
#> b 0.2721655 0.1543033 .         0.1666667 .         .         0.3086067
#> c .         0.3086067 0.4330127 .         .         .         0.1543033
#> d 0.1054093 .         0.2236068 0.1290994 0.2108185 0.1825742 .        
#> e 0.1490712 .         .         .         .         .         0.1690309
#> f .         0.2672612 .         0.1443376 0.2357023 0.2041241 .        
#> g 0.1490712 0.3380617 .         .         0.1490712 0.2581989 .        
#> h .         .         0.2041241 .         0.1924501 .         0.2182179
#> i 0.1360828 .         .         0.1666667 0.2721655 .         0.1543033
#> j 0.1490712 .         0.1581139 0.1825742 0.1490712 .         .        
#>                                
#> a .         .         0.1443376
#> b .         .         .        
#> c .         .         .        
#> d .         0.2236068 0.2236068
#> e .         0.3162278 0.3162278
#> f 0.3535534 .         0.1250000
#> g .         .         .        
#> h .         .         .        
#> i .         .         0.1443376
#> j .         .         0.1581139
```
