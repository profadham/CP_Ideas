# note that in FFT we can use distance matching
### for example lets see this string and pattern:
acabergxa  |  axa

lets match each character alone

for a:

##### string:  array a =  1 0 1 0 0 0 0 0 1

##### pattern: array b =  1 0 1

then use distance matching

lets say, for a given k, I want to find:

$$
\sum_{i=1}^{n} a[i]b[i+k]
$$

### FFT works only to calculate this convolution

what does this represent??

for k = 0 this means I want to see how many (a) characters match when pattern is standing under i = 0 of the string

for k = 1 this means we shifted the pattern

### so now we can know for each character for every index how many times this character matches if pattern is put at that index


