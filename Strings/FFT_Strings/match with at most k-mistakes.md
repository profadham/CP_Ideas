# note that in FFT we can use distance matching
### for example lets see this string and pattern:
acabergxa  |  axa

lets match each character alone

for a:

##### string:  array a =  1 0 1 0 0 0 0 0 1

##### pattern: array b =  1 0 1

then use distance matching

lets say I want to find:

$$
\sum_{i=1}^{n} a[i]b[i]
$$

