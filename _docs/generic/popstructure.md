---
title: Population Structure
permalink: popstructure.html
sidebar: generic
#tags: [getting-started, formatting]
product: Generic
---

## Test of post

<!--<p><strong>Test of post</strong></p>-->

Add a table:

| Col1 | Col2 |
|----|-----|
| row 1 | row 1 |
| row 2 | row 2 |

## Chunks with theme

<!--<p><strong>Chunks with theme</strong></p>-->

Python chunk:

```python
def function(x):
    if x == True:
        print('Hello World!')
    else:
	print('Goodbye World?')
```

Bash chunck:

```bash
for i in `seq 1 10`
do
    echo Hello world $i
done
```

R chunck:

```r
library(tidyverse)
x <- read.csv('/home/agomez/dataset.csv')
x %>% filter(level == 'High')
```

## Highlight chunks

<!--<p><strong>Highlight chunks</strong></p>-->

<p>Python</p>
{% highlight python %}
def function(x):
    if x == True:
        print('Hello World!')
    else:
	print('Goodbye World?')
{% endhighlight %}


<p>Bash</p>
{% highlight bash %}
for i in `seq 1 10`
do
    echo Hello world $i
done
{% endhighlight %}

<p>R</p>
{% highlight R %}
library(tidyverse)
x <- read.csv('/home/agomez/dataset.csv')
x %>% filter(level == 'High')
{% endhighlight %}

