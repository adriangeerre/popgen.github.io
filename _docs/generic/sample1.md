---
title: Sample 1
permalink: sample1.html
sidebar: generic
tags: [getting-started, formatting]
product: Generic
---

<p><strong>Test of post</strong></p>

Add a table:

| sample | table |
|----|-----|
| row 1 | row 1 |
| row 2 | row 2|

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


R chunck:

```r
library(tidyverse)
x <- read.csv('/home/agomez/dataset.csv')
x %>% filter(level == 'High')
```
