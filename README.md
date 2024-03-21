# pnormGO

## About

Rewrite the pnorm function of R in Golang.

https://github.com/wch/r-source/blob/tags/R-3-3-0/src/nmath/pnorm.c

## Demo

```
package main

import pnormGO "github.com/tallcode/pnormGO"

func main() {
	println(pnormGO.Pnorm(30, 0, 1, false, false))
}
```