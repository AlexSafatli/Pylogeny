#!/bin/sh
wget http://libpll.org/Downloads/libpll-1.0.2-sse3-64.tar.gz
tar -xvf libpll-1.0.2-sse3-64.tar.gz
cd libpll-1.0.2-sse3-64 && cp libpll* /usr/local/lib/ && cp -r include/pll /usr/local/include/
