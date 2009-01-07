#!/bin/sh

for i in autoconf; do
    echo "$i"
    $i
    if [ $? -ne 0 ]; then
	echo "Error $? in $i"
	exit 1
    fi
done

echo "./configure --enable-autogen --enable-flex $@"
./configure --enable-autogen --enable-flex $@
if [ $? -ne 0 ]; then
    echo "Error $? in ./configure"
    exit 1
fi
