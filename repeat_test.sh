#!/usr/bin/env bash

ulimit -c unlimited

rm -rf core

for i in {1..100}; do
	if ! ./build/test_packing3d < ../static/request.json; then
		break
	fi
done
