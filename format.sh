#!/usr/bin/env sh

find $(pwd)/src -iname '*.h' -o -iname '*.cpp' | xargs clang-format -i
