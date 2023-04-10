#!/bin/sh
debug_flag=''
while getopts 'dr' flag; do
    case "${flag}" in
        d) debug_flag='Debug' ;;
        r) debug_flag='Release' ;;
        *) print_usage
            exit 1 ;;
    esac
done

echo "$debug_flag Build"
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE="$debug_flag" ..
make
