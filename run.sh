#!/bin/bash

gcc -Wall -O2 -Werror -std=gnu11 -g -ftrapv -Wformat-security -Wignored-qualifiers -Winit-self -Wswitch-default -Wshadow -Wpointer-arith -Wtype-limits -Wempty-body -Wlogical-op -Wstrict-prototypes -Wold-style-declaration -Wold-style-definition -Wmissing-parameter-type -Wmissing-field-initializers -Wnested-externs -Wno-pointer-sign -lm -fsanitize=address,undefined -o a.out main.c
