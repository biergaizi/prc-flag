CC = gcc
#CC = clang

INCLUDE = -I/usr/include/cairo -lcairo -lm
CFLAGS = -O2 -std=c99 -Wall -Wextra

all:
	${CC} prc_flag.c ${INCLUDE} ${CFLAGS} -o prc_flag
