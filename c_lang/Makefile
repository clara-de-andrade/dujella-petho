CC = gcc
CFLAGS = -Wall -Wextra 
INCLUDES = -I./
LIBS = -lm -lgmp -lmpfr

SRCS = 	findroot.c 	\
       	polynom.c  	\
		nthterm.c   \
		dptest.c
OBJS = $(SRCS:.c=.o)

.PHONY: build clean

test_1: $(SRCS) test_1.c
	$(CC) $(CFLAGS) $(INCLUDES) $(SRCS) test_1.c -o test_1 $(LIBS)
	./test_1

test_2: $(SRCS) test_2.c
	$(CC) $(CFLAGS) $(INCLUDES) $(SRCS) test_2.c -o test_2 $(LIBS)
	./test_2

clean:
	rm -f test_1 test_2
