CC = gcc
CFLAGS = -Wall -Wextra 
INCLUDES = -I./
LIBS = -lm -lgmp -lmpfr

SRCS = 	findroot.c 	\
       	polynom.c  	\
		dptest.c 	\
		reqseq.c    \
		main.c
OBJS = $(SRCS:.c=.o)
EXEC = main

.PHONY: build clean

$(EXEC): $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) $(SRCS) -o $(EXEC) $(LIBS)

run: $(EXEC)
	./$(EXEC)

clean:
	rm -f $(EXEC)
