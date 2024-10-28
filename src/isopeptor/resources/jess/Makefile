# File: src/isopeptor/resources/jess/Makefile
CC=gcc
SRCS=$(wildcard *.c)
OBJS=$(SRCS:.c=.o)
TARGET=jess

all: $(TARGET)

# First compile all .c files
$(OBJS): $(SRCS)
	$(CC) -c $(SRCS)

# Then link all .o files with math library
$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) -lm

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean