CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99
LDFLAGS = -lm

SRC = src/main.c src/dist.c src/random.c src/combi.c src/bayes.c src/stats.c
OBJ = $(SRC:.c=.o)
TARGET = probcalc

.PHONY: all clean install

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(TARGET)
	@echo "Clean complete"

install:
	@echo "Installing probcalc..."
	cp $(TARGET) ~/.local/bin/
	@echo "Installed to ~/.local/bin"