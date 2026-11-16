CC      := gcc
AR      := ar
CFLAGS  := -Wall -Wextra -O2

TARGET  := libpolynomials.a
INC_DIR := include
SRC_DIR := source
OBJ_DIR := build

SRCS    := $(wildcard $(SRC_DIR)/*.c)
OBJS    := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRCS))
DEPS    := $(OBJS:.o=.d)

.PHONY: all clean

all: $(TARGET)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)


$(TARGET): $(OBJS)
	$(AR) rcs $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -MMD -MP -c $< -o $@

-include $(DEPS)