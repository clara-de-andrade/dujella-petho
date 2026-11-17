CC        := gcc
AR        := ar
REL_FLAGS := -Wall -Wextra -O2
DBG_FLAGS := -Wall -Wextra -g -O0 -fsanitize=address,undefined
CFLAGS    ?= $(REL_FLAGS)

POLY_DIR  := polynomials
POLY_INC  := $(POLY_DIR)/include
POLY_LIB  := $(POLY_DIR)/libpolynomials.a

RSEQ_DIR  := sequences
RSEQ_INC  := $(RSEQ_DIR)/include
RSEQ_LIB  := $(RSEQ_DIR)/libsequences.a

LDLIBS    := -lmpfr -lgmp
LIBS      := -L$(POLY_DIR) -L$(RSEQ_DIR) -lpolynomials -lsequences
INCLUDES  := -I./ -I$(POLY_INC) -I$(RSEQ_INC)

SOURCES   := dujella_petho.c
HEADERS   := dujella_petho.h

TEST_DIR  := tests
TEST_BLD  := $(TEST_DIR)/build
TEST_SRCS := $(wildcard $(TEST_DIR)/*.c)
TEST_BINS := $(patsubst $(TEST_DIR)/%.c,$(TEST_BLD)/%,$(TEST_SRCS))

.PHONY: all test debug clean $(POLY_LIB) $(RSEQ_LIB)

all: $(POLY_LIB) $(RSEQ_LIB) test

# Release build and run all tests
test: $(TEST_BINS)

# Clean rebuild of libraries and tests with debug symbols and sanitizers
debug:
	$(MAKE) clean
	$(MAKE) test CFLAGS="$(DBG_FLAGS)"

$(TARGET): $(SOURCES) $(HEADERS) $(POLY_LIB) $(RSEQ_LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $(SOURCES) $(LIBS) $(LDLIBS) -o $@

# Compile each test directly into tests/build/ with no intermediate .o or .d files
$(TEST_BLD)/%: $(TEST_DIR)/%.c dujella_petho.c $(HEADERS) $(POLY_LIB) $(RSEQ_LIB)
	@mkdir -p $(TEST_BLD)
	$(CC) $(CFLAGS) $(INCLUDES) $< dujella_petho.c $(LIBS) $(LDLIBS) -o $@

$(POLY_LIB):
	$(MAKE) -C $(POLY_DIR) CFLAGS="$(CFLAGS)"

$(RSEQ_LIB):
	$(MAKE) -C $(RSEQ_DIR) CFLAGS="$(CFLAGS)"

clean:
	rm -rf $(TARGET) $(TARGET).exe $(TEST_BLD)
	$(MAKE) clean -C $(POLY_DIR)
	$(MAKE) clean -C $(RSEQ_DIR)