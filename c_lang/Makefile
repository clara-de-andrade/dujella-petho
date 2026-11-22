CC        := gcc
AR        := ar
REL_FLAGS := -Wall -Wextra -O2
DBG_FLAGS := -Wall -Wextra -g -O0 -fsanitize=address,undefined
CFLAGS    ?= $(REL_FLAGS)

POLY_DIR  := polynomials
POLY_INC  := $(POLY_DIR)/include
POLY_LIB  := $(POLY_DIR)/libpolynomials.a

LDLIBS    := -lmpfr -lgmp
LIBS      := -L$(POLY_DIR) -lpolynomials
INCLUDES  := -I./ -I$(POLY_INC)

SOURCES   := dujella_petho.c matveev.c
HEADERS   := dujella_petho.h matveev.h

TEST_DIR  := tests
TEST_BLD  := $(TEST_DIR)/build
TEST_SRCS := $(wildcard $(TEST_DIR)/*.c)
TEST_BINS := $(patsubst $(TEST_DIR)/%.c,$(TEST_BLD)/%,$(TEST_SRCS))

.PHONY: all test debug clean $(POLY_LIB) $(RSEQ_LIB)

all: $(POLY_LIB) $(RSEQ_LIB) $(TEST_BINS)

debug:
	$(MAKE) clean
	$(MAKE) test CFLAGS="$(DBG_FLAGS)"

$(TARGET): $(SOURCES) $(HEADERS) $(POLY_LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $(SOURCES) $(LIBS) $(LDLIBS) -o $@

$(TEST_BLD)/%: $(TEST_DIR)/%.c $(SOURCES) $(HEADERS) $(POLY_LIB)
	@mkdir -p $(TEST_BLD)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(SOURCES) $(LIBS) $(LDLIBS) -o $@

$(POLY_LIB):
	$(MAKE) -C $(POLY_DIR) CFLAGS="$(CFLAGS)"

clean:
	rm -rf $(TARGET) $(TARGET).exe $(TEST_BLD)
	$(MAKE) clean -C $(POLY_DIR)
