WCSLIBDIR := wcslib-$(WCSLIB_VERSION)

CC := gcc
CFLAGS ?= -Wall -Wextra -O2 -fopenmp
LDFLAGS ?= -lm -fopenmp
MEMFLAG := -fsanitize=undefined -fsanitize=bounds -fsanitize=address -fsanitize-undefined-trap-on-error -fstack-protector-all
INCLUDE := -I$(WCSLIBDIR)/C -I$(WCSLIBDIR)
WCSLIB := $(WCSLIBDIR)/C/libwcs-$(WCSLIB_VERSION).a

ifneq ($(OS),Windows_NT)
  CFLAGS += $(MEMFLAG)
  LDFLAGS += $(MEMFLAG)
endif

TARGET := test_wcs_threads
SRC := test_wcs_threads.c
OBJS := $(SRC:.c=.o)

all: $(TARGET) $(TARGET)-bypass

$(TARGET): $(OBJS) $(WCSLIB)
	$(CC) $^ $(LDFLAGS) -o $@

%.o: %.c Makefile
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

$(TARGET)-bypass: $(SRC)
	$(CC) -DUSE_FLAG_TO_BYPASS $(INCLUDE) $(CFLAGS) $< $(WCSLIB) $(LDFLAGS) -o $@

.PHONY: clean
clean:
	rm -f $(TARGET) $(TARGET)-bypass
