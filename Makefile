WCSLIBDIR := wcslib-$(WCSLIB_VERSION)

CC := gcc
CFLAGS := -Wall -Wextra -O2 -fopenmp -I$(WCSLIBDIR)/C -I$(WCSLIBDIR)
LDFLAGS := -lm -fopenmp
WCSLIB := $(WCSLIBDIR)/C/libwcs-$(WCSLIB_VERSION).a

TARGET := test_wcs_threads
SRC := test_wcs_threads.c

all: $(TARGET) $(TARGET)-bypass

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(LDFLAGS) $< $(WCSLIB) -o $@

$(TARGET)-bypass: $(SRC)
	$(CC) -DUSE_FLAG_TO_BYPASS $(CFLAGS) $(LDFLAGS) $< $(WCSLIB) -o $@

clean:
	rm -f $(TARGET) $(TARGET)-bypass
