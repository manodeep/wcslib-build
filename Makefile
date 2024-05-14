WCSLIBDIR := wcslib-$(WCSLIB_VERSION)

CC := gcc
CFLAGS += -Wall -Wextra -O2 -fopenmp -I$(WCSLIBDIR)/C -I$(WCSLIBDIR)
LDFLAGS += -lm -fopenmp
WCSLIB := $(WCSLIBDIR)/C/libwcs-$(WCSLIB_VERSION).a

TARGET := test_wcs_threads
SRC := test_wcs_threads.c

all: $(TARGET) $(TARGET)-bypass

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -fsanitize=undefined -fsanitize=bounds -fsanitize=address -fsanitize-undefined-trap-on-error -fstack-protector-all $< $(WCSLIB) $(LDFLAGS) -o $@

$(TARGET)-bypass: $(SRC)
	$(CC) -DUSE_FLAG_TO_BYPASS -fsanitize=undefined -fsanitize=bounds -fsanitize=address -fsanitize-undefined-trap-on-error -fstack-protector-all $(CFLAGS) $< $(WCSLIB)  $(LDFLAGS) -o $@

clean:
	rm -f $(TARGET) $(TARGET)-bypass
