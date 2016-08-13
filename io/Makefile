ROOT_DIR := ..
UTILS_DIR := $(ROOT_DIR)/utils

include $(ROOT_DIR)/common.mk
SRC   := ftread.c io.c
OBJS  := $(SRC:.c=.o)
INCL  := ftread.h io.h 

all: $(OBJS) $(INCL) Makefile 

%.o: %.c $(INCL) Makefile $(ROOT_DIR)/common.mk 
	$(CC) $(CFLAGS) -I$(UTILS_DIR) -c $< -o $@
clean:
	$(RM) $(OBJS)

clena: clean
celan: celan


