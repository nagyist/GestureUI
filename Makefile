#
#	wiiuse Makefile
#

#
# Change this to your GCC version.
#
CC = g++-4.4
debug=0

####################################################
#
# You should not need to edit below this line.
#
####################################################

#
# Universal cflags
#
#CFLAGS = -Wall -pipe -fPIC -funroll-loops
CFLAGS = -Wall -pipe -fPIC -funroll-loops

ifeq ($(debug),1)
OBJ_PREFIX = debug
CFLAGS += -g -pg -DWITH_WIIUSE_DEBUG
else
OBJ_PREFIX = release
CFLAGS += -O2
endif

#OBJ_DIR = $(OBJ_PREFIX)-$(shell $(CC) -v 2>&1|grep ^Target:|cut -d' ' -f2)
OBJ_DIR = .

#
# Linking flags
#
#LDFLAGS = -L../src/$(OBJ_DIR) -lm -lwiiuse
LDFLAGS = -L./lib -L$(OBJ_DIR) -lm -lwiiuse

#
# Target binaries (always created as BIN)
#
BIN = ./gesture3d

#
# Inclusion paths.
#
INCLUDES = -I./libwiiuse -I./seqdb -I./math

#
# Generate a list of object files
#
OBJS = $(OBJ_DIR)/gesture3d.o ./seqdb/database.o

###############################
#
# Build targets.
#
###############################

all: $(BIN)

clean:
	@-rm $(OBJS) 2> /dev/null

distclean:
	@-rm -r debug-* release-* 2> /dev/null

install:
	@if [ -e $(BIN) ]; then \
	cp -v $(BIN) /usr/bin ; \
	fi


$(BIN): $(OBJS) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) -o $@

$(OBJ_DIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

mkdir:
	@if [ ! -d $(OBJ_DIR) ]; then \
	mkdir $(OBJ_DIR); \
	fi

run: all
	LD_LIBRARY_PATH=`pwd`/../src/$(OBJ_DIR):./seqdb:$(LD_LIBRARY_PATH) $(BIN)

