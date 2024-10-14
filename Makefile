#
# FLAGS
#

PROGRAM_NAME = correct-barcodes

# PATHS
SRC_DIR = src
BUILD_DIR = build

SRC_CPP := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_CPP))

CPPFLAGS = -pedantic -Wall -Wextra -Wwrite-strings -Wstrict-overflow -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong -fno-strict-aliasing -fPIC $(INCLUDES)
CXXFLAGS = -std=c++11 -O3 -g $(CPPFLAGS)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CXXFLAGS += -D LINUX $(INCLUDES) -O3 -g
endif
ifeq ($(UNAME_S),Darwin)  # Mac with Homebrew
	CXXFLAGS += -D OSX -O3
endif

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
	CXXFLAGS += -D AMD64
endif
ifneq ($(filter %86,$(UNAME_P)),)
	CXXFLAGS += -D IA32
endif
ifneq ($(filter arm%,$(UNAME_P)),)
	CXXFLAGS += -D ARM
endif

ifdef BOOST_TAGGED
	BOOST_LIBS = -lboost_filesystem-mt -lboost_iostreams-mt -lboost_system-mt -lboost_chrono-mt -lz
else
	BOOST_LIBS = -lboost_filesystem -lboost_iostreams -lboost_system -lboost_chrono -lz
endif

ifdef BOOST_ROOT
	CPPFLAGS += -I$(BOOST_ROOT)/include
	LDFLAGS += -L$(BOOST_ROOT)/lib
else
	ifdef BOOST_INCLUDE
		CPPFLAGS += -I$(BOOST_INCLUDE)
	endif

	ifdef BOOST_LIB
		LDFLAGS += -L$(BOOST_LIB)
	endif
endif

LDLIBS = $(BOOST_LIBS)

PREFIX = /usr/local

#
# TARGETS
#

.PHONY: all checkdirs clean install

all: $(BUILD_DIR)/$(PROGRAM_NAME)

checkdirs: $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

install: $(BUILD_DIR)/$(PROGRAM_NAME)
	install -d -m 0755 $(PREFIX)/bin
	install -m 0755 $(BUILD_DIR)/$(PROGRAM_NAME) $(PREFIX)/bin/$(PROGRAM_NAME)

$(BUILD_DIR):
	@mkdir -p $@

$(BUILD_DIR)/$(PROGRAM_NAME): checkdirs $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<