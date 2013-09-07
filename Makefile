TARGET 	= heh1p
SOURCES = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%.o,$(SOURCES))

RM 	:= rm
CXX := g++
CC 	:= gcc

CXX_DEBUG_FLAGS		=	-g -O0 -DDEBUG
CXX_RELEASE_FLAGS	=	-O2
UNAME := $(shell uname -a)
#target
ifeq ($(firstword $(filter Linux,$(UNAME))),Linux)
# for Linux
CXX_RELEASE_FLAGS += -march=native
endif

ifeq ($(firstword $(filter Darwin,$(UNAME))),Darwin)
# for MacOSX
CXX_RELEASE_FLAGS += -march=nocona
endif
CPPFLAGS = -Wall -std=c++0x -Wno-unused-but-set-variable
LDFLAGS  = -lm
VPATH	 = 
.DEFAULT_GOAL := release

# debug
.PHONY	: debug
debug 	: CXXFLAGS+=$(CXX_DEBUG_FLAGS)
debug 	: all
# release
.PHONY	: release
release	: CXXFLAGS+=$(CXX_RELEASE_FLAGS)
release	: all

all : $(TARGET)

${TARGET}:${OBJS}
	${CXX} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}
%.o: %.c
	${CXX} -c ${SOURCES} ${CXXFLAGS} 
# make clean
.PHONY: clean
clean:
	@for FILE in ${OBJS} ; do \
		if [ -f $$FILE ] ; then \
			echo "rm $$FILE." ; \
			${RM} -f $$FILE ; \
		fi ; \
	done
	@for FILE in ${TARGET} ; do \
		if [ -f $$FILE ] ; then \
			echo "rm $$FILE." ; \
			${RM} -f $$FILE ; \
		fi ; \
	done
