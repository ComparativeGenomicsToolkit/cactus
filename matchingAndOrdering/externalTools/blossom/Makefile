DIRS := . MinCost GEOM 

SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

#CFLAGS := -O5 -D_NDEBUG -ggdb
CFLAGS := -D_NDEBUG -ggdb
CC := c++
LIBS :=  -lrt 
INCLUDES := 
LIBDIR := 

# Mac OS X specific stuff
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
	CFLAGS += -DPM_TIMER_GETRUSAGE
	LIBS := 
endif

binPath = ../../bin

all: ${binPath}/blossom5 ${binPath}/blossom_randGraph.py

${binPath}/blossom_randGraph.py: blossom_randGraph.py
	cp blossom_randGraph.py ${binPath}
	chmod +x ${binPath}/blossom_randGraph.py

${binPath}/blossom5: ${OBJS}
	${CC} ${CFLAGS} ${LIBDIR} -o $@ ${OBJS} ${LIBS}
	chmod +x ${binPath}/blossom5

.cpp.o:
	$(CC) $(CFLAGS) ${INCLUDES} $< -c -o $@

clean:
	rm -f ${OBJS}
	rm -f ${binPath}/blossom5 ${binPath}/blossom_randGraph.py
