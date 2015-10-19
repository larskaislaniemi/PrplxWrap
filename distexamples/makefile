CC=gcc
CFLAGS=-c -Wall -g -I.. -std=c99
#CFLAGS=-c -Wall -O3 -I.. -std=c99
LDFLAGS=-L.. -L$(PERPLEXSRC) -lperplex -lmeemum -lgfortran -g

SOURCES=gridding.c

OBJECTS=$(SOURCES:.c=.o)
WRAPPEROBJECTS=../perplex.o
PERPLEXSRC=../../perplex672
PERPLEXOBJECTS=$(PERPLEXSRC)/meemum.o

EXECUTABLE=gridding

all: $(SOURCES) $(EXECUTABLE) $(OBJECTS) $(WRAPPEROBJECTS) $(PERPLEXOBJECTS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

obj: $(SOURCES) $(OBJECTS) $(WRAPPEROBJECTS) $(PERPLEXOBJECTS)
	$(CC) $(CFLAGS) $(SOURCES)

$(EXECUTABLE): $(OBJECTS) $(WRAPPEROBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $(EXECUTABLE)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

