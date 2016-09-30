# simple make file
# September 2004
SOURCES=main.f derivs.f
PRODUCT=a.out




GFFLAGS=-O2

GF=gfortran-fsf-5


all: $(PRODUCT)

$(PRODUCT) : $(SOURCES)
	$(GF) $(GFFLAGS) -o $(PRODUCT) $(SOURCES)
