COLOR_DIR=/home/dotter/science/color


FC=gfortran
#FWARN = -Wimplicit-interface -Wunderflow -Wunused-parameter
FWARN = 
#FDBG=-g
FDBG=
FFLAGS = -fPIC -pedantic -Wall -O2 -ffixed-line-length-132
COMPILE = $(FC) $(FWARN) $(FDBG) $(FFLAGS)

COLOR_SRC_DIR = $(COLOR_DIR)/src

LIB  = libcolor.a
LIB_OBJS = color_lib.o phx_color.o kur_color.o color_def.o

$(LIB): $(LIB_OBJS)
	ar crs $(LIB) $(LIB_OBJS)

color_def.o: $(COLOR_SRC_DIR)/color_def.f
	$(COMPILE) -c $<

phx_color.o: $(COLOR_SRC_DIR)/phx_color.f color_def.o
	$(COMPILE) -c $<

kur_color.o: $(COLOR_SRC_DIR)/kur_color.f color_def.o
	$(COMPILE) -c $<

color_lib.o: $(COLOR_SRC_DIR)/color_lib.f kur_color.o phx_color.o color_def.o
	$(COMPILE) -c $<

clean:
	rm -f *~ *.o *.mod *genmod.f90 $(LIB)