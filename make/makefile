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
LIB_OBJS = all_color.o phx_color.o kur_color.o

#$(TEST): color_test.f $(LIB)
#	$(FC) $(FFLAGS) -o $(TEST) color_test.f -L$(COLOR_LIB) -lcolor

$(LIB): $(LIB_OBJS)
	ar crs $(LIB) $(LIB_OBJS)

phx_color.o: $(COLOR_SRC_DIR)/phx_color.f
	$(COMPILE) -c $<

kur_color.o: $(COLOR_SRC_DIR)/kur_color.f
	$(COMPILE) -c $<

all_color.o: $(COLOR_SRC_DIR)/all_color.f kur_color.o phx_color.o
	$(COMPILE) -c $<

clean:
	rm -f *~ *.o *.mod *genmod.f90 $(LIB)