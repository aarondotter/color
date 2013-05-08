COLOR_DIR=/home/dotter/science/color

FC=gfortran
#FWARN = -Wimplicit-interface -Wunderflow -Wunused-parameter
FWARN = 
#FDBG=-g
FDBG=
FFLAGS = -fPIC -pedantic -Wall -O2 -ffixed-line-length-132
COMPILE = $(FC) $(FWARN) $(FDBG) $(FFLAGS)

COLOR_TEST_DIR = $(COLOR_DIR)/test
COLOR_LIB_DIR = $(COLOR_DIR)/make
LOAD_COLOR= -L$(COLOR_LIB_DIR) -lcolor
INCLUDE_COLOR = -I$(COLOR_LIB_DIR)
TEST=color_test

$(TEST): color_test.f $(COLOR_LIB_DIR)/$(LIB)
	$(COMPILE) -o $(TEST) color_test.f $(LOAD_COLOR) $(INCLUDE_COLOR)

clean:
	rm -f $(TEST) *~ *.o *.mod *genmod.f90 $(LIB)
