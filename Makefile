CXX=icpc

# Flags passed to the C++ compiler.
CXXFLAGS =-m64 -fast -xCORE-AVX-I -unroll -Wall -std=c++11 -std=c++0x -pthread

TESTS_DIR = tests
INCLUDE_DIR = ./include
SRC_DIR = src
LIB_DIR = lib
BIN_DIR = bin


MKLPATH   = $(MKLROOT)/lib/intel64
MKLINCLUDE = $(MKLROOT)/include

#dynamic linking
LAPACKLIB = -ldl -L$(MKLPATH) -I$(MKLINCLUDE) -mkl=sequential

all: susceptibility susceptibility_sc susceptibility_sc_square susceptibility_square rho_lieb rho_square test 

susceptibility: common.o susceptibility.o	
	xiar rcs $(LIB_DIR)/libsusceptibility.a $(SRC_DIR)/common.o $(SRC_DIR)/susceptibility.o	
	$(CXX) $(CXXFLAGS) -static $(SRC_DIR)/X_q_Lieb.cc -I$(INCLUDE_DIR) -L$(LIB_DIR) -lsusceptibility -o $(BIN_DIR)/susceptibility

susceptibility_sc: common.o susceptibility_sc.o	
	xiar rcs $(LIB_DIR)/libsusceptibility_sc.a $(SRC_DIR)/common.o $(SRC_DIR)/susceptibility_sc.o	
	$(CXX) $(CXXFLAGS) -static $(SRC_DIR)/X_q_Lieb_sc.cc -I$(INCLUDE_DIR) -L$(LIB_DIR) -lsusceptibility_sc -o $(BIN_DIR)/susceptibility_sc

susceptibility_square: common.o susceptibility_square.o	
	xiar rcs $(LIB_DIR)/libsusceptibility_square.a $(SRC_DIR)/common.o $(SRC_DIR)/susceptibility_square.o	
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/X_q_square.cc -I$(INCLUDE_DIR) -L$(LIB_DIR) -lsusceptibility_square -o $(BIN_DIR)/susceptibility_square

susceptibility_sc_square: common.o susceptibility_sc_square.o	
	xiar rcs $(LIB_DIR)/libsusceptibility_sc_square.a $(SRC_DIR)/common.o $(SRC_DIR)/susceptibility_sc_square.o	
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/X_q_square_sc.cc -I$(INCLUDE_DIR) -L$(LIB_DIR) -lsusceptibility_sc_square -o $(BIN_DIR)/susceptibility_sc_square

test: 	
	(cd $(TESTS_DIR);$(MAKE))

rho_lieb: rho_lieb.o common.o $(INCLUDE_DIR)/common.h	
	xiar rcs $(LIB_DIR)/librho.a $(SRC_DIR)/rho.o $(SRC_DIR)/common.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) $(SRC_DIR)/rho_lieb.cc -L$(LIB_DIR) -lrho -o $(BIN_DIR)/rho_lieb

rho_square: rho_square.o common.o $(INCLUDE_DIR)/common.h	
	xiar rcs $(LIB_DIR)/librho_square.a $(SRC_DIR)/rho_square.o $(SRC_DIR)/common.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) $(SRC_DIR)/rho_square_generate.cc -L$(LIB_DIR) -lrho_square -o $(BIN_DIR)/rho_square


mean_field_phase_diagram: mean_field_phase_diagram.o mean_field.o common.o lapack_wrapper.o
	xiar rcs $(LIB_DIR)/libmean_field.a $(SRC_DIR)/common.o $(SRC_DIR)/mean_field.o $(SRC_DIR)/mean_field_phase_diagram.o $(SRC_DIR)/lapack_wrapper.o
	$(CXX) $(CXXFLAGS) -static -I$(INCLUDE_DIR) $(SRC_DIR)/mean_field_phase_diagram.cc -L$(LIB_DIR) $(LAPACKLIB) -lmean_field -o $(BIN_DIR)/create_mean_field_phase_diagram

common.o: $(SRC_DIR)/common.cc $(INCLUDE_DIR)/common.h
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(SRC_DIR)/common.cc -o $(SRC_DIR)/common.o

lapack_wrapper.o: $(SRC_DIR)/lapack_wrapper.cc $(INCLUDE_DIR)/lapack_wrapper.h
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(SRC_DIR)/lapack_wrapper.cc -o $(SRC_DIR)/lapack_wrapper.o

susceptibility.o: common.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(SRC_DIR)/susceptibility.cc -o $(SRC_DIR)/susceptibility.o

rho_lieb.o: common.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(SRC_DIR)/rho.cc -o $(SRC_DIR)/rho.o

rho_square.o: common.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(SRC_DIR)/rho_square.cc -o $(SRC_DIR)/rho_square.o

susceptibility_sc.o: common.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(LAPACKLIB) $(SRC_DIR)/susceptibility_sc.cc -o $(SRC_DIR)/susceptibility_sc.o

susceptibility_sc_square.o: common.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(LAPACKLIB) $(SRC_DIR)/susceptibility_sc_square.cc -o $(SRC_DIR)/susceptibility_sc_square.o

susceptibility_square.o: common.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(LAPACKLIB) $(SRC_DIR)/susceptibility_square.cc -o $(SRC_DIR)/susceptibility_square.o

mean_field.o:
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(LAPACKLIB) $(SRC_DIR)/mean_field.cc -o $(SRC_DIR)/mean_field.o

mean_field_phase_diagram.o: $(SRC_DIR)/mean_field_phase_diagram.cc $(INCLUDE_DIR)/mean_field_phase_diagram.h mean_field.o common.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE_DIR) $(LAPACKLIB) $(SRC_DIR)/mean_field_phase_diagram.cc -o $(SRC_DIR)/mean_field_phase_diagram.o
